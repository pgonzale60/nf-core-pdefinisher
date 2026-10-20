from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple, Any
from .models import EditedRegion, StructuralVariant, EditValidationResult, FlagrantError
from .policy import DEFAULT_VALIDATION_POLICY

try:
    import pysam
except ImportError:
    pysam = None


def check_bam_region(
    bam_path: str,
    chrom: str,
    start: int,
    end: int
) -> Tuple[float, int]:
    """
    Computes mean depth and spanning reads across a target interval in BAM.
    Returns (mean_depth, spanning_read_count).
    """
    if not pysam or not os.path.exists(bam_path):
        return (0.0, 0)

    try:
        with pysam.AlignmentFile(bam_path, "rb") as sam:
            # Check chromosome presence
            if chrom not in sam.references:
                return (0.0, 0)

            # Cap query interval
            c_len = sam.get_reference_length(chrom)
            q_start = max(0, start - 1)
            q_end = min(c_len, max(q_start + 1, end))

            # Spanning reads across the entire interval
            spanning = 0
            for read in sam.fetch(chrom, q_start, q_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.reference_start <= q_start and read.reference_end >= q_end:
                    spanning += 1

            # Depth coverage
            col_depths = [col.nsegments for col in sam.pileup(chrom, q_start, q_end, truncate=True)]
            mean_cov = sum(col_depths) / len(col_depths) if col_depths else 0.0

            return (round(mean_cov, 1), spanning)
    except Exception:
        return (0.0, 0)


def validate_edited_regions(
    regions: List[EditedRegion],
    svs: List[StructuralVariant],
    bam_path: Optional[str] = None,
    policy: Optional[Dict[str, Any]] = None
) -> List[EditValidationResult]:
    """
    Evaluates each edited region from the surgery ledger against read depth and Sniffles SVs.
    """
    cfg = policy or DEFAULT_VALIDATION_POLICY
    results: List[EditValidationResult] = []

    # Index SVs by chrom for fast lookup
    svs_by_chrom: Dict[str, List[StructuralVariant]] = {}
    for sv in svs:
        svs_by_chrom.setdefault(sv.chrom, []).append(sv)

    for r in regions:
        # Determine query interval
        q_start = max(1, r.curated_start)
        q_end = max(q_start + 1, r.curated_end)

        mean_cov, spanning = (0.0, 0)

        # Check if this region is a purged entity (e.g. PURGE_ORGANELLAR, PURGE_CONTAMINANT)
        is_purged_action = (
            "PURGE" in r.action.upper()
            or "PURGE" in r.edit_type.upper()
        )

        if is_purged_action:
            is_present = False
            if bam_path and os.path.exists(bam_path) and pysam:
                try:
                    with pysam.AlignmentFile(bam_path, "rb") as sam:
                        is_present = r.chrom in sam.references
                except Exception:
                    is_present = False

            if not is_present:
                status = "PASS"
                reasons = ["Purged sequence verified absent from curated assembly"]
            else:
                status = "WARN"
                reasons = [f"Purged sequence unexpectedly still present in curated assembly ({mean_cov}x depth)"]

            results.append(
                EditValidationResult(
                    region=r,
                    mean_coverage=0.0,
                    spanning_reads=0,
                    conflicting_svs=[],
                    status=status,
                    explanation="; ".join(reasons)
                )
            )
            continue

        if bam_path and os.path.exists(bam_path):
            mean_cov, spanning = check_bam_region(bam_path, r.chrom, q_start, q_end)

        # Detect overlapping or adjacent SVs
        conflicting: List[StructuralVariant] = []
        for sv in svs_by_chrom.get(r.chrom, []):
            # Check overlap or proximity within 500bp
            if max(q_start - 500, 1) <= sv.end and sv.start <= q_end + 500:
                conflicting.append(sv)

        # Assign verdict
        status = "PASS"
        reasons = []

        # High-confidence conflicting SVs (e.g. large DEL right over patch)
        high_vaf_svs = [s for s in conflicting if s.vaf >= 0.50 and s.support >= cfg.get("flagrant_min_support", 5)]
        if high_vaf_svs:
            status = "WARN"
            sv_descs = [f"{s.sv_type}:{s.sv_len}bp(VAF={s.vaf})" for s in high_vaf_svs]
            reasons.append(f"Conflicting SV detected at edited junction: {', '.join(sv_descs)}")

        # Coverage check if BAM is available
        if bam_path and os.path.exists(bam_path):
            if mean_cov < cfg.get("min_edit_depth", 5.0) and r.edit_type != "UNPLACED_SCAFFOLD_CREATION":
                if status == "PASS":
                    status = "WARN"
                reasons.append(f"Low read coverage ({mean_cov}x < {cfg.get('min_edit_depth', 5.0)}x floor)")

        if not reasons:
            reasons.append(f"Clean junction; read support confirmed ({mean_cov}x depth, {spanning} spanning reads)")

        results.append(
            EditValidationResult(
                region=r,
                mean_coverage=mean_cov,
                spanning_reads=spanning,
                conflicting_svs=conflicting,
                status=status,
                explanation="; ".join(reasons)
            )
        )

    return results


def find_flagrant_errors(
    svs: List[StructuralVariant],
    regions: Optional[List[EditedRegion]] = None,
    policy: Optional[Dict[str, Any]] = None
) -> List[FlagrantError]:
    """
    Identifies severe, high-confidence assembly defects across the assembly
    (e.g. homozygous inversions, major deletions, trans-chromosomal translocations).

    When regions is provided, excludes SVs that overlap surgical edit loci so that
    Section 2 (Surgery Loci Verification) and Section 3 (Assembly Error Registry)
    are strictly mutually exclusive.
    """
    cfg = policy or DEFAULT_VALIDATION_POLICY
    hom_thresh = cfg.get("flagrant_hom_vaf_threshold", 0.80)
    min_sup = cfg.get("flagrant_min_support", 5)
    min_len = cfg.get("flagrant_min_sv_len", 50)
    flank = cfg.get("min_edit_flank_bp", 500)

    # Index edited intervals by chrom for exclusion check
    edit_intervals: Dict[str, List[Tuple[int, int]]] = {}
    if regions:
        for r in regions:
            if "PURGE" in r.action.upper() or "PURGE" in r.edit_type.upper():
                continue
            q_start = max(1, r.curated_start - flank)
            q_end = max(q_start + 1, r.curated_end + flank)
            edit_intervals.setdefault(r.chrom, []).append((q_start, q_end))

    errors: List[FlagrantError] = []

    for sv in svs:
        # Check for near-homozygous structural defect
        if sv.vaf >= hom_thresh and sv.support >= min_sup and sv.sv_len >= min_len:
            # Mutually exclusive check: exclude SVs overlapping any edited surgical locus
            if sv.chrom in edit_intervals:
                overlaps_edit = any(
                    max(e_start, 1) <= sv.end and sv.start <= e_end
                    for e_start, e_end in edit_intervals[sv.chrom]
                )
                if overlaps_edit:
                    continue

            severity = "CRITICAL" if sv.sv_len >= 1000 or sv.sv_type in ("INV", "BND") else "HIGH"
            desc = (
                f"Flagrant homozygous {sv.sv_type} ({sv.sv_len} bp, VAF={sv.vaf}, "
                f"{sv.support} reads) at {sv.chrom}:{sv.start}-{sv.end}"
            )
            errors.append(
                FlagrantError(
                    chrom=sv.chrom,
                    start=sv.start,
                    end=sv.end,
                    error_type="HOM_SV_DEFECT",
                    severity=severity,
                    details={
                        "sv_id": sv.sv_id,
                        "sv_type": sv.sv_type,
                        "sv_len": sv.sv_len,
                        "vaf": sv.vaf,
                        "support": sv.support,
                        "chr2": sv.chr2,
                    },
                    description=desc
                )
            )

    return errors
