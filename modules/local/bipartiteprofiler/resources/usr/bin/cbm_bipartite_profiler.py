#!/usr/bin/env python3
"""Stage 5: Soft-Erasure Multi-Motif Discovery, SpaMo Spacing Validation & Bipartite Profiler."""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from bipartite.aux_grid import (  # noqa: E402
    is_core_echo,
    run_localized_auxiliary_meme,
)

__all__ = ["build_higher_order_markov_bg", "is_core_echo", "main", "parse_args"]
from bipartite.masker import (  # noqa: E402
    load_fasta_dict,
    load_soma_len_map,
    mask_core_motifs,
    run_fimo,
    write_sliced_fasta,
)
from bipartite.reporter import generate_cbm_discovery_report  # noqa: E402
from bipartite.spacing import (  # noqa: E402
    detect_homotypic_bipartite,
    profile_bipartite_spacing,
    write_spacing_tsv,
)
from cbm_defaults import SSOT  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Soft-erasure multi-motif discovery, SpaMo spacing validation, and bipartite profiler"
    )
    parser.add_argument("--break-windows", required=True, help="Path to break_windows_300bp.fasta")
    parser.add_argument("--core-meme", required=True, help="Path to best_core_motif.meme")
    parser.add_argument("--census-json", required=True, help="Path to coverage_and_chromosome_census.json")
    parser.add_argument("--core-metrics-json", required=True, help="Path to core_motif_metrics.json")
    parser.add_argument("--background-bfile", default=None, help="Winning Markov background from Stage 4")
    parser.add_argument("--sample-id", default="sample", help="Sample ID")
    parser.add_argument("--output-aux-meme", default="best_auxiliary_motif.meme", help="Output aux MEME")
    parser.add_argument("--output-spacing-tsv", default="bipartite_spacing_summary.tsv", help="Spacing TSV")
    parser.add_argument("--output-report", default="cbm_discovery_report.md", help="Output markdown report")
    parser.add_argument("--output-masked-fasta", default="break_sequences.core_masked.fasta", help="Masked FASTA")
    parser.add_argument("--threads", type=int, default=1, help="Threads for MEME/FIMO/SpaMo")
    parser.add_argument("--telomere-repeat", default=SSOT["telomere_repeat"], help="Canonical telomere repeat for aux echo filter")
    return parser.parse_args()


def build_higher_order_markov_bg(fasta_file: str | Path, out_bg: Path) -> Path:
    """Build higher-order Markov background from distal flanking arms outside [-excl_r..+excl_r]."""
    fgm_bin = shutil.which("fasta-get-markov")
    if not fgm_bin:
        raise RuntimeError("Required binary 'fasta-get-markov' not found in PATH")

    seqs = load_fasta_dict(fasta_file)
    soma_lens = load_soma_len_map(fasta_file)
    distal_fa = out_bg.parent / "distal_arms_bg.fasta"
    excl_r = SSOT["distal_bg_exclusion_radius"]

    total_len = 0
    with open(distal_fa, "w", encoding="utf-8") as f_out:
        for idx, (hdr, s) in enumerate(seqs.items()):
            s = s.upper()
            pos0 = soma_lens.get(hdr, len(s) // 2)
            left_arm = s[: max(0, pos0 - excl_r)]
            right_arm = s[pos0 + excl_r :] if len(s) > pos0 + excl_r else ""
            if left_arm:
                f_out.write(f">distal_L_{idx}\n{left_arm}\n")
                total_len += len(left_arm)
            if right_arm:
                f_out.write(f">distal_R_{idx}\n{right_arm}\n")
                total_len += len(right_arm)

    order = 3 if total_len >= 3600 else (2 if total_len >= 1440 else 1)
    subprocess.run(
        [fgm_bin, "-dna", "-m", str(order), str(distal_fa), str(out_bg)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    return out_bg


def main():
    args = parse_args()
    # Check required binaries fail loudly
    for b in ["meme", "fimo"]:
        if not shutil.which(b):
            raise RuntimeError(f"Required binary '{b}' not found in PATH")

    raw_seqs = load_fasta_dict(args.break_windows)
    soma_lens = load_soma_len_map(args.break_windows)
    seq_names = list(raw_seqs.keys())

    work_dir = Path("bipartite_work")
    work_dir.mkdir(exist_ok=True)

    # Handle negative control / empty break windows / single sequence (< 2 sequences)
    if len(seq_names) < 2:
        with open(args.output_aux_meme, "w", encoding="utf-8") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n")
        with open(args.output_masked_fasta, "w", encoding="utf-8") as f:
            pass
        spacing_rows, spacing_summary = profile_bipartite_spacing(
            seq_names=seq_names,
            core_hits_by_seq={},
            aux_hits_by_seq=None,
            pos0_coord={},
            architecture="SIMPLE",
        )
        write_spacing_tsv(spacing_rows, args.output_spacing_tsv)
        generate_cbm_discovery_report(
            sample_id=args.sample_id,
            census_json_path=args.census_json,
            core_metrics_json_path=args.core_metrics_json,
            best_aux_motif=None,
            bipartite_summary=spacing_summary,
            output_report_md=args.output_report,
        )
        print(
            f"[BIPARTITEPROFILER] Completed for {args.sample_id} "
            f"({len(seq_names)} break window(s) found - insufficient sequences for statistical motif discovery)."
        )
        return

    # 1. Parse core metrics JSON to retrieve winning radius and parameters
    core_metrics = {}
    if args.core_metrics_json and os.path.exists(args.core_metrics_json):
        try:
            with open(args.core_metrics_json, encoding="utf-8") as f_cm:
                core_metrics = json.load(f_cm)
        except Exception:
            pass

    best_radius = (
        core_metrics.get("best_radius")
        or core_metrics.get("radius")
        or 50
    )
    try:
        best_radius = int(best_radius)
    except (ValueError, TypeError):
        best_radius = 50

    best_mod = (
        core_metrics.get("best_model")
        or core_metrics.get("mod")
        or "zoops"
    )
    if best_mod not in ("zoops", "anr", "oops"):
        best_mod = "zoops"
    core_consensus = core_metrics.get("consensus", "")

    # 2. Slice input FASTA to optimal spatial radius R*
    sliced_fa = work_dir / f"windows_radius_{best_radius}bp.fasta"
    sliced_seqs, sliced_soma_lens = write_sliced_fasta(
        raw_seqs=raw_seqs,
        soma_lens=soma_lens,
        radius=best_radius,
        output_fasta=sliced_fa,
    )

    # 3. Core motif instance mapping via FIMO and core-masking
    _masked_seqs, core_hits = mask_core_motifs(
        core_meme=args.core_meme,
        fasta_file=sliced_fa,
        output_masked_fasta=args.output_masked_fasta,
        fimo_out_dir=work_dir / "fimo_core_hits",
        threshold=SSOT["fimo_core_threshold"],
    )

    core_hits_by_seq = defaultdict(list)
    for h in core_hits:
        core_hits_by_seq[h["seq_name"]].append(h)

    # 4. Step 1: Homotypic Core Partner Scan
    is_homotypic, homotypic_info = detect_homotypic_bipartite(
        core_hits=core_hits,
        seq_names=seq_names,
        min_fraction=0.30,
        max_spacer_iqr=25.0,
    )

    best_aux = None
    if is_homotypic:
        architecture = "HOMOTYPIC_BIPARTITE"
        best_aux = {
            "consensus": core_consensus,
            "width": core_metrics.get("width", len(core_consensus)),
            "evalue": core_metrics.get("evalue", 0.0),
            "evalue_str": str(core_metrics.get("evalue_str", core_metrics.get("evalue", "N/A"))),
            "mean_ic": core_metrics.get("mean_ic", 0.0),
            "architecture": "HOMOTYPIC_BIPARTITE",
            "partner_type": "Secondary_Core",
            "median_spacer_bp": homotypic_info.get("median_spacer", 0.0),
            "iqr_spacer_bp": homotypic_info.get("iqr_spacer", 0.0),
        }
        # In homotypic bipartite, partner motif is the core motif itself
        shutil.copy(args.core_meme, args.output_aux_meme)

        spacing_rows, spacing_summary = profile_bipartite_spacing(
            seq_names=seq_names,
            core_hits_by_seq=core_hits_by_seq,
            aux_hits_by_seq=None,
            pos0_coord=sliced_soma_lens,
            architecture=architecture,
        )
        spacing_summary["architecture"] = architecture
    else:
        # 5. Step 2: Localized Multi-Motif Re-run on sliced window
        if args.background_bfile and os.path.exists(args.background_bfile):
            bg_file = Path(args.background_bfile)
        else:
            bg_file = build_higher_order_markov_bg(
                args.break_windows, work_dir / "bg_bipartite.txt"
            )

        meme_aux_dir = work_dir / "meme_aux_out"
        best_candidate, _aux_candidates = run_localized_auxiliary_meme(
            sliced_fasta=sliced_fa,
            out_dir=meme_aux_dir,
            bfile=str(bg_file) if bg_file else None,
            nmotifs=3,
            minw=6,
            maxw=30,
            mod=best_mod,
            revcomp=True,
            core_consensus=core_consensus,
            total_sequences=len(seq_names),
            telomere_repeat=args.telomere_repeat,
        )

        if best_candidate and best_candidate.get("meme_txt_path") and os.path.exists(best_candidate["meme_txt_path"]):
            architecture = "HETEROTYPIC_BIPARTITE"
            best_aux = best_candidate
            best_aux["architecture"] = architecture
            shutil.copy(best_candidate["meme_txt_path"], args.output_aux_meme)

            aux_fimo_hits = run_fimo(
                meme_file=args.output_aux_meme,
                fasta_file=sliced_fa,
                out_dir=work_dir / "fimo_aux_hits",
                threshold=SSOT["fimo_aux_threshold"],
            )
            aux_hits_by_seq = defaultdict(list)
            for h in aux_fimo_hits:
                aux_hits_by_seq[h["seq_name"]].append(h)

            spacing_rows, spacing_summary = profile_bipartite_spacing(
                seq_names=seq_names,
                core_hits_by_seq=core_hits_by_seq,
                aux_hits_by_seq=aux_hits_by_seq,
                pos0_coord=sliced_soma_lens,
                architecture=architecture,
            )
            spacing_summary["architecture"] = architecture
        else:
            architecture = "SIMPLE"
            with open(args.output_aux_meme, "w", encoding="utf-8") as f:
                f.write("MEME version 4\n\nALPHABET= ACGT\n")

            spacing_rows, spacing_summary = profile_bipartite_spacing(
                seq_names=seq_names,
                core_hits_by_seq=core_hits_by_seq,
                aux_hits_by_seq=None,
                pos0_coord=sliced_soma_lens,
                architecture=architecture,
            )
            spacing_summary["architecture"] = architecture

    # 6. Write spacing TSV and Markdown report
    write_spacing_tsv(spacing_rows, args.output_spacing_tsv)

    generate_cbm_discovery_report(
        sample_id=args.sample_id,
        census_json_path=args.census_json,
        core_metrics_json_path=args.core_metrics_json,
        best_aux_motif=best_aux,
        bipartite_summary=spacing_summary,
        output_report_md=args.output_report,
    )

    partner_desc = (
        f"Homotypic Core ({best_aux['consensus']})"
        if architecture == "HOMOTYPIC_BIPARTITE"
        else (f"Auxiliary ({best_aux['consensus']})" if best_aux else "None")
    )
    print(
        f"[BIPARTITEPROFILER] Completed for {args.sample_id}. "
        f"Architecture: {architecture}. Partner motif: {partner_desc}. "
        f"Paired loci: {spacing_summary.get('n_paired_occurrences', 0)}."
    )


if __name__ == "__main__":
    main()
