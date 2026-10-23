#!/usr/bin/env python3
"""
pde_curation_summary.py - Cohort-Level Curation Synthesis & Executive Scorecard

Aggregates individual sample curation scorecards (*.curation_scorecard.json)
produced by CURATIONEVALUATOR across an entire cohort or batch into a single,
comprehensive summary.

Outputs:
  1. ${prefix}.curation_summary.md  - Executive Markdown report with KPI cards,
                                     tier classifications, league table, and triage.
  2. ${prefix}.curation_summary.tsv - Flat tabular matrix of all cohort metrics.
  3. ${prefix}.curation_summary.json- Complete structured JSON payload.
"""

import argparse
import hashlib
import json
import os
import sys
from collections import defaultdict

__version__ = "1.3.0"


def get_script_fingerprint():
    try:
        with open(__file__, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except Exception:
        return "unknown"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cohort-Level Post-Surgery Curation Summary and Scorecard Aggregator"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__} (sha256:{get_script_fingerprint()})",
    )
    parser.add_argument(
        "--scorecards",
        nargs="+",
        required=True,
        help="List of sample curation scorecard JSON files (*.curation_scorecard.json)",
    )
    parser.add_argument(
        "--validations",
        nargs="*",
        default=[],
        help="Optional list of sample validation JSON files (*.curation_validation.json)",
    )
    parser.add_argument(
        "--prefix",
        default="cohort",
        help="Prefix for output summary files (default: cohort)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: current directory)",
    )
    return parser.parse_args()


def load_scorecard(filepath):
    """Load and validate a single sample curation scorecard JSON file."""
    if not os.path.exists(filepath):
        print(f"WARNING: Scorecard file not found: {filepath}", file=sys.stderr)
        return None
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        print(f"WARNING: Failed to parse scorecard {filepath}: {exc}", file=sys.stderr)
        return None


def classify_tier(sample_data):
    """
    Classify sample into an executive curation status tier:
      - Tier 1: Gold Standard (100% T2T AND Non-rDNA Gapless, permitting 100bp rDNA sizing gap)
      - Tier 2: Silver Standard (100% T2T with gaps OR Non-rDNA Gapless near-T2T >=85%)
      - Tier 3: Bronze Standard (Substantially Improved)
      - Tier 4: Minimal / Attention Needed
    """
    is_t2t = sample_data["is_t2t"]
    is_gapless = sample_data["is_gapless_non_rdna"]
    post_t2t_pct = sample_data["post_t2t_pct"]
    pre_t2t_pct = sample_data["pre_t2t_pct"]
    gaps_resolved = sample_data["gaps_resolved"]

    if is_t2t and is_gapless:
        return "Tier 1: Gold (T2T Gapless)"
    elif is_t2t and not is_gapless:
        return "Tier 2: Silver (T2T Scaffolding)"
    elif not is_t2t and is_gapless and post_t2t_pct >= 85.0:
        return "Tier 2: Silver (Gapless Near-T2T)"
    elif post_t2t_pct > pre_t2t_pct or gaps_resolved > 0:
        return "Tier 3: Bronze (Substantially Improved)"
    else:
        return "Tier 4: Attention Needed"


def parse_sample_record(sc):
    """Extract and normalise metrics from a single sample scorecard."""
    sample = sc.get("sample", "unknown")
    contig = sc.get("contiguity", {})
    pre_c = contig.get("pre_curation", {})
    post_c = contig.get("post_curation", {})
    deltas = contig.get("deltas", {})

    t2t = sc.get("t2t_completeness", {})
    gap = sc.get("gap_resolution", {})
    trans = sc.get("transchromosomal_surgery", {})
    purge = sc.get("purge_ledger", {})
    recon = sc.get("reconciliation", {})

    total_termini = t2t.get("total_chromosome_termini", 0)
    chrom_count = total_termini // 2 if total_termini > 0 else 0
    pre_capped = t2t.get("pre_capped_termini", 0)
    post_capped = t2t.get("post_capped_termini", 0)
    pre_t2t_pct = t2t.get("pre_t2t_score_pct", 0.0)
    post_t2t_pct = t2t.get("post_t2t_score_pct", 0.0)
    is_t2t = (post_capped == total_termini) and (total_termini > 0)
    uncapped_ends = t2t.get("uncapped_ends", [])
    net_telo_bp = t2t.get("net_telomeric_bp_added", 0)
    term_breakdown = t2t.get("terminus_status_breakdown", {})

    baseline_gaps = gap.get("baseline_gaps", pre_c.get("total_gaps", 0))
    remaining_gaps = gap.get("remaining_gaps", post_c.get("total_gaps", 0))
    is_gapless = (remaining_gaps == 0)

    # Distinguish non-rDNA gaps from rDNA array sizing gaps (100 bp N placeholder)
    rdna_gaps = gap.get("rdna_gaps", 0)
    remaining_non_rdna = gap.get("remaining_non_rdna_gaps", 0)
    is_gapless_non_rdna = gap.get("is_gapless_non_rdna", (remaining_non_rdna == 0))

    # Size deltas: overall delta vs non-rDNA delta
    total_delta = deltas.get("length_delta_bp", recon.get("measured_delta_bp", 0))
    rdna_delta = deltas.get("rdna_length_delta_bp", recon.get("rdna_delta_bp", 0))
    non_rdna_delta = deltas.get("non_rdna_length_delta_bp", recon.get("non_rdna_delta_bp", total_delta - rdna_delta))

    gaps_excised = gap.get("gaps_excised_count", 0)
    gaps_patched = gap.get("gaps_patched_count", 0)
    gaps_resolved = max(0, baseline_gaps - remaining_gaps)
    gap_resolution_pct = (
        round((gaps_resolved / baseline_gaps) * 100.0, 1)
        if baseline_gaps > 0
        else (100.0 if remaining_gaps == 0 else 0.0)
    )

    pde = sc.get("karyotype_compartments", {})
    core_info = pde.get("putative_core", {})
    grs_info = pde.get("putative_grs", {})
    core_span_bp = core_info.get("span_bp", 0)
    core_median_cov = core_info.get("median_coverage", 0.0)
    grs_span_bp = grs_info.get("span_bp", 0)
    grs_median_cov = grs_info.get("median_coverage", 0.0)
    grs_thresh = pde.get("classification_threshold", 0.0)
    tot_pde_len = core_span_bp + grs_span_bp
    grs_frac_pct = round((grs_span_bp / tot_pde_len * 100.0), 2) if tot_pde_len > 0 else 0.0

    elim_range = pde.get("elimination_range", {})
    canon_grs_bp = elim_range.get("canon_grs_bp", grs_span_bp)
    unplaced_grs_bp = elim_range.get("unplaced_grs_bp", 0)
    unplaced_uncertain_bp = elim_range.get("unplaced_uncertain_bp", 0)
    unplaced_core_bp = elim_range.get("unplaced_core_bp", 0)
    min_elim_bp = elim_range.get("min_elim_bp", grs_span_bp)
    max_elim_bp = elim_range.get("max_elim_bp", grs_span_bp)
    est_elim_bp = elim_range.get("est_elim_bp", grs_span_bp)
    min_elim_pct = round((min_elim_bp / tot_pde_len * 100.0), 2) if tot_pde_len > 0 else 0.0
    max_elim_pct = round((max_elim_bp / tot_pde_len * 100.0), 2) if tot_pde_len > 0 else 0.0
    est_elim_pct = round((est_elim_bp / tot_pde_len * 100.0), 2) if tot_pde_len > 0 else 0.0

    data = {
        "sample": sample,
        "chromosomes_count": chrom_count,
        "total_termini": total_termini,
        "pre_capped": pre_capped,
        "post_capped": post_capped,
        "pre_t2t_pct": pre_t2t_pct,
        "post_t2t_pct": post_t2t_pct,
        "t2t_gain_pct": round(post_t2t_pct - pre_t2t_pct, 1),
        "is_t2t": is_t2t,
        "uncapped_ends": uncapped_ends,
        "net_telo_bp_added": net_telo_bp,
        "terminus_status_breakdown": term_breakdown,
        "baseline_gaps": baseline_gaps,
        "remaining_gaps": remaining_gaps,
        "rdna_gaps": rdna_gaps,
        "remaining_non_rdna_gaps": remaining_non_rdna,
        "gaps_excised": gaps_excised,
        "gaps_patched": gaps_patched,
        "gaps_resolved": gaps_resolved,
        "gap_resolution_pct": gap_resolution_pct,
        "is_gapless": is_gapless,
        "is_gapless_non_rdna": is_gapless_non_rdna,
        "patched_bp": gap.get("total_patched_bp", 0),
        "excised_bp": gap.get("total_excised_bp", 0),
        "pre_total_bp": pre_c.get("total_bp", 0),
        "post_total_bp": post_c.get("total_bp", 0),
        "delta_bp": total_delta,
        "rdna_delta_bp": rdna_delta,
        "non_rdna_delta_bp": non_rdna_delta,
        "pre_scaffolds": pre_c.get("scaffold_count", 0),
        "post_scaffolds": post_c.get("scaffold_count", 0),
        "scaffolds_delta": deltas.get("scaffolds_delta", 0),
        "pre_n50": pre_c.get("n50", 0),
        "post_n50": post_c.get("n50", 0),
        "pre_chr_n50": pre_c.get("chr_n50", 0),
        "post_chr_n50": post_c.get("chr_n50", 0),
        "pre_n_bp": pre_c.get("total_n_bp", 0),
        "post_n_bp": post_c.get("total_n_bp", 0),
        "transfers_count": trans.get("transfers_count", 0),
        "fissions_count": trans.get("fissions_count", 0),
        "transferred_bp": trans.get("transferred_bp", 0),
        "severed_bp": trans.get("severed_bp", 0),
        "purged_haplotigs_count": purge.get("purged_haplotigs_count", 0),
        "purged_haplotigs_bp": purge.get("purged_haplotigs_bp", 0),
        "purged_rdna_count": purge.get("purged_rdna_scaffolds_count", 0),
        "purged_rdna_bp": purge.get("purged_rdna_scaffolds_bp", 0),
        "primary_purity_pct": purge.get("primary_release_purity_pct", 0.0),
        "mathematically_conserved": recon.get("mathematically_conserved", True),
        "discrepancy_bp": recon.get("discrepancy_bp", 0),
        "core_span_bp": core_span_bp,
        "core_median_cov": core_median_cov,
        "grs_span_bp": grs_span_bp,
        "grs_median_cov": grs_median_cov,
        "grs_threshold": grs_thresh,
        "grs_fraction_pct": grs_frac_pct,
        "canon_grs_bp": canon_grs_bp,
        "unplaced_grs_bp": unplaced_grs_bp,
        "unplaced_uncertain_bp": unplaced_uncertain_bp,
        "unplaced_core_bp": unplaced_core_bp,
        "min_elim_bp": min_elim_bp,
        "max_elim_bp": max_elim_bp,
        "min_elim_pct": min_elim_pct,
        "max_elim_pct": max_elim_pct,
        "est_elim_bp": est_elim_bp,
        "est_elim_pct": est_elim_pct,
    }

    telo_u = sc.get("telomere_unitigs", {})
    data["telo_unitigs_assembled"] = telo_u.get("total_assembled", 0)
    data["telo_unitigs_assigned"] = telo_u.get("assigned_count", 0)
    data["telo_unitigs_unassigned"] = telo_u.get("unassigned_count", 0)

    data["tier"] = classify_tier(data)
    data["is_finished_gold"] = (data["is_t2t"] and data["is_gapless_non_rdna"])
    return data


def compute_cohort_aggregates(samples_data):
    """Aggregate cohort-wide statistics across all processed samples."""
    total_genomes = len(samples_data)
    if total_genomes == 0:
        return {}

    t2t_pre = sum(1 for s in samples_data if s["pre_t2t_pct"] == 100.0)
    t2t_post = sum(1 for s in samples_data if s["is_t2t"])
    near_t2t_post = sum(
        1 for s in samples_data if s["post_t2t_pct"] >= 85.0 and not s["is_t2t"]
    )

    gapless_pre = sum(1 for s in samples_data if s["baseline_gaps"] == 0)
    gapless_post = sum(1 for s in samples_data if s["is_gapless"])
    gapless_non_rdna_post = sum(1 for s in samples_data if s["is_gapless_non_rdna"])

    gold_finished = sum(1 for s in samples_data if s["is_finished_gold"])

    total_termini = sum(s["total_termini"] for s in samples_data)
    pre_capped_termini = sum(s["pre_capped"] for s in samples_data)
    post_capped_termini = sum(s["post_capped"] for s in samples_data)
    net_telo_added = post_capped_termini - pre_capped_termini
    net_telo_bp = sum(s["net_telo_bp_added"] for s in samples_data)

    total_u_assembled = sum(s.get("telo_unitigs_assembled", 0) for s in samples_data)
    total_u_assigned = sum(s.get("telo_unitigs_assigned", 0) for s in samples_data)
    total_u_unassigned = sum(s.get("telo_unitigs_unassigned", 0) for s in samples_data)
    u_assignment_pct = (
        round((total_u_assigned / total_u_assembled) * 100.0, 1)
        if total_u_assembled > 0
        else 0.0
    )

    baseline_gaps = sum(s["baseline_gaps"] for s in samples_data)
    remaining_gaps = sum(s["remaining_gaps"] for s in samples_data)
    gaps_resolved = sum(s["gaps_resolved"] for s in samples_data)
    gaps_excised = sum(s["gaps_excised"] for s in samples_data)
    gaps_patched = sum(s["gaps_patched"] for s in samples_data)
    patched_bp = sum(s["patched_bp"] for s in samples_data)
    excised_bp = sum(s["excised_bp"] for s in samples_data)

    gap_reduction_pct = (
        round((gaps_resolved / baseline_gaps) * 100.0, 1)
        if baseline_gaps > 0
        else (100.0 if remaining_gaps == 0 else 0.0)
    )

    unplaced_scaffolds_absorbed = sum(
        -s["scaffolds_delta"] for s in samples_data if s["scaffolds_delta"] < 0
    )
    purged_haplotigs_bp = sum(s["purged_haplotigs_bp"] for s in samples_data)
    purged_rdna_bp = sum(s["purged_rdna_bp"] for s in samples_data)
    net_length_delta_bp = sum(s["delta_bp"] for s in samples_data)
    net_rdna_delta_bp = sum(s["rdna_delta_bp"] for s in samples_data)
    net_non_rdna_delta_bp = sum(s["non_rdna_delta_bp"] for s in samples_data)

    conserved_count = sum(1 for s in samples_data if s["mathematically_conserved"])

    tot_core_bp = sum(s.get("core_span_bp", 0) for s in samples_data)
    tot_grs_bp = sum(s.get("grs_span_bp", 0) for s in samples_data)
    tot_compartment_bp = tot_core_bp + tot_grs_bp
    cohort_grs_fraction_pct = round((tot_grs_bp / tot_compartment_bp) * 100.0, 2) if tot_compartment_bp > 0 else 0.0
    weighted_core_cov = sum(s.get("core_span_bp", 0) * s.get("core_median_cov", 0.0) for s in samples_data)
    weighted_grs_cov = sum(s.get("grs_span_bp", 0) * s.get("grs_median_cov", 0.0) for s in samples_data)
    cohort_core_cov = round(weighted_core_cov / tot_core_bp, 1) if tot_core_bp > 0 else 0.0
    cohort_grs_cov = round(weighted_grs_cov / tot_grs_bp, 1) if tot_grs_bp > 0 else 0.0

    tot_canon_grs_bp = sum(s.get("canon_grs_bp", s.get("grs_span_bp", 0)) for s in samples_data)
    tot_unplaced_grs_bp = sum(s.get("unplaced_grs_bp", 0) for s in samples_data)
    tot_unplaced_uncertain_bp = sum(s.get("unplaced_uncertain_bp", 0) for s in samples_data)
    tot_unplaced_core_bp = sum(s.get("unplaced_core_bp", 0) for s in samples_data)
    tot_min_elim_bp = sum(s.get("min_elim_bp", s.get("grs_span_bp", 0)) for s in samples_data)
    tot_max_elim_bp = sum(s.get("max_elim_bp", s.get("grs_span_bp", 0)) for s in samples_data)
    tot_est_elim_bp = sum(s.get("est_elim_bp", s.get("grs_span_bp", 0)) for s in samples_data)
    cohort_min_elim_pct = round((tot_min_elim_bp / tot_compartment_bp) * 100.0, 2) if tot_compartment_bp > 0 else 0.0
    cohort_max_elim_pct = round((tot_max_elim_bp / tot_compartment_bp) * 100.0, 2) if tot_compartment_bp > 0 else 0.0
    cohort_est_elim_pct = round((tot_est_elim_bp / tot_compartment_bp) * 100.0, 2) if tot_compartment_bp > 0 else 0.0

    tier_counts = defaultdict(int)
    for s in samples_data:
        tier_counts[s["tier"]] += 1

    return {
        "total_genomes": total_genomes,
        "gold_finished_genomes": gold_finished,
        "gold_finished_pct": round((gold_finished / total_genomes) * 100.0, 1),
        "t2t_genomes_pre": t2t_pre,
        "t2t_genomes_post": t2t_post,
        "t2t_genomes_post_pct": round((t2t_post / total_genomes) * 100.0, 1),
        "near_t2t_genomes_post": near_t2t_post,
        "gapless_genomes_pre": gapless_pre,
        "gapless_genomes_post": gapless_post,
        "gapless_genomes_post_pct": round((gapless_post / total_genomes) * 100.0, 1),
        "gapless_non_rdna_post": gapless_non_rdna_post,
        "gapless_non_rdna_post_pct": round((gapless_non_rdna_post / total_genomes) * 100.0, 1),
        "total_termini": total_termini,
        "pre_capped_termini": pre_capped_termini,
        "post_capped_termini": post_capped_termini,
        "net_telo_added": net_telo_added,
        "cohort_t2t_completeness_pct": (
            round((post_capped_termini / total_termini) * 100.0, 1)
            if total_termini > 0
            else 0.0
        ),
        "net_telo_bp": net_telo_bp,
        "baseline_gaps": baseline_gaps,
        "remaining_gaps": remaining_gaps,
        "gaps_resolved": gaps_resolved,
        "gaps_excised": gaps_excised,
        "gaps_patched": gaps_patched,
        "patched_bp": patched_bp,
        "excised_bp": excised_bp,
        "gap_reduction_pct": gap_reduction_pct,
        "unplaced_scaffolds_absorbed": unplaced_scaffolds_absorbed,
        "purged_haplotigs_bp": purged_haplotigs_bp,
        "purged_rdna_bp": purged_rdna_bp,
        "net_length_delta_bp": net_length_delta_bp,
        "net_rdna_delta_bp": net_rdna_delta_bp,
        "net_non_rdna_delta_bp": net_non_rdna_delta_bp,
        "tot_core_bp": tot_core_bp,
        "tot_grs_bp": tot_grs_bp,
        "cohort_grs_fraction_pct": cohort_grs_fraction_pct,
        "cohort_core_cov": cohort_core_cov,
        "cohort_grs_cov": cohort_grs_cov,
        "tot_canon_grs_bp": tot_canon_grs_bp,
        "tot_unplaced_grs_bp": tot_unplaced_grs_bp,
        "tot_unplaced_uncertain_bp": tot_unplaced_uncertain_bp,
        "tot_unplaced_core_bp": tot_unplaced_core_bp,
        "tot_min_elim_bp": tot_min_elim_bp,
        "tot_max_elim_bp": tot_max_elim_bp,
        "tot_est_elim_bp": tot_est_elim_bp,
        "cohort_min_elim_pct": cohort_min_elim_pct,
        "cohort_max_elim_pct": cohort_max_elim_pct,
        "cohort_est_elim_pct": cohort_est_elim_pct,
        "total_unitigs_assembled": total_u_assembled,
        "total_unitigs_assigned": total_u_assigned,
        "total_unitigs_unassigned": total_u_unassigned,
        "unitig_assignment_pct": u_assignment_pct,
        "mathematically_conserved_count": conserved_count,
        "mathematically_conserved_all": (conserved_count == total_genomes),
        "tier_distribution": dict(tier_counts),
    }


def generate_markdown_report(prefix, agg, samples_data):
    """Generate executive Markdown curation summary report."""
    md = []
    md.append(f"# Automated Curation Cohort Summary: {prefix}\n")
    md.append(
        "Executive synthesis and quality scorecard for automated chromosome-scale assembly "
        "curation, gap resolution, and telomere capping across the cohort.\n"
    )

    # 1. KPI Highlight Cards
    md.append("## Executive Key Performance Indicators\n")
    md.append("| Metric | Baseline | Curated Outcome | Cohort Impact / Gain |")
    md.append("| :--- | :--- | :--- | :--- |")
    md.append(
        f"| **Finished Genomes (T2T + Non-rDNA Gapless)** | "
        f"{agg.get('gold_finished_genomes', 0)} / {agg['total_genomes']} | "
        f"**{agg.get('gold_finished_genomes', 0)} / {agg['total_genomes']} ({agg.get('gold_finished_pct', 0.0)}%)** | "
        f"Gold standard assemblies ready for NCBI/ENA release (excluding expected rDNA sizing gaps) |"
    )
    md.append(
        f"| **Telomere-to-Telomere Genomes** | "
        f"{agg['t2t_genomes_pre']} / {agg['total_genomes']} | "
        f"**{agg['t2t_genomes_post']} / {agg['total_genomes']} ({agg['t2t_genomes_post_pct']}%)** | "
        f"+{agg['t2t_genomes_post'] - agg['t2t_genomes_pre']} assemblies ({agg.get('near_t2t_genomes_post', 0)} near-T2T) |"
    )
    md.append(
        f"| **Non-rDNA Gapless Genomes** | "
        f"{agg['gapless_genomes_pre']} / {agg['total_genomes']} | "
        f"**{agg['gapless_non_rdna_post']} / {agg['total_genomes']} ({agg['gapless_non_rdna_post_pct']}%)** | "
        f"+{agg['gapless_non_rdna_post'] - agg['gapless_genomes_pre']} non-rDNA gapless assemblies ({agg['gapless_genomes_post']} strictly 0 Ns) |"
    )
    md.append(
        f"| **Chromosome Termini Capped** | "
        f"{agg['pre_capped_termini']} / {agg['total_termini']} | "
        f"**{agg['post_capped_termini']} / {agg['total_termini']} ({agg['cohort_t2t_completeness_pct']}%)** | "
        f"+{agg['net_telo_added']} ends (+{agg['net_telo_bp']:,} bp telomeric tracts) |"
    )
    if agg.get("total_unitigs_assembled", 0) > 0:
        md.append(
            f"| **Telomeric Unitigs Assigned** | "
            f"{agg['total_unitigs_assembled']} assembled | "
            f"**{agg['total_unitigs_assigned']} / {agg['total_unitigs_assembled']} assigned ({agg['unitig_assignment_pct']}%)** | "
            f"{agg['total_unitigs_unassigned']} unassigned unitigs across cohort |"
        )
    md.append(
        f"| **Assembly Gaps Resolved** | "
        f"{agg['baseline_gaps']} baseline gaps | "
        f"**{agg['remaining_gaps']} remaining gaps** | "
        f"**{agg['gaps_resolved']} gaps resolved ({agg['gap_reduction_pct']}%)** |"
    )
    md.append(
        f"| **Gap Resolution Breakdown** | "
        f"Excised: {agg['gaps_excised']} | "
        f"Patched: {agg['gaps_patched']} (+{agg['patched_bp']:,} bp) | "
        f"Excised Ns: {agg['excised_bp']:,} bp |"
    )
    unplaced_str = f"{agg['unplaced_scaffolds_absorbed']} scaffolds assimilated / retired"
    md.append(
        f"| **Scaffold Assimilation & Purity** | "
        f"Baseline unplaced | "
        f"{unplaced_str} | "
        f"Net non-rDNA delta: {agg['net_non_rdna_delta_bp']:+,} bp (Total net delta: {agg['net_length_delta_bp']:+,} bp) |"
    )
    conserved_badge = "100% Verified (0 Discrepancy)" if agg["mathematically_conserved_all"] else f"WARNING: {agg['total_genomes'] - agg['mathematically_conserved_count']} Mismatches"
    md.append(
        f"| **Mathematical Conservation** | "
        f"Audit reconciliation | "
        f"**{conserved_badge}** | "
        f"{agg['mathematically_conserved_count']}/{agg['total_genomes']} assemblies strictly conserved |"
    )
    if agg.get("tot_core_bp", 0) > 0 or agg.get("tot_grs_bp", 0) > 0:
        tot_comp_mb = (agg["tot_core_bp"] + agg["tot_grs_bp"]) / 1e6
        core_mb = agg["tot_core_bp"] / 1e6
        grs_mb = agg["tot_grs_bp"] / 1e6
        min_elim_mb = agg.get("tot_min_elim_bp", agg["tot_grs_bp"]) / 1e6
        max_elim_mb = agg.get("tot_max_elim_bp", agg["tot_grs_bp"]) / 1e6
        min_pct = agg.get("cohort_min_elim_pct", agg["cohort_grs_fraction_pct"])
        max_pct = agg.get("cohort_max_elim_pct", agg["cohort_grs_fraction_pct"])
        if min_elim_mb != max_elim_mb:
            elim_str = f"**{min_elim_mb:.2f} – {max_elim_mb:.2f} Mb ({min_pct:.1f}%–{max_pct:.1f}%)**"
        else:
            elim_str = f"**{grs_mb:.2f} Mb ({agg['cohort_grs_fraction_pct']:.1f}%)**"
        md.append(
            f"| **Programmed DNA Elimination (PDE)** | "
            f"Core: {core_mb:.2f} Mb (`{agg['cohort_core_cov']:.1f}×`) | "
            f"Eliminated GRS: {elim_str} (`{agg['cohort_grs_cov']:.1f}×`) | "
            f"Total classified: {tot_comp_mb:.2f} Mb across {agg['total_genomes']} assemblies |"
        )
    md.append("")

    # 2. Curation Tiers Summary
    md.append("## Cohort Curation League Table\n")
    md.append("Assemblies are categorised by their final curation outcome:")
    for tier_name, count in sorted(agg.get("tier_distribution", {}).items()):
        pct = round((count / agg["total_genomes"]) * 100.0, 1)
        md.append(f"- **{tier_name}**: `{count}` assemblies ({pct}%)")
    md.append("")

    # 3. Master Scorecard Matrix Table
    md.append("### Per-Genome Curation Matrix\n")
    md.append(
        "| Sample | Tier | Chrs | Gaps (Pre &rarr; Post) | Gaps Res. (%) | Gapless? | T2T Termini | T2T Score (%) | Finished? | Uncapped Ends | Non-rDNA Delta (bp) | Total Delta (bp) | Conserved? |"
    )
    md.append(
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :---: | :---: | :---: |"
    )

    for s in sorted(samples_data, key=lambda x: (x["tier"], -x["post_t2t_pct"], x["remaining_non_rdna_gaps"])):
        gaps_str = f"{s['baseline_gaps']} &rarr; {s['remaining_gaps']}"
        t2t_str = f"{s['post_capped']}/{s['total_termini']}"
        t2t_score_str = f"{s['post_t2t_pct']}%"
        if s["t2t_gain_pct"] > 0:
            t2t_score_str += f" (+{s['t2t_gain_pct']}%)"
        if s["is_gapless"]:
            gapless_badge = "Yes"
        elif s["is_gapless_non_rdna"]:
            gapless_badge = f"Yes* ({s['rdna_gaps']} rDNA)"
        else:
            gapless_badge = f"No ({s['remaining_non_rdna_gaps']})"
        finished_badge = "**Gold**" if s["is_finished_gold"] else ("Silver" if s["is_t2t"] or s["is_gapless_non_rdna"] else "No")
        uncapped_str = ", ".join(s["uncapped_ends"]) if s["uncapped_ends"] else "None"
        conserved_str = "Yes" if s["mathematically_conserved"] else "No"
        non_rdna_delta_str = f"{s['non_rdna_delta_bp']:+,}"
        delta_str = f"{s['delta_bp']:+,}"

        md.append(
            f"| `{s['sample']}` | {s['tier']} | {s['chromosomes_count']} | "
            f"{gaps_str} | {s['gap_resolution_pct']}% | {gapless_badge} | "
            f"{t2t_str} | {t2t_score_str} | {finished_badge} | "
            f"`{uncapped_str}` | {non_rdna_delta_str} | {delta_str} | {conserved_str} |"
        )
    md.append("")
    md.append("*\\*Note: `Yes* (1 rDNA)` indicates that all non-rDNA sequence is completely gapless, with only an expected 100 bp INSDC/ENA rDNA array sizing placeholder remaining.*\n")

    # 4. Actionable Uncapped Termini Diagnostic / Triage Guide
    uncapped_samples = [s for s in samples_data if not s["is_t2t"]]
    if uncapped_samples:
        md.append("## Uncapped Termini Triage Guide\n")
        md.append(
            "The following chromosome ends remain uncapped after automated telomere walking. "
            "Inspect the biological context and recommended triage actions below:\n"
        )
        md.append("| Sample | Chromosome End | Baseline Status | Post Status | Actionable Triage Recommendation |")
        md.append("| :--- | :--- | :--- | :--- | :--- |")

        for s in uncapped_samples:
            for end in s["uncapped_ends"]:
                term_status = s["terminus_status_breakdown"].get(end, "UNCAPPED")
                recommendation = (
                    "Inspect satellite / rDNA junction repeat; sequence across terminal gap with targeted ultralong reads"
                    if "SATELLITE" in term_status or "RDNA" in term_status
                    else "Check PretextView for transchromosomal chimeric junction or distal telomeric cluster"
                )
                md.append(
                    f"| `{s['sample']}` | `{end}` | UNCAPPED | `{term_status}` | {recommendation} |"
                )
        md.append("")
    else:
        md.append("## Uncapped Termini Triage Guide\n")
        md.append("All chromosome ends across the entire cohort achieved complete telomere capping (100% T2T).\n")

    # 5. Gap Resolution Ledger
    remaining_gap_samples = [s for s in samples_data if s["remaining_non_rdna_gaps"] > 0]
    if remaining_gap_samples:
        md.append("## Remaining Gaps Diagnostic\n")
        md.append(
            "The following assemblies retain unresolved intra-chromosomal non-rDNA gaps:\n"
        )
        md.append("| Sample | Baseline Gaps | Patched | Excised | Remaining Non-rDNA | Remaining Total | Biological Context |")
        md.append("| :--- | :---: | :---: | :---: | :---: | :---: | :--- |")
        for s in remaining_gap_samples:
            md.append(
                f"| `{s['sample']}` | {s['baseline_gaps']} | {s['gaps_patched']} | {s['gaps_excised']} | "
                f"**{s['remaining_non_rdna_gaps']}** | {s['remaining_gaps']} | "
                f"Complex repetitive non-rDNA region / segmental duplication |"
            )
        md.append("")
        md.append("*Note: Assemblies where the only remaining gap is an expected 100 bp tandem rDNA array sizing gap are omitted from this diagnostic as non-rDNA regions are fully resolved.*\n")
    else:
        md.append("## Remaining Gaps Diagnostic\n")
        md.append("All assemblies across the entire cohort are 100% non-rDNA gapless (0 unresolved non-rDNA Ns).\n")

    # 6. Structural Surgery & Repeat Purge Summary
    md.append("## Structural Surgery & Purge Summary\n")
    md.append("| Sample | Contig Transfers | Contig Fissions | Purged Haplotigs (bp) | Purged rDNA Scaffolds (bp) | Primary Purity (%) |")
    md.append("| :--- | :---: | :---: | :---: | :---: | :---: |")
    for s in samples_data:
        md.append(
            f"| `{s['sample']}` | {s['transfers_count']} | {s['fissions_count']} | "
            f"{s['purged_haplotigs_bp']:,} bp | {s['purged_rdna_bp']:,} bp | "
            f"{s['primary_purity_pct']}% |"
        )
    md.append("")

    # 7. Programmed DNA Elimination (PDE) Compartments
    if any(s.get("core_span_bp", 0) > 0 or s.get("grs_span_bp", 0) > 0 for s in samples_data):
        md.append("## Programmed DNA Elimination (PDE) Compartment Allocation\n")
        md.append(
            "Genome-wide partition between retained somatic core chromosomes and germline-restricted sequences (GRS) "
            "eliminated during developmental diminution, classified via EM mixture modelling and neo-telomeric addition breakpoints:\n"
        )
        md.append(
            "| Sample | Threshold ($T_{\\text{GRS}}$) | Core Span (bp) | Core Cov | Canon GRS (bp) | Unplaced GRS (bp) | Uncertain (bp) | Min Elim (bp) | Max Elim (bp) | Est Elim (bp) | Elim Range (%) |"
        )
        md.append(
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |"
        )
        for s in sorted(samples_data, key=lambda x: x["sample"]):
            c_bp = s.get("core_span_bp", 0)
            c_cov = s.get("core_median_cov", 0.0)
            t_grs = s.get("grs_threshold", 0.0)
            canon_bp = s.get("canon_grs_bp", s.get("grs_span_bp", 0))
            unpl_grs_bp = s.get("unplaced_grs_bp", 0)
            uncert_bp = s.get("unplaced_uncertain_bp", 0)
            min_bp = s.get("min_elim_bp", s.get("grs_span_bp", 0))
            max_bp = s.get("max_elim_bp", s.get("grs_span_bp", 0))
            est_bp = s.get("est_elim_bp", s.get("grs_span_bp", 0))
            min_pct = s.get("min_elim_pct", s.get("grs_fraction_pct", 0.0))
            max_pct = s.get("max_elim_pct", s.get("grs_fraction_pct", 0.0))
            range_pct_str = f"**{min_pct:.1f}%–{max_pct:.1f}%**" if min_pct != max_pct else f"**{min_pct:.1f}%**"
            md.append(
                f"| `{s['sample']}` | `{t_grs:.2f}×` | {c_bp:,d} | `{c_cov:.1f}×` | {canon_bp:,d} | {unpl_grs_bp:,d} | {uncert_bp:,d} | {min_bp:,d} | {max_bp:,d} | {est_bp:,d} | {range_pct_str} |"
            )
        tot_c_bp = agg.get("tot_core_bp", 0)
        cohort_c_cov = agg.get("cohort_core_cov", 0.0)
        tot_canon = agg.get("tot_canon_grs_bp", agg.get("tot_grs_bp", 0))
        tot_unpl = agg.get("tot_unplaced_grs_bp", 0)
        tot_uncert = agg.get("tot_unplaced_uncertain_bp", 0)
        tot_min = agg.get("tot_min_elim_bp", agg.get("tot_grs_bp", 0))
        tot_max = agg.get("tot_max_elim_bp", agg.get("tot_grs_bp", 0))
        tot_est = agg.get("tot_est_elim_bp", agg.get("tot_grs_bp", 0))
        cohort_min_pct = agg.get("cohort_min_elim_pct", agg.get("cohort_grs_fraction_pct", 0.0))
        cohort_max_pct = agg.get("cohort_max_elim_pct", agg.get("cohort_grs_fraction_pct", 0.0))
        cohort_range_str = f"**{cohort_min_pct:.1f}%–{cohort_max_pct:.1f}%**" if cohort_min_pct != cohort_max_pct else f"**{cohort_min_pct:.1f}%**"
        md.append(
            f"| **Cohort Total** | — | **{tot_c_bp:,d}** | **`{cohort_c_cov:.1f}×`** | **{tot_canon:,d}** | **{tot_unpl:,d}** | **{tot_uncert:,d}** | **{tot_min:,d}** | **{tot_max:,d}** | **{tot_est:,d}** | {cohort_range_str} |\n"
        )

    return "\n".join(md)


def generate_tsv_table(samples_data):
    """Generate TSV tabular matrix of per-sample metrics."""
    fields = [
        "sample",
        "tier",
        "chromosomes_count",
        "total_termini",
        "pre_capped_termini",
        "post_capped_termini",
        "telo_unitigs_assembled",
        "telo_unitigs_assigned",
        "telo_unitigs_unassigned",
        "pre_t2t_pct",
        "post_t2t_pct",
        "t2t_gain_pct",
        "is_t2t",
        "uncapped_ends",
        "net_telo_bp_added",
        "baseline_gaps",
        "remaining_gaps",
        "rdna_gaps",
        "remaining_non_rdna_gaps",
        "gaps_resolved",
        "gaps_excised",
        "gaps_patched",
        "gap_resolution_pct",
        "is_gapless",
        "is_gapless_non_rdna",
        "is_finished_gold",
        "patched_bp",
        "excised_bp",
        "pre_total_bp",
        "post_total_bp",
        "net_length_delta_bp",
        "non_rdna_length_delta_bp",
        "rdna_length_delta_bp",
        "pre_scaffolds",
        "post_scaffolds",
        "scaffolds_delta",
        "pre_n50",
        "post_n50",
        "pre_chr_n50",
        "post_chr_n50",
        "pre_n_bp",
        "post_n_bp",
        "transfers_count",
        "fissions_count",
        "transferred_bp",
        "severed_bp",
        "purged_haplotigs_count",
        "purged_haplotigs_bp",
        "purged_rdna_count",
        "purged_rdna_bp",
        "primary_purity_pct",
        "mathematically_conserved",
        "discrepancy_bp",
        "core_span_bp",
        "core_median_cov",
        "grs_span_bp",
        "grs_median_cov",
        "grs_threshold",
        "grs_fraction_pct",
        "canon_grs_bp",
        "unplaced_grs_bp",
        "unplaced_uncertain_bp",
        "unplaced_core_bp",
        "min_elim_bp",
        "max_elim_bp",
        "min_elim_pct",
        "max_elim_pct",
        "est_elim_bp",
        "est_elim_pct",
    ]

    lines = ["\t".join(fields)]
    for s in sorted(samples_data, key=lambda x: x["sample"]):
        row = []
        for f in fields:
            val = s.get(f, "")
            if f == "uncapped_ends":
                row.append(",".join(val) if val else "None")
            elif f == "net_length_delta_bp":
                row.append(str(s.get("delta_bp", 0)))
            elif f == "non_rdna_length_delta_bp":
                row.append(str(s.get("non_rdna_delta_bp", 0)))
            elif f == "rdna_length_delta_bp":
                row.append(str(s.get("rdna_delta_bp", 0)))
            elif isinstance(val, bool):
                row.append("TRUE" if val else "FALSE")
            else:
                row.append(str(val))
        lines.append("\t".join(row))

    return "\n".join(lines) + "\n"


def load_validation(filepath):
    """Load a single sample curation validation JSON file."""
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        print(f"WARNING: Failed to parse validation {filepath}: {exc}", file=sys.stderr)
        return None


def render_validation_summary_section(validations_map):
    """Renders cohort-level validation section for executive Markdown report."""
    if not validations_map:
        lines = []
        lines.append("\n## 7. Post-Curation Mapping & Structural Variant Validation\n\n")
        lines.append("*Post-curation read re-mapping and Sniffles structural variant validation was skipped (`--skip_validation`).*\n")
        return "".join(lines)

    lines = []
    lines.append("\n## 7. Post-Curation Mapping & Structural Variant Validation\n")
    lines.append("*Independent re-mapping of HiFi reads against curated assemblies and Sniffles SV calling across surgical loci:*\n")
    lines.append("| Sample | Surgical Edits Evaluated | Concordant (Pass) | Warnings | Flagrant Errors ($VAF \\ge 0.8$) | SV Concordance |")
    lines.append("| :--- | :---: | :---: | :---: | :---: | :---: |")

    tot_edits = 0
    tot_pass = 0
    tot_warn = 0
    tot_errors = 0

    for sample in sorted(validations_map.keys()):
        val = validations_map[sample]
        s = val.get("summary", {})
        e_eval = s.get("total_edits_evaluated", 0)
        e_pass = s.get("edits_passed", 0)
        e_warn = s.get("edits_warned", 0)
        e_err = s.get("flagrant_errors_count", 0)
        conc = s.get("edit_concordance_pct", 100.0)

        tot_edits += e_eval
        tot_pass += e_pass
        tot_warn += e_warn
        tot_errors += e_err

        status_flag = f"**{conc}%**" if conc >= 95.0 else f"<span style='color:orange;'>**{conc}%**</span>"
        err_str = f"**{e_err}**" if e_err > 0 else "0"
        lines.append(f"| `{sample}` | {e_eval} | {e_pass} | {e_warn} | {err_str} | {status_flag} |")

    cohort_conc = round((tot_pass / tot_edits * 100.0), 1) if tot_edits > 0 else 100.0
    lines.append(f"| **Cohort Total** | **{tot_edits}** | **{tot_pass}** | **{tot_warn}** | **{tot_errors}** | **{cohort_conc}%** |\n")
    return "\n".join(lines)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    samples_data = []
    for sc_file in args.scorecards:
        sc = load_scorecard(sc_file)
        if sc:
            samples_data.append(parse_sample_record(sc))

    if not samples_data:
        print("ERROR: No valid sample scorecards were loaded.", file=sys.stderr)
        sys.exit(1)

    print(
        f"[{args.prefix}] Processing cohort summary across {len(samples_data)} sample scorecards...",
        flush=True,
    )

    validations_map = {}
    for val_file in getattr(args, "validations", []):
        val = load_validation(val_file)
        if val and "sample" in val:
            validations_map[val["sample"]] = val

    agg = compute_cohort_aggregates(samples_data)

    # 1. Emit Markdown Report
    md_content = generate_markdown_report(args.prefix, agg, samples_data)
    md_content += render_validation_summary_section(validations_map)

    md_path = os.path.join(args.outdir, f"{args.prefix}.curation_summary.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)
    print(f"[{args.prefix}] Wrote executive Markdown report: {md_path}", flush=True)

    # 2. Emit TSV Table
    tsv_content = generate_tsv_table(samples_data)
    tsv_path = os.path.join(args.outdir, f"{args.prefix}.curation_summary.tsv")
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write(tsv_content)
    print(f"[{args.prefix}] Wrote TSV matrix: {tsv_path}", flush=True)

    # 3. Emit Comprehensive JSON Payload
    json_payload = {
        "provenance": {
            "software": "pde_curation_summary.py",
            "software_version": __version__,
            "script_fingerprint": get_script_fingerprint(),
        },
        "cohort_prefix": args.prefix,
        "cohort_aggregates": agg,
        "samples": samples_data,
        "validations": validations_map,
    }
    json_path = os.path.join(args.outdir, f"{args.prefix}.curation_summary.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_payload, f, indent=2, sort_keys=True)
    print(f"[{args.prefix}] Wrote JSON payload: {json_path}", flush=True)


if __name__ == "__main__":
    main()
