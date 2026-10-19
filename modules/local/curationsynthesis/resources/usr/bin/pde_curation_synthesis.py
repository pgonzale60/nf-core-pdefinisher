#!/usr/bin/env python3
"""pde_curation_synthesis.py

Module 2 Curation Readiness & Synthesis Engine (Thin Orchestrator).
Orchestrates specialised sub-engines from `resources/usr/lib/curationsynthesis`:
  - `policy`: Policy loading & deterministic package fingerprinting
  - `models`: Core Open Port Graph domain models (`Port`, `UnplacedSegment`)
  - `evidence_io`: Pure Module 1 reconnaissance loaders and interval utilities
  - `port_valency`: Molecule deduplication, scarcity-aware unitig arbitration & port graph
  - `internal_resolver`: Internal gap, BND/INV misjoin, macro SV & rDNA array card synthesis
  - `terminus_resolver`: 7-level chromosome terminus resolution hierarchy & chimera pruning
  - `unplaced_triage`: Redundancy triage & chromosomal port extension candidacy
  - `reporting`: Markdown dossier and JSON manifest rendering

British spelling standard enforced across all documentation and comments.
"""

import argparse
import gc
import json
import os
from pathlib import Path
import sys

# Add resources/usr/lib to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

from curationsynthesis import (
    __version__,
    build_manifest_json,
    build_unitig_assignment,
    evaluate_extension_candidacy,
    get_script_fingerprint as _pkg_fingerprint,
    load_delfies_points,
    load_flye_info,
    load_grs_summary,
    load_mosdepth_tip_depths,
    load_policy as _load_policy,
    load_rdna_bed,
    load_red_repeats,
    load_telo_junctions,
    load_terminal_softclips,
    load_terminal_summary,
    match_open_termini_to_scaffolds,
    parse_telo_paf,
    parse_unplaced_paf,
    partition_chromosomes,
    reconstruct_rdna_operons,
    render_markdown_blueprint,
    survey_telomere_unitigs,
    synthesize_curation_cards,
    synthesize_rdna_cards,
    synthesize_telomere_cards,
    triage_unplaced_scaffolds,
)


def get_script_fingerprint():
    return _pkg_fingerprint(__file__)


def load_policy(policy_path):
    return _load_policy(policy_path, caller_file=__file__)


def parse_args():
    parser = argparse.ArgumentParser(description="Compile Module 2 Curation Synthesis Blueprint.")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__} (sha256:{get_script_fingerprint()})",
    )
    parser.add_argument("--sample", required=True, help="Sample ID (e.g. nxOscSper1_DF5120)")
    parser.add_argument("--ref", required=True, help="Reference FASTA path")
    parser.add_argument("--grs-summary", required=True, help="Path to grs_summary.tsv from identifyputativegrs")
    parser.add_argument("--mosdepth-regions", required=True, help="Path to mosdepth regions.bed.gz")
    parser.add_argument("--terminal-summary", required=True, help="Path to terminal_ends_summary.tsv from softclipends")
    parser.add_argument("--terminal-softclips", required=True, help="Path to terminal_softclips.tsv from softclipends")
    parser.add_argument("--delfies-bed", required=True, help="Path to breakpoint_locations.bed from delfies")
    parser.add_argument("--rdna-bed", required=True, help="Path to rdna_regions.bed from extractrdnareads")
    parser.add_argument("--barrnap-gff", required=True, help="Path to rRNA.gff from barrnap")
    parser.add_argument("--flye-info", required=True, help="Path to assembly_info.txt from flye")
    parser.add_argument("--unplaced-paf", required=True, help="Path to unplaced contigs alignment PAF")
    parser.add_argument("--red-bed", required=True, help="Path to repeat masking BED from red")
    parser.add_argument("--telo-paf", required=True, help="Path to telomeric unitigs alignment PAF from minimap2")
    parser.add_argument("--gap-evidence-json", default=None, help="Optional path to gap_evidence.json")
    parser.add_argument("--output-md", default=None, help="Optional path to write markdown summary")
    parser.add_argument("--output-json", default=None, help="Optional path to write curation_manifest.json")
    parser.add_argument("--policy", default=None, help="Path to curation_policy.json (single source of truth)")
    parser.add_argument("--telo-motif", default=None, help="Canonical telomeric motif (optional override)")
    parser.add_argument("--telo-junctions", default=None, help="Path to telomeric junction catalog JSON")
    parser.add_argument("--telo-unitigs", default=None, help="Optional path to telomeric unitigs FASTA")
    parser.add_argument("--chromosomes", default=None, help="Comma-separated list of canonical nuclear chromosome names")
    return parser.parse_args()


def main():
    args = parse_args()
    policy = load_policy(args.policy)
    if args.telo_motif:
        policy.setdefault("telomere", {})["motif"] = args.telo_motif

    # 1. Load baseline karyotype, Core/GRS segmentation & terminal coverage
    (
        chrom_lens,
        core_segments,
        grs_segments,
        grs_arm_depths,
        local_depth_lookup,
        grs_threshold,
    ) = load_grs_summary(args.grs_summary)
    tip_depths_1kb = load_mosdepth_tip_depths(args.mosdepth_regions, chrom_lens)
    sorted_chroms, nuclear_chroms, unplaced_chroms = partition_chromosomes(
        chrom_lens, args.chromosomes, args.unplaced_paf
    )

    tot_core_len = sum(sum(seg[2] for seg in core_segments[c]) for c in sorted_chroms)
    core_weighted_depth = sum(sum(seg[2] * seg[3] for seg in core_segments[c]) for c in sorted_chroms)
    cohort_core_d = (core_weighted_depth / tot_core_len) if tot_core_len > 0 else 100.0

    # 2. Synthesise Domain 1 Internal Action Cards (Gaps & SVs)
    rdna_intervals, rdna_major_regions = load_rdna_bed(args.rdna_bed, nuclear_chroms)
    internal_cards = []
    gap_provenance = None
    gap_parse_warning = None
    if args.gap_evidence_json and os.path.exists(args.gap_evidence_json):
        try:
            with open(args.gap_evidence_json) as f:
                evidence_data = json.load(f)
            internal_cards = synthesize_curation_cards(
                evidence_data,
                policy=policy,
                local_depth_lookup=local_depth_lookup,
                genome_core_depth=cohort_core_d,
                chrom_lens=chrom_lens,
                rdna_intervals=rdna_intervals,
            )
            gap_provenance = evidence_data.get("provenance", {})
            del evidence_data
            gc.collect()
        except Exception as e:
            gap_parse_warning = str(e)

    # 3. Synthesise Domain 2 Chromosome Terminus Action Cards
    baseline_status, tip_tract_bp = load_terminal_summary(args.terminal_summary)
    raw_softclips, single_read_ext = load_terminal_softclips(args.terminal_softclips, policy)
    unitig_extensions = parse_telo_paf(args.telo_paf, nuclear_chroms, policy)
    delfies_points = load_delfies_points(args.delfies_bed, nuclear_chroms)
    telo_junctions = load_telo_junctions(args.telo_junctions)

    unitig_assignment = build_unitig_assignment(
        nuclear_chroms, chrom_lens, grs_segments, core_segments, delfies_points, unitig_extensions, policy
    )
    telomere_unitig_inventory = survey_telomere_unitigs(
        telo_paf=args.telo_paf,
        nuclear_chroms=nuclear_chroms,
        chrom_lens=chrom_lens,
        unitig_assignment=unitig_assignment,
        unitig_extensions=unitig_extensions,
        telo_junctions=telo_junctions,
        telo_unitigs_fasta=args.telo_unitigs,
        policy=policy,
    )

    _scaff_align, scaff_primary_records, scaff_any_records = parse_unplaced_paf(args.unplaced_paf)
    min_ext_bp = policy.get("telomere", {}).get("min_extension_bp", 200)
    rdna_dist_bp = policy.get("repeat_and_rdna", {}).get("rdna_dist_from_end", 25000)

    open_termini = []
    for c in nuclear_chroms:
        c_len = chrom_lens.get(c, 0)
        for e in ("5p", "3p"):
            if baseline_status.get((c, e), "UNKNOWN") == "ALREADY_CAPPED":
                continue
            if any(
                (e == "5p" and r_st <= rdna_dist_bp) or (e == "3p" and r_en >= c_len - rdna_dist_bp)
                for r_st, r_en in rdna_intervals.get(c, [])
            ):
                continue
            best_u = unitig_assignment.get((c, e))
            if best_u and best_u[0] >= min_ext_bp and best_u[1] != "None":
                continue
            open_termini.append((c, e))

    terminal_matches, scaffold_attributions = match_open_termini_to_scaffolds(
        nuclear_chroms,
        chrom_lens,
        unplaced_chroms,
        scaff_primary_records,
        telo_junctions,
        open_termini,
        unitig_assignment=unitig_assignment,
        policy=policy,
        grs_segments=grs_segments,
    )

    telo_cards = synthesize_telomere_cards(
        nuclear_chroms,
        chrom_lens,
        grs_segments,
        core_segments,
        delfies_points,
        rdna_intervals,
        baseline_status,
        unitig_extensions,
        single_read_ext,
        tip_tract_bp=tip_tract_bp,
        policy=policy,
        start_card_idx=1,
        telo_junctions=telo_junctions,
        terminal_matches=terminal_matches,
        raw_softclips=raw_softclips,
        internal_cards=internal_cards,
    )

    # 4. Synthesise Domain 3 Major rDNA Locus Audit & Replacement Cards
    scaff_rdna_operons, scaff_rdna_operon_bp, all_r_units = reconstruct_rdna_operons(
        args.barrnap_gff, chrom_lens, max_spacer_bp=3500, terminal_padding_bp=1000
    )
    flye_data = load_flye_info(args.flye_info)
    rdna_cards, rdna_locus_audits = synthesize_rdna_cards(
        rdna_major_regions, chrom_lens, all_r_units, flye_data, policy
    )

    # 5. Synthesise Domain 4 Unplaced Scaffolds Triage & Port Integration Cards
    scaff_repeats = load_red_repeats(args.red_bed)
    triage_result = triage_unplaced_scaffolds(
        unplaced_chroms,
        chrom_lens,
        scaff_repeats,
        scaff_any_records,
        scaff_rdna_operons,
        scaff_rdna_operon_bp,
    )
    scaffold_cards, scaffold_dispositions, extension_by_seq = evaluate_extension_candidacy(
        triage_result,
        scaff_primary_records,
        scaff_rdna_operon_bp,
        scaffold_attributions,
        telo_cards,
        rdna_cards,
        nuclear_chroms=nuclear_chroms,
        chrom_lens=chrom_lens,
        scaff_any_records=scaff_any_records,
        internal_cards=internal_cards,
        open_termini=open_termini,
        policy=policy,
    )

    # 6. Render Markdown Blueprint and Machine-Executable JSON Manifest
    full_report, karyotype_stats = render_markdown_blueprint(
        sample_id=args.sample,
        sorted_chroms=sorted_chroms,
        nuclear_chroms=nuclear_chroms,
        unplaced_chroms=unplaced_chroms,
        chrom_lens=chrom_lens,
        core_segments=core_segments,
        grs_segments=grs_segments,
        grs_threshold=grs_threshold,
        grs_arm_depths=grs_arm_depths,
        tip_depths_1kb=tip_depths_1kb,
        baseline_status=baseline_status,
        unitig_assignment=unitig_assignment,
        single_read_ext=single_read_ext,
        rdna_intervals=rdna_intervals,
        rdna_dist_bp=rdna_dist_bp,
        terminal_matches=terminal_matches,
        internal_cards=internal_cards,
        gap_provenance=gap_provenance,
        gap_json_path=args.gap_evidence_json,
        gap_parse_warning=gap_parse_warning,
        telo_cards=telo_cards,
        flye_info_exists=bool(args.flye_info and os.path.exists(args.flye_info)),
        flye_data=flye_data,
        rdna_locus_audits=rdna_locus_audits,
        rdna_cards=rdna_cards,
        triage_result=triage_result,
        scaffold_cards=scaffold_cards,
        scaffold_dispositions=scaffold_dispositions,
        extension_by_seq=extension_by_seq,
        telomere_unitig_inventory=telomere_unitig_inventory,
    )
    print(full_report)

    if args.output_md:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_md)), exist_ok=True)
        with open(args.output_md, "w") as out:
            out.write(full_report)
        print(f"\n[SUCCESS] Wrote curation reconnaissance blueprint to: {args.output_md}")

    if args.output_json:
        manifest_data = build_manifest_json(
            sample_id=args.sample,
            script_fingerprint=get_script_fingerprint(),
            nuclear_chroms=nuclear_chroms,
            unplaced_chroms=unplaced_chroms,
            grs_threshold=grs_threshold,
            karyotype_stats=karyotype_stats,
            policy=policy,
            internal_cards=internal_cards,
            telo_cards=telo_cards,
            rdna_cards=rdna_cards,
            scaffold_cards=scaffold_cards,
            scaffold_dispositions=scaffold_dispositions,
            telomere_unitig_inventory=telomere_unitig_inventory,
        )
        os.makedirs(os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True)
        with open(args.output_json, "w") as out:
            json.dump(manifest_data, out, indent=2, sort_keys=True)
        print(f"[SUCCESS] Wrote machine-executable curation manifest to: {args.output_json}")


if __name__ == "__main__":
    main()
