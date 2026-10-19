"""Module 2 Curation Readiness & Synthesis Engine package (`curationsynthesis`).

Organised around the Open Port Graph & Valency Conservation philosophy:
  - `policy`: Single Source of Truth policy loading & deterministic package fingerprinting
  - `models`: Core `Port` and `UnplacedSegment` domain models
  - `evidence_io`: Pure I/O loaders and interval utilities
  - `port_valency`: Molecule deduplication, scarcity-aware unitig arbitration & open port graph
  - `internal_resolver`: Internal gap, BND/INV misjoin, macro SV & rDNA array card synthesis
  - `terminus_resolver`: 7-level chromosome terminus resolution hierarchy & chimera pruning
  - `unplaced_triage`: Redundancy triage & chromosomal port extension candidacy
  - `reporting`: Markdown dossier and JSON manifest rendering
"""

from .evidence_io import (
    compute_query_represented_intervals,
    load_delfies_points,
    load_flye_info,
    load_grs_summary,
    load_mosdepth_tip_depths,
    load_rdna_bed,
    load_red_repeats,
    load_telo_junctions,
    load_terminal_softclips,
    load_terminal_summary,
    merge_intervals,
    merge_query_spans,
    parse_bnd_target,
    parse_telo_paf,
    parse_unplaced_paf,
    partition_chromosomes,
    reconstruct_rdna_operons,
    survey_telomere_unitigs,
)
from .internal_resolver import synthesize_curation_cards, synthesize_rdna_cards
from .models import Port, UnplacedSegment
from .policy import DEFAULT_POLICY, __version__, get_script_fingerprint, load_policy
from .port_valency import (
    assign_unitigs_to_termini,
    build_molecule_ids,
    build_open_port_graph,
    build_unitig_assignment,
    candidate_preference,
    eligible_unitig_claims,
    match_open_termini_to_scaffolds,
    match_unplaced_ports_by_valency,
    terminus_state,
)
from .reporting import build_manifest_json, render_markdown_blueprint
from .terminus_resolver import detect_chimeric_terminus, synthesize_telomere_cards
from .unplaced_triage import evaluate_extension_candidacy, triage_unplaced_scaffolds

__all__ = [
    "DEFAULT_POLICY",
    "Port",
    "UnplacedSegment",
    "__version__",
    "assign_unitigs_to_termini",
    "build_manifest_json",
    "build_molecule_ids",
    "build_open_port_graph",
    "build_unitig_assignment",
    "candidate_preference",
    "compute_query_represented_intervals",
    "detect_chimeric_terminus",
    "eligible_unitig_claims",
    "evaluate_extension_candidacy",
    "get_script_fingerprint",
    "load_delfies_points",
    "load_flye_info",
    "load_grs_summary",
    "load_mosdepth_tip_depths",
    "load_policy",
    "load_rdna_bed",
    "load_red_repeats",
    "load_telo_junctions",
    "load_terminal_softclips",
    "load_terminal_summary",
    "match_open_termini_to_scaffolds",
    "match_unplaced_ports_by_valency",
    "merge_intervals",
    "merge_query_spans",
    "parse_bnd_target",
    "parse_telo_paf",
    "parse_unplaced_paf",
    "partition_chromosomes",
    "reconstruct_rdna_operons",
    "render_markdown_blueprint",
    "synthesize_curation_cards",
    "synthesize_rdna_cards",
    "synthesize_telomere_cards",
    "survey_telomere_unitigs",
    "terminus_state",
    "triage_unplaced_scaffolds",
]

