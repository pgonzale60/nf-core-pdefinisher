"""Module 3 Engine A: Targeted Gap Assembly library (`gapassembly`).

Submodules:
  - `flank_io`: Flank extraction, palindromic anchor reduction, and BAM read baiting
  - `consensus_bridge`: Spanning-read abPOA consensus patching
  - `unitig_valency`: Genome-wide unitig valency conservation & PAF bridge evaluation
"""

from .consensus_bridge import (
    claimed_insertion_bp,
    collect_spanning_read_candidates,
    decide_consensus_patch,
    median_patch_support,
    resolve_spanned_gap,
    run_abpoa,
    trim_flank_residue,
    write_consensus_input,
)
from .flank_io import (
    CONSENSUS_ACTIONS,
    FLANK_RECORDS,
    GAP_ACTIONS,
    bait_reads_for_gap,
    choose_anchor,
    extract_fasta_region,
    gfa_to_fasta,
    longest_inverted_repeat,
    parse_gaps_from_input,
    read_fasta_records,
    read_fastq_records,
    require_indexable_reference,
    revcomp,
    reverse_complement,
    write_flank_reference,
)
from .unitig_valency import (
    evaluate_gap_paf,
    parse_cigar_ops,
    project_target_interval_via_cigar,
    survey_repeat_anchor_valency,
    survey_unitig_valencies_genome_wide,
)

__all__ = [
    "CONSENSUS_ACTIONS",
    "FLANK_RECORDS",
    "GAP_ACTIONS",
    "bait_reads_for_gap",
    "choose_anchor",
    "claimed_insertion_bp",
    "collect_spanning_read_candidates",
    "decide_consensus_patch",
    "evaluate_gap_paf",
    "extract_fasta_region",
    "gfa_to_fasta",
    "longest_inverted_repeat",
    "median_patch_support",
    "parse_cigar_ops",
    "parse_gaps_from_input",
    "project_target_interval_via_cigar",
    "read_fasta_records",
    "read_fastq_records",
    "require_indexable_reference",
    "resolve_spanned_gap",
    "revcomp",
    "reverse_complement",
    "run_abpoa",
    "survey_repeat_anchor_valency",
    "survey_unitig_valencies_genome_wide",
    "trim_flank_residue",
    "write_consensus_input",
    "write_flank_reference",
]
