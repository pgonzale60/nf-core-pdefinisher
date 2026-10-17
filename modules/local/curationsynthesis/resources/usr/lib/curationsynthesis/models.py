"""Core domain models for the Open Port Graph & Valency Conservation philosophy.

Every unfinished boundary in a draft assembly is represented as an Open Port (`Port`):
  - Uncapped chromosome termini (`TERM_5P`, `TERM_3P`)
  - Internal gap flanks (`GAP_LEFT`, `GAP_RIGHT`)
  - Unplaced segment perimeter ends (`CONTIG_5P`, `CONTIG_3P`)

Every curation operation is a deterministic port-pairing or port-capping transformation
governed by global Anchor Valency Conservation (in-degree m, out-degree n) and the
Principle of Exclusion.

British spelling standard enforced across all documentation and comments.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class Port:
    """An open boundary port on a nuclear chromosome or unplaced segment."""

    port_id: str
    molecule: str
    coord: int
    orientation: str  # '5p' / 'left' or '3p' / 'right'
    port_type: str    # 'TERM_5P', 'TERM_3P', 'GAP_LEFT', 'GAP_RIGHT', 'CONTIG_5P', 'CONTIG_3P'
    is_capped: bool = False
    is_germline_tip: bool = False
    is_soma_boundary: bool = False
    d_soma: int = 0
    repeat_family: Optional[str] = None
    target_chroms: Set[str] = field(default_factory=set)
    valency_in: int = 0
    valency_out: int = 0

    @property
    def is_monovalent(self) -> bool:
        """True when this port's repeat or sequence signature maps to a single locus."""
        return self.valency_in <= 1 and self.valency_out <= 1 and len(self.target_chroms) <= 1

    @property
    def is_chrom_exclusive(self) -> bool:
        """True when all alignments of this port are confined to a single nuclear chromosome."""
        return len(self.target_chroms) == 1


@dataclass
class UnplacedSegment:
    """An unplaced contig or scaffold evaluated for redundancy triage and port integration."""

    scaffold_id: str
    length: int
    chrom_aligned_bp: int = 0
    chrom_aligned_pct: float = 0.0
    rdna_operon_bp: int = 0
    rdna_operon_pct: float = 0.0
    novel_bp: int = 0
    novel_pct: float = 0.0
    repeat_bp: int = 0
    repeat_pct: float = 0.0
    triage_category: str = "Unknown"
    recommended_action: str = "REPORT_ONLY"
    port_5p: Optional[Port] = None
    port_3p: Optional[Port] = None
    primary_records: List[ dict ] = field(default_factory=list)
    any_records: List[ dict ] = field(default_factory=list)
    target_chromosomes: Set[str] = field(default_factory=set)
    terminal_repeat_ports: Dict[str, List[Tuple[str, str, dict]]] = field(default_factory=dict)
