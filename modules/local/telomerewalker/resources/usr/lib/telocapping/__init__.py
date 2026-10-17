"""Telomeric capping sub-package: unitig materialisation and single-read overhang capping."""

from .overhang_capper import (
    extract_read_by_qname,
    inspect_terminal_reads,
    maximal_tip_overhang,
    named_read_overhang,
    outward_tail,
    overhang_of_read,
    terminal_read_candidates,
    terminal_window_lines,
)
from .unitig_capper import (
    authorized_unitig_materialisation,
    extract_unitig_by_id,
    outboard_extension,
    region_extension,
)

__all__ = [
    "outboard_extension",
    "region_extension",
    "extract_unitig_by_id",
    "authorized_unitig_materialisation",
    "overhang_of_read",
    "terminal_window_lines",
    "maximal_tip_overhang",
    "named_read_overhang",
    "terminal_read_candidates",
    "inspect_terminal_reads",
    "outward_tail",
    "extract_read_by_qname",
]
