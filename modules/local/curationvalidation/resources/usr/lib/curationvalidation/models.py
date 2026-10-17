from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List


@dataclass
class EditedRegion:
    """Represents a surgically modified region in the curated assembly."""
    sample: str
    chrom: str
    card_id: str
    action: str
    edit_type: str
    curated_start: int
    curated_end: int
    orig_start: int
    orig_end: int
    delta_bp: int
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class StructuralVariant:
    """Represents a structural variant parsed from Sniffles VCF."""
    sv_id: str
    chrom: str
    start: int
    end: int
    sv_type: str
    sv_len: int
    support: int
    vaf: float
    filter_status: str
    qual: float
    chr2: str
    info: Dict[str, Any] = field(default_factory=dict)


@dataclass
class EditValidationResult:
    """Validation outcome for a specific edited locus."""
    region: EditedRegion
    mean_coverage: float
    spanning_reads: int
    conflicting_svs: List[StructuralVariant] = field(default_factory=list)
    status: str = "PASS"  # PASS, WARN, FAIL
    explanation: str = ""


@dataclass
class FlagrantError:
    """Anomalous structural defect or severe coverage lesion across the assembly."""
    chrom: str
    start: int
    end: int
    error_type: str  # HOM_SV_DEFECT, COVERAGE_DROPOUT, CHIMERIC_TRANSLOCATION
    severity: str    # CRITICAL, HIGH, MEDIUM
    details: Dict[str, Any] = field(default_factory=dict)
    description: str = ""
