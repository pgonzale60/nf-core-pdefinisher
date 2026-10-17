"""Package init for curationvalidation module."""

from .models import EditedRegion, StructuralVariant, EditValidationResult, FlagrantError
from .policy import DEFAULT_VALIDATION_POLICY, __version__, get_script_fingerprint
from .loaders import parse_agp_regions, parse_curation_audit, parse_sniffles_vcf
from .verifier import validate_edited_regions, find_flagrant_errors
from .reporting import generate_validation_json, render_validation_markdown

__all__ = [
    "EditedRegion",
    "StructuralVariant",
    "EditValidationResult",
    "FlagrantError",
    "DEFAULT_VALIDATION_POLICY",
    "__version__",
    "get_script_fingerprint",
    "parse_agp_regions",
    "parse_curation_audit",
    "parse_sniffles_vcf",
    "validate_edited_regions",
    "find_flagrant_errors",
    "generate_validation_json",
    "render_validation_markdown",
]
