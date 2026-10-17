"""Satellite repeat scaffolding sub-package: bridges satellite-terminated chromosome ends to telomeres."""

from .repeat_scaffolder import find_sequence_overlap, materialise_satellite_scaffold

__all__ = ["find_sequence_overlap", "materialise_satellite_scaffold"]
