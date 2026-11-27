# Phased telomere evidence for gfapdefinisher

`nf-core/pdefinisher` writes `evidence/phased_end_evidence.json` before running `gfapdefinisher finish`. It passes `--end-evidence` when that option appears in the installed executable's help. Older installations continue to emit graph candidates, and the workflow attaches the family matches during reconciliation.

## JSON contract version 1

| Field | Content |
| --- | --- |
| `version` | Integer `1`. Reject unsupported versions. |
| `coordinate_system` | `zero-based half-open on forward sequences`. All reference, graph, and read coordinates are on the original forward record. |
| `graph_sha256` | Digest of the exact GFA input bytes, keyed by `primary` and, when supplied, `raw`. Unitig IDs are local to these graphs. |
| `families` | Array of records keyed by `cluster_id`. |
| `families[].members` | Array with `end_id`, `molecule_id`, `read_end` (`5p` or `3p`), `tract_start` and `tract_end`. A molecule can contribute more than one observed end. |
| `families[].haplotype_patterns` | Array with `pattern_id`, `read_count`, and the archive’s `alleles` string. These are local aggregate patterns. Version 1 carries no assignments of individual molecules to patterns. |
| `consensus_placements` | Oriented placements keyed by `cluster_id`, `source` (`primary` or `raw`) and `record`. Includes `strand`, `start`, `end`, `target_bp`, `consensus_start`, `consensus_end`, `matches`, `alignment_bp`, `mapq`, `family_molecules`, and `patterns` (comma-separated IDs). Multiple placements for a family can remain. |

Validate the primary digest before using its placements. Use raw placements only when the exact raw graph is available. Keep family matches as candidate links, and retain competing placements. For sequence closure, require evidence that a molecule connects the uniquely anchored path through to the claimed continuation. A family’s molecule count describes the archive and does not count reads spanning a particular graph junction. Include the relevant family and molecule IDs in the evidence behind a candidate when available.
