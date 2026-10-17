#!/usr/bin/env python3
"""
pde_validate_curation.py - Post-Curation Mapping & Sniffles Validation Engine

Evaluates physical read mapping and Sniffles SV calls against the curated assembly,
focusing on:
  1. "Edited" regions (gap patches, telomere extensions, rDNA insertions, unplaced severings)
  2. Flagrant whole-genome assembly defects (homozygous inversions, large deletions, translocations)

Outputs:
  - ${sample}.curation_validation.json
  - ${sample}.curation_validation.md
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Add resources/usr/lib to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

from curationvalidation import (
    __version__,
    get_script_fingerprint,
    parse_agp_regions,
    parse_curation_audit,
    parse_sniffles_vcf,
    validate_edited_regions,
    find_flagrant_errors,
    generate_validation_json,
    render_validation_markdown,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-Curation Mapping and Sniffles Structural Variant Validation Engine"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__} (sha256:{get_script_fingerprint()})",
    )
    parser.add_argument(
        "--audit",
        required=True,
        help="Curation audit TSV from CURATIONAGPGENERATOR (*.curation_audit.tsv)",
    )
    parser.add_argument(
        "--agp",
        required=False,
        default=None,
        help="Curated AGP file (*.curated.agp)",
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Sniffles VCF called against curated assembly (*.vcf.gz or *.vcf)",
    )
    parser.add_argument(
        "--bam",
        required=False,
        default=None,
        help="HiFi alignments mapped to curated assembly (*.bam)",
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Parse AGP & Audit
    agp_entries = parse_agp_regions(args.agp) if args.agp else {}
    regions = parse_curation_audit(args.audit, agp_entries)

    # 2. Parse Sniffles VCF
    svs = parse_sniffles_vcf(args.vcf)

    # 3. Validate Edited Regions
    edit_results = validate_edited_regions(
        regions=regions,
        svs=svs,
        bam_path=args.bam,
    )

    # 4. Discover Flagrant Assembly Errors
    flagrant_errors = find_flagrant_errors(svs)

    # 5. Build JSON Payload
    json_payload = generate_validation_json(
        sample=args.sample,
        edit_results=edit_results,
        flagrant_errors=flagrant_errors,
        total_svs_called=len(svs),
    )

    json_path = os.path.join(args.outdir, f"{args.sample}.curation_validation.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_payload, f, indent=2)
    print(f"[{args.sample}] Emitted validation JSON: {json_path}")

    # 6. Render Markdown Dossier
    md_content = render_validation_markdown(json_payload)
    md_path = os.path.join(args.outdir, f"{args.sample}.curation_validation.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)
    print(f"[{args.sample}] Emitted validation Markdown: {md_path}")


if __name__ == "__main__":
    main()
