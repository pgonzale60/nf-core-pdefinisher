#!/usr/bin/env python3
"""
Unit test for curationvalidation library and CLI engine.
"""

import os
import tempfile
import unittest
from pathlib import Path

# Add resources/usr/lib to path
base_dir = Path(__file__).resolve().parent.parent
lib_dir = base_dir / "resources" / "usr" / "lib"
import sys
sys.path.insert(0, str(lib_dir))

from curationvalidation import (
    EditedRegion,
    StructuralVariant,
    parse_agp_regions,
    parse_curation_audit,
    validate_edited_regions,
    find_flagrant_errors,
    generate_validation_json,
    render_validation_markdown,
)


class TestCurationValidation(unittest.TestCase):

    def test_clean_patch_validation(self):
        region = EditedRegion(
            sample="test_sample",
            chrom="I",
            card_id="GAP_001",
            action="PATCH_GAP",
            edit_type="GAP_PATCH",
            curated_start=1000,
            curated_end=1500,
            orig_start=1000,
            orig_end=1200,
            delta_bp=300
        )
        # No conflicting SVs
        svs = [
            StructuralVariant(
                sv_id="SV_001",
                chrom="I",
                start=50000,
                end=50200,
                sv_type="DEL",
                sv_len=200,
                support=15,
                vaf=0.5,
                filter_status="PASS",
                qual=60.0,
                chr2="I"
            )
        ]

        results = validate_edited_regions([region], svs, bam_path=None)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].status, "PASS")
        self.assertEqual(len(results[0].conflicting_svs), 0)

    def test_conflicting_sv_warning(self):
        region = EditedRegion(
            sample="test_sample",
            chrom="I",
            card_id="GAP_001",
            action="PATCH_GAP",
            edit_type="GAP_PATCH",
            curated_start=1000,
            curated_end=1500,
            orig_start=1000,
            orig_end=1200,
            delta_bp=300
        )
        # Conflicting deletion right across the patch
        svs = [
            StructuralVariant(
                sv_id="SV_CONFLICT",
                chrom="I",
                start=1100,
                end=1400,
                sv_type="DEL",
                sv_len=300,
                support=12,
                vaf=0.85,
                filter_status="PASS",
                qual=60.0,
                chr2="I"
            )
        ]

        results = validate_edited_regions([region], svs, bam_path=None)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].status, "WARN")
        self.assertEqual(len(results[0].conflicting_svs), 1)

    def test_flagrant_error_detection(self):
        svs = [
            # High VAF large deletion (flagrant defect)
            StructuralVariant(
                sv_id="SV_DEFECT",
                chrom="II",
                start=200000,
                end=205000,
                sv_type="DEL",
                sv_len=5000,
                support=20,
                vaf=0.95,
                filter_status="PASS",
                qual=60.0,
                chr2="II"
            ),
            # Low VAF heterozygous variant (normal somatic/heterozygosity)
            StructuralVariant(
                sv_id="SV_HET",
                chrom="II",
                start=300000,
                end=300100,
                sv_type="INS",
                sv_len=100,
                support=5,
                vaf=0.30,
                filter_status="PASS",
                qual=30.0,
                chr2="II"
            )
        ]

        errors = find_flagrant_errors(svs)
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0].chrom, "II")
        self.assertEqual(errors[0].severity, "CRITICAL")
        self.assertEqual(errors[0].details["sv_type"], "DEL")

    def test_report_and_json_generation(self):
        region = EditedRegion(
            sample="test_sample",
            chrom="I",
            card_id="TELO_001",
            action="CAP_TERMINUS",
            edit_type="TELOMERE_CAPPING_5P",
            curated_start=1,
            curated_end=2500,
            orig_start=0,
            orig_end=0,
            delta_bp=2500
        )
        res = validate_edited_regions([region], [], bam_path=None)
        payload = generate_validation_json("test_sample", res, [], 0)
        self.assertEqual(payload["sample"], "test_sample")
        self.assertEqual(payload["summary"]["edits_passed"], 1)

        md = render_validation_markdown(payload)
        self.assertIn("# Post-Curation Mapping & SV Validation Dossier: `test_sample`", md)
        self.assertIn("TELOMERE_CAPPING_5P", md)

    def test_loaders_parse_agp_and_audit(self):
        with tempfile.TemporaryDirectory() as td:
            agp_path = os.path.join(td, "test.agp")
            with open(agp_path, "w") as f:
                f.write("chr1\t1\t1000\t1\tW\tcontig_1\t1\t1000\t+\n")
                f.write("chr1\t1001\t1100\t2\tN\t100\tscaffold\tyes\tpaired-ends\n")
            
            agp_res = parse_agp_regions(agp_path)
            self.assertIn("chr1", agp_res)
            self.assertEqual(len(agp_res["chr1"]), 2)

            audit_path = os.path.join(td, "test.audit.tsv")
            with open(audit_path, "w") as f:
                f.write("sample\tchrom\tcard_id\taction\tedit_type\torig_start\torig_end\tdelta_bp\n")
                f.write("sample1\tchr1\tCARD-01\tPATCH_GAP\tGAP_PATCH\t1000\t1100\t500\n")

            regions = parse_curation_audit(audit_path, agp_entries=agp_res)
            self.assertEqual(len(regions), 1)
            self.assertEqual(regions[0].chrom, "chr1")
            self.assertEqual(regions[0].card_id, "CARD-01")


if __name__ == "__main__":
    unittest.main()
