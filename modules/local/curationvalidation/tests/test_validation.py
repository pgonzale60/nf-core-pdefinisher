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

    def test_purged_sequence_validation(self):
        purged_region = EditedRegion(
            sample="test_sample",
            chrom="scaffold_MT",
            card_id="SCAFF-MT",
            action="PURGE_ORGANELLAR",
            edit_type="PURGE_ORGANELLAR",
            curated_start=0,
            curated_end=15000,
            orig_start=0,
            orig_end=15000,
            delta_bp=-15000
        )
        res = validate_edited_regions([purged_region], [], bam_path=None)
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0].status, "PASS")
        self.assertIn("Purged/retired sequence", res[0].explanation)

        # Test RETIRE_ABSORBED action
        absorbed_region = EditedRegion(
            sample="test_sample",
            chrom="I_unloc_1",
            card_id="SCAFF-I_unloc_1",
            action="RETIRE_ABSORBED",
            edit_type="RETIRE_ABSORBED",
            curated_start=0,
            curated_end=8500,
            orig_start=0,
            orig_end=8500,
            delta_bp=-8500
        )
        res_abs = validate_edited_regions([absorbed_region], [], bam_path=None)
        self.assertEqual(len(res_abs), 1)
        self.assertEqual(res_abs[0].status, "PASS")
        self.assertIn("Purged/retired sequence verified absent", res_abs[0].explanation)

    def test_compensating_tandem_repeat_shift(self):
        # DEL 58 bp at X:3841 and INS 59 bp at X:3688 (distance 153 bp, diff 1 bp <= 5 bp)
        svs = [
            StructuralVariant(
                sv_id="SV_DEL_58",
                chrom="X",
                start=3841,
                end=3899,
                sv_type="DEL",
                sv_len=58,
                support=10,
                vaf=0.833,
                filter_status="PASS",
                qual=60.0,
                chr2="X"
            ),
            StructuralVariant(
                sv_id="SV_INS_59",
                chrom="X",
                start=3688,
                end=3688,
                sv_type="INS",
                sv_len=59,
                support=9,
                vaf=0.692,
                filter_status="PASS",
                qual=60.0,
                chr2="X"
            )
        ]
        errors = find_flagrant_errors(svs)
        self.assertEqual(len(errors), 2)
        for err in errors:
            self.assertEqual(err.error_type, "TANDEM_REPEAT_ALIGNMENT_SHIFT")
            self.assertEqual(err.severity, "LOW")

    def test_agp_coordinate_projection(self):
        with tempfile.TemporaryDirectory() as td:
            agp_path = os.path.join(td, "test.agp")
            with open(agp_path, "w") as f:
                # 5' telomere extension of 5000 bp, then draft contig slice of 50000 bp
                f.write("chr1\t1\t5000\t1\tW\tchr1_5p_telo_ext\t1\t5000\t+\n")
                f.write("chr1\t5001\t55000\t2\tW\tcontig_1\t1\t50000\t+\n")
                f.write("chr1\t55001\t58000\t3\tW\tchr1_3p_telo_ext\t1\t3000\t+\n")

            agp_res = parse_agp_regions(agp_path)

            audit_path = os.path.join(td, "test.audit.tsv")
            with open(audit_path, "w") as f:
                f.write("sample\tchrom\tcard_id\taction\tedit_type\torig_start\torig_end\tdelta_bp\n")
                f.write("sample1\tchr1\tTELO-chr1-5p\tEXTEND_TELOMERE\tTELO_5P\t0\t0\t5000\n")
                f.write("sample1\tchr1\tTELO-chr1-3p\tEXTEND_TELOMERE\tTELO_3P\t50000\t50000\t3000\n")

            regions = parse_curation_audit(audit_path, agp_entries=agp_res)
            self.assertEqual(len(regions), 2)
            # 5' extension boundary: [5000 - 250, 5000 + 250] = [4750, 5250]
            self.assertEqual(regions[0].curated_start, 4750)
            self.assertEqual(regions[0].curated_end, 5250)
            # 3' extension boundary: [55001 - 250, 55001 + 250] = [54751, 55251]
            self.assertEqual(regions[1].curated_start, 54751)
            self.assertEqual(regions[1].curated_end, 55251)

    def test_mutually_exclusive_flagrant_errors(self):
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
        svs = [
            # SV 1: Inside edited locus (+/- 500bp flank) -> belongs to Section 2
            StructuralVariant(
                sv_id="SV_INSIDE_EDIT",
                chrom="I",
                start=1200,
                end=1300,
                sv_type="DEL",
                sv_len=100,
                support=20,
                vaf=0.95,
                filter_status="PASS",
                qual=60.0,
                chr2="I"
            ),
            # SV 2: Far outside edited locus -> belongs to Section 3
            StructuralVariant(
                sv_id="SV_OUTSIDE_EDIT",
                chrom="I",
                start=50000,
                end=50200,
                sv_type="DEL",
                sv_len=200,
                support=20,
                vaf=0.95,
                filter_status="PASS",
                qual=60.0,
                chr2="I"
            )
        ]
        # When regions is provided, SV_INSIDE_EDIT should be excluded from Section 3
        errors = find_flagrant_errors(svs, regions=[region])
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0].details["sv_id"], "SV_OUTSIDE_EDIT")

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
