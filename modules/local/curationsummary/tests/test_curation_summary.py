#!/usr/bin/env python3
"""Unit tests for curationsummary cohort aggregation and elimination range propagation."""

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(BIN_DIR))

from pde_curation_summary import (
    load_scorecard,
    parse_sample_record,
    compute_cohort_aggregates,
    generate_markdown_report,
    generate_tsv_table,
)


class TestCurationSummary(unittest.TestCase):

    def setUp(self):
        self.sample_sc_with_elim = {
            "sample": "sample_A",
            "t2t_completeness": {
                "total_chromosome_termini": 12,
                "pre_capped_termini": 10,
                "post_capped_termini": 12,
                "pre_t2t_score_pct": 83.3,
                "post_t2t_score_pct": 100.0,
                "uncapped_ends": [],
                "net_telomeric_bp_added": 4500,
                "terminus_status_breakdown": {},
            },
            "gap_resolution": {
                "baseline_gaps": 2,
                "remaining_gaps": 0,
                "rdna_gaps": 0,
                "remaining_non_rdna_gaps": 0,
                "is_gapless_non_rdna": True,
                "gaps_excised_count": 1,
                "gaps_patched_count": 1,
                "total_patched_bp": 1200,
                "total_excised_bp": 200,
            },
            "contiguity": {
                "pre_curation": {"total_bp": 100000000, "scaffold_count": 10, "n50": 15000000, "chr_n50": 16000000, "total_n_bp": 200, "total_gaps": 2},
                "post_curation": {"total_bp": 100001000, "scaffold_count": 6, "n50": 16000000, "chr_n50": 16000000, "total_n_bp": 0, "total_gaps": 0},
                "deltas": {"length_delta_bp": 1000, "rdna_length_delta_bp": 0, "non_rdna_length_delta_bp": 1000, "scaffolds_delta": -4},
            },
            "karyotype_compartments": {
                "classification_threshold": 25.0,
                "putative_core": {
                    "span_bp": 80000000,
                    "span_mb": 80.0,
                    "span_pct": 80.0,
                    "median_coverage": 65.0,
                },
                "putative_grs": {
                    "span_bp": 20000000,
                    "span_mb": 20.0,
                    "span_pct": 20.0,
                    "median_coverage": 12.0,
                },
                "elimination_range": {
                    "canon_grs_bp": 18000000,
                    "unplaced_grs_bp": 1500000,
                    "unplaced_uncertain_bp": 500000,
                    "unplaced_core_bp": 200000,
                    "min_elim_bp": 19500000,
                    "max_elim_bp": 20000000,
                    "est_elim_bp": 19750000,
                },
            },
            "reconciliation": {
                "measured_delta_bp": 1000,
                "rdna_delta_bp": 0,
                "non_rdna_delta_bp": 1000,
                "mathematically_conserved": True,
                "discrepancy_bp": 0,
            },
        }

    def test_load_scorecard_with_elimination_range(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sc_path = os.path.join(tmpdir, "sample_A.curation_scorecard.json")
            with open(sc_path, "w") as f:
                json.dump(self.sample_sc_with_elim, f)

            sc = load_scorecard(sc_path)
            self.assertIsNotNone(sc)
            sample_data = parse_sample_record(sc)
            self.assertEqual(sample_data["sample"], "sample_A")
            self.assertEqual(sample_data["canon_grs_bp"], 18000000)
            self.assertEqual(sample_data["unplaced_grs_bp"], 1500000)
            self.assertEqual(sample_data["unplaced_uncertain_bp"], 500000)
            self.assertEqual(sample_data["unplaced_core_bp"], 200000)
            self.assertEqual(sample_data["min_elim_bp"], 19500000)
            self.assertEqual(sample_data["max_elim_bp"], 20000000)
            self.assertEqual(sample_data["est_elim_bp"], 19750000)
            self.assertEqual(sample_data["min_elim_pct"], 19.5)
            self.assertEqual(sample_data["max_elim_pct"], 20.0)
            self.assertEqual(sample_data["est_elim_pct"], 19.75)

    def test_load_scorecard_legacy_fallback(self):
        sc_legacy = json.loads(json.dumps(self.sample_sc_with_elim))
        del sc_legacy["karyotype_compartments"]["elimination_range"]

        with tempfile.TemporaryDirectory() as tmpdir:
            sc_path = os.path.join(tmpdir, "sample_legacy.curation_scorecard.json")
            with open(sc_path, "w") as f:
                json.dump(sc_legacy, f)

            sc = load_scorecard(sc_path)
            sample_data = parse_sample_record(sc)
            self.assertIsNotNone(sample_data)
            self.assertEqual(sample_data["min_elim_bp"], 20000000)
            self.assertEqual(sample_data["max_elim_bp"], 20000000)
            self.assertEqual(sample_data["est_elim_bp"], 20000000)
            self.assertEqual(sample_data["canon_grs_bp"], 20000000)
            self.assertEqual(sample_data["unplaced_grs_bp"], 0)
            self.assertEqual(sample_data["unplaced_uncertain_bp"], 0)

    def test_compute_cohort_aggregates_and_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sc_path = os.path.join(tmpdir, "sample_A.curation_scorecard.json")
            with open(sc_path, "w") as f:
                json.dump(self.sample_sc_with_elim, f)

            sc = load_scorecard(sc_path)
            s_a = parse_sample_record(sc)
            samples = [s_a]
            agg = compute_cohort_aggregates(samples)

            self.assertEqual(agg["tot_canon_grs_bp"], 18000000)
            self.assertEqual(agg["tot_unplaced_grs_bp"], 1500000)
            self.assertEqual(agg["tot_unplaced_uncertain_bp"], 500000)
            self.assertEqual(agg["tot_min_elim_bp"], 19500000)
            self.assertEqual(agg["tot_max_elim_bp"], 20000000)
            self.assertEqual(agg["tot_est_elim_bp"], 19750000)
            self.assertEqual(agg["cohort_min_elim_pct"], 19.5)
            self.assertEqual(agg["cohort_max_elim_pct"], 20.0)

            # TSV table checks
            tsv = generate_tsv_table(samples)
            headers = tsv.splitlines()[0].split("\t")
            for col in [
                "canon_grs_bp",
                "unplaced_grs_bp",
                "unplaced_uncertain_bp",
                "unplaced_core_bp",
                "min_elim_bp",
                "max_elim_bp",
                "min_elim_pct",
                "max_elim_pct",
                "est_elim_bp",
                "est_elim_pct",
            ]:
                self.assertIn(col, headers)

            row = tsv.splitlines()[1].split("\t")
            min_idx = headers.index("min_elim_bp")
            max_idx = headers.index("max_elim_bp")
            self.assertEqual(row[min_idx], "19500000")
            self.assertEqual(row[max_idx], "20000000")

            # Markdown report checks
            md = generate_markdown_report("TEST_COHORT", agg, samples)
            self.assertIn("19.50 – 20.00 Mb (19.5%–20.0%)", md)
            self.assertIn("## Programmed DNA Elimination (PDE) Compartment Allocation", md)
            self.assertIn("| Sample | Threshold ($T_{\\text{GRS}}$) | Core Span (bp) | Core Cov | Canon GRS (bp) | Unplaced GRS (bp) | Uncertain (bp) | Min Elim (bp) | Max Elim (bp) | Est Elim (bp) | Elim Range (%) |", md)
            self.assertIn("19.5%–20.0%", md)


if __name__ == "__main__":
    unittest.main()
