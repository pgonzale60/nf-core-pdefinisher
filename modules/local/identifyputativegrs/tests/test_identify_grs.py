"""Unit tests for identifyputativegrs helper functions and SSOT policy calibration.

British spelling standard enforced across all documentation and comments.
"""

import os
from pathlib import Path
import sys
import unittest

bin_dir = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(bin_dir))

from pde_identify_grs import (
    build_chromosome_segments,
    cluster_breakpoints_by_strand,
    find_curation_policy,
    fit_em_2comp,
    load_curation_policy,
    run_pushdown_automaton,
)


class TestIdentifyPutativeGRS(unittest.TestCase):

    def test_load_curation_policy(self):
        policy_path = find_curation_policy()
        self.assertIsNotNone(policy_path)
        self.assertTrue(os.path.exists(policy_path))
        policy = load_curation_policy(policy_path)
        self.assertIn("pde_segmentation", policy)
        pde_policy = policy["pde_segmentation"]
        self.assertAlmostEqual(pde_policy["max_grs_to_core_ratio"], 0.40)
        self.assertAlmostEqual(pde_policy["min_transition_contrast"], 2.0)
        self.assertAlmostEqual(pde_policy["default_grs_coverage_ratio"], 0.35)
        self.assertEqual(pde_policy["max_terminal_arm_bp"], 1500000)

    def test_cluster_breakpoints_by_strand(self):
        # 3 breakpoints on '+' strand within 3 bp -> cluster together
        raw_bps = [
            (100, "+", 10.0, {"name": "bp1"}),
            (102, "+", 20.0, {"name": "bp2"}),
            (103, "+", 5.0, {"name": "bp3"}),
            (500, "-", 15.0, {"name": "bp4"}),
            (501, "-", 12.0, {"name": "bp5"}),
        ]
        clustered = cluster_breakpoints_by_strand(raw_bps, max_dist=3, min_support=5.0)
        self.assertEqual(len(clustered), 2)
        # First cluster is '+' at 102 with total support 35.0
        self.assertEqual(clustered[0][0], 102)
        self.assertEqual(clustered[0][1], "+")
        self.assertEqual(clustered[0][2], 35.0)
        # Second cluster is '-' at 500 with total support 27.0
        self.assertEqual(clustered[1][0], 500)
        self.assertEqual(clustered[1][1], "-")
        self.assertEqual(clustered[1][2], 27.0)

    def test_run_pushdown_automaton(self):
        # Synthesize coverage windows for Chr I:
        # 0 to 200 kb: 10x (GRS)
        # 200 kb to 4 Mb: 60x (Core)
        # 4 Mb to 4.02 Mb (20 kb): 8x (Internal GRS)
        # 4.02 Mb to 10 Mb: 60x (Core)
        # 10 Mb to 10.5 Mb: 12x (3' GRS)
        cov_windows = {"chr1": []}
        for pos in range(0, 10_500_000, 10_000):
            if pos < 200_000:
                cov = 10.0
            elif 4_000_000 <= pos < 4_020_000:
                cov = 8.0
            elif pos >= 10_000_000:
                cov = 12.0
            else:
                cov = 60.0
            cov_windows["chr1"].append((pos, pos + 10_000, cov))

        clustered_bps = [
            (200_000, "-", 25.0, {}),
            (4_000_000, "+", 30.0, {}),
            (4_020_000, "-", 30.0, {}),
            (10_000_000, "+", 20.0, {}),
        ]
        t_grs = 25.0
        ok, grs_intervals, frontier_bps = run_pushdown_automaton(
            "chr1",
            10_500_000,
            clustered_bps,
            cov_windows,
            t_grs,
            max_terminal_arm_bp=1_500_000,
            max_internal_grs_bp=500_000,
            min_transition_contrast=2.0,
        )
        self.assertTrue(ok)
        self.assertEqual(len(grs_intervals), 3)
        self.assertEqual(grs_intervals[0], (0, 200_000))
        self.assertEqual(grs_intervals[1], (4_000_000, 4_020_000))
        self.assertEqual(grs_intervals[2], (10_000_000, 10_500_000))

    def test_run_pushdown_automaton_rejects_huge_terminal_arm(self):
        # A '-' breakpoint at 4 Mb (violates max_terminal_arm_bp=1.5 Mb)
        cov_windows = {"chr1": [(p, p + 10_000, 50.0) for p in range(0, 10_000_000, 10_000)]}
        clustered_bps = [(4_000_000, "-", 40.0, {})]
        ok, grs_intervals, _ = run_pushdown_automaton(
            "chr1",
            10_000_000,
            clustered_bps,
            cov_windows,
            25.0,
            max_terminal_arm_bp=1_500_000,
        )
        self.assertFalse(ok)
        self.assertEqual(len(grs_intervals), 0)

    def test_fit_em_2comp_unimodal_rejection(self):
        # Unimodal normal-like distribution around 65x (e.g. 50x to 80x)
        # Should yield component ratio > 0.45, confirming unimodality
        import random
        random.seed(42)
        vals = [max(1.0, random.gauss(65.0, 10.0)) for _ in range(5000)]
        em_res = fit_em_2comp(vals)
        mu_ratio = em_res["mu_linear"][0] / max(em_res["mu_linear"][1], 1e-3)
        # Ratio should be high (> 0.60) since it's a single unimodal distribution
        self.assertGreater(mu_ratio, 0.45)

    def test_fit_em_2comp_massive_elimination_bimodal(self):
        # Massive elimination genome (e.g. APS4): 65% GRS at 12x, 35% Core at 65x
        import random
        random.seed(42)
        grs_vals = [max(1.0, random.gauss(12.0, 3.0)) for _ in range(6500)]
        core_vals = [max(1.0, random.gauss(65.0, 10.0)) for _ in range(3500)]
        vals = grs_vals + core_vals
        em_res = fit_em_2comp(vals)
        mu_ratio = em_res["mu_linear"][0] / max(em_res["mu_linear"][1], 1e-3)
        # Should be bimodal with separation ratio <= 0.45
        self.assertLessEqual(mu_ratio, 0.45)
        # Somatic core (higher component) must be near 65x
        self.assertGreater(em_res["mu_linear"][1], 50.0)

    def test_build_chromosome_segments_low_depth_flanks_labeled_grs(self):
        # Scenario: SUPER_X with 206 kb candidate GRS in the middle, but flanks are all low depth (~37x-41x <= t_grs=89.5x)
        # Flanks MUST NOT be mislabeled putative_core, and contiguous GRS intervals must merge.
        cov_windows = {
            "SUPER_X": [
                (0, 10_000_000, 37.0),
                (10_000_000, 10_200_000, 41.0),
                (10_200_000, 26_000_000, 41.0),
            ]
        }
        all_grs = [(10_000_000, 10_200_000, 41.0)]
        t_grs = 89.5
        segs = build_chromosome_segments("SUPER_X", 26_000_000, all_grs, cov_windows, t_grs)
        self.assertEqual(len(segs), 1)
        self.assertEqual(segs[0][0], "SUPER_X")
        self.assertEqual(segs[0][1], 0)
        self.assertEqual(segs[0][2], 26_000_000)
        self.assertEqual(segs[0][3], "putative_GRS")

    def test_build_chromosome_segments_high_depth_flanks_labeled_core(self):
        # Scenario: Normal chromosome with internal GRS between true core regions
        cov_windows = {
            "chr1": [
                (0, 5_000_000, 240.0),
                (5_000_000, 6_000_000, 30.0),
                (6_000_000, 10_000_000, 240.0),
            ]
        }
        all_grs = [(5_000_000, 6_000_000, 30.0)]
        t_grs = 89.5
        segs = build_chromosome_segments("chr1", 10_000_000, all_grs, cov_windows, t_grs)
        self.assertEqual(len(segs), 3)
        self.assertEqual(segs[0], ("chr1", 0, 5_000_000, "putative_core", 240.0))
        self.assertEqual(segs[1], ("chr1", 5_000_000, 6_000_000, "putative_GRS", 30.0))
        self.assertEqual(segs[2], ("chr1", 6_000_000, 10_000_000, "putative_core", 240.0))


if __name__ == "__main__":
    unittest.main()


