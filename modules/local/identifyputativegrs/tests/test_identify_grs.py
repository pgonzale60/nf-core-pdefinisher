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
    build_chromosomal_segments,
    classify_unplaced,
    cluster_breakpoints_by_strand,
    delfies_bounded_grs,
    drop_overlap_bp,
    find_curation_policy,
    fit_em_2comp,
    load_curation_policy,
    quantize_drops,
)

RULES = {
    "drop_overlap_fraction": 0.15,
    "micro_grs_max_bp": 1000,
    "anchor_proximity_bp": 2000,
    "min_closed_grs_bp": 1000,
    "unplaced_min_nonrepeat_bp": 1000,
    "unplaced_grs_drop_fraction": 0.5,
}


def bp(pos, strand, end=None):
    return (pos, strand, 10.0, {"end": pos + 1 if end is None else end})


def windows(chrom, length, drops, core=60.0, low=5.0, step=10_000):
    out = []
    for p in range(0, length, step):
        in_drop = any(ds <= p and p + step <= de for ds, de in drops)
        out.append((p, p + step, low if in_drop else core))
    return {chrom: out}


def grs(length, bps, drops, chrom="chr1"):
    cov = windows(chrom, length, drops)
    raw, merged = quantize_drops(cov, chrom, 20.0, 1000)
    return delfies_bounded_grs(length, bps, raw, merged, RULES)


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
        self.assertIsNone(pde_policy.get("max_terminal_arm_bp"))
        self.assertIsNone(pde_policy.get("max_internal_grs_bp"))

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


    def test_policy_defines_delfies_bounded_rules(self):
        pde = load_curation_policy(find_curation_policy())["pde_segmentation"]
        for key in RULES:
            self.assertIn(key, pde)

    def test_quantize_and_overlap(self):
        cov = windows("chr1", 100_000, [(20_000, 40_000), (40_500, 50_000)], step=500)
        raw, merged = quantize_drops(cov, "chr1", 20.0, 1000)
        self.assertEqual(merged, [(20_000, 50_000)])
        self.assertEqual(drop_overlap_bp(raw, 0, 30_000), 10_000)

    def test_low_coverage_repeat_without_breakpoint_is_core(self):
        # MAPQ-filtered somatic repeat shows as a drop but has no DELFIES anchor
        out = grs(1_000_000, [], [(300_000, 500_000)])
        self.assertEqual(out, [])
        segs = build_chromosomal_segments("chr1", 1_000_000, out, windows("chr1", 1_000_000, []))
        self.assertEqual([s[3] for s in segs], ["putative_core"])

    def test_paired_macro_grs_requires_drop_support(self):
        bps = [bp(200_000, "+"), bp(400_000, "-")]
        supported = grs(1_000_000, bps, [(200_000, 400_000)])
        self.assertEqual([(s, e, t) for s, e, t, _ in supported], [(200_000, 400_001, "Macro_GRS_Paired")])
        self.assertEqual(grs(1_000_000, bps, []), [])

    def test_micro_grs_with_reciprocal_cut_needs_no_drop(self):
        out = grs(1_000_000, [bp(500_000, "+", 500_001), bp(500_300, "-", 500_301)], [])
        self.assertEqual(out[0][2], "Micro_GRS_Paired")

    def test_distal_grs_5p_and_3p(self):
        out = grs(1_000_000, [bp(100_000, "-"), bp(900_000, "+")], [(0, 100_000), (900_000, 1_000_000)])
        self.assertEqual([t for _, _, t, _ in out], ["5_prime_distal_GRS", "3_prime_distal_GRS"])
        self.assertEqual((out[0][0], out[1][1]), (0, 1_000_000))

    def test_orphan_plus_closed_by_coverage(self):
        out = grs(1_000_000, [bp(300_000, "+")], [(300_000, 380_000)])
        self.assertEqual(out[0][:3], (300_000, 380_000, "Start_DELFIES_Anchored_End_Coverage_Closed"))

    def test_orphan_minus_opened_by_coverage(self):
        out = grs(1_000_000, [bp(100_000, "+"), bp(920_000, "-", 920_001)], [(840_000, 920_000)])
        # the pair spans 820 kb with <15 % drop, so it is unsupported and the - opens via coverage
        self.assertIn("Start_Coverage_Opened_End_DELFIES_Anchored", [t for _, _, t, _ in out])

    def test_nested_variable_grs_reports_outer_and_inner(self):
        bps = [bp(200_000, "+"), bp(210_000, "+"), bp(300_000, "-", 300_001), bp(320_000, "-", 320_001)]
        out = grs(1_000_000, bps, [(200_000, 320_000)])
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0][:3], (200_000, 320_001, "Nested_Variable_GRS_Outer"))
        self.assertEqual(out[0][3], "inner=210000-300001")

    def test_chromosomal_tiling_flanks_are_always_core(self):
        cov = windows("chr1", 1_000_000, [(0, 1_000_000)])  # whole chromosome low depth
        segs = build_chromosomal_segments("chr1", 1_000_000, [(300_000, 400_000, "Macro_GRS_Paired", "")], cov)
        self.assertEqual([s[3] for s in segs], ["putative_core", "putative_GRS", "putative_core"])

    def test_unplaced_repeat_aware_classes(self):
        drops_all = [(0, 50_000)]
        # mostly repeat: < 1 kb unique sequence -> uncertain, whatever the coverage
        self.assertEqual(classify_unplaced("u1", 50_000, drops_all, {"u1": [(0, 49_500)]}, RULES)[0],
                         "putative_GRS_uncertain")
        # enough unique sequence: deterministic on drop fraction
        self.assertEqual(classify_unplaced("u1", 50_000, drops_all, {}, RULES)[0], "putative_GRS")
        self.assertEqual(classify_unplaced("u1", 50_000, [(0, 10_000)], {}, RULES)[0], "putative_core")

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



    def test_chromosomes_list_filtering(self):
        import subprocess, tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = os.path.join(tmpdir, "a.fa")
            with open(fa, "w") as f:
                f.write(">chrA\n" + "A" * 10 + "\n>chrB\n" + "A" * 10 + "\n>u1\n" + "A" * 10 + "\n")
            with open(fa + ".fai", "w") as f:
                f.write("chrA\t2000000\t0\t0\t0\nchrB\t500000\t0\t0\t0\nu1\t2000000\t0\t0\t0\n")
            cov = os.path.join(tmpdir, "cov.bed")
            with open(cov, "w") as f:
                for p in range(0, 2000000, 10000):
                    f.write(f"chrA\t{p}\t{p+10000}\t60\n")
                for p in range(0, 500000, 10000):
                    d = 4 if p < 100000 else 60
                    f.write(f"chrB\t{p}\t{p+10000}\t{d}\n")
                for p in range(0, 2000000, 10000):
                    f.write(f"u1\t{p}\t{p+10000}\t4\n")
            bp = os.path.join(tmpdir, "bp.bed")
            with open(bp, "w") as f:
                f.write("chrB\t100000\t100001\tb1\t12\t-\n")
            red = os.path.join(tmpdir, "red.bed")
            with open(red, "w") as f:
                f.write("u1\t0\t1999500\n")

            script = str(bin_dir / "pde_identify_grs.py")
            cmd = [
                "python3", script,
                "--assembly", fa,
                "--breakpoints-bed", bp,
                "--regions-bed", cov,
                "--red-bed", red,
                "--chromosomes", "chrA,chrB",
                "--outdir", tmpdir,
                "--prefix", "test"
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"Error: {res.stderr}")
            bed_out = os.path.join(tmpdir, "test.pde_segments.bed")
            with open(bed_out) as bf:
                rows = [l.split("\t") for l in bf.read().strip().split("\n")]
            # chrB (500 kb) was recognized as chromosomal despite being < 1 Mb
            self.assertTrue(any(r[0] == "chrB" and r[3] == "putative_GRS" for r in rows))
            # u1 (2 Mb) was recognized as unplaced despite being >= 1 Mb, labeled putative_GRS_uncertain
            u1_rows = [r for r in rows if r[0] == "u1"]
            self.assertEqual(u1_rows[0][3], "putative_GRS_uncertain")

if __name__ == "__main__":
    unittest.main()
