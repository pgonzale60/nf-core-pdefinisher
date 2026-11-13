"""Unit tests for pde_carve_core_islands: core-island carving from GRS intervals.

British spelling standard enforced across all documentation and comments.
"""

from pathlib import Path
import sys
import unittest

bin_dir = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(bin_dir))

from pde_carve_core_islands import carve_row

POLICY = {
    "core_subtraction_min_grs_bp": 30000,
    "core_subtraction_switch_penalty_bp": 5000,
    "core_subtraction_masked_repeat_fraction": 0.9,
}
BASELINE = 30.0
SUBTYPE = "delfies_bounded"
GRS = ("chr1", 0, 100000, "putative_GRS", 12.0, SUBTYPE)


def quant_of(runs):
    """Quantized runs [(s, e, is_low)] in the structure read_quantized builds."""
    runs = sorted(runs)
    return ([r[0] for r in runs], runs)


def red_of(ivs):
    """RED intervals in the structure read_red builds (already merged, sorted)."""
    ivs = sorted(ivs)
    return ([s for s, _ in ivs], ivs)


def carve(row, runs, gaps=None, red=None):
    quant = {"chr1": quant_of(runs)}
    gap_map = {"chr1": sorted(gaps or [])}
    red_map = {"chr1": red_of(red)} if red is not None else {}
    return carve_row(row, quant, gap_map, red_map, BASELINE, POLICY)


class CarveCoreIslandsTest(unittest.TestCase):

    def test_island_with_interior_gap_carved(self):
        runs = [(0, 40000, True), (40000, 60000, False), (60000, 100000, True)]
        pieces, n, carved = carve(GRS, runs, gaps=[(49000, 50000)])
        self.assertEqual(
            pieces,
            [
                ("chr1", 0, 40000, "putative_GRS", 12.0, SUBTYPE),
                ("chr1", 40000, 60000, "putative_core", BASELINE, "core_island_gap"),
                ("chr1", 60000, 100000, "putative_GRS", 12.0, SUBTYPE),
            ],
        )
        self.assertEqual((n, carved), (1, 20000))

    def test_high_run_without_gap_not_carved(self):
        runs = [(0, 40000, True), (40000, 60000, False), (60000, 100000, True)]
        row, n, _ = carve(GRS, runs)
        self.assertEqual(row, [GRS])
        self.assertEqual(n, 0)

    def test_gap_at_edge_counts(self):
        # The gap follows the HIGH run and is carried forward into it, so the island
        # span ends at the gap's far edge.
        runs = [(0, 40000, True), (40000, 60000, False), (60000, 100000, True)]
        pieces, n, _ = carve(GRS, runs, gaps=[(60000, 61000)])
        self.assertEqual(
            pieces[1],
            ("chr1", 40000, 61000, "putative_core", BASELINE, "core_island_gap"),
        )
        self.assertEqual(n, 1)

    def test_short_grs_untouched(self):
        short = ("chr1", 0, 20000, "putative_GRS", 12.0, SUBTYPE)
        runs = [(0, 10000, True), (10000, 20000, False)]
        row, n, _ = carve(short, runs, gaps=[(15000, 16000)])
        self.assertEqual(row, [short])
        self.assertEqual(n, 0)

    def test_gap_inside_low_run_not_carved(self):
        runs = [(0, 100000, True)]
        row, n, _ = carve(GRS, runs, gaps=[(49000, 50000)])
        self.assertEqual(row, [GRS])
        self.assertEqual(n, 0)

    def test_short_high_blip_absorbed(self):
        # 2 kb HIGH blip costs 2000 as LOW; a round trip costs 10000 in switches.
        runs = [(0, 49000, True), (49000, 51000, False), (51000, 100000, True)]
        row, n, _ = carve(GRS, runs, gaps=[(49500, 50000)])
        self.assertEqual(row, [GRS])
        self.assertEqual(n, 0)

    def test_masked_repeat_run_imputed(self):
        # LOW 45000-65000 is fully RED-covered: masked (free), so the HIGH stretch
        # 30000-70000 costs only two switches (10000) against 20000 for LOW-everywhere.
        # If the masked run were counted as LOW evidence, no island would form.
        runs = [
            (0, 30000, True),
            (30000, 45000, False),
            (45000, 65000, True),
            (65000, 70000, False),
            (70000, 100000, True),
        ]
        pieces, n, carved = carve(
            GRS, runs, gaps=[(67000, 68000)], red=[(45000, 65000)]
        )
        self.assertEqual(
            pieces[1],
            ("chr1", 30000, 70000, "putative_core", BASELINE, "core_island_gap"),
        )
        self.assertEqual((n, carved), (1, 40000))

    def test_partial_repeat_run_not_masked(self):
        # Same layout, RED covers 50% of the LOW run: it stays informative. LOW-everywhere
        # costs 19000; any HIGH path must pay the 20000 LOW-run mismatch or switches, so
        # the stretch is not carved.
        runs = [
            (0, 30000, True),
            (30000, 45000, False),
            (45000, 65000, True),
            (65000, 69000, False),
            (69000, 100000, True),
        ]
        row, n, _ = carve(
            GRS, runs, gaps=[(66000, 67000)], red=[(45000, 55000)]
        )
        self.assertEqual(row, [GRS])
        self.assertEqual(n, 0)

    def test_unplaced_subtype_untouched(self):
        unplaced = ("chr1", 0, 100000, "putative_GRS", 12.0, "unplaced_coverage_drop")
        runs = [(0, 40000, True), (40000, 60000, False), (60000, 100000, True)]
        row, n, _ = carve(unplaced, runs, gaps=[(49000, 50000)])
        self.assertEqual(row, [unplaced])
        self.assertEqual(n, 0)


if __name__ == "__main__":
    unittest.main()
