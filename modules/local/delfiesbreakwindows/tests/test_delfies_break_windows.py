#!/usr/bin/env python3
"""Unit tests for pde_delfies_break_windows.py."""

import importlib.util
import pathlib
import tempfile
import unittest

SCRIPT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "resources" / "usr" / "bin" / "pde_delfies_break_windows.py"
)
spec = importlib.util.spec_from_file_location("pde_delfies_break_windows", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestBalance(unittest.TestCase):
    def test_both_strands_within_half_of_max_is_balanced(self):
        self.assertTrue(mod.is_balanced(10, 6))

    def test_disparity_over_half_is_not_balanced(self):
        self.assertFalse(mod.is_balanced(10, 4))

    def test_single_strand_is_not_balanced(self):
        self.assertFalse(mod.is_balanced(12, 0))


class TestBorders(unittest.TestCase):
    def test_core_grs_contiguous_pair_gives_border(self):
        segs = {"c": [(0, 1000, "putative_core"), (1000, 2000, "putative_GRS")]}
        self.assertEqual(mod.find_borders(segs), {"c": [1000]})

    def test_non_contiguous_pair_is_not_border(self):
        segs = {"c": [(0, 1000, "putative_core"), (1100, 2000, "putative_GRS")]}
        self.assertEqual(mod.find_borders(segs), {"c": []})

    def test_near_border_respects_tolerance(self):
        borders = {"c": [1000]}
        self.assertTrue(mod.near_border(borders, "c", 1100))
        self.assertFalse(mod.near_border(borders, "c", 1101))


class TestCallSites(unittest.TestCase):
    def test_floor_is_half_median_of_top_three(self):
        sites = {("c", 1): (10, 0), ("c", 2): (8, 0), ("c", 3): (6, 0), ("c", 4): (1, 0)}
        self.assertEqual(mod.support_floor(sites), 0.5 * 8)

    def test_balanced_adjacent_pair_is_fission_anywhere(self):
        sites = {("c", 100): (10, 8), ("c", 101): (9, 7)}
        calls, _ = mod.call_sites(sites, {})
        self.assertEqual(calls["c"][100], "fission")

    def test_unbalanced_adjacent_pair_is_not_fission(self):
        sites = {("c", 100): (10, 8), ("c", 101): (20, 0)}
        calls, _ = mod.call_sites(sites, {"c": [100]})
        self.assertEqual(calls["c"][100], "break")

    def test_site_off_border_is_dropped(self):
        sites = {("c", 500): (10, 0), ("c", 900): (8, 2)}
        calls, _ = mod.call_sites(sites, {"c": [900]})
        self.assertNotIn(500, calls["c"])
        self.assertEqual(calls["c"][900], "break")


class TestWindows(unittest.TestCase):
    def test_close_sites_merge_and_window_is_centred(self):
        groups = mod.cluster([100, 105, 110])
        self.assertEqual(groups, [[100, 105, 110]])
        start, end = mod.window_for(groups[0], 10000)
        self.assertEqual(end - start, 150)
        self.assertEqual(start, 105 - 75)

    def test_window_clipped_at_contig_start(self):
        self.assertEqual(mod.window_for([10], 1000), (0, 150))


class TestBuild(unittest.TestCase):
    def test_end_to_end_break_at_border_and_fission(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            tsv = tmp / "s2g.tsv"
            tsv.write_text(
                "contig\tstart\tend\tread_depth\tbreakpoint_type\tsample_name\t"
                "num_supporting_reads__forward\tnum_supporting_reads__reverse\n"
                "c1\t500\t501\t30\tS2G\ts\t10\t0\n"
                "c1\t900\t901\t30\tS2G\ts\t8\t2\n"
                "c1\t1500\t1501\t30\tS2G\ts\t10\t8\n"
                "c1\t1501\t1502\t30\tS2G\ts\t9\t7\n"
            )
            segs = tmp / "pde_segments.bed"
            segs.write_text(
                "c1\t0\t1000\tputative_core\n"
                "c1\t1000\t2000\tputative_GRS\n"
            )
            fa = tmp / "asm.fa"
            fa.write_text(">c1\n" + "A" * 2000 + "\n")
            prefix = str(tmp / "out")
            mod.build(str(tsv), str(segs), str(fa), prefix)
            bed = pathlib.Path(prefix + ".breakwindows.bed").read_text().splitlines()
            kinds = sorted(line.split("\t")[6] for line in bed)
            self.assertEqual(kinds, ["break", "fission"])
            fasta = pathlib.Path(prefix + ".breakwindows.fa").read_text()
            self.assertEqual(fasta.count(">"), 2)


if __name__ == "__main__":
    unittest.main()
