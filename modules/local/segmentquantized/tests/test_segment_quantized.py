"""Unit tests for pde_segment_quantized: final three-class call from quantized depth.

British spelling standard enforced across all documentation and comments.
"""

from pathlib import Path
import os
import sys
import tempfile
import unittest

bin_dir = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(bin_dir))

from pde_segment_quantized import (
    call_chromosome,
    check_tiling,
    elimination_line,
    main,
)

POLICY = {
    "anchor_proximity_bp": 2000,
    "quantized_min_low_run_bp": 1,
    "quantized_max_bridge_bp": 5000,
    "quantized_masked_repeat_fraction": 0.9,
    "chromosomal_min_bp": 1000000,
}
LENGTH = 300000
CHROM = "chr1"


def run(s, e, is_low):
    return (s, e, is_low)


def bp(pos, strand, end=None):
    return (pos, strand, pos if end is None else end)


def call(runs, gaps=None, red=None, bps=None, length=LENGTH):
    pieces = call_chromosome(CHROM, length, runs, gaps or [], red, bps or [], POLICY)
    check_tiling(CHROM, pieces, length)
    return pieces


BASE_RUNS = [run(0, 100000, False), run(100000, 150000, True), run(150000, LENGTH, False)]


class TestQuantizedCall(unittest.TestCase):
    def test_low_run_with_breakpoint_at_left_border_is_grs(self):
        pieces = call(BASE_RUNS, bps=[bp(100500, "+")])
        self.assertEqual(pieces, [
            (0, 100000, "putative_core", "-"),
            (100000, 150000, "putative_GRS", "quantized_delfies_anchored"),
            (150000, LENGTH, "putative_core", "-"),
        ])

    def test_low_run_without_breakpoint_is_uncertain(self):
        pieces = call(BASE_RUNS)
        self.assertIn((100000, 150000, "putative_uncertain", "no_delfies_support"), pieces)

    def test_low_run_with_minus_breakpoint_at_right_border_is_grs(self):
        runs = [run(0, 50000, True), run(50000, LENGTH, False)]
        pieces = call(runs, bps=[bp(0, "-", end=50300)])
        self.assertEqual(pieces[0], (0, 50000, "putative_GRS", "quantized_delfies_anchored"))

    def test_breakpoint_near_right_border_only_does_not_anchor_plus(self):
        pieces = call(BASE_RUNS, bps=[bp(150100, "+")])
        self.assertIn((100000, 150000, "putative_uncertain", "no_delfies_support"), pieces)

    def test_gap_inside_low_run_is_uncertain_and_split_from_grs(self):
        pieces = call(BASE_RUNS, gaps=[(120000, 121000)], bps=[bp(100000, "+")])
        self.assertEqual(pieces[1:4], [
            (100000, 120000, "putative_GRS", "quantized_delfies_anchored"),
            (120000, 121000, "putative_uncertain", "gap_spanning"),
            (121000, 150000, "putative_GRS", "quantized_delfies_anchored"),
        ])

    def test_short_high_bridge_carries_anchor_to_next_low_run(self):
        runs = [run(0, 100000, False), run(100000, 150000, True),
                run(150000, 152000, False), run(152000, 160000, True),
                run(160000, LENGTH, False)]
        pieces = call(runs, bps=[bp(100500, "+")])
        self.assertEqual(pieces[1:4], [
            (100000, 150000, "putative_GRS", "quantized_delfies_anchored"),
            (150000, 152000, "putative_GRS", "quantized_bridged"),
            (152000, 160000, "putative_GRS", "quantized_delfies_anchored"),
        ])

    def test_long_high_bridge_breaks_block(self):
        runs = [run(0, 100000, False), run(100000, 150000, True),
                run(150000, 160000, False), run(160000, 170000, True),
                run(170000, LENGTH, False)]
        pieces = call(runs, bps=[bp(100500, "+")])
        self.assertEqual(pieces[1:4], [
            (100000, 150000, "putative_GRS", "quantized_delfies_anchored"),
            (150000, 160000, "putative_core", "-"),
            (160000, 170000, "putative_uncertain", "no_delfies_support"),
        ])

    def test_unanchored_block_keeps_bridge_as_core(self):
        runs = [run(0, 100000, False), run(100000, 150000, True),
                run(150000, 152000, False), run(152000, 160000, True),
                run(160000, LENGTH, False)]
        pieces = call(runs)
        self.assertEqual(pieces[1:4], [
            (100000, 150000, "putative_uncertain", "no_delfies_support"),
            (150000, 152000, "putative_core", "-"),
            (152000, 160000, "putative_uncertain", "no_delfies_support"),
        ])

    def test_minus_breakpoint_at_block_end_anchors_whole_block(self):
        runs = [run(0, 100000, False), run(100000, 150000, True),
                run(150000, 152000, False), run(152000, 160000, True),
                run(160000, LENGTH, False)]
        pieces = call(runs, bps=[bp(0, "-", end=159500)])
        self.assertEqual(pieces[1:4], [
            (100000, 150000, "putative_GRS", "quantized_delfies_anchored"),
            (150000, 152000, "putative_GRS", "quantized_bridged"),
            (152000, 160000, "putative_GRS", "quantized_delfies_anchored"),
        ])

    def test_fully_repeat_low_run_is_uncertain_without_breakpoint_and_grs_with(self):
        red = [(100000, 150000)]
        self.assertIn((100000, 150000, "putative_uncertain", "repeat_masked"),
                      call(BASE_RUNS, red=red))
        self.assertIn((100000, 150000, "putative_GRS", "quantized_delfies_anchored"),
                      call(BASE_RUNS, red=red, bps=[bp(100000, "+")]))

    def test_low_run_below_minimum_length_is_core(self):
        policy = dict(POLICY, quantized_min_low_run_bp=3000)
        runs = [run(0, 100000, False), run(100000, 102000, True), run(102000, LENGTH, False)]
        bps = [bp(100000, "+")]
        pieces = call_chromosome(CHROM, LENGTH, runs, [], None, bps, policy)
        self.assertEqual(pieces, [(0, LENGTH, "putative_core", "-")])

    def test_uncovered_bases_are_uncertain_and_tiling_holds(self):
        runs = [run(0, 40000, False), run(50000, LENGTH, False)]
        pieces = call(runs)
        self.assertIn((40000, 50000, "putative_uncertain", "no_quantized_coverage"), pieces)

    def test_check_tiling_rejects_one_base_hole(self):
        with self.assertRaises(ValueError):
            check_tiling(CHROM, [(0, 10, "putative_core", "-"), (11, 20, "putative_core", "-")], 20)

    def test_adjacent_identical_pieces_are_merged(self):
        runs = [run(0, 100000, False), run(100000, 110000, False), run(110000, LENGTH, False)]
        pieces = call(runs)
        self.assertEqual(pieces, [(0, LENGTH, "putative_core", "-")])


class TestEliminationLine(unittest.TestCase):
    def test_grs_bp_and_uncertain_reported_separately(self):
        segments = [
            ("chr1", 0, 1000, "putative_GRS", 5.0, "quantized_delfies_anchored"),
            ("chr1", 1000, 3000, "putative_GRS", 5.0, "quantized_delfies_anchored"),
            ("chr1", 3000, 4000, "putative_uncertain", 5.0, "no_delfies_support"),
            ("u1", 0, 500, "putative_uncertain", 1.0, "fully_repetitive_unplaced"),
        ]
        line = elimination_line(segments, canonical={"chr1"})
        self.assertIn("canon_grs_bp=3000", line)
        self.assertIn("unplaced_uncertain_bp=500", line)
        self.assertIn("min_elim_bp=3000\tmax_elim_bp=3000\test_elim_bp=3000", line)
        self.assertIn("canon_uncertain_bp=1000", line)


class TestMainFixture(unittest.TestCase):
    def test_cli_writes_bed_and_summary(self):
        with tempfile.TemporaryDirectory() as d:
            def w(name, text):
                path = os.path.join(d, name)
                with open(path, "w") as f:
                    f.write(text)
                return path

            summary = w("in.tsv",
                        "# SOFTWARE=pde_identify_grs.py\tGRS_THRESHOLD=10.00\n"
                        "# ELIMINATION_RANGE\tcanon_grs_bp=0\n"
                        "chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype\n"
                        "chr1\t0\t300000\t300000\tputative_core\t30.00\t-\n"
                        "u1\t0\t5000\t5000\tputative_uncertain\t1.00\tfully_repetitive_unplaced\n")
            quant = w("q.bed", "chr1\t0\t100000\t10:inf\t0\nchr1\t100000\t150000\t0:10\t0\n"
                               "chr1\t150000\t300000\t10:inf\t0\n")
            bps = w("bp.tsv", "chrom\tpos\tend\tstrand\tsupport\nchr1\t100500\t100500\t+\t8.0\n")
            outdir = os.path.join(d, "out")
            os.makedirs(outdir)
            main(["--summary-tsv", summary, "--quantized-bed", quant, "--breakpoints-tsv", bps,
                  "--prefix", "t", "--outdir", outdir, "--chromosomes", "chr1"])
            with open(os.path.join(outdir, "t.pde_segments.bed")) as f:
                bed = f.read().splitlines()
            with open(os.path.join(outdir, "t.grs_summary.tsv")) as f:
                tsv = f.read().splitlines()
        self.assertEqual(bed, [
            "chr1\t0\t100000\tputative_core",
            "chr1\t100000\t150000\tputative_GRS",
            "chr1\t150000\t300000\tputative_core",
            "u1\t0\t5000\tputative_uncertain",
        ])
        self.assertTrue(tsv[1].startswith("# ELIMINATION_RANGE\tcanon_grs_bp=50000\t"))
        self.assertIn("unplaced_uncertain_bp=5000", tsv[1])
        self.assertIn("min_elim_bp=50000\tmax_elim_bp=50000\test_elim_bp=50000", tsv[1])


if __name__ == "__main__":
    unittest.main()
