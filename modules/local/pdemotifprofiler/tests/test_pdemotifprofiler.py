#!/usr/bin/env python3
"""Unit tests for PDEMOTIFPROFILER component logic."""

import os
import tempfile
import unittest
from unittest import mock
from pathlib import Path
import sys

BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(BIN_DIR))

import pde_motif_profiler


class TestPdeMotifProfiler(unittest.TestCase):

    def test_parse_breakpoints(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bp_bed = os.path.join(tmpdir, "breakpoints.bed")
            with open(bp_bed, "w") as f:
                f.write("chr1\t1000\t1001\tbp1\t100\t+\n")
                f.write("chr1\t5000\t5001\tbp2\t100\t-\n")

            bps = pde_motif_profiler.parse_breakpoints(bp_bed, flank_tolerance=500)
            self.assertIn("chr1", bps)
            self.assertEqual(len(bps["chr1"]), 2)
            self.assertEqual(bps["chr1"][0]["zone_start"], 500)
            self.assertEqual(bps["chr1"][0]["zone_end"], 1501)

    def test_parse_segments(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            seg_bed = os.path.join(tmpdir, "segments.bed")
            with open(seg_bed, "w") as f:
                f.write("chr1\t0\t1000\tputative_GRS\n")
                f.write("chr1\t1000\t5000\tputative_core\n")

            segs = pde_motif_profiler.parse_segments(seg_bed)
            self.assertIn("chr1", segs)
            self.assertEqual(len(segs["chr1"]), 2)
            self.assertEqual(segs["chr1"][1]["label"], "putative_core")

    def test_hit_in_uncertain_segment_is_not_core(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            seg_bed = os.path.join(tmpdir, "segments.bed")
            with open(seg_bed, "w") as f:
                f.write("chr1\t0\t1000\tputative_core\n")
                f.write("chr1\t1000\t2000\tuncertain\n")
            bp_bed = os.path.join(tmpdir, "breakpoints.bed")
            open(bp_bed, "w").close()
            fimo = os.path.join(tmpdir, "fimo.tsv")
            with open(fimo, "w") as f:
                f.write("motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n")
                f.write("M1\t\tchr1\t1200\t1210\t+\t10.0\t1e-5\t\tACGT\n")
            outdir = os.path.join(tmpdir, "out")
            argv = ["prog", "--fimo-tsv", fimo, "--breakpoints-bed", bp_bed,
                    "--pde-segments-bed", seg_bed, "--prefix", "t", "--outdir", outdir]
            with mock.patch.object(sys, "argv", argv):
                pde_motif_profiler.main()
            with open(os.path.join(outdir, "t.motif_distribution_summary.tsv")) as f:
                metrics = dict(line.rstrip("\n").split("\t") for line in f.readlines()[1:])
            self.assertEqual(metrics["significant_hits_in_core"], "0")
            self.assertEqual(metrics["significant_hits_in_uncertain"], "1")


if __name__ == "__main__":
    unittest.main()
