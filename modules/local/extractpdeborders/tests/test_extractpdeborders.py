#!/usr/bin/env python3
"""Unit tests for EXTRACTPDEBORDERS component logic."""

import os
import tempfile
import unittest
from pathlib import Path
import sys

BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(BIN_DIR))

import pde_extract_borders


class TestExtractPdeBorders(unittest.TestCase):

    def test_revcomp(self):
        self.assertEqual(pde_extract_borders.revcomp("ACGTN"), "NACGT")

    def test_parse_segments_and_derive_frontiers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            seg_bed = os.path.join(tmpdir, "pde_segments.bed")
            with open(seg_bed, "w") as f:
                f.write("chr1\t0\t1000\tputative_GRS\n")
                f.write("chr1\t1000\t5000\tputative_core\n")
                f.write("chr1\t5000\t6000\tputative_GRS\n")

            segs = pde_extract_borders.parse_segments(seg_bed)
            self.assertIn("chr1", segs)
            self.assertEqual(len(segs["chr1"]), 3)

            # Test frontier transition logic
            seg_list = segs["chr1"]
            frontiers = []
            for i in range(len(seg_list) - 1):
                s1, e1, lab1 = seg_list[i]
                s2, e2, lab2 = seg_list[i + 1]
                if e1 == s2:
                    if "GRS" in lab1 and "core" in lab2:
                        frontiers.append((e1, "-"))
                    elif "core" in lab1 and "GRS" in lab2:
                        frontiers.append((e1, "+"))

            self.assertEqual(len(frontiers), 2)
            self.assertEqual(frontiers[0], (1000, "-"))
            self.assertEqual(frontiers[1], (5000, "+"))

    def test_compute_census_metrics_full(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            grs_tsv = os.path.join(tmpdir, "grs_summary.tsv")
            with open(grs_tsv, "w") as f:
                f.write("# SOMATIC_CORE_BASELINE=50.0\n")
                f.write("chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tratio\n")
                f.write("chr1\t0\t1000\t1000\tputative_GRS\t5.0\t0.1\n")
                f.write("chr1\t1000\t5000\t4000\tputative_core\t50.0\t1.0\n")
                f.write("chr1\t5000\t6000\t1000\tputative_GRS\t7.0\t0.14\n")

            ends_tsv = os.path.join(tmpdir, "terminal_ends_summary.tsv")
            with open(ends_tsv, "w") as f:
                f.write("chrom_end\tchrom_len\tstatus\tn_telo_reads\tn_subtelo_reads\tproposed_action\tterminal_telomere_tract_length\n")
                f.write("chr1_5p\t6000\tEXTEND_AND_CAP\t15\t5\tEXTEND_AND_CAP\t250\n")
                f.write("chr1_3p\t6000\tALREADY_CAPPED\t20\t2\tALREADY_CAPPED\t400\n")

            segs_by_chrom = {
                "chr1": [
                    (0, 1000, "putative_GRS"),
                    (1000, 5000, "putative_core"),
                    (5000, 6000, "putative_GRS"),
                ]
            }
            frontier_bps = [
                {"chrom": "chr1", "start": 1000, "end": 1001, "name": "bp1", "score": "100", "strand": "-"},
                {"chrom": "chr1", "start": 5000, "end": 5001, "name": "bp2", "score": "100", "strand": "+"},
            ]

            metrics = pde_extract_borders.compute_census_metrics(
                frontier_breakpoints=frontier_bps,
                sample_id="test_sample",
                segs_by_chrom=segs_by_chrom,
                grs_summary_tsv=grs_tsv,
                terminal_ends_summary_tsv=ends_tsv,
            )

            self.assertEqual(metrics["sample_id"], "test_sample")
            self.assertEqual(metrics["n_breaks"], 2)
            self.assertEqual(metrics["median_c_soma"], 50.0)
            self.assertEqual(metrics["median_c_germ"], 6.0)
            self.assertEqual(metrics["germline_fraction"], 0.12)
            self.assertEqual(metrics["n_germ_telo"], 2)
            self.assertEqual(metrics["estimated_germline_chromosomes"], 1)

    def test_compute_census_metrics_empty(self):
        metrics = pde_extract_borders.compute_census_metrics(
            frontier_breakpoints=[],
            sample_id="test_sample",
            segs_by_chrom={},
            grs_summary_tsv=None,
            terminal_ends_summary_tsv=None,
        )
        self.assertEqual(metrics["sample_id"], "test_sample")
        self.assertEqual(metrics["n_breaks"], 0)
        self.assertEqual(metrics["median_c_soma"], "N/A")
        self.assertEqual(metrics["median_c_germ"], "N/A")
        self.assertIsNone(metrics["germline_fraction"])
        self.assertEqual(metrics["n_germ_telo"], "N/A")
        self.assertEqual(metrics["estimated_germline_chromosomes"], "N/A")


if __name__ == "__main__":
    unittest.main()
