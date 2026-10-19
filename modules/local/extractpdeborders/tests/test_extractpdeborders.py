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


if __name__ == "__main__":
    unittest.main()
