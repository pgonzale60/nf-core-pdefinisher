#!/usr/bin/env python3
"""Unit tests for TELOMERICJUNCTIONPROFILER component logic."""

import shutil
import sys
import unittest
from pathlib import Path

# Add resources/usr/bin to sys.path
BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(BIN_DIR))

import pde_telomeric_junction_profiler


class TestTelomericJunctionProfiler(unittest.TestCase):

    def setUp(self):
        self.has_vsearch = shutil.which("vsearch") is not None

    def test_vsearch_clustering_and_consensus(self):
        """Verifies vsearch clusters high-identity junction flanks and generates consensus sequences."""
        if not self.has_vsearch:
            self.skipTest("vsearch binary not available in test environment")

        flank_a1 = "GATC" * 15 + "AAAA" * 10
        flank_a2 = "GATC" * 15 + "AAAA" * 9 + "AAAG"
        flank_b1 = "TGCA" * 15 + "CCCC" * 10

        obs1 = ("read_1", "5p", 120, 120, flank_a1, "TTAGGC")
        obs2 = ("read_2", "5p", 150, 120, flank_a2, "TTAGGC")
        obs3 = ("read_3", "3p", 80, 50, flank_b1, "TTAGGC")

        clusters = pde_telomeric_junction_profiler.cluster_junction_flanks(
            [obs1, obs2, obs3], threads=1, id_threshold=0.97, core_flank_len=60
        )
        self.assertEqual(len(clusters), 2, f"Expected 2 clusters, got {len(clusters)}")
        sizes = sorted([len(c["reads"]) for c in clusters])
        self.assertEqual(sizes, [1, 2])

        # Verify consensus extensions
        for c in clusters:
            self.assertGreaterEqual(len(c["rep_ext"]), 60)
            self.assertEqual(len(c["rep_core"]), 60)
            self.assertEqual(c["rep_core"], c["rep_ext"][:60])

    def test_empty_observations_return_empty_clusters(self):
        """Verifies empty observation list yields empty cluster set."""
        clusters = pde_telomeric_junction_profiler.cluster_junction_flanks([], threads=1)
        self.assertEqual(clusters, [])


if __name__ == "__main__":
    unittest.main()
