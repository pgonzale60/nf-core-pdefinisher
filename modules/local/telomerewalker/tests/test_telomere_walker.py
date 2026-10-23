#!/usr/bin/env python3
"""Unit tests for TELOMEREWALKER component logic."""

import json
import os
import tempfile
import unittest
from pathlib import Path

# Add resources/usr/bin and lib to sys.path
BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
LIB_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "lib"
import sys
sys.path.insert(0, str(BIN_DIR))
sys.path.insert(0, str(LIB_DIR))

import pde_telomere_walker


class TestTelomereWalker(unittest.TestCase):

    def test_parse_blueprint_preserves_terminus_card_over_scaffold_cards(self):
        """Verifies that terminus cards are not clobbered by subsequent unplaced scaffold cards."""
        blueprint_data = {
            "action_cards": [
                {
                    "card_id": "CARD-TELO-03",
                    "domain": "DOMAIN_TERMINUS",
                    "chrom": "II",
                    "terminus": "5p",
                    "action": "SCAFFOLD_REPEAT_TO_TELOMERE",
                    "classification": "SATELLITE_LINKED_TELOMERE",
                    "evidence": {
                        "associated_scaffolds": ["II_unloc_3", "II_unloc_1"],
                        "telomeric_tract_bp": 10408,
                    }
                },
                {
                    "card_id": "CARD-SCAFF-INT-01",
                    "domain": "DOMAIN_SCAFFOLD",
                    "chrom": "II",
                    "terminus": "5p",
                    "scaffold": "II_unloc_1",
                    "action": "INTEGRATE_INTO_CHROMOSOME",
                    "classification": "TERMINAL_EXTENSION",
                },
                {
                    "card_id": "CARD-SCAFF-INT-02",
                    "domain": "DOMAIN_SCAFFOLD",
                    "chrom": "II",
                    "terminus": "5p",
                    "scaffold": "II_unloc_3",
                    "action": "INTEGRATE_INTO_CHROMOSOME",
                    "classification": "TERMINAL_EXTENSION",
                }
            ]
        }

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
            json.dump(blueprint_data, tf)
            tf_path = tf.name

        try:
            term_cards, _, _ = pde_telomere_walker.parse_blueprint(tf_path)
            self.assertIn(("II", "5p"), term_cards)
            card = term_cards[("II", "5p")]
            self.assertEqual(card["card_id"], "CARD-TELO-03")
            self.assertEqual(card["classification"], "SATELLITE_LINKED_TELOMERE")
            self.assertEqual(card["action"], "SCAFFOLD_REPEAT_TO_TELOMERE")
        finally:
            if os.path.exists(tf_path):
                os.unlink(tf_path)

    def test_run_iterative_walker_zero_unmapped_status(self):
        """Verifies that empty unmapped read pool produces ZERO_EXTENDING_READS status."""
        from telowalking import run_iterative_walker
        class DummyEmptyIndex:
            def __len__(self):
                return 0

        res_germ = run_iterative_walker("dummy.bam", "ACGT" * 100, "TTAGGC", "GCCTAA",
                                        max_iters=100, unmapped=DummyEmptyIndex(), tip_is_germline=True)
        self.assertEqual(res_germ["status"], "ZERO_EXTENDING_READS")
        self.assertIn("Zero unmapped reads", res_germ["notes"])
        self.assertIn("distal GRS", res_germ["notes"])

        res_soma = run_iterative_walker("dummy.bam", "ACGT" * 100, "TTAGGC", "GCCTAA",
                                        max_iters=100, unmapped=DummyEmptyIndex(), is_somatic_boundary=True)
        self.assertEqual(res_soma["status"], "ZERO_EXTENDING_READS")
        self.assertIn("somatic diminution breaksite", res_soma["notes"])

    def test_abpoa_overhang_consensus(self):
        """Verifies abpoa computes SIMD-POA consensus from multiple soft-clipped overhangs."""
        import shutil
        from telocapping import overhang_capper
        if not shutil.which("abpoa"):
            self.skipTest("abpoa not available in test environment")

        mock_lines = [
            "read1\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15 + "AAAA") + "\t*",
            "read2\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15 + "AAAT") + "\t*",
        ]
        # Monkey patch terminal_window_lines
        orig_twl = overhang_capper.terminal_window_lines
        overhang_capper.terminal_window_lines = lambda *args, **kwargs: mock_lines
        try:
            res = overhang_capper.maximal_tip_overhang("dummy.bam", "chrI", "3p", 10000, window=1000, chrom_start=0)
            self.assertIsNotNone(res)
            label, clip_len, cons_seq = res
            self.assertIn("abpoa_consensus", label)
            self.assertIn("TTAGGC" * 10, cons_seq)
            self.assertGreaterEqual(clip_len, 90)
        finally:
            overhang_capper.terminal_window_lines = orig_twl

    def test_extract_cluster_reads_from_card(self):
        """Verifies cluster reads are only extracted from genuine junction cluster cards."""
        cluster_card = {
            "mechanism": "SCAFFOLD_REPEAT_TO_TELOMERE",
            "evidence": {
                "junction_cluster_id": "CL_01",
                "supporting_reads": ["readA", "readB"],
            }
        }
        self.assertEqual(
            pde_telomere_walker.extract_cluster_reads_from_card(cluster_card),
            ["readA", "readB"]
        )

        read_cap_card = {
            "mechanism": "CAP_WITH_READ",
            "authorized_read_qname": "readC",
            "evidence": {
                "read_count": 1,
                "supporting_reads": ["readC"],
            }
        }
        self.assertIsNone(pde_telomere_walker.extract_cluster_reads_from_card(read_cap_card))

        legacy_read_card = {
            "mechanism": "CAP_WITH_READ",
            "evidence": {
                "supporting_reads": 1,
            }
        }
        self.assertIsNone(pde_telomere_walker.extract_cluster_reads_from_card(legacy_read_card))

    def test_iterative_walker_handles_scalar_or_invalid_cluster_reads(self):
        """Verifies that non-iterable cluster_reads (e.g. legacy int count) does not raise TypeError."""
        from telowalking import run_iterative_walker
        class DummyEmptyIndex:
            def __len__(self):
                return 0

        # Should handle scalar int without raising TypeError: 'int' object is not iterable
        res = run_iterative_walker(
            "dummy.bam", "ACGT" * 100, "TTAGGC", "GCCTAA",
            max_iters=10, unmapped=DummyEmptyIndex(),
            tip_is_germline=True, cluster_reads=1
        )
        self.assertEqual(res["status"], "ZERO_EXTENDING_READS")


    def test_align_seed_to_unitigs(self):
        """Verifies local unitig bridging aligns outward tip seed to telomeric unitigs."""
        import shutil
        from telocapping import align_seed_to_unitigs
        if not shutil.which("minimap2"):
            self.skipTest("minimap2 not available in test environment")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Chromosome seed: non-palindromic sequence ending at tip
            seed = ("GACTACGATAGCTAGGACATTAGCTA" * 20)
            # Unitig: seed followed by 300 bp canonical telomere TTAGGC
            unitig_seq = seed + ("TTAGGC" * 50)
            u_path = os.path.join(tmpdir, "unitigs.fa")
            with open(u_path, "w") as uf:
                uf.write(f">utg001\n{unitig_seq}\n")

            res = align_seed_to_unitigs(seed, u_path, "TTAGGC", "GCCTAA", min_overlap=200, min_telo_repeats=15)
            self.assertIsNotNone(res)
            self.assertEqual(res["unitig_id"], "utg001")
            self.assertGreaterEqual(res["added_bp"], 300)
            self.assertGreaterEqual(res["telo_bp"], 300)
            self.assertTrue(res["extension_seq"].startswith("TTAGGC"))

if __name__ == "__main__":
    unittest.main()
