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


if __name__ == "__main__":
    unittest.main()

