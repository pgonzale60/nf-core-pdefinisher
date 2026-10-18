#!/usr/bin/env python3
"""Unit tests for CURATIONAGPGENERATOR terminal N-gap stripping and AGP block integrity."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# Add resources/usr/bin to sys.path
BIN_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(BIN_DIR))

from pde_generate_agp import AGPBlock, strip_terminal_gaps


class TestCurationAGPGenerator(unittest.TestCase):

    def test_strip_terminal_gaps_leading_and_trailing(self):
        comp_pool = {
            "c_lead": "NNNNATCGATCG",
            "c_mid": "GGCCAATT",
            "c_trail": "ATCGATCGNNNN",
        }
        blocks = [
            AGPBlock(type="N", component_id="gap_start", length=100, orientation=None, source_seq="N" * 100),
            AGPBlock(type="W", component_id="c_lead", length=12, orientation="+", source_seq=comp_pool["c_lead"]),
            AGPBlock(type="N", component_id="gap_mid", length=50, orientation=None, source_seq="N" * 50),
            AGPBlock(type="W", component_id="c_mid", length=8, orientation="+", source_seq=comp_pool["c_mid"]),
            AGPBlock(type="W", component_id="c_trail", length=12, orientation="+", source_seq=comp_pool["c_trail"]),
            AGPBlock(type="N", component_id="gap_end", length=200, orientation=None, source_seq="N" * 200),
        ]

        stripped = strip_terminal_gaps(blocks, comp_pool)
        self.assertEqual(len(stripped), 4)
        self.assertEqual(stripped[0].type, "W")
        self.assertEqual(stripped[0].component_id, "c_lead")
        self.assertEqual(stripped[0].length, 8)
        self.assertEqual(comp_pool["c_lead"], "ATCGATCG")

        # Internal gap must be preserved
        self.assertEqual(stripped[1].type, "N")
        self.assertEqual(stripped[1].component_id, "gap_mid")
        self.assertEqual(stripped[2].type, "W")
        self.assertEqual(stripped[2].component_id, "c_mid")

        # Trailing gap stripped and trailing W trimmed
        self.assertEqual(stripped[3].type, "W")
        self.assertEqual(stripped[3].component_id, "c_trail")
        self.assertEqual(stripped[3].length, 8)
        self.assertEqual(comp_pool["c_trail"], "ATCGATCG")

    def test_strip_terminal_gaps_pure_n_blocks(self):
        comp_pool = {
            "c_all_n": "NNNNNNNN",
            "c_mid": "ACGTACGT",
        }
        blocks = [
            AGPBlock(type="W", component_id="c_all_n", length=8, orientation="+", source_seq=comp_pool["c_all_n"]),
            AGPBlock(type="N", component_id="gap_start", length=50, orientation=None, source_seq="N" * 50),
            AGPBlock(type="W", component_id="c_mid", length=8, orientation="+", source_seq=comp_pool["c_mid"]),
        ]
        stripped = strip_terminal_gaps(blocks, comp_pool)
        self.assertEqual(len(stripped), 1)
        self.assertEqual(stripped[0].component_id, "c_mid")
        self.assertNotIn("c_all_n", comp_pool)

    def test_e2e_agp_generation_with_terminal_gap_stripping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_path = os.path.join(tmpdir, "ref.fa")
            with open(ref_path, "w") as f:
                # Chromosome with an internal gap: [1000 bp][100 bp Ns][1000 bp]
                f.write(">chr1\n" + ("A" * 1000) + ("N" * 100) + ("C" * 1000) + "\n")

            blueprint_path = os.path.join(tmpdir, "blueprint.json")
            # Cleave chr1 at pos 1050 (inside the scaffolding N-gap)
            blueprint = {
                "action_cards": [
                    {
                        "card_id": "CARD-SPLIT-01",
                        "chrom": "chr1",
                        "domain": "DOMAIN_INTERNAL",
                        "action": "EXCISE_GAP",
                        "start_pos": 1000,
                        "end_pos": 1100,
                        "classification": "EXCISE_SCAFFOLDING_GAP"
                    }
                ]
            }
            with open(blueprint_path, "w") as f:
                json.dump(blueprint, f)

            outdir = os.path.join(tmpdir, "out")
            cmd = [
                sys.executable,
                str(BIN_DIR / "pde_generate_agp.py"),
                "--blueprint", blueprint_path,
                "--ref", ref_path,
                "--sample", "TEST_SAMPLE",
                "--outdir", outdir,
                "--synthesize-fasta"
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}")

            agp_file = os.path.join(outdir, "TEST_SAMPLE.curated.agp")
            self.assertTrue(os.path.exists(agp_file))
            with open(agp_file) as f:
                lines = [line.strip().split("\t") for line in f if line.strip() and not line.startswith("#")]

            # Verify that neither first nor last component is an N-gap
            self.assertTrue(len(lines) > 0)
            self.assertNotEqual(lines[0][4], "N", "First AGP record must not be an N-gap")
            self.assertNotEqual(lines[-1][4], "N", "Last AGP record must not be an N-gap")


if __name__ == "__main__":
    unittest.main()
