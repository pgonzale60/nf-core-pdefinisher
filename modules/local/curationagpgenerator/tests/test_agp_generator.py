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

from pde_generate_agp import AGPBlock, assign_chromosome_pde_segments, merge_curated_pde_segments, strip_terminal_gaps


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

    def test_e2e_curated_bed_generation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_path = os.path.join(tmpdir, "ref.fa")
            with open(ref_path, "w") as f:
                # 2000 bp chromosome
                f.write(">chr1\n" + ("A" * 500) + ("N" * 100) + ("C" * 1400) + "\n")

            blueprint_path = os.path.join(tmpdir, "blueprint.json")
            blueprint = {
                "action_cards": [
                    {
                        "card_id": "CARD-EXCISE-01",
                        "chrom": "chr1",
                        "domain": "DOMAIN_INTERNAL",
                        "action": "EXCISE_GAP",
                        "start_pos": 500,
                        "end_pos": 600,
                        "classification": "EXCISE_SCAFFOLDING_GAP"
                    }
                ]
            }
            with open(blueprint_path, "w") as f:
                json.dump(blueprint, f)

            pde_bed = os.path.join(tmpdir, "pde_segments.bed")
            with open(pde_bed, "w") as f:
                f.write("chr1\t0\t500\tCORE\n")
                f.write("chr1\t500\t600\tGRS\n")
                f.write("chr1\t600\t2000\tCORE\n")

            red_bed = os.path.join(tmpdir, "red_repeats.bed")
            with open(red_bed, "w") as f:
                f.write("chr1\t100\t200\trepeat_A\n")
                f.write("chr1\t520\t580\trepeat_in_gap\n")
                f.write("chr1\t700\t800\trepeat_B\n")

            outdir = os.path.join(tmpdir, "out")
            cmd = [
                sys.executable,
                str(BIN_DIR / "pde_generate_agp.py"),
                "--blueprint", blueprint_path,
                "--ref", ref_path,
                "--sample", "TEST_SAMPLE",
                "--outdir", outdir,
                "--pde-segments", pde_bed,
                "--red-bed", red_bed
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}")

            curated_pde = os.path.join(outdir, "TEST_SAMPLE.curated.pde_segments.bed")
            curated_rep = os.path.join(outdir, "TEST_SAMPLE.curated.repeats.bed")

            # Check that both curated BED files exist
            self.assertTrue(os.path.exists(curated_pde), "Curated PDE segments BED must exist")
            self.assertTrue(os.path.exists(curated_rep), "Curated repeats BED must exist")

            # Verify NO TSV or input coordinate files were generated for segments/repeats
            out_files = os.listdir(outdir)
            self.assertNotIn("TEST_SAMPLE.curated.pde_segments.tsv", out_files)
            self.assertNotIn("TEST_SAMPLE.curated.repeats.tsv", out_files)
            self.assertNotIn("TEST_SAMPLE.input.pde_segments.bed", out_files)
            self.assertNotIn("TEST_SAMPLE.input.repeats.bed", out_files)

            # Validate PDE segment liftover coordinates (adjacent CORE regions merge after GRS excision)
            with open(curated_pde) as f:
                pde_lines = [l.strip().split("\t") for l in f if l.strip()]
            self.assertEqual(len(pde_lines), 1)
            self.assertEqual(pde_lines[0], ["chr1", "0", "1900", "CORE"])

            # Validate repeat liftover coordinates
            with open(curated_rep) as f:
                rep_lines = [l.strip().split("\t") for l in f if l.strip()]
            self.assertEqual(len(rep_lines), 2)
            self.assertEqual(rep_lines[0], ["chr1", "100", "200", "repeat_A"])
            self.assertEqual(rep_lines[1], ["chr1", "600", "700", "repeat_B"])

    def test_e2e_curated_bed_generation_with_whitespace_fasta_headers(self):
        """BED files containing space-separated descriptions (e.g. from RED or NCBI headers like DF5020) must match accessions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_path = os.path.join(tmpdir, "ref.fa")
            with open(ref_path, "w") as f:
                # FASTA header containing space-separated organism description (e.g. DF5020 / NCBI standard)
                f.write(">OZ059902.1 Oscheius myriophilus genome assembly, chromosome: I\n" + ("A" * 1000) + ("N" * 100) + ("C" * 1000) + "\n")

            blueprint_path = os.path.join(tmpdir, "blueprint.json")
            blueprint = {
                "action_cards": [
                    {
                        "card_id": "CARD-EXCISE-01",
                        "chrom": "OZ059902.1",
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

            pde_bed = os.path.join(tmpdir, "pde_segments.bed")
            with open(pde_bed, "w") as f:
                f.write("OZ059902.1 Oscheius myriophilus genome assembly, chromosome: I\t0\t2100\tputative_core\n")

            # RED outputs headers with full descriptions when given FASTA with description lines
            red_bed = os.path.join(tmpdir, "red_repeats.bed")
            with open(red_bed, "w") as f:
                f.write("OZ059902.1 Oscheius myriophilus genome assembly, chromosome: I\t100\t250\trepeat_A\n")
                f.write("OZ059902.1 Oscheius myriophilus genome assembly, chromosome: I\t1500\t1700\trepeat_B\n")

            outdir = os.path.join(tmpdir, "out")
            cmd = [
                sys.executable,
                str(BIN_DIR / "pde_generate_agp.py"),
                "--blueprint", blueprint_path,
                "--ref", ref_path,
                "--sample", "DF5020_TEST",
                "--outdir", outdir,
                "--pde-segments", pde_bed,
                "--red-bed", red_bed
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}")

            curated_pde = os.path.join(outdir, "DF5020_TEST.curated.pde_segments.bed")
            curated_rep = os.path.join(outdir, "DF5020_TEST.curated.repeats.bed")

            self.assertTrue(os.path.exists(curated_pde))
            self.assertTrue(os.path.exists(curated_rep))

            with open(curated_rep) as f:
                rep_lines = [l.strip().split("\t") for l in f if l.strip()]

            # Must have lifted over both repeats rather than producing an empty (0 byte) file
            self.assertEqual(len(rep_lines), 2)
            self.assertEqual(rep_lines[0], ["OZ059902.1", "100", "250", "repeat_A"])
            # Second repeat was at 1500-1700, after 100 bp gap excision at 1000-1100 -> 1400-1600
            self.assertEqual(rep_lines[1], ["OZ059902.1", "1400", "1600", "repeat_B"])

    def test_internal_gap_and_patch_within_grs_inherits_grs_and_merges(self):
        """Internal gaps and unitig patches between GRS blocks must stay GRS and merge continuously."""
        chrom = "SUPER_X"
        blocks = [
            AGPBlock(type="W", component_id="SUPER_X_0_1000", length=1000, orientation="+", source_seq="A" * 1000),
            AGPBlock(type="N", component_id="gap_1", length=100, orientation=None, source_seq="N" * 100),
            AGPBlock(type="W", component_id="CARD-42_patch", length=500, orientation="+", source_seq="C" * 500),
            AGPBlock(type="N", component_id="gap_2", length=100, orientation=None, source_seq="N" * 100),
            AGPBlock(type="W", component_id="SUPER_X_1200_2000", length=800, orientation="+", source_seq="G" * 800),
        ]
        input_segs = {
            "SUPER_X": [(0, 2000, "putative_GRS")]
        }
        repeats = {}
        raw_segs, _ = assign_chromosome_pde_segments(chrom, blocks, input_segs, repeats)
        merged = merge_curated_pde_segments(raw_segs, sum(b.length for b in blocks))
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0], (0, 2500, "putative_GRS"))

    def test_transition_gap_assigned_to_core_side(self):
        """When an internal gap straddles GRS and Core, the gap must be assigned to the Core side."""
        chrom = "SUPER_1"
        blocks = [
            AGPBlock(type="W", component_id="SUPER_1_0_1000", length=1000, orientation="+", source_seq="A" * 1000),
            AGPBlock(type="N", component_id="transition_gap", length=200, orientation=None, source_seq="N" * 200),
            AGPBlock(type="W", component_id="SUPER_1_1200_3000", length=1800, orientation="+", source_seq="C" * 1800),
        ]
        input_segs = {
            "SUPER_1": [
                (0, 1000, "putative_GRS"),
                (1200, 3000, "putative_core")
            ]
        }
        repeats = {}
        raw_segs, _ = assign_chromosome_pde_segments(chrom, blocks, input_segs, repeats)
        merged = merge_curated_pde_segments(raw_segs, sum(b.length for b in blocks))
        self.assertEqual(len(merged), 2)
        # GRS segment ends at 1000; the 200 bp gap is assigned to core side (1000-1200), merging with 1200-3000
        self.assertEqual(merged[0], (0, 1000, "putative_GRS"))
        self.assertEqual(merged[1], (1000, 3000, "putative_core"))

    def test_telomere_extension_inherits_adjacent_terminus_status(self):
        """A 5' telomere extension attached to a GRS arm must inherit GRS, not core."""
        chrom = "SUPER_1"
        blocks = [
            AGPBlock(type="W", component_id="SUPER_1_5p_telo_ext", length=500, orientation="+", source_seq="TTAGGG" * 83),
            AGPBlock(type="W", component_id="SUPER_1_0_2000", length=2000, orientation="+", source_seq="A" * 2000),
            AGPBlock(type="W", component_id="SUPER_1_3p_telo_ext", length=300, orientation="+", source_seq="CCCTAA" * 50),
        ]
        input_segs = {
            "SUPER_1": [
                (0, 1500, "putative_GRS"),
                (1500, 2000, "putative_core")
            ]
        }
        repeats = {}
        raw_segs, _ = assign_chromosome_pde_segments(chrom, blocks, input_segs, repeats)
        merged = merge_curated_pde_segments(raw_segs, sum(b.length for b in blocks))
        self.assertEqual(len(merged), 2)
        # 5' telo ext (500 bp) inherits GRS and merges with 0-1500 GRS slice -> [0, 2000, putative_GRS]
        self.assertEqual(merged[0], (0, 2000, "putative_GRS"))
        # 3' telo ext (300 bp) inherits Core from adjacent 1500-2000 slice -> [2000, 2800, putative_core]
        self.assertEqual(merged[1], (2000, 2800, "putative_core"))

    def test_rdna_cassette_always_assigned_to_core(self):
        """An rDNA cassette must always be classified as putative_core even if flanked by GRS."""
        chrom = "SUPER_1"
        blocks = [
            AGPBlock(type="W", component_id="SUPER_1_0_1000", length=1000, orientation="+", source_seq="A" * 1000),
            AGPBlock(type="W", component_id="CARD-RDNA-1_rdna_cassette", length=1500, orientation="+", source_seq="G" * 1500),
            AGPBlock(type="W", component_id="SUPER_1_1000_2000", length=1000, orientation="+", source_seq="T" * 1000),
        ]
        input_segs = {
            "SUPER_1": [(0, 2000, "putative_GRS")]
        }
        repeats = {}
        raw_segs, _ = assign_chromosome_pde_segments(chrom, blocks, input_segs, repeats)
        merged = merge_curated_pde_segments(raw_segs, sum(b.length for b in blocks))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged[0], (0, 1000, "putative_GRS"))
        self.assertEqual(merged[1], (1000, 2500, "putative_core"))
        self.assertEqual(merged[2], (2500, 3500, "putative_GRS"))

    def test_scaffold_addition_retains_prior_status(self):
        """An integrated unplaced scaffold must retain its prior status from input segmentation."""
        chrom = "SUPER_2"
        blocks = [
            AGPBlock(type="W", component_id="SUPER_2_0_1000", length=1000, orientation="+", source_seq="A" * 1000),
            AGPBlock(type="W", component_id="scaffold_grs_candidate", length=800, orientation="+", source_seq="C" * 800),
            AGPBlock(type="W", component_id="SUPER_2_1000_2000", length=1000, orientation="+", source_seq="G" * 1000),
        ]
        input_segs = {
            "SUPER_2": [(0, 2000, "putative_core")],
            "scaffold_grs_candidate": [(0, 800, "putative_GRS")]
        }
        repeats = {}
        raw_segs, _ = assign_chromosome_pde_segments(chrom, blocks, input_segs, repeats)
        merged = merge_curated_pde_segments(raw_segs, sum(b.length for b in blocks))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged[0], (0, 1000, "putative_core"))
        self.assertEqual(merged[1], (1000, 1800, "putative_GRS"))
        self.assertEqual(merged[2], (1800, 2800, "putative_core"))


if __name__ == "__main__":
    unittest.main()

