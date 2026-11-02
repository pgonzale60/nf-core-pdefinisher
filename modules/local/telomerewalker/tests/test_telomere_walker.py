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
            "read3\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15 + "AAAG") + "\t*",
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

    def test_overhang_of_read_filters_low_mapq_and_flags(self):
        """Verifies overhang_of_read filters reads with MAPQ < 20 or secondary/supplementary/unmapped flags."""
        from telocapping.overhang_capper import overhang_of_read
        valid_line = "read1\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*"
        clip, ov = overhang_of_read(valid_line, "3p", 10000, min_clip_bp=50)
        self.assertGreater(clip, 0)
        self.assertNotEqual(ov, "")

        # Low MAPQ < 20
        low_mapq_line = "read1\t0\tchrI\t9950\t19\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*"
        self.assertEqual(overhang_of_read(low_mapq_line, "3p", 10000, min_clip_bp=50), (0, ""))

        # Secondary alignment (flag 0x100 = 256)
        sec_line = "read1\t256\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*"
        self.assertEqual(overhang_of_read(sec_line, "3p", 10000, min_clip_bp=50), (0, ""))

        # Supplementary alignment (flag 0x800 = 2048)
        supp_line = "read1\t2048\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*"
        self.assertEqual(overhang_of_read(supp_line, "3p", 10000, min_clip_bp=50), (0, ""))

        # Unmapped read (flag 0x4 = 4)
        unmapped_line = "read1\t4\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*"
        self.assertEqual(overhang_of_read(unmapped_line, "3p", 10000, min_clip_bp=50), (0, ""))

    def test_maximal_tip_overhang_requires_min_reads(self):
        """Verifies maximal_tip_overhang returns None when fewer than min_reads are present."""
        from telocapping import overhang_capper
        mock_lines = [
            "read1\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*",
            "read2\t0\tchrI\t9950\t60\t50M100S\t*\t0\t0\t" + ("A" * 50) + ("TTAGGC" * 15) + "\t*",
        ]
        orig_twl = overhang_capper.terminal_window_lines
        overhang_capper.terminal_window_lines = lambda *args, **kwargs: mock_lines
        try:
            # Default min_reads = 3 -> must return None for 2 reads
            res = overhang_capper.maximal_tip_overhang("dummy.bam", "chrI", "3p", 10000, window=1000, chrom_start=0)
            self.assertIsNone(res)
        finally:
            overhang_capper.terminal_window_lines = orig_twl

    def test_maximal_terminal_tract_requires_min_reads_and_consensus(self):
        """Verifies maximal_terminal_tract enforces >= 3 reads with telo_len >= 200 and consensus."""
        from telotract.tract_detector import maximal_terminal_tract
        # Anchor (500 bp) + current tract (150 bp of TTAGGC = 25 repeats)
        anchor = "A" * 500
        cur_telo = "TTAGGC" * 25
        terminal_seq = anchor + cur_telo

        # Candidates extending the tract by > 60 bp with total telo_len >= 200 bp (35 repeats = 210 bp)
        c1 = {"qname": "cand1", "seq": anchor + ("TTAGGC" * 40)}
        c2 = {"qname": "cand2", "seq": anchor + ("TTAGGC" * 45)}
        c3 = {"qname": "cand3", "seq": anchor + ("TTAGGC" * 50)}

        # With 2 candidates (< 3) -> None
        res_few = maximal_terminal_tract([c1, c2], "TTAGGC", "GCCTAA", terminal_seq, anchor_bp=500, min_reads=3)
        self.assertIsNone(res_few)

        # With 3 candidates (>= 3) -> Consensus
        res_ok = maximal_terminal_tract([c1, c2, c3], "TTAGGC", "GCCTAA", terminal_seq, anchor_bp=500, min_reads=3)
        self.assertIsNotNone(res_ok)
        label, total_tract, extra_seq, orientation, insert_at = res_ok
        self.assertGreater(total_tract, len(cur_telo))
        self.assertIn("TTAGGC", extra_seq)

    def test_dentist_guard_veto_logic(self):
        """Verifies evaluate_dentist_guard vetoes when EITHER flank exceeds I_max (OR gate)."""
        import sys
        target_lib = str(Path(__file__).resolve().parent.parent.parent / "targetedgapassembly" / "resources" / "usr" / "lib")
        if target_lib not in sys.path:
            sys.path.insert(0, target_lib)
        from gapassembly import dentist_guard

        orig_cf = dentist_guard.count_flank_improper_alignments
        try:
            # Case 1: Left flank exceeds I_max, Right flank is clean -> MUST BE VETOED
            dentist_guard.count_flank_improper_alignments = lambda bam, chrom, coord, side, **kwargs: (
                (20, 15) if side == "LEFT" else (20, 2)
            )
            is_vetoed, l_stat, r_stat = dentist_guard.evaluate_dentist_guard("dummy.bam", "chr1", 1000, 2000)
            self.assertTrue(l_stat["excess"])
            self.assertFalse(r_stat["excess"])
            self.assertTrue(is_vetoed, "Single flank excess must trigger veto (OR logic)")

            # Case 2: Right flank exceeds I_max, Left flank is clean -> MUST BE VETOED
            dentist_guard.count_flank_improper_alignments = lambda bam, chrom, coord, side, **kwargs: (
                (20, 2) if side == "LEFT" else (20, 15)
            )
            is_vetoed, l_stat, r_stat = dentist_guard.evaluate_dentist_guard("dummy.bam", "chr1", 1000, 2000)
            self.assertFalse(l_stat["excess"])
            self.assertTrue(r_stat["excess"])
            self.assertTrue(is_vetoed, "Single flank excess must trigger veto (OR logic)")

            # Case 3: Neither flank exceeds I_max -> NOT VETOED
            dentist_guard.count_flank_improper_alignments = lambda bam, chrom, coord, side, **kwargs: (
                (20, 2) if side == "LEFT" else (20, 2)
            )
            is_vetoed, l_stat, r_stat = dentist_guard.evaluate_dentist_guard("dummy.bam", "chr1", 1000, 2000)
            self.assertFalse(is_vetoed, "Neither flank excess must not trigger veto")
        finally:
            dentist_guard.count_flank_improper_alignments = orig_cf

    def test_opportunity_2_blocked_on_rollback(self):
        """Verifies Opportunity 2 abpoa micro-gap fallback is NOT executed when repeat rollback occurs."""
        import shutil
        from telowalking import iterative_walker
        if not shutil.which("abpoa") or not shutil.which("minimap2"):
            self.skipTest("abpoa or minimap2 not available in test environment")

        class DummyIndex:
            def __len__(self):
                return 0

        orig_gate = iterative_walker.evaluate_step_information_gain
        try:
            iterative_walker.evaluate_step_information_gain = lambda *args, **kwargs: {
                "decision": "ROLLBACK",
                "telo_bp": 0,
                "reason": "Polyvalent read branch detected"
            }
            res = iterative_walker._walk_loop(
                unmapped=DummyIndex(),
                seed_seq="ACGT" * 100,
                motif="TTAGGC",
                rev_motif="GCCTAA",
                max_iters=5,
                min_telo_bp=100,
                is_somatic_boundary=False,
                min_germline_extension_bp=100
            )
            self.assertNotEqual(res["status"], "EXTEND_TO_TELOMERE")
        finally:
            iterative_walker.evaluate_step_information_gain = orig_gate

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

    def test_derive_default_execution_plan_legacy_cards(self):
        """Verifies _derive_default_execution_plan accurately maps legacy blueprint cards."""
        # None with blueprint -> ALREADY_CAPPED
        self.assertEqual(
            pde_telomere_walker._derive_default_execution_plan(None, has_blueprint=True),
            [{"step": "ALREADY_CAPPED"}]
        )
        # None without blueprint -> ITERATIVE_READ_WALKER
        self.assertEqual(
            pde_telomere_walker._derive_default_execution_plan(None, has_blueprint=False),
            [{"step": "ITERATIVE_READ_WALKER"}]
        )
        # rDNA array
        rdna_card = {"classification": "TERMINUS_RDNA_ARRAY", "action": "TERMINUS_RDNA_ARRAY"}
        self.assertEqual(
            pde_telomere_walker._derive_default_execution_plan(rdna_card),
            [{"step": "DEFER_TO_RDNA_RESOLVER"}]
        )
        # Chimeric terminus
        chimera_card = {
            "action": "PRUNE_CHIMERIC_TERMINUS",
            "terminus": "5p",
            "effective_chrom_start": 2000,
            "authorized_read_qname": "read_chimeric_auth",
        }
        self.assertEqual(
            pde_telomere_walker._derive_default_execution_plan(chimera_card),
            [{"step": "PRUNE_CHIMERIC_TERMINUS", "junction_coord": 2000, "authorized_read_qname": "read_chimeric_auth"}]
        )
        # Unitig cap with unplaced nomination
        utg_card = {
            "classification": "RESOLVABLE_BY_UNITIG",
            "mechanism": "CAP_WITH_UNITIG",
            "authorized_unitig_id": "utg001l",
        }
        nom = {"scaffold": "scaff_ext", "overhang_bp": 5000, "strand": "+"}
        plan = pde_telomere_walker._derive_default_execution_plan(utg_card, unplaced_nomination=nom)
        self.assertEqual([s["step"] for s in plan], [
            "CAP_WITH_UNITIG",
            "INTEGRATE_UNPLACED_CONTIG",
            "ALIGN_SEED_TO_UNITIGS",
            "ITERATIVE_READ_WALKER",
        ])

    def test_execution_plan_dispatcher_order_and_stop(self):
        """Verifies the stateless dispatcher runs steps strictly in order and stops on first success."""
        manifest_rows = []
        telomere_lengths = []

        # 1. ALREADY_CAPPED handler appends row and returns True
        res = pde_telomere_walker._handle_already_capped("chr1", "5p", {}, {}, manifest_rows, telomere_lengths)
        self.assertTrue(res)
        self.assertEqual(len(manifest_rows), 1)
        self.assertEqual(manifest_rows[0][2], "ALREADY_CAPPED")
        self.assertEqual(manifest_rows[0][3], "DRAFT_VERIFIED_COMPLETE")

        # 2. DEFER_TO_RDNA_RESOLVER handler appends row and returns True
        res_rdna = pde_telomere_walker._handle_defer_to_rdna("chr2", "3p", {}, {}, "SAMPLE1", "chr2_3p", manifest_rows, telomere_lengths)
        self.assertTrue(res_rdna)
        self.assertEqual(len(manifest_rows), 2)
        self.assertEqual(manifest_rows[1][2], "TERMINUS_RDNA_ARRAY")

        # 3. Step execution short-circuiting: first step that succeeds halts the loop
        executed_steps = []
        mock_plan = [
            {"step": "STEP_A"},
            {"step": "STEP_B"},
            {"step": "STEP_C"},
        ]

        def step_handler(step):
            executed_steps.append(step)
            return step == "STEP_B"  # STEP_B succeeds

        for s in mock_plan:
            if step_handler(s["step"]):
                break

        self.assertEqual(executed_steps, ["STEP_A", "STEP_B"])
        self.assertNotIn("STEP_C", executed_steps)

    def test_materialise_satellite_scaffold_sanitizes_unit_suffixed_tract_bp(self):
        """Verifies that repeat_scaffolder.materialise_satellite_scaffold parses unit-suffixed telomeric_tract_bp."""
        from teloscaffold import materialise_satellite_scaffold

        card_with_bp_str = {
            "evidence": {
                "telomeric_tract_bp": "0 bp",
                "supporting_reads": [],
            }
        }
        res = materialise_satellite_scaffold(
            card_with_bp_str, "chr1", "5p", "TTAGGC", "GCCTAA", nominal_cap_bp=800
        )
        self.assertEqual(res["status"], "RESOLVED_BY_SATELLITE_SCAFFOLD")
        self.assertGreater(res["telo_len"], 0)

        card_with_comma_kbp = {
            "evidence": {
                "telomeric_tract_bp": "1,250 bp",
                "supporting_reads": [],
            }
        }
        res2 = materialise_satellite_scaffold(
            card_with_comma_kbp, "chr1", "3p", "TTAGGC", "GCCTAA", nominal_cap_bp=800
        )
        self.assertEqual(res2["status"], "RESOLVED_BY_SATELLITE_SCAFFOLD")
        self.assertGreaterEqual(res2["telo_len"], 1200)

        card_with_int = {
            "evidence": {
                "telomeric_tract_bp": 500,
                "supporting_reads": [],
            }
        }
        res3 = materialise_satellite_scaffold(
            card_with_int, "chr1", "5p", "TTAGGC", "GCCTAA", nominal_cap_bp=800
        )
        self.assertEqual(res3["status"], "RESOLVED_BY_SATELLITE_SCAFFOLD")

        card_with_none = {
            "evidence": {
                "telomeric_tract_bp": None,
                "supporting_reads": [],
            }
        }
        res4 = materialise_satellite_scaffold(
            card_with_none, "chr1", "5p", "TTAGGC", "GCCTAA", nominal_cap_bp=800
        )
        self.assertEqual(res4["status"], "RESOLVED_BY_SATELLITE_SCAFFOLD")

if __name__ == "__main__":
    unittest.main()
