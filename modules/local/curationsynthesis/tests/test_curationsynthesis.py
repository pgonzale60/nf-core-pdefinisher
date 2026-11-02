#!/usr/bin/env python3
"""Unit tests for curationsynthesis telomere unitig surveying logic."""

import os
from pathlib import Path
import sys
import tempfile
import unittest

LIB_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "lib"
sys.path.insert(0, str(LIB_DIR))

from curationsynthesis import (
    build_execution_plan,
    detect_chimeric_terminus,
    load_grs_summary,
    load_terminal_softclips,
    load_terminal_summary,
    reconcile_unplaced_contigs_into_terminus_cards,
    survey_telomere_unitigs,
    synthesize_telomere_cards,
)


class TestCurationsynthesis(unittest.TestCase):

    def test_survey_telomere_unitigs_empty(self):
        res = survey_telomere_unitigs(
            telo_paf=None,
            nuclear_chroms=["I", "II"],
            chrom_lens={"I": 1000000, "II": 1000000},
            unitig_assignment={},
        )
        self.assertEqual(res["total_assembled"], 0)
        self.assertEqual(res["assigned_count"], 0)
        self.assertEqual(res["unassigned_count"], 0)
        self.assertEqual(len(res["termini_lacking_unitigs"]), 4)

    def test_survey_telomere_unitigs_paf_and_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paf_file = os.path.join(tmpdir, "telo.paf")
            with open(paf_file, "w") as f:
                # utg000013l: 34948 bp, maps to Chr III 5' (0 to 12284), MAPQ 60, strand +
                # Overhang at 5' = qstart (22664) if strand +? Wait: qlen=34948, qstart=22664, qend=34948
                f.write("utg000013l\t34948\t22664\t34948\t+\tIII\t10000000\t0\t12284\t12200\t12284\t60\n")
                # ptg000004l: 33055 bp, maps to Chr II 5' (0 to 6762), MAPQ 60, strand +
                f.write("ptg000004l\t33055\t26293\t33055\t+\tII\t10000000\t0\t6762\t6700\t6762\t60\n")
                # utg000099l: 15000 bp, maps internally to Chr I (500000 to 515000), MAPQ 60, strand +
                f.write("utg000099l\t15000\t0\t15000\t+\tI\t10000000\t500000\t515000\t15000\t15000\t60\n")

            telo_junctions = {
                "clusters": [
                    {
                        "cluster_id": "cl_12",
                        "read_count": 18,
                        "mean_telomere_tract_bp": 850,
                        "target_chrom": "III",
                        "target_arm": "5p",
                    },
                    {
                        "cluster_id": "cl_17",
                        "read_count": 5,
                        "mean_telomere_tract_bp": 720,
                        "target_chrom": "IV",
                        "target_arm": "5p",
                    },
                ]
            }

            # utg000013l assigned to III 5p; ptg000004l lost arbitration at II 5p
            unitig_assignment = {
                ("III", "5p"): (22664, "utg000013l", 600, 0),
                ("II", "5p"): (4000, "utg000007l", 400, 0),
            }

            res = survey_telomere_unitigs(
                telo_paf=paf_file,
                nuclear_chroms=["I", "II", "III", "IV"],
                chrom_lens={"I": 10000000, "II": 10000000, "III": 10000000, "IV": 10000000},
                unitig_assignment=unitig_assignment,
                telo_junctions=telo_junctions,
            )

            self.assertEqual(res["total_assembled"], 3)
            self.assertEqual(res["assigned_count"], 1)
            self.assertEqual(res["unassigned_count"], 2)

            u13 = next(u for u in res["unitigs"] if u["unitig_id"] == "utg000013l")
            self.assertEqual(u13["status"], "ASSIGNED")
            self.assertEqual(u13["assigned_to"], "Chr III 5p")
            self.assertIn("18 reads (cl_12)", u13["coverage_str"])
            self.assertEqual(u13["hits"][0]["overhang_bp"], 22664)

            p04 = next(u for u in res["unitigs"] if u["unitig_id"] == "ptg000004l")
            self.assertEqual(p04["status"], "COMPETED_OUT")
            self.assertIn("Competed out at Chr II 5p", p04["disposition"])

            u99 = next(u for u in res["unitigs"] if u["unitig_id"] == "utg000099l")
            self.assertEqual(u99["status"], "INTERNAL_ALIGNMENT")
            self.assertIn("Internal alignment on Chr I:500,000–515,000", u99["disposition"])

            # Check that Chr IV 5p is recognized as lacking unitigs and notes cluster cl_17
            self.assertTrue(any("Chr IV 5p (supported by cluster cl_17 with 5 reads)" in s for s in res["termini_lacking_unitigs"]))

    def test_load_grs_summary_elimination_range(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            grs_tsv = os.path.join(tmpdir, "test.grs_summary.tsv")
            with open(grs_tsv, "w") as f:
                f.write("# SOFTWARE=pde_identify_grs.py\tVERSION=2.0.0\tFINGERPRINT=abc123\tGRS_THRESHOLD=25.00\tSOMATIC_CORE_BASELINE=60.00\n")
                f.write("# ELIMINATION_RANGE\tcanon_grs_bp=1500000\tunplaced_grs_bp=250000\tunplaced_uncertain_bp=80000\tunplaced_core_bp=50000\tmin_elim_bp=1750000\tmax_elim_bp=1830000\test_elim_bp=1790000\n")
                f.write("chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype\n")
                f.write("chr1\t0\t500000\t500000\tputative_GRS\t12.50\tDistal_5p\n")
                f.write("chr1\t500000\t2000000\t1500000\tputative_core\t62.00\tSomatic_Core\n")

            (
                chrom_lens,
                core_segments,
                grs_segments,
                grs_arm_depths,
                local_depth_lookup,
                parsed_grs_threshold,
                elimination_range,
            ) = load_grs_summary(grs_tsv)

            self.assertEqual(parsed_grs_threshold, 25.0)
            self.assertEqual(elimination_range["canon_grs_bp"], 1500000)
            self.assertEqual(elimination_range["unplaced_grs_bp"], 250000)
            self.assertEqual(elimination_range["unplaced_uncertain_bp"], 80000)
            self.assertEqual(elimination_range["unplaced_core_bp"], 50000)
            self.assertEqual(elimination_range["min_elim_bp"], 1750000)
            self.assertEqual(elimination_range["max_elim_bp"], 1830000)
            self.assertEqual(elimination_range["est_elim_bp"], 1790000)

    def test_detect_chimeric_terminus_guardrails(self):
        # 3p Arm Test
        chrom = "chr1"
        chrom_len = 1000000
        policy = {
            "telomere": {
                "min_repeats": 15,
                "chimeric_terminus_min_dist_bp": 1000,
                "chimeric_terminus_max_dist_bp": 50000,
                "subtelomeric_min_overlap_bp": 500,
                "min_chimeric_clip_support": 3,
                "chimeric_cluster_window_bp": 50,
            },
            "mapping_and_anchors": {
                "mapq_anchored": 30,
                "mapq_evidence_floor": 20,
                "anchor_min_length_bp": 1000,
            },
            "genomic_architecture": {"germline_buffer_bp": 100},
        }

        # 1. Reject when < 3 clips (1 or 2 reads)
        cands_3p = [
            {
                "qname": f"read_3p_{i}",
                "align_start": 965000,
                "align_end": 980000 + (i * 10),
                "clip_len": 500,
                "overhang_bp": 0,
                "telo_cnt": 20,
                "mapq": 60,
                "collinear_spanning_reads": 0,
            }
            for i in range(2)
        ]
        self.assertIsNone(detect_chimeric_terminus(chrom, "3p", chrom_len, cands_3p, 0, [], policy))

        # 2. Accept when >= 3 clips clustered near breakpoint and 0 collinear spanning reads
        cands_3p.append(
            {
                "qname": "read_3p_2",
                "align_start": 965000,
                "align_end": 980005,
                "clip_len": 520,
                "overhang_bp": 0,
                "telo_cnt": 25,
                "mapq": 60,
                "collinear_spanning_reads": 0,
            }
        )
        res = detect_chimeric_terminus(chrom, "3p", chrom_len, cands_3p, 0, [], policy)
        self.assertIsNotNone(res)
        self.assertEqual(res["chimeric_start"], 980005)
        self.assertEqual(res["chimeric_end"], chrom_len)
        self.assertEqual(res["authorized_read_qname"], "read_3p_2")
        self.assertEqual(res["supporting_read_count"], 3)
        self.assertEqual(len(res["supporting_reads"]), 3)

        # 3. Reject when collinear_spanning_reads > 0 on candidate record
        cands_with_collinear = [dict(c) for c in cands_3p]
        cands_with_collinear[0]["collinear_spanning_reads"] = 2
        self.assertIsNone(detect_chimeric_terminus(chrom, "3p", chrom_len, cands_with_collinear, 0, [], policy))

        # 4. Reject when an unclipped spanning read is in clips
        cands_with_clip_span = list(cands_3p) + [
            {
                "qname": "unclipped_span",
                "align_start": 975000,
                "align_end": 985000,
                "clip_len": 0,
                "overhang_bp": 0,
                "telo_cnt": 0,
                "mapq": 60,
                "collinear_spanning_reads": 0,
            }
        ]
        self.assertIsNone(detect_chimeric_terminus(chrom, "3p", chrom_len, cands_with_clip_span, 0, [], policy))

        # 5p Arm Test
        cands_5p = [
            {
                "qname": f"read_5p_{i}",
                "align_start": 20000 + (i * 10),
                "align_end": 35000,
                "clip_len": 500,
                "overhang_bp": 0,
                "telo_cnt": 20,
                "mapq": 60,
                "collinear_spanning_reads": 0,
            }
            for i in range(3)
        ]
        res_5p = detect_chimeric_terminus(chrom, "5p", chrom_len, cands_5p, 0, [], policy)
        self.assertIsNotNone(res_5p)
        self.assertEqual(res_5p["chimeric_start"], 0)
        self.assertEqual(res_5p["chimeric_end"], 20020)
        self.assertEqual(res_5p["supporting_read_count"], 3)

        # Reject 5p when collinear_spanning_reads > 0
        cands_5p_coll = [dict(c) for c in cands_5p]
        cands_5p_coll[1]["collinear_spanning_reads"] = 1
        self.assertIsNone(detect_chimeric_terminus(chrom, "5p", chrom_len, cands_5p_coll, 0, [], policy))

    def test_level_4_unitig_priority_over_longer_read(self):
        # Even if a single read overhang is significantly longer (+2500 bp vs +500 bp unitig),
        # assembly-first hierarchy must strictly prefer the unitig.
        chrom_lens = {"chr1": 1000000}
        nuclear_chroms = ["chr1"]
        unitig_extensions = {
            ("chr1", "3p"): [
                (500, "utg000010l", 600, 0, 0, 10000, "+", 10000, 990000, 1000000)
            ]
        }
        single_read_ext = {
            ("chr1", "3p"): (2500, 5, "m64000_long_read", 25, ["r1", "r2", "r3", "r4", "r5"])
        }

        # Case 1: ALREADY_CAPPED
        cards_capped = synthesize_telomere_cards(
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            grs_segments={},
            core_segments={},
            delfies_points={},
            rdna_intervals={},
            baseline_status={("chr1", "3p"): "ALREADY_CAPPED"},
            unitig_extensions=unitig_extensions,
            single_read_ext=single_read_ext,
            tip_tract_bp={("chr1", "3p"): 500},
        )
        card_3p = next(c for c in cards_capped if c["chrom"] == "chr1" and c["terminus"] == "3p")
        self.assertEqual(card_3p["action"], "EXTEND_TELOMERE_TRACT")
        self.assertEqual(card_3p["authorized_unitig_id"], "utg000010l")
        self.assertIsNone(card_3p["authorized_read_qname"])
        self.assertIsInstance(card_3p["evidence"]["deep_spanning_overhang_bp"], int)
        self.assertEqual(card_3p["evidence"]["deep_spanning_overhang_bp"], 2500)

        # Case 2: TELOMERE_IN_OVERHANG
        cards_overhang = synthesize_telomere_cards(
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            grs_segments={},
            core_segments={},
            delfies_points={},
            rdna_intervals={},
            baseline_status={("chr1", "3p"): "TELOMERE_IN_OVERHANG"},
            unitig_extensions=unitig_extensions,
            single_read_ext=single_read_ext,
            tip_tract_bp={("chr1", "3p"): 0},
        )
        card_oh = next(c for c in cards_overhang if c["chrom"] == "chr1" and c["terminus"] == "3p")
        self.assertEqual(card_oh["action"], "EXTEND_FROM_OVERHANG")
        self.assertEqual(card_oh["authorized_unitig_id"], "utg000010l")
        self.assertIsNone(card_oh["authorized_read_qname"])
        self.assertIsInstance(card_oh["evidence"]["deep_spanning_overhang_bp"], int)
        self.assertEqual(card_oh["evidence"]["deep_spanning_overhang_bp"], 2500)
        self.assertIsInstance(card_oh["evidence"]["terminal_tract_bp"], int)
        self.assertEqual(card_oh["evidence"]["terminal_tract_bp"], 0)

    def test_level_5c_cap_with_read_multi_read_support(self):
        chrom_lens = {"chr1": 1000000}
        nuclear_chroms = ["chr1"]

        # Case A: < 3 supporting reads (e.g. 2 reads) -> must NOT trigger RESOLVABLE_BY_READ / CAP_WITH_READ
        single_read_insufficient = {
            ("chr1", "3p"): (500, 2, "read_capper_2x", 20, ["r1", "r2"])
        }
        cards_insufficient = synthesize_telomere_cards(
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            grs_segments={},
            core_segments={},
            delfies_points={},
            rdna_intervals={},
            baseline_status={("chr1", "3p"): "NO_CLIPS"},
            unitig_extensions={},
            single_read_ext=single_read_insufficient,
            tip_tract_bp={("chr1", "3p"): 0},
        )
        card_rej = next(c for c in cards_insufficient if c["chrom"] == "chr1" and c["terminus"] == "3p")
        self.assertNotEqual(card_rej.get("mechanism"), "CAP_WITH_READ")
        self.assertEqual(card_rej["action"], "EXTEND_TO_TELOMERE")
        self.assertEqual(card_rej["classification"], "UNCAPPED_SUBTELOMERE")

        # Case B: >= 3 supporting reads (e.g. 3 reads) -> successfully authorizes CAP_WITH_READ
        single_read_sufficient = {
            ("chr1", "3p"): (500, 3, "read_capper_3x", 20, ["r1", "r2", "r3"])
        }
        cards_sufficient = synthesize_telomere_cards(
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            grs_segments={},
            core_segments={},
            delfies_points={},
            rdna_intervals={},
            baseline_status={("chr1", "3p"): "NO_CLIPS"},
            unitig_extensions={},
            single_read_ext=single_read_sufficient,
            tip_tract_bp={("chr1", "3p"): 0},
        )
        card_acc = next(c for c in cards_sufficient if c["chrom"] == "chr1" and c["terminus"] == "3p")
        self.assertEqual(card_acc["classification"], "RESOLVABLE_BY_READ")
        self.assertEqual(card_acc["mechanism"], "CAP_WITH_READ")
        self.assertEqual(card_acc["authorized_read_qname"], "read_capper_3x")
        self.assertEqual(card_acc["evidence"]["read_count"], 3)
        self.assertEqual(card_acc["evidence"]["supporting_reads"], ["r1", "r2", "r3"])
        self.assertIsInstance(card_acc["evidence"]["read_extension_bp"], int)
        self.assertEqual(card_acc["evidence"]["read_extension_bp"], 500)

    def test_load_terminal_softclips_collinear_and_multi_read(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = os.path.join(tmpdir, "test.terminal_softclips.tsv")
            with open(tsv_path, "w") as f:
                f.write("read_id\tchrom\tside\talign_start\talign_end\tclip_len\toverhang_bp\ttelo_count\tmapq\tcollinear_spanning_reads\n")
                # 10-column record with collinear_spanning_reads = 4
                f.write("read1\tchr1\t3p\t980000\t995000\t500\t200\t25\t60\t4\n")
                # 9-column legacy record (backward compatibility)
                f.write("read2\tchr1\t3p\t980000\t995000\t450\t180\t20\t50\n")

            raw_clips, single_ext = load_terminal_softclips(tsv_path)
            clips_3p = raw_clips[("chr1", "3p")]
            self.assertEqual(len(clips_3p), 2)
            self.assertEqual(clips_3p[0]["collinear_spanning_reads"], 4)
            self.assertEqual(clips_3p[1]["collinear_spanning_reads"], 0)

            self.assertIn(("chr1", "3p"), single_ext)
            ext_info = single_ext[("chr1", "3p")]
            # ext_info format: (overhang_bp, read_count, best_qname, best_telo_cnt, [qnames])
            self.assertEqual(ext_info[0], 200)
            self.assertEqual(ext_info[1], 2)
            self.assertEqual(ext_info[2], "read1")
            self.assertEqual(ext_info[3], 25)
            self.assertEqual(ext_info[4], ["read1", "read2"])

    def test_execution_plan_assembly_first_priority(self):
        """Verifies curationsynthesis emits execution_plan in strict assembly-first priority order."""
        # 1. Major rDNA array precedence -> DEFER_TO_RDNA_RESOLVER
        rdna_card = {
            "classification": "TERMINUS_RDNA_ARRAY",
            "action": "TERMINUS_RDNA_ARRAY",
            "chrom": "chr1",
            "terminus": "5p",
        }
        plan_rdna = build_execution_plan(rdna_card)
        self.assertEqual(plan_rdna, [{"step": "DEFER_TO_RDNA_RESOLVER"}])

        # 2. Chimeric terminus pruning -> PRUNE_CHIMERIC_TERMINUS
        chimera_card = {
            "classification": "RESOLVABLE_BY_READ",
            "action": "PRUNE_CHIMERIC_TERMINUS",
            "chrom": "chr1",
            "terminus": "3p",
            "start": 980000,
            "end": 1000000,
            "effective_chrom_len": 980000,
            "authorized_read_qname": "read_chimera",
        }
        plan_chimera = build_execution_plan(chimera_card)
        self.assertEqual(len(plan_chimera), 1)
        self.assertEqual(plan_chimera[0]["step"], "PRUNE_CHIMERIC_TERMINUS")
        self.assertEqual(plan_chimera[0]["junction_coord"], 980000)
        self.assertEqual(plan_chimera[0]["authorized_read_qname"], "read_chimera")

        # 3. Plain ALREADY_CAPPED without deeper extension -> ALREADY_CAPPED
        capped_card = {
            "classification": "ALREADY_CAPPED",
            "action": "SKIP_WALKER",
            "chrom": "chr1",
            "terminus": "5p",
        }
        plan_capped = build_execution_plan(capped_card)
        self.assertEqual(plan_capped, [{"step": "ALREADY_CAPPED"}])

        # 4. Terminus with unitig, satellite, unplaced contig, and read evidence -> strict priority
        complex_card = {
            "classification": "RESOLVABLE_BY_UNITIG",
            "action": "EXTEND_TO_TELOMERE",
            "mechanism": "CAP_WITH_UNITIG",
            "authorized_unitig_id": "utg000010l",
            "authorized_read_qname": "read_backup",
            "chrom": "chr1",
            "terminus": "5p",
            "evidence": {
                "unitig_extension_region": {"q_start": 0, "q_end": 5000, "strand": "+", "unitig_len": 10000},
                "distal_cluster_id": "cl_distal_1",
                "associated_scaffolds": ["scaff_repeat_1"],
            },
        }
        nom = {
            "scaffold": "scaff_unplaced_1",
            "overhang_bp": 3500,
            "strand": "+",
            "ext_qstart": 0,
            "ext_qend": 3500,
        }
        plan_complex = build_execution_plan(complex_card, unplaced_nomination=nom)
        step_names = [s["step"] for s in plan_complex]
        expected_order = [
            "CAP_WITH_UNITIG",
            "INTEGRATE_UNPLACED_CONTIG",
            "SCAFFOLD_REPEAT_TO_TELOMERE",
            "ALIGN_SEED_TO_UNITIGS",
            "CAP_WITH_READ",
            "ITERATIVE_READ_WALKER",
        ]
        self.assertEqual(step_names, expected_order)
        self.assertEqual(plan_complex[0]["unitig_id"], "utg000010l")
        self.assertEqual(plan_complex[1]["scaffold"], "scaff_unplaced_1")
        self.assertEqual(plan_complex[2]["cluster_id"], "cl_distal_1")
        self.assertEqual(plan_complex[4]["authorized_read_qname"], "read_backup")

    def test_reconcile_unplaced_contigs_into_terminus_cards(self):
        """Verifies reconcile_unplaced_contigs_into_terminus_cards incorporates nominations in assembly-first order."""
        telo_cards = [
            {
                "card_id": "CARD-TELO-01",
                "chrom": "chr1",
                "terminus": "5p",
                "action": "EXTEND_TO_TELOMERE",
                "classification": "RESOLVABLE_BY_UNITIG",
                "authorized_unitig_id": "utg001",
                "evidence": {},
                "execution_plan": [
                    {"step": "CAP_WITH_UNITIG", "unitig_id": "utg001"},
                    {"step": "ALIGN_SEED_TO_UNITIGS"},
                    {"step": "ITERATIVE_READ_WALKER"},
                ],
            },
            {
                "card_id": "CARD-TELO-02",
                "chrom": "chr1",
                "terminus": "3p",
                "action": "EXTEND_TO_TELOMERE",
                "classification": "UNCAPPED_SUBTELOMERE",
                "evidence": {},
                "execution_plan": [
                    {"step": "ALIGN_SEED_TO_UNITIGS"},
                    {"step": "ITERATIVE_READ_WALKER"},
                ],
            },
        ]

        scaffold_cards = [
            {
                "card_id": "CARD-SCAFF-INT-01",
                "chrom": "chr1",
                "terminus": "5p",
                "scaffold": "unplaced_scaff_A",
                "length": 15000,
                "ext_qstart": 0,
                "ext_qend": 8000,
                "orientation": "+",
                "action": "INTEGRATE_INTO_CHROMOSOME",
                "evidence": {"extension_beyond_terminus": "8,000 bp"},
            },
            {
                "card_id": "CARD-SCAFF-INT-02",
                "chrom": "chr1",
                "terminus": "3p",
                "scaffold": "unplaced_scaff_B",
                "length": 10000,
                "ext_qstart": 1000,
                "ext_qend": 5000,
                "orientation": "-",
                "action": "INTEGRATE_INTO_CHROMOSOME",
                "evidence": {"extension_beyond_terminus": "4,000 bp"},
            },
        ]

        reconcile_unplaced_contigs_into_terminus_cards(telo_cards, scaffold_cards=scaffold_cards)

        # chr1 5p: nomination reconciled into evidence and placed after CAP_WITH_UNITIG
        c5 = telo_cards[0]
        self.assertIn("nominated_unplaced_contig", c5["evidence"])
        self.assertEqual(c5["evidence"]["nominated_unplaced_contig"]["scaffold"], "unplaced_scaff_A")
        plan_5p = c5["execution_plan"]
        self.assertEqual([s["step"] for s in plan_5p], [
            "CAP_WITH_UNITIG",
            "INTEGRATE_UNPLACED_CONTIG",
            "ALIGN_SEED_TO_UNITIGS",
            "ITERATIVE_READ_WALKER",
        ])
        self.assertEqual(plan_5p[1]["scaffold"], "unplaced_scaff_A")

        # chr1 3p: nomination placed at index 0 (no unitig present)
        c3 = telo_cards[1]
        self.assertIn("nominated_unplaced_contig", c3["evidence"])
        self.assertEqual(c3["evidence"]["nominated_unplaced_contig"]["scaffold"], "unplaced_scaff_B")
        plan_3p = c3["execution_plan"]
        self.assertEqual([s["step"] for s in plan_3p], [
            "INTEGRATE_UNPLACED_CONTIG",
            "ALIGN_SEED_TO_UNITIGS",
            "ITERATIVE_READ_WALKER",
        ])
        self.assertEqual(plan_3p[0]["scaffold"], "unplaced_scaff_B")

    def test_level_5a_unitig_capping_emits_raw_integers_in_evidence(self):
        """Verifies that Level 5a Phase 1 Unitig Capping emits raw integers in evidence_dict."""
        chrom_lens = {"chr1": 1000000}
        nuclear_chroms = ["chr1"]
        # unitig_extensions format: list of (overhang_bp, unitig_id, telomeric_tract_bp, subtelo_gap_bp, q_start, q_end, strand, unitig_len)
        unitig_extensions = {
            ("chr1", "5p"): [(12000, "utg000005l", 850, 400, 0, 12000, "+", 15000)]
        }
        cards = synthesize_telomere_cards(
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            grs_segments={},
            core_segments={},
            delfies_points={},
            rdna_intervals={},
            baseline_status={("chr1", "5p"): "UNCAPPED"},
            unitig_extensions=unitig_extensions,
            single_read_ext={},
            tip_tract_bp={("chr1", "5p"): 0},
        )
        card = next(c for c in cards if c["chrom"] == "chr1" and c["terminus"] == "5p")
        self.assertEqual(card["classification"], "RESOLVABLE_BY_UNITIG")
        evidence = card["evidence"]
        self.assertEqual(evidence["unitig_id"], "utg000005l")
        self.assertIsInstance(evidence["unitig_extension_bp"], int)
        self.assertEqual(evidence["unitig_extension_bp"], 12000)
        self.assertIsInstance(evidence["subtelomeric_grs_bp"], int)
        self.assertEqual(evidence["subtelomeric_grs_bp"], 400)
        self.assertIsInstance(evidence["telomeric_tract_bp"], int)
        self.assertEqual(evidence["telomeric_tract_bp"], 850)

    def test_load_terminal_summary_scaffold_with_underscores(self):
        """Verifies that scaffold names with underscores (e.g. scaffold_1_5p) are split correctly using rsplit."""
        with tempfile.NamedTemporaryFile("w", delete=False) as tf:
            tf.write("chrom_end\tchrom_length\tn_softclipped_reads\tn_telomeric_reads\tmax_telo_count\tproposed_status\tterminal_tract_bp\n")
            tf.write("scaffold_1_5p\t1000000\t15\t12\t45\tALREADY_CAPPED\t250\n")
            tf.write("chr_2_sub_3p\t2000000\t8\t6\t30\tUNCAPPED\t180\n")
            tf_path = tf.name

        try:
            baseline_status, tip_tract_bp = load_terminal_summary(tf_path)
            self.assertIn(("scaffold_1", "5p"), baseline_status)
            self.assertEqual(baseline_status[("scaffold_1", "5p")], "ALREADY_CAPPED")
            self.assertEqual(tip_tract_bp[("scaffold_1", "5p")], 250)

            self.assertIn(("chr_2_sub", "3p"), baseline_status)
            self.assertEqual(baseline_status[("chr_2_sub", "3p")], "UNCAPPED")
            self.assertEqual(tip_tract_bp[("chr_2_sub", "3p")], 180)

            # Verify that synthesize_telomere_cards correctly consumes baseline_status for scaffold_1
            cards = synthesize_telomere_cards(
                nuclear_chroms=["scaffold_1"],
                chrom_lens={"scaffold_1": 1000000},
                grs_segments={},
                core_segments={},
                delfies_points={},
                rdna_intervals={},
                baseline_status=baseline_status,
                unitig_extensions={},
                single_read_ext={},
                tip_tract_bp=tip_tract_bp,
            )
            scaff_5p_card = next(c for c in cards if c["chrom"] == "scaffold_1" and c["terminus"] == "5p")
            self.assertEqual(scaff_5p_card["classification"], "ALREADY_CAPPED")
        finally:
            os.unlink(tf_path)


if __name__ == "__main__":
    unittest.main()
