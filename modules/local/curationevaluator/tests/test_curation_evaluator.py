"""Unit tests for curationevaluator helper functions and repeat calculation."""

import tempfile
import unittest
from pathlib import Path

import sys
# Add usr/bin to path
bin_dir = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(bin_dir))

from pde_evaluate_curation import (
    compute_red_repeat_pct,
    analyze_terminal_repeats,
    classify_unresolved_gap,
    render_telomere_unitig_table,
)


class TestCurationEvaluator(unittest.TestCase):

    def test_compute_red_repeat_pct_no_file(self):
        self.assertIsNone(compute_red_repeat_pct(None, "I", 0, 5000))
        self.assertIsNone(compute_red_repeat_pct("/nonexistent/path.bed", "I", 0, 5000))

    def test_compute_red_repeat_pct_disjoint_and_overlapping(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as tf:
            # Chromosome I: repeat from 0 to 1000, 500 to 1500 (overlapping), and 2000 to 3000
            tf.write("I\t0\t1000\n")
            tf.write("I\t500\t1500\n")
            tf.write("I\t2000\t3000\n")
            tf.write("II\t0\t5000\n")
            tf_path = tf.name

        try:
            # For Chromosome I in window [0, 5000]:
            # Merged intervals: [0, 1500] (1500 bp) and [2000, 3000] (1000 bp) = 2500 bp / 5000 bp = 50.0%
            pct = compute_red_repeat_pct(tf_path, "I", 0, 5000)
            self.assertAlmostEqual(pct, 50.0)

            # For Chromosome II in window [0, 5000]:
            # [0, 5000] = 5000 bp / 5000 bp = 100.0%
            pct_ii = compute_red_repeat_pct(tf_path, "II", 0, 5000)
            self.assertAlmostEqual(pct_ii, 100.0)

            # For Chromosome III (absent from BED):
            pct_iii = compute_red_repeat_pct(tf_path, "III", 0, 5000)
            self.assertEqual(pct_iii, 0.0)
        finally:
            Path(tf_path).unlink(missing_ok=True)

    def test_analyze_terminal_repeats_tip_stretch(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as tf:
            # Chromosome SUPER_2: massive repeat from 0 to 763004 bp
            tf.write("SUPER_2\t0\t763004\n")
            # Chromosome SUPER_3: repeat at 3' tip: length 10,000,000, repeat from 9,900,000 to 10,000,000
            tf.write("SUPER_3\t9900000\t10000000\n")
            # Chromosome SUPER_4: repeat not at tip, starts at 500 to 1500
            tf.write("SUPER_4\t500\t1500\n")
            tf_path = tf.name

        try:
            # 5' tip test on SUPER_2:
            res_s2 = analyze_terminal_repeats(tf_path, "SUPER_2", "5p", 25000000, window_bp=5000)
            self.assertIsNotNone(res_s2)
            self.assertAlmostEqual(res_s2["window_pct"], 100.0)
            self.assertEqual(res_s2["window_masked_bp"], 5000)
            self.assertEqual(res_s2["tip_repeat_stretch_bp"], 763004)
            self.assertIn("contiguous 763,004 bp (763.0 kb) repeat tract directly abutting tip", res_s2["summary_str"])

            # 3' tip test on SUPER_3:
            res_s3 = analyze_terminal_repeats(tf_path, "SUPER_3", "3p", 10000000, window_bp=5000)
            self.assertIsNotNone(res_s3)
            self.assertAlmostEqual(res_s3["window_pct"], 100.0)
            self.assertEqual(res_s3["window_masked_bp"], 5000)
            self.assertEqual(res_s3["tip_repeat_stretch_bp"], 100000)
            self.assertIn("contiguous 100,000 bp (100.0 kb) repeat tract directly abutting tip", res_s3["summary_str"])

            # 5' tip test on SUPER_4 (repeat starts at 500 bp, not abutting tip):
            res_s4 = analyze_terminal_repeats(tf_path, "SUPER_4", "5p", 10000000, window_bp=5000)
            self.assertIsNotNone(res_s4)
            self.assertEqual(res_s4["tip_repeat_stretch_bp"], 0)
            self.assertEqual(res_s4["nearest_repeat_dist_bp"], 500)
            self.assertEqual(res_s4["first_repeat_len_bp"], 1000)
            self.assertIn("nearest repeat 500 bp from tip spanning 1,000 bp", res_s4["summary_str"])
        finally:
            Path(tf_path).unlink(missing_ok=True)

    def test_classify_unresolved_gap(self):
        # 1. Unbridged repeat array
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "Unbridged repeat array; contigs terminate at repeat boundary without spanning; gap preserved to protect unplaced intervening elements"}
        )
        self.assertEqual(cls, "Unbridged Repeat Array")
        self.assertIn("terminate at the repeat boundary", takeaway)

        # 2. Inversion breakpoint
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "Inversion breakpoint detected at gap boundary (ptg000001l discordant strands: left +, right -); gap preserved"}
        )
        self.assertEqual(cls, "Inversion Breakpoint Discontinuity")
        self.assertIn("discordant strand orientations", takeaway)

        # 3. Ambiguous multi-mapping anchor
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "No unitigs anchored with sufficient length and MAPQ >= 50"}
        )
        self.assertEqual(cls, "Ambiguous Multi-Mapping Anchor")
        self.assertIn("MAPQ", takeaway)

        # 4. De novo assembly failure
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "Hifiasm emitted zero contigs from 45 baited reads"}
        )
        self.assertEqual(cls, "De Novo Assembly Failure")
        self.assertIn("produced zero contigs", takeaway)

        # 5. Read coverage drop
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "Insufficient baited reads (1 reads, required >= 3)"},
            d_5p=0.5, d_3p=0.8
        )
        self.assertEqual(cls, "Read Coverage Drop to Zero")
        self.assertIn("evidence floor", takeaway)

        # 6. Consensus bridge rejected
        cls, note, takeaway = classify_unresolved_gap(
            {"notes": "Consensus bridge rejected: spanning reads do not agree on a patch"}
        )
        self.assertEqual(cls, "Consensus Bridge Rejected")
        self.assertIn("consensus bridging", takeaway)

        # 7. Fallback: spanning reads > 0, 0 unitigs
        cls, note, takeaway = classify_unresolved_gap(None, span_reads=3, n_unitigs=0)
        self.assertEqual(cls, "Unassembled Spanning Reads")
        self.assertIn("3 raw spanning read(s)", takeaway)

        # 8. Fallback: 0 spanning reads, 0 unitigs
        cls, note, takeaway = classify_unresolved_gap(None, span_reads=0, n_unitigs=0)
        self.assertEqual(cls, "Physical Sequencing Gap")
        self.assertIn("Zero spanning HiFi reads", takeaway)

    def test_render_telomere_unitig_table(self):
        import io

        # 1. Empty inventory -> no output
        buf = io.StringIO()
        render_telomere_unitig_table(buf, {})
        self.assertEqual(buf.getvalue(), "")

        # 2. Populated inventory -> Markdown table with hits, overhangs, and disposition
        inv = {
            "total_assembled": 2,
            "assigned_count": 1,
            "unassigned_count": 1,
            "termini_lacking_unitigs": ["Chr IV 5p (supported by cluster cl_17 with 5 reads)"],
            "unitigs": [
                {
                    "unitig_id": "utg000013l",
                    "length_bp": 34948,
                    "coverage_str": "18.5× | 18 reads (cl_12)",
                    "status": "ASSIGNED",
                    "assigned_to": "Chr III 5p",
                    "disposition": "Assigned to Chr III 5p (+22,664 bp extension).",
                    "hits": [
                        {
                            "target_chrom": "III",
                            "target_arm": "5p",
                            "target_start": 0,
                            "target_end": 12284,
                            "strand": "+",
                            "mapq": 60,
                            "overhang_bp": 22664,
                        }
                    ],
                },
                {
                    "unitig_id": "ptg000004l",
                    "length_bp": 33055,
                    "coverage_str": "N/A",
                    "status": "COMPETED_OUT",
                    "assigned_to": None,
                    "disposition": "Competed out at Chr II 5p (+26,293 bp) by `utg000007l` (+4,000 bp).",
                    "hits": [
                        {
                            "target_chrom": "II",
                            "target_arm": "5p",
                            "target_start": 0,
                            "target_end": 6762,
                            "strand": "+",
                            "mapq": 60,
                            "overhang_bp": 26293,
                        }
                    ],
                },
            ],
        }

        buf = io.StringIO()
        render_telomere_unitig_table(buf, inv)
        output = buf.getvalue()
        self.assertIn("### Telomeric Unitig Inventory & Hit Disposition", output)
        self.assertIn("`utg000013l`", output)
        self.assertIn("34,948 bp", output)
        self.assertIn("18.5× | 18 reads (cl_12)", output)
        self.assertIn("Chr III (5p): 0–12,284 (+, MAPQ 60)", output)
        self.assertIn("+22,664 bp", output)
        self.assertIn("**ASSIGNED**", output)
        self.assertIn("`ptg000004l`", output)
        self.assertIn("**COMPETED_OUT**", output)
        self.assertIn("Chromosome Termini Lacking Assembled Unitigs", output)
        self.assertIn("Chr IV 5p (supported by cluster cl_17 with 5 reads)", output)


if __name__ == "__main__":
    unittest.main()
