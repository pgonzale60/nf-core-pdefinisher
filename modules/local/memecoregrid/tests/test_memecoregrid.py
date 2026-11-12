#!/usr/bin/env python3
"""Unit tests for memecoregrid module."""

import pathlib
import sys
import unittest

# Ensure module lib is importable
CUR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = CUR
while REPO_ROOT.parent != REPO_ROOT and not (REPO_ROOT / "modules").is_dir():
    REPO_ROOT = REPO_ROOT.parent

LIB_DIR = (
    REPO_ROOT
    / "modules"
    / "local"
    / "memecoregrid"
    / "resources"
    / "usr"
    / "lib"
)
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from memegrid.evaluator import (  # noqa: E402
    calculate_consensus_entropy,
    detect_tandem_repeats,
    evaluate_meme_motif,
    is_tandem_repeat_motif,
    parse_meme_txt,
)
from memegrid.grid_runner import generate_grid_conditions  # noqa: E402
from memegrid.slicer import slice_window_at_radius  # noqa: E402


class TestMemeCoreGrid(unittest.TestCase):
    def test_window_slicing(self):
        # 300 bp test sequence: 150 'A's + 150 'C's
        seq = "A" * 150 + "C" * 150
        # Radius 30 should yield 30 'A's + 30 'C's = 60 bp
        sliced_30 = slice_window_at_radius(seq, radius=30, center_soma_len=150)
        self.assertEqual(len(sliced_30), 60)
        self.assertEqual(sliced_30, "A" * 30 + "C" * 30)

        # Radius 50 should yield 50 'A's + 50 'C's = 100 bp
        sliced_50 = slice_window_at_radius(seq, radius=50, center_soma_len=150)
        self.assertEqual(len(sliced_50), 100)
        self.assertEqual(sliced_50, "A" * 50 + "C" * 50)

    def test_consensus_entropy(self):
        # Homopolymer should have 0 entropy
        self.assertEqual(calculate_consensus_entropy("AAAAAAAAAA"), 0.0)
        # 4 bases balanced should have 2.0 bits
        self.assertAlmostEqual(calculate_consensus_entropy("ACGT"), 2.0, places=2)
        # Realistic motif TGTTCA
        ent = calculate_consensus_entropy("TGTTCA")
        self.assertGreater(ent, 1.0)

    def test_generate_grid_conditions(self):
        from memegrid.grid_runner import RADII

        self.assertIn(300, RADII)
        conds = generate_grid_conditions(radii=[30, 50])
        # 2 radii * 3 models * 2 strands * 3 width ranges = 36 conditions
        self.assertEqual(len(conds), 36)
        self.assertTrue(any(c["id"] == "R30_zoops_revcomp_w6_15" for c in conds))

    def test_parse_meme_txt(self):
        # Create a mock MEME output text
        mock_meme = """MEME version 4

MOTIF TGTTCA 1 width = 6 sites = 10 llr = 50 E-value = 1.2e-005

letter-probability matrix: alength= 4 w= 6 nsites= 10 E= 1.2e-005
  0.0 0.0 0.0 1.0
  0.0 0.0 1.0 0.0
  0.0 0.0 0.0 1.0
  0.0 0.0 0.0 1.0
  0.0 1.0 0.0 0.0
  1.0 0.0 0.0 0.0

Multilevel TGTTCA

sites sorted by position p-value
seq1  +  25  0.0001
seq2  +  28  0.0001
seq3  +  26  0.0001
seq4  +  27  0.0001
"""
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            tf.write(mock_meme)
            tf_path = tf.name

        metrics = parse_meme_txt(tf_path, radius=30, total_sequences=4)
        self.assertIsNotNone(metrics)
        self.assertEqual(metrics["width"], 6)
        self.assertEqual(metrics["consensus"], "TGTTCA")
        # Pure discrete matrix should give maximum IC = 2.0 bits per position
        self.assertAlmostEqual(metrics["mean_ic"], 2.0, places=1)
        self.assertEqual(metrics["occupancy"], 1.0)
        self.assertGreater(metrics["composite_score"], 0.5)

    def test_core_30bp_beats_short_6bp_motif(self):
        """Verify high-information monomer is not defeated by repetitive concatemer."""
        import tempfile

        # Condition 1: 6 bp, mean_ic ~ 1.85, evalue 1e-6
        # matrix of 6 rows
        m6 = """MEME version 4
MOTIF SHORT 1 width = 6 sites = 20 llr = 50 E-value = 1e-006
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 1e-006
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.0 0.9 0.05 0.05
  0.0 0.05 0.9 0.05
  0.05 0.05 0.0 0.9
Multilevel AAACGT
sites sorted by position p-value
s1 + 25 1e-5
s2 + 25 1e-5
"""
        # Condition 2: 30 bp, mean_ic ~ 1.55, evalue 1e-40 (repetitive AC concatemer)
        row_30 = "\n".join(["  0.8 0.1 0.1 0.0" if i % 2 == 0 else "  0.0 0.8 0.1 0.1" for i in range(30)])
        m30 = f"""MEME version 4
MOTIF LONGCORE 1 width = 30 sites = 20 llr = 200 E-value = 1e-040
letter-probability matrix: alength= 4 w= 30 nsites= 20 E= 1e-040
{row_30}
Multilevel {"AC" * 15}
sites sorted by position p-value
s1 + 25 1e-5
s2 + 25 1e-5
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f6:
            f6.write(m6)
            p6 = f6.name
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f30:
            f30.write(m30)
            p30 = f30.name

        met6 = parse_meme_txt(p6, radius=30, total_sequences=2)
        met30 = parse_meme_txt(p30, radius=30, total_sequences=2)

        self.assertIsNotNone(met6)
        self.assertIsNotNone(met30)
        # Condition 1 (high-IC 6 bp monomer) must beat Condition 2 (30 bp repetitive concatemer)
        self.assertGreater(met6["composite_score"], met30["composite_score"])

    def test_anr_site_multiplicity_credit(self):
        """Verify ANR site multiplicity credits site density and boosts composite score."""
        import tempfile

        def make_meme(nsites: int, n_seqs: int) -> str:
            rows = "\n".join(["  0.8 0.05 0.1 0.05"] * 8)
            sites_lines = "\n".join([f"s{i % n_seqs} + 25 1e-5" for i in range(nsites)])
            return f"""MEME version 4
MOTIF TEST 1 width = 8 sites = {nsites} llr = 100 E-value = 1e-10
letter-probability matrix: alength= 4 w= 8 nsites= {nsites} E= 1e-10
{rows}
Multilevel AAAAAAAA
sites sorted by position p-value
{sites_lines}
"""

        # 55 sites across 31 sequences (ANR multi-site recurrence)
        m_anr = make_meme(55, 31)
        # 31 single sites across 31 sequences (ZOOPS single-site)
        m_zoops = make_meme(31, 31)

        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f_anr, tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f_zoops:
            f_anr.write(m_anr)
            f_anr.flush()
            f_zoops.write(m_zoops)
            f_zoops.flush()

            res_anr = parse_meme_txt(f_anr.name, radius=30, total_sequences=31)
            res_zoops = parse_meme_txt(f_zoops.name, radius=30, total_sequences=31)

            self.assertIsNotNone(res_anr)
            self.assertIsNotNone(res_zoops)
            self.assertAlmostEqual(res_zoops["site_density_factor"], 1.0, places=2)
            self.assertGreater(res_anr["site_density_factor"], 1.5)
            self.assertGreater(res_anr["composite_score"], res_zoops["composite_score"])

    def test_phase2_joint_model_width_exploration(self):
        """Verify Phase 2 grid evaluates both zoops and anr across all width ranges."""
        import inspect
        from memegrid import grid_runner

        src = inspect.getsource(grid_runner.run_meme_grid)
        self.assertIn('for mod in ["zoops", "anr"]:', src)
        self.assertIn("P2_R{best_radius}_m{best_order}_{mod}_{st_lbl}_{w['label']}", src)

    def test_insignificant_evalue_filtered_out(self):
        import tempfile

        # Motif with evalue = 2.5
        m_insig = """MEME version 4
MOTIF NOISE 1 width = 6 sites = 5 llr = 10 E-value = 2.5
letter-probability matrix: alength= 4 w= 6 nsites= 5 E= 2.5
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.0 0.9 0.05 0.05
  0.0 0.05 0.9 0.05
  0.05 0.05 0.0 0.9
Multilevel AAACGT
sites sorted by position p-value
s1 + 25 0.1
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            tf.write(m_insig)
            p_tf = tf.name

        met = parse_meme_txt(p_tf, radius=30, total_sequences=1)
        self.assertIsNotNone(met)
        # Must have score 0 or be marked insignificant
        self.assertEqual(met["composite_score"], 0.0)

    def test_parallel_grid_runner_no_mpi_flag(self):
        import inspect
        from memegrid import grid_runner

        src = inspect.getsource(grid_runner)
        # Verify -p is NOT passed to meme command
        self.assertNotIn('"-p"', src)
        # Verify ProcessPoolExecutor is used
        self.assertIn("ProcessPoolExecutor", src)

    def test_dynamic_soma_len_normalization(self):
        import tempfile
        from memegrid.slicer import extract_soma_len_from_header, write_sliced_fasta

        # Header with short Illumina somatic flank (e.g. 93 bp)
        hdr = "anchor_5 class=PDE_BREAK_SITE soma_len=93"
        self.assertEqual(extract_soma_len_from_header(hdr), 93)

        # 243 bp sequence: 93 bp soma + 150 bp germline
        seq = "A" * 93 + "C" * 150
        with tempfile.NamedTemporaryFile("w", suffix=".fasta") as tf_fa:
            write_sliced_fasta([(hdr, seq)], radius=150, output_path=tf_fa.name)
            sliced_content = pathlib.Path(tf_fa.name).read_text()
            self.assertIn("anchor_5_r150_p0_93", sliced_content)

        # Mock MEME output where site is at start 92 (1 bp before cleavage)
        mock_meme = """MEME version 4
MOTIF TEST 1 width = 6 sites = 1 llr = 10 E-value = 1e-5
letter-probability matrix: alength= 4 w= 6 nsites= 1 E= 1e-5
  1.0 0.0 0.0 0.0
  0.0 1.0 0.0 0.0
  0.0 0.0 1.0 0.0
  0.0 0.0 0.0 1.0
  1.0 0.0 0.0 0.0
  0.0 1.0 0.0 0.0
Multilevel ACGTAC
sites sorted by position p-value
anchor_5_r150_p0_93  +  92  1e-5
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as tf_txt:
            tf_txt.write(mock_meme)
            tf_txt.flush()
            met = parse_meme_txt(tf_txt.name, radius=150, total_sequences=1)
            # Offset from motif center to pos0 must be (92 + 2.5) - 93 = 1.5 bp
            self.assertEqual(met["median_offset_to_pos0"], 1.5)

    def test_telomeric_echo_rejection_and_residual_masking(self):
        import tempfile
        from memegrid.evaluator import is_telomeric_repeat_echo
        from memegrid.slicer import mask_residual_telomeric_tracts

        # 1. Exact and 1-mismatch telomeric repeat motifs must be rejected
        is_echo, frac = is_telomeric_repeat_echo("TTAGGCTTAGGCTTAGGCTTA")
        self.assertTrue(is_echo)
        self.assertGreaterEqual(frac, 0.5)

        is_echo_deg, _ = is_telomeric_repeat_echo("TTAGGCTTAGACTAAGGCTTA")
        self.assertTrue(is_echo_deg)

        # Authentic Oscheius CBM must NOT be rejected
        is_echo_cbm, _ = is_telomeric_repeat_echo("TTCTGGTTTAGGCTCAGCCTATGCCATAA")
        self.assertFalse(is_echo_cbm)

        # Authentic Caenorhabditis core motifs must NOT be rejected
        is_echo_auri, _ = is_telomeric_repeat_echo("TTTTTCGCTTAGCC")
        self.assertFalse(is_echo_auri)

        is_echo_mono, _ = is_telomeric_repeat_echo("TTTCGCTTAGGCC")
        self.assertFalse(is_echo_mono)

        is_echo_parv, _ = is_telomeric_repeat_echo("TATCAGACTCGCATTTTGCTTGCT")
        self.assertFalse(is_echo_parv)

        # Spurious long telomeric arrays must be rejected
        is_echo_spurious, _ = is_telomeric_repeat_echo("AGTCTTAGGCTTAGGCTTAGGCTCATGCG")
        self.assertTrue(is_echo_spurious)

        # 2. parse_meme_txt must hard-reject telomeric echo with composite_score = 0.0
        row_12 = "\n".join(["  1.0 0.0 0.0 0.0"] * 12)
        telo_meme = f"""MEME version 4
MOTIF TELO 1 width = 12 sites = 10 llr = 100 E-value = 1e-30
letter-probability matrix: alength= 4 w= 12 nsites= 10 E= 1e-30
{row_12}
Multilevel TTAGGCTTAGGC
sites sorted by position p-value
s1 + 25 1e-8
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as tf_txt:
            tf_txt.write(telo_meme)
            tf_txt.flush()
            met = parse_meme_txt(tf_txt.name, radius=30, total_sequences=1)
            self.assertEqual(met["composite_score"], 0.0)
            self.assertEqual(met["rejection_flag"], "REJECTED_TELOMERIC_ECHO")

        # 3. Pre-MEME masking masks distal (TTAGGC)2+ tracts while preserving [-6..+6] around pos0
        win = "A" * 50 + "TTAGGCTTAGGCTTAGGC" + "C" * 76 + "TTAGGCTTAGGC" + "G" * 144
        self.assertEqual(len(win), 300)
        masked = mask_residual_telomeric_tracts(win, center_soma_len=150)
        # Distal tract at 50..68 must be masked with N
        self.assertEqual(masked[50:68], "N" * 18)
        # Cleavage junction [-6..+6] (144..156) must remain unmasked
        self.assertEqual(masked[144:156], win[144:156])

    def test_single_and_zero_sequence_graceful_handling(self):
        import json
        import tempfile
        from memegrid.grid_runner import export_grid_outputs, run_meme_grid

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = pathlib.Path(tmp_dir)
            single_fa = tmp_p / "single_seq.fasta"
            single_fa.write_text(">anchor_1 soma_len=150\n" + "A" * 150 + "C" * 150 + "\n")

            out_dir = tmp_p / "meme_grid_work"
            best_motif, grid_results = run_meme_grid(single_fa, out_dir)
            self.assertIsNone(best_motif)
            self.assertEqual(grid_results, [])

            out_meme = tmp_p / "best_core_motif.meme"
            out_tsv = tmp_p / "core_grid_summary.tsv"
            out_json = tmp_p / "core_motif_metrics.json"
            out_bg = tmp_p / "best_background.bfile"

            export_grid_outputs(
                best_motif=best_motif,
                grid_results=grid_results,
                output_meme=out_meme,
                summary_tsv=out_tsv,
                metrics_json=out_json,
                output_background=out_bg,
            )

            self.assertTrue(out_meme.exists())
            self.assertIn("MEME version 4", out_meme.read_text())

            self.assertTrue(out_tsv.exists())
            tsv_lines = out_tsv.read_text().strip().splitlines()
            self.assertEqual(len(tsv_lines), 1)  # Only header
            self.assertIn("composite_score", tsv_lines[0])

            self.assertTrue(out_json.exists())
            self.assertEqual(json.loads(out_json.read_text()), {})

            self.assertTrue(out_bg.exists())
            self.assertIn("# 0-order Markov background", out_bg.read_text())

            # Test 0 sequences
            zero_fa = tmp_p / "zero_seq.fasta"
            zero_fa.write_text("")
            best_0, res_0 = run_meme_grid(zero_fa, tmp_p / "meme_grid_work_0")
            self.assertIsNone(best_0)
            self.assertEqual(res_0, [])

    def test_subprocess_surfaces_stderr_on_failure(self):
        from unittest.mock import patch
        import subprocess
        from memegrid.grid_runner import _run_single_condition

        cond = {
            "id": "P1_test_fail",
            "radius": 30,
            "mod": "zoops",
            "revcomp": True,
            "w_args": ["-minw", "8", "-maxw", "30"],
            "w_label": "w8_30",
        }

        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                returncode=1,
                cmd=["meme"],
                stderr="Your (primary) dataset must contain at least 2 sequences",
            )
            with self.assertRaises(RuntimeError) as ctx:
                _run_single_condition(
                    meme_bin="/usr/bin/meme",
                    cond=cond,
                    r_fa_str="/tmp/test.fasta",
                    cond_dir_str="/tmp/cond_dir",
                    total_seqs=1,
                )
            self.assertIn("Your (primary) dataset must contain at least 2 sequences", str(ctx.exception))
            self.assertIn("P1_test_fail", str(ctx.exception))

    def test_is_tandem_repeat_motif(self):
        """Verify tandem satellite motifs are rejected while authentic CBMs pass."""
        import tempfile

        # 0. Direct helper tests for detect_tandem_repeats
        found, p, mon = detect_tandem_repeats("CTCAGGCTCAGG")
        self.assertTrue(found)
        self.assertEqual(p, 6)
        self.assertEqual(mon, "CTCAGG")

        # 1. Direct function tests for tandem arrays
        self.assertTrue(is_tandem_repeat_motif("CTCAGGCTCAGG"))
        self.assertTrue(is_tandem_repeat_motif("CTCAGGCTCAGGCTC"))
        self.assertTrue(is_tandem_repeat_motif("CAGCAGCAGCAG"))
        self.assertTrue(is_tandem_repeat_motif("AC" * 15))
        self.assertTrue(is_tandem_repeat_motif("ATG" * 4))

        # 2. Authentic CBMs must not be flagged as tandem satellite repeats
        self.assertFalse(is_tandem_repeat_motif("ACATTAGTGCT"))
        self.assertFalse(is_tandem_repeat_motif("TGTTCA"))
        self.assertFalse(is_tandem_repeat_motif("TTCTGGTTTAGGCTCAGCCTATGCCATAA"))
        self.assertFalse(is_tandem_repeat_motif("TTTTTCGCTTAGCC"))

        # 3. MEME parser hard-rejects tandem satellite motifs with composite_score = 0.0
        row_12 = "\n".join(["  0.8 0.1 0.1 0.0"] * 12)
        sat_meme = f"""MEME version 4
MOTIF SATELLITE 1 width = 12 sites = 20 llr = 150 E-value = 1e-50
letter-probability matrix: alength= 4 w= 12 nsites= 20 E= 1e-50
{row_12}
Multilevel CTCAGGCTCAGG
sites sorted by position p-value
s1 + 25 1e-8
s2 + 25 1e-8
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as tf_txt:
            tf_txt.write(sat_meme)
            tf_txt.flush()
            met = parse_meme_txt(tf_txt.name, radius=300, total_sequences=2)
            self.assertEqual(met["composite_score"], 0.0)
            self.assertEqual(met["rejection_flag"], "REJECTED_TANDEM_SATELLITE")

    def test_proximal_cleavage_pinned_cbm_beats_diffuse_satellite_at_r300(self):
        """Verify proximal cleavage-pinned CBM (-40 bp, IQR 15 bp) beats diffuse distal satellite at R=300."""
        import tempfile

        # Proximal cleavage-pinned CBM (e.g. Stevens benchmark -40 bp, tight IQR 15 bp)
        # Position 0 is at 300, so site at 260 -> offset -40 bp
        row_11 = "\n".join(["  0.8 0.1 0.05 0.05"] * 11)
        # 30 sites clustered around index 260 (IQR ~ 15 bp)
        prox_sites = "\n".join([f"s{i} + {255 + (i % 10)} 1e-5" for i in range(30)])
        prox_meme = f"""MEME version 4
MOTIF PROXIMAL_CBM 1 width = 11 sites = 30 llr = 200 E-value = 1e-15
letter-probability matrix: alength= 4 w= 11 nsites= 30 E= 1e-15
{row_11}
Multilevel ACATTAGTGCT
sites sorted by position p-value
{prox_sites}
"""

        # Distal diffuse tandem satellite motif at R=300 (offset -220 bp, diffuse spread)
        # Position 0 is at 300, so site at 80 -> offset -220 bp
        row_12 = "\n".join(["  0.85 0.05 0.05 0.05"] * 12)
        dist_sites = "\n".join([f"s{i % 30} + {50 + (i * 2)} 1e-20" for i in range(120)])
        dist_meme = f"""MEME version 4
MOTIF DISTAL_SAT 1 width = 12 sites = 120 llr = 600 E-value = 1e-100
letter-probability matrix: alength= 4 w= 12 nsites= 120 E= 1e-100
{row_12}
Multilevel CTCAGGCTCAGG
sites sorted by position p-value
{dist_sites}
"""

        with tempfile.NamedTemporaryFile("w", suffix=".txt") as f_prox, tempfile.NamedTemporaryFile("w", suffix=".txt") as f_dist:
            f_prox.write(prox_meme)
            f_prox.flush()
            f_dist.write(dist_meme)
            f_dist.flush()

            res_prox = parse_meme_txt(f_prox.name, radius=300, total_sequences=30)
            res_dist = parse_meme_txt(f_dist.name, radius=300, total_sequences=30)

            # 1. Distal satellite motif must be rejected
            self.assertEqual(res_dist["rejection_flag"], "REJECTED_TANDEM_SATELLITE")
            self.assertEqual(res_dist["composite_score"], 0.0)

            # 2. Proximal CBM passes and scores high
            self.assertEqual(res_prox["rejection_flag"], "PASS")
            self.assertGreater(res_prox["composite_score"], 0.5)
            self.assertGreater(res_prox["composite_score"], res_dist["composite_score"])

        # 3. Direct evaluator scoring verification without tandem rejection:
        # Verify that catalytic proximity factor and absolute IQR prioritize proximal over distal
        score_prox = evaluate_meme_motif(
            consensus="ACATTAGTGCT",
            width=11,
            mean_ic=1.8,
            occupancy=1.0,
            pos_iqr=15.0,
            median_offset=-40.0,
            evalue=1e-15,
            site_density=1.0,
        )
        score_dist_non_sat = evaluate_meme_motif(
            consensus="ACATTAGTGCT",
            width=11,
            mean_ic=1.8,
            occupancy=1.0,
            pos_iqr=60.0,
            median_offset=-220.0,
            evalue=1e-15,
            site_density=4.0,  # Hyper-dense ANR
        )
        self.assertGreater(score_prox["proximity_factor"], score_dist_non_sat["proximity_factor"])
        self.assertGreater(score_prox["iqr_penalty"], score_dist_non_sat["iqr_penalty"])
        self.assertGreater(score_prox["composite_score"], score_dist_non_sat["composite_score"])


if __name__ == "__main__":
    unittest.main()
