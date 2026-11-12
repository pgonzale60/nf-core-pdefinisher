#!/usr/bin/env python3
"""Unit tests for bipartiteprofiler module."""

import pathlib
import sys
import tempfile
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
    / "bipartiteprofiler"
    / "resources"
    / "usr"
    / "lib"
)
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from bipartite.aux_grid import (  # noqa: E402
    detect_tandem_repeats,
    extract_motif_n_from_meme,
    is_core_echo,
    run_localized_auxiliary_meme,
    score_auxiliary_candidate,
)
from bipartite.masker import (  # noqa: E402
    load_soma_len_map,
    mask_sequences_with_hits,
    slice_window_at_radius,
    write_sliced_fasta,
)
from bipartite.reporter import generate_cbm_discovery_report  # noqa: E402
from bipartite.spacing import (  # noqa: E402
    calculate_edge_to_edge_spacer,
    detect_homotypic_bipartite,
    profile_bipartite_spacing,
)


class TestBipartiteProfiler(unittest.TestCase):
    def test_edge_to_edge_spacer(self):
        # Case 1: Embedded
        # Core: 100..120. Aux: 105..115
        sp_emb, typ_emb = calculate_edge_to_edge_spacer(100, 120, 105, 115)
        self.assertEqual(sp_emb, 0)
        self.assertEqual(typ_emb, "Embedded")

        # Case 2: Upstream
        # Core: 100..120. Aux: 70..80 -> spacer = 100 - 80 = 20 bp
        sp_up, typ_up = calculate_edge_to_edge_spacer(100, 120, 70, 80)
        self.assertEqual(sp_up, 20)
        self.assertEqual(typ_up, "Flank_Upstream")

        # Case 3: Downstream
        # Core: 100..120. Aux: 135..145 -> spacer = 135 - 120 = 15 bp
        sp_dn, typ_dn = calculate_edge_to_edge_spacer(100, 120, 135, 145)
        self.assertEqual(sp_dn, 15)
        self.assertEqual(typ_dn, "Flank_Downstream")

    def test_mask_sequences_with_hits(self):
        seqs = {"seq1": "AAAAACCCCCGGGGGTTTTT"}
        hits = [
            {"seq_name": "seq1", "start": 6, "stop": 10}  # CCCCC -> NNNNN
        ]
        masked = mask_sequences_with_hits(seqs, hits)
        self.assertEqual(masked["seq1"], "AAAAANNNNNGGGGGTTTTT")

    def test_bipartite_spacing_profiling(self):
        seq_names = ["border_1", "border_2"]
        # Core at 100..120 (+), pos0 at 150 -> offset = 150 - 120 = 30 bp
        core_hits = {
            "border_1": [{"start": 100, "stop": 120, "strand": "+"}],
            "border_2": [{"start": 100, "stop": 120, "strand": "+"}],
        }
        # Aux at 70..80 on same strand (+) -> spacer = 20 bp
        aux_hits = {
            "border_1": [{"start": 70, "stop": 80, "strand": "+"}],
            "border_2": [{"start": 105, "stop": 115, "strand": "+"}],  # Embedded
        }

        rows, summary = profile_bipartite_spacing(
            seq_names, core_hits, aux_hits, pos0_coord=150
        )
        self.assertEqual(len(rows), 2)
        self.assertEqual(summary["n_paired_occurrences"], 2)
        self.assertEqual(summary["pct_embedded_0bp"], 50.0)
        self.assertEqual(summary["pct_spacer_le30bp"], 100.0)
        self.assertEqual(summary["pct_same_strand"], 100.0)
        self.assertEqual(summary["pct_retained_domain"], 100.0)
        self.assertEqual(summary["pct_eliminated_domain"], 0.0)
        self.assertEqual(rows[0]["partner_domain"], "Retained")

        # Test dynamic pos0 mapping (e.g. border_1 with pos0=93)
        pos0_map = {"border_1": 93, "border_2": 150}
        rows_dyn, _summary_dyn = profile_bipartite_spacing(
            seq_names, core_hits, aux_hits, pos0_coord=pos0_map
        )
        # In border_1 with pos0=93: Core at 100..120 (+) has offset 93 - 120 = -27 bp
        self.assertEqual(rows_dyn[0]["core_offset_to_pos0"], -27)
        self.assertEqual(rows[0]["core_offset_to_pos0"], 30)

        # Test load_soma_len_map from FASTA headers
        with tempfile.NamedTemporaryFile("w", suffix=".fasta") as tf_fa:
            tf_fa.write(">border_1 soma_len=93\nACGT\n>border_2\nACGT\n")
            tf_fa.flush()
            m_lens = load_soma_len_map(tf_fa.name, default=150)
            self.assertEqual(m_lens["border_1"], 93)
            self.assertEqual(m_lens["border_2"], 150)

    def test_report_generation(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json") as f_census, \
             tempfile.NamedTemporaryFile("w", suffix=".json") as f_core, \
             tempfile.NamedTemporaryFile("w", suffix=".md") as f_rep:

            f_census.write('{"n_breaks": 12, "median_c_soma": 25.0, "median_c_germ": 5.0, "germline_fraction": 0.167}')
            f_census.flush()

            f_core.write('{"consensus": "TGTTCA", "width": 6, "mean_ic": 1.95, "evalue_str": "1e-5", "occupancy": 0.92}')
            f_core.flush()

            bipartite_summary = {
                "n_paired_occurrences": 10,
                "median_spacer_bp": 15.0,
                "iqr_spacer_bp": 5.0,
                "pct_embedded_0bp": 20.0,
                "pct_spacer_le30bp": 90.0,
                "pct_spacer_le50bp": 100.0,
                "pct_same_strand": 90.0,
                "pct_retained_domain": 20.0,
                "pct_eliminated_domain": 80.0,
            }

            best_aux = {
                "consensus": "WWTWATGAACA",
                "width": 11,
                "mean_ic": 1.45,
                "evalue_str": "1e-4",
                "proximal_sites_6_50bp": 8,
            }

            generate_cbm_discovery_report(
                sample_id="test_species",
                census_json_path=f_census.name,
                core_metrics_json_path=f_core.name,
                best_aux_motif=best_aux,
                bipartite_summary=bipartite_summary,
                output_report_md=f_rep.name,
            )

            report_text = pathlib.Path(f_rep.name).read_text(encoding="utf-8")
            self.assertIn("CBM Discovery & PDE Junction Report: `test_species`", report_text)
            self.assertIn("TGTTCA", report_text)
            self.assertIn("WWTWATGAACA", report_text)
            self.assertIn("Median Edge-to-Edge Spacer", report_text)
            self.assertIn("Located in Retained Region", report_text)
            self.assertIn("Located in Eliminated Region", report_text)

            # Test Single-Motif Architecture (core found, no aux motif)
            generate_cbm_discovery_report(
                sample_id="test_species_single",
                census_json_path=f_census.name,
                core_metrics_json_path=f_core.name,
                best_aux_motif=None,
                bipartite_summary={},
                output_report_md=f_rep.name,
            )
            single_rep_text = pathlib.Path(f_rep.name).read_text(encoding="utf-8")
            self.assertIn("### Single-Motif Architecture", single_rep_text)
            self.assertIn("- No statistically significant auxiliary motif detected after core motif masking.", single_rep_text)
            self.assertNotIn("The cleavage system operates", single_rep_text)

            # Test Non-Motif Architecture (no core motif found, no aux motif)
            with tempfile.NamedTemporaryFile("w", suffix=".json") as f_empty_core:
                f_empty_core.write("{}")
                f_empty_core.flush()
                generate_cbm_discovery_report(
                    sample_id="test_species_nomotif",
                    census_json_path=f_census.name,
                    core_metrics_json_path=f_empty_core.name,
                    best_aux_motif=None,
                    bipartite_summary={},
                    output_report_md=f_rep.name,
                )
                nomotif_rep_text = pathlib.Path(f_rep.name).read_text(encoding="utf-8")
                self.assertIn("### Non-Motif Architecture", nomotif_rep_text)
                self.assertIn("- No statistically significant core or auxiliary recognition motif detected.", nomotif_rep_text)
                self.assertNotIn("The cleavage system operates", nomotif_rep_text)

    def test_score_auxiliary_candidate_favors_longer_palindrome_over_core_hexamer(self):
        # 11 bp imperfect palindrome (e.g. ATAGCGCTAT) across 15 sequences
        pal_11bp = {
            "width": 11,
            "nsites": 15,
            "evalue": 1.5e-4,
            "total_ic": 16.5,
            "mean_ic": 1.50,
            "consensus_entropy": 1.75,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 15,
        }
        score_11bp = score_auxiliary_candidate(pal_11bp, total_sequences=15)

        # 6 bp central core (AGCAGC) with higher per-column IC and multiple ANR hits
        core_6bp = {
            "width": 6,
            "nsites": 22,
            "evalue": 9.2e-6,
            "total_ic": 10.8,
            "mean_ic": 1.80,
            "consensus_entropy": 1.45,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 18,
        }
        score_6bp = score_auxiliary_candidate(core_6bp, total_sequences=15)

        self.assertGreater(score_11bp, score_6bp)
        self.assertGreater(score_11bp, 14.0)

        # Auxiliary scoring must NOT reward proximity to cleavage site
        pal_distal = dict(pal_11bp, proximal_sites_6_50bp=0)
        self.assertEqual(
            score_auxiliary_candidate(pal_distal, total_sequences=15),
            score_11bp,
        )

        # Telomeric repeat echo in auxiliary candidate must be hard-rejected (score 0.0)
        telo_aux = {
            "consensus": "TTAGGCTTAGGC",
            "width": 12,
            "nsites": 15,
            "evalue": 1e-12,
            "total_ic": 20.0,
            "mean_ic": 1.67,
            "consensus_entropy": 1.85,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 15,
        }
        self.assertEqual(score_auxiliary_candidate(telo_aux, total_sequences=15), 0.0)
        self.assertEqual(telo_aux["rejection_flag"], "REJECTED_TELOMERIC_ECHO")

    def test_aux_grid_surfaces_stderr_on_failure(self):
        from unittest.mock import patch
        import subprocess
        from bipartite.aux_grid import _run_single_aux_condition

        cond = {
            "id": "A03_test_fail",
            "mod": "zoops",
            "revcomp": True,
            "w_args": ["-w", "6"],
        }

        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                returncode=1,
                cmd=["meme"],
                stderr="Your (primary) dataset must contain at least 2 sequences",
            )
            with self.assertRaises(RuntimeError) as ctx:
                _run_single_aux_condition(
                    meme_bin="/usr/bin/meme",
                    cond=cond,
                    masked_fa_str="/tmp/test_masked.fasta",
                    cond_dir_str="/tmp/aux_cond_dir",
                )
            self.assertIn("Your (primary) dataset must contain at least 2 sequences", str(ctx.exception))
            self.assertIn("A03_test_fail", str(ctx.exception))

    def test_homotypic_bipartite_detection_mesorhabditis_belari_archetype(self):
        """Test detection and spacing profiling of dual-core homotypic bipartite architecture (M. belari archetype)."""
        # M. belari archetype: CCATTAATGC (10 bp) at -42 bp and +29 bp relative to pos0=100
        # Site 1 (somatic flank): start=49, stop=58 (offset to pos0 = 100 - 58 = 42 bp)
        # Site 2 (germline flank): start=129, stop=138 (offset to pos0 = 100 - 138 = -38 bp -> 29 bp into germline)
        # Edge-to-edge spacer = 129 - 58 = 71 bp (or 62 bp)
        seq_names = [f"nxMesBela1_break_{i}" for i in range(1, 11)]  # 10 loci
        core_hits = []
        core_hits_by_seq = {}
        for s in seq_names:
            h1 = {
                "seq_name": s,
                "motif_id": "CCATTAATGC",
                "start": 49,
                "stop": 58,
                "strand": "+",
                "score": 15.2,
                "p_value": 1e-6,
            }
            h2 = {
                "seq_name": s,
                "motif_id": "CCATTAATGC",
                "start": 129,
                "stop": 138,
                "strand": "+",
                "score": 14.8,
                "p_value": 2e-6,
            }
            core_hits.extend([h1, h2])
            core_hits_by_seq[s] = [h1, h2]

        is_homo, info = detect_homotypic_bipartite(core_hits, seq_names, min_fraction=0.30, max_spacer_iqr=25.0)
        self.assertTrue(is_homo)
        self.assertEqual(info["fraction_multi"], 1.0)
        self.assertEqual(info["n_multi_seqs"], 10)
        self.assertEqual(info["median_spacer"], 71.0)
        self.assertEqual(info["iqr_spacer"], 0.0)

        rows, summary = profile_bipartite_spacing(
            seq_names=seq_names,
            core_hits_by_seq=core_hits_by_seq,
            aux_hits_by_seq=None,
            pos0_coord=100,
            architecture="HOMOTYPIC_BIPARTITE",
        )
        # Uniform reporting: exactly 1 row per locus
        self.assertEqual(len(rows), 10)
        self.assertEqual(summary["n_paired_occurrences"], 10)
        self.assertEqual(summary["median_spacer_bp"], 71.0)
        self.assertEqual(summary["iqr_spacer_bp"], 0.0)
        self.assertEqual(summary["pct_same_strand"], 100.0)
        self.assertEqual(summary["pct_retained_domain"], 100.0)
        for r in rows:
            self.assertTrue(r["has_partner"])
            self.assertEqual(r["partner_type"], "Secondary_Core")
            self.assertEqual(r["spacer_bp"], 71)
            self.assertEqual(r["relative_strand"], "Same")

        # Negative control: only 1 sequence has a secondary hit (fraction 0.10 < 0.30)
        single_hits = [h for h in core_hits if h["start"] == 49]
        single_hits.append(core_hits_by_seq[seq_names[0]][1])  # add second hit for only 1 sequence
        is_homo_neg, info_neg = detect_homotypic_bipartite(single_hits, seq_names, min_fraction=0.30)
        self.assertFalse(is_homo_neg)
        self.assertEqual(info_neg["fraction_multi"], 0.1)

    def test_heterotypic_bipartite_resolution_df5120_archetype(self):
        """Test resolution of heterotypic bipartite architecture (Oscheius DF5120 archetype: 29 bp core + 11 bp palindrome with 11-13 bp spacer)."""
        seq_names = [f"nxOscSper1_break_{i}" for i in range(1, 11)]
        # Core at 50..78 (29 bp), pos0 at 100
        # Aux palindrome at 91..101 (11 bp) -> edge-to-edge spacer = 91 - 78 = 13 bp
        core_hits_by_seq = {
            s: [{"start": 50, "stop": 78, "strand": "+", "score": 28.5, "p_value": 1e-12}]
            for s in seq_names
        }
        aux_hits_by_seq = {
            s: [{"start": 91, "stop": 101, "strand": "+", "score": 14.2, "p_value": 1e-5}]
            for s in seq_names
        }

        rows, summary = profile_bipartite_spacing(
            seq_names=seq_names,
            core_hits_by_seq=core_hits_by_seq,
            aux_hits_by_seq=aux_hits_by_seq,
            pos0_coord=100,
            architecture="HETEROTYPIC_BIPARTITE",
        )
        self.assertEqual(len(rows), 10)
        self.assertEqual(summary["n_paired_occurrences"], 10)
        self.assertEqual(summary["median_spacer_bp"], 13.0)
        self.assertEqual(summary["pct_spacer_le30bp"], 100.0)
        self.assertEqual(summary["pct_same_strand"], 100.0)
        for r in rows:
            self.assertTrue(r["has_partner"])
            self.assertEqual(r["partner_type"], "Auxiliary")
            self.assertEqual(r["spacer_bp"], 13)

        # Verify candidate scoring and localized MEME extraction
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = pathlib.Path(tmp_d)
            meme_txt = tmp_p / "meme.txt"
            # Synthetic MEME text with Motif 1 (29 bp core) and Motif 2 (11 bp palindrome ATAGCWGCTAT)
            meme_txt.write_text(
                "MEME version 5\n\nALPHABET= ACGT\n\n"
                "MOTIF 1 MEME-1 width=29 nsites=10 E= 1e-25\n"
                "letter-probability matrix: alength= 4 w= 29 nsites= 10 E= 1e-25\n"
                + "0.85 0.05 0.05 0.05\n" * 29 +
                "Multilevel TTTTGGTTTTGGCAACACCAAAACCAATA\n\n"
                "MOTIF 2 MEME-2 width=11 nsites=10 E= 1e-5\n"
                "letter-probability matrix: alength= 4 w= 11 nsites= 10 E= 1e-5\n"
                + "0.85 0.05 0.05 0.05\n" * 11 +
                "Multilevel ATAGCGCTAT\n"
            )
            out_aux = tmp_p / "aux_motif2.meme"
            parsed = extract_motif_n_from_meme(meme_txt, 2, out_aux)
            self.assertIsNotNone(parsed)
            self.assertEqual(parsed["consensus"], "ATAGCGCTAT")
            self.assertEqual(parsed["width"], 11)

            # Echo and repeat checks
            self.assertFalse(is_core_echo("ATAGCGCTAT", "TTTTGGTTTTGGCAACACCAAAACCAATA"))
            self.assertFalse(detect_tandem_repeats("ATAGCGCTAT")[0])

            # Localized runner extracts and selects Motif 2
            best_cand, _cands = run_localized_auxiliary_meme(
                sliced_fasta=tmp_p / "fake.fa",
                out_dir=tmp_p,
                core_consensus="TTTTGGTTTTGGCAACACCAAAACCAATA",
                total_sequences=10,
            )
            self.assertIsNotNone(best_cand)
            self.assertEqual(best_cand["consensus"], "ATAGCGCTAT")
            self.assertEqual(best_cand["rejection_flag"], "PASS")
            self.assertGreater(best_cand["composite_score"], 0.0)

    def test_window_slicing_and_tandem_repeat_filter(self):
        """Test slicing 600 bp windows to optimal radius and filtering tandem repeats."""
        # Test 1: Slicing 600 bp window with pos0=300 to radius=50
        seq_600 = "A" * 300 + "C" * 300
        sliced_seq, new_pos0 = slice_window_at_radius(seq_600, radius=50, pos0=300)
        self.assertEqual(len(sliced_seq), 100)
        self.assertEqual(new_pos0, 50)
        self.assertEqual(sliced_seq[:50], "A" * 50)
        self.assertEqual(sliced_seq[50:], "C" * 50)

        # Test 2: Already sliced window (100 bp, pos0=50, radius=50) is unchanged
        s_unchanged, p_unchanged = slice_window_at_radius(sliced_seq, radius=50, pos0=50)
        self.assertEqual(s_unchanged, sliced_seq)
        self.assertEqual(p_unchanged, 50)

        # Test 3: write_sliced_fasta and load_soma_len_map
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_fa = pathlib.Path(tmp_d) / "sliced_out.fasta"
            raw_seqs = {"seqA": seq_600, "seqB": "G" * 600}
            soma_lens = {"seqA": 300, "seqB": 300}
            sliced_dict, new_lens = write_sliced_fasta(raw_seqs, soma_lens, radius=50, output_fasta=tmp_fa)
            self.assertEqual(len(sliced_dict), 2)
            self.assertEqual(new_lens["seqA"], 50)
            loaded_lens = load_soma_len_map(tmp_fa)
            self.assertEqual(loaded_lens["seqA"], 50)
            self.assertEqual(loaded_lens["seqB"], 50)

        # Test 4: Tandem satellite repeat rejection (M. belari flank hexamer CTCAGGCTCAGG)
        sat_hex = "CTCAGGCTCAGG"
        is_tr, period, frac = detect_tandem_repeats(sat_hex)
        self.assertTrue(is_tr)
        self.assertEqual(period, 6)
        self.assertEqual(frac, 1.0)

        cand_sat = {
            "consensus": sat_hex,
            "width": 12,
            "nsites": 10,
            "evalue": 1e-8,
            "total_ic": 18.0,
            "mean_ic": 1.5,
            "consensus_entropy": 1.4,
            "unique_sequences_hit": 10,
        }
        score_sat = score_auxiliary_candidate(cand_sat, total_sequences=10)
        self.assertEqual(score_sat, 0.0)
        self.assertEqual(cand_sat["rejection_flag"], "REJECTED_TANDEM_REPEAT")


if __name__ == "__main__":
    unittest.main()
