#!/usr/bin/env python3
"""
Unit test for policy-driven targeted gap assembly in repetitive and asymmetric flank regions.
Ensures SSOT curation_policy.json is strictly enforced and that evaluate_gap_paf handles:
1. Standard high-MAPQ dual-anchored bridges
2. Asymmetric anchor bridges (e.g. 20 kb left flank, 507 bp right flank, MAPQ 60)
3. Repetitive low-MAPQ flank bridges with high sequence identity
4. Proper rejection of insufficient anchors and low-identity candidates
"""

import os
import sys
import json
import tempfile
import importlib.util

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
POLICY_PATH = os.path.join(REPO_ROOT, "assets", "curation_policy.json")
TGA_SCRIPT = os.path.join(
    REPO_ROOT,
    "modules", "local", "targetedgapassembly",
    "resources", "usr", "bin", "pde_targeted_gap_assembly.py"
)


def load_module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_gap_assembly_policy():
    print("Testing gap assembly policy enforcement...")

    # 1. Verify SSOT policy parameters exist in assets/curation_policy.json
    assert os.path.exists(POLICY_PATH), f"Policy file not found: {POLICY_PATH}"
    with open(POLICY_PATH) as f:
        policy = json.load(f)

    gp = policy.get("gap_patch", {})
    ma = policy.get("mapping_and_anchors", {})

    min_mapq = ma.get("mapq_high_confidence", 50)
    min_aln_len = ma.get("anchor_min_length_bp", 1000)
    min_sec_anchor = gp.get("min_secondary_anchor_bp")
    min_identity = gp.get("min_flank_identity")
    rep_mapq_floor = gp.get("mapq_repetitive_gap_floor")
    require_unique_anchor = gp.get("require_unique_anchor_for_repetitive_gap")
    min_rep_anchor = gp.get("min_repetitive_anchor_bp")
    valency_enabled = gp.get("anchor_valency_survey_enabled")
    valency_kmer = gp.get("valency_kmer_size")
    valency_offtarget = gp.get("max_offtarget_repeat_fraction")
    allow_mono_overhang = gp.get("allow_monovalent_overhang_patch")
    min_mono_overhang = gp.get("min_monovalent_overhang_bp")
    mm_preset = gp.get("minimap2_preset")

    assert min_mapq == 50, f"Expected mapq_high_confidence == 50, got {min_mapq}"
    assert min_aln_len == 1000, f"Expected anchor_min_length_bp == 1000, got {min_aln_len}"
    assert min_sec_anchor == 300, f"Expected min_secondary_anchor_bp == 300, got {min_sec_anchor}"
    assert min_identity == 0.98, f"Expected min_flank_identity == 0.98, got {min_identity}"
    assert rep_mapq_floor == 0, f"Expected mapq_repetitive_gap_floor == 0, got {rep_mapq_floor}"
    assert require_unique_anchor is True, f"Expected require_unique_anchor_for_repetitive_gap is True, got {require_unique_anchor}"
    assert min_rep_anchor == 1000, f"Expected min_repetitive_anchor_bp == 1000, got {min_rep_anchor}"
    assert valency_enabled is True, f"Expected anchor_valency_survey_enabled is True, got {valency_enabled}"
    assert valency_kmer == 21, f"Expected valency_kmer_size == 21, got {valency_kmer}"
    assert valency_offtarget == 0.02, f"Expected max_offtarget_repeat_fraction == 0.02, got {valency_offtarget}"
    assert allow_mono_overhang is True, f"Expected allow_monovalent_overhang_patch is True, got {allow_mono_overhang}"
    assert min_mono_overhang == 200, f"Expected min_monovalent_overhang_bp == 200, got {min_mono_overhang}"
    assert mm_preset == "asm20", f"Expected minimap2_preset == 'asm20', got {mm_preset}"
    print("1. SSOT policy parameters in assets/curation_policy.json verified.")

    # 2. Test evaluate_gap_paf behavior
    tga = load_module_from_path("tga", TGA_SCRIPT)

    # Local gap coordinates: left flank 20000 bp, gap 200 bp, right flank 20000 bp (total 40200 bp)
    left_len = 20000
    gap_len = 200
    right_len = 20000
    gap_start = 20000
    gap_end = 20200

    def make_paf_line(qname, qlen, qstart, qend, strand, tstart, tend, matches, aln_len, mapq):
        return f"{qname}\t{qlen}\t{qstart}\t{qend}\t{strand}\tloc_ref\t40200\t{tstart}\t{tend}\t{matches}\t{aln_len}\t{mapq}\n"

    # Case A: Standard high-MAPQ bridge (both anchors >= 1000 bp, MAPQ 60)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_std", 5000, 0, 5000, "+", 18500, 21500, 3000, 3000, 60))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor
        )
        assert res["outcome"] == "PATCH_GAP_INSERTION", res
        assert "standard" in res["notes"]
    os.unlink(tf.name)
    print("2. Standard high-MAPQ bridge correctly accepted.")

    # Case B: Asymmetric anchor bridge (e.g. CARD-01 in PS1017: 20 kb left anchor, 507 bp right anchor, MAPQ 60, id 0.999)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("ptg000001l", 36000, 15400, 36000, "+", 0, 20707, 20500, 20510, 60))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor
        )
        assert res["outcome"] == "PATCH_GAP_INSERTION", res
        assert "asymmetric" in res["notes"]
        assert res["left_anchor"] == 20000
        assert res["right_anchor"] == 507
    os.unlink(tf.name)
    print("3. Asymmetric anchor bridge (PS1017 CARD-01 case) correctly accepted.")

    # Case C: Repetitive flank bridge with verified unique anchor (left anchor 1500 bp >= 1000 bp, right anchor 600 bp, MAPQ 0, identity 0.99)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_rep", 3000, 0, 3000, "+", 18500, 20800, 2277, 2300, 0))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor,
            require_unique_anchor=require_unique_anchor,
            min_repetitive_anchor=min_rep_anchor
        )
        assert res["outcome"] == "PATCH_GAP_INSERTION", res
        assert "repetitive flank" in res["notes"]
    os.unlink(tf.name)
    print("4. Repetitive low-MAPQ flank bridge with verified unique anchor correctly accepted.")

    # Case C2: Rejection of ambiguous repetitive candidate lacking unique anchor (left anchor 600 bp, right anchor 600 bp, MAPQ 0)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_ambig", 2000, 0, 2000, "+", 19400, 20800, 1386, 1400, 0))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor,
            require_unique_anchor=require_unique_anchor,
            min_repetitive_anchor=min_rep_anchor
        )
        assert res["outcome"] == "UNRESOLVED_GAP", res
        assert "Ambiguous repetitive bridge rejected" in res["notes"]
    os.unlink(tf.name)
    print("5. Ambiguous multi-locus repetitive bridge correctly rejected (gap preserved).")

    # Case C2b: Acceptance of monovalent candidate by Anchor Valency Conservation (left anchor 600 bp, right anchor 600 bp, MAPQ 0)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_ambig", 2000, 0, 2000, "+", 19400, 20800, 1386, 1400, 0))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor,
            require_unique_anchor=require_unique_anchor,
            min_repetitive_anchor=min_rep_anchor,
            monovalent_unitigs={"utg_ambig"}
        )
        assert res["outcome"] == "PATCH_GAP_INSERTION", res
        assert "monovalent anchor valency (P=1.0)" in res["notes"]
    os.unlink(tf.name)
    print("6. Monovalent repetitive bridge accepted by Anchor Valency Conservation (P=1.0).")

    # Case C2c: Monovalent repeat array overhang bridging unbridged repeat wall (e.g. CARD-03)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_overhang", 20000, 5000, 15000, "-", 10000, 20000, 10000, 10000, 60))
        tf.write(make_paf_line("utg_r", 20000, 0, 10000, "+", 20200, 30200, 10000, 10000, 60))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            monovalent_unitigs={"utg_overhang"},
            allow_monovalent_overhang_patch=allow_mono_overhang,
            min_monovalent_overhang=min_mono_overhang
        )
        assert res["outcome"] == "PATCH_GAP_INSERTION", res
        assert res["patch_len"] == 5000, res
        assert "Monovalent repeat array extension" in res["notes"]
    os.unlink(tf.name)
    print("6b. Monovalent repeat array overhang bridging unbridged repeat wall correctly accepted.")

    # Case C3: Rejection of inversion breakpoint candidate (discordant strands: left +, right -)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_inv", 4000, 0, 1500, "+", 18500, 20000, 1500, 1500, 60))
        tf.write(make_paf_line("utg_inv", 4000, 1500, 3000, "-", 20200, 21700, 1500, 1500, 60))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity
        )
        assert res["outcome"] == "UNRESOLVED_GAP", res
        assert "Inversion breakpoint detected" in res["notes"]
    os.unlink(tf.name)
    print("7. Inversion breakpoint candidate correctly detected and collinear bridge rejected.")

    # Case D: Rejection of low identity bridge in repetitive context (identity 0.90 < 0.98, MAPQ 0)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_bad_id", 2000, 0, 2000, "+", 19400, 20800, 1260, 1400, 0))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor
        )
        assert res["outcome"] == "UNRESOLVED_GAP", res
    os.unlink(tf.name)
    print("8. Divergent / low-identity candidate correctly rejected.")

    # Case E: Rejection of insufficient secondary anchor (< 500 bp)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_short_sec", 30000, 10000, 30000, "+", 0, 20300, 20300, 20300, 60))
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor
        )
        # Only 100 bp right anchor past 20200 -> must not be accepted as full bridge
        assert res["outcome"] != "PATCH_GAP_INSERTION", res
    os.unlink(tf.name)
    print("9. Insufficient secondary anchor (< 500 bp) correctly rejected from full bridge.")

    # Case F: False duplication overlap artefact excision (DF5120 Chr V case: 13,542 bp overlap + 200 bp N-gap)
    max_overlap = gp.get("max_overlap_excision_bp", 50000)
    assert max_overlap == 50000, f"Expected max_overlap_excision_bp == 50000, got {max_overlap}"
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        # CIGAR: 6459=13742D20000=
        cigar_paf = (
            "ptg000001l\t64808\t19985\t46444\t+\tV:1938669-1978869\t40201\t0\t40201\t"
            "26459\t40001\t60\tcg:Z:6459=13742D20000=\n"
        )
        tf.write(cigar_paf)
        tf.flush()
        res = tga.evaluate_gap_paf(
            tf.name, left_flank_len=20001, gap_len=200, right_flank_len=20000,
            g_start=1958669, g_end=1958869,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            repetitive_mapq_floor=rep_mapq_floor,
            max_overlap_excision_bp=max_overlap
        )
        assert res["outcome"] == "EXCISE_OVERLAP_ARTIFACT", res
        assert res["patch_len"] == 0, res
        assert res["excised_len"] == 13742, res
        assert res["excise_start"] == 1945127, res
        assert res["excise_end"] == 1958869, res
        assert "excise_range=1945127-1958869" in res["notes"], res
    os.unlink(tf.name)
    print("10. False duplication overlap artefact correctly identified and calculated.")

    # 11. Test AGP generator excision parsing
    AGP_SCRIPT = os.path.join(
        REPO_ROOT,
        "modules", "local", "curationagpgenerator",
        "resources", "usr", "bin", "pde_generate_agp.py"
    )
    pde_agp = load_module_from_path("pde_agp", AGP_SCRIPT)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write("card_id\tchrom\tstart\tend\tgap_len\toutcome\tleft_anchor\tright_anchor\tpatch_len\texcised_len\tunitig_id\tpatch_fasta\tnotes\n")
        tf.write("CARD-03\tV\t1958669\t1958869\t200\tEXCISE_OVERLAP_ARTIFACT\t6459\t20000\t0\t13742\tptg000001l\tNONE\tUnitig ptg000001l proves 13542 bp false duplication overlap artefact across gap; excise_range=1945127-1958869\n")
        tf.flush()
        manifest_rows = pde_agp.parse_tsv_manifest(tf.name)
        assert len(manifest_rows) == 1
        row = manifest_rows[0]
        assert row["outcome"] == "EXCISE_OVERLAP_ARTIFACT"
        import re
        m = re.search(r"excise_range=(\d+)-(\d+)", row["notes"])
        assert m is not None, "excise_range not found in manifest row notes"
        e_start = int(m.group(1))
        e_end = int(m.group(2))
        assert e_start == 1945127
        assert e_end == 1958869
        assert e_end - e_start == 13742
    os.unlink(tf.name)
    print("11. AGP generator manifest loading and coordinate excision parsing verified.")

    # 12. CIGAR-exact query coordinate projection across flank indels (+ and - strand)
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        # 50 bp deletion in left flank (10000=50D9950=), 200D + 1500I at gap, 20000= on right flank
        tf.write(
            "utg_cig\t41450\t0\t41450\t+\tloc_ref\t40200\t0\t40200\t"
            "39950\t40200\t60\tcg:Z:10000=50D9950=200D1500I20000=\n"
        )
        tf.flush()
        res_cig = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity
        )
        assert res_cig["outcome"] == "PATCH_GAP_INSERTION", res_cig
        assert (res_cig["q_patch_start"], res_cig["q_patch_end"], res_cig["patch_len"]) == (19950, 21450, 1500), res_cig
    os.unlink(tf.name)
    print("12. CIGAR-exact query coordinate projection across flank indels verified.")

    # 13. Priority ordering: minor partial flank extension must NOT mask a monovalent repeat overhang or inversion
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        # utg_minor extends 20 bp past gap_start_local (20020); utg_overhang has 5000 bp monovalent overhang
        tf.write(make_paf_line("utg_minor", 15000, 0, 15000, "+", 5020, 20020, 15000, 15000, 60))
        tf.write(make_paf_line("utg_overhang", 20000, 5000, 15000, "-", 10000, 20000, 10000, 10000, 60))
        tf.write(make_paf_line("utg_r", 20000, 0, 10000, "+", 20200, 30200, 10000, 10000, 60))
        tf.flush()
        res_prio = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            monovalent_unitigs={"utg_overhang"},
            allow_monovalent_overhang_patch=True,
            min_monovalent_overhang=200
        )
        assert res_prio["outcome"] == "PATCH_GAP_INSERTION" and res_prio["unitig_id"] == "utg_overhang", res_prio
    os.unlink(tf.name)
    print("13. Priority ordering verified: minor flank extension does not mask monovalent repeat overhang.")

    # 14. Dual monovalent repeat overhangs overlapping across a repeat gap
    with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
        tf.write(make_paf_line("utg_L", 2500, 0, 2000, "+", 18000, 20000, 2000, 2000, 60))
        tf.write(make_paf_line("utg_R", 2500, 500, 2500, "+", 20200, 22200, 2000, 2000, 60))
        tf.flush()
        shared_ovl = "ACGTACGTACGTACGTACGT" * 10  # 200 bp overlap
        utg_seqs = {
            "utg_L": ("A" * 2000) + ("C" * 300) + shared_ovl,
            "utg_R": shared_ovl + ("G" * 300) + ("T" * 2000),
        }
        res_dual = tga.evaluate_gap_paf(
            tf.name, left_len, gap_len, right_len,
            min_mapq=min_mapq, min_aln_len=min_aln_len,
            min_secondary_anchor=min_sec_anchor, min_identity=min_identity,
            monovalent_unitigs={"utg_L", "utg_R"},
            allow_monovalent_overhang_patch=True,
            min_monovalent_overhang=200,
            unitig_seqs=utg_seqs
        )
        assert res_dual["outcome"] == "PATCH_GAP_INSERTION", res_dual
        assert res_dual["unitig_id"] == "utg_L+utg_R" and res_dual["patch_len"] == 800, res_dual
    os.unlink(tf.name)
    print("14. Dual monovalent repeat overhang overlap bridge verified.")

    print("\nALL GAP ASSEMBLY POLICY TESTS PASSED!")
    return 0


if __name__ == "__main__":
    sys.exit(test_gap_assembly_policy())
