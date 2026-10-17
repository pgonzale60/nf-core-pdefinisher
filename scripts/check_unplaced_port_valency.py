#!/usr/bin/env python3
"""Regression and invariant suite for `unplaced_triage.py`, `port_valency.py`,
and Conservation of Sequence in `pde_generate_agp.py`.

Usage: python3 scripts/check_unplaced_port_valency.py
"""

import json
import pathlib
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[1]
SYN_LIB = REPO / "modules/local/curationsynthesis/resources/usr/lib"
if str(SYN_LIB) not in sys.path:
    sys.path.insert(0, str(SYN_LIB))

import curationsynthesis as syn

AGP_GEN = REPO / "modules/local/curationagpgenerator/resources/usr/bin/pde_generate_agp.py"


def main():

    # --- 1. Allelic Haplotig vs Diverged Paralogous Gene-Island Triage ---
    unplaced_chroms = ["scaffold_gene_island", "scaffold_true_haplotig", "scaffold_sat"]
    chrom_lens = {
        "SUPER_X": 10000000,
        "scaffold_gene_island": 45000,
        "scaffold_true_haplotig": 45000,
        "scaffold_sat": 20000,
    }
    scaff_repeats = {
        "scaffold_gene_island": 5000,    # 11.1% repeat density (gene-bearing contig)
        "scaffold_true_haplotig": 4000,  # 8.9% repeat density
        "scaffold_sat": 17000,           # 85.0% repeat density (satellite array)
    }
    scaff_any_records = {
        # Diverged paralogous gene hits at 89% identity across 44.5 kb: must NOT count as allelic haplotig!
        "scaffold_gene_island": [
            {
                "qstart": 0, "qend": 44500, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 1000000, "tend": 1044500, "mapq": 0, "pid": 89.0, "alen": 44500,
                "matches": 39605, "strand": "+", "tp": "S",
            }
        ],
        # True collinear allelic haplotig at 99.4% identity across 44.6 kb: must be purged as Haplotig
        "scaffold_true_haplotig": [
            {
                "qstart": 0, "qend": 44600, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 2000000, "tend": 2044600, "mapq": 60, "pid": 99.4, "alen": 44600,
                "matches": 44332, "strand": "+", "tp": "P",
            }
        ],
        # Redundant satellite array at 90% identity with 85% repeat density: must be purged as Satellite
        "scaffold_sat": [
            {
                "qstart": 0, "qend": 19500, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 0, "tend": 19500, "mapq": 0, "pid": 90.0, "alen": 19500,
                "matches": 17550, "strand": "+", "tp": "S",
            }
        ],
    }

    triage = syn.triage_unplaced_scaffolds(
        unplaced_chroms=unplaced_chroms,
        chrom_lens=chrom_lens,
        scaff_repeats=scaff_repeats,
        scaff_any_records=scaff_any_records,
        scaff_rdna_operons={},
        scaff_rdna_operon_bp={},
    )
    assert triage["disposition"]["scaffold_gene_island"][0] == "Genuine_Contig", triage["disposition"]
    assert triage["scaff_novel_bp"]["scaffold_gene_island"] == 45000, triage["scaff_novel_bp"]
    assert triage["disposition"]["scaffold_true_haplotig"] == ("Haplotig", "PURGE_HAPLOTIG")
    assert triage["disposition"]["scaffold_sat"] == ("Satellite", "PURGE_HAPLOTIG")
    print("1 ok: diverged paralogous gene island retained as Genuine_Contig; true haplotig & satellite purged")

    # --- 2. Directional Outward Overhang Geometry + Case 2 Internal Gap Island Bridge + Low-MAPQ Exclusive Terminus ---
    # Test sc_inward_trap: aligns at 5p on '+' strand with qstart=200, qend=5000 (slen=20000).
    # Outward overhang at 5p is qstart=200 (<1000), NOT slen-qend=15000 (which points inward!).
    # Test ptg000262l: 50 kb gene-bearing island whose 5' (0..2000) and 3' (48000..50000) repeat ports
    # bridge internal gap CARD-INT-01 on SUPER_X (500000..500200) at MAPQ 0, 99.2% identity.
    # Test scaffold_x_tip: 25 kb contig whose 5' end (0..4000) aligns at MAPQ 0, 99.1% identity
    # exclusively to open terminus SUPER_X 3p (9996000..10000000) on '+' strand (outward overhang = 21000 bp).
    chrom_lens_2 = {
        "SUPER_X": 10000000,
        "SUPER_1": 12000000,
        "sc_inward_trap": 20000,
        "ptg000262l": 50000,
        "scaffold_x_tip": 25000,
    }
    any_recs_2 = {
        "sc_inward_trap": [
            {
                "qstart": 200, "qend": 5000, "tname": "SUPER_1", "tlen": 12000000,
                "tstart": 100, "tend": 4900, "mapq": 60, "pid": 99.8, "alen": 4800,
                "matches": 4790, "strand": "+", "tp": "P",
            }
        ],
        "ptg000262l": [
            {
                "qstart": 0, "qend": 2000, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 498000, "tend": 500000, "mapq": 0, "pid": 99.2, "alen": 2000,
                "matches": 1984, "strand": "+", "tp": "S",
            },
            {
                "qstart": 48000, "qend": 50000, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 500200, "tend": 502200, "mapq": 0, "pid": 99.4, "alen": 2000,
                "matches": 1988, "strand": "+", "tp": "S",
            },
        ],
        "scaffold_x_tip": [
            {
                "qstart": 0, "qend": 4000, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 9996000, "tend": 10000000, "mapq": 0, "pid": 99.1, "alen": 4000,
                "matches": 3964, "strand": "+", "tp": "S",
            },
            {
                "qstart": 0, "qend": 2000, "tname": "SUPER_X", "tlen": 10000000,
                "tstart": 9992000, "tend": 9994000, "mapq": 0, "pid": 98.5, "alen": 2000,
                "matches": 1970, "strand": "+", "tp": "S",
            },
        ],
    }
    prim_recs_2 = {
        "sc_inward_trap": any_recs_2["sc_inward_trap"],
        "ptg000262l": [],
        "scaffold_x_tip": [],
    }
    internal_cards = [
        {
            "card_id": "CARD-INT-01",
            "domain": "DOMAIN_INTERNAL",
            "chrom": "SUPER_X",
            "start": 500000,
            "end": 500200,
            "action": "REQUIRES_UNITIG_PATCH",
        }
    ]
    triage_2 = syn.triage_unplaced_scaffolds(
        unplaced_chroms=["sc_inward_trap", "ptg000262l", "scaffold_x_tip"],
        chrom_lens=chrom_lens_2,
        scaff_repeats={},
        scaff_any_records=any_recs_2,
        scaff_rdna_operons={},
        scaff_rdna_operon_bp={},
    )
    scaff_cards, scaff_disp, ext_by_seq = syn.evaluate_extension_candidacy(
        triage_result=triage_2,
        scaff_primary_records=prim_recs_2,
        scaff_rdna_operon_bp={},
        scaffold_attributions={},
        telo_cards=[],
        rdna_cards=[],
        nuclear_chroms=["SUPER_X", "SUPER_1"],
        chrom_lens=chrom_lens_2,
        scaff_any_records=any_recs_2,
        internal_cards=internal_cards,
        open_termini=[("SUPER_X", "3p"), ("SUPER_1", "5p")],
        policy=syn.DEFAULT_POLICY,
    )

    assert "sc_inward_trap" not in ext_by_seq, (
        "inward overhang must not be mistaken for an outward terminal extension"
    )
    assert ext_by_seq["ptg000262l"]["action"] == "INTEGRATE_UNPLACED_REPEAT", ext_by_seq["ptg000262l"]
    assert ext_by_seq["ptg000262l"]["target_gap_id"] == "CARD-INT-01"
    assert (ext_by_seq["ptg000262l"]["insert_qstart"], ext_by_seq["ptg000262l"]["insert_qend"]) == (2000, 48000)

    assert ext_by_seq["scaffold_x_tip"]["action"] == "INTEGRATE_INTO_CHROMOSOME", ext_by_seq["scaffold_x_tip"]
    assert ext_by_seq["scaffold_x_tip"]["valency_mode"] == "CHROMOSOME_EXCLUSIVE_PORT"
    assert (ext_by_seq["scaffold_x_tip"]["ext_qstart"], ext_by_seq["scaffold_x_tip"]["ext_qend"]) == (4000, 25000)
    print("2 ok: directional overhang rejected inward trap; ptg000262l bridged CARD-INT-01; scaffold_x_tip extended SUPER_X 3p")

    # --- 3. SUPER_X_unloc_1 Prefix Extraction & Core Contig Protection ---
    _, attributions = syn.match_open_termini_to_scaffolds(
        nuclear_chroms=["SUPER_X"],
        chrom_lens={"SUPER_X": 10000000, "SUPER_X_unloc_1": 60000, "SUPER_X_unloc_sat": 15000},
        unplaced_chroms=["SUPER_X_unloc_1", "SUPER_X_unloc_sat"],
        scaff_records={
            "SUPER_X_unloc_1": [
                {
                    "tname": "SUPER_X", "tlen": 10000000, "tstart": 3000000, "tend": 3005000,
                    "qstart": 0, "qend": 5000, "mapq": 60, "strand": "+",
                }
            ],
            "SUPER_X_unloc_sat": [],
        },
        telo_junctions={"clusters": []},
        open_termini=[("SUPER_X", "5p")],
        policy=syn.DEFAULT_POLICY,
    )
    assert "SUPER_X_unloc_1" not in attributions, (
        "unloc contig with internal chromosomal alignment must not be absorbed into terminal cassette"
    )
    assert "SUPER_X_unloc_sat" in attributions and attributions["SUPER_X_unloc_sat"]["chromosome"] == "SUPER_X"
    print("3 ok: SUPER_X_unloc_1 protected from terminal cassette absorption; SUPER_X prefix parsed accurately")

    # --- 4. End-to-End AGP Materialisation & Conservation of Sequence Guard ---
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)
        ref_fa = tmp / "ref.fa"
        # SUPER_X: 1,000 bp before gap (0..1000), 200 bp N-gap (1000..1200), 1,000 bp after gap (1200..2200)
        super_x_seq = ("A" * 1000) + ("N" * 200) + ("C" * 1000)
        # ptg000262l: 100 bp left repeat + 500 bp unique gene island ('G'*500) + 100 bp right repeat
        ptg_seq = ("A" * 100) + ("G" * 500) + ("C" * 100)
        # scaffold_x_tip: 100 bp overlap ('C'*100) + 300 bp extension ('T'*300)
        tip_seq = ("C" * 100) + ("T" * 300)
        # scaffold_unmaterialized: marked RETIRED_ABSORBED for 5p, but walker did NOT extend 5p -> must be RETAINED!
        unmat_seq = "ACGT" * 200

        ref_fa.write_text(
            f">SUPER_X\n{super_x_seq}\n"
            f">scaffold_ptg000262l\n{ptg_seq}\n"
            f">scaffold_x_tip\n{tip_seq}\n"
            f">scaffold_unmaterialized\n{unmat_seq}\n"
        )

        blueprint = {
            "action_cards": [
                {
                    "card_id": "CARD-INT-01",
                    "domain": "DOMAIN_INTERNAL",
                    "chrom": "SUPER_X",
                    "start": 1000,
                    "end": 1200,
                    "action": "REQUIRES_UNITIG_PATCH",
                },
                {
                    "card_id": "CARD-SCAFF-INT-01",
                    "domain": "DOMAIN_SCAFFOLD",
                    "scaffold": "scaffold_ptg000262l",
                    "chrom": "SUPER_X",
                    "start": 1000,
                    "end": 1200,
                    "target_gap_id": "CARD-INT-01",
                    "insert_qstart": 100,
                    "insert_qend": 600,
                    "left_gap_bp": 0,
                    "right_gap_bp": 0,
                    "orientation": "+",
                    "action": "INTEGRATE_UNPLACED_REPEAT",
                },
                {
                    "card_id": "CARD-SCAFF-INT-02",
                    "domain": "DOMAIN_SCAFFOLD",
                    "scaffold": "scaffold_x_tip",
                    "chrom": "SUPER_X",
                    "terminus": "3p",
                    "ext_qstart": 100,
                    "ext_qend": 400,
                    "orientation": "+",
                    "action": "INTEGRATE_INTO_CHROMOSOME",
                },
            ],
            "scaffold_dispositions": [
                {
                    "scaffold": "scaffold_ptg000262l",
                    "action": "INTEGRATE_UNPLACED_REPEAT",
                },
                {
                    "scaffold": "scaffold_x_tip",
                    "action": "INTEGRATE_INTO_CHROMOSOME",
                },
                {
                    "scaffold": "scaffold_unmaterialized",
                    "chromosome": "SUPER_X",
                    "terminus": "5p",
                    "action": "RETIRED_ABSORBED",
                },
            ],
        }
        bp_json = tmp / "blueprint.json"
        bp_json.write_text(json.dumps(blueprint, indent=2))

        subprocess.run(
            [
                sys.executable,
                str(AGP_GEN),
                "--blueprint", str(bp_json),
                "--ref", str(ref_fa),
                "--sample", "test_valency",
                "--outdir", str(tmp),
                "--synthesize-fasta",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        curated_fa = (tmp / "test_valency.curated.fasta").read_text()
        # Parse FASTA
        seqs = {}
        curr = None
        for line in curated_fa.splitlines():
            if line.startswith(">"):
                curr = line[1:].strip()
                seqs[curr] = ""
            elif curr:
                seqs[curr] += line.strip()

        expected_super_x = ("A" * 1000) + ("G" * 500) + ("C" * 1000) + ("T" * 300)
        assert seqs["SUPER_X"] == expected_super_x, (
            f"SUPER_X mismatch: len={len(seqs['SUPER_X'])}, expected={len(expected_super_x)}"
        )
        assert "scaffold_ptg000262l" not in seqs, "materialised gap island should be retired from unplaced"
        assert "scaffold_x_tip" not in seqs, "materialised terminal extension should be retired from unplaced"
        assert "scaffold_unmaterialized" in seqs and seqs["scaffold_unmaterialized"] == unmat_seq, (
            "unmaterialised RETIRED_ABSORBED scaffold MUST be retained by Conservation of Sequence guard!"
        )
    print("4 ok: pde_generate_agp.py splices gap island + terminal extension and retains unmaterialised scaffold")
    print("all unplaced port valency checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
