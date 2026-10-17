#!/usr/bin/env python3
"""Checks for the terminus-extension rules in pde_curation_synthesis.py.

Two rules that decide what gets appended to a chromosome end, both of which were wrong in the first
end-to-end run:

1. A Hifiasm unitig (`utg…`) is a maximal non-branching stretch of the assembly graph, while a contig
   (`ptg…`, `atg…`) is a chain of unitigs joined across junctions, so a contig's outboard tail can
   carry sequence belonging to a different genomic context. DF5120 III 5p was capped with
   `ptg000004l`, whose alignment covered only its first 1,893 bp while the other 33,055 bp were
   claimed as extension.
2. One assembled sequence is one molecule. The same sequence was authorised at two ends --
   `ptg000004l` at II 5p and III 5p, `atg000002l` at I 3p and X 3p -- which inserts it twice.

Usage: python3 scripts/check_unitig_assignment.py
"""
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
for p in (
    REPO / "modules/local/curationsynthesis/resources/usr/lib",
    REPO / "modules/local/telomerewalker/resources/usr/lib",
):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import curationsynthesis as syn
import telotract


def candidate(overhang, unitig_id, tract=0, subgap=0):
    # (overhang_bp, unitig_id, telomeric_tract_bp, subtelo_gap_bp, q_start, q_end, strand, unitig_len)
    return (overhang, unitig_id, tract, subgap, 0, overhang, "+", overhang)


def main():

    # --- A: a unitig is preferred over a longer contig for the same terminus.
    contig = candidate(33055, "ptg000004l")
    unitig = candidate(9000, "utg000013l", tract=600)
    ranking = syn.eligible_unitig_claims([contig, unitig], {"is_soma_boundary": False}, 0)
    assert [c[1] for c in ranking] == ["utg000013l", "ptg000004l"], ranking
    print("A ok: unitig preferred over a contig adding more sequence ({} before {})".format(
        ranking[0][1], ranking[1][1]))

    # --- B: a contig remains the fallback when no unitig reaches the end.
    ranking = syn.eligible_unitig_claims([contig], {"is_soma_boundary": False}, 0)
    assert [c[1] for c in ranking] == ["ptg000004l"], ranking
    print("B ok: contig still available when it is the only candidate")

    # --- C: a somatic stub still requires a unitig reaching past the break site.
    ranking = syn.eligible_unitig_claims([candidate(5000, "utg1", subgap=10),
                                          candidate(4000, "utg2", subgap=500)],
                                         {"is_soma_boundary": True}, germline_buf_bp=100)
    assert [c[1] for c in ranking] == ["utg2"], ranking
    print("C ok: somatic-stub buffer still applied ({} survives)".format(ranking[0][1]))

    # --- D: the double-authorisation seen in DF5120 is broken up, best extension first.
    claims = [
        ("II", "5p", candidate(26293, "ptg000004l")),
        ("III", "5p", candidate(33055, "ptg000004l")),
        ("X", "3p", candidate(17604, "atg000002l")),
        ("I", "3p", candidate(10878, "atg000002l")),
    ]
    assignment = syn.assign_unitigs_to_termini(claims)
    used = [c[1] for c in assignment.values()]
    assert len(used) == len(set(used)), "a sequence was assigned to more than one terminus: {}".format(used)
    assert sorted(used) == ["atg000002l", "ptg000004l"], used
    assert assignment[("III", "5p")][0] == 33055, "the larger extension should keep ptg000004l"
    assert assignment[("X", "3p")][0] == 17604, "the larger extension should keep atg000002l"
    assert ("II", "5p") not in assignment and ("I", "3p") not in assignment, assignment
    print("D ok: each sequence assigned once, to its largest extension ({})".format(
        ", ".join("{} {} <- {}".format(k[0], k[1], v[1]) for k, v in sorted(assignment.items()))))

    # --- E: a terminus that loses its first choice falls back to its next candidate.
    claims = [
        ("II", "5p", candidate(26293, "ptg000004l")),
        ("III", "5p", candidate(33055, "ptg000004l")),
        ("II", "5p", candidate(4000, "utg000007l")),
    ]
    assignment = syn.assign_unitigs_to_termini(claims)
    assert assignment[("III", "5p")][1] == "ptg000004l", assignment
    assert assignment[("II", "5p")][1] == "utg000007l", \
        "II 5p should fall back to its next candidate, got {}".format(assignment.get(("II", "5p")))
    print("E ok: a terminus losing its first choice falls back to {}".format(assignment[("II", "5p")][1]))

    # --- F: a terminus with no candidate left is absent, so the read/walker route resolves it.
    claims = [("II", "5p", candidate(1000, "utg1")), ("III", "5p", candidate(1000, "utg1"))]
    assignment = syn.assign_unitigs_to_termini(claims)
    assert len(assignment) == 1, assignment
    print("F ok: exactly one terminus keeps {} and the other is left to the read path".format(
        list(assignment.values())[0][1]))

    # --- G: molecule identity. Hifiasm names one molecule twice -- as a unitig and as the contig built
    #     from it -- and name-based arbitration still let DF5120 insert two slices of the same molecule:
    #     26,293 bp at II 5p and 33,055 bp at III 5p of one 34,948 bp sequence.
    def placed(overhang, ident, length, t_start, t_end):
        # (overhang, id, tract, subgap, q_start, q_end, strand, length, t_start, t_end)
        return (overhang, ident, 0, 0, 0, overhang, "-", length, t_start, t_end)

    lengths = {"ptg000004l": 34948, "utg000007l": 34948, "utg000001l": 34948}
    placements = {
        "ptg000004l": {("III", "-", 600, 1200)},   # the contig is placed here
        "utg000007l": {("III", "-", 600, 1200)},   # so is the unitig it is built from
        "utg000001l": {("V", "+", 40, 80)},        # same length, never shares a placement
    }
    molecule_ids = syn.build_molecule_ids(lengths, placements)
    assert molecule_ids["ptg000004l"] == molecule_ids["utg000007l"], molecule_ids
    assert molecule_ids["utg000001l"] != molecule_ids["ptg000004l"], \
        "equal-length molecules that never share a placement are distinct"
    claims = [("II", "5p", placed(26293, "ptg000004l", 34948, 0, 0)),
              ("III", "5p", placed(33055, "utg000007l", 34948, 0, 0))]
    same = syn.assign_unitigs_to_termini(
        claims, lambda chrom, candidate: molecule_ids.get(candidate[1]))
    assert len(same) == 1, "one molecule must be used at one terminus, got {}".format(same)
    assert list(same)[0] == ("III", "5p"), "the larger extension keeps the molecule, got {}".format(same)
    print("G ok: ptg000004l and utg000007l are one molecule; only III 5p keeps it (larger extension)")

    # --- H: extending the telomeric tract itself. The walker's terminal seed is oriented so that
    #     outward is its 3' end, so the tract is the seed's trailing run and the anchor is the
    #     non-telomeric sequence immediately before it. A single read is enough (N=1) to say how far
    #     the telomere runs; only the repeat bases beyond what is assembled get appended, never the
    #     read, which is how an "8,313 bp added / 654 bp telomere" extension came to be written.
    walker = telotract
    motif, rev = "TTAGGC", "GCCTAA"
    anchor = "".join("ACGT"[(i * 5 + 2) % 4] for i in range(500))
    seed = anchor + motif * 100                        # 600 bp of tract at the seed's 3' end
    assert walker.leading_telomere_bp(motif * 30, motif, rev) == 180
    assert walker.trailing_telomere_bp(anchor + motif * 30, motif, rev) == 180
    assert walker.trailing_telomere_bp(motif * 30 + anchor, motif, rev) == 0, \
        "a tract that is not terminal is not this end's tract"

    read_plus = "GGGG" + anchor + motif * 400            # same orientation, longer tract
    read_minus = walker.revcomp("GGGG" + anchor + motif * 200 + "TTTT")   # opposite orientation
    elsewhere = "GGGG" + motif * 400 + "TTTT"            # telomere, but not this tract
    best = walker.maximal_terminal_tract(
        [("short", read_minus), ("long", read_plus), ("elsewhere", elsewhere)],
        motif, rev, seed, anchor_bp=500, min_extra_bp=60)
    assert best is not None, "a read following the tract must be found"
    name, total, extra, orientation, insert_at = best
    assert (name, orientation) == ("long", "+"), (name, orientation)
    assert len(extra) == 1800 and total == 2400, (len(extra), total)
    assert insert_at == len(seed), insert_at
    assert walker.terminal_tract_region(seed, motif, rev, "3p", min_repeats=15)[0] == 600
    # A tract is allowed to be imperfect and to be followed by a short remnant: neither may read as 0.
    imperfect = anchor + motif * 100 + "TTTT"
    assert 500 <= walker.terminal_tract_region(imperfect, motif, rev, "3p", min_repeats=15)[0] <= 600
    assert walker.terminal_tract_region(imperfect, motif, rev, "3p", min_repeats=15)[2] == len(imperfect) - 4
    print("H ok: tract extended 600 -> {} bp from read {} ({} bp of repeat appended)".format(
        total, name, len(extra)))

    # A read that only looks telomeric elsewhere must not be mistaken for one continuing this tract.
    assert walker.maximal_terminal_tract([("elsewhere", elsewhere)], motif, rev, seed,
                                         anchor_bp=500, min_extra_bp=60) is None
    # Nor may an extension below the floor be reported as one.
    assert walker.maximal_terminal_tract([("tiny", "AAAA" + anchor + motif * 101)], motif, rev, seed,
                                         anchor_bp=500, min_extra_bp=90) is None
    print("I ok: a run elsewhere in the read, and an extension below the floor, are both refused")

    # --- J: terminus scarcity overrides overhang length.
    # When a terminus has only 1 candidate while another has alternative fallbacks, the scarce
    # terminus takes priority even if its overhang is shorter, allowing both ends to resolve.
    claims_scarcity = [
        ("II", "5p", candidate(26293, "utg000007l")),
        ("II", "5p", candidate(6798, "utg000006l")),
        ("III", "5p", candidate(22917, "utg000007l")),
    ]
    assignment_scarcity = syn.assign_unitigs_to_termini(claims_scarcity)
    assert assignment_scarcity[("III", "5p")][1] == "utg000007l", assignment_scarcity
    assert assignment_scarcity[("II", "5p")][1] == "utg000006l", assignment_scarcity
    print("J ok: scarce terminus III 5p (1 candidate) keeps utg000007l over II 5p (2 candidates, larger overhang)")

    # --- K: parse_telo_paf intra-chromosomal low-MAPQ exclusivity rule.
    import tempfile
    with tempfile.NamedTemporaryFile("w+", suffix=".paf") as tf:
        # utg_intra: 2 alignments, both on chr IV, MAPQ 0 (intra-chromosomal repeat)
        tf.write("utg_intra\t6000\t0\t5000\t-\tIV\t1000000\t0\t5000\t5000\t5000\t0\n")
        tf.write("utg_intra\t6000\t0\t5000\t-\tIV\t1000000\t2000\t7000\t5000\t5000\t0\n")
        # utg_inter: alignments on chr I and chr II with MAPQ 0 (inter-chromosomal ambiguity -> reject)
        tf.write("utg_inter\t6000\t0\t5000\t-\tI\t1000000\t0\t5000\t5000\t5000\t0\n")
        tf.write("utg_inter\t6000\t0\t5000\t-\tII\t1000000\t0\t5000\t5000\t5000\t0\n")
        # utg_high: MAPQ 60 on chr X -> accept
        tf.write("utg_high\t6000\t0\t5000\t-\tX\t1000000\t0\t5000\t5000\t5000\t60\n")
        tf.flush()

        exts = syn.parse_telo_paf(tf.name, ["I", "II", "IV", "X"])
        assert ("IV", "5p") in exts, "utg_intra must be accepted at IV 5p due to intra-chromosomal exclusivity"
        assert exts[("IV", "5p")][0][1] == "utg_intra"
        assert exts[("IV", "5p")][0][0] == 1000  # overhang 6000 - 5000
        assert ("I", "5p") not in exts, "utg_inter with MAPQ 0 across chromosomes must be rejected"
        assert ("II", "5p") not in exts, "utg_inter with MAPQ 0 across chromosomes must be rejected"
        assert ("X", "5p") in exts, "utg_high with MAPQ 60 must be accepted"
        assert exts[("X", "5p")][0][1] == "utg_high"
    print("K ok: parse_telo_paf admits low-MAPQ unitigs only when exclusive to a single chromosome")

    # --- L: synthesize_telomere_cards satellite repeat scaffolding from telo_junctions
    dummy_telo_junc = {
        "sample": "test_sample",
        "clusters": [
            {
                "cluster_id": "cl_sat_01",
                "read_count": 12,
                "meets_policy_floor": True,
                "classification": "GERMLINE_JUNCTION",
                "subtype": "SATELLITE_REPEAT",
                "target_chrom": "II",
                "target_arm": "5p",
                "target_pos": 0,
                "consensus_flank_60bp": "TTGATGAAGCTGTGAAGTTTCTGAAATTAGTGAATCCACTTATTATACTTAAAATTGAAT",
                "supporting_read_ids": [f"read_{i}" for i in range(12)],
                "mean_telomere_tract_bp": 950.0,
                "dominant_motif": "GCCTAA"
            }
        ]
    }
    cards = syn.synthesize_telomere_cards(
        ["II"], {"II": 15000000}, {}, {},
        {}, {},
        {("II", "5p"): "UNCAPPED", ("II", "3p"): "ALREADY_CAPPED"},
        {}, {},
        policy=syn.DEFAULT_POLICY,
        telo_junctions=dummy_telo_junc
    )
    card_5p = next(c for c in cards if c["chrom"] == "II" and c["terminus"] == "5p")
    assert card_5p["action"] == "SCAFFOLD_REPEAT_TO_TELOMERE", card_5p["action"]
    assert card_5p["classification"] == "SATELLITE_LINKED_TELOMERE"
    assert card_5p["gap_size"] == 100
    assert card_5p["evidence"]["read_count"] == 12
    assert card_5p["evidence"]["junction_cluster_id"] == "cl_sat_01"
    print("L ok: uncapped satellite arm II 5p emits SCAFFOLD_REPEAT_TO_TELOMERE action card with 100 bp N gap")

    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())