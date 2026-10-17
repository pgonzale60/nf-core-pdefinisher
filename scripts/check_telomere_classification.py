#!/usr/bin/env python3
"""Regression check for the germline-tip vs somatic-stub classification.

Two defects motivated this, both found on the real nxOscSper1_DF5120 output:

  1. A terminus whose distal sequence is GRS was labelled SOMA_TRUNCATED. That asserts the somatic
     core ends there, but the core actually ends d_soma away, *inside* the assembly: the tip is the
     germline end of the arm. For DF5120 III, the tip is at coordinate 0 and d_soma is 61,996.
  2. The same flag pre-filtered germline unitigs (requiring them to "bridge" the breaksite), so the
     unitig sitting at the tip -- ptg000004l at III 253-11,756, MAPQ 60 -- was discarded. The end
     then fell through to a read cap, which is why DF5120 had zero RESOLVED_BY_UNITIG rows.

Run: python3 scripts/check_telomere_classification.py
"""
import importlib.util
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
SYN_LIB = REPO / "modules/local/curationsynthesis/resources/usr/lib"
WALKER_LIB = REPO / "modules/local/telomerewalker/resources/usr/lib"
WALKER = REPO / "modules/local/telomerewalker/resources/usr/bin/pde_telomere_walker.py"
for p in (SYN_LIB, WALKER_LIB):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import curationsynthesis as syn
import telocapping


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for attr in ("outboard_extension", "region_extension"):
        setattr(module, attr, getattr(telocapping, attr))
    return module


def terminus_card(cards, chrom, terminus):
    for card in cards:
        if card.get("chrom") == chrom and card.get("terminus") == terminus:
            return card
    raise AssertionError("no card for {} {}".format(chrom, terminus))


def main():
    walk = load(WALKER, "pde_telomere_walker")
    policy = syn.DEFAULT_POLICY

    # --- Case A: DF5120 III 5p, the real geometry: GRS at the tip, core inward, S2G at 61,996,
    #     and a MAPQ-60 telomeric unitig sitting at the tip.
    cards = syn.synthesize_telomere_cards(
        nuclear_chroms=["III"],
        chrom_lens={"III": 7104590},
        grs_segments={"III": [(0, 61996, 61996, 16.19), (7052312, 7104590, 52278, 9.99)]},
        core_segments={"III": [(61996, 7052312, 6990316, 65.65)]},
        delfies_points={"III": [61996, 7052312]},
        rdna_intervals={},
        baseline_status={("III", "5p"): "UNCAPPED", ("III", "3p"): "UNCAPPED"},
        unitig_extensions={("III", "5p"): [(10976, "ptg000004l", 5000, 0)]},
        single_read_ext={("III", "5p"): (0, 0, "None", 0)},
        policy=policy,
    )
    card = terminus_card(cards, "III", "5p")
    assert card["classification"] == "RESOLVABLE_BY_UNITIG", \
        "germline tip should accept its unitig, got {}".format(card["classification"])
    assert card.get("mechanism") == "CAP_WITH_UNITIG"
    assert card.get("terminus_is_germline") is True
    assert "somatic_diminution_resolved" in card["evidence"], \
        "the unitig bridges past an existing diminution breakpoint; that evidence must be recorded"
    print("A ok: III 5p germline tip -> {} by {} (bridges Dsoma)".format(
        card["classification"], card["evidence"].get("unitig_id")))

    # --- Case B: control. A genuine somatic stub: no distal GRS, core stops short, breakpoint
    #     present. Must still resolve to the somatic branch, unchanged.
    cards = syn.synthesize_telomere_cards(
        nuclear_chroms=["IV"],
        chrom_lens={"IV": 900000},
        grs_segments={"IV": [(1000, 5000, 4000, 10.0)]},
        core_segments={"IV": [(0, 500000, 500000, 60.0)]},
        delfies_points={"IV": [500000]},
        rdna_intervals={},
        baseline_status={("IV", "5p"): "UNCAPPED", ("IV", "3p"): "UNCAPPED"},
        unitig_extensions={},
        single_read_ext={},
        policy=policy,
    )
    card = terminus_card(cards, "IV", "3p")
    assert card["classification"] == "SOMATIC_DIMINUTION_BOUNDARY", \
        "a stub whose core reaches the tip must stay somatic, got {}".format(card["classification"])
    assert card.get("is_somatic_boundary") is True
    print("B ok: IV 3p core-reaching stub -> SOMATIC_DIMINUTION_BOUNDARY (unchanged)")

    # --- Case C: the walker reports ZERO_EXTENDING_READS with distinct germline vs somatic context notes
    #     when the unmapped read index is empty.
    class EmptyIndex:
        def __len__(self):
            return 0

    common = dict(seed_seq="ACGT" * 100, motif="TTAGGC", rev_motif="GCCTAA",
                  unmapped=EmptyIndex(), max_iters=1)
    germline = walk.run_iterative_walker("unused.bam", tip_is_germline=True, **common)
    assert germline["status"] == "ZERO_EXTENDING_READS" and "distal GRS arm" in germline["notes"], germline
    somatic = walk.run_iterative_walker("unused.bam", is_somatic_boundary=True, **common)
    assert somatic["status"] == "ZERO_EXTENDING_READS" and "somatic diminution breaksite" in somatic["notes"], somatic
    print("C ok: walker -> {} ({}) vs {} ({})".format(
        germline["status"], germline["notes"], somatic["status"], somatic["notes"]))

    # --- Case D: the net extension of a unitig cap is the outboard remainder, not the unitig.
    #     Real geometry from DF5120 III 5p: ptg000004l is 34,948 bp, aligning to III 253-11,756 at
    #     unitig coordinates 22,917-34,948. Only the 22,664 bp before the tip is new sequence.
    unitig = "A" * 22664 + "C" * (34948 - 22664)
    hit = {"q_aln_start": 22917, "q_aln_end": 34948, "t_start": 253, "t_end": 11756, "aln_strand": "+"}
    net = walk.outboard_extension(unitig, hit, "5p", 7052312)
    assert len(net) == 22664, "5p net extension should be 22,664 bp, got {}".format(len(net))
    net3 = walk.outboard_extension(unitig, hit, "3p", 11756)
    # The unitig's tail *is* the aligned part, so a 3p terminus at the alignment end gains nothing.
    assert len(net3) == 0, "expected no outboard sequence at the alignment end, got {}".format(len(net3))
    print("D ok: unitig 34,948 bp -> net 5p extension {} bp (was reporting the whole unitig)".format(len(net)))

    # --- Case E: telomeric bp are counted by scanning the whole added sequence, including arrays
    #     that sit in the middle of it, and in either orientation.
    mixed = "A" * 100 + "TTAGGC" * 50 + "G" * 100 + "GCCTAA" * 20 + "T" * 50
    assert walk.count_telomeric_bp(mixed, "TTAGGC", "GCCTAA") == 50 * 6 + 20 * 6, \
        "internal arrays must be counted: {}".format(walk.count_telomeric_bp(mixed, "TTAGGC", "GCCTAA"))
    assert walk.count_telomeric_bp("TTAGAC" * 10, "TTAGGC", "GCCTAA") == 60, "TTAGAC is canonical too"
    assert walk.count_telomeric_bp("ACGT" * 500, "TTAGGC", "GCCTAA") == 0, "non-telomeric must be 0"
    print("E ok: internal and reverse-orientation arrays counted (420 bp from two internal arrays)")

    # --- Case F: the metrics reported for an extension come from the sequence written.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        path, added, telo, non_telo = walk.write_extension_fasta(
            tmp, "III_5p", "unitig_cap", "unitig_id=ptg000004l",
            "N" * 1000 + "TTAGGC" * 100, "TTAGGC", "GCCTAA")
        assert (added, telo, non_telo) == (1600, 600, 1000), (added, telo, non_telo)
        assert pathlib.Path(path).exists()
    print("F ok: extension metrics = {} added / {} telomeric / {} non-telomeric".format(added, telo, non_telo))

    # --- Case G: the extension is the interval Module 2 authorised, in chromosome orientation.
    #     III 5p's card reports ptg000004l (34,948 bp) adding 33,055 bp, so the overlap is 1,893 bp.
    #     A forward alignment puts the extension at the unitig's head; a reverse alignment puts it
    #     at the tail and it must be complemented to appear in the chromosome.
    real_unitig = "".join("ACGT"[(i // 7) % 4] for i in range(34948))
    fwd = walk.region_extension(real_unitig, {"q_start": 0, "q_end": 33055, "strand": "+",
                                              "unitig_len": 34948})
    rev = walk.region_extension(real_unitig, {"q_start": 1893, "q_end": 34948, "strand": "-",
                                              "unitig_len": 34948})
    assert len(fwd) == 33055 and fwd == real_unitig[:33055], "forward interval mis-sliced"
    assert len(rev) == 33055 and rev == walk.revcomp(real_unitig[1893:34948]), \
        "reverse interval must be complemented into chromosome orientation"
    assert len(fwd) == len(rev) == 33055, "both orientations must yield the card's 33,055 bp"
    assert walk.region_extension(real_unitig, {}) == "", "absent region must yield nothing"
    assert walk.region_extension(real_unitig, {"q_start": 10, "q_end": 999999}) == "", \
        "an interval beyond the unitig must yield nothing rather than a partial extension"
    print("G ok: interval 33,055 bp both orientations (overlap {} bp) and not the 34,948 unitig"
          .format(34948 - 33055))

    # --- Case H: a card with no interval and no usable profile must not fall back to the unitig.
    #     This is the regression that produced 34,948: the row must report nothing added.
    no_region = {"unitig_extension_bp": "+33,055 bp", "telomeric_tract_bp": "0 bp"}
    assert walk.safe_int(no_region["unitig_extension_bp"]) == 33055, "card extension must parse"
    assert walk.region_extension(real_unitig, no_region.get("unitig_extension_region") or {}) == "", \
        "a card without an interval yields no extension, never the whole unitig"
    print("H ok: card without unitig_extension_region yields no extension (not 34,948 bp)")

    # --- Case I: the reported defect. Reads carry telomere at the tip, but the assembled tip holds no
    #     tract: the draft is NOT capped and the action is an overhang extension. DF5120 IV 5p and X 5p
    #     are this case -- baseline ALREADY_CAPPED from read clips (MAPQ 1-14) with terminal_tract_bp 0
    #     -- and they were classified DRAFT_HAS_TELOMERE, so the walker looked for a tract that is not
    #     there and finished UNRESOLVED while the blueprint endorsed a candidate for each.
    def _one_end_cards(baseline, tract, read_ext, unitigs=None):
        return syn.synthesize_telomere_cards(
            nuclear_chroms=["IV"],
            chrom_lens={"IV": 8418973},
            grs_segments={"IV": [(0, 3000, 3000, 12.0)]},
            core_segments={"IV": [(3000, 8418973, 8415973, 65.0)]},
            delfies_points={"IV": [3000]},
            rdna_intervals={},
            baseline_status=baseline,
            unitig_extensions=unitigs or {},
            single_read_ext=read_ext,
            tip_tract_bp=tract,
            policy=policy,
        )

    card = terminus_card(_one_end_cards({("IV", "5p"): "TELOMERE_IN_OVERHANG", ("IV", "3p"): "UNCAPPED"},
                                        {("IV", "5p"): 0},
                                        {("IV", "5p"): (0, 0, "None", 0)}), "IV", "5p")
    assert card["classification"] == "TELOMERE_IN_OVERHANG", card["classification"]
    assert card["action"] == "EXTEND_FROM_OVERHANG", card["action"]
    assert card["evidence"]["terminal_tract_bp"] == "0 bp", card["evidence"]
    print("I ok: telomere in the reads but not in the tip -> {} / {}".format(
        card["classification"], card["action"]))

    # --- Case I2: the same end WITH a tract in the tip keeps its tract-extension card, so the split
    #     relabels only the ends it should.
    card = terminus_card(_one_end_cards({("IV", "5p"): "ALREADY_CAPPED", ("IV", "3p"): "UNCAPPED"},
                                        {("IV", "5p"): 960},
                                        {("IV", "5p"): (1500, 1, "read_a", 150)}), "IV", "5p")
    assert card["classification"] == "DRAFT_HAS_TELOMERE", card["classification"]
    assert card["action"] == "EXTEND_TELOMERE_TRACT", card["action"]
    assert card["authorized_read_qname"] == "read_a", card
    print("I2 ok: tip holds 960 bp of tract -> {} / {}".format(
        card["classification"], card["action"]))

    # --- Case I3: a summary predating the tract column (no value for this end) must keep the old
    #     single-status behaviour rather than being guessed into the overhang branch.
    card = terminus_card(_one_end_cards({("IV", "5p"): "ALREADY_CAPPED", ("IV", "3p"): "UNCAPPED"},
                                        {},
                                        {("IV", "5p"): (1500, 1, "read_a", 150)}), "IV", "5p")
    assert card["classification"] == "DRAFT_HAS_TELOMERE", \
        "an unmeasured tip tract must not be reinterpreted: {}".format(card["classification"])
    print("I3 ok: no tract measurement available -> unchanged {} / {}".format(
        card["classification"], card["action"]))

    # --- Case K: the two tract measurements answer different questions and must not be conflated.
    #     Module 1 asks the taxonomy's question for this end -- does the tip contain canonical repeats
    #     at all (EXTEND_TELOMERE_TRACT) or is it unique sequence (EXTEND_TO_TELOMERE) -- so its
    #     threshold is presence. Module 3's terminal_tract_region applies the policy's min_repeats,
    #     because it measures a tract it is about to extend. With the thresholds matched they must
    #     agree; with the policy threshold a short run is legitimately 0 and the end stays A.
    sc = load(REPO / "modules/local/softclipends/resources/usr/bin/pde_terminal_softclips.py",
              "pde_terminal_softclips")
    motif, rev = "TTAGGC", "GCCTAA"
    for units in (2, 5, 15, 40):
        tip = "A" * 200 + motif * units
        m1 = sc.terminal_tract_bp(tip, motif, rev)
        m3 = walk.terminal_tract_region(tip, motif, rev, "3p", min_repeats=2)[0]
        assert m1 > 0 and m3 > 0, "both must see {} units: {} vs {}".format(units, m1, m3)
        assert m1 <= m3 + 12, "Module 1 must not claim more tract than Module 3: {} vs {}".format(m1, m3)
    short = "A" * 200 + motif * 5
    assert sc.terminal_tract_bp(short, motif, rev) > 0, "5 units is a tract for the A/B question"
    assert walk.terminal_tract_region(short, motif, rev, "3p")[0] == 0, \
        "5 units is below the policy's min_repeats, so Module 3 does not call it a tract to extend"
    bare = "A" * 200 + motif          # one unit is not a tract
    assert sc.terminal_tract_bp(bare, motif, rev) == 0, "a single unit is not a tract"
    assert sc.terminal_tract_bp("ACGT" * 300, motif, rev) == 0, "unique sequence holds no tract"
    # A tip stopping inside a repeat unit is still a tract (assemblies rarely stop on the boundary).
    partial = "A" * 200 + motif * 20 + motif[:4]
    assert sc.terminal_tract_bp(partial, motif, rev) >= 19 * 6, (
        "a tip inside a unit must still measure the run: {}".format(sc.terminal_tract_bp(partial, motif, rev))
    )
    print("K ok: tip-tract measurement agrees with Module 3 at matched thresholds, and the "
          "presence-vs-min_repeats difference is intentional")

    # --- Case L: 5p extensions in chromosome orientation place the true chromosome tip at base 1
    #     (index 0 / 5' end), whilst the 3' end is the internal junction joining the chromosome.
    #     terminal_tract_region and cap_status_for must measure the tract at the outward 5' tip for
    #     5p ends, and at the 3' tip for 3p ends.
    ext_5p = (rev * 150) + ("A" * 400)  # 900 bp telomere at 5' tip, 400 bp non-telo junction at 3' end
    t5_bp, _, insert_5p = walk.terminal_tract_region(ext_5p, motif, rev, min_repeats=15, terminus="5p")
    assert t5_bp == 900, "5p tract should measure 900 bp at the 5' tip, got {}".format(t5_bp)
    assert insert_5p == 0, "insert_at for 5p should be 0, got {}".format(insert_5p)
    status_5p, tract_5p = walk.cap_status_for("RESOLVED_BY_UNITIG", ext_5p, motif, rev, 15, terminus="5p")
    assert status_5p == "RESOLVED_BY_UNITIG" and tract_5p == 900, (status_5p, tract_5p)

    # The same 5p extension evaluated as 3p must report 0 (the 3' end is unique junction sequence).
    t5_as_3p = walk.terminal_tract_region(ext_5p, motif, rev, min_repeats=15, terminus="3p")[0]
    assert t5_as_3p == 0, "5p extension must not be mistaken for capped when evaluated as 3p"

    ext_3p = ("A" * 400) + (motif * 150)  # 400 bp non-telo junction at 5' end, 900 bp telomere at 3' tip
    t3_bp, _, insert_3p = walk.terminal_tract_region(ext_3p, motif, rev, min_repeats=15, terminus="3p")
    assert t3_bp == 900, "3p tract should measure 900 bp at the 3' tip, got {}".format(t3_bp)
    status_3p, tract_3p = walk.cap_status_for("RESOLVED_BY_UNITIG", ext_3p, motif, rev, 15, terminus="3p")
    assert status_3p == "RESOLVED_BY_UNITIG" and tract_3p == 900, (status_3p, tract_3p)

    # The same 3p extension evaluated as 5p must report 0 (the 5' end is unique junction sequence).
    t3_as_5p = walk.terminal_tract_region(ext_3p, motif, rev, min_repeats=15, terminus="5p")[0]
    assert t3_as_5p == 0, "3p extension must not be mistaken for capped when evaluated as 5p"

    # Subtelomeric extension with internal repeats (like DF5120 III 5p) must still report 0 at 5p.
    ext_subtelo = ("A" * 32800) + (motif * 36)
    t_sub = walk.terminal_tract_region(ext_subtelo, motif, rev, min_repeats=15, terminus="5p")[0]
    assert t_sub == 0, "internal repeats in a 5p extension must not be reported as a terminal tract"
    status_sub, tract_sub = walk.cap_status_for("RESOLVED_BY_UNITIG", ext_subtelo, motif, rev, 15, terminus="5p")
    assert status_sub == "EXTENDED_SUBTELOMERE_ONLY" and tract_sub == 0, (status_sub, tract_sub)
    print("L ok: 5p vs 3p tract orientation measured at respective outward chromosome tips")

    # --- Case M: PAF alignments for unitigs must respect policy mapq_anchored (30).
    #     For DF5120 III 5p, utg000007l had an authentic MAPQ 60 alignment (+ strand, q=22917..34948,
    #     overhang 22,917 bp) carrying the 846 bp terminal telomere, and a secondary repeat alignment
    #     (- strand, q=939..1893, MAPQ 1, overhang 33,055 bp). Without MAPQ filtering, the 33,055 bp
    #     bogus alignment won and excised the telomere. Enforcing MAPQ >= 30 ensures the authentic
    #     alignment is selected.
    with tempfile.NamedTemporaryFile("w", suffix=".paf") as tf:
        # Authentic MAPQ 60 alignment (+ strand)
        tf.write("utg000007l\t34948\t22917\t34948\t+\tIII\t7104590\t253\t11756\t10976\t12168\t60\n")
        # Bogus MAPQ 1 alignment (- strand)
        tf.write("utg000007l\t34948\t939\t1893\t-\tIII\t7104590\t4001\t4955\t941\t954\t1\n")
        tf.flush()

        unitig_exts = {}
        t_win = policy.get("telomere", {}).get("terminal_window_bp", 10000)
        min_ext = policy.get("telomere", {}).get("min_extension_bp", 200)
        mapq_anch = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30)
        with open(tf.name) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 12:
                    mapq = int(p[11])
                    if mapq < mapq_anch:
                        continue
                    qname, qlen, qstart, qend = p[0], int(p[1]), int(p[2]), int(p[3])
                    strand, tname, tlen, tstart, tend = p[4], p[5], int(p[6]), int(p[7]), int(p[8])
                    if tstart <= t_win:
                        oh = qstart if strand == "+" else (qlen - qend)
                        if oh >= min_ext:
                            unitig_exts.setdefault((tname, "5p"), []).append(
                                (oh, qname, 0, 0, qstart, qend, strand, qlen, tstart, tend))
        assert len(unitig_exts[("III", "5p")]) == 1, "only MAPQ >= 30 alignment should be accepted"
        best_cand = unitig_exts[("III", "5p")][0]
        assert best_cand[0] == 22917 and best_cand[6] == "+", "authentic + strand 22,917 bp alignment must win"
    print("M ok: PAF unitig candidate parser enforces mapq_anchored (rejects MAPQ 1 33,055 bp overhang)")

    # --- Case N: INVERT_AND_STITCH distinguishes near_coord (gap) from far_coord (distal end of
    #     inverted contig) in both directions (s_start near gap vs target_coord near gap), and
    #     pde_generate_agp.py excises the dummy N-gap while inverting only [inv_start, inv_end].
    from curationsynthesis.internal_resolver import synthesize_curation_cards
    import os
    import subprocess

    ev_inv = {
        "gaps": [
            {"chrom": "I", "start": 1000, "end": 1200, "gap_length": 200, "spanning_reads": []},
            {"chrom": "II", "start": 5000, "end": 5200, "gap_length": 200, "spanning_reads": []},
        ],
        "structural_variants": [
            # On chrom I: s_start (1200) is at the gap, end/target_coord (4200) is the distal end of the 3,000 bp inverted contig to the right
            {"id": "Sniffles2.INV.1", "chrom": "I", "start": 1200, "end": 4200, "svtype": "INV", "svlen": 3000, "support": 25, "vaf": 0.85, "alt": "<INV>"},
            # On chrom II: s_start (2000) is the distal end to the left, end/target_coord (5000) is at the gap
            {"id": "Sniffles2.INV.2", "chrom": "II", "start": 2000, "end": 5000, "svtype": "INV", "svlen": 3000, "support": 25, "vaf": 0.85, "alt": "<INV>"},
            # Homozygous novel insertion (INSERT_SEQUENCE) on chrom I at 4800
            {"id": "Sniffles2.INS.1", "chrom": "I", "start": 4800, "end": 4800, "svtype": "INS", "svlen": 800, "support": 30, "vaf": 0.95, "alt": "<INS>"},
            # Unverified DUP must NOT emit a phantom EXPAND_TANDEM_REPEAT card
            {"id": "Sniffles2.DUP.1", "chrom": "I", "start": 500, "end": 2500, "svtype": "DUP", "svlen": 2000, "support": 30, "vaf": 0.80, "alt": "<DUP>"},
        ],
    }
    int_cards = synthesize_curation_cards(ev_inv, policy=policy, chrom_lens={"I": 6000, "II": 8000})
    actions_by_chrom = {(c["chrom"], c["action"]): c for c in int_cards}
    assert ("I", "EXPAND_TANDEM_REPEAT") not in actions_by_chrom, "EXPAND_TANDEM_REPEAT phantom card must not be emitted"
    c_inv1 = actions_by_chrom[("I", "INVERT_AND_STITCH")]
    assert c_inv1["start"] == 1000 and c_inv1["end"] == 4200
    assert c_inv1["inv_start"] == 1200 and c_inv1["inv_end"] == 4200
    c_inv2 = actions_by_chrom[("II", "INVERT_AND_STITCH")]
    assert c_inv2["start"] == 2000 and c_inv2["end"] == 5200
    assert c_inv2["inv_start"] == 2000 and c_inv2["inv_end"] == 5000
    c_ins = actions_by_chrom[("I", "INSERT_SEQUENCE")]
    assert c_ins["start"] == 4800 and c_ins["end"] == 4800
    assert c_ins["evidence"]["max_spanning_insertion"] == "800 bp"

    # Run pde_generate_agp.py to verify INVERT_AND_STITCH excises the 200 bp N-gap and inverts only [1200, 4200],
    # and INSERT_SEQUENCE splices the 800 bp patch at 4800..4800.
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_fa = os.path.join(tmpdir, "ref.fa")
        left_flank = "ACGT" * 250           # 1000 bp (0..1000)
        dummy_gap = "N" * 200               # 200 bp (1000..1200)
        inv_block = ("AAACCC" * 500)        # 3000 bp (1200..4200)
        mid_flank = "TGCA" * 150            # 600 bp (4200..4800)
        right_flank = "GCTA" * 300          # 1200 bp (4800..6000)
        chrom_i_seq = left_flank + dummy_gap + inv_block + mid_flank + right_flank
        with open(ref_fa, "w") as fh:
            fh.write(">I\n" + chrom_i_seq + "\n")
        ins_patch_fa = os.path.join(tmpdir, "ins_patch.fa")
        ins_seq = "TTAA" * 200              # 800 bp
        with open(ins_patch_fa, "w") as fh:
            fh.write(">{}\n{}\n".format(c_ins["card_id"], ins_seq))
        patches_tsv = os.path.join(tmpdir, "patches.tsv")
        with open(patches_tsv, "w") as fh:
            fh.write("gap_id\tchrom\tstart\tend\torig_len\toutcome\tleft_anchor\tright_anchor\tpatch_len\texcised_len\tunitig_id\tpatch_fasta\tnotes\n")
            fh.write("{}\tI\t4800\t4800\t800\tPATCH_GAP_INSERTION\t2000\t2000\t800\t0\tabpoa_cons\t{}\tInserted\n".format(
                c_ins["card_id"], ins_patch_fa
            ))
        bp_json = os.path.join(tmpdir, "blueprint.json")
        import json
        with open(bp_json, "w") as fh:
            json.dump({"curation_cards": [c_inv1, c_ins], "telomere_cards": []}, fh)
        agp_script = str(REPO / "modules/local/curationagpgenerator/resources/usr/bin/pde_generate_agp.py")
        subprocess.run(
            [sys.executable, agp_script, "--sample", "TEST_INV", "--ref", ref_fa, "--blueprint", bp_json, "--gap-manifest", patches_tsv, "--outdir", tmpdir],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        comp_fa = os.path.join(tmpdir, "TEST_INV.components.fasta")
        comp_seqs = {}
        cur_k = None
        with open(comp_fa) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    cur_k = line[1:].split()[0]
                    comp_seqs[cur_k] = ""
                elif cur_k:
                    comp_seqs[cur_k] += line
        assert "I_1200_4200_inv" in comp_seqs, "inverted component must slice [1200, 4200] excluding dummy gap"
        assert "N" not in comp_seqs["I_1200_4200_inv"], "inverted component must not contain dummy gap Ns"
        assert len(comp_seqs["I_1200_4200_inv"]) == 3000
        assert comp_seqs[c_ins["card_id"] + "_patch"] == ins_seq
    print("N ok: INVERT_AND_STITCH excises dummy gap and inverts [inv_start, inv_end]; INSERT_SEQUENCE splices cleanly")

    # --- Case O: EXTEND_TELOMERE_TRACT / EXTEND_FROM_OVERHANG arbitrates unitig against valid_read_ext_bp
    #     (so an invalid non-telomeric soft-clip read cannot block authorized_unitig_id) and includes
    #     start, end, and unitig_extension_region.
    from curationsynthesis.terminus_resolver import synthesize_telomere_cards, detect_chimeric_terminus
    from telocapping.overhang_capper import overhang_of_read

    telo_cards_o = synthesize_telomere_cards(
        nuclear_chroms=["I"],
        chrom_lens={"I": 100000},
        grs_segments={"I": []},
        core_segments={"I": [(1, 100000, 100000, 100.0)]},
        delfies_points={},
        rdna_intervals={},
        baseline_status={("I", "5p"): "ALREADY_CAPPED", ("I", "3p"): "ALREADY_CAPPED"},
        unitig_extensions={("I", "3p"): [(1500, "utg000099l", 900, 600, 0, 8500, "+", 10000)]},
        # Raw read_ext_bp is 3000 (> 1500), but read_telo_cnt is 0 (invalid germline read!)
        single_read_ext={("I", "3p"): (3000, 1, "read_non_telo", 0)},
        tip_tract_bp={("I", "5p"): 600, ("I", "3p"): 600},
        policy=policy,
    )
    card_o = terminus_card(telo_cards_o, "I", "3p")
    assert card_o["action"] == "EXTEND_TELOMERE_TRACT"
    assert card_o["authorized_unitig_id"] == "utg000099l", "valid unitig must win over invalid non-telomeric read"
    assert card_o["authorized_read_qname"] is None
    assert card_o["start"] == 100000 and card_o["end"] == 100000
    assert card_o["evidence"]["unitig_extension_region"] == {"q_start": 8500, "q_end": 10000, "strand": "+", "unitig_len": 10000}
    print("O ok: EXTEND_TELOMERE_TRACT arbitrates unitig vs valid_read_ext_bp and carries unitig_extension_region")

    # --- Case P: detect_chimeric_terminus enforces min_repeats on 5p, and overhang_of_read respects chrom_start

    clips_5p = [
        {"qname": "r_12_reps", "align_start": 2000, "align_end": 8000, "clip_len": 900, "telo_cnt": 12, "mapq": 60},
    ]
    assert detect_chimeric_terminus("I", "5p", 100000, clips_5p, d_soma=0, delfies_pts=[], policy=policy) is None, (
        "5p read with 12 repeats (< min_repeats=15) must not trigger chimeric terminus pruning"
    )
    clips_5p[0]["telo_cnt"] = 20
    chim_5p = detect_chimeric_terminus("I", "5p", 100000, clips_5p, d_soma=0, delfies_pts=[], policy=policy)
    assert chim_5p is not None and chim_5p["junction_coord"] == 2000

    # Verify overhang_of_read at 5p with chrom_start=2000 (pos=2001, 1-based -> gap=0 relative to pruned junction 2000)
    sam_line_5p = "r_12_reps\t0\tI\t2001\t60\t300S1000M\t*\t0\t0\t" + ("G" * 300) + ("A" * 1000) + "\t*"
    oh_unpruned, _ = overhang_of_read(sam_line_5p, "5p", 100000, chrom_start=0)
    assert oh_unpruned == 0, "without chrom_start=2000, a 300 bp clip at pos=2001 cannot reach coordinate 1"
    oh_pruned, oh_seq = overhang_of_read(sam_line_5p, "5p", 100000, chrom_start=2000)
    assert oh_pruned == 300 and oh_seq == ("G" * 300), "with chrom_start=2000, all 300 bp overhang the pruned 5p junction"
    print("P ok: detect_chimeric_terminus enforces min_repeats on 5p and overhang_of_read honours chrom_start")

    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
