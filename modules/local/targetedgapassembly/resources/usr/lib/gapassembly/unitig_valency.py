"""Unitig Anchor Valency Conservation & Gap PAF Evaluation Engine for Targeted Gap Assembly.

Implements:
  - Genome-wide unitig valency survey (`survey_unitig_valencies_genome_wide`, `survey_repeat_anchor_valency`)
  - Deterministic PAF evaluation (`evaluate_gap_paf`) across standard bridges, asymmetric bridges,
    repetitive unique-anchored bridges, monovalent exclusion bridges, and monovalent overhang patches.

British spelling standard enforced across all documentation and comments.
"""

from collections import defaultdict
import os
import re
import subprocess
import sys
import tempfile


def parse_cigar_ops(cigar_str):
    """Parse CIGAR string into list of (length, op)."""
    if not cigar_str:
        return []
    return [(int(m.group(1)), m.group(2)) for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_str)]


def _interval_union_len(intervals):
    """Compute merged length of 0-based half-open intervals [(start, end), ...]."""
    clean = sorted((min(s, e), max(s, e)) for s, e in intervals if e > s)
    if not clean:
        return 0
    total = 0
    cur_s, cur_e = clean[0]
    for s, e in clean[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return total


def project_target_interval_via_cigar(a, t_left, t_right):
    """Map target interval [t_left, t_right] through a PAF alignment's CIGAR string
    to exact query coordinates (qs, qe). Falls back to linear anchor offsets if no CIGAR is present.
    """
    cigar_str = a.get("cigar")
    strand = a["strand"]
    tstart = a["tstart"]
    tend = a["tend"]
    qstart = a["qstart"]
    qend = a["qend"]

    left_anchor_len = max(0, t_left - tstart)
    right_anchor_len = max(0, tend - t_right)

    if not cigar_str:
        if strand == "+":
            return qstart + left_anchor_len, qend - right_anchor_len
        return qstart + right_anchor_len, qend - left_anchor_len

    ops = parse_cigar_ops(cigar_str)
    if not ops:
        if strand == "+":
            return qstart + left_anchor_len, qend - right_anchor_len
        return qstart + right_anchor_len, qend - left_anchor_len

    tpos = tstart
    q_off = 0
    q_off_left = None
    q_off_right = None

    if t_left <= tpos:
        q_off_left = 0
    if t_right <= tpos:
        q_off_right = 0

    for length, op in ops:
        if op in ("M", "=", "X"):
            next_t = tpos + length
            if q_off_left is None and t_left <= next_t:
                q_off_left = q_off + max(0, t_left - tpos)
            if q_off_right is None and t_right <= next_t:
                q_off_right = q_off + max(0, t_right - tpos)
            tpos = next_t
            q_off += length
        elif op in ("I", "S"):
            # An insertion right at t_right (gap_end_local) belongs inside the gap patch
            if q_off_right is not None and tpos == t_right:
                q_off_right += length
            q_off += length
        elif op in ("D", "N"):
            next_t = tpos + length
            if q_off_left is None and t_left <= next_t:
                q_off_left = q_off
            if q_off_right is None and t_right <= next_t:
                q_off_right = q_off
            tpos = next_t

    if q_off_left is None:
        q_off_left = q_off
    if q_off_right is None:
        q_off_right = q_off

    if strand == "+":
        qs = qstart + q_off_left
        qe = qstart + q_off_right
    else:
        qs = qend - q_off_right
        qe = qend - q_off_left

    return max(qstart, min(qs, qe)), min(qend, max(qs, qe))


def _find_suffix_prefix_overlap(left_seq, right_seq, min_ovl=100, max_ovl=10000, min_identity=0.98):
    """Find the longest suffix of `left_seq` matching a prefix of `right_seq` at >= min_identity."""
    if not left_seq or not right_seq:
        return 0
    limit = min(len(left_seq), len(right_seq), max_ovl)
    if limit < min_ovl:
        return 0
    # Use a 31-mer seed from the prefix of right_seq to test candidate overlap lengths fast
    seed_k = min(31, min_ovl)
    seed = right_seq[:seed_k]
    search_start = max(0, len(left_seq) - limit)
    pos = left_seq.find(seed, search_start)
    best_ovl = 0
    while pos != -1:
        ovl_len = len(left_seq) - pos
        if ovl_len >= min_ovl:
            s1 = left_seq[pos:]
            s2 = right_seq[:ovl_len]
            matches = sum(1 for a, b in zip(s1, s2) if a == b)
            if matches / ovl_len >= min_identity:
                return ovl_len
        pos = left_seq.find(seed, pos + 1)
    return best_ovl


def survey_unitig_valencies_genome_wide(
    ref_fasta,
    target_chrom,
    gap_start,
    gap_end,
    unitigs_fa,
    threads=4,
    max_offtarget_fraction=0.02,
    min_anchor_bp=1000,
    window_buffer=50000,
    min_identity=0.90,
):
    """Survey unitig alignments across the entire reference genome using minimap2 to evaluate
    whether unitigs are strictly monovalent to the gap locus (Principle of Anchor Valency Conservation).

    Uses merged interval unions (`_interval_union_len`) and overlap geometry (`tend >= win_start`
    and `tstart <= win_end`) so stacked tandem-repeat self-alignments at the gap cannot inflate
    on-target coverage and long unitigs extending past `window_buffer` are not miscounted as off-target.
    """
    monovalent = set()
    if not os.path.exists(unitigs_fa) or os.path.getsize(unitigs_fa) == 0:
        return monovalent

    cmd = [
        "minimap2",
        "-x",
        "asm20",
        "-c",
        "--secondary=yes",
        "-N",
        "50",
        "-t",
        str(threads),
        ref_fasta,
        unitigs_fa,
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except Exception as e:
        print(f"Valency survey warning: minimap2 whole-genome alignment failed: {e}", file=sys.stderr)
        return monovalent

    win_start = max(0, gap_start - window_buffer)
    win_end = gap_end + window_buffer

    on_target_q_intervals = defaultdict(list)
    off_target_t_intervals = defaultdict(lambda: defaultdict(list))

    for line in res.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        qname = parts[0]
        qstart = int(parts[2])
        qend = int(parts[3])
        tname = parts[5]
        tstart = int(parts[7])
        tend = int(parts[8])
        matches = int(parts[9])
        aln_len = int(parts[10])

        if aln_len < 150 or (matches / max(1, aln_len)) < min_identity:
            continue

        if tname == target_chrom and tend >= win_start and tstart <= win_end:
            on_target_q_intervals[qname].append((qstart, qend))
        else:
            off_target_t_intervals[qname][tname].append((tstart, tend))

    for qname, q_ivs in on_target_q_intervals.items():
        on_len = _interval_union_len(q_ivs)
        off_len = sum(
            _interval_union_len(t_ivs) for t_ivs in off_target_t_intervals[qname].values()
        )
        total_len = on_len + off_len
        off_frac = off_len / total_len if total_len > 0 else 1.0
        if on_len >= min_anchor_bp and off_frac <= max_offtarget_fraction:
            monovalent.add(qname)

    return monovalent


def survey_repeat_anchor_valency(
    ref_fasta,
    target_chrom,
    gap_start,
    gap_end,
    candidate_seq,
    kmer_size=21,
    max_offtarget_fraction=0.02,
    min_repeat_copies=2,
    threads=4,
    min_anchor_bp=1000,
):
    """Evaluate whether a candidate repeat sequence is monovalent to the target gap locus."""
    if not candidate_seq:
        return False, "sequence_empty"
    if not os.path.exists(ref_fasta):
        return False, "ref_not_found"
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".fa") as tf:
        tf.write(f">candidate\n{candidate_seq}\n")
        tf.flush()
        mono_set = survey_unitig_valencies_genome_wide(
            ref_fasta,
            target_chrom,
            gap_start,
            gap_end,
            tf.name,
            threads=threads,
            max_offtarget_fraction=max_offtarget_fraction,
            min_anchor_bp=min_anchor_bp,
        )
    os.unlink(tf.name)
    if "candidate" in mono_set:
        return True, "monovalent_solitary_locus (in_degree=1, out_degree=1)"
    return False, "polyvalent_or_unanchored"


def evaluate_gap_paf(
    paf_path,
    left_flank_len,
    gap_len,
    right_flank_len,
    g_start=None,
    g_end=None,
    min_mapq=50,
    min_aln_len=1000,
    min_secondary_anchor=500,
    min_identity=0.98,
    repetitive_mapq_floor=0,
    max_overlap_excision_bp=50000,
    require_unique_anchor=True,
    mapq_anchored=30,
    min_repetitive_anchor=1000,
    monovalent_unitigs=None,
    allow_monovalent_overhang_patch=True,
    min_monovalent_overhang=200,
    unitig_seqs=None,
):
    """Evaluate alignment of unitigs against [left_flank][gap][right_flank]."""
    if g_start is None:
        g_start = left_flank_len
    if g_end is None:
        g_end = left_flank_len + gap_len

    if not os.path.exists(paf_path) or os.path.getsize(paf_path) == 0:
        return {
            "outcome": "UNRESOLVED_GAP",
            "left_anchor": 0,
            "right_anchor": 0,
            "patch_len": 0,
            "excised_len": 0,
            "unitig_id": None,
            "notes": "No PAF alignments generated",
        }

    gap_start_local = left_flank_len
    gap_end_local = left_flank_len + gap_len

    if monovalent_unitigs is None:
        monovalent_unitigs = set()

    ambiguous_repeats_seen = []
    inversion_candidates = []

    unitigs = defaultdict(list)
    with open(paf_path) as paf:
        for line in paf:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            cigar = next((p[5:] for p in parts[12:] if p.startswith("cg:Z:")), None)
            unitigs[qname].append(
                {
                    "qlen": int(parts[1]),
                    "qstart": int(parts[2]),
                    "qend": int(parts[3]),
                    "strand": parts[4],
                    "tstart": int(parts[7]),
                    "tend": int(parts[8]),
                    "matches": int(parts[9]),
                    "aln_len": int(parts[10]),
                    "mapq": int(parts[11]),
                    "cigar": cigar,
                }
            )

    for qname, alns in unitigs.items():
        for a in alns:
            # 1. CIGAR deletion spanning across the gap (Overlap Artifact)
            if a.get("cigar"):
                tpos = a["tstart"]
                for c_len, c_op in parse_cigar_ops(a["cigar"]):
                    if c_op in ("M", "=", "X"):
                        tpos += c_len
                    elif c_op in ("D", "N"):
                        d_start = tpos
                        d_end = tpos + c_len
                        tpos += c_len
                        if gap_len > 0 and d_start <= gap_start_local and d_end >= gap_end_local and c_len > gap_len:
                            qs_del, qe_del = project_target_interval_via_cigar(a, d_start, d_end)
                            if (qe_del - qs_del) == 0:
                                l_anch = d_start - a["tstart"]
                                r_anch = a["tend"] - d_end
                                has_anchors = (
                                    (l_anch >= min_aln_len and r_anch >= min_aln_len and a["mapq"] >= min_mapq)
                                    or (
                                        l_anch >= 2 * min_aln_len
                                        and r_anch >= min_secondary_anchor
                                        and a["mapq"] >= min_mapq
                                    )
                                    or (
                                        r_anch >= 2 * min_aln_len
                                        and l_anch >= min_secondary_anchor
                                        and a["mapq"] >= min_mapq
                                    )
                                    or (
                                        l_anch >= min_secondary_anchor
                                        and r_anch >= min_secondary_anchor
                                        and a["mapq"] >= repetitive_mapq_floor
                                    )
                                )
                                if has_anchors and c_len <= max_overlap_excision_bp:
                                    del_left_bp = gap_start_local - d_start
                                    del_right_bp = d_end - gap_end_local
                                    excise_start = g_start - del_left_bp
                                    excise_end = g_end + del_right_bp
                                    excised_len = excise_end - excise_start
                                    return {
                                        "outcome": "EXCISE_OVERLAP_ARTIFACT",
                                    "left_anchor": l_anch,
                                    "right_anchor": r_anch,
                                    "patch_len": 0,
                                    "excised_len": excised_len,
                                    "unitig_id": qname,
                                    "excise_start": excise_start,
                                    "excise_end": excise_end,
                                    "notes": f"Unitig {qname} proves {excised_len - gap_len} bp false duplication overlap artefact across gap; excise_range={excise_start}-{excise_end}",
                                }

            left_anchor_len = gap_start_local - a["tstart"] if a["tstart"] < gap_start_local else 0
            right_anchor_len = a["tend"] - gap_end_local if a["tend"] > gap_end_local else 0
            identity = a["matches"] / max(1, a["aln_len"])
            mapq = a["mapq"]

            standard_bridge = (
                left_anchor_len >= min_aln_len and right_anchor_len >= min_aln_len and mapq >= min_mapq
            )
            asymmetric_bridge = (
                (
                    (left_anchor_len >= 2 * min_aln_len and right_anchor_len >= min_secondary_anchor)
                    or (right_anchor_len >= 2 * min_aln_len and left_anchor_len >= min_secondary_anchor)
                )
                and mapq >= mapq_anchored
                and identity >= min_identity
            )

            is_monovalent = qname in monovalent_unitigs
            has_unique_anchor = (
                mapq >= mapq_anchored
                or left_anchor_len >= min_aln_len
                or right_anchor_len >= min_aln_len
                or is_monovalent
            )
            meets_rep_anchors = (
                left_anchor_len >= min_secondary_anchor
                and right_anchor_len >= min_secondary_anchor
                and identity >= min_identity
                and mapq >= repetitive_mapq_floor
            )

            if (
                meets_rep_anchors
                and require_unique_anchor
                and not is_monovalent
                and not (has_unique_anchor and max(left_anchor_len, right_anchor_len) >= min_repetitive_anchor)
            ):
                ambiguous_repeats_seen.append((qname, left_anchor_len, right_anchor_len, mapq, identity))

            repetitive_bridge = meets_rep_anchors and (
                is_monovalent
                or not require_unique_anchor
                or (has_unique_anchor and max(left_anchor_len, right_anchor_len) >= min_repetitive_anchor)
            )

            if standard_bridge or asymmetric_bridge or repetitive_bridge:
                strand = a["strand"]
                qs, qe = project_target_interval_via_cigar(a, gap_start_local, gap_end_local)
                patch_len = max(0, qe - qs)
                bridge_mode = (
                    "standard"
                    if standard_bridge
                    else (
                        "asymmetric"
                        if asymmetric_bridge
                        else ("monovalent anchor valency (P=1.0)" if is_monovalent else "repetitive flank")
                    )
                )
                return {
                    "outcome": "PATCH_GAP_INSERTION",
                    "left_anchor": left_anchor_len,
                    "right_anchor": right_anchor_len,
                    "patch_len": patch_len,
                    "excised_len": gap_len,
                    "unitig_id": qname,
                    "q_patch_start": qs,
                    "q_patch_end": qe,
                    "strand": strand,
                    "notes": f"Unitig {qname} continuously bridges both flanks ({bridge_mode}, MAPQ={mapq}, id={identity:.3f})",
                }

        # Check split/collinear alignments from same unitig
        left_hits = [
            a
            for a in alns
            if (gap_start_local - a["tstart"] >= min_secondary_anchor)
            and a["tend"] <= gap_end_local + gap_len
            and (a["matches"] / max(1, a["aln_len"]) >= min_identity or a["mapq"] >= mapq_anchored)
            and a["mapq"] >= repetitive_mapq_floor
        ]
        right_hits = [
            a
            for a in alns
            if a["tstart"] >= gap_start_local - gap_len
            and (a["tend"] - gap_end_local >= min_secondary_anchor)
            and (a["matches"] / max(1, a["aln_len"]) >= min_identity or a["mapq"] >= mapq_anchored)
            and a["mapq"] >= repetitive_mapq_floor
        ]

        if left_hits and right_hits:
            best_l = max(left_hits, key=lambda x: x["aln_len"])
            best_r = max(right_hits, key=lambda x: x["aln_len"])
            l_anchor = gap_start_local - best_l["tstart"]
            r_anchor = best_r["tend"] - gap_end_local

            if best_l["strand"] != best_r["strand"]:
                inversion_candidates.append((qname, best_l["strand"], best_r["strand"]))
                continue

            has_anchors = (
                (l_anchor >= min_aln_len and r_anchor >= min_aln_len)
                or (l_anchor >= 2 * min_aln_len and r_anchor >= min_secondary_anchor)
                or (r_anchor >= 2 * min_aln_len and l_anchor >= min_secondary_anchor)
                or (
                    l_anchor >= min_secondary_anchor
                    and r_anchor >= min_secondary_anchor
                    and (best_l["matches"] / max(1, best_l["aln_len"]) >= min_identity)
                    and (best_r["matches"] / max(1, best_r["aln_len"]) >= min_identity)
                )
            )

            is_monovalent = qname in monovalent_unitigs
            has_unique_split = (
                best_l["mapq"] >= mapq_anchored
                or best_r["mapq"] >= mapq_anchored
                or l_anchor >= min_aln_len
                or r_anchor >= min_aln_len
                or is_monovalent
            )

            if (
                has_anchors
                and require_unique_anchor
                and not is_monovalent
                and not (has_unique_split and max(l_anchor, r_anchor) >= min_repetitive_anchor)
            ):
                ambiguous_repeats_seen.append(
                    (
                        qname,
                        l_anchor,
                        r_anchor,
                        max(best_l["mapq"], best_r["mapq"]),
                        min(
                            best_l["matches"] / max(1, best_l["aln_len"]),
                            best_r["matches"] / max(1, best_r["aln_len"]),
                        ),
                    )
                )
                continue

            if has_anchors and (
                is_monovalent
                or not require_unique_anchor
                or (has_unique_split and max(l_anchor, r_anchor) >= min_repetitive_anchor)
            ):
                strand = best_l["strand"]
                if strand == "+":
                    qs, qe = best_l["qend"], best_r["qstart"]
                else:
                    qs, qe = best_r["qend"], best_l["qstart"]
                patch_span = max(0, qe - qs)
                if (strand == "+" and best_l["qend"] <= best_r["qstart"]) or (
                    strand == "-" and best_r["qend"] <= best_l["qstart"]
                ):
                    val_str = " (monovalent anchor valency, P=1.0)" if is_monovalent else ""
                    return {
                        "outcome": "PATCH_GAP_INSERTION",
                        "left_anchor": l_anchor,
                        "right_anchor": r_anchor,
                        "patch_len": patch_span,
                        "excised_len": gap_len,
                        "unitig_id": qname,
                        "q_patch_start": qs,
                        "q_patch_end": qe,
                        "strand": strand,
                        "notes": f"Unitig {qname} collinear pair bridges gap with {patch_span} bp novel insertion{val_str}",
                    }

    if inversion_candidates:
        qname, sl, sr = inversion_candidates[0]
        return {
            "outcome": "UNRESOLVED_GAP",
            "left_anchor": 0,
            "right_anchor": 0,
            "patch_len": 0,
            "excised_len": 0,
            "unitig_id": qname,
            "notes": f"Inversion breakpoint detected at gap boundary ({qname} discordant strands: left {sl}, right {sr}); gap preserved",
        }

    if ambiguous_repeats_seen:
        qname, la, ra, mq, _idt = ambiguous_repeats_seen[0]
        return {
            "outcome": "UNRESOLVED_GAP",
            "left_anchor": la,
            "right_anchor": ra,
            "patch_len": 0,
            "excised_len": 0,
            "unitig_id": qname,
            "notes": f"Ambiguous repetitive bridge rejected ({qname} MAPQ={mq} < {mapq_anchored}, anchors L={la} R={ra} lack unique sequence >= {min_aln_len} bp); repeat may be multi-locus; gap preserved",
        }

    # Check if unitigs terminate at repeat boundaries without a single unitig spanning both walls
    left_wall = any(
        abs(gap_start_local - a["tend"]) <= 200
        for alns in unitigs.values()
        for a in alns
        if a["tstart"] < gap_start_local
    )
    right_wall = any(
        abs(a["tstart"] - gap_end_local) <= 200
        for alns in unitigs.values()
        for a in alns
        if a["tend"] > gap_end_local
    )
    if left_wall and right_wall:
        if allow_monovalent_overhang_patch and monovalent_unitigs:
            best_left_oh = None
            best_right_oh = None
            for qname in monovalent_unitigs:
                if qname in unitigs:
                    for a in unitigs[qname]:
                        if (
                            abs(a["tend"] - gap_start_local) <= 200
                            and a["aln_len"] >= min_repetitive_anchor
                            and a["tstart"] < gap_start_local
                        ):
                            overhang = a["qstart"] if a["strand"] == "-" else (a["qlen"] - a["qend"])
                            if overhang >= min_monovalent_overhang:
                                qs = 0 if a["strand"] == "-" else a["qend"]
                                qe = a["qstart"] if a["strand"] == "-" else a["qlen"]
                                cand = {
                                    "outcome": "PATCH_GAP_INSERTION",
                                    "left_anchor": a["aln_len"],
                                    "right_anchor": 0,
                                    "patch_len": overhang,
                                    "excised_len": gap_len,
                                    "unitig_id": qname,
                                    "q_patch_start": qs,
                                    "q_patch_end": qe,
                                    "strand": a["strand"],
                                    "notes": f"Monovalent repeat array extension ({qname}, MAPQ={a['mapq']}, anchor={a['aln_len']} bp) bridges unbridged repeat wall into gap by {overhang} bp (P=1.0)",
                                }
                                if not best_left_oh or overhang > best_left_oh["patch_len"]:
                                    best_left_oh = cand
                        if (
                            abs(a["tstart"] - gap_end_local) <= 200
                            and a["aln_len"] >= min_repetitive_anchor
                            and a["tend"] > gap_end_local
                        ):
                            overhang = (a["qlen"] - a["qend"]) if a["strand"] == "-" else a["qstart"]
                            if overhang >= min_monovalent_overhang:
                                qs = a["qend"] if a["strand"] == "-" else 0
                                qe = a["qlen"] if a["strand"] == "-" else a["qstart"]
                                cand = {
                                    "outcome": "PATCH_GAP_INSERTION",
                                    "left_anchor": 0,
                                    "right_anchor": a["aln_len"],
                                    "patch_len": overhang,
                                    "excised_len": gap_len,
                                    "unitig_id": qname,
                                    "q_patch_start": qs,
                                    "q_patch_end": qe,
                                    "strand": a["strand"],
                                    "notes": f"Monovalent repeat array extension ({qname}, MAPQ={a['mapq']}, anchor={a['aln_len']} bp) bridges unbridged repeat wall backwards into gap by {overhang} bp (P=1.0)",
                                }
                                if not best_right_oh or overhang > best_right_oh["patch_len"]:
                                    best_right_oh = cand

            # If both left and right monovalent overhangs exist and sequences are provided,
            # check whether their overhangs overlap in '+' chromosome orientation to form a complete dual bridge!
            if best_left_oh and best_right_oh and unitig_seqs:
                comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

                def _rc(s):
                    return s.translate(comp)[::-1]

                u_l = best_left_oh["unitig_id"]
                u_r = best_right_oh["unitig_id"]
                if u_l in unitig_seqs and u_r in unitig_seqs:
                    raw_l = unitig_seqs[u_l][best_left_oh["q_patch_start"] : best_left_oh["q_patch_end"]]
                    seq_l = _rc(raw_l) if best_left_oh["strand"] == "-" else raw_l
                    raw_r = unitig_seqs[u_r][best_right_oh["q_patch_start"] : best_right_oh["q_patch_end"]]
                    seq_r = _rc(raw_r) if best_right_oh["strand"] == "-" else raw_r
                    ovl_bp = _find_suffix_prefix_overlap(seq_l, seq_r, min_ovl=100, min_identity=min_identity)
                    if ovl_bp >= 100:
                        merged_seq = seq_l + seq_r[ovl_bp:]
                        return {
                            "outcome": "PATCH_GAP_INSERTION",
                            "left_anchor": best_left_oh["left_anchor"],
                            "right_anchor": best_right_oh["right_anchor"],
                            "patch_len": len(merged_seq),
                            "excised_len": gap_len,
                            "unitig_id": f"{u_l}+{u_r}",
                            "patch_seq": merged_seq,
                            "strand": "+",
                            "notes": (
                                f"Dual monovalent repeat overhangs ({u_l} L={len(seq_l)} bp, {u_r} R={len(seq_r)} bp) "
                                f"overlap by {ovl_bp} bp to bridge repeat array ({len(merged_seq)} bp net, P=1.0)"
                            ),
                        }

            best_candidate = None
            if best_left_oh and best_right_oh:
                best_candidate = (
                    best_left_oh
                    if best_left_oh["patch_len"] >= best_right_oh["patch_len"]
                    else best_right_oh
                )
            else:
                best_candidate = best_left_oh or best_right_oh
            if best_candidate:
                return best_candidate

    # Test for Partial Extensions (only after inversion, ambiguous repeat, and monovalent repeat wall checks)
    left_ext = None
    right_ext = None
    for qname, alns in unitigs.items():
        for a in alns:
            if a["mapq"] >= mapq_anchored and a["aln_len"] >= min_aln_len:
                if a["tstart"] <= gap_start_local - min_aln_len and a["tend"] > gap_start_local:
                    ext_bp = a["tend"] - gap_start_local
                    if not left_ext or ext_bp > left_ext["ext_bp"]:
                        left_ext = {"unitig_id": qname, "ext_bp": ext_bp, "aln": a}
                if a["tend"] >= gap_end_local + min_aln_len and a["tstart"] < gap_end_local:
                    ext_bp = gap_end_local - a["tstart"]
                    if not right_ext or ext_bp > right_ext["ext_bp"]:
                        right_ext = {"unitig_id": qname, "ext_bp": ext_bp, "aln": a}

    if left_ext and right_ext:
        return {
            "outcome": "EXTEND_BOTH_FLANKS",
            "left_anchor": left_ext["ext_bp"],
            "right_anchor": right_ext["ext_bp"],
            "patch_len": left_ext["ext_bp"] + right_ext["ext_bp"],
            "excised_len": 0,
            "unitig_id": f"{left_ext['unitig_id']};{right_ext['unitig_id']}",
            "notes": f"Left extended by {left_ext['ext_bp']} bp; right extended by {right_ext['ext_bp']} bp",
        }
    if left_ext:
        return {
            "outcome": "EXTEND_LEFT_FLANK",
            "left_anchor": left_ext["ext_bp"],
            "right_anchor": 0,
            "patch_len": left_ext["ext_bp"],
            "excised_len": 0,
            "unitig_id": left_ext["unitig_id"],
            "notes": f"Left flank extended inward by {left_ext['ext_bp']} bp",
        }
    if right_ext:
        return {
            "outcome": "EXTEND_RIGHT_FLANK",
            "left_anchor": 0,
            "right_anchor": right_ext["ext_bp"],
            "patch_len": right_ext["ext_bp"],
            "excised_len": 0,
            "unitig_id": right_ext["unitig_id"],
            "notes": f"Right flank extended inward by {right_ext['ext_bp']} bp",
        }

    if left_wall and right_wall:
        return {
            "outcome": "UNRESOLVED_GAP",
            "left_anchor": 0,
            "right_anchor": 0,
            "patch_len": 0,
            "excised_len": 0,
            "unitig_id": None,
            "notes": "Unbridged repeat array; contigs terminate at repeat boundary without spanning; gap preserved to protect unplaced intervening elements",
        }

    return {
        "outcome": "UNRESOLVED_GAP",
        "left_anchor": 0,
        "right_anchor": 0,
        "patch_len": 0,
        "excised_len": 0,
        "unitig_id": None,
        "notes": "No unitigs anchored with sufficient length and MAPQ >= 50",
    }
