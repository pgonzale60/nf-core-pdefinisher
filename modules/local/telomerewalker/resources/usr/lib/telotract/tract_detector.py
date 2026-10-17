"""Tract detector sub-module.

Locates and quantifies canonical and variant telomeric repeat tracts in terminal sequences.
British spelling standard enforced.
"""

import re


def revcomp(seq):
    """Return the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return (seq or "").translate(trans)[::-1]


def count_telomeric_bp(seq, motif, rev_motif):
    """Bases of `seq` that form canonical tandem telomere arrays.

    Counted by scanning the whole sequence in both orientations, not by inspecting a terminal
    window and not by multiplying a repeat count by six.
    TTAGAC and GTCTAA are accepted alongside the given motif because the taxonomy names both as canonical.
    """
    variants = {motif, rev_motif, "TTAGAC", "GTCTAA"}
    pattern = re.compile("|".join("(?:%s)+" % v for v in sorted(variants) if v))
    return sum(len(m.group(0)) for m in pattern.finditer((seq or "").upper()))


def leading_telomere_bp(seq, motif, rev_motif):
    """Canonical telomeric bases at the start of `seq`: the tract as one read carries it."""
    variants = {motif, rev_motif, "TTAGAC", "GTCTAA"}
    pattern = re.compile("^(?:%s)+" % "|".join(sorted(v for v in variants if v)))
    match = pattern.match((seq or "").upper())
    return len(match.group(0)) if match else 0


def trailing_telomere_bp(seq, motif, rev_motif):
    """Canonical telomeric bases at the end of `seq`: the tract as the assembly currently has it."""
    return leading_telomere_bp(revcomp(seq or ""), motif, rev_motif)


def terminal_tract_region(terminal_seq, motif, rev_motif, terminus, min_repeats=15,
                          tip_tolerance_bp=100, gap_tolerance_bp=30):
    """Locate this end's telomeric tract in a sequence.

    For a 3' terminus (or when `terminal_seq` is oriented with the tip at its 3' end), the tip is at
    index `len(terminal_seq)`.
    For a 5' extension in chromosome orientation (prepended before coordinate 1), the true chromosome tip
    is at the 5' end (base 1 / index 0), whilst the 3' end is the internal junction that joins the
    chromosome.

    Returns (tract_bp, tract_sequence, insert_at).
    The tract is the run of repeat units closest to the tip in which no two units are further
    apart than `gap_tolerance_bp`.
    """
    variants = [v for v in dict.fromkeys([motif, rev_motif, "TTAGAC", "GTCTAA"]) if v]
    hits = []
    for variant in variants:
        start = terminal_seq.find(variant)
        while start >= 0:
            hits.append((start, start + len(variant)))
            start = terminal_seq.find(variant, start + 1)
    if not hits:
        return 0, "", (0 if terminus == "5p" else len(terminal_seq))
    hits.sort()

    runs = [list(hits[0])]
    for start, stop in hits[1:]:
        if start - runs[-1][1] <= gap_tolerance_bp:
            runs[-1][1] = max(runs[-1][1], stop)
        else:
            runs.append([start, stop])

    qualifying = []
    for start, stop in runs:
        tract_seq = terminal_seq[start:stop]
        if len(tract_seq) >= min_repeats * len(motif) and count_telomeric_bp(tract_seq, motif, rev_motif) >= 0.5 * len(tract_seq):
            qualifying.append((start, stop, tract_seq))

    if not qualifying:
        return 0, "", (0 if terminus == "5p" else len(terminal_seq))

    if terminus == "5p":
        start, stop, tract_seq = qualifying[0]
        if start > tip_tolerance_bp + gap_tolerance_bp:
            return 0, "", 0
        return len(tract_seq), tract_seq, start
    else:
        start, stop, tract_seq = qualifying[-1]
        if len(terminal_seq) - stop > tip_tolerance_bp + gap_tolerance_bp:
            return 0, "", len(terminal_seq)
        return len(tract_seq), tract_seq, stop


def outward_tail(read_seq, anchor, anchor_bp):
    """The read's sequence outward of the tract boundary, oriented for appending."""
    for orientation, probe in (("+", anchor), ("-", revcomp(anchor))):
        start = (read_seq or "").find(probe)
        if start < 0:
            continue
        if orientation == "+":
            return orientation, read_seq[start + anchor_bp:]
        return orientation, revcomp(read_seq[:start])
    return None, ""


def maximal_terminal_tract(candidates, motif, rev_motif, terminal_seq,
                           anchor_bp=500, min_extra_bp=60, min_repeats=15, tip_tolerance_bp=100):
    """Whichever read follows the telomeric tract furthest out, and how long the tract then is.

    Returns (read_name, total_tract_bp, extra_sequence, orientation, insert_at), or None when no read
    continues the tract.
    """
    current_tract, _tract_seq, insert_at = terminal_tract_region(
        terminal_seq, motif, rev_motif, "3p", min_repeats=min_repeats, tip_tolerance_bp=tip_tolerance_bp)
    if current_tract <= 0:
        return None
    anchor_end = insert_at - current_tract
    anchor = terminal_seq[max(0, anchor_end - anchor_bp):anchor_end]
    if len(anchor) < anchor_bp:
        return None

    best = None
    for cand in candidates:
        if isinstance(cand, dict):
            name, seq = cand["qname"], cand["seq"]
        else:
            name, seq = cand[0], cand[1]
        orientation, tail = outward_tail(seq, anchor, anchor_bp)
        if not tail:
            continue
        read_tract, read_tract_seq, _read_insert = terminal_tract_region(
            tail, motif, rev_motif, "3p", min_repeats=min_repeats, tip_tolerance_bp=tip_tolerance_bp)
        extra = read_tract - current_tract
        if extra >= min_extra_bp and (best is None or extra > best[1]):
            best = (name, extra, read_tract_seq[current_tract:current_tract + extra], orientation)
    if best is None:
        return None
    return best[0], current_tract + best[1], best[2], best[3], insert_at


def cap_status_for(base_status, extension, motif, rev_motif, min_repeats, terminus):
    """Differentiate true telomeric capping from subtelomeric extension."""
    tract_bp, _, _ = terminal_tract_region(
        extension, motif, rev_motif, terminus=terminus, min_repeats=min_repeats)
    if tract_bp > 0:
        return base_status, tract_bp
    return "EXTENDED_SUBTELOMERE_ONLY", 0
