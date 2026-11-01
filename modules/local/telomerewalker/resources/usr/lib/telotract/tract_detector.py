"""Tract detector sub-module.

Locates and quantifies canonical and variant telomeric repeat tracts in terminal sequences.
British spelling standard enforced.
"""

import os
import re
import shutil
import subprocess
import tempfile


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
                           anchor_bp=500, min_extra_bp=60, min_repeats=15, tip_tolerance_bp=100, min_reads=3):
    """Compute multi-read consensus of telomeric tract extension across supporting reads.

    Requires at least `min_reads` supporting reads matching the anchor with telomeric repeats
    (telo_len >= 200), and computes their abpoa consensus (falling back to the median supported
    tract if abpoa is unavailable in unit tests).
    Returns (read_name, total_tract_bp, extra_sequence, orientation, insert_at), or None.
    """
    current_tract, _tract_seq, insert_at = terminal_tract_region(
        terminal_seq, motif, rev_motif, "3p", min_repeats=min_repeats, tip_tolerance_bp=tip_tolerance_bp)
    if current_tract <= 0:
        return None
    anchor_end = insert_at - current_tract
    anchor = terminal_seq[max(0, anchor_end - anchor_bp):anchor_end]
    if len(anchor) < anchor_bp:
        return None

    qualifying_reads = []
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
        telo_len = count_telomeric_bp(read_tract_seq, motif, rev_motif)
        extra = read_tract - current_tract
        if extra >= min_extra_bp and (telo_len >= 200 or read_tract >= 200):
            qualifying_reads.append({
                "name": name,
                "extra": extra,
                "read_tract": read_tract,
                "extra_seq": read_tract_seq[current_tract:current_tract + extra],
                "orientation": orientation,
            })

    if len(qualifying_reads) < min_reads:
        return None

    # Compute abpoa consensus across qualifying reads
    abpoa_bin = shutil.which("abpoa")
    if abpoa_bin:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp_fa:
            tmp_fa_path = tmp_fa.name
            for idx, q in enumerate(qualifying_reads):
                tmp_fa.write(f">seq_{idx}_{q['name']}\n{q['extra_seq']}\n")

        try:
            cmd = [abpoa_bin, tmp_fa_path, "-r", "0"]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            cons_lines = [line.strip() for line in res.stdout.splitlines() if line.strip() and not line.startswith(">")]
            consensus_seq = "".join(cons_lines)
            if consensus_seq and len(consensus_seq) >= min_extra_bp:
                top_qname = max(qualifying_reads, key=lambda x: x["extra"])["name"]
                orientation = qualifying_reads[0]["orientation"]
                label = f"abpoa_consensus_{len(qualifying_reads)}reads_{top_qname}"
                return (label, current_tract + len(consensus_seq), consensus_seq, orientation, insert_at)
        except Exception:
            pass
        finally:
            if os.path.exists(tmp_fa_path):
                os.remove(tmp_fa_path)

    # Fallback to median supported tract when abpoa is unavailable
    sorted_reads = sorted(qualifying_reads, key=lambda x: x["extra"])
    chosen = sorted_reads[len(sorted_reads) // 2]
    label = f"median_consensus_{len(qualifying_reads)}reads_{chosen['name']}"
    return (label, current_tract + chosen["extra"], chosen["extra_seq"], chosen["orientation"], insert_at)


def cap_status_for(base_status, extension, motif, rev_motif, min_repeats, terminus):
    """Differentiate true telomeric capping from subtelomeric extension."""
    tract_bp, _, _ = terminal_tract_region(
        extension, motif, rev_motif, terminus=terminus, min_repeats=min_repeats)
    if tract_bp > 0:
        return base_status, tract_bp
    return "EXTENDED_SUBTELOMERE_ONLY", 0
