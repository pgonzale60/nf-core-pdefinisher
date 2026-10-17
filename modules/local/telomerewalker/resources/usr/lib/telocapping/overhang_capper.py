"""Overhang capping sub-module.

Handles single-read overhang extraction, soft-clip validation, and terminal candidate screening.
British spelling standard enforced.
"""

import re
import subprocess
from telotract.tract_detector import count_telomeric_bp, revcomp


def overhang_of_read(line, terminus, chrom_len, min_clip_bp=50, chrom_start=0):
    """The outward overhang of one SAM record: (bp, sequence in chromosome orientation) or (0, "").

    A read only overhangs the tip when the clip on the outward side is at least as long as the gap
    between its alignment and the tip (or `chrom_start` when a 5' chimeric prefix has been pruned).
    """
    parts = line.split("\t")
    if len(parts) < 11 or parts[5] == "*":
        return 0, ""
    try:
        flag = int(parts[1])
        pos = int(parts[3])
    except ValueError:
        return 0, ""
    if flag & 0x100 or flag & 0x800:
        return 0, ""
    ops = re.findall(r"(\d+)([MIDNSHP=X])", parts[5])
    if not ops:
        return 0, ""
    seq = parts[9]
    lead = int(ops[0][0]) if ops[0][1] in "SH" else 0
    tail = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
    ref_span = sum(int(n) for n, op in ops if op in "MDN=X")
    # In SAM format, SEQ (parts[9]) and CIGAR ops are stored oriented along the
    # reference from 5' to 3'. ops[0] (lead) is always the 5' clip on the reference
    # strand, and ops[-1] (tail) is always the 3' clip on the reference strand,
    # regardless of whether the query aligned to the forward or reverse strand.
    if terminus == "5p":
        clip = lead
        gap = max(0, pos - 1 - chrom_start)
        if clip < min_clip_bp or clip < gap:
            return 0, ""
        # Outward overhang extending beyond the effective 5' boundary
        overhang = seq[:lead - gap] if gap > 0 else seq[:lead]
    else:
        clip = tail
        align_end = pos - 1 + ref_span
        gap = max(0, chrom_len - align_end)
        if clip < min_clip_bp or clip < gap:
            return 0, ""
        # Outward overhang extending past chrom_len
        overhang = seq[len(seq) - (tail - gap):] if gap > 0 else seq[-tail:]

    if not overhang or len(overhang) < min_clip_bp:
        return 0, ""
    return len(overhang), overhang


def terminal_window_lines(bam_path, chrom, terminus, chrom_len, window, chrom_start=0):
    """SAM lines for the terminal window of a chromosome end."""
    if terminus == "5p":
        region = "{}:{}-{}".format(chrom, max(1, chrom_start + 1), chrom_start + window)
    else:
        region = "{}:{}-{}".format(chrom, max(1, chrom_len - window + 1), chrom_len)
    res = subprocess.run(["samtools", "view", bam_path, region], capture_output=True, text=True)
    if res.returncode != 0:
        return []
    return res.stdout.splitlines()


def maximal_tip_overhang(bam_path, chrom, terminus, chrom_len, window=10000, chrom_start=0):
    """The read whose alignment stops at this tip with the longest outward overhang."""
    best = None
    for line in terminal_window_lines(bam_path, chrom, terminus, chrom_len, window, chrom_start=chrom_start):
        clip, overhang = overhang_of_read(line, terminus, chrom_len, chrom_start=chrom_start)
        if clip and (best is None or clip > best[1]):
            best = (line.split("\t")[0], clip, overhang)
    return best


def named_read_overhang(bam_path, qname, chrom, terminus, chrom_len, window=10000, chrom_start=0):
    """The outward overhang of one nominated read, or (0, "")."""
    if not qname or qname == "None":
        return 0, ""
    for line in terminal_window_lines(bam_path, chrom, terminus, chrom_len, window, chrom_start=chrom_start):
        if line.split("\t")[0] == qname:
            return overhang_of_read(line, terminus, chrom_len, chrom_start=chrom_start)
    return 0, ""


def terminal_read_candidates(bam_path, chrom, terminus, chrom_len, motif, rev_motif,
                             window=10000, min_copies=5, chrom_start=0):
    """Reads near a terminus carrying any canonical telomeric repeat, for tract extension."""
    if terminus == "5p":
        region = "{}:{}-{}".format(chrom, max(1, chrom_start + 1), chrom_start + window)
    else:
        region = "{}:{}-{}".format(chrom, max(1, chrom_len - window), chrom_len)
    res = subprocess.run(["samtools", "view", bam_path, region], capture_output=True, text=True)
    if res.returncode != 0:
        return []
    candidates = []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        seq = parts[9]
        if max(seq.count(motif), seq.count(rev_motif)) >= min_copies:
            candidates.append({"qname": parts[0], "seq": seq})
    return candidates


def inspect_terminal_reads(bam_path, chrom, end_type, chrom_len, motif, rev_motif, window=10000, min_repeats=15, chrom_start=0):
    """Inspect terminal BAM reads for soft-clipped canonical telomere repeats."""
    if end_type == "5p":
        region = f"{chrom}:{max(1, chrom_start + 1)}-{chrom_start + window}"
    else:
        region = f"{chrom}:{max(1, chrom_len - window)}-{chrom_len}"

    cmd = ["samtools", "view", bam_path, region]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return None

    best_candidate = None
    max_telo_bp = 0

    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        cigar = parts[5]
        if "S" not in cigar and "H" not in cigar:
            continue
        seq = parts[9]
        telo_bp = count_telomeric_bp(seq, motif, rev_motif)
        if telo_bp > max_telo_bp and telo_bp >= (min_repeats * len(motif)):
            max_telo_bp = telo_bp
            best_candidate = {"qname": parts[0], "seq": seq, "telo_bp": telo_bp}

    return best_candidate


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


def extract_read_by_qname(bam_path, qname, chrom=None, end_type=None, chrom_len=None, window=100000, chrom_start=0):
    """Extract full query sequence of a read by QNAME from BAM using an indexed window."""
    if not qname or qname == "None":
        return None

    if chrom and end_type and chrom_len:
        if end_type == "5p":
            region = f"{chrom}:{max(1, chrom_start + 1)}-{chrom_start + window}"
        else:
            region = f"{chrom}:{max(1, chrom_len - window)}-{chrom_len}"
        cmd = ["samtools", "view", bam_path, region]
        res = subprocess.run(cmd, capture_output=True, text=True)
        for line in res.stdout.splitlines():
            p = line.split("\t")
            if len(p) >= 11 and p[0] == qname:
                return p[9]

    return None
