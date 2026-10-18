"""Dentist Improper-Alignment Guard for Targeted Gap Assembly.

Calculates Dentist's improper-alignment threshold I_max and audits flank boundaries
for excessive diverging soft-clipped alignments to prevent ectopic repeat chimeras.
Formula:
    I_max = floor(0.5 * x + exp(0.1875 * (8 - x)))
where x is the local read depth at the boundary.
"""

import math
import re
import subprocess


def dentist_max_improper(read_depth):
    """Dentist's maximum allowable improper alignments at flank boundary."""
    x = max(0.0, float(read_depth))
    return math.floor(0.5 * x + math.exp(0.1875 * (8.0 - x)))


def count_flank_improper_alignments(bam_path, chrom, flank_coord, side, min_mapq=20, min_clip_bp=100, window_bp=500):
    """Count total read depth and diverging soft-clipped reads at a gap boundary.

    side='LEFT': boundary is at flank_coord (read extends leftwards into flank; tail soft-clip into gap).
    side='RIGHT': boundary is at flank_coord (read extends rightwards into flank; lead soft-clip into gap).
    """
    if side == "LEFT":
        region = f"{chrom}:{max(1, flank_coord - window_bp)}-{flank_coord}"
    else:
        region = f"{chrom}:{flank_coord}-{flank_coord + window_bp}"

    cmd = ["samtools", "view", "-F", "0x904", "-q", str(min_mapq), bam_path, region]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return 0, 0

    depth = 0
    improper = 0

    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        cigar = parts[5]
        if cigar == "*":
            continue

        depth += 1
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if not ops:
            continue

        # In SAM format:
        # ops[0] is 5' clip on reference strand (leftmost)
        # ops[-1] is 3' clip on reference strand (rightmost)
        if side == "LEFT":
            # Tail clip extending past the left flank into the gap
            tail_clip = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
            if tail_clip >= min_clip_bp:
                improper += 1
        else:
            # Lead clip extending before the right flank into the gap
            lead_clip = int(ops[0][0]) if ops[0][1] in "SH" else 0
            if lead_clip >= min_clip_bp:
                improper += 1

    return depth, improper


def evaluate_dentist_guard(bam_path, chrom, g_start, g_end, min_mapq=20, min_clip_bp=100):
    """Audit both flank boundaries against Dentist's I_max threshold.

    Returns (is_vetoed, left_stats, right_stats).
    """
    left_depth, left_improper = count_flank_improper_alignments(
        bam_path, chrom, g_start, side="LEFT", min_mapq=min_mapq, min_clip_bp=min_clip_bp
    )
    right_depth, right_improper = count_flank_improper_alignments(
        bam_path, chrom, g_end, side="RIGHT", min_mapq=min_mapq, min_clip_bp=min_clip_bp
    )

    left_imax = dentist_max_improper(left_depth)
    right_imax = dentist_max_improper(right_depth)

    left_excess = left_improper > left_imax
    right_excess = right_improper > right_imax

    left_stats = {
        "depth": left_depth,
        "improper": left_improper,
        "i_max": left_imax,
        "excess": left_excess,
    }
    right_stats = {
        "depth": right_depth,
        "improper": right_improper,
        "i_max": right_imax,
        "excess": right_excess,
    }

    is_vetoed = left_excess and right_excess
    return is_vetoed, left_stats, right_stats
