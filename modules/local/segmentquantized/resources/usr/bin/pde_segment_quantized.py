#!/usr/bin/env python3
"""Final three-class segmentation from base-resolution quantized depth.

Canonical chromosomes are tiled from mosdepth --quantize runs: HIGH runs are putative_core;
LOW runs are putative_GRS when a clustered DELFIES breakpoint sits at a border (within
anchor_proximity_bp), otherwise putative_uncertain. Assembly gaps and repeat-dominated LOW
stretches without breakpoint support are putative_uncertain. Non-canonical scaffolds keep the
call from pde_identify_grs.py. Every base of a canonical chromosome gets exactly one label.
"""

import argparse
import gzip
import json
import os
import statistics
import sys
from pathlib import Path

__version__ = "1.0.0"

SUMMARY_HEADER = "chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype"


def find_curation_policy(policy_path=None):
    """Locate assets/curation_policy.json, the Single Source of Truth."""
    if policy_path and os.path.exists(policy_path):
        return os.path.abspath(policy_path)
    if "CURATION_POLICY" in os.environ and os.path.exists(os.environ["CURATION_POLICY"]):
        return os.path.abspath(os.environ["CURATION_POLICY"])
    current = Path(__file__).resolve()
    for parent in [current] + list(current.parents):
        candidate = parent / "assets" / "curation_policy.json"
        if candidate.is_file():
            return str(candidate)
        candidate_direct = parent / "curation_policy.json"
        if candidate_direct.is_file():
            return str(candidate_direct)
    return None


def load_curation_policy(policy_path=None):
    """Load curation policy parameters strictly from assets/curation_policy.json."""
    resolved_path = find_curation_policy(policy_path)
    if not resolved_path:
        raise FileNotFoundError(
            "Could not locate assets/curation_policy.json (the Single Source of Truth). "
            "Policy parameters must not be hardcoded or duplicated across components."
        )
    with open(resolved_path, "r", encoding="utf-8") as f:
        return json.load(f)


def read_lines(path):
    """Return lines of an optional input; a missing or empty file is no evidence."""
    if not path or not os.path.exists(path):
        return []
    open_fn = gzip.open if path.endswith(".gz") else open
    with open_fn(path, "rt") as f:
        return f.read().splitlines()


def merge_intervals(ivs):
    """Sort and merge (start, end) intervals."""
    ivs = sorted(ivs)
    merged = []
    for s, e in ivs:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def read_red(path):
    """RED repeat intervals per chromosome, merged: {chrom: [(s, e)]}."""
    repeats = {}
    for line in read_lines(path):
        cols = line.split()
        if len(cols) >= 3 and not line.startswith("#"):
            try:
                repeats.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))
            except ValueError:
                continue
    return {chrom: merge_intervals(ivs) for chrom, ivs in repeats.items()}


def covered_bp(s, e, red_ivs):
    """Bases of [s, e) covered by merged, sorted RED intervals."""
    total = 0
    for a, b in red_ivs or []:
        if a >= e:
            break
        lo, hi = max(a, s), min(b, e)
        if hi > lo:
            total += hi - lo
    return total


def read_quantized(path):
    """Quantized depth runs per chromosome: {chrom: [(s, e, is_low)]}, sorted.

    mosdepth --quantize 0:T: writes LOW runs labelled "0:T" and HIGH runs "T:inf".
    Column 5 (repeat flag) is ignored: repeats are judged on RED coverage instead.
    """
    runs = {}
    for line in read_lines(path):
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        runs.setdefault(cols[0], []).append((int(cols[1]), int(cols[2]), cols[3].startswith("0:")))
    for rs in runs.values():
        rs.sort()
    return runs


def read_gaps(path):
    """Assembly gap intervals per chromosome, sorted: {chrom: [(s, e)]}."""
    gaps = {}
    for line in read_lines(path):
        cols = line.split("\t")
        if len(cols) >= 3:
            gaps.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))
    for ivs in gaps.values():
        ivs.sort()
    return gaps


def read_breakpoints(path):
    """Clustered DELFIES breakpoints per chromosome: {chrom: [(pos, strand, end)]}."""
    bps = {}
    for line in read_lines(path):
        cols = line.split("\t")
        if len(cols) < 4 or cols[0] == "chrom":
            continue
        bps.setdefault(cols[0], []).append((int(cols[1]), cols[3], int(cols[2])))
    return bps


def read_windows(path):
    """Mosdepth windows per chromosome: {chrom: [(s, e, depth)]}, sorted."""
    wins = {}
    for line in read_lines(path):
        cols = line.split()
        if len(cols) < 4 or line.startswith("#"):
            continue
        try:
            wins.setdefault(cols[0], []).append((int(cols[1]), int(cols[2]), float(cols[3])))
        except ValueError:
            continue
    for ws in wins.values():
        ws.sort()
    return wins


def median_depth(windows, chrom, s, e):
    """Median depth of mosdepth windows overlapping [s, e); 0.0 when none overlap."""
    vals = [v for ws, we, v in windows.get(chrom, []) if ws < e and we > s]
    return statistics.median(vals) if vals else 0.0


def subtract_gaps(lo, hi, gaps):
    """Split [lo, hi) into (start, end, is_gap) pieces around sorted gap intervals."""
    out = []
    cur = lo
    for gs, ge in gaps:
        if ge <= cur:
            continue
        if gs >= hi:
            break
        a, b = max(gs, cur), min(ge, hi)
        if a > cur:
            out.append((cur, a, False))
        out.append((a, b, True))
        cur = b
    if cur < hi:
        out.append((cur, hi, False))
    return out


def call_chromosome(chrom, length, runs, gaps, red, bps, policy):
    """Tile [0, length) with (start, end, label, subtype) pieces from quantized runs."""
    anchor = policy["anchor_proximity_bp"]
    min_low = policy["quantized_min_low_run_bp"]
    max_rep = policy["quantized_masked_repeat_fraction"]
    plus_starts = [pos for pos, strand, _end in bps if strand == "+"]
    minus_ends = [end for _pos, strand, end in bps if strand == "-"]

    pieces = []

    def add(s, e, label, subtype):
        if e > s:
            pieces.append((s, e, label, subtype))

    cursor = 0
    for rs, re_, is_low in runs:
        rs, re_ = max(rs, 0), min(re_, length)
        if re_ <= rs:
            continue
        if rs > cursor:
            add(cursor, rs, "putative_uncertain", "no_quantized_coverage")
        cursor = max(cursor, re_)
        if not is_low or (re_ - rs) < min_low:
            add(rs, re_, "putative_core", "-")
            continue
        for a, b, is_gap in subtract_gaps(rs, re_, gaps):
            if is_gap:
                add(a, b, "putative_uncertain", "gap_spanning")
                continue
            supported = any(abs(p - a) <= anchor for p in plus_starts) or any(
                abs(q - b) <= anchor for q in minus_ends
            )
            if supported:
                add(a, b, "putative_GRS", "quantized_delfies_anchored")
            elif covered_bp(a, b, red) / (b - a) > max_rep:
                add(a, b, "putative_uncertain", "repeat_masked")
            else:
                add(a, b, "putative_uncertain", "no_delfies_support")
    if cursor < length:
        add(cursor, length, "putative_uncertain", "no_quantized_coverage")
    return merge_pieces(pieces)


def merge_pieces(pieces):
    """Join adjacent pieces that share label and subtype."""
    out = []
    for s, e, label, subtype in pieces:
        if out and out[-1][1] == s and out[-1][2] == label and out[-1][3] == subtype:
            out[-1] = (out[-1][0], e, label, subtype)
        else:
            out.append((s, e, label, subtype))
    return out


def check_tiling(chrom, pieces, length):
    """Raise unless pieces are contiguous from 0 to length."""
    cur = 0
    for s, e, _label, _subtype in pieces:
        if s != cur or e <= s:
            raise ValueError(f"[Tiling] {chrom}: pieces do not tile 0-{length}")
        cur = e
    if cur != length:
        raise ValueError(f"[Tiling] {chrom}: pieces do not tile 0-{length}")


def read_summary(path):
    """Return (metadata lines, data rows) of an identify grs_summary.tsv."""
    meta, rows = [], []
    for line in read_lines(path):
        if line.startswith("#"):
            meta.append(line)
        elif line.startswith("chromosome\t"):
            continue
        elif line.strip():
            cols = line.split("\t")
            rows.append((cols[0], int(cols[1]), int(cols[2]), cols[4], float(cols[5]), cols[6]))
    return meta, rows


def parse_chromosomes(value):
    """Comma- or whitespace-separated chromosome list; empty when unset."""
    if not value:
        return []
    return [c for c in value.replace(",", " ").split() if c]


def build_segments(rows, quant, gaps, red, bps, windows, policy, chromosomes):
    """Ordered (chrom, start, end, label, median, subtype) rows for the final summary."""
    lengths = {}
    for chrom, _s, e, _lab, _med, _sub in rows:
        lengths[chrom] = max(lengths.get(chrom, 0), e)

    def is_canonical(chrom):
        if chromosomes:
            return chrom in chromosomes
        return lengths[chrom] >= policy["chromosomal_min_bp"]

    out = []
    done = set()
    for chrom, s, e, label, med, subtype in rows:
        if not is_canonical(chrom):
            out.append((chrom, s, e, label, med, subtype))
            continue
        if chrom in done:
            continue
        done.add(chrom)
        length = lengths[chrom]
        pieces = call_chromosome(
            chrom, length, quant.get(chrom, []), gaps.get(chrom, []),
            red.get(chrom), bps.get(chrom, []), policy,
        )
        check_tiling(chrom, pieces, length)
        for ps, pe, plabel, psub in pieces:
            out.append((chrom, ps, pe, plabel, median_depth(windows, chrom, ps, pe), psub))
    return out, done


def elimination_line(segments, canonical):
    """Summarise bp per class; canonical and unplaced are counted separately."""
    canon_grs = canon_unc = unpl_grs = unpl_unc = unpl_core = 0
    for chrom, s, e, label, _med, _sub in segments:
        length = e - s
        if chrom in canonical:
            if label == "putative_GRS":
                canon_grs += length
            elif label == "putative_uncertain":
                canon_unc += length
        elif label == "putative_GRS":
            unpl_grs += length
        elif label == "putative_uncertain":
            unpl_unc += length
        else:
            unpl_core += length
    min_elim = canon_grs + unpl_grs
    return (
        f"# ELIMINATION_RANGE\tcanon_grs_bp={canon_grs}\tunplaced_grs_bp={unpl_grs}"
        f"\tunplaced_uncertain_bp={unpl_unc}\tunplaced_core_bp={unpl_core}"
        f"\tmin_elim_bp={min_elim}\tmax_elim_bp={min_elim}\test_elim_bp={min_elim}"
        f"\tcanon_uncertain_bp={canon_unc}\n"
    )


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Final three-class segmentation of quantized depth runs"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--quantized-bed", required=True)
    parser.add_argument("--gaps-bed", default=None)
    parser.add_argument("--red-bed", default=None)
    parser.add_argument("--regions-bed", default=None)
    parser.add_argument("--breakpoints-tsv", default=None)
    parser.add_argument("--chromosomes", default=None)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--policy", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    policy = load_curation_policy(args.policy)["pde_segmentation"]
    quant = read_quantized(args.quantized_bed)
    gaps = read_gaps(args.gaps_bed)
    red = read_red(args.red_bed)
    bps = read_breakpoints(args.breakpoints_tsv)
    windows = read_windows(args.regions_bed)
    meta, rows = read_summary(args.summary_tsv)
    chromosomes = parse_chromosomes(args.chromosomes)

    segments, canonical = build_segments(rows, quant, gaps, red, bps, windows, policy, chromosomes)

    bed_path = os.path.join(args.outdir, f"{args.prefix}.pde_segments.bed")
    tsv_path = os.path.join(args.outdir, f"{args.prefix}.grs_summary.tsv")
    with open(bed_path, "w") as out_bed, open(tsv_path, "w") as out_tsv:
        for line in meta:
            if not line.startswith("# ELIMINATION_RANGE"):
                out_tsv.write(line + "\n")
        out_tsv.write(elimination_line(segments, canonical))
        out_tsv.write(SUMMARY_HEADER + "\n")
        for chrom, s, e, label, med, subtype in segments:
            out_bed.write(f"{chrom}\t{s}\t{e}\t{label}\n")
            out_tsv.write(f"{chrom}\t{s}\t{e}\t{e - s}\t{label}\t{med:.2f}\t{subtype}\n")


if __name__ == "__main__":
    main()
