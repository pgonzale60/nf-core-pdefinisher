#!/usr/bin/env python3
"""Carve high-confidence core islands out of long putative_GRS intervals.

A long GRS interval is re-segmented with a two-state (LOW/HIGH) Viterbi model over
quantized coverage runs. Repeat-dominated runs and assembly gaps carry no evidence:
they take the state of the preceding evidence (carry-forward). A maximal HIGH stretch
that touches an assembly gap is a core island and is re-labelled putative_core.
"""

import argparse
import bisect
import json
import os
import sys
from pathlib import Path

__version__ = "1.0.0"

LOW = 0
HIGH = 1
EXCLUDED_SUBTYPES = {"fully_repetitive_unplaced", "unplaced_coverage_drop"}
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
    with open(path) as f:
        return f.read().splitlines()


def read_red(path):
    """RED repeat intervals per chromosome, merged: {chrom: (starts, intervals)}."""
    repeats = {}
    for line in read_lines(path):
        cols = line.split()
        if len(cols) >= 3 and not line.startswith("#"):
            try:
                repeats.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))
            except ValueError:
                continue
    out = {}
    for chrom, ivs in repeats.items():
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        merged = [(s, e) for s, e in merged]
        out[chrom] = ([s for s, _ in merged], merged)
    return out


def covered_bp(s, e, red):
    """Bases of [s, e) covered by merged, sorted RED intervals (starts, intervals)."""
    if red is None:
        return 0
    starts, ivs = red
    i = max(bisect.bisect_right(starts, s) - 1, 0)
    total = 0
    for a, b in ivs[i:]:
        if a >= e:
            break
        lo, hi = max(a, s), min(b, e)
        if hi > lo:
            total += hi - lo
    return total


def read_quantized(path):
    """Quantized depth runs per chromosome: {chrom: (starts, [(s, e, is_low)])}.

    mosdepth --quantize 0:T: writes LOW runs labelled "0:T" and HIGH runs "T:inf".
    Column 5 (repeat flag) is ignored: repeats are judged on RED coverage instead.
    """
    runs = {}
    for line in read_lines(path):
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        runs.setdefault(cols[0], []).append(
            (int(cols[1]), int(cols[2]), cols[3].startswith("0:"))
        )
    out = {}
    for chrom, rs in runs.items():
        rs.sort()
        out[chrom] = ([r[0] for r in rs], rs)
    return out


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


def read_baseline(path):
    """Somatic core baseline depth: column 2 of row 2 in core_coverage.tsv."""
    lines = read_lines(path)
    if len(lines) < 2:
        raise ValueError(f"core coverage table has no data row: {path}")
    return float(lines[1].split("\t")[1])


def subtract_gaps(lo, hi, gaps):
    """Split [lo, hi) around sorted gap intervals; returns the non-gap sub-intervals."""
    out = []
    cur = lo
    for gs, ge in gaps:
        if ge <= cur:
            continue
        if gs >= hi:
            break
        if gs > cur:
            out.append((cur, gs))
        cur = max(cur, ge)
        if cur >= hi:
            break
    if cur < hi:
        out.append((cur, hi))
    return out


def build_observations(s, e, quant, gaps, red, repeat_fraction):
    """Ordered, non-overlapping observations over [s, e).

    Each is (start, end, kind, is_low) with kind in {info, mask, gap}. Repeat-dominated
    runs are 'mask' (weight 0); gaps are 'gap' (weight 0) and take precedence.
    """
    starts, runs = quant
    gaps_in = [(max(gs, s), min(ge, e)) for gs, ge in gaps if ge > s and gs < e]
    obs = []
    i = max(bisect.bisect_right(starts, s) - 1, 0)
    for a, b, is_low in runs[i:]:
        if a >= e:
            break
        lo, hi = max(a, s), min(b, e)
        if hi <= lo:
            continue
        kind = "mask" if covered_bp(lo, hi, red) / (hi - lo) > repeat_fraction else "info"
        for slo, shi in subtract_gaps(lo, hi, gaps_in):
            obs.append((slo, shi, kind, is_low))
    obs.extend((gs, ge, "gap", None) for gs, ge in gaps_in)
    obs.sort(key=lambda o: o[0])
    return obs


def viterbi(obs, switch_penalty):
    """Two-state path over observations; returns a state per observation.

    Informative runs cost their length when the state disagrees with them. A state change
    costs switch_penalty and is charged only when entering an informative observation:
    masked and gap observations carry the previous state forward at zero cost and can
    never initiate a switch (carry-forward). Ties keep the previous state. Ties at the
    end resolve to LOW, so unsupported stretches never become core.
    """
    if not obs:
        return []
    values = [0.0, 0.0]
    back = []
    for idx, (lo, hi, kind, is_low) in enumerate(obs):
        emit = [0.0, 0.0]
        if kind == "info":
            length = hi - lo
            emit[LOW] = 0.0 if is_low else float(length)
            emit[HIGH] = float(length) if is_low else 0.0
        ptr = [LOW, HIGH]
        new = [0.0, 0.0]
        for state in (LOW, HIGH):
            if idx == 0:
                best = 0.0
                ptr[state] = state
            elif kind != "info":
                best, ptr[state] = values[state], state
            else:
                stay = values[state]
                switch = values[1 - state] + switch_penalty
                if stay <= switch:
                    best, ptr[state] = stay, state
                else:
                    best, ptr[state] = switch, 1 - state
            new[state] = best + emit[state]
        values = new
        back.append(ptr)
    state = HIGH if values[HIGH] < values[LOW] else LOW
    path = [0] * len(obs)
    for idx in range(len(obs) - 1, -1, -1):
        path[idx] = state
        state = back[idx][state]
    return path


def find_islands(obs, path, gaps_clip):
    """Maximal HIGH stretches [a, b) whose span touches an assembly gap."""
    islands = []
    i = 0
    while i < len(obs):
        if path[i] != HIGH:
            i += 1
            continue
        j = i
        while j + 1 < len(obs) and path[j + 1] == HIGH:
            j += 1
        a, b = obs[i][0], obs[j][1]
        has_evidence = any(obs[k][2] == "info" and not obs[k][3] for k in range(i, j + 1))
        if has_evidence and any(gs <= b and ge >= a for gs, ge in gaps_clip):
            islands.append((a, b))
        i = j + 1
    return islands


def carve_row(row, quant, gaps, red, baseline, policy):
    """Return (pieces, n_islands, carved_bp) for one summary row.

    row = (chrom, start, end, classification, median_depth, subtype).
    """
    chrom, s, e, label, med, subtype = row
    if (
        label != "putative_GRS"
        or e - s <= policy["core_subtraction_min_grs_bp"]
        or subtype in EXCLUDED_SUBTYPES
    ):
        return [row], 0, 0
    gaps_chrom = gaps.get(chrom, [])
    gaps_clip = [(gs, ge) for gs, ge in gaps_chrom if ge >= s and gs <= e]
    if not gaps_clip or chrom not in quant:
        return [row], 0, 0
    obs = build_observations(
        s, e, quant[chrom], gaps_chrom, red.get(chrom),
        policy["core_subtraction_masked_repeat_fraction"],
    )
    path = viterbi(obs, policy["core_subtraction_switch_penalty_bp"])
    islands = find_islands(obs, path, gaps_clip)
    if not islands:
        return [row], 0, 0
    pieces = []
    cur = s
    for a, b in islands:
        if a > cur:
            pieces.append((chrom, cur, a, label, med, subtype))
        pieces.append((chrom, a, b, "putative_core", baseline, "core_island_gap"))
        cur = b
    if e > cur:
        pieces.append((chrom, cur, e, label, med, subtype))
    carved = sum(b - a for a, b in islands)
    return pieces, len(islands), carved


def read_summary(path):
    """Return (metadata lines, data rows) of a grs_summary.tsv."""
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


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Carve core islands out of long putative_GRS intervals"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--quantized-bed", required=True)
    parser.add_argument("--gaps-bed", required=True)
    parser.add_argument("--red-bed", default=None)
    parser.add_argument("--core-coverage-tsv", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--policy", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    policy = load_curation_policy(args.policy)["pde_segmentation"]
    baseline = read_baseline(args.core_coverage_tsv)
    quant = read_quantized(args.quantized_bed)
    gaps = read_gaps(args.gaps_bed)
    red = read_red(args.red_bed)
    meta, rows = read_summary(args.summary_tsv)

    out_rows = []
    n_islands = 0
    carved_bp = 0
    for row in rows:
        pieces, n, bp = carve_row(row, quant, gaps, red, baseline, policy)
        out_rows.extend(pieces)
        n_islands += n
        carved_bp += bp

    bed_path = os.path.join(args.outdir, f"{args.prefix}.pde_segments.bed")
    tsv_path = os.path.join(args.outdir, f"{args.prefix}.grs_summary.tsv")
    with open(bed_path, "w") as out_bed, open(tsv_path, "w") as out_tsv:
        for line in meta:
            out_tsv.write(line + "\n")
        out_tsv.write(f"# CORE_ISLANDS\tn_islands={n_islands}\tcarved_bp={carved_bp}\n")
        out_tsv.write(SUMMARY_HEADER + "\n")
        for chrom, s, e, label, med, subtype in out_rows:
            out_bed.write(f"{chrom}\t{s}\t{e}\t{label}\n")
            out_tsv.write(f"{chrom}\t{s}\t{e}\t{e - s}\t{label}\t{med:.2f}\t{subtype}\n")


if __name__ == "__main__":
    main()
