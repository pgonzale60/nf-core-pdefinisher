#!/usr/bin/env python3
"""Build MEME break windows from DELFIES S2G foci and the genome classification.

Inputs
  --s2g-tsv         breakpoint_foci__S2G.tsv (contig, start,
                    num_supporting_reads__forward, num_supporting_reads__reverse).
                    Duplicate (contig, start) rows are summed.
  --pde-segments-bed  pde_segments.bed (contig, start, end, label). Core/GRS
                    borders are contiguous putative_core / putative_GRS pairs only
                    (uncertain never forms a border); the border sits at the
                    shared coordinate, as in pde_extract_borders.py.

Rules
  Site score      = max(fwd, rev).
  Support floor   = 0.5 * median(top 3 site scores) of the sample; sites below it
                    are dropped.
  Balanced site   = min(fwd, rev) >= 0.5 * max(fwd, rev), i.e. both strands
                    represented with disparity <= 50% of the max.
  Fission         = adjacent positions p, p+1, both balanced and above the floor.
                    Kept anywhere on the genome. Reported at p with kind "fission".
  Break site      = any other site above the floor within BORDER_TOL_BP of a
                    core/GRS border on the same contig. Kind "break".
  Clusters        = sorted kept positions less than MERGE_BP apart are merged.
                    Each cluster is centred on its midpoint and widened to
                    WINDOW_BP, clipped to the contig.

Outputs
  <prefix>.breakwindows.bed  BED6 + kind + n_sites (0-based half-open)
  <prefix>.breakwindows.fa   one record per BED line, header = BED name
"""

import argparse
import csv
import gzip
import json
import statistics
import sys
from collections import defaultdict

WINDOW_BP = 150
MERGE_BP = 15
BORDER_TOL_BP = 100
BALANCE = 0.5
TOP_N = 3


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def load_sites(s2g_tsv):
    """Return {(contig, pos): (fwd, rev)} with duplicate rows summed."""
    totals = defaultdict(lambda: [0, 0])
    with open_text(s2g_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                pos = int(row["start"])
                fwd = int(float(row.get("num_supporting_reads__forward") or 0))
                rev = int(float(row.get("num_supporting_reads__reverse") or 0))
            except (KeyError, ValueError):
                continue
            if pos < 0:
                continue
            cell = totals[(row["contig"], pos)]
            cell[0] += fwd
            cell[1] += rev
    return {key: (f, r) for key, (f, r) in totals.items()}


def load_segments(segments_bed):
    """Return {contig: sorted [(start, end, label)]}."""
    segs = defaultdict(list)
    with open_text(segments_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            segs[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
    for chrom in segs:
        segs[chrom].sort()
    return segs


def _is_grs(label):
    return label == "putative_GRS"


def _is_core(label):
    return label == "putative_core"


def find_borders(segs):
    """Return {contig: sorted border positions} where core and GRS meet."""
    borders = {}
    for chrom, seg_list in segs.items():
        pos_list = []
        for (_s1, e1, lab1), (s2, _e2, lab2) in zip(seg_list, seg_list[1:]):
            if e1 != s2:
                continue
            if (_is_grs(lab1) and _is_core(lab2)) or (_is_core(lab1) and _is_grs(lab2)):
                pos_list.append(e1)
        borders[chrom] = sorted(pos_list)
    return borders


def near_border(borders, contig, pos):
    """True if pos lies within BORDER_TOL_BP of a border on contig."""
    positions = borders.get(contig, [])
    lo = 0
    hi = len(positions)
    while lo < hi:
        mid = (lo + hi) // 2
        if positions[mid] < pos:
            lo = mid + 1
        else:
            hi = mid
    for idx in (lo - 1, lo):
        if 0 <= idx < len(positions) and abs(positions[idx] - pos) <= BORDER_TOL_BP:
            return True
    return False


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def site_score(fwd, rev):
    return max(fwd, rev)


def is_balanced(fwd, rev):
    hi = max(fwd, rev)
    return hi > 0 and min(fwd, rev) >= BALANCE * hi


def support_floor(sites):
    scores = sorted((site_score(f, r) for f, r in sites.values()), reverse=True)[:TOP_N]
    return 0.5 * statistics.median(scores) if scores else float("inf")


def call_sites(sites, borders):
    """Return ({contig: {pos: kind}}, kept_sites) for fissions and border break sites."""
    floor = support_floor(sites)
    kept = {key: fr for key, fr in sites.items() if site_score(*fr) >= floor}
    calls = defaultdict(dict)
    for (contig, pos), (fwd, rev) in kept.items():
        nxt = kept.get((contig, pos + 1))
        if is_balanced(fwd, rev) and nxt is not None and is_balanced(*nxt):
            calls[contig][pos] = "fission"
        elif near_border(borders, contig, pos):
            calls[contig][pos] = "break"
    return calls, kept


def cluster(positions):
    """Group sorted positions where consecutive gaps are < MERGE_BP."""
    groups = []
    for pos in sorted(positions):
        if groups and pos - groups[-1][-1] < MERGE_BP:
            groups[-1].append(pos)
        else:
            groups.append([pos])
    return groups


def window_for(group, contig_len):
    mid = (group[0] + group[-1]) // 2
    start = max(0, mid - WINDOW_BP // 2)
    end = start + WINDOW_BP
    if end > contig_len:
        end = contig_len
        start = max(0, end - WINDOW_BP)
    return start, end


def build(s2g_tsv, segments_bed, assembly, prefix):
    sites = load_sites(s2g_tsv)
    borders = find_borders(load_segments(segments_bed))
    calls, kept = call_sites(sites, borders)
    seqs = read_fasta(assembly)

    bed_lines = []
    fa_lines = []
    n_break = 0
    n_fission = 0
    idx = 0
    for contig in sorted(calls):
        if contig not in seqs:
            print(f"[WARN] contig {contig} not in assembly; skipped", file=sys.stderr)
            continue
        contig_len = len(seqs[contig])
        kinds = calls[contig]
        for group in cluster(kinds):
            idx += 1
            name = f"{prefix}_bw{idx}"
            start, end = window_for(group, contig_len)
            kind = "fission" if any(kinds[p] == "fission" for p in group) else "break"
            if kind == "fission":
                n_fission += 1
            else:
                n_break += 1
            score = min(1000, sum(site_score(*kept[(contig, p)]) for p in group))
            bed_lines.append(f"{contig}\t{start}\t{end}\t{name}\t{score}\t.\t{kind}\t{len(group)}")
            fa_lines.append(f">{name}\n{seqs[contig][start:end]}\n")

    with open(f"{prefix}.breakwindows.bed", "w", encoding="utf-8") as fh:
        fh.write("".join(line + "\n" for line in bed_lines))
    with open(f"{prefix}.breakwindows.fa", "w", encoding="utf-8") as fh:
        fh.write("".join(fa_lines))
    census = {
        "n_breaks": n_break + n_fission,
        "n_break_windows": n_break,
        "n_fission_windows": n_fission,
        "n_sites_kept": len(kept),
    }
    with open(f"{prefix}.census.json", "w", encoding="utf-8") as fh:
        json.dump(census, fh, indent=2)
        fh.write("\n")


def main():
    global WINDOW_BP, MERGE_BP, BORDER_TOL_BP
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--s2g-tsv", required=True)
    parser.add_argument("--pde-segments-bed", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--window-bp", type=int, default=WINDOW_BP)
    parser.add_argument("--merge-bp", type=int, default=MERGE_BP)
    parser.add_argument("--border-tol-bp", type=int, default=BORDER_TOL_BP)
    args = parser.parse_args()
    WINDOW_BP = args.window_bp
    MERGE_BP = args.merge_bp
    BORDER_TOL_BP = args.border_tol_bp
    build(args.s2g_tsv, args.pde_segments_bed, args.assembly, args.prefix)


if __name__ == "__main__":
    main()
