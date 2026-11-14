#!/usr/bin/env python3
"""Segmentation of a sequence into putative GRS, core, variable and uncertain segments.

Implements R1-R11 of PDE_origin docs/segmentation_definition.md:
break sites are classified from DELFIES breakpoints (echo/weak/site/two_strand);
two-strand pairs within site_merge_bp are fissions (sidecar output, cuts only);
span state (inside/outside/undetermined) is read from the nearest break sites;
depth is measured on 500 bp mosdepth windows relative to a core baseline (E/I/R);
the label is set by a depth x span x border table; unplaced sequences use depth only;
rDNA 45S arrays (BARRNAP) mask depth and mark all-array sequences. Every base of every
sequence gets exactly one label.
"""

import argparse
import bisect
import gzip
import json
import os
import statistics
import sys
from collections import namedtuple
from pathlib import Path

__version__ = "2.0.0"

SUMMARY_HEADER = (
    "chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype\trepeat_bp\trdna_bp"
)
FISSION_HEADER = "chromosome\tplus_pos\tminus_pos\tdistance_bp\tcut_pos"

GRS = "putative_GRS"
CORE = "putative_core"
VARIABLE = "putative_variable"
UNCERTAIN = "uncertain"
UNPLACED_ELIM = "unplaced_eliminated"
UNPLACED_RDNA = "unplaced_rdna"

Site = namedtuple("Site", "pos end strand reads")


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
    """Bases of [s, e) covered by merged, sorted intervals."""
    total = 0
    for a, b in red_ivs or []:
        if a >= e:
            break
        lo, hi = max(a, s), min(b, e)
        if hi > lo:
            total += hi - lo
    return total


class Intervals:
    """Merged, sorted intervals with O(log n) coverage queries."""

    def __init__(self, ivs):
        self.ivs = merge_intervals(ivs)
        self.starts = [s for s, _ in self.ivs]
        self.ends = [e for _, e in self.ivs]

    def covered(self, s, e):
        i = bisect.bisect_right(self.ends, s)
        j = bisect.bisect_left(self.starts, e)
        return covered_bp(s, e, self.ivs[i:j])


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


def read_site_rows(path):
    """Clustered DELFIES breakpoints: {chrom: [(pos, end, strand, support)]}."""
    rows = {}
    for line in read_lines(path):
        cols = line.split("\t")
        if len(cols) < 5 or cols[0] == "chrom":
            continue
        rows.setdefault(cols[0], []).append((int(cols[1]), int(cols[2]), cols[3], float(cols[4])))
    return rows


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


class WindowSet:
    """Mosdepth windows with midpoint-based interval lookup."""

    def __init__(self, wins):
        self.wins = wins
        self.starts = {c: [w[0] for w in ws] for c, ws in wins.items()}
        self.maxlen = max((e - s for ws in wins.values() for s, e, _ in ws), default=0)

    def in_interval(self, chrom, s, e):
        """Windows whose midpoint lies in [s, e)."""
        ws = self.wins.get(chrom, [])
        if not ws or e <= s:
            return []
        i = bisect.bisect_left(self.starts[chrom], s - self.maxlen)
        out = []
        while i < len(ws) and ws[i][0] < e:
            a, b, _ = ws[i]
            if 2 * s <= a + b < 2 * e:
                out.append(ws[i])
            i += 1
        return out

    def usable(self, chrom, s, e, mask):
        """Windows with midpoint in [s, e) that do not overlap the rDNA mask."""
        return [w for w in self.in_interval(chrom, s, e) if mask.covered(w[0], w[1]) == 0]


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


def format_subranges(subranges):
    """Render subtype subranges as '<subtype>:<bp>,<subtype>:<bp>'."""
    return ",".join(f"{subtype}:{bp}" for subtype, bp in subranges)


def check_tiling(chrom, pieces, length):
    """Raise unless pieces are contiguous from 0 to length."""
    cur = 0
    for s, e, _label, _subtype in pieces:
        if s != cur or e <= s:
            raise ValueError(f"[Tiling] {chrom}: pieces do not tile 0-{length}")
        cur = e
    if cur != length:
        raise ValueError(f"[Tiling] {chrom}: pieces do not tile 0-{length}")


# --- R2-R4: break sites, fissions, span state -------------------------------------------------

def read_sites(rows, policy):
    """Merged break sites per chromosome: {chrom: [Site]} sorted by pos (R2)."""
    merge_bp = policy["site_merge_bp"]
    min_reads = policy["site_min_reads"]
    sites = {}
    for chrom, items in rows.items():
        out = []
        for strand in ("+", "-"):
            group = sorted((r for r in items if r[2] == strand), key=lambda r: r[0])
            clusters = []
            for pos, end, _, support in group:
                if clusters and pos - clusters[-1][-1][0] <= merge_bp:
                    clusters[-1].append((pos, end, support))
                else:
                    clusters.append([(pos, end, support)])
            for members in clusters:
                total = sum(m[2] for m in members)
                best = max(members, key=lambda m: m[2])
                if total >= min_reads:
                    out.append(Site(best[0], best[1], strand, total))
        sites[chrom] = sorted(out)
    return sites


def top_support(sites, n):
    """Median of the n largest merged site supports across all sequences."""
    vals = sorted((s.reads for ss in sites.values() for s in ss), reverse=True)[:n]
    return statistics.median(vals) if vals else 0.0


def site_kind(reads, partner, floor, balance):
    """R3 classification of one merged site against its opposite-strand partner."""
    if reads < balance * partner:
        return "echo"
    if reads < floor:
        return "weak"
    if partner >= max(floor, balance * reads):
        return "two_strand"
    return "site"


def classify_sites(sites, policy):
    """{chrom: [(Site, kind)]} with R3 kinds."""
    merge_bp = policy["site_merge_bp"]
    balance = policy["site_balance"]
    floor = policy["site_floor_fraction"] * top_support(sites, policy["site_top_n"])
    out = {}
    for chrom, ss in sites.items():
        kinded = []
        for s in ss:
            partner = max(
                (o.reads for o in ss if o.strand != s.strand and abs(o.pos - s.pos) <= merge_bp),
                default=0,
            )
            kinded.append((s, site_kind(s.reads, partner, floor, balance)))
        out[chrom] = kinded
    return out


def cut_coord(site):
    """Cut coordinate of a non-fission site: pos for '+', end for '-'."""
    return site.pos if site.strand == "+" else site.end


def find_fissions(kinded, merge_bp):
    """Pair consecutive opposite-strand two_strand sites into fissions (R4).

    Returns (fissions, bounds): fissions are (plus_pos, minus_pos, distance_bp, cut_pos);
    bounds are (coord, kind) with kind '+', '-' for sites and 'F' for fissions.
    """
    two = [s for s, k in kinded if k == "two_strand"]
    fissions = []
    bounds = [(cut_coord(s), s.strand) for s, k in kinded if k == "site"]
    i = 0
    while i < len(two):
        a = two[i]
        b = two[i + 1] if i + 1 < len(two) else None
        if b is not None and a.strand != b.strand and b.pos - a.pos <= merge_bp:
            plus, minus = (a, b) if a.strand == "+" else (b, a)
            cut = plus.pos + (1 if minus.pos >= plus.pos else -1)
            fissions.append((plus.pos, minus.pos, minus.pos - plus.pos, cut))
            bounds.append((cut, "F"))
            i += 2
        else:
            bounds.append((cut_coord(a), a.strand))
            i += 1
    return fissions, sorted(bounds)


def span_state(bounds, x):
    """Span state at x from sorted (coord, kind) bounds, as collect_evidence.span_state."""
    if not bounds:
        return "no_break_site"
    coords = [c for c, _ in bounds]
    i = bisect.bisect_right(coords, x)
    left = None if i == 0 else ("grs" if bounds[i - 1][1] == "+" else "core")
    right = None if i == len(bounds) else ("grs" if bounds[i][1] == "-" else "core")
    if left is None or right is None or left == right:
        return "inside_span" if (left or right) == "grs" else "outside_span"
    return "undetermined"


# --- R5/R6: pieces and depth -------------------------------------------------------------------

def make_units(length, coords, gap_ivs):
    """Gap-free pieces and gap units tiling [0, length) with their border kinds.

    Each piece dict: gap False, start, end, lb, rb in {site, gap, end}, left_cut.
    Cuts strictly inside a gap never occur because pieces are split only inside
    non-gap segments.
    """
    units = []
    for s, e, is_gap in subtract_gaps(0, length, gap_ivs.ivs):
        if is_gap:
            units.append({"gap": True, "start": s, "end": e})
            continue
        cuts = sorted({c for c in coords if s < c < e})
        edges = [s, *cuts, e]
        for a, b in zip(edges, edges[1:]):
            units.append({"gap": False, "start": a, "end": b})
    for i, u in enumerate(units):
        if u["gap"]:
            continue
        for side, j in (("lb", i - 1), ("rb", i + 1)):
            if j < 0 or j >= len(units):
                u[side] = "end"
            elif units[j]["gap"]:
                u[side] = "gap"
            else:
                u[side] = "site"
        u["left_cut"] = u["lb"] == "site"
    return units


class Seq:
    """Per-sequence evidence used by the labelling rules."""

    def __init__(self, chrom, length, canonical, bounds, gaps, red, mask, annotated, windows, policy):
        self.chrom = chrom
        self.length = length
        self.canonical = canonical
        self.bounds = bounds
        self.gaps = gaps
        self.red = red
        self.mask = mask
        self.annotated = annotated
        self.windows = windows
        self.policy = policy
        self.baseline = None


def state_of(seq, x):
    """Span state at x; no break sites read as outside_span (R2/R7)."""
    st = span_state(seq.bounds, x)
    return "outside_span" if st == "no_break_site" else st


def window_mode(rel, policy):
    if rel < policy["depth_cut_eliminated"]:
        return "E"
    if rel < policy["depth_cut_retained"]:
        return "I"
    return "R"


def measure(seq, s, e):
    """(mode, median) of [s, e) against the sequence baseline; (None, 0.0) if unmeasurable."""
    if e - s < seq.policy["min_depth_piece_bp"] or not seq.baseline or seq.baseline <= 0:
        return None, 0.0
    vals = [v for _, _, v in seq.windows.usable(seq.chrom, s, e, seq.mask)]
    if not vals:
        return None, 0.0
    med = statistics.median(vals)
    return window_mode(med / seq.baseline, seq.policy), med


def mode_runs(seq, s, e):
    """Maximal contiguous runs [mode, start, end] over usable windows of [s, e)."""
    runs = []
    for ws, we, v in seq.windows.usable(seq.chrom, s, e, seq.mask):
        m = window_mode(v / seq.baseline, seq.policy)
        if runs and runs[-1][0] == m and runs[-1][2] == ws:
            runs[-1][2] = we
        else:
            runs.append([m, ws, we])
    return runs


def depth_boundary_cuts(seq, s, e):
    """R8: cut coordinates where a sustained E run meets a sustained R run."""
    if not seq.baseline:
        return []
    min_run = seq.policy["depth_boundary_min_run_bp"]
    runs = mode_runs(seq, s, e)
    cuts = []
    for (m1, a1, b1), (m2, a2, b2) in zip(runs, runs[1:]):
        if {m1, m2} == {"E", "R"} and b1 - a1 >= min_run and b2 - a2 >= min_run and b1 == a2:
            cuts.append(b1)
    return cuts


# --- R7 table, R8 and R10 labelling ------------------------------------------------------------

def table_label(state, lb, rb, mode):
    """R7 table: (label, subtype) for a gap-free piece with span state, borders and mode."""
    gap = lb == "gap" or rb == "gap"
    gaps_both = lb == "gap" and rb == "gap"
    one_end = gap and not gaps_both
    if mode is None:
        if gap:
            return UNCERTAIN, "no_usable_depth"
        if state == "inside_span":
            return GRS, "short_piece_in_span"
        if state == "outside_span":
            return CORE, "short_piece_outside_span"
        return UNCERTAIN, "no_usable_depth"
    elim_gap = "eliminated_depth_one_end_only" if one_end else "eliminated_depth_no_break_site"
    if state == "inside_span":
        if not gap:
            return {
                "E": (GRS, "span_eliminated_depth"),
                "I": (GRS, "span_intermediate_depth"),
                "R": (UNCERTAIN, "retained_depth_inside_span"),
            }[mode]
        return {
            "E": (GRS, "span_eliminated_depth"),
            "I": (UNCERTAIN, "intermediate_depth"),
            "R": (CORE, "span_misjoin_suspected") if gaps_both
            else (UNCERTAIN, "retained_depth_inside_span"),
        }[mode]
    if state == "outside_span":
        if not gap:
            return {
                "E": (UNCERTAIN, "eliminated_depth_no_break_site"),
                "I": (CORE, "retained_between_sites"),
                "R": (CORE, "retained_depth"),
            }[mode]
        return {
            "E": (UNCERTAIN, elim_gap),
            "I": (UNCERTAIN, "intermediate_depth"),
            "R": (CORE, "retained_depth"),
        }[mode]
    # undetermined
    if not gap:
        return {
            "E": (GRS, "undetermined_eliminated_depth"),
            "I": (VARIABLE, "undetermined_intermediate_depth"),
            "R": (CORE, "undetermined_retained_depth"),
        }[mode]
    return {
        "E": (UNCERTAIN, elim_gap),
        "I": (UNCERTAIN, "intermediate_depth"),
        "R": (CORE, "retained_depth"),
    }[mode]


def labelled(unit, label, subtype):
    out = dict(unit)
    out["label"] = label
    out["sub"] = subtype
    return out


def table_piece(seq, unit):
    state = state_of(seq, unit["start"])
    mode, _ = measure(seq, unit["start"], unit["end"])
    return labelled(unit, *table_label(state, unit["lb"], unit["rb"], mode))


def depth_only_label(mode):
    """Label of an R8 sub-piece from its own median mode."""
    if mode is None:
        return UNCERTAIN, "no_usable_depth"
    return {
        "R": (CORE, "depth_only_boundary"),
        "E": (UNCERTAIN, "depth_only_low_side"),
        "I": (UNCERTAIN, "intermediate_depth"),
    }[mode]


def split_piece(seq, unit, cuts):
    """R8: split an undetermined piece at qualifying cuts; relabel each sub-piece by depth."""
    edges = [unit["start"], *cuts, unit["end"]]
    out = []
    for i, (a, b) in enumerate(zip(edges, edges[1:])):
        mode, _ = measure(seq, a, b)
        sub = {"gap": False, "start": a, "end": b,
               "left_cut": unit["left_cut"] if i == 0 else True}
        out.append(labelled(sub, *depth_only_label(mode)))
    return out


def label_unplaced(seq, units):
    """R10: non-canonical sequence without break sites."""
    total = seq.length
    if seq.mask.covered(0, total) / total >= seq.policy["rdna_all_array_fraction"]:
        whole = {"gap": False, "start": 0, "end": total, "left_cut": False, "lb": "end", "rb": "end"}
        return [labelled(whole, UNPLACED_RDNA, "all_array")]
    out = []
    for u in units:
        if u["gap"]:
            out.append(u)
            continue
        mode, _ = measure(seq, u["start"], u["end"])
        if mode is None:
            out.append(labelled(u, UNCERTAIN, "no_usable_depth"))
        elif mode == "E":
            out.append(labelled(u, UNPLACED_ELIM, "eliminated_depth"))
        elif mode == "I":
            out.append(labelled(u, UNCERTAIN, "intermediate_depth"))
        else:
            out.append(labelled(u, CORE, "retained_depth"))
    return out


def label_sequence(seq, units):
    """Labelled units (pieces and gap units) tiling the sequence."""
    if not seq.canonical and not seq.bounds:
        return label_unplaced(seq, units)
    out = []
    for u in units:
        if u["gap"]:
            out.append(u)
            continue
        if state_of(seq, u["start"]) == "undetermined":
            cuts = depth_boundary_cuts(seq, u["start"], u["end"])
            if cuts:
                out.extend(split_piece(seq, u, cuts))
                continue
        out.append(table_piece(seq, u))
    return out


# --- R5: segments -------------------------------------------------------------------------------

def add_run(runs, subtype, bp):
    if runs and runs[-1][0] == subtype:
        runs[-1][1] += bp
    else:
        runs.append([subtype, bp])


def assemble(items):
    """Attach gap units to neighbours, then merge same-label pieces without a cut between them."""
    pieces = []
    pending_start = None
    for it in items:
        if it["gap"]:
            if pieces:
                pieces[-1]["end"] = it["end"]
            else:
                pending_start = it["start"]
            continue
        if pending_start is not None:
            it["start"] = pending_start
            pending_start = None
        pieces.append(it)
    segs = []
    for pc in pieces:
        if segs and segs[-1]["label"] == pc["label"] and not pc["left_cut"]:
            seg = segs[-1]
            seg["end"] = pc["end"]
            add_run(seg["runs"], pc["sub"], pc["end"] - pc["start"])
        else:
            segs.append({"start": pc["start"], "end": pc["end"], "label": pc["label"],
                         "runs": [[pc["sub"], pc["end"] - pc["start"]]]})
    return segs


def median_of(seq, s, e):
    vals = [v for _, _, v in seq.windows.usable(seq.chrom, s, e, seq.mask)]
    return statistics.median(vals) if vals else 0.0


def core_depths(seq):
    """Depths of windows that are outside spans, unmasked, gap-free and not repeat-dominated."""
    frac = seq.policy["core_depth_max_repeat_fraction"]
    vals = []
    for ws, we, v in seq.windows.wins.get(seq.chrom, []):
        if state_of(seq, ws) != "outside_span":
            continue
        if seq.gaps.covered(ws, we) or seq.mask.covered(ws, we):
            continue
        if seq.red.covered(ws, we) / (we - ws) > frac:
            continue
        vals.append(v)
    return vals


# --- R11: rDNA 45S ------------------------------------------------------------------------------

def read_45s_hits(path):
    """45S rRNA hits per chromosome: {chrom: [(start, end, strand)]}, 0-based half-open."""
    hits = {}
    for line in read_lines(path):
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9 or cols[2] != "rRNA":
            continue
        attr = cols[8]
        if any(x in attr for x in ("16S", "23S")):
            continue
        if not any(x in attr for x in ("18S", "5.8S", "5_8S", "28S")):
            continue
        hits.setdefault(cols[0], []).append((int(cols[3]) - 1, int(cols[4]), cols[6]))
    return hits


def group_regions(hits, max_sep):
    """Same-strand 45S hits within max_sep of each other form regions: [(start, end)]."""
    regions = []
    for strand in ("+", "-"):
        cur = None
        for s, e, st in sorted(h for h in hits if h[2] == strand):
            if cur and s - cur[1] <= max_sep:
                cur[1] = max(cur[1], e)
            else:
                if cur:
                    regions.append((cur[0], cur[1]))
                cur = [s, e]
        if cur:
            regions.append((cur[0], cur[1]))
    return regions


def group_loci(regions, merge_bp):
    """Regions within merge_bp form loci: [{start, end, size}] where size is the 45S bp."""
    loci = []
    for s, e in sorted(regions):
        if loci and s - loci[-1]["end"] <= merge_bp:
            loci[-1]["end"] = max(loci[-1]["end"], e)
            loci[-1]["size"] += e - s
        else:
            loci.append({"start": s, "end": e, "size": e - s})
    return loci


def read_45s_loci(path, policy):
    """{chrom: [locus]} from BARRNAP GFF3."""
    hits = read_45s_hits(path)
    return {
        chrom: group_loci(group_regions(hs, policy["rdna_45s_max_separation_bp"]),
                          policy["rdna_locus_merge_bp"])
        for chrom, hs in hits.items()
    }


def rdna_masks(loci, lengths, canonical, policy):
    """mask and annotated spans per chromosome (R11).

    Canonical: depth mask = all loci unless the main (largest) array is internal.
    Non-canonical: always masked. Annotated: all loci regardless of masking.
    """
    mask, annotated = {}, {}
    for chrom, ls in loci.items():
        spans = [(l["start"], l["end"]) for l in ls]
        annotated[chrom] = spans
        if chrom not in canonical:
            mask[chrom] = spans
            continue
        main = max(ls, key=lambda l: l["size"])
        internal = min(main["start"], lengths[chrom] - main["end"]) > policy["rdna_internal_min_dist_bp"]
        mask[chrom] = [] if internal else spans
    return mask, annotated


# --- output -------------------------------------------------------------------------------------

def elimination_line(segments, canonical):
    """bp per label, counted explicitly; canonical and unplaced are separate."""
    keys = [
        "canon_grs_bp", "canon_uncertain_bp", "canon_variable_bp",
        "unplaced_grs_bp", "unplaced_uncertain_bp", "unplaced_core_bp", "unplaced_variable_bp",
        "unplaced_eliminated_bp", "unplaced_rdna_bp",
    ]
    tot = dict.fromkeys(keys, 0)
    canon_map = {GRS: "canon_grs_bp", UNCERTAIN: "canon_uncertain_bp", VARIABLE: "canon_variable_bp"}
    unpl_map = {
        GRS: "unplaced_grs_bp", UNCERTAIN: "unplaced_uncertain_bp", CORE: "unplaced_core_bp",
        VARIABLE: "unplaced_variable_bp", UNPLACED_ELIM: "unplaced_eliminated_bp",
        UNPLACED_RDNA: "unplaced_rdna_bp",
    }
    for seg in segments:
        chrom, s, e, label = seg[0], seg[1], seg[2], seg[3]
        key = (canon_map if chrom in canonical else unpl_map).get(label)
        if key:
            tot[key] += e - s
    min_elim = tot["canon_grs_bp"] + tot["unplaced_grs_bp"]
    return (
        "# ELIMINATION_RANGE"
        f"\tcanon_grs_bp={tot['canon_grs_bp']}\tunplaced_grs_bp={tot['unplaced_grs_bp']}"
        f"\tunplaced_uncertain_bp={tot['unplaced_uncertain_bp']}"
        f"\tunplaced_core_bp={tot['unplaced_core_bp']}"
        f"\tmin_elim_bp={min_elim}\tmax_elim_bp={min_elim}\test_elim_bp={min_elim}"
        f"\tcanon_uncertain_bp={tot['canon_uncertain_bp']}"
        f"\tcanon_variable_bp={tot['canon_variable_bp']}"
        f"\tunplaced_variable_bp={tot['unplaced_variable_bp']}"
        f"\tunplaced_eliminated_bp={tot['unplaced_eliminated_bp']}"
        f"\tunplaced_rdna_bp={tot['unplaced_rdna_bp']}\n"
    )


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Segmentation of sequences into putative GRS, core, variable and uncertain"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--barrnap-gff", required=True)
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
    gaps_raw = read_gaps(args.gaps_bed)
    red_raw = read_red(args.red_bed)
    windows = WindowSet(read_windows(args.regions_bed))
    meta, rows = read_summary(args.summary_tsv)
    chromosome_list = parse_chromosomes(args.chromosomes)

    order = list(dict.fromkeys(r[0] for r in rows))
    lengths = {}
    for chrom, _s, e, *_ in rows:
        lengths[chrom] = max(lengths.get(chrom, 0), e)

    def is_canonical(chrom):
        if chromosome_list:
            return chrom in chromosome_list
        return lengths[chrom] >= policy["chromosomal_min_bp"]

    canonical = {c for c in order if is_canonical(c)}

    kinded = classify_sites(read_sites(read_site_rows(args.breakpoints_tsv), policy), policy)
    fissions, bounds = {}, {}
    for chrom in order:
        fissions[chrom], bounds[chrom] = find_fissions(kinded.get(chrom, []), policy["site_merge_bp"])

    loci = read_45s_loci(args.barrnap_gff, policy)
    mask_raw, annotated_raw = rdna_masks(loci, lengths, canonical, policy)

    seqs = {}
    for chrom in order:
        seqs[chrom] = Seq(
            chrom, lengths[chrom], chrom in canonical, bounds[chrom],
            Intervals(gaps_raw.get(chrom, [])), Intervals(red_raw.get(chrom, [])),
            Intervals(mask_raw.get(chrom, [])), Intervals(annotated_raw.get(chrom, [])),
            windows, policy,
        )

    core = {chrom: core_depths(seqs[chrom]) for chrom in order if chrom in canonical}
    pooled = [v for vals in core.values() for v in vals]
    taxon_baseline = statistics.median(pooled) if pooled else None
    for chrom in order:
        seqs[chrom].baseline = statistics.median(core[chrom]) if core.get(chrom) else taxon_baseline

    fission_path = os.path.join(args.outdir, f"{args.prefix}.fissions.tsv")
    bed_path = os.path.join(args.outdir, f"{args.prefix}.pde_segments.bed")
    tsv_path = os.path.join(args.outdir, f"{args.prefix}.grs_summary.tsv")
    with open(fission_path, "w") as out_fis:
        out_fis.write(FISSION_HEADER + "\n")
        for chrom in order:
            for plus, minus, dist, cut in fissions[chrom]:
                out_fis.write(f"{chrom}\t{plus}\t{minus}\t{dist}\t{cut}\n")

    segments = []
    for chrom in order:
        seq = seqs[chrom]
        units = make_units(seq.length, [c for c, _ in seq.bounds], seq.gaps)
        segs = assemble(label_sequence(seq, units))
        check_tiling(chrom, [(g["start"], g["end"], g["label"], None) for g in segs], seq.length)
        for g in segs:
            s, e = g["start"], g["end"]
            segments.append((
                chrom, s, e, g["label"], median_of(seq, s, e),
                format_subranges(g["runs"]), seq.red.covered(s, e), seq.annotated.covered(s, e),
            ))

    with open(bed_path, "w") as out_bed, open(tsv_path, "w") as out_tsv:
        for line in meta:
            if not line.startswith("# ELIMINATION_RANGE"):
                out_tsv.write(line + "\n")
        out_tsv.write(elimination_line(segments, canonical))
        out_tsv.write(SUMMARY_HEADER + "\n")
        for chrom, s, e, label, med, subtype, rep, rdna in segments:
            out_bed.write(f"{chrom}\t{s}\t{e}\t{label}\n")
            out_tsv.write(
                f"{chrom}\t{s}\t{e}\t{e - s}\t{label}\t{med:.2f}\t{subtype}\t{rep}\t{rdna}\n"
            )


if __name__ == "__main__":
    main()
