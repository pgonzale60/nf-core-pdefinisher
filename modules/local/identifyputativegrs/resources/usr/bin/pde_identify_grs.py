#!/usr/bin/env python3
"""
pde_identify_grs.py
Identify putative GRS (eliminated) and Core segments using:
1. Clustering of DELFIES neo-telomere breakpoints by strand within 3 bp.
2. A coverage threshold T_GRS (2-component EM for bimodal genomes; DELFIES-calibrated otherwise).
3. Chromosomal scaffolds: a GRS interval must be anchored by >= 1 DELFIES breakpoint;
   mosdepth coverage drops (depth <= T_GRS) only validate or close it. Repeats are never
   consulted and everything else is putative_core.
4. Non-chromosomal scaffolds: repeat-aware (RED) call with an explicit
   putative_GRS_uncertain class and a [min_elim_bp, max_elim_bp] range.
"""

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import random
import sys

__version__ = "1.3.1"


def get_script_fingerprint():
    try:
        with open(__file__, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except Exception:
        return "unknown"


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify assembly into putative_core and putative_GRS"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__} (sha256:{get_script_fingerprint()})",
    )
    parser.add_argument("--assembly", required=True, help="Assembly FASTA file")
    parser.add_argument(
        "--breakpoints-bed",
        required=False,
        default=None,
        help="Delfies breakpoints BED file",
    )
    parser.add_argument(
        "--regions-bed",
        required=False,
        default=None,
        help="Mosdepth window coverage BED file",
    )
    parser.add_argument(
        "--red-bed",
        required=False,
        default=None,
        help="RED repeat intervals BED; used only for non-chromosomal scaffolds",
    )
    parser.add_argument(
        "--chromosomes",
        required=False,
        default=None,
        help="Optional comma-separated list of canonical chromosome names",
    )
    parser.add_argument(
        "--policy",
        required=False,
        default=None,
        help="Path to assets/curation_policy.json (single source of truth)",
    )
    parser.add_argument(
        "--min-support",
        required=False,
        type=float,
        default=5.0,
        help="Minimum clustered read support for automaton consideration",
    )
    parser.add_argument(
        "--prefix", required=False, default="sample", help="Output prefix"
    )
    parser.add_argument(
        "--outdir", required=False, default=".", help="Output directory"
    )
    return parser.parse_args()


def get_chromosome_lengths(assembly_file):
    fai_file = f"{assembly_file}.fai"
    chr_lens = {}
    if os.path.exists(fai_file):
        try:
            with open(fai_file) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        chr_lens[parts[0]] = int(parts[1])
        except Exception:
            pass

    if not chr_lens:
        open_fn = gzip.open if assembly_file.endswith(".gz") else open
        curr_h = None
        curr_len = 0
        with open_fn(assembly_file, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    if curr_h:
                        chr_lens[curr_h] = curr_len
                    curr_h = line[1:].strip().split()[0]
                    curr_len = 0
                else:
                    curr_len += len(line.strip())
            if curr_h:
                chr_lens[curr_h] = curr_len
    return chr_lens


def parse_breakpoints(bed_file):
    breakpoints = {}
    bp_records_by_coord = {}
    if not bed_file or not os.path.exists(bed_file):
        return breakpoints, bp_records_by_coord

    try:
        with open(bed_file) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    chrom = parts[0]
                    pos = max(0, int(parts[1]))
                    strand = parts[5]
                    score = float(parts[4]) if len(parts) > 4 else 1.0
                    name = parts[3] if len(parts) > 3 else f"{chrom}_{pos}"
                    extra = parts[6:] if len(parts) > 6 else []
                    rec = {
                        "chrom": chrom,
                        "start": pos,
                        "end": int(parts[2]),
                        "name": name,
                        "score": score,
                        "strand": strand,
                        "extra": extra,
                    }
                    breakpoints.setdefault(chrom, []).append((pos, strand, score, rec))
                    bp_records_by_coord[(chrom, pos)] = rec
    except Exception:
        pass
    return breakpoints, bp_records_by_coord


def parse_coverage_windows(mosdepth_file):
    cov_windows = {}
    if not mosdepth_file or not os.path.exists(mosdepth_file):
        return cov_windows

    try:
        open_fn = gzip.open if mosdepth_file.endswith(".gz") else open
        with open_fn(mosdepth_file, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    c = parts[0]
                    s, e = int(parts[1]), int(parts[2])
                    val = float(parts[3])
                    cov_windows.setdefault(c, []).append((s, e, val))
    except Exception:
        pass
    return cov_windows




def get_median_depth(cov_windows, chrom, start, end):
    if chrom not in cov_windows:
        return 0.0
    vals = [w[2] for w in cov_windows[chrom] if w[0] >= start and w[1] <= end]
    if not vals:
        return 0.0
    vals.sort()
    mid = len(vals) // 2
    return vals[mid] if len(vals) % 2 != 0 else (vals[mid - 1] + vals[mid]) / 2.0


def cluster_breakpoints_by_strand(bp_list, max_dist=3, min_support=5.0):
    """
    Groups DELFIES signals by strand where consecutive signals are <= max_dist bp apart.
    Aggregates total read support in the cluster.
    Picks the position with the highest read support as representative.
    Filters out clusters with total support < min_support for macro-automaton consideration.
    """
    by_strand = {"+": [], "-": []}
    for b in bp_list:
        by_strand[b[1]].append(b)

    clustered = []
    for strand, items in by_strand.items():
        if not items:
            continue
        items.sort(key=lambda x: x[0])
        curr_cluster = [items[0]]
        for it in items[1:]:
            if it[0] - curr_cluster[-1][0] <= max_dist:
                curr_cluster.append(it)
            else:
                total_support = sum(x[2] for x in curr_cluster)
                if total_support >= min_support:
                    rep = max(curr_cluster, key=lambda x: x[2])
                    # Update cluster representative with total aggregated support
                    rep_updated = (rep[0], rep[1], total_support, rep[3])
                    clustered.append(rep_updated)
                curr_cluster = [it]
        if curr_cluster:
            total_support = sum(x[2] for x in curr_cluster)
            if total_support >= min_support:
                rep = max(curr_cluster, key=lambda x: x[2])
                rep_updated = (rep[0], rep[1], total_support, rep[3])
                clustered.append(rep_updated)

    clustered.sort(key=lambda x: x[0])
    return clustered


def fit_em_2comp(vals, max_iter=100, tol=1e-4):
    """
    Fit a 2-component Gaussian mixture model:
      Component 0: GRS / Eliminated (diploid/low coverage)
      Component 1: Somatic Core (retained diploid core coverage)
    """
    data = [math.log2(v + 1.0) for v in vals if v > 0]
    if not data:
        return {
            "mu_linear": [15.0, 150.0],
            "t_linear": 40.0,
            "autosome_core": 150.0,
        }
    n = len(data)
    if n > 25000:
        random.seed(42)
        sample = random.sample(data, 25000)
    else:
        sample = list(data)
    n_s = len(sample)
    s_data = sorted(sample)
    mu = [
        s_data[int(0.10 * n_s)],
        s_data[int(0.85 * n_s)],
    ]
    var = [0.5, 0.5]
    pi = [0.25, 0.75]
    for _ in range(max_iter):
        resp = []
        for x in sample:
            dens = []
            for k in range(2):
                sd = math.sqrt(max(var[k], 1e-4))
                p = (1.0 / (sd * math.sqrt(2 * math.pi))) * math.exp(
                    -0.5 * ((x - mu[k]) / sd) ** 2
                )
                dens.append(pi[k] * max(p, 1e-12))
            sdens = sum(dens)
            resp.append([d / sdens for d in dens])
        new_mu, new_var, new_pi = [], [], []
        for k in range(2):
            Nk = sum(resp[i][k] for i in range(n_s))
            if Nk < 1e-6:
                new_pi.append(pi[k])
                new_mu.append(mu[k])
                new_var.append(var[k])
                continue
            mk = sum(resp[i][k] * sample[i] for i in range(n_s)) / Nk
            vk = sum(resp[i][k] * ((sample[i] - mk) ** 2) for i in range(n_s)) / Nk
            new_mu.append(mk)
            new_var.append(max(vk, 0.01))
            new_pi.append(Nk / n_s)
        diff = sum(abs(new_mu[k] - mu[k]) for k in range(2))
        mu, var, pi = new_mu, new_var, new_pi
        if diff < tol:
            break
    order = sorted(range(2), key=lambda k: mu[k])
    sorted_mu = [mu[k] for k in order]

    mu_grs = sorted_mu[0]
    mu_core = sorted_mu[1]
    t_log2 = (mu_grs + mu_core) / 2.0
    t_linear = 2 ** t_log2 - 1.0

    return {
        "mu_linear": [2**m - 1.0 for m in sorted_mu],
        "t_linear": t_linear,
        "autosome_core": 2 ** sorted_mu[1] - 1.0,
    }


def parse_red_bed(red_file):
    """Load repeat intervals (BED) merged per scaffold; empty when no mask is available."""
    repeats = {}
    if not red_file or not os.path.exists(red_file):
        return repeats
    with open(red_file) as f:
        for line in f:
            cols = line.split()
            if len(cols) >= 3 and not line.startswith("#"):
                try:
                    repeats.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))
                except ValueError:
                    continue
    for c, ivs in repeats.items():
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        repeats[c] = [(s, e) for s, e in merged]
    return repeats


def quantize_drops(cov_windows, chrom, t_grs, merge_gap_bp):
    """Quantize coverage windows into low-coverage drops (depth <= t_grs).

    Returns (raw_drops, merged_drops); merged_drops joins drops separated by
    <= merge_gap_bp and is used only to trace where a single-anchored GRS closes.
    """
    raw = []
    for s, e, v in sorted(cov_windows.get(chrom, [])):
        if v <= t_grs:
            if raw and s <= raw[-1][1]:
                raw[-1] = (raw[-1][0], max(raw[-1][1], e))
            else:
                raw.append((s, e))
    merged = []
    for s, e in raw:
        if merged and s - merged[-1][1] <= merge_gap_bp:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return raw, merged


def drop_overlap_bp(drops, start, end):
    total = 0
    for ds, de in drops:
        if de <= start:
            continue
        if ds >= end:
            break
        total += max(0, min(end, de) - max(start, ds))
    return total


def delfies_bounded_grs(length, clustered_bps, raw_drops, merged_drops, rules):
    """GRS intervals on a chromosomal scaffold. Every interval is anchored by >= 1 DELFIES breakpoint.

    Repeats play no part here. Coverage only validates or closes an interval that a
    breakpoint has already opened. Returns merged [(start, end, subtype, note)].
    """
    frac = rules["drop_overlap_fraction"]
    micro = rules["micro_grs_max_bp"]
    prox = rules["anchor_proximity_bp"]
    min_closed = rules["min_closed_grs_bp"]

    def bp_end(b):
        return max(b[0], int(b[3].get("end", b[0])))

    def supported(s, e):
        return drop_overlap_bp(raw_drops, s, e) / max(1, e - s) >= frac

    bps = sorted(clustered_bps, key=lambda b: b[0])
    n = len(bps)
    found = []

    if n and bps[0][1] == "-":
        first_neg = bp_end(bps[0])
        if first_neg < micro or supported(0, first_neg):
            found.append((0, first_neg, "5_prime_distal_GRS", ""))

    i = 0
    while i < n:
        pos, strand = bps[i][0], bps[i][1]
        if strand == "+":
            # Concentric nested pair (+1, +2, -1, -2): the outer span is the GRS
            j = i
            while j < n and bps[j][1] == "+":
                j += 1
            k = j
            while k < n and bps[k][1] == "-":
                k += 1
            if j - i >= 2 and k - j >= 2:
                outer_s, outer_e = bps[i][0], bp_end(bps[k - 1])
                if supported(outer_s, outer_e):
                    inner = f"inner={bps[j - 1][0]}-{bp_end(bps[j])}"
                    found.append((outer_s, outer_e, "Nested_Variable_GRS_Outer", inner))
                    i = k
                    continue
            if i + 1 < n and bps[i + 1][1] == "-":
                pair_e = bp_end(bps[i + 1])
                span = pair_e - pos
                if span < micro:
                    found.append((pos, pair_e, "Micro_GRS_Paired", ""))
                    i += 2
                    continue
                if supported(pos, pair_e):
                    found.append((pos, pair_e, "Macro_GRS_Paired", ""))
                    i += 2
                    continue
            if i == n - 1:
                to_end = length - pos
                if to_end < micro or supported(pos, length):
                    found.append((pos, length, "3_prime_distal_GRS", ""))
                    i += 1
                    continue
            for ds, de in merged_drops:
                if abs(ds - pos) <= prox or ds <= pos <= de:
                    if de - pos >= min_closed:
                        found.append((pos, de, "Start_DELFIES_Anchored_End_Coverage_Closed", ""))
                    break
        elif i > 0:
            e_pos = bp_end(bps[i])
            for ds, de in merged_drops:
                if abs(de - e_pos) <= prox or ds <= e_pos <= de:
                    if e_pos - ds >= min_closed:
                        found.append((ds, e_pos, "Start_Coverage_Opened_End_DELFIES_Anchored", ""))
                    break
        i += 1

    if not found:
        return []
    found.sort(key=lambda x: x[0])
    merged = [list(found[0])]
    for s, e, t, note in found[1:]:
        if s <= merged[-1][1] + 100:
            merged[-1][1] = max(merged[-1][1], e)
            if "distal" not in merged[-1][2]:
                merged[-1][2] = t
            merged[-1][3] = merged[-1][3] or note
        else:
            merged.append([s, e, t, note])
    return [tuple(m) for m in merged]


def build_chromosomal_segments(chrom, length, grs_intervals, cov_windows):
    """Tile a chromosomal scaffold: DELFIES-anchored GRS intervals, everything else putative_core."""
    segs = []
    pos = 0
    for s, e, subtype, note in grs_intervals:
        s, e = max(0, s), min(length, e)
        if e <= s or s < pos:
            continue
        if s > pos:
            segs.append((chrom, pos, s, "putative_core", get_median_depth(cov_windows, chrom, pos, s), "-"))
        sub = f"{subtype}({note})" if note else subtype
        segs.append((chrom, s, e, "putative_GRS", get_median_depth(cov_windows, chrom, s, e), sub))
        pos = e
    if pos < length:
        segs.append((chrom, pos, length, "putative_core", get_median_depth(cov_windows, chrom, pos, length), "-"))
    return segs


def classify_unplaced(chrom, length, raw_drops, repeats, rules):
    """Non-chromosomal scaffold: repeat-aware call with an explicit uncertain class."""
    repeat_bp = sum(min(length, e) - max(0, s) for s, e in repeats.get(chrom, []) if e > 0 and s < length)
    non_rep_bp = max(0, length - repeat_bp)
    drop_bp = sum(e - s for s, e in raw_drops)
    if non_rep_bp < rules["unplaced_min_nonrepeat_bp"]:
        return "putative_GRS_uncertain", "fully_repetitive_unplaced"
    if drop_bp / max(1, length) >= rules["unplaced_grs_drop_fraction"]:
        return "putative_GRS", "unplaced_coverage_drop"
    return "putative_core", "-"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load curation policy (single source of truth)
    policy = load_curation_policy(args.policy)
    pde_policy = policy.get("pde_segmentation", {})
    max_grs_to_core_ratio = float(pde_policy.get("max_grs_to_core_ratio", 0.40))
    min_transition_contrast = float(pde_policy.get("min_transition_contrast", 2.0))
    default_grs_coverage_ratio = float(pde_policy.get("default_grs_coverage_ratio", 0.35))
    min_em_component_separation = float(pde_policy.get("min_em_component_separation", 0.45))
    rules = {
        "chromosomal_min_bp": int(pde_policy.get("chromosomal_min_bp", 1_000_000)),
        "drop_overlap_fraction": float(pde_policy.get("drop_overlap_fraction", 0.15)),
        "micro_grs_max_bp": int(pde_policy.get("micro_grs_max_bp", 1000)),
        "drop_merge_gap_bp": int(pde_policy.get("drop_merge_gap_bp", 1000)),
        "anchor_proximity_bp": int(pde_policy.get("anchor_proximity_bp", 2000)),
        "min_closed_grs_bp": int(pde_policy.get("min_closed_grs_bp", 1000)),
        "unplaced_min_nonrepeat_bp": int(pde_policy.get("unplaced_min_nonrepeat_bp", 1000)),
        "unplaced_grs_drop_fraction": float(pde_policy.get("unplaced_grs_drop_fraction", 0.50)),
    }
    raw_term = pde_policy.get("max_terminal_arm_bp")
    max_terminal_arm_bp = int(raw_term) if raw_term is not None and int(raw_term) > 0 else None
    raw_internal = pde_policy.get("max_internal_grs_bp")
    max_internal_grs_bp = int(raw_internal) if raw_internal is not None and int(raw_internal) > 0 else None

    chr_lens = get_chromosome_lengths(args.assembly)
    canon_chrom_set = None
    if args.chromosomes:
        canon_chrom_set = {c.strip() for c in args.chromosomes.split(",") if c.strip()}
        sys.stderr.write(
            f"[Chromosomes] Restricting canonical chromosomal evaluation to {len(canon_chrom_set)} chromosomes: "
            f"{', '.join(sorted(canon_chrom_set))}\n"
        )

    def is_canonical_chromosome(name, length):
        if canon_chrom_set is not None:
            return name in canon_chrom_set
        return length >= rules["chromosomal_min_bp"]

    breakpoints, bp_records_by_coord = parse_breakpoints(args.breakpoints_bed)
    cov_windows = parse_coverage_windows(args.regions_bed)

    # 2. Gather whole-genome window depths across major chromosomes (>= 1 Mb)
    all_core_vals = []
    for c, l in chr_lens.items():
        if is_canonical_chromosome(c, l) and c in cov_windows:
            all_core_vals.extend([w[2] for w in cov_windows[c]])

    # 3. Cluster breakpoints by strand within 3 bp
    clustered_bps_by_chr = {}
    for chrom in chr_lens:
        raw_bps = breakpoints.get(chrom, [])
        clustered_bps_by_chr[chrom] = cluster_breakpoints_by_strand(
            raw_bps, max_dist=3, min_support=args.min_support
        )

    # 4. Check whether coverage distribution is bimodal (e.g. APS4 with massive >60% GRS)
    em_fit = fit_em_2comp(all_core_vals) if all_core_vals else None
    is_bimodal = False
    if em_fit:
        mu_ratio = em_fit["mu_linear"][0] / max(em_fit["mu_linear"][1], 1e-3)
        if mu_ratio <= min_em_component_separation:
            is_bimodal = True
            # In a bimodal genome (e.g. APS4), Component 1 is ALWAYS somatic core (high peak)
            # and Component 0 is GRS (low peak), regardless of which component is more abundant.
            somatic_core_baseline = float(em_fit["mu_linear"][1])
            empirical_grs_depth = float(em_fit["mu_linear"][0])
            max_grs_threshold = max_grs_to_core_ratio * somatic_core_baseline
            t_grs = min(float(em_fit["t_linear"]), max_grs_threshold)
            sys.stderr.write(
                f"[Bimodal-EM] High peak (somatic core) = {somatic_core_baseline:.2f}x, "
                f"Low peak (GRS) = {empirical_grs_depth:.2f}x, separation ratio = {mu_ratio:.2f} <= {min_em_component_separation:.2f}; "
                f"GRS threshold = {t_grs:.2f}x (ceiling = {max_grs_threshold:.2f}x)\n"
            )

    if not is_bimodal:
        # Unimodal somatic core (e.g. PS1017 where GRS is only ~2-3%)
        if all_core_vals:
            s_vals = sorted(all_core_vals)
            somatic_core_baseline = float(s_vals[len(s_vals) // 2])
        else:
            somatic_core_baseline = 100.0
        max_grs_threshold = max_grs_to_core_ratio * somatic_core_baseline

        # Sample empirical GRS depth from DELFIES candidate intervals or small scaffolds
        candidate_grs_depths = []
        for c, l in chr_lens.items():
            c_bps = clustered_bps_by_chr.get(c, [])
            if c_bps:
                if c_bps[0][1] == "-" and (max_terminal_arm_bp is None or c_bps[0][0] <= max_terminal_arm_bp):
                    d = get_median_depth(cov_windows, c, 0, c_bps[0][0])
                    if 0 < d <= max_grs_threshold:
                        candidate_grs_depths.append(d)
                if c_bps[-1][1] == "+" and (max_terminal_arm_bp is None or (l - c_bps[-1][0]) <= max_terminal_arm_bp):
                    d = get_median_depth(cov_windows, c, c_bps[-1][0], l)
                    if 0 < d <= max_grs_threshold:
                        candidate_grs_depths.append(d)
                for i in range(len(c_bps) - 1):
                    if c_bps[i][1] == "+" and c_bps[i + 1][1] == "-":
                        s, e = c_bps[i][0], c_bps[i + 1][0]
                        if e > s and (max_internal_grs_bp is None or (e - s) <= max_internal_grs_bp):
                            d = get_median_depth(cov_windows, c, s, e)
                            if 0 < d <= max_grs_threshold:
                                candidate_grs_depths.append(d)
            elif (not is_canonical_chromosome(c, l)) and c in cov_windows:
                d = get_median_depth(cov_windows, c, 0, l)
                if 0 < d <= max_grs_threshold:
                    candidate_grs_depths.append(d)

        if candidate_grs_depths:
            cand_sorted = sorted(candidate_grs_depths)
            empirical_grs_depth = float(cand_sorted[len(cand_sorted) // 2])
            t_geom = (
                2
                ** (
                    (
                        math.log2(empirical_grs_depth + 1.0)
                        + math.log2(somatic_core_baseline + 1.0)
                    )
                    / 2.0
                )
                - 1.0
            )
            t_grs = min(t_geom, max_grs_threshold)
            sys.stderr.write(
                f"[Unimodal-Calibrated] Empirical GRS depth = {empirical_grs_depth:.2f}x across {len(candidate_grs_depths)} candidate intervals; "
                f"Somatic baseline = {somatic_core_baseline:.2f}x, GRS threshold = {t_grs:.2f}x (ceiling = {max_grs_threshold:.2f}x)\n"
            )
        else:
            t_grs = default_grs_coverage_ratio * somatic_core_baseline
            sys.stderr.write(
                f"[Policy-Default] Enforcing policy default GRS threshold = {t_grs:.2f}x "
                f"({default_grs_coverage_ratio * 100:.0f}% of somatic core baseline {somatic_core_baseline:.2f}x)\n"
            )

    segments = []
    repeats = parse_red_bed(args.red_bed)
    canon_grs_bp = unplaced_grs_bp = unplaced_uncertain_bp = unplaced_core_bp = 0

    for chrom, length in chr_lens.items():
        raw_drops, merged_drops = quantize_drops(
            cov_windows, chrom, t_grs, rules["drop_merge_gap_bp"]
        )
        if is_canonical_chromosome(chrom, length):
            # Chromosomal: DELFIES-bounded only; repeats are never consulted
            grs_ivs = delfies_bounded_grs(
                length, clustered_bps_by_chr.get(chrom, []), raw_drops, merged_drops, rules
            )
            chr_segs = build_chromosomal_segments(chrom, length, grs_ivs, cov_windows)
            canon_grs_bp += sum(e - s for _c, s, e, lab, _m, _t in chr_segs if lab == "putative_GRS")
            segments.extend(chr_segs)
        else:
            label, subtype = classify_unplaced(chrom, length, raw_drops, repeats, rules)
            med = get_median_depth(cov_windows, chrom, 0, length)
            segments.append((chrom, 0, length, label, med, subtype))
            if label == "putative_GRS_uncertain":
                unplaced_uncertain_bp += length
            elif label == "putative_GRS":
                unplaced_grs_bp += length
            else:
                unplaced_core_bp += length

    min_elim_bp = canon_grs_bp + unplaced_grs_bp
    max_elim_bp = min_elim_bp + unplaced_uncertain_bp
    est_elim_bp = (min_elim_bp + max_elim_bp) // 2

    # Write outputs
    prefix = args.prefix
    outdir = args.outdir
    bed_path = os.path.join(outdir, f"{prefix}.pde_segments.bed")
    tsv_path = os.path.join(outdir, f"{prefix}.grs_summary.tsv")
    cov_tsv_path = os.path.join(outdir, f"{prefix}.core_coverage.tsv")
    cov_txt_path = os.path.join(outdir, f"{prefix}.core_coverage.txt")

    with open(bed_path, "w") as out_bed, open(tsv_path, "w") as out_tsv, open(
        cov_tsv_path, "w"
    ) as out_cov_tsv, open(cov_txt_path, "w") as out_cov:
        out_cov_tsv.write("sample\tsomatic_core_coverage\tgrs_threshold\n")
        out_cov_tsv.write(f"{prefix}\t{somatic_core_baseline:.2f}\t{t_grs:.2f}\n")
        out_cov.write(f"{somatic_core_baseline:.2f}\n")
        out_tsv.write(
            f"# SOFTWARE=pde_identify_grs.py\tVERSION={__version__}\tFINGERPRINT={get_script_fingerprint()}\tGRS_THRESHOLD={t_grs:.2f}\tSOMATIC_CORE_BASELINE={somatic_core_baseline:.2f}\n"
        )
        out_tsv.write(
            f"# ELIMINATION_RANGE\tcanon_grs_bp={canon_grs_bp}\tunplaced_grs_bp={unplaced_grs_bp}"
            f"\tunplaced_uncertain_bp={unplaced_uncertain_bp}\tunplaced_core_bp={unplaced_core_bp}"
            f"\tmin_elim_bp={min_elim_bp}\tmax_elim_bp={max_elim_bp}\test_elim_bp={est_elim_bp}\n"
        )
        out_tsv.write(
            "chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype\n"
        )
        for chrom, s, e, label, med, subtype in segments:
            out_bed.write(f"{chrom}\t{s}\t{e}\t{label}\n")
            out_tsv.write(f"{chrom}\t{s}\t{e}\t{e-s}\t{label}\t{med:.2f}\t{subtype}\n")


if __name__ == "__main__":
    main()
