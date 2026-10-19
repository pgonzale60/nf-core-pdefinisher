#!/usr/bin/env python3
"""
pde_identify_grs.py
Identify putative GRS (eliminated) and Core segments using:
1. DELFIES neo-telomeric addition breakpoints parsed via Pushdown Automaton (PDA).
2. Clustering of breakpoints by strand within 3 bp.
3. Fallback 2-component Gaussian Mixture Model (GRS vs Somatic Core) for scaffolds
   without clean bracket matching or lacking DELFIES signals.
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
        help="Red repeat intervals BED file (deprecated; GRS is predominantly repetitive)",
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


def run_pushdown_automaton(
    chrom,
    length,
    clustered_bps,
    cov_windows,
    t_grs,
    max_terminal_arm_bp=None,
    max_internal_grs_bp=None,
    min_transition_contrast=2.0,
):
    """
    Models breakpoint pairing as a Pushdown Automaton (Dyck language):
      '+' represents open parenthesis '(' (exits core into GRS)
      '-' represents close parenthesis ')' (exits GRS into core)

    Validates that terminal boundaries have true somatic core contrast on their
    retained flank and do not exceed biologically realistic arm lengths if configured.

    Returns:
      (is_neat_match, grs_intervals, frontier_bps)
    """
    if not clustered_bps:
        return False, [], []

    grs_intervals = []
    frontier_bps = []
    stack = []
    idx = 0

    # Optional 5' terminal GRS: starts with '-' (closing bracket from 0)
    # MUST have somatic core contrast downstream (> t_grs) and satisfy max arm length if specified
    if clustered_bps[0][1] == "-":
        bp_pos = clustered_bps[0][0]
        if max_terminal_arm_bp is not None and bp_pos > max_terminal_arm_bp:
            return False, [], []
        med_up = get_median_depth(cov_windows, chrom, 0, bp_pos)
        downstream_end = min(length, bp_pos + 50000)
        med_downstream = get_median_depth(cov_windows, chrom, bp_pos, downstream_end)
        if (
            med_downstream > t_grs
            and (
                med_up <= t_grs
                or med_downstream >= min_transition_contrast * max(med_up, 0.1)
            )
        ):
            grs_intervals.append((0, bp_pos))
            frontier_bps.append(clustered_bps[0])
            idx = 1
        else:
            return False, [], []

    for i in range(idx, len(clustered_bps)):
        b = clustered_bps[i]
        pos, strand = b[0], b[1]
        if strand == "+":
            stack.append(b)
        elif strand == "-":
            if not stack:
                # Unmatched closing breakpoint in the interior (e.g. isolated point scar).
                # Skip rather than failing the whole chromosome.
                continue
            open_b = stack.pop()
            seg_len = pos - open_b[0]
            if max_internal_grs_bp is not None and seg_len > max_internal_grs_bp:
                # Exceeds max internal GRS length
                continue
            # Validate coverage drop inside the candidate internal segment
            med_internal = get_median_depth(cov_windows, chrom, open_b[0], pos)
            up_start = max(0, open_b[0] - 50000)
            med_up = get_median_depth(cov_windows, chrom, up_start, open_b[0])
            down_end = min(length, pos + 50000)
            med_down = get_median_depth(cov_windows, chrom, pos, down_end)
            flank_min = min(med_up, med_down)
            if med_internal <= t_grs or med_internal <= 0.45 * flank_min or seg_len < 1000:
                grs_intervals.append((open_b[0], pos))
                frontier_bps.extend([open_b, b])

    # Optional 3' terminal GRS: look for unclosed '+' within max_terminal_arm_bp from chromosome end (if configured)
    candidates_3p = [b for b in stack if (max_terminal_arm_bp is None or (length - b[0]) <= max_terminal_arm_bp)]
    if candidates_3p:
        # Pick the leftmost valid 3' candidate (most comprehensive terminal range)
        for open_b in candidates_3p:
            bp_pos = open_b[0]
            upstream_start = max(0, bp_pos - 50000)
            med_upstream = get_median_depth(cov_windows, chrom, upstream_start, bp_pos)
            med_down = get_median_depth(cov_windows, chrom, bp_pos, length)
            if (
                med_upstream > t_grs
                and (
                    med_down <= t_grs
                    or med_upstream >= min_transition_contrast * max(med_down, 0.1)
                )
            ):
                grs_intervals.append((bp_pos, length))
                frontier_bps.append(open_b)
                break

    grs_intervals.sort(key=lambda x: x[0])
    # Deduplicate and ensure intervals do not overlap
    merged = []
    for it in grs_intervals:
        if not merged:
            merged.append(it)
        else:
            prev = merged[-1]
            if it[0] < prev[1]:
                # Overlap: merge into union
                merged[-1] = (prev[0], max(prev[1], it[1]))
            else:
                merged.append(it)

    if not merged:
        return False, [], []

    return True, merged, frontier_bps


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


def build_chromosome_segments(chrom, length, all_grs, cov_windows, t_grs):
    """Build and classify segments for a chromosome, checking depth of all intervals.

    Never assumes flanks around candidate GRS intervals are somatic core. Contiguous
    adjacent segments sharing the same classification label are merged.
    """
    if not all_grs:
        chr_med = get_median_depth(cov_windows, chrom, 0, length)
        label = "putative_GRS" if chr_med <= t_grs else "putative_core"
        return [(chrom, 0, length, label, chr_med)]

    chr_segments = []
    curr_pos = 0
    for grs_s, grs_e, grs_med in all_grs:
        if grs_s > curr_pos:
            flank_med = get_median_depth(cov_windows, chrom, curr_pos, grs_s)
            flank_label = "putative_core" if flank_med > t_grs else "putative_GRS"
            chr_segments.append((chrom, curr_pos, grs_s, flank_label, flank_med))
        chr_segments.append((chrom, grs_s, grs_e, "putative_GRS", grs_med))
        curr_pos = grs_e

    if curr_pos < length:
        flank_med = get_median_depth(cov_windows, chrom, curr_pos, length)
        flank_label = "putative_core" if flank_med > t_grs else "putative_GRS"
        chr_segments.append((chrom, curr_pos, length, flank_label, flank_med))

    merged_chr_segments = []
    for c, s, e, lab, med in chr_segments:
        if (
            merged_chr_segments
            and merged_chr_segments[-1][3] == lab
            and merged_chr_segments[-1][2] == s
        ):
            prev_c, prev_s, _prev_e, prev_lab, _prev_med = merged_chr_segments[-1]
            comb_med = get_median_depth(cov_windows, prev_c, prev_s, e)
            merged_chr_segments[-1] = (prev_c, prev_s, e, prev_lab, comb_med)
        else:
            merged_chr_segments.append((c, s, e, lab, med))

    return merged_chr_segments


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
    raw_term = pde_policy.get("max_terminal_arm_bp")
    max_terminal_arm_bp = int(raw_term) if raw_term is not None and int(raw_term) > 0 else None
    raw_internal = pde_policy.get("max_internal_grs_bp")
    max_internal_grs_bp = int(raw_internal) if raw_internal is not None and int(raw_internal) > 0 else None

    chr_lens = get_chromosome_lengths(args.assembly)
    breakpoints, bp_records_by_coord = parse_breakpoints(args.breakpoints_bed)
    cov_windows = parse_coverage_windows(args.regions_bed)

    # 2. Gather whole-genome window depths across major chromosomes (>= 1 Mb)
    all_core_vals = []
    for c, l in chr_lens.items():
        if l >= 1_000_000 and c in cov_windows:
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
            elif l < 1_000_000 and c in cov_windows:
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
    frontier_breakpoints = []

    for chrom, length in chr_lens.items():
        raw_bps = breakpoints.get(chrom, [])
        chr_med = get_median_depth(cov_windows, chrom, 0, length)
        clustered_bps = clustered_bps_by_chr.get(chrom, [])

        if not raw_bps:
            # Scaffolds with 0 breakpoints: classify based on calibrated threshold
            label = "putative_GRS" if chr_med <= t_grs else "putative_core"
            segments.append((chrom, 0, length, label, chr_med))
            continue

        # 5. Attempt Pushdown Automaton (bracket matching)
        pda_success, pda_grs_intervals, pda_frontier = run_pushdown_automaton(
            chrom,
            length,
            clustered_bps,
            cov_windows,
            t_grs,
            max_terminal_arm_bp=max_terminal_arm_bp,
            max_internal_grs_bp=max_internal_grs_bp,
            min_transition_contrast=min_transition_contrast,
        )

        all_grs = []
        if pda_success:
            sys.stderr.write(
                f"[PDA] Chromosome {chrom}: Successfully resolved {len(pda_grs_intervals)} GRS regions via matching parentheses.\n"
            )
            for s, e in pda_grs_intervals:
                med = get_median_depth(cov_windows, chrom, s, e)
                all_grs.append((s, e, med))
            for b in pda_frontier:
                if (chrom, b[0]) in bp_records_by_coord:
                    frontier_breakpoints.append(bp_records_by_coord[(chrom, b[0])])
        else:
            # 6. Fallback: Coverage-guided segmentation using policy ceilings & local fold-change
            sys.stderr.write(
                f"[Fallback] Chromosome {chrom}: PDA did not resolve neatly. Using coverage-guided segmentation.\n"
            )
            bps = sorted(clustered_bps, key=lambda x: x[0])

            # 5' Terminal GRS
            chosen_5p_grs = None
            chosen_5p_bp = None
            minus_bps = [
                bp
                for bp in bps
                if bp[1] == "-"
                and (max_terminal_arm_bp is None or bp[0] <= max_terminal_arm_bp)
            ]
            for bp in reversed(minus_bps):
                pos = bp[0]
                if pos <= 0 or pos >= length:
                    continue
                med_grs = get_median_depth(cov_windows, chrom, 0, pos)
                adj_end = min(length, pos + 50000)
                med_adj = get_median_depth(cov_windows, chrom, pos, adj_end)

                if med_adj > t_grs and (
                    med_grs <= t_grs
                    or med_adj >= min_transition_contrast * max(med_grs, 0.1)
                ):
                    chosen_5p_grs = (0, pos, med_grs)
                    chosen_5p_bp = pos
                    break

            # 3' Terminal GRS
            chosen_3p_grs = None
            chosen_3p_bp = None
            plus_bps = [
                bp
                for bp in bps
                if bp[1] == "+"
                and (
                    max_terminal_arm_bp is None
                    or (length - bp[0]) <= max_terminal_arm_bp
                )
            ]
            for bp in plus_bps:
                pos = bp[0]
                if pos <= 0 or pos >= length:
                    continue
                if chosen_5p_grs and pos <= chosen_5p_grs[1]:
                    continue
                med_grs = get_median_depth(cov_windows, chrom, pos, length)
                adj_start = max(0, pos - 50000)
                med_adj = get_median_depth(cov_windows, chrom, adj_start, pos)

                if med_adj > t_grs and (
                    med_grs <= t_grs
                    or med_adj >= min_transition_contrast * max(med_grs, 0.1)
                ):
                    chosen_3p_grs = (pos, length, med_grs)
                    chosen_3p_bp = pos
                    break

            # Internal GRS with relative fold-drop
            internal_min = chosen_5p_grs[1] if chosen_5p_grs else 0
            internal_max = chosen_3p_grs[0] if chosen_3p_grs else length

            candidate_internal_grs = []
            internal_plus = [
                bp
                for bp in bps
                if bp[1] == "+" and internal_min <= bp[0] < internal_max
            ]
            internal_minus = [
                bp
                for bp in bps
                if bp[1] == "-" and internal_min < bp[0] <= internal_max
            ]

            for p_left in internal_plus:
                for p_right in internal_minus:
                    if p_right[0] > p_left[0] and (
                        max_internal_grs_bp is None
                        or (p_right[0] - p_left[0]) <= max_internal_grs_bp
                    ):
                        med_grs = get_median_depth(
                            cov_windows, chrom, p_left[0], p_right[0]
                        )
                        seg_len = p_right[0] - p_left[0]
                        up_start = max(internal_min, p_left[0] - 50000)
                        med_up = get_median_depth(
                            cov_windows, chrom, up_start, p_left[0]
                        )
                        down_end = min(internal_max, p_right[0] + 50000)
                        med_down = get_median_depth(
                            cov_windows, chrom, p_right[0], down_end
                        )

                        flank_min = min(med_up, med_down)
                        if (
                            (
                                med_grs <= t_grs
                                or med_grs <= 0.45 * flank_min
                                or seg_len < 1000
                            )
                            and med_up > t_grs
                            and med_down > t_grs
                        ):
                            candidate_internal_grs.append(
                                (
                                    p_left[0],
                                    p_right[0],
                                    seg_len,
                                    med_grs,
                                    p_left[0],
                                    p_right[0],
                                )
                            )

            candidate_internal_grs.sort(key=lambda x: x[2], reverse=True)
            chosen_internal_grs = []
            chosen_internal_bps = []
            for s, e, l, med, b_l, b_r in candidate_internal_grs:
                overlaps = False
                for cs, ce, _ in chosen_internal_grs:
                    if not (e <= cs or s >= ce):
                        overlaps = True
                        break
                if not overlaps:
                    chosen_internal_grs.append((s, e, med))
                    chosen_internal_bps.extend([b_l, b_r])

            chosen_internal_grs.sort(key=lambda x: x[0])

            if chosen_5p_grs:
                all_grs.append((chosen_5p_grs[0], chosen_5p_grs[1], chosen_5p_grs[2]))
                if (chrom, chosen_5p_bp) in bp_records_by_coord:
                    frontier_breakpoints.append(bp_records_by_coord[(chrom, chosen_5p_bp)])

            for s, e, med in chosen_internal_grs:
                all_grs.append((s, e, med))
            for b_pos in chosen_internal_bps:
                if (chrom, b_pos) in bp_records_by_coord:
                    frontier_breakpoints.append(bp_records_by_coord[(chrom, b_pos)])

            if chosen_3p_grs:
                all_grs.append((chosen_3p_grs[0], chosen_3p_grs[1], chosen_3p_grs[2]))
                if (chrom, chosen_3p_bp) in bp_records_by_coord:
                    frontier_breakpoints.append(bp_records_by_coord[(chrom, chosen_3p_bp)])

        chr_segs = build_chromosome_segments(
            chrom, length, all_grs, cov_windows, t_grs
        )
        segments.extend(chr_segs)

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
            "chromosome\tstart\tend\tlength\tclassification\tmedian_depth\n"
        )
        for chrom, s, e, label, med in segments:
            out_bed.write(f"{chrom}\t{s}\t{e}\t{label}\n")
            out_tsv.write(f"{chrom}\t{s}\t{e}\t{e-s}\t{label}\t{med:.2f}\n")


if __name__ == "__main__":
    main()
