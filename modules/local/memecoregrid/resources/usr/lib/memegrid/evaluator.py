"""MEME output parser, Information Content calculator, and motif scoring evaluator."""

from collections import Counter
import math
import os
from pathlib import Path
import re
from typing import Any

from cbm_defaults import SSOT


def calculate_consensus_entropy(consensus: str) -> float:
    """Calculate Shannon entropy (in bits) across the characters of a consensus sequence."""
    if not consensus:
        return 0.0
    cons_upper = consensus.upper()
    n = len(cons_upper)
    counts = Counter(cons_upper)
    entropy = 0.0
    for count in counts.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def is_telomeric_repeat_echo(
    consensus: str, repeat: str = SSOT["telomere_repeat"]
) -> tuple[bool, float]:
    """Check if motif consensus is a telomeric repeat echo (>=50% exact 6-mers or >=70% in-phase 1-mismatch coverage)."""
    if not consensus or not repeat:
        return False, 0.0
    cons = consensus.upper().strip()
    rep = repeat.upper().strip()
    k = len(rep)
    n = len(cons)
    if n < k:
        return False, 0.0

    rc_rep = _revcomp(rep)
    fwd_perms = {(rep * 2)[i : i + k] for i in range(k)}
    rev_perms = {(rc_rep * 2)[i : i + k] for i in range(k)}

    n_kmers = n - k + 1
    fwd_exact = sum(1 for i in range(n_kmers) if cons[i : i + k] in fwd_perms)
    rev_exact = sum(1 for i in range(n_kmers) if cons[i : i + k] in rev_perms)
    exact_frac = max(fwd_exact, rev_exact) / n_kmers

    # Check in-phase 1-mismatch coverage across all cyclic phases of forward and reverse repeat
    best_cov_1mm = 0.0
    reps_needed = (n // k) + 3
    for base_rep in (rep, rc_rep):
        full_periodic = base_rep * reps_needed
        for phase in range(k):
            tmpl = full_periodic[phase : phase + n]
            covered = [False] * n
            mm1_count = 0
            for i in range(n_kmers):
                mm = sum(
                    1
                    for a, b in zip(cons[i : i + k], tmpl[i : i + k], strict=False)
                    if a != b
                )
                if mm <= 1:
                    mm1_count += 1
                    for j in range(i, i + k):
                        covered[j] = True
            if mm1_count >= 2:
                best_cov_1mm = max(best_cov_1mm, sum(covered) / n)

    is_echo = (exact_frac > SSOT["telomeric_echo_exact_kmer_cutoff"]) or (
        n >= 18 and best_cov_1mm >= SSOT["telomeric_echo_1mm_cov_cutoff"]
    )
    return is_echo, round(max(exact_frac, best_cov_1mm), 4)


def detect_tandem_repeats(
    seq: str,
    min_period: int = 2,
    max_period: int = 10,
    min_copies: float = 2.0,
    min_identity: float = 0.85,
    target_unit: str | None = None,
) -> tuple[bool, int | None, str | None]:
    """Detects tandem periodic repeats via lagged self-comparison.

    If target_unit is provided, restricts to that monomer's period and cyclic permutations.
    """
    if not seq:
        return False, None, None
    s = seq.upper()
    n = len(s)

    if target_unit:
        tu = target_unit.upper()
        p = len(tu)
        rc_tu = tu.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        valid_cyclic_kmers = {(tu * 2)[i : i + p] for i in range(p)} | {
            (rc_tu * 2)[i : i + p] for i in range(p)
        }
        periods = [p]
    else:
        valid_cyclic_kmers = None
        periods = range(min_period, min(max_period, n // int(min_copies)) + 1)

    for p in periods:
        lag_len = int((min_copies - 1) * p)
        if n < p + lag_len:
            continue
        for i in range(n - p - lag_len + 1):
            s1 = s[i : i + lag_len]
            s2 = s[i + p : i + p + lag_len]
            matches = sum(1 for a, b in zip(s1, s2, strict=False) if a == b)
            if matches / lag_len >= min_identity:
                monomer = s[i : i + p]
                if valid_cyclic_kmers and monomer not in valid_cyclic_kmers:
                    continue
                if not valid_cyclic_kmers and len(set(monomer)) <= 1:
                    continue
                return True, p, monomer

    return False, None, None


def is_tandem_repeat_motif(
    consensus: str,
    min_period: int = 2,
    max_period: int = 10,
    min_copies: float = 2.0,
    min_identity: float = 0.85,
) -> bool:
    """Check whether consensus sequence has tandem repeating periodicity (e.g. 2-10 bp repeats >= 2 copies)."""
    found, _, _ = detect_tandem_repeats(
        consensus,
        min_period=min_period,
        max_period=max_period,
        min_copies=min_copies,
        min_identity=min_identity,
    )
    return found


def evaluate_meme_motif(
    consensus: str,
    width: int,
    mean_ic: float,
    occupancy: float,
    pos_iqr: float,
    median_offset: float,
    evalue: float,
    site_density: float = 1.0,
    telomere_repeat: str = SSOT["telomere_repeat"],
    proximity_decay_bp: float | None = None,
    max_density_credit: float | None = None,
) -> dict[str, Any]:
    """Evaluate motif quality, rejection criteria, and compute composite score."""
    if proximity_decay_bp is None:
        proximity_decay_bp = float(SSOT.get("cbm_catalytic_proximity_decay_bp", 75.0))
    if max_density_credit is None:
        max_density_credit = float(SSOT.get("max_anr_site_density_credit", 2.0))

    effective_density = min(max_density_credit, site_density)
    site_density_factor = 1.0 + math.log(max(1.0, effective_density))

    cons_entropy = calculate_consensus_entropy(consensus)
    homopolymer_factor = (
        1.0
        if cons_entropy >= 1.2
        else (cons_entropy / 1.2 if cons_entropy >= 0.6 else 0.05)
    )

    iqr_penalty = 1.0 / (1.0 + math.log1p(pos_iqr / 15.0))
    proximity_factor = math.exp(-abs(median_offset) / proximity_decay_bp)

    is_telo_echo, telo_echo_frac = is_telomeric_repeat_echo(
        consensus, repeat=telomere_repeat
    )
    is_tandem_rep = is_tandem_repeat_motif(consensus)

    if is_telo_echo:
        composite_score = 0.0
        rejection_flag = "REJECTED_TELOMERIC_ECHO"
    elif is_tandem_rep:
        composite_score = 0.0
        rejection_flag = "REJECTED_TANDEM_SATELLITE"
    elif evalue > SSOT["core_evalue_cutoff"]:
        composite_score = 0.0
        rejection_flag = "REJECTED_INSIGNIFICANT_EVALUE"
    else:
        width_factor = math.pow(width, 0.2)
        composite_score = (
            (mean_ic * width_factor)
            * site_density_factor
            * occupancy
            * iqr_penalty
            * proximity_factor
            * homopolymer_factor
        )
        rejection_flag = "PASS"

    return {
        "composite_score": round(composite_score, 4),
        "rejection_flag": rejection_flag,
        "site_density_factor": round(site_density_factor, 4),
        "proximity_factor": round(proximity_factor, 4),
        "iqr_penalty": round(iqr_penalty, 4),
        "telomeric_echo_fraction": telo_echo_frac,
        "is_tandem_repeat": is_tandem_rep,
        "consensus_entropy": cons_entropy,
    }


def parse_meme_txt(
    txt_path: str | Path,
    radius: int,
    total_sequences: int = 1,
    telomere_repeat: str = SSOT["telomere_repeat"],
) -> dict[str, Any] | None:
    """Parse MEME text output and compute comprehensive motif metrics."""
    if not os.path.exists(txt_path):
        return None

    with open(txt_path, encoding="utf-8") as f:
        text = f.read()

    # Unwrap backslash-wrapped continuation lines if any
    unwrapped = re.sub(r"\\\s*\n\s*", "", text)

    # 1. Header regex for width, sites, E-value
    m = re.search(
        r"letter-probability matrix:\s+alength=\s*(\d+)\s+w=\s*(\d+)\s+nsites=\s*(\d+)\s+E=\s*(\S+)",
        unwrapped,
    )
    if not m:
        return None

    w = int(m.group(2))
    nsites = int(m.group(3))
    evalue_str = m.group(4)
    try:
        evalue = float(evalue_str)
    except ValueError:
        evalue = 1.0

    # 2. Consensus sequence
    cons_m = re.search(r"Multilevel\s+([A-Za-z]+)", unwrapped)
    consensus = cons_m.group(1).upper() if cons_m else ""

    # 3. Parse letter-probability matrix to compute Information Content
    idx = unwrapped.find(m.group(0)) + len(m.group(0))
    matrix_lines = unwrapped[idx:].strip().splitlines()
    ic_per_col = []
    for line in matrix_lines:
        line = line.strip()
        if not line or line.startswith("URL") or line.startswith("-"):
            break
        parts = [float(x) for x in line.split()]
        if len(parts) == 4:
            # IC = 2 - sum(-p * log2(p))
            h = sum(-p * math.log2(p) for p in parts if p > 0)
            ic_per_col.append(max(0.0, 2.0 - h))
        if len(ic_per_col) == w:
            break

    total_ic = sum(ic_per_col)
    mean_ic = total_ic / w if w > 0 else 0.0

    # 4. Parse individual site positions relative to Position 0
    # In sliced window of radius R, Position 0 is at 1-based index R
    in_sites = False
    site_records = []
    for line in unwrapped.splitlines():
        if "sites sorted by position p-value" in line:
            in_sites = True
            continue
        if in_sites:
            if line.startswith("---"):
                continue
            if "block diagrams" in line or line.startswith("MOTIF"):
                break
            parts = line.split()
            if len(parts) >= 4 and ("+" in parts or "-" in parts):
                seq_name = parts[0]
                strand = parts[1]
                try:
                    start_pos = int(parts[2])
                    # Distance relative to Position 0 (dynamic pos0 from sequence header, fallback to radius)
                    p0_match = re.search(r"_p0_(\d+)", seq_name)
                    pos0 = int(p0_match.group(1)) if p0_match else radius
                    # Distance from motif center to Position 0
                    motif_center = start_pos + (w - 1) / 2.0
                    dist_to_pos0 = motif_center - pos0
                    site_records.append({
                        "seq": seq_name,
                        "strand": strand,
                        "start": start_pos,
                        "dist_to_pos0": dist_to_pos0,
                    })
                except ValueError:
                    continue

    unique_seqs = len({s["seq"] for s in site_records})
    occupancy = unique_seqs / total_sequences if total_sequences > 0 else 0.0

    site_density = (nsites / unique_seqs) if unique_seqs > 0 else 1.0

    # 5. Positional IQR
    dists = sorted(s["dist_to_pos0"] for s in site_records)
    if len(dists) >= 4:
        q25 = dists[len(dists) // 4]
        q75 = dists[(3 * len(dists)) // 4]
        pos_iqr = float(q75 - q25)
        median_offset = float(dists[len(dists) // 2])
    elif dists:
        pos_iqr = float(max(dists) - min(dists))
        median_offset = float(dists[0])
    else:
        pos_iqr = 50.0
        median_offset = 0.0

    eval_res = evaluate_meme_motif(
        consensus=consensus,
        width=w,
        mean_ic=mean_ic,
        occupancy=occupancy,
        pos_iqr=pos_iqr,
        median_offset=median_offset,
        evalue=evalue,
        site_density=site_density,
        telomere_repeat=telomere_repeat,
    )

    return {
        "width": w,
        "nsites": nsites,
        "evalue": evalue,
        "evalue_str": evalue_str,
        "consensus": consensus,
        "total_ic": round(total_ic, 2),
        "mean_ic": round(mean_ic, 3),
        "consensus_entropy": eval_res["consensus_entropy"],
        "unique_occupied_sites": unique_seqs,
        "occupancy": round(occupancy, 4),
        "site_density_factor": eval_res["site_density_factor"],
        "median_offset_to_pos0": round(median_offset, 2),
        "positional_iqr": round(pos_iqr, 2),
        "telomeric_echo_fraction": eval_res["telomeric_echo_fraction"],
        "rejection_flag": eval_res["rejection_flag"],
        "composite_score": eval_res["composite_score"],
    }
