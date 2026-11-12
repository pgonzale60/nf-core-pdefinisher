"""Tier-2 Auxiliary motif grid search on core-masked break windows."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

from cbm_defaults import SSOT

AUX_GRID_CONDITIONS = [
    # Fixed width 6 (hexamer)
    {"id": "A01_anr_w6_revcomp", "mod": "anr", "w_args": ["-w", "6"], "revcomp": True},
    {"id": "A02_anr_w6_fwd", "mod": "anr", "w_args": ["-w", "6"], "revcomp": False},
    {"id": "A03_zoops_w6_revcomp", "mod": "zoops", "w_args": ["-w", "6"], "revcomp": True},
    {"id": "A04_zoops_w6_fwd", "mod": "zoops", "w_args": ["-w", "6"], "revcomp": False},
    # Fixed width 8 (octamer)
    {"id": "A05_anr_w8_revcomp", "mod": "anr", "w_args": ["-w", "8"], "revcomp": True},
    {"id": "A06_anr_w8_fwd", "mod": "anr", "w_args": ["-w", "8"], "revcomp": False},
    {"id": "A07_zoops_w8_revcomp", "mod": "zoops", "w_args": ["-w", "8"], "revcomp": True},
    {"id": "A08_zoops_w8_fwd", "mod": "zoops", "w_args": ["-w", "8"], "revcomp": False},
    # Fixed width 11
    {"id": "A09_anr_w11_revcomp", "mod": "anr", "w_args": ["-w", "11"], "revcomp": True},
    {"id": "A10_anr_w11_fwd", "mod": "anr", "w_args": ["-w", "11"], "revcomp": False},
    {"id": "A11_zoops_w11_revcomp", "mod": "zoops", "w_args": ["-w", "11"], "revcomp": True},
    {"id": "A12_zoops_w11_fwd", "mod": "zoops", "w_args": ["-w", "11"], "revcomp": False},
    # Flexible width 6 to 12 bp
    {"id": "A13_anr_w6_12_revcomp", "mod": "anr", "w_args": ["-minw", "6", "-maxw", "12"], "revcomp": True},
    {"id": "A14_anr_w6_12_fwd", "mod": "anr", "w_args": ["-minw", "6", "-maxw", "12"], "revcomp": False},
    {"id": "A15_zoops_w6_12_revcomp", "mod": "zoops", "w_args": ["-minw", "6", "-maxw", "12"], "revcomp": True},
    {"id": "A16_zoops_w6_12_fwd", "mod": "zoops", "w_args": ["-minw", "6", "-maxw", "12"], "revcomp": False},
]


def parse_aux_meme_results(
    txt_path: str | Path, center_pos: int = SSOT["max_soma_flank_len"]
) -> dict[str, Any] | None:
    """Parse secondary MEME search results."""
    if not os.path.exists(txt_path):
        return None
    with open(txt_path, encoding="utf-8") as f:
        text = f.read()

    m = re.search(
        r"letter-probability matrix:\s+alength=\s*(\d+)\s+w=\s*(\d+)\s+nsites=\s*(\d+)\s+E=\s*(\S+)",
        text,
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

    cons_m = re.search(r"Multilevel\s+([A-Za-z]+)", text)
    consensus = cons_m.group(1).upper() if cons_m else ""

    idx = text.find(m.group(0)) + len(m.group(0))
    matrix_lines = text[idx:].strip().splitlines()
    ic_list = []
    for line in matrix_lines:
        line = line.strip()
        if not line or line.startswith("URL") or line.startswith("-"):
            break
        parts = [float(x) for x in line.split()]
        if len(parts) == 4:
            h = sum(-p * math.log2(p) for p in parts if p > 0)
            ic_list.append(max(0.0, 2.0 - h))
        if len(ic_list) == w:
            break

    total_ic = sum(ic_list)
    mean_ic = total_ic / w if w > 0 else 0.0

    # Parse individual sites
    in_sites = False
    site_records = []
    for line in text.splitlines():
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
                    start = int(parts[2])
                    dist_to_center = abs(start - center_pos)
                    site_records.append({
                        "seq": seq_name,
                        "strand": strand,
                        "start": start,
                        "dist": dist_to_center,
                    })
                except ValueError:
                    continue

    unique_seqs = len({s["seq"] for s in site_records})
    proximal_sites = sum(1 for s in site_records if 6 <= s["dist"] <= 50)

    # Entropy of consensus
    cons_counts = Counter(consensus)
    c_ent = 0.0
    for cnt in cons_counts.values():
        p = cnt / len(consensus) if len(consensus) > 0 else 0
        if p > 0:
            c_ent -= p * math.log2(p)

    return {
        "width": w,
        "nsites": nsites,
        "evalue": evalue,
        "evalue_str": evalue_str,
        "consensus": consensus,
        "total_ic": round(total_ic, 2),
        "mean_ic": round(mean_ic, 3),
        "consensus_entropy": round(c_ent, 3),
        "unique_sequences_hit": unique_seqs,
        "proximal_sites_6_50bp": proximal_sites,
    }


def _run_single_aux_condition(
    meme_bin: str | None,
    cond: dict[str, Any],
    masked_fa_str: str,
    cond_dir_str: str,
    bfile: str | None = None,
) -> dict[str, Any] | None:
    """Run MEME for a single auxiliary condition and parse metrics."""
    cond_dir = Path(cond_dir_str)
    meme_txt = cond_dir / "meme.txt"
    cid = cond["id"]

    if not meme_txt.exists():
        if not meme_bin:
            raise RuntimeError("Required binary 'meme' not found in PATH")
        cmd = [
            meme_bin,
            masked_fa_str,
            "-dna",
            "-nmotifs",
            "1",
            "-mod",
            cond["mod"],
            "-oc",
            cond_dir_str,
        ] + cond["w_args"]
        if bfile and os.path.exists(bfile):
            cmd.extend(["-bfile", str(bfile)])
        if cond["revcomp"]:
            cmd.append("-revcomp")
        try:
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr.strip() if e.stderr else "Unknown error"
            raise RuntimeError(
                f"Auxiliary MEME condition '{cid}' failed with exit code {e.returncode}: {err_msg}"
            ) from e

    metrics = parse_aux_meme_results(str(meme_txt))
    if metrics:
        metrics["id"] = cid
        metrics["mod"] = cond["mod"]
        metrics["revcomp"] = cond["revcomp"]
        metrics["meme_txt_path"] = str(meme_txt)
        return metrics
    return None


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def is_core_echo(aux_consensus: str, core_consensus: str) -> bool:
    """Return True if aux_consensus is a substring/superstring of core_consensus or its revcomp."""
    if not aux_consensus or not core_consensus:
        return False
    a = aux_consensus.upper().strip()
    c = core_consensus.upper().strip()
    rc_c = _revcomp(c)
    return (a in c) or (c in a) or (a in rc_c) or (rc_c in a)


def detect_tandem_repeats(
    seq: str,
    min_period: int = 2,
    max_period: int = 10,
    min_copies: float = 2.0,
    min_identity: float = 0.85,
) -> tuple[bool, int, float]:
    """Detect if sequence is a tandem repeat using lagged self-comparison."""
    if not seq:
        return False, 0, 0.0
    s = seq.upper().strip()
    n = len(s)
    upper_p = min(max_period, int(n // min_copies))
    for p in range(min_period, upper_p + 1):
        if n < p * min_copies:
            continue
        lag_len = n - p
        matches = sum(1 for a, b in zip(s[:-p], s[p:], strict=False) if a == b)
        frac = matches / lag_len
        if frac >= min_identity:
            return True, p, round(frac, 4)
    return False, 0, 0.0


def extract_motif_n_from_meme(
    meme_txt_path: Path | str, motif_num: int, output_aux_meme: Path | str
) -> dict[str, Any] | None:
    """Extract secondary or tertiary motif from multi-motif MEME text output."""
    p_in = Path(meme_txt_path)
    if not p_in.exists():
        return None
    text = p_in.read_text(encoding="utf-8")
    m_split = re.split(rf"\n(?=MOTIF\s+{motif_num}\b|\bMEME-{motif_num}\b)", text)
    if len(m_split) < 2:
        return None

    header = text.split("MOTIF")[0].strip() + "\n\n"
    next_num = motif_num + 1
    motif_content = re.split(
        rf"\n(?=MOTIF\s+{next_num}\b|\bMEME-{next_num}\b)", m_split[1]
    )[0]
    out_p = Path(output_aux_meme)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text(header + motif_content, encoding="utf-8")

    parsed = parse_aux_meme_results(out_p)
    if parsed:
        parsed["id"] = f"MEME_motif_{motif_num}"
        parsed["meme_txt_path"] = str(out_p)
        return parsed

    m_e = re.search(r"E=\s*(\S+)", motif_content)
    m_w = re.search(r"w=\s*(\d+)", motif_content)
    m_ns = re.search(r"nsites=\s*(\d+)", motif_content)
    cons_m = re.search(r"Multilevel\s+([A-Za-z]+)", motif_content)
    w = int(m_w.group(1)) if m_w else 8
    return {
        "id": f"MEME_motif_{motif_num}",
        "consensus": cons_m.group(1).upper() if cons_m else "",
        "width": w,
        "nsites": int(m_ns.group(1)) if m_ns else 0,
        "evalue": float(m_e.group(1)) if m_e else 1.0,
        "total_ic": w * 1.5,
        "mean_ic": 1.5,
        "consensus_entropy": 1.0,
        "meme_txt_path": str(out_p),
    }


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


def score_auxiliary_candidate(
    cand: dict[str, Any] | None,
    total_sequences: int = 1,
    telomere_repeat: str = SSOT["telomere_repeat"],
) -> float:
    """Compute composite score for an auxiliary motif candidate without proximity bias."""
    if not cand:
        return 0.0
    is_telo, telo_frac = is_telomeric_repeat_echo(
        cand.get("consensus", ""), repeat=telomere_repeat
    )
    cand["telomeric_echo_fraction"] = telo_frac
    if is_telo:
        cand["rejection_flag"] = "REJECTED_TELOMERIC_ECHO"
        cand["composite_score"] = 0.0
        return 0.0
    is_tr, _tr_p, _tr_frac = detect_tandem_repeats(cand.get("consensus", ""))
    if is_tr:
        cand["rejection_flag"] = "REJECTED_TANDEM_REPEAT"
        cand["composite_score"] = 0.0
        return 0.0
    if cand.get("evalue", 1.0) > SSOT["aux_evalue_cutoff"]:
        cand["rejection_flag"] = "REJECTED_INSIGNIFICANT_EVALUE"
        return 0.0
    cand["rejection_flag"] = "PASS"
    total_ic = cand.get("total_ic", 0.0)
    if total_ic <= 0.0:
        w = cand.get("width", 0)
        mean_ic = cand.get("mean_ic", 0.0)
        total_ic = w * mean_ic
    if total_ic <= 0.0:
        return 0.0

    unique_seqs = cand.get("unique_sequences_hit", 0)
    nsites = max(1, cand.get("nsites", 1))
    c_ent = cand.get("consensus_entropy", 1.0)

    tot_seq = max(1, total_sequences)
    occupancy = (
        min(1.0, unique_seqs / tot_seq)
        if unique_seqs > 0
        else min(1.0, nsites / tot_seq)
    )
    homopolymer_factor = (
        1.0 if c_ent >= 1.2 else (c_ent / 1.2 if c_ent >= 0.6 else 0.05)
    )

    score = total_ic * (0.3 + 0.7 * occupancy) * homopolymer_factor
    return round(score, 4)


def run_localized_auxiliary_meme(
    sliced_fasta: str | Path,
    out_dir: str | Path,
    bfile: str | None = None,
    nmotifs: int = 3,
    minw: int = 6,
    maxw: int = 30,
    mod: str = "zoops",
    revcomp: bool = True,
    core_consensus: str = "",
    total_sequences: int = 1,
    telomere_repeat: str = SSOT["telomere_repeat"],
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Execute single targeted localized MEME command on sliced FASTA and evaluate auxiliary candidate motifs."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    meme_txt = out_path / "meme.txt"

    if not meme_txt.exists():
        meme_bin = shutil.which("meme")
        if not meme_bin:
            raise RuntimeError("Required binary 'meme' not found in PATH")
        cmd = [
            meme_bin,
            str(sliced_fasta),
            "-dna",
            "-mod",
            mod,
            "-nmotifs",
            str(nmotifs),
            "-minw",
            str(minw),
            "-maxw",
            str(maxw),
            "-oc",
            str(out_path),
        ]
        if revcomp:
            cmd.append("-revcomp")
        if bfile and os.path.exists(bfile):
            cmd.extend(["-bfile", str(bfile)])

        try:
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr.strip() if e.stderr else "Unknown error"
            raise RuntimeError(
                f"Localized auxiliary MEME discovery failed with exit code {e.returncode}: {err_msg}"
            ) from e

    candidates = []
    for m_idx in range(2, nmotifs + 1):
        sub_meme_path = out_path / f"motif_{m_idx}.meme"
        cand = extract_motif_n_from_meme(meme_txt, m_idx, sub_meme_path)
        if not cand:
            continue

        cons = cand.get("consensus", "")
        if core_consensus and is_core_echo(cons, core_consensus):
            cand["rejection_flag"] = "REJECTED_CORE_ECHO"
            cand["composite_score"] = 0.0
            candidates.append(cand)
            continue

        is_telo, telo_frac = is_telomeric_repeat_echo(cons, repeat=telomere_repeat)
        cand["telomeric_echo_fraction"] = telo_frac
        if is_telo:
            cand["rejection_flag"] = "REJECTED_TELOMERIC_ECHO"
            cand["composite_score"] = 0.0
            candidates.append(cand)
            continue

        is_tr, _tr_p, _tr_frac = detect_tandem_repeats(cons)
        if is_tr:
            cand["rejection_flag"] = "REJECTED_TANDEM_REPEAT"
            cand["composite_score"] = 0.0
            candidates.append(cand)
            continue

        score = score_auxiliary_candidate(
            cand, total_sequences=total_sequences, telomere_repeat=telomere_repeat
        )
        cand["composite_score"] = score
        candidates.append(cand)

    valid_candidates = [
        c
        for c in candidates
        if c.get("rejection_flag") == "PASS" and c.get("composite_score", 0.0) > 0
    ]
    if valid_candidates:
        valid_candidates.sort(
            key=lambda x: (
                x.get("composite_score", 0.0),
                x.get("total_ic", 0.0),
                -math.log10(max(x.get("evalue", 1.0), 1e-100)),
            ),
            reverse=True,
        )
        return valid_candidates[0], candidates

    return None, candidates


def run_auxiliary_grid(
    masked_fasta: str | Path,
    out_dir: str | Path,
    bfile: str | None = None,
    threads: int = 1,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Execute auxiliary grid search across 16 conditions on core-masked FASTA."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    meme_bin = shutil.which("meme")
    if not meme_bin and not any(out_path.glob("*/meme.txt")):
        raise RuntimeError("Required binary 'meme' not found in PATH")

    n_seqs = 0
    with open(masked_fasta, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                n_seqs += 1
    total_sequences = max(1, n_seqs)

    tasks = []
    for cond in AUX_GRID_CONDITIONS:
        cid = cond["id"]
        cond_dir = out_path / cid
        cond_dir.mkdir(exist_ok=True)
        tasks.append((cond, str(masked_fasta), str(cond_dir), bfile))

    results = []
    max_workers = max(1, threads)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_run_single_aux_condition, meme_bin, cond, m_fa, c_dir, bf)
            for cond, m_fa, c_dir, bf in tasks
        ]
        for f in as_completed(futures):
            res = f.result()
            if res:
                res["composite_score"] = score_auxiliary_candidate(res, total_sequences)
                results.append(res)

    # Check for significant auxiliary motif ranked by composite score
    sig_candidates = [
        r
        for r in results
        if r["evalue"] <= 0.05
        and r["mean_ic"] >= 0.50
        and r["consensus_entropy"] >= 0.80
        and r.get("composite_score", 0.0) > 0
    ]
    if sig_candidates:
        sig_candidates.sort(
            key=lambda x: (
                x["composite_score"],
                x["total_ic"],
                -math.log10(max(x["evalue"], 1e-100)),
            ),
            reverse=True,
        )
        return sig_candidates[0], results

    return None, results
