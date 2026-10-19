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
    txt_path: str | Path, center_pos: int = 150
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
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    metrics = parse_aux_meme_results(str(meme_txt))
    if metrics:
        metrics["id"] = cid
        metrics["mod"] = cond["mod"]
        metrics["revcomp"] = cond["revcomp"]
        metrics["meme_txt_path"] = str(meme_txt)
        return metrics
    return None


def score_auxiliary_candidate(
    cand: dict[str, Any] | None, total_sequences: int = 1
) -> float:
    """Compute composite score for an auxiliary motif candidate."""
    if not cand or cand.get("evalue", 1.0) > 0.05:
        return 0.0
    total_ic = cand.get("total_ic", 0.0)
    if total_ic <= 0.0:
        w = cand.get("width", 0)
        mean_ic = cand.get("mean_ic", 0.0)
        total_ic = w * mean_ic
    if total_ic <= 0.0:
        return 0.0

    unique_seqs = cand.get("unique_sequences_hit", 0)
    nsites = max(1, cand.get("nsites", 1))
    proximal = cand.get("proximal_sites_6_50bp", 0)
    c_ent = cand.get("consensus_entropy", 1.0)

    tot_seq = max(1, total_sequences)
    occupancy = min(1.0, unique_seqs / tot_seq) if unique_seqs > 0 else min(1.0, nsites / tot_seq)
    spatial_spec = min(1.0, proximal / nsites) if proximal > 0 else 0.75
    homopolymer_factor = (
        1.0 if c_ent >= 1.2 else (c_ent / 1.2 if c_ent >= 0.6 else 0.05)
    )

    score = (
        total_ic
        * (0.3 + 0.7 * occupancy)
        * (0.4 + 0.6 * spatial_spec)
        * homopolymer_factor
    )
    return round(score, 4)


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
