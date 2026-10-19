"""MEME parameter grid search runner with distal flank Markov background and 3-phase search."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any

from .evaluator import parse_meme_txt
from .slicer import load_fasta, write_sliced_fasta

RADII = [30, 50, 100, 150]

WIDTH_RANGES = [
    {"label": "w6_15", "args": ["-minw", "6", "-maxw", "15"]},
    {"label": "w12_30", "args": ["-minw", "12", "-maxw", "30"]},
    {"label": "w15_50", "args": ["-minw", "15", "-maxw", "50"]},
]

MODELS = ["zoops", "anr", "oops"]


def generate_grid_conditions(radii: list[int] | None = None) -> list[dict[str, Any]]:
    """Generate parameter combinations across radius, model, strand, and width ranges."""
    selected_radii = radii or RADII
    conditions = []
    idx = 1
    for r in selected_radii:
        for mod in MODELS:
            for revcomp in [True, False]:
                for w in WIDTH_RANGES:
                    strand_lbl = "revcomp" if revcomp else "fwd"
                    conditions.append({
                        "id": f"R{r}_{mod}_{strand_lbl}_{w['label']}",
                        "cond_idx": idx, "radius": r, "mod": mod,
                        "revcomp": revcomp, "w_args": w["args"], "w_label": w["label"],
                    })
                    idx += 1
    return conditions


def build_distal_flank_markov_backgrounds(
    break_windows_fasta: str | Path,
    output_dir: Path,
) -> tuple[list[int], dict[int, str]]:
    """Build candidate Markov backgrounds (m=0,1,2,3) from distal flanking arms [-150..-60] & [+60..+150]."""
    fgm_bin = shutil.which("fasta-get-markov")
    if not fgm_bin:
        raise RuntimeError("Required binary 'fasta-get-markov' not found in PATH")

    records = load_fasta(break_windows_fasta)
    distal_fa = output_dir / "distal_flanks_bg.fasta"

    total_len = 0
    with open(distal_fa, "w", encoding="utf-8") as f_out:
        for idx, (_hdr, seq) in enumerate(records):
            s = seq.upper()
            # Distal 5' arm [-150 .. -60] (first 90 bp)
            left_arm = s[:90] if len(s) >= 90 else s
            # Distal 3' arm [+60 .. +150] (last 90 bp from 210 to 300)
            right_arm = s[210:] if len(s) >= 210 else ""
            if left_arm:
                f_out.write(f">distal_L_{idx}\n{left_arm}\n")
                total_len += len(left_arm)
            if right_arm:
                f_out.write(f">distal_R_{idx}\n{right_arm}\n")
                total_len += len(right_arm)

    if total_len >= 3600:
        candidate_orders = [2, 3]
    elif total_len >= 1500:
        candidate_orders = [1, 2]
    else:
        candidate_orders = [0, 1]

    bg_files: dict[int, str] = {}
    for m in candidate_orders:
        bg_path = output_dir / f"bg_m{m}.txt"
        subprocess.run(
            [fgm_bin, "-dna", "-m", str(m), str(distal_fa), str(bg_path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        bg_files[m] = str(bg_path)

    return candidate_orders, bg_files


def _run_single_condition(
    meme_bin: str,
    cond: dict[str, Any],
    r_fa_str: str,
    cond_dir_str: str,
    total_seqs: int,
) -> dict[str, Any] | None:
    """Run MEME for a single condition and parse metrics."""
    cond_dir = Path(cond_dir_str)
    meme_txt = cond_dir / "meme.txt"
    r = cond["radius"]
    cid = cond["id"]

    if not meme_txt.exists():
        cmd = [
            meme_bin,
            r_fa_str,
            "-dna",
            "-nmotifs",
            "1",
            "-mod",
            cond["mod"],
            "-oc",
            cond_dir_str,
        ] + cond["w_args"]
        if cond.get("bfile"):
            cmd.extend(["-bfile", cond["bfile"]])
        if cond.get("pal"):
            cmd.append("-pal")
        if cond["revcomp"]:
            cmd.append("-revcomp")

        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    if meme_txt.exists():
        metrics = parse_meme_txt(
            str(meme_txt), radius=r, total_sequences=total_seqs
        )
        if metrics:
            metrics["condition_id"] = cid
            metrics["phase"] = cond.get("phase", 1)
            metrics["radius"] = r
            metrics["markov_order"] = cond.get("markov_order", 0)
            metrics["bfile"] = cond.get("bfile", "")
            metrics["mod"] = cond["mod"]
            metrics["revcomp"] = cond["revcomp"]
            metrics["pal"] = cond.get("pal", False)
            metrics["w_label"] = cond["w_label"]
            metrics["meme_txt_path"] = str(meme_txt)
            return metrics
    return None


def run_meme_grid(
    break_windows_fasta: str | Path,
    output_dir: str | Path,
    radii: list[int] | None = None,
    threads: int = 1,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Execute 3-phase hierarchical MEME parameter grid search."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    records = load_fasta(break_windows_fasta)
    total_seqs = len(records)
    if total_seqs == 0:
        return None, []

    meme_bin = shutil.which("meme")
    if not meme_bin:
        raise RuntimeError("Required binary 'meme' not found in PATH")

    # Step 0: Markov backgrounds
    bg_dir = out_path / "backgrounds"
    bg_dir.mkdir(exist_ok=True)
    candidate_orders, bg_files = build_distal_flank_markov_backgrounds(
        break_windows_fasta, bg_dir
    )

    # Slice FASTA for each radius
    sliced_dir = out_path / "sliced_windows"
    sliced_dir.mkdir(exist_ok=True)
    selected_radii = radii or RADII
    sliced_fastas = {}
    for r in selected_radii:
        r_fa = sliced_dir / f"windows_radius_{r}bp.fasta"
        write_sliced_fasta(records, r, r_fa, center_soma_len=150)
        sliced_fastas[r] = r_fa

    runs_dir = out_path / "runs"
    runs_dir.mkdir(exist_ok=True)
    all_results: list[dict[str, Any]] = []

    def execute_batch(conditions: list[dict[str, Any]]) -> list[dict[str, Any]]:
        batch_tasks = []
        for cond in conditions:
            cid = cond["id"]
            cond_dir = runs_dir / cid
            cond_dir.mkdir(exist_ok=True)
            r = cond["radius"]
            r_fa = sliced_fastas[r]
            batch_tasks.append((cond, str(r_fa), str(cond_dir)))

        results = []
        with ProcessPoolExecutor(max_workers=max(1, threads)) as executor:
            futures = [
                executor.submit(
                    _run_single_condition,
                    meme_bin,
                    cond,
                    r_fa_str,
                    c_dir_str,
                    total_seqs,
                )
                for cond, r_fa_str, c_dir_str in batch_tasks
            ]
            for f in as_completed(futures):
                res = f.result()
                if res:
                    results.append(res)
        return results

    # Phase 1: Spatial Radius & Markov Order Localization
    phase1_conds = []
    for r in selected_radii:
        for m in candidate_orders:
            phase1_conds.append({
                "id": f"P1_R{r}_m{m}_zoops_revcomp_w8_30", "phase": 1,
                "radius": r, "markov_order": m, "bfile": bg_files[m],
                "mod": "zoops", "revcomp": True, "w_args": ["-minw", "8", "-maxw", "30"],
                "w_label": "w8_30",
            })
    p1_results = execute_batch(phase1_conds)
    all_results.extend(p1_results)

    if not p1_results:
        return None, []

    p1_best = max(p1_results, key=lambda x: x["composite_score"])
    best_radius = p1_best["radius"]
    best_order = p1_best["markov_order"]

    # Phase 2: Width Bracket & Strand Refinement
    phase2_conds = []
    for w in WIDTH_RANGES:
        for revcomp in [True, False]:
            st_lbl = "revcomp" if revcomp else "fwd"
            phase2_conds.append({
                "id": f"P2_R{best_radius}_m{best_order}_zoops_{st_lbl}_{w['label']}",
                "phase": 2, "radius": best_radius, "markov_order": best_order,
                "bfile": bg_files[best_order], "mod": "zoops", "revcomp": revcomp,
                "w_args": w["args"], "w_label": w["label"],
            })
    p2_results = execute_batch(phase2_conds)
    all_results.extend(p2_results)

    p2_best = max(
        p2_results if p2_results else p1_results,
        key=lambda x: x["composite_score"],
    )
    best_w_args = next(
        w["args"] for w in WIDTH_RANGES if w["label"] == p2_best["w_label"]
    )
    best_w_label = p2_best["w_label"]
    best_revcomp = p2_best["revcomp"]

    # Phase 3: Site Distribution Model & Palindrome Check
    phase3_conds = []
    st_lbl = "revcomp" if best_revcomp else "fwd"
    for mod in ["zoops", "oops", "anr"]:
        phase3_conds.append({
            "id": f"P3_R{best_radius}_m{best_order}_{mod}_{st_lbl}_{best_w_label}",
            "phase": 3, "radius": best_radius, "markov_order": best_order,
            "bfile": bg_files[best_order], "mod": mod, "revcomp": best_revcomp,
            "pal": False, "w_args": best_w_args, "w_label": best_w_label,
        })
    phase3_conds.append({
        "id": f"P3_R{best_radius}_m{best_order}_zoops_{st_lbl}_pal_{best_w_label}",
        "phase": 3, "radius": best_radius, "markov_order": best_order,
        "bfile": bg_files[best_order], "mod": "zoops", "revcomp": best_revcomp,
        "pal": True, "w_args": best_w_args, "w_label": best_w_label,
    })
    p3_results = execute_batch(phase3_conds)
    all_results.extend(p3_results)

    all_results.sort(key=lambda r: r["composite_score"], reverse=True)
    best_motif = all_results[0] if all_results else None

    # Save winning Markov background
    if best_motif and best_motif.get("bfile"):
        shutil.copy(best_motif["bfile"], out_path / "best_background.bfile")

    return best_motif, all_results


def export_grid_outputs(
    best_motif: dict[str, Any] | None,
    grid_results: list[dict[str, Any]],
    output_meme: str | Path,
    summary_tsv: str | Path,
    metrics_json: str | Path,
    output_background: str | Path | None = None,
):
    """Save best_core_motif.meme, core_grid_summary.tsv, core_motif_metrics.json, and best_background.bfile."""
    if best_motif and os.path.exists(best_motif["meme_txt_path"]):
        shutil.copy(best_motif["meme_txt_path"], output_meme)
    else:
        with open(output_meme, "w", encoding="utf-8") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n")

    if (
        best_motif
        and best_motif.get("bfile")
        and os.path.exists(best_motif["bfile"])
    ):
        target_bg = output_background or (
            Path(output_meme).parent / "best_background.bfile"
        )
        shutil.copy(best_motif["bfile"], target_bg)
    elif output_background:
        with open(output_background, "w", encoding="utf-8") as f:
            f.write("# 0-order Markov background\nA 0.25\nC 0.25\nG 0.25\nT 0.25\n")

    fieldnames = [
        "phase", "condition_id", "radius", "markov_order", "mod", "revcomp",
        "pal", "w_label", "width", "nsites", "evalue", "consensus", "total_ic",
        "mean_ic", "consensus_entropy", "occupancy", "median_offset_to_pos0",
        "positional_iqr", "composite_score",
    ]
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(grid_results)

    with open(metrics_json, "w", encoding="utf-8") as f:
        json.dump(best_motif or {}, f, indent=2)
