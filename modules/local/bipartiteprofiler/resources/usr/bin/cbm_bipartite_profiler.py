#!/usr/bin/env python3
"""Stage 5: Soft-Erasure Multi-Motif Discovery, SpaMo Spacing Validation & Bipartite Profiler."""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from bipartite.aux_grid import (  # noqa: E402
    parse_aux_meme_results,
    run_auxiliary_grid,
    score_auxiliary_candidate,
)
from bipartite.masker import (  # noqa: E402
    load_fasta_dict,
    load_soma_len_map,
    mask_core_motifs,
    run_fimo,
)
from bipartite.reporter import generate_cbm_discovery_report  # noqa: E402
from bipartite.spacing import (  # noqa: E402
    profile_bipartite_spacing,
    write_spacing_tsv,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Soft-erasure multi-motif discovery, SpaMo spacing validation, and bipartite profiler"
    )
    parser.add_argument("--break-windows", required=True, help="Path to break_windows_300bp.fasta")
    parser.add_argument("--core-meme", required=True, help="Path to best_core_motif.meme")
    parser.add_argument("--census-json", required=True, help="Path to coverage_and_chromosome_census.json")
    parser.add_argument("--core-metrics-json", required=True, help="Path to core_motif_metrics.json")
    parser.add_argument("--background-bfile", default=None, help="Winning Markov background from Stage 4")
    parser.add_argument("--sample-id", default="sample", help="Sample ID")
    parser.add_argument("--output-aux-meme", default="best_auxiliary_motif.meme", help="Output aux MEME")
    parser.add_argument("--output-spacing-tsv", default="bipartite_spacing_summary.tsv", help="Spacing TSV")
    parser.add_argument("--output-report", default="cbm_discovery_report.md", help="Output markdown report")
    parser.add_argument("--output-masked-fasta", default="break_sequences.core_masked.fasta", help="Masked FASTA")
    parser.add_argument("--threads", type=int, default=1, help="Threads for MEME/FIMO/SpaMo")
    return parser.parse_args()


def build_higher_order_markov_bg(fasta_file: str | Path, out_bg: Path) -> Path:
    """Build higher-order Markov background from distal flanking arms [-150..-60] & [+60..+150]."""
    fgm_bin = shutil.which("fasta-get-markov")
    if not fgm_bin:
        raise RuntimeError("Required binary 'fasta-get-markov' not found in PATH")

    seqs = load_fasta_dict(fasta_file)
    distal_fa = out_bg.parent / "distal_arms_bg.fasta"

    total_len = 0
    with open(distal_fa, "w", encoding="utf-8") as f_out:
        for idx, (_hdr, s) in enumerate(seqs.items()):
            s = s.upper()
            left_arm = s[:90] if len(s) >= 90 else s
            right_arm = s[210:] if len(s) >= 210 else ""
            if left_arm:
                f_out.write(f">distal_L_{idx}\n{left_arm}\n")
                total_len += len(left_arm)
            if right_arm:
                f_out.write(f">distal_R_{idx}\n{right_arm}\n")
                total_len += len(right_arm)

    order = 3 if total_len >= 3600 else (2 if total_len >= 1440 else 1)
    subprocess.run(
        [fgm_bin, "-dna", "-m", str(order), str(distal_fa), str(out_bg)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    return out_bg


def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def is_core_echo(aux_consensus: str, core_consensus: str) -> bool:
    """Return True if aux_consensus is a substring/superstring of core_consensus or its revcomp."""
    if not aux_consensus or not core_consensus:
        return False
    a = aux_consensus.upper()
    c = core_consensus.upper()
    rc_c = revcomp(c)
    return (a in c) or (c in a) or (a in rc_c) or (rc_c in a)


def extract_motif_n_from_meme(
    meme_txt_path: Path, motif_num: int, output_aux_meme: Path
) -> dict[str, Any] | None:
    """Extract secondary or tertiary motif from multi-motif MEME text output."""
    if not meme_txt_path.exists():
        return None
    text = meme_txt_path.read_text(encoding="utf-8")
    m_split = re.split(rf"\n(?=MOTIF\s+{motif_num}\b|\bMEME-{motif_num}\b)", text)
    if len(m_split) < 2:
        return None

    header = text.split("MOTIF")[0].strip() + "\n\n"
    next_num = motif_num + 1
    motif_content = re.split(
        rf"\n(?=MOTIF\s+{next_num}\b|\bMEME-{next_num}\b)", m_split[1]
    )[0]
    output_aux_meme.write_text(header + motif_content, encoding="utf-8")

    parsed = parse_aux_meme_results(output_aux_meme)
    if parsed:
        parsed["id"] = f"MEME_soft_erasure_motif_{motif_num}"
        parsed["meme_txt_path"] = str(output_aux_meme)
        return parsed

    m_e = re.search(r"E=\s*(\S+)", motif_content)
    m_w = re.search(r"w=\s*(\d+)", motif_content)
    m_ns = re.search(r"nsites=\s*(\d+)", motif_content)
    cons_m = re.search(r"Multilevel\s+([A-Za-z]+)", motif_content)
    w = int(m_w.group(1)) if m_w else 8
    return {
        "id": f"MEME_soft_erasure_motif_{motif_num}",
        "consensus": cons_m.group(1).upper() if cons_m else "",
        "width": w,
        "nsites": int(m_ns.group(1)) if m_ns else 0,
        "evalue": float(m_e.group(1)) if m_e else 1.0,
        "total_ic": w * 1.5,
        "mean_ic": 1.5,
        "consensus_entropy": 1.0,
        "meme_txt_path": str(output_aux_meme),
    }


def main():
    args = parse_args()
    # Check all required binaries fail loudly
    for b in ["meme", "fasta-get-markov", "fimo", "spamo"]:
        if not shutil.which(b):
            raise RuntimeError(f"Required binary '{b}' not found in PATH")

    raw_seqs = load_fasta_dict(args.break_windows)
    seq_names = list(raw_seqs.keys())

    work_dir = Path("bipartite_work")
    work_dir.mkdir(exist_ok=True)

    # 1. Higher-order Markov background
    if args.background_bfile and os.path.exists(args.background_bfile):
        bg_file = Path(args.background_bfile)
    else:
        bg_file = build_higher_order_markov_bg(
            args.break_windows, work_dir / "bg_bipartite.txt"
        )

    # 2. Native Multi-Motif Soft-Erasure Discovery (MEME -nmotifs 3)
    meme_out_dir = work_dir / "meme_soft_erasure_out"
    subprocess.run(
        [
            "meme",
            args.break_windows,
            "-dna",
            "-mod",
            "zoops",
            "-nmotifs",
            "3",
            "-minw",
            "6",
            "-maxw",
            "25",
            "-revcomp",
            "-bfile",
            str(bg_file),
            "-oc",
            str(meme_out_dir),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )

    # 3. Spatial Spacing Validation via SpaMo
    meme_xml = meme_out_dir / "meme.xml"
    spamo_out_tsv = work_dir / "spamo.tsv"
    if meme_xml.exists():
        with open(spamo_out_tsv, "w", encoding="utf-8") as f_spamo:
            subprocess.run(
                [
                    "spamo",
                    "-text",
                    "-margin",
                    "100",
                    "-range",
                    "100",
                    args.break_windows,
                    str(meme_xml),
                    str(meme_xml),
                ],
                stdout=f_spamo,
                stderr=subprocess.DEVNULL,
                check=False,
            )

    # 4. Mask core motifs with N and run masked grid as complementary auxiliary search
    _masked_seqs, core_hits = mask_core_motifs(
        core_meme=args.core_meme,
        fasta_file=args.break_windows,
        output_masked_fasta=args.output_masked_fasta,
        fimo_out_dir="fimo_core_hits",
        threshold=1e-4,
    )

    core_hits_by_seq = defaultdict(list)
    for h in core_hits:
        core_hits_by_seq[h["seq_name"]].append(h)

    best_aux_grid, _grid_results = run_auxiliary_grid(
        masked_fasta=args.output_masked_fasta,
        out_dir="aux_grid_runs",
        bfile=str(bg_file),
        threads=args.threads,
    )

    core_consensus = ""
    if args.core_metrics_json and os.path.exists(args.core_metrics_json):
        try:
            with open(args.core_metrics_json, encoding="utf-8") as f_cm:
                core_consensus = json.load(f_cm).get("consensus", "")
        except Exception:
            pass

    # Compare soft-erasure Motif 2 / Motif 3 vs masked-grid auxiliary motif
    candidate_soft = None
    m2 = extract_motif_n_from_meme(
        meme_out_dir / "meme.txt", 2, work_dir / "soft_motif2.meme"
    )
    if m2 and not is_core_echo(m2["consensus"], core_consensus):
        candidate_soft = m2
    else:
        m3 = extract_motif_n_from_meme(
            meme_out_dir / "meme.txt", 3, work_dir / "soft_motif3.meme"
        )
        if m3 and not is_core_echo(m3["consensus"], core_consensus):
            candidate_soft = m3

    best_aux = None
    soft_score = (
        score_auxiliary_candidate(candidate_soft, len(seq_names))
        if candidate_soft
        else 0.0
    )
    grid_score = (
        score_auxiliary_candidate(best_aux_grid, len(seq_names))
        if best_aux_grid
        else 0.0
    )

    if candidate_soft and soft_score > 0 and soft_score >= grid_score:
        best_aux = candidate_soft
        best_aux["composite_score"] = soft_score
        shutil.copy(candidate_soft["meme_txt_path"], args.output_aux_meme)
    elif best_aux_grid and grid_score > 0 and best_aux_grid.get("meme_txt_path"):
        best_aux = best_aux_grid
        best_aux["composite_score"] = grid_score
        shutil.copy(best_aux["meme_txt_path"], args.output_aux_meme)
    else:
        with open(args.output_aux_meme, "w", encoding="utf-8") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n")

    # 5. Profile bipartite spacing
    aux_hits_by_seq = defaultdict(list)
    if best_aux and os.path.exists(args.output_aux_meme):
        aux_fimo_hits = run_fimo(
            meme_file=args.output_aux_meme,
            fasta_file=args.break_windows,
            out_dir="fimo_aux_hits",
            threshold=1e-3,
        )
        for h in aux_fimo_hits:
            aux_hits_by_seq[h["seq_name"]].append(h)

    soma_lens = load_soma_len_map(args.break_windows)
    spacing_rows, spacing_summary = profile_bipartite_spacing(
        seq_names=seq_names,
        core_hits_by_seq=core_hits_by_seq,
        aux_hits_by_seq=aux_hits_by_seq if best_aux else None,
        pos0_coord=soma_lens,
    )
    write_spacing_tsv(spacing_rows, args.output_spacing_tsv)

    # 6. Generate CBM Discovery Markdown Report
    generate_cbm_discovery_report(
        sample_id=args.sample_id,
        census_json_path=args.census_json,
        core_metrics_json_path=args.core_metrics_json,
        best_aux_motif=best_aux,
        bipartite_summary=spacing_summary,
        output_report_md=args.output_report,
    )

    print(
        f"[BIPARTITEPROFILER] Completed for {args.sample_id}. "
        f"Auxiliary motif: {'Found (' + best_aux['consensus'] + ')' if best_aux else 'None'}. "
        f"Paired loci: {spacing_summary['n_paired_occurrences']}."
    )


if __name__ == "__main__":
    main()
