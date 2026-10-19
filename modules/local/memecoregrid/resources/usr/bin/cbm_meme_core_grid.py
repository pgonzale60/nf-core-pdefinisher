#!/usr/bin/env python3
"""Stage 4: Tier-1 Core Motif Parameter Grid Search."""

import argparse
from pathlib import Path
import sys

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from memegrid.grid_runner import (  # noqa: E402
    export_grid_outputs,
    run_meme_grid,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run multi-radius MEME parameter grid search to select best core motif"
    )
    parser.add_argument(
        "--break-windows",
        required=True,
        help="Path to break_windows_300bp.fasta",
    )
    parser.add_argument(
        "--output-dir", default="meme_grid_work", help="Working directory for grid runs"
    )
    parser.add_argument(
        "--output-meme",
        default="best_core_motif.meme",
        help="Output best core motif MEME file",
    )
    parser.add_argument(
        "--output-tsv",
        default="core_grid_summary.tsv",
        help="Output TSV summary of grid search conditions",
    )
    parser.add_argument(
        "--output-metrics",
        default="core_motif_metrics.json",
        help="Output JSON of best core motif metrics",
    )
    parser.add_argument(
        "--output-background",
        default="best_background.bfile",
        help="Output winning Markov background file",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads for MEME"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    best_motif, grid_results = run_meme_grid(
        break_windows_fasta=args.break_windows,
        output_dir=args.output_dir,
        threads=args.threads,
    )

    export_grid_outputs(
        best_motif=best_motif,
        grid_results=grid_results,
        output_meme=args.output_meme,
        summary_tsv=args.output_tsv,
        metrics_json=args.output_metrics,
        output_background=args.output_background,
    )

    if best_motif:
        print(
            f"[MEMECOREGRID] Selected best core motif: Condition={best_motif['condition_id']} "
            f"Consensus={best_motif['consensus']} Width={best_motif['width']} "
            f"IC={best_motif['mean_ic']} E-value={best_motif['evalue']} "
            f"Occupancy={best_motif['occupancy']} Score={best_motif['composite_score']}"
        )
    else:
        print("[MEMECOREGRID] Warning: No valid motif identified in grid search.")


if __name__ == "__main__":
    main()
