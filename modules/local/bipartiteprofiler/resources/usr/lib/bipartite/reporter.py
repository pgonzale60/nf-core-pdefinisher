"""Markdown report generator for CBM Discovery pipeline."""

import json
import os
from pathlib import Path
from typing import Any


def generate_cbm_discovery_report(
    sample_id: str,
    census_json_path: str | Path,
    core_metrics_json_path: str | Path,
    best_aux_motif: dict[str, Any] | None,
    bipartite_summary: dict[str, Any],
    output_report_md: str | Path,
):
    """Generate comprehensive cbm_discovery_report.md markdown summary."""
    census = {}
    if os.path.exists(census_json_path):
        with open(census_json_path, encoding="utf-8") as f:
            census = json.load(f)

    core_metrics = {}
    if os.path.exists(core_metrics_json_path):
        with open(core_metrics_json_path, encoding="utf-8") as f:
            core_metrics = json.load(f)

    md = [
        f"# CBM Discovery & PDE Junction Report: `{sample_id}`\n",
        "## 1. Executive Summary\n",
        f"- **Sample ID**: `{sample_id}`",
        f"- **PDE Breakpoints Identified ($N_{{breaks}}$)**: `{census.get('n_breaks', 'N/A')}`",
        f"- **Germline Native Telomeres ($N_{{germ_telo}}$)**: `{census.get('n_germ_telo', 'N/A')}`",
        f"- **Estimated Germline Chromosome Count**: `{census.get('estimated_germline_chromosomes', 'N/A')}`",
        f"- **Median Somatic Coverage ($\\tilde{{C}}_{{soma}}$)**: `{census.get('median_c_soma', 'N/A')}×`",
        f"- **Median Germline Coverage ($\\tilde{{C}}_{{germ}}$)**: `{census.get('median_c_germ', 'N/A')}×`",
        f"- **Estimated Germline Fraction**: `{round(float(census.get('germline_fraction', 0)) * 100, 1)}%`",
        "\n---\n",
        "## 2. Tier-1 Core Motif Architecture\n",
    ]

    if core_metrics:
        md.extend([
            f"- **Core Consensus**: `{core_metrics.get('consensus', 'N/A')}`",
            f"- **Motif Width**: `{core_metrics.get('width', 'N/A')} bp`",
            f"- **Mean Information Content**: `{core_metrics.get('mean_ic', 'N/A')} bits/pos`",
            f"- **E-value**: `{core_metrics.get('evalue_str', core_metrics.get('evalue', 'N/A'))}`",
            f"- **Break Site Occupancy**: `{round(float(core_metrics.get('occupancy', 0)) * 100, 1)}%` ({core_metrics.get('unique_occupied_sites', 0)} sites)",
            f"- **Median Cleavage Offset (to Pos 0)**: `{core_metrics.get('median_offset_to_pos0', 'N/A')} bp` (IQR: `{core_metrics.get('positional_iqr', 'N/A')} bp`)",
            f"- **Grid Search Condition**: `{core_metrics.get('condition_id', 'N/A')}`",
        ])
    else:
        md.append("- No statistically significant core motif detected.\n")

    md.extend(["\n---\n", "## 3. Tier-2 Auxiliary & Bipartite Architecture\n"])

    if best_aux_motif:
        md.extend([
            "### Tier-2 Auxiliary Motif Detected\n",
            f"- **Auxiliary Consensus**: `{best_aux_motif.get('consensus', 'N/A')}`",
            f"- **Width**: `{best_aux_motif.get('width', 'N/A')} bp`",
            f"- **E-value**: `{best_aux_motif.get('evalue_str', best_aux_motif.get('evalue', 'N/A'))}`",
            f"- **Mean Information Content**: `{best_aux_motif.get('mean_ic', 'N/A')} bits/pos`",
            f"- **Proximal Sites (6–50 bp from cleavage)**: `{best_aux_motif.get('proximal_sites_6_50bp', 0)}`",
            "\n### Bipartite Spacing & Orientation Metrics\n",
            f"- **Paired Core-Auxiliary Loci**: `{bipartite_summary.get('n_paired_occurrences', 0)}`",
            f"- **Median Edge-to-Edge Spacer**: `{bipartite_summary.get('median_spacer_bp', 0)} bp` (IQR: `{bipartite_summary.get('iqr_spacer_bp', 0)} bp`)",
            f"- **Embedded Occurrences (0 bp)**: `{bipartite_summary.get('pct_embedded_0bp', 0)}%`",
            f"- **Tightly Spaced (≤ 30 bp)**: `{bipartite_summary.get('pct_spacer_le30bp', 0)}%`",
            f"- **Flanking (≤ 50 bp)**: `{bipartite_summary.get('pct_spacer_le50bp', 0)}%`",
            f"- **Strand Concordance (Same Strand)**: `{bipartite_summary.get('pct_same_strand', 0)}%`",
            f"- **Located in Retained Region (≤ Pos 0)**: `{bipartite_summary.get('pct_retained_domain', 0.0)}%`",
            f"- **Located in Eliminated Region (> Pos 0)**: `{bipartite_summary.get('pct_eliminated_domain', 0.0)}%`",
        ])
        if bipartite_summary.get("pct_spanning_domain", 0.0) > 0:
            md.append(f"- **Spanning Cleavage Junction**: `{bipartite_summary.get('pct_spanning_domain')}%`")
    else:
        md.extend([
            "### Simple (Single Motif) Architecture\n",
            "- No statistically significant auxiliary motif was detected after core masking.",
            "- The cleavage system operates primarily via the single core recognition element.",
        ])

    md.append("\n---\n")

    with open(output_report_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
