"""Bipartite motif spacing, relative strand orientation, and cleavage offset profiler."""

import csv
from pathlib import Path
import statistics
from typing import Any


def calculate_edge_to_edge_spacer(
    c_start: int, c_stop: int, p_start: int, p_stop: int
) -> tuple[int, str]:
    """Calculate edge-to-edge spacer distance and location type."""
    if p_start >= c_start and p_stop <= c_stop:
        return 0, "Embedded"
    if p_stop < c_start:
        return c_start - p_stop, "Flank_Upstream"
    if p_start > c_stop:
        return p_start - c_stop, "Flank_Downstream"
    # Partial overlap
    return 0, "Partial_Overlap"


def profile_bipartite_spacing(
    seq_names: list[str],
    core_hits_by_seq: dict[str, list[dict[str, Any]]],
    aux_hits_by_seq: dict[str, list[dict[str, Any]]] | None = None,
    pos0_coord: int | dict[str, int] = 150,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Profile inter-motif edge-to-edge spacers, relative strand, and cleavage offsets."""
    rows = []
    spacers = []
    core_offsets = []
    partner_offsets = []

    aux_dict = aux_hits_by_seq or {}
    pos0_map = pos0_coord if isinstance(pos0_coord, dict) else {}
    default_pos0 = pos0_coord if isinstance(pos0_coord, int) else 150

    for s_name in seq_names:
        c_hits = core_hits_by_seq.get(s_name, [])
        a_hits = aux_dict.get(s_name, [])

        if not c_hits:
            continue

        pos0 = pos0_map.get(s_name, default_pos0)

        for c in c_hits:
            c_start = c["start"]
            c_stop = c["stop"]
            c_strand = c["strand"]
            # Distance from Core edge to Position 0
            c_dist_to_pos0 = (
                (pos0 - c_stop)
                if c_strand == "+"
                else (c_start - pos0)
            )
            core_offsets.append(c_dist_to_pos0)

            # Find partner: prefer auxiliary hits, fallback to second core hit
            partners = a_hits if a_hits else [h for h in c_hits if h != c]
            if not partners:
                rows.append({
                    "sequence_name": s_name,
                    "core_start": c_start,
                    "core_stop": c_stop,
                    "core_strand": c_strand,
                    "core_offset_to_pos0": c_dist_to_pos0,
                    "has_partner": False,
                    "partner_type": "None",
                    "partner_domain": "None",
                    "partner_start": "",
                    "partner_stop": "",
                    "partner_strand": "",
                    "partner_offset_to_pos0": "",
                    "spacer_bp": "",
                    "spacer_type": "None",
                    "relative_strand": "None",
                })
                continue

            # Pick nearest partner
            best_p = None
            min_d = 999999
            for p in partners:
                d, _typ = calculate_edge_to_edge_spacer(
                    c_start, c_stop, p["start"], p["stop"]
                )
                if d < min_d:
                    min_d = d
                    best_p = p

            if best_p is not None:
                p_start = best_p["start"]
                p_stop = best_p["stop"]
                p_strand = best_p["strand"]
                spacer_bp, spacer_type = calculate_edge_to_edge_spacer(
                    c_start, c_stop, p_start, p_stop
                )
                p_dist_to_pos0 = (
                    (pos0 - p_stop)
                    if p_strand == "+"
                    else (p_start - pos0)
                )
                p_domain = (
                    "Retained"
                    if p_stop <= pos0
                    else ("Eliminated" if p_start > pos0 else "Spanning")
                )
                rel_strand = "Same" if c_strand == p_strand else "Opposite"
                p_type = "Auxiliary" if a_hits else "Secondary_Core"

                spacers.append(spacer_bp)
                partner_offsets.append(p_dist_to_pos0)

                rows.append({
                    "sequence_name": s_name,
                    "core_start": c_start,
                    "core_stop": c_stop,
                    "core_strand": c_strand,
                    "core_offset_to_pos0": c_dist_to_pos0,
                    "has_partner": True,
                    "partner_type": p_type,
                    "partner_domain": p_domain,
                    "partner_start": p_start,
                    "partner_stop": p_stop,
                    "partner_strand": p_strand,
                    "partner_offset_to_pos0": p_dist_to_pos0,
                    "spacer_bp": spacer_bp,
                    "spacer_type": spacer_type,
                    "relative_strand": rel_strand,
                })

    # Summary statistics
    n_paired = len(spacers)
    if n_paired > 0:
        med_sp = float(statistics.median(spacers))
        sp_sorted = sorted(spacers)
        q25 = sp_sorted[len(sp_sorted) // 4]
        q75 = sp_sorted[(3 * len(sp_sorted)) // 4]
        iqr_sp = float(q75 - q25)
        pct_embedded = round(
            100.0 * sum(1 for s in spacers if s == 0) / n_paired, 1
        )
        pct_le30 = round(
            100.0 * sum(1 for s in spacers if s <= 30) / n_paired, 1
        )
        pct_le50 = round(
            100.0 * sum(1 for s in spacers if s <= 50) / n_paired, 1
        )
        same_count = sum(
            1 for r in rows if r.get("relative_strand") == "Same"
        )
        pct_same = round(100.0 * same_count / n_paired, 1)
        n_ret = sum(1 for r in rows if r.get("partner_domain") == "Retained")
        n_elim = sum(1 for r in rows if r.get("partner_domain") == "Eliminated")
        n_span = sum(1 for r in rows if r.get("partner_domain") == "Spanning")
        pct_retained = round(100.0 * n_ret / n_paired, 1)
        pct_eliminated = round(100.0 * n_elim / n_paired, 1)
        pct_spanning = round(100.0 * n_span / n_paired, 1)
    else:
        med_sp = 0.0
        iqr_sp = 0.0
        pct_embedded = 0.0
        pct_le30 = 0.0
        pct_le50 = 0.0
        pct_same = 0.0
        pct_retained = 0.0
        pct_eliminated = 0.0
        pct_spanning = 0.0

    med_core_off = (
        float(statistics.median(core_offsets)) if core_offsets else 0.0
    )
    med_partner_off = (
        float(statistics.median(partner_offsets)) if partner_offsets else 0.0
    )

    summary = {
        "n_paired_occurrences": n_paired,
        "median_spacer_bp": med_sp,
        "iqr_spacer_bp": iqr_sp,
        "pct_embedded_0bp": pct_embedded,
        "pct_spacer_le30bp": pct_le30,
        "pct_spacer_le50bp": pct_le50,
        "pct_same_strand": pct_same,
        "pct_retained_domain": pct_retained,
        "pct_eliminated_domain": pct_eliminated,
        "pct_spanning_domain": pct_spanning,
        "median_core_offset_to_pos0": med_core_off,
        "median_partner_offset_to_pos0": med_partner_off,
    }

    return rows, summary


def write_spacing_tsv(rows: list[dict[str, Any]], output_tsv: str | Path):
    """Write bipartite spacing details to TSV."""
    fieldnames = [
        "sequence_name",
        "core_start",
        "core_stop",
        "core_strand",
        "core_offset_to_pos0",
        "has_partner",
        "partner_type",
        "partner_domain",
        "partner_start",
        "partner_stop",
        "partner_strand",
        "partner_offset_to_pos0",
        "spacer_bp",
        "spacer_type",
        "relative_strand",
    ]
    with open(output_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
