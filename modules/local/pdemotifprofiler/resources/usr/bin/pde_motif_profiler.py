#!/usr/bin/env python3
"""PDE Motif Genome-Wide Profiler & Spatial Distribution Analyzer.

Intersects genome-wide FIMO motif occurrences with:
  1. Breakpoint locations (high-confidence cleavage junctions)
  2. Putative Core vs Putative Eliminated (GRS) segments (from IDENTIFYPUTATIVEGRS)

Answers the fundamental question:
Is the discovered elimination motif strictly localized to break sites, or
broadly distributed across satellite cores / transposable element bodies?

Outputs:
  - <prefix>.fimo_annotated_hits.tsv: hit table annotated with genomic compartment
  - <prefix>.motif_distribution_summary.tsv: tabular metrics of compartment enrichment & periodicity
  - <prefix>.motif_distribution_report.md: human-readable markdown assessment
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_breakpoints(bed_path, flank_tolerance=500):
    """Load breakpoints and build tolerance intervals."""
    bps = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return bps
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else f"{chrom}_{start}"
                bps[chrom].append(
                    {
                        "start": start,
                        "end": end,
                        "name": name,
                        "zone_start": max(0, start - flank_tolerance),
                        "zone_end": end + flank_tolerance,
                    }
                )
    return bps


def parse_segments(bed_path):
    """Load PDE core vs GRS segmentation."""
    segs = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return segs
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                label = parts[3]
                segs[chrom].append({"start": start, "end": end, "label": label})
    return segs


def main():
    parser = argparse.ArgumentParser(description="Profile genome-wide FIMO motif occurrences across PDE compartments")
    parser.add_argument("--fimo-tsv", required=True, help="Path to fimo_genome_hits.tsv file")
    parser.add_argument(
        "--breakpoints-bed",
        required=True,
        help="Path to breakpoint_locations.bed file",
    )
    parser.add_argument(
        "--pde-segments-bed",
        required=True,
        help="Path to pde_segments.bed file from IDENTIFYPUTATIVEGRS",
    )
    parser.add_argument(
        "--flank-tolerance",
        type=int,
        default=500,
        help="Window around breakpoints to classify as 'Break Site' (default: 500 bp)",
    )
    parser.add_argument("--prefix", default="sample", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    bps = parse_breakpoints(args.breakpoints_bed, args.flank_tolerance)
    segs = parse_segments(args.pde_segments_bed)

    sys.stderr.write(
        f"Loaded {sum(len(v) for v in bps.values())} breakpoints and {sum(len(v) for v in segs.values())} segments.\n"
    )

    annotated_hits = []
    hits_per_chrom = defaultdict(list)

    if os.path.exists(args.fimo_tsv):
        with open(args.fimo_tsv, encoding="utf-8") as f:
            header_line = f.readline().strip()
            # If commented with '#'
            if header_line.startswith("#"):
                header_line = header_line.lstrip("#").strip()
            headers = header_line.split("\t")

            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < len(headers):
                    continue
                row = dict(zip(headers, parts))
                chrom = row.get("sequence_name", "")
                try:
                    start = int(row.get("start", 0))
                    stop = int(row.get("stop", 0))
                except ValueError:
                    continue

                # Determine compartment
                compartment = "Retained Genome (Core)"
                min_dist = float("inf")
                nearest_bp_id = "NA"

                # Check proximity to breakpoints
                if chrom in bps:
                    for bp in bps[chrom]:
                        dist = min(abs(start - bp["start"]), abs(stop - bp["start"]))
                        if dist < min_dist:
                            min_dist = dist
                            nearest_bp_id = bp["name"]

                        if start <= bp["zone_end"] and stop >= bp["zone_start"]:
                            compartment = "Break Site"

                # If not at break site, check PDE segment
                if compartment != "Break Site" and chrom in segs:
                    for s in segs[chrom]:
                        if start >= s["start"] and stop <= s["end"]:
                            if "GRS" in s["label"] or "eliminated" in s["label"].lower():
                                compartment = "Eliminated Region (GRS)"
                            else:
                                compartment = "Retained Genome (Core)"
                            break

                try:
                    score_val = float(row.get("score", 0.0))
                except ValueError:
                    score_val = 0.0

                row["genomic_compartment"] = compartment
                row["dist_to_nearest_breakpoint_bp"] = str(min_dist) if min_dist != float("inf") else "NA"
                row["nearest_breakpoint_id"] = nearest_bp_id

                annotated_hits.append(row)
                hits_per_chrom[chrom].append((start, stop, compartment, score_val))

    total_hits = len(annotated_hits)
    counts = defaultdict(int)
    for h in annotated_hits:
        counts[h["genomic_compartment"]] += 1

    break_count = counts["Break Site"]
    grs_count = counts["Eliminated Region (GRS)"]
    core_count = counts["Retained Genome (Core)"]

    break_pct = (break_count / total_hits * 100.0) if total_hits > 0 else 0.0
    grs_pct = (grs_count / total_hits * 100.0) if total_hits > 0 else 0.0
    core_pct = (core_count / total_hits * 100.0) if total_hits > 0 else 0.0

    # Spacing / periodicity analysis (satellite detection)
    inter_hit_distances = []
    for chrom, hit_list in hits_per_chrom.items():
        hit_list.sort(key=lambda x: x[0])
        for i in range(len(hit_list) - 1):
            d = hit_list[i + 1][0] - hit_list[i][0]
            if 0 < d <= 50000:
                inter_hit_distances.append(d)

    median_inter_distance = sorted(inter_hit_distances)[len(inter_hit_distances) // 2] if inter_hit_distances else 0

    # Diagnosis
    if total_hits == 0:
        diagnosis = "NO_MOTIF_HITS_DETECTED"
    elif break_pct >= 60.0:
        diagnosis = "STRICTLY_LOCALIZED_TO_BREAK_SITES"
    elif (break_pct + grs_pct) >= 70.0 and median_inter_distance > 0 and median_inter_distance <= 10000:
        diagnosis = "ENRICHED_AT_BREAKS_AND_GRS_SATELLITE_CORES"
    elif grs_pct >= 60.0:
        diagnosis = "ENRICHED_IN_ELIMINATED_GRS_REPEATS"
    else:
        diagnosis = "BROADLY_DISTRIBUTED_GENOME_WIDE"

    # Write annotated hits TSV
    out_annot = os.path.join(args.outdir, f"{args.prefix}.fimo_annotated_hits.tsv")
    if annotated_hits:
        fieldnames = list(annotated_hits[0].keys())
        with open(out_annot, "w", encoding="utf-8") as f:
            f.write("\t".join(fieldnames) + "\n")
            for h in annotated_hits:
                f.write("\t".join([str(h.get(k, "")) for k in fieldnames]) + "\n")
    else:
        with open(out_annot, "w", encoding="utf-8") as f:
            f.write(
                "motif_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\tgenomic_compartment\tdist_to_nearest_breakpoint_bp\tnearest_breakpoint_id\n"
            )

    # -------------------------------------------------------------------------
    # Border-Threshold Specificity Analysis:
    # 1. Identify the top-scoring motif hit within the flank tolerance of each border.
    # 2. Determine the minimum score among these best-per-border hits.
    # 3. Filter genome-wide hits for score >= min_border_score, and report how many
    #    occur in GRS vs Retained Core.
    # -------------------------------------------------------------------------
    border_top_scores = []
    for chrom, bp_list in bps.items():
        for bp in bp_list:
            b_s = bp["zone_start"]
            b_e = bp["zone_end"]
            matching_scores = []
            if chrom in hits_per_chrom:
                for h_start, h_stop, _, h_score in hits_per_chrom[chrom]:
                    if not (h_stop < b_s or h_start > b_e):
                        matching_scores.append(h_score)
            if matching_scores:
                border_top_scores.append(max(matching_scores))
            else:
                border_top_scores.append(None)

    valid_border_scores = [s for s in border_top_scores if s is not None]
    min_border_score = min(valid_border_scores) if valid_border_scores else 0.0

    # Filter genome-wide hits at or above min_border_score
    thresh_hits = [h for h in annotated_hits if float(h.get("score", 0.0)) >= min_border_score]
    thresh_counts = defaultdict(int)
    for h in thresh_hits:
        thresh_counts[h["genomic_compartment"]] += 1

    thresh_break_count = thresh_counts["Break Site"]
    thresh_grs_count = thresh_counts["Eliminated Region (GRS)"]
    thresh_core_count = thresh_counts["Retained Genome (Core)"]

    # Write summary TSV
    out_summary = os.path.join(args.outdir, f"{args.prefix}.motif_distribution_summary.tsv")
    with open(out_summary, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"sample_id\t{args.prefix}\n")
        f.write(f"total_borders\t{len(border_top_scores)}\n")
        f.write(f"lowest_score_at_border\t{min_border_score:.4f}\n")
        f.write(f"significant_hits_at_border\t{thresh_break_count}\n")
        f.write(f"significant_hits_in_grs\t{thresh_grs_count}\n")
        f.write(f"significant_hits_in_core\t{thresh_core_count}\n")
        f.write(f"significant_hits_total\t{len(thresh_hits)}\n")

    # Write Markdown Report
    out_report = os.path.join(args.outdir, f"{args.prefix}.motif_distribution_report.md")
    with open(out_report, "w", encoding="utf-8") as f:
        f.write(f"# PDE Motif Specificity & Border Threshold Report: {args.prefix}\n\n")
        f.write("## Border-Level Specificity (Lowest Top-Scoring Border Hit)\n\n")
        f.write(f"- **Total PDE Borders:** {len(border_top_scores)}\n")
        f.write(f"- **Lowest Score at Border (Threshold):** `{min_border_score:.4f}`\n\n")
        f.write(f"When filtering genome-wide for hits with **Score $\\ge$ {min_border_score:.4f}**:\n")
        f.write(f"- **Significant Hits at Border:** {thresh_break_count:,}\n")
        f.write(f"- **Significant Hits in Eliminated Region (GRS):** {thresh_grs_count:,}\n")
        f.write(f"- **Significant Hits in Retained Genome (Core):** {thresh_core_count:,}\n")
        f.write(f"- **Total Significant Hits Genome-Wide:** {len(thresh_hits):,}\n")

    sys.stderr.write(
        f"Border-level threshold: min_border_score={min_border_score:.2f} -> "
        f"Total={len(thresh_hits)}, Borders={thresh_break_count}, GRS={thresh_grs_count}, Core={thresh_core_count}\n"
    )


if __name__ == "__main__":
    main()
