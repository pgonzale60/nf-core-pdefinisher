#!/usr/bin/env python3
"""Extract PDE Border Flank Sequences directly from Classified Segments & Genome Assembly.

Identifies the exact boundaries where `putative_GRS` meets `putative_core` from pde_segments.bed,
orients them relative to cleavage:
  - Transition from GRS (upstream) to Core (downstream): 5' end of core -> strand '-'
  - Transition from Core (upstream) to GRS (downstream): 3' end of core -> strand '+'

Extracts uniform flanking windows (±150 bp, ±350 bp) directly from the genome assembly FASTA,
orienting each sequence so that the core/cleavage polarity is uniform.

Outputs:
  - <prefix>.pde_border_breakpoints.bed
  - <prefix>.pde_border_flanks_150bp.fasta
  - <prefix>.pde_border_flanks_350bp.fasta
  - <prefix>.pde_border_flanks.bed
"""

import argparse
import gzip
import json
import os
import re
import statistics
import sys


def revcomp(seq):
    """Return reverse complement of DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def load_fasta(fasta_path):
    """Load FASTA into dictionary of chromosome sequences."""
    genome = {}
    open_fn = gzip.open if fasta_path.endswith(".gz") else open
    curr_chrom = None
    chunks = []
    with open_fn(fasta_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if curr_chrom:
                    genome[curr_chrom] = "".join(chunks).upper()
                curr_chrom = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if curr_chrom:
            genome[curr_chrom] = "".join(chunks).upper()
    return genome


def parse_breakpoints(bed_path):
    """Parse candidate breakpoints to retain original names and read support scores."""
    bps = {}
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
                name = parts[3] if len(parts) > 3 else f"{chrom}_{start}"
                score = parts[4] if len(parts) > 4 else "0"
                strand = parts[5] if len(parts) > 5 else "+"
                bps.setdefault(chrom, []).append({
                    "start": start,
                    "end": int(parts[2]),
                    "name": name,
                    "score": score,
                    "strand": strand,
                })
    return bps


def parse_segments(bed_path):
    """Parse classified segments into ordered per-chromosome intervals."""
    segs = {}
    if not os.path.exists(bed_path):
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
                segs.setdefault(chrom, []).append((start, end, label))
    for chrom in segs:
        segs[chrom].sort(key=lambda x: x[0])
    return segs


def compute_census_metrics(
    frontier_breakpoints,
    sample_id,
    segs_by_chrom,
    grs_summary_tsv=None,
    terminal_ends_summary_tsv=None,
):
    """Compute comprehensive coverage, breakpoint, and chromosome census."""
    n_breaks = len(frontier_breakpoints)
    median_c_soma = None
    median_c_germ = None
    germline_fraction = None
    n_germ_telo = None
    estimated_germline_chromosomes = None

    # 1. Coverage Metrics from grs_summary.tsv
    if grs_summary_tsv and os.path.exists(grs_summary_tsv):
        grs_depths = []
        core_depths = []
        try:
            with open(grs_summary_tsv, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith("#"):
                        for part in line.split("\t"):
                            if part.startswith("SOMATIC_CORE_BASELINE="):
                                try:
                                    median_c_soma = float(part.split("=")[1])
                                except ValueError:
                                    pass
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 6 and parts[0] != "chromosome":
                        classification = parts[4]
                        try:
                            depth = float(parts[5])
                            if "GRS" in classification:
                                grs_depths.append(depth)
                            elif "core" in classification.lower():
                                core_depths.append(depth)
                        except ValueError:
                            pass

            if median_c_soma is None and core_depths:
                median_c_soma = statistics.median(core_depths)

            if grs_depths:
                median_c_germ = statistics.median(grs_depths)
            elif median_c_soma is not None:
                median_c_germ = 0.0

            if median_c_soma is not None and median_c_germ is not None:
                if median_c_soma > 0:
                    germline_fraction = round(median_c_germ / median_c_soma, 4)
                else:
                    germline_fraction = 0.0
        except Exception as e:
            sys.stderr.write(f"Warning: Failed to parse grs_summary_tsv {grs_summary_tsv}: {e}\n")

    # 2. Native Germline Telomeres from terminal_ends_summary.tsv and segs_by_chrom
    if terminal_ends_summary_tsv and os.path.exists(terminal_ends_summary_tsv):
        try:
            germline_telomeres = []
            with open(terminal_ends_summary_tsv, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("chrom_end"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 7:
                        chrom_end = parts[0]
                        chrom_len = int(parts[1]) if parts[1].isdigit() else 0
                        n_telo_reads = int(parts[3]) if parts[3].isdigit() else 0
                        prop_status = parts[5]
                        term_tract = int(parts[6]) if parts[6].isdigit() else 0

                        if "_" in chrom_end:
                            chrom, arm = chrom_end.rsplit("_", 1)
                        else:
                            continue

                        # Determine if this arm terminates at a native germline end or a diminution break
                        is_germline_arm = False
                        if chrom in segs_by_chrom and segs_by_chrom[chrom]:
                            chrom_segs = segs_by_chrom[chrom]
                            if arm == "5p":
                                tip_seg = chrom_segs[0]
                                if "GRS" in tip_seg[2] or len(chrom_segs) == 1:
                                    is_germline_arm = True
                                elif tip_seg[0] == 0 and not any("GRS" in s[2] for s in chrom_segs if s[0] < tip_seg[1]):
                                    is_germline_arm = True
                            elif arm == "3p":
                                tip_seg = chrom_segs[-1]
                                if "GRS" in tip_seg[2] or len(chrom_segs) == 1:
                                    is_germline_arm = True
                                elif tip_seg[1] >= chrom_len:
                                    is_germline_arm = True
                        else:
                            is_germline_arm = True

                        has_telomere = (
                            prop_status in ("ALREADY_CAPPED", "EXTEND_AND_CAP", "TELOMERE_IN_OVERHANG")
                            or term_tract > 0
                            or n_telo_reads >= 3
                        )
                        if is_germline_arm and has_telomere:
                            germline_telomeres.append(chrom_end)

            n_germ_telo = len(germline_telomeres)
            if n_germ_telo > 0:
                estimated_germline_chromosomes = (n_germ_telo + 1) // 2
        except Exception as e:
            sys.stderr.write(f"Warning: Failed to parse terminal_ends_summary_tsv {terminal_ends_summary_tsv}: {e}\n")

    return {
        "sample_id": sample_id,
        "n_breaks": n_breaks,
        "n_germ_telo": n_germ_telo if n_germ_telo is not None else "N/A",
        "estimated_germline_chromosomes": estimated_germline_chromosomes if estimated_germline_chromosomes is not None else "N/A",
        "median_c_soma": round(median_c_soma, 1) if median_c_soma is not None else "N/A",
        "median_c_germ": round(median_c_germ, 1) if median_c_germ is not None else "N/A",
        "germline_fraction": germline_fraction if germline_fraction is not None else None,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract PDE border sequences directly from classified segments and assembly FASTA"
    )
    parser.add_argument("--pde-segments-bed", required=True, help="Path to pde_segments.bed")
    parser.add_argument("--assembly", required=True, help="Path to genome FASTA assembly")
    parser.add_argument("--breakpoints-bed", default=None, help="Optional path to filtered_breakpoints.bed for metadata")
    parser.add_argument("--grs-summary", default=None, help="Optional path to grs_summary.tsv")
    parser.add_argument("--terminal-ends-summary", default=None, help="Optional path to terminal_ends_summary.tsv")
    parser.add_argument("--window-small", type=int, default=150, help="Small half-window in bp (default: 150)")
    parser.add_argument("--window-large", type=int, default=350, help="Large half-window in bp (default: 350)")
    parser.add_argument("--prefix", default="sample", help="Output prefix")
    parser.add_argument("--outdir", default=".", help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sys.stderr.write(f"Loading assembly from {args.assembly}...\n")
    genome = load_fasta(args.assembly)
    sys.stderr.write(f"Loaded {len(genome)} sequences from assembly.\n")

    segs_by_chrom = parse_segments(args.pde_segments_bed)
    raw_bps_by_chrom = parse_breakpoints(args.breakpoints_bed)

    # Derive GRS/Core frontier boundaries
    frontier_breakpoints = []
    for chrom, seg_list in segs_by_chrom.items():
        if len(seg_list) < 2:
            continue
        for i in range(len(seg_list) - 1):
            s1, e1, lab1 = seg_list[i]
            s2, e2, lab2 = seg_list[i + 1]
            if e1 != s2:
                continue

            is_transition = False
            strand = "+"
            if ("GRS" in lab1 or "eliminated" in lab1.lower()) and ("core" in lab2.lower()):
                is_transition = True
                strand = "-"  # Cleavage at 5' end of core
            elif ("core" in lab1.lower()) and ("GRS" in lab2 or "eliminated" in lab2.lower()):
                is_transition = True
                strand = "+"  # Cleavage at 3' end of core

            if is_transition:
                pos = e1
                # Check if matching candidate breakpoint exists within 100 bp
                matched_name = f"{chrom}_{pos}"
                matched_score = "100"
                if chrom in raw_bps_by_chrom:
                    candidates = [
                        bp for bp in raw_bps_by_chrom[chrom]
                        if abs(bp["start"] - pos) <= 100
                    ]
                    if candidates:
                        candidates.sort(key=lambda b: abs(b["start"] - pos))
                        matched_name = candidates[0]["name"]
                        matched_score = candidates[0]["score"]

                frontier_breakpoints.append({
                    "chrom": chrom,
                    "start": pos,
                    "end": pos + 1,
                    "name": matched_name,
                    "score": matched_score,
                    "strand": strand,
                })

    frontier_breakpoints.sort(key=lambda x: (x["chrom"], x["start"]))
    sys.stderr.write(
        f"Identified {len(frontier_breakpoints)} genuine GRS/Core frontier boundaries across genome.\n"
    )

    # 1. Output frontier breakpoints BED
    out_frontier_bed = os.path.join(args.outdir, f"{args.prefix}.pde_border_breakpoints.bed")
    with open(out_frontier_bed, "w", encoding="utf-8") as f:
        for bp in frontier_breakpoints:
            f.write(f"{bp['chrom']}\t{bp['start']}\t{bp['end']}\t{bp['name']}\t{bp['score']}\t{bp['strand']}\n")

    # 2. Extract and orient flank sequences
    windows = [
        ("150bp", args.window_small),
        ("350bp", args.window_large),
    ]

    bed_intervals = []
    for label, win in windows:
        out_fasta = os.path.join(args.outdir, f"{args.prefix}.pde_border_flanks_{label}.fasta")
        fasta_records = []

        for bp in frontier_breakpoints:
            chrom = bp["chrom"]
            if chrom not in genome:
                continue

            chr_len = len(genome[chrom])
            pos = bp["start"]
            strand = bp["strand"]

            flank_start = max(0, pos - win)
            flank_end = min(chr_len, pos + win)
            subseq = genome[chrom][flank_start:flank_end]

            clean_bp_name = re.sub(r"[^A-Za-z0-9_-]+", "_", bp["name"]).strip("_")
            strand_str = "minus" if strand == "-" else "plus"
            seq_id = f"{chrom}_{pos + 1}_{strand_str}_win{win}bp_oriented"
            header = f"{seq_id} soma_len={win}"

            if strand == "-":
                subseq_oriented = revcomp(subseq)
            else:
                subseq_oriented = subseq

            fasta_records.append(f">{header}\n{subseq_oriented}")

            if win == args.window_small:
                bed_intervals.append(f"{chrom}\t{flank_start}\t{flank_end}\t{header.split()[0]}\t{bp['score']}\t{strand}\n")

        with open(out_fasta, "w", encoding="utf-8") as f:
            if fasta_records:
                f.write("\n".join(fasta_records) + "\n")
        sys.stderr.write(f"Wrote {len(fasta_records)} oriented sequences to {out_fasta}\n")

        if win == 150 or label == "150bp":
            out_300bp = os.path.join(args.outdir, f"{args.prefix}.break_windows_300bp.fasta")
            with open(out_300bp, "w", encoding="utf-8") as f:
                if fasta_records:
                    f.write("\n".join(fasta_records) + "\n")

    out_bed_intervals = os.path.join(args.outdir, f"{args.prefix}.pde_border_flanks.bed")
    with open(out_bed_intervals, "w", encoding="utf-8") as f:
        f.writelines(bed_intervals)

    census_data = compute_census_metrics(
        frontier_breakpoints=frontier_breakpoints,
        sample_id=args.prefix,
        segs_by_chrom=segs_by_chrom,
        grs_summary_tsv=args.grs_summary,
        terminal_ends_summary_tsv=args.terminal_ends_summary,
    )
    out_census = os.path.join(args.outdir, f"{args.prefix}.census.json")
    with open(out_census, "w", encoding="utf-8") as f:
        json.dump(census_data, f, indent=2)

    sys.stderr.write("Border flank extraction completed successfully.\n")


if __name__ == "__main__":
    main()
