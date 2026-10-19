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
import os
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


def main():
    parser = argparse.ArgumentParser(
        description="Extract PDE border sequences directly from classified segments and assembly FASTA"
    )
    parser.add_argument("--pde-segments-bed", required=True, help="Path to pde_segments.bed")
    parser.add_argument("--assembly", required=True, help="Path to genome FASTA assembly")
    parser.add_argument("--breakpoints-bed", default=None, help="Optional path to filtered_breakpoints.bed for metadata")
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
    import json

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

            if strand == "-":
                subseq_oriented = revcomp(subseq)
                header = f"{bp['name']}_{chrom}_{pos + 1}_minus_win{win}bp_oriented soma_len={win}"
            else:
                subseq_oriented = subseq
                header = f"{bp['name']}_{chrom}_{pos + 1}_plus_win{win}bp_oriented soma_len={win}"

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

    out_census = os.path.join(args.outdir, f"{args.prefix}.census.json")
    with open(out_census, "w", encoding="utf-8") as f:
        json.dump({
            "n_breaks": len(frontier_breakpoints),
            "sample_id": args.prefix,
        }, f, indent=2)

    sys.stderr.write("Border flank extraction completed successfully.\n")


if __name__ == "__main__":
    main()
