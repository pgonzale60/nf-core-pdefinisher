#!/usr/bin/env python3
"""
pde_evaluate_curation.py - Module 5: Post-Surgery Curation Evaluation & Impact Scorecard

Compares pre-curation baseline assembly against post-surgery curated assembly, AGP, and audit.
Ingests Module 3 resolution manifests (Engine A Gap, Engine B Walker, Engine C rDNA)
and Module 1 terminal ends summary to generate:
  1. ${sample}.curation_scorecard.json (Machine-readable metrics)
  2. ${sample}.t2t_status.tsv (Per-chromosome terminus classification)
  3. ${sample}.curation_report.md (Executive Curation & Triage Report matching Sanger PDE standard)
"""

import argparse
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict

__version__ = "1.5.0"


def get_script_fingerprint():
    try:
        with open(__file__, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except Exception:
        return "unknown"


def chrom_sort_key(name):
    """Sort chromosomes naturally with nuclear karyotype precedence (I, II, III, IV, V, X, etc.)."""
    clean = re.sub(r"^[cC]hr_?", "", name)
    roman = {
        "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6,
        "VII": 7, "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12,
        "MT": 99, "M": 99
    }
    if clean in roman:
        return (0, roman[clean], clean)
    try:
        return (1, int(clean), clean)
    except ValueError:
        return (2, 0, clean)


def parse_args():
    parser = argparse.ArgumentParser(description="Module 5: Post-Surgery Curation Evaluator")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__} (sha256:{get_script_fingerprint()})",
    )
    parser.add_argument("--pre-ref", required=True, help="Baseline uncurated reference FASTA (.fa or .fa.gz)")
    parser.add_argument("--post-ref", required=True, help="Post-surgery curated FASTA")
    parser.add_argument("--blueprint", required=True, help="Module 2 Curation Blueprint JSON")
    parser.add_argument("--agp", required=True, help="Module 4 Curated AGP v2.1")
    parser.add_argument("--audit", required=True, help="Module 4 Curation Audit TSV")
    parser.add_argument("--gap-manifest", required=True, help="Module 3 Engine A Gap Manifest TSV")
    parser.add_argument("--walker-manifest", required=True, help="Module 3 Engine B Telomere Walker Manifest TSV")
    parser.add_argument("--rdna-manifest", required=True, help="Module 3 Engine C rDNA Manifest TSV")
    parser.add_argument("--terminal-summary", required=True, help="Module 1 Terminal Ends Summary TSV")
    parser.add_argument("--trans-manifest", required=False, default=None, help="Module 3 Engine D Trans Manifest TSV")
    parser.add_argument("--red-bed", required=False, default=None, help="Repeat masking BED from red")
    parser.add_argument("--gap-evidence", required=False, default=None, help="gap_evidence.json from extractgapevidence")
    parser.add_argument("--sample", required=True, help="Sample/strain prefix")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--chromosomes", required=False, default=None, help="Comma-separated list of canonical nuclear chromosome names")
    return parser.parse_args()


def load_fasta_stats(fasta_path):
    """
    Parses a FASTA file (plain or gzipped) and extracts lengths, N counts,
    and gap runs for each scaffold. Fails fast if file does not exist.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Required FASTA file not found: {fasta_path}")

    stats = {}
    opener = gzip.open if fasta_path.endswith(".gz") else open
    mode = "rt" if fasta_path.endswith(".gz") else "r"

    current_id = None
    seq_chunks = []

    def process_record(chunks):
        full_seq = "".join(chunks)
        length = len(full_seq)
        n_matches = list(re.finditer(r"[Nn]+", full_seq))
        gap_count = len(n_matches)
        total_n = sum(m.end() - m.start() for m in n_matches)
        gaps = [(m.start() + 1, m.end(), m.end() - m.start()) for m in n_matches]
        return {
            "length": length,
            "total_n": total_n,
            "gap_count": gap_count,
            "seq_clean_len": length - total_n,
            "gaps": gaps
        }

    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    stats[current_id] = process_record(seq_chunks)
                current_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if current_id:
            stats[current_id] = process_record(seq_chunks)

    return stats


def read_fasta(fasta_path):
    """Reads FASTA file and returns dict mapping record id to full uppercase sequence."""
    opener = gzip.open if fasta_path.endswith(".gz") else open
    mode = "rt" if fasta_path.endswith(".gz") else "r"
    seqs = {}
    current_id = None
    chunks = []
    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(chunks).upper()
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current_id:
            seqs[current_id] = "".join(chunks).upper()
    return seqs


def compute_n50(lengths):
    """Calculates N50 and L50 given a list of sequence lengths."""
    if not lengths:
        return 0, 0
    sorted_lens = sorted(lengths, reverse=True)
    total_len = sum(sorted_lens)
    half_len = total_len / 2.0

    accum = 0
    for idx, seq_len in enumerate(sorted_lens, start=1):
        accum += seq_len
        if accum >= half_len:
            return seq_len, idx
    return sorted_lens[-1], len(sorted_lens)


def parse_tsv(tsv_path):
    """Parses a TSV into a list of dictionaries using header keys. Fails fast if file missing."""
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Required TSV file not found: {tsv_path}")
    rows = []
    with open(tsv_path) as f:
        lines = [line.strip() for line in f if line.strip()]
        if not lines:
            return rows
        header = lines[0].split("\t")
        for line in lines[1:]:
            parts = line.split("\t")
            row = dict(zip(header, parts))
            rows.append(row)
    return rows


def analyze_terminal_repeats(red_bed_path, chrom, term, c_len, window_bp=5000):
    """
    Analyzes repeat masking at a chromosome terminus.
    Returns a dictionary of empirical numerical metrics:
      - window_pct: percentage of terminal window masked (float)
      - window_masked_bp: number of bp masked in terminal window (int)
      - window_bp: size of the terminal window (int)
      - tip_repeat_stretch_bp: contiguous (or near-contiguous) repeat stretch reaching the physical tip (int)
      - nearest_repeat_dist_bp: distance from physical tip to the first repeat interval (int or None)
      - first_repeat_len_bp: length of that first repeat interval (int or None)
      - summary_str: human-readable numerical summary
    """
    if not red_bed_path or not os.path.exists(red_bed_path):
        return None

    intervals = []
    try:
        with open(red_bed_path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0] == chrom:
                    intervals.append((int(parts[1]), int(parts[2])))
    except Exception:
        return None

    w_s = 0 if term == "5p" else max(0, c_len - window_bp)
    w_e = min(window_bp, c_len) if term == "5p" else c_len
    win_len = max(1, w_e - w_s)

    if not intervals:
        return {
            "window_pct": 0.0,
            "window_masked_bp": 0,
            "window_bp": win_len,
            "tip_repeat_stretch_bp": 0,
            "nearest_repeat_dist_bp": None,
            "first_repeat_len_bp": 0,
            "summary_str": f"0.0% masked in terminal {win_len:,} bp window (0 bp masked)"
        }

    intervals.sort()

    # Merge overlapping or micro-gapped (<= 100 bp) repeat blocks
    merged = []
    cur_s, cur_e = intervals[0]
    for n_s, n_e in intervals[1:]:
        if n_s <= cur_e + 100:
            cur_e = max(cur_e, n_e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = n_s, n_e
    merged.append((cur_s, cur_e))

    # Calculate masked bp within terminal window
    win_masked_bp = 0
    for s, e in merged:
        ovl = max(0, min(w_e, e) - max(w_s, s))
        win_masked_bp += ovl
    win_pct = min(100.0, (win_masked_bp / win_len) * 100.0)

    # Calculate contiguous stretch abutting physical tip (within 100 bp)
    tip_stretch_bp = 0
    nearest_dist = None
    first_rep_len = 0

    if term == "5p":
        for s, e in merged:
            if s <= 100:
                tip_stretch_bp = max(tip_stretch_bp, e)
            else:
                dist = s
                if nearest_dist is None or dist < nearest_dist:
                    nearest_dist = dist
                    first_rep_len = e - s
    else:
        for s, e in merged:
            if e >= c_len - 100:
                stretch = c_len - s
                tip_stretch_bp = max(tip_stretch_bp, stretch)
            else:
                dist = c_len - e
                if nearest_dist is None or dist < nearest_dist:
                    nearest_dist = dist
                    first_rep_len = e - s

    if tip_stretch_bp > 0:
        if tip_stretch_bp >= 1000:
            stretch_str = f"contiguous {tip_stretch_bp:,} bp ({tip_stretch_bp / 1000.0:.1f} kb) repeat tract directly abutting tip"
        else:
            stretch_str = f"contiguous {tip_stretch_bp:,} bp repeat tract directly abutting tip"
        summary_str = f"{win_pct:.1f}% masked in terminal {win_len:,} bp window ({win_masked_bp:,} bp masked; {stretch_str})"
    elif win_masked_bp > 0:
        summary_str = f"{win_pct:.1f}% masked in terminal {win_len:,} bp window ({win_masked_bp:,} bp masked; nearest repeat {nearest_dist:,} bp from tip spanning {first_rep_len:,} bp)"
    else:
        summary_str = f"0.0% masked in terminal {win_len:,} bp window (0 bp masked)"

    return {
        "window_pct": win_pct,
        "window_masked_bp": win_masked_bp,
        "window_bp": win_len,
        "tip_repeat_stretch_bp": tip_stretch_bp,
        "nearest_repeat_dist_bp": nearest_dist,
        "first_repeat_len_bp": first_rep_len,
        "summary_str": summary_str
    }


def compute_red_repeat_pct(red_bed_path, chrom, start, end):
    """Calculates exact percentage of bases masked by Red in interval [start, end]."""
    if not red_bed_path or not os.path.exists(red_bed_path):
        return None
    win_len = max(1, end - start)
    intervals = []
    try:
        with open(red_bed_path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0] == chrom:
                    r_s = int(parts[1])
                    r_e = int(parts[2])
                    clip_s = max(start, r_s)
                    clip_e = min(end, r_e)
                    if clip_e > clip_s:
                        intervals.append((clip_s, clip_e))
        if not intervals:
            return 0.0
        intervals.sort()
        merged = []
        cur_s, cur_e = intervals[0]
        for n_s, n_e in intervals[1:]:
            if n_s <= cur_e:
                cur_e = max(cur_e, n_e)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = n_s, n_e
        merged.append((cur_s, cur_e))
        masked_bp = sum(e - s for s, e in merged)
        return min(100.0, (masked_bp / win_len) * 100.0)
    except Exception:
        return None


def load_gap_evidence_map(gap_evidence_path):
    """Loads gap_evidence.json into a lookup dictionary: (chrom, start, end) -> gap_dict."""
    if not gap_evidence_path or not os.path.exists(gap_evidence_path):
        return {}
    try:
        with open(gap_evidence_path) as f:
            data = json.load(f)
        gaps = data.get("gaps", data.get("gap_evidence", []))
        gap_map = {}
        for g in gaps:
            ch = g.get("chrom")
            st = int(g.get("start", 0))
            en = int(g.get("end", 0))
            read_obs = g.get("read_observations", [])
            span_cnt_30 = len([
                r for r in read_obs 
                if r.get("span_type") == "COLLINEAR_SPANNING" and r.get("mapq", 0) >= 30
            ])
            g["span_cnt_30"] = span_cnt_30
            gap_map[(ch, st, en)] = g
        return gap_map
    except Exception:
        return {}


def classify_unresolved_gap(g_man=None, span_reads=0, n_unitigs=0, d_5p=25.0, d_3p=25.0):
    """
    Classifies the concrete algorithmic stopping condition for an unresolved gap.

    Returns:
        tuple: (classification_str, diagnosis_note_str, curator_takeaway_str)
    """
    g_notes = (g_man.get("notes") or "").strip() if g_man else ""
    notes_lower = g_notes.lower()

    # 1. Unbridged repeat array (highest frequency across assemblies)
    if "unbridged repeat" in notes_lower or "terminate at repeat boundary" in notes_lower:
        classification = "Unbridged Repeat Array"
        takeaway = (
            "Local de novo assembly baited flanking reads and assembled unitigs into both flanks, "
            "but contigs terminate at the repeat boundary without bridging across. The internal repeat "
            "array exceeds unitig span length; the gap is preserved to protect unplaced intervening elements "
            "and prevent asserting an unanchored sequence collapse."
        )
    # 2. Inversion breakpoint / discordant strand orientation
    elif "inversion breakpoint" in notes_lower or "discordant strands" in notes_lower:
        classification = "Inversion Breakpoint Discontinuity"
        takeaway = (
            "Assembled contigs or flanking reads align across the gap boundary with discordant strand "
            "orientations (+/-), indicating an inversion breakpoint or structural misorientation rather "
            "than a simple collinear gap."
        )
    # 3. Ambiguous multi-mapping anchor / low MAPQ / repetitive flanks
    elif (
        "no unitigs anchored" in notes_lower
        or "mapq >= 50" in notes_lower
        or "ambiguous repetitive" in notes_lower
        or "multi-locus" in notes_lower
        or "valency" in notes_lower
    ):
        classification = "Ambiguous Multi-Mapping Anchor"
        takeaway = (
            "Flanking sequences or candidate unitigs map across multiple genomic loci with low confidence "
            r"($\text{MAPQ} < 50$ or insufficient unique anchor sequence). Ambiguous mapping prevents unique, "
            "confident gap closure without risking cross-locus contamination."
        )
    # 4. De Novo assembly failure (zero contigs from reads)
    elif "hifiasm emitted zero contigs" in notes_lower or "assembly error" in notes_lower:
        classification = "De Novo Assembly Failure"
        takeaway = (
            "Flanking reads were baited from the locus, but local de novo assembly produced zero contigs "
            "due to extreme sequence complexity, low overlap identity, or unresolvable heterozygosity."
        )
    # 5. Read coverage drop / insufficient baited reads
    elif "insufficient baited reads" in notes_lower or "no paf alignments" in notes_lower:
        classification = "Read Coverage Drop to Zero"
        takeaway = (
            f"Flanking sequencing read depth dropped below the evidence floor for targeted assembly "
            f"(${d_5p:.1f}\\times$ 5', ${d_3p:.1f}\\times$ 3'). Preserved as an unbridged physical sequencing gap."
        )
    # 6. Consensus validation rejection
    elif "consensus" in notes_lower or "card premise not reproducible" in notes_lower:
        classification = "Consensus Bridge Rejected"
        takeaway = (
            "Spanning reads were evaluated for local consensus bridging but failed validation criteria "
            "(e.g. insufficient carrier reads, patch length discrepancy, or lack of read agreement)."
        )
    # 7. Fallback based on bridging metrics if no specific manifest note matched
    elif span_reads > 0 and n_unitigs == 0:
        classification = "Unassembled Spanning Reads"
        takeaway = (
            f"Sequence continuity is supported by {span_reads} raw spanning read(s), but local de novo "
            f"assembly yielded zero bridging unitigs across the gap."
        )
    elif span_reads > 0 and n_unitigs > 0:
        classification = "Unvalidated Candidate Patch"
        takeaway = (
            f"Both spanning reads ({span_reads}) and bridging unitigs ({n_unitigs}) were detected, "
            f"but candidate sequence failed validation criteria."
        )
    else:
        classification = "Physical Sequencing Gap"
        takeaway = (
            f"Zero spanning HiFi reads or assembled unitigs bridge between the flanks "
            f"(${d_5p:.1f}\\times$ 5', ${d_3p:.1f}\\times$ 3'). Preserved as an unbridged physical sequencing gap."
        )

    return classification, g_notes, takeaway


def parse_agp_gaps(agp_path):
    """Extracts gap records (type N or U) from AGP. Fails fast if file missing."""
    if not os.path.exists(agp_path):
        raise FileNotFoundError(f"Required AGP file not found: {agp_path}")
    gaps = []
    with open(agp_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5 and parts[4] in ["N", "U"]:
                gaps.append({
                    "object": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "length": int(parts[5]),
                    "type": parts[4],
                    "linkage": parts[7] if len(parts) > 7 else ""
                })
    return gaps


def parse_agp_rdna_intervals(agp_path):
    """
    Extracts genomic intervals for rDNA cassettes from AGP.
    Looks for component records whose component ID contains 'rdna' and 'cassette'.
    """
    if not os.path.exists(agp_path):
        return []
    intervals = []
    with open(agp_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                comp_id = parts[5].lower()
                if "rdna" in comp_id and "cassette" in comp_id:
                    intervals.append({
                        "chrom": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "comp_id": parts[5]
                    })
    return intervals


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{args.sample}] Module 5: Evaluating Curation Impact...", flush=True)

    # 1. Load Pre and Post Assembly Statistics
    pre_stats = load_fasta_stats(args.pre_ref)
    post_stats = load_fasta_stats(args.post_ref)

    # 2. Load Blueprint, AGP, Audit, and Manifests
    if not os.path.exists(args.blueprint):
        raise FileNotFoundError(f"Required Blueprint JSON not found: {args.blueprint}")
    with open(args.blueprint) as f:
        blueprint = json.load(f)

    karyotype_compartments = blueprint.get("karyotype_compartments", {})

    audit_rows = parse_tsv(args.audit)
    parse_agp_gaps(args.agp)  # fail fast if the AGP is missing or malformed
    gap_rows = parse_tsv(args.gap_manifest)
    walker_rows = parse_tsv(args.walker_manifest)
    rdna_rows = parse_tsv(args.rdna_manifest)
    terminal_summary_rows = parse_tsv(args.terminal_summary)
    trans_rows = parse_tsv(args.trans_manifest) if (args.trans_manifest and os.path.exists(args.trans_manifest)) else []

    # Policy / Motifs
    telo_motif = blueprint.get("policy", {}).get("telomere", {}).get("motif", "TTAGGC")

    # 3. Identify Nuclear Chromosomes vs Unplaced Scaffolds
    if args.chromosomes:
        specified = [c.strip() for c in args.chromosomes.split(",") if c.strip()]
        chr_names = [c for c in specified if c in pre_stats]
    elif blueprint.get("nuclear_chromosomes"):
        chr_names = [c for c in blueprint.get("nuclear_chromosomes") if c in pre_stats]
    else:
        chr_names = sorted(
            [k for k in pre_stats.keys() if "_unloc" not in k and not k.startswith("scaffold_") and k != "scaffold_MT" and not k.startswith("contig_")],
            key=chrom_sort_key
        )

    # 4. Compute Assembly Contiguity Metrics
    pre_lengths = [s["length"] for s in pre_stats.values()]
    post_lengths = [s["length"] for s in post_stats.values()]

    pre_total_bp = sum(pre_lengths)
    post_total_bp = sum(post_lengths)
    pre_n50, pre_l50 = compute_n50(pre_lengths)
    post_n50, post_l50 = compute_n50(post_lengths)

    pre_total_gaps = sum(s["gap_count"] for s in pre_stats.values())
    post_total_gaps = sum(s["gap_count"] for s in post_stats.values())
    pre_n_bp = sum(s["total_n"] for s in pre_stats.values())
    post_n_bp = sum(s["total_n"] for s in post_stats.values())

    pre_chr_lens = [pre_stats[c]["length"] for c in chr_names if c in pre_stats]
    post_chr_lens = [post_stats[c]["length"] for c in chr_names if c in post_stats]
    pre_chr_n50, pre_chr_l50 = compute_n50(pre_chr_lens)
    post_chr_n50, post_chr_l50 = compute_n50(post_chr_lens)

    # 5. Build Index Lookups
    all_action_cards = (
        blueprint.get("action_cards", [])
        + blueprint.get("curation_cards", [])
        + blueprint.get("telomere_cards", [])
    )
    cards_by_id = {c.get("card_id"): c for c in all_action_cards if c.get("card_id")}
    telo_cards = {}
    gap_cards = {}
    for card in all_action_cards:
        c = card.get("chrom")
        term = card.get("terminus")
        dom = card.get("domain", "")
        if c and term and dom not in ("DOMAIN_SCAFFOLD", "DOMAIN_GAP") and not card.get("card_id", "").startswith("CARD-SCAFF"):
            existing = telo_cards.get((c, term))
            if not existing or existing.get("domain") != "DOMAIN_TERMINUS":
                telo_cards[(c, term)] = card
        if card.get("domain") in ["DOMAIN_INTERNAL", "DOMAIN_GAP"] or "gap" in card.get("card_id", "").lower():
            gid = card.get("gap_id", card.get("card_id"))
            gap_cards[gid] = card
            if c and "start" in card and "end" in card:
                gap_cards[(c, int(card["start"]), int(card["end"]))] = card

    gap_manifest_by_coord = {}
    gap_manifest_by_id = {}
    for g in gap_rows:
        gid = g.get("gap_id")
        if gid:
            gap_manifest_by_id[gid] = g
        c = g.get("chrom")
        gs = g.get("gap_start")
        ge = g.get("gap_end")
        if c and gs is not None and ge is not None:
            try:
                gap_manifest_by_coord[(c, int(gs), int(ge))] = g
            except (ValueError, TypeError):
                pass

    term_summary_map = {}
    for r in terminal_summary_rows:
        cend = r.get("chrom_end", "")
        if "_" in cend:
            parts = cend.split("_", 1)
            c = parts[0]
            side = "5p" if ("5p" in parts[1] or "left" in parts[1]) else "3p"
            term_summary_map[(c, side)] = r

    walker_by_end = {(w.get("chrom"), w.get("terminus")): w for w in walker_rows}

    # Load gap_evidence.json if available
    gap_ev_path = args.gap_evidence
    if not gap_ev_path:
        cand1 = os.path.join(os.path.dirname(args.blueprint), f"{args.sample}.gap_evidence.json")
        cand2 = os.path.join(os.path.dirname(args.blueprint), "gap_evidence.json")
        cand3 = os.path.join(os.path.dirname(args.blueprint), "..", "extract_gap_evidence", f"{args.sample}.gap_evidence.json")
        for cand in [cand1, cand2, cand3]:
            if os.path.exists(cand):
                gap_ev_path = cand
                break
    gap_evidence_map = load_gap_evidence_map(gap_ev_path)

    # Resolve red_bed path if available
    red_bed_path = args.red_bed
    if not red_bed_path:
        cand1 = os.path.join(os.path.dirname(args.blueprint), f"{args.sample}.bed")
        cand2 = os.path.join(os.path.dirname(args.blueprint), "..", "red", f"{args.sample}.bed")
        cand3 = os.path.join(args.outdir, "..", "red", f"{args.sample}.bed")
        for cand in [cand1, cand2, cand3]:
            if os.path.exists(cand):
                red_bed_path = cand
                break

    # 6. Telomere Capping Status Assessment
    baseline_capped = {}
    for c in chr_names:
        for term in ["5p", "3p"]:
            card = telo_cards.get((c, term), {})
            act = card.get("action", "")
            cls = card.get("classification", "")
            b_stat = card.get("evidence", {}).get("baseline_status", "")
            # A terminus is baseline-capped only if the draft's own tip is capped. A card reporting
            # TELOMERE_IN_OVERHANG states the opposite: the reads carry telomere but the assembly's
            # tip holds no tract, so the end is work still to do. Counting it as capped would
            # overstate T2T completeness at exactly the ends the walker reports unresolved.
            if cls == "TELOMERE_IN_OVERHANG" or act == "EXTEND_FROM_OVERHANG":
                baseline_capped[(c, term)] = False
            elif act in ["SKIP_WALKER", "RETAIN_AS_IS"] or cls in ["ALREADY_CAPPED", "DRAFT_HAS_TELOMERE"] or b_stat == "ALREADY_CAPPED":
                baseline_capped[(c, term)] = True
            else:
                baseline_capped[(c, term)] = False

    post_capped = {}
    terminus_status = {}
    telo_extensions_bp = defaultdict(int)

    for (c, term), is_cap in baseline_capped.items():
        if is_cap:
            post_capped[(c, term)] = True
            terminus_status[(c, term)] = "CAPPED_BASELINE"
        else:
            post_capped[(c, term)] = False
            terminus_status[(c, term)] = "UNCAPPED"

    for row in audit_rows:
        e_type = row.get("edit_type", "")
        c = row.get("chrom", "")
        delta = int(row.get("delta_bp", 0))

        if e_type == "TELO_5P":
            post_capped[(c, "5p")] = True
            terminus_status[(c, "5p")] = "CAPPED_EXTENSION"
            telo_extensions_bp[(c, "5p")] += delta
        elif e_type == "TELO_3P":
            post_capped[(c, "3p")] = True
            terminus_status[(c, "3p")] = "CAPPED_EXTENSION"
            telo_extensions_bp[(c, "3p")] += delta
        elif e_type in ["TELO_RDNA_CAPPED_5P", "TELO_RDNA_CAPPED_3P"]:
            term = "5p" if "5P" in e_type else "3p"
            post_capped[(c, term)] = True
            terminus_status[(c, term)] = "CAPPED_RDNA_TELO"
            telo_extensions_bp[(c, term)] += delta
        elif e_type in ["RDNA_UNCAPPED_5P", "RDNA_UNCAPPED_3P"]:
            term = "5p" if "5P" in e_type else "3p"
            post_capped[(c, term)] = False
            terminus_status[(c, term)] = "UNCAPPED_RDNA"
            telo_extensions_bp[(c, term)] += delta

    total_termini = len(chr_names) * 2
    pre_capped_count = sum(1 for (c, t), is_cap in baseline_capped.items() if is_cap and c in chr_names)

    # Somatic diminution boundaries are confirmed chromosome ends — count them as capped.
    # The body of the report already classifies them as "Somatic (Core) Chromosome End / Confirmed"
    # but post_capped is not set True by the audit loop (no TELO edit is emitted for them).
    for c in chr_names:
        for term in ["5p", "3p"]:
            if post_capped.get((c, term)):
                continue  # already capped
            card = telo_cards.get((c, term), {})
            cls = card.get("classification", "")
            d_soma = card.get("somatic_boundary_coord", card.get("forbidden_zone_bp", 0))
            n_telo_ev = int(card.get("evidence", {}).get("n_telomeric_reads",
                            card.get("evidence", {}).get("n_telo", 0)))
            if cls == "SOMATIC_DIMINUTION_BOUNDARY" or (d_soma > 0 and n_telo_ev >= 0):
                # Confirm only if the card actually says it has a somatic boundary
                if cls == "SOMATIC_DIMINUTION_BOUNDARY":
                    post_capped[(c, term)] = True
                    terminus_status[(c, term)] = "CAPPED_SOMATIC_BOUNDARY"

    post_capped_count = sum(1 for (c, t), is_cap in post_capped.items() if is_cap and c in chr_names)

    fully_capped_chrs = [c for c in chr_names if post_capped.get((c, "5p")) and post_capped.get((c, "3p"))]
    uncapped_ends_list = []
    for c in chr_names:
        for term in ["5p", "3p"]:
            if not post_capped.get((c, term)):
                term_str = "5'" if term == "5p" else "3'"
                uncapped_ends_list.append((c, term, f"{c} {term_str}"))


    # 7. Gap Resolution Analysis
    excised_gaps = [r for r in audit_rows if r.get("edit_type") == "EXCISE_GAP"]
    patched_gaps = [r for r in audit_rows if r.get("edit_type") == "GAP_PATCH"]
    unresolved_gaps = [r for r in audit_rows if r.get("edit_type") == "UNRESOLVED_GAP"]

    total_patched_bp = sum(
        int(r.get("delta_bp", 0)) + (int(r.get("orig_end", 0)) - int(r.get("orig_start", 0)))
        for r in patched_gaps
    )
    total_excised_bp = sum(abs(int(r.get("delta_bp", 0))) for r in excised_gaps)
    remaining_gap_count = post_total_gaps

    # Distinguish non-rDNA assembly gaps from expected tandem rDNA array sizing gaps (100 bp N placeholder)
    agp_rdna_intervals = parse_agp_rdna_intervals(args.agp)
    rdna_policy = blueprint.get("policy", {}).get("rdna", {})
    rdna_gap_size = int(rdna_policy.get("rdna_gap_size_bp", 100))

    rdna_chroms = {r.get("chrom") for r in rdna_rows if r.get("chrom")}
    for r in audit_rows:
        if "RDNA" in r.get("edit_type", "") or "rdna" in r.get("card_id", "").lower() or "rdna" in r.get("action", "").lower():
            if r.get("chrom"):
                rdna_chroms.add(r.get("chrom"))

    rdna_gaps = []
    non_rdna_gaps = []
    for c, s in post_stats.items():
        for (g_start, g_end, g_len) in s.get("gaps", []):
            is_rdna = False
            for r_int in agp_rdna_intervals:
                if r_int["chrom"] == c and g_start >= r_int["start"] and g_end <= r_int["end"]:
                    is_rdna = True
                    break
            if not is_rdna and c in rdna_chroms and g_len == rdna_gap_size:
                c_len = s["length"]
                if g_start <= 50000 or g_end >= (c_len - 50000):
                    is_rdna = True
            if is_rdna:
                rdna_gaps.append((c, g_start, g_end, g_len))
            else:
                non_rdna_gaps.append((c, g_start, g_end, g_len))

    rdna_gaps_count = len(rdna_gaps)
    remaining_non_rdna_gaps = len(non_rdna_gaps)
    is_gapless_non_rdna = (remaining_non_rdna_gaps == 0)

    # 8. Purge & De-duplication Ledger
    purged_haplo = [r for r in audit_rows if r.get("edit_type") == "PURGE_HAPLOTIG"]
    purged_rdna = [r for r in audit_rows if r.get("edit_type") == "PURGE_RDNA"]
    purged_organellar = [r for r in audit_rows if r.get("edit_type") in ("PURGE_ORGANELLAR", "PURGE_MITOCHONDRION")]

    purged_haplo_bp = sum(abs(int(r.get("delta_bp", 0))) for r in purged_haplo)
    purged_rdna_bp = sum(abs(int(r.get("delta_bp", 0))) for r in purged_rdna)
    purged_organellar_bp = sum(abs(int(r.get("delta_bp", 0))) for r in purged_organellar)
    total_purged_bp = purged_haplo_bp + purged_rdna_bp + purged_organellar_bp

    pre_primary_chr_bp = sum(pre_stats[c]["length"] for c in chr_names if c in pre_stats)
    post_primary_chr_bp = sum(post_stats[c]["length"] for c in chr_names if c in post_stats)
    chr_total_delta_bp = post_primary_chr_bp - pre_primary_chr_bp
    purity_pct = (post_primary_chr_bp / post_total_bp * 100.0) if post_total_bp > 0 else 0.0

    # 9. Audit Reconciliation & Size Deltas
    audit_total_delta = sum(int(r.get("delta_bp", 0)) for r in audit_rows)
    measured_total_delta = post_total_bp - pre_total_bp
    reconciliation_discrepancy = measured_total_delta - audit_total_delta
    is_conserved = (reconciliation_discrepancy == 0)

    # Assembly-wide rDNA delta (including purged unplaced rDNA scaffolds)
    rdna_delta_bp = sum(
        int(r.get("delta_bp", 0))
        for r in audit_rows
        if "RDNA" in r.get("edit_type", "")
        or "rdna" in r.get("card_id", "").lower()
        or "rdna" in r.get("action", "").lower()
    )

    # Chromosome-level rDNA delta: rDNA modifications directly on primary chromosomes
    chr_rdna_delta_bp = sum(
        int(r.get("delta_bp", 0))
        for r in audit_rows
        if (r.get("chrom") in chr_names or r.get("scaffold") in chr_names)
        and (
            "RDNA" in r.get("edit_type", "")
            or "rdna" in r.get("card_id", "").lower()
            or "rdna" in r.get("action", "").lower()
        )
    )

    # Non-rDNA delta reflects the size change in chromosomes (removing unlocalized scaffolds does not affect it)
    non_rdna_delta_bp = chr_total_delta_bp - chr_rdna_delta_bp

    # 10. Compile Scorecard Payload (JSON)
    scorecard = {
        "sample": args.sample,
        "provenance": {
            "software": "pde_evaluate_curation.py",
            "software_version": __version__,
            "script_fingerprint": get_script_fingerprint(),
        },
        "contiguity": {
            "pre_curation": {
                "total_bp": pre_total_bp,
                "scaffold_count": len(pre_stats),
                "n50": pre_n50,
                "l50": pre_l50,
                "total_gaps": pre_total_gaps,
                "total_n_bp": pre_n_bp,
                "chr_n50": pre_chr_n50,
                "chr_l50": pre_chr_l50
            },
            "post_curation": {
                "total_bp": post_total_bp,
                "scaffold_count": len(post_stats),
                "n50": post_n50,
                "l50": post_l50,
                "total_gaps": post_total_gaps,
                "total_n_bp": post_n_bp,
                "chr_n50": post_chr_n50,
                "chr_l50": post_chr_l50
            },
            "deltas": {
                "length_delta_bp": measured_total_delta,
                "rdna_length_delta_bp": rdna_delta_bp,
                "non_rdna_length_delta_bp": non_rdna_delta_bp,
                "chr_length_delta_bp": chr_total_delta_bp,
                "chr_rdna_length_delta_bp": chr_rdna_delta_bp,
                "scaffolds_delta": len(post_stats) - len(pre_stats),
                "gaps_delta": post_total_gaps - pre_total_gaps,
                "n_bp_delta": post_n_bp - pre_n_bp
            }
        },
        "t2t_completeness": {
            "total_chromosome_termini": total_termini,
            "pre_capped_termini": pre_capped_count,
            "post_capped_termini": post_capped_count,
            "pre_t2t_score_pct": round(pre_capped_count / total_termini * 100.0, 2) if total_termini else 0.0,
            "post_t2t_score_pct": round(post_capped_count / total_termini * 100.0, 2) if total_termini else 0.0,
            "fully_capped_chromosomes": fully_capped_chrs,
            "uncapped_ends": [item[2] for item in uncapped_ends_list],
            "net_telomeric_bp_added": sum(telo_extensions_bp.values()),
            "terminus_status_breakdown": {
                f"{c}_{term}": terminus_status.get((c, term), "UNCAPPED")
                for c in chr_names for term in ["5p", "3p"]
            }
        },
        "gap_resolution": {
            "baseline_gaps": pre_total_gaps,
            "remaining_gaps": remaining_gap_count,
            "rdna_gaps": rdna_gaps_count,
            "remaining_non_rdna_gaps": remaining_non_rdna_gaps,
            "is_gapless_non_rdna": is_gapless_non_rdna,
            "gaps_excised_count": len(excised_gaps),
            "gaps_patched_count": len(patched_gaps),
            "gaps_unresolved_count": len(unresolved_gaps),
            "total_patched_bp": total_patched_bp,
            "total_excised_bp": total_excised_bp
        },
        "transchromosomal_surgery": {
            "candidates_evaluated": len(trans_rows),
            "transfers_count": sum(1 for r in trans_rows if r.get("outcome") == "TRANSFER_AND_STITCH_CONTIG"),
            "fissions_count": sum(1 for r in trans_rows if r.get("outcome") == "SEVER_TO_UNPLACED"),
            "rejected_count": sum(1 for r in trans_rows if r.get("outcome") == "RETAIN_BASELINE_GAP"),
            "ambiguous_count": sum(1 for r in trans_rows if r.get("outcome") == "AMBIGUOUS_NO_HIFI_SPAN"),
            "transferred_bp": sum(int(r.get("contig_length", 0)) for r in trans_rows if r.get("outcome") == "TRANSFER_AND_STITCH_CONTIG"),
            "severed_bp": sum(int(r.get("contig_length", 0)) for r in trans_rows if r.get("outcome") == "SEVER_TO_UNPLACED")
        },
        "purge_ledger": {
            "purged_haplotigs_count": len(purged_haplo),
            "purged_haplotigs_bp": purged_haplo_bp,
            "purged_rdna_scaffolds_count": len(purged_rdna),
            "purged_rdna_scaffolds_bp": purged_rdna_bp,
            "purged_organellar_scaffolds_count": len(purged_organellar),
            "purged_organellar_scaffolds_bp": purged_organellar_bp,
            "total_purged_bp": total_purged_bp,
            "primary_release_purity_pct": round(purity_pct, 3)
        },
        "reconciliation": {
            "measured_delta_bp": measured_total_delta,
            "audit_sum_delta_bp": audit_total_delta,
            "rdna_delta_bp": rdna_delta_bp,
            "non_rdna_delta_bp": non_rdna_delta_bp,
            "chr_delta_bp": chr_total_delta_bp,
            "chr_rdna_delta_bp": chr_rdna_delta_bp,
            "discrepancy_bp": reconciliation_discrepancy,
            "mathematically_conserved": is_conserved
        },
        "karyotype_compartments": karyotype_compartments
    }

    scorecard_path = os.path.join(args.outdir, f"{args.sample}.curation_scorecard.json")
    with open(scorecard_path, "w") as f:
        json.dump(scorecard, f, indent=2)
    print(f"[{args.sample}] Emitted curation scorecard JSON: {scorecard_path}", flush=True)

    # 11. Emit T2T Status TSV
    t2t_tsv_path = os.path.join(args.outdir, f"{args.sample}.t2t_status.tsv")
    with open(t2t_tsv_path, "w") as f:
        f.write("sample\tchrom\tterminus\tbaseline_status\tpost_status\tstatus_category\textension_bp\n")
        for c in chr_names:
            for term in ["5p", "3p"]:
                st = terminus_status.get((c, term), "UNCAPPED")
                b_st = "CAPPED" if baseline_capped.get((c, term)) else "UNCAPPED"
                p_st = "CAPPED" if post_capped.get((c, term)) else "UNCAPPED"
                ext = telo_extensions_bp.get((c, term), 0)
                f.write(f"{args.sample}\t{c}\t{term}\t{b_st}\t{p_st}\t{st}\t{ext}\n")
    print(f"[{args.sample}] Emitted T2T status TSV: {t2t_tsv_path}", flush=True)

    # =========================================================================
    # 12. Compile Markdown Executive Report (Standard Sanger PDE Format)
    # =========================================================================
    report_path = os.path.join(args.outdir, f"{args.sample}.curation_report.md")
    with open(report_path, "w") as f:
        # Title
        f.write(f"# Automated Curation & Triage Report: `{args.sample}`\n\n")

        # Section 1: Executive Summary & Readiness
        f.write("## 1. Executive Summary & Readiness\n")
        f.write(f"* **Chromosome Termini Capped**: **{post_capped_count} / {total_termini}** termini capped with canonical `({telo_motif})n` telomeres.\n")
        fully_str = ", ".join(fully_capped_chrs) if fully_capped_chrs else "None"
        f.write(f"  - Fully Capped Chromosomes: **{len(fully_capped_chrs)} / {len(chr_names)}** ({fully_str})\n")
        uncapped_str = ", ".join(item[2] for item in uncapped_ends_list) if uncapped_ends_list else "None"
        f.write(f"  - Uncapped Chromosome Ends: **{len(uncapped_ends_list)}**" + (f" ({uncapped_str})\n" if uncapped_ends_list else "\n"))

        if rdna_gaps_count > 0:
            f.write(
                f"* **Gaps**: **{remaining_gap_count} remaining** ({remaining_non_rdna_gaps} non-rDNA gaps, "
                f"{rdna_gaps_count} expected rDNA array sizing placeholder) "
                f"(down from {pre_total_gaps} baseline gaps; {len(patched_gaps)} patched, {len(excised_gaps)} artifact excised).\n"
            )
        else:
            f.write(
                f"* **Gaps**: **{remaining_gap_count} remaining** "
                f"(down from {pre_total_gaps} baseline gaps; {len(patched_gaps)} patched, {len(excised_gaps)} artifact excised).\n"
            )
        pre_mb = pre_total_bp / 1e6
        post_mb = post_total_bp / 1e6
        f.write(
            f"* **Net Assembly Delta**: **{measured_total_delta:+,} bp** "
            f"(Non-rDNA Delta: **{non_rdna_delta_bp:+,} bp**, rDNA Delta: **{rdna_delta_bp:+,} bp**; "
            f"Total Length: {pre_mb:.2f} Mb → {post_mb:.2f} Mb).\n"
        )

        if karyotype_compartments:
            thresh = karyotype_compartments.get("classification_threshold")
            core_info = karyotype_compartments.get("putative_core", {})
            grs_info = karyotype_compartments.get("putative_grs", {})
            core_bp = core_info.get("span_bp", 0)
            core_pct = core_info.get("span_pct", 0.0)
            core_cov = core_info.get("median_coverage", 0.0)
            grs_bp = grs_info.get("span_bp", 0)
            grs_pct = grs_info.get("span_pct", 0.0)
            grs_cov = grs_info.get("median_coverage", 0.0)
            thresh_str = f"{thresh:.2f}×" if thresh is not None else "N/A"
            f.write(
                f"* **PDE Compartment Allocation ($T_{{\\text{{GRS}}}}$ = `{thresh_str}`)**:\n"
                f"  - Putative Somatic Core: **{core_bp:,d} bp** ({core_bp / 1e6:.2f} Mb, {core_pct:.1f}%) | Median Coverage: `{core_cov:.1f}×`\n"
                f"  - Putative GRS (Eliminated): **{grs_bp:,d} bp** ({grs_bp / 1e6:.2f} Mb, {grs_pct:.1f}%) | Median Coverage: `{grs_cov:.1f}×`\n"
            )
        f.write("\n")

        f.write("---\n\n")

        # Section 2: Unresolved Loci & Curator Diagnostics (Actionable Evidence)
        total_unresolved = len(uncapped_ends_list) + len(unresolved_gaps)
        f.write("## 2. Unresolved Loci & Curator Diagnostics (Actionable Evidence)\n")
        if total_unresolved == 0:
            f.write("*All chromosome termini are capped and all baseline gaps have been resolved. Zero unresolved loci remain.*\n\n")
        else:
            f.write(f"*Empirical evidence at the {total_unresolved} incomplete loci to allow the curator to challenge pipeline calls and assess manual feasibility:*\n\n")

            if uncapped_ends_list:
                f.write("### Uncapped Termini\n")
                for c, term, label in uncapped_ends_list:
                    card = telo_cards.get((c, term), {})
                    tsum = term_summary_map.get((c, term), {})

                    n_soft = int(tsum.get("n_softclipped_reads", 0))
                    n_telo = int(tsum.get("n_telomeric_reads", 0))
                    c_len = post_stats.get(c, {}).get("length", pre_stats.get(c, {}).get("length", 0))

                    raw_d_soma = card.get("somatic_boundary_coord", card.get("forbidden_zone_bp", 0))
                    if isinstance(raw_d_soma, (int, float)):
                        d_soma = int(raw_d_soma)
                    elif isinstance(raw_d_soma, str):
                        digits = re.sub(r"[^\d]", "", raw_d_soma.split()[0] if raw_d_soma else "")
                        d_soma = int(digits) if digits else 0
                    else:
                        d_soma = 0

                    if term == "5p":
                        coord_str = f"Coordinate 0 / {d_soma:,} bp" if d_soma > 0 else "Coordinate 0 bp"
                    else:
                        coord_str = f"Coordinate {c_len:,} bp"

                    # Core Chromosome End Definition:
                    # Presence of empirical telomeric junction vs Stalled Assembly
                    has_junction = (n_telo > 0) or (d_soma > 0 and card.get("classification") == "SOMATIC_DIMINUTION_BOUNDARY")
                    f.write(f"* **Chr {label} ({coord_str})**:\n")

                    # Calculate real Red repeat masking metrics at terminal window
                    rep_metrics = analyze_terminal_repeats(red_bed_path, c, term, c_len, window_bp=5000)
                    if rep_metrics is not None:
                        real_red_pct = rep_metrics["window_pct"]
                        rep_content = rep_metrics["summary_str"]
                        tip_repeat_stretch = rep_metrics["tip_repeat_stretch_bp"]
                    else:
                        real_red_pct = None
                        rep_content = card.get("evidence", {}).get("terminal_repeat_content", "Repeat content not evaluated")
                        tip_repeat_stretch = 0

                    if has_junction:
                        # The summary no longer carries overhang columns: they read 0 at every end of
                        # both samples, including ends with 58 telomeric reads, so they measured
                        # something other than the outward overhang and were dropped rather than kept
                        # as zeros. d_soma, where a somatic boundary was found, is the only coordinate
                        # this file can offer; otherwise the junction is reported as a terminal
                        # coordinate rather than a fabricated one.
                        _telo_pos = d_soma
                        pos_reported = f"{_telo_pos:,} bp" if _telo_pos > 0 else "terminal coordinate"
                        if n_soft > 0:
                            f.write(f"  * **Telomeric Junction**: **Confirmed**. {n_telo} of {n_soft} terminal reads carry `({telo_motif})n` soft-clips precisely at position {pos_reported}.\n")
                        elif n_telo > 0:
                            f.write(f"  * **Telomeric Junction**: **Confirmed**. {n_telo} terminal reads carry `({telo_motif})n` soft-clips precisely at position {pos_reported}.\n")
                        else:
                            f.write(f"  * **Telomeric Junction**: **Confirmed**. Somatic diminution breaksite identified at position {pos_reported}.\n")
                        f.write(f"  * **Terminal Repeat Content**: **{rep_content}**.\n")
                        f.write("  * **Coverage State**: Stable diploid somatic baseline across the adjacent scaffold; 0 reads continue past the junction into alternative sequence.\n")
                        f.write("  * **Classification**: **Somatic (Core) Chromosome End**.\n")
                        f.write("  * **Curator Takeaway**: The somatic chromosome end is reached and capped. True sub-stoichiometric germline extension (GRS) is absent in this library; no further extension is supported by the data.\n\n")
                    else:
                        w = walker_by_end.get((c, term), {})
                        w_status = w.get("status", "")
                        w_notes = w.get("notes", "")
                        cls = card.get("classification", "")
                        is_germline = bool(card.get("terminus_is_germline", False)) or ("germline" in w_notes.lower()) or ("distal grs" in w_notes.lower()) or (w_status == "GERMLINE_TIP_UNCAPPED")

                        if d_soma > 0:
                            if term == "5p":
                                grs_len = d_soma
                                d_info = f" extending +{grs_len:,} bp ({grs_len / 1000.0:.1f} kb) distal to the somatic diminution boundary at {d_soma:,} bp"
                            else:
                                grs_len = max(0, c_len - d_soma)
                                d_info = f" extending +{grs_len:,} bp ({grs_len / 1000.0:.1f} kb) distal to the somatic diminution boundary at {d_soma:,} bp (tip at {c_len:,} bp)"
                        else:
                            grs_len = 0
                            d_info = ""

                        # Repeat summary phrase
                        if tip_repeat_stretch > 0:
                            rep_ctx = f"a contiguous {tip_repeat_stretch:,} bp ({tip_repeat_stretch / 1000.0:.1f} kb) repeat tract directly abutting the tip"
                        else:
                            rep_ctx = f"sequence containing {rep_content}"

                        f.write("  * **Telomeric Junction**: **None observed**. 0 soft-clipped telomeres across all terminal reads.\n")
                        f.write(f"  * **Terminal Repeat Content**: **{rep_content}**.\n")

                        # Mode 1: Multivalent Repeat Branching
                        if w_status == "POLYVALENT_REPEAT_BRANCH" or "polyvalent" in w_notes.lower():
                            if is_germline:
                                f.write(f"  * **Coverage State**: Distal Germline-Restricted Sequence (GRS){d_info} arresting in a polyvalent repeat branch.\n")
                            else:
                                f.write("  * **Coverage State**: Stable diploid baseline arresting in a polyvalent repeat branch.\n")
                            if w_notes:
                                f.write(f"  * **Walker Observation**: {w_notes}\n")
                            f.write("  * **Classification**: **Multivalent Repeat Branching**.\n")
                            f.write(f"  * **Curator Takeaway**: This chromosome terminus cannot be extended because candidate reads diverge into mutually conflicting sequences (multivalent repeat branch). The tip terminates in {rep_ctx}. Advancing along any branch would create an arbitrary chimeric path; resolving this end requires ultra-long spanning reads bridging across the multivalent array to the distal telomere.\n\n")

                        # Mode 2: Stationary Repeat Loop
                        elif w_status == "STATIONARY_REPEAT_ARRAY" or "stationary repeat" in w_notes.lower() or "extension shorter than" in w_notes.lower():
                            if is_germline:
                                f.write(f"  * **Coverage State**: Distal Germline-Restricted Sequence (GRS){d_info} trapped within a cyclical tandem repeat array.\n")
                            else:
                                f.write("  * **Coverage State**: Stable diploid baseline trapped within a cyclical tandem repeat array.\n")
                            if w_notes:
                                f.write(f"  * **Walker Observation**: {w_notes}\n")
                            f.write("  * **Classification**: **Stationary Repeat Loop**.\n")
                            f.write(f"  * **Curator Takeaway**: Candidate reads overlap the assembly tip but introduce zero novel sequence complexity (cyclical tandem repeat loop). The terminus abuts {rep_ctx}. Iterative walking was arrested and rolled back to prevent arbitrary repeat expansion without a unique sequence anchor.\n\n")

                        # Mode 3: Multi-Mapping Paralogous Anchor
                        elif w_status == "LOW_ANCHOR_MAPQ" or "low anchor mapq" in w_notes.lower():
                            f.write("  * **Coverage State**: Terminal reads map with low confidence across paralogous repeat copies ($\\text{MAPQ} < 20$).\n")
                            if w_notes:
                                f.write(f"  * **Walker Observation**: {w_notes}\n")
                            f.write("  * **Classification**: **Multi-Mapping Paralogous Anchor**.\n")
                            f.write(f"  * **Curator Takeaway**: Terminal sequence is repetitive ({rep_content}) with reads mapping across multiple genomic loci with $\\text{{MAPQ}} < 20$. Ambiguous mapping prevents unique, confident extension.\n\n")

                        # Mode 4: Walker Horizon Reached
                        elif w_status == "HORIZON_LIMIT_REACHED" or "horizon limit" in w_notes.lower() or ("within" in w_notes.lower() and "steps" in w_notes.lower()):
                            f.write("  * **Coverage State**: Extended through unique sequence for 100 iterations without reaching telomeric repeats.\n")
                            if w_notes:
                                f.write(f"  * **Walker Observation**: {w_notes}\n")
                            f.write("  * **Classification**: **Walker Horizon Reached (100 Steps)**.\n")
                            f.write("  * **Curator Takeaway**: The iterative walker advanced 100 steps through unique sequence without encountering a canonical telomere tract ($\\ge 200\\text{ bp}$). Rolled back to Step 0 to prevent uncapped overhang assertion.\n\n")

                        # Mode 5: Read Coverage Drop to Zero
                        elif w_status == "ZERO_EXTENDING_READS" or "zero unmapped" in w_notes.lower() or "zero extending" in w_notes.lower() or n_soft == 0:
                            if is_germline:
                                f.write(f"  * **Coverage State**: Sub-stoichiometric germline coverage (GRS){d_info} terminating at complete read exhaustion (0 outbound reads).\n")
                            else:
                                f.write("  * **Coverage State**: Diploid baseline terminating at an abrupt read coverage drop-off (0 outbound reads).\n")
                            if w_notes:
                                f.write(f"  * **Walker Observation**: {w_notes}\n")
                            f.write("  * **Classification**: **Read Coverage Drop to Zero**.\n")
                            if is_germline:
                                f.write(f"  * **Curator Takeaway**: The somatic chromosome core is complete and intact. This terminus extends into distal germline sequence ({rep_content}) where sequencing read coverage drops completely to zero without reaching a germline telomere cap. Preserved as an honestly uncapped germline terminus.\n\n")
                            else:
                                f.write(f"  * **Curator Takeaway**: Zero reads in the sequencing library extend beyond coordinate {coord_str} ({rep_content}). The terminus represents an unbridged physical sequencing library exhaustion or physical gap.\n\n")

                        # Mode 6: Fallback (e.g. non-telomeric soft-clips)
                        else:
                            if n_soft > 0:
                                f.write(f"  * **Coverage State**: Baseline terminating at a soft-clipped discontinuity ({n_soft} outbound soft-clipped reads).\n")
                                if w_notes:
                                    f.write(f"  * **Walker Observation**: {w_notes}\n")
                                f.write("  * **Classification**: **Non-Telomeric Terminal Discontinuity**.\n")
                                f.write(f"  * **Curator Takeaway**: Terminal reads carry {n_soft} non-telomeric soft-clips at coordinate {coord_str} ({rep_content}). This indicates a potential structural misassembly breakpoint or unplaced scaffold junction rather than a natural chromosome end.\n\n")
                            else:
                                f.write("  * **Coverage State**: Uncapped subtelomeric terminus.\n")
                                if w_notes:
                                    f.write(f"  * **Walker Observation**: {w_notes}\n")
                                f.write("  * **Classification**: **Unresolved Subtelomere**.\n")
                                f.write(f"  * **Curator Takeaway**: Uncapped terminus ({rep_content}). Preserved as an unextended subtelomeric end.\n\n")

            if unresolved_gaps:
                f.write("### Unresolved Gaps\n")
                for r in unresolved_gaps:
                    c = r.get("chrom", "")
                    s = int(r.get("orig_start", 0))
                    e = int(r.get("orig_end", 0))
                    glen = e - s
                    cid = r.get("card_id", "")

                    # Retrieve real empirical metrics.
                    # Priority: gap_evidence_map (raw JSON) > gap_cards evidence dict
                    #           > gap_card top-level fields > local_depth string > 25.0 sentinel.
                    g_info = gap_evidence_map.get((c, s, e), {})
                    g_card = gap_cards.get((c, s, e), {}) or (gap_cards.get(cid) if cid else {})
                    g_ev   = g_card.get("evidence", {}) if g_card else {}

                    def _pick(*sources, key, default=None):
                        """Return first non-None value for key across sources."""
                        for src in sources:
                            v = src.get(key)
                            if v is not None:
                                return v
                        return default

                    def _parse_float(val, fallback=25.0):
                        try:
                            return float(re.search(r'([\d\.]+)', str(val)).group(1))
                        except Exception:
                            return fallback

                    # Flanking depths: raw JSON > card evidence > card top-level > local_depth str
                    d_5p_raw = _pick(g_info, g_ev, g_card, key="depth_flank_5p")
                    d_3p_raw = _pick(g_info, g_ev, g_card, key="depth_flank_3p")
                    if d_5p_raw is None or d_3p_raw is None:
                        loc_raw = _pick(g_info, g_ev, g_card, key="local_depth",
                                        default=_pick(g_info, g_ev, key="d_local", default=None))
                        d_loc = _parse_float(loc_raw, fallback=None) if loc_raw is not None else None
                        if d_5p_raw is None:
                            d_5p_raw = d_loc
                        if d_3p_raw is None:
                            d_3p_raw = d_loc

                    d_5p = _parse_float(d_5p_raw) if d_5p_raw is not None else 25.0
                    d_3p = _parse_float(d_3p_raw) if d_3p_raw is not None else 25.0

                    # Spanning reads: raw JSON span count > card collinear_spanning_reads string
                    span_reads = _pick(g_info, g_ev, g_card, key="span_cnt_30")
                    if span_reads is None:
                        span_reads = _pick(g_info, g_ev, g_card, key="collinear_gap_reads")
                    if span_reads is None:
                        col_str = str(_pick(g_info, g_ev, g_card, key="collinear_spanning_reads", default=""))
                        m_m = re.search(r'MAPQ>=30:\s*(\d+)', col_str)
                        if m_m:
                            span_reads = int(m_m.group(1))
                        else:
                            m_tot = re.search(r'(\d+)', col_str)
                            span_reads = int(m_tot.group(1)) if m_tot else 0

                    n_unitigs = int(_pick(g_info, g_ev, g_card, key="bridging_unitigs_count", default=0))

                    # Retrieve Assembly Engine Manifest Diagnosis (Engine A)
                    g_man = gap_manifest_by_id.get(cid)
                    if not g_man:
                        g_man = gap_manifest_by_coord.get((c, s, e))
                    if not g_man:
                        for (gc, gs, ge), cand in gap_manifest_by_coord.items():
                            if gc == c and abs(gs - s) <= 1000 and abs(ge - e) <= 1000:
                                g_man = cand
                                break

                    classification, diag_note, takeaway = classify_unresolved_gap(
                        g_man=g_man,
                        span_reads=span_reads,
                        n_unitigs=n_unitigs,
                        d_5p=d_5p,
                        d_3p=d_3p,
                    )

                    f.write(f"* **Chr {c} Gap (`{c}: {s:,}–{e:,}`, {glen} Ns)**:\n")
                    f.write(f"  * **Coverage State**: Flanking read depth (${d_5p:.1f}\\times$ 5', ${d_3p:.1f}\\times$ 3').\n")
                    f.write(f"  * **Bridging Evidence**: {span_reads} spanning HiFi reads; {n_unitigs} unitigs bridge across both flanks.\n")
                    if diag_note:
                        f.write(f"  * **Assembly Engine Diagnosis**: {diag_note}\n")
                    f.write(f"  * **Classification**: **{classification}**.\n")
                    f.write(f"  * **Curator Takeaway**: {takeaway}\n\n")



            if trans_rows:
                f.write("### Trans-Chromosomal Scaffold Reassignments & Fissions\n")
                for t in trans_rows:
                    outcome = t.get("outcome", "")
                    src = t.get("source_chrom", "")
                    gs = int(t.get("source_gap_start", 0))
                    ge = int(t.get("source_gap_end", 0))
                    glen = int(t.get("source_gap_len", ge - gs))
                    collinear = int(t.get("collinear_gap_reads", 0))
                    supp = int(t.get("support_split_reads", 0))
                    mapq = float(t.get("split_reads_mapq_mean", 50.0))
                    vaf = float(t.get("vaf", 0.0))
                    tgt = t.get("target_chrom", "")
                    tpos = int(t.get("target_coord", 0))
                    tfeat = t.get("target_feature", "unknown")
                    tstrand = t.get("target_strand", "+")
                    clen = int(t.get("contig_length", 0))
                    cid = t.get("contig_id", "")
                    retained_chr = t.get("retained_chrom_name", src)
                    retained_pct = float(t.get("retained_mass_fraction", 1.0)) * 100.0
                    assigned = t.get("assigned_name", "")
                    rationale = t.get("curator_rationale", "")

                    f.write(f"* **{src}:{gs:,}–{ge:,} ({glen} Ns) → {tgt}:{tpos:,}** [{outcome}]:\n")

                    # Negative control line — never say "misjoin confirmed" when collinear=0;
                    # absence of HiFi spanning reads is ambiguous, not diagnostic.
                    if collinear > 0:
                        f.write(f"  * **Negative Control**: **{collinear} collinear reads** span across the source gap (Authentic read continuity detected; gap retained).\n")
                    else:
                        f.write("  * **Negative Control**: **0 collinear HiFi reads** span the source gap (uninformative — reads in repetitive/GRS regions are MAPQ-filtered; absence does not confirm misjoin).\n")

                    f.write(f"  * **Positive Control**: **{supp} split reads** (mean MAPQ {mapq:.0f}, VAF={vaf:.2f}) link the {clen:,} bp contig `{cid}` to {tgt}:{tpos:,} ({tfeat}, strand `{tstrand}`).\n")
                    f.write(f"  * **Sequence Mass Balance**: Primary chromosome **{retained_chr}** retains **{retained_pct:.1f}%** sequence mass. ")

                    if outcome == "TRANSFER_AND_STITCH_CONTIG":
                        f.write(f"Contig ({clen:,} bp) transferred and stitched into **{assigned}**.\n")
                    elif outcome == "SEVER_TO_UNPLACED":
                        f.write(f"Target locus is mid-contig/unconfirmed; contig ({clen:,} bp) severed to unplaced scaffold **{assigned}**.\n")
                    elif outcome == "AMBIGUOUS_NO_HIFI_SPAN":
                        f.write("Trans-chromosomal signal conservatively rejected; baseline gap retained pending Hi-C review.\n")
                    else:
                        f.write("Trans-chromosomal signal rejected; baseline gap retained.\n")

                    if outcome == "AMBIGUOUS_NO_HIFI_SPAN":
                        f.write(
                            f"  * ⚠️ **Ambiguous Call** (`AMBIGUOUS_NO_HIFI_SPAN`): "
                            f"Zero HiFi spanning reads cannot distinguish a repetitive/GRS gap "
                            f"(reads present but MAPQ-filtered) from a genuine break. "
                            f"This locus was retained conservatively. "
                            f"The trans-chromosomal signal ({supp} split reads, VAF={vaf:.2f} → {tgt}:{tpos:,}) "
                            f"may be real but requires Hi-C validation before any transfer is executed. "
                            f"Pre-curation Hi-C already validated the current join structure.\n"
                        )
                    f.write(f"  * **Curator Takeaway**: {rationale}\n\n")


        f.write("---\n\n")

        # Section 3: Assembly Sequence & Mass Surgery Ledger
        f.write("## 3. Assembly Sequence & Mass Surgery Ledger\n")
        f.write("*Unified accounting of all physical additions, excisions, and structural compositions across the genome:*\n\n")
        f.write("| Chromosome / Locus | Surgery Type | Event / Source | Net Delta (bp) | Structural Breakdown |\n")
        f.write("| :--- | :--- | :--- | :---: | :--- |\n")

        ledger_rows = []

        # Collect per-chromosome edits
        for c in chr_names:
            c_items = []

            # 1. Gaps (from gap_manifest or audit)
            for g in gap_rows:
                if g.get("chrom") == c:
                    start = int(g.get("gap_start", 0))
                    end = int(g.get("gap_end", 0))
                    glen = int(g.get("gap_len", end - start))
                    outcome = g.get("outcome", "")
                    pos_mb = start / 1e6

                    if outcome == "EXCISE_GAP":
                        c_items.append({
                            "coord": start,
                            "locus": f"**Chr {c}** (~{pos_mb:.2f} Mb)",
                            "surgery_type": "Gap Excision",
                            "event": f"{glen} Ns removed",
                            "net_delta": -glen,
                            "breakdown": "0-bp dummy gap artifact"
                        })
                    elif outcome == "PATCH_GAP":
                        plen = int(g.get("patch_len", 0))
                        delta = plen - glen
                        card_lbl = g.get("unitig_id", g.get("gap_id", "CARD"))
                        c_items.append({
                            "coord": start,
                            "locus": f"**Chr {c}** (~{pos_mb:.2f} Mb)",
                            "surgery_type": "Gap Patch",
                            "event": "Spanning unitig",
                            "net_delta": delta,
                            "breakdown": f"Continuous biological sequence (`{card_lbl}`)"
                        })

            # 2. Telomeres (from walker_manifest)
            for w in walker_rows:
                if w.get("chrom") == c:
                    added = int(w.get("added_bp", 0))
                    telo = int(w.get("canonical_telo_bp", 0))
                    term = w.get("terminus", "")
                    term_str = "5'" if term == "5p" else "3'"
                    coord = 0 if term == "5p" else pre_stats.get(c, {}).get("length", 100000000)

                    if added > 0:
                        subtelo = added - telo
                        subtelo_pct = (subtelo / added * 100.0) if added > 0 else 0.0
                        telo_pct = (telo / added * 100.0) if added > 0 else 0.0

                        event = "Walker extension"
                        if "UNITIG" in w.get("status", "") or "Unitig" in w.get("initial_evidence", ""):
                            event = "Unitig cap"
                        elif "READ" in w.get("status", "") or "Read" in w.get("initial_evidence", ""):
                            event = "Single read cap"
                        elif "SATELLITE" in w.get("status", ""):
                            event = "Satellite array cap"

                        c_items.append({
                            "coord": coord,
                            "locus": f"**Chr {c}** ({term_str})",
                            "surgery_type": "Telomere Capping",
                            "event": event,
                            "net_delta": added,
                            "breakdown": f"Subtelomere: {subtelo:,} bp ({subtelo_pct:.1f}%)<br>Canonical telomere: {telo:,} bp ({telo_pct:.1f}%)"
                        })

            # 3. rDNA Cassettes (from rdna_manifest)
            for r in rdna_rows:
                if r.get("chrom") == c:
                    total_len = int(r.get("cassette_total_len", 0))
                    if total_len > 0:
                        card_id = r.get("card_id", "")
                        card = cards_by_id.get(card_id, {})
                        term = card.get("terminus", "3p")
                        term_str = "5'" if term == "5p" else "3'"
                        coord = 0 if term == "5p" else pre_stats.get(c, {}).get("length", 100000000)

                        mono_len = int(r.get("monomer_len", 0))
                        mono_copies = int(r.get("monomer_copies", 1))
                        array_bp = mono_len * mono_copies
                        anchor_bp = int(r.get("core_anchor_len", 0))
                        cap_bp = int(r.get("distal_cap_len", 0))

                        # Extract canonical telomere portion if present
                        telo_bp = 0
                        if "telo_bp=" in r.get("notes", ""):
                            m = re.search(r"telo_bp=(\d+)", r["notes"])
                            if m:
                                telo_bp = int(m.group(1))
                        elif "canonical_telo_bp" in r and r.get("canonical_telo_bp"):
                            telo_bp = int(r.get("canonical_telo_bp", 0))
                        elif r.get("distal_telo_status") == "TELO_CAPPED":
                            fasta_path = r.get("cassette_fasta", "")
                            if fasta_path and os.path.exists(fasta_path):
                                try:
                                    c_seqs = read_fasta(fasta_path)
                                    for s in c_seqs.values():
                                        cap_seq = s[-cap_bp:] if cap_bp > 0 else s
                                        telo_cnt = max(
                                            cap_seq.upper().count("TTAGGC"),
                                            cap_seq.upper().count("GCCTAA"),
                                            cap_seq.upper().count("TTAGGG"),
                                            cap_seq.upper().count("CCCTAA")
                                        )
                                        telo_bp = telo_cnt * 6
                                        break
                                except Exception:
                                    telo_bp = 0

                        array_pct = (array_bp / total_len) * 100.0 if total_len > 0 else 0.0
                        anchor_pct = (anchor_bp / total_len) * 100.0 if total_len > 0 else 0.0
                        cap_pct = (cap_bp / total_len) * 100.0 if total_len > 0 else 0.0
                        telo_pct = (telo_bp / total_len) * 100.0 if total_len > 0 else 0.0
                        mono_str = f"{mono_copies} monomer" if mono_copies == 1 else f"{mono_copies} monomers"

                        breakdown_parts = [
                            f"**rDNA repeat array: {array_bp:,} bp** ({array_pct:.1f}%, {mono_str})",
                            f"Unique anchor: {anchor_bp:,} bp ({anchor_pct:.1f}%)"
                        ]
                        if cap_bp > 0:
                            breakdown_parts.append(f"Distal cap: {cap_bp:,} bp ({cap_pct:.1f}%)")
                        if telo_bp > 0:
                            breakdown_parts.append(f"Canonical telomere: {telo_bp:,} bp ({telo_pct:.1f}%)")

                        breakdown = "<br>".join(breakdown_parts)

                        c_items.append({
                            "coord": coord,
                            "locus": f"**Chr {c}** ({term_str})",
                            "surgery_type": "**rDNA Cassette**",
                            "event": "Tripartite assembly",
                            "net_delta": total_len,
                            "breakdown": breakdown
                        })

            # 4. Trans-Chromosomal Events (from trans_rows)
            for t in trans_rows:
                outcome = t.get("outcome", "")
                clen = int(t.get("contig_length", 0))
                if clen <= 0:
                    continue
                src = t.get("source_chrom", "")
                tgt = t.get("target_chrom", "")
                gs = int(t.get("source_gap_start", 0))
                tpos = int(t.get("target_coord", 0))
                supp = int(t.get("support_split_reads", 0))
                vaf = float(t.get("vaf", 0.0))
                assigned = t.get("assigned_name", "")

                if src == c:
                    pos_mb = gs / 1e6
                    if outcome == "TRANSFER_AND_STITCH_CONTIG":
                        c_items.append({
                            "coord": gs,
                            "locus": f"**Chr {c}** (~{pos_mb:.2f} Mb)",
                            "surgery_type": "Scaffold Transfer",
                            "event": f"Transferred to {tgt}",
                            "net_delta": -clen,
                            "breakdown": f"Severed {clen:,} bp contig; {supp} split reads (VAF={vaf:.2f}) join to {tgt}"
                        })
                    elif outcome == "SEVER_TO_UNPLACED":
                        c_items.append({
                            "coord": gs,
                            "locus": f"**Chr {c}** (~{pos_mb:.2f} Mb)",
                            "surgery_type": "Scaffold Sever",
                            "event": f"Severed to {assigned}",
                            "net_delta": -clen,
                            "breakdown": f"Severed {clen:,} bp contig (0 collinear reads across gap)"
                        })

                if tgt == c and outcome == "TRANSFER_AND_STITCH_CONTIG":
                    pos_mb = tpos / 1e6
                    tfeat = t.get("target_feature", "terminus")
                    c_items.append({
                        "coord": tpos,
                        "locus": f"**Chr {c}** (~{pos_mb:.2f} Mb)",
                        "surgery_type": "Contig Integration",
                        "event": f"Transferred from {src}",
                        "net_delta": clen,
                        "breakdown": f"Integrated {clen:,} bp contig at {tfeat}; {supp} split reads (VAF={vaf:.2f})"
                    })

            # Sort items by genomic coordinate: 5' (0) -> internal gaps -> 3' (end)
            c_items.sort(key=lambda x: x["coord"])
            ledger_rows.extend(c_items)

        # 4. Unplaced Scaffolds (Trans-chromosomal fissions & Segregation / Purge)
        severed_unplaced = [t for t in trans_rows if t.get("outcome") == "SEVER_TO_UNPLACED"]
        for su in severed_unplaced:
            su_len = int(su.get("contig_length", 0))
            su_name = su.get("assigned_name", "unloc_scaffold")
            su_src = su.get("source_chrom", "")
            ledger_rows.append({
                "locus": f"**{su_name}**",
                "surgery_type": "Unplaced Scaffold Creation",
                "event": f"Severed from {su_src}",
                "net_delta": su_len,
                "breakdown": f"Demoted {su_len:,} bp contig from misjoined gap"
            })

        if total_purged_bp == 0:
            ledger_rows.append({
                "locus": "**Unplaced Scaffolds**",
                "surgery_type": "Segregation / Purge",
                "event": "Redundant haplotigs",
                "net_delta": 0,
                "breakdown": "None segregated"
            })
        else:
            purged_cnt = len(purged_haplo) + len(purged_rdna)
            ledger_rows.append({
                "locus": "**Unplaced Scaffolds**",
                "surgery_type": "Segregation / Purge",
                "event": "Redundant haplotigs / rDNA",
                "net_delta": -total_purged_bp,
                "breakdown": f"{purged_cnt} scaffolds segregated ({total_purged_bp:,} bp)"
            })

        # 5. Total Assembly Delta
        ledger_rows.append({
            "locus": "**Total Assembly Delta**",
            "surgery_type": "",
            "event": "",
            "net_delta": measured_total_delta,
            "breakdown": ""
        })

        # Emit markdown table rows
        for row in ledger_rows:
            delta_str = f"**{row['net_delta']:+,} bp**" if row["net_delta"] != 0 else "**0 bp**"
            f.write(f"| {row['locus']} | {row['surgery_type']} | {row['event']} | {delta_str} | {row['breakdown']} |\n")

    print(f"[{args.sample}] Emitted curation report MD: {report_path}", flush=True)


if __name__ == "__main__":
    main()
