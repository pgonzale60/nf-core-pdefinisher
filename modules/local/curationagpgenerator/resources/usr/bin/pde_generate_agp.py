#!/usr/bin/env python3
"""
pde_generate_agp.py - Module 4: Master Curated AGP Generator

Compiles the Curation Blueprint (Module 2) and physical patches (Module 3 Engines A, B, C)
into an authoritative, NCBI/INSDC-compliant AGP v2.1 file and component FASTA bundle.

Uses the Stateless Block Accumulator pattern:
- Decomposes each chromosome into an ordered sequence of discrete building blocks:
  (Baseline slices, Gap patches, Telomere extensions, rDNA cassettes, Unresolved gaps)
- Accumulates block lengths strictly sequentially:
  Eliminates dynamic coordinate shifting math, off-by-one errors, and negative lengths.
- Segregates unplaced scaffolds (PURGE_RDNA_FRAGMENT, PURGE_HAPLOTIG).
- Emits:
  1. ${sample}.curated.agp
  2. ${sample}.components.fasta
  3. ${sample}.haplotigs.fasta
  4. ${sample}.purged_rdna.fasta
  5. ${sample}.curation_audit.tsv
  (Optional: ${sample}.curated.fasta directly via --synthesize-fasta)
"""

import argparse
import gzip
import json
import os
import re
import sys
from collections import defaultdict, namedtuple

# Named tuple representing an invariant assembly building block
AGPBlock = namedtuple("AGPBlock", ["type", "component_id", "length", "orientation", "source_seq"])


def safe_int(val, default=0):
    if val is None:
        return default
    if isinstance(val, (int, float)):
        return int(val)
    cleaned = re.sub(r"[^\d\-]", "", str(val))
    return int(cleaned) if cleaned else default


def strip_terminal_gaps(blocks, components_pool):
    """
    Strips leading and trailing N blocks (type == "N") and trims terminal N characters
    from boundary W blocks to guarantee no chromosome or scaffold begins or ends
    with an N-gap (e.g. when cleaved at PDE breakpoints or scaffolding gaps).
    """
    changed = True
    while changed and blocks:
        changed = False
        # Strip leading N blocks
        while blocks and blocks[0].type == "N":
            blocks.pop(0)
            changed = True
        # Trim leading Ns from first W block
        if blocks and blocks[0].type == "W" and blocks[0].source_seq:
            seq = blocks[0].source_seq
            trimmed_seq = seq.lstrip("Nn")
            if len(trimmed_seq) < len(seq):
                cid = blocks[0].component_id
                if not trimmed_seq:
                    blocks.pop(0)
                    components_pool.pop(cid, None)
                else:
                    blocks[0] = AGPBlock(
                        type="W", component_id=cid, length=len(trimmed_seq),
                        orientation=blocks[0].orientation, source_seq=trimmed_seq
                    )
                    components_pool[cid] = trimmed_seq
                changed = True

        # Strip trailing N blocks
        while blocks and blocks[-1].type == "N":
            blocks.pop()
            changed = True
        # Trim trailing Ns from last W block
        if blocks and blocks[-1].type == "W" and blocks[-1].source_seq:
            seq = blocks[-1].source_seq
            trimmed_seq = seq.rstrip("Nn")
            if len(trimmed_seq) < len(seq):
                cid = blocks[-1].component_id
                if not trimmed_seq:
                    blocks.pop()
                    components_pool.pop(cid, None)
                else:
                    blocks[-1] = AGPBlock(
                        type="W", component_id=cid, length=len(trimmed_seq),
                        orientation=blocks[-1].orientation, source_seq=trimmed_seq
                    )
                    components_pool[cid] = trimmed_seq
                changed = True
    return blocks


def run_self_test():
    """
    Validates terminal N-gap stripping and AGP block integrity under breakpoint cleavage.
    """
    comp_pool = {
        "c_lead_trail_n": "NNNNATCGATCGNNNN",
        "c_mid": "GGCCAATT",
        "c_pure_n": "NNNNNNNN",
    }
    test_blocks = [
        AGPBlock(type="N", component_id="gap_start", length=100, orientation=None, source_seq="N" * 100),
        AGPBlock(type="W", component_id="c_lead_trail_n", length=16, orientation="+", source_seq=comp_pool["c_lead_trail_n"]),
        AGPBlock(type="N", component_id="gap_mid", length=50, orientation=None, source_seq="N" * 50),
        AGPBlock(type="W", component_id="c_mid", length=8, orientation="+", source_seq=comp_pool["c_mid"]),
        AGPBlock(type="N", component_id="gap_end", length=200, orientation=None, source_seq="N" * 200),
    ]

    stripped = strip_terminal_gaps(test_blocks, comp_pool)
    assert len(stripped) == 3, f"Expected 3 blocks, got {len(stripped)}"
    assert stripped[0].type == "W" and stripped[0].component_id == "c_lead_trail_n"
    assert stripped[0].length == 12, f"Expected length 12 after leading N trim, got {stripped[0].length}"
    assert comp_pool["c_lead_trail_n"] == "ATCGATCGNNNN", f"Expected 'ATCGATCGNNNN', got {comp_pool['c_lead_trail_n']}"
    assert stripped[1].type == "N" and stripped[1].component_id == "gap_mid"
    assert stripped[2].type == "W" and stripped[2].component_id == "c_mid"

    # Test trailing W block with trailing Ns followed by pure N block
    comp_pool2 = {
        "c_mid2": "GGCCAATT",
        "c_trail_n": "ATCGNNNN",
        "c_pure_n": "NNNNNNNN",
    }
    test_blocks2 = [
        AGPBlock(type="W", component_id="c_mid2", length=8, orientation="+", source_seq=comp_pool2["c_mid2"]),
        AGPBlock(type="W", component_id="c_trail_n", length=8, orientation="+", source_seq=comp_pool2["c_trail_n"]),
        AGPBlock(type="W", component_id="c_pure_n", length=8, orientation="+", source_seq=comp_pool2["c_pure_n"]),
        AGPBlock(type="N", component_id="gap_end2", length=100, orientation=None, source_seq="N" * 100),
    ]
    stripped2 = strip_terminal_gaps(test_blocks2, comp_pool2)
    assert len(stripped2) == 2, f"Expected 2 blocks, got {len(stripped2)}"
    assert stripped2[0].component_id == "c_mid2"
    assert stripped2[1].component_id == "c_trail_n"
    assert stripped2[1].length == 4, f"Expected length 4, got {stripped2[1].length}"
    assert comp_pool2["c_trail_n"] == "ATCG"
    assert "c_pure_n" not in comp_pool2

    print("CURATIONAGPGENERATOR self-test PASS")
    sys.exit(0)


def parse_args():
    parser = argparse.ArgumentParser(description="Module 4: Master Curated AGP Generator")
    parser.add_argument("--blueprint", required=False, help="Module 2 Curation Blueprint JSON")
    parser.add_argument("--ref", required=False, help="Reference FASTA (chromosomes + unplaced scaffolds)")
    parser.add_argument("--gap-manifest", required=False, default=None, help="Engine A patches manifest TSV")
    parser.add_argument("--walker-manifest", required=False, default=None, help="Engine B telomere manifest TSV")
    parser.add_argument("--rdna-manifest", required=False, default=None, help="Engine C rDNA manifest TSV")
    parser.add_argument("--trans-manifest", required=False, default=None, help="Engine D trans-chromosomal manifest TSV")
    parser.add_argument("--pde-segments", required=False, default=None, help="Input assembly pde_segments.bed")
    parser.add_argument("--red-bed", required=False, default=None, help="Input assembly RED repeat mask BED")
    parser.add_argument("--sample", required=False, help="Sample/strain prefix")
    parser.add_argument("--outdir", required=False, help="Output directory")
    parser.add_argument("--synthesize-fasta", action="store_true", help="Synthesize curated.fasta directly using pysam/internal engine")
    parser.add_argument("--self-test", action="store_true", help="Run internal self-test and exit")
    args = parser.parse_args()
    if not args.self_test:
        if not (args.blueprint and args.ref and args.sample and args.outdir):
            parser.error("--blueprint, --ref, --sample, and --outdir are required when not running --self-test")
    return args


def load_bed_intervals(bed_path, default_label=""):
    """Loads intervals from a BED file into {chrom: [(start, end, label), ...]}."""
    intervals = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return intervals
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                c = parts[0]
                try:
                    s = int(parts[1])
                    e = int(parts[2])
                except ValueError:
                    continue
                lbl = parts[3] if len(parts) >= 4 else default_label
                intervals[c].append((s, e, lbl))
    for c in intervals:
        intervals[c].sort(key=lambda x: (x[0], x[1]))
    return intervals


def merge_curated_repeats(raw_intervals):
    """Sorts and merges overlapping/adjacent repeat intervals."""
    if not raw_intervals:
        return []
    raw_intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []
    for s, e, lbl in raw_intervals:
        if s >= e:
            continue
        if not merged:
            merged.append([s, e, lbl])
        elif s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e, lbl])
    return [(s, e, lbl) for s, e, lbl in merged]


def merge_curated_pde_segments(raw_intervals, total_len):
    """Sorts, resolves, and merges contiguous GRS/Core segments, filling novel/uncovered regions with putative_core."""
    if not raw_intervals:
        return [(0, total_len, "putative_core")] if total_len > 0 else []
    raw_intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []
    for s, e, cls in raw_intervals:
        if s >= e:
            continue
        if merged and merged[-1][2] == cls and merged[-1][1] == s:
            merged[-1][1] = e
        elif merged and merged[-1][1] > s:
            if e > merged[-1][1]:
                if merged[-1][2] == cls:
                    merged[-1][1] = e
                else:
                    merged.append([merged[-1][1], e, cls])
        else:
            merged.append([s, e, cls])
    if not merged and total_len > 0:
        merged.append([0, total_len, "putative_core"])
    return [(s, e, cls) for s, e, cls in merged]


def load_fasta(fasta_path):
    """Loads a FASTA file (plain or gzipped) into a dictionary of {name: seq}."""
    seqs = {}
    if not fasta_path or not os.path.exists(fasta_path):
        return seqs
    opener = gzip.open if fasta_path.endswith(".gz") else open
    mode = "rt" if fasta_path.endswith(".gz") else "r"
    current_name = None
    current_seq = []
    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    seqs[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            seqs[current_name] = "".join(current_seq)
    return seqs


def revcomp(seq):
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def parse_tsv_manifest(tsv_path):
    """Parses a TSV manifest with headers into a list of dicts."""
    rows = []
    if not tsv_path or not os.path.exists(tsv_path):
        return rows
    with open(tsv_path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
            else:
                row = dict(zip(header, parts))
                rows.append(row)
    return rows


def main():
    args = parse_args()
    if args.self_test:
        run_self_test()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load Reference Genome
    ref_seqs = load_fasta(args.ref)
    if not ref_seqs:
        sys.stderr.write(f"ERROR: Failed to load reference FASTA: {args.ref}\n")
        sys.exit(1)

    # 2. Load Module 2 Blueprint
    if not os.path.exists(args.blueprint):
        sys.stderr.write(f"ERROR: Blueprint not found: {args.blueprint}\n")
        sys.exit(1)
    with open(args.blueprint) as f:
        blueprint = json.load(f)

    # Load input segmentation (GRS vs Core)
    input_segs_by_chr = load_bed_intervals(args.pde_segments, default_label="putative_core")
    if not input_segs_by_chr and "input_segments" in blueprint:
        for seg in blueprint["input_segments"]:
            c = seg.get("chrom", "")
            s = safe_int(seg.get("start", 0))
            e = safe_int(seg.get("end", 0))
            cls = seg.get("classification", "putative_core")
            if c:
                input_segs_by_chr[c].append((s, e, cls))
    if not input_segs_by_chr:
        for c in sorted(ref_seqs.keys()):
            clen = len(ref_seqs[c])
            cls = "putative_GRS" if "grs" in c.lower() else "putative_core"
            input_segs_by_chr[c].append((0, clen, cls))

    # Load RED repeat mask
    red_repeats_by_chr = load_bed_intervals(args.red_bed, default_label="repeat")

    curated_segments = []
    curated_repeats = []

    # 3. Load Module 3 Manifests
    gap_rows = parse_tsv_manifest(args.gap_manifest)
    walker_rows = parse_tsv_manifest(args.walker_manifest)
    rdna_rows = parse_tsv_manifest(args.rdna_manifest)
    trans_rows = parse_tsv_manifest(args.trans_manifest)

    # Index Engine A Patches: gap_id -> dict
    patches_by_id = {}
    gap_dir = os.path.dirname(os.path.abspath(args.gap_manifest)) if args.gap_manifest else ""
    for r in gap_rows:
        gid = r.get("gap_id", "")
        if gid:
            fa_p = r.get("patch_fasta", "")
            if fa_p and fa_p != "NONE" and not os.path.isabs(fa_p) and not os.path.exists(fa_p):
                alt_p = os.path.join(gap_dir, fa_p)
                alt_p_base = os.path.join(gap_dir, os.path.basename(fa_p))
                alt_p_cwd = os.path.basename(fa_p)
                if os.path.exists(alt_p):
                    r["patch_fasta"] = alt_p
                elif os.path.exists(alt_p_base):
                    r["patch_fasta"] = alt_p_base
                elif os.path.exists(alt_p_cwd):
                    r["patch_fasta"] = alt_p_cwd
            patches_by_id[gid] = r

    # Index Engine B Telomeres: (chrom, terminus) -> dict
    walker_by_end = {}
    walker_dir = os.path.dirname(os.path.abspath(args.walker_manifest)) if args.walker_manifest else ""
    for r in walker_rows:
        c = r.get("chrom", "")
        term = r.get("terminus", "")
        if c and term:
            fa_p = r.get("extension_fasta", "")
            if fa_p and fa_p != "NONE" and not os.path.isabs(fa_p) and not os.path.exists(fa_p):
                alt_p = os.path.join(walker_dir, fa_p)
                alt_p_base = os.path.join(walker_dir, os.path.basename(fa_p))
                alt_p_cwd = os.path.basename(fa_p)
                if os.path.exists(alt_p):
                    r["extension_fasta"] = alt_p
                elif os.path.exists(alt_p_base):
                    r["extension_fasta"] = alt_p_base
                elif os.path.exists(alt_p_cwd):
                    r["extension_fasta"] = alt_p_cwd
            walker_by_end[(c, term)] = r

    # Index Engine C rDNA: (chrom, card_id) -> dict
    rdna_by_card = {}
    rdna_dir = os.path.dirname(os.path.abspath(args.rdna_manifest)) if args.rdna_manifest else ""
    for r in rdna_rows:
        cid = r.get("card_id", "")
        c = r.get("chrom", "")
        if cid:
            fa_p = r.get("cassette_fasta", "")
            if fa_p and fa_p != "NONE" and not os.path.isabs(fa_p) and not os.path.exists(fa_p):
                alt_p = os.path.join(rdna_dir, fa_p)
                alt_p_base = os.path.join(rdna_dir, os.path.basename(fa_p))
                alt_p_cwd = os.path.basename(fa_p)
                if os.path.exists(alt_p):
                    r["cassette_fasta"] = alt_p
                elif os.path.exists(alt_p_base):
                    r["cassette_fasta"] = alt_p_base
                elif os.path.exists(alt_p_cwd):
                    r["cassette_fasta"] = alt_p_cwd
            rdna_by_card[cid] = r
            rdna_by_card[(c, cid)] = r

    components_pool = {}  # component_id -> seq

    # Index Engine D Trans-Chromosomal Links:
    trans_by_src = defaultdict(list)
    trans_by_tgt = defaultdict(list)
    unplaced_severed = []

    for r in trans_rows:
        outcome = r.get("outcome", "")
        src = r.get("source_chrom", "")
        tgt = r.get("target_chrom", "")
        cstart = safe_int(r.get("contig_start", 0))
        cend = safe_int(r.get("contig_end", 0))
        cid = r.get("contig_id", f"{src}_{cstart}_{cend}")
        clen = safe_int(r.get("contig_length", cend - cstart))

        if outcome in ["TRANSFER_AND_STITCH_CONTIG", "SEVER_TO_UNPLACED"]:
            if outcome == "TRANSFER_AND_STITCH_CONTIG" and (src == tgt or not tgt):
                # Safety Guard: Never execute self-targeting contig transfer
                continue
            if src in ref_seqs and clen > 0:
                components_pool[cid] = ref_seqs[src][cstart:cend]
            trans_by_src[src].append(r)
            if outcome == "TRANSFER_AND_STITCH_CONTIG" and tgt:
                trans_by_tgt[tgt].append(r)
            elif outcome == "SEVER_TO_UNPLACED":
                unplaced_severed.append(r)

    # Collect Blueprint Actions by Domain
    cards_by_chrom = defaultdict(list)
    scaffold_actions = {}
    scaffold_gap_cards_by_id = {}
    scaffold_term_cards_by_end = defaultdict(list)
    cassette_scaffolds_by_end = defaultdict(list)
    materialized_scaffolds = set()

    all_cards = (
        blueprint.get("action_cards", []) +
        blueprint.get("curation_cards", []) +
        blueprint.get("telomere_cards", [])
    )

    seen_cards = set()
    for card in all_cards:
        cid = card.get("card_id", "")
        if cid in seen_cards:
            continue
        seen_cards.add(cid)

        dom = card.get("domain", "")
        act = card.get("action", "")
        c = card.get("chrom", "")
        scaff = card.get("scaffold", "")

        if dom == "DOMAIN_SCAFFOLD" or act.startswith("PURGE_"):
            target = scaff or c
            if target:
                scaffold_actions[target] = act
            if act == "INTEGRATE_UNPLACED_REPEAT" and card.get("target_gap_id"):
                scaffold_gap_cards_by_id[card["target_gap_id"]] = card
            elif act == "INTEGRATE_INTO_CHROMOSOME" and c and card.get("terminus"):
                if "ext_qstart" in card and "ext_qend" in card:
                    scaffold_term_cards_by_end[(c, card["terminus"])].append(card)
                elif scaff:
                    cassette_scaffolds_by_end[(c, card["terminus"])].append(scaff)
        elif c:
            cards_by_chrom[c].append(card)
            t_arm = card.get("terminus")
            if t_arm:
                for a_scaff in card.get("evidence", {}).get("associated_scaffolds", []):
                    if a_scaff not in cassette_scaffolds_by_end[(c, t_arm)]:
                        cassette_scaffolds_by_end[(c, t_arm)].append(a_scaff)

    for row in blueprint.get("scaffold_dispositions", []) or []:
        scaff = row.get("scaffold")
        act = row.get("action")
        if scaff and act:
            scaffold_actions[scaff] = act
            if act == "RETIRED_ABSORBED" and row.get("chromosome") and row.get("terminus"):
                cassette_scaffolds_by_end[(row["chromosome"], row["terminus"])].append(scaff)

    # File paths
    agp_path = os.path.join(args.outdir, f"{args.sample}.curated.agp")
    comp_path = os.path.join(args.outdir, f"{args.sample}.components.fasta")
    haplo_path = os.path.join(args.outdir, f"{args.sample}.haplotigs.fasta")
    purged_rdna_path = os.path.join(args.outdir, f"{args.sample}.purged_rdna.fasta")
    audit_path = os.path.join(args.sample and args.outdir, f"{args.sample}.curation_audit.tsv")

    audit_rows = []
    agp_lines = []

    # Write AGP Header
    agp_lines.append("## AGP-version 2.1")
    agp_lines.append(f"# Organism: {args.sample}")
    agp_lines.append("# Pipeline: Automated PDE Curation Pipeline (Module 4 Surgery)")

    # Identify Chromosomes vs Unplaced Scaffolds
    # Primary chromosomes typically don't have '_unloc' or 'scaffold'
    chrom_names = [k for k in ref_seqs.keys() if "_unloc" not in k and not k.startswith("scaffold_") and k != "scaffold_MT"]
    unplaced_names = [k for k in ref_seqs.keys() if k not in chrom_names]

    # --- PROCESS CHROMOSOMES ---
    for chrom in sorted(chrom_names):
        chrom_seq = ref_seqs[chrom]
        chrom_len = len(chrom_seq)
        chrom_cards = cards_by_chrom.get(chrom, [])

        # Build list of Edits for this chromosome
        edits = []

        # 1. 5' Telomere Extension
        w_5p = walker_by_end.get((chrom, "5p"))
        extended_5p = False
        if w_5p:
            status = w_5p.get("status", "")
            added = safe_int(w_5p.get("added_bp", 0))
            fasta_p = w_5p.get("extension_fasta", "")
            if added > 0 and fasta_p and fasta_p != "NONE" and os.path.exists(fasta_p) and os.path.getsize(fasta_p) > 0:
                ext_seqs = load_fasta(fasta_p)
                ext_seq = list(ext_seqs.values())[0] if ext_seqs else ""
                if ext_seq:
                    comp_id = f"{chrom}_5p_telo_ext"
                    components_pool[comp_id] = ext_seq
                    edits.append({
                        "type": "TELO_5P",
                        "start": 0, "end": 0,
                        "comp_id": comp_id,
                        "length": len(ext_seq),
                        "card_id": f"TELO-{chrom}-5p",
                        "action": status
                    })
                    extended_5p = True
                    for c_scaff in cassette_scaffolds_by_end.get((chrom, "5p"), []):
                        materialized_scaffolds.add(c_scaff)
            elif status == "ALREADY_CAPPED" or w_5p.get("action") == "SKIP_WALKER":
                audit_rows.append({
                    "sample": args.sample, "chrom": chrom, "card_id": f"TELO-{chrom}-5p",
                    "action": "SKIP_WALKER", "edit_type": "ALREADY_CAPPED",
                    "orig_start": 0, "orig_end": 0, "delta_bp": 0
                })
        if not extended_5p and scaffold_term_cards_by_end.get((chrom, "5p")):
            t_card = scaffold_term_cards_by_end[(chrom, "5p")][0]
            sc_name = t_card.get("scaffold", "")
            if sc_name in ref_seqs and sc_name not in materialized_scaffolds:
                qs = safe_int(t_card.get("ext_qstart", 0))
                qe = safe_int(t_card.get("ext_qend", 0))
                ori = t_card.get("orientation", "+")
                raw_slice = ref_seqs[sc_name][qs:qe]
                if raw_slice:
                    ext_seq = revcomp(raw_slice) if ori == "-" else raw_slice
                    comp_id = f"{chrom}_5p_{sc_name}_ext"
                    components_pool[comp_id] = ext_seq
                    edits.append({
                        "type": "TELO_5P",
                        "start": 0, "end": 0,
                        "comp_id": comp_id,
                        "length": len(ext_seq),
                        "card_id": t_card.get("card_id", f"SCAFF-EXT-{chrom}-5p"),
                        "action": "INTEGRATE_INTO_CHROMOSOME"
                    })
                    materialized_scaffolds.add(sc_name)

        # 2. Trans-Chromosomal Severs (Source Chromosome Excision)
        for tr in trans_by_src.get(chrom, []):
            gs = safe_int(tr.get("source_gap_start", 0))
            ge = safe_int(tr.get("source_gap_end", 0))
            cstart = safe_int(tr.get("contig_start", 0))
            cend = safe_int(tr.get("contig_end", 0))
            cid = tr.get("contig_id", f"{chrom}_{cstart}_{cend}")
            excise_start = min(cstart, gs)
            excise_end = max(cend, ge)
            edits.append({
                "type": "EXCISE_SEVERED_CONTIG",
                "start": excise_start, "end": excise_end,
                "comp_id": cid,
                "length": 0,
                "card_id": tr.get("card_id", "TRANS"),
                "action": tr.get("outcome", "SEVER")
            })

        # 3. Trans-Chromosomal Integrations (Target Chromosome Insertion)
        for tr in trans_by_tgt.get(chrom, []):
            tpos = safe_int(tr.get("target_coord", 0))
            tstrand = tr.get("target_strand", "+")
            cid = tr.get("contig_id", "")
            clen = safe_int(tr.get("contig_length", 0))
            tfeat = tr.get("target_feature", "")
            if clen > 0 and cid:
                if tpos <= 10000 or tfeat == "5_prime_terminus":
                    edits.append({
                        "type": "TRANS_INSERT_5P",
                        "start": 0, "end": 0,
                        "comp_id": cid,
                        "length": clen,
                        "orientation": tstrand,
                        "card_id": tr.get("card_id", "TRANS"),
                        "action": "TRANSFER_AND_STITCH_CONTIG"
                    })
                elif abs(tpos - chrom_len) <= 10000 or tfeat == "3_prime_terminus":
                    edits.append({
                        "type": "TRANS_INSERT_3P",
                        "start": chrom_len, "end": chrom_len,
                        "comp_id": cid,
                        "length": clen,
                        "orientation": tstrand,
                        "card_id": tr.get("card_id", "TRANS"),
                        "action": "TRANSFER_AND_STITCH_CONTIG"
                    })
                else:
                    edits.append({
                        "type": "TRANS_INSERT_GAP",
                        "start": tpos, "end": tpos,
                        "comp_id": cid,
                        "length": clen,
                        "orientation": tstrand,
                        "card_id": tr.get("card_id", "TRANS"),
                        "action": "TRANSFER_AND_STITCH_CONTIG"
                    })

        # 4. 3' Telomere Extension
        w_3p = walker_by_end.get((chrom, "3p"))
        extended_3p = False
        if w_3p:
            status = w_3p.get("status", "")
            added = safe_int(w_3p.get("added_bp", 0))
            fasta_p = w_3p.get("extension_fasta", "")
            if added > 0 and fasta_p and fasta_p != "NONE" and os.path.exists(fasta_p) and os.path.getsize(fasta_p) > 0:
                ext_seqs = load_fasta(fasta_p)
                ext_seq = list(ext_seqs.values())[0] if ext_seqs else ""
                if ext_seq:
                    comp_id = f"{chrom}_3p_telo_ext"
                    components_pool[comp_id] = ext_seq
                    edits.append({
                        "type": "TELO_3P",
                        "start": chrom_len, "end": chrom_len,
                        "comp_id": comp_id,
                        "length": len(ext_seq),
                        "card_id": f"TELO-{chrom}-3p",
                        "action": status
                    })
                    extended_3p = True
                    for c_scaff in cassette_scaffolds_by_end.get((chrom, "3p"), []):
                        materialized_scaffolds.add(c_scaff)
            elif status == "ALREADY_CAPPED" or w_3p.get("action") == "SKIP_WALKER":
                audit_rows.append({
                    "sample": args.sample, "chrom": chrom, "card_id": f"TELO-{chrom}-3p",
                    "action": "SKIP_WALKER", "edit_type": "ALREADY_CAPPED",
                    "orig_start": chrom_len, "orig_end": chrom_len, "delta_bp": 0
                })
        if not extended_3p and scaffold_term_cards_by_end.get((chrom, "3p")):
            t_card = scaffold_term_cards_by_end[(chrom, "3p")][0]
            sc_name = t_card.get("scaffold", "")
            if sc_name in ref_seqs and sc_name not in materialized_scaffolds:
                qs = safe_int(t_card.get("ext_qstart", 0))
                qe = safe_int(t_card.get("ext_qend", 0))
                ori = t_card.get("orientation", "+")
                raw_slice = ref_seqs[sc_name][qs:qe]
                if raw_slice:
                    ext_seq = revcomp(raw_slice) if ori == "-" else raw_slice
                    comp_id = f"{chrom}_3p_{sc_name}_ext"
                    components_pool[comp_id] = ext_seq
                    edits.append({
                        "type": "TELO_3P",
                        "start": chrom_len, "end": chrom_len,
                        "comp_id": comp_id,
                        "length": len(ext_seq),
                        "card_id": t_card.get("card_id", f"SCAFF-EXT-{chrom}-3p"),
                        "action": "INTEGRATE_INTO_CHROMOSOME"
                    })
                    materialized_scaffolds.add(sc_name)

        # 5. Gap Cards and Structural Actions
        for card in chrom_cards:
            cid = card.get("card_id", "")
            act = card.get("action", "")
            start = card.get("start") if card.get("start") is not None else card.get("start_pos")
            end = card.get("end") if card.get("end") is not None else card.get("end_pos")

            if start is None or end is None:
                term = card.get("terminus")
                if "TERMINUS" in act or "RDNA" in act or card.get("domain") == "DOMAIN_RDNA":
                    if term == "5p":
                        start, end = 0, 0
                    elif term == "3p":
                        start, end = chrom_len, chrom_len
                    else:
                        continue
                else:
                    continue

            start = safe_int(start)
            end = safe_int(end)

            if act in ("EXCISE_DUMMY_GAP", "EXCISE_GAP"):
                edits.append({
                    "type": "EXCISE_GAP",
                    "start": start, "end": end,
                    "comp_id": None,
                    "length": 0,
                    "card_id": cid,
                    "action": act
                })
            elif act == "INVERT_AND_STITCH":
                inv_start = safe_int(card.get("inv_start", start))
                inv_end = safe_int(card.get("inv_end", end))
                edits.append({
                    "type": "INVERT",
                    "start": start, "end": end,
                    "inv_start": inv_start, "inv_end": inv_end,
                    "comp_id": f"{chrom}_{inv_start}_{inv_end}_inv",
                    "length": inv_end - inv_start,
                    "card_id": cid,
                    "action": act
                })
            elif act == "CLOSE_GAP_COLINEAR":
                # Reads bridge this gap collinearly with zero missing sequence: stitch the flanks
                # directly. Physically the same edit as EXCISE_DUMMY_GAP, and it needs no patch
                # from Module 3. Without this branch these gaps were left as Ns, because no engine
                # emits a patch for them and nothing else recognised the action.
                edits.append({
                    "type": "EXCISE_GAP",
                    "start": start, "end": end,
                    "comp_id": None,
                    "length": 0,
                    "card_id": cid,
                    "action": act
                })
            elif act == "PRUNE_CHIMERIC_TERMINUS":
                edits.append({
                    "type": "PRUNE_CHIMERIC_TERMINUS",
                    "start": start, "end": end,
                    "comp_id": None,
                    "length": 0,
                    "card_id": cid,
                    "action": act
                })
            elif act in ("REQUIRES_UNITIG_PATCH", "EXCISE_OVERLAP_ARTIFACT", "INSERT_SEQUENCE") or "PATCH_GAP" in act:
                patch_info = patches_by_id.get(cid)
                if patch_info and patch_info.get("outcome") in ["EXCISE_OVERLAP_ARTIFACT", "EXCISE_DUMMY_GAP"]:
                    notes = patch_info.get("notes", "")
                    m = re.search(r"excise_range=(\d+)-(\d+)", notes)
                    if m:
                        e_start = int(m.group(1))
                        e_end = int(m.group(2))
                    else:
                        e_start = start
                        e_end = end
                    edits.append({
                        "type": "EXCISE_GAP",
                        "start": e_start,
                        "end": e_end,
                        "comp_id": None,
                        "length": 0,
                        "card_id": cid,
                        "action": "EXCISE_OVERLAP_ARTIFACT"
                    })
                    continue
                elif patch_info and patch_info.get("outcome") == "PATCH_GAP_INSERTION":
                    p_fasta = patch_info.get("patch_fasta", "")
                    if p_fasta and p_fasta != "NONE" and os.path.exists(p_fasta) and os.path.getsize(p_fasta) > 0:
                        p_seqs = load_fasta(p_fasta)
                        p_seq = p_seqs.get(cid) or p_seqs.get(patch_info.get("unitig_id")) or (list(p_seqs.values())[0] if p_seqs else "")
                        if p_seq:
                            comp_id = f"{cid}_patch"
                            components_pool[comp_id] = p_seq
                            edits.append({
                                "type": "GAP_PATCH",
                                "start": start, "end": end,
                                "comp_id": comp_id,
                                "length": len(p_seq),
                                "card_id": cid,
                                "action": act if act == "INSERT_SEQUENCE" else "PATCH_GAP_INSERTION"
                            })
                            u_id = patch_info.get("unitig_id", "")
                            if u_id in ref_seqs:
                                materialized_scaffolds.add(u_id)
                            continue
                # Fallback to Domain 4 Port-Valency Unplaced Contig Integration (`INTEGRATE_UNPLACED_REPEAT`)
                if cid in scaffold_gap_cards_by_id:
                    g_card = scaffold_gap_cards_by_id[cid]
                    sc_name = g_card.get("scaffold", "")
                    if sc_name in ref_seqs and sc_name not in materialized_scaffolds:
                        qs = safe_int(g_card.get("insert_qstart", 0))
                        qe = safe_int(g_card.get("insert_qend", len(ref_seqs[sc_name])))
                        ori = g_card.get("orientation", "+")
                        l_gap = safe_int(g_card.get("left_gap_bp", 0))
                        r_gap = safe_int(g_card.get("right_gap_bp", 0))
                        body = ref_seqs[sc_name][qs:qe]
                        if body:
                            if ori == "-":
                                body = revcomp(body)
                            p_seq = ("N" * l_gap) + body + ("N" * r_gap)
                            comp_id = f"{cid}_{sc_name}_patch"
                            components_pool[comp_id] = p_seq
                            edits.append({
                                "type": "GAP_PATCH",
                                "start": start, "end": end,
                                "comp_id": comp_id,
                                "length": len(p_seq),
                                "card_id": cid,
                                "action": "INTEGRATE_UNPLACED_REPEAT"
                            })
                            materialized_scaffolds.add(sc_name)
                            continue
                if end > start:
                    edits.append({
                        "type": "UNRESOLVED_GAP",
                        "start": start, "end": end,
                        "comp_id": None,
                        "length": end - start,
                        "card_id": cid,
                        "action": "UNRESOLVED_GAP"
                    })
            elif act == "SUBSUMED_BY_RDNA_ARRAY":
                edits.append({
                    "type": "SUBSUMED_BY_RDNA_ARRAY",
                    "start": start, "end": end,
                    "comp_id": None,
                    "length": 0,
                    "card_id": cid,
                    "action": act
                })
            elif "RDNA" in act or card.get("domain") == "DOMAIN_RDNA":
                rdna_info = rdna_by_card.get(cid) or rdna_by_card.get((chrom, cid))
                if rdna_info and rdna_info.get("gatekeeper_result") == "PASS":
                    c_fasta = rdna_info.get("cassette_fasta", "")
                    if c_fasta and c_fasta != "NONE" and os.path.exists(c_fasta) and os.path.getsize(c_fasta) > 0:
                        c_seqs = load_fasta(c_fasta)
                        c_seq = list(c_seqs.values())[0] if c_seqs else ""
                        if c_seq:
                            comp_id = f"{cid}_rdna_cassette"
                            components_pool[comp_id] = c_seq
                            status = rdna_info.get("status", "")
                            distal_telo = rdna_info.get("distal_telo_status", "")
                            term = card.get("terminus", "")
                            if term:
                                is_telo = (status == "CASSETTE_RESOLVED_CAPPED" or distal_telo == "TELO_CAPPED")
                                edit_type = f"TELO_RDNA_CAPPED_{term.upper()}" if is_telo else f"RDNA_UNCAPPED_{term.upper()}"
                            else:
                                edit_type = "RDNA_CASSETTE"

                            orientation = "-" if term == "5p" else "+"

                            edits.append({
                                "type": edit_type,
                                "start": start, "end": end,
                                "comp_id": comp_id,
                                "length": len(c_seq),
                                "orientation": orientation,
                                "card_id": cid,
                                "action": status or "RESOLVE_RDNA_CASSETTE"
                            })

        # Sort edits strictly by start coordinate
        edits.sort(key=lambda x: (x["start"], x["end"]))

        # Build Ordered Block List (Stateless Accumulator Pattern)
        blocks = []
        curr_idx = 0

        for e in edits:
            e["applied_delta"] = 0
            e["orig_start_pos"] = e["start"]
            e["orig_end_pos"] = e["end"]

        # Check 5' Telomere Extension
        telo_5p_edits = [e for e in edits if e["type"] == "TELO_5P"]
        if telo_5p_edits:
            e = telo_5p_edits[0]
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation="+", source_seq=components_pool[e["comp_id"]]
            ))
            e["applied_delta"] = e["length"]
            for extra in telo_5p_edits[1:]:
                extra["applied_delta"] = 0
                extra["note"] = "SKIPPED_DUPLICATE_TELO_5P"

        # Check 5' Translocated Contig Insertions
        trans_5p_edits = [e for e in edits if e["type"] == "TRANS_INSERT_5P"]
        for e in trans_5p_edits:
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation=e.get("orientation", "+"), source_seq=components_pool.get(e["comp_id"], "")
            ))
            e["applied_delta"] = e["length"]

        # Check 5' Terminal rDNA Cassette
        rdna_5p_edits = [e for e in edits if e["type"] in ["TELO_RDNA_CAPPED_5P", "RDNA_UNCAPPED_5P"]]
        if rdna_5p_edits:
            e = rdna_5p_edits[0]
            seq = components_pool[e["comp_id"]]
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation=e.get("orientation", "-"),
                source_seq=revcomp(seq) if e.get("orientation") == "-" else seq
            ))
            consumed_5p = min(e["end"], chrom_len)
            e["applied_delta"] = e["length"] - consumed_5p
            if consumed_5p > curr_idx:
                curr_idx = consumed_5p
            for extra in rdna_5p_edits[1:]:
                extra["applied_delta"] = 0
                extra["note"] = "SKIPPED_DUPLICATE_RDNA_5P"

        internal_edits = [e for e in edits if e["type"] not in [
            "TELO_5P", "TELO_3P", "TRANS_INSERT_5P", "TRANS_INSERT_3P",
            "TELO_RDNA_CAPPED_5P", "RDNA_UNCAPPED_5P", "TELO_RDNA_CAPPED_3P", "RDNA_UNCAPPED_3P",
            "SUBSUMED_BY_RDNA_ARRAY"
        ]]
        for e in internal_edits:
            # Skip edit if already enclosed within an earlier excised range
            if e["end"] < curr_idx or (e["end"] == curr_idx and e["start"] < e["end"]):
                e["applied_delta"] = 0
                e["note"] = "SKIPPED_ENCLOSED_OR_OVERLAPPING"
                continue
            if e["start"] < curr_idx:
                e["start"] = curr_idx

            # Baseline contig slice before edit
            if e["start"] > curr_idx:
                slice_len = e["start"] - curr_idx
                slice_seq = chrom_seq[curr_idx:e["start"]]
                slice_id = f"{chrom}_{curr_idx}_{e['start']}"
                components_pool[slice_id] = slice_seq
                blocks.append(AGPBlock(
                    type="W", component_id=slice_id, length=slice_len,
                    orientation="+", source_seq=slice_seq
                ))

            # Process Edit
            consumed = e["end"] - e["start"]
            if e["type"] in ["EXCISE_GAP", "EXCISE_SEVERED_CONTIG", "PRUNE_CHIMERIC_TERMINUS"]:
                e["applied_delta"] = -consumed
            elif e["type"] == "INVERT":
                inv_start = max(e["start"], e.get("inv_start", e["start"]))
                inv_end = min(e["end"], e.get("inv_end", e["end"]))
                inv_seq = chrom_seq[inv_start:inv_end]
                components_pool[e["comp_id"]] = inv_seq
                blocks.append(AGPBlock(
                    type="W", component_id=e["comp_id"], length=len(inv_seq),
                    orientation="-", source_seq=revcomp(inv_seq)
                ))
                e["applied_delta"] = len(inv_seq) - consumed
            elif e["type"] in ["GAP_PATCH", "RDNA_CASSETTE"] or "RDNA" in e["type"] or e["type"] == "TRANS_INSERT_GAP":
                seq = components_pool[e["comp_id"]]
                gap_match = re.search(r'N{10,}', seq)
                if gap_match:
                    g_start, g_end = gap_match.span()
                    part1 = seq[:g_start]
                    gap_len = g_end - g_start
                    part2 = seq[g_end:]
                    p1_id = f"{e['comp_id']}_p1"
                    p2_id = f"{e['comp_id']}_p2"
                    if part1:
                        components_pool[p1_id] = part1
                        blocks.append(AGPBlock(
                            type="W", component_id=p1_id, length=len(part1),
                            orientation=e.get("orientation", "+"), source_seq=part1
                        ))
                    blocks.append(AGPBlock(
                        type="N", component_id=f"{e['comp_id']}_gap", length=gap_len,
                        orientation=None, source_seq="N" * gap_len
                    ))
                    if part2:
                        components_pool[p2_id] = part2
                        blocks.append(AGPBlock(
                            type="W", component_id=p2_id, length=len(part2),
                            orientation=e.get("orientation", "+"), source_seq=part2
                        ))
                    components_pool.pop(e["comp_id"], None)
                else:
                    blocks.append(AGPBlock(
                        type="W", component_id=e["comp_id"], length=e["length"],
                        orientation=e.get("orientation", "+"), source_seq=seq
                    ))
                e["applied_delta"] = e["length"] - consumed
            elif e["type"] == "UNRESOLVED_GAP":
                blocks.append(AGPBlock(
                    type="N", component_id=None, length=e["length"],
                    orientation=None, source_seq="N" * e["length"]
                ))
                e["applied_delta"] = e["length"] - consumed

            curr_idx = e["end"]

        # Trailing baseline slice
        if curr_idx < chrom_len:
            slice_end = chrom_len
            rdna_3p_edits = [e for e in edits if e["type"] in ["TELO_RDNA_CAPPED_3P", "RDNA_UNCAPPED_3P"]]
            if rdna_3p_edits:
                slice_end = min(slice_end, max(curr_idx, rdna_3p_edits[0]["start"]))
            if curr_idx < slice_end:
                slice_len = slice_end - curr_idx
                slice_seq = chrom_seq[curr_idx:slice_end]
                slice_id = f"{chrom}_{curr_idx}_{slice_end}"
                components_pool[slice_id] = slice_seq
                blocks.append(AGPBlock(
                    type="W", component_id=slice_id, length=slice_len,
                    orientation="+", source_seq=slice_seq
                ))
                curr_idx = slice_end

        # Check 3' Translocated Contig Insertions
        trans_3p_edits = [e for e in edits if e["type"] == "TRANS_INSERT_3P"]
        for e in trans_3p_edits:
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation=e.get("orientation", "+"), source_seq=components_pool.get(e["comp_id"], "")
            ))
            e["applied_delta"] = e["length"]

        # Check 3' Terminal rDNA Cassette
        rdna_3p_edits = [e for e in edits if e["type"] in ["TELO_RDNA_CAPPED_3P", "RDNA_UNCAPPED_3P"]]
        if rdna_3p_edits:
            e = rdna_3p_edits[0]
            seq = components_pool[e["comp_id"]]
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation=e.get("orientation", "+"), source_seq=seq
            ))
            consumed_3p = max(0, chrom_len - max(curr_idx, e["start"]))
            e["applied_delta"] = e["length"] - consumed_3p
            curr_idx = chrom_len
            for extra in rdna_3p_edits[1:]:
                extra["applied_delta"] = 0
                extra["note"] = "SKIPPED_DUPLICATE_RDNA_3P"

        # Check 3' Telomere Extension
        telo_3p_edits = [e for e in edits if e["type"] == "TELO_3P"]
        if telo_3p_edits:
            e = telo_3p_edits[0]
            blocks.append(AGPBlock(
                type="W", component_id=e["comp_id"], length=e["length"],
                orientation="+", source_seq=components_pool[e["comp_id"]]
            ))
            e["applied_delta"] = e["length"]
            for extra in telo_3p_edits[1:]:
                extra["applied_delta"] = 0
                extra["note"] = "SKIPPED_DUPLICATE_TELO_3P"

        # Strip terminal N-gap blocks and trim terminal Ns from boundary W blocks
        blocks = strip_terminal_gaps(blocks, components_pool)

        # Emit AGP Lines for Chromosome (Stateless Accumulator)
        chr_curated_segs_raw = []
        chr_curated_repeats_raw = []
        curr_pos = 1
        valid_part = 1
        for b in blocks:
            if b.length <= 0:
                sys.stderr.write(f"WARNING: Dropping zero/negative length block {b.component_id}: {b.length}\n")
                continue
            beg = curr_pos
            end = curr_pos + b.length - 1
            if b.type == "W":
                agp_lines.append(f"{chrom}\t{beg}\t{end}\t{valid_part}\tW\t{b.component_id}\t1\t{b.length}\t{b.orientation}")
                m_inv = re.match(rf"^{re.escape(chrom)}_inv_(\d+)_(\d+)$", b.component_id) if b.component_id else None
                m_slice = re.match(rf"^{re.escape(chrom)}_(\d+)_(\d+)$", b.component_id) if b.component_id else None
                m = m_inv or m_slice
                if m:
                    src_s, src_e = int(m.group(1)), int(m.group(2))
                    # Map GRS/Core segments
                    for seg_s, seg_e, cls in input_segs_by_chr.get(chrom, []):
                        ovl_s = max(seg_s, src_s)
                        ovl_e = min(seg_e, src_e)
                        if ovl_s < ovl_e:
                            if b.orientation == "-":
                                f_s = (beg - 1) + (src_e - ovl_e)
                                f_e = (beg - 1) + (src_e - ovl_s)
                            else:
                                f_s = (beg - 1) + (ovl_s - src_s)
                                f_e = (beg - 1) + (ovl_e - src_s)
                            chr_curated_segs_raw.append((f_s, f_e, cls))
                    # Map RED repeat mask
                    for rep_s, rep_e, lbl in red_repeats_by_chr.get(chrom, []):
                        ovl_s = max(rep_s, src_s)
                        ovl_e = min(rep_e, src_e)
                        if ovl_s < ovl_e:
                            if b.orientation == "-":
                                f_s = (beg - 1) + (src_e - ovl_e)
                                f_e = (beg - 1) + (src_e - ovl_s)
                            else:
                                f_s = (beg - 1) + (ovl_s - src_s)
                                f_e = (beg - 1) + (ovl_e - src_s)
                            chr_curated_repeats_raw.append((f_s, f_e, lbl))
                else:
                    chr_curated_segs_raw.append((beg - 1, end, "putative_core"))
            else:
                gap_type = "repeat" if (b.component_id and "rdna" in str(b.component_id).lower()) else "scaffold"
                gap_evidence = "unspecified" if gap_type == "repeat" else "paired-ends"
                agp_lines.append(f"{chrom}\t{beg}\t{end}\t{valid_part}\tN\t{b.length}\t{gap_type}\tyes\t{gap_evidence}")
                chr_curated_segs_raw.append((beg - 1, end, "putative_core"))

            curr_pos = end + 1
            valid_part += 1

        total_chr_len = curr_pos - 1
        for s, e, cls in merge_curated_pde_segments(chr_curated_segs_raw, total_chr_len):
            curated_segments.append((chrom, s, e, cls))
        for s, e, lbl in merge_curated_repeats(chr_curated_repeats_raw):
            curated_repeats.append((chrom, s, e, lbl))

        # Record Audit Info
        for e in edits:
            if e["type"] == "SUBSUMED_BY_RDNA_ARRAY":
                audit_rows.append({
                    "sample": args.sample, "chrom": chrom, "card_id": e["card_id"],
                    "action": e["action"], "edit_type": e["type"],
                    "orig_start": e.get("orig_start_pos", e["start"]), "orig_end": e.get("orig_end_pos", e["end"]),
                    "delta_bp": 0
                })
                continue
            row = {
                "sample": args.sample, "chrom": chrom, "card_id": e["card_id"],
                "action": e["action"], "edit_type": e["type"],
                "orig_start": e.get("orig_start_pos", e["start"]), "orig_end": e.get("orig_end_pos", e["end"]),
                "delta_bp": e.get("applied_delta", 0)
            }
            if "note" in e:
                row["note"] = e["note"]
            audit_rows.append(row)

    # --- PROCESS UNPLACED SCAFFOLDS ---
    haplotigs = {}
    purged_rdna = {}

    for scaff in unplaced_names:
        scaff_seq = ref_seqs[scaff]
        scaff_len = len(scaff_seq)
        act = scaffold_actions.get(scaff, "")

        if act == "PURGE_RDNA_FRAGMENT" or "RDNA" in act or act == "ROUTE_TO_RDNA_CASSETTE":
            purged_rdna[scaff] = scaff_seq
            audit_rows.append({
                "sample": args.sample, "chrom": scaff, "card_id": f"SCAFF-{scaff}",
                "action": act, "edit_type": "PURGE_RDNA",
                "orig_start": 0, "orig_end": scaff_len, "delta_bp": -scaff_len
            })
        elif (
            act in ("PURGE_HAPLOTIG", "PURGE_SATELLITE", "RETIRE_REDUNDANT_SATELLITE", "Tag as Unplaced Satellite Array")
            or "HAPLOTIG" in act
            or "SATELLITE" in act
        ):
            haplotigs[scaff] = scaff_seq
            audit_rows.append({
                "sample": args.sample, "chrom": scaff, "card_id": f"SCAFF-{scaff}",
                "action": act, "edit_type": "PURGE_HAPLOTIG",
                "orig_start": 0, "orig_end": scaff_len, "delta_bp": -scaff_len
            })
        elif (
            act in ("PURGE_ORGANELLAR", "ISOLATE_MITOCHONDRION", "Isolate into Mitochondrial Release")
            or "MT" in scaff
            or "mitochondrion" in scaff.lower()
        ):
            haplotigs[scaff] = scaff_seq
            audit_rows.append({
                "sample": args.sample, "chrom": scaff, "card_id": f"SCAFF-{scaff}",
                "action": act, "edit_type": "PURGE_ORGANELLAR",
                "orig_start": 0, "orig_end": scaff_len, "delta_bp": -scaff_len
            })
        elif act in ("INTEGRATE_INTO_CHROMOSOME", "INTEGRATE_UNPLACED_REPEAT", "RETIRE_ABSORBED", "RETIRED_ABSORBED"):
            if scaff in materialized_scaffolds:
                # Absorbed into chromosome arm or internal gap; do not emit as standalone unplaced scaffold
                audit_rows.append({
                    "sample": args.sample, "chrom": scaff, "card_id": f"SCAFF-{scaff}",
                    "action": act, "edit_type": "RETIRE_ABSORBED",
                    "orig_start": 0, "orig_end": scaff_len, "delta_bp": -scaff_len
                })
            else:
                # Conservation of Sequence guard: target port was not materialised; retain scaffold!
                components_pool[scaff] = scaff_seq
                agp_lines.append(f"{scaff}\t1\t{scaff_len}\t1\tW\t{scaff}\t1\t{scaff_len}\t+")
                scaff_segs = input_segs_by_chr.get(scaff, [])
                if scaff_segs:
                    for s, e, cls in scaff_segs:
                        curated_segments.append((scaff, s, e, cls))
                else:
                    cls = "putative_GRS" if "grs" in scaff.lower() else "putative_core"
                    curated_segments.append((scaff, 0, scaff_len, cls))
                for s, e, lbl in red_repeats_by_chr.get(scaff, []):
                    curated_repeats.append((scaff, s, e, lbl))
                audit_rows.append({
                    "sample": args.sample, "chrom": scaff, "card_id": f"SCAFF-{scaff}",
                    "action": "RETAIN_UNMATERIALIZED", "edit_type": "RETENTION_GUARD_UNMATERIALIZED",
                    "orig_start": 0, "orig_end": scaff_len, "delta_bp": 0
                })
        else:
            # Retain in primary AGP release
            components_pool[scaff] = scaff_seq
            agp_lines.append(f"{scaff}\t1\t{scaff_len}\t1\tW\t{scaff}\t1\t{scaff_len}\t+")
            scaff_segs = input_segs_by_chr.get(scaff, [])
            if scaff_segs:
                for s, e, cls in scaff_segs:
                    curated_segments.append((scaff, s, e, cls))
            else:
                cls = "putative_GRS" if "grs" in scaff.lower() else "putative_core"
                curated_segments.append((scaff, 0, scaff_len, cls))
            for s, e, lbl in red_repeats_by_chr.get(scaff, []):
                curated_repeats.append((scaff, s, e, lbl))


    # Process Severed Unplaced Scaffolds from Trans-Chromosomal Fissions
    for r in unplaced_severed:
        assigned = r.get("assigned_name", "")
        cid = r.get("contig_id", "")
        clen = safe_int(r.get("contig_length", 0))
        if assigned and clen > 0 and cid:
            if cid not in components_pool:
                components_pool[cid] = "N" * clen
            seq = components_pool[cid]
            trimmed_seq = seq.strip("Nn")
            if len(trimmed_seq) > 0:
                components_pool[cid] = trimmed_seq
                tlen = len(trimmed_seq)
                agp_lines.append(f"{assigned}\t1\t{tlen}\t1\tW\t{cid}\t1\t{tlen}\t+")
                curated_segments.append((assigned, 0, tlen, "putative_GRS" if "grs" in assigned.lower() else "putative_core"))
                for s, e, lbl in red_repeats_by_chr.get(cid, []):
                    if s < tlen:
                        curated_repeats.append((assigned, max(0, s), min(tlen, e), lbl))
                audit_rows.append({
                    "sample": args.sample, "chrom": assigned, "card_id": r.get("card_id", "TRANS"),
                    "action": "SEVER_TO_UNPLACED", "edit_type": "UNPLACED_SCAFFOLD_CREATION",
                    "orig_start": 0, "orig_end": clen, "delta_bp": tlen
                })

    # 4. Write Output AGP
    with open(agp_path, "w") as f:
        f.write("\n".join(agp_lines) + "\n")
    print(f"[{args.sample}] Emitted AGP ({len(agp_lines)} lines): {agp_path}")

    # 5. Write Components FASTA
    with open(comp_path, "w") as f:
        for cid, seq in components_pool.items():
            f.write(f">{cid}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    print(f"[{args.sample}] Emitted components FASTA ({len(components_pool)} records): {comp_path}")

    # 6. Write Segregated FASTAs
    if haplotigs:
        with open(haplo_path, "w") as f:
            for hid, seq in haplotigs.items():
                f.write(f">{hid}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        print(f"[{args.sample}] Emitted haplotigs FASTA ({len(haplotigs)} records): {haplo_path}")
    else:
        open(haplo_path, "w").close()

    if purged_rdna:
        with open(purged_rdna_path, "w") as f:
            for rid, seq in purged_rdna.items():
                f.write(f">{rid}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        print(f"[{args.sample}] Emitted purged rDNA FASTA ({len(purged_rdna)} records): {purged_rdna_path}")
    else:
        # No rDNA fragment was classified for purging. The file is deliberately not created:
        # an empty file is indistinguishable from a purge step that failed, and the absence
        # is already recorded in the audit (and therefore in the scorecard), so nothing is
        # lost by omitting it.
        print(f"[{args.sample}] No rDNA fragments to purge; {os.path.basename(purged_rdna_path)} not emitted.")

    # 7. Write Audit TSV
    audit_header = ["sample", "chrom", "card_id", "action", "edit_type", "orig_start", "orig_end", "delta_bp"]
    with open(audit_path, "w") as f:
        f.write("\t".join(audit_header) + "\n")
        for r in audit_rows:
            f.write("\t".join(str(r.get(h, "")) for h in audit_header) + "\n")
    print(f"[{args.sample}] Emitted curation audit TSV: {audit_path}")

    # 8. Optional: Synthesize Curated FASTA Directly
    if args.synthesize_fasta:
        curated_fasta_path = os.path.join(args.outdir, f"{args.sample}.curated.fasta")
        # Build sequences per object from AGP lines
        agp_by_obj = defaultdict(list)
        for line in agp_lines:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            agp_by_obj[parts[0]].append(parts)

        with open(curated_fasta_path, "w") as f:
            for obj, rows in agp_by_obj.items():
                obj_seq_parts = []
                for r in rows:
                    comp_type = r[4]
                    if comp_type == "W":
                        cid = r[5]
                        ori = r[8]
                        seq = components_pool.get(cid, "")
                        if ori == "-":
                            seq = revcomp(seq)
                        obj_seq_parts.append(seq)
                    elif comp_type == "N":
                        gap_len = safe_int(r[5])
                        obj_seq_parts.append("N" * gap_len)
                full_obj_seq = "".join(obj_seq_parts)
                f.write(f">{obj}\n")
                for i in range(0, len(full_obj_seq), 80):
                    f.write(full_obj_seq[i:i+80] + "\n")
        print(f"[{args.sample}] Synthesized curated FASTA directly: {curated_fasta_path}")

    # 9. Emit Final Curated BED Files (exclusively final assembly coordinates)
    curated_seg_bed = os.path.join(args.outdir, f"{args.sample}.curated.pde_segments.bed")
    with open(curated_seg_bed, "w", encoding="utf-8") as f:
        for c, s, e, cls in curated_segments:
            f.write(f"{c}\t{s}\t{e}\t{cls}\n")
    print(f"[{args.sample}] Emitted curated PDE segmentation BED: {curated_seg_bed}")

    curated_rep_bed = os.path.join(args.outdir, f"{args.sample}.curated.repeats.bed")
    with open(curated_rep_bed, "w", encoding="utf-8") as f:
        for c, s, e, lbl in curated_repeats:
            f.write(f"{c}\t{s}\t{e}\t{lbl}\n")
    print(f"[{args.sample}] Emitted curated repeat mask BED: {curated_rep_bed}")


if __name__ == "__main__":
    main()
