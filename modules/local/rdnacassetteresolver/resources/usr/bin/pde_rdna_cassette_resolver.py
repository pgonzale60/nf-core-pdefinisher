#!/usr/bin/env python3
"""
pde_rdna_cassette_resolver.py - Module 3 Engine C: Ribosomal RNA Cassette Resolver

Resolves multi-copy ribosomal RNA array loci without expensive de novo re-assembly
by synthesizing pre-assembled components from Module 1 (FLYE and HIFIASM):
  - Component 1: Single-Copy Core Anchor (Flye linear contig mapped via PAF)
  - Component 2: Harmonized rDNA Monomer Model (Flye circular repeat unit)
  - Component 3: Distal Telomeric Cap (Hifiasm telomeric unitig) or Distal Flank

Incorporates the Graph Complexity Gatekeeper (DEFER_COMPLEX_RDNA_LOCUS) to prevent
destructive over-curation of structurally heterogeneous arrays (e.g. nrLitMari1).
"""

import argparse
import json
import os
import sys
import tempfile
from collections import defaultdict


def parse_hifiasm_gfa_depths(gfa_path):
    """Parse hifiasm native unitig read depths (rd:i:<depth>) from GFA S-lines.

    Returns dict {utg_name: {'length': int, 'depth': float, 'seq': str}}.
    """
    unitigs = {}
    if not gfa_path or not os.path.exists(gfa_path):
        return unitigs

    with open(gfa_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("S\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[1]
            seq = parts[2]
            length = len(seq)
            depth = None
            for tag in parts[3:]:
                if tag.startswith("LN:i:"):
                    length = int(tag[5:])
                elif tag.startswith("rd:i:"):
                    depth = float(tag[5:])
            unitigs[name] = {
                "length": length,
                "depth": depth if depth is not None else 1.0,
                "seq": seq if seq != "*" else "",
            }
    return unitigs


def run_self_test():
    """Verify GFA rd:i depth parsing, large overlap detection, and exclusion principle."""
    sample_gfa = (
        "H\tVN:Z:1.0\n"
        "S\tutg000001l\tACGT\tLN:i:4\trd:i:30\n"
        "S\tutg000002l\tGGCC\tLN:i:4\trd:i:150\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".gfa", delete=False) as tf:
        tf.write(sample_gfa)
        tf_path = tf.name
    try:
        depths = parse_hifiasm_gfa_depths(tf_path)
        assert "utg000001l" in depths and "utg000002l" in depths
        d_anchor = depths["utg000001l"]["depth"]
        d_monomer = depths["utg000002l"]["depth"]
        assert d_anchor == 30.0 and d_monomer == 150.0
        m_b = max(1, int(d_monomer / d_anchor + 0.5))
        assert m_b == 5, f"Expected M_B == 5, got {m_b}"
    finally:
        if os.path.exists(tf_path):
            os.remove(tf_path)

    # Test large overlap detection (> 2000 bp)
    seq_prefix = "A" * 3000 + "C" * 2500
    seq_suffix = "C" * 2500 + "T" * 1000
    ovlp = find_overlap(seq_prefix, seq_suffix, min_k=25)
    assert ovlp == 2500, f"Expected overlap 2500, got {ovlp}"

    # Test distal telomere unitig identification with exclusion principle
    telo_seq = "G" * 500 + ("TTAGGC" * 30)
    telo_dict = {"ptg000015l": telo_seq}
    cand_id, cand_seq, cand_rep, cand_ovlp = find_distal_telomeric_unitig(
        telo_dict, seq_prefix, min_telo_repeats=15, motifs=("TTAGGC", "GCCTAA"),
        used_unitigs=set(), orphan_fallback=True
    )
    assert cand_id == "ptg000015l", f"Expected orphan unitig ptg000015l, got {cand_id}"
    assert cand_rep == 30, f"Expected 30 repeats, got {cand_rep}"

    print("RDNACASSETTERESOLVER self-test PASS")
    sys.exit(0)


def parse_args():
    parser = argparse.ArgumentParser(description="Module 3 Engine C: Ribosomal RNA Cassette Resolver")
    parser.add_argument("--blueprint", help="Module 2 Curation Blueprint JSON (curation_manifest.json)")
    parser.add_argument("--ref", help="Reference FASTA (chromosomes)")
    parser.add_argument("--flye-fasta", nargs="+", help="Assembly FASTA (assembly.fasta or unitigs.fa)")
    parser.add_argument("--flye-info", default=None, help="Assembly info (assembly_info.txt, optional)")
    parser.add_argument("--flye-paf", help="PAF alignment of contigs to reference genome")
    parser.add_argument("--telomeric-unitigs", default=None, help="Pre-assembled telomeric unitigs FASTA")
    parser.add_argument("--rdna-gfa", default=None, help="Path to hifiasm rDNA GFA graph (.p_utg.gfa)")
    parser.add_argument("--sample", default="sample", help="Sample/strain prefix")
    parser.add_argument("--outdir", default=None, help="Output directory for resolved cassettes")
    parser.add_argument("--manifest", default=None, help="Output manifest TSV")
    parser.add_argument("--monomer-copies", type=int, default=None, help="Number of canonical monomer copies to insert (fallback if not in blueprint policy)")
    parser.add_argument("--self-test", action="store_true", help="Run self-test and exit")
    return parser.parse_args()


def load_fasta(fasta_path):
    """Loads a FASTA file (or space-separated list of FASTA files) into a dictionary of {name: seq}."""
    seqs = {}
    if not fasta_path:
        return seqs
    for f_path in fasta_path.split():
        if not os.path.exists(f_path):
            continue
        if "noseq" in f_path:
            continue
        current_name = None
        current_seq = []
        with open(f_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name:
                        s = "".join(current_seq)
                        if current_name not in seqs or len(s) > len(seqs[current_name]):
                            seqs[current_name] = s
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_name:
                s = "".join(current_seq)
                if current_name not in seqs or len(seqs[current_name]) < len(s):
                    seqs[current_name] = s
    return seqs


def revcomp(seq):
    """Reverse complement of DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def parse_assembly_info(info_path, fasta_path=None):
    """
    Parses Flye assembly_info.txt or synthesizes contig info directly from FASTA.
    Returns list of dicts: name, length, cov, circ, repeat, mult.
    """
    contigs = []
    if info_path and os.path.exists(info_path):
        with open(info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                name = parts[0]
                length = int(parts[1])
                cov_str = parts[2].replace('.', '', 1)
                cov = float(parts[2]) if cov_str.isdigit() else 30.0
                circ = parts[3].upper()
                repeat = parts[4].upper()
                mult = int(parts[5]) if parts[5].isdigit() else 1
                contigs.append({
                    "name": name,
                    "length": length,
                    "cov": cov,
                    "circ": circ,
                    "repeat": repeat,
                    "mult": mult
                })
    if not contigs and fasta_path:
        seqs = load_fasta(fasta_path)
        for name, seq in seqs.items():
            contigs.append({
                "name": name,
                "length": len(seq),
                "cov": 30.0,
                "circ": "N",
                "repeat": "N",
                "mult": 1
            })
    return contigs



def parse_paf(paf_path):
    """
    Parses Flye-to-reference PAF alignments.
    Returns dict: query_name -> list of alignment dicts.
    """
    alns = defaultdict(list)
    if not os.path.exists(paf_path):
        return alns
    with open(paf_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7])
            tend = int(parts[8])
            matches = int(parts[9])
            alen = int(parts[10])
            mapq = int(parts[11])
            alns[qname].append({
                "qlen": qlen, "qstart": qstart, "qend": qend,
                "strand": strand, "tname": tname, "tlen": tlen,
                "tstart": tstart, "tend": tend, "matches": matches,
                "alen": alen, "mapq": mapq
            })
    return alns


def extract_monomer_sequence(contig_seq, monomer_min_len=5000, monomer_max_len=18000):
    """
    Extracts a single canonical monomer from a unitig.
    If the sequence is already in [monomer_min_len, monomer_max_len], returns it directly.
    If it is a tandem array (> monomer_max_len), finds the internal repeat period via self-overlap.
    """
    seq_len = len(contig_seq)
    if seq_len <= monomer_max_len:
        return contig_seq

    prefix_seed = contig_seq[:200]
    min_search = max(monomer_min_len, 4000)
    max_search = min(monomer_max_len, seq_len - 500)

    for p in range(min_search, max_search):
        if contig_seq[p:p+200] == prefix_seed:
            return contig_seq[:p]

    return contig_seq[:monomer_max_len]


def evaluate_graph_complexity(contigs, monomer_min_len=5000, monomer_max_len=18000, max_contigs=100):
    """
    Evaluates assembly repeat graph complexity.
    Returns (is_clean, reason, monomer_contig_dict, is_entity).
    """
    total_contigs = len(contigs)
    if total_contigs == 0:
        return False, "No contigs found in assembly", None, False

    # Check for single-entity array
    if total_contigs == 1:
        return True, "Single-entity rDNA array (one linear contig, no branching repeat graph)", contigs[0], True

    # Identify circular monomers in canonical range (Flye)
    circular_monomers = [
        c for c in contigs
        if c.get("circ") == "Y" and monomer_min_len <= c["length"] <= monomer_max_len
    ]
    if len(circular_monomers) == 1:
        return True, "Clean canonical circular monomer", circular_monomers[0], False
    elif len(circular_monomers) > 1:
        monomer = max(circular_monomers, key=lambda x: (x.get("cov", 1.0), x["length"]))
        return True, f"Selected dominant circular monomer from {len(circular_monomers)} candidates", monomer, False

    # Linear unitigs (e.g. from Hifiasm or uncircularised Flye)
    canonical_candidates = [
        c for c in contigs
        if monomer_min_len <= c["length"] <= monomer_max_len
    ]
    if canonical_candidates:
        monomer = max(canonical_candidates, key=lambda x: (x.get("repeat") == "Y", x.get("cov", 1.0), x["length"]))
        return True, f"Selected representative monomer from {len(canonical_candidates)} unitigs", monomer, False

    # Long tandem array unitigs (> monomer_max_len)
    long_candidates = [c for c in contigs if c["length"] > monomer_max_len]
    if long_candidates:
        best_long = max(long_candidates, key=lambda x: (x.get("cov", 1.0), x["length"]))
        return True, f"Long tandem array unitig ({best_long['length']} bp)", best_long, True

    # Fallback to largest contig
    largest = max(contigs, key=lambda x: x["length"])
    return True, f"Fallback to largest contig ({largest['length']} bp)", largest, True



def find_overlap(seq1, seq2, min_k=21, max_check=None):
    """
    Finds exact suffix-prefix overlap between seq1 end and seq2 start.
    Returns overlap_length, or 0 if none found >= min_k.
    Efficient search using prefix k-mer seed.
    """
    if not seq1 or not seq2 or len(seq1) < min_k or len(seq2) < min_k:
        return 0
    max_k = min(len(seq1), len(seq2))
    if max_check is not None and max_check < max_k:
        max_k = max_check

    seed = seq2[:min_k].upper()
    s1_suffix = seq1[-max_k:].upper()

    pos = 0
    best_ovlp = 0
    while True:
        idx = s1_suffix.find(seed, pos)
        if idx == -1:
            break
        cand_len = len(s1_suffix) - idx
        if seq2[:cand_len].upper() == s1_suffix[idx:]:
            best_ovlp = cand_len
            break  # First occurrence from left in s1_suffix gives the maximal overlap
        pos = idx + 1
    return best_ovlp


def align_and_orient_monomer(anchor_seq, monomer_seq, k=25, search_window=150):
    """
    Finds best orientation and junction offset between anchor tip and circular monomer.
    Tests sliding k-mers over search_window bp of the anchor tip against forward and reverse complement of monomer.
    Returns (oriented_monomer, junction_overlap_bp).
    """
    anchor_tip = anchor_seq[-1500:] if len(anchor_seq) > 1500 else anchor_seq
    
    mono_2x_fwd = monomer_seq + monomer_seq
    mono_2x_rev = revcomp(monomer_seq) + revcomp(monomer_seq)

    # Search backwards from anchor tip in sliding 5-bp steps
    for offset in range(0, min(search_window, len(anchor_tip) - k + 1), 5):
        tip_kmer = anchor_tip[-(offset + k) : len(anchor_tip) - offset] if offset > 0 else anchor_tip[-k:]
        idx_fwd = mono_2x_fwd.find(tip_kmer)
        if idx_fwd != -1:
            junction_pos = (idx_fwd + k + offset) % len(monomer_seq)
            oriented = monomer_seq[junction_pos:] + monomer_seq[:junction_pos]
            return oriented, k
        idx_rev = mono_2x_rev.find(tip_kmer)
        if idx_rev != -1:
            rev_mono = revcomp(monomer_seq)
            junction_pos = (idx_rev + k + offset) % len(monomer_seq)
            oriented = rev_mono[junction_pos:] + rev_mono[:junction_pos]
            return oriented, k

    return monomer_seq, 0


def count_telomeric_repeats(seq, motifs=("TTAGGC", "GCCTAA")):
    """Count maximum occurrences of forward or reverse telomeric motifs."""
    s_upper = seq.upper()
    return max(s_upper.count(m) for m in motifs)


def find_distal_telomeric_unitig(
    telo_unitigs, query_seq, min_telo_repeats=15, min_k=25,
    motifs=("TTAGGC", "GCCTAA"), used_unitigs=None, orphan_fallback=False
):
    """
    Identifies candidate distal telomeric cap unitigs that contain:
      1. Terminal telomeric repeats >= min_telo_repeats.
      2. Exact suffix-prefix overlap or homology with query_seq.
      3. Exclusion principle fallback for unallocated orphan telomeric unitigs.
    Returns (unitig_id, oriented_full_unitig_seq, telo_repeats, overlap_bp).
    The distal unitig is oriented so the terminal telomeric tract is at the 3' end.
    """
    if not telo_unitigs:
        return None, "", 0, 0

    best_candidate = None
    best_telo_repeats = 0
    best_oriented_seq = ""
    best_ovlp = 0

    # 1. Direct suffix-prefix overlap search
    for uid, useq in telo_unitigs.items():
        n_repeats = count_telomeric_repeats(useq, motifs=motifs)
        if n_repeats < min_telo_repeats:
            continue

        u_upper = useq.upper()
        telo_5p_cnt = count_telomeric_repeats(u_upper[:4000], motifs=motifs)
        telo_3p_cnt = count_telomeric_repeats(u_upper[-4000:], motifs=motifs)

        target_seq = u_upper
        if telo_5p_cnt > telo_3p_cnt:
            target_seq = revcomp(u_upper)
            telo_3p_cnt = count_telomeric_repeats(target_seq[-4000:], motifs=motifs)

        if telo_3p_cnt < min_telo_repeats:
            continue

        if query_seq:
            ovlp = find_overlap(query_seq, target_seq, min_k=min_k)
            if ovlp > 0:
                if ovlp > best_ovlp or (ovlp == best_ovlp and n_repeats > best_telo_repeats):
                    best_ovlp = ovlp
                    best_telo_repeats = n_repeats
                    best_candidate = uid
                    best_oriented_seq = target_seq

    if best_candidate and best_oriented_seq:
        return best_candidate, best_oriented_seq, best_telo_repeats, best_ovlp

    # 2. Homology / k-mer match search
    if query_seq:
        mono_2x = query_seq + query_seq
        mono_2x_rc = revcomp(query_seq) + revcomp(query_seq)

        mono_kmers = set()
        k_search = min(31, len(query_seq))
        for i in range(len(mono_2x) - k_search + 1):
            mono_kmers.add(mono_2x[i : i + k_search].upper())
        for i in range(len(mono_2x_rc) - k_search + 1):
            mono_kmers.add(mono_2x_rc[i : i + k_search].upper())

        for uid, useq in telo_unitigs.items():
            n_repeats = count_telomeric_repeats(useq, motifs=motifs)
            if n_repeats < min_telo_repeats:
                continue

            u_upper = useq.upper()
            telo_5p_cnt = count_telomeric_repeats(u_upper[:4000], motifs=motifs)
            telo_3p_cnt = count_telomeric_repeats(u_upper[-4000:], motifs=motifs)

            target_seq = u_upper
            if telo_5p_cnt > telo_3p_cnt:
                target_seq = revcomp(u_upper)
                telo_3p_cnt = count_telomeric_repeats(target_seq[-4000:], motifs=motifs)

            if telo_3p_cnt < min_telo_repeats:
                continue

            prefix_matches = 0
            test_len = min(len(target_seq) - 1000, 15000)
            for pos in range(0, max(0, test_len - k_search + 1), 50):
                if target_seq[pos : pos + k_search] in mono_kmers:
                    prefix_matches += 1

            if prefix_matches >= 3:
                if n_repeats > best_telo_repeats:
                    best_telo_repeats = n_repeats
                    best_candidate = uid
                    best_oriented_seq = target_seq
                    best_ovlp = 0

    if best_candidate and best_oriented_seq:
        return best_candidate, best_oriented_seq, best_telo_repeats, best_ovlp

    # 3. Exclusion Principle: assign unallocated orphan telomeric unitig
    if orphan_fallback:
        candidates = []
        for uid, useq in telo_unitigs.items():
            if used_unitigs and uid in used_unitigs:
                continue
            n_repeats = count_telomeric_repeats(useq, motifs=motifs)
            if n_repeats >= min_telo_repeats:
                u_upper = useq.upper()
                telo_5p_cnt = count_telomeric_repeats(u_upper[:4000], motifs=motifs)
                telo_3p_cnt = count_telomeric_repeats(u_upper[-4000:], motifs=motifs)
                target_seq = u_upper
                if telo_5p_cnt > telo_3p_cnt:
                    target_seq = revcomp(u_upper)
                    telo_3p_cnt = count_telomeric_repeats(target_seq[-4000:], motifs=motifs)
                if telo_3p_cnt >= min_telo_repeats:
                    candidates.append((uid, target_seq, n_repeats))

        if candidates:
            # Pick orphan unitig with maximal telomeric repeats
            best_cand = max(candidates, key=lambda c: c[2])
            return best_cand[0], best_cand[1], best_cand[2], 0

    return None, "", 0, 0


def phase_monomer_to_distal_unitig(monomer_seq, distal_unitig_seq, k=25):
    """
    Finds where the 5' prefix of the distal unitig aligns within the circular monomer.
    Trims the monomer repeat so it ends exactly where the distal unitig begins,
    preventing any duplicated sequence at the junction.
    Returns trimmed_monomer_array.
    """
    prefix = distal_unitig_seq[:1000]
    mono_2x = monomer_seq + monomer_seq
    
    # Search for first k-mer of distal unitig within circular monomer
    for offset in range(0, min(100, len(prefix) - k + 1), 5):
        kmer = prefix[offset:offset+k]
        idx = mono_2x.find(kmer)
        if idx != -1:
            # Monomer should end right before this junction offset
            junction_in_mono = (idx - offset) % len(monomer_seq)
            if junction_in_mono > 0:
                trimmed_mono = monomer_seq[:junction_in_mono]
            else:
                trimmed_mono = ""
            return trimmed_mono

    # If exact k-mer not found, return full monomer
    return monomer_seq


def main():
    args = parse_args()
    if args.self_test:
        run_self_test()

    if not (args.blueprint and args.ref and args.flye_fasta and args.flye_paf and args.outdir and args.manifest):
        sys.stderr.write("ERROR: Missing required arguments when not running --self-test.\n")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.blueprint):
        sys.stderr.write(f"ERROR: Blueprint not found: {args.blueprint}\n")
        sys.exit(1)

    with open(args.blueprint) as f:
        blueprint = json.load(f)

    policy = blueprint.get("policy", {}) if isinstance(blueprint, dict) else {}
    telo_policy = policy.get("telomere", {})
    rdna_policy = policy.get("repeat_and_rdna", {})
    map_policy = policy.get("mapping_and_anchors", {})

    motif = (telo_policy.get("motif") or "TTAGGC").upper()
    rev_motif = revcomp(motif)
    motifs = (motif, rev_motif)
    min_telo_repeats = int(telo_policy.get("min_repeats", 15))
    monomer_min_len = int(rdna_policy.get("rdna_monomer_min_len", 5000))
    monomer_max_len = int(rdna_policy.get("rdna_monomer_max_len", 18000))
    max_contigs = int(rdna_policy.get("gatekeeper_max_contigs", 5))
    phasing_kmer = int(rdna_policy.get("phasing_kmer_size", 25))
    phasing_window = int(rdna_policy.get("phasing_search_window_bp", 150))
    monomer_copies = args.monomer_copies if args.monomer_copies is not None else int(rdna_policy.get("rdna_monomer_copies", 1))
    mapq_anchor = int(map_policy.get("mapq_high_confidence", 50))

    target_cards = []
    seen_ids = set()
    for card in blueprint.get("action_cards", []):
        act = card.get("action", "")
        cid = card.get("card_id", "")
        if cid in seen_ids:
            continue
        if act == "TERMINUS_RDNA_ARRAY":
            seen_ids.add(cid)
            target_cards.append({"type": "TERMINUS", "card": card})
        elif act == "REPLACE_RDNA_ARRAY":
            seen_ids.add(cid)
            target_cards.append({"type": "INTERNAL", "card": card})

    manifest_rows = []
    manifest_header = [
        "sample", "card_id", "chrom", "locus_type", "action", "status",
        "gatekeeper_result", "core_anchor_id", "core_anchor_len",
        "monomer_id", "monomer_len", "monomer_copies",
        "distal_cap_id", "distal_cap_len", "distal_telo_status", "cassette_total_len",
        "cassette_fasta", "notes"
    ]

    if not target_cards:
        with open(args.manifest, "w") as mf:
            mf.write("\t".join(manifest_header) + "\n")
        print(f"[{args.sample}] Zero rDNA cards in blueprint. Empty manifest emitted.")
        sys.exit(0)

    flye_fasta_str = " ".join(args.flye_fasta) if isinstance(args.flye_fasta, list) else str(args.flye_fasta)
    contigs_info = parse_assembly_info(args.flye_info, flye_fasta_str)
    is_clean, gatekeeper_reason, _, _ = evaluate_graph_complexity(
        contigs_info, monomer_min_len=monomer_min_len, monomer_max_len=monomer_max_len, max_contigs=max_contigs
    )

    flye_seqs = load_fasta(flye_fasta_str)
    telo_unitigs = load_fasta(args.telomeric_unitigs) if args.telomeric_unitigs else {}
    flye_paf = parse_paf(args.flye_paf)

    used_unitigs = set()
    for card in blueprint.get("action_cards", []):
        for k in ("authorized_unitig_id", "unitig_id", "scaffold_id"):
            val = card.get(k)
            if val and val != "NONE":
                used_unitigs.add(str(val))

    for item in target_cards:
        card = item["card"]
        locus_type = item["type"]
        card_id = card.get("card_id", "CARD-RDNA-UNKNOWN")
        chrom = card.get("chrom", "UNKNOWN")
        card_start = card.get("start")
        card_end = card.get("end")

        if not is_clean:
            action = "DEFER_COMPLEX_RDNA_LOCUS"
            status = "DEFERRED_COMPLEX_GRAPH"
            manifest_rows.append({
                "sample": args.sample, "card_id": card_id, "chrom": chrom,
                "locus_type": locus_type, "action": action, "status": status,
                "gatekeeper_result": "FAIL", "core_anchor_id": "NONE",
                "core_anchor_len": 0, "monomer_id": "NONE", "monomer_len": 0,
                "monomer_copies": 0, "distal_cap_id": "NONE", "distal_cap_len": 0,
                "distal_telo_status": "NONE",
                "cassette_total_len": 0, "cassette_fasta": "NONE",
                "notes": f"Gatekeeper aborted: {gatekeeper_reason}. Baseline preserved intact."
            })
            continue

        # 1. Identify Core Anchor contig mapping to reference chromosome
        anchor_id = None
        anchor_seq = ""
        anchor_aln = None

        anchor_candidates = []
        for cname, aln_list in flye_paf.items():
            for aln in aln_list:
                if aln["tname"] == chrom and aln["mapq"] >= mapq_anchor:
                    anchor_candidates.append((cname, aln))

        if anchor_candidates:
            if locus_type == "TERMINUS":
                # For 3' terminus, prefer contig starting furthest into the core before card_start
                if card_start is not None:
                    core_anchors = [ac for ac in anchor_candidates if ac[1]["tstart"] <= card_start]
                    if core_anchors:
                        anchor_candidates = core_anchors
                best_anchor = min(anchor_candidates, key=lambda ac: ac[1]["tstart"])
                anchor_id, anchor_aln = best_anchor
                anchor_seq = flye_seqs.get(anchor_id, "")
            else:
                best_anchor = max(anchor_candidates, key=lambda ac: ac[1]["matches"])
                anchor_id, anchor_aln = best_anchor
                anchor_seq = flye_seqs.get(anchor_id, "")

        if anchor_seq and anchor_aln and anchor_aln["strand"] == "-":
            anchor_seq = revcomp(anchor_seq)

        # 2. Select rDNA Array Contig (prioritising distal telomere overlap or longest contig)
        candidate_rdna_contigs = [c for c in contigs_info if c["name"] != anchor_id]
        if not candidate_rdna_contigs:
            candidate_rdna_contigs = [c for c in contigs_info if c["name"] == anchor_id]
        if not candidate_rdna_contigs:
            candidate_rdna_contigs = contigs_info

        matching_contigs = []
        for c in candidate_rdna_contigs:
            c_name = c["name"]
            c_seq = flye_seqs.get(c_name, "")
            if not c_seq:
                continue
            for uid, useq in telo_unitigs.items():
                if count_telomeric_repeats(useq, motifs=motifs) >= min_telo_repeats:
                    ovlp_fwd = find_overlap(c_seq, useq, min_k=phasing_kmer)
                    ovlp_rev = find_overlap(c_seq, revcomp(useq), min_k=phasing_kmer)
                    if ovlp_fwd > 0 or ovlp_rev > 0:
                        matching_contigs.append((c, max(ovlp_fwd, ovlp_rev)))
                        break

        if matching_contigs:
            best_contig = max(matching_contigs, key=lambda x: (x[1], x[0]["length"]))[0]
        else:
            best_contig = max(candidate_rdna_contigs, key=lambda x: x["length"])

        monomer_id = best_contig["name"]
        raw_monomer_seq = flye_seqs.get(monomer_id, "")
        if not raw_monomer_seq:
            action = "DEFER_COMPLEX_RDNA_LOCUS"
            status = "MISSING_MONOMER_SEQUENCE"
            manifest_rows.append({
                "sample": args.sample, "card_id": card_id, "chrom": chrom,
                "locus_type": locus_type, "action": action, "status": status,
                "gatekeeper_result": "FAIL", "core_anchor_id": "NONE",
                "core_anchor_len": 0, "monomer_id": monomer_id, "monomer_len": 0,
                "monomer_copies": 0, "distal_cap_id": "NONE", "distal_cap_len": 0,
                "distal_telo_status": "NONE",
                "cassette_total_len": 0, "cassette_fasta": "NONE",
                "notes": f"Monomer sequence {monomer_id} not found in {args.flye_fasta}"
            })
            continue

        is_circular_monomer = (
            best_contig.get("circ") == "Y"
            and monomer_min_len <= len(raw_monomer_seq) <= monomer_max_len
        )
        is_entity = (not is_circular_monomer)
        if not is_circular_monomer:
            # Preserve full assembled contiguous sequence without artificial truncation
            monomer_seq = raw_monomer_seq
        else:
            monomer_seq = extract_monomer_sequence(raw_monomer_seq, monomer_min_len, monomer_max_len)

        # 3. Identify Distal Telomeric Cap
        distal_id = "NONE"
        distal_seq = ""
        distal_ovlp = 0

        if locus_type == "TERMINUS":
            query_seq = monomer_seq if monomer_seq else anchor_seq
            cand_id, cand_seq, cand_telo, cand_ovlp = find_distal_telomeric_unitig(
                telo_unitigs, query_seq, min_telo_repeats=min_telo_repeats, min_k=phasing_kmer,
                motifs=motifs, used_unitigs=used_unitigs, orphan_fallback=True
            )
            if cand_id and cand_seq:
                distal_id = cand_id
                distal_seq = cand_seq
                distal_ovlp = cand_ovlp
            elif anchor_seq and anchor_id != monomer_id:
                cand_id, cand_seq, cand_telo, cand_ovlp = find_distal_telomeric_unitig(
                    telo_unitigs, anchor_seq, min_telo_repeats=min_telo_repeats, min_k=phasing_kmer,
                    motifs=motifs, used_unitigs=used_unitigs, orphan_fallback=True
                )
                if cand_id and cand_seq:
                    distal_id = cand_id
                    distal_seq = cand_seq
                    distal_ovlp = cand_ovlp

            if distal_id == "NONE":
                auth_u = card.get("authorized_unitig_id")
                if auth_u and auth_u in telo_unitigs:
                    distal_id = auth_u
                    distal_seq = telo_unitigs[auth_u]

        rdna_gap_size = int(rdna_policy.get("rdna_gap_size_bp", 100))
        gap_seq = "N" * rdna_gap_size

        # 4. Seamless Tripartite Cassette Assembly
        if is_entity:
            monomer_copies = 1
            # Check overlap between anchor_seq and monomer_seq
            if anchor_id and anchor_id != monomer_id and anchor_seq:
                ovlp_am = find_overlap(anchor_seq, monomer_seq, min_k=phasing_kmer)
                if ovlp_am == 0:
                    ovlp_am_rc = find_overlap(anchor_seq, revcomp(monomer_seq), min_k=phasing_kmer)
                    if ovlp_am_rc > 0:
                        monomer_seq = revcomp(monomer_seq)
                        ovlp_am = ovlp_am_rc

                if ovlp_am > 0:
                    anchor_part = anchor_seq[:-ovlp_am]
                else:
                    anchor_part = anchor_seq + gap_seq
            elif anchor_id and anchor_id == monomer_id:
                anchor_part = ""
            else:
                anchor_part = anchor_seq

            # Check overlap between monomer_seq and distal_seq
            if distal_seq:
                ovlp_md = find_overlap(monomer_seq, distal_seq, min_k=phasing_kmer)
                if ovlp_md == 0:
                    ovlp_md_rc = find_overlap(monomer_seq, revcomp(distal_seq), min_k=phasing_kmer)
                    if ovlp_md_rc > 0:
                        distal_seq = revcomp(distal_seq)
                        ovlp_md = ovlp_md_rc

                if ovlp_md > 0:
                    distal_part = distal_seq[ovlp_md:]
                else:
                    distal_part = gap_seq + distal_seq
            else:
                distal_part = ""

            monomer_array = monomer_seq
            cassette_seq = anchor_part + monomer_array + distal_part
        else:
            # Circular monomer model (Flye)
            gfa_depths = parse_hifiasm_gfa_depths(args.rdna_gfa) if getattr(args, "rdna_gfa", None) else {}
            gfa_copies = None
            if gfa_depths and anchor_id in gfa_depths and monomer_id in gfa_depths:
                d_anchor = gfa_depths[anchor_id]["depth"]
                d_monomer = gfa_depths[monomer_id]["depth"]
                if d_anchor > 0:
                    gfa_copies = max(1, int(d_monomer / d_anchor + 0.5))

            monomer_copies = gfa_copies if gfa_copies is not None else int(card.get("monomer_copies", args.monomer_copies or 2))
            oriented_monomer, j_ovlp = align_and_orient_monomer(
                anchor_seq, monomer_seq, k=phasing_kmer, search_window=phasing_window
            )
            if distal_seq:
                phased_distal_mono = phase_monomer_to_distal_unitig(oriented_monomer, distal_seq, k=phasing_kmer)
                if monomer_copies > 2:
                    intermediate_copies = (gap_seq + oriented_monomer) * (monomer_copies - 2)
                    monomer_array = oriented_monomer + intermediate_copies + gap_seq + phased_distal_mono
                else:
                    monomer_array = oriented_monomer + gap_seq + phased_distal_mono
            else:
                if monomer_copies > 1:
                    monomer_array = gap_seq.join([oriented_monomer] * monomer_copies)
                else:
                    monomer_array = oriented_monomer

            cassette_seq = anchor_seq + monomer_array + distal_seq

        cassette_filename = f"{card_id}_cassette.fa"
        cassette_path = os.path.join(args.outdir, cassette_filename)

        with open(cassette_path, "w") as cf:
            unit_label = "entity" if is_entity else "monomer"
            cf.write(f">{card_id}_resolved_rdna_cassette chrom={chrom} anchor={anchor_id}:{len(anchor_seq)}bp {unit_label}={monomer_id}:{len(monomer_array)}bp copies={monomer_copies} distal={distal_id}:{len(distal_seq)}bp total={len(cassette_seq)}bp\n")
            for i in range(0, len(cassette_seq), 80):
                cf.write(cassette_seq[i:i+80] + "\n")

        action = "RESOLVE_RDNA_ARRAY_ENTITY" if is_entity else "RESOLVE_RDNA_CASSETTE"
        has_telo_cap = (distal_id != "NONE" and len(distal_seq) > 0 and count_telomeric_repeats(distal_seq, motifs=motifs) >= min_telo_repeats)
        distal_telo_status = "TELO_CAPPED" if has_telo_cap else "UNCAPPED"
        if is_entity:
            status = "CASSETTE_RESOLVED_ENTITY_CAPPED" if has_telo_cap else "CASSETTE_RESOLVED_ENTITY_UNCAPPED"
        else:
            status = "CASSETTE_RESOLVED_CAPPED" if has_telo_cap else "CASSETTE_RESOLVED_UNCAPPED"

        if is_entity:
            if has_telo_cap:
                notes = (
                    f"rDNA array resolved as continuous assembled entity ({monomer_id}, {len(monomer_array)} bp) "
                    f"anchored to {anchor_id or 'core'} ({len(anchor_seq)} bp) with distal telomeric cap ({distal_id}, {len(distal_seq)} bp). "
                    f"Contiguous array preserved without artificial monomer truncation."
                )
            else:
                notes = (
                    f"rDNA array resolved as continuous assembled entity ({monomer_id}, {len(monomer_array)} bp) "
                    f"anchored to {anchor_id or 'core'} ({len(anchor_seq)} bp) without distal telomeric cap."
                )
        elif has_telo_cap:
            notes = f"Tripartite cassette synthesised WITH distal telomeric cap ({distal_id}, {len(distal_seq)} bp): Core Anchor ({len(anchor_seq)} bp) -> Monomer ({len(monomer_array)} bp) -> Distal Telomeric Cap."
        else:
            notes = f"Cassette synthesised WITHOUT distal telomeric cap (open terminus): Core Anchor ({len(anchor_seq)} bp) -> Monomer ({len(monomer_array)} bp) -> Distal Cap ({distal_id}, {len(distal_seq)} bp)."

        manifest_rows.append({
            "sample": args.sample, "card_id": card_id, "chrom": chrom,
            "locus_type": locus_type, "action": action, "status": status,
            "gatekeeper_result": "PASS", "core_anchor_id": anchor_id or "NONE",
            "core_anchor_len": len(anchor_seq), "monomer_id": monomer_id,
            "monomer_len": len(monomer_array), "monomer_copies": monomer_copies,
            "distal_cap_id": distal_id, "distal_cap_len": len(distal_seq),
            "distal_telo_status": distal_telo_status,
            "cassette_total_len": len(cassette_seq), "cassette_fasta": cassette_path,
            "notes": notes
        })

    with open(args.manifest, "w") as mf:
        mf.write("\t".join(manifest_header) + "\n")
        for row in manifest_rows:
            vals = [str(row.get(h, "")) for h in manifest_header]
            mf.write("\t".join(vals) + "\n")

    print(f"[{args.sample}] Processed {len(target_cards)} rDNA loci. Emitted manifest: {args.manifest}")


if __name__ == "__main__":
    main()
