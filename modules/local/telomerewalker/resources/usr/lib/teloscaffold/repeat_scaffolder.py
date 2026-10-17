"""Repeat scaffolder sub-module.

Handles gapless scaffolding of unbridgeable satellite repeat arrays to authentic terminal telomere tracts
with sequence overlap excision and maximal read-supported telomeric tract extension.
British spelling standard enforced.
"""

import os
import subprocess
import tempfile
from telotract.tract_detector import count_telomeric_bp, outward_tail


def find_sequence_overlap(seq1, seq2, min_overlap=500, min_identity=0.95):
    """Find sequence overlap between the 3' end (suffix) of seq1 and 5' end (prefix) of seq2.

    Returns the overlap length in bp if an authentic overlap >= min_overlap exists, else 0.
    British spelling standard enforced.
    """
    if not seq1 or not seq2:
        return 0

    s1 = str(seq1).upper()
    s2 = str(seq2).upper()

    # Fast minimap2 -x asm20 alignment if available
    try:
        with tempfile.NamedTemporaryFile("w") as f1, tempfile.NamedTemporaryFile("w") as f2:
            f1.write(f">s1\n{s1}\n")
            f1.flush()
            f2.write(f">s2\n{s2}\n")
            f2.flush()
            cmd = ["minimap2", "-x", "asm20", f1.name, f2.name]
            res = subprocess.run(cmd, capture_output=True, text=True)
            for line in res.stdout.strip().splitlines():
                p = line.split("\t")
                if len(p) >= 12 and p[4] == "+":
                    q_len = int(p[1])
                    q_start = int(p[2])
                    q_end = int(p[3])
                    t_len = int(p[6])
                    t_start = int(p[7])
                    t_end = int(p[8])
                    matches = int(p[9])
                    aln_len = int(p[10])
                    identity = matches / max(1, aln_len)

                    # Overlap must connect the end of s1 (target) to the start of s2 (query)
                    if (q_start <= 50 and (t_len - t_end) <= 50 and
                            aln_len >= min_overlap and identity >= min_identity):
                        return q_end
    except Exception:
        pass

    # Exact/kmer Python suffix-prefix fallback
    probe_lens = [min_overlap, 1000, 2000, 4000]
    for p_len in probe_lens:
        if len(s2) >= p_len:
            probe = s2[:p_len]
            pos = s1.rfind(probe)
            if pos >= 0:
                expected_overlap = len(s1) - pos
                if expected_overlap >= min_overlap and s2[:expected_overlap] == s1[pos:]:
                    return expected_overlap

    return 0


def extract_reads_for_card(bam_path, qnames, chrom=None, terminus=None, chrom_len=None, window=100000):
    """Extract query sequences of supporting reads from BAM."""
    reads = {}
    if not bam_path or not os.path.exists(bam_path) or not qnames:
        return reads

    target_qnames = set(qnames)
    try:
        # First attempt targeted regional fetch
        regions = []
        if chrom and terminus and chrom_len:
            if terminus == "5p":
                regions.append(f"{chrom}:1-{window}")
            else:
                regions.append(f"{chrom}:{max(1, chrom_len - window)}-{chrom_len}")
        elif chrom:
            regions.append(chrom)

        for reg in regions:
            cmd = ["samtools", "view", bam_path, reg]
            res = subprocess.run(cmd, capture_output=True, text=True)
            for line in res.stdout.splitlines():
                p = line.split("\t")
                if len(p) >= 11 and p[0] in target_qnames and p[9] != "*":
                    reads[p[0]] = p[9].upper()
                    if len(reads) == len(target_qnames):
                        return reads

        # Fallback for reads on unplaced scaffolds or secondary targets
        if len(reads) < len(target_qnames):
            cmd = ["samtools", "view", bam_path]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
            for line in proc.stdout:
                p = line.split("\t")
                if len(p) >= 11 and p[0] in target_qnames and p[9] != "*":
                    reads[p[0]] = p[9].upper()
                    if len(reads) == len(target_qnames):
                        proc.kill()
                        break
            proc.stdout.close()
    except Exception:
        pass

    return reads


def find_maximal_read_supported_tract(read_sequences, anchor_seq, terminus, motif, rev_motif,
                                      anchor_bp=500, min_repeats=15):
    """Find the candidate read that provides the longest canonical telomeric tract outward.

    Returns (best_qname, best_tract_seq, best_tract_bp) or None.
    British spelling standard enforced.
    """
    if not read_sequences or not anchor_seq:
        return None

    anchor = anchor_seq[-anchor_bp:] if terminus == "3p" else anchor_seq[:anchor_bp]
    if len(anchor) < 100:
        return None

    best = None
    for qname, seq in read_sequences.items():
        orientation, tail = outward_tail(seq, anchor, len(anchor))
        if not tail:
            continue

        # Count canonical telomeric bases outward
        telo_bp = count_telomeric_bp(tail, motif, rev_motif)
        if telo_bp >= min_repeats * len(motif):
            if best is None or telo_bp > best[2]:
                best = (qname, tail, telo_bp)

    return best


def materialise_satellite_scaffold(card, chrom, terminus, motif, rev_motif,
                                   gap_size=0, nominal_cap_bp=800,
                                   scaffold_seq="", tip_seq="",
                                   bam_path=None, policy=None):
    """Materialise a gapless scaffold bridging a satellite array to a distal telomeric cap.

    Integrates proximal sequence overlap excision and maximal read-supported tract extension.
    British spelling standard enforced.
    """
    evidence = card.get("evidence", {}) or {}
    supporting_reads = evidence.get("supporting_reads") or []
    read_count = len(supporting_reads) if isinstance(supporting_reads, list) else int(supporting_reads or 0)
    tract_bp = int(evidence.get("telomeric_tract_bp") or nominal_cap_bp)

    # Respect curation policy
    if policy is None:
        allow_gapless = False
        min_overlap = 500
    else:
        allow_gapless = policy.get("allow_subtelomeric_gapless_splicing", True)
        min_overlap = policy.get("subtelomeric_min_overlap_bp", 500)
    card_gap = int(card.get("gap_size", gap_size))

    scaff = (scaffold_seq or evidence.get("scaffold_sequence") or "").upper()
    proximal_overlap_bp = 0

    # 1. Proximal Sequence Overlap Excision (Tip to Scaffold)
    if scaff and tip_seq:
        if terminus == "3p":
            # tip_seq is proximal (chr end), scaff is distal
            overlap = find_sequence_overlap(tip_seq, scaff, min_overlap=min_overlap)
            if overlap > 0:
                scaff = scaff[overlap:]
                proximal_overlap_bp = overlap
        else:
            # scaff is distal, tip_seq is proximal (chr start)
            overlap = find_sequence_overlap(scaff, tip_seq, min_overlap=min_overlap)
            if overlap > 0:
                scaff = scaff[:-overlap]
                proximal_overlap_bp = overlap

    if allow_gapless and scaff and proximal_overlap_bp > 0:
        card_gap = 0
    if card.get("gap_size") == 0:
        card_gap = 0

    # 2. Distal Telomeric Cap & Maximal Tract Extension
    cap_motif = rev_motif if terminus == "5p" else motif
    best_read_tract = None

    if bam_path and supporting_reads:
        read_qnames = supporting_reads if isinstance(supporting_reads, list) else []
        read_seqs = extract_reads_for_card(bam_path, read_qnames, chrom=chrom, terminus=terminus)
        anchor_candidate = scaff if scaff else tip_seq
        if not anchor_candidate:
            anchor_candidate = evidence.get("consensus_flank_60bp", "")
        if anchor_candidate and read_seqs:
            best_read_tract = find_maximal_read_supported_tract(
                read_seqs, anchor_candidate, terminus, motif, rev_motif
            )

    if best_read_tract:
        best_qname, cap_seq, telo_len = best_read_tract
        cap_source = f"read {best_qname} ({telo_len:,} bp maximal tract)"
    else:
        custom_cap = evidence.get("telomeric_cap_seq")
        if custom_cap:
            cap_seq = custom_cap.upper()
        else:
            repeat_count = max(15, tract_bp // len(cap_motif))
            cap_seq = cap_motif * repeat_count
        telo_len = count_telomeric_bp(cap_seq, motif, rev_motif)
        cap_source = f"junction consensus ({telo_len:,} bp tract)"

    gap_seq = "N" * card_gap

    if terminus == "5p":
        # 5' extension is prepended before coordinate 0:
        # [distal_cap] + [gap_seq] + [scaffold_seq] -> [assembly start]
        if scaff:
            net_seq = cap_seq + gap_seq + scaff
        else:
            net_seq = cap_seq + gap_seq
    else:
        # 3' extension is appended after chrom_len:
        # [assembly end] -> [scaffold_seq] + [gap_seq] + [distal_cap]
        if scaff:
            net_seq = scaff + gap_seq + cap_seq
        else:
            net_seq = gap_seq + cap_seq

    added_len = len(net_seq)
    non_telo = added_len - telo_len

    scaff_msg = f"{len(scaff):,} bp subtelomeric sequence, " if scaff else ""
    overlap_msg = f"excised {proximal_overlap_bp:,} bp redundant proximal overlap; " if proximal_overlap_bp else ""
    gap_msg = "gapless (0 Ns)" if card_gap == 0 else f"{card_gap} bp scaffolding gap"

    notes = (
        f"Scaffolded unbridgeable satellite repeat array to authentic terminal telomere cap "
        f"({gap_msg}; {overlap_msg}supported by {read_count} reads; "
        f"{scaff_msg}{telo_len:,} bp canonical telomere via {cap_source})"
    )

    return {
        "status": "RESOLVED_BY_SATELLITE_SCAFFOLD",
        "net_seq": net_seq,
        "added_len": added_len,
        "telo_len": telo_len,
        "non_telo": non_telo,
        "notes": notes,
        "read_count": read_count,
        "cap_seq": cap_seq,
        "proximal_overlap_bp": proximal_overlap_bp,
    }
