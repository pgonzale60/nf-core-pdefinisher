"""Spanning-read consensus gap bridging (abPOA) for Targeted Gap Assembly.

British spelling standard enforced across all documentation and comments.
"""

from collections import defaultdict
import os
import re
import subprocess

from .flank_io import (
    FLANK_RECORDS,
    bait_reads_for_gap,
    longest_inverted_repeat,
    read_fasta_records,
    read_fastq_records,
    reverse_complement,
    write_flank_reference,
)


def collect_spanning_read_candidates(
    paf_path,
    reads,
    min_mapq=30,
    min_anchor_bp=1000,
    min_secondary_anchor_bp=None,
    max_boundary_slack=500,
):
    """Reads that span the gap across separate LEFT and RIGHT flank records, registered to the gap boundary."""
    if min_secondary_anchor_bp is None:
        min_secondary_anchor_bp = min_anchor_bp

    per_read = defaultdict(lambda: defaultdict(list))
    with open(paf_path) as paf:
        for line in paf:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12 or parts[5] not in FLANK_RECORDS:
                continue
            per_read[parts[0]][parts[5]].append(
                {
                    "qstart": int(parts[2]),
                    "qend": int(parts[3]),
                    "strand": parts[4],
                    "tlen": int(parts[6]),
                    "tstart": int(parts[7]),
                    "tend": int(parts[8]),
                    "aln_len": int(parts[10]),
                    "mapq": int(parts[11]),
                }
            )

    candidates = []
    for qname in sorted(per_read):
        if qname not in reads:
            raise RuntimeError(f"read {qname} aligns to the flanks but is not in the baited reads")
        hits = per_read[qname]
        best = None
        for left in hits.get("LEFT", []):
            left_slack = max(0, left["tlen"] - left["tend"])
            if left_slack > max_boundary_slack:
                continue
            for right in hits.get("RIGHT", []):
                right_slack = max(0, right["tstart"])
                if right_slack > max_boundary_slack:
                    continue
                if left["strand"] != right["strand"]:
                    continue
                if left["mapq"] < min_mapq or right["mapq"] < min_mapq:
                    continue
                left_len = left["aln_len"]
                right_len = right["aln_len"]
                if (
                    max(left_len, right_len) < min_anchor_bp
                    or min(left_len, right_len) < min_secondary_anchor_bp
                ):
                    continue
                if left["strand"] == "+":
                    begin = left["qend"] + left_slack
                    end = right["qstart"] - right_slack
                else:
                    begin = right["qend"] + right_slack
                    end = left["qstart"] - left_slack
                if end < begin:
                    continue
                candidate = {
                    "qname": qname,
                    "strand": left["strand"],
                    "q_from": begin,
                    "q_to": end,
                    "left_anchor": left["aln_len"],
                    "right_anchor": right["aln_len"],
                    "left_slack": left_slack,
                    "right_slack": right_slack,
                }
                if best is None or (candidate["left_anchor"] + candidate["right_anchor"]) > (
                    best["left_anchor"] + best["right_anchor"]
                ):
                    best = candidate
        if best is None:
            continue
        span = reads[qname][0][best["q_from"] : best["q_to"]]
        if best["strand"] == "-":
            span = reverse_complement(span)
        best["seq"] = span
        candidates.append(best)

    return {"candidates": candidates, "n_spanning": len(candidates)}


def write_consensus_input(reads, spanning, out_fq):
    """Write the gap-interval sequence of every spanning read as FASTQ, name-sorted."""
    written = 0
    with open(out_fq, "w") as fh:
        for candidate in sorted(spanning["candidates"], key=lambda c: c["qname"]):
            if not candidate["seq"]:
                continue
            _, qualities = reads[candidate["qname"]]
            span_qual = qualities[candidate["q_from"] : candidate["q_to"]]
            if candidate["strand"] == "-":
                span_qual = span_qual[::-1]
            if len(span_qual) != len(candidate["seq"]):
                raise RuntimeError(f"quality length mismatch for {candidate['qname']}")
            fh.write(f"@{candidate['qname']}\n{candidate['seq']}\n+\n{span_qual}\n")
            written += 1
    return written


def trim_flank_residue(consensus, flank_left, flank_right, max_trim=600, min_identity=0.9):
    """Trim the flank bases that alignment slack leaves at the consensus ends."""

    def identity(a, b):
        return sum(1 for x, y in zip(a, b) if x == y) / float(len(a)) if a else 0.0

    start = 0
    limit = min(len(consensus), len(flank_left), max_trim)
    for length in range(limit, 0, -1):
        if identity(consensus[:length], flank_left[-length:]) >= min_identity:
            start = length
            break

    end = len(consensus)
    limit = min(len(consensus) - start, len(flank_right), max_trim)
    for length in range(limit, 0, -1):
        if identity(consensus[len(consensus) - length :], flank_right[:length]) >= min_identity:
            end = len(consensus) - length
            break

    return start, end


def median_patch_support(carriers, patch, kmer=31):
    """Median fraction of the patch's k-mers present in each read's own candidate sequence."""
    if len(patch) < kmer:
        return 0.0
    patch_kmers = {patch[i : i + kmer] for i in range(len(patch) - kmer + 1)}
    fractions = []
    for carrier in carriers:
        sequence = carrier["seq"]
        kmers = {sequence[i : i + kmer] for i in range(max(0, len(sequence) - kmer + 1))}
        fractions.append(len(patch_kmers & kmers) / float(len(patch_kmers)))
    fractions.sort()
    return fractions[len(fractions) // 2]


def run_abpoa(reads_fq, out_fa, max_num_cons):
    """Consensus of the spanning reads with abPOA, in the order it emits them."""
    result = subprocess.run(
        ["abpoa", "-Q", "-d", str(max_num_cons), "-o", out_fa, reads_fq],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        err = (result.stderr or "").strip()[:300]
        raise RuntimeError(f"abpoa failed (rc={result.returncode}): {err}")
    records = read_fasta_records(out_fa)
    if not records:
        raise RuntimeError(f"abpoa produced no consensus sequence for {reads_fq}")
    return list(records.values())


def claimed_insertion_bp(evidence):
    """The insertion length the card asserts, e.g. {'max_spanning_insertion': '184 bp'}."""
    match = re.search(r"\d+", str((evidence or {}).get("max_spanning_insertion", "")))
    return int(match.group(0)) if match else None


def decide_consensus_patch(
    evidence,
    spanning,
    consensus_seqs,
    flanks,
    min_spanning_reads,
    min_insertion_bp,
    max_excess_ratio,
    min_patch_support,
    gap_len=0,
):
    """Execute the card's decision, or refuse it by name."""
    n_spanning = spanning["n_spanning"]
    if n_spanning < min_spanning_reads:
        return (
            "UNRESOLVED_GAP",
            None,
            f"card premise not reproducible: {n_spanning} reads span the gap, policy floor is {min_spanning_reads}",
        )
    carriers = [c for c in spanning["candidates"] if c["seq"]]
    if len(carriers) < min_spanning_reads:
        return (
            "UNRESOLVED_GAP",
            None,
            f"card claims an inserted sequence but only {len(carriers)} of {n_spanning} spanning reads carry bases for the gap interval",
        )

    pooled = consensus_seqs[0]
    start, end = trim_flank_residue(pooled, flanks["LEFT"], flanks["RIGHT"])
    if end <= start:
        return (
            "UNRESOLVED_GAP",
            None,
            f"the consensus is entirely flank sequence ({len(pooled)} bp), so no missing sequence was recovered",
        )
    patch = pooled[start:end]

    support = median_patch_support(carriers, patch)
    if support < min_patch_support:
        return (
            "UNRESOLVED_GAP",
            None,
            f"spanning reads do not agree on a patch: median k-mer support {support:.2f} is below the policy floor {min_patch_support:.2f}",
        )

    if len(patch) < min_insertion_bp:
        return (
            "UNRESOLVED_GAP",
            None,
            f"card claims an insertion but the trimmed consensus is {len(patch)} bp, below the policy floor of {min_insertion_bp} bp",
        )

    claimed = claimed_insertion_bp(evidence)
    expected_span = (gap_len or 0) + (claimed or 0)
    if expected_span > 0 and len(patch) > expected_span * max_excess_ratio:
        return (
            "UNRESOLVED_GAP",
            None,
            f"card asserts an expected gap-spanning length of {expected_span} bp (gap {gap_len} bp + asserted insertion {claimed} bp) but the consensus is {len(patch)} bp, beyond the policy ratio of {max_excess_ratio}",
        )

    claimed_str = f"{claimed} bp" if claimed else "an unparsed length"
    notes = (
        f"consensus (abpoa) of {n_spanning} spanning reads, {len(carriers)} carrying sequence; "
        f"support {support:.2f}; trimmed {len(pooled) - len(patch)} bp of flank residue; "
        f"{len(consensus_seqs)} consensus sequence(s); card asserted {claimed_str}"
    )
    if len(consensus_seqs) > 1:
        notes += "; a second consensus passed the allele-frequency floor, so the gap may be heterozygous"
    return ("PATCH_GAP_INSERTION", patch, notes)


def resolve_spanned_gap(gap, gap_work, outdir, bam, ref, threads, policy, sample_id, _index, _total):
    """Recover the missing sequence of a spanned gap as the consensus of the reads spanning it."""
    gid = gap["gap_id"]
    chrom = gap["chrom"]
    g_start, g_end, g_len = gap["start"], gap["end"], gap["length"]

    flank_len = policy.get("genomic_architecture", {}).get("gap_assembly_flank_bp", 20000)
    min_mapq = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30)
    min_anchor_bp = policy.get("mapping_and_anchors", {}).get("anchor_min_length_bp", 1000)
    min_spanning_reads = policy.get("evidence_floors", {}).get("min_read_support_floor", 3)
    min_baited_reads = min_spanning_reads
    patch_policy = policy.get("gap_patch", {})
    min_insertion_bp = policy.get("evidence_floors", {}).get("min_insertion_bp", 50)
    max_num_cons = patch_policy.get("consensus_max_num", 2)
    max_excess_ratio = patch_policy.get("premise_max_excess_ratio", 2.0)
    min_patch_support = patch_policy.get("min_patch_support_fraction", 0.9)
    max_inverted_repeat_bp = patch_policy.get("max_anchor_inverted_repeat_bp", 500)

    def row(outcome, left_anchor, right_anchor, patch_len, excised, source, patch_file, notes):
        return [
            gid,
            chrom,
            g_start,
            g_end,
            g_len,
            outcome,
            left_anchor,
            right_anchor,
            patch_len,
            excised,
            source,
            patch_file,
            notes,
        ]

    reads_fq = os.path.join(gap_work, "spanning_reads.fq")
    n_baited = bait_reads_for_gap(bam, chrom, g_start, g_end, flank_len, reads_fq, fmt="fastq")
    if n_baited < min_baited_reads:
        return row(
            "UNRESOLVED_GAP",
            0,
            0,
            0,
            0,
            "NONE",
            "NONE",
            f"Insufficient baited reads ({n_baited} reads, required >= {min_baited_reads})",
        )

    flank_ref_fa = os.path.join(gap_work, "flank_ref.fa")
    anchor_seqs = write_flank_reference(
        ref,
        chrom,
        g_start,
        g_end,
        flank_len,
        flank_ref_fa,
        min_anchor_bp=min_anchor_bp,
        max_inverted_repeat_bp=max_inverted_repeat_bp,
    )
    print(
        f"[{sample_id}] {gid}: anchor windows LEFT={len(anchor_seqs['LEFT'])} bp (inverted repeat {longest_inverted_repeat(anchor_seqs['LEFT'])}), "
        f"RIGHT={len(anchor_seqs['RIGHT'])} bp (inverted repeat {longest_inverted_repeat(anchor_seqs['RIGHT'])})",
        flush=True,
    )

    read_paf = os.path.join(gap_work, "reads_to_flanks.paf")
    with open(read_paf, "w") as paf_fh:
        subprocess.run(
            ["minimap2", "-x", "map-hifi", "--secondary=no", "-t", str(threads), flank_ref_fa, reads_fq],
            stdout=paf_fh,
            stderr=subprocess.DEVNULL,
        )

    reads = read_fastq_records(reads_fq)
    min_secondary_anchor_bp = patch_policy.get("min_secondary_anchor_bp", 300)
    spanning = collect_spanning_read_candidates(
        read_paf,
        reads,
        min_mapq=min_mapq,
        min_anchor_bp=min_anchor_bp,
        min_secondary_anchor_bp=min_secondary_anchor_bp,
    )

    consensus_input_fq = os.path.join(gap_work, "consensus_input.fq")
    n_written = write_consensus_input(reads, spanning, consensus_input_fq)
    if n_written < min_spanning_reads:
        return row(
            "UNRESOLVED_GAP",
            0,
            0,
            0,
            0,
            "NONE",
            "NONE",
            f"card premise not reproducible: {spanning['n_spanning']} reads span the gap but only {n_written} carry sequence for it (policy floor {min_spanning_reads})",
        )

    abpoa_fa = os.path.join(gap_work, "abpoa_consensus.fa")
    consensus_seqs = run_abpoa(consensus_input_fq, abpoa_fa, max_num_cons)

    flanks = read_fasta_records(flank_ref_fa)
    outcome, patch, notes = decide_consensus_patch(
        gap.get("evidence", {}),
        spanning,
        consensus_seqs,
        flanks,
        min_spanning_reads,
        min_insertion_bp,
        max_excess_ratio,
        min_patch_support,
        gap_len=g_len,
    )
    if outcome != "PATCH_GAP_INSERTION":
        return row(outcome, 0, 0, 0, 0, "NONE", "NONE", notes)

    carriers = [c for c in spanning["candidates"] if c["seq"]]
    left_anchor = max(c["left_anchor"] for c in carriers)
    right_anchor = max(c["right_anchor"] for c in carriers)
    patch_file = os.path.join(outdir, f"{gid}_patch.fa")
    with open(patch_file, "w") as fh:
        fh.write(f">{gid}\n{patch}\n")
    source = f"CONSENSUS(abpoa;n={len(carriers)};len={len(patch)})"
    print(f"[{sample_id}] {gid}: {len(patch)} bp consensus from {len(carriers)} spanning reads", flush=True)
    return row("PATCH_GAP_INSERTION", left_anchor, right_anchor, len(patch), g_len, source, patch_file, notes)
