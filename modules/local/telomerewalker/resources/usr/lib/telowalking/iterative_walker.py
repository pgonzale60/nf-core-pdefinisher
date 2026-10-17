"""Iterative walker sub-module.

Implements prefix-suffix minimap2 overlap graph traversal outward from terminal seeds into unmapped reads.
Enforces the somatic diminution exclusion zone. British spelling standard enforced.
"""

import os
import shutil
import subprocess
import tempfile
from telotract.tract_detector import revcomp


class UnmappedReadIndex:
    """Disk-backed index over a library's unmapped reads.

    Streams unmapped reads straight to disk to maintain constant memory footprint.
    """

    def __init__(self, bam_path, workdir):
        self.bam_path = bam_path
        self.fa_path = os.path.join(workdir, "unmapped_reads.fa")
        self._offsets = {}
        self._cache = {}
        self._build()
        self._index()

    def _build(self):
        """Stream `samtools view -b BAM '*' | samtools fasta -` to disk."""
        cmd_view = ["samtools", "view", "-b", self.bam_path, "*"]
        cmd_fasta = ["samtools", "fasta", "-"]
        with open(self.fa_path, "wb") as fa_fh:
            p1 = subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen(cmd_fasta, stdin=p1.stdout, stdout=fa_fh, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.wait()
            p1.wait()

    def _index(self):
        """Single streaming pass recording each record's sequence byte range."""
        if not os.path.exists(self.fa_path):
            return
        with open(self.fa_path, "rb") as fh:
            cur_id = None
            seq_start = 0
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if cur_id is not None:
                        self._offsets[cur_id] = (seq_start, pos - seq_start)
                    cur_id = line[1:].split()[0].decode()
                    seq_start = fh.tell()
            if cur_id is not None:
                self._offsets[cur_id] = (seq_start, fh.tell() - seq_start)

    def __len__(self):
        return len(self._offsets)

    def get_seq(self, qname):
        """Return one read's full uppercased sequence, read back from disk."""
        if qname in self._cache:
            return self._cache[qname]
        span = self._offsets.get(qname)
        if span is None:
            return None
        start, length = span
        with open(self.fa_path, "rb") as fh:
            fh.seek(start)
            raw = fh.read(length)
        seq = "".join(raw.decode().split()).upper()
        self._cache[qname] = seq
        return seq


def _kmer_set(seq, k=21):
    """Extract the set of k-mers from a sequence."""
    if not seq or len(seq) < k:
        return set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def evaluate_step_information_gain(
    consensus_seq,
    ext_seq,
    best_cand,
    all_candidates,
    motif,
    rev_motif,
    min_telo_bp=200,
    is_somatic_boundary=False,
    total_added_before=0,
    min_germline_extension_bp=100,
    kmer_size=21,
):
    """Evaluate whether a single walker step improves the information state or must roll back.

    Three possible outcomes after taking 1 step (`ext_seq`):
      1. `TELOMERE_REACHED`: The new tail carries >= `min_telo_bp` of canonical telomeric repeats
         beyond any somatic exclusion zone. Walk succeeds immediately.
      2. `UNIQUE_EXIT`: The step is uniquely anchored (`mapq >= 20`), non-branching across candidate
         reads, and introduces genuinely novel sequence complexity (`novel_kmer_fraction >= 0.20`
         and internal `unique_kmer_ratio >= 0.35`), proving we are not cycling inside a tandem repeat.
      3. `ROLLBACK`: The step has `mapq < 20` (internal repeat paralogy), conflicts with competing
         read branches (polyvalent repeat), or fails to add novel k-mer complexity (stationary
         satellite array). Roll back immediately to Step 0.
    """
    k = len(motif)
    total_after = total_added_before + len(ext_seq)
    combined_tail = (consensus_seq + ext_seq)[-1000:]
    fwd_hits = combined_tail.count(motif)
    rev_hits = combined_tail.count(rev_motif)
    telo_bp = max(fwd_hits, rev_hits) * k

    if telo_bp >= min_telo_bp:
        if not (is_somatic_boundary and total_after < min_germline_extension_bp):
            return {
                "decision": "TELOMERE_REACHED",
                "telo_bp": telo_bp,
                "novel_kmer_fraction": 1.0,
                "unique_kmer_ratio": 1.0,
                "reason": f"Reached canonical telomere ({telo_bp} bp repeats)",
            }

    # If telomere was NOT reached on this step, check whether information state improved:
    if len(ext_seq) < kmer_size:
        return {
            "decision": "ROLLBACK",
            "telo_bp": telo_bp,
            "novel_kmer_fraction": 0.0,
            "unique_kmer_ratio": 0.0,
            "reason": "Extension shorter than k-mer window without reaching telomere",
        }

    # 1. Internal repeat paralogy check (MAPQ < 20 against target seed indicates tandem repeat array)
    if best_cand.get("mapq", 0) < 20:
        return {
            "decision": "ROLLBACK",
            "telo_bp": telo_bp,
            "novel_kmer_fraction": 0.0,
            "unique_kmer_ratio": 0.0,
            "reason": f"Low anchor MAPQ ({best_cand.get('mapq', 0)} < 20) without reaching telomere (repetitive array)",
        }

    # 2. Polyvalent branching check across competing candidate reads
    best_prefix = ext_seq[: min(250, len(ext_seq))]
    best_prefix_kmers = _kmer_set(best_prefix, kmer_size)
    if best_prefix_kmers:
        for other in all_candidates:
            if other["id"] == best_cand["id"] or other.get("mapq", 0) < 20:
                continue
            other_prefix = other["ext"][: min(250, len(other["ext"]))]
            other_kmers = _kmer_set(other_prefix, kmer_size)
            if not other_kmers:
                continue
            union_n = len(best_prefix_kmers | other_kmers)
            jaccard = len(best_prefix_kmers & other_kmers) / union_n if union_n > 0 else 1.0
            if jaccard < 0.50:
                return {
                    "decision": "ROLLBACK",
                    "telo_bp": telo_bp,
                    "novel_kmer_fraction": 0.0,
                    "unique_kmer_ratio": 0.0,
                    "reason": f"Polyvalent read branch detected ({best_cand['id']} vs {other['id']}, Jaccard={jaccard:.2f})",
                }

    # 3. Novel k-mer information gain & internal sequence complexity
    history_kmers = _kmer_set(consensus_seq, kmer_size)
    total_windows = len(ext_seq) - kmer_size + 1
    ext_kmers = set()
    novel_windows = 0
    for i in range(total_windows):
        km = ext_seq[i : i + kmer_size]
        ext_kmers.add(km)
        if km not in history_kmers:
            novel_windows += 1

    novel_fraction = novel_windows / max(1, total_windows)
    unique_ratio = len(ext_kmers) / max(1, total_windows)

    if novel_fraction < 0.20 or unique_ratio < 0.35:
        return {
            "decision": "ROLLBACK",
            "telo_bp": telo_bp,
            "novel_kmer_fraction": round(novel_fraction, 3),
            "unique_kmer_ratio": round(unique_ratio, 3),
            "reason": (
                f"Stationary repeat array: novel_kmer_fraction={novel_fraction:.2f}, "
                f"unique_kmer_ratio={unique_ratio:.2f}; rolling back to Step 0"
            ),
        }

    return {
        "decision": "UNIQUE_EXIT",
        "telo_bp": telo_bp,
        "novel_kmer_fraction": round(novel_fraction, 3),
        "unique_kmer_ratio": round(unique_ratio, 3),
        "reason": f"Information state improved (novel_kmer_fraction={novel_fraction:.2f}, unique_kmer_ratio={unique_ratio:.2f})",
    }


def run_iterative_walker(bam_path, seed_seq, motif, rev_motif, max_iters, kmer_size=21,
                         min_telo_bp=200, is_somatic_boundary=False, min_germline_extension_bp=100,
                         unmapped=None, threads=1, tip_is_germline=False):
    """Iterative consensus extension outward into unmapped/clipped reads using minimap2 prefix-suffix overlaps
    gated by 1-step information-gain and valency conservation.
    """
    if unmapped is None:
        with tempfile.TemporaryDirectory() as tmpdir:
            return _walk_loop(UnmappedReadIndex(bam_path, tmpdir), seed_seq, motif, rev_motif,
                              max_iters, min_telo_bp, is_somatic_boundary, min_germline_extension_bp,
                              threads, tip_is_germline, kmer_size=kmer_size)

    if len(unmapped) == 0:
        notes = "Zero unmapped reads in library"
        if tip_is_germline:
            notes += " (distal GRS arm)"
        elif is_somatic_boundary:
            notes += " to extend past somatic diminution breaksite"
        return {
            "status": "ZERO_EXTENDING_READS",
            "added_bp": 0,
            "telo_bp": 0,
            "extension_seq": "",
            "notes": notes,
            "failure_reason": "Zero unmapped reads available in BAM library"
        }

    return _walk_loop(unmapped, seed_seq, motif, rev_motif, max_iters, min_telo_bp,
                      is_somatic_boundary, min_germline_extension_bp, threads, tip_is_germline,
                      kmer_size=kmer_size)


def _walk_loop(unmapped, seed_seq, motif, rev_motif, max_iters, min_telo_bp,
               is_somatic_boundary, min_germline_extension_bp, threads=1, tip_is_germline=False,
               kmer_size=21):
    consensus = seed_seq
    total_added = 0
    k = len(motif)
    used_read_ids = set()
    rollback_reason = None

    if not shutil.which("minimap2"):
        raise RuntimeError("minimap2 is required for iterative telomere walking but not found in PATH")

    effective_max_iters = max(1, max_iters)

    with tempfile.TemporaryDirectory() as tmpdir:
        candidates_fa = unmapped.fa_path

        for iteration in range(1, effective_max_iters + 1):
            tail = consensus[-1000:]
            fwd_hits = tail.count(motif)
            rev_hits = tail.count(rev_motif)
            if max(fwd_hits, rev_hits) * k >= min_telo_bp:
                if not (is_somatic_boundary and total_added < min_germline_extension_bp):
                    return {
                        "status": "EXTEND_TO_TELOMERE",
                        "added_bp": total_added,
                        "telo_bp": max(fwd_hits, rev_hits) * k,
                        "extension_seq": consensus[len(seed_seq):],
                        "notes": f"Reached canonical telomere ({max(fwd_hits, rev_hits) * k} bp repeats) at iteration {iteration}"
                    }

            valid_candidates = []
            win_len = min(len(consensus), 4000)
            target_seq = consensus[-win_len:]
            target_fa = os.path.join(tmpdir, f"target_iter_{iteration}.fa")
            with open(target_fa, "w") as tf:
                tf.write(f">target\n{target_seq}\n")

            cmd_mm2 = [
                "minimap2", "-x", "map-hifi", "-c", "--secondary=yes",
                "-t", str(threads),
                "-N", "10", target_fa, candidates_fa
            ]
            mm2_res = subprocess.run(cmd_mm2, capture_output=True, text=True)

            if mm2_res.returncode == 0 and mm2_res.stdout.strip():
                for line in mm2_res.stdout.splitlines():
                    parts = line.strip().split("\t")
                    if len(parts) < 12:
                        continue
                    qname = parts[0]
                    if qname in used_read_ids:
                        continue
                    qlen = int(parts[1])
                    qstart = int(parts[2])
                    qend = int(parts[3])
                    strand = parts[4]
                    tlen = int(parts[6])
                    tend = int(parts[8])
                    nmatch = int(parts[9])
                    alen = int(parts[10])
                    mapq = int(parts[11])

                    # Require high-identity HiFi overlap (>=500 bp at >=98% identity, or >=150 bp at >=95% if read is short)
                    min_req_alen = 500 if tlen >= 1000 else 150
                    min_req_id = 0.98 if tlen >= 1000 else 0.95
                    if alen < min_req_alen or (nmatch / max(1, alen)) < min_req_id:
                        continue

                    if tend < (tlen - 50):
                        continue

                    ext_seq = ""
                    if strand == "+":
                        if (qlen - qend) >= 50:
                            read_seq = unmapped.get_seq(qname)
                            ext_seq = read_seq[qend:] if read_seq else ""
                    else:
                        if qstart >= 50:
                            read_seq = unmapped.get_seq(qname)
                            ext_seq = revcomp(read_seq[:qstart]) if read_seq else ""

                    if not ext_seq or len(ext_seq) < 50:
                        continue

                    if is_somatic_boundary and total_added < min_germline_extension_bp:
                        rem = min_germline_extension_bp - total_added
                        check_zone = ext_seq[:rem + 100]
                        zone_telo = max(check_zone.count(motif), check_zone.count(rev_motif)) * k
                        if zone_telo >= min_germline_extension_bp:
                            continue

                    ext_tail = ext_seq[-1000:]
                    ext_telo_bp = max(ext_tail.count(motif), ext_tail.count(rev_motif)) * k
                    valid_candidates.append({
                        "id": qname,
                        "ext": ext_seq,
                        "mapq": mapq,
                        "len": len(ext_seq),
                        "telo_bp": ext_telo_bp,
                    })

            if not valid_candidates:
                if total_added > 0 and rollback_reason is None:
                    rollback_reason = (
                        f"Walk of {total_added} bp ran out of extending reads at iteration {iteration} "
                        "without reaching canonical telomere; rolled back to Step 0"
                    )
                break

            # Prefer candidates that directly reach a telomere, then higher MAPQ, then length
            best = max(
                valid_candidates,
                key=lambda x: (1 if x["telo_bp"] >= min_telo_bp else 0, x["mapq"], x["len"]),
            )

            gate = evaluate_step_information_gain(
                consensus_seq=consensus,
                ext_seq=best["ext"],
                best_cand=best,
                all_candidates=valid_candidates,
                motif=motif,
                rev_motif=rev_motif,
                min_telo_bp=min_telo_bp,
                is_somatic_boundary=is_somatic_boundary,
                total_added_before=total_added,
                min_germline_extension_bp=min_germline_extension_bp,
                kmer_size=kmer_size,
            )

            if gate["decision"] == "TELOMERE_REACHED":
                consensus += best["ext"]
                total_added += best["len"]
                return {
                    "status": "EXTEND_TO_TELOMERE",
                    "added_bp": total_added,
                    "telo_bp": gate["telo_bp"],
                    "extension_seq": consensus[len(seed_seq):],
                    "notes": f"{gate['reason']} at iteration {iteration}"
                }

            if gate["decision"] == "ROLLBACK":
                # Information state did not improve: roll back to Step 0 immediately!
                rollback_reason = f"Step {iteration} rolled back to Step 0: {gate['reason']}"
                consensus = seed_seq
                total_added = 0
                break

            # UNIQUE_EXIT: advance one step and re-evaluate on next iteration
            used_read_ids.add(best["id"])
            consensus += best["ext"]
            total_added += best["len"]

    # If the loop finished without reaching a canonical telomere, roll back any partial uncapped walk to Step 0!
    if total_added > 0 and rollback_reason is None:
        rollback_reason = (
            f"Walk of {total_added} bp reached horizon limit of {effective_max_iters} steps "
            "without reaching canonical telomere; rolled back to Step 0"
        )
        consensus = seed_seq
        total_added = 0

    # Determine structured algorithmic failure status
    if rollback_reason:
        if "Polyvalent" in rollback_reason:
            status = "POLYVALENT_REPEAT_BRANCH"
        elif "Stationary repeat array" in rollback_reason or "Extension shorter than" in rollback_reason:
            status = "STATIONARY_REPEAT_ARRAY"
        elif "Low anchor MAPQ" in rollback_reason:
            status = "LOW_ANCHOR_MAPQ"
        elif "horizon limit" in rollback_reason or "within" in rollback_reason:
            status = "HORIZON_LIMIT_REACHED"
        elif "ran out of extending reads" in rollback_reason:
            status = "ZERO_EXTENDING_READS"
        else:
            status = "UNRESOLVED_WALK"
        diag_reason = rollback_reason
    else:
        status = "ZERO_EXTENDING_READS"
        diag_reason = "Zero unmapped or extending reads overlap tip beyond seed"

    context_parts = []
    if tip_is_germline:
        context_parts.append("Distal GRS")
    elif is_somatic_boundary:
        context_parts.append("Somatic breaksite")
    prefix = f"[{'; '.join(context_parts)}] " if context_parts else ""

    return {
        "status": status,
        "added_bp": 0,
        "telo_bp": 0,
        "extension_seq": "",
        "notes": f"{prefix}{diag_reason}",
        "failure_reason": diag_reason,
        "tip_is_germline": tip_is_germline
    }

