"""Unitig capping sub-module.

Handles deterministic materialisation of authorized unitig caps, net outboard sequence slicing,
and region intervals. British spelling standard enforced.
"""

import gzip
import os
import subprocess
from telotract.tract_detector import revcomp


def safe_int(val, default=0):
    """Safely cast value to integer."""
    try:
        if val is None or val == "" or str(val).strip().upper() == "NA":
            return default
        return int(float(val))
    except (ValueError, TypeError):
        return default


def outboard_extension(unitig_seq, hit, terminus, chrom_len):
    """The part of a unitig that lies beyond the assembly tip, i.e. what actually gets added.

    The alignment coordinates in the telomere profile are extrapolated 1:1 to the tip to find the boundary.
    """
    if not unitig_seq:
        return ""
    q_start = safe_int(hit.get("q_aln_start", 0))
    q_end = safe_int(hit.get("q_aln_end", 0))
    t_start = safe_int(hit.get("t_start", hit.get("t_aln_start", 0)))
    t_end = safe_int(hit.get("t_end", hit.get("t_aln_end", 0)))
    if t_end <= t_start or q_end <= q_start:
        return unitig_seq

    seq = unitig_seq if hit.get("aln_strand", "+") == "+" else revcomp(unitig_seq)
    q_len = len(seq)
    if hit.get("aln_strand", "+") != "+":
        q_start, q_end = q_len - q_end, q_len - q_start

    if terminus == "5p":
        q_at_tip = q_start - t_start
        return seq[:max(0, min(q_len, q_at_tip))]
    q_at_tip = q_end + (chrom_len - t_end)
    return seq[max(0, min(q_len, q_at_tip)):]


def region_extension(unitig_seq, region):
    """The exact unitig interval Module 2 authorised, oriented as it appears in the chromosome."""
    if not unitig_seq or not region:
        return ""
    q_start = safe_int(region.get("q_start", 0))
    q_end = safe_int(region.get("q_end", 0))
    if q_end <= q_start or q_start < 0 or q_end > len(unitig_seq):
        return ""
    sub = unitig_seq[q_start:q_end]
    if str(region.get("strand", "+")) == "-":
        sub = revcomp(sub)
    return sub


def extract_unitig_by_id(unitigs_fasta, unitig_id):
    """Extract full sequence of a named unitig from FASTA, FASTA.gz, or space-separated list of FASTAs."""
    if not unitigs_fasta or not unitig_id or str(unitig_id) in ("None", "NONE", ""):
        return None

    if isinstance(unitigs_fasta, str):
        fasta_paths = [p.strip() for p in unitigs_fasta.split() if p.strip()]
    elif isinstance(unitigs_fasta, (list, tuple)):
        fasta_paths = list(unitigs_fasta)
    else:
        fasta_paths = [unitigs_fasta]

    best_seq = ""
    for f_path in fasta_paths:
        if not os.path.exists(f_path) or "noseq" in f_path:
            continue

        # Try samtools faidx first if index exists
        fai_path = f_path + ".fai"
        if os.path.exists(fai_path):
            cmd = ["samtools", "faidx", f_path, unitig_id]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode == 0 and res.stdout.strip():
                lines = res.stdout.strip().split("\n")
                cand = "".join(lines[1:]).upper()
                if len(cand) > len(best_seq):
                    best_seq = cand

        # Streaming linear fallback
        open_fn = gzip.open if f_path.endswith(".gz") else open
        mode = "rt" if f_path.endswith(".gz") else "r"
        try:
            with open_fn(f_path, mode) as f:
                target = False
                seq_parts = []
                for line in f:
                    if line.startswith(">"):
                        name = line[1:].strip().split()[0]
                        if name == unitig_id or name.startswith(unitig_id):
                            target = True
                        elif target:
                            break
                    elif target:
                        seq_parts.append(line.strip().upper())
                if seq_parts:
                    cand = "".join(seq_parts)
                    if len(cand) > len(best_seq):
                        best_seq = cand
        except Exception:
            continue

    return best_seq if len(best_seq) > 1 else None


def authorized_unitig_materialisation(unitigs_path, card, chrom, terminus, c_len,
                                      profile_hits=None):
    """The outboard sequence of the unitig Module 2 endorsed for this end, or (None, why not)."""
    unitig_id = card.get("authorized_unitig_id")
    if not unitig_id:
        return None, "the card endorses no unitig for this end"
    region = (card.get("evidence", {}) or {}).get("unitig_extension_region") or {}
    unitig_seq = extract_unitig_by_id(unitigs_path, unitig_id) if unitigs_path else None
    if not unitig_seq:
        return None, "unitig sequence not found in the provided FASTA"
    if region:
        net_seq = region_extension(unitig_seq, region)
        if not net_seq:
            return None, "the recorded extension interval lies outside the unitig"
    else:
        cap_hit = None
        if profile_hits is not None:
            hits = [h for h in profile_hits.get(chrom, []) if h.get("unitig_id") == unitig_id]
            if hits:
                if terminus == "5p":
                    cap_hit = min(hits, key=lambda h: h.get("t_start", 0))
                else:
                    cap_hit = max(hits, key=lambda h: h.get("t_end", 0))
        if not cap_hit:
            return None, ("the card records no unitig_extension_region and no telomere profile hit "
                          "locates the unitig's alignment, so the extension cannot be measured "
                          "(the unitig's own length is not the extension)")
        net_seq = outboard_extension(unitig_seq, cap_hit, terminus, c_len)
        if not net_seq:
            return None, "the unitig lies entirely within the assembly, so it adds no new sequence"

    return {
        "unitig_id": unitig_id,
        "net_seq": net_seq,
        "unitig_seq": unitig_seq,
        "unitig_len": len(unitig_seq),
    }, None


def align_seed_to_unitigs(seed_seq, unitigs_path, motif, rev_motif, min_overlap=300, min_telo_repeats=15, threads=1):
    """Align outward-oriented tip seed sequence against unitigs FASTA and return extending sequence if it carries telomere."""
    if not unitigs_path or not os.path.exists(unitigs_path) or not seed_seq:
        return None

    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as qf:
        qf_path = qf.name
        qf.write(f">seed\n{seed_seq}\n")

    try:
        cmd = ["minimap2", "-x", "map-hifi", "-c", "--secondary=no", "-t", str(threads), unitigs_path, qf_path]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not res.stdout.strip():
            return None

        best_cand = None
        for line in res.stdout.splitlines():
            p = line.strip().split("\t")
            if len(p) < 12:
                continue
            try:
                qlen, qstart, qend = int(p[1]), int(p[2]), int(p[3])
                strand = p[4]
                tname, tlen, tstart, tend = p[5], int(p[6]), int(p[7]), int(p[8])
                mapq = int(p[11])
                aln_len = int(p[10])
            except (ValueError, IndexError):
                continue

            if aln_len < min_overlap or mapq < 20:
                continue

            # Ensure alignment reaches near the outward tip of seed (3' end of seed_seq)
            if (qlen - qend) > 50:
                continue

            unitig_seq = extract_unitig_by_id(unitigs_path, tname)
            if not unitig_seq:
                continue

            if strand == "+":
                ext = unitig_seq[tend:]
            else:
                ext = revcomp(unitig_seq[:tstart])

            if not ext or len(ext) < 50:
                continue

            k = len(motif)
            fwd_h = ext[-1000:].count(motif)
            rev_h = ext[-1000:].count(rev_motif)
            telo_bp = max(fwd_h, rev_h) * k

            if telo_bp >= min_telo_repeats * k:
                if best_cand is None or telo_bp > best_cand["telo_bp"]:
                    best_cand = {
                        "unitig_id": tname,
                        "extension_seq": ext,
                        "added_bp": len(ext),
                        "telo_bp": telo_bp,
                        "mapq": mapq,
                        "strand": strand
                    }
        return best_cand
    finally:
        if os.path.exists(qf_path):
            try:
                os.remove(qf_path)
            except OSError:
                pass
