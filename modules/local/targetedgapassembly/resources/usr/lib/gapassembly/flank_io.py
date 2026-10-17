"""Flank extraction, inverted-repeat anchor selection, and BAM read baiting for Targeted Gap Assembly.

British spelling standard enforced across all documentation and comments.
"""

import json
import os
import subprocess

GAP_ACTIONS = ("PATCH_GAP_INSERTION", "REQUIRES_UNITIG_PATCH", "INSERT_SEQUENCE")
CONSENSUS_ACTIONS = ("PATCH_GAP_INSERTION", "INSERT_SEQUENCE")
FLANK_RECORDS = ("LEFT", "RIGHT")

_INDEXABLE_REF_CACHE = {}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


revcomp = reverse_complement


def require_indexable_reference(ref_fasta):
    """Return the reference, having proved that samtools can resolve regions from it."""
    if not ref_fasta or not os.path.exists(ref_fasta):
        raise RuntimeError(f"reference {ref_fasta} does not exist")
    if ref_fasta in _INDEXABLE_REF_CACHE:
        return _INDEXABLE_REF_CACHE[ref_fasta]
    try:
        res = subprocess.run(["samtools", "faidx", ref_fasta], capture_output=True, text=True)
    except OSError as e:
        raise RuntimeError(f"cannot run samtools to index {ref_fasta}: {e}")
    if res.returncode != 0 or not os.path.exists(ref_fasta + ".fai"):
        err = (res.stderr or "").strip().replace("\n", " ")
        raise RuntimeError(
            f"reference {ref_fasta} is not indexable, so no region can be read from it. This engine requires "
            f"the bgzipped assembly produced by BGZIPFASTA, with its .gzi staged alongside it. "
            f"samtools said: {err}"
        )
    _INDEXABLE_REF_CACHE[ref_fasta] = ref_fasta
    return ref_fasta


def read_fasta_records(path):
    """Read FASTA into {name: sequence} in file order."""
    records = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise RuntimeError(f"{path} does not start with a FASTA header")
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def read_fastq_records(path):
    """Read FASTQ into {name: (sequence, qualities)} in file order."""
    records = {}
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise RuntimeError(f"malformed FASTQ in {path}: {header.strip()[:80]}")
            seq = fh.readline().rstrip("\n")
            separator = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not separator.startswith("+") or len(qual) != len(seq):
                raise RuntimeError(
                    f"malformed FASTQ record {header.strip()[:80]}: {len(seq)} bases, {len(qual)} qualities"
                )
            records[header[1:].split()[0]] = (seq, qual)
    return records


def parse_gaps_from_input(blueprint_json=None, gaps_bed=None):
    """Extract list of gap dicts strictly conforming to Module 2 curation blueprint cards."""
    gaps = []
    if blueprint_json and os.path.exists(blueprint_json):
        with open(blueprint_json) as f:
            data = json.load(f)
        cards = data if isinstance(data, list) else data.get("action_cards", [])
        for c in cards:
            action = c.get("action")
            if action not in GAP_ACTIONS:
                continue
            gaps.append(
                {
                    "gap_id": c["card_id"],
                    "chrom": c["chrom"],
                    "start": int(c["start"]),
                    "end": int(c["end"]),
                    "length": int(c["length"]),
                    "action": action,
                    "evidence": c.get("evidence", {}),
                    "action_source": "card",
                }
            )
    elif gaps_bed and os.path.exists(gaps_bed):
        with open(gaps_bed) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                gaps.append(
                    {
                        "gap_id": f"gap_{chrom}_{start}_{end}",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "length": end - start,
                        "action": "REQUIRES_UNITIG_PATCH",
                        "evidence": {},
                        "action_source": "default_no_card",
                    }
                )
    return gaps


def extract_fasta_region(ref_fasta, chrom, start, end, out_fa):
    """Extract a genomic region with samtools faidx, or fail loudly."""
    ref_fasta = require_indexable_reference(ref_fasta)
    start_1 = max(1, start)
    region = f"{chrom}:{start_1}-{end}"
    res = subprocess.run(["samtools", "faidx", ref_fasta, region], capture_output=True, text=True)
    if res.returncode != 0 or not res.stdout.strip():
        err = (res.stderr or "").strip().replace("\n", " ")
        raise RuntimeError(f"could not read flank {region} from {ref_fasta}: {err}")
    with open(out_fa, "w") as out:
        out.write(res.stdout)
    return True


def longest_inverted_repeat(window, max_k=None):
    """Longest k for which the window holds a k-mer that also occurs in its own reverse complement."""
    if not window:
        return 0
    needle = window.upper()
    complement = reverse_complement(window).upper()
    lo, hi = 0, (max_k or len(window) // 2)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        step = max(1, mid // 4)
        if any(needle[i : i + mid] in complement for i in range(0, len(needle) - mid + 1, step)):
            lo = mid
        else:
            hi = mid - 1
    return lo


def choose_anchor(window, gap_at_end, min_anchor_bp=1000, max_inverted_repeat_bp=500):
    """Shrink flank window toward the gap until its self-inverted repeat cannot mislead alignment."""
    n = len(window)
    while n > min_anchor_bp:
        anchor = window[-n:] if gap_at_end else window[:n]
        if longest_inverted_repeat(anchor) <= max_inverted_repeat_bp:
            break
        n = max(min_anchor_bp, n // 2)
    return window[-n:] if gap_at_end else window[:n]


def write_flank_reference(
    ref_fasta, chrom, gap_start, gap_end, flank_len, out_fa, min_anchor_bp=1000, max_inverted_repeat_bp=500
):
    """Write the two flank anchors as separate records, LEFT and RIGHT, and return them."""
    ref_fasta = require_indexable_reference(ref_fasta)
    windows = (
        ("LEFT", max(1, gap_start - flank_len + 1), gap_start, True),
        ("RIGHT", gap_end + 1, gap_end + flank_len, False),
    )
    anchors = {}
    with open(out_fa, "w") as out:
        for name, start, end, gap_at_end in windows:
            res = subprocess.run(
                ["samtools", "faidx", ref_fasta, f"{chrom}:{start}-{end}"],
                capture_output=True,
                text=True,
            )
            if res.returncode != 0 or not res.stdout.strip():
                err = (res.stderr or "").strip().replace("\n", " ")
                raise RuntimeError(f"could not read flank {chrom}:{start}-{end} from {ref_fasta} -- {err}")
            window = "".join(line.strip() for line in res.stdout.splitlines() if not line.startswith(">"))
            anchor = choose_anchor(
                window, gap_at_end, min_anchor_bp=min_anchor_bp, max_inverted_repeat_bp=max_inverted_repeat_bp
            )
            anchors[name] = anchor
            out.write(
                f">{name} flank_bp={len(window)} anchor_bp={len(anchor)} inverted_repeat_bp={longest_inverted_repeat(anchor)}\n{anchor}\n"
            )
    return anchors


def bait_reads_for_gap(bam_path, chrom, gap_start, gap_end, flank, out_reads, fmt="fasta"):
    """Extract reads overlapping the gap flanks into a deduplicated FASTA or FASTQ file."""
    if fmt not in ("fasta", "fastq"):
        raise RuntimeError(f"unsupported bait format: {fmt}")
    r1 = f"{chrom}:{max(0, gap_start - flank)}-{gap_start}"
    r2 = f"{chrom}:{gap_end}-{gap_end + flank}"

    view = subprocess.Popen(
        ["samtools", "view", "-F", "0x900", "-b", bam_path, r1, r2],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    convert = subprocess.Popen(
        ["samtools", fmt, "-"],
        stdin=view.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    view.stdout.close()
    out, convert_err = convert.communicate()
    view_err = view.stderr.read().decode() if view.stderr else ""
    if view.wait() != 0 or convert.returncode != 0:
        err = (view_err or convert_err or "").strip()[:300]
        raise RuntimeError(f"samtools could not extract reads around {chrom}:{gap_start} -- {err}")
    if not out.strip():
        return 0

    seen = set()
    with open(out_reads, "w") as fh:
        if fmt == "fasta":
            qname = None
            chunks = []
            for line in out.splitlines():
                if line.startswith(">"):
                    if qname and qname not in seen:
                        seen.add(qname)
                        fh.write(f">{qname}\n{''.join(chunks)}\n")
                    qname = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            if qname and qname not in seen:
                seen.add(qname)
                fh.write(f">{qname}\n{''.join(chunks)}\n")
        else:
            lines = out.splitlines()
            if len(lines) % 4 != 0:
                raise RuntimeError(f"samtools fastq returned {len(lines)} lines, not a multiple of 4")
            for i in range(0, len(lines), 4):
                header, sequence, separator, qualities = (
                    lines[i],
                    lines[i + 1],
                    lines[i + 2],
                    lines[i + 3],
                )
                if not header.startswith("@"):
                    raise RuntimeError(f"malformed FASTQ record: {header[:60]}")
                qname = header[1:].split()[0]
                if qname in seen:
                    continue
                seen.add(qname)
                fh.write(f"{header}\n{sequence}\n{separator}\n{qualities}\n")
    return len(seen)


def gfa_to_fasta(gfa_path, out_fa):
    """Extract segment sequences from GFA to FASTA."""
    if not os.path.exists(gfa_path):
        return 0
    count = 0
    with open(gfa_path) as gfa, open(out_fa, "w") as fa:
        for line in gfa:
            if line.startswith("S\t"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    fa.write(f">{parts[1]}\n{parts[2]}\n")
                    count += 1
    return count
