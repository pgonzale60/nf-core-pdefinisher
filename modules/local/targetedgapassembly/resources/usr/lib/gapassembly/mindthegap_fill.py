"""MindTheGap fill integration for Tier-2 multi-read gap bridging.

Performs bounded de Bruijn graph path traversal between left and right flank anchors
when single-read spanning coverage is zero (n_spanning == 0).
"""

import os
import shutil
import subprocess
from .flank_io import read_fasta_records


def run_mindthegap_fill(reads_fq, flank_left_seq, flank_right_seq, gap_work, kmer_size=63, max_nodes=500):
    """Run MindTheGap fill between left and right flank breakpoint anchors.

    Returns the assembled bridge sequence or None if unbridged.
    """
    mtg_bin = shutil.which("mindthegap")
    if not mtg_bin:
        raise RuntimeError("mindthegap binary not found in PATH. mindthegap is required for Tier-2 gap filling.")

    bkpt_fa = os.path.join(gap_work, "gap_breakpoints.fa")
    with open(bkpt_fa, "w", encoding="utf-8") as f:
        f.write(f">LEFT\n{flank_left_seq}\n>RIGHT\n{flank_right_seq}\n")

    out_prefix = os.path.join(gap_work, "mtg_fill")
    cmd = [
        mtg_bin,
        "fill",
        "-in", reads_fq,
        "-bkpt", bkpt_fa,
        "-kmer-size", str(kmer_size),
        "-max-nodes", str(max_nodes),
        "-out", out_prefix,
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return None

    # MindTheGap outputs filled insertions to <prefix>.insertions.fasta or <prefix>.fasta
    possible_outputs = [
        f"{out_prefix}.insertions.fasta",
        f"{out_prefix}.insertions.fa",
        f"{out_prefix}.fasta",
    ]
    for out_fa in possible_outputs:
        if os.path.exists(out_fa):
            records = read_fasta_records(out_fa)
            for header, seq in records.items():
                if seq and len(seq) > 0:
                    return seq

    return None
