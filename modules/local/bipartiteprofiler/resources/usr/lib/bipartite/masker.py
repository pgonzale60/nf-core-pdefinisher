"""Core motif instance locator (via FIMO) and sequence masker."""

from collections import defaultdict
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any


def load_soma_len_map(fasta_path: str | Path, default: int = 150) -> dict[str, int]:
    """Map clean sequence name to soma_len from FASTA header."""
    soma_map = {}
    if not os.path.exists(fasta_path):
        return soma_map
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split()[0].replace(":", "_")
                m = re.search(r"\bsoma_len=(\d+)", line)
                soma_map[seq_id] = int(m.group(1)) if m else default
    return soma_map


def load_fasta_dict(fasta_path: str | Path) -> dict[str, str]:
    """Load FASTA into dict: header -> sequence."""
    seqs = {}
    if not os.path.exists(fasta_path):
        return seqs
    with open(fasta_path, encoding="utf-8") as f:
        cur_h = None
        cur_s = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_h:
                    seqs[cur_h] = "".join(cur_s)
                # Keep ID clean (e.g. up to space)
                cur_h = line[1:].split()[0].replace(":", "_")
                cur_s = []
            else:
                cur_s.append(line.upper())
        if cur_h:
            seqs[cur_h] = "".join(cur_s)
    return seqs


def run_fimo(
    meme_file: str | Path,
    fasta_file: str | Path,
    out_dir: str | Path,
    threshold: float = 1e-4,
) -> list[dict[str, Any]]:
    """Execute FIMO to find motif occurrences passing threshold."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fimo_tsv = out_path / "fimo.tsv"

    if not fimo_tsv.exists():
        fimo_bin = shutil.which("fimo")
        if not fimo_bin:
            raise RuntimeError("Required binary 'fimo' not found in PATH")
        cmd = [
            fimo_bin,
            "--thresh",
            str(threshold),
            "--oc",
            str(out_path),
            str(meme_file),
            str(fasta_file),
        ]
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    hits = []
    if fimo_tsv.exists():
        with open(fimo_tsv, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                s = r.get("sequence_name")
                if not s or s.startswith("#"):
                    continue
                try:
                    hits.append({
                        "seq_name": s.replace(":", "_"),
                        "motif_id": r.get("motif_id", ""),
                        "start": int(r["start"]),  # 1-based inclusive
                        "stop": int(r["stop"]),  # 1-based inclusive
                        "strand": r.get("strand", "+"),
                        "score": float(r.get("score", 0.0)),
                        "p_value": float(r.get("p-value", 1.0)),
                        "matched_seq": r.get("matched_sequence", ""),
                    })
                except (ValueError, KeyError):
                    continue

    return hits


def mask_sequences_with_hits(
    seq_dict: dict[str, str], hits: list[dict[str, Any]]
) -> dict[str, str]:
    """Mask hits with 'N' in each sequence."""
    hits_by_seq = defaultdict(list)
    for h in hits:
        hits_by_seq[h["seq_name"]].append(h)

    masked_seqs = {}
    for name, seq in seq_dict.items():
        char_list = list(seq)
        for h in hits_by_seq.get(name, []):
            start_0 = max(0, h["start"] - 1)
            stop_1 = min(len(char_list), h["stop"])
            for i in range(start_0, stop_1):
                char_list[i] = "N"
        masked_seqs[name] = "".join(char_list)

    return masked_seqs


def mask_core_motifs(
    core_meme: str | Path,
    fasta_file: str | Path,
    output_masked_fasta: str | Path,
    fimo_out_dir: str | Path,
    threshold: float = 1e-4,
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    """Run FIMO with core motif on FASTA, mask instances with N, and save masked FASTA."""
    seq_dict = load_fasta_dict(fasta_file)
    hits = run_fimo(
        core_meme, fasta_file, out_dir=fimo_out_dir, threshold=threshold
    )
    masked_dict = mask_sequences_with_hits(seq_dict, hits)

    with open(output_masked_fasta, "w", encoding="utf-8") as f:
        for name, mseq in masked_dict.items():
            f.write(f">{name}\n{mseq}\n")

    return masked_dict, hits
