#!/usr/bin/env python3
"""Measure what the reads and the assembly itself say about each chromosome end.

Two independent signals are reported per end, and keeping them apart is the point of this module:

  1. `n_telomeric_reads` / `max_telo_count` -- telomeric repeats in *reads* whose soft clip sits at
     the tip. This is what the reads say. It says nothing about the draft's own tip sequence: a read
     can carry a deep telomeric tail while the assembled tip ends in unique subtelomeric sequence.
     That is the DRAFT_HAS_TELOMERE versus TELOMERE_IN_OVERHANG distinction Module 2 has to make and
     cannot make from reads alone.

  2. `terminal_tract_bp` -- telomeric repeat in the *assembly's own* terminal sequence. This is what
     the draft says. It is measured here because this is the only stage that has the assembly:
     Module 2 classifies an end but receives evidence, never sequence.

Conflating the two is a real defect. DF5120 IV 5p and X 5p carry 84-150 hexamers in MAPQ 1-14 read
clips that stop *inside* the assembly (overhang_bp = 0, so no read extends the tip), and an
unqualified read-based status called both ends ALREADY_CAPPED. The walker then found no tract at
either tip, and both ends finished UNRESOLVED while Module 2's own blueprint endorsed a candidate
for each. The status now names the situation instead of asserting a cap:

    EXTEND_AND_CAP        telomeric reads at the tip, at least one extending past it
    ALREADY_CAPPED        telomeric reads at the tip AND a canonical tract in the assembly's tip
    TELOMERE_IN_OVERHANG  telomeric reads at the tip, but no tract in the assembly's tip: the draft
                          is not capped, the telomere is in the reads' overhang tails, and the end
                          needs an overhang extension rather than a tract extension
    UNCAPPED_EXTENDABLE   no telomeric reads, but a soft clip extends past the tip
    NO_CLIPS              nothing to extend with

`proposed_status` is consumed by pde_curation_synthesis.py, which decides the action. This module
measures; it does not decide.
"""
import argparse
import sys

SCHEMA_VERSION = "1.1.0"

# Reads are collected over 50 kb because a read whose clip reaches the tip may start far from it.
WINDOW = 50000

# The tract scan sees the terminal 10 kb, matching the policy's terminal_window_bp.
TRACT_WINDOW = 10000

# An assembly rarely stops exactly on a repeat boundary, so the scan allows the final partial unit
# to be up to this many bp into the array before the tandem run is counted. Requiring the final base
# to complete a unit would report 0 bp for most genuinely capped ends.
TIP_OFFSET_MAX_BP = 5

# TELOMERE_IN_OVERHANG is the conservative label: it claims only what the reads show and never
# asserts a cap, so a mislabelled end loses nothing but a card's wording.
#
# A telomeric read carries at least this many canonical units in its clip. Two units (12 bp) is
# deliberately permissive -- this module measures, Module 2 weighs -- and it is why the MAPQ column
# is kept in the soft-clip table: DF5120 IV 5p's clips are MAPQ 1-14, and a caller that wants to
# discount them can.
MIN_TELO_UNITS_IN_READ = 2

# A clip shorter than this is noise unless it extends past the tip.
MIN_CLIP_BP = 30

# An overhang of at least this many bp marks a read as genuinely extending the end.
MIN_OVERHANG_BP = 50

def revcomp(seq, table):
    return seq.translate(table)[::-1]


def terminal_tract_bp(terminal_seq, motif, rev_motif, tip_offset_max_bp=TIP_OFFSET_MAX_BP):
    """Length of the canonical tandem run anchored at the 3' end of `terminal_seq`.

    The sequence must arrive in chromosome orientation with the tip last. Every canonical variant is
    in frame with every other, so the run is read backwards in whole units, allowing the tip to stop
    `tip_offset_max_bp` inside a unit. Scattered substitutions are not tolerated: this is the
    conservative test that decides whether the draft can be said to hold a tract at all, and
    TELOMERE_IN_OVERHANG (the answer when it cannot) asserts nothing.

    This answers the taxonomy's question for a terminus -- does the tip contain canonical repeats
    (EXTEND_TELOMERE_TRACT) or is it unique sequence (EXTEND_TO_TELOMERE) -- so the threshold is
    presence, not length. Module 3's terminal_tract_region measures a tract it is about to extend and
    applies the policy's min_repeats instead, which is why a five-unit tip measures 12 bp here and 0
    there. Both are right; check_telomere_classification.py holds them to that at matched thresholds.
    """
    variants = [v.upper() for v in dict.fromkeys([motif, rev_motif, "TTAGAC", "GTCTAA"]) if v]
    unit = len(motif)
    seq = terminal_seq.upper()

    best = 0
    for offset in range(0, min(tip_offset_max_bp, len(seq)) + 1):
        end = len(seq) - offset
        count = 0
        while end - unit >= 0 and seq[end - unit:end] in variants:
            count += 1
            end -= unit
        if count >= 2:
            if offset == 0:
                # In frame from the very tip: unambiguous, and no offset can beat it.
                return count * unit
            best = max(best, count * unit)
    return best


class FastaReader:
    """Wrapper around pysam.FastaFile with transparent fallback for non-bgzipped gzipped FASTA."""

    def __init__(self, path):
        import pysam
        self._fasta = None
        self._dict = None
        try:
            self._fasta = pysam.FastaFile(path)
        except Exception:
            self._load(path)

    def _load(self, path):
        import gzip
        opener = gzip.open if str(path).endswith(".gz") else open
        self._dict = {}
        curr_name = None
        chunks = []
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    if curr_name is not None:
                        self._dict[curr_name] = "".join(chunks)
                    curr_name = line[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            if curr_name is not None:
                self._dict[curr_name] = "".join(chunks)

    def fetch(self, ref, start, end):
        if self._fasta is not None:
            return self._fasta.fetch(ref, start, end)
        seq = self._dict.get(ref, "")
        return seq[start:end]

    def close(self):
        if self._fasta is not None:
            self._fasta.close()
        self._dict = None


def tip_sequence(fasta, ref, length, side="3p", window=TRACT_WINDOW):
    """The terminal `window` bp of `ref` in chromosome orientation, tip last."""
    if side == "5p":
        try:
            # For 5p, index 0 is the physical tip. Reverse it so the tip is last,
            # matching terminal_tract_bp's expectation of reading backwards from the tip.
            seq = fasta.fetch(ref, 0, min(window, length)).upper()
            return seq[::-1]
        except Exception as exc:  # pragma: no cover - depends on the FASTA
            sys.stderr.write(f"Warning: could not read 5p tip sequence for {ref}: {exc}\n")
            return ""
    else:
        start = max(0, length - window)
        try:
            return fasta.fetch(ref, start, length).upper()
        except Exception as exc:  # pragma: no cover - depends on the FASTA
            sys.stderr.write(f"Warning: could not read 3p tip sequence for {ref}: {exc}\n")
            return ""


def scan_end(bam, ref, length, side, motif, rev_motif, fasta):
    """Reads at one end, the tract in the assembly's own tip, and the status that follows."""
    if side == "5p":
        fetch_start, fetch_end = 0, min(WINDOW, length)
    else:
        fetch_start, fetch_end = max(0, length - WINDOW), length

    records = []
    try:
        for r in bam.fetch(ref, fetch_start, fetch_end):
            # Secondary and supplementary alignments are excluded: an end must be established by
            # primary alignments. What that means for overhang evidence is the walker's business.
            if r.is_secondary or r.is_supplementary or not r.cigartuples:
                continue
            if side == "5p":
                cigar_terminal = r.cigartuples[0]
                if cigar_terminal[0] != 4:
                    continue
                clip_len = cigar_terminal[1]
                align_start = r.reference_start
                overhang_bp = max(0, clip_len - align_start)
                clip_seq = (r.query_sequence or "")[:clip_len].upper()
                records.append((r.query_name, ref, "5p", align_start, r.reference_end,
                                clip_len, overhang_bp, clip_seq, r.mapping_quality))
            else:
                cigar_terminal = r.cigartuples[-1]
                if cigar_terminal[0] != 4:
                    continue
                clip_len = cigar_terminal[1]
                align_end = r.reference_end
                overhang_bp = max(0, clip_len - (length - align_end))
                clip_seq = (r.query_sequence or "")[-clip_len:].upper()
                records.append((r.query_name, ref, "3p", r.reference_start, align_end,
                                clip_len, overhang_bp, clip_seq, r.mapping_quality))
    except Exception as exc:
        sys.stderr.write(f"Warning: fetch error on {ref} {side}: {exc}\n")

    kept = []
    for rec in records:
        clip_len, overhang_bp, clip_seq = rec[5], rec[6], rec[7]
        telo_count = clip_seq.count(motif) + clip_seq.count(rev_motif)
        if clip_len >= MIN_CLIP_BP or overhang_bp > 0:
            kept.append((rec[0], rec[1], rec[2], rec[3], rec[4], clip_len, overhang_bp,
                         telo_count, rec[8]))

    n_telo = sum(1 for x in kept if x[7] >= MIN_TELO_UNITS_IN_READ)
    max_oh = max([x[6] for x in kept], default=0)
    max_telo_cnt = max([x[7] for x in kept], default=0)
    tract_bp = terminal_tract_bp(tip_sequence(fasta, ref, length, side), motif, rev_motif)

    status = "NO_CLIPS"
    if n_telo >= 1 and max_oh >= MIN_OVERHANG_BP:
        status = "EXTEND_AND_CAP"
    elif n_telo >= 1 and tract_bp > 0:
        # Reads carry telomere and the draft's own tip ends in a tandem run: the end holds a tract,
        # and the action is to extend that tract.
        status = "ALREADY_CAPPED"
    elif n_telo >= 1:
        # Reads carry telomere but the draft's tip does not. The draft is NOT capped: the telomere
        # is in the reads' overhang tails, so the end needs an overhang extension, not a tract
        # extension. Calling this ALREADY_CAPPED is what let the walker conclude there was nothing
        # to retain and still finish unresolved.
        status = "TELOMERE_IN_OVERHANG"
    elif max_oh >= MIN_OVERHANG_BP:
        status = "UNCAPPED_EXTENDABLE"

    summary = (f"{ref}_{side}", length, len(kept), n_telo, max_telo_cnt, status, tract_bp)
    return kept, summary
def write_outputs(prefix, clip_records, end_summaries):
    with open(f"{prefix}.terminal_softclips.tsv", "w") as f:
        f.write("read_id\tchrom\tside\talign_start\talign_end\tclip_len"
                "\toverhang_bp\ttelo_count\tmapq\n")
        for rec in clip_records:
            f.write("\t".join(map(str, rec)) + "\n")

    with open(f"{prefix}.terminal_ends_summary.tsv", "w") as f:
        # terminal_tract_bp is APPENDED, never inserted: pde_curation_synthesis.py parses this table
        # positionally (it takes the first six fields), so a column in the middle would shift every
        # field after it by one and silently corrupt every baseline status.
        f.write("chrom_end\tchrom_length\tn_softclipped_reads\tn_telomeric_reads\tmax_telo_count"
                "\tproposed_status\tterminal_tract_bp\n")
        for s in end_summaries:
            f.write("\t".join(map(str, s)) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument("--bam", required=True, help="Coordinate-sorted HiFi BAM")
    parser.add_argument("--assembly", required=True,
                        help="The assembly the BAM was aligned to; its tip sequence is measured")
    parser.add_argument("--motif", default="TTAGGC", help="Canonical telomeric repeat unit")
    parser.add_argument("--prefix", required=True, help="Output filename prefix")
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCHEMA_VERSION}")
    args = parser.parse_args(argv)

    motif = args.motif.upper()
    rc_table = str.maketrans("ACGTN", "TGCAN")
    rev_motif = revcomp(motif, rc_table)

    # Imported here rather than at module scope so the measurement functions above can be imported
    # and checked without the bio-stack on the path -- which is what the repository's harness does,
    # and why no other engine here imports pysam either.
    import pysam

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = FastaReader(args.assembly)

    clip_records = []
    end_summaries = []
    for ref in bam.references:
        length = bam.get_reference_length(ref)
        if length < 1000:
            continue
        for side in ("5p", "3p"):
            records, summary = scan_end(bam, ref, length, side, motif, rev_motif, fasta)
            clip_records.extend(records)
            end_summaries.append(summary)

    bam.close()
    fasta.close()
    write_outputs(args.prefix, clip_records, end_summaries)
    return 0


if __name__ == "__main__":
    sys.exit(main())
