#!/usr/bin/env python3
"""Audit the ends of a curated assembly using aligned reads and a satellite family panel."""

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")
MATCH_OPERATIONS = {0, 7, 8}


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as source:
        name, chunks = None, []
        for line in source:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_annotations(path, lengths):
    intervals = defaultdict(list)
    with open_text(path) as source:
        for line in source:
            if line.startswith("#") or not line.strip():
                continue
            chromosome, _, feature, start, end, _, _, _, attributes = line.rstrip().split("\t")
            if feature != "satellite_array":
                continue
            start, end = int(start) - 1, int(end)
            if chromosome not in lengths or not 0 <= start < end <= lengths[chromosome]:
                raise ValueError(f"Satellite annotation outside assembly: {chromosome}:{start}-{end}")
            labels = dict(field.split("=", 1) for field in attributes.split(";") if "=" in field)
            intervals[chromosome].append((start, end, labels["Family"]))
    return intervals


def read_breaksites(path, lengths):
    sites = defaultdict(list)
    if path is not None:
        with open_text(path) as source:
            for line in source:
                if line.startswith("#") or not line.strip():
                    continue
                chromosome, start, end, label, _, strand, *_ = line.split()
                if "Type:S2G" not in label or strand not in {"+", "-"}:
                    continue
                start, end = int(start), int(end)
                # DELFIES represents a cut at the physical boundary as [-1, 0)
                # on the left or [length, length + 1) on the right.
                position = start if strand == "+" else end
                if chromosome not in lengths or start >= end or not 0 <= position <= lengths[chromosome]:
                    raise ValueError(f"Breakpoint outside assembly: {chromosome}:{start}-{end}")
                sites[chromosome, strand].append(position)
    return sites


def candidate_sequences(assembly, chromosomes, graphs):
    for name in sorted(assembly.keys() - set(chromosomes)):
        yield "assembly", name, assembly[name]
    for graph in graphs:
        with open_text(graph) as source:
            for line in source:
                if line.startswith("S\t"):
                    _, name, sequence, *_ = line.rstrip().split("\t")
                    if sequence != "*":
                        yield graph.name, name, sequence.upper()


def marker_panel(monomers, size):
    """Use cyclic monomer k-mers exact in one annotated family, on either strand."""
    families = defaultdict(set)
    for name, sequence in monomers.items():
        cyclic = sequence * ((len(sequence) + size - 1) // len(sequence) + 1)
        for i in range(len(sequence)):
            mer = cyclic[i:i + size]
            if set(mer) <= {"A", "C", "G", "T"}:
                families[mer].add(name)
                families[reverse_complement(mer)].add(name)
    return {mer: next(iter(names)) for mer, names in families.items() if len(names) == 1}


def repeat_runs(sequence, markers, periods, size, offset=0):
    """Yield repeat marker footprints, counting each observed base at most once per run."""
    runs = {}
    for pos in range(offset, len(sequence) - size + 1):
        family = markers.get(sequence[pos:pos + size])
        if family is None:
            continue
        start, stop, covered = runs.get(family, (pos, pos, 0))
        if pos - stop > max(size, periods[family]):
            if covered >= max(100, 2 * periods[family]):
                yield family, start, stop, covered
            start, stop, covered = pos, pos, 0
        end = pos + size
        runs[family] = start, end, covered + end - max(pos, stop)
    for family, (start, stop, covered) in runs.items():
        if covered >= max(100, 2 * periods[family]):
            yield family, start, stop, covered


def anchor_intervals(start, end, masked):
    current = start
    for left, right in sorted(masked):
        if left > current:
            yield current, min(left, end)
        current = max(current, right)
        if current >= end:
            return
    if current < end:
        yield current, end


def query_anchor(read, intervals, inward_from_right, minimum):
    """Locate an inward anchor in one unmasked reference interval using aligned query bases."""
    aligned = []
    reference, query = read.reference_start, 0
    for operation, size in read.cigartuples:
        if operation in MATCH_OPERATIONS:
            aligned.append((reference, reference + size, query))
        if operation in MATCH_OPERATIONS or operation in {2, 3}:
            reference += size
        if operation in MATCH_OPERATIONS or operation in {1, 4}:
            query += size
    ordered = reversed(intervals) if inward_from_right else intervals
    for start, end in ordered:
        blocks = reversed(aligned) if inward_from_right else aligned
        remaining = minimum
        for ref_start, ref_end, query_start in blocks:
            left, right = max(start, ref_start), min(end, ref_end)
            if left >= right:
                continue
            if remaining <= right - left:
                if inward_from_right:
                    return len(read.query_sequence) - (query_start + right - ref_start - remaining)
                return query_start + left - ref_start + remaining
            remaining -= right - left
    return None


def terminal_pattern(motif):
    rotations = {strand[i:] + strand[:i] for strand in (motif, reverse_complement(motif))
                 for i in range(len(motif))}
    return re.compile("|".join(f"(?:{word}){{3,}}" for word in sorted(rotations)))


def terminal_tract(sequence, pattern):
    # Search only the physical tip, allowing up to 25 bases after the exact tract.
    start = max(0, len(sequence) - 20000)
    return max((match.end() - match.start() for match in pattern.finditer(sequence, start)
                if len(sequence) - match.end() <= 25), default=0)


def write_tsv(path, columns, rows):
    with path.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def audit(args):
    import pysam

    assembly = dict(fasta_records(args.assembly))
    lengths = {name: len(sequence) for name, sequence in assembly.items()}
    chromosomes = args.chromosomes.split(",")
    if len(set(chromosomes)) != len(chromosomes) or set(chromosomes) - lengths.keys():
        raise ValueError("Chromosomes must be distinct FASTA record names")
    monomers = dict(fasta_records(args.monomers))
    annotations = read_annotations(args.arrays, lengths)
    breaksites = read_breaksites(args.breakpoints, lengths)
    markers = marker_panel(monomers, args.kmer)
    periods = {name: len(sequence) for name, sequence in monomers.items()}
    pattern = terminal_pattern(args.telomere.upper())
    ends, read_rows, terminal_rows, sequence_rows = [], [], [], []
    support = defaultdict(dict)
    candidates = defaultdict(set)
    assessable = 0

    with pysam.AlignmentFile(args.bam, index_filename=args.bai, require_index=True) as bam:
        if dict(zip(bam.references, bam.lengths)) != lengths:
            raise ValueError("BAM sequence dictionary does not match assembly names and lengths")
        md5_verified = 0
        for record in bam.header.to_dict()["SQ"]:
            if "M5" in record:
                digest = hashlib.md5(assembly[record["SN"]].encode()).hexdigest()
                if digest != record["M5"].lower():
                    raise ValueError(f"BAM reference checksum mismatch for {record['SN']}")
                md5_verified += 1

        for chromosome in chromosomes:
            sequence = assembly[chromosome]
            for side in ("5p", "3p"):
                left = side == "5p"
                end_id = f"{chromosome}:{side}"
                window = min(args.window, len(sequence) // 2)
                start, end = (0, window) if left else (len(sequence) - window, len(sequence))
                arrays = [(max(start, s), min(end, e), f) for s, e, f in annotations[chromosome]
                          if s < end and e > start]
                gaps = [(start + match.start(), start + match.end())
                        for match in re.finditer("N+", sequence[start:end])]
                anchors = list(anchor_intervals(start, end, [(s, e) for s, e, _ in arrays] + gaps))
                oriented = reverse_complement(sequence[start:end]) if left else sequence[start:end]
                placed_cap = terminal_tract(oriented, pattern)
                sites = breaksites[chromosome, "-" if left else "+"]
                anchored, cap_reads, extend_reads = set(), set(), set()
                for read in bam.fetch(chromosome, start, end):
                    if (read.is_secondary or read.is_supplementary or read.is_duplicate
                            or read.is_qcfail or read.mapping_quality < args.min_mapq
                            or not read.query_sequence or not read.cigartuples):
                        continue
                    offset = query_anchor(read, anchors, left, args.min_anchor)
                    if offset is None:
                        continue
                    molecule = read.query_name
                    anchored.add(molecule)
                    query = reverse_complement(read.query_sequence) if left else read.query_sequence
                    tract = terminal_tract(query, pattern)
                    operation, clip = read.cigartuples[0 if left else -1]
                    outward_clip = clip if operation == 4 else 0
                    distance = read.reference_start if left else len(sequence) - read.reference_end
                    reaches_tip = tract and distance <= 25 and outward_clip <= 25
                    endpoint = read.reference_start if left else read.reference_end
                    nearest = min((abs(site - endpoint) for site in sites), default=None) if tract else None
                    near_break = nearest is not None and nearest <= 50
                    extends_tip = tract and distance <= 100 and outward_clip - distance >= 100 and not near_break
                    if reaches_tip:
                        cap_reads.add(molecule)
                    if extends_tip:
                        extend_reads.add(molecule)
                    if tract:
                        terminal_rows.append((end_id, molecule, read.mapping_quality,
                                              read.reference_start, read.reference_end,
                                              tract, distance, outward_clip, nearest,
                                              int(bool(reaches_tip)), int(bool(extends_tip))))
                    for family, first, stop, covered in repeat_runs(query, markers, periods, args.kmer, offset):
                        group = support[family, end_id]
                        group[molecule] = max(covered, group.get(molecule, 0))
                        # A distal tract is attributed to this molecule only if it follows the repeat footprint.
                        distal_tract = tract if len(query) - tract >= stop else 0
                        read_rows.append((end_id, family, molecule, read.mapping_quality,
                                          offset, first, stop, covered, distal_tract))
                status = "unresolved"
                if anchored:
                    assessable += 1
                if placed_cap and cap_reads and not gaps:
                    status = "terminal_tract_supported"
                if extend_reads or gaps:
                    status = "review"
                if not anchored:
                    status = "insufficient_anchor"
                ends.append((end_id, len(sequence), start, end, placed_cap, len(anchored),
                             len(cap_reads), len(extend_reads), len(gaps), status,
                             ",".join(sorted({family for _, _, family in arrays}))))

        # Independent candidates are matched by repeat family. No placement is inferred here.
        for source, name, sequence in candidate_sequences(assembly, chromosomes, args.candidate_gfa):
            window = min(args.window, len(sequence))
            for side, oriented in (("5p", reverse_complement(sequence[:window])),
                                   ("3p", sequence[-window:])):
                tract = terminal_tract(oriented, pattern)
                if not tract:
                    continue
                for family, first, stop, covered in repeat_runs(oriented, markers, periods, args.kmer):
                    if stop <= len(oriented) - tract:
                        candidate_id = f"{source}:{name}:{side}"
                        candidates[family].add(candidate_id)
                        sequence_rows.append((candidate_id, family, covered, tract, len(oriented) - stop - tract))

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "ends.tsv",
              ("end", "sequence_bp", "window_start", "window_end", "assembly_exact_tip_repeat_bp",
               "anchored_molecules", "tip_tract_molecules", "outward_telo_overhang_molecules", "n_runs",
               "status", "annotated_array_families"), ends)
    write_tsv(args.outdir / "read_repeat_evidence.tsv",
              ("end", "family", "read_id", "mapq", "outward_anchor_end", "repeat_start", "repeat_end",
               "diagnostic_repeat_bp", "distal_exact_telo_bp"), sorted(read_rows))
    write_tsv(args.outdir / "terminal_read_evidence.tsv",
              ("end", "read_id", "mapq", "alignment_start", "alignment_end", "exact_telo_bp",
               "bp_to_assembly_tip", "outward_softclip_bp", "bp_to_oriented_s2g",
               "reaches_tip", "extends_tip"),
              sorted(terminal_rows))
    write_tsv(args.outdir / "unplaced_ending_candidates.tsv",
              ("candidate", "family", "diagnostic_repeat_bp", "exact_tip_repeat_bp", "repeat_to_telo_bp"),
              sorted(sequence_rows))
    end_repeat_rows = []
    count_by_family = Counter()
    for (family, end_id), molecules in sorted(support.items()):
        count = len(molecules)
        counted = count >= args.min_molecules
        count_by_family[family] += counted
        spans = sorted(molecules.values(), reverse=True)
        supported_bp = spans[args.min_molecules - 1] if counted else 0
        end_repeat_rows.append((family, end_id, count, spans[0], supported_bp, int(counted)))
    write_tsv(args.outdir / "end_repeat_links.tsv",
              ("family", "end", "molecules", "max_observed_diagnostic_bp",
               "replicated_diagnostic_bp_lower_bound", "counted"), end_repeat_rows)
    marker_counts = Counter(markers.values())
    write_tsv(args.outdir / "repeat_families.tsv",
              ("family", "monomer_bp", "diagnostic_kmers", "counted_ends", "assessable_ends",
               "surveyed_ends", "compatible_unplaced_endings"),
              ((family, periods[family], marker_counts[family], count_by_family[family],
                assessable, len(ends), ",".join(sorted(candidates[family])))
               for family in sorted(monomers)))
    metadata = vars(args).copy()
    with open(args.assembly, "rb") as source:
        sha = hashlib.sha256()
        for block in iter(lambda: source.read(1024 * 1024), b""):
            sha.update(block)
    metadata.update(assembly_file_sha256=sha.hexdigest(), bam_reference_md5_verified_records=md5_verified,
                    pysam_version=pysam.__version__)
    (args.outdir / "audit.json").write_text(json.dumps(metadata, default=str, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--bai", required=True)
    parser.add_argument("--arrays", type=Path, required=True, help="SatelliteProfiler arrays GFF3")
    parser.add_argument("--monomers", type=Path, required=True, help="SatelliteProfiler family monomers FASTA")
    parser.add_argument("--breakpoints", type=Path, help="DELFIES breakpoint_locations.bed")
    parser.add_argument("--candidate-gfa", type=Path, action="append", default=[],
                        help="Optional sequence-bearing GFA1 to search for compatible distal candidates")
    parser.add_argument("--chromosomes", required=True, help="Comma-separated assembly record names")
    parser.add_argument("--telomere", required=True)
    parser.add_argument("--window", type=int, default=200000)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-anchor", type=int, default=1500)
    parser.add_argument("--min-molecules", type=int, default=2)
    parser.add_argument("--kmer", type=int, default=15)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch("[ACGTacgt]+", args.telomere) or min(
            args.window, args.min_anchor, args.min_molecules, args.kmer) <= 0 or not 0 <= args.min_mapq <= 60:
        parser.error("Invalid motif, MAPQ, or evidence window/threshold")
    audit(args)


if __name__ == "__main__":
    main()
