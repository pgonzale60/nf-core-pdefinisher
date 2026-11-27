#!/usr/bin/env python3
"""Run the curated-assembly graph finishing pilot and publish its evidence and assembly."""

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from package_reviewed_assembly import package, read_fasta, sha256


GAP = re.compile(r"N+")
PORT = re.compile(r"([^,+\-\[\]]+)([+-])(?:\[[^\]]*\])?")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ARCHIVE_FILES = (
    "consensus.fasta", "families.tsv", "haplotype_patterns.tsv.gz", "membership.tsv.gz",
    "end_records.tsv.gz", "locus_support.tsv", "MANIFEST.yaml",
)


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_rows(path, fields, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def export_graph(gfa, fasta):
    opener = gzip.open if gfa.suffix == ".gz" else open
    with opener(gfa, "rt") as handle, fasta.open("w") as dest:
        for line in handle:
            if line.startswith("S\t"):
                _, name, sequence, *_ = line.rstrip().split("\t")
                if sequence == "*":
                    raise ValueError(f"Graph segment {name} has no sequence")
                dest.write(f">{name}\n{sequence}\n")


def align(target, query, paf, preset, threads):
    with paf.open("w") as handle:
        subprocess.run([
            "minimap2", "-c", "-x", preset, "-t", str(threads), "-N", "50", "--secondary=yes", target, query,
        ], stdout=handle, check=True)


def paf_rows(path):
    with path.open() as handle:
        for line in handle:
            query, qlen, qstart, qend, strand, target, tlen, start, end, matches, span, mapq, *_ = line.split()
            yield (query, int(qlen), int(qstart), int(qend), strand, target, int(tlen),
                   int(start), int(end), int(matches), int(span), int(mapq))


def families_on_graph(archive, targets, out, threads):
    """Project family consensuses and keep all passing placements with molecule provenance."""
    consensus = archive / "consensus.fasta"
    families = {row["cluster_id"]: row for row in rows(archive / "families.tsv")}
    if set(families) != {name for name, _ in read_fasta(consensus)}:
        raise ValueError("Family table and consensus FASTA must contain the same cluster IDs")
    patterns = defaultdict(list)
    for row in rows(archive / "haplotype_patterns.tsv.gz"):
        patterns[row["cluster_id"]].append(row)
    observations = {row["end_id"]: row["read_id"] for row in rows(archive / "end_records.tsv.gz")}
    molecules = defaultdict(set)
    for row in rows(archive / "membership.tsv.gz"):
        molecules[row["cluster_id"]].add(observations[row["end_id"]])
    if (set(patterns) | set(molecules)) - set(families):
        raise ValueError("Haplotype or membership table refers to an unknown family")
    config = json.loads((archive / "MANIFEST.yaml").read_text())["parameters"]
    min_bp, min_identity = config["mapping_min_bp"], config["mapping_identity"]
    placements = []
    for label, target in targets.items():
        paf = out / f"families_to_{label}.paf"
        align(target, consensus, paf, "asm20", threads)
        for family, _, qstart, qend, strand, record, _, start, end, matches, span, mapq in paf_rows(paf):
            if span >= min_bp and matches / span >= min_identity:
                placements.append(dict(
                    cluster_id=family, source=label, record=record, strand=strand, start=start, end=end,
                    consensus_start=qstart, consensus_end=qend, matches=matches, alignment_bp=span, mapq=mapq,
                    family_molecules=len(molecules[family]),
                    patterns=",".join(sorted(row["pattern_id"] for row in patterns[family])),
                ))
    for name in ARCHIVE_FILES:
        shutil.copyfile(archive / name, out / name)
    return placements, config


def report_candidates(catalog, placements, curated_paf, config, out):
    by_unitig = defaultdict(list)
    for hit in placements:
        if hit["source"] == "primary":
            by_unitig[hit["record"]].append(hit)
    alignments = defaultdict(list)
    for hit in paf_rows(curated_paf):
        query, _, _, _, _, _, _, _, _, matches, span, _ = hit
        if span >= config["mapping_min_bp"] and matches / span >= config["mapping_identity"]:
            alignments[query].append(hit)

    annotations = []
    by_curated_end = defaultdict(list)
    for end in rows(catalog):
        name, side = end["object_id"], end["side"]
        candidates = []
        for _, qlen, qstart, qend, strand, chrom, clen, start, stop, _, _, mapq in alignments[name]:
            flip = strand == "-"
            is_left = (side == "5p") != flip
            target_side = "5p" if is_left else "3p"
            distance = start if is_left else clen - stop
            qdistance = qstart if side == "5p" else qlen - qend
            if distance <= config["reference_end_window_bp"]:
                candidates.append((chrom, target_side, strand, distance, qdistance, mapq))
                by_curated_end[chrom, target_side].append(end["end_id"])
        ports = {unitig for unitig, _ in PORT.findall(",".join((end["horizon_port"], end["candidate_endings"])))}
        support = sorted({(p["record"], p["cluster_id"]) for unitig in ports for p in by_unitig[unitig]})
        annotations.append({
            **end,
            "curated_end_candidates": ";".join(":".join(map(str, c)) for c in sorted(set(candidates))),
            "family_unitig_matches": ";".join(f"{unitig}:{family}" for unitig, family in support),
            "decision": "review_required" if candidates else "unplaced_graph_end",
        })
    with catalog.open() as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    write_rows(out, fields + ["curated_end_candidates", "family_unitig_matches", "decision"], annotations)
    return by_curated_end


def report_delivery(delivery, candidates, recipe, out):
    gaps = defaultdict(list)
    for gap in rows(delivery / "sequence_gaps.tsv"):
        gaps[gap["object"]].append(gap)
    records = []
    for end in rows(delivery / "ends.tsv"):
        chrom, side = end["chromosome"], end["side"]
        chrom_gaps = gaps[chrom]
        gap = (chrom_gaps[0] if side == "5p" else chrom_gaps[-1]) if chrom_gaps else None
        records.append({
            **end,
            "distal_cap_in_layout": "reviewed" if recipe else "unverified",
            "chromosome_sequence": "contains_placeholders" if gap else "gapless_in_layout",
            "nearest_sequence_gap": f"{gap['start']}-{gap['end']}" if gap else "",
            "nearest_gap_class": gap["gap_class"] if gap else "",
            "nearest_gap_linkage": gap["linkage_status"] if gap else "",
            "graph_end_candidates": ",".join(sorted(set(candidates[chrom, side]))),
        })
    write_rows(out, list(records[0]), records)


def source_records(path, ids):
    """Use seqkit to fetch requested records from FASTA or FASTQ without indexing all reads."""
    with tempfile.TemporaryDirectory() as work:
        keys, fasta = Path(work) / "keys", Path(work) / "source.fasta"
        keys.write_text("\n".join(sorted(ids)) + "\n")
        with subprocess.Popen(["seqkit", "grep", "-f", keys, path], stdout=subprocess.PIPE) as grep:
            with fasta.open("w") as output:
                subprocess.run(["seqkit", "fq2fa"], stdin=grep.stdout, stdout=output, check=True)
            if grep.wait():
                raise subprocess.CalledProcessError(grep.returncode, grep.args)
        records = dict(read_fasta(fasta))
        if ids - records.keys():
            raise ValueError(f"Missing source records: {sorted(ids - records.keys())}")
        return records


def trace_recipe(recipe_path, fasta, sources, evidence_sources, out, finish_agp=None):
    """Bind reviewed recipe slices to the exact FASTA emitted by gfapdefinisher."""
    recipe = json.loads(recipe_path.read_text())
    required = defaultdict(set)
    for chrom in recipe["chromosomes"]:
        for part in chrom["parts"]:
            if "source" in part:
                required[part["source"]].add(part["record"])
    sequences = {key: source_records(sources[key], ids) for key, ids in required.items()}
    result = dict(read_fasta(fasta))
    objects = {}
    gap_bp = recipe["gap_bp"]
    rdna_intervals = set()
    if finish_agp:
        with finish_agp.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if columns[4] in ("N", "U") and columns[6:9] == ["repeat", "yes", "rna_tangle"]:
                    rdna_intervals.add((columns[0], int(columns[1]) - 1, int(columns[2])))
    for chrom in recipe["chromosomes"]:
        name = chrom["name"]
        pieces = []
        parts = []
        offset = 0
        for piece in chrom["parts"]:
            source_length = 0
            if "gap" in piece:
                seq, start = "N" * gap_bp, 0
            else:
                seq = sequences[piece["source"]][piece["record"]].upper()
                source_length = len(seq)
                if piece.get("strand", "+") == "-":
                    seq = seq.translate(COMPLEMENT)[::-1]
                start = piece.get("start", 0)
                seq = seq[start:piece.get("end", len(seq))]
            pos = 0
            for match in GAP.finditer(seq):
                left, right = match.span()
                if left > pos:
                    offset = append_source(parts, pieces, piece, seq[pos:left], start + pos, offset)
                if not (parts and "gap" in parts[-1] and left == 0):
                    kind = piece.get("gap", "source_gap")
                    if piece.get("source") == "finish":
                        a, b = start + left, start + right
                        if piece.get("strand") == "-":
                            a, b = source_length - b, source_length - a
                        if (piece["record"], a, b) in rdna_intervals:
                            kind = "rdna"
                    parts.append(dict(start=offset, end=offset + gap_bp, gap=kind))
                    pieces.append("N" * gap_bp)
                    offset += gap_bp
                pos = right
            if pos < len(seq):
                offset = append_source(parts, pieces, piece, seq[pos:], start + pos, offset)
        if "".join(pieces) != result[name]:
            raise ValueError(f"Recipe source coordinates fail to reconstruct {name}")
        objects[name] = {"parts": parts, "caps": {side: chrom[f"cap_{side}"] for side in ("5p", "3p")}}
    out.write_text(json.dumps({
        "version": 1, "coordinate_system": "zero-based half-open in the oriented source",
        "fasta_sha256": sha256(fasta), "recipe_sha256": sha256(recipe_path),
        "sources": {key: sha256(path) for key, path in sorted(evidence_sources.items()) if path is not None},
        "objects": objects,
    }, indent=2) + "\n")


def append_source(parts, pieces, part, seq, start, offset):
    stop = offset + len(seq)
    parts.append(dict(start=offset, end=stop, source=part["source"], record=part["record"],
                      strand=part.get("strand", "+"), source_start=start, source_end=start + len(seq)))
    pieces.append(seq)
    return stop


def graph_args(command, **values):
    for name, value in values.items():
        if value is not None:
            command.extend(("--" + name.replace("_", "-"), value))
    subprocess.run(command, check=True)


def run_pilot(args, out):
    evidence = out / "evidence"
    evidence.mkdir()
    with tempfile.TemporaryDirectory(dir=out) as workspace:
        graphs = {}
        for label, path in (("primary", args.gfa), ("raw", args.raw_gfa)):
            if path:
                target = Path(workspace) / f"{label}.fa"
                export_graph(path, target)
                graphs[label] = target
        placements, config = families_on_graph(
            args.end_archive, {**graphs, "curated": args.curated}, evidence, args.threads,
        )
        fields = ("cluster_id source record strand start end consensus_start consensus_end matches alignment_bp mapq "
                  "family_molecules patterns").split()
        write_rows(evidence / "family_placements.tsv", fields, sorted(placements, key=lambda hit: tuple(hit[f] for f in fields)))

        graph = out / "graph"
        ont_paf = args.ont_paf
        if args.ont_reads and not ont_paf:
            ont_paf = evidence / "ont_to_primary.paf"
            align(graphs["primary"], args.ont_reads, ont_paf, "map-ont", args.threads)
        graph_args([args.gfapdefinisher, "finish", "--nominate-super"], gfa=args.gfa, tiara=args.tiara,
                   barrnap=args.barrnap, ont_paf=ont_paf, ont_reads=args.ont_reads, out_dir=graph)
        scaffolds = graph / "scaffolds.fasta"
        scaffold_paf = evidence / "scaffolds_to_curated.paf"
        align(args.curated, scaffolds, scaffold_paf, "asm20", args.threads)
        candidates = report_candidates(
            graph / "end_catalog.tsv", placements, scaffold_paf, config, evidence / "graph_end_proposals.tsv",
        )

        assembly, plan = args.curated, None
        if args.recipe:
            assembly, plan = Path(workspace) / "reviewed.fa", Path(workspace) / "assembly_plan.json"
            graph_args([args.gfapdefinisher, "curate"], recipe=args.recipe, reference=args.curated, gfa=args.gfa,
                       raw_gfa=args.raw_gfa, scaffolds=scaffolds, ont_reads=args.ont_reads,
                       hifi_reads=args.hifi_reads, out=assembly)
            sources = dict(primary=args.gfa, raw=args.raw_gfa, finish=scaffolds,
                           reference=args.curated, ont=args.ont_reads, hifi=args.hifi_reads)
            finish_agp = graph / "scaffolds.agp"
            trace_recipe(args.recipe, assembly, {**sources, **graphs}, dict(sources, finish_agp=finish_agp),
                         plan, finish_agp)
        chromosomes = set(args.chromosomes.split(","))
        delivery = out / "assembly"
        package(assembly, delivery, chromosomes, args.telomere.upper(), 3, args.recipe, plan)
        report_delivery(delivery, candidates, args.recipe, evidence / "curated_end_review.tsv")
        inputs = dict(curated=args.curated, gfa=args.gfa, raw_gfa=args.raw_gfa, recipe=args.recipe,
                      tiara=args.tiara, barrnap=args.barrnap, ont_paf=ont_paf,
                      ont_reads=args.ont_reads, hifi_reads=args.hifi_reads)
        (out / "inputs.json").write_text(json.dumps({
            "mode": "reviewed_replay" if args.recipe else "discovery",
            "sha256": {key: sha256(path) for key, path in inputs.items() if path is not None},
            "end_archive_sha256": {name: sha256(args.end_archive / name) for name in ARCHIVE_FILES},
        }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("curated", "gfa", "end-archive", "outdir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name in ("raw-gfa", "ont-reads", "ont-paf", "hifi-reads", "tiara", "barrnap", "recipe"):
        parser.add_argument(f"--{name}", type=Path)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument("--telomere", default="TTAGGC")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--gfapdefinisher", default="gfapdefinisher")
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error(f"Output already exists: {args.outdir}")
    if args.threads < 1 or not re.fullmatch("[ACGT]+", args.telomere.upper()):
        parser.error("Supply positive threads and a DNA telomere motif")
    args.outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.outdir.parent) as tmp:
        out = Path(tmp) / "results"
        out.mkdir()
        run_pilot(args, out)
        out.rename(args.outdir)


if __name__ == "__main__":
    main()
