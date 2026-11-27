#!/usr/bin/env python3
"""Package an accepted FASTA as AGP, component FASTA, and an end/gap inventory."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path


GAP = re.compile("N+", re.IGNORECASE)
MITOCHONDRIAL_IDS = {"chrmt", "chrm", "mt", "scaffold_mt", "mitochondrion", "mitochondria"}


def read_fasta(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as handle:
        name, parts = None, []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"Missing FASTA record ID in {path}")
                name, parts = header[0], []
            elif name is None:
                raise ValueError(f"Expected a FASTA header in {path}")
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def terminal_repeats(sequence, motif, side):
    if side == "3p":
        sequence = sequence[::-1]
        motif = motif[::-1]
    rotations = (motif[i:] + motif[:i] for i in range(len(motif)))
    patterns = [re.compile(f"(?:{repeat})+") for repeat in rotations]
    return max((match.end() - match.start()
                for shift in range(min(25, len(sequence)) + 1)
                for pattern in patterns
                if (match := pattern.match(sequence, shift))), default=0)


def assembly_intervals(seq, plan):
    if plan is not None:
        offset = 0
        for part in plan:
            start, end = part["start"], part["end"]
            gap = part.get("gap")
            if start != offset or not start < end <= len(seq):
                raise ValueError("Assembly plan parts must tile the FASTA record")
            bases = seq[start:end].upper()
            if (gap and bases.strip("N")) or (not gap and "N" in bases):
                raise ValueError("Assembly plan gap disagrees with FASTA")
            if not gap and (
                not part["source"] or not part["record"] or part["strand"] not in ("+", "-")
                or part["source_start"] < 0 or part["source_end"] - part["source_start"] != end - start
            ):
                raise ValueError("Invalid source interval in assembly plan")
            yield start, end, gap
            offset = end
        if offset != len(seq):
            raise ValueError("Assembly plan does not cover the FASTA record")
        return
    offset = 0
    for gap in GAP.finditer(seq):
        if offset < gap.start():
            yield offset, gap.start(), None
        yield gap.start(), gap.end(), "unknown"
        offset = gap.end()
    if offset < len(seq):
        yield offset, len(seq), None


def write_layout(fasta, out, chromosomes, motif, min_repeats, plan=None):
    reverse_motif = motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seen = set()
    with (
        (out / "assembly.agp").open("w") as agp,
        (out / "components.fasta").open("w") as components,
        (out / "ends.tsv").open("w") as ends,
        (out / "sequence_gaps.tsv").open("w") as gaps,
        (out / "unresolved_joins.tsv").open("w") as joins,
    ):
        agp.write("##agp-version 2.1\n")
        agp.write("# Components are slices of the supplied reviewed FASTA.\n")
        if plan is not None:
            agp.write("# Source intervals are recorded in assembly_plan.json.\n")
        end_rows = csv.writer(ends, delimiter="\t", lineterminator="\n")
        gap_rows = csv.writer(gaps, delimiter="\t", lineterminator="\n")
        join_rows = csv.writer(joins, delimiter="\t", lineterminator="\n")
        end_rows.writerow(("chromosome", "side", "canonical_repeat_bp", "physical_tip_status"))
        gap_header = ("object", "start", "end", "placeholder_bp", "placement", "agp_part",
                      "gap_class", "linkage_status", "sequence_status")
        gap_rows.writerow(gap_header)
        join_rows.writerow(gap_header)

        for record_number, (name, seq) in enumerate(read_fasta(fasta), 1):
            if name in seen or not seq or not re.fullmatch("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", seq):
                raise ValueError(f"Duplicate, empty, or invalid FASTA record: {name}")
            seen.add(name)
            if name in chromosomes:
                for side, end_motif in (("5p", reverse_motif), ("3p", motif)):
                    bp = terminal_repeats(seq.upper(), end_motif, side)
                    status = "canonical_tract" if bp >= min_repeats * len(motif) else "open_or_unconfirmed"
                    end_rows.writerow((name, side, bp, status))

            if plan is not None and name not in plan:
                raise ValueError(f"Assembly plan has no path for {name}")
            intervals = list(assembly_intervals(seq, plan[name]["parts"] if plan is not None else None))
            if not any(not gap for _, _, gap in intervals):
                raise ValueError(f"FASTA record contains only gap characters: {name}")

            for part, (start, stop, gap) in enumerate(intervals, 1):
                length = stop - start
                if gap:
                    terminal = start == 0 or stop == len(seq)
                    gap_type = "telomere" if terminal and name in chromosomes else "contig"
                    fields = ("N", length, gap_type, "no", "na")
                    gap_class, linkage, sequence_status = "unclassified", "unassessed", "sequence_unresolved"
                    if gap == "rdna":
                        gap_class, linkage, sequence_status = "rdna_array", "flanking_connection_supported", "sequence_and_copy_order_unresolved"
                        fields = ("N", length, "repeat", "yes", "unspecified")
                    elif gap == "satellite_port":
                        gap_class, linkage = "tandem_repeat", "distal_placement_unresolved"
                    row = (name, start, stop, length, "terminal" if terminal else "internal", part,
                           gap_class, linkage, sequence_status)
                    gap_rows.writerow(row)
                    if gap != "rdna":
                        join_rows.writerow(row)
                else:
                    component = f"component_{record_number}_{part}"
                    fields = ("W", component, 1, length, "+")
                    components.write(f">{component}\n")
                    for i in range(start, stop, 60):
                        components.write(seq[i:min(i + 60, stop)] + "\n")
                agp.write("\t".join(map(str, (name, start + 1, stop, part, *fields))) + "\n")

    if plan is not None and set(plan) != seen:
        raise ValueError("The assembly plan and FASTA must contain the same records")
    if not seen or chromosomes - seen:
        raise ValueError(f"Assembly is empty or is missing chromosome IDs: {sorted(chromosomes - seen)}")


def validate_layout(fasta, out):
    """Round-trip the written AGP and components against each original FASTA record."""
    components = iter(read_fasta(out / "components.fasta"))
    with (out / "assembly.agp").open() as handle:
        rows = (line.rstrip().split("\t") for line in handle if not line.startswith("#"))
        for name, seq in read_fasta(fasta):
            offset, part = 0, 1
            while offset < len(seq):
                row = next(rows, None)
                if row is None or row[0] != name or int(row[1]) != offset + 1 or int(row[3]) != part:
                    raise ValueError(f"AGP order or coordinates disagree with FASTA for {name}")
                stop = int(row[2])
                if row[4] == "W":
                    component = next(components, None)
                    if component is None or component[0] != row[5] or row[6:] != ["1", str(stop - offset), "+"]:
                        raise ValueError(f"AGP component coordinates disagree with FASTA for {name}")
                    assembled = component[1].upper()
                else:
                    assembled = "N" * int(row[5])
                if stop <= offset or stop > len(seq) or assembled != seq[offset:stop].upper():
                    raise ValueError(f"AGP does not reconstruct FASTA at {name}:{offset + 1}-{stop}")
                offset, part = stop, part + 1
        if next(rows, None) is not None or next(components, None) is not None:
            raise ValueError("AGP or component FASTA has extra rows")


def split_mitochondrial_fasta(fasta, staging, chromosomes):
    """Preserve original FASTA record bytes when splitting recognized mitochondrial IDs."""
    nuclear_path = staging / "assembly.fasta"
    mito_path = staging / "mitochondrial.fasta"
    mitochondrial_names = set()
    opener = gzip.open if fasta.suffix == ".gz" else open
    with opener(fasta, "rb") as source, nuclear_path.open("wb") as nuclear, mito_path.open("wb") as mito:
        destination = nuclear
        for line in source:
            if line.startswith(b">"):
                header = line[1:].split()
                if not header:
                    raise ValueError(f"Missing FASTA record ID in {fasta}")
                name = header[0].decode("ascii")
                if name not in chromosomes and name.casefold() in MITOCHONDRIAL_IDS:
                    if name in mitochondrial_names:
                        raise ValueError(f"Duplicate mitochondrial FASTA record: {name}")
                    mitochondrial_names.add(name)
                    destination = mito
                else:
                    destination = nuclear
            destination.write(line)
    if mitochondrial_names:
        mitochondrial_records = dict(read_fasta(mito_path))
        if any(not re.fullmatch("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", seq) for seq in mitochondrial_records.values()):
            raise ValueError("Empty or invalid mitochondrial FASTA record")
    else:
        if fasta.suffix == ".gz":
            nuclear_path.unlink()
        mito_path.unlink()
    return mitochondrial_names


def package(fasta, outdir, chromosomes, motif, min_repeats, reviewed_recipe=None, assembly_plan=None, reviewed_base=None):
    if outdir.exists():
        raise FileExistsError(outdir)
    original_digest = sha256(fasta)
    if reviewed_base and not (reviewed_recipe and assembly_plan):
        raise ValueError("A separate reviewed base requires its recipe and the final assembly plan")
    reviewed_digest = sha256(reviewed_base) if reviewed_base else original_digest
    if reviewed_recipe:
        recipe = json.loads(reviewed_recipe.read_text())
        if recipe.get("expected_sha256") != reviewed_digest:
            raise ValueError("Reviewed recipe checksum does not match the reviewed FASTA")
    plan = None
    if assembly_plan:
        plan = json.loads(assembly_plan.read_text())
        if plan["version"] != 1 or plan["fasta_sha256"] != original_digest:
            raise ValueError("Assembly plan version or FASTA checksum does not match")
        if reviewed_recipe and plan["recipe_sha256"] != sha256(reviewed_recipe):
            raise ValueError("The assembly plan was produced by a different recipe")
        if reviewed_base and plan.get("reviewed_fasta_sha256", plan.get("base_fasta_sha256")) != reviewed_digest:
            raise ValueError("The final assembly plan was produced from a different reviewed FASTA")

    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=outdir.parent) as temp:
        staging = Path(temp) / "package"
        staging.mkdir()
        mitochondrial_names = split_mitochondrial_fasta(fasta, staging, chromosomes)
        if mitochondrial_names:
            staged_fasta = staging / "assembly.fasta"
            if plan:
                plan["mitochondrial_objects"] = {
                    record: plan["objects"].pop(record) for record in sorted(mitochondrial_names)
                }
                plan["combined_fasta_sha256"] = original_digest
                plan["fasta_sha256"] = sha256(staged_fasta)
                plan["mitochondrial_fasta_sha256"] = sha256(staging / "mitochondrial.fasta")
        else:
            name = "assembly.fasta.gz" if fasta.suffix == ".gz" else "assembly.fasta"
            staged_fasta = staging / name
            if fasta.suffix == ".gz":
                shutil.copyfile(fasta, staged_fasta)
            if sha256(staged_fasta) != original_digest:
                raise ValueError("FASTA changed during packaging")
        write_layout(staged_fasta, staging, chromosomes, motif, min_repeats, plan["objects"] if plan else None)
        validate_layout(staged_fasta, staging)
        if reviewed_recipe:
            shutil.copyfile(reviewed_recipe, staging / "reviewed_recipe.json")
        if assembly_plan:
            if mitochondrial_names:
                (staging / "assembly_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
            else:
                shutil.copyfile(assembly_plan, staging / "assembly_plan.json")
        metadata = {
            "input_fasta": staged_fasta.name if plan else str(fasta.resolve()),
            "chromosomes": sorted(chromosomes),
            "telomere_repeat": motif,
            "min_telomere_repeats": min_repeats,
            "agp_provenance": "components are slices of the accepted FASTA",
            "gap_coordinates": "zero-based half-open",
            "checksums_sha256": {item.name: sha256(item) for item in sorted(staging.iterdir())},
        }
        if mitochondrial_names:
            metadata["mitochondrial_records"] = sorted(mitochondrial_names)
            metadata["combined_fasta_sha256"] = original_digest
        (staging / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
        staging.rename(outdir)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--chromosomes", required=True, help="Comma-separated FASTA record IDs")
    parser.add_argument("--telomere", required=True, help="Forward telomere motif")
    parser.add_argument("--min-telomere-repeats", type=int, default=3)
    parser.add_argument("--reviewed-recipe", type=Path, help="Optional provenance with the accepted FASTA SHA-256")
    parser.add_argument("--assembly-plan", type=Path, help="Reviewed component coordinates from gfapdefinisher")
    args = parser.parse_args()
    chromosomes = set(args.chromosomes.split(","))
    motif = args.telomere.upper()
    if "" in chromosomes or not re.fullmatch("[ACGT]+", motif) or args.min_telomere_repeats <= 0:
        parser.error("Supply nonempty chromosome IDs, a DNA motif and a positive repeat count")
    package(args.assembly, args.outdir, chromosomes, motif, args.min_telomere_repeats,
            args.reviewed_recipe, args.assembly_plan)


if __name__ == "__main__":
    main()
