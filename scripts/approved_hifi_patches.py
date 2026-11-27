"""Apply approved BAM-resolver sequence patches to uniquely located assembly gaps."""

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from package_reviewed_assembly import assembly_intervals, read_fasta, sha256


DNA = re.compile(r"[ACGT]+")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def candidate_patches(manifests, payloads):
    manifest = next((p for p in manifests if p.name.endswith("_patches_manifest.tsv")), None)
    if manifest is None:
        return {}
    if payloads is None:
        raise ValueError("HiFi gap manifest requires its patch FASTA directory")
    candidates = {}
    with manifest.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["outcome"] != "PATCH_GAP_INSERTION":
                continue
            name = row["gap_id"]
            if not name or name in candidates:
                raise ValueError(f"Missing or duplicate HiFi gap ID: {name!r}")
            patch = payloads / Path(row["patch_fasta"]).name
            records = list(read_fasta(patch)) if patch.is_file() else []
            seq = records[0][1].upper() if len(records) == 1 else ""
            supported = bool(DNA.fullmatch(seq)) and len(seq) == int(row["patch_len"])
            candidates[name] = {
                "gap_id": name, "record": row["chrom"], "start": int(row["gap_start"]), "end": int(row["gap_end"]),
                "left_anchor_len": int(row["left_anchor_len"]), "right_anchor_len": int(row["right_anchor_len"]),
                "patch_sha256": hashlib.sha256(seq.encode("ascii")).hexdigest() if supported else None,
                "patch_bp": int(row["patch_len"]), "sequence": seq if supported else None,
            }
    return candidates


def read_approvals(path, reference_digest):
    if path is None:
        return {}
    document = json.loads(path.read_text())
    if document["version"] != 1 or document["reference_sha256"] != reference_digest:
        raise ValueError("HiFi patch approvals do not match the curated reference or approval format")
    approvals = {}
    for entry in document["patches"]:
        name, digest = entry["gap_id"], entry["patch_sha256"]
        if not name or name in approvals or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"Invalid or repeated HiFi patch approval: {name!r}")
        approvals[name] = digest
    return approvals


def unique_anchor(records, query):
    hits = []
    for strand, pattern in (("+", query), ("-", reverse_complement(query))):
        for record, target in records.items():
            pos = target.find(pattern)
            if pos >= 0:
                if hits or target.find(pattern, pos + 1) >= 0:
                    return None
                hits.append((record, strand, pos))
    return hits[0] if hits else None


def lift_patch(patch, reference, targets, plans):
    chrom, start, stop = (patch[key] for key in ("record", "start", "end"))
    seq = reference.get(chrom, "")
    left_bp, right_bp = patch["left_anchor_len"], patch["right_anchor_len"]
    if (
        not 0 < start < stop < len(seq) or seq[start:stop] != "N" * (stop - start)
        or seq[start - 1] == "N" or seq[stop] == "N"
    ):
        return {}, "not_an_internal_reference_gap"
    if not 0 < left_bp <= start or not 0 < right_bp <= len(seq) - stop:
        return {}, "missing_reference_flanks"
    left, right = seq[start - left_bp:start], seq[stop:stop + right_bp]
    if not DNA.fullmatch(left) or not DNA.fullmatch(right):
        return {}, "ambiguous_reference_flanks"
    ref_left, ref_right = unique_anchor(reference, left), unique_anchor(reference, right)
    if not (ref_left == (chrom, "+", start - left_bp) or ref_right == (chrom, "+", stop)):
        return {}, "no_unique_anchor_in_curated_reference"
    a, b = unique_anchor(targets, left), unique_anchor(targets, right)
    if a is None and b is None:
        return {}, "no_unique_anchor_in_accepted_assembly"
    anchored = a or b
    if (a and b and a[:2] != b[:2]) or anchored[0] != chrom:
        return {}, "flanks_on_different_records_or_strands"
    strand = anchored[1]
    first, last, replacement = (left, right, patch["sequence"]) if strand == "+" else (
        reverse_complement(right), reverse_complement(left), reverse_complement(patch["sequence"]),
    )
    # A repeat flank can be disambiguated by adjacency to the gap and the unique anchor.
    pattern = first + f"(N+|{replacement})" + last
    match = re.search(pattern, targets[chrom])
    if not match:
        return {}, "accepted_sequence_between_flanks_or_flanks_changed"
    left_bound, right_bound = match.span(1)
    placement = {"record": chrom, "start": left_bound, "end": right_bound, "strand": strand}
    current = match.group(1)
    if current == replacement:
        return placement, "sequence_already_present"
    gap = next((part.get("gap") for part in plans[chrom]["parts"]
                if part["start"] == left_bound and part["end"] == right_bound), None)
    if gap not in ("unknown", "source_gap"):
        return placement, "protected_or_inexact_gap"
    return {**placement, "sequence": replacement}, "compatible_gap"


def reference_layout(sequences, reference_digest):
    objects = {}
    for name, seq in sequences.items():
        parts = []
        for start, end, gap in assembly_intervals(seq, None):
            part = {"start": start, "end": end}
            if gap:
                part["gap"] = gap
            else:
                part.update(source="reference", record=name, strand="+", source_start=start, source_end=end)
            parts.append(part)
        objects[name] = {"parts": parts}
    return {"version": 1, "coordinate_system": "zero-based half-open in the oriented source",
            "sources": {"reference": reference_digest}, "objects": objects}


def materialize(base_records, layout, placements, final_fasta):
    replacements = {(patch["placement"]["record"], patch["placement"]["start"]): patch for patch in placements}
    with final_fasta.open("w") as handle:
        for record, seq in base_records.items():
            offset = 0
            pieces, parts = [], []
            for part in layout["objects"][record]["parts"]:
                start, stop = part["start"], part["end"]
                patch = replacements.get((record, start))
                if patch:
                    placement = patch["placement"]
                    piece = placement.pop("sequence")
                    part = {"source": "hifi_gap_patch", "record": patch["gap_id"], "strand": placement["strand"],
                            "source_start": 0, "source_end": len(piece), "source_sha256": patch["patch_sha256"]}
                    patch["final_interval"] = {"record": record, "start": offset, "end": offset + len(piece)}
                else:
                    piece = seq[start:stop]
                parts.append({**part, "start": offset, "end": offset + len(piece)})
                pieces.append(piece)
                offset += len(piece)
            layout["objects"][record]["parts"] = parts
            handle.write(f">{record}\n")
            final_seq = "".join(pieces)
            for start in range(0, len(final_seq), 60):
                handle.write(final_seq[start:start + 60] + "\n")


def apply_approved_patches(reference, base, plan_path, manifests, payloads, approval_path, work, evidence):
    candidates = candidate_patches(manifests, payloads)
    reference_digest = sha256(reference)
    approvals = read_approvals(approval_path, reference_digest)
    if not candidates and not approvals and not manifests:
        return base, plan_path, []
    if candidates:
        with (evidence / "hifi_patch_sequences.fasta").open("w") as handle:
            for name, patch in sorted(candidates.items()):
                if patch["sequence"]:
                    handle.write(f">{name}\n{patch['sequence']}\n")
    inventory = {"version": 1, "reference_sha256": reference_digest, "patches": []}
    compatible = []
    if approvals:
        reference_records = {name: seq.upper() for name, seq in read_fasta(reference)}
        base_records = dict(read_fasta(base))
        targets = {name: seq.upper() for name, seq in base_records.items()}
        layout = json.loads(plan_path.read_text()) if plan_path else reference_layout(base_records, reference_digest)
    for name in sorted(candidates.keys() | approvals.keys()):
        candidate = candidates.get(name, {"gap_id": name})
        patch = {key: value for key, value in candidate.items() if key != "sequence"}
        digest = approvals.get(name)
        if digest is None:
            patch.update(status="review_required", reason="awaiting_approval" if patch["patch_sha256"] else "invalid_patch_fasta")
        elif digest != patch.get("patch_sha256"):
            patch.update(status="withheld", reason="approved_sequence_not_in_resolver_output", approved_patch_sha256=digest)
        else:
            location, reason = lift_patch(candidate, reference_records, targets, layout["objects"])
            patch.update(status="withheld", reason=reason, placement=location)
            if reason == "compatible_gap":
                compatible.append(patch)
            elif reason == "sequence_already_present":
                patch["status"] = "already_present"
        inventory["patches"].append(patch)
    counts = Counter((p["placement"]["record"], p["placement"]["start"]) for p in compatible)
    applied = []
    for patch in compatible:
        placement = patch["placement"]
        if counts[placement["record"], placement["start"]] > 1:
            placement.pop("sequence")
            patch.update(status="withheld", reason="multiple_approved_patches_for_gap")
        else:
            patch.update(status="applied", reason="approved_patch_unique_flanks")
            applied.append(patch)
    if applied:
        selected = work / "patched_assembly.fa"
        materialize(base_records, layout, applied, selected)
        layout["base_fasta_sha256"] = sha256(base)
        layout["fasta_sha256"] = sha256(selected)
        layout["sources"]["hifi_gap_patch"] = sha256(evidence / "hifi_patch_sequences.fasta")
        patched_plan = work / "patched_assembly_plan.json"
        patched_plan.write_text(json.dumps(layout, indent=2) + "\n")
    else:
        selected, patched_plan = base, plan_path
    (evidence / "hifi_gap_patches.json").write_text(json.dumps(inventory, indent=2) + "\n")
    return selected, patched_plan, inventory["patches"]
