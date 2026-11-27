"""Normalize graph and BAM observations into a reference-bound curation decision plan."""

import csv
import json

from package_reviewed_assembly import read_fasta, sha256


def interval(lengths, record, start, end, side=None):
    if record not in lengths:
        return None
    if start in (None, "") or end in (None, ""):
        if side not in ("5p", "3p"):
            return None
        start = end = 0 if side == "5p" else lengths[record]
    start, end = int(start), int(end)
    if not 0 <= start <= end <= lengths[record]:
        raise ValueError(f"Edit interval {record}:{start}-{end} is outside the curated reference")
    location = {"record": record, "start": start, "end": end}
    if side in ("5p", "3p"):
        location["side"] = side
    return location


def read_proposals(lengths, artifacts):
    files = {suffix: path for path in artifacts for suffix in (
        "curation_manifest.json", "patches_manifest.tsv", "walker_manifest.tsv", "rdna_manifest.tsv", "trans_manifest.tsv"
    ) if path.name.endswith(suffix)}
    if not artifacts:
        return []
    if len(files) != len(artifacts) or "curation_manifest.json" not in files:
        raise ValueError("Expected one blueprint and distinct BAM resolver manifests")
    blueprint = json.loads(files.pop("curation_manifest.json").read_text())
    cards = {}
    proposals = []
    for key in ("action_cards", "curation_cards", "telomere_cards"):
        for entry in blueprint.get(key, []):
            card_id = entry["card_id"]
            if card_id in cards:
                continue
            cards[card_id] = entry
            span = card_interval(lengths, entry)
            proposals.append(proposal("bam", f"blueprint:{card_id}", entry.get("action", "unclassified"),
                                      [span] if span else [], {"card_id": card_id, "observations": entry.get("evidence", {})}))
    for suffix, path in sorted(files.items()):
        with path.open() as handle:
            for index, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
                card_id = row.get("card_id", row.get("gap_id", ""))
                if suffix == "walker_manifest.tsv":
                    targets = [interval(lengths, row["chrom"], None, None, row["terminus"])]
                elif suffix == "patches_manifest.tsv":
                    targets = [interval(lengths, row["chrom"], row["gap_start"], row["gap_end"])]
                elif suffix == "trans_manifest.tsv":
                    targets = transfer_intervals(lengths, row)
                else:
                    targets = [card_interval(lengths, cards.get(card_id, {}))]
                observations = {key: value for key, value in row.items() if key and not key.endswith("_fasta")}
                proposals.append(proposal("bam", f"{suffix}:{index}", row.get("outcome") or row.get("status", "unclassified"),
                                          [target for target in targets if target],
                                          {"manifest": path.name, "card_id": card_id, "observations": observations}))
    return proposals


def card_interval(lengths, card):
    start = card.get("start") if card.get("start") is not None else card.get("start_pos")
    end = card.get("end") if card.get("end") is not None else card.get("end_pos")
    return interval(lengths, card.get("chrom"), start, end, card.get("terminus"))


def transfer_intervals(lengths, row):
    sources = [interval(lengths, row["source_chrom"], row.get(f"{prefix}_start"), row.get(f"{prefix}_end"))
               for prefix in ("source_gap", "contig")]
    source = [span for span in sources if span]
    if len(source) == 2:
        source = [interval(lengths, row["source_chrom"], min(s["start"] for s in source), max(s["end"] for s in source))]
    target = row.get("target_chrom")
    if target not in lengths or not row.get("target_coord"):
        return source
    site = int(row["target_coord"])
    feature = row.get("target_feature")
    if site <= 10000 or feature == "5_prime_terminus":
        source.append(interval(lengths, target, None, None, "5p"))
    elif abs(site - lengths[target]) <= 10000 or feature == "3_prime_terminus":
        source.append(interval(lengths, target, None, None, "3p"))
    else:
        source.append(interval(lengths, target, site, site))
    return source


def proposal(source, identifier, kind, targets, evidence):
    return {"id": f"{source}:{identifier}", "source": source, "kind": kind, "targets": targets, "evidence": evidence}


def catalog_candidates(row):
    """Tie the catalog's positional evidence columns to their candidate ports."""
    endings = row["candidate_endings"].split(";") if row["candidate_endings"] else []
    groups = []
    for name in ("candidate_families", "candidate_molecules"):
        cells = row[name].split(";") if row.get(name) else [""] * len(endings)
        if len(cells) != len(endings):
            raise ValueError(f"{row['end_id']}: {name} does not match the candidate ends")
        groups.append(cells)
    candidates = []
    for ending, families, molecules in zip(endings, *groups, strict=True):
        port, _, tags = ending.partition("[")
        candidates.append({
            "port": port,
            "evidence_tags": list(filter(None, tags.removesuffix("]").split(","))),
            "family_ids": families.split(",") if families else [],
            "reported_molecule_ids": molecules.split(",") if molecules else [],
        })
    return candidates


def graph_proposals(lengths, candidate_rows):
    proposals = []
    for row in candidate_rows:
        locations = {}
        alignments = []
        for hit in filter(None, row["curated_end_candidates"].split(";")):
            record, side, strand, ref_dist, graph_dist, mapq = hit.rsplit(":", 5)
            location = interval(lengths, record, None, None, side)
            if location:
                locations[record, side] = location
                alignments.append({"record": record, "side": side, "strand": strand, "mapq": int(mapq),
                                   "reference_tip_distance": int(ref_dist), "graph_tip_distance": int(graph_dist)})
        proposals.append(proposal("graph", row["end_id"], "graph_end_candidate", list(locations.values()), {
            "graph_category": row["category"], "horizon_port": row["horizon_port"],
            "competing_graph_ends": row["candidate_endings"], "family_unitig_matches": row["family_unitig_matches"],
            "candidate_ends": catalog_candidates(row),
            "certain_ports": row.get("certain_ports", ""), "certain_joins": row.get("certain_joins", ""),
            "reference_alignments": alignments,
        }))
    return proposals


def overlaps(left, right):
    if left["record"] != right["record"]:
        return False
    if left["start"] == left["end"] or right["start"] == right["end"]:
        return left["start"] <= right["end"] and right["start"] <= left["end"]
    return left["start"] < right["end"] and right["start"] < left["end"]


def review_reason(candidate):
    if candidate["source"] == "graph":
        if len(candidate["targets"]) > 1:
            return "competing_curated_end_placements"
        category = candidate["evidence"]["graph_category"]
        if category == "shared":
            return "shared_terminal_repeat"
        if category in ("bounded", "open"):
            return "germline_linkage_unresolved"
        return "reference_splice_not_validated"
    return "bam_edit_requires_review"


def decide(reference, graph_rows, bam_artifacts, reviewed, reviewed_recipe, output, nuclear_chromosomes=None):
    """Choose the reference or an explicit replay and record every candidate's disposition."""
    reference_records = dict(read_fasta(reference))
    lengths = {name: len(seq) for name, seq in reference_records.items()}
    reference_digest = sha256(reference)
    candidates = graph_proposals(lengths, graph_rows) + read_proposals(lengths, bam_artifacts)
    accepted = None
    selected = reference
    if reviewed:
        recipe = json.loads(reviewed_recipe.read_text())
        if recipe.get("reference_sha256", reference_digest) != reference_digest:
            raise ValueError("Reviewed recipe describes a different curated reference")
        reviewed_digest = sha256(reviewed)
        if reviewed_digest != recipe.get("expected_sha256"):
            raise ValueError("Reviewed replay requires a matching expected FASTA checksum")
        reviewed_records = dict(read_fasta(reviewed))
        required_records = lengths.keys() if nuclear_chromosomes is None else nuclear_chromosomes
        if missing := sorted(required_records - reviewed_records.keys()):
            raise ValueError(f"Reviewed replay omits curated records: {', '.join(missing)}")
        changed = [interval(lengths, record, 0, lengths[record]) for record in sorted(lengths.keys() & reviewed_records.keys())
                   if reviewed_records[record].upper() != reference_records[record].upper()]
        accepted = proposal("graph", "reviewed_replay", "reviewed_assembly", changed, {
            "recipe_sha256": sha256(reviewed_recipe), "fasta_sha256": reviewed_digest,
            "added_records": sorted(reviewed_records.keys() - lengths.keys()),
            "preserved_reference_records": sorted(lengths.keys() - reviewed_records.keys()),
        })
        candidates.append(accepted)
        selected = reviewed

    decisions = []
    for candidate in sorted(candidates, key=lambda item: item["id"]):
        if candidate is accepted:
            state, reason = "accepted", "checksum_verified_reviewed_recipe"
        elif not candidate["targets"]:
            state, reason = "unplaced", "no_curated_reference_interval"
        elif accepted and any(overlaps(a, b) for a in candidate["targets"] for b in accepted["targets"]):
            state, reason = "review_required", "within_reviewed_replay_scope"
        else:
            state, reason = "review_required", review_reason(candidate)
        decisions.append({**candidate, "decision": state, "reason": reason})
    output.write_text(json.dumps({
        "version": 1, "coordinate_system": "curated FASTA, zero-based half-open",
        "reference_sha256": reference_digest, "delivery_sha256": sha256(selected) if accepted else reference_digest,
        "accepted_proposal": accepted["id"] if accepted else None, "proposals": decisions,
    }, sort_keys=True, indent=2) + "\n")
    return selected
