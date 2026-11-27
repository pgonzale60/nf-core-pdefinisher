import hashlib
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from approved_hifi_patches import apply_approved_patches, reference_layout, reverse_complement
from curation_edit_plan import decide
from gfapdefinisher_handoff import retain_reference_records
from package_reviewed_assembly import package, read_fasta, sha256


def fixture(tmp_path, *, reverse=False, resolved=False, protected=False):
    rng = random.Random(42)
    anchored = lambda n: "".join(rng.choices("ACGT", k=n))
    left, right, caps, patch = anchored(40), anchored(40), anchored(32), anchored(14)
    reference = tmp_path / "reference.fa"
    reference.write_text(f">I\n{left}NNNN{right}\n")
    body = left + ("GGGG" if resolved else "NNNN") + right
    sequence = caps + (reverse_complement(body) if reverse else body)
    base = tmp_path / "graph.fa"
    base.write_text(f">I\n{sequence}\n")
    plan = tmp_path / "graph_plan.json"
    layout = reference_layout({"I": sequence}, sha256(reference))
    layout["fasta_sha256"] = sha256(base)
    layout["objects"]["I"]["caps"] = {"5p": {"source": "primary", "record": "cap"}}
    if protected:
        layout["objects"]["I"]["parts"][1]["gap"] = "rdna"
    recipe = tmp_path / "recipe.json"
    recipe.write_text(json.dumps({"expected_sha256": sha256(base)}))
    layout["recipe_sha256"] = sha256(recipe)
    plan.write_text(json.dumps(layout))
    manifests = tmp_path / "I_patches_manifest.tsv"
    manifests.write_text("gap_id\tchrom\tgap_start\tgap_end\toutcome\tleft_anchor_len\tright_anchor_len\tpatch_len\tpatch_fasta\n"
                         "gap1\tI\t40\t44\tPATCH_GAP_INSERTION\t20\t20\t14\t/from/nextflow/work/gap1.fa\n")
    payloads = tmp_path / "payloads"
    payloads.mkdir()
    (payloads / "gap1.fa").write_text(f">gap1\n{patch}\n")
    approvals = tmp_path / "approvals.json"
    approvals.write_text(json.dumps({"version": 1, "reference_sha256": sha256(reference), "patches": [
        {"gap_id": "gap1", "patch_sha256": hashlib.sha256(patch.encode("ascii")).hexdigest()},
    ]}))
    return reference, base, plan, [manifests], payloads, approvals, recipe, sequence, patch


@pytest.mark.parametrize("reverse", [False, True])
def test_patch_round_trips_alongside_reviewed_graph_cap(tmp_path, reverse):
    reference, base, plan, manifests, payloads, approvals, recipe, sequence, patch = fixture(tmp_path, reverse=reverse)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    fasta, final_plan, results = apply_approved_patches(reference, base, plan, manifests, payloads, approvals, tmp_path, evidence)
    assert dict(read_fasta(fasta))["I"] == sequence.replace("NNNN", reverse_complement(patch) if reverse else patch)
    assert results[0]["status"] == "applied"
    assert results[0]["final_interval"] == {"record": "I", "start": 72, "end": 86}
    package(fasta, tmp_path / "delivery", {"I"}, "TTAGGC", 3, recipe, final_plan, reviewed_base=base)
    layout = json.loads((tmp_path / "delivery/assembly_plan.json").read_text())
    assert layout["objects"]["I"]["caps"]["5p"]["record"] == "cap"
    parts = layout["objects"]["I"]["parts"]
    assert parts[1]["source"] == "hifi_gap_patch"
    assert parts[1]["strand"] == ("-" if reverse else "+")
    assert parts[2]["start"] == 86
    assert layout["base_fasta_sha256"] == sha256(base)


@pytest.mark.parametrize(("options", "reason"), [
    ({"protected": True}, "protected_or_inexact_gap"),
    ({"resolved": True}, "accepted_sequence_between_flanks_or_flanks_changed"),
    ({}, "approved_sequence_not_in_resolver_output"),
])
def test_conflicting_patches_keep_the_accepted_graph_assembly(tmp_path, options, reason):
    reference, base, plan, manifests, payloads, approvals, *_ = fixture(tmp_path, **options)
    if not options:
        (payloads / "gap1.fa").write_text(">gap1\nACGTACGTACGTAC\n")
    assert apply_approved_patches(reference, base, plan, manifests, payloads, approvals, tmp_path, tmp_path)[:2] == (base, plan)
    row = json.loads((tmp_path / "hifi_gap_patches.json").read_text())["patches"][0]
    assert (row["status"], row["reason"]) == ("withheld", reason)


def test_reviewed_nuclear_recipe_and_gap_patch_retain_other_reference_records(tmp_path):
    reference, reviewed, plan, manifests, payloads, approvals, recipe, sequence, patch = fixture(tmp_path)
    mt_sequence = "cggaatagtctacgacgagatccccgggatt"
    with reference.open("a") as handle:
        handle.write(f">scaffold_MT\n{mt_sequence}\n")
    document = json.loads(approvals.read_text())
    document["reference_sha256"] = sha256(reference)
    approvals.write_text(json.dumps(document))
    decision = tmp_path / "decision.json"
    decide(reference, [], [], reviewed, recipe, decision, nuclear_chromosomes={"I"})
    base, plan = retain_reference_records(reference, reviewed, plan, tmp_path)
    fasta, final_plan, results = apply_approved_patches(reference, base, plan, manifests, payloads, approvals,
                                                       tmp_path, tmp_path)
    assert results[0]["status"] == "applied"
    assert dict(read_fasta(fasta)) == {"I": sequence.replace("NNNN", patch), "scaffold_MT": mt_sequence}
    package(fasta, tmp_path / "delivery", {"I"}, "TTAGGC", 3, recipe, final_plan, reviewed_base=reviewed)
    layout = json.loads(final_plan.read_text())
    assert layout["objects"]["scaffold_MT"]["parts"][0]["source"] == "reference"
    assert layout["reviewed_fasta_sha256"] == sha256(reviewed)
    assert layout["base_fasta_sha256"] == sha256(base)
    delivery = tmp_path / "delivery"
    assert dict(read_fasta(delivery / "assembly.fasta")) == {"I": sequence.replace("NNNN", patch)}
    assert dict(read_fasta(delivery / "mitochondrial.fasta")) == {"scaffold_MT": mt_sequence}
    assert "scaffold_MT" not in (delivery / "assembly.agp").read_text()
    packaged = json.loads((delivery / "assembly_plan.json").read_text())
    assert set(packaged["objects"]) == {"I"}
    assert set(packaged["mitochondrial_objects"]) == {"scaffold_MT"}
    assert packaged["fasta_sha256"] == sha256(delivery / "assembly.fasta")
