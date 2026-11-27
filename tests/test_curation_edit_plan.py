import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from curation_edit_plan import catalog_candidates, decide
from package_reviewed_assembly import sha256


def test_graph_and_bam_proposals_share_a_curated_coordinate_plan(tmp_path):
    reference = tmp_path / "curated.fa"
    reference.write_text(">I\nAACCGGTT\n>II\nAACCGGTT\n")
    blueprint = tmp_path / "test.curation_manifest.json"
    blueprint.write_text(json.dumps({"action_cards": [
        {"card_id": "card-1", "action": "EXCISE_GAP", "chrom": "I", "start": 2, "end": 4},
    ]}))
    transfers = tmp_path / "test_trans_manifest.tsv"
    transfers.write_text(
        "card_id\toutcome\tsource_chrom\tsource_gap_start\tsource_gap_end\tcontig_start\tcontig_end\t"
        "target_chrom\ttarget_coord\ttarget_feature\n"
        "t1\tTRANSFER_AND_STITCH_CONTIG\tII\t2\t3\t3\t6\tI\t0\t5_prime_terminus\n"
    )
    candidates = [{"end_id": "S:5p", "curated_end_candidates": "I:5p:+:0:2:50;II:5p:+:0:2:0",
                   "category": "shared", "horizon_port": "u+", "candidate_endings": "v-;w+",
                   "family_unitig_matches": "u:C000001"}]
    plan = tmp_path / "plan.json"
    assert decide(reference, candidates, [blueprint, transfers], None, None, plan) == reference
    discovery = json.loads(plan.read_text())
    assert discovery["delivery_sha256"] == sha256(reference)
    assert {p["source"] for p in discovery["proposals"]} == {"graph", "bam"}
    assert all(p["decision"] == "review_required" for p in discovery["proposals"])
    assert any(p["reason"] == "competing_curated_end_placements" for p in discovery["proposals"])

    reviewed = tmp_path / "reviewed.fa"
    reviewed.write_text(">I\nTTAACCGGTT\n>II\nAACCGGTT\n")
    recipe = tmp_path / "recipe.json"
    recipe.write_text(json.dumps({"expected_sha256": sha256(reviewed), "reference_sha256": sha256(reference)}))
    assert decide(reference, candidates, [blueprint, transfers], reviewed, recipe, plan) == reviewed
    replay = json.loads(plan.read_text())
    assert replay["accepted_proposal"] == "graph:reviewed_replay"
    transfer = next(p for p in replay["proposals"] if p["kind"] == "TRANSFER_AND_STITCH_CONTIG")
    assert transfer["targets"] == [{"record": "II", "start": 2, "end": 6},
                                   {"record": "I", "side": "5p", "start": 0, "end": 0}]
    assert transfer["reason"] == "within_reviewed_replay_scope"

    reviewed.write_text(">I\nTTAACCGGTT\n")
    recipe.write_text(json.dumps({"expected_sha256": sha256(reviewed)}))
    with pytest.raises(ValueError, match="omits curated records: II"):
        decide(reference, candidates, [blueprint], reviewed, recipe, plan)


def test_native_catalog_evidence_keeps_candidate_order_and_empty_positions():
    row = {
        "end_id": "S:5p",
        "candidate_endings": "A+[reads:1,family:C1];B-[graph];C+[family:C2]",
        "candidate_families": "C1;;C2",
        "candidate_molecules": "r1,r2;;r2",
    }
    assert catalog_candidates(row) == [
        {"port": "A+", "family_ids": ["C1"], "reported_molecule_ids": ["r1", "r2"],
         "evidence_tags": ["reads:1", "family:C1"]},
        {"port": "B-", "family_ids": [], "reported_molecule_ids": [], "evidence_tags": ["graph"]},
        {"port": "C+", "family_ids": ["C2"], "reported_molecule_ids": ["r2"], "evidence_tags": ["family:C2"]},
    ]
    with pytest.raises(ValueError, match="candidate_families does not match"):
        catalog_candidates({**row, "candidate_families": "C1;C2"})
