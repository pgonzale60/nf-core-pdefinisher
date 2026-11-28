import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import gfapdefinisher_handoff
from package_reviewed_assembly import package, sha256


def test_reviewed_path_round_trips_with_reverse_source_and_unknown_join(tmp_path, monkeypatch):
    reference = tmp_path / "sources.fa"
    reference.write_text(">forward\nACCGNN\n>reverse\nTTGA\n")
    reviewed = tmp_path / "reviewed.fa"
    reviewed.write_text(">I\nACCGNNNCAA\n")
    recipe = tmp_path / "recipe.json"
    recipe.write_text(json.dumps({
        "gap_bp": 3,
        "expected_sha256": sha256(reviewed),
        "chromosomes": [{
            "name": "I", "cap_5p": {"record": "forward"}, "cap_3p": {"record": "reverse"},
            "parts": [
                {"source": "primary", "record": "forward"},
                {"gap": "satellite_port"},
                {"source": "primary", "record": "reverse", "strand": "-", "start": 1},
            ],
        }],
    }))
    monkeypatch.setattr(gfapdefinisher_handoff, "source_records", lambda path, _: dict(gfapdefinisher_handoff.read_fasta(path)))
    plan = tmp_path / "plan.json"
    gfapdefinisher_handoff.trace_recipe(recipe, reviewed, {"primary": reference}, {"primary": reference}, plan)
    parts = json.loads(plan.read_text())["objects"]["I"]["parts"]
    assert parts[1] == {"start": 4, "end": 7, "gap": "source_gap"}
    assert parts[2]["strand"] == "-"
    assert (parts[2]["source_start"], parts[2]["source_end"]) == (1, 4)
    delivery = tmp_path / "delivery"
    package(reviewed, delivery, {"I"}, "TTAGGC", 3, recipe, plan)
    assert reviewed.read_bytes() == (delivery / "assembly.fasta").read_bytes()
    assert "I\t5\t7\t2\tN\t3\tcontig\tno\tna" in (delivery / "assembly.agp").read_text()

    reviewed.write_text(">I\nACCGNNNCCC\n")
    with pytest.raises(ValueError, match="fail to reconstruct"):
        gfapdefinisher_handoff.trace_recipe(recipe, reviewed, {"primary": reference}, {"primary": reference}, plan)


def test_placed_rdna_with_unresolved_copy_order_is_not_an_unresolved_join(tmp_path, monkeypatch):
    reference = tmp_path / "sources.fa"
    reference.write_text(">scaffold\nAAANNGGGG\n>satellite\nACGT\n")
    graph_agp = tmp_path / "graph.agp"
    graph_agp.write_text("scaffold\t4\t5\t2\tU\t2\trepeat\tyes\trna_tangle\n")
    reviewed = tmp_path / "reviewed.fa"
    reviewed.write_text(">II\nCCCCNNNTTT\n>III\nACGTNNNACGT\n")
    recipe = tmp_path / "recipe.json"
    satellite = {"source": "finish", "record": "satellite"}
    recipe.write_text(json.dumps({
        "gap_bp": 3,
        "expected_sha256": sha256(reviewed),
        "chromosomes": [
            {"name": "II", "cap_5p": {}, "cap_3p": {},
             "parts": [{"source": "finish", "record": "scaffold", "strand": "-"}]},
            {"name": "III", "cap_5p": {}, "cap_3p": {},
             "parts": [satellite, {"gap": "satellite_port"}, satellite]},
        ],
    }))
    monkeypatch.setattr(gfapdefinisher_handoff, "source_records", lambda path, _: dict(gfapdefinisher_handoff.read_fasta(path)))
    plan = tmp_path / "plan.json"
    gfapdefinisher_handoff.trace_recipe(recipe, reviewed, {"finish": reference}, {"finish": reference}, plan, graph_agp)
    delivery = tmp_path / "delivery"
    package(reviewed, delivery, {"II", "III"}, "TTAGGC", 3, recipe, plan)
    gaps = list(gfapdefinisher_handoff.rows(delivery / "sequence_gaps.tsv"))
    assert [(r["object"], r["linkage_status"], r["sequence_status"]) for r in gaps] == [
        ("II", "flanking_connection_supported", "sequence_and_copy_order_unresolved"),
        ("III", "distal_placement_unresolved", "sequence_unresolved"),
    ]
    assert list(gfapdefinisher_handoff.rows(delivery / "unresolved_joins.tsv")) == gaps[1:]
    assert "II\t5\t7\t2\tN\t3\trepeat\tyes\tunspecified" in (delivery / "assembly.agp").read_text()


def test_reference_alignment_is_linked_to_its_proximal_scaffold_end(tmp_path):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "end_id\tobject_id\tside\thorizon_port\tcandidate_endings\n"
        "S:5p\tS\t5p\tu+\t\n"
        "S:3p\tS\t3p\tv-\t\n"
    )
    paf = tmp_path / "scaffolds.paf"
    paf.write_text("S\t1000\t10\t190\t+\tI\t180\t0\t180\t180\t180\t60\n")
    config = {"mapping_min_bp": 100, "mapping_identity": 0.98, "reference_end_window_bp": 20}
    outfile = tmp_path / "proposals.tsv"
    matches = gfapdefinisher_handoff.report_candidates(catalog, [], paf, config, outfile)
    assert matches["I", "5p"] == ["S:5p"]
    assert matches["I", "3p"] == []
    assert [row["decision"] for row in gfapdefinisher_handoff.rows(outfile)] == [
        "review_required", "unplaced_graph_end",
    ]


def test_large_end_catalog_field_survives_proposal_round_trip(tmp_path):
    catalog = tmp_path / "end_catalog.tsv"
    candidate_endings = ",".join(f"utg{i:06d}+" for i in range(15_000))
    catalog.write_text(
        "end_id\tobject_id\tside\thorizon_port\tcandidate_endings\n"
        f"S:5p\tS\t5p\tu+\t{candidate_endings}\n"
    )
    paf = tmp_path / "scaffolds.paf"
    paf.touch()
    config = {"mapping_min_bp": 100, "mapping_identity": 0.98, "reference_end_window_bp": 20}
    placements = [{"source": "primary", "record": "utg014999", "cluster_id": "family"}]
    proposals = tmp_path / "graph_end_proposals.tsv"

    gfapdefinisher_handoff.report_candidates(catalog, placements, paf, config, proposals)

    [proposal] = gfapdefinisher_handoff.rows(proposals)
    assert proposal["candidate_endings"] == candidate_endings
    assert proposal["family_unitig_matches"] == "utg014999:family"
