import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import graph_finish
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
    monkeypatch.setattr(graph_finish, "source_records", lambda path, _: dict(graph_finish.read_fasta(path)))
    plan = tmp_path / "plan.json"
    graph_finish.trace_recipe(recipe, reviewed, {"primary": reference}, {"primary": reference}, plan)
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
        graph_finish.trace_recipe(recipe, reviewed, {"primary": reference}, {"primary": reference}, plan)
