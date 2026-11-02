#!/usr/bin/env python3
"""Verify inter-module JSON/TSV data contracts and AST/topology documentation against live code.

Validates:
  1. `assets/schemas/curation_blueprint.schema.json` and `assets/schemas/inter_module_tsv_contracts.json`
     are valid JSON and contain all required domains, step enums, and TSV contracts.
  2. Live output from `synthesize_telomere_cards` + `reconcile_unplaced_contigs_into_terminus_cards`
     matches the allowed `status` (`classification`), `recommended_action` (`action`), and
     `execution_plan[].step` enums in `assets/schemas/curation_blueprint.schema.json`.
  3. TSV headers and status gates in `pde_terminal_softclips.py`, `pde_telomere_walker.py`,
     `pde_targeted_gap_assembly.py`, `pde_rdna_cassette_resolver.py`,
     `pde_transchromosomal_link_resolver.py`, and `pde_generate_agp.py` match
     `assets/schemas/inter_module_tsv_contracts.json` (specifically asserting
     `EXTENDED_SUBTELOMERE_ONLY` is in `non_spliceable_statuses` and rejected in `pde_generate_agp.py`).
  4. `scripts/extract_code_contracts.py --check` passes (zero drift between code/schemas and
     `docs/CODE_CONTRACTS.md`).
"""

import json
import pathlib
import subprocess
import sys
import tempfile


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
BLUEPRINT_SCHEMA_PATH = REPO_ROOT / "assets" / "schemas" / "curation_blueprint.schema.json"
TSV_CONTRACTS_PATH = REPO_ROOT / "assets" / "schemas" / "inter_module_tsv_contracts.json"

# Add curationsynthesis lib to sys.path for live synthesis verification
sys.path.insert(
    0,
    str(REPO_ROOT / "modules" / "local" / "curationsynthesis" / "resources" / "usr" / "lib"),
)

from curationsynthesis import (  # noqa: E402
    DEFAULT_POLICY,
    build_manifest_json,
    load_terminal_softclips,
    load_terminal_summary,
    reconcile_unplaced_contigs_into_terminus_cards,
    synthesize_telomere_cards,
)


def verify_schemas_structure() -> tuple[dict, dict]:
    assert BLUEPRINT_SCHEMA_PATH.is_file(), f"Missing {BLUEPRINT_SCHEMA_PATH}"
    assert TSV_CONTRACTS_PATH.is_file(), f"Missing {TSV_CONTRACTS_PATH}"

    bp_schema = json.loads(BLUEPRINT_SCHEMA_PATH.read_text(encoding="utf-8"))
    tsv_contracts = json.loads(TSV_CONTRACTS_PATH.read_text(encoding="utf-8"))

    assert bp_schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema"
    props = bp_schema.get("properties", {})
    for key in (
        "sample_id",
        "provenance",
        "summary",
        "internal_gaps",
        "chromosome_termini",
        "trans_chromosomal_links",
        "unplaced_scaffolds",
        "action_cards",
    ):
        assert key in props, f"Missing top-level property '{key}' in curation_blueprint.schema.json"

    defs = bp_schema.get("$defs", {})
    step_enum = set(defs.get("terminusExecutionStepName", {}).get("enum", []))
    required_steps = {
        "DEFER_TO_RDNA_RESOLVER",
        "PRUNE_CHIMERIC_TERMINUS",
        "ALREADY_CAPPED",
        "CAP_WITH_UNITIG",
        "INTEGRATE_UNPLACED_CONTIG",
        "SCAFFOLD_REPEAT_TO_TELOMERE",
        "ALIGN_SEED_TO_UNITIGS",
        "EXTEND_TELOMERE_TRACT",
        "CAP_WITH_READ",
        "ITERATIVE_READ_WALKER",
        "RETAIN_AS_IS",
    }
    assert step_enum == required_steps, f"Step enum mismatch: {step_enum ^ required_steps}"

    contracts = tsv_contracts.get("contracts", {})
    required_tsvs = {
        "terminal_softclips.tsv",
        "terminal_softclips_summary.tsv",
        "telomere_extension_manifest.tsv",
        "targeted_gap_manifest.tsv",
        "rdna_cassette_manifest.tsv",
        "trans_links_manifest.tsv",
        "curation_audit.tsv",
    }
    missing = required_tsvs - set(contracts.keys())
    assert not missing, f"Missing required TSV contracts: {missing}"

    return bp_schema, tsv_contracts


def verify_live_terminus_cards_against_schema(bp_schema: dict) -> None:
    defs = bp_schema["$defs"]
    allowed_steps = set(defs["terminusExecutionStepName"]["enum"])
    allowed_statuses = set(defs["terminusStatusEnum"]["enum"])
    allowed_actions = set(defs["terminusActionEnum"]["enum"])

    nuclear_chroms = ["I", "II", "III", "IV", "V", "X"]
    chrom_lens = {c: 1_000_000 for c in nuclear_chroms}
    grs_segments = {
        "I": [(0, 50_000, 50_000, 100.0)],
        "II": [],
        "III": [(0, 40_000, 40_000, 100.0), (960_000, 1_000_000, 40_000, 100.0)],
        "IV": [],
        "V": [(950_000, 1_000_000, 50_000, 100.0)],
        "X": [],
    }
    core_segments = {
        c: [(50_000, 950_000, 900_000, 30.0)] for c in nuclear_chroms
    }
    delfies_points = {"II": [50_000], "IV": [950_000]}
    rdna_intervals = {"I": [(985_000, 1_000_000)]}
    baseline_status = {
        ("I", "5p"): "ALREADY_CAPPED",
        ("I", "3p"): "UNCAPPED_EXTENDABLE",
        ("II", "5p"): "NO_CLIPS",
        ("II", "3p"): "TELOMERE_IN_OVERHANG",
        ("III", "5p"): "EXTEND_AND_CAP",
        ("III", "3p"): "ALREADY_CAPPED",
        ("IV", "5p"): "EXTEND_AND_CAP",
        ("IV", "3p"): "NO_CLIPS",
        ("V", "5p"): "NO_CLIPS",
        ("V", "3p"): "UNCAPPED_EXTENDABLE",
        ("X", "5p"): "NO_CLIPS",
        ("X", "3p"): "NO_CLIPS",
    }
    tip_tract_bp = {
        ("I", "5p"): 300,
        ("II", "3p"): 0,
        ("III", "3p"): 450,
    }
    unitig_extensions = {
        ("III", "5p"): [(1500, "utg0001l", 600, 900, 0, 1500, "+", 10_000, 0, 8500)],
        ("III", "3p"): [(800, "utg0002l", 800, 0, 4200, 5000, "+", 5000, 995_000, 1_000_000)],
    }
    single_read_ext = {
        ("II", "3p"): (600, 4, "read_overhang_1", 25, ["read_overhang_1", "r2", "r3", "r4"]),
        ("IV", "5p"): (900, 5, "read_cap_1", 35, ["read_cap_1", "r2", "r3", "r4", "r5"]),
    }
    terminal_matches = {
        ("V", "5p"): {
            "cluster": {
                "cluster_id": "cl_0",
                "classification": "GERMLINE_JUNCTION",
                "subtype": "SATELLITE_JUNCTION",
                "read_count": 6,
                "mean_telomere_tract_bp": 750,
                "dominant_motif": "GCCTAA",
                "supporting_read_ids": ["r_sat_1", "r_sat_2", "r_sat_3"],
            },
            "scaffolds": ["scaffold_10"],
            "orientations": {"scaffold_10": "+"},
        }
    }
    # Chimeric terminus candidate at X 3p
    raw_softclips = {
        ("X", "3p"): [
            {
                "qname": f"chim_read_{i}",
                "align_start": 980_000,
                "align_end": 990_000,
                "clip_len": 1200,
                "overhang_bp": 0,
                "telo_cnt": 30,
                "mapq": 60,
                "collinear_spanning_reads": 0,
            }
            for i in range(4)
        ]
    }

    cards = synthesize_telomere_cards(
        nuclear_chroms=nuclear_chroms,
        chrom_lens=chrom_lens,
        grs_segments=grs_segments,
        core_segments=core_segments,
        delfies_points=delfies_points,
        rdna_intervals=rdna_intervals,
        baseline_status=baseline_status,
        unitig_extensions=unitig_extensions,
        single_read_ext=single_read_ext,
        tip_tract_bp=tip_tract_bp,
        policy=DEFAULT_POLICY,
        terminal_matches=terminal_matches,
        raw_softclips=raw_softclips,
    )

    scaffold_dispositions = [
        {
            "scaffold": "scaffold_22",
            "length": 15_000,
            "chromosome": "V",
            "terminus": "3p",
            "action": "INTEGRATE_INTO_CHROMOSOME",
            "strand": "+",
            "ext_qstart": 3000,
            "ext_qend": 15_000,
        }
    ]
    reconcile_unplaced_contigs_into_terminus_cards(
        cards,
        scaffold_dispositions=scaffold_dispositions,
        min_overhang_bp=200,
    )

    assert len(cards) == 12, f"Expected 12 terminus cards, got {len(cards)}"
    seen_steps = set()
    for card in cards:
        cls = card.get("classification")
        act = card.get("action")
        assert cls in allowed_statuses, f"Card {card['card_id']} classification '{cls}' not in schema {allowed_statuses}"
        assert act in allowed_actions, f"Card {card['card_id']} action '{act}' not in schema {allowed_actions}"
        assert isinstance(card.get("evidence"), dict), f"Card {card['card_id']} missing evidence dict"
        plan = card.get("execution_plan")
        assert isinstance(plan, list) and len(plan) >= 1, f"Card {card['card_id']} has empty execution_plan"
        for step_obj in plan:
            step_name = step_obj.get("step")
            assert step_name in allowed_steps, f"Step '{step_name}' not in schema {allowed_steps}"
            seen_steps.add(step_name)

    # Confirm all 11 step types were exercised across the synthetic karyotype
    assert seen_steps == allowed_steps, f"Synthetic karyotype did not exercise all steps: missing {allowed_steps - seen_steps}"

    # Verify manifest top-level keys match schema required keys
    manifest = build_manifest_json(
        sample_id="test_sample",
        script_fingerprint="0123456789abcdef",
        nuclear_chroms=nuclear_chroms,
        unplaced_chroms=["scaffold_10", "scaffold_22"],
        grs_threshold=50.0,
        karyotype_stats={
            "tot_core_len": 5_400_000,
            "core_pct": 90.0,
            "cohort_core_d": 30.0,
            "tot_grs_len": 600_000,
            "grs_pct": 10.0,
            "cohort_grs_d": 100.0,
        },
        policy=DEFAULT_POLICY,
        internal_cards=[],
        telo_cards=cards,
        rdna_cards=[],
        scaffold_cards=[],
        scaffold_dispositions=scaffold_dispositions,
    )
    for req_key in bp_schema["required"]:
        assert req_key in manifest, f"Required manifest key '{req_key}' missing from build_manifest_json output"


def verify_tsv_headers_and_status_gates(tsv_contracts: dict) -> None:
    contracts = tsv_contracts["contracts"]

    # 1. SOFTCLIPENDS: terminal_softclips.tsv (10-col & 9-col legacy) and terminal_softclips_summary.tsv
    softclips_src = (
        REPO_ROOT / "modules" / "local" / "softclipends" / "resources" / "usr" / "bin" / "pde_terminal_softclips.py"
    ).read_text(encoding="utf-8")

    expected_sc_cols = contracts["terminal_softclips.tsv"]["columns"]
    assert "\\t".join(expected_sc_cols[:6]) in softclips_src
    assert "\\t".join(expected_sc_cols[6:]) in softclips_src

    expected_sum_cols = contracts["terminal_softclips_summary.tsv"]["columns"]
    assert "\\t".join(expected_sum_cols[:5]) in softclips_src
    assert "\\t".join(expected_sum_cols[5:]) in softclips_src

    # Verify consumer load_terminal_softclips handles both 9-col and 10-col TSVs and load_terminal_summary handles 7-col
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)
        sc_10 = tmp / "sc10.tsv"
        sc_10.write_text(
            "\t".join(expected_sc_cols) + "\n"
            "r1\tI\t5p\t0\t5000\t600\t600\t25\t60\t2\n",
            encoding="utf-8",
        )
        raw_10, _ = load_terminal_softclips(str(sc_10))
        assert raw_10[("I", "5p")][0]["collinear_spanning_reads"] == 2

        sc_9 = tmp / "sc9.tsv"
        sc_9.write_text(
            "\t".join(contracts["terminal_softclips.tsv"]["legacy_9col_columns"]) + "\n"
            "r1\tI\t5p\t0\t5000\t600\t600\t25\t60\n",
            encoding="utf-8",
        )
        raw_9, _ = load_terminal_softclips(str(sc_9))
        assert raw_9[("I", "5p")][0]["collinear_spanning_reads"] == 0

        sum_tsv = tmp / "sum.tsv"
        sum_tsv.write_text(
            "\t".join(expected_sum_cols) + "\n"
            "I_5p\t100000\t5\t3\t30\tALREADY_CAPPED\t180\n",
            encoding="utf-8",
        )
        b_stat, t_tract = load_terminal_summary(str(sum_tsv))
        assert b_stat[("I", "5p")] == "ALREADY_CAPPED"
        assert t_tract[("I", "5p")] == 180

    # 2. TELOMEREWALKER: telomere_extension_manifest.tsv
    walker_src = (
        REPO_ROOT / "modules" / "local" / "telomerewalker" / "resources" / "usr" / "bin" / "pde_telomere_walker.py"
    ).read_text(encoding="utf-8")
    walker_spec = contracts["telomere_extension_manifest.tsv"]
    for col in walker_spec["columns"]:
        assert f'"{col}"' in walker_src, f"Column {col} missing in pde_telomere_walker.py"
    assert "EXTENDED_SUBTELOMERE_ONLY" in walker_spec["non_spliceable_statuses"]
    assert "EXTENDED_SUBTELOMERE_ONLY" not in walker_spec["spliceable_statuses"]

    # 3. TARGETEDGAPASSEMBLY: targeted_gap_manifest.tsv
    gap_src = (
        REPO_ROOT
        / "modules"
        / "local"
        / "targetedgapassembly"
        / "resources"
        / "usr"
        / "bin"
        / "pde_targeted_gap_assembly.py"
    ).read_text(encoding="utf-8")
    for col in contracts["targeted_gap_manifest.tsv"]["columns"]:
        assert f'"{col}"' in gap_src, f"Column {col} missing in pde_targeted_gap_assembly.py"

    # 4. RDNACASSETTERESOLVER: rdna_cassette_manifest.tsv
    rdna_src = (
        REPO_ROOT
        / "modules"
        / "local"
        / "rdnacassetteresolver"
        / "resources"
        / "usr"
        / "bin"
        / "pde_rdna_cassette_resolver.py"
    ).read_text(encoding="utf-8")
    for col in contracts["rdna_cassette_manifest.tsv"]["columns"]:
        assert f'"{col}"' in rdna_src, f"Column {col} missing in pde_rdna_cassette_resolver.py"

    # 5. TRANSCHROMOSOMALLINKRESOLVER: trans_links_manifest.tsv
    trans_src = (
        REPO_ROOT
        / "modules"
        / "local"
        / "transchromosomallinkresolver"
        / "resources"
        / "usr"
        / "bin"
        / "pde_transchromosomal_link_resolver.py"
    ).read_text(encoding="utf-8")
    for col in contracts["trans_links_manifest.tsv"]["columns"]:
        assert f'"{col}"' in trans_src, f"Column {col} missing in pde_transchromosomal_link_resolver.py"

    # 6. CURATIONAGPGENERATOR: curation_audit.tsv and EXTENDED_SUBTELOMERE_ONLY rejection gate
    agp_src = (
        REPO_ROOT
        / "modules"
        / "local"
        / "curationagpgenerator"
        / "resources"
        / "usr"
        / "bin"
        / "pde_generate_agp.py"
    ).read_text(encoding="utf-8")
    for col in contracts["curation_audit.tsv"]["columns"]:
        assert f'"{col}"' in agp_src, f"Column {col} missing in pde_generate_agp.py"
    assert 'status != "EXTENDED_SUBTELOMERE_ONLY"' in agp_src, (
        "pde_generate_agp.py must explicitly reject EXTENDED_SUBTELOMERE_ONLY from splicing"
    )
    assert "REJECT_SUBTELOMERE_ONLY" in agp_src


def verify_code_contracts_doc_zero_drift() -> None:
    res = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "extract_code_contracts.py"), "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, (
        f"extract_code_contracts.py --check failed:\nstdout: {res.stdout}\nstderr: {res.stderr}"
    )


def main() -> int:
    bp_schema, tsv_contracts = verify_schemas_structure()
    verify_live_terminus_cards_against_schema(bp_schema)
    verify_tsv_headers_and_status_gates(tsv_contracts)
    verify_code_contracts_doc_zero_drift()
    print("check_data_contracts: PASS (all JSON/TSV schemas, live synthesis enums, and docs/CODE_CONTRACTS.md verified)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
