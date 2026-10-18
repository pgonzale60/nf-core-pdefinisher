#!/usr/bin/env python3
"""pde_targeted_gap_assembly.py

Module 3 Engine A: Targeted Resolution of Assembly Gaps (Thin Orchestrator).
Orchestrates specialised sub-engines from `resources/usr/lib/gapassembly`:
  - `flank_io`: Flank extraction, palindromic anchor reduction, and BAM read baiting
  - `consensus_bridge`: Spanning-read abPOA consensus patching
  - `unitig_valency`: Genome-wide unitig valency conservation & PAF bridge evaluation

Authoritative Manifest Schema Version: 1.0.0
British spelling standard enforced across all documentation and comments.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

# Add resources/usr/lib to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

from gapassembly import (
    CONSENSUS_ACTIONS,
    bait_reads_for_gap,
    choose_anchor,
    collect_spanning_read_candidates,
    evaluate_gap_paf,
    extract_fasta_region,
    gfa_to_fasta,
    longest_inverted_repeat,
    parse_gaps_from_input,
    project_target_interval_via_cigar,
    read_fasta_records,
    resolve_spanned_gap,
    revcomp,
    reverse_complement,
    survey_unitig_valencies_genome_wide,
)

SCHEMA_VERSION = "1.0.0"


def get_script_fingerprint():
    try:
        hasher = hashlib.sha256()
        with open(__file__, "rb") as f:
            hasher.update(f.read())
        pkg_dir = Path(__file__).resolve().parent.parent / "lib" / "gapassembly"
        for py_file in sorted(pkg_dir.glob("*.py")):
            with open(py_file, "rb") as f:
                hasher.update(f.read())
        return hasher.hexdigest()[:16]
    except Exception:
        return "unknown"


def run_self_test():
    """Verify internal functions and dependencies."""
    sha = get_script_fingerprint()
    print(f"pde_targeted_gap_assembly self-test: schema={SCHEMA_VERSION} sha256={sha}")
    test_eval = evaluate_gap_paf("/nonexistent/path.paf", 10000, 200, 10000)
    assert test_eval["outcome"] == "UNRESOLVED_GAP", "Empty PAF test failed"

    qname = "read1"
    seq = "A" * 20000
    reads = {qname: (seq, "I" * 20000)}
    with tempfile.TemporaryDirectory() as tmp:
        paf = os.path.join(tmp, "reads_to_flanks.paf")

        def paf_line(qs, qe, strand, target, aln, mapq, tstart=None, tend=None, tlen=20000):
            if tstart is None or tend is None:
                if target == "LEFT":
                    tstart = tlen - aln
                    tend = tlen
                else:
                    tstart = 0
                    tend = aln
            return (
                "\t".join(
                    [
                        qname,
                        "20000",
                        str(qs),
                        str(qe),
                        strand,
                        target,
                        str(tlen),
                        str(tstart),
                        str(tend),
                        str(aln),
                        str(aln),
                        str(mapq),
                    ]
                )
                + "\n"
            )

        with open(paf, "w") as fh:
            fh.write(paf_line(12, 11018, "-", "LEFT", 11010, 60))
            fh.write(paf_line(6448, 7748, "+", "LEFT", 1300, 60))
            fh.write(paf_line(11018, 14518, "+", "RIGHT", 3500, 60))
        spanning = collect_spanning_read_candidates(paf, reads, min_mapq=30, min_anchor_bp=1000)
        assert spanning["n_spanning"] == 1
        cand = spanning["candidates"][0]
        assert (cand["strand"], cand["q_from"], cand["q_to"]) == ("+", 7748, 11018), cand

        # Test slack registration & internal flank rejection:
        # 1) LEFT alignment stopping 5,000 bp before boundary (tend=15000 < 20000-500) must be rejected
        with open(paf, "w") as fh:
            fh.write(paf_line(5000, 6500, "+", "LEFT", 1500, 60, tstart=13500, tend=15000))
            fh.write(paf_line(7000, 8500, "+", "RIGHT", 1500, 60, tstart=0, tend=1500))
        spanning_far = collect_spanning_read_candidates(paf, reads, min_mapq=30, min_anchor_bp=1000)
        assert spanning_far["n_spanning"] == 0, "Internal flank hit stopping 5kb before gap must be rejected"

        # 2) LEFT alignment stopping 40 bp short of boundary (tend=19960) and RIGHT starting at tstart=30
        # must register q_from = 6500 + 40 = 6540 and q_to = 7000 - 30 = 6970
        with open(paf, "w") as fh:
            fh.write(paf_line(5000, 6500, "+", "LEFT", 1500, 60, tstart=18460, tend=19960))
            fh.write(paf_line(7000, 7350, "+", "RIGHT", 350, 60, tstart=30, tend=380))
        spanning = collect_spanning_read_candidates(
            paf, reads, min_mapq=30, min_anchor_bp=1000, min_secondary_anchor_bp=300
        )
        assert spanning["n_spanning"] == 1
        assert (spanning["candidates"][0]["q_from"], spanning["candidates"][0]["q_to"]) == (6540, 6970)

    # Test CIGAR-exact query coordinate projection across flank indels
    aln_with_indel = {
        "qstart": 0,
        "qend": 41450,
        "strand": "+",
        "tstart": 0,
        "tend": 40200,
        "cigar": "10000=50D9950=200D1500I20000=",
    }
    qs_cig, qe_cig = project_target_interval_via_cigar(aln_with_indel, 20000, 20200)
    assert (qs_cig, qe_cig) == (19950, 21450), f"CIGAR projection mismatch: {(qs_cig, qe_cig)}"

    tail = "".join("ACGT"[(i * 3 + 1) % 4] for i in range(5000))
    window = reverse_complement(tail) + tail + tail
    anchor = choose_anchor(window, gap_at_end=True, min_anchor_bp=1000, max_inverted_repeat_bp=500)
    assert 1000 <= len(anchor) < len(window) and longest_inverted_repeat(anchor) <= 500

    with tempfile.TemporaryDirectory() as tmp:
        paf_ambig = os.path.join(tmp, "ambiguous_rep.paf")
        with open(paf_ambig, "w") as fh:
            fh.write("utg_ambig\t2000\t0\t2000\t+\tloc_ref\t40200\t19400\t20800\t1386\t1400\t0\n")
        res_ambig = evaluate_gap_paf(
            paf_ambig,
            left_flank_len=20000,
            gap_len=200,
            right_flank_len=20000,
            min_mapq=50,
            min_aln_len=1000,
            min_secondary_anchor=500,
            min_identity=0.98,
            repetitive_mapq_floor=0,
            require_unique_anchor=True,
            mapq_anchored=30,
            min_repetitive_anchor=1000,
        )
        assert res_ambig["outcome"] == "UNRESOLVED_GAP"

        res_mono = evaluate_gap_paf(
            paf_ambig,
            left_flank_len=20000,
            gap_len=200,
            right_flank_len=20000,
            min_mapq=50,
            min_aln_len=1000,
            min_secondary_anchor=500,
            min_identity=0.98,
            repetitive_mapq_floor=0,
            require_unique_anchor=True,
            mapq_anchored=30,
            min_repetitive_anchor=1000,
            monovalent_unitigs={"utg_ambig"},
        )
        assert res_mono["outcome"] == "PATCH_GAP_INSERTION"
        assert "monovalent anchor valency (P=1.0)" in res_mono["notes"]

        # Test dual monovalent repeat overhangs overlapping across a repeat gap
        paf_dual = os.path.join(tmp, "dual_overhang.paf")
        with open(paf_dual, "w") as fh:
            fh.write("utg_L\t2500\t0\t2000\t+\tloc_ref\t40200\t18000\t20000\t2000\t2000\t60\n")
            fh.write("utg_R\t2500\t500\t2500\t+\tloc_ref\t40200\t20200\t22200\t2000\t2000\t60\n")
        shared_ovl = "ACGTACGTACGTACGTACGT" * 10  # 200 bp overlap
        utg_seqs = {
            "utg_L": ("A" * 2000) + ("C" * 300) + shared_ovl,
            "utg_R": shared_ovl + ("G" * 300) + ("T" * 2000),
        }
        res_dual = evaluate_gap_paf(
            paf_dual,
            left_flank_len=20000,
            gap_len=200,
            right_flank_len=20000,
            min_mapq=50,
            min_aln_len=1000,
            min_secondary_anchor=500,
            min_identity=0.98,
            monovalent_unitigs={"utg_L", "utg_R"},
            allow_monovalent_overhang_patch=True,
            min_monovalent_overhang=200,
            unitig_seqs=utg_seqs,
        )
        assert res_dual["outcome"] == "PATCH_GAP_INSERTION"
        assert res_dual["patch_len"] == 800, f"Expected 300 + 200 + 300 = 800 bp merged patch, got {res_dual}"
        assert res_dual["patch_seq"] == ("C" * 300) + shared_ovl + ("G" * 300)

    print("Self-test PASS")
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(description="Targeted Local Gap Assembly using Hifiasm CCS-only.")
    parser.add_argument("--blueprint", help="Path to Module 2 blueprint JSON or curation cards")
    parser.add_argument("--gaps-bed", help="Optional path to gaps.bed")
    parser.add_argument("--bam", help="Path to PacBio HiFi sorted BAM")
    parser.add_argument("--ref", help="Path to reference FASTA")
    parser.add_argument("--sample", default="SAMPLE", help="Sample ID")
    parser.add_argument("--outdir", default=None, help="Directory to save gap patches and unitigs")
    parser.add_argument("--manifest", default=None, help="Output manifest TSV file")
    parser.add_argument("--policy", default=None, help="Path to curation_policy.json (SSOT)")
    parser.add_argument("--flank", type=int, default=20000, help="Flanking distance around gap in bp")
    parser.add_argument("--threads", type=int, default=4, help="Threads for hifiasm and minimap2")
    parser.add_argument("--self-test", action="store_true", help="Run self-test and exit")
    parser.add_argument("--version", action="store_true", help="Print version and exit")
    args = parser.parse_args()

    if args.version:
        print(f"pde_targeted_gap_assembly schema={SCHEMA_VERSION} sha256={get_script_fingerprint()}")
        sys.exit(0)

    if args.self_test:
        run_self_test()

    if not (args.bam and args.ref and (args.blueprint or args.gaps_bed)):
        parser.error("--bam, --ref, and either --blueprint or --gaps-bed are required.")

    sample_id = args.sample
    outdir = args.outdir or f"{sample_id}_gap_patches"
    manifest_path = args.manifest or f"{sample_id}_patches_manifest.tsv"
    os.makedirs(outdir, exist_ok=True)

    policy = {}
    if getattr(args, "policy", None) and os.path.exists(args.policy):
        with open(args.policy) as pf:
            policy = json.load(pf)
    elif args.blueprint and os.path.exists(args.blueprint):
        with open(args.blueprint) as f:
            bp_obj = json.load(f)
            if isinstance(bp_obj, dict):
                policy = bp_obj.get("policy", {})

    flank_len = policy.get("genomic_architecture", {}).get("gap_assembly_flank_bp", args.flank or 20000)
    min_mapq = policy.get("mapping_and_anchors", {}).get("mapq_high_confidence", 50)
    mapq_anchored = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30)
    min_aln_len = policy.get("mapping_and_anchors", {}).get("anchor_min_length_bp", 1000)
    min_baited_reads = policy.get("evidence_floors", {}).get("min_read_support_floor", 3)
    gp = policy.get("gap_patch", {})
    min_secondary_anchor_bp = gp.get("min_secondary_anchor_bp", 500)
    min_flank_identity = gp.get("min_flank_identity", 0.98)
    mapq_repetitive_gap_floor = gp.get("mapq_repetitive_gap_floor", 0)
    max_overlap_excision_bp = gp.get("max_overlap_excision_bp", 50000)
    require_unique_anchor = gp.get("require_unique_anchor_for_repetitive_gap", True)
    min_repetitive_anchor = gp.get("min_repetitive_anchor_bp", 1000)
    anchor_valency_survey_enabled = gp.get("anchor_valency_survey_enabled", True)
    max_offtarget_frac = gp.get("max_offtarget_repeat_fraction", 0.02)
    allow_monovalent_overhang = gp.get("allow_monovalent_overhang_patch", True)
    min_monovalent_overhang = gp.get("min_monovalent_overhang_bp", 200)
    minimap2_preset = gp.get("minimap2_preset", "asm20")

    gaps = parse_gaps_from_input(args.blueprint, args.gaps_bed)
    print(f"[{sample_id}] Identified {len(gaps)} gap targets for targeted assembly.", flush=True)

    manifest_rows = []
    header = [
        "gap_id",
        "chrom",
        "gap_start",
        "gap_end",
        "gap_len",
        "outcome",
        "left_anchor_len",
        "right_anchor_len",
        "patch_len",
        "excised_gap_len",
        "unitig_id",
        "patch_fasta",
        "notes",
    ]

    scratch_dir = os.path.join(os.path.dirname(outdir) or ".", f".{os.path.basename(outdir)}_scratch")
    os.makedirs(scratch_dir, exist_ok=True)

    try:
        for idx, g in enumerate(gaps, start=1):
            gid, chrom, g_start, g_end, g_len = g["gap_id"], g["chrom"], g["start"], g["end"], g["length"]
            gap_work = os.path.join(scratch_dir, gid)
            os.makedirs(gap_work, exist_ok=True)
            print(
                f"[{sample_id}] [{idx}/{len(gaps)}] Processing {gid} ({chrom}:{g_start:,}-{g_end:,}, {g_len} bp Ns) as {g['action']}...",
                flush=True,
            )

            cons_fallback_note = ""
            if g["action"] in CONSENSUS_ACTIONS:
                cons_row = resolve_spanned_gap(
                    g, gap_work, outdir, args.bam, args.ref, args.threads, policy, sample_id, idx, len(gaps)
                )
                if cons_row[5] == "PATCH_GAP_INSERTION":
                    manifest_rows.append(cons_row)
                    shutil.rmtree(gap_work, ignore_errors=True)
                    continue
                cons_fallback_note = f" [abpoa fallback: {cons_row[12]}]"
                print(
                    f"[{sample_id}] {gid}: Direct read consensus did not resolve gap ({cons_row[12]}); "
                    f"falling back to local hifiasm unitig assembly + valency evaluation...",
                    flush=True,
                )

            reads_fa = os.path.join(gap_work, "baited_reads.fa")
            n_reads = bait_reads_for_gap(args.bam, chrom, g_start, g_end, flank_len, reads_fa)
            if n_reads < min_baited_reads:
                manifest_rows.append(
                    [
                        gid,
                        chrom,
                        g_start,
                        g_end,
                        g_len,
                        "UNRESOLVED_GAP",
                        0,
                        0,
                        0,
                        0,
                        "NONE",
                        "NONE",
                        f"Insufficient baited reads ({n_reads} reads, required >= {min_baited_reads}){cons_fallback_note}",
                    ]
                )
                shutil.rmtree(gap_work, ignore_errors=True)
                continue

            hifi_prefix = os.path.join(gap_work, "hifi_asm")
            hifi_cmd = ["hifiasm", "-o", hifi_prefix, "-t", str(args.threads), "-f0", "-n", "1", reads_fa]
            try:
                subprocess.run(hifi_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            except Exception as e:
                manifest_rows.append(
                    [
                        gid,
                        chrom,
                        g_start,
                        g_end,
                        g_len,
                        "UNRESOLVED_GAP",
                        0,
                        0,
                        0,
                        0,
                        "NONE",
                        "NONE",
                        f"Hifiasm execution error: {e}{cons_fallback_note}",
                    ]
                )
                shutil.rmtree(gap_work, ignore_errors=True)
                continue

            unitigs_fa = os.path.join(gap_work, "unitigs.fa")
            if gfa_to_fasta(f"{hifi_prefix}.bp.p_ctg.gfa", unitigs_fa) == 0:
                manifest_rows.append(
                    [
                        gid,
                        chrom,
                        g_start,
                        g_end,
                        g_len,
                        "UNRESOLVED_GAP",
                        0,
                        0,
                        0,
                        0,
                        "NONE",
                        "NONE",
                        f"Hifiasm emitted zero contigs{cons_fallback_note}",
                    ]
                )
                shutil.rmtree(gap_work, ignore_errors=True)
                continue

        flank_ref_fa = os.path.join(gap_work, "flank_ref.fa")
        left_coord = max(1, g_start - flank_len + 1)
        right_coord = g_end + flank_len
        left_len = g_start - left_coord + 1
        right_len = right_coord - g_end
        extract_fasta_region(args.ref, chrom, left_coord, right_coord, flank_ref_fa)

        paf_out = os.path.join(gap_work, "unitigs_to_flanks.paf")
        mm_cmd = [
            "minimap2",
            "-x",
            minimap2_preset,
            "-c",
            "--secondary=yes",
            "-N",
            "50",
            "-p",
            "0.1",
            "-t",
            str(args.threads),
            flank_ref_fa,
            unitigs_fa,
        ]
        with open(paf_out, "w") as paf_f:
            subprocess.run(mm_cmd, stdout=paf_f, stderr=subprocess.DEVNULL)

        monovalent_unitigs = set()
        if anchor_valency_survey_enabled and os.path.exists(unitigs_fa):
            monovalent_unitigs = survey_unitig_valencies_genome_wide(
                args.ref,
                chrom,
                g_start,
                g_end,
                unitigs_fa,
                threads=args.threads,
                max_offtarget_fraction=max_offtarget_frac,
                min_anchor_bp=min_aln_len,
            )
            for u_id in monovalent_unitigs:
                print(
                    f"[{sample_id}] {gid}: Unitig {u_id} certified as MONOVALENT across genome (P=1.0)",
                    flush=True,
                )

        utg_records = read_fasta_records(unitigs_fa) if os.path.exists(unitigs_fa) else {}
        ref_gap_len = max(0, g_end - g_start) if g["action"] == "INSERT_SEQUENCE" else g_len
        res = evaluate_gap_paf(
            paf_out,
            left_len,
            ref_gap_len,
            right_len,
            g_start=g_start,
            g_end=g_end,
            min_mapq=min_mapq,
            min_aln_len=min_aln_len,
            min_secondary_anchor=min_secondary_anchor_bp,
            min_identity=min_flank_identity,
            repetitive_mapq_floor=mapq_repetitive_gap_floor,
            max_overlap_excision_bp=max_overlap_excision_bp,
            require_unique_anchor=require_unique_anchor,
            mapq_anchored=mapq_anchored,
            min_repetitive_anchor=min_repetitive_anchor,
            monovalent_unitigs=monovalent_unitigs,
            allow_monovalent_overhang_patch=allow_monovalent_overhang,
            min_monovalent_overhang=min_monovalent_overhang,
            unitig_seqs=utg_records,
        )
        patch_file = "NONE"
        if res["outcome"] == "PATCH_GAP_INSERTION" and os.path.exists(unitigs_fa):
            patch_file = os.path.join(outdir, f"{gid}_patch.fa")
            utg_seq = utg_records.get(res.get("unitig_id"), "")
            if res.get("patch_seq"):
                with open(patch_file, "w") as pf_out:
                    pf_out.write(f">{gid}\n{res['patch_seq']}\n")
            elif utg_seq and "q_patch_start" in res and "q_patch_end" in res:
                qs, qe = res["q_patch_start"], res["q_patch_end"]
                raw_patch = utg_seq[qs:qe] if qe > qs else ""
                patch_seq = revcomp(raw_patch) if res.get("strand") == "-" else raw_patch
                with open(patch_file, "w") as pf_out:
                    pf_out.write(f">{gid}\n{patch_seq}\n")
            elif res.get("unitig_id") in utg_records:
                with open(patch_file, "w") as pf_out:
                    pf_out.write(f">{gid}\n{utg_records[res['unitig_id']]}\n")
            else:
                shutil.copyfile(unitigs_fa, patch_file)

        manifest_rows.append(
            [
                gid,
                chrom,
                g_start,
                g_end,
                g_len,
                res["outcome"],
                res["left_anchor"],
                res["right_anchor"],
                res["patch_len"],
                res["excised_len"],
                res["unitig_id"] or "NONE",
                patch_file,
                f"{res['notes']}{cons_fallback_note}",
            ]
        )
        shutil.rmtree(gap_work, ignore_errors=True)
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)

    with open(manifest_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in manifest_rows:
            out.write("\t".join(map(str, r)) + "\n")

    print(f"[{sample_id}] Completed targeted gap assembly. Manifest written to: {manifest_path}", flush=True)


if __name__ == "__main__":
    main()
