#!/usr/bin/env python3
"""pde_telomere_walker.py

Module 3 Engine B: Modular Two-Phase Telomere Extension, Capping, Scaffolding & Walker Engine.
Orchestrates four specialised sub-engines from resources/usr/lib:
  - telotract: Repeat tract quantification and terminal boundary analysis
  - telocapping: Unitig interval materialisation and single-read overhang capping
  - telowalking: Prefix-suffix graph traversal with somatic exclusion enforcement
  - teloscaffold: Nominal 100 bp N-gap scaffolding for unbridgeable satellite repeat arrays

Authoritative Manifest Schema Version: 1.0.0
British spelling standard enforced across all documentation and comments.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# Add resources/usr/lib to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

from telocapping import (
    authorized_unitig_materialisation,
    extract_read_by_qname,
    inspect_terminal_reads,
    maximal_tip_overhang,
    named_read_overhang,
    terminal_read_candidates,
)
from teloscaffold import find_sequence_overlap, materialise_satellite_scaffold
from telotract import (
    cap_status_for,
    count_telomeric_bp,
    maximal_terminal_tract,
    revcomp,
    terminal_tract_region,
)
from telowalking import (
    UnmappedReadIndex,
    evaluate_step_information_gain,
    run_iterative_walker,
)

SCHEMA_VERSION = "1.0.0"


def safe_int(val, default=0):
    """Safely convert value to integer."""
    if val is None:
        return default
    if isinstance(val, (int, float)):
        return int(val)
    cleaned = re.sub(r"[^\d\-]", "", str(val))
    return int(cleaned) if cleaned else default


def get_script_fingerprint():
    """Compute sha256 checksum fingerprint of the driver script and lib sub-packages."""
    try:
        hasher = hashlib.sha256()
        with open(__file__, "rb") as f:
            hasher.update(f.read())
        lib_root = Path(__file__).resolve().parent.parent / "lib"
        if lib_root.is_dir():
            for py_file in sorted(lib_root.rglob("*.py")):
                with open(py_file, "rb") as f:
                    hasher.update(f.read())
        return hasher.hexdigest()[:16]
    except Exception:
        return "unknown"


def parse_chrom_lens_from_bam(bam_path):
    """Parse chromosome names and lengths from BAM header @SQ tags."""
    cmd = ["samtools", "view", "-H", bam_path]
    res = subprocess.run(cmd, capture_output=True, text=True)
    chrom_lens = {}
    for line in res.stdout.splitlines():
        if line.startswith("@SQ\t"):
            parts = line.split("\t")
            sn = None
            ln = None
            for p in parts[1:]:
                if p.startswith("SN:"):
                    sn = p[3:]
                elif p.startswith("LN:"):
                    ln = int(p[3:])
            if sn and ln:
                chrom_lens[sn] = ln
    return chrom_lens


def parse_telomere_profile(tsv_path, min_repeats=15, min_mapq=30):
    """Parse PROFILE_TELOMERES TSV:

    unitig_id, unitig_len, chrom, chrom_len, aln_strand, q_aln_start, q_aln_end,
    t_aln_start, t_aln_end, mapq, telo_start, telo_end, telo_strand, telo_repeats, telo_tract_bp, subtelo_gap_bp
    """
    profile_hits = defaultdict(list)
    if not tsv_path or not os.path.exists(tsv_path):
        return profile_hits

    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
    expected = ["unitig_id", "unitig_len", "chrom", "chrom_len", "aln_strand", "q_aln_start"]
    if header[:2] == ["chrom_end", "chrom_length"]:
        print(f"[telomerewalker] {tsv_path} is the softclipends terminal summary, not a PROFILE_TELOMERES "
              "table, so it carries no unitig alignment coordinates. Unitig extensions are taken "
              "from the blueprint's unitig_extension_region.", flush=True)
        return profile_hits
    if header[:len(expected)] != expected:
        raise SystemExit(
            f"--telomere-profile {tsv_path} is not a PROFILE_TELOMERES table (header: {header[:6]}; expected the "
            f"first six columns to be {expected})")

    with open(tsv_path) as f:
        header = None
        for line in f:
            parts = line.strip().split("\t")
            if not header:
                header = parts
                continue
            if len(parts) < 14:
                continue
            chrom = parts[2]
            t_start = int(parts[7])
            t_end = int(parts[8])
            mapq = int(parts[9])
            telo_strand = parts[12]
            telo_repeats = int(parts[13])
            telo_tract_bp = int(parts[14]) if len(parts) > 14 else 0
            subtelo_gap = int(parts[15]) if len(parts) > 15 else 0

            if telo_strand != "NONE" and telo_repeats >= min_repeats and mapq >= min_mapq:
                profile_hits[chrom].append({
                    "unitig_id": parts[0],
                    "unitig_len": int(parts[1]) if len(parts) > 1 else 0,
                    "t_start": t_start,
                    "t_end": t_end,
                    "q_aln_start": int(parts[5]) if len(parts) > 5 else 0,
                    "q_aln_end": int(parts[6]) if len(parts) > 6 else 0,
                    "aln_strand": parts[4] if len(parts) > 4 else "+",
                    "mapq": mapq,
                    "telo_strand": telo_strand,
                    "telo_repeats": telo_repeats,
                    "telo_tract_bp": telo_tract_bp,
                    "subtelo_gap": subtelo_gap
                })
    return profile_hits


def write_extension_fasta(outdir, term_id, suffix, note, seq, motif, rev_motif):
    """Write an extension and return (path, added_bp, telo_bp, non_telo_bp)."""
    path = os.path.join(outdir, f"{term_id}_{suffix}.fa")
    with open(path, "w") as f:
        f.write(f">{term_id} {note}\n{seq}\n")
    telo = count_telomeric_bp(seq, motif, rev_motif)
    return path, len(seq), telo, len(seq) - telo


_INDEXABLE_REF_CACHE = {}


def require_indexable_reference(ref_fasta):
    """Return the reference, having proved that samtools can resolve regions from it."""
    if not ref_fasta or not os.path.exists(ref_fasta):
        raise RuntimeError(f"reference {ref_fasta} does not exist")
    if ref_fasta in _INDEXABLE_REF_CACHE:
        return _INDEXABLE_REF_CACHE[ref_fasta]
    try:
        res = subprocess.run(["samtools", "faidx", ref_fasta], capture_output=True, text=True)
    except OSError as e:
        raise RuntimeError(f"cannot run samtools to index {ref_fasta}: {e}")
    if res.returncode != 0 or not os.path.exists(ref_fasta + ".fai"):
        raise RuntimeError(
            f"reference {ref_fasta} is not indexable, so no region can be read from it. This engine requires "
            "the bgzipped assembly produced by BGZIPFASTA, with its .gzi staged alongside it. "
            f"samtools said: {(res.stderr or '').strip().replace(chr(10), ' ')}")
    _INDEXABLE_REF_CACHE[ref_fasta] = ref_fasta
    return ref_fasta


def extract_seed_sequence(ref_fasta, chrom, end_type, chrom_len, seed_len=2000, chrom_start=0):
    """Extract terminal genomic sequence to seed the walker."""
    ref_fasta = require_indexable_reference(ref_fasta)
    if end_type == "5p":
        r_start = max(1, chrom_start + 1)
        r_end = min(chrom_len, chrom_start + seed_len)
        direction = "forward"
    else:
        r_start = max(1, chrom_len - seed_len + 1)
        r_end = chrom_len
        direction = "backward"

    cmd = ["samtools", "faidx", ref_fasta, f"{chrom}:{r_start}-{r_end}"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not res.stdout.strip():
        raise RuntimeError(f"could not read the seed region {chrom}:{r_start}-{r_end} from {ref_fasta}: "
                           f"{(res.stderr or '').strip().replace(chr(10), ' ')}")
    lines = [line.strip() for line in res.stdout.splitlines() if line.strip() and not line.startswith(">")]
    seed_seq = "".join(lines).upper()

    if direction == "forward":
        return revcomp(seed_seq)
    return seed_seq


def extract_sequence_by_name(ref_fasta, name):
    """Full sequence of one record from the reference, or a loud failure."""
    ref = require_indexable_reference(ref_fasta)
    try:
        res = subprocess.run(["samtools", "faidx", ref, name], capture_output=True, text=True)
    except OSError as e:
        raise RuntimeError(f"cannot run samtools to read {name} from {ref}: {e}")
    if res.returncode != 0 or not res.stdout.strip():
        raise RuntimeError(f"record {name} is not present in {ref}: {(res.stderr or '').strip().replace(chr(10), ' ')}")
    body = "".join(l.strip() for l in res.stdout.splitlines() if not l.startswith(">"))
    return body.upper()


def parse_unplaced_extensions(paf_path, dispositions, min_overhang_bp=200, min_mapq=30, terminal_window_bp=10000):
    """Map chromosome termini to unplaced sequences that CONTINUE them.

    Honours pre-computed directional slice coordinates (`ext_qstart`, `ext_qend`) from Module 2's
    `scaffold_dispositions` when present, and falls back to directional outward overhang geometry
    on `unplaced.paf`. British spelling standard enforced.
    """
    nominations = {}
    if not dispositions:
        return nominations

    eligible_by_name = {}
    for d in dispositions:
        if d.get("action") == "INTEGRATE_INTO_CHROMOSOME" and d.get("scaffold"):
            eligible_by_name[d["scaffold"]] = d
            if (
                d.get("chromosome")
                and d.get("terminus")
                and "ext_qstart" in d
                and "ext_qend" in d
                and int(d["ext_qend"]) - int(d["ext_qstart"]) >= min_overhang_bp
            ):
                key = (d["chromosome"], d["terminus"])
                oh = int(d["ext_qend"]) - int(d["ext_qstart"])
                prev = nominations.get(key)
                if prev is None or oh > prev["overhang_bp"]:
                    nominations[key] = {
                        "scaffold": d["scaffold"],
                        "scaffold_length": int(d.get("length", oh)),
                        "target_chrom": d["chromosome"],
                        "terminus": d["terminus"],
                        "overhang_bp": oh,
                        "mapq": int(d.get("mapq", 0)),
                        "strand": d.get("strand", "+"),
                        "ext_qstart": int(d["ext_qstart"]),
                        "ext_qend": int(d["ext_qend"]),
                        "nominating_action": "INTEGRATE_INTO_CHROMOSOME",
                    }

    if not eligible_by_name or not paf_path or not os.path.exists(paf_path):
        return nominations

    with open(paf_path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 12:
                continue
            qname = p[0]
            if qname not in eligible_by_name:
                continue
            disp = eligible_by_name[qname]
            saw_tp = next((t for t in p[12:] if t.startswith("tp:A:")), None)
            if saw_tp and saw_tp != "tp:A:P" and disp.get("valency_mode") != "CHROMOSOME_EXCLUSIVE_PORT":
                continue
            try:
                qlen = int(p[1])
                qstart = int(p[2])
                qend = int(p[3])
                strand = p[4]
                tname = p[5]
                tlen = int(p[6])
                tstart = int(p[7])
                tend = int(p[8])
                mapq = int(p[11])
            except ValueError:
                continue
            if tend >= tlen - terminal_window_bp:
                term = "3p"
                tip_gap = max(0, tlen - tend)
            elif tstart <= terminal_window_bp:
                term = "5p"
                tip_gap = tstart
            else:
                continue
            if (term == "5p" and strand == "+") or (term == "3p" and strand == "-"):
                overhang, ext_s, ext_e = qstart, 0, qstart
            else:
                overhang, ext_s, ext_e = max(0, qlen - qend), qend, qlen
            net_beyond = overhang - tip_gap
            if overhang < min_overhang_bp or net_beyond <= 0:
                continue
            key = (tname, term)
            best = nominations.get(key)
            if best is None or overhang > best["overhang_bp"]:
                nominations[key] = {
                    "scaffold": qname, "scaffold_length": qlen, "target_chrom": tname,
                    "terminus": term, "overhang_bp": overhang, "mapq": mapq, "strand": strand,
                    "aligned_qstart": qstart, "aligned_qend": qend,
                    "ext_qstart": ext_s, "ext_qend": ext_e,
                    "target_tlen": tlen, "target_tstart": tstart, "target_tend": tend,
                    "nominating_action": "INTEGRATE_INTO_CHROMOSOME",
                }
    return nominations


def materialise_unplaced_extension(nomination, ref_fasta):
    """Return (extension_sequence, note) for a nomination in chromosome '+' orientation."""
    seq = extract_sequence_by_name(ref_fasta, nomination["scaffold"])
    if "ext_qstart" in nomination and "ext_qend" in nomination:
        ext_s = int(nomination["ext_qstart"])
        ext_e = int(nomination["ext_qend"])
    else:
        qstart = int(nomination["aligned_qstart"])
        qend = int(nomination["aligned_qend"])
        if (nomination["terminus"] == "5p" and nomination["strand"] == "+") or (
            nomination["terminus"] == "3p" and nomination["strand"] == "-"
        ):
            ext_s, ext_e = 0, qstart
        else:
            ext_s, ext_e = qend, len(seq)

    raw_slice = seq[ext_s:ext_e]
    if nomination.get("strand") == "-":
        ext = revcomp(raw_slice)
        note = f"contig slice {ext_s}..{ext_e}, reverse-complemented into chromosome orientation"
    else:
        ext = raw_slice
        note = f"contig slice {ext_s}..{ext_e} in forward chromosome orientation"
    if len(ext) < 1:
        raise RuntimeError(
            f"nomination for {nomination['scaffold']} on {nomination['target_chrom']} "
            f"{nomination['terminus']} leaves no sequence beyond its anchored alignment")
    return ext, note



class BlueprintParseResult(tuple):
    """Result of parse_blueprint: behaves as a 3-tuple (term_cards, policy, dispositions)
    for backwards compatibility with existing callers, while exposing
    terminal_excisions_5p and terminal_excisions_3p as attributes.
    """

    def __new__(
        cls,
        term_cards,
        policy,
        dispositions,
        terminal_excisions_5p=None,
        terminal_excisions_3p=None,
    ):
        instance = super().__new__(cls, (term_cards, policy, dispositions))
        instance.term_cards = term_cards
        instance.policy = policy
        instance.dispositions = dispositions
        instance.terminal_excisions_5p = terminal_excisions_5p or {}
        instance.terminal_excisions_3p = terminal_excisions_3p or {}
        return instance


def parse_blueprint(blueprint_path, chrom_lens=None):
    """Parse curation action cards for Domain 2 (Terminus) and Domain 3 (Terminal rDNA),

    the policy, Module 2's per-sequence scaffold dispositions, and dynamic terminal excisions
    from curation_manifest.json.
    """
    term_cards = {}
    policy = {}
    dispositions = []
    terminal_excisions_5p = {}
    terminal_excisions_3p = {}
    if not blueprint_path or not os.path.exists(blueprint_path):
        return BlueprintParseResult(term_cards, policy, dispositions, terminal_excisions_5p, terminal_excisions_3p)

    if chrom_lens is None:
        chrom_lens = {}

    try:
        with open(blueprint_path) as f:
            data = json.load(f)

        cards = []
        all_cards = []
        if isinstance(data, list):
            cards = data
            all_cards = data
        elif isinstance(data, dict):
            policy = data.get("policy", {})
            dispositions = data.get("scaffold_dispositions", []) or []
            doms = data.get("domains", {})
            if isinstance(doms, dict):
                cards.extend(doms.get("DOMAIN_TERMINUS", []))
                cards.extend(doms.get("DOMAIN_RDNA", []))
                for d_cards in doms.values():
                    if isinstance(d_cards, list):
                        all_cards.extend(d_cards)
            if not cards:
                cards = data.get("action_cards", data.get("telomere_cards", data.get("curation_cards", [])))
            if not all_cards:
                all_cards = cards

        for c in cards:
            dom = c.get("domain", "")
            if dom in ("DOMAIN_SCAFFOLD", "DOMAIN_GAP"):
                continue
            if c.get("card_id", "").startswith("CARD-SCAFF"):
                continue
            if c.get("scaffold") and dom != "DOMAIN_TERMINUS":
                continue

            chrom = c.get("chrom")
            terminus = c.get("terminus")
            if chrom and terminus:
                existing = term_cards.get((chrom, terminus))
                if existing and (existing.get("domain") == "DOMAIN_TERMINUS" or existing.get("card_id", "").startswith("CARD-TELO")):
                    continue
                term_cards[(chrom, terminus)] = c

        # Dynamic terminal excision tracking across all domains
        for ic in all_cards:
            act = ic.get("action", "")
            c_target = ic.get("chrom", "")
            s = safe_int(ic.get("start", ic.get("start_pos", 0)))
            e_pos = safe_int(ic.get("end", ic.get("end_pos", 0)))
            c_l = chrom_lens.get(c_target, 0)
            if act in ("EXCISE_DUMMY_GAP", "EXCISE_GAP", "PRUNE_CHIMERIC_TERMINUS"):
                if s == 0 and e_pos > 0:
                    terminal_excisions_5p[c_target] = max(terminal_excisions_5p.get(c_target, 0), e_pos)
                if c_l > 0 and e_pos >= c_l and s < c_l:
                    terminal_excisions_3p[c_target] = min(terminal_excisions_3p.get(c_target, c_l), s)
            if "effective_chrom_start" in ic:
                terminal_excisions_5p[c_target] = max(terminal_excisions_5p.get(c_target, 0), safe_int(ic["effective_chrom_start"]))
            if "effective_chrom_len" in ic and safe_int(ic["effective_chrom_len"]) > 0:
                terminal_excisions_3p[c_target] = min(terminal_excisions_3p.get(c_target, c_l or 10**12), safe_int(ic["effective_chrom_len"]))
    except Exception as e:
        sys.stderr.write(f"Warning: Failed to parse blueprint JSON {blueprint_path}: {e}\n")

    return BlueprintParseResult(term_cards, policy, dispositions, terminal_excisions_5p, terminal_excisions_3p)


def extract_cluster_reads_from_card(card):
    """Extract verified list of cluster read QNAMEs from action card if and only if it represents a cluster."""
    if not card or not isinstance(card, dict):
        return None
    evidence = card.get("evidence", {})
    if not isinstance(evidence, dict):
        return None
    is_cluster_card = (
        "junction_cluster_id" in evidence
        or card.get("mechanism") == "SCAFFOLD_REPEAT_TO_TELOMERE"
        or card.get("classification") == "SATELLITE_LINKED_TELOMERE"
    )
    if is_cluster_card:
        supp = evidence.get("supporting_reads")
        if isinstance(supp, (list, tuple, set)):
            return [str(r) for r in supp if r]
    return None


def run_self_test():
    """Verify internal modular packages and dependency contracts."""
    sha = get_script_fingerprint()
    print(f"pde_telomere_walker self-test: schema={SCHEMA_VERSION} sha256={sha}")
    assert revcomp("AAGCTA") == "TAGCTT", "revcomp test failed"
    assert count_telomeric_bp("TTAGGCTTAGGC", "TTAGGC", "GCCTAA") == 12, "count_telomeric_bp test failed"

    # Test satellite scaffolding sub-module
    dummy_card = {
        "gap_size": 100,
        "evidence": {
            "telomeric_tract_bp": 90,
            "supporting_reads": ["r1", "r2", "r3"]
        }
    }
    scaff_res = materialise_satellite_scaffold(dummy_card, "chr1", "5p", "TTAGGC", "GCCTAA", gap_size=100)
    assert scaff_res["status"] == "RESOLVED_BY_SATELLITE_SCAFFOLD"
    assert scaff_res["net_seq"].endswith("N" * 100)
    assert scaff_res["read_count"] == 3

    # Test 1-step walker information-gain / valency gate
    seed_sat = ("AACCGGTTAACCGGTT" * 125)  # 2,000 bp low-complexity satellite
    # 1. Stationary satellite extension -> ROLLBACK
    sat_cand = {"id": "r_sat", "ext": "AACCGGTTAACCGGTT" * 50, "mapq": 60}
    gate_sat = evaluate_step_information_gain(
        consensus_seq=seed_sat,
        ext_seq=sat_cand["ext"],
        best_cand=sat_cand,
        all_candidates=[sat_cand],
        motif="TTAGGC",
        rev_motif="GCCTAA",
        min_telo_bp=60,
    )
    assert gate_sat["decision"] == "ROLLBACK" and "Stationary repeat array" in gate_sat["reason"], (
        f"Expected stationary repeat rollback, got {gate_sat}"
    )

    # 2. Telomere reached -> TELOMERE_REACHED (even if mapq=0 in subtelomeric repeat)
    telo_cand = {"id": "r_telo", "ext": ("AACCGGTT" * 10) + ("TTAGGC" * 20), "mapq": 0}
    gate_telo = evaluate_step_information_gain(
        consensus_seq=seed_sat,
        ext_seq=telo_cand["ext"],
        best_cand=telo_cand,
        all_candidates=[telo_cand],
        motif="TTAGGC",
        rev_motif="GCCTAA",
        min_telo_bp=60,
    )
    assert gate_telo["decision"] == "TELOMERE_REACHED", f"Expected TELOMERE_REACHED, got {gate_telo}"

    # 3. Polyvalent branching -> ROLLBACK
    import random
    rng = random.Random(42)
    uniq_ext_a = "".join(rng.choice("ACGT") for _ in range(500))
    uniq_ext_b = "".join(rng.choice("ACGT") for _ in range(500))
    cand_a = {"id": "r_a", "ext": uniq_ext_a, "mapq": 60}
    cand_b = {"id": "r_b", "ext": uniq_ext_b, "mapq": 60}
    gate_branch = evaluate_step_information_gain(
        consensus_seq=seed_sat,
        ext_seq=uniq_ext_a,
        best_cand=cand_a,
        all_candidates=[cand_a, cand_b],
        motif="TTAGGC",
        rev_motif="GCCTAA",
        min_telo_bp=60,
    )
    assert gate_branch["decision"] == "ROLLBACK" and "Polyvalent" in gate_branch["reason"], (
        f"Expected Polyvalent branch rollback, got {gate_branch}"
    )

    # 4. Unique exit -> UNIQUE_EXIT
    gate_uniq = evaluate_step_information_gain(
        consensus_seq=seed_sat,
        ext_seq=uniq_ext_a,
        best_cand=cand_a,
        all_candidates=[cand_a],
        motif="TTAGGC",
        rev_motif="GCCTAA",
        min_telo_bp=60,
    )
    assert gate_uniq["decision"] == "UNIQUE_EXIT", f"Expected UNIQUE_EXIT, got {gate_uniq}"

    # 5. Unplaced extension orientation across all 4 (terminus, strand) pairs
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_fa = os.path.join(tmpdir, "ref.fa")
        # Contig scaff1: LEFT='AAAAACCCCC' (0..10), CORE='GGGGGGGGGG' (10..20), RIGHT='TTTTTAAAAA' (20..30)
        # Contig scaff2: LEFT='AAACCCGGG' (0..9), CORE='GGGGGGGGGG' (9..19), RIGHT='TTTCCCAAA' (19..28)
        with open(ref_fa, "w") as f:
            f.write(">scaff1\nAAAAACCCCCGGGGGGGGGGTTTTTAAAAA\n>scaff2\nAAACCCGGGGGGGGGGGGGTTTCCCAAA\n")

        # 5p, + strand: outward overhang is [0..10] ('AAAAACCCCC') in forward chromosome orientation
        ext_5p_plus, _ = materialise_unplaced_extension(
            {"scaffold": "scaff1", "target_chrom": "chr1", "terminus": "5p", "strand": "+", "aligned_qstart": 10, "aligned_qend": 20},
            ref_fa,
        )
        assert ext_5p_plus == "AAAAACCCCC", f"5p/+ mismatch: {ext_5p_plus}"

        # 5p, - strand: outward overhang is [19..28] ('TTTCCCAAA'), reverse-complemented -> 'TTTGGGAAA'
        ext_5p_minus, _ = materialise_unplaced_extension(
            {"scaffold": "scaff2", "target_chrom": "chr1", "terminus": "5p", "strand": "-", "aligned_qstart": 9, "aligned_qend": 19},
            ref_fa,
        )
        assert ext_5p_minus == revcomp("TTTCCCAAA"), f"5p/- mismatch: {ext_5p_minus}"

        # 3p, + strand: outward overhang is [19..28] ('TTTCCCAAA') in forward orientation
        ext_3p_plus, _ = materialise_unplaced_extension(
            {"scaffold": "scaff2", "target_chrom": "chr1", "terminus": "3p", "strand": "+", "aligned_qstart": 9, "aligned_qend": 19},
            ref_fa,
        )
        assert ext_3p_plus == "TTTCCCAAA", f"3p/+ mismatch: {ext_3p_plus}"

        # 3p, - strand: outward overhang is [0..9] ('AAACCCGGG'), reverse-complemented -> 'CCCGGGTTT'
        ext_3p_minus, _ = materialise_unplaced_extension(
            {"scaffold": "scaff2", "target_chrom": "chr1", "terminus": "3p", "strand": "-", "aligned_qstart": 9, "aligned_qend": 19},
            ref_fa,
        )
        assert ext_3p_minus == revcomp("AAACCCGGG"), f"3p/- mismatch: {ext_3p_minus}"

    print("Self-test PASS")
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(description="Modular Telomere Extension, Capping, Scaffolding & Walker Engine.")
    parser.add_argument("--ref", help="Reference FASTA")
    parser.add_argument("--bam", help="Sorted HiFi BAM")
    parser.add_argument("--blueprint", default=None, help="Path to Module 2 curation_manifest.json")
    parser.add_argument("--telomere-profile", default=None,
                        help="Optional PROFILE_TELOMERES TSV (unitig alignments with query and target coordinates)")
    parser.add_argument("--unitigs", default=None, help="Path to assembled telomeric unitigs FASTA")
    parser.add_argument("--telo-motif", default=None, help="Canonical telomere repeat motif")
    parser.add_argument("--chromosomes", default=None, help="Comma-separated chromosome whitelist")
    parser.add_argument("--sample", default="SAMPLE", help="Sample ID")
    parser.add_argument("--outdir", default=None, help="Output directory for extensions")
    parser.add_argument("--manifest", default=None, help="Output manifest TSV")
    parser.add_argument("--unplaced-paf", default=None, help="PAF of non-chromosomal sequences against chromosomes")
    parser.add_argument("--threads", type=int, default=1, help="Threads for minimap2")
    parser.add_argument("--self-test", action="store_true", help="Run self-test and exit")
    parser.add_argument("--version", action="store_true", help="Print version and exit")
    args = parser.parse_args()

    if args.version:
        print(f"pde_telomere_walker schema={SCHEMA_VERSION} sha256={get_script_fingerprint()}")
        sys.exit(0)

    if args.self_test:
        run_self_test()

    if not (args.ref and args.bam):
        parser.error("--ref and --bam are required when not running --self-test or --version")

    sample_id = args.sample
    outdir = args.outdir or f"{sample_id}_telomere_extensions"
    manifest_path = args.manifest or f"{sample_id}_walker_manifest.tsv"
    os.makedirs(outdir, exist_ok=True)

    chrom_lens = parse_chrom_lens_from_bam(args.bam)
    target_chroms = [c.strip() for c in args.chromosomes.split(",")] if args.chromosomes else list(chrom_lens.keys())
    target_chroms = [c for c in target_chroms if c in chrom_lens]

    if args.blueprint:
        bp_res = parse_blueprint(args.blueprint, chrom_lens=chrom_lens)
        term_cards, policy, scaffold_dispositions = bp_res
        terminal_excisions_5p = bp_res.terminal_excisions_5p
        terminal_excisions_3p = bp_res.terminal_excisions_3p
    else:
        term_cards, policy, scaffold_dispositions = {}, {}, []
        terminal_excisions_5p, terminal_excisions_3p = {}, {}
    has_blueprint = bool(term_cards or args.blueprint)

    telo_policy = policy.get("telomere", {})
    arch_policy = policy.get("genomic_architecture", {})
    map_policy = policy.get("mapping_and_anchors", {})

    motif = (telo_policy.get("motif") or args.telo_motif or "TTAGGC").upper()
    rev_motif = revcomp(motif)
    min_telo_repeats = safe_int(telo_policy.get("min_repeats", 15), default=15)
    min_extension_bp = safe_int(telo_policy.get("min_extension_bp", 200), default=200)
    terminal_window_bp = safe_int(telo_policy.get("terminal_window_bp", 10000), default=10000)
    seed_len_bp = safe_int(telo_policy.get("seed_len_bp", 2000), default=2000)
    max_iters = safe_int(telo_policy.get("walker_max_iters"))
    germline_buffer_bp = safe_int(arch_policy.get("germline_buffer_bp", 100), default=100)
    mapq_anchored = safe_int(map_policy.get("mapq_anchored", 30), default=30)

    print(f"[{sample_id}] Evaluating {len(target_chroms)} chromosomes (2N = {2 * len(target_chroms)} termini)...", flush=True)

    walk_cache = {"index": None, "tmpdir": None}

    def shared_unmapped():
        if walk_cache["index"] is None:
            walk_cache["tmpdir"] = tempfile.TemporaryDirectory()
            walk_cache["index"] = UnmappedReadIndex(args.bam, walk_cache["tmpdir"].name)
        return walk_cache["index"]

    unplaced_nominations = parse_unplaced_extensions(
        args.unplaced_paf, scaffold_dispositions,
        min_overhang_bp=min_extension_bp,
        min_mapq=mapq_anchored,
        terminal_window_bp=terminal_window_bp
    )
    if unplaced_nominations:
        print(f"[{sample_id}] Module 2 nominated unplaced sequences for {len(unplaced_nominations)} terminus/termini", flush=True)

    profile_hits = parse_telomere_profile(args.telomere_profile, min_repeats=min_telo_repeats, min_mapq=mapq_anchored)

    manifest_rows = []
    telomere_lengths = []
    header = [
        "chrom", "terminus", "status", "initial_evidence", "added_bp",
        "canonical_telo_bp", "non_telo_added_bp", "extension_fasta", "notes"
    ]

    for chrom in target_chroms:
        for terminus in ["5p", "3p"]:
            orig_len = chrom_lens[chrom]
            c_start = terminal_excisions_5p.get(chrom, 0)
            c_len = terminal_excisions_3p.get(chrom, orig_len)
            term_id = f"{chrom}_{terminus}"
            print(f"[{sample_id}] Evaluating terminus: {term_id}...", flush=True)

            if has_blueprint:
                card = term_cards.get((chrom, terminus))
                classification = card.get("classification", "") if card else ""
                action = card.get("action", "") if card else ""
                action_subcode = card.get("action_subcode", "") if card else ""
                mechanism = card.get("mechanism", "") if card else ""
                tip_is_germline = bool(card.get("terminus_is_germline", False)) if card else False

                if card:
                    if terminus == "5p" and "effective_chrom_start" in card:
                        c_start = max(c_start, safe_int(card["effective_chrom_start"]))
                    elif terminus == "3p" and "effective_chrom_len" in card:
                        c_len = min(c_len, safe_int(card["effective_chrom_len"]))

                # 1. Baseline Capped: SKIP_WALKER or ALREADY_CAPPED
                if card is None or action == "SKIP_WALKER" or classification == "ALREADY_CAPPED":
                    manifest_rows.append([
                        chrom, terminus, "ALREADY_CAPPED", "DRAFT_VERIFIED_COMPLETE",
                        0, 0, 0, "NONE",
                        "Designated ALREADY_CAPPED / SKIP_WALKER by Module 2. Retained as-is in baseline assembly."
                    ])
                    continue

                # 2. Major rDNA Array Precedence: Bypass Walker for Module 4 Cassette
                if classification == "TERMINUS_RDNA_ARRAY" or action in ("TERMINUS_RDNA_ARRAY", "BYPASS_WALKER"):
                    print(f"[{sample_id}] Terminus {term_id} is designated TERMINUS_RDNA_ARRAY. Bypassing Telomere Walker for Module 4 Cassette.", flush=True)
                    manifest_rows.append([
                        chrom, terminus, "TERMINUS_RDNA_ARRAY", "MAJOR_RDNA_CLUSTER",
                        0, 0, 0, "NONE",
                        "Directly abuts multi-copy rDNA array. Bypassed by Telomere Walker for Module 4 tripartite cassette resolution."
                    ])
                    continue

                # 2b. Satellite Repeat Scaffolding to Distal Telomeric Cap
                if classification == "SATELLITE_LINKED_TELOMERE" or action == "SCAFFOLD_REPEAT_TO_TELOMERE":
                    evidence = card.get("evidence", {}) or {}
                    assoc_scaffs = evidence.get("associated_scaffolds")
                    if not assoc_scaffs and evidence.get("associated_scaffold"):
                        assoc_scaffs = [evidence.get("associated_scaffold")]
                    combined_scaff_seq = ""
                    if assoc_scaffs and args.ref and os.path.exists(args.ref):
                        seq_parts = []
                        for s_name in assoc_scaffs:
                            try:
                                s_seq = extract_sequence_by_name(args.ref, s_name)
                                strand = evidence.get("scaffold_orientations", {}).get(s_name, "+")
                                if strand == "-":
                                    s_seq = revcomp(s_seq)
                                seq_parts.append(s_seq)
                            except Exception as ex:
                                sys.stderr.write(f"Warning: could not extract scaffold {s_name} from {args.ref}: {ex}\n")
                        if seq_parts:
                            min_scaff_overlap = safe_int(telo_policy.get("subtelomeric_min_overlap_bp", 500), default=500)
                            combined_scaff_seq = seq_parts[0]
                            for next_part in seq_parts[1:]:
                                ov = find_sequence_overlap(combined_scaff_seq, next_part, min_overlap=min_scaff_overlap)
                                if ov > 0:
                                    combined_scaff_seq += next_part[ov:]
                                else:
                                    combined_scaff_seq += next_part

                    tip_seq = extract_seed_sequence(args.ref, chrom, terminus, c_len, seed_len=terminal_window_bp, chrom_start=c_start) if (args.ref and os.path.exists(args.ref)) else ""
                    scaffold_res = materialise_satellite_scaffold(
                        card, chrom, terminus, motif, rev_motif,
                        gap_size=safe_int(card.get("gap_size", 0)),
                        scaffold_seq=combined_scaff_seq,
                        tip_seq=tip_seq,
                        bam_path=args.bam,
                        policy=telo_policy
                    )
                    ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                        outdir, term_id, "satellite_scaffold_ext",
                        f"status={scaffold_res['status']} reads={scaffold_res['read_count']}",
                        scaffold_res["net_seq"], motif, rev_motif)
                    telomere_lengths.append((chrom, terminus, scaffold_res["status"], "SATELLITE_LINKED",
                                             scaffold_res["telo_len"], scaffold_res["added_len"]))
                    manifest_rows.append([
                        chrom, terminus, scaffold_res["status"], f"SATELLITE_ARRAY_LINKED_{scaffold_res['read_count']}_READS",
                        added_len, telo_len, non_telo, ext_fa, scaffold_res["notes"]
                    ])
                    continue

                # 3. Somatic Diminution Boundary: Enforce Somatic Exclusion Zone (Dsoma + germline_buffer_bp)
                if classification == "SOMATIC_DIMINUTION_BOUNDARY" or action_subcode == "WALK_PAST_BREAKSITE":
                    d_soma = card.get("somatic_boundary_coord", card.get("forbidden_zone_bp", 0))
                    min_buf = card.get("min_germline_extension_bp", germline_buffer_bp)
                    term_max_iters = safe_int(card.get("walker_max_iters") if card else None, default=max_iters)
                    seed_seq = extract_seed_sequence(args.ref, chrom, terminus, c_len, seed_len=seed_len_bp, chrom_start=c_start)
                    walk_res = run_iterative_walker(
                        args.bam, seed_seq, motif, rev_motif,
                        max_iters=term_max_iters,
                        min_telo_bp=min_extension_bp,
                        is_somatic_boundary=True,
                        min_germline_extension_bp=min_buf,
                        unmapped=shared_unmapped(),
                        threads=args.threads,
                        tip_is_germline=tip_is_germline,
                        chrom=chrom,
                        terminus=terminus,
                        c_len=c_len,
                        c_start=c_start,
                        cluster_reads=extract_cluster_reads_from_card(card),
                    )
                    ext_fa, added_len, telo_len, non_telo = "NONE", 0, 0, 0
                    if walk_res["added_bp"] > 0:
                        ext_seq_chrom = (
                            revcomp(walk_res["extension_seq"])
                            if terminus == "5p"
                            else walk_res["extension_seq"]
                        )
                        ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                            outdir, term_id, "walker_ext", f"status={walk_res['status']}",
                            ext_seq_chrom, motif, rev_motif)

                    manifest_rows.append([
                        chrom, terminus, walk_res["status"], f"SOMATIC_BOUNDARY_{d_soma}",
                        added_len, telo_len, non_telo, ext_fa,
                        f"{walk_res['notes'] or ''} | {added_len} bp added, {telo_len} bp canonical telomere, {non_telo} bp non-telomeric"
                    ])

                    continue

                # 4. Telomere at the draft tip
                if action in ("EXTEND_TELOMERE_TRACT", "EXTEND_FROM_OVERHANG") or classification in (
                        "DRAFT_HAS_TELOMERE", "TELOMERE_IN_OVERHANG"):
                    terminal_seq = extract_seed_sequence(args.ref, chrom, terminus, c_len,
                                                         seed_len=terminal_window_bp, chrom_start=c_start)
                    candidates = terminal_read_candidates(args.bam, chrom, terminus, c_len, motif,
                                                          rev_motif, window=terminal_window_bp, chrom_start=c_start)
                    best_tract = maximal_terminal_tract(
                        candidates, motif, rev_motif, terminal_seq,
                        anchor_bp=safe_int(telo_policy.get("tract_anchor_bp", 500), default=500),
                        min_extra_bp=min_telo_repeats * len(motif))
                    if best_tract:
                        qname, total_tract, extra_seq, orientation, _insert_at = best_tract
                        ext_seq_chrom = revcomp(extra_seq) if terminus == "5p" else extra_seq
                        ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                            outdir, term_id, "telo_ext",
                            f"qname={qname} strand={orientation} tract_bp={total_tract}",
                            ext_seq_chrom, motif, rev_motif)
                        telomere_lengths.append((chrom, terminus, "EXTEND_TELOMERE_TRACT", qname,
                                                 total_tract, added_len))
                        manifest_rows.append([
                            chrom, terminus, "EXTEND_TELOMERE_TRACT", f"Read_{qname}",
                            added_len, telo_len, non_telo, ext_fa,
                            f"Telomeric tract extended to {total_tract} bp using read {qname} ({added_len} bp of canonical repeat "
                            f"appended; the tract was {total_tract - added_len} bp)"
                        ])
                        continue
                    else:
                        existing = terminal_tract_region(terminal_seq, motif, rev_motif, "3p",
                                                         min_repeats=min_telo_repeats)[0]
                        if existing > 0:
                            telomere_lengths.append((chrom, terminus, "RETAIN_AS_IS", "NONE", existing, 0))
                            manifest_rows.append([
                                chrom, terminus, "RETAIN_AS_IS", "DRAFT_TELOMERE_CONFIRMED",
                                0, 0, 0, "NONE",
                                f"Draft telomeric tract is {existing} bp; no read follows it further out. Retained as-is."
                            ])
                            continue

                        print(f"[{sample_id}] {term_id}: no telomeric tract at the draft tip; trying endorsed candidates.", flush=True)
                        tried = []

                        # (a) The unitig Module 2 endorsed
                        materialised, why = authorized_unitig_materialisation(
                            args.unitigs, card, chrom, terminus, c_len,
                            profile_hits=profile_hits)
                        if materialised:
                            net_seq = materialised["net_seq"]
                            card_evidence = card.get("evidence", {}) or {}
                            distal_tract_bp = safe_int(card_evidence.get("distal_telomere_tract_bp", 0))
                            distal_cl_id = card_evidence.get("distal_cluster_id")
                            distal_reads = safe_int(card_evidence.get("distal_cluster_reads", 0))
                            is_unresolved_repeat = bool(card_evidence.get("unresolved_repeat_gap", False))

                            cap_motif = rev_motif if terminus == "5p" else motif
                            existing_tract = terminal_tract_region(
                                net_seq, motif, rev_motif, terminus, min_repeats=min_telo_repeats)[0]

                            added_overhang = 0
                            if distal_tract_bp > 0:
                                overhang_needed = max(0, distal_tract_bp - existing_tract)
                                if overhang_needed > 0:
                                    repeats_to_add = max(15, overhang_needed // len(cap_motif))
                                    overhang_seq = cap_motif * repeats_to_add
                                    added_overhang = len(overhang_seq)
                                    if is_unresolved_repeat:
                                        gap_seq = "N" * 100
                                        if terminus == "5p":
                                            net_seq = overhang_seq + gap_seq + net_seq
                                        else:
                                            net_seq = net_seq + gap_seq + overhang_seq
                                        added_overhang += 100
                                    else:
                                        if terminus == "5p":
                                            net_seq = overhang_seq + net_seq
                                        else:
                                            net_seq = net_seq + overhang_seq

                            ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                                outdir, term_id, "tip_unitig_ext",
                                f"unitig_id={materialised['unitig_id']} from=EXTEND_FROM_OVERHANG",
                                net_seq, motif, rev_motif)
                            tip_status, tract_bp = cap_status_for(
                                "RESOLVED_BY_UNITIG", net_seq, motif, rev_motif, min_telo_repeats,
                                terminus=terminus)
                            telomere_lengths.append((chrom, terminus, tip_status,
                                                     materialised["unitig_id"], tract_bp, telo_len))
                            notes_cluster = f" + distal cap {distal_cl_id} ({added_overhang} bp overhang, {distal_reads} reads)" if distal_tract_bp > 0 and added_overhang > 0 else ""
                            manifest_rows.append([
                                chrom, terminus, tip_status,
                                f"Unitig_{materialised['unitig_id']}",
                                added_len, telo_len, non_telo, ext_fa,
                                f"EXTEND_FROM_OVERHANG: used unitig Module 2 endorsed for this end ({added_len} bp of "
                                f"{materialised['unitig_len']} bp; {telo_len} bp canonical telomere; terminal tract {tract_bp} bp{notes_cluster})"
                            ])
                            continue
                        if why:
                            tried.append(f"endorsed unitig: {why}")

                        # (b) The read Module 2 endorsed
                        overhang = None
                        endorsed_qname = card.get("authorized_read_qname")
                        if endorsed_qname:
                            clip_bp, overhang_seq = named_read_overhang(
                                args.bam, endorsed_qname, chrom, terminus, c_len,
                                window=terminal_window_bp, chrom_start=c_start)
                            if overhang_seq:
                                overhang = (endorsed_qname, clip_bp, overhang_seq)
                            else:
                                tried.append(f"endorsed read {endorsed_qname} does not reach past the tip")

                        # (c) Maximal tip overhang
                        if overhang is None:
                            overhang = maximal_tip_overhang(args.bam, chrom, terminus, c_len,
                                                            window=terminal_window_bp, chrom_start=c_start)
                        if overhang:
                            qname, clip_bp, overhang_seq = overhang
                            ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                                outdir, term_id, "tip_overhang",
                                f"qname={qname} overhang_bp={clip_bp}",
                                overhang_seq, motif, rev_motif)
                            tip_status, tract_bp = cap_status_for(
                                "EXTEND_TO_TELOMERE", overhang_seq, motif, rev_motif, min_telo_repeats,
                                terminus=terminus)
                            telomere_lengths.append((chrom, terminus, tip_status, qname, tract_bp, telo_len))
                            manifest_rows.append([
                                chrom, terminus, tip_status, f"Read_{qname}",
                                added_len, telo_len, non_telo, ext_fa,
                                f"Extended from the draft tip to the end of read {qname} ({clip_bp} bp overhang; "
                                f"{telo_len} bp canonical telomere; terminal tract {tract_bp} bp)"
                            ])
                            continue

                        # (d) Fallback to satellite repeat scaffolding if designated cluster or scaffolds exist
                        card_evidence = (card.get("evidence", {}) or {}) if card else {}
                        distal_cl_id = card_evidence.get("distal_cluster_id")
                        assoc_scaffs = card_evidence.get("associated_scaffolds") or (
                            [card_evidence.get("associated_scaffold")] if card_evidence.get("associated_scaffold") else None
                        )
                        if distal_cl_id or assoc_scaffs:
                            print(f"[{sample_id}] {term_id}: no single-read overhang; falling back to satellite scaffolding for cluster {distal_cl_id or assoc_scaffs}.", flush=True)
                            combined_scaff_seq = ""
                            if assoc_scaffs and args.ref and os.path.exists(args.ref):
                                seq_parts = []
                                for s_name in assoc_scaffs:
                                    try:
                                        s_seq = extract_sequence_by_name(args.ref, s_name)
                                        strand = card_evidence.get("scaffold_orientations", {}).get(s_name, "+")
                                        if strand == "-":
                                            s_seq = revcomp(s_seq)
                                        seq_parts.append(s_seq)
                                    except Exception as ex:
                                        sys.stderr.write(f"Warning: could not extract scaffold {s_name} from {args.ref}: {ex}\n")
                                if seq_parts:
                                    min_scaff_overlap = safe_int(telo_policy.get("subtelomeric_min_overlap_bp", 500), default=500)
                                    combined_scaff_seq = seq_parts[0]
                                    for next_part in seq_parts[1:]:
                                        ov = find_sequence_overlap(combined_scaff_seq, next_part, min_overlap=min_scaff_overlap)
                                        if ov > 0:
                                            combined_scaff_seq += next_part[ov:]
                                        else:
                                            combined_scaff_seq += next_part

                            tip_seq = extract_seed_sequence(args.ref, chrom, terminus, c_len, seed_len=terminal_window_bp, chrom_start=c_start) if (args.ref and os.path.exists(args.ref)) else ""
                            scaffold_res = materialise_satellite_scaffold(
                                card, chrom, terminus, motif, rev_motif,
                                gap_size=safe_int(card.get("gap_size", 0)),
                                scaffold_seq=combined_scaff_seq,
                                tip_seq=tip_seq,
                                bam_path=args.bam,
                                policy=telo_policy,
                            )
                            ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                                outdir, term_id, "satellite_scaffold_ext",
                                f"status={scaffold_res['status']} reads={scaffold_res['read_count']}",
                                scaffold_res["net_seq"], motif, rev_motif)
                            telomere_lengths.append((chrom, terminus, scaffold_res["status"], "SATELLITE_LINKED",
                                                     scaffold_res["telo_len"], scaffold_res["added_len"]))
                            manifest_rows.append([
                                chrom, terminus, scaffold_res["status"], f"SATELLITE_ARRAY_LINKED_{scaffold_res['read_count']}_READS",
                                added_len, telo_len, non_telo, ext_fa, scaffold_res["notes"]
                            ])
                            continue

                        print(f"[{sample_id}] {term_id}: handing end to iterative walker ({'; '.join(tried) or 'no candidates'})", flush=True)

                # 5. Phase 1 Deterministic Unitig Capping
                if mechanism == "CAP_WITH_UNITIG" or action == "CAP_WITH_UNITIG" or classification == "RESOLVABLE_BY_UNITIG":
                    unitig_id = card.get("authorized_unitig_id", "UNKNOWN_UNITIG")
                    card_evidence = card.get("evidence", {}) or {}
                    ext_bp = safe_int(card_evidence.get("unitig_extension_bp", 0))
                    telo_bp = safe_int(card_evidence.get("telomeric_tract_bp", ext_bp))
                    region = card_evidence.get("unitig_extension_region") or {}

                    materialised, why = authorized_unitig_materialisation(
                        args.unitigs, card, chrom, terminus, c_len,
                        profile_hits=profile_hits)
                    ext_fa = "NONE"
                    net_seq = materialised["net_seq"] if materialised else ""
                    unitig_seq = materialised["unitig_seq"] if materialised else ""
                    if net_seq:
                        distal_tract_bp = safe_int(card_evidence.get("distal_telomere_tract_bp", 0))
                        distal_cl_id = card_evidence.get("distal_cluster_id")
                        distal_reads = safe_int(card_evidence.get("distal_cluster_reads", 0))
                        is_unresolved_repeat = bool(card_evidence.get("unresolved_repeat_gap", False))

                        cap_motif = rev_motif if terminus == "5p" else motif
                        existing_tract = terminal_tract_region(
                            net_seq, motif, rev_motif, terminus, min_repeats=min_telo_repeats)[0]

                        added_overhang = 0
                        if distal_tract_bp > 0:
                            overhang_needed = max(0, distal_tract_bp - existing_tract)
                            if overhang_needed > 0:
                                repeats_to_add = max(15, overhang_needed // len(cap_motif))
                                overhang_seq = cap_motif * repeats_to_add
                                added_overhang = len(overhang_seq)
                                if is_unresolved_repeat:
                                    gap_seq = "N" * 100
                                    if terminus == "5p":
                                        net_seq = overhang_seq + gap_seq + net_seq
                                    else:
                                        net_seq = net_seq + gap_seq + overhang_seq
                                    added_overhang += 100
                                else:
                                    if terminus == "5p":
                                        net_seq = overhang_seq + net_seq
                                    else:
                                        net_seq = net_seq + overhang_seq

                        ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                            outdir, term_id, "unitig_cap",
                            f"unitig_id={unitig_id} unitig_len={len(unitig_seq)} region={region.get('q_start', '?')}-{region.get('q_end', '?')}",
                            net_seq, motif, rev_motif)
                        cap_status, tract_bp = cap_status_for("RESOLVED_BY_UNITIG", net_seq, motif,
                                                              rev_motif, min_telo_repeats,
                                                              terminus=terminus)
                        telomere_lengths.append((chrom, terminus, cap_status, unitig_id, tract_bp, telo_len))
                        notes_cluster = f" + distal cap {distal_cl_id} ({added_overhang} bp overhang, {distal_reads} reads)" if distal_tract_bp > 0 and added_overhang > 0 else ""
                        manifest_rows.append([
                            chrom, terminus, cap_status, f"Unitig_{unitig_id}",
                            added_len, telo_len, non_telo, ext_fa,
                            f"Capped: outboard {added_len} bp of unitig {unitig_id} ({telo_len} bp canonical telomere, "
                            f"{non_telo} bp non-telomeric; terminal tract {tract_bp} bp{notes_cluster})"
                        ])
                        continue
                    else:
                        if not why:
                            why = "unitig interval yields no sequence"
                        sys.stderr.write(
                            f"Warning: Authorized unitig {unitig_id} at {chrom} {terminus} not materialised: {why} "
                            f"(card reported {ext_bp} bp / {telo_bp} bp telomere); falling back to unplaced contigs and iterative consensus walker.\n"
                        )

                # 6. Phase 1 Deterministic Read Capping
                if mechanism == "CAP_WITH_READ" or action == "CAP_WITH_READ" or classification == "RESOLVABLE_BY_READ":
                    auth_qname = card.get("authorized_read_qname")
                    read_seq = extract_read_by_qname(args.bam, auth_qname, chrom, terminus, c_len, window=terminal_window_bp, chrom_start=c_start) if auth_qname else None
                    if not read_seq:
                        read_cap = inspect_terminal_reads(args.bam, chrom, terminus, c_len, motif, rev_motif, window=terminal_window_bp, min_repeats=min_telo_repeats, chrom_start=c_start)
                        if read_cap:
                            read_seq = read_cap["seq"]
                            auth_qname = read_cap["qname"]

                    if read_seq:
                        clip_bp, overhang_seq = named_read_overhang(args.bam, auth_qname, chrom, terminus, c_len, window=terminal_window_bp, chrom_start=c_start)
                        if overhang_seq:
                            ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                                outdir, term_id, "read_cap", f"qname={auth_qname} overhang_bp={clip_bp}",
                                overhang_seq, motif, rev_motif)
                            cap_status, tract_bp = cap_status_for("RESOLVED_BY_READ", overhang_seq, motif,
                                                                  rev_motif, min_telo_repeats, terminus=terminus)
                            telomere_lengths.append((chrom, terminus, cap_status, auth_qname, tract_bp, telo_len))
                            manifest_rows.append([
                                chrom, terminus, cap_status, f"Read_{auth_qname}",
                                added_len, telo_len, non_telo, ext_fa,
                                f"Terminal HiFi read {auth_qname} ({added_len} bp added, {telo_len} bp canonical telomere, {non_telo} bp non-telo)"
                            ])
                            continue
                        else:
                            sys.stderr.write(
                                f"Warning: Authorized read {auth_qname} at {chrom} {terminus} does not reach past the tip; "
                                f"falling back to unplaced contigs and iterative consensus walker.\n"
                            )
                    else:
                        sys.stderr.write(
                            f"Warning: Authorized read {auth_qname} at {chrom} {terminus} could not be extracted from BAM; "
                            f"falling back to unplaced contigs and iterative consensus walker.\n"
                        )

                # 7. Unplaced Contigs & Phase 2 Iterative Consensus Walker
                contig_row = None
                nom = unplaced_nominations.get((chrom, terminus))
                if nom:
                    ext, note = materialise_unplaced_extension(nom, args.ref)
                    ext_fa = os.path.join(outdir, f"{term_id}_unplaced_contig_ext.fa")
                    with open(ext_fa, "w") as f:
                        f.write(f">{term_id}_unplaced_contig_ext scaffold={nom['scaffold']} strand={nom['strand']} overhang={nom['overhang_bp']}bp\n{ext}\n")
                    telo_len = count_telomeric_bp(ext, motif, rev_motif)
                    contig_row = [
                        chrom, terminus, "EXTENDED_BY_UNPLACED_CONTIG",
                        f"UnplacedContig_{nom['scaffold']}", len(ext), telo_len, len(ext) - telo_len, ext_fa,
                        f"Terminus continued by unplaced sequence {nom['scaffold']} ({note})"
                    ]

                if contig_row is not None:
                    manifest_rows.append(contig_row)
                    continue

                term_max_iters = safe_int(card.get("walker_max_iters") if card else None, default=max_iters)
                seed_seq = extract_seed_sequence(args.ref, chrom, terminus, c_len, seed_len=seed_len_bp, chrom_start=c_start)
                walk_cluster_reads = extract_cluster_reads_from_card(card)
                walk_res = run_iterative_walker(
                    args.bam, seed_seq, motif, rev_motif,
                    max_iters=term_max_iters,
                    min_telo_bp=min_extension_bp,
                    unmapped=shared_unmapped(),
                    threads=args.threads,
                    tip_is_germline=tip_is_germline,
                    chrom=chrom,
                    terminus=terminus,
                    c_len=c_len,
                    c_start=c_start,
                    cluster_reads=walk_cluster_reads,
                )
                ext_fa, added_len, telo_len, non_telo = "NONE", 0, 0, 0
                walk_status = walk_res["status"]
                if walk_res["added_bp"] > 0:
                    ext_seq_chrom = (
                        revcomp(walk_res["extension_seq"])
                        if terminus == "5p"
                        else walk_res["extension_seq"]
                    )
                    ext_fa, added_len, telo_len, non_telo = write_extension_fasta(
                        outdir, term_id, "walker_ext", f"status={walk_res['status']}",
                        ext_seq_chrom, motif, rev_motif)
                    walk_status, _ = cap_status_for(walk_res["status"], ext_seq_chrom,
                                                    motif, rev_motif, min_telo_repeats, terminus=terminus)


                manifest_rows.append([
                    chrom, terminus, walk_status, "IterativeWalker",
                    added_len, telo_len, non_telo, ext_fa,
                    f"{walk_res['notes'] or ''} | {added_len} bp added, {telo_len} bp canonical telomere, {non_telo} bp non-telomeric"
                ])

    # Write Manifest TSV
    with open(manifest_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in manifest_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"[{sample_id}] Finished evaluating {len(target_chroms)} chromosomes. Manifest written to {manifest_path}.", flush=True)


if __name__ == "__main__":
    main()
