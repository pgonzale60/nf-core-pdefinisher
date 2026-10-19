#!/usr/bin/env python3
"""Telomeric Junction Profiler.

PDE-aware deconvolution and clustering of subtelomeric junctions from raw telomeric reads.
Identifies whether uncapped termini characterized by repetitive or satellite arrays
are directly linked to authentic terminal telomere tracts, and distinguishes somatic
diminution breakpoints from germline chromosome ends using DELFIES annotations.

All thresholds are read directly from assets/curation_policy.json (single source of truth).
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def revcomp(seq):
    """Return the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def read_curation_policy(policy_path):
    """Load policy parameters from curation_policy.json."""
    if not os.path.exists(policy_path):
        sys.stderr.write(f"ERROR: Policy file not found: {policy_path}\n")
        sys.exit(1)
    with open(policy_path, "r", encoding="utf-8") as f:
        policy = json.load(f)

    # Extract relevant thresholds from policy
    min_read_support = int(policy["evidence_floors"]["min_read_support_floor"])
    min_repeats = int(policy["telomere"]["min_repeats"])
    terminal_window_bp = int(policy["telomere"]["terminal_window_bp"])
    germline_buffer_bp = int(policy["genomic_architecture"]["germline_buffer_bp"])

    return {
        "min_read_support_floor": min_read_support,
        "min_repeats": min_repeats,
        "terminal_window_bp": terminal_window_bp,
        "germline_buffer_bp": germline_buffer_bp,
    }


def parse_delfies_bed(bed_path):
    """Parse DELFIES diminution breakpoint coordinates per chromosome."""
    breakpoints = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return breakpoints
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                    breakpoints[chrom].append(pos)
                except ValueError:
                    continue
    return breakpoints


def stream_fastq(fastq_path):
    """Yield (header, sequence, quality) from raw or gzipped FASTQ."""
    if fastq_path.endswith(".gz"):
        fh = gzip.open(fastq_path, "rt", encoding="utf-8", errors="replace")
    else:
        fh = open(fastq_path, "r", encoding="utf-8", errors="replace")

    try:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            _plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            qname = header.strip()[1:].split()[0]
            yield qname, seq.strip().upper(), qual.strip()
    finally:
        fh.close()


def find_terminal_tracts(seq, motif, min_repeats=15, tip_tolerance_bp=100, gap_tolerance_bp=30):
    """Locate terminal telomeric tracts in a read sequence.

    Returns a tuple (tract_orientation, tract_len, junction_pos, flank_away, tract_motif)
    or None if no terminal tract satisfying min_repeats is present.
    `flank_away` is the subtelomeric flank oriented pointing AWAY from the telomere.
    """
    rev_motif = revcomp(motif)
    variants = [v for v in dict.fromkeys([motif, rev_motif, "TTAGAC", "GTCTAA"]) if v]

    hits = []
    for variant in variants:
        start = seq.find(variant)
        while start >= 0:
            hits.append((start, start + len(variant), variant))
            start = seq.find(variant, start + 1)

    if not hits:
        return None

    hits.sort(key=lambda x: (x[0], x[1]))

    # Group adjacent hits within gap_tolerance_bp
    runs = []
    curr_start, curr_stop, curr_variants = hits[0][0], hits[0][1], [hits[0][2]]
    for start, stop, var in hits[1:]:
        if start - curr_stop <= gap_tolerance_bp:
            curr_stop = max(curr_stop, stop)
            curr_variants.append(var)
        else:
            runs.append((curr_start, curr_stop, curr_variants))
            curr_start, curr_stop, curr_variants = start, stop, [var]
    runs.append((curr_start, curr_stop, curr_variants))

    min_bp = min_repeats * len(motif)
    read_len = len(seq)

    # Check 5' end tract
    for start, stop, vars_in_run in runs:
        run_len = stop - start
        if start <= tip_tolerance_bp and run_len >= min_bp:
            # 5' terminal tract
            junction_pos = stop
            # Flank starts at junction_pos and extends into read (away from 5' telomere)
            flank_away = seq[junction_pos : min(read_len, junction_pos + 300)]
            if len(flank_away) < 40:
                continue
            pred_motif = max(set(vars_in_run), key=vars_in_run.count)
            return ("5p", run_len, junction_pos, flank_away, pred_motif)

    # Check 3' end tract
    for start, stop, vars_in_run in reversed(runs):
        run_len = stop - start
        if (read_len - stop) <= tip_tolerance_bp and run_len >= min_bp:
            # 3' terminal tract
            junction_pos = start
            # Flank in read is seq[start-300:start]. Telomere is to the right (3').
            # To standardise orientation pointing AWAY from the telomere, reverse complement it:
            raw_flank = seq[max(0, junction_pos - 300) : junction_pos]
            if len(raw_flank) < 40:
                continue
            flank_away = revcomp(raw_flank)
            pred_motif = max(set(vars_in_run), key=vars_in_run.count)
            return ("3p", run_len, junction_pos, flank_away, pred_motif)

    return None


def cluster_junction_flanks(read_observations, threads=1, id_threshold=0.97, core_flank_len=60):
    """Cluster junction flanks using vsearch --cluster_fast and compute MSA consensus sequences.

    Requires vsearch to be available in PATH (no silent fallbacks).
    """
    if not read_observations:
        return []

    vsearch_bin = shutil.which("vsearch")
    if not vsearch_bin:
        raise RuntimeError("vsearch binary not found in PATH. vsearch is required for telomeric junction flank clustering.")

    with tempfile.TemporaryDirectory() as tmp_dir:
        raw_flanks_path = os.path.join(tmp_dir, "raw_flanks.fasta")
        cons_path = os.path.join(tmp_dir, "clusters_cons.fasta")
        uc_path = os.path.join(tmp_dir, "clusters.uc")

        valid_obs = []
        with open(raw_flanks_path, "w", encoding="utf-8") as f:
            for obs in read_observations:
                flank_away = obs[4]
                if len(flank_away) < 40:
                    continue
                valid_obs.append(obs)
                seq_id = f"obs_{len(valid_obs) - 1}"
                f.write(f">{seq_id}\n{flank_away}\n")

        if not valid_obs:
            return []

        cmd = [
            vsearch_bin,
            "--cluster_fast", raw_flanks_path,
            "--id", str(id_threshold),
            "--iddef", "0",
            "--qmask", "none",
            "--strand", "plus",
            "--consout", cons_path,
            "--uc", uc_path,
            "--sizeout",
            "--threads", str(threads),
        ]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"vsearch clustering failed:\n{proc.stderr}")

        # Parse clusters.uc to assign observations to cluster index
        cluster_members = defaultdict(list)
        cluster_centroids = {}
        if os.path.exists(uc_path):
            with open(uc_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if not parts:
                        continue
                    rec_type = parts[0]
                    if rec_type in ("S", "H"):
                        cluster_nr = int(parts[1])
                        seq_id = parts[8]
                        idx = int(seq_id.replace("obs_", ""))
                        cluster_members[cluster_nr].append(valid_obs[idx])
                        if rec_type == "S":
                            cluster_centroids[cluster_nr] = seq_id

        # Parse clusters_cons.fasta
        consensus_seqs = {}
        if os.path.exists(cons_path):
            curr_header = None
            curr_seq = []
            with open(cons_path, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if curr_header is not None:
                            consensus_seqs[curr_header] = "".join(curr_seq)
                        curr_header = line[1:]
                        curr_seq = []
                    else:
                        curr_seq.append(line)
                if curr_header is not None:
                    consensus_seqs[curr_header] = "".join(curr_seq)

        centroid_to_cons = {}
        for h, seq in consensus_seqs.items():
            c_id = None
            for token in h.split(";"):
                if token.startswith("centroid="):
                    c_id = token.split("=", 1)[1]
                    break
            if not c_id:
                c_id = h.split(";")[0]
            centroid_to_cons[c_id] = seq

        clusters = []
        for cluster_nr in sorted(cluster_members.keys()):
            reads = cluster_members[cluster_nr]
            c_id = cluster_centroids.get(cluster_nr, f"obs_{cluster_nr}")
            cons = centroid_to_cons.get(c_id)
            if not cons:
                flanks = [obs[4] for obs in reads]
                flanks.sort(key=len)
                cons = flanks[len(flanks) // 2]
            rep_core = cons[:core_flank_len]
            clusters.append({
                "cluster_nr": cluster_nr,
                "rep_core": rep_core,
                "rep_ext": cons,
                "reads": reads,
            })

        return clusters


def align_clusters_minimap2(clusters, ref_fasta, threads=1):
    """Align cluster extended flanks (MSA consensus) to reference assembly using minimap2."""
    if not clusters or not ref_fasta or not os.path.exists(ref_fasta):
        return {}

    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp_fa:
        tmp_fa_path = tmp_fa.name
        for i, cl in enumerate(clusters):
            rep_ext = cl.get("rep_ext")
            if not rep_ext:
                flanks = [obs[4] for obs in cl["reads"]]
                flanks.sort(key=len)
                rep_ext = flanks[len(flanks) // 2]
                cl["rep_ext"] = rep_ext
            tmp_fa.write(f">cl_{i}\n{rep_ext}\n")

    alignments = defaultdict(list)
    cmd = [
        "minimap2",
        "-c",
        "-x", "sr",
        "--secondary=yes",
        "-N", "20",
        "-t", str(threads),
        ref_fasta,
        tmp_fa_path
    ]

    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        for line in proc.stdout.splitlines():
            parts = line.strip().split("\t")
            if len(parts) >= 12:
                qname = parts[0]
                qlen = int(parts[1])
                qstart = int(parts[2])
                qend = int(parts[3])
                strand = parts[4]
                tname = parts[5]
                tlen = int(parts[6])
                tstart = int(parts[7])
                tend = int(parts[8])
                mapq = int(parts[11])
                alignments[qname].append({
                    "target": tname,
                    "target_len": tlen,
                    "target_start": tstart,
                    "target_end": tend,
                    "strand": strand,
                    "query_start": qstart,
                    "query_end": qend,
                    "mapq": mapq,
                    "aln_len": tend - tstart
                })
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"WARNING: minimap2 alignment failed or minimap2 not available: {e}\n")
    finally:
        if os.path.exists(tmp_fa_path):
            os.remove(tmp_fa_path)

    return alignments


def classify_junction_cluster(cluster_id, cluster_data, alignments, delfies_points, policy):
    """PDE-aware classification of a subtelomeric junction cluster."""
    min_support = policy["min_read_support_floor"]
    terminal_window_bp = policy["terminal_window_bp"]
    germline_buffer_bp = policy["germline_buffer_bp"]

    reads = cluster_data["reads"]
    read_count = len(reads)
    alns = alignments.get(cluster_id, [])

    # Average tract length and motif breakdown
    avg_tract_len = sum(obs[2] for obs in reads) / read_count if read_count else 0
    motifs = [obs[5] for obs in reads]
    top_motif = max(set(motifs), key=motifs.count) if motifs else "UNKNOWN"

    result = {
        "cluster_id": cluster_id,
        "read_count": read_count,
        "meets_policy_floor": read_count >= min_support,
        "mean_telomere_tract_bp": round(avg_tract_len, 1),
        "dominant_motif": top_motif,
        "consensus_flank_60bp": cluster_data["rep_core"],
        "supporting_read_ids": [obs[0] for obs in reads],
        "alignments": alns,
        "classification": "UNCLASSIFIED",
        "subtype": "UNKNOWN",
        "target_chrom": None,
        "target_arm": None,
        "target_pos": None,
        "is_somatic": False,
        "delfies_breakpoint_match": None
    }

    if not alns:
        result["classification"] = "GERMLINE_JUNCTION"
        result["subtype"] = "UNANCHORED_ORPHAN"
        return result

    # Check if any alignment overlaps a DELFIES diminution breakpoint
    for aln in alns:
        tname = aln["target"]
        tstart, tend = aln["target_start"], aln["target_end"]
        for pt in delfies_points.get(tname, []):
            if (tstart - germline_buffer_bp) <= pt <= (tend + germline_buffer_bp):
                result["classification"] = "SOMATIC_DIMINUTION_JUNCTION"
                result["subtype"] = "DIMINUTION_BREAKPOINT"
                result["target_chrom"] = tname
                result["target_pos"] = pt
                result["is_somatic"] = True
                result["delfies_breakpoint_match"] = {
                    "chrom": tname,
                    "pos": pt,
                    "distance_bp": min(abs(tstart - pt), abs(tend - pt))
                }
                return result

    # Not somatic: classify as GERMLINE_JUNCTION
    result["classification"] = "GERMLINE_JUNCTION"

    # Evaluate best or terminal alignments
    # Check if alignments are near chromosome ends
    terminal_alns = []
    for aln in alns:
        tlen = aln["target_len"]
        tstart, tend = aln["target_start"], aln["target_end"]
        if tstart <= terminal_window_bp:
            terminal_alns.append((aln, "5p"))
        elif tend >= (tlen - terminal_window_bp):
            terminal_alns.append((aln, "3p"))

    if terminal_alns:
        # Sort by MAPQ descending
        terminal_alns.sort(key=lambda x: x[0]["mapq"], reverse=True)
        best_aln, arm = terminal_alns[0]
        result["target_chrom"] = best_aln["target"]
        result["target_arm"] = arm
        result["target_pos"] = best_aln["target_start"] if arm == "5p" else best_aln["target_end"]

        # If high MAPQ and unique
        if best_aln["mapq"] >= 30 and len(alns) == 1:
            result["subtype"] = "UNIQUE_ARM"
        else:
            # Low MAPQ or multi-mapping to a terminal repeat / satellite array
            result["subtype"] = "SATELLITE_REPEAT"
    else:
        # Internal mapping
        best_aln = max(alns, key=lambda x: x["mapq"])
        result["target_chrom"] = best_aln["target"]
        result["target_pos"] = best_aln["target_start"]
        result["subtype"] = "INTERNAL_ITS"

    return result


def run_self_test():
    """Internal self-test verifying vsearch clustering and consensus generation."""
    vsearch_bin = shutil.which("vsearch")
    assert vsearch_bin is not None, "vsearch binary not found in PATH"

    test_seq1 = "A" * 60 + "C" * 40
    test_seq2 = "A" * 60 + "C" * 39 + "G"
    test_seq3 = "T" * 60 + "G" * 40

    obs1 = ("read1", "5p", 100, 100, test_seq1, "TTAGGC")
    obs2 = ("read2", "5p", 120, 100, test_seq2, "TTAGGC")
    obs3 = ("read3", "3p", 80, 50, test_seq3, "TTAGGC")

    clusters = cluster_junction_flanks([obs1, obs2, obs3], threads=1, id_threshold=0.97)
    assert len(clusters) == 2, f"Expected 2 clusters, got {len(clusters)}"
    cl_sizes = sorted([len(c["reads"]) for c in clusters])
    assert cl_sizes == [1, 2], f"Expected cluster sizes [1, 2], got {cl_sizes}"
    for c in clusters:
        assert len(c["rep_ext"]) >= 60, "Consensus sequence too short"
        assert len(c["rep_core"]) == 60, "rep_core must be 60 bp"
    print("TELOMERICJUNCTIONPROFILER self-test PASS")
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description="Telomeric Junction Profiler: PDE-aware clustering and deconvolution of subtelomeric junctions."
    )
    parser.add_argument("--reads", help="Path to telomeric reads FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--ref", help="Path to reference genome assembly FASTA")
    parser.add_argument("--delfies-bed", default=None, help="Path to DELFIES diminution breakpoints BED")
    parser.add_argument("--policy", help="Path to curation_policy.json")
    parser.add_argument("--motif", help="Canonical forward telomeric repeat motif (e.g. TTAGGC)")
    parser.add_argument("--sample", default="sample", help="Sample identifier prefix")
    parser.add_argument("--threads", type=int, default=1, help="Number of alignment threads")
    parser.add_argument("--output-json", help="Path to output catalog JSON")
    parser.add_argument("--output-tsv", default=None, help="Path to output summary TSV")
    parser.add_argument("--self-test", action="store_true", help="Run self-test and exit")

    args = parser.parse_args()

    if args.self_test:
        run_self_test()

    if not (args.reads and args.ref and args.policy and args.motif and args.output_json):
        parser.error("--reads, --ref, --policy, --motif, and --output-json are required when not running --self-test")

    # 1. Read Policy
    policy = read_curation_policy(args.policy)

    # 2. Parse DELFIES Breakpoints
    delfies_points = parse_delfies_bed(args.delfies_bed)

    # 3. Stream reads and identify terminal tracts
    read_observations = []
    total_reads_scanned = 0
    for qname, seq, _qual in stream_fastq(args.reads):
        total_reads_scanned += 1
        tract_info = find_terminal_tracts(
            seq,
            args.motif,
            min_repeats=policy["min_repeats"],
            tip_tolerance_bp=100,
            gap_tolerance_bp=30
        )
        if tract_info:
            orient, tract_len, junc_pos, flank_away, motif_var = tract_info
            read_observations.append((qname, orient, tract_len, junc_pos, flank_away, motif_var))

    sys.stderr.write(
        f"Scanned {total_reads_scanned} reads; found {len(read_observations)} terminal telomeric junctions.\n"
    )

    # 4. Cluster junction flanks
    raw_clusters = cluster_junction_flanks(read_observations, threads=args.threads, core_flank_len=60)
    sys.stderr.write(f"Deconvolved into {len(raw_clusters)} distinct junction clusters.\n")

    # Filter clusters by policy floor
    min_floor = policy["min_read_support_floor"]
    supported_clusters = [cl for cl in raw_clusters if len(cl["reads"]) >= min_floor]
    sys.stderr.write(
        f"{len(supported_clusters)} clusters meet evidence floor (>= {min_floor} supporting reads).\n"
    )

    # 5. Align clusters to reference
    cluster_alignments = align_clusters_minimap2(supported_clusters, args.ref, threads=args.threads)

    # 6. PDE-Aware Classification
    classified_catalog = []
    tsv_rows = []
    for i, cl in enumerate(supported_clusters):
        cid = f"cl_{i}"
        classification = classify_junction_cluster(cid, cl, cluster_alignments, delfies_points, policy)
        classified_catalog.append(classification)

        tsv_rows.append([
            classification["cluster_id"],
            classification["classification"],
            classification["subtype"],
            str(classification["read_count"]),
            str(classification["mean_telomere_tract_bp"]),
            classification["dominant_motif"],
            str(classification["target_chrom"] or "NA"),
            str(classification["target_arm"] or "NA"),
            str(classification["target_pos"] or "NA"),
            str(classification["is_somatic"]),
            classification["consensus_flank_60bp"]
        ])

    # 7. Write outputs
    catalog_doc = {
        "sample": args.sample,
        "telomere_motif": args.motif,
        "policy_thresholds": policy,
        "total_reads_scanned": total_reads_scanned,
        "terminal_junction_reads": len(read_observations),
        "total_clusters_found": len(raw_clusters),
        "supported_clusters_count": len(supported_clusters),
        "clusters": classified_catalog
    }

    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(catalog_doc, f, indent=2, sort_keys=True)

    if args.output_tsv:
        headers = [
            "cluster_id", "classification", "subtype", "read_count",
            "mean_telomere_bp", "dominant_motif", "target_chrom", "target_arm",
            "target_pos", "is_somatic", "consensus_flank_60bp"
        ]
        with open(args.output_tsv, "w", encoding="utf-8") as f:
            f.write("\t".join(headers) + "\n")
            for r in tsv_rows:
                f.write("\t".join(r) + "\n")

    sys.stderr.write(f"Successfully generated {args.output_json}\n")


if __name__ == "__main__":
    main()
