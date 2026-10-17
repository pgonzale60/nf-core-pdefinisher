"""Pure I/O loaders and interval utilities for Module 2 Curation Synthesis.

Speaks strictly the internal schemas produced by Module 1 reconnaissance tools.
British spelling standard enforced across all documentation and comments.
"""

from collections import defaultdict
import gzip
import json
import os
import re
import sys

from .policy import DEFAULT_POLICY


def parse_bnd_target(alt):
    """Extract target chromosome and 1-based coordinate from VCF BND ALT field."""
    if not alt:
        return None, None
    m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def merge_intervals(intervals):
    """Merge overlapping and abutting half-open intervals [s, e)."""
    if not intervals:
        return []
    sorted_ivs = sorted([(int(s), int(e)) for s, e in intervals if e > s], key=lambda x: x[0])
    if not sorted_ivs:
        return []
    merged = [sorted_ivs[0]]
    for s, e in sorted_ivs[1:]:
        prev_s, prev_e = merged[-1]
        if s <= prev_e:
            merged[-1] = (prev_s, max(prev_e, e))
        else:
            merged.append((s, e))
    return merged


def merge_query_spans(records):
    """Compute the non-overlapping cumulative query coverage in bp."""
    if not records:
        return 0
    intervals = [(r["qstart"], r["qend"]) for r in records if r.get("qend", 0) > r.get("qstart", 0)]
    merged = merge_intervals(intervals)
    return sum(e - s for s, e in merged)


def compute_query_represented_intervals(records, min_pid=85.0, min_alen=100):
    """Compute merged non-overlapping query intervals aligned against chromosomes."""
    if not records:
        return []
    intervals = []
    for r in records:
        alen = r.get("alen", 0)
        matches = r.get("matches", 0)
        pid = (matches / alen * 100.0) if alen > 0 else 0.0
        if alen >= min_alen and pid >= min_pid:
            intervals.append((r["qstart"], r["qend"]))
    return merge_intervals(intervals)


def reconstruct_rdna_operons(barrnap_gff, seq_lens, max_spacer_bp=3500, terminal_padding_bp=1000):
    """Reconstruct reasoned rDNA operon models from raw Barrnap subunit annotations."""
    all_raw_units = []
    features_by_seq = defaultdict(list)
    if barrnap_gff and os.path.exists(barrnap_gff):
        with open(barrnap_gff) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.strip().split("\t")
                if len(p) >= 9 and p[2] == "rRNA":
                    c = p[0]
                    st_1 = int(p[3])
                    en_1 = int(p[4])
                    strand = p[6]
                    info = p[8]
                    st_0 = max(0, st_1 - 1)
                    en_0 = en_1
                    features_by_seq[c].append((st_0, en_0, strand, info))
                    all_raw_units.append((c, st_1, en_1, strand, info))

    scaff_rdna_operons = defaultdict(list)
    scaff_rdna_operon_bp = defaultdict(int)

    for c, feats in features_by_seq.items():
        slen = seq_lens.get(c, 0)
        sorted_feats = sorted(feats, key=lambda x: (x[0], x[1]))
        if not sorted_feats:
            continue

        clusters = []
        cur_st, cur_en = sorted_feats[0][0], sorted_feats[0][1]
        for st, en, _strand, _info in sorted_feats[1:]:
            if st <= cur_en + max_spacer_bp:
                cur_en = max(cur_en, en)
            else:
                clusters.append((cur_st, cur_en))
                cur_st, cur_en = st, en
        clusters.append((cur_st, cur_en))

        padded_clusters = []
        for c_st, c_en in clusters:
            p_st = max(0, c_st - terminal_padding_bp)
            p_en = min(slen, c_en + terminal_padding_bp) if slen > 0 else (c_en + terminal_padding_bp)
            padded_clusters.append((p_st, p_en))

        merged_op = merge_intervals(padded_clusters)
        scaff_rdna_operons[c] = merged_op
        scaff_rdna_operon_bp[c] = sum(e - s for s, e in merged_op)

    return scaff_rdna_operons, scaff_rdna_operon_bp, all_raw_units


def parse_unplaced_paf(unplaced_paf):
    """Parse alignments of unplaced contigs against chromosomes."""
    scaff_align = defaultdict(int)
    scaff_primary_records = defaultdict(list)
    scaff_any_records = defaultdict(list)
    if not unplaced_paf or not os.path.exists(unplaced_paf):
        return scaff_align, scaff_primary_records, scaff_any_records
    primary_records = defaultdict(list)
    any_records = defaultdict(list)
    saw_tp_tag = False
    with open(unplaced_paf) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 12:
                continue
            try:
                rec = {
                    "qlen": int(p[1]),
                    "qstart": int(p[2]),
                    "qend": int(p[3]),
                    "strand": p[4],
                    "tname": p[5],
                    "tlen": int(p[6]),
                    "tstart": int(p[7]),
                    "tend": int(p[8]),
                    "matches": int(p[9]) if len(p) > 9 else (int(p[3]) - int(p[2])),
                    "alen": int(p[10]) if len(p) > 10 else (int(p[3]) - int(p[2])),
                    "mapq": int(p[11]),
                }
            except ValueError:
                continue
            qname = p[0]
            any_records[qname].append(rec)
            for tag in p[12:]:
                if tag.startswith("tp:A:"):
                    saw_tp_tag = True
                    if tag[5:] == "P":
                        primary_records[qname].append(rec)
                    break

    primary_spans = defaultdict(int)
    any_spans = defaultdict(int)
    for q, recs in any_records.items():
        any_spans[q] = merge_query_spans(recs)
    if saw_tp_tag:
        for q, recs in primary_records.items():
            primary_spans[q] = merge_query_spans(recs)
        for q in any_records:
            if q not in primary_spans or primary_spans[q] == 0:
                primary_spans[q] = any_spans[q]
            if q not in primary_records or not primary_records[q]:
                primary_records[q] = any_records[q]
        return primary_spans, primary_records, any_records
    return any_spans, any_records, any_records


def parse_telo_paf(telo_paf, nuclear_chroms, policy=None):
    """Parse telomeric unitigs alignment PAF, extracting candidate extensions for chromosome termini.

    Admits an alignment if:
      1. mapq >= mapq_anchored (high-confidence anchor), OR
      2. All alignments of the unitig are exclusive to a single nuclear chromosome (chrom-exclusive
         valency) and reach within terminal_window_bp of a terminus.
    """
    if policy is None:
        policy = DEFAULT_POLICY
    t_win = policy.get("telomere", {}).get("terminal_window_bp", 10000)
    min_ext = policy.get("telomere", {}).get("min_extension_bp", 200)
    mapq_anch = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30)

    unitig_extensions = defaultdict(list)
    if not telo_paf or not os.path.exists(telo_paf):
        return unitig_extensions

    paf_by_query = defaultdict(list)
    with open(telo_paf) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = line.split("\t")
            if len(p) >= 12:
                paf_by_query[p[0]].append(p)

    for qname, records in paf_by_query.items():
        target_chroms = {rec[5] for rec in records}
        is_single_chrom = len(target_chroms) == 1

        best_cand_by_term = {}
        for p in records:
            mapq = int(p[11])
            qlen, qstart, qend = int(p[1]), int(p[2]), int(p[3])
            strand, tname, tlen, tstart, tend = p[4], p[5], int(p[6]), int(p[7]), int(p[8])

            if tname not in nuclear_chroms:
                continue

            is_near_terminal = (tstart <= t_win) or (tend >= tlen - t_win)
            if mapq < mapq_anch and not (is_single_chrom and is_near_terminal):
                continue

            if tstart <= t_win:
                overhang = qstart if strand == "+" else (qlen - qend)
                if overhang >= min_ext:
                    cand = (overhang, qname, 0, 0, qstart, qend, strand, qlen, tstart, tend)
                    key = (tname, "5p")
                    if key not in best_cand_by_term or overhang > best_cand_by_term[key][0]:
                        best_cand_by_term[key] = cand

            if tend >= tlen - t_win:
                overhang = (qlen - qend) if strand == "+" else qstart
                if overhang >= min_ext:
                    cand = (overhang, qname, 0, 0, qstart, qend, strand, qlen, tstart, tend)
                    key = (tname, "3p")
                    if key not in best_cand_by_term or overhang > best_cand_by_term[key][0]:
                        best_cand_by_term[key] = cand

        for key, cand in best_cand_by_term.items():
            unitig_extensions[key].append(cand)

    return unitig_extensions


def load_grs_summary(grs_file):
    """Load GRS/Core segmentation and threshold from identifyputativegrs TSV."""
    chrom_lens = {}
    core_segments = defaultdict(list)
    grs_segments = defaultdict(list)
    grs_arm_depths = {}
    local_depth_lookup = defaultdict(dict)
    grs_rows = []
    parsed_grs_threshold = None

    if os.path.exists(grs_file):
        with open(grs_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    if "GRS_THRESHOLD=" in line:
                        for token in line.split():
                            if token.startswith("GRS_THRESHOLD="):
                                try:
                                    parsed_grs_threshold = float(token.split("=")[1])
                                except ValueError:
                                    pass
                    continue
                parts = line.split("\t")
                if parts[0] == "chromosome":
                    continue
                if len(parts) >= 6:
                    c, s, e, length, cls, depth = (
                        parts[0],
                        int(parts[1]),
                        int(parts[2]),
                        int(parts[3]),
                        parts[4],
                        float(parts[5]),
                    )
                    chrom_lens[c] = max(chrom_lens.get(c, 0), e)
                    grs_rows.append((c, s, e, length, cls, depth))
                    local_depth_lookup[c][(s, e)] = depth
                    if cls == "putative_core":
                        core_segments[c].append((s, e, length, depth))
                    elif cls == "putative_GRS":
                        grs_segments[c].append((s, e, length, depth))
                        if s == 0:
                            grs_arm_depths[(c, "5p")] = (length, depth)

        for c, s, e, length, cls, depth in grs_rows:
            if cls == "putative_GRS" and e == chrom_lens.get(c, 0):
                grs_arm_depths[(c, "3p")] = (length, depth)

    if parsed_grs_threshold is None:
        raise ValueError(
            f"Missing required metadata '# GRS_THRESHOLD=...' in GRS summary file: {grs_file}"
        )

    return (
        chrom_lens,
        core_segments,
        grs_segments,
        grs_arm_depths,
        local_depth_lookup,
        parsed_grs_threshold,
    )


def load_mosdepth_tip_depths(mosdepth_regions, chrom_lens):
    """Extract immediate 1 kb terminal window depths from mosdepth regions.bed.gz."""
    tip_depths_1kb = {}
    if mosdepth_regions and os.path.exists(mosdepth_regions):
        with gzip.open(mosdepth_regions, "rt") as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 4:
                    c, s, e, d = p[0], int(p[1]), int(p[2]), float(p[3])
                    if s == 0:
                        tip_depths_1kb[(c, "5p")] = d
                    if c in chrom_lens and e >= chrom_lens[c] - 1000:
                        tip_depths_1kb[(c, "3p")] = d
    return tip_depths_1kb


def partition_chromosomes(chrom_lens, chromosomes_arg=None, unplaced_paf=None):
    """Partition assembly sequences into canonical nuclear chromosomes and unplaced scaffolds."""
    sorted_chroms = sorted(
        chrom_lens.keys(),
        key=lambda x: (not x.startswith("I") and not x.startswith("V") and not x.startswith("X"), x),
    )
    if chromosomes_arg:
        specified = [c.strip() for c in chromosomes_arg.split(",") if c.strip()]
        nuclear_chroms = [c for c in specified if c in chrom_lens]
        unplaced_chroms = [c for c in sorted_chroms if c not in nuclear_chroms]
        return sorted_chroms, nuclear_chroms, unplaced_chroms

    paf_targets = set()
    paf_queries = set()
    if unplaced_paf and os.path.exists(unplaced_paf):
        try:
            with open(unplaced_paf) as pf:
                for line in pf:
                    parts = line.split("\t")
                    if len(parts) >= 6:
                        paf_queries.add(parts[0])
                        paf_targets.add(parts[5])
        except Exception:
            pass

    nuclear_chroms = []
    unplaced_chroms = []
    has_megabase_chroms = len([x for x in sorted_chroms if chrom_lens.get(x, 0) > 1000000]) >= 3
    for c in sorted_chroms:
        c_len = chrom_lens.get(c, 0)
        is_unplaced_name = (
            c.startswith("scaffold")
            or c.startswith("contig")
            or "MT" in c
            or "unloc" in c.lower()
            or c in paf_queries
        )
        if c in paf_targets and not is_unplaced_name:
            nuclear_chroms.append(c)
        elif is_unplaced_name:
            unplaced_chroms.append(c)
        elif re.match(r"^[A-Z]{4,6}\d{8,}", c):
            unplaced_chroms.append(c)
        elif c_len < 500000 and has_megabase_chroms:
            unplaced_chroms.append(c)
        else:
            nuclear_chroms.append(c)

    return sorted_chroms, nuclear_chroms, unplaced_chroms


def load_rdna_bed(rdna_bed, nuclear_chroms):
    """Load rDNA intervals from rdna_regions.bed."""
    rdna_intervals = defaultdict(list)
    rdna_major_regions = []
    if rdna_bed and os.path.exists(rdna_bed):
        with open(rdna_bed) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = line.split("\t")
                if len(p) >= 3:
                    c, s, en = p[0], int(p[1]), int(p[2])
                    rdna_intervals[c].append((s, en))
                    if c in nuclear_chroms:
                        rdna_major_regions.append((c, s, en))
    return rdna_intervals, rdna_major_regions


def load_terminal_summary(terminal_summary):
    """Load baseline terminal status and measured terminal_tract_bp from softclipends."""
    baseline_status = {}
    tip_tract_bp = {}
    if terminal_summary and os.path.exists(terminal_summary):
        with open(terminal_summary) as f:
            f.readline()  # advance past header
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 6:
                    cend, _clen, _n_clips, _n_telo, _max_tc, p_stat = p[:6]
                    if "_" in cend:
                        c, e = cend.split("_", 1)
                        baseline_status[(c, e)] = p_stat
                        if len(p) >= 7:
                            try:
                                tip_tract_bp[(c, e)] = int(float(p[6]))
                            except ValueError:
                                pass
    return baseline_status, tip_tract_bp


def load_terminal_softclips(terminal_softclips, policy=None):
    """Load per-read terminal soft-clip records and select best single-read extension per end."""
    if policy is None:
        policy = DEFAULT_POLICY
    min_t_reps = policy.get("telomere", {}).get("min_repeats", 15)

    softclip_candidates = defaultdict(list)
    raw_softclips = defaultdict(list)
    if terminal_softclips and os.path.exists(terminal_softclips):
        with open(terminal_softclips) as f:
            for line in f:
                if line.startswith("query_name") or line.startswith("qname") or "overhang_bp" in line:
                    continue
                p = line.strip().split("\t")
                if len(p) >= 9:
                    qname, ref, end_type = p[0], p[1], p[2]
                    try:
                        align_start = int(p[3])
                        align_end = int(p[4])
                        clip_len = int(p[5])
                        overhang_bp = int(p[6])
                        telo_cnt = int(p[7])
                        mapq = int(p[8])
                    except ValueError:
                        continue
                    rec = {
                        "qname": qname,
                        "align_start": align_start,
                        "align_end": align_end,
                        "clip_len": clip_len,
                        "overhang_bp": overhang_bp,
                        "telo_cnt": telo_cnt,
                        "mapq": mapq,
                    }
                    raw_softclips[(ref, end_type)].append(rec)
                    if mapq >= 20 and telo_cnt >= min_t_reps and overhang_bp >= 50:
                        softclip_candidates[(ref, end_type)].append(
                            (overhang_bp, telo_cnt, qname, mapq)
                        )

    single_read_ext = {}
    for (c, e), cand_list in softclip_candidates.items():
        if cand_list:
            cand_list.sort(key=lambda x: (x[0] // 10, x[1], x[0]), reverse=True)
            best = cand_list[0]
            single_read_ext[(c, e)] = (best[0], len(cand_list), best[2], best[1])

    return raw_softclips, single_read_ext


def load_delfies_points(delfies_bed, nuclear_chroms):
    """Load somatic diminution breakpoint coordinates per nuclear chromosome."""
    delfies_points = defaultdict(list)
    if delfies_bed and os.path.exists(delfies_bed):
        with open(delfies_bed) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 3:
                    c, s = p[0], int(p[1])
                    if c in nuclear_chroms:
                        delfies_points[c].append(s)
    return delfies_points


def load_telo_junctions(telo_junctions_path):
    """Load telomeric junction clusters JSON from TELOMERICJUNCTIONPROFILER."""
    if telo_junctions_path and os.path.exists(telo_junctions_path):
        try:
            with open(telo_junctions_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            sys.stderr.write(f"Warning: Failed to load telo-junctions JSON {telo_junctions_path}: {e}\n")
    return None


def load_flye_info(flye_info_path):
    """Parse Flye repeat-graph assembly_info.txt for circular rDNA monomer resolution."""
    flye_rows = []
    circ_monomer_name = None
    circ_monomer_len = 0
    circ_cov = "0"
    anchor_name = None

    if flye_info_path and os.path.exists(flye_info_path):
        with open(flye_info_path) as f:
            f.readline()  # header
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 6:
                    cname, clen, ccov, ccirc, crep, cmult = (
                        p[0],
                        int(p[1]),
                        p[2],
                        p[3],
                        p[4],
                        p[5],
                    )
                    flye_rows.append((cname, clen, ccov, ccirc, crep, cmult))
                    if ccirc == "Y":
                        circ_monomer_name = cname
                        circ_monomer_len = clen
                        circ_cov = ccov
                    elif not anchor_name and clen > 20000:
                        anchor_name = cname

    return {
        "rows": flye_rows,
        "circ_monomer_name": circ_monomer_name,
        "circ_monomer_len": circ_monomer_len,
        "circ_cov": circ_cov,
        "anchor_name": anchor_name,
    }


def load_red_repeats(red_bed):
    """Compute cumulative repeat-masked bp per sequence from Red BED."""
    scaff_repeats = defaultdict(int)
    if red_bed and os.path.exists(red_bed):
        with open(red_bed) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 3:
                    sc = p[0]
                    span = int(p[2]) - int(p[1])
                    scaff_repeats[sc] += span
    return scaff_repeats
