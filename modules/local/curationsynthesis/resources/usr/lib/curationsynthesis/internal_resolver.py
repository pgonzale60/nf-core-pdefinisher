"""Internal Gap, Structural Variant, and Major rDNA Locus Resolution Engine (Domains 1 & 3).

Synthesises deterministic action cards for internal chromosomal gaps, BND/INV scaffolding
misjoins, novel insertions, collapsed tandem duplications, and internal inverted rDNA arrays.

British spelling standard enforced across all documentation and comments.
"""

from collections import defaultdict
from .evidence_io import parse_bnd_target
from .policy import DEFAULT_POLICY


def _index_bnd_svs(all_svs):
    """Pre-index BND and INV structural variants by source and target chromosome."""
    bnd_svs_by_chrom = defaultdict(list)
    bnd_svs_by_target = defaultdict(list)
    for s in all_svs:
        if s.get("svtype") in ("INV", "BND"):
            bnd_svs_by_chrom[s["chrom"]].append(s)
            s_chr2 = s.get("chr2")
            t_chr, _ = parse_bnd_target(s.get("alt", ""))
            tgt = s_chr2 or t_chr
            if tgt:
                bnd_svs_by_target[tgt].append(s)
    return bnd_svs_by_chrom, bnd_svs_by_target


def _measure_gap_reads(gap, glen, mapq_anch, mapq_floor, min_ins_bp):
    """Aggregate per-read micro-observations across an internal assembly gap."""
    read_obs = gap.get("read_observations") or []
    span_cnt = 0
    span_cnt_30 = 0
    del_ops_30 = []
    ins_ops_30 = []
    left_cnt = 0
    left_cnt_20 = 0
    right_cnt = 0
    right_cnt_20 = 0
    for r in read_obs:
        stype = r.get("span_type")
        mapq = r.get("mapq", 0)
        if stype == "COLLINEAR_SPANNING":
            span_cnt += 1
            if mapq >= mapq_anch:
                span_cnt_30 += 1
                del_op = r.get("deletion_across_gap", 0)
                if abs(del_op - glen) <= 25:
                    del_ops_30.append(del_op)
                ins_op = r.get("insertion_across_gap", 0)
                if ins_op >= min_ins_bp:
                    ins_ops_30.append(ins_op)
        elif stype == "LEFT_ABUTTING":
            left_cnt += 1
            if mapq >= mapq_floor:
                left_cnt_20 += 1
        elif stype == "RIGHT_ABUTTING":
            right_cnt += 1
            if mapq >= mapq_floor:
                right_cnt_20 += 1

    return {
        "span_cnt": span_cnt,
        "span_cnt_30": span_cnt_30,
        "span_del_cnt": len(del_ops_30),
        "ins_bp": max(ins_ops_30) if ins_ops_30 else 0,
        "left_cnt": left_cnt,
        "left_cnt_20": left_cnt_20,
        "right_cnt": right_cnt,
        "right_cnt_20": right_cnt_20,
    }


def _find_best_bnd_candidate(
    chrom, start, end, d_local, bnd_svs_by_chrom, bnd_svs_by_target, gap_link_dist, min_floor
):
    """Identify the highest-scoring BND/INV candidate associated with an unspanned gap."""
    seen_sv_ids = set()
    relevant_svs = []
    for s in bnd_svs_by_chrom.get(chrom, []) + bnd_svs_by_target.get(chrom, []):
        sid = s.get("id")
        if sid not in seen_sv_ids:
            seen_sv_ids.add(sid)
            relevant_svs.append(s)

    bnd_candidates_by_target = defaultdict(list)
    for s in relevant_svs:
        if s.get("svtype") not in ("INV", "BND"):
            continue
        s_chrom = s["chrom"]
        s_start = s["start"]
        s_end = s.get("end", s_start)
        s_chr2 = s.get("chr2")
        s_alt = s.get("alt", "")
        t_chr, t_pos = parse_bnd_target(s_alt)
        if not s_chr2 or s_chr2 == s_chrom:
            s_chr2 = t_chr
        target_coord = t_pos if t_pos is not None else s_end

        if s_chrom == chrom:
            target = s_chr2 if s_chr2 else chrom
            if target == chrom:
                d_s = min(abs(s_start - start), abs(s_start - end))
                d_t = min(abs(target_coord - start), abs(target_coord - end))
                dist = min(d_s, d_t)
                if d_s <= d_t:
                    near_coord = s_start
                    far_coord = target_coord
                else:
                    near_coord = target_coord
                    far_coord = s_start
            else:
                dist = min(abs(s_start - start), abs(s_start - end), abs(s_end - start), abs(s_end - end))
                near_coord = s_start
                far_coord = target_coord
            if dist <= gap_link_dist:
                bnd_candidates_by_target[target].append(
                    {
                        "direction": "forward",
                        "sv": s,
                        "dist": dist,
                        "breakend_coord": near_coord,
                        "target_chrom": target,
                        "target_coord": far_coord,
                        "support": s.get("support", 0),
                        "vaf": s.get("vaf", 0.0),
                        "alt": s_alt,
                    }
                )
        elif s_chr2 == chrom or t_chr == chrom:
            src_pos = target_coord if (t_chr == chrom or s_chr2 == chrom) else s_end
            dist = min(abs(src_pos - start), abs(src_pos - end))
            if dist <= gap_link_dist:
                target = s_chrom
                bnd_candidates_by_target[target].append(
                    {
                        "direction": "reciprocal",
                        "sv": s,
                        "dist": dist,
                        "breakend_coord": src_pos,
                        "target_chrom": target,
                        "target_coord": s_start,
                        "support": s.get("support", 0),
                        "vaf": s.get("vaf", 0.0),
                        "alt": s_alt,
                    }
                )

    for _target, bnd_list in list(bnd_candidates_by_target.items()):
        for b in list(bnd_list):
            split_obs = b["sv"].get("split_read_observations")
            if split_obs:
                sa_targets = defaultdict(int)
                sa_coords = {}
                for r in split_obs:
                    sa = r.get("sa_tag", "")
                    for entry in sa.split(";"):
                        if entry.strip():
                            parts = entry.split(",")
                            if len(parts) >= 3 and parts[0] != chrom:
                                sa_chr = parts[0]
                                try:
                                    sa_pos = int(parts[1])
                                    sa_targets[sa_chr] += 1
                                    sa_coords[sa_chr] = sa_pos
                                except ValueError:
                                    pass
                for sa_chr, count in sa_targets.items():
                    if count >= min_floor:
                        bnd_candidates_by_target[sa_chr].append(
                            {
                                "direction": "forward_sa",
                                "sv": b["sv"],
                                "dist": b["dist"],
                                "breakend_coord": b["breakend_coord"],
                                "target_chrom": sa_chr,
                                "target_coord": sa_coords.get(sa_chr, 0),
                                "support": count,
                                "vaf": count / d_local,
                                "alt": f"[{sa_chr}:{sa_coords.get(sa_chr, 0)}[",
                            }
                        )

    scored_candidates = []
    for tgt, b_list in bnd_candidates_by_target.items():
        if not tgt:
            continue
        tot_supp = sum(b["support"] for b in b_list)
        max_vaf = max(b["vaf"] for b in b_list)
        has_fwd = any("forward" in b["direction"] for b in b_list)
        has_recip = any(b["direction"] == "reciprocal" for b in b_list)
        recip_confirmed = has_fwd and has_recip
        score = tot_supp * 3.0 if recip_confirmed else tot_supp
        rep_bnd = max(b_list, key=lambda x: x["support"])
        scored_candidates.append(
            {
                "score": score,
                "target_chrom": tgt,
                "target_coord": rep_bnd["target_coord"],
                "breakend_coord": rep_bnd["breakend_coord"],
                "total_support": tot_supp,
                "max_vaf": max_vaf,
                "reciprocal_confirmed": recip_confirmed,
                "alt": rep_bnd["alt"],
                "sv": rep_bnd["sv"],
            }
        )

    if scored_candidates:
        scored_candidates.sort(key=lambda x: x["score"], reverse=True)
        best_cand = scored_candidates[0]
        if best_cand["total_support"] >= min_floor:
            return best_cand
    return None


def synthesize_curation_cards(
    evidence_data,
    policy=None,
    local_depth_lookup=None,
    genome_core_depth=100.0,
    chrom_lens=None,
    rdna_intervals=None,
):
    """Synthesise Domain 1 internal gap and structural variant action cards."""
    if policy is None:
        policy = DEFAULT_POLICY

    min_floor = policy.get("evidence_floors", {}).get("min_read_support_floor", 3)
    min_support_fraction = policy.get("evidence_floors", {}).get("min_support_fraction", 0.15)
    min_ins_vaf = policy.get("evidence_floors", {}).get("insertion_min_vaf", 0.70)
    min_ins_bp = policy.get("evidence_floors", {}).get("min_insertion_bp", 50)
    min_sv_len = policy.get("evidence_floors", {}).get("structural_variant_min_len_bp", 500)
    mapq_anch = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30)
    mapq_floor = policy.get("mapping_and_anchors", {}).get("mapq_evidence_floor", 20)
    gap_link_dist = policy.get("genomic_architecture", {}).get("gap_sv_link_dist_bp", 1500)
    t_win = policy.get("telomere", {}).get("terminal_window_bp", 10000)

    cards = []
    card_idx = 1

    raw_gaps = evidence_data.get("gaps", [])
    all_svs = evidence_data.get("structural_variants", [])
    bnd_svs_by_chrom, bnd_svs_by_target = _index_bnd_svs(all_svs)
    handled_sv_ids = set()

    def get_local_depth(chrom, pos):
        if local_depth_lookup and chrom in local_depth_lookup:
            for (s, e), d in local_depth_lookup[chrom].items():
                if s <= pos <= e:
                    return max(d, 1.0)
        return max(genome_core_depth, 1.0)

    seen_gap_coords = set()
    gaps = []
    for g in raw_gaps:
        key = (g.get("chrom"), g.get("start"), g.get("end"))
        if key not in seen_gap_coords:
            seen_gap_coords.add(key)
            gaps.append(g)

    for g in gaps:
        chrom = g["chrom"]
        start = g["start"]
        end = g["end"]
        glen = g.get("gap_length", g.get("length", end - start))
        d_local = get_local_depth(chrom, start)

        # 1. rDNA Array Locus Check
        if rdna_intervals and chrom in rdna_intervals:
            matching_rdna_span = next(
                ((r_st, r_en) for r_st, r_en in rdna_intervals[chrom] if max(start, r_st) < min(end, r_en)),
                None,
            )
            if matching_rdna_span:
                cards.append(
                    {
                        "card_id": f"CARD-{card_idx:02d}",
                        "domain": "DOMAIN_INTERNAL",
                        "category": "Assembly Gap: Subsumed by rDNA Array",
                        "action": "SUBSUMED_BY_RDNA_ARRAY",
                        "title": f"Gap within rDNA Array at {chrom}:{start:,}-{end:,} ({glen} bp)",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "length": glen,
                        "confidence": "HIGH",
                        "evidence": {
                            "locus": f"{chrom}:{start}-{end}",
                            "rdna_array": f"{chrom}:{matching_rdna_span[0]}-{matching_rdna_span[1]}",
                        },
                        "instructions": [
                            f"Gap at {chrom}:{start:,}-{end:,} is located within a multi-copy rDNA array ({chrom}:{matching_rdna_span[0]:,}-{matching_rdna_span[1]:,}).",
                            "Omit from local targeted gap assembly; the entire array is replaced by the Module 4 tripartite rDNA cassette.",
                        ],
                    }
                )
                card_idx += 1
                continue

        m = _measure_gap_reads(g, glen, mapq_anch, mapq_floor, min_ins_bp)
        span_cnt = m["span_cnt"]
        span_cnt_30 = m["span_cnt_30"]
        span_del_cnt = m["span_del_cnt"]
        ins_bp = m["ins_bp"]
        left_cnt = m["left_cnt"]
        left_cnt_20 = m["left_cnt_20"]
        right_cnt = m["right_cnt"]
        right_cnt_20 = m["right_cnt_20"]

        collinear_confirmed = span_cnt_30 >= min_floor
        best_bnd = None
        if not collinear_confirmed:
            best_bnd = _find_best_bnd_candidate(
                chrom, start, end, d_local, bnd_svs_by_chrom, bnd_svs_by_target, gap_link_dist, min_floor
            )
            if best_bnd:
                handled_sv_ids.add(best_bnd["sv"]["id"])

        # 2. Trans-chromosomal Misjoin or Intra-chromosomal Inversion
        if best_bnd:
            target_chr = best_bnd["target_chrom"]
            target_pos = best_bnd["target_coord"]
            near_pos = best_bnd["breakend_coord"]
            supp_all = best_bnd["total_support"]
            max_vaf = best_bnd["max_vaf"]
            recip_conf = best_bnd["reciprocal_confirmed"]
            alt = best_bnd["alt"]
            sv = best_bnd["sv"]

            if target_chr != chrom:
                cards.append(
                    {
                        "card_id": f"CARD-{card_idx:02d}",
                        "domain": "DOMAIN_INTERNAL",
                        "category": "Scaffolding Error: Trans-chromosomal Misjoin",
                        "action": "TRANSLOCATE_CONTIG_CANDIDATE",
                        "title": f"Trans-chromosomal Link Candidate: {chrom}:{start:,}-{end:,} -> {target_chr}:{target_pos:,}",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "gap_start": start,
                        "gap_end": end,
                        "breakend_coord": near_pos,
                        "target_chrom": target_chr,
                        "target_coord": target_pos,
                        "alt": alt,
                        "confidence": "HIGH" if recip_conf or max_vaf >= 0.30 else "MEDIUM",
                        "evidence": {
                            "sniffles_id": sv["id"],
                            "svtype": sv.get("svtype", "BND"),
                            "chr2": target_chr,
                            "target_pos": target_pos,
                            "reciprocal_confirmed": recip_conf,
                            "split_reads_supporting": f"{supp_all} total (Reciprocal: {recip_conf})",
                            "local_depth": f"{d_local:.1f}× (VAF: {max_vaf * 100:.1f}%)",
                            "associated_gap": f"{chrom}:{start:,}-{end:,} ({glen} Ns)",
                        },
                        "instructions": [
                            "Delegate to Module 3 Engine D (TRANSCHROMOSOMAL_LINK_RESOLVER)",
                            f"Evaluate negative control (zero collinear reads across {chrom}:{start}-{end})",
                            f"Evaluate positive control (reciprocal split-reads connecting to {target_chr}:{target_pos})",
                            f"Check target locus integrity: if {target_chr}:{target_pos} is unbroken, sever contig to unplaced",
                        ],
                    }
                )
                card_idx += 1
                continue
            else:
                far_pos = target_pos
                if far_pos < start:
                    inv_block_start = far_pos
                    inv_block_end = start
                    card_start = far_pos
                    card_end = end
                elif far_pos > end:
                    inv_block_start = end
                    inv_block_end = far_pos
                    card_start = start
                    card_end = far_pos
                else:
                    inv_block_start = min(near_pos, far_pos)
                    inv_block_end = max(near_pos, far_pos)
                    card_start = min(inv_block_start, start)
                    card_end = max(inv_block_end, end)
                inv_len = inv_block_end - inv_block_start
                if inv_len >= 500:
                    cards.append(
                        {
                            "card_id": f"CARD-{card_idx:02d}",
                            "domain": "DOMAIN_INTERNAL",
                            "category": "Scaffolding Error: Inverted Contig",
                            "action": "INVERT_AND_STITCH",
                            "title": f"Inverted Contig Block: {chrom}:{inv_block_start:,}-{inv_block_end:,}",
                            "chrom": chrom,
                            "start": card_start,
                            "end": card_end,
                            "inv_start": inv_block_start,
                            "inv_end": inv_block_end,
                            "gap_start": start,
                            "gap_end": end,
                            "length": inv_len,
                            "confidence": "HIGH" if max_vaf >= 0.30 else "MEDIUM",
                            "evidence": {
                                "sniffles_id": sv["id"],
                                "svlen": f"{sv.get('svlen', inv_len):,} bp",
                                "inverted_contig_span": f"{chrom}:{inv_block_start:,}-{inv_block_end:,} ({inv_len:,} bp)",
                                "split_reads_supporting": f"{supp_all} total",
                                "local_depth": f"{d_local:.1f}× (VAF: {max_vaf * 100:.1f}%)",
                                "associated_gap": f"{chrom}:{start:,}-{end:,} ({glen} Ns)",
                            },
                            "instructions": [
                                f"Slice left collinear flank: 1..{card_start - 1}",
                                f"Reverse-complement inverted contig block: {inv_block_start}..{inv_block_end}",
                                f"Excise dummy scaffolding gap: {start}..{end} ({glen} Ns)",
                                f"Stitch directly to right collinear flank: {card_end}..end",
                            ],
                        }
                    )
                    card_idx += 1
                    continue

        # 3. Subtelomeric False-Tail Scaffolder Gap
        c_len = chrom_lens.get(chrom, 0) if chrom_lens else 0
        if c_len > 0:
            dist_to_3p = c_len - end
            if dist_to_3p <= t_win and span_cnt_30 <= 1 and right_cnt <= 1 and left_cnt >= min_floor:
                d_3p = g.get("depth_flank_3p", right_cnt)
                if d_3p <= min_floor:
                    cards.append(
                        {
                            "card_id": f"CARD-{card_idx:02d}",
                            "domain": "DOMAIN_INTERNAL",
                            "category": "Assembly Gap: False Terminal Scaffolding Tail",
                            "action": "EXCISE_DUMMY_GAP",
                            "title": f"Terminal Scaffolding Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                            "chrom": chrom,
                            "start": start,
                            "end": c_len,
                            "length": c_len - start,
                            "confidence": "HIGH",
                            "evidence": {
                                "collinear_spanning_reads": f"{span_cnt} (MAPQ>=30: {span_cnt_30})",
                                "left_boundary_reads": f"{left_cnt} (MAPQ>=20: {left_cnt_20})",
                                "right_boundary_reads": f"{right_cnt} (MAPQ>=20: {right_cnt_20})",
                                "terminal_tail_length_bp": dist_to_3p,
                                "distal_flank_depth": d_3p,
                            },
                            "instructions": [
                                f"Excise {c_len - start} bp comprising {glen} Ns dummy gap and {dist_to_3p} bp unanchored false terminal scaffolding tail from {chrom}:{start:,} to {c_len:,}.",
                                f"Expose authentic upstream chromosome boundary at {chrom}:{start:,} to telomere walk/capping.",
                            ],
                        }
                    )
                    card_idx += 1
                    continue

            dist_to_5p = start
            if dist_to_5p <= t_win and span_cnt_30 <= 1 and left_cnt <= 1 and right_cnt >= min_floor:
                d_5p = g.get("depth_flank_5p", left_cnt)
                if d_5p <= min_floor:
                    cards.append(
                        {
                            "card_id": f"CARD-{card_idx:02d}",
                            "domain": "DOMAIN_INTERNAL",
                            "category": "Assembly Gap: False Terminal Scaffolding Tail",
                            "action": "EXCISE_DUMMY_GAP",
                            "title": f"Terminal Scaffolding Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                            "chrom": chrom,
                            "start": 0,
                            "end": end,
                            "length": end,
                            "confidence": "HIGH",
                            "evidence": {
                                "collinear_spanning_reads": f"{span_cnt} (MAPQ>=30: {span_cnt_30})",
                                "left_boundary_reads": f"{left_cnt} (MAPQ>=20: {left_cnt_20})",
                                "right_boundary_reads": f"{right_cnt} (MAPQ>=20: {right_cnt_20})",
                                "terminal_tail_length_bp": dist_to_5p,
                                "distal_flank_depth": d_5p,
                            },
                            "instructions": [
                                f"Excise {end} bp comprising {dist_to_5p} bp unanchored false terminal scaffolding tail and {glen} Ns dummy gap from {chrom}:0 to {end:,}.",
                                f"Expose authentic downstream chromosome boundary at {chrom}:{end:,} to telomere walk/capping.",
                            ],
                        }
                    )
                    card_idx += 1
                    continue

        # 4. Standard Assembly Gap Resolution
        span_frac = span_cnt_30 / d_local
        if span_cnt_30 >= min_floor and span_del_cnt >= min_floor and (span_del_cnt / d_local) >= min_support_fraction:
            cards.append(
                {
                    "card_id": f"CARD-{card_idx:02d}",
                    "domain": "DOMAIN_INTERNAL",
                    "category": "Assembly Gap: Spurious Scaffolder Gap",
                    "action": "EXCISE_DUMMY_GAP",
                    "title": f"Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": glen,
                    "confidence": "HIGH" if (span_del_cnt / d_local) >= 0.30 else "MEDIUM",
                    "evidence": {
                        "collinear_spanning_reads": f"{span_cnt} total (MAPQ>=30: {span_cnt_30})",
                        "spanning_deletions_matching_gap": f"{span_del_cnt} reads",
                        "local_depth": f"{d_local:.1f}× (Support Ratio: {span_frac * 100:.1f}%)",
                        "left_boundary_reads": f"{left_cnt} (MAPQ>=20: {left_cnt_20})",
                        "right_boundary_reads": f"{right_cnt} (MAPQ>=20: {right_cnt_20})",
                    },
                    "instructions": [
                        f"Excise {glen} Ns and stitch {chrom}:{start:,} directly to {chrom}:{end:,} (confirmed by {span_del_cnt} high-quality spanning reads deleting the dummy gap)"
                    ],
                }
            )
        elif span_cnt_30 >= min_floor and ins_bp > min_ins_bp and span_frac >= min_support_fraction:
            cards.append(
                {
                    "card_id": f"CARD-{card_idx:02d}",
                    "domain": "DOMAIN_INTERNAL",
                    "category": "Assembly Gap: Spanned by Reads",
                    "action": "PATCH_GAP_INSERTION",
                    "title": f"Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": glen,
                    "confidence": "HIGH",
                    "evidence": {
                        "collinear_spanning_reads": f"{span_cnt} total (MAPQ>=30: {span_cnt_30})",
                        "max_spanning_insertion": f"{ins_bp} bp",
                        "local_depth": f"{d_local:.1f}× (Support Ratio: {span_frac * 100:.1f}%)",
                        "left_boundary_reads": f"{left_cnt}",
                        "right_boundary_reads": f"{right_cnt}",
                    },
                    "instructions": [
                        f"Replace {glen} Ns with consensus sequence ({ins_bp} bp) derived from {span_cnt_30} high-quality spanning reads"
                    ],
                }
            )
        elif span_cnt_30 >= min_floor and span_frac >= min_support_fraction:
            cards.append(
                {
                    "card_id": f"CARD-{card_idx:02d}",
                    "domain": "DOMAIN_INTERNAL",
                    "category": "Assembly Gap: Spanned by Reads",
                    "action": "CLOSE_GAP_COLINEAR",
                    "title": f"Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": glen,
                    "confidence": "HIGH",
                    "evidence": {
                        "collinear_spanning_reads": f"{span_cnt} total (MAPQ>=30: {span_cnt_30})",
                        "local_depth": f"{d_local:.1f}× (Support Ratio: {span_frac * 100:.1f}%)",
                        "left_boundary_reads": f"{left_cnt}",
                        "right_boundary_reads": f"{right_cnt}",
                    },
                    "instructions": [
                        f"Gap bridged collinearly by {span_cnt_30} high-quality reads. Stitch flanks directly."
                    ],
                }
            )
        else:
            cards.append(
                {
                    "card_id": f"CARD-{card_idx:02d}",
                    "domain": "DOMAIN_INTERNAL",
                    "category": "Assembly Gap: Unresolved Repeat",
                    "action": "REQUIRES_UNITIG_PATCH",
                    "title": f"Gap at {chrom}:{start:,}-{end:,} ({glen} bp)",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": glen,
                    "confidence": "HIGH" if (left_cnt >= min_floor or right_cnt >= min_floor) else "MEDIUM",
                    "evidence": {
                        "collinear_spanning_reads": f"{span_cnt} (MAPQ>=30: {span_cnt_30})",
                        "local_depth": f"{d_local:.1f}×",
                        "left_boundary_reads": f"{left_cnt} (MAPQ>=20: {left_cnt_20})",
                        "right_boundary_reads": f"{right_cnt} (MAPQ>=20: {right_cnt_20})",
                    },
                    "instructions": [
                        f"Gap flanked by {left_cnt} left-abutting reads and {right_cnt} right-abutting reads. Missing sequence exceeds single-read linear alignment; patch using local unitig or Flye repeat graph."
                    ],
                }
            )
        card_idx += 1

    # 5. Macro Novel Sequence Insertions
    for s in all_svs:
        if s["id"] in handled_sv_ids:
            continue
        d_local = get_local_depth(s["chrom"], s["start"])
        supp = s.get("support", 0)

        if (
            s["svtype"] == "INS"
            and s["svlen"] >= min_sv_len
            and supp >= min_floor
            and s.get("vaf", 0) >= min_ins_vaf
        ):
            cards.append(
                {
                    "card_id": f"CARD-{card_idx:02d}",
                    "domain": "DOMAIN_INTERNAL",
                    "category": "Missing Genomic Sequence: Novel Insertion",
                    "action": "INSERT_SEQUENCE",
                    "title": f"Unassembled Sequence at {s['chrom']}:{s['start']:,} ({s['svlen']:,} bp)",
                    "chrom": s["chrom"],
                    "start": s["start"],
                    "end": s["start"],
                    "length": s["svlen"],
                    "confidence": "HIGH" if (s.get("vaf", 0) >= 0.9 and supp >= min_floor) else "MEDIUM",
                    "evidence": {
                        "sniffles_id": s["id"],
                        "svlen": f"{s['svlen']:,} bp",
                        "max_spanning_insertion": f"{s['svlen']} bp",
                        "supporting_reads": f"{supp} reads (Local Depth: {d_local:.1f}×, VAF: {s.get('vaf', 0) * 100:.1f}%)",
                        "reconstructed_sequence_available": "YES" if s.get("has_explicit_alt_seq") else "NO",
                    },
                    "instructions": [
                        f"Insert {s['svlen']:,} bp novel sequence into {s['chrom']} at position {s['start']:,}",
                        f"Sequence source: Sniffles ALT consensus ({s.get('alt_seq_len', 0)} bp)"
                        if s.get("has_explicit_alt_seq")
                        else "Extract spanning HiFi reads to generate insertion consensus",
                    ],
                }
            )
            card_idx += 1

    return cards


def synthesize_rdna_cards(rdna_major_regions, chrom_lens, all_r_units, flye_data, policy=None):
    """Synthesise Domain 3 action cards and locus audit records for major rDNA arrays."""
    if policy is None:
        policy = DEFAULT_POLICY
    rdna_dist_bp = policy.get("repeat_and_rdna", {}).get("rdna_dist_from_end", 25000)

    circ_monomer_name = flye_data.get("circ_monomer_name")
    circ_monomer_len = flye_data.get("circ_monomer_len", 0)
    circ_cov = flye_data.get("circ_cov", "0")

    rdna_cards = []
    locus_audits = []
    rdna_card_idx = 1

    for r_chrom, r_st, r_en in rdna_major_regions:
        c_len = chrom_lens.get(r_chrom, 0)
        is_term = (r_st <= rdna_dist_bp) or (r_en >= c_len - rdna_dist_bp)
        r_units_region = [
            u for u in all_r_units if u[0] == r_chrom and u[1] >= r_st - 5000 and u[2] <= r_en + 5000
        ]
        strands_region = [u[3] for u in r_units_region]
        strand_inversions = sum(
            1 for i in range(1, len(strands_region)) if strands_region[i] != strands_region[i - 1]
        )

        locus_audits.append(
            {
                "chrom": r_chrom,
                "start": r_st,
                "end": r_en,
                "is_terminal": is_term,
                "unit_count": len(r_units_region),
                "strand_inversions": strand_inversions,
            }
        )

        if not is_term and (strand_inversions > 0 or circ_monomer_name is not None):
            rdna_cards.append(
                {
                    "card_id": f"CARD-RDNA-{rdna_card_idx:02d}",
                    "domain": "DOMAIN_RDNA",
                    "category": "Ribosomal RNA: Major Locus Replacement",
                    "title": f"Major rDNA Array on {r_chrom}:{r_st:,}-{r_en:,} ({r_en - r_st:,} bp)",
                    "chrom": r_chrom,
                    "start": r_st,
                    "end": r_en,
                    "action": "REPLACE_RDNA_ARRAY",
                    "classification": "COLLAPSED_INVERTED_RDNA",
                    "confidence": "HIGH",
                    "flye_monomer_id": circ_monomer_name,
                    "flye_monomer_length": circ_monomer_len,
                    "evidence": {
                        "host_chromosome": r_chrom,
                        "array_span": f"{r_en - r_st:,} bp",
                        "assembled_units": len(r_units_region),
                        "strand_inversions": strand_inversions,
                        "flye_monomer": f"{circ_monomer_name} ({circ_monomer_len:,} bp, {circ_cov}×)"
                        if circ_monomer_name
                        else "Repeat-graph pending",
                    },
                    "instructions": [
                        f"Excise baseline erroneous rDNA region on {r_chrom}:{r_st:,}-{r_en:,} bp.",
                        f"Insert Flye-resolved circular monomer unit ({circ_monomer_name}, {circ_monomer_len:,} bp) anchored to flanking chromosome core.",
                    ],
                }
            )
            rdna_card_idx += 1

    return rdna_cards, locus_audits
