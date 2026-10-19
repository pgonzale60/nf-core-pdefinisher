"""Chromosome Terminus Resolution Hierarchy Engine (Domain 2).

Implements the 7-level terminus resolution hierarchy:
  - Level 1: Terminal rDNA Array Precedence (`TERMINUS_RDNA_ARRAY`)
  - Level 2: Somatic Diminution Boundary vs. Germline Tip State (`terminus_state`)
  - Level 3: Arbitrated Unitig & Read Candidate Selection
  - Level 4: Baseline Tract vs. Overhang Extension (`DRAFT_HAS_TELOMERE` / `TELOMERE_IN_OVERHANG`)
  - Level 5a: Phase 1 Unitig Capping (`RESOLVABLE_BY_UNITIG`)
  - Level 5b: Satellite-Linked Distal Telomere Scaffolding (`SATELLITE_LINKED_TELOMERE`)
  - Level 5c: Phase 1 Single-Read Capping & Chimeric Terminus Pruning (`RESOLVABLE_BY_READ`)
  - Level 6: Somatic Diminution Boundary Walk (`SOMATIC_DIMINUTION_BOUNDARY`)
  - Level 7: Plain Uncapped Subtelomere Walk (`UNCAPPED_SUBTELOMERE`)

British spelling standard enforced across all documentation and comments.
"""

from .policy import DEFAULT_POLICY
from .port_valency import build_unitig_assignment, terminus_state


def safe_int(val, default=0):
    """Safely convert value to integer."""
    if val is None:
        return default
    if isinstance(val, (int, float)):
        return int(val)
    try:
        return int(str(val).replace(",", "").strip())
    except (ValueError, TypeError):
        return default


def detect_chimeric_terminus(chrom, terminus, chrom_len, clips, d_soma, delfies_pts, policy):
    """Detect whether an uncapped chromosome end holds an internal authentic telomeric junction
    flanked by a distal chimeric misassembly.
    """
    if not clips:
        return None

    germline_buf = policy.get("genomic_architecture", {}).get("germline_buffer_bp", 100) if policy else 100
    mapq_min = policy.get("mapping_and_anchors", {}).get("mapq_anchored", 30) if policy else 30
    mapq_floor = policy.get("mapping_and_anchors", {}).get("mapq_evidence_floor", 20) if policy else 20
    min_repeats = policy.get("telomere", {}).get("min_repeats", 15) if policy else 15
    chimeric_min_dist = policy.get("telomere", {}).get("chimeric_terminus_min_dist_bp", 1000) if policy else 1000
    chimeric_max_dist = policy.get("telomere", {}).get("chimeric_terminus_max_dist_bp", 50000) if policy else 50000
    subtelo_overlap = policy.get("telomere", {}).get("subtelomeric_min_overlap_bp", 500) if policy else 500
    anchor_len = policy.get("mapping_and_anchors", {}).get("anchor_min_length_bp", 1000) if policy else 1000

    if terminus == "3p":
        cands = []
        for r in clips:
            coord = r["align_end"]
            if r["mapq"] >= mapq_min and r["telo_cnt"] >= min_repeats:
                dist_from_tip = chrom_len - coord
                if chimeric_min_dist <= dist_from_tip <= chimeric_max_dist:
                    if d_soma and d_soma > 0 and coord <= (d_soma + subtelo_overlap):
                        continue
                    if any(abs(coord - pt) <= germline_buf for pt in (delfies_pts or [])):
                        continue
                    cands.append(r)

        if not cands:
            return None

        cands.sort(
            key=lambda x: (
                (x["align_end"] - x["align_start"]) >= (subtelo_overlap * 6),
                x["telo_cnt"],
                -x["align_end"],
            ),
            reverse=True,
        )

        for best_cand in cands:
            junc_coord = best_cand["align_end"]
            spanning_reads = [
                r
                for r in clips
                if r["align_start"] < (junc_coord - subtelo_overlap)
                and r["align_end"] > (junc_coord + anchor_len)
                and r["mapq"] >= mapq_floor
            ]
            if len(spanning_reads) == 0:
                return {
                    "junction_coord": junc_coord,
                    "chimeric_start": junc_coord,
                    "chimeric_end": chrom_len,
                    "chimeric_length": chrom_len - junc_coord,
                    "authorized_read_qname": best_cand["qname"],
                    "clip_len": best_cand["clip_len"],
                    "telo_cnt": best_cand["telo_cnt"],
                    "mapq": best_cand["mapq"],
                }

    elif terminus == "5p":
        cands = []
        for r in clips:
            coord = r["align_start"]
            if r["mapq"] >= mapq_min and r["telo_cnt"] >= min_repeats:
                dist_from_tip = coord
                if chimeric_min_dist <= dist_from_tip <= chimeric_max_dist:
                    if d_soma and d_soma > 0 and coord >= (d_soma - subtelo_overlap):
                        continue
                    if any(abs(coord - pt) <= germline_buf for pt in (delfies_pts or [])):
                        continue
                    cands.append(r)

        if not cands:
            return None

        cands.sort(
            key=lambda x: (
                (x["align_end"] - x["align_start"]) >= (subtelo_overlap * 6),
                x["telo_cnt"],
                x["align_start"],
            ),
            reverse=True,
        )

        for best_cand in cands:
            junc_coord = best_cand["align_start"]
            spanning_reads = [
                r
                for r in clips
                if r["align_start"] < (junc_coord - anchor_len)
                and r["align_end"] > (junc_coord + subtelo_overlap)
                and r["mapq"] >= mapq_floor
            ]
            if len(spanning_reads) == 0:
                return {
                    "junction_coord": junc_coord,
                    "chimeric_start": 0,
                    "chimeric_end": junc_coord,
                    "chimeric_length": junc_coord,
                    "authorized_read_qname": best_cand["qname"],
                    "clip_len": best_cand["clip_len"],
                    "telo_cnt": best_cand["telo_cnt"],
                    "mapq": best_cand["mapq"],
                }

    return None


def synthesize_telomere_cards(
    nuclear_chroms,
    chrom_lens,
    grs_segments,
    core_segments,
    delfies_points,
    rdna_intervals,
    baseline_status,
    unitig_extensions,
    single_read_ext,
    tip_tract_bp=None,
    policy=None,
    start_card_idx=1,
    telo_junctions=None,
    terminal_matches=None,
    raw_softclips=None,
    internal_cards=None,
):
    """Synthesise Domain 2 (Terminus) and Domain 3 (Terminal rDNA) Curation Action Cards."""
    if policy is None:
        policy = DEFAULT_POLICY

    min_ext_bp = policy.get("telomere", {}).get("min_extension_bp", 200)
    min_t_reps = policy.get("telomere", {}).get("min_repeats", 15)
    rdna_dist_bp = policy.get("repeat_and_rdna", {}).get("rdna_dist_from_end", 25000)
    germline_buf_bp = policy.get("genomic_architecture", {}).get("germline_buffer_bp", 100)
    motif_len = len(policy.get("telomere", {}).get("motif", "TTAGGC"))
    subtelo_overlap = policy.get("telomere", {}).get("subtelomeric_min_overlap_bp", 500)
    subtelo_scaff_win = policy.get("telomere", {}).get("subtelomeric_scaffold_window_bp", 25000)
    walker_max_iters = policy.get("telomere", {}).get("walker_max_iters", 100)

    effective_5p = {}  # chrom -> max excised end coordinate
    effective_3p = {}  # chrom -> min excised start coordinate
    if internal_cards:
        for ic in internal_cards:
            act = ic.get("action", "")
            c_target = ic.get("chrom", "")
            s = safe_int(ic.get("start", ic.get("start_pos", 0)))
            e_pos = safe_int(ic.get("end", ic.get("end_pos", 0)))
            c_l = chrom_lens.get(c_target, 0)
            if act in ("EXCISE_DUMMY_GAP", "EXCISE_GAP", "PRUNE_CHIMERIC_TERMINUS"):
                if s == 0 and e_pos > 0:
                    effective_5p[c_target] = max(effective_5p.get(c_target, 0), e_pos)
                if c_l > 0 and e_pos >= c_l and s < c_l:
                    effective_3p[c_target] = min(effective_3p.get(c_target, c_l), s)

    unitig_assignment = build_unitig_assignment(
        nuclear_chroms, chrom_lens, grs_segments, core_segments, delfies_points, unitig_extensions, policy
    )

    cards = []
    card_idx = start_card_idx

    for c in nuclear_chroms:
        orig_c_len = chrom_lens.get(c, 0)
        eff_start = effective_5p.get(c, 0)
        eff_len = effective_3p.get(c, orig_c_len)
        c_len = eff_len
        c_fully_covered_by_rdna = False
        for e in ["5p", "3p"]:
            if c_fully_covered_by_rdna:
                continue
            term_id = f"{c}_{e}"
            b_stat = baseline_status.get((c, e), "UNKNOWN")

            # LEVEL 1: Major rDNA Array Locus Precedence
            is_rdna_term = False
            matching_rdna = None
            if c in rdna_intervals:
                for r_start, r_end in rdna_intervals[c]:
                    if e == "5p" and r_start <= rdna_dist_bp:
                        is_rdna_term = True
                        clamped_end = min(r_end, c_len)
                        matching_rdna = (r_start, clamped_end)
                        if clamped_end >= c_len:
                            c_fully_covered_by_rdna = True
                        break
                    elif e == "3p" and r_end >= c_len - rdna_dist_bp:
                        is_rdna_term = True
                        matching_rdna = (max(0, r_start), c_len)
                        break

            if is_rdna_term:
                r_start, r_end = matching_rdna
                dist_to_tip = (r_start - eff_start) if e == "5p" else (eff_len - r_end)
                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_RDNA",
                        "category": "Ribosomal RNA: Terminal Array Precedence",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": eff_start if e == "5p" else r_start,
                        "end": r_end if e == "5p" else eff_len,
                        "effective_chrom_start": eff_start,
                        "effective_chrom_len": eff_len,
                        "action": "TERMINUS_RDNA_ARRAY",
                        "classification": "TERMINUS_RDNA_ARRAY",
                        "confidence": "HIGH",
                        "target_cassette": "MODULE4_RDNA_CASSETTE",
                        "evidence": {
                            "host_chromosome": c,
                            "array_coordinates": f"{c}:{r_start:,}-{r_end:,}",
                            "distance_to_tip": f"{dist_to_tip:,} bp",
                        },
                        "instructions": [
                            f"Chromosome end {term_id} abuts multi-copy rDNA array at {c}:{r_start:,}-{r_end:,} bp.",
                            "Bypass generic Telomere Walker; resolve as tripartite cassette in Module 4 (Core Anchor -> Flye rDNA Array Model -> Distal Cap Contig).",
                        ],
                    }
                )
                card_idx += 1
                continue

            # LEVEL 2: Somatic Diminution Breaksite Analysis
            state = terminus_state(c, e, chrom_lens, grs_segments, core_segments, delfies_points, policy=policy)
            somatic_boundary_found = state["somatic_boundary_found"]
            d_soma = state["d_soma"]
            tip_is_germline = state["tip_is_germline"]
            is_soma_boundary = state["is_soma_boundary"]

            # LEVEL 3: Identify Candidate Germline Unitigs & Reads
            best_u = unitig_assignment.get((c, e))
            best_u_bp = best_u[0] if best_u else 0
            best_u_name = best_u[1] if best_u else "None"
            best_u_tract = best_u[2] if best_u else 0
            best_u_subgap = best_u[3] if best_u and len(best_u) > 3 else 0

            best_u_region = None
            if best_u and len(best_u) >= 8:
                u_q_start, u_q_end = best_u[4], best_u[5]
                u_strand, u_len = best_u[6], best_u[7]
                if e == "5p":
                    u_out_start, u_out_end = (0, u_q_start) if u_strand == "+" else (u_q_end, u_len)
                else:
                    u_out_start, u_out_end = (u_q_end, u_len) if u_strand == "+" else (0, u_q_start)
                best_u_region = {
                    "q_start": u_out_start,
                    "q_end": u_out_end,
                    "strand": u_strand,
                    "unitig_len": u_len,
                }

            read_info = single_read_ext.get((c, e), (0, 0, "None", 0))
            read_ext_bp = read_info[0]
            read_cnt = read_info[1]
            read_qname = read_info[2] if len(read_info) > 2 else "None"
            read_telo_cnt = read_info[3] if len(read_info) > 3 else 0

            is_valid_germline_read = (
                read_ext_bp >= min_ext_bp and read_qname != "None" and read_telo_cnt >= min_t_reps
            )
            if is_soma_boundary and is_valid_germline_read:
                if (read_ext_bp - (read_telo_cnt * motif_len)) < germline_buf_bp and read_ext_bp < subtelo_overlap:
                    is_valid_germline_read = False

            valid_read_ext_bp = read_ext_bp if is_valid_germline_read else 0
            deeper_ext = max(best_u_bp, valid_read_ext_bp)

            # LEVEL 4: Baseline Capped vs. Overhang Extension
            if b_stat in ("ALREADY_CAPPED", "TELOMERE_IN_OVERHANG"):
                measured_tract = (tip_tract_bp or {}).get((c, e))
                if b_stat == "TELOMERE_IN_OVERHANG" or measured_tract == 0:
                    t_match = (terminal_matches or {}).get((c, e))
                    junc_match = t_match.get("cluster") if t_match else None
                    if not junc_match and telo_junctions:
                        for cl in telo_junctions.get("clusters", []):
                            if (
                                cl.get("classification") == "GERMLINE_JUNCTION"
                                and cl.get("target_chrom") == c
                                and cl.get("target_arm") == e
                            ):
                                junc_match = cl
                                break

                    # Opportunity 3: If neither unitig nor single read provides a verified extending anchor (deeper_ext < min_ext_bp),
                    # but a distal telomeric junction cluster or scaffold match exists, fall through to Level 5b
                    # (SATELLITE_LINKED_TELOMERE) rather than stranding this end in a dead-end overhang card.
                    if deeper_ext < min_ext_bp and (junc_match or (t_match and t_match.get("scaffolds"))):
                        pass
                    else:
                        evidence_dict = {
                            "baseline_status": b_stat,
                            "terminal_tract_bp": f"{measured_tract} bp",
                            "deep_spanning_overhang_bp": f"+{deeper_ext:,} bp",
                            "candidate_read": read_qname if is_valid_germline_read else "None",
                            "candidate_unitig": best_u_name,
                            **({"unitig_extension_region": best_u_region} if best_u_region else {}),
                        }
                        if junc_match:
                            r_cnt = junc_match.get("read_count", 0)
                            t_bp = int(junc_match.get("mean_telomere_tract_bp", 800))
                            motif_var = junc_match.get("dominant_motif", "GCCTAA" if e == "5p" else "TTAGGC")
                            evidence_dict["distal_cluster_id"] = junc_match.get("cluster_id")
                            evidence_dict["distal_telomere_tract_bp"] = t_bp
                            evidence_dict["distal_cluster_reads"] = r_cnt
                            evidence_dict["distal_motif"] = motif_var
                            evidence_dict["unresolved_repeat_gap"] = False

                        cards.append(
                            {
                                "card_id": f"CARD-TELO-{card_idx:02d}",
                                "domain": "DOMAIN_TERMINUS",
                                "category": "Chromosome Terminus: Extend From Read Overhang",
                                "title": f"Terminus at {term_id}",
                                "chrom": c,
                                "terminus": e,
                                "start": eff_start if e == "5p" else eff_len,
                                "end": eff_start if e == "5p" else eff_len,
                                "action": "EXTEND_FROM_OVERHANG",
                                "classification": "TELOMERE_IN_OVERHANG",
                                "confidence": "HIGH" if deeper_ext >= min_ext_bp else "MEDIUM",
                                "authorized_read_qname": read_qname
                                if (valid_read_ext_bp >= best_u_bp and is_valid_germline_read)
                                else None,
                                "authorized_unitig_id": best_u_name
                                if (best_u_bp > valid_read_ext_bp and best_u_name != "None")
                                else None,
                                "effective_chrom_start": eff_start,
                                "effective_chrom_len": eff_len,
                                "evidence": evidence_dict,
                                "instructions": [
                                    "Assembled tip ends in unique subtelomeric sequence: the draft holds no canonical tract and is not capped.",
                                    "Reads overlapping this tip carry canonical repeats in their soft-clipped overhang tails.",
                                    "Extend the tip by the overhang reaching furthest into the repeat array, then confirm the appended sequence ends in a tract.",
                                    "If no candidate reaches past the tip, hand the end to the iterative walker rather than reporting it resolved.",
                                ],
                            }
                        )
                        card_idx += 1
                        continue

                if deeper_ext >= min_ext_bp:
                    cards.append(
                        {
                            "card_id": f"CARD-TELO-{card_idx:02d}",
                            "domain": "DOMAIN_TERMINUS",
                            "category": "Chromosome Terminus: Extend Telomere Tract",
                            "title": f"Terminus at {term_id}",
                            "chrom": c,
                            "terminus": e,
                            "start": eff_start if e == "5p" else eff_len,
                            "end": eff_start if e == "5p" else eff_len,
                            "action": "EXTEND_TELOMERE_TRACT",
                            "classification": "DRAFT_HAS_TELOMERE",
                            "confidence": "HIGH",
                            "authorized_read_qname": read_qname
                            if (valid_read_ext_bp >= best_u_bp and is_valid_germline_read)
                            else None,
                            "authorized_unitig_id": best_u_name
                            if (best_u_bp > valid_read_ext_bp and best_u_name != "None")
                            else None,
                            "effective_chrom_start": eff_start,
                            "effective_chrom_len": eff_len,
                            "evidence": {
                                "baseline_status": "ALREADY_CAPPED",
                                "deep_spanning_overhang_bp": f"+{deeper_ext:,} bp",
                                "candidate_read": read_qname if is_valid_germline_read else "None",
                                "candidate_unitig": best_u_name,
                                **({"unitig_extension_region": best_u_region} if best_u_region else {}),
                            },
                            "instructions": [
                                "Draft tip contains canonical telomeric repeats in baseline assembly.",
                                "Spanning sequence extends array deeper beyond baseline.",
                            ],
                        }
                    )
                else:
                    cards.append(
                        {
                            "card_id": f"CARD-TELO-{card_idx:02d}",
                            "domain": "DOMAIN_TERMINUS",
                            "category": "Chromosome Terminus: Baseline Capped (Retain As-Is)",
                            "title": f"Terminus at {term_id}",
                            "chrom": c,
                            "terminus": e,
                            "start": eff_start if e == "5p" else eff_len,
                            "end": eff_start if e == "5p" else eff_len,
                            "action": "SKIP_WALKER",
                            "classification": "ALREADY_CAPPED",
                            "confidence": "HIGH",
                            "effective_chrom_start": eff_start,
                            "effective_chrom_len": eff_len,
                            "evidence": {
                                "baseline_status": "ALREADY_CAPPED",
                                "deep_spanning_overhang_bp": f"+{deeper_ext:,} bp",
                            },
                            "instructions": [
                                "Draft tip contains canonical telomeric repeats in baseline assembly.",
                                "No deeper spanning extension observed. Retain draft baseline intact (SKIP_WALKER).",
                            ],
                        }
                    )
                card_idx += 1
                continue

            # LEVEL 5a: Phase 1 Deterministic Unitig Capping
            if best_u_bp >= min_ext_bp and best_u_name != "None":
                t_match = (terminal_matches or {}).get((c, e))
                junc_match = t_match.get("cluster") if t_match else None

                evidence_dict = {
                    "unitig_id": best_u_name,
                    "unitig_extension_bp": f"+{best_u_bp:,} bp",
                    "subtelomeric_grs_bp": f"+{best_u_subgap:,} bp",
                    "telomeric_tract_bp": f"{best_u_tract:,} bp",
                }
                if best_u_region:
                    evidence_dict["unitig_extension_region"] = best_u_region
                if somatic_boundary_found:
                    evidence_dict["somatic_diminution_resolved"] = (
                        f"Bridges past Dsoma ({d_soma:,} bp) into germline telomere"
                    )

                instructions = [
                    f"Cap terminus directly with authorized germline telomeric unitig {best_u_name} (+{best_u_bp:,} bp)."
                ]
                if junc_match:
                    r_cnt = junc_match.get("read_count", 0)
                    t_bp = int(junc_match.get("mean_telomere_tract_bp", 800))
                    motif_var = junc_match.get("dominant_motif", "GCCTAA" if e == "5p" else "TTAGGC")
                    evidence_dict["distal_cluster_id"] = junc_match.get("cluster_id")
                    evidence_dict["distal_telomere_tract_bp"] = t_bp
                    evidence_dict["distal_cluster_reads"] = r_cnt
                    evidence_dict["distal_motif"] = motif_var
                    evidence_dict["unresolved_repeat_gap"] = False
                    instructions.append(
                        f"Extend distal tip with authentic telomeric repeat array from cluster {junc_match.get('cluster_id')} "
                        f"({t_bp:,} bp {motif_var}, supported by {r_cnt} reads; add outboard overhang only)."
                    )

                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_TERMINUS",
                        "category": "Chromosome Terminus: Capped by Unitig",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": eff_start if e == "5p" else eff_len,
                        "end": eff_start if e == "5p" else eff_len,
                        "terminus_is_germline": tip_is_germline,
                        "action": "EXTEND_TO_TELOMERE",
                        "classification": "RESOLVABLE_BY_UNITIG",
                        "mechanism": "CAP_WITH_UNITIG",
                        "authorized_unitig_id": best_u_name,
                        "confidence": "HIGH",
                        "effective_chrom_start": eff_start,
                        "effective_chrom_len": eff_len,
                        "evidence": evidence_dict,
                        "instructions": instructions,
                    }
                )
                card_idx += 1
                continue

            # LEVEL 5b: Satellite Repeat Linked to Distal Telomere Cap
            t_match = (terminal_matches or {}).get((c, e))
            junc_match = t_match.get("cluster") if t_match else None
            matched_scaffs = t_match.get("scaffolds", []) if t_match else []
            scaff_orients = t_match.get("orientations", {}) if t_match else {}

            if not junc_match and telo_junctions:
                for cl in telo_junctions.get("clusters", []):
                    if (
                        cl.get("classification") == "GERMLINE_JUNCTION"
                        and cl.get("target_chrom") == c
                        and cl.get("target_arm") == e
                    ):
                        junc_match = cl
                        break

            alns = junc_match.get("alignments", []) if junc_match else []
            chr_alns = [a for a in alns if a.get("target") == c]
            dist_to_tip = None
            if chr_alns:
                if e == "5p":
                    dist_to_tip = max(0, min(a.get("target_start", 0) for a in chr_alns) - eff_start)
                else:
                    dist_to_tip = max(0, eff_len - max(a.get("target_end", 0) for a in chr_alns))

            has_anchored_cluster = (
                junc_match is not None
                and junc_match.get("read_count", 0) >= 3
                and (dist_to_tip is None or dist_to_tip <= subtelo_scaff_win)
            )

            if has_anchored_cluster:
                r_cnt = junc_match.get("read_count", 0)
                default_tract = policy.get("telomere", {}).get("tract_anchor_bp", 500)
                t_bp = int(junc_match.get("mean_telomere_tract_bp", default_tract))
                motif_var = junc_match.get("dominant_motif", "GCCTAA" if e == "5p" else "TTAGGC")
                allow_gapless = policy.get("telomere", {}).get("allow_subtelomeric_gapless_splicing", True)
                subtelo_gap = (
                    0
                    if (allow_gapless and matched_scaffs)
                    else policy.get("telomere", {}).get("default_scaffold_gap_bp", 100)
                )
                evidence_dict = {
                    "junction_cluster_id": junc_match.get("cluster_id"),
                    "subtype": junc_match.get("subtype"),
                    "read_count": r_cnt,
                    "supporting_reads": junc_match.get("supporting_read_ids", []),
                    "telomeric_tract_bp": t_bp,
                    "predominant_motif": motif_var,
                    "consensus_flank_60bp": junc_match.get("consensus_flank_60bp"),
                    "target_chrom": c,
                    "target_arm": e,
                    "scaffold_gap_bp": subtelo_gap,
                }
                if dist_to_tip is not None:
                    evidence_dict["distal_anchor_dist_to_tip_bp"] = dist_to_tip
                if matched_scaffs:
                    evidence_dict["associated_scaffolds"] = matched_scaffs
                    evidence_dict["scaffold_orientations"] = scaff_orients

                instructions = [
                    "Terminus ends in extensive unbridgeable satellite array linking directly to distal telomere cap.",
                    f"Identified {r_cnt} reads spanning the satellite-to-telomere transition.",
                ]
                if dist_to_tip is not None:
                    instructions.append(
                        f"Cluster {junc_match.get('cluster_id')} distal alignment anchor is just {dist_to_tip:,} bp from the physical assembly tip."
                    )
                if matched_scaffs:
                    instructions.append(
                        f"Incorporate unplaced subtelomeric scaffold cassette ({', '.join(matched_scaffs)}) into {c} {e}."
                    )
                if subtelo_gap == 0:
                    instructions.append(
                        f"Scaffold assembly tip to distal telomeric cap (~{t_bp} bp {motif_var}) gaplessly (0 Ns) using spanning read evidence."
                    )
                else:
                    instructions.append(
                        f"Scaffold assembly tip to distal telomeric cap (~{t_bp} bp {motif_var}) across a nominal {subtelo_gap} bp N gap."
                    )

                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_TERMINUS",
                        "category": "Chromosome Terminus: Satellite Array Linked to Telomere",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": eff_start if e == "5p" else eff_len,
                        "end": eff_start if e == "5p" else eff_len,
                        "terminus_is_germline": tip_is_germline,
                        "action": "SCAFFOLD_REPEAT_TO_TELOMERE",
                        "classification": "SATELLITE_LINKED_TELOMERE",
                        "mechanism": "SCAFFOLD_REPEAT_TO_TELOMERE",
                        "gap_size": subtelo_gap,
                        "confidence": "HIGH",
                        "effective_chrom_start": eff_start,
                        "effective_chrom_len": eff_len,
                        "evidence": evidence_dict,
                        "instructions": instructions,
                    }
                )
                card_idx += 1
                continue

            # LEVEL 5c: Phase 1 Deterministic Read Capping
            if is_valid_germline_read:
                evidence_dict = {
                    "authorized_read_qname": read_qname,
                    "read_extension_bp": f"+{read_ext_bp:,} bp",
                    "telo_repeats": read_telo_cnt,
                    "read_count": read_cnt,
                    "supporting_reads": [read_qname] if read_qname and read_qname != "None" else [],
                }
                if somatic_boundary_found:
                    evidence_dict["somatic_diminution_resolved"] = (
                        f"Bridges past Dsoma ({d_soma:,} bp) into germline telomere"
                    )

                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_TERMINUS",
                        "category": "Chromosome Terminus: Capped by Read",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": eff_start if e == "5p" else eff_len,
                        "end": eff_start if e == "5p" else eff_len,
                        "terminus_is_germline": tip_is_germline,
                        "action": "EXTEND_TO_TELOMERE",
                        "classification": "RESOLVABLE_BY_READ",
                        "mechanism": "CAP_WITH_READ",
                        "authorized_read_qname": read_qname,
                        "confidence": "HIGH",
                        "effective_chrom_start": eff_start,
                        "effective_chrom_len": eff_len,
                        "evidence": evidence_dict,
                        "instructions": [
                            f"Cap terminus directly with terminal soft-clipped read {read_qname} extending +{read_ext_bp:,} bp ({read_telo_cnt} repeats)."
                        ],
                    }
                )
                card_idx += 1
                continue

            # LEVEL 5d: Internal Authentic Telomere with Chimeric Terminal Misassembly
            chimera = detect_chimeric_terminus(
                c, e, c_len, (raw_softclips or {}).get((c, e), []), d_soma, delfies_points.get(c, []), policy
            )
            if chimera:
                junc = chimera["junction_coord"]
                qname = chimera["authorized_read_qname"]
                telo_reps = chimera["telo_cnt"]
                clip_bp = chimera["clip_len"]
                c_start = chimera["chimeric_start"]
                c_end = chimera["chimeric_end"]
                c_len_pruned = chimera["chimeric_length"]
                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_TERMINUS",
                        "category": "Chromosome Terminus: Prune Terminal Chimera and Cap",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": c_start,
                        "end": c_end,
                        "terminus_is_germline": tip_is_germline,
                        "action": "PRUNE_CHIMERIC_TERMINUS",
                        "classification": "RESOLVABLE_BY_READ",
                        "mechanism": "CAP_WITH_READ",
                        "authorized_read_qname": qname,
                        "confidence": "HIGH",
                        "effective_chrom_len": junc if e == "3p" else c_len,
                        "effective_chrom_start": 0 if e == "3p" else junc,
                        "pruned_terminal_bp": c_len_pruned,
                        "evidence": {
                            "authentic_telomere_coord": f"{c}:{junc:,} bp",
                            "chimeric_segment_pruned": f"{c}:{c_start:,}-{c_end:,} bp ({c_len_pruned:,} bp)",
                            "authorized_read": qname,
                            "telomeric_repeats": telo_reps,
                            "softclip_bp": clip_bp,
                            "mapq": chimera["mapq"],
                        },
                        "instructions": [
                            f"Terminal sequence {c}:{c_start:,}-{c_end:,} bp ({c_len_pruned:,} bp) represents a chimeric misassembly scaffolded beyond the authentic telomere.",
                            f"Prune chimeric segment {c}:{c_start:,}-{c_end:,} bp.",
                            f"Cap authentic terminus at {junc:,} bp with authorized telomeric read {qname} (+{clip_bp} bp, {telo_reps} repeats).",
                        ],
                    }
                )
                card_idx += 1
                continue

            # LEVEL 6: Unresolved Somatic Diminution Boundary -> Iterative Walker
            if is_soma_boundary and d_soma > 0:
                cards.append(
                    {
                        "card_id": f"CARD-TELO-{card_idx:02d}",
                        "domain": "DOMAIN_TERMINUS",
                        "category": "Chromosome Terminus: Somatic Diminution Boundary",
                        "title": f"Terminus at {term_id}",
                        "chrom": c,
                        "terminus": e,
                        "start": eff_start if e == "5p" else eff_len,
                        "end": eff_start if e == "5p" else eff_len,
                        "action": "EXTEND_TO_TELOMERE",
                        "classification": "SOMATIC_DIMINUTION_BOUNDARY",
                        "action_subcode": "WALK_PAST_BREAKSITE",
                        "confidence": "HIGH",
                        "effective_chrom_start": eff_start,
                        "effective_chrom_len": eff_len,
                        "is_somatic_boundary": True,
                        "somatic_boundary_coord": d_soma,
                        "forbidden_zone_bp": d_soma,
                        "min_germline_extension_bp": 100,
                        "walker_max_iters": walker_max_iters,
                        "evidence": {
                            "diminution_breakpoint": f"{d_soma:,} bp",
                            "somatic_exclusion_zone": f"<= {d_soma:,} bp",
                            "enforced_min_germline_distance": f"> {d_soma:,} + 100 bp",
                            "core_segment_status": "Terminates at core boundary (zero germline unitigs/reads bridge into telomere)",
                        },
                        "instructions": [
                            f"Tip terminates at putative_core boundary (somatic diminution breaksite at {d_soma:,} bp).",
                            "Somatic telomere unitigs/reads disqualified from Phase 1 direct capping.",
                            "Iterative walker must extend >= Dsoma + 100 bp into GRS before accepting canonical telomere cap.",
                            "If zero reads extend past Dsoma, preserve as unclosed stub (SOMA_TRUNCATED).",
                        ],
                    }
                )
                card_idx += 1
                continue

            # LEVEL 7: Plain Uncapped Subtelomere -> Phase 2 Iterative Consensus Walker
            evidence_dict = {"baseline_status": b_stat}
            if somatic_boundary_found:
                evidence_dict["somatic_boundary_coord"] = f"{d_soma:,} bp (internal; terminus lies distal to it)"
                evidence_dict["terminus_is_germline"] = str(tip_is_germline)
            cards.append(
                {
                    "card_id": f"CARD-TELO-{card_idx:02d}",
                    "domain": "DOMAIN_TERMINUS",
                    "category": "Chromosome Terminus: Uncapped Subtelomere",
                    "title": f"Terminus at {term_id}",
                    "chrom": c,
                    "terminus": e,
                    "start": eff_start if e == "5p" else eff_len,
                    "end": eff_start if e == "5p" else eff_len,
                    "terminus_is_germline": tip_is_germline,
                    "action": "EXTEND_TO_TELOMERE",
                    "classification": "UNCAPPED_SUBTELOMERE",
                    "mechanism": "WALK_TO_TELOMERE",
                    "confidence": "MEDIUM",
                    "effective_chrom_start": eff_start,
                    "effective_chrom_len": eff_len,
                    "walker_max_iters": walker_max_iters,
                    "evidence": evidence_dict,
                    "instructions": [
                        "Terminus ends in unique subtelomeric sequence with zero pre-assembled caps.",
                        "Perform iterative consensus walking outward using candidate telomeric and unmapped reads.",
                        "If canonical array reached -> EXTEND_TO_TELOMERE. If iteration limit reached -> preserve partial sequence as EXTEND_SUBTELOMERE / INCOMPLETE.",
                        "If nothing extends the tip at all -> GERMLINE_TIP_UNCAPPED (no sequence lies beyond the terminus), never a somatic stub.",
                    ],
                }
            )
            card_idx += 1

    return cards
