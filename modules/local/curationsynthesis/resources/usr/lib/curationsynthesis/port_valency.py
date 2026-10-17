"""Global Open Port Graph, Molecule Identity Arbitration & Valency Conservation.

Implements:
  - Physical molecule deduplication (`build_molecule_ids`) so Hifiasm unitig/contig aliases
    cannot be inserted at two termini simultaneously.
  - Scarcity-aware, valency-constrained unitig-to-terminus assignment (`assign_unitigs_to_termini`).
  - Open terminus to unplaced scaffold & telomeric junction matching (`match_open_termini_to_scaffolds`).
  - Open port graph survey (`survey_open_ports`) unifying internal gaps, chromosome termini,
    and unplaced contig perimeter ports.

British spelling standard enforced across all documentation and comments.
"""

from collections import defaultdict
from .models import Port
from .policy import DEFAULT_POLICY


def candidate_preference(candidate):
    """Sort key for one terminus's eligible extension candidates.

    A `utg…` sequence is a Hifiasm unitig: a maximal non-branching stretch of the assembly graph.
    `ptg…` and `atg…` are contigs, chains of unitigs joined across junctions, so their outboard tail
    can include sequence that belongs to a different genomic context rather than continuing this
    chromosome end. A unitig that anchors a tip and continues outward is therefore preferred,
    with contigs kept as the fallback for ends no unitig reaches.
    """
    unitig_id = str(candidate[1])
    return (0 if unitig_id.startswith("utg") else 1, -candidate[0], -candidate[2])


def eligible_unitig_claims(candidates, state, germline_buf_bp=0):
    """The candidates a terminus may be capped with, in preference order.

    A genuine somatic stub requires a unitig that reaches past the break site; a germline tip lies
    beyond d_soma already and is not held to that rule.
    """
    claims = []
    for candidate in candidates:
        subgap = candidate[3] if len(candidate) > 3 else 0
        if state["is_soma_boundary"] and subgap < germline_buf_bp:
            continue
        claims.append(candidate)
    return sorted(claims, key=candidate_preference)


def terminus_state(c, e, chrom_lens, grs_segments, core_segments, delfies_points, policy=None):
    """Diminution state of one terminus: whether a break boundary exists, where, and what the tip is.

    A distal GRS segment means the assembly already extends past the breakpoint, so the tip is
    germline and d_soma is an internal calibration coordinate rather than a label for the end;
    only a core that itself reaches the tip makes the terminus a somatic stub.
    """
    if policy is None:
        policy = DEFAULT_POLICY
    t_win = policy.get("telomere", {}).get("terminal_window_bp", 10000)
    germline_buf_bp = policy.get("genomic_architecture", {}).get("germline_buffer_bp", 100)

    c_len = chrom_lens.get(c, 0)
    somatic_boundary_found = False
    has_distal_grs = False
    d_soma = 0
    if e == "5p":
        has_distal_grs = any(s == 0 for s, _end, _l, _d in grs_segments.get(c, []))
        core_starts = [s for s, _end, _l, _d in core_segments.get(c, [])]
        min_core_start = min(core_starts) if core_starts else 0
        pts = [pt for pt in delfies_points.get(c, []) if pt < min_core_start + t_win]
        if has_distal_grs and pts:
            somatic_boundary_found = True
            d_soma = min(pts)
        elif min_core_start > t_win and pts:
            somatic_boundary_found = True
            d_soma = min(pts)
    else:
        has_distal_grs = any(end >= c_len - germline_buf_bp for _s, end, _l, _d in grs_segments.get(c, []))
        core_ends = [end for _s, end, _l, _d in core_segments.get(c, [])]
        max_core_end = max(core_ends) if core_ends else c_len
        pts = [pt for pt in delfies_points.get(c, []) if pt > max_core_end - t_win]
        if has_distal_grs and pts:
            somatic_boundary_found = True
            d_soma = max(pts)
        elif max_core_end < c_len - t_win and pts:
            somatic_boundary_found = True
            d_soma = max(pts)

    tip_is_germline = has_distal_grs
    return {
        "c_len": c_len,
        "somatic_boundary_found": somatic_boundary_found,
        "has_distal_grs": has_distal_grs,
        "d_soma": d_soma,
        "tip_is_germline": tip_is_germline,
        "is_soma_boundary": somatic_boundary_found and not tip_is_germline,
    }


def assign_unitigs_to_termini(claims, molecule_of=None):
    """Assign each assembled molecule to the single terminus it extends best.

    Claims are arbitrated by:
      1. Terminus scarcity (fewer candidate molecules takes priority so constrained ends resolve first).
      2. Capping status (claims providing an authentic terminal telomeric tract take priority).
      3. Unitig preference (`utg...` preferred over `ptg...`/`atg...`).
      4. Overhang length.
    """
    if molecule_of is None:
        molecule_of = lambda chrom, candidate: candidate[1]

    distinct_molecules = defaultdict(set)
    for chrom, terminus, candidate in claims:
        distinct_molecules[(chrom, terminus)].add(molecule_of(chrom, candidate))
    scarcity = {term: len(mols) for term, mols in distinct_molecules.items()}

    ordered = sorted(
        claims,
        key=lambda x: (
            scarcity.get((x[0], x[1]), 0),
            -(1 if len(x[2]) > 2 and x[2][2] > 0 else 0),
            -x[2][2] if len(x[2]) > 2 else 0,
            0 if str(x[2][1]).startswith("utg") else 1,
            -x[2][0],
            x[0],
            x[1],
        ),
    )
    assignment = {}
    taken = set()
    for chrom, terminus, candidate in ordered:
        if (chrom, terminus) in assignment:
            continue
        molecule = molecule_of(chrom, candidate)
        if molecule in taken:
            continue
        assignment[(chrom, terminus)] = candidate
        taken.add(molecule)
    return assignment


def build_molecule_ids(lengths, placements):
    """Identity for each query's physical molecule: its length plus a shared chromosome placement."""
    groups = []
    molecule_ids = {}
    for query in sorted(lengths):
        for index, (group_len, group_placements, key) in enumerate(groups):
            if group_len == lengths[query] and group_placements & placements.get(query, set()):
                groups[index] = (group_len, group_placements | placements[query], key)
                molecule_ids[query] = key
                break
        else:
            key = f"molecule_{len(groups) + 1}"
            groups.append((lengths[query], set(placements.get(query, set())), key))
            molecule_ids[query] = key
    return molecule_ids


def build_unitig_assignment(
    nuclear_chroms, chrom_lens, grs_segments, core_segments, delfies_points, unitig_extensions, policy=None
):
    """Build arbitrated 1-to-1 molecule-to-terminus assignment across all nuclear chromosome ends."""
    if policy is None:
        policy = DEFAULT_POLICY
    germline_buf_bp = policy.get("genomic_architecture", {}).get("germline_buffer_bp", 100)

    def field(candidate, index, default=0):
        return int(candidate[index]) if len(candidate) > index else default

    lengths = {}
    placements = defaultdict(set)
    for (chrom, _terminus), candidates in unitig_extensions.items():
        for candidate in candidates:
            query = candidate[1]
            lengths[query] = field(candidate, 7)
            placements[query].add(
                (
                    chrom,
                    candidate[6] if len(candidate) > 6 else "+",
                    round(field(candidate, 8) / 50),
                    round(field(candidate, 9) / 50),
                )
            )
    molecule_ids = build_molecule_ids(lengths, placements)

    def molecule_of(chrom, candidate):
        return molecule_ids.get(candidate[1], candidate[1])

    claims = []
    for c in nuclear_chroms:
        for e in ("5p", "3p"):
            state = terminus_state(c, e, chrom_lens, grs_segments, core_segments, delfies_points, policy=policy)
            for candidate in eligible_unitig_claims(unitig_extensions.get((c, e), []), state, germline_buf_bp):
                claims.append((c, e, candidate))
    return assign_unitigs_to_termini(claims, molecule_of)


def _outward_overhang(rec, slen, term):
    """Compute the directional outward overhang (and slice coordinates) of a query on a terminus.

    At a 5' chromosome terminus:
      - '+' alignment places the outward extension in the query's 5' prefix [0, qstart].
      - '-' alignment places the outward extension in the query's 3' suffix [qend, slen].
    At a 3' chromosome terminus:
      - '+' alignment places the outward extension in the query's 3' suffix [qend, slen].
      - '-' alignment places the outward extension in the query's 5' prefix [0, qstart].
    """
    qstart = int(rec.get("qstart", 0))
    qend = int(rec.get("qend", 0))
    strand = rec.get("strand", "+")
    if (term == "5p" and strand == "+") or (term == "3p" and strand == "-"):
        return qstart, 0, qstart
    return max(0, slen - qend), qend, slen


def _extract_chrom_prefix(sc):
    """Extract the chromosome prefix from an unlocalised scaffold name (e.g. SUPER_X_unloc_1 -> SUPER_X)."""
    lower = sc.lower()
    if "_unloc" in lower:
        idx = lower.find("_unloc")
        return sc[:idx]
    return sc.split("_")[0] if "_" in sc else ""


def match_open_termini_to_scaffolds(
    nuclear_chroms,
    chrom_lens,
    unplaced_chroms,
    scaff_records,
    telo_junctions,
    open_termini,
    unitig_assignment=None,
    policy=None,
    grs_segments=None,
):
    """Perform constrained bipartite matching between open termini, unplaced scaffolds,
    and distal telomeric junction clusters.
    """
    if policy is None:
        policy = DEFAULT_POLICY

    subtelo_scaff_win = policy.get("telomere", {}).get("subtelomeric_scaffold_window_bp", 25000)
    min_ext_bp = policy.get("telomere", {}).get("min_extension_bp", 200)

    terminal_matches = {}
    scaffold_attributions = {}

    clusters = (telo_junctions or {}).get("clusters", [])
    germline_clusters = [cl for cl in clusters if cl.get("classification") == "GERMLINE_JUNCTION"]

    cl_by_target = {}
    orphan_clusters = []
    for cl in germline_clusters:
        t_chrom = cl.get("target_chrom")
        t_arm = cl.get("target_arm")
        if t_chrom in unplaced_chroms:
            cl_by_target[t_chrom] = cl
        elif t_chrom in nuclear_chroms and t_arm:
            cl_by_target[(t_chrom, t_arm)] = cl
        elif not t_chrom or cl.get("subtype") == "UNANCHORED_ORPHAN":
            orphan_clusters.append(cl)

    candidate_scaffolds_per_term = defaultdict(list)
    for sc in unplaced_chroms:
        recs = scaff_records.get(sc, [])
        slen = chrom_lens.get(sc, 0)
        for rec in recs:
            tname = rec.get("tname")
            tlen = chrom_lens.get(tname, rec.get("tlen", 0))
            if tname not in nuclear_chroms:
                continue
            term = None
            if rec.get("tstart", 999999999) <= subtelo_scaff_win:
                term = "5p"
            elif rec.get("tend", 0) >= tlen - subtelo_scaff_win:
                term = "3p"
            if term and (tname, term) in open_termini:
                oh, ext_s, ext_e = _outward_overhang(rec, slen, term)
                matches = rec.get("matches", rec.get("qend", 0) - rec.get("qstart", 0))
                mapq = rec.get("mapq", 0)
                score = matches * (1 + mapq)
                candidate_scaffolds_per_term[(tname, term)].append(
                    {
                        "scaffold": sc,
                        "score": score,
                        "strand": rec.get("strand", "+"),
                        "overhang": oh,
                        "ext_qstart": ext_s,
                        "ext_qend": ext_e,
                        "rec": rec,
                    }
                )

    assigned_scaffolds = set()
    assigned_cluster_ids = set()

    for c, e in open_termini:
        c_len = chrom_lens.get(c, 0)
        matched_scaffs = []
        orientations = {}
        matched_cluster = None

        # 1. Direct PAF candidate matches
        candidates = candidate_scaffolds_per_term.get((c, e), [])
        candidates.sort(key=lambda x: -x["score"])
        for cand in candidates:
            sc = cand["scaffold"]
            if sc not in assigned_scaffolds:
                matched_scaffs.append(sc)
                orientations[sc] = cand["strand"]
                assigned_scaffolds.add(sc)
                if sc in cl_by_target and matched_cluster is None:
                    matched_cluster = cl_by_target[sc]

        # 2. Companion scaffolds with matching chromosome prefix (e.g. SUPER_X_unloc_* or II_unloc_*)
        #    Only absorb a prefix-named scaffold into a terminal cassette if it carries a telomeric
        #    junction cluster OR has no non-terminal alignments (satellite/cassette fragment).
        for sc in unplaced_chroms:
            if sc in assigned_scaffolds:
                continue
            sc_prefix = _extract_chrom_prefix(sc)
            if sc_prefix and (sc_prefix == c or sc_prefix == f"chr{c}" or f"chr{sc_prefix}" == c):
                recs = scaff_records.get(sc, [])
                has_internal_core_hit = any(
                    r.get("tname") in nuclear_chroms
                    and r.get("tstart", 0) > subtelo_scaff_win
                    and r.get("tend", 0) < chrom_lens.get(r.get("tname"), r.get("tlen", 0)) - subtelo_scaff_win
                    and r.get("mapq", 0) >= 30
                    for r in recs
                )
                if has_internal_core_hit and sc not in cl_by_target:
                    continue
                matched_scaffs.append(sc)
                orientations[sc] = "+"
                assigned_scaffolds.add(sc)
                if sc in cl_by_target and matched_cluster is None:
                    matched_cluster = cl_by_target[sc]

        # 3. Associate distal telomeric cluster
        if matched_cluster is None:
            if (c, e) in cl_by_target:
                matched_cluster = cl_by_target[(c, e)]
            else:
                expected_motif = "GCCTAA" if e == "5p" else "TTAGGC"
                matching_orphans = [
                    cl
                    for cl in orphan_clusters
                    if cl.get("cluster_id") not in assigned_cluster_ids
                    and cl.get("dominant_motif") == expected_motif
                ]
                if matching_orphans:
                    best_orphan = max(
                        matching_orphans,
                        key=lambda x: (x.get("read_count", 0), x.get("mean_telomere_tract_bp", 0)),
                    )
                    matched_cluster = best_orphan
                    assigned_cluster_ids.add(best_orphan["cluster_id"])

        # 4. Satellite flank-homology fallback for terminal GRS arrays
        if matched_cluster is None and grs_segments is not None:
            germline_buf_bp = policy.get("genomic_architecture", {}).get("germline_buffer_bp", 100)
            if e == "5p":
                has_distal_grs = any(s == 0 for s, _end, _l, _d in grs_segments.get(c, []))
            else:
                has_distal_grs = any(
                    end >= c_len - germline_buf_bp for _s, end, _l, _d in grs_segments.get(c, [])
                )
            if has_distal_grs:
                expected_motif = "GCCTAA" if e == "5p" else "TTAGGC"
                satellite_cls = [
                    cl
                    for cl in germline_clusters
                    if cl.get("subtype") == "SATELLITE_REPEAT" and cl.get("dominant_motif") == expected_motif
                ]
                if satellite_cls:
                    matched_cluster = max(
                        satellite_cls,
                        key=lambda x: (x.get("read_count", 0), x.get("mean_telomere_tract_bp", 0)),
                    )

        if matched_scaffs or matched_cluster:
            terminal_matches[(c, e)] = {
                "cluster": matched_cluster,
                "scaffolds": matched_scaffs,
                "orientations": orientations,
            }
            if matched_cluster:
                assigned_cluster_ids.add(matched_cluster["cluster_id"])
            for sc in matched_scaffs:
                scaffold_attributions[sc] = {
                    "chromosome": c,
                    "terminus": e,
                    "action": "RETIRED_ABSORBED",
                    "note": f"absorbed into {c} {e} as terminal repeat cassette",
                }

    # 5. Generalised cluster assignment to all nuclear termini
    for c in nuclear_chroms:
        for e in ["5p", "3p"]:
            tm = terminal_matches.setdefault(
                (c, e), {"cluster": None, "scaffolds": [], "orientations": {}}
            )
            if tm.get("cluster"):
                continue

            arm_cls = [
                cl
                for cl in germline_clusters
                if cl.get("target_chrom") == c
                and cl.get("target_arm") == e
                and cl.get("cluster_id") not in assigned_cluster_ids
            ]
            if arm_cls:
                best_cl = max(
                    arm_cls, key=lambda x: (x.get("read_count", 0), x.get("mean_telomere_tract_bp", 0))
                )
                tm["cluster"] = best_cl
                for cl in arm_cls:
                    assigned_cluster_ids.add(cl["cluster_id"])
                continue

            if unitig_assignment:
                u_info = unitig_assignment.get((c, e))
                if u_info and u_info[0] >= min_ext_bp and u_info[1] != "None":
                    orphan_cls = [
                        cl for cl in orphan_clusters if cl.get("cluster_id") not in assigned_cluster_ids
                    ]
                    expected_motif = "GCCTAA" if e == "5p" else "TTAGGC"
                    matching = [cl for cl in orphan_cls if cl.get("dominant_motif") == expected_motif]
                    if matching:
                        best_orphan = max(
                            matching,
                            key=lambda x: (
                                x.get("dominant_motif") == expected_motif,
                                x.get("read_count", 0),
                                x.get("mean_telomere_tract_bp", 0),
                            ),
                        )
                        tm["cluster"] = best_orphan
                        for m in matching:
                            assigned_cluster_ids.add(m["cluster_id"])

    return terminal_matches, scaffold_attributions


def build_open_port_graph(
    nuclear_chroms,
    chrom_lens,
    open_termini,
    internal_cards,
    unplaced_chroms,
    scaff_any_records,
    grs_segments,
    core_segments,
    delfies_points,
    policy=None,
):
    """Build the unified Open Port Graph across chromosome termini, unresolved internal gaps,
    and unplaced segment perimeters.

    Returns a dictionary with:
      - 'chrom_ports': list of Port objects for open termini and unresolved internal gaps
      - 'contig_ports': dict[scaffold_id -> (Port_5p, Port_3p)]
      - 'gap_ports_by_id': dict[card_id -> (Port_left, Port_right, card_dict)]
    """
    if policy is None:
        policy = DEFAULT_POLICY

    min_pid = policy.get("gap_patch", {}).get("min_valency_anchor_identity", 0.98) * 100.0
    min_alen = policy.get("gap_patch", {}).get("min_secondary_anchor_bp", 300)

    chrom_ports = []
    for c, e in open_termini:
        state = terminus_state(c, e, chrom_lens, grs_segments, core_segments, delfies_points, policy=policy)
        coord = 0 if e == "5p" else chrom_lens.get(c, 0)
        chrom_ports.append(
            Port(
                port_id=f"{c}:{e}",
                molecule=c,
                coord=coord,
                orientation=e,
                port_type=f"TERM_{e.upper()}",
                is_capped=False,
                is_germline_tip=state["tip_is_germline"],
                is_soma_boundary=state["is_soma_boundary"],
                d_soma=state["d_soma"],
                target_chroms={c},
            )
        )

    gap_ports_by_id = {}
    for card in internal_cards or []:
        if card.get("action") == "REQUIRES_UNITIG_PATCH":
            c = card["chrom"]
            st, en = int(card["start"]), int(card["end"])
            cid = card.get("card_id", f"GAP-{c}-{st}")
            p_left = Port(
                port_id=f"{c}:{st}:GAP_LEFT",
                molecule=c,
                coord=st,
                orientation="3p",
                port_type="GAP_LEFT",
                target_chroms={c},
            )
            p_right = Port(
                port_id=f"{c}:{en}:GAP_RIGHT",
                molecule=c,
                coord=en,
                orientation="5p",
                port_type="GAP_RIGHT",
                target_chroms={c},
            )
            chrom_ports.append(p_left)
            chrom_ports.append(p_right)
            gap_ports_by_id[cid] = (p_left, p_right, card)

    contig_ports = {}
    for sc in unplaced_chroms:
        slen = chrom_lens.get(sc, 0)
        port_win = min(5000, max(1000, slen // 3)) if slen > 0 else 5000
        recs = [
            r
            for r in scaff_any_records.get(sc, [])
            if r.get("tname") in nuclear_chroms
            and r.get("pid", 0.0) >= min_pid
            and r.get("alen", 0) >= min_alen
        ]
        recs_5p = [r for r in recs if r.get("qstart", slen) <= port_win]
        recs_3p = [r for r in recs if r.get("qend", 0) >= max(0, slen - port_win)]
        chroms_5p = {r["tname"] for r in recs_5p}
        chroms_3p = {r["tname"] for r in recs_3p}

        p5 = Port(
            port_id=f"{sc}:5p",
            molecule=sc,
            coord=0,
            orientation="5p",
            port_type="CONTIG_5P",
            target_chroms=chroms_5p,
            valency_in=len(recs_5p),
            valency_out=len(chroms_5p),
        )
        p3 = Port(
            port_id=f"{sc}:3p",
            molecule=sc,
            coord=slen,
            orientation="3p",
            port_type="CONTIG_3P",
            target_chroms=chroms_3p,
            valency_in=len(recs_3p),
            valency_out=len(chroms_3p),
        )
        contig_ports[sc] = (p5, p3)

    return {
        "chrom_ports": chrom_ports,
        "contig_ports": contig_ports,
        "gap_ports_by_id": gap_ports_by_id,
    }


def match_unplaced_ports_by_valency(
    unplaced_seqs,
    nuclear_chroms,
    chrom_lens,
    scaff_any_records,
    internal_cards=None,
    open_termini=None,
    policy=None,
):
    """Match unplaced contig 5' and 3' perimeter ports against open chromosomal gap ports
    and open chromosome termini using Anchor Valency Conservation.

    Handles:
      1. Case A (Internal Gap Island Bridge / Monovalent Repeat Dock -> `INTEGRATE_UNPLACED_REPEAT`):
         - An unplaced contig's perimeter ports (5' and/or 3') align at high identity (>=98%)
           to the flanking regions of an unresolved internal gap (`REQUIRES_UNITIG_PATCH`),
           even when MAPQ = 0 due to local tandem repeats.
         - Accepted if either:
           (i) Dual-flank bridge: 5' and 3' perimeter ports bridge `GAP_LEFT` and `GAP_RIGHT`
               of a single gap in consistent orientation, OR
           (ii) Monovalent repeat locus (m=1, n=1): 100% of the contig's high-identity perimeter
                alignments across the entire nuclear genome land within that single gap's flanks.
      2. Case B (Monovalent / Chromosome-Exclusive Terminal Extension -> `INTEGRATE_INTO_CHROMOSOME`):
         - A perimeter port aligns within `terminal_window_bp` of an open chromosome terminus with
           outward overhang >= 1,000 bp at high identity (>=98%, alen >= 500 bp), even if MAPQ < 50,
           provided its perimeter alignments are exclusive to that chromosome and match a single
           open terminus (`len(candidate_termini) == 1`).

    Returns:
      - `gap_integrations`: dict[scaffold_id -> integration_dict]
      - `valency_terminal_extensions`: dict[scaffold_id -> extension_dict]
    """
    if policy is None:
        policy = DEFAULT_POLICY

    min_pid = policy.get("gap_patch", {}).get("min_valency_anchor_identity", 0.98) * 100.0
    min_alen = policy.get("gap_patch", {}).get("min_secondary_anchor_bp", 300)
    flank_win = policy.get("genomic_architecture", {}).get("gap_assembly_flank_bp", 20000)
    term_win = policy.get("telomere", {}).get("terminal_window_bp", 10000)
    default_gap_bp = policy.get("telomere", {}).get("default_scaffold_gap_bp", 100)

    nuclear_set = set(nuclear_chroms or [])
    open_term_set = set(open_termini or [])

    gap_targets = []
    for card in internal_cards or []:
        if card.get("action") == "REQUIRES_UNITIG_PATCH" and card.get("chrom") in nuclear_set:
            gap_targets.append(
                {
                    "card_id": card.get("card_id", f"GAP-{card['chrom']}-{card['start']}"),
                    "chrom": card["chrom"],
                    "start": int(card["start"]),
                    "end": int(card["end"]),
                }
            )

    gap_integrations = {}
    valency_terminal_extensions = {}
    claimed_gaps = set()

    for sc in sorted(unplaced_seqs):
        slen = unplaced_seqs[sc]
        if slen <= 0:
            continue
        port_win = min(5000, max(1000, slen // 3))

        # High-identity perimeter alignments across all nuclear chromosomes (including MAPQ = 0)
        perim_recs = [
            r
            for r in scaff_any_records.get(sc, [])
            if r.get("tname") in nuclear_set
            and r.get("pid", 0.0) >= min_pid
            and r.get("alen", 0) >= min_alen
            and (r.get("qstart", slen) <= port_win or r.get("qend", 0) >= max(0, slen - port_win))
        ]
        if not perim_recs:
            continue

        # --- 1. Evaluate Internal Gap Port Valency (Case A) ---
        if gap_targets:
            matching_gaps = {}
            off_gap_loci = 0
            for r in perim_recs:
                tname = r["tname"]
                ts, te = int(r["tstart"]), int(r["tend"])
                qs, qe = int(r["qstart"]), int(r["qend"])
                is_5p = qs <= port_win
                is_3p = qe >= max(0, slen - port_win)
                hit_any_gap = False
                for g in gap_targets:
                    if g["chrom"] != tname:
                        continue
                    g_left_s = max(0, g["start"] - flank_win)
                    g_right_e = min(chrom_lens.get(tname, te), g["end"] + flank_win)
                    hits_left = ts <= g["start"] and te >= g_left_s
                    hits_right = ts <= g_right_e and te >= g["end"]
                    if hits_left or hits_right:
                        hit_any_gap = True
                        info = matching_gaps.setdefault(
                            g["card_id"],
                            {
                                "gap": g,
                                "left_5p": [],
                                "left_3p": [],
                                "right_5p": [],
                                "right_3p": [],
                                "all_recs": [],
                            },
                        )
                        info["all_recs"].append(r)
                        if hits_left and is_5p:
                            info["left_5p"].append(r)
                        if hits_left and is_3p:
                            info["left_3p"].append(r)
                        if hits_right and is_5p:
                            info["right_5p"].append(r)
                        if hits_right and is_3p:
                            info["right_3p"].append(r)
                if not hit_any_gap:
                    off_gap_loci += 1

            # Check for consistent dual-flank bridges first, then monovalent single-locus docks
            bridged_candidates = []
            monovalent_candidates = []
            for gid, info in matching_gaps.items():
                if gid in claimed_gaps:
                    continue
                g = info["gap"]
                # '+' orientation bridge: 5p hits GAP_LEFT (+) and 3p hits GAP_RIGHT (+)
                plus_l = [r for r in info["left_5p"] if r.get("strand") == "+"]
                plus_r = [r for r in info["right_3p"] if r.get("strand") == "+"]
                # '-' orientation bridge: 3p hits GAP_LEFT (-) and 5p hits GAP_RIGHT (-)
                minus_l = [r for r in info["left_3p"] if r.get("strand") == "-"]
                minus_r = [r for r in info["right_5p"] if r.get("strand") == "-"]

                if plus_l and plus_r:
                    best_l = max(plus_l, key=lambda r: (r.get("alen", 0), r.get("pid", 0)))
                    best_r = max(plus_r, key=lambda r: (r.get("alen", 0), r.get("pid", 0)))
                    if best_l["qend"] < best_r["qstart"]:
                        bridged_candidates.append((gid, "+", best_l, best_r, info))
                elif minus_l and minus_r:
                    best_l = max(minus_l, key=lambda r: (r.get("alen", 0), r.get("pid", 0)))
                    best_r = max(minus_r, key=lambda r: (r.get("alen", 0), r.get("pid", 0)))
                    if best_r["qend"] < best_l["qstart"]:
                        bridged_candidates.append((gid, "-", best_l, best_r, info))
                elif len(matching_gaps) == 1 and off_gap_loci == 0:
                    # Strictly monovalent across the entire assembly (m=1, n=1)
                    best_r = max(info["all_recs"], key=lambda r: (r.get("alen", 0), r.get("pid", 0)))
                    monovalent_candidates.append((gid, best_r.get("strand", "+"), best_r, info))

            if len(bridged_candidates) == 1:
                gid, strand, best_l, best_r, info = bridged_candidates[0]
                g = info["gap"]
                if strand == "+":
                    ins_s, ins_e = int(best_l["qend"]), int(best_r["qstart"])
                else:
                    ins_s, ins_e = int(best_r["qend"]), int(best_l["qstart"])
                left_tight = abs(g["start"] - int(best_l["tend"])) <= 500
                right_tight = abs(int(best_r["tstart"]) - g["end"]) <= 500
                gap_integrations[sc] = {
                    "scaffold": sc,
                    "length": slen,
                    "chromosome": g["chrom"],
                    "target_gap_id": gid,
                    "gap_start": g["start"],
                    "gap_end": g["end"],
                    "strand": strand,
                    "bridge_mode": "DUAL_FLANK_ISLAND_BRIDGE",
                    "insert_qstart": ins_s if (left_tight and right_tight) else 0,
                    "insert_qend": ins_e if (left_tight and right_tight) else slen,
                    "left_gap_bp": 0 if left_tight else default_gap_bp,
                    "right_gap_bp": 0 if right_tight else default_gap_bp,
                    "mapq": max(best_l.get("mapq", 0), best_r.get("mapq", 0)),
                    "pid": round(min(best_l.get("pid", 100.0), best_r.get("pid", 100.0)), 2),
                    "action": "INTEGRATE_UNPLACED_REPEAT",
                    "note": (
                        f"bridges internal gap {gid} ({g['chrom']}:{g['start']:,}-{g['end']:,}) "
                        f"via dual perimeter repeat ports"
                    ),
                }
                claimed_gaps.add(gid)
                continue

            if len(monovalent_candidates) == 1:
                gid, strand, best_r, info = monovalent_candidates[0]
                g = info["gap"]
                gap_integrations[sc] = {
                    "scaffold": sc,
                    "length": slen,
                    "chromosome": g["chrom"],
                    "target_gap_id": gid,
                    "gap_start": g["start"],
                    "gap_end": g["end"],
                    "strand": strand,
                    "bridge_mode": "MONOVALENT_REPEAT_PORT",
                    "insert_qstart": 0,
                    "insert_qend": slen,
                    "left_gap_bp": default_gap_bp,
                    "right_gap_bp": default_gap_bp,
                    "mapq": best_r.get("mapq", 0),
                    "pid": round(best_r.get("pid", 100.0), 2),
                    "action": "INTEGRATE_UNPLACED_REPEAT",
                    "note": (
                        f"monovalent repeat port (m=1, n=1) docks exclusively into {gid} "
                        f"({g['chrom']}:{g['start']:,}-{g['end']:,})"
                    ),
                }
                claimed_gaps.add(gid)
                continue

        # --- 2. Evaluate Chromosome-Exclusive Terminal Port Valency (Case B) ---
        target_chroms = {r["tname"] for r in perim_recs}
        if len(target_chroms) == 1:
            only_chrom = next(iter(target_chroms))
            tlen = chrom_lens.get(only_chrom, 0)
            term_hits = {}
            for r in perim_recs:
                if r.get("alen", 0) < 500:
                    continue
                term = None
                tip_gap = 0
                if int(r["tstart"]) <= term_win:
                    term = "5p"
                    tip_gap = int(r["tstart"])
                elif int(r["tend"]) >= tlen - term_win:
                    term = "3p"
                    tip_gap = max(0, tlen - int(r["tend"]))
                if term and (not open_term_set or (only_chrom, term) in open_term_set):
                    oh, ext_s, ext_e = _outward_overhang(r, slen, term)
                    net_beyond = oh - tip_gap
                    if oh >= 1000 and net_beyond > 0:
                        rank_key = (net_beyond, -tip_gap, oh)
                        prev = term_hits.get(term)
                        if prev is None or rank_key > prev["rank_key"]:
                            term_hits[term] = {
                                "term": term,
                                "overhang": oh,
                                "ext_qstart": ext_s,
                                "ext_qend": ext_e,
                                "rec": r,
                                "rank_key": rank_key,
                            }

            if len(term_hits) == 1:
                best_t = next(iter(term_hits.values()))
                valency_terminal_extensions[sc] = {
                    "chromosome": only_chrom,
                    "terminus": best_t["term"],
                    "overhang": best_t["overhang"],
                    "ext_qstart": best_t["ext_qstart"],
                    "ext_qend": best_t["ext_qend"],
                    "rec": best_t["rec"],
                    "valency_mode": "CHROMOSOME_EXCLUSIVE_PORT",
                }

    return gap_integrations, valency_terminal_extensions

