"""Markdown Dossier and JSON Manifest Renderer for Module 2 Curation Synthesis.

Separates human/machine presentation formatting from biological decision logic.
British spelling standard enforced across all documentation and comments.
"""

import os
from .policy import __version__


def render_markdown_blueprint(
    sample_id,
    sorted_chroms,
    nuclear_chroms,
    unplaced_chroms,
    chrom_lens,
    core_segments,
    grs_segments,
    grs_threshold,
    grs_arm_depths,
    tip_depths_1kb,
    baseline_status,
    unitig_assignment,
    single_read_ext,
    rdna_intervals,
    rdna_dist_bp,
    terminal_matches,
    internal_cards,
    gap_provenance,
    gap_json_path,
    gap_parse_warning,
    telo_cards,
    flye_info_exists,
    flye_data,
    rdna_locus_audits,
    rdna_cards,
    triage_result,
    scaffold_cards,
    scaffold_dispositions,
    extension_by_seq,
    telomere_unitig_inventory=None,
):
    """Render the full Markdown Curation Reconnaissance Blueprint."""
    report_lines = [
        f"# Curation Reconnaissance Blueprint: `{sample_id}`\n",
        "**Synthesized from Module 1 Reconnaissance Outputs**\n",
        "## 1. Ground-Truth Baseline Karyotype & Dual Coverage (Core vs. GRS)\n",
    ]

    tot_len = 0
    tot_core_len = 0
    tot_grs_len = 0
    core_weighted_depth = 0.0
    grs_weighted_depth = 0.0

    for c in sorted_chroms:
        tlen = chrom_lens[c]
        tot_len += tlen
        c_len = sum(seg[2] for seg in core_segments[c])
        c_depth = (sum(seg[2] * seg[3] for seg in core_segments[c]) / c_len) if c_len > 0 else 0.0
        tot_core_len += c_len
        core_weighted_depth += c_len * c_depth

        g_len = sum(seg[2] for seg in grs_segments[c])
        g_depth = (sum(seg[2] * seg[3] for seg in grs_segments[c]) / g_len) if g_len > 0 else 0.0
        tot_grs_len += g_len
        grs_weighted_depth += g_len * g_depth

    cohort_core_d = (core_weighted_depth / tot_core_len) if tot_core_len > 0 else 0.0
    cohort_grs_d = (grs_weighted_depth / tot_grs_len) if tot_grs_len > 0 else 0.0
    core_pct = (tot_core_len / tot_len * 100.0) if tot_len > 0 else 0.0
    grs_pct = (tot_grs_len / tot_len * 100.0) if tot_len > 0 else 0.0

    report_lines.append(
        f"* **Classification Threshold ($T_{{\\text{{GRS}}}}$)**: `{grs_threshold:.2f}×` (depth threshold distinguishing putative GRS from somatic core)"
    )
    report_lines.append(
        f"* **Putative Somatic Core**: **{tot_core_len:,d} bp** ({tot_core_len / 1e6:.2f} Mb, {core_pct:.1f}%) | **Median Coverage**: `{cohort_core_d:.1f}×`"
    )
    report_lines.append(
        f"* **Putative GRS (Eliminated)**: **{tot_grs_len:,d} bp** ({tot_grs_len / 1e6:.2f} Mb, {grs_pct:.1f}%) | **Median Coverage**: `{cohort_grs_d:.1f}×`\n"
    )

    report_lines.append(
        "| Chromosome | Total Length (bp) | Core Length (bp) | Core Median Depth | GRS Length (bp) | GRS Median Depth | GRS Fraction (%) |"
    )
    report_lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

    for c in sorted_chroms:
        tlen = chrom_lens[c]
        c_len = sum(seg[2] for seg in core_segments[c])
        c_depth = (sum(seg[2] * seg[3] for seg in core_segments[c]) / c_len) if c_len > 0 else 0.0
        g_len = sum(seg[2] for seg in grs_segments[c])
        g_depth = (sum(seg[2] * seg[3] for seg in grs_segments[c]) / g_len) if g_len > 0 else 0.0
        chr_grs_pct = (g_len / tlen * 100.0) if tlen > 0 else 0.0
        report_lines.append(
            f"| **Chr {c}** | {tlen:,d} | {c_len:,d} | {c_depth:.1f}× | {g_len:,d} | {g_depth:.1f}× | {chr_grs_pct:.1f}% |"
        )

    report_lines.append(
        f"| **Cohort Total** | **{tot_len:,d}** | **{tot_core_len:,d}** | **{cohort_core_d:.1f}×** | **{tot_grs_len:,d}** | **{cohort_grs_d:.1f}×** | **{grs_pct:.1f}%** |\n"
    )

    unplaced_bp = sum(chrom_lens[c] for c in unplaced_chroms)
    report_lines.append(f"* **Unplaced Scaffolds Count:** {len(unplaced_chroms)}")
    report_lines.append(f"* **Unplaced Total Megabases:** {unplaced_bp / 1e6:.3f} Mb ({unplaced_bp:,d} bp)")
    report_lines.append(
        f"* **Identified Scaffolds:** {', '.join(f'`{c}`' for c in unplaced_chroms) if unplaced_chroms else 'None (100% Chromosome-Scale)'}\n"
    )

    # Step 2: Internal Cards
    report_lines.append("## 2. Curation Action Cards (Assembly Gaps & Structural Fixes)\n")
    if gap_provenance is not None and gap_json_path:
        report_lines.append(
            f"> **Evidence Source:** `{os.path.basename(gap_json_path)}` (Fingerprint: `{gap_provenance.get('script_fingerprint', 'N/A')}`, Schema: `{gap_provenance.get('schema_version', '1.0.0')}`)\n"
        )
    elif gap_parse_warning:
        report_lines.append(f"> *Warning: Could not parse gap evidence JSON ({gap_parse_warning})*\n")

    if internal_cards:
        for c in internal_cards:
            report_lines.append(f"### [{c['card_id']}] {c['category']}")
            report_lines.append(f"- **Target:** `{c['title']}`")
            report_lines.append(f"- **Action Required:** `{c['action']}`")
            report_lines.append(f"- **Confidence:** **{c['confidence']}**")
            report_lines.append("- **Evidence Breakdown:**")
            for k, v in c["evidence"].items():
                report_lines.append(f"  - {k.replace('_', ' ').title()}: `{v}`")
            report_lines.append("- **Deterministic Curation Recipe:**")
            for step in c["instructions"]:
                report_lines.append(f"  1. {step}")
            report_lines.append("")
    else:
        report_lines.append(
            "No actionable scaffolding errors, assembly gaps, novel insertions, or collapsed repeats detected. Assembly is contig-grade.\n"
        )

    # Step 3: Telomere Matrix & Cards
    tot_ends = len(nuclear_chroms) * 2
    report_lines.append(f"## 3. Dynamic {tot_ends}-End Telomere Evidence Matrix & Evidence Hierarchy\n")
    if telomere_unitig_inventory:
        tot_u = telomere_unitig_inventory.get("total_assembled", 0)
        as_u = telomere_unitig_inventory.get("assigned_count", 0)
        unas_u = telomere_unitig_inventory.get("unassigned_count", 0)
        report_lines.append(
            f"* **Independently Assembled Telomeric Unitigs**: **{tot_u} unitigs** ({as_u} assigned to chromosome ends, {unas_u} unassigned)\n"
        )
    report_lines.append(
        "| End | Baseline Tip Status | Distal GRS Depth | Tip 1kb Depth | Single Read Overhang | Best Unitig Overhang | Evidence Decision | Accidental Somatic Warning |"
    )
    report_lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

    for c in nuclear_chroms:
        for e in ("5p", "3p"):
            b_stat = baseline_status.get((c, e), "UNKNOWN")
            if (c, e) in grs_arm_depths:
                g_len, g_d = grs_arm_depths[(c, e)]
                grs_str = f"{g_d:.1f}× ({g_len:,d} bp)"
            else:
                grs_str = "None (Core to Tip)"

            tip_1k = tip_depths_1kb.get((c, e), 0.0)
            tip_1k_str = f"{tip_1k:.1f}×" if tip_1k > 0 else "N/A"

            best_u = unitig_assignment.get((c, e))
            best_u_bp, best_u_name = (best_u[0], best_u[1]) if best_u else (0, "None")

            read_ext_bp = 0
            read_ext_str = "0 bp"
            if (c, e) in single_read_ext:
                r_info = single_read_ext[(c, e)]
                read_ext_bp, read_cnt = r_info[0], r_info[1]
                read_ext_str = f"+{read_ext_bp:,d} bp ({read_cnt} reads)"

            is_rdna_term = False
            if c in rdna_intervals:
                for r_start, r_end in rdna_intervals[c]:
                    if (e == "5p" and r_start <= rdna_dist_bp) or (
                        e == "3p" and r_end >= chrom_lens.get(c, 0) - rdna_dist_bp
                    ):
                        is_rdna_term = True
                        break

            t_match = terminal_matches.get((c, e))
            junc_match = t_match.get("cluster") if t_match else None

            if is_rdna_term:
                decision = "**TERMINUS_RDNA_ARRAY (MODULE 4 CASSETTE)**"
            elif best_u_bp > 0 and best_u_bp >= read_ext_bp:
                decision = f"**EXTEND_VIA_UNITIG** (+{best_u_bp:,d} bp)"
            elif read_ext_bp > best_u_bp:
                decision = f"**EXTEND_VIA_SINGLE_READ** (+{read_ext_bp:,d} bp)"
            elif best_u_bp > 0:
                decision = f"**EXTEND_VIA_UNITIG** (+{best_u_bp:,d} bp)"
            elif b_stat == "ALREADY_CAPPED":
                decision = "**MAINTAIN_CAPPED**"
            elif junc_match:
                decision = (
                    f"**SCAFFOLD_SATELLITE_TO_TELOMERE** ({junc_match.get('read_count', 0)} reads, 100 bp N gap)"
                )
            else:
                decision = "**SEARCH_UNMAPPED_READS**"

            somatic_warning = "CLEAN (>50 kb)"
            if (c, e) in grs_arm_depths:
                g_len, _g_d = grs_arm_depths[(c, e)]
                if g_len < 20000:
                    somatic_warning = f"PROXIMAL ({g_len:,d} bp)"

            u_str = f"`{best_u_name}` (+{best_u_bp:,d} bp)" if best_u_bp > 0 else "None"
            report_lines.append(
                f"| **{c}_{e}** | `{b_stat}` | {grs_str} | {tip_1k_str} | {read_ext_str} | {u_str} | {decision} | {somatic_warning} |"
            )

    report_lines.append("\n### Domain 2 Action Cards (Chromosome Termini)\n")
    if telo_cards:
        for c in telo_cards:
            report_lines.append(f"#### [{c['card_id']}] {c['category']}")
            report_lines.append(f"- **Target:** `{c['title']}`")
            report_lines.append(
                f"- **Action Required:** `{c['action']}`"
                + (f" (`{c['action_subcode']}`)" if "action_subcode" in c else "")
            )
            report_lines.append(f"- **Classification:** `{c['classification']}`")
            report_lines.append(f"- **Confidence:** **{c['confidence']}**")
            report_lines.append("- **Evidence Breakdown:**")
            for k, v in c["evidence"].items():
                report_lines.append(f"  - {k.replace('_', ' ').title()}: `{v}`")
            report_lines.append("- **Curation Protocol:**")
            for step in c.get("instructions", []):
                report_lines.append(f"  1. {step}")
            report_lines.append("")
    else:
        report_lines.append(
            "Zero active terminus action cards required. All chromosome ends are verified complete in baseline assembly.\n"
        )

    if telomere_unitig_inventory and telomere_unitig_inventory.get("unitigs"):
        report_lines.append("### Telomeric Unitig Inventory & Hit Disposition\n")
        report_lines.append(
            "| Unitig ID | Length | Coverage / Reads | Alignment Hits (Target, Coordinates, MAPQ) | Terminal Overhang | Status | Curator Disposition |"
        )
        report_lines.append("| :--- | :---: | :---: | :--- | :---: | :---: | :--- |")
        for u in telomere_unitig_inventory["unitigs"]:
            uid = u["unitig_id"]
            ulen = f"{u['length_bp']:,d} bp"
            cov = u.get("coverage_str", "N/A")
            status = f"**{u['status']}**"
            disp = u.get("disposition", "")
            hit_strs = []
            max_ovh = 0
            for h in u.get("hits", []):
                t_arm = f" ({h['target_arm']})" if h.get("target_arm") != "internal" else ""
                hit_strs.append(
                    f"Chr {h['target_chrom']}{t_arm}: {h['target_start']:,}–{h['target_end']:,} ({h['strand']}, MAPQ {h['mapq']})"
                )
                if h.get("overhang_bp", 0) > max_ovh:
                    max_ovh = h["overhang_bp"]
            hits_cell = "<br>".join(hit_strs) if hit_strs else "Zero alignments"
            ovh_cell = f"+{max_ovh:,d} bp" if max_ovh > 0 else "0 bp"
            report_lines.append(f"| `{uid}` | {ulen} | {cov} | {hits_cell} | {ovh_cell} | {status} | {disp} |")
        report_lines.append("")
        if telomere_unitig_inventory.get("termini_lacking_unitigs"):
            lacking_str = ", ".join(telomere_unitig_inventory["termini_lacking_unitigs"])
            report_lines.append(f"* **Chromosome Termini Lacking Assembled Unitigs**: {lacking_str}\n")

    # Step 4: rDNA Audit
    report_lines.append("## 4. Major rDNA Locus Audit & Repeat-Graph Resolution\n")
    if flye_info_exists:
        report_lines.append("\n**Flye Repeat-Graph Resolution:**")
        for cname, clen, ccov, ccirc, _crep, cmult in flye_data.get("rows", []):
            circ_tag = "Circular Monomer" if ccirc == "Y" else "Linear Contig"
            report_lines.append(f"- `{cname}`: {clen:,d} bp, {ccov}× coverage ({circ_tag}, multiplicity {cmult})")
        report_lines.append(
            "- **Curation Recommendation:** Replace baseline mis-oriented array with resolved circular monomer array anchored to flanking chromosome core."
        )
    else:
        report_lines.append("\n*Flye repeat-graph output pending or not found.*")

    for audit in rdna_locus_audits:
        r_chrom, r_st, r_en = audit["chrom"], audit["start"], audit["end"]
        locus_type = "Terminal Array" if audit["is_terminal"] else "Internal Array"
        inv_cnt = audit["strand_inversions"]
        report_lines.append(
            f"* **Host Locus:** `{r_chrom}:{r_st:,d}-{r_en:,d}` ({locus_type}, Total Span: {r_en - r_st:,d} bp)"
        )
        report_lines.append(f"  - Assembled rRNA Units (Barrnap): {audit['unit_count']} units")
        report_lines.append(
            f"  - Monomer Strand Harmonization: {'**CLEAN (Unidirectional Tandem)**' if inv_cnt == 0 else f'**WARNING: {inv_cnt} STRAND INVERSIONS DETECTED**'}"
        )

    report_lines.append("\n### Domain 3 Action Cards (Ribosomal RNA)\n")
    if rdna_cards:
        for c in rdna_cards:
            report_lines.append(f"#### [{c['card_id']}] {c['category']}")
            report_lines.append(f"- **Target:** `{c['title']}`")
            report_lines.append(f"- **Action Required:** `{c['action']}`")
            report_lines.append(f"- **Confidence:** **{c['confidence']}**")
            report_lines.append("- **Evidence Breakdown:**")
            for k, v in c["evidence"].items():
                report_lines.append(f"  - {k.replace('_', ' ').title()}: `{v}`")
            report_lines.append("- **Deterministic Curation Recipe:**")
            for step in c["instructions"]:
                report_lines.append(f"  1. {step}")
            report_lines.append("")
    else:
        report_lines.append("Zero active rDNA locus replacement cards required.\n")

    # Step 5: Unplaced Scaffolds Triage
    report_lines.append("## 5. Unplaced Scaffolds Triage\n")
    report_lines.append(
        "| Scaffold | Length (bp) | Core Alignment Span | Reasoned rDNA Span | Novel Sequence Span | Repeat Density (Red %) | Triage Category | Recommended Action |"
    )
    report_lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

    table_rows = triage_result["table_rows"]
    if table_rows:
        for r in table_rows:
            report_lines.append(
                f"| `{r['scaffold']}` | {r['length']:,d} | {r['chrom_aligned_pct']:.1f}% ({r['chrom_aligned_bp']:,d} bp) | "
                f"{r['rdna_operon_pct']:.1f}% ({r['rdna_operon_bp']:,d} bp) | {r['novel_pct']:.1f}% ({r['novel_bp']:,d} bp) | "
                f"{r['rep_pct']:.1f}% | `{r['cat']}` | {r['action']} |"
            )
    else:
        report_lines.append(
            "| None | 0 | 0.0% | 0.0% | 0.0% | 0.0% | `N/A` | No unplaced scaffolds in assembly |"
        )

    if extension_by_seq:
        report_lines.append("\n**Chromosomal Port Integration Candidacy (Termini & Internal Gaps):**")
        report_lines.append("")
        report_lines.append(
            "| Scaffold | Length (bp) | Target Port | Port Mode | Inserted / Overhang (bp) | Anchored / Identity | MAPQ | Disposition |"
        )
        report_lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
        for row in scaffold_dispositions:
            if "terminus" in row:
                report_lines.append(
                    f"| `{row['scaffold']}` | {row['length']:,d} | {row['chromosome']} {row['terminus']} | "
                    f"{row['terminus_owner']} | {row['overhang_bp']:,d} | {row['anchored_span_pct']:.1f}% | "
                    f"{row['mapq']} | `{row['action']}` |"
                )
            elif "target_gap_id" in row:
                ins_bp = int(row.get("insert_qend", row["length"])) - int(row.get("insert_qstart", 0))
                report_lines.append(
                    f"| `{row['scaffold']}` | {row['length']:,d} | {row['chromosome']} `{row['target_gap_id']}` | "
                    f"{row.get('bridge_mode', 'GAP_PORT')} | {ins_bp:,d} | {row.get('pid', 100.0):.1f}% | "
                    f"{row.get('mapq', 0)} | `{row['action']}` |"
                )
        report_lines.append("")


    triage_counts = triage_result["triage_counts"]
    triage_bp = triage_result["triage_bp"]
    report_lines.append("\n**Triage Category Summary:**")
    report_lines.append(
        f"- **True Collinear Haplotigs:** {triage_counts['Haplotig']} scaffolds ({triage_bp['Haplotig'] / 1e6:.3f} Mb)"
    )
    report_lines.append(
        f"- **Orphan rDNA Fragments:** {triage_counts['rDNA_Fragment']} scaffolds ({triage_bp['rDNA_Fragment'] / 1e6:.3f} Mb)"
    )
    report_lines.append(
        f"- **Unplaced Satellite Arrays:** {triage_counts['Satellite']} scaffolds ({triage_bp['Satellite'] / 1e6:.3f} Mb)"
    )
    report_lines.append(
        f"- **Genuine Non-Redundant Contigs:** {triage_counts['Genuine_Contig']} scaffolds ({triage_bp['Genuine_Contig'] / 1e6:.3f} Mb)"
    )
    report_lines.append(
        f"- **Organellar (Mitochondrial):** {triage_counts['Organellar']} scaffolds ({triage_bp['Organellar'] / 1e6:.3f} Mb)\n"
    )

    report_lines.append("### Domain 4 Action Cards (Unplaced Scaffolds)\n")
    if scaffold_cards:
        for c in scaffold_cards:
            report_lines.append(f"#### [{c['card_id']}] {c['category']}")
            report_lines.append(f"- **Target:** `{c['title']}`")
            report_lines.append(f"- **Action Required:** `{c['action']}`")
            report_lines.append(f"- **Confidence:** **{c['confidence']}**")
            report_lines.append("- **Evidence Breakdown:**")
            for k, v in c["evidence"].items():
                report_lines.append(f"  - {k.replace('_', ' ').title()}: `{v}`")
            report_lines.append("- **Deterministic Curation Recipe:**")
            for step in c["instructions"]:
                report_lines.append(f"  1. {step}")
            report_lines.append("")
    else:
        report_lines.append(
            "Zero active unplaced scaffold action cards required. Non-redundant scaffolds retained as-is.\n"
        )

    return "\n".join(report_lines), {
        "tot_core_len": tot_core_len,
        "core_pct": core_pct,
        "cohort_core_d": cohort_core_d,
        "tot_grs_len": tot_grs_len,
        "grs_pct": grs_pct,
        "cohort_grs_d": cohort_grs_d,
    }


def build_manifest_json(
    sample_id,
    script_fingerprint,
    nuclear_chroms,
    unplaced_chroms,
    grs_threshold,
    karyotype_stats,
    policy,
    internal_cards,
    telo_cards,
    rdna_cards,
    scaffold_cards,
    scaffold_dispositions,
    telomere_unitig_inventory=None,
):
    """Assemble the machine-executable curation_manifest.json payload."""
    all_action_cards = internal_cards + telo_cards + rdna_cards + scaffold_cards
    tot_core_len = karyotype_stats["tot_core_len"]
    tot_grs_len = karyotype_stats["tot_grs_len"]
    return {
        "provenance": {
            "sample": sample_id,
            "software": "pde_curation_synthesis.py",
            "software_version": __version__,
            "schema_version": "2.0.0",
            "script_fingerprint": script_fingerprint,
        },
        "nuclear_chromosomes": nuclear_chroms,
        "unplaced_scaffolds": unplaced_chroms,
        "karyotype_compartments": {
            "classification_threshold": round(grs_threshold, 2),
            "putative_core": {
                "span_bp": tot_core_len,
                "span_mb": round(tot_core_len / 1e6, 3),
                "span_pct": round(karyotype_stats["core_pct"], 2),
                "median_coverage": round(karyotype_stats["cohort_core_d"], 2),
            },
            "putative_grs": {
                "span_bp": tot_grs_len,
                "span_mb": round(tot_grs_len / 1e6, 3),
                "span_pct": round(karyotype_stats["grs_pct"], 2),
                "median_coverage": round(karyotype_stats["cohort_grs_d"], 2),
            },
        },
        "policy": policy,
        "summary": {
            "total_action_cards": len(all_action_cards),
            "domain_internal": len(internal_cards),
            "domain_terminus": len(telo_cards),
            "domain_rdna": len(rdna_cards),
            "domain_scaffold": len(scaffold_cards),
        },
        "scaffold_dispositions": scaffold_dispositions,
        "telomere_unitig_inventory": telomere_unitig_inventory or {},
        "action_cards": all_action_cards,
    }
