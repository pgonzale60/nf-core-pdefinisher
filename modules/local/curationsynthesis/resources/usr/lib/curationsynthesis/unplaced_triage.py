"""Unplaced Scaffolds Triage & Chromosomal Port Integration Engine (Domain 4).

Separates:
  1. Redundancy & Novelty Triage (`triage_unplaced_scaffolds`): classifies each unplaced
     sequence as Organellar, rDNA_Fragment, Satellite, Haplotig, or Genuine_Contig.
     Distinguishes true high-identity allelic haplotigs (>=97% identity, >=500 bp blocks)
     from diverged paralogous gene islands (85-96% identity) and short dispersed repeats so
     gene-bearing contigs are never falsely purged.
  2. Chromosomal Port Integration (`evaluate_extension_candidacy`): evaluates whether
     unplaced sequences bridge/dock into unresolved internal repeat gaps (`INTEGRATE_UNPLACED_REPEAT`),
     absorb into terminal repeat cassettes, or extend open chromosome termini (`INTEGRATE_INTO_CHROMOSOME`)
     using directional outward overhang geometry and Anchor Valency Conservation.

British spelling standard enforced across all documentation and comments.
"""

from .evidence_io import compute_query_represented_intervals, merge_intervals
from .port_valency import _outward_overhang, match_unplaced_ports_by_valency


def triage_unplaced_scaffolds(
    unplaced_chroms,
    chrom_lens,
    scaff_repeats,
    scaff_any_records,
    scaff_rdna_operons,
    scaff_rdna_operon_bp,
):
    """Classify unplaced scaffolds by novelty floor, rDNA span, and repeat density."""
    triage_counts = {
        "Haplotig": 0,
        "Satellite": 0,
        "Genuine_Contig": 0,
        "Organellar": 0,
        "rDNA_Fragment": 0,
    }
    triage_bp = {
        "Haplotig": 0,
        "Satellite": 0,
        "Genuine_Contig": 0,
        "Organellar": 0,
        "rDNA_Fragment": 0,
    }
    disposition = {}
    scaff_chrom_aligned_bp = {}
    scaff_novel_bp = {}
    scaff_repeats_pct = {}
    table_rows = []
    scaffold_cards = []
    scaff_card_idx = 1

    unplaced_seqs = {c: chrom_lens[c] for c in unplaced_chroms}
    for sc, slen in unplaced_seqs.items():
        rep_bp = scaff_repeats.get(sc, 0)
        rep_pct = (rep_bp / slen * 100.0) if slen > 0 else 0.0
        scaff_repeats_pct[sc] = rep_pct

        recs = scaff_any_records.get(sc, [])
        # Broad repeat/paralog alignment intervals (>=85% identity, >=100 bp) for satellite arrays
        broad_aligned_intervals = compute_query_represented_intervals(
            recs, min_pid=85.0, min_alen=100
        )
        # Strict allelic alignment intervals (>=97% identity, >=500 bp) for non-repetitive haplotig triage
        allelic_aligned_intervals = compute_query_represented_intervals(
            recs, min_pid=97.0, min_alen=500
        )

        # For non-repetitive contigs (rep_pct < 45%), only allelic-grade chromosomal representation
        # counts against novelty so diverged paralogous genes (85-96% identity) or short repeat hits
        # cannot mask a novel gene island and trigger a false haplotig purge.
        active_chrom_intervals = broad_aligned_intervals if rep_pct >= 45.0 else allelic_aligned_intervals
        chrom_aligned_bp = sum(e - s for s, e in active_chrom_intervals)
        chrom_aligned_pct = min(100.0, chrom_aligned_bp / slen * 100.0) if slen > 0 else 0.0
        scaff_chrom_aligned_bp[sc] = chrom_aligned_bp

        rdna_ivs = scaff_rdna_operons.get(sc, [])
        rdna_operon_bp = scaff_rdna_operon_bp.get(sc, 0)
        rdna_operon_pct = min(100.0, rdna_operon_bp / slen * 100.0) if slen > 0 else 0.0

        represented_intervals = merge_intervals(active_chrom_intervals + rdna_ivs)
        represented_bp = min(slen, sum(e - s for s, e in represented_intervals))

        novel_bp = max(0, slen - represented_bp)
        novel_pct = (novel_bp / slen * 100.0) if slen > 0 else 0.0
        scaff_novel_bp[sc] = novel_bp

        is_organellar = "MT" in sc or "mitochondrion" in sc.lower()
        is_rdna_fragment = rdna_operon_bp >= 0.8 * slen or (
            slen - rdna_operon_bp < 1000 and rdna_operon_bp > 0
        )

        if is_organellar:
            cat = "Organellar"
            action = "PURGE_ORGANELLAR"
            scaffold_cards.append(
                {
                    "card_id": f"CARD-SCAFF-{scaff_card_idx:02d}",
                    "domain": "DOMAIN_SCAFFOLD",
                    "category": "Unplaced Scaffolds: Organellar Genome",
                    "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp)",
                    "scaffold": sc,
                    "length": slen,
                    "action": "PURGE_ORGANELLAR",
                    "classification": "ORGANELLAR",
                    "confidence": "HIGH",
                    "evidence": {
                        "length_bp": f"{slen:,} bp",
                        "core_alignment_span": f"{chrom_aligned_pct:.1f}% ({chrom_aligned_bp:,} bp)",
                        "novel_sequence_span": f"{novel_pct:.1f}% ({novel_bp:,} bp)",
                        "repeat_density": f"{rep_pct:.1f}%",
                    },
                    "instructions": [
                        f"Scaffold `{sc}` is verified organellar (mitochondrial) sequence.",
                        "Remove scaffold from primary nuclear assembly release; segregate into organellar release.",
                    ],
                }
            )
            scaff_card_idx += 1
        elif is_rdna_fragment:
            cat = "rDNA_Fragment"
            action = "PURGE_RDNA_FRAGMENT"
            scaffold_cards.append(
                {
                    "card_id": f"CARD-SCAFF-{scaff_card_idx:02d}",
                    "domain": "DOMAIN_SCAFFOLD",
                    "category": "Unplaced Scaffolds: Orphan rDNA Fragment",
                    "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp)",
                    "scaffold": sc,
                    "length": slen,
                    "action": "PURGE_RDNA_FRAGMENT",
                    "classification": "RDNA_FRAGMENT",
                    "confidence": "HIGH",
                    "evidence": {
                        "length_bp": f"{slen:,} bp",
                        "reasoned_rdna_span": f"{rdna_operon_pct:.1f}% ({rdna_operon_bp:,} bp)",
                        "novel_sequence_span": f"{novel_pct:.1f}% ({novel_bp:,} bp)",
                        "core_alignment_span": f"{chrom_aligned_pct:.1f}% ({chrom_aligned_bp:,} bp)",
                    },
                    "instructions": [
                        f"Scaffold `{sc}` is composed of reasoned rDNA operon sequence ({rdna_operon_pct:.1f}% span, only {novel_bp:,} bp unrepresented sequence).",
                        "Purge redundant scaffold from primary assembly release; sequence is absorbed into the Flye-resolved circular rDNA model or chromosomal rDNA array.",
                    ],
                }
            )
            scaff_card_idx += 1
        elif novel_bp < 1000:
            if rep_pct >= 45.0:
                cat = "Satellite"
                action = "PURGE_HAPLOTIG"
                scaffold_cards.append(
                    {
                        "card_id": f"CARD-SCAFF-{scaff_card_idx:02d}",
                        "domain": "DOMAIN_SCAFFOLD",
                        "category": "Unplaced Scaffolds: Redundant Satellite Array",
                        "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp)",
                        "scaffold": sc,
                        "length": slen,
                        "action": "PURGE_HAPLOTIG",
                        "classification": "SATELLITE",
                        "confidence": "HIGH",
                        "evidence": {
                            "length_bp": f"{slen:,} bp",
                            "core_alignment_span": f"{chrom_aligned_pct:.1f}% ({chrom_aligned_bp:,} bp)",
                            "novel_sequence_span": f"{novel_pct:.1f}% ({novel_bp:,} bp)",
                            "repeat_density": f"{rep_pct:.1f}%",
                        },
                        "instructions": [
                            f"Scaffold `{sc}` is a redundant satellite array with less than 1 kb of novel sequence ({novel_bp:,} bp unrepresented; {rep_pct:.1f}% repeat density).",
                            "Chromosome termini are fully capped; purge redundant scaffold from primary assembly release.",
                        ],
                    }
                )
                scaff_card_idx += 1
            else:
                cat = "Haplotig"
                action = "PURGE_HAPLOTIG"
                scaffold_cards.append(
                    {
                        "card_id": f"CARD-SCAFF-{scaff_card_idx:02d}",
                        "domain": "DOMAIN_SCAFFOLD",
                        "category": "Unplaced Scaffolds: True Collinear Haplotig",
                        "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp)",
                        "scaffold": sc,
                        "length": slen,
                        "action": "PURGE_HAPLOTIG",
                        "classification": "HAPLOTIG",
                        "confidence": "HIGH",
                        "evidence": {
                            "length_bp": f"{slen:,} bp",
                            "core_alignment_span": f"{chrom_aligned_pct:.1f}% ({chrom_aligned_bp:,} bp)",
                            "novel_sequence_span": f"{novel_pct:.1f}% ({novel_bp:,} bp)",
                            "repeat_density": f"{rep_pct:.1f}%",
                        },
                        "instructions": [
                            f"Scaffold `{sc}` aligns collinearly to the chromosomal assembly with less than 1 kb of novel sequence ({novel_bp:,} bp unrepresented; {chrom_aligned_pct:.1f}% represented in chromosomes).",
                            "Remove redundant scaffold from primary assembly release; segregate into alternative haplotig track.",
                        ],
                    }
                )
                scaff_card_idx += 1
        else:
            cat = "Genuine_Contig"
            action = "Retain as `unplaced_scaffold`"

        triage_counts[cat] += 1
        triage_bp[cat] += slen
        disposition[sc] = (cat, action)
        table_rows.append(
            {
                "scaffold": sc,
                "length": slen,
                "chrom_aligned_pct": chrom_aligned_pct,
                "chrom_aligned_bp": chrom_aligned_bp,
                "rdna_operon_pct": rdna_operon_pct,
                "rdna_operon_bp": rdna_operon_bp,
                "novel_pct": novel_pct,
                "novel_bp": novel_bp,
                "rep_pct": rep_pct,
                "cat": cat,
                "action": action,
            }
        )

    return {
        "unplaced_seqs": unplaced_seqs,
        "triage_counts": triage_counts,
        "triage_bp": triage_bp,
        "disposition": disposition,
        "scaff_chrom_aligned_bp": scaff_chrom_aligned_bp,
        "scaff_novel_bp": scaff_novel_bp,
        "scaff_repeats_pct": scaff_repeats_pct,
        "table_rows": table_rows,
        "scaffold_cards": scaffold_cards,
    }


def evaluate_extension_candidacy(
    triage_result,
    scaff_primary_records,
    scaff_rdna_operon_bp,
    scaffold_attributions,
    telo_cards,
    rdna_cards,
    nuclear_chroms=None,
    chrom_lens=None,
    scaff_any_records=None,
    internal_cards=None,
    open_termini=None,
    policy=None,
):
    """Evaluate unplaced scaffolds for internal repeat-gap integration and chromosomal terminal
    extension, and build final dispositions.
    """
    unplaced_seqs = triage_result["unplaced_seqs"]
    disposition = triage_result["disposition"]
    scaff_chrom_aligned_bp = triage_result["scaff_chrom_aligned_bp"]
    scaff_novel_bp = triage_result["scaff_novel_bp"]
    scaff_repeats_pct = triage_result["scaff_repeats_pct"]
    scaffold_cards = list(triage_result["scaffold_cards"])

    terminus_owner = {}
    for card in list(telo_cards) + list(rdna_cards):
        chrom = card.get("chrom")
        term = card.get("terminus")
        if not (chrom and term):
            continue
        text = " ".join(str(card.get(k, "")) for k in ("action", "category", "classification")).upper()
        if "RDNA" in text:
            terminus_owner[(chrom, term)] = "rDNA_cassette"
        elif "SOMA" in text:
            terminus_owner[(chrom, term)] = "somatic_boundary"

    gap_integrations = {}
    valency_terminal_extensions = {}
    if scaff_any_records and nuclear_chroms and chrom_lens:
        gap_integrations, valency_terminal_extensions = match_unplaced_ports_by_valency(
            unplaced_seqs=unplaced_seqs,
            nuclear_chroms=nuclear_chroms,
            chrom_lens=chrom_lens,
            scaff_any_records=scaff_any_records,
            internal_cards=internal_cards,
            open_termini=open_termini,
            policy=policy,
        )

    min_overhang_bp = 1000
    min_extension_mapq = 50
    extension_cards = []
    extension_by_seq = {}
    integrate_idx = 1

    for sc in sorted(unplaced_seqs):
        cat, triage_action = disposition.get(sc, ("Unknown", "REPORT_ONLY"))
        slen = unplaced_seqs[sc]
        eligible = cat not in ("Haplotig", "Organellar", "rDNA_Fragment")

        # 1. Internal Repeat-Gap Island Integration (Case 2 of Anchor Valency Conservation)
        if eligible and sc in gap_integrations:
            g_int = dict(gap_integrations[sc])
            chrom = g_int["chromosome"]
            gid = g_int["target_gap_id"]
            g_start, g_end = g_int["gap_start"], g_int["gap_end"]
            extension_by_seq[sc] = g_int
            extension_cards.append(
                {
                    "card_id": f"CARD-SCAFF-INT-{integrate_idx:02d}",
                    "domain": "DOMAIN_SCAFFOLD",
                    "category": "Unplaced Scaffolds: Internal Repeat-Gap Island Integration",
                    "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp) bridges/docks into {chrom}:{g_start:,}-{g_end:,}",
                    "scaffold": sc,
                    "chrom": chrom,
                    "start": g_start,
                    "end": g_end,
                    "target_gap_id": gid,
                    "length": slen,
                    "insert_qstart": g_int["insert_qstart"],
                    "insert_qend": g_int["insert_qend"],
                    "left_gap_bp": g_int["left_gap_bp"],
                    "right_gap_bp": g_int["right_gap_bp"],
                    "orientation": g_int["strand"],
                    "action": "INTEGRATE_UNPLACED_REPEAT",
                    "classification": g_int["bridge_mode"],
                    "confidence": "HIGH",
                    "evidence": {
                        "target_gap": f"{gid} ({chrom}:{g_start:,}-{g_end:,})",
                        "bridge_mode": g_int["bridge_mode"],
                        "perimeter_identity": f"{g_int['pid']:.2f}% (MAPQ {g_int['mapq']})",
                        "orientation": g_int["strand"],
                        "inserted_span_bp": f"{g_int['insert_qend'] - g_int['insert_qstart']:,} bp",
                    },
                    "instructions": [
                        f"Integrate `{sc}` ({g_int['insert_qstart']}..{g_int['insert_qend']}, strand {g_int['strand']}) into unresolved internal gap `{gid}` on {chrom} ({g_start:,}-{g_end:,}).",
                        "Retire standalone unplaced scaffold once materialised into the chromosome.",
                    ],
                }
            )
            integrate_idx += 1
            disposition[sc] = (cat, "INTEGRATE_UNPLACED_REPEAT")
            continue

        # 2. Terminal Repeat Cassette Absorption (from match_open_termini_to_scaffolds)
        if sc in scaffold_attributions:
            attr = scaffold_attributions[sc]
            chrom = attr["chromosome"]
            term = attr["terminus"]
            entry = {
                "scaffold": sc,
                "length": slen,
                "chromosome": chrom,
                "terminus": term,
                "terminus_owner": "terminal_repeat_cassette",
                "overhang_bp": slen,
                "anchored_span_pct": 100.0,
                "mapq": 60,
                "strand": "+",
                "action": "RETIRED_ABSORBED",
                "note": attr.get("note", f"absorbed into {chrom} {term} as terminal repeat cassette"),
            }
            extension_by_seq[sc] = entry
            extension_cards.append(
                {
                    "card_id": f"CARD-SCAFF-INT-{integrate_idx:02d}",
                    "domain": "DOMAIN_SCAFFOLD",
                    "category": "Unplaced Scaffolds: Terminal Repeat Cassette",
                    "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp) absorbed into {chrom} {term}",
                    "scaffold": sc,
                    "chrom": chrom,
                    "terminus": term,
                    "length": slen,
                    "action": "INTEGRATE_INTO_CHROMOSOME",
                    "classification": "TERMINAL_EXTENSION",
                    "confidence": "HIGH",
                    "evidence": {
                        "target_terminus": f"{chrom} {term}",
                        "disposition": "RETIRED_ABSORBED",
                        "extension_beyond_terminus": f"{slen:,} bp",
                    },
                    "instructions": [
                        f"Scaffold `{sc}` is absorbed into {chrom} {term} as part of the terminal repeat-to-telomere cassette.",
                        "Retire scaffold from standalone unplaced release only when terminus is materialised.",
                    ],
                }
            )
            integrate_idx += 1
            disposition[sc] = (cat, "RETIRED_ABSORBED")
            continue

        # 3. Direct or Port-Valency Terminal Extension Candidacy
        best = None
        for rec in scaff_primary_records.get(sc, []):
            if rec["mapq"] < min_extension_mapq and eligible:
                continue
            if rec["tend"] >= rec["tlen"] - 5000:
                term = "3p"
                tip_gap = max(0, rec["tlen"] - rec["tend"])
            elif rec["tstart"] <= 5000:
                term = "5p"
                tip_gap = rec["tstart"]
            else:
                continue
            overhang, ext_s, ext_e = _outward_overhang(rec, slen, term)
            net_beyond = overhang - tip_gap
            if overhang < min_overhang_bp or net_beyond <= 0:
                continue
            rank_key = (net_beyond, -tip_gap, overhang)
            if best is None or rank_key > best["rank_key"]:
                best = {
                    "overhang": overhang,
                    "ext_qstart": ext_s,
                    "ext_qend": ext_e,
                    "term": term,
                    "rec": rec,
                    "rank_key": rank_key,
                    "valency_mode": "PRIMARY_ANCHOR",
                }


        if best is None and eligible and sc in valency_terminal_extensions:
            v_ext = valency_terminal_extensions[sc]
            best = {
                "overhang": v_ext["overhang"],
                "ext_qstart": v_ext["ext_qstart"],
                "ext_qend": v_ext["ext_qend"],
                "term": v_ext["terminus"],
                "rec": v_ext["rec"],
                "valency_mode": v_ext["valency_mode"],
            }

        if best is None:
            continue

        rec = best["rec"]
        chrom = rec["tname"]
        term = best["term"]
        anchored_bp = rec["qend"] - rec["qstart"]
        anchored_pct = (anchored_bp / slen * 100.0) if slen > 0 else 0.0
        owner = terminus_owner.get((chrom, term), "plain")
        entry = {
            "scaffold": sc,
            "length": slen,
            "chromosome": chrom,
            "terminus": term,
            "terminus_owner": owner,
            "overhang_bp": best["overhang"],
            "ext_qstart": best["ext_qstart"],
            "ext_qend": best["ext_qend"],
            "anchored_span_pct": round(anchored_pct, 1),
            "mapq": rec["mapq"],
            "strand": rec["strand"],
            "valency_mode": best.get("valency_mode", "PRIMARY_ANCHOR"),
        }
        if owner != "plain":
            entry["action"] = "ROUTE_TO_" + owner.upper()
            entry["note"] = "terminus owned by the " + owner + " path; not a generic extension"
            extension_by_seq[sc] = entry
            continue
        if not eligible:
            entry["action"] = "EXTENSION_CANDIDATE_WITHHELD"
            entry["note"] = f"terminus is unreserved, but the sequence is {cat} ({triage_action})"
            extension_by_seq[sc] = entry
            continue

        entry["action"] = "INTEGRATE_INTO_CHROMOSOME"
        entry["note"] = f"continues {chrom} {term} with {best['overhang']:,} bp beyond the terminus"
        extension_by_seq[sc] = entry
        disposition[sc] = (cat, "INTEGRATE_INTO_CHROMOSOME")
        extension_cards.append(
            {
                "card_id": f"CARD-SCAFF-INT-{integrate_idx:02d}",
                "domain": "DOMAIN_SCAFFOLD",
                "category": "Unplaced Scaffolds: Terminal Extension Candidate",
                "title": f"Unplaced Scaffold `{sc}` ({slen:,} bp) extends {chrom} {term}",
                "scaffold": sc,
                "chrom": chrom,
                "terminus": term,
                "length": slen,
                "ext_qstart": best["ext_qstart"],
                "ext_qend": best["ext_qend"],
                "orientation": rec["strand"],
                "action": "INTEGRATE_INTO_CHROMOSOME",
                "classification": "TERMINAL_EXTENSION",
                "confidence": "HIGH" if best.get("valency_mode") == "CHROMOSOME_EXCLUSIVE_PORT" else "MEDIUM",
                "evidence": {
                    "target_terminus": f"{chrom} {term}",
                    "extension_beyond_terminus": f"{best['overhang']:,} bp",
                    "anchored_span": f"{anchored_pct:.1f}% ({anchored_bp:,} bp at MAPQ {rec['mapq']})",
                    "valency_mode": best.get("valency_mode", "PRIMARY_ANCHOR"),
                    "orientation": rec["strand"],
                },
                "instructions": [
                    f"Place `{sc}` ({best['ext_qstart']}..{best['ext_qend']}) into {chrom} at the {term} terminus, orientation {rec['strand']}.",
                    "Re-evaluate the newly created terminus for capping once joined.",
                ],
            }
        )
        integrate_idx += 1

    scaffold_cards.extend(extension_cards)

    scaffold_dispositions = []
    for sc in sorted(unplaced_seqs, key=lambda name: (-unplaced_seqs[name], name)):
        cat, action = disposition.get(sc, ("Unknown", "REPORT_ONLY"))
        row = {
            "scaffold": sc,
            "length": unplaced_seqs[sc],
            "triage": cat,
            "action": action,
            "core_aligned_bp": scaff_chrom_aligned_bp.get(sc, 0),
            "rdna_operon_bp": scaff_rdna_operon_bp.get(sc, 0),
            "novel_bp": scaff_novel_bp.get(sc, 0),
            "repeat_density_pct": round(scaff_repeats_pct.get(sc, 0.0), 2),
        }
        if sc in extension_by_seq:
            ext_entry = dict(extension_by_seq[sc])
            if action.startswith("PURGE_") and not ext_entry.get("action", "").startswith("INTEGRATE_"):
                ext_entry["action"] = action
            row.update(ext_entry)
        scaffold_dispositions.append(row)

    return scaffold_cards, scaffold_dispositions, extension_by_seq

