from __future__ import annotations

import csv
import gzip
import os
from typing import Any, Dict, List, Optional
from .models import EditedRegion, StructuralVariant


def parse_agp_regions(agp_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Parses AGP v2.1 file and indexes components and gap coordinates."""
    objects: Dict[str, List[Dict[str, Any]]] = {}
    if not os.path.exists(agp_path):
        return objects

    opener = gzip.open if agp_path.endswith(".gz") else open
    with opener(agp_path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            obj = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            part_num = int(parts[3])
            comp_type = parts[4]

            entry = {
                "object": obj,
                "start": start,
                "end": end,
                "part_num": part_num,
                "type": comp_type,
            }
            if comp_type == "W":
                entry["comp_id"] = parts[5]
                entry["comp_start"] = int(parts[6])
                entry["comp_end"] = int(parts[7])
                entry["orientation"] = parts[8] if len(parts) > 8 else "+"
            elif comp_type == "N":
                entry["gap_len"] = int(parts[5])
                entry["gap_type"] = parts[6]
                entry["linkage"] = parts[7]
            objects.setdefault(obj, []).append(entry)

    return objects


def parse_curation_audit(audit_path: str, agp_entries: Optional[Dict[str, List[Dict[str, Any]]]] = None) -> List[EditedRegion]:
    """
    Parses curation audit TSV (from CURATIONAGPGENERATOR) and maps edited regions
    to curated assembly coordinates using AGP where available.
    """
    regions: List[EditedRegion] = []
    if not os.path.exists(audit_path):
        return regions

    with open(audit_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get("sample", "")
            chrom = row.get("chrom", "")
            card_id = row.get("card_id", "")
            action = row.get("action", "")
            edit_type = row.get("edit_type", "")
            try:
                orig_start = int(row.get("orig_start", 0))
                orig_end = int(row.get("orig_end", 0))
                delta_bp = int(row.get("delta_bp", 0))
            except ValueError:
                orig_start = orig_end = delta_bp = 0

            cur_start = orig_start
            cur_end = orig_end

            agp_list = agp_entries.get(chrom, []) if agp_entries else []
            e_upper = edit_type.upper()
            act_upper = action.upper()

            # 1. 5' Telomere extension
            if ("5P" in e_upper or "5P" in act_upper) and "TELO" in (e_upper + act_upper):
                if agp_list and delta_bp > 0:
                    w_first = next((e for e in agp_list if e.get("type") == "W"), None)
                    ext_len = (w_first["end"] - w_first["start"] + 1) if w_first else delta_bp
                    cur_start = max(1, ext_len - 250)
                    cur_end = ext_len + 250
                elif delta_bp > 0:
                    cur_start = max(1, delta_bp - 250)
                    cur_end = delta_bp + 250
            # 2. 3' Telomere extension
            elif ("3P" in e_upper or "3P" in act_upper) and "TELO" in (e_upper + act_upper):
                if agp_list and delta_bp > 0:
                    w_last = next((e for e in reversed(agp_list) if e.get("type") == "W"), None)
                    j_pos = w_last["start"] if w_last else orig_end
                    cur_start = max(1, j_pos - 250)
                    cur_end = j_pos + 250
                elif delta_bp > 0:
                    cur_start = max(1, orig_end - 250)
                    cur_end = orig_end + delta_bp + 250
            # 3. Gap patch
            elif "PATCH" in e_upper or "PATCH" in act_upper:
                matched_entry = None
                if agp_list:
                    for e in agp_list:
                        cid = e.get("comp_id", "")
                        if card_id and card_id in cid:
                            matched_entry = e
                            break
                        if "_patch" in cid and chrom in cid:
                            if abs(e["start"] - orig_start) <= 500000 or (e["start"] <= orig_end + 100000 and e["end"] >= orig_start - 100000):
                                matched_entry = e
                                break
                if matched_entry:
                    cur_start = max(1, matched_entry["start"] - 200)
                    cur_end = matched_entry["end"] + 200
                else:
                    cur_start = max(1, orig_start - 200)
                    cur_end = max(cur_start + 1, orig_end + delta_bp + 200)
            # 4. Gap excision
            elif "EXCISE" in e_upper or "EXCISE" in act_upper:
                matched_junction = None
                if agp_list:
                    for i in range(len(agp_list) - 1):
                        e1, e2 = agp_list[i], agp_list[i + 1]
                        if e1.get("type") == "W" and e2.get("type") == "W":
                            s1 = e1.get("comp_end", 0)
                            s2 = e2.get("comp_start", 0)
                            if s1 <= orig_end + 500 and s2 >= orig_start - 500:
                                matched_junction = e1["end"]
                                break
                if matched_junction:
                    cur_start = max(1, matched_junction - 250)
                    cur_end = matched_junction + 250
                else:
                    cur_start = max(1, orig_start - 200)
                    cur_end = max(cur_start + 1, orig_end + delta_bp + 200)
            # 5. rDNA insertion / sequence insertion
            elif "RDNA" in e_upper or "INSERT" in e_upper:
                matched_entry = None
                if agp_list:
                    for e in agp_list:
                        cid = e.get("comp_id", "")
                        if card_id and card_id in cid:
                            matched_entry = e
                            break
                        if "rdna" in cid.lower():
                            matched_entry = e
                            break
                if matched_entry:
                    cur_start = max(1, matched_entry["start"] - 200)
                    cur_end = matched_entry["end"] + 200
                else:
                    cur_start = max(1, orig_start - 200)
                    cur_end = max(cur_start + 1, orig_end + delta_bp + 200)
            # 6. Legacy/default handling
            elif edit_type == "TELOMERE_CAPPING_5P":
                cur_start = max(1, delta_bp - 250)
                cur_end = delta_bp + 250
            elif edit_type in ("GAP_PATCH", "GAP_EXCISION"):
                cur_start = max(1, orig_start - 200)
                cur_end = max(cur_start + 1, orig_end + delta_bp + 200)

            regions.append(
                EditedRegion(
                    sample=sample,
                    chrom=chrom,
                    card_id=card_id,
                    action=action,
                    edit_type=edit_type,
                    curated_start=cur_start,
                    curated_end=cur_end,
                    orig_start=orig_start,
                    orig_end=orig_end,
                    delta_bp=delta_bp,
                    extra=dict(row)
                )
            )
    return regions


def parse_sniffles_vcf(vcf_path: str) -> List[StructuralVariant]:
    """Parses structural variants called by Sniffles against curated assembly."""
    svs: List[StructuralVariant] = []
    if not os.path.exists(vcf_path):
        return svs

    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            sv_id = parts[2]
            qual = float(parts[5]) if parts[5] != "." else 0.0
            filter_status = parts[6]

            info_items = parts[7].split(";")
            info = {}
            for item in info_items:
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = True

            sv_type = info.get("SVTYPE", "SV")
            end_pos = int(info.get("END", pos))
            try:
                sv_len = abs(int(info.get("SVLEN", end_pos - pos)))
            except ValueError:
                sv_len = abs(end_pos - pos)

            try:
                support = int(info.get("SUPPORT", 0))
            except ValueError:
                support = 0

            try:
                vaf = float(info.get("VAF", 0.0))
            except ValueError:
                vaf = 0.0

            chr2 = info.get("CHR2", chrom)

            svs.append(
                StructuralVariant(
                    sv_id=sv_id,
                    chrom=chrom,
                    start=pos,
                    end=end_pos,
                    sv_type=sv_type,
                    sv_len=sv_len,
                    support=support,
                    vaf=vaf,
                    filter_status=filter_status,
                    qual=qual,
                    chr2=chr2,
                    info=info
                )
            )
    return svs
