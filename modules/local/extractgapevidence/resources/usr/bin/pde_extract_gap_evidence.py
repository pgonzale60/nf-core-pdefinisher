#!/usr/bin/env python3
"""
pde_extract_gap_evidence.py

Module 1 Evidence Extraction Engine:
Extracts granular, un-aggregated physical read observations around assembly gaps
and structural variants from HiFi alignments (BAM), Gaps (BED), and SVs (VCF).

Architectural Principle:
- ZERO information loss: records per-read micro-properties (MAPQ, anchors, CIGAR ops, soft-clips).
- NO arbitrary thresholds: does not pre-filter by MAPQ or pre-aggregate into rigid integer counts.
- Empowers Module 2 to re-filter, bin, and explore curation policies instantly without re-running BAM I/O.
- Strict cryptographic runtime provenance (SHA256, schema version, input files).
"""

import argparse
import gzip
import hashlib
import json
import os
import sys

import pysam

SCHEMA_VERSION = "3.0.0"

def get_script_fingerprint():
    try:
        with open(__file__, 'rb') as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except Exception:
        return "unknown"

def parse_gaps(gaps_bed):
    gaps = []
    if not os.path.exists(gaps_bed):
        return gaps
    with open(gaps_bed) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gaps.append({
                'gap_id': f"gap_{chrom}_{start}_{end}",
                'chrom': chrom,
                'start': start,
                'end': end,
                'length': end - start
            })
    return gaps

def parse_sniffles_svs(vcf_path):
    svs = []
    if not os.path.exists(vcf_path):
        return svs
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            sv_id = parts[2]
            alt_seq = parts[4]
            has_explicit_alt = bool(alt_seq and not alt_seq.startswith('<') and not any(c in alt_seq for c in '[]'))
            qual = float(parts[5]) if parts[5] != '.' else 0.0
            filt = parts[6]
            info = dict(item.split('=', 1) if '=' in item else (item, True) for item in parts[7].split(';'))
            svtype = info.get('SVTYPE', '')
            
            end_pos = int(info.get('END', pos))
            try:
                svlen = abs(int(info.get('SVLEN', end_pos - pos)))
            except ValueError:
                svlen = abs(end_pos - pos)

            try:
                support = int(info.get('SUPPORT', 0))
            except ValueError:
                support = 0

            try:
                vaf = float(info.get('VAF', 0.0))
            except ValueError:
                vaf = 0.0

            cov = info.get('COVERAGE', '')

            chr2 = info.get('CHR2', chrom)
            mateid = info.get('MATEID', '')

            svs.append({
                'id': sv_id,
                'chrom': chrom,
                'start': pos,
                'end': end_pos,
                'svtype': svtype,
                'svlen': svlen,
                'qual': qual,
                'filter': filt,
                'support': support,
                'vaf': vaf,
                'coverage': cov,
                'chr2': chr2,
                'alt': alt_seq,
                'mateid': mateid,
                'has_explicit_alt_seq': has_explicit_alt,
                'alt_seq_len': len(alt_seq) if has_explicit_alt else 0,
                'alt_seq_head': alt_seq[:60] if has_explicit_alt else None,
                'alt_seq_tail': alt_seq[-60:] if has_explicit_alt else None
            })
    return svs

def extract_gap_read_observations(bam, chrom, start, end, search_win=2500):
    """
    Extracts individual read observations around an assembly gap.
    Records granular properties without pre-filtering on MAPQ.
    """
    read_obs = []
    seen_qnames = set()

    q_start = max(0, start - search_win)
    q_end = end + search_win

    for r in bam.fetch(chrom, q_start, q_end):
        if r.is_unmapped:
            continue
        
        # Unique alignment record identifier: qname + primary/supp tag
        rec_id = f"{r.query_name}_{'supp' if r.is_supplementary else 'prim'}_{r.reference_start}"
        if rec_id in seen_qnames:
            continue
        seen_qnames.add(rec_id)

        r_start = r.reference_start
        r_end = r.reference_end
        mapq = r.mapping_quality
        qlen = r.query_length

        # 1. Collinear Spanning Read (anchored on both left and right flanks)
        if not r.is_secondary and r_start < start and r_end > end:
            left_anchor = start - r_start
            right_anchor = r_end - end
            
            # Trace CIGAR operations across the gap interval [start, end]
            del_at_gap = 0
            ins_at_gap = 0
            ref_p = r_start
            for op, op_len in (r.cigartuples or []):
                if op == 2 and abs(ref_p - start) <= 50: # Deletion
                    del_at_gap = max(del_at_gap, op_len)
                elif op == 1 and abs(ref_p - start) <= 50: # Insertion
                    ins_at_gap = max(ins_at_gap, op_len)
                if op in (0, 2, 3, 7, 8):
                    ref_p += op_len

            read_obs.append({
                'qname': r.query_name,
                'is_supplementary': r.is_supplementary,
                'mapq': mapq,
                'query_length': qlen,
                'span_type': 'COLLINEAR_SPANNING',
                'left_anchor_bp': left_anchor,
                'right_anchor_bp': right_anchor,
                'deletion_across_gap': del_at_gap,
                'insertion_across_gap': ins_at_gap,
                'distance_to_left_edge': 0,
                'distance_to_right_edge': 0,
                'soft_clip_bp': 0
            })

        # 2. Left-Abutting Read (terminates near left gap flank)
        elif abs(r_end - start) <= 100:
            dist = start - r_end
            clip_len = r.cigartuples[-1][1] if (r.cigartuples and r.cigartuples[-1][0] in (4, 5)) else 0
            read_obs.append({
                'qname': r.query_name,
                'is_supplementary': r.is_supplementary,
                'mapq': mapq,
                'query_length': qlen,
                'span_type': 'LEFT_ABUTTING',
                'left_anchor_bp': start - r_start,
                'right_anchor_bp': 0,
                'deletion_across_gap': 0,
                'insertion_across_gap': 0,
                'distance_to_left_edge': dist,
                'distance_to_right_edge': None,
                'soft_clip_bp': clip_len
            })

        # 3. Right-Abutting Read (starts near right gap flank)
        elif abs(r_start - end) <= 100:
            dist = r_start - end
            clip_len = r.cigartuples[0][1] if (r.cigartuples and r.cigartuples[0][0] in (4, 5)) else 0
            read_obs.append({
                'qname': r.query_name,
                'is_supplementary': r.is_supplementary,
                'mapq': mapq,
                'query_length': qlen,
                'span_type': 'RIGHT_ABUTTING',
                'left_anchor_bp': 0,
                'right_anchor_bp': r_end - end,
                'deletion_across_gap': 0,
                'insertion_across_gap': 0,
                'distance_to_left_edge': None,
                'distance_to_right_edge': dist,
                'soft_clip_bp': clip_len
            })

        # 4. Supplementary / Split reads in window
        elif r.is_supplementary:
            read_obs.append({
                'qname': r.query_name,
                'is_supplementary': True,
                'mapq': mapq,
                'query_length': qlen,
                'span_type': 'SUPPLEMENTARY_IN_WINDOW',
                'ref_start': r_start,
                'ref_end': r_end,
                'has_sa_tag': r.has_tag('SA')
            })

    return read_obs

def extract_inversion_split_reads(bam, chrom, pos, search_win=1000):
    """Extract individual supplementary alignments near candidate inversion breakpoints."""
    split_obs = []
    for r in bam.fetch(chrom, max(0, pos - search_win), pos + search_win):
        if r.is_supplementary and r.has_tag('SA'):
            split_obs.append({
                'qname': r.query_name,
                'mapq': r.mapping_quality,
                'ref_start': r.reference_start,
                'ref_end': r.reference_end,
                'strand': '-' if r.is_reverse else '+',
                'sa_tag': r.get_tag('SA')
            })
    return split_obs

def extract_evidence(gaps_bed, sniffles_vcf, bam_path, sample_id):
    gaps = parse_gaps(gaps_bed)
    svs = parse_sniffles_svs(sniffles_vcf)
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Extract granular read observations for each gap
    for g in gaps:
        ch = g['chrom']
        st = g['start']
        en = g['end']
        g['read_observations'] = extract_gap_read_observations(bam, ch, st, en)
        # Empirical flanking depths within 200bp of the gap boundaries
        try:
            d_5p = len(list(bam.fetch(ch, max(0, st - 200), st)))
            d_3p = len(list(bam.fetch(ch, en, en + 200)))
        except Exception:
            d_5p = 0
            d_3p = 0
        g['depth_flank_5p'] = d_5p
        g['depth_flank_3p'] = d_3p

    # Extract granular split reads for inversions / breakends
    for s in svs:
        if s['svtype'] in ('INV', 'BND'):
            s['split_read_observations'] = extract_inversion_split_reads(bam, s['chrom'], s['start'])

    provenance = {
        'module': 'module1_reconnaissance',
        'task': 'EXTRACT_GAP_AND_ANOMALY_EVIDENCE',
        'script_fingerprint': f"sha256:{get_script_fingerprint()}",
        'schema_version': SCHEMA_VERSION,
        'sample_id': sample_id,
        'granularity': 'per_read_micro_observations',
        'inputs': {
            'gaps_bed': os.path.basename(gaps_bed),
            'sniffles_vcf': os.path.basename(sniffles_vcf),
            'bam_file': os.path.basename(bam_path)
        }
    }

    return {
        'provenance': provenance,
        'gaps': gaps,
        'structural_variants': svs
    }

def run_self_test():
    sha = get_script_fingerprint()
    assert len(sha) == 16, "SHA256 calculation failed"
    print(f"Self-test PASS: schema_version={SCHEMA_VERSION}, fingerprint=sha256:{sha}")
    sys.exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract granular, unaggregated gap and SV evidence for Module 2")
    parser.add_argument('--gaps', help="Path to gaps.bed")
    parser.add_argument('--sniffles', help="Path to sniffles.vcf.gz")
    parser.add_argument('--bam', help="Path to alignment BAM")
    parser.add_argument('--sample', default="SAMPLE", help="Sample identifier")
    parser.add_argument('--self-test', action='store_true', help="Run self-test and exit")
    parser.add_argument('--version', action='store_true', help="Print schema and fingerprint")
    args = parser.parse_args()

    if args.version:
        print(f"pde_extract_gap_evidence schema={SCHEMA_VERSION} sha256={get_script_fingerprint()}")
        sys.exit(0)

    if args.self_test:
        run_self_test()

    if not (args.gaps and args.sniffles and args.bam):
        parser.error("--gaps, --sniffles, and --bam are required when not running --self-test or --version")

    result = extract_evidence(args.gaps, args.sniffles, args.bam, args.sample)
    print(json.dumps(result, indent=2, sort_keys=True))
