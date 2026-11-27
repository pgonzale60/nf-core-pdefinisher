#!/usr/bin/env python3
"""Sequence-family analysis of inward-facing telomere flanks; families may merge ends.

Ported from PDE_origin/scripts/telomeric_read_fishing/end_clusters.py for
self-contained execution inside nf-core/pdefinisher.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
from array import array
from collections import Counter, defaultdict
from pathlib import Path


def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def fasta(path):
    with op(path) as handle:
        name, chunks = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def table(path):
    with op(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_table(path, rows, fields):
    with op(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def terminal_groups(hits, length, cfg):
    groups = []
    for strand in ("-", "+"):
        for h in sorted((h for h in hits if h["strand"] == strand), key=lambda h: h["start"]):
            if groups and groups[-1]["strand"] == strand and h["start"] <= groups[-1]["end"] + cfg["tract_join_gap_bp"]:
                groups[-1]["end"] = max(groups[-1]["end"], h["end"])
                groups[-1]["hits"].append(h)
            else:
                groups.append(dict(start=h["start"], end=h["end"], strand=strand, hits=[h]))
    selected = {}
    left = [g for g in groups if g["strand"] == "-" and g["start"] <= cfg["terminal_gap_bp"]]
    right = [g for g in groups if g["strand"] == "+" and length - g["end"] <= cfg["terminal_gap_bp"]]
    if left:
        selected["5p"] = min(left, key=lambda g: (g["start"], -g["end"]))
    if right:
        selected["3p"] = min(right, key=lambda g: (length - g["end"], g["start"]))
    return selected


END_RECORD_FIELDS = [
    "end_id",
    "read_id",
    "movie",
    "read_end",
    "read_length",
    "tract_start",
    "tract_end",
    "tract_span",
    "tract_parts",
    "end_gap",
    "strict_gap",
    "flank_length",
    "prefix_length",
    "status",
]

READ_FATE_FIELDS = [
    "read_id",
    "movie",
    "read_length",
    "ncrf_hits",
    "terminal_ends",
    "status",
]


def extract(args, cfg):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    hits = defaultdict(list)
    for row in table(args.ncrf):
        row.update({key: int(row[key]) for key in ("start", "end", "seqLen", "m", "mm", "i", "d")})
        hits[row["seq"]].append(row)
    records, fates, internal = [], [], []
    with gzip.open(out / "inward_full.fasta.gz", "wt") as full, (out / "prefixes.fasta").open("w") as prefixes:
        for name, seq in fasta(args.reads):
            read_hits = hits.pop(name, [])
            selected = terminal_groups(read_hits, len(seq), cfg)
            used = {id(h) for group in selected.values() for h in group["hits"]}
            for h in read_hits:
                if id(h) not in used:
                    internal.append(
                        dict(
                            read_id=name,
                            start=h["start"],
                            end=h["end"],
                            strand=h["strand"],
                            read_length=len(seq),
                            match_ratio=h["mRatio"],
                            status="nonterminal_or_other_orientation",
                        )
                    )
            fate = "no_ncrf_hit" if not read_hits else "nonterminal_only" if not selected else "terminal"
            fates.append(
                dict(
                    read_id=name,
                    movie=name.split("/")[0],
                    read_length=len(seq),
                    ncrf_hits=len(read_hits),
                    terminal_ends=len(selected),
                    status=fate,
                )
            )
            for end, group in selected.items():
                inner_left = selected["5p"]["end"] if "5p" in selected else 0
                inner_right = selected["3p"]["start"] if "3p" in selected else len(seq)
                core = seq[inner_left:inner_right] if inner_left < inner_right else ""
                flank = core if end == "5p" else revcomp(core)
                end_id = f"E{len(records) + 1:08d}"
                gap = group["start"] if end == "5p" else len(seq) - group["end"]
                records.append(
                    dict(
                        end_id=end_id,
                        read_id=name,
                        movie=name.split("/")[0],
                        read_end=end,
                        read_length=len(seq),
                        tract_start=group["start"],
                        tract_end=group["end"],
                        tract_span=group["end"] - group["start"],
                        tract_parts=len(group["hits"]),
                        end_gap=gap,
                        strict_gap=int(gap <= cfg["strict_terminal_gap_bp"]),
                        flank_length=len(flank),
                        prefix_length=min(len(flank), cfg["prefix_bp"]),
                        status="clusterable" if len(flank) >= cfg["min_flank_bp"] else "short_or_no_flank",
                    )
                )
                if len(flank) >= cfg["min_flank_bp"]:
                    full.write(f">{end_id}\n{flank}\n")
                    prefixes.write(f">{end_id}\n{flank[: cfg['prefix_bp']]}\n")
    if hits:
        raise ValueError("NCRF contains read IDs absent from FASTA")
    write_table(out / "end_records.tsv.gz", records, END_RECORD_FIELDS)
    write_table(out / "read_fate.tsv.gz", fates, READ_FATE_FIELDS)
    write_table(
        out / "nonterminal_hits.tsv.gz",
        internal,
        ["read_id", "start", "end", "strand", "read_length", "match_ratio", "status"],
    )
    # Recover per-movie input exposure only when header/log cardinality is unambiguous.
    exposure = []
    if args.headers and args.fishing_log:
        headers_path = Path(args.headers)
        log_path = Path(args.fishing_log)
        if headers_path.exists() and log_path.exists():
            movies = []
            for line in headers_path.read_text().splitlines():
                if line.startswith("@RG"):
                    for value in line.split("\t"):
                        if value.startswith("PU:") and value[3:] not in movies:
                            movies.append(value[3:])
            counts = re.findall(r"\[M::bam2fq_mainloop\] processed (\d+) reads", log_path.read_text())
            if len(movies) == len(counts):
                exposure = [{"movie": movie, "input_reads": count} for movie, count in zip(movies, counts)]
    write_table(out / "movie_exposure.tsv", exposure, ["movie", "input_reads"])
    print("EXTRACT", len(fates), "reads", len(records), "end observations", Counter(r["status"] for r in records))


def paf(path):
    with op(path) as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12 or f[5] == "*":
                continue
            yield dict(
                q=f[0],
                ql=int(f[1]),
                qs=int(f[2]),
                qe=int(f[3]),
                strand=f[4],
                t=f[5],
                tl=int(f[6]),
                ts=int(f[7]),
                te=int(f[8]),
                matches=int(f[9]),
                block=int(f[10]),
                mapq=int(f[11]),
                tags={tag[:2]: tag[5:] for tag in f[12:]},
                line=line,
            )


def overlap_ok(p, cfg, identity=None):
    return (
        p["q"] != p["t"]
        and p["strand"] == "+"
        and p["block"] > 0
        and p["matches"] / p["block"] >= (identity or cfg["cluster_identity"])
        and min(p["qe"] - p["qs"], p["te"] - p["ts"]) >= cfg["min_flank_bp"] * cfg["cluster_query_coverage"]
        and min(p["qe"] - p["qs"], p["te"] - p["ts"]) / min(p["ql"], p["tl"]) >= cfg["cluster_query_coverage"]
        and max(p["qs"], p["ts"]) <= cfg["anchor_offset_bp"]
    )


class UnionFind:
    def __init__(self, names):
        self.parent = {name: name for name in names}

    def root(self, name):
        while self.parent[name] != name:
            self.parent[name] = self.parent[self.parent[name]]
            name = self.parent[name]
        return name

    def join(self, a, b):
        a, b = self.root(a), self.root(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def cluster(args, cfg):
    out = Path(args.out)
    seqs = dict(fasta(out / "prefixes.fasta"))
    uf = UnionFind(seqs)
    degree = Counter()
    sensitivity = {x: UnionFind(seqs) for x in cfg["identity_sensitivity"]}
    edges = 0
    for p in paf(args.overlaps):
        for identity, obj in sensitivity.items():
            if overlap_ok(p, cfg, identity):
                obj.join(p["q"], p["t"])
        if overlap_ok(p, cfg):
            uf.join(p["q"], p["t"])
            degree.update((p["q"], p["t"]))
            edges += 1
    groups = defaultdict(list)
    for name in sorted(seqs):
        groups[uf.root(name)].append(name)
    memberships = []
    with (out / "representatives.fasta").open("w") as handle:
        for i, members in enumerate(sorted(groups.values(), key=lambda m: m[0]), 1):
            cid = f"C{i:06d}"
            seed = max(members, key=lambda name: (len(seqs[name]), degree[name], name))
            handle.write(f">{cid}\n{seqs[seed]}\n")
            for name in members:
                memberships.append(
                    dict(
                        end_id=name,
                        cluster_id=cid,
                        seed_end_id=seed,
                        prefix_length=len(seqs[name]),
                        overlap_degree=degree[name],
                    )
                )
    write_table(
        out / "membership.tsv.gz",
        memberships,
        ["end_id", "cluster_id", "seed_end_id", "prefix_length", "overlap_degree"],
    )
    write_table(
        out / "identity_sensitivity.tsv",
        [{"identity": x, "sequence_families": len({u.root(n) for n in seqs})} for x, u in sensitivity.items()],
        ["identity", "sequence_families"],
    )
    print("CLUSTER", len(seqs), "observations", len(groups), "sequence families", edges, "accepted directed overlaps")


def assigned(args, cfg):
    out = Path(args.out)
    assignment = {r["end_id"]: r["cluster_id"] for r in table(out / "membership.tsv.gz")}
    best = {}
    for p in paf(args.alignments):
        if assignment.get(p["q"]) != p["t"] or p["strand"] != "+":
            continue
        if (p["qe"] - p["qs"]) / p["ql"] < cfg["cluster_query_coverage"]:
            continue
        key = (p["matches"], -p["block"], p["mapq"])
        if p["q"] not in best or key > best[p["q"]][0]:
            best[p["q"]] = (key, p["line"])
    with (out / "assigned.paf").open("w") as handle:
        for name in sorted(best):
            handle.write(best[name][1])
    print("ASSIGNED", len(best), "of", len(assignment), "end observations mapped to their family representative")


def cs_calls(p):
    """Return observed substitutions and deletion intervals in target coordinates."""
    pos = p["ts"]
    subs = {}
    deletions = []
    for match in re.finditer(
        r":(\d+)|\*([a-zA-Z])([a-zA-Z])|\+([a-zA-Z]+)|-([a-zA-Z]+)|=([a-zA-Z]+)", p["tags"].get("cs", "")
    ):
        n, ref, alt, insertion, deletion, equal = match.groups()
        if n:
            pos += int(n)
        elif ref:
            subs[pos] = alt.upper()
            pos += 1
        elif deletion:
            deletions.append((pos, pos + len(deletion)))
            pos += len(deletion)
        elif equal:
            pos += len(equal)
    if pos != p["te"]:
        raise ValueError(f"cs span differs from PAF span for {p['q']}")
    return subs, deletions


FAMILY_FIELDS = [
    "cluster_id",
    "end_observations",
    "unique_reads",
    "strict_gap_reads",
    "movies",
    "rdna_proximal_observations",
    "rdna_any_observations",
    "candidate_reference_loci",
    "uniquely_supported_loci",
    "multiplicity_status",
    "profiled_observations",
    "callable_sites",
    "mean_pairwise_diversity",
    "informative_sites",
    "supported_patterns",
    "mean_prefix_depth",
]


def summarize(args, cfg):
    out = Path(args.out)
    records = {r["end_id"]: r for r in table(out / "end_records.tsv.gz")}
    members = defaultdict(list)
    assignment = {}
    for r in table(out / "membership.tsv.gz"):
        assignment[r["end_id"]] = r["cluster_id"]
        members[r["cluster_id"]].append(r["end_id"])
    refs = dict(fasta(out / "representatives.fasta"))
    exposure_path = out / "movie_exposure.tsv"
    exposure = {r["movie"]: int(r["input_reads"]) for r in table(exposure_path)} if exposure_path.exists() else {}
    mapped = defaultdict(list)
    diff = {c: array("i", [0]) * (len(s) + 1) for c, s in refs.items()}
    alt_counts = defaultdict(Counter)
    assigned_path = Path(args.assigned) if args.assigned else (out / "assigned.paf")
    for p in paf(assigned_path):
        subs, deletions = cs_calls(p)
        cid = p["t"]
        mapped[cid].append((p, subs, deletions))
        diff[cid][p["ts"]] += 1
        diff[cid][p["te"]] -= 1
        for a, b in deletions:
            diff[cid][a] -= 1
            diff[cid][b] += 1
        for pos, base in subs.items():
            alt_counts[cid, pos][base] += 1
    variant_rows, hap_rows, diversity = [], [], {}
    for cid, seq in refs.items():
        depth = []
        running = 0
        for change in diff[cid][:-1]:
            running += change
            depth.append(running)
        sites = []
        pi = 0.0
        callable_sites = 0
        for pos, n in enumerate(depth):
            if n < cfg["variant_min_depth"]:
                continue
            alleles = alt_counts.get((cid, pos), Counter()).copy()
            alleles[seq[pos]] += n - sum(alleles.values())
            counts = sorted(alleles.values(), reverse=True)
            pi += (1 - sum((v / n) ** 2 for v in alleles.values())) * n / (n - 1)
            callable_sites += 1
            if (
                len(counts) > 1
                and counts[1] >= cfg["variant_min_minor_count"]
                and counts[1] / n >= cfg["variant_min_minor_fraction"]
            ):
                sites.append(pos)
                variant_rows.append(
                    dict(
                        cluster_id=cid,
                        position=pos,
                        depth=n,
                        alleles=json.dumps(dict(alleles), sort_keys=True),
                        minor_fraction=counts[1] / n,
                    )
                )
        pattern_reads = defaultdict(set)
        for p, subs, deletions in mapped[cid]:
            if sites and all(p["ts"] <= pos < p["te"] and not any(a <= pos < b for a, b in deletions) for pos in sites):
                pattern = "".join(subs.get(pos, seq[pos]) for pos in sites)
                pattern_reads[pattern].add(records[p["q"]]["read_id"])
        patterns = {pattern: len(reads) for pattern, reads in pattern_reads.items()}
        for pattern, count in patterns.items():
            hap_rows.append(
                dict(
                    cluster_id=cid,
                    pattern_id=hashlib.sha256(pattern.encode()).hexdigest()[:16],
                    reads=count,
                    alleles=pattern,
                )
            )
        diversity[cid] = dict(
            profiled_observations=len(mapped[cid]),
            callable_sites=callable_sites,
            mean_pairwise_diversity=pi / callable_sites if callable_sites else "",
            informative_sites=len(sites),
            supported_patterns=sum(n >= cfg["haplotype_min_reads"] for n in patterns.values()),
            mean_prefix_depth=sum(depth) / len(depth) if depth else 0,
        )
    # Full inward flanks may reach unique sequence beyond the shared clustering prefix.
    placements = []
    locus_reads = defaultdict(set)
    locus_unique = defaultdict(set)
    for p in paf(args.genome):
        cid = assignment[p["q"]]
        if (
            p["block"] < cfg["mapping_min_bp"]
            or p["matches"] / p["block"] < cfg["mapping_identity"]
            or p["qs"] > cfg["anchor_offset_bp"]
        ):
            continue
        coordinate = p["ts"] - p["qs"] if p["strand"] == "+" else p["te"] + p["qs"]
        coordinate = max(0, min(p["tl"], coordinate))
        side = (
            "left"
            if p["strand"] == "+" and coordinate <= cfg["reference_end_window_bp"]
            else "right"
            if p["strand"] == "-" and p["tl"] - coordinate <= cfg["reference_end_window_bp"]
            else "internal"
        )
        read = records[p["q"]]["read_id"]
        placements.append(
            dict(
                end_id=p["q"],
                read_id=read,
                cluster_id=cid,
                reference=p["t"],
                strand=p["strand"],
                coordinate=coordinate,
                native_side=side,
                mapq=p["mapq"],
                alignment_identity=p["matches"] / p["block"],
                query_start=p["qs"],
                query_end=p["qe"],
                primary=int(p["tags"].get("tp") == "P"),
            )
        )
    grouped = defaultdict(list)
    for placement in placements:
        grouped[placement["cluster_id"], placement["reference"], placement["strand"]].append(placement)
    locus_groups = {}
    for key, group in grouped.items():
        locus_id, previous = 0, None
        for placement in sorted(group, key=lambda p: p["coordinate"]):
            if previous is None or placement["coordinate"] - previous > cfg["locus_bin_bp"]:
                locus_id += 1
            previous = placement["coordinate"]
            locus = (*key, locus_id)
            locus_groups.setdefault(locus, []).append(placement)
            locus_reads[locus].add(placement["read_id"])
            if placement["mapq"] >= cfg["unique_mapq"] and placement["primary"]:
                locus_unique[locus].add(placement["read_id"])
    locus_rows = []
    for locus, group in locus_groups.items():
        coordinates = sorted(p["coordinate"] for p in group)
        locus_rows.append(
            dict(
                cluster_id=locus[0],
                reference=locus[1],
                strand=locus[2],
                locus_group=locus[3],
                median_coordinate=coordinates[len(coordinates) // 2],
                native_sides=",".join(sorted({p["native_side"] for p in group})),
                all_candidate_reads=len(locus_reads[locus]),
                unique_mapping_reads=len(locus_unique[locus]),
                ambiguity_status="unique_support" if len(locus_unique[locus]) >= 2 else "unresolved_mapping",
            )
        )
    # rDNA sequence identity is an annotation, not a requirement for inclusion.
    rdna_any = set()
    rdna_tip = set()
    if args.rdna and Path(args.rdna).exists():
        for p in paf(args.rdna):
            if p["block"] >= cfg["rdna_min_alignment_bp"] and p["matches"] / p["block"] >= cfg["rdna_min_identity"]:
                rdna_any.add(p["q"])
                if p["qs"] <= cfg["anchor_offset_bp"]:
                    rdna_tip.add(p["q"])
    rows = []
    library_rows = []
    for cid, names in members.items():
        read_set = {records[n]["read_id"] for n in names}
        by_movie = defaultdict(set)
        for name in names:
            by_movie[records[name]["movie"]].add(records[name]["read_id"])
        for movie, read_ids in by_movie.items():
            library_rows.append(
                dict(
                    cluster_id=cid,
                    movie=movie,
                    unique_reads=len(read_ids),
                    input_reads=exposure.get(movie, ""),
                    support_per_million_input_reads=len(read_ids) * 1e6 / exposure[movie]
                    if exposure.get(movie)
                    else "",
                )
            )
        loci = [key for key in locus_reads if key[0] == cid]
        rows.append(
            dict(
                cluster_id=cid,
                end_observations=len(names),
                unique_reads=len(read_set),
                strict_gap_reads=len({records[n]["read_id"] for n in names if int(records[n]["strict_gap"])}),
                movies=len(by_movie),
                rdna_proximal_observations=sum(n in rdna_tip for n in names),
                rdna_any_observations=sum(n in rdna_any for n in names),
                candidate_reference_loci=len(loci),
                uniquely_supported_loci=sum(len(locus_unique[k]) >= 2 for k in loci),
                multiplicity_status="unresolved_sequence_family",
                **diversity[cid],
            )
        )
    fai = Path(args.reference + ".fai")
    if fai.exists():
        with fai.open() as handle:
            lengths = {f[0]: int(f[1]) for f in (line.split() for line in handle)}
    else:
        lengths = {name: len(seq) for name, seq in fasta(args.reference)}
    inventory = []
    for ref, length in lengths.items():
        for side in ("left", "right"):
            supported = [p for p in placements if p["reference"] == ref and p["native_side"] == side]
            inventory.append(
                dict(
                    reference=ref,
                    side=side,
                    sequence_length=length,
                    candidate_reads=len({p["read_id"] for p in supported}),
                    unique_mapping_reads=len(
                        {p["read_id"] for p in supported if p["mapq"] >= cfg["unique_mapq"] and p["primary"]}
                    ),
                    sequence_families=",".join(sorted({p["cluster_id"] for p in supported})),
                    status="supported" if supported else "no_assigned_terminal_read",
                )
            )
    write_table(out / "families.tsv", rows, FAMILY_FIELDS)
    write_table(
        out / "library_support.tsv",
        library_rows,
        ["cluster_id", "movie", "unique_reads", "input_reads", "support_per_million_input_reads"],
    )
    write_table(
        out / "variant_sites.tsv.gz", variant_rows, ["cluster_id", "position", "depth", "alleles", "minor_fraction"]
    )
    write_table(out / "haplotype_patterns.tsv.gz", hap_rows, ["cluster_id", "pattern_id", "reads", "alleles"])
    write_table(
        out / "placements.tsv.gz",
        placements,
        [
            "end_id",
            "read_id",
            "cluster_id",
            "reference",
            "strand",
            "coordinate",
            "native_side",
            "mapq",
            "alignment_identity",
            "query_start",
            "query_end",
            "primary",
        ],
    )
    write_table(
        out / "reference_ends.tsv",
        inventory,
        [
            "reference",
            "side",
            "sequence_length",
            "candidate_reads",
            "unique_mapping_reads",
            "sequence_families",
            "status",
        ],
    )
    write_table(
        out / "locus_support.tsv",
        locus_rows,
        [
            "cluster_id",
            "reference",
            "strand",
            "locus_group",
            "median_coordinate",
            "native_sides",
            "all_candidate_reads",
            "unique_mapping_reads",
            "ambiguity_status",
        ],
    )
    manifest = {
        "parameters": cfg,
        "input_reads": args.reads,
        "ncrf_table": args.ncrf,
        "reference": args.reference,
        "interpretation": "Clusters are sequence families, not physical-end counts. Consensus does not replace original reads.",
    }
    if args.rdna_reference:
        manifest["rdna_monomer"] = args.rdna_reference
    Path(out / "MANIFEST.yaml").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        "SUMMARY",
        len(rows),
        "families",
        len(placements),
        "candidate placements",
        len(inventory),
        "assembly sequence ends",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["extract", "cluster", "assigned", "summarize"])
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    for name in (
        "reads",
        "ncrf",
        "headers",
        "fishing-log",
        "overlaps",
        "alignments",
        "assigned",
        "genome",
        "rdna",
        "reference",
        "rdna-reference",
    ):
        parser.add_argument("--" + name)
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    {"extract": extract, "cluster": cluster, "assigned": assigned, "summarize": summarize}[args.action](args, cfg)


if __name__ == "__main__":
    main()
