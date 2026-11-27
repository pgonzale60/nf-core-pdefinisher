import gzip
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import gfapdefinisher_handoff
import telomere_end_clusters


def _write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_telomere_end_clusters_builds_handoff_compatible_archive(tmp_path):
    cfg = json.loads((Path(__file__).resolve().parents[1] / "assets" / "telomeric_end_clusters.json").read_text())
    cfg["min_flank_bp"] = 40
    cfg["mapping_min_bp"] = 40
    cfg["prefix_bp"] = 60
    cfg["variant_min_depth"] = 2
    cfg["variant_min_minor_count"] = 1
    cfg["haplotype_min_reads"] = 1

    flank_a = "ACGTTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    flank_b = "ACGTTAGCATCGATCGATCGATCAATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    telo_minus = "GCCTAA" * 6
    read1_seq = telo_minus + flank_a
    read2_seq = telo_minus + flank_b

    reads_gz = tmp_path / "reads.fasta.gz"
    _write_gz(reads_gz, f">m1/1/ccs\n{read1_seq}\n>m1/2/ccs\n{read2_seq}\n")

    ncrf_gz = tmp_path / "ncrf_summary.tsv.gz"
    _write_gz(
        ncrf_gz,
        "line\tmotif\tseq\tstart\tend\tstrand\tseqLen\tquerybp\tmRatio\tm\tmm\ti\td\n"
        f"1\tTTAGGC\tm1/1/ccs\t0\t{len(telo_minus)}\t-\t{len(read1_seq)}\t36\t1.0\t36\t0\t0\t0\n"
        f"2\tTTAGGC\tm1/2/ccs\t0\t{len(telo_minus)}\t-\t{len(read2_seq)}\t36\t1.0\t36\t0\t0\t0\n",
    )

    archive = tmp_path / "end_archive"
    telomere_end_clusters.extract(
        SimpleNamespace(out=str(archive), reads=str(reads_gz), ncrf=str(ncrf_gz), headers=None, fishing_log=None),
        cfg,
    )

    overlaps = tmp_path / "overlaps.paf.gz"
    _write_gz(
        overlaps,
        f"E00000001\t{len(flank_a)}\t0\t{len(flank_a)}\t+\tE00000002\t{len(flank_b)}\t0\t{len(flank_b)}\t59\t60\t60\ttp:A:P\tcs:Z::23*ag:36\n",
    )
    telomere_end_clusters.cluster(SimpleNamespace(out=str(archive), overlaps=str(overlaps)), cfg)

    members = tmp_path / "members.paf.gz"
    _write_gz(
        members,
        f"E00000001\t{len(flank_a)}\t0\t{len(flank_a)}\t+\tC000001\t{len(flank_b)}\t0\t{len(flank_b)}\t59\t60\t60\ttp:A:P\tcs:Z::23*ag:36\n"
        f"E00000002\t{len(flank_b)}\t0\t{len(flank_b)}\t+\tC000001\t{len(flank_b)}\t0\t{len(flank_b)}\t60\t60\t60\ttp:A:P\tcs:Z::60\n",
    )
    telomere_end_clusters.assigned(SimpleNamespace(out=str(archive), alignments=str(members)), cfg)

    # Simulate Racon consensus output matching representatives
    (archive / "consensus.fasta").write_text((archive / "representatives.fasta").read_text())

    reference = tmp_path / "ref.fasta"
    reference.write_text(f">I\n{flank_a}{'N' * 100}\n")
    genome_paf = tmp_path / "genome.paf.gz"
    _write_gz(
        genome_paf,
        f"E00000001\t{len(flank_a)}\t0\t{len(flank_a)}\t+\tI\t{len(flank_a) + 100}\t0\t{len(flank_a)}\t60\t60\t60\ttp:A:P\n"
        f"E00000002\t{len(flank_b)}\t0\t{len(flank_b)}\t+\tI\t{len(flank_a) + 100}\t0\t{len(flank_b)}\t59\t60\t60\ttp:A:P\n",
    )

    telomere_end_clusters.summarize(
        SimpleNamespace(
            out=str(archive),
            assigned=str(archive / "assigned.paf"),
            genome=str(genome_paf),
            rdna=None,
            reference=str(reference),
            rdna_reference=None,
            reads=str(reads_gz),
            ncrf=str(ncrf_gz),
        ),
        cfg,
    )

    for required_file in gfapdefinisher_handoff.ARCHIVE_FILES:
        assert (archive / required_file).exists(), f"Missing archive file: {required_file}"

    families = list(gfapdefinisher_handoff.rows(archive / "families.tsv"))
    assert len(families) == 1
    assert families[0]["cluster_id"] == "C000001"
    assert int(families[0]["unique_reads"]) == 2
    assert int(families[0]["rdna_any_observations"]) == 0

    haps = list(gfapdefinisher_handoff.rows(archive / "haplotype_patterns.tsv.gz"))
    assert len(haps) == 2


def test_telomere_end_clusters_handles_zero_reads(tmp_path):
    cfg = json.loads((Path(__file__).resolve().parents[1] / "assets" / "telomeric_end_clusters.json").read_text())
    reads_gz = tmp_path / "empty.fasta.gz"
    _write_gz(reads_gz, "")
    ncrf_gz = tmp_path / "empty.ncrf.tsv.gz"
    _write_gz(ncrf_gz, "line\tmotif\tseq\tstart\tend\tstrand\tseqLen\tquerybp\tmRatio\tm\tmm\ti\td\n")

    archive = tmp_path / "empty_archive"
    telomere_end_clusters.extract(
        SimpleNamespace(out=str(archive), reads=str(reads_gz), ncrf=str(ncrf_gz), headers=None, fishing_log=None),
        cfg,
    )
    empty_paf = tmp_path / "empty.paf.gz"
    _write_gz(empty_paf, "")
    telomere_end_clusters.cluster(SimpleNamespace(out=str(archive), overlaps=str(empty_paf)), cfg)
    telomere_end_clusters.assigned(SimpleNamespace(out=str(archive), alignments=str(empty_paf)), cfg)
    (archive / "consensus.fasta").write_text("")
    reference = tmp_path / "ref.fasta"
    reference.write_text(">I\nACGT\n")
    telomere_end_clusters.summarize(
        SimpleNamespace(
            out=str(archive),
            assigned=str(archive / "assigned.paf"),
            genome=str(empty_paf),
            rdna=None,
            reference=str(reference),
            rdna_reference=None,
            reads=str(reads_gz),
            ncrf=str(ncrf_gz),
        ),
        cfg,
    )
    for required_file in gfapdefinisher_handoff.ARCHIVE_FILES:
        assert (archive / required_file).exists(), f"Missing archive file on empty input: {required_file}"
