"""Tests for pde_segment_quantized (v2.0.0): R1-R11 segmentation via the CLI.

Synthetic fixtures use 500 bp mosdepth windows. Expected values follow the
R1-R11 plan (Verification section).
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

bin_dir = Path(__file__).resolve().parent.parent / "resources" / "usr" / "bin"
SCRIPT = bin_dir / "pde_segment_quantized.py"
WIN = 500


def write(path, lines):
    path.write_text("".join(line + "\n" for line in lines))


def depth_windows(chrom, length, depth_fn):
    """Mosdepth-format regions: chrom, start, end, depth."""
    return [f"{chrom}\t{s}\t{s + WIN}\t{depth_fn(s):g}" for s in range(0, length, WIN)]


def run_case(length, depth_fn, sites=(), gaps=(), chromosomes="chr1", barrnap=(), summary_chroms=None):
    """Run the CLI on one synthetic sequence set. Returns (outdir, bed_rows, tsv_lines, fission_lines)."""
    tmp = Path(tempfile.mkdtemp())
    chroms = summary_chroms or {"chr1": length}
    summary = ["# GRS_THRESHOLD=0.5", "# ELIMINATION_RANGE\tcanon_uncertain_bp=0"]
    summary.append("chromosome\tstart\tend\tlength\tclassification\tmedian_depth\tsubtype")
    for c, n in chroms.items():
        summary.append(f"{c}\t0\t{n}\t{n}\tcore\t0\tx")
    write(tmp / "summary.tsv", summary)
    rows = ["chrom\tpos\tend\tstrand\tsupport"] + list(sites)
    write(tmp / "breaks.tsv", rows)
    write(tmp / "gaps.bed", list(gaps))
    write(tmp / "red.bed", [])
    regions = []
    for c, n in chroms.items():
        regions += depth_windows(c, n, lambda s, c=c: depth_fn(c, s))
    write(tmp / "regions.bed", regions)
    write(tmp / "barrnap.gff", list(barrnap))
    cmd = [
        sys.executable, str(SCRIPT),
        "--summary-tsv", str(tmp / "summary.tsv"),
        "--barrnap-gff", str(tmp / "barrnap.gff"),
        "--gaps-bed", str(tmp / "gaps.bed"),
        "--red-bed", str(tmp / "red.bed"),
        "--regions-bed", str(tmp / "regions.bed"),
        "--breakpoints-tsv", str(tmp / "breaks.tsv"),
        "--chromosomes", chromosomes,
        "--prefix", "t", "--outdir", str(tmp),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    bed = [l.split("\t") for l in (tmp / "t.pde_segments.bed").read_text().splitlines() if l]
    tsv = (tmp / "t.final.grs_summary.tsv").read_text().splitlines()
    fis = (tmp / "t.fissions.tsv").read_text().splitlines()
    return tmp, bed, tsv, fis


def labels(bed):
    return [(int(r[1]), int(r[2]), r[3]) for r in bed]


class TestFission(unittest.TestCase):
    """Case 2: opposing +/- sites that are close form a fission, not a GRS."""

    def test_plus_minus_pair_cuts_one_base_after_plus(self):
        _, bed, _, fis = run_case(
            5000, lambda c, s: 100,
            sites=["chr1\t1000\t1000\t+\t100", "chr1\t1003\t1003\t-\t80"],
        )
        self.assertIn("chr1\t1000\t1003\t3\t1001", fis)
        segs = labels(bed)
        self.assertEqual([x[2] for x in segs], ["putative_core", "putative_core"])
        self.assertEqual(segs[0][1], 1001)

    def test_minus_before_plus_cuts_one_base_before_plus(self):
        _, bed, _, fis = run_case(
            5000, lambda c, s: 100,
            sites=["chr1\t1000\t1000\t+\t100", "chr1\t997\t997\t-\t80"],
        )
        # Cut = plus - 1 when the minus site is upstream of the plus site.
        self.assertIn("chr1\t1000\t997\t-3\t999", fis)
        self.assertNotIn("putative_GRS", [x[2] for x in labels(bed)])


class TestSpan(unittest.TestCase):
    """Cases 3 and 4: a +/- span with depth inside it is eliminated, gaps change the call."""

    def span_sites(self):
        return ["chr1\t10000\t10000\t+\t100", "chr1\t30000\t30000\t-\t100"]

    def test_depth_inside_span_is_grs(self):
        _, bed, tsv, _ = run_case(
            40000, lambda c, s: 20 if 10000 <= s < 30000 else 100, sites=self.span_sites(),
        )
        segs = [x for x in labels(bed) if x[2] == "putative_GRS"]
        self.assertEqual(segs, [(10000, 30000, "putative_GRS")])
        self.assertTrue(any("span_eliminated_depth:20000" in l for l in tsv))

    def test_gap_at_span_start_with_retained_depth_is_uncertain(self):
        _, bed, tsv, _ = run_case(
            40000, lambda c, s: 100, sites=self.span_sites(), gaps=["chr1\t10000\t10500"],
        )
        self.assertIn("uncertain", [x[2] for x in labels(bed)])
        self.assertTrue(any("retained_depth_inside_span" in l for l in tsv))

    def test_gaps_at_both_span_edges_with_retained_depth_is_misjoin(self):
        _, bed, _, _ = run_case(
            40000, lambda c, s: 100, sites=self.span_sites(),
            gaps=["chr1\t10000\t10500", "chr1\t29500\t30000"],
        )
        self.assertIn("putative_core", [x[2] for x in labels(bed)])


class TestVariable(unittest.TestCase):
    """Case 6: half-depth span is putative_variable, not GRS."""

    def test_half_depth_span_is_variable(self):
        sites = ["chr1\t10000\t10000\t+\t100", "chr1\t20000\t20000\t+\t100", "chr1\t30000\t30000\t-\t100"]
        _, bed, _, _ = run_case(
            40000, lambda c, s: 50 if 10000 <= s < 20000 else 100, sites=sites,
        )
        self.assertIn("putative_variable", [x[2] for x in labels(bed)])


class TestUnplaced(unittest.TestCase):
    """Case 9: a non-canonical sequence with low depth and no sites is unplaced_eliminated."""

    def test_low_depth_non_canonical_is_unplaced_eliminated(self):
        _, bed, _, _ = run_case(
            40000, lambda c, s: 20 if c == "chr2" else 100, chromosomes="chr1",
            summary_chroms={"chr1": 40000, "chr2": 40000},
        )
        chr2 = [r[3] for r in bed if r[0] == "chr2"]
        self.assertEqual(chr2, ["unplaced_eliminated"])


class TestOutput(unittest.TestCase):
    """Case 10: tiling, header and output contract."""

    def test_header_and_tiling(self):
        _, bed, tsv, _ = run_case(40000, lambda c, s: 100)
        header = [l for l in tsv if l.startswith("chromosome\t")]
        self.assertEqual(len(header), 1)
        self.assertEqual(header[0].split("\t")[:6], ["chromosome", "start", "end", "length", "classification", "median_depth"])
        self.assertEqual(sum(int(r[2]) - int(r[1]) for r in bed), 40000)


if __name__ == "__main__":
    unittest.main()
