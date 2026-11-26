import csv
import gzip
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


RESOLVER = Path(__file__).resolve().parents[1] / "resources/usr/bin/pde_rdna_cassette_resolver.py"


class RdnaEvidenceTests(unittest.TestCase):
    def test_reads_compressed_primary_graph(self):
        parse_gfa = runpy.run_path(str(RESOLVER))["parse_hifiasm_gfa_depths"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.p_utg.gfa.gz"
            with gzip.open(path, "wt") as graph:
                graph.write("S\tutg1\tACGT\tLN:i:4\trd:i:23\n")
            self.assertEqual(parse_gfa(path)["utg1"]["depth"], 23)

    def test_missing_rDNA_evidence_retains_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            blueprint, reference, manifest = (
                directory / "blueprint.json", directory / "reference.fa", directory / "manifest.tsv"
            )
            blueprint.write_text(json.dumps({
                "action_cards": [{"action": "TERMINUS_RDNA_ARRAY", "card_id": "card1", "chrom": "II"}]
            }))
            reference.write_text(">II\nACGT\n")
            subprocess.run([
                sys.executable, str(RESOLVER), "--blueprint", str(blueprint), "--ref", str(reference),
                "--sample", "test", "--outdir", str(directory / "cassettes"), "--manifest", str(manifest),
            ], check=True, capture_output=True)
            with manifest.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["status"], "DEFERRED_NO_ASSEMBLY_EVIDENCE")
            self.assertEqual(rows[0]["cassette_fasta"], "NONE")


if __name__ == "__main__":
    unittest.main()
