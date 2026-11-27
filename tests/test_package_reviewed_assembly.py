import gzip
import json
import tempfile
import unittest
from pathlib import Path

from scripts.package_reviewed_assembly import package, sha256


class ReviewedAssemblyPackageTests(unittest.TestCase):
    def test_preserves_fasta_and_gap_with_matching_recipe(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, recipe, out = root / "input.fa.gz", root / "recipe.json", root / "package"
            with gzip.open(fasta, "wt") as handle:
                handle.write(">I\n" + "GCCTAA" * 4 + "ACG" + "N" * 200 + "CCA" + "TTAGGC" * 5 + "\n")
            recipe.write_text(json.dumps({"expected_sha256": sha256(fasta)}))
            package(fasta, out, {"I"}, "TTAGGC", 3, recipe)
            self.assertEqual(sha256(fasta), sha256(out / "assembly.fasta.gz"))
            self.assertIn("I\t28\t227\t2\tN\t200\tcontig\tno\tna", (out / "assembly.agp").read_text())
            self.assertIn("I\t27\t227\t200\tinternal\t2", (out / "unresolved_joins.tsv").read_text())
            self.assertIn("I\t3p\t30\tcanonical_tract", (out / "ends.tsv").read_text())

    def test_wrong_recipe_leaves_no_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, recipe, out = root / "input.fa", root / "recipe.json", root / "package"
            fasta.write_text(">I\nACGT\n")
            recipe.write_text('{"expected_sha256": "invalid"}')
            with self.assertRaisesRegex(ValueError, "checksum"):
                package(fasta, out, {"I"}, "TTAGGC", 3, recipe)
            self.assertFalse(out.exists())

    def test_terminal_repeats_bridges_tip_indels_and_variants(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, out = root / "input.fa", root / "package"
            # 5p: 4x GCCTAA (24 bp) + 1-bp A insertion + 20x GCCTAA (120 bp) -> 144 canonical bp
            # 3p: 15x TTAGGC (90 bp) + 1-bp deletion in TTAGG + 7x TTAGGC (42 bp) -> 132 canonical bp
            seq_5p = "CTAAGC" * 4 + "A" + "GCCTAA" * 20
            seq_3p = "TTAGGC" * 15 + "TTAGG" + "TTAGGC" * 7
            fasta.write_text(f">I\n{seq_5p}{'N' * 100}{seq_3p}\n")
            package(fasta, out, {"I"}, "TTAGGC", 3)
            ends = (out / "ends.tsv").read_text()
            self.assertIn("I\t5p\t144\tcanonical_tract", ends)
            self.assertIn("I\t3p\t132\tcanonical_tract", ends)


if __name__ == "__main__":
    unittest.main()
