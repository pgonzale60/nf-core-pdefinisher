import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from scripts.audit_chromosome_ends import (
    marker_panel,
    query_anchor,
    read_breaksites,
    repeat_runs,
    reverse_complement,
)


class ChromosomeEndAuditTests(unittest.TestCase):
    def test_breaks_at_both_physical_boundaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "breakpoints.bed"
            path.write_text("I\t-1\t0\tType:S2G\t2\t-\nI\t100\t101\tType:S2G\t2\t+\n")
            self.assertEqual(read_breaksites(path, {"I": 100}), {("I", "-"): [0], ("I", "+"): [100]})

    def test_footprint_ignores_shared_markers_and_overlapping_bases(self):
        monomer = "GATTACA"
        sequence = reverse_complement(monomer * 20)
        self.assertFalse(marker_panel({"A": monomer, "B": monomer}, 5))
        markers = marker_panel({"A": monomer}, 5)
        self.assertEqual(list(repeat_runs(sequence, markers, {"A": 7}, 5)), [("A", 0, 140, 140)])
        self.assertEqual(list(repeat_runs(sequence, markers, {"A": 7}, 5, offset=60)), [])

    def test_anchor_uses_inward_orientation_and_excludes_softclip(self):
        read = SimpleNamespace(reference_start=100, query_sequence="A" * 120,
                               cigartuples=[(4, 10), (0, 100), (4, 10)])
        self.assertEqual(query_anchor(read, [(120, 180)], False, 20), 50)
        self.assertEqual(query_anchor(read, [(120, 180)], True, 20), 50)
        self.assertIsNone(query_anchor(read, [(120, 180)], False, 61))


if __name__ == "__main__":
    unittest.main()
