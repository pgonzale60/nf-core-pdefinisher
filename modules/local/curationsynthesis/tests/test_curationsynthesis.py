#!/usr/bin/env python3
"""Unit tests for curationsynthesis telomere unitig surveying logic."""

import os
from pathlib import Path
import sys
import tempfile
import unittest

LIB_DIR = Path(__file__).resolve().parent.parent / "resources" / "usr" / "lib"
sys.path.insert(0, str(LIB_DIR))

from curationsynthesis import survey_telomere_unitigs


class TestCurationsynthesis(unittest.TestCase):

    def test_survey_telomere_unitigs_empty(self):
        res = survey_telomere_unitigs(
            telo_paf=None,
            nuclear_chroms=["I", "II"],
            chrom_lens={"I": 1000000, "II": 1000000},
            unitig_assignment={},
        )
        self.assertEqual(res["total_assembled"], 0)
        self.assertEqual(res["assigned_count"], 0)
        self.assertEqual(res["unassigned_count"], 0)
        self.assertEqual(len(res["termini_lacking_unitigs"]), 4)

    def test_survey_telomere_unitigs_paf_and_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paf_file = os.path.join(tmpdir, "telo.paf")
            with open(paf_file, "w") as f:
                # utg000013l: 34948 bp, maps to Chr III 5' (0 to 12284), MAPQ 60, strand +
                # Overhang at 5' = qstart (22664) if strand +? Wait: qlen=34948, qstart=22664, qend=34948
                f.write("utg000013l\t34948\t22664\t34948\t+\tIII\t10000000\t0\t12284\t12200\t12284\t60\n")
                # ptg000004l: 33055 bp, maps to Chr II 5' (0 to 6762), MAPQ 60, strand +
                f.write("ptg000004l\t33055\t26293\t33055\t+\tII\t10000000\t0\t6762\t6700\t6762\t60\n")
                # utg000099l: 15000 bp, maps internally to Chr I (500000 to 515000), MAPQ 60, strand +
                f.write("utg000099l\t15000\t0\t15000\t+\tI\t10000000\t500000\t515000\t15000\t15000\t60\n")

            telo_junctions = {
                "clusters": [
                    {
                        "cluster_id": "cl_12",
                        "read_count": 18,
                        "mean_telomere_tract_bp": 850,
                        "target_chrom": "III",
                        "target_arm": "5p",
                    },
                    {
                        "cluster_id": "cl_17",
                        "read_count": 5,
                        "mean_telomere_tract_bp": 720,
                        "target_chrom": "IV",
                        "target_arm": "5p",
                    },
                ]
            }

            # utg000013l assigned to III 5p; ptg000004l lost arbitration at II 5p
            unitig_assignment = {
                ("III", "5p"): (22664, "utg000013l", 600, 0),
                ("II", "5p"): (4000, "utg000007l", 400, 0),
            }

            res = survey_telomere_unitigs(
                telo_paf=paf_file,
                nuclear_chroms=["I", "II", "III", "IV"],
                chrom_lens={"I": 10000000, "II": 10000000, "III": 10000000, "IV": 10000000},
                unitig_assignment=unitig_assignment,
                telo_junctions=telo_junctions,
            )

            self.assertEqual(res["total_assembled"], 3)
            self.assertEqual(res["assigned_count"], 1)
            self.assertEqual(res["unassigned_count"], 2)

            u13 = next(u for u in res["unitigs"] if u["unitig_id"] == "utg000013l")
            self.assertEqual(u13["status"], "ASSIGNED")
            self.assertEqual(u13["assigned_to"], "Chr III 5p")
            self.assertIn("18 reads (cl_12)", u13["coverage_str"])
            self.assertEqual(u13["hits"][0]["overhang_bp"], 22664)

            p04 = next(u for u in res["unitigs"] if u["unitig_id"] == "ptg000004l")
            self.assertEqual(p04["status"], "COMPETED_OUT")
            self.assertIn("Competed out at Chr II 5p", p04["disposition"])

            u99 = next(u for u in res["unitigs"] if u["unitig_id"] == "utg000099l")
            self.assertEqual(u99["status"], "INTERNAL_ALIGNMENT")
            self.assertIn("Internal alignment on Chr I:500,000–515,000", u99["disposition"])

            # Check that Chr IV 5p is recognized as lacking unitigs and notes cluster cl_17
            self.assertTrue(any("Chr IV 5p (supported by cluster cl_17 with 5 reads)" in s for s in res["termini_lacking_unitigs"]))


if __name__ == "__main__":
    unittest.main()
