#!/usr/bin/env python3
"""Container matrix and URL verification test.

Ensures that every container URL across the v2-tool-first architectural upgrades
is non-fabricated, matches the ARCHITECT_SPEC matrix, and returns HTTP 200 when reachable.
"""

import json
import re
import unittest
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Verified image tags specified in ARCHITECT_SPEC.md for v2-tool-first
UPGRADE_CONTAINER_MATRIX = {
    "modules/local/extracttelomericreads/main.nf": "samtools_bbmap_pigz:c8f695696758ff77",
    "modules/local/telomericjunctionprofiler/main.nf": "samtools_minimap2_vsearch_python:1aa370c14225a820",
    "modules/local/teloclip/extract/main.nf": "teloclip_samtools:2e39bcbf9e597d36",
    "modules/local/targetedgapassembly/main.nf": "hifiasm_abpoa_minimap2_samtools_mindthegap_python:4e99aa7a05e907bf",
    "modules/local/telomerewalker/main.nf": "hifiasm_abpoa_minimap2_samtools_gfatools_python:b88b84d24c308348",
    "modules/local/ragtag/agpcheck/main.nf": "ragtag:2.1.0--6125126fe9d553ab",
}


def get_wave_token(image_name: str) -> str:
    url = f"https://cerbero.seqera.io/auth/token?service=community.wave.seqera.io&scope=repository:library/{image_name}:pull"
    req = urllib.request.Request(url, headers={"User-Agent": "curl/7.68.0"})
    with urllib.request.urlopen(req, timeout=10) as resp:
        data = json.loads(resp.read().decode())
        return data["token"]


def check_wave_manifest(image_name: str, tag: str) -> int:
    token = get_wave_token(image_name)
    url = f"https://community.wave.seqera.io/v2/library/{image_name}/manifests/{tag}"
    req = urllib.request.Request(
        url,
        headers={
            "Authorization": f"Bearer {token}",
            "Accept": "application/vnd.oci.image.manifest.v1+json, application/vnd.docker.distribution.manifest.v2+json",
            "User-Agent": "curl/7.68.0",
        },
    )
    with urllib.request.urlopen(req, timeout=10) as resp:
        return resp.status


class TestContainerMatrix(unittest.TestCase):

    def test_upgrade_modules_match_architect_spec(self):
        """Verifies that all modified/new modules in v2-tool-first match ARCHITECT_SPEC tags."""
        for rel_path, expected_img in UPGRADE_CONTAINER_MATRIX.items():
            nf_path = REPO_ROOT / rel_path
            self.assertTrue(nf_path.exists(), f"Module file {rel_path} does not exist!")
            content = nf_path.read_text()
            matches = re.findall(r"community\.wave\.seqera\.io/library/([a-zA-Z0-9_\-\.]+:[a-zA-Z0-9_\-\.]+)", content)
            self.assertTrue(len(matches) > 0, f"No Wave container found in {rel_path}")
            self.assertEqual(
                matches[0],
                expected_img,
                f"Container in {rel_path} ({matches[0]}) does not match expected {expected_img}",
            )

    def test_live_wave_container_manifests(self):
        """Verifies that all 5 upgrade containers return HTTP 200 from Seqera Wave registry."""
        upgrades = [
            ("samtools_bbmap_pigz", "c8f695696758ff77"),
            ("samtools_minimap2_vsearch_python", "1aa370c14225a820"),
            ("teloclip_samtools", "2e39bcbf9e597d36"),
            ("hifiasm_abpoa_minimap2_samtools_mindthegap_python", "4e99aa7a05e907bf"),
            ("hifiasm_abpoa_minimap2_samtools_gfatools_python", "b88b84d24c308348"),
        ]
        for name, tag in upgrades:
            status = check_wave_manifest(name, tag)
            self.assertEqual(status, 200, f"Wave container {name}:{tag} failed with status {status}")


if __name__ == "__main__":
    unittest.main()
