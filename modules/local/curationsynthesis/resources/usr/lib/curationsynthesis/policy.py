"""Policy definitions and package fingerprinting for Module 2 Curation Synthesis.

British spelling standard enforced across all documentation and comments.
"""

import hashlib
import json
import os
from pathlib import Path

__version__ = "2.2.0"


def find_curation_policy(policy_path=None, start_path=None):
    """Locate assets/curation_policy.json, the Single Source of Truth."""
    if policy_path and os.path.exists(policy_path):
        return os.path.abspath(policy_path)
    if "CURATION_POLICY" in os.environ and os.path.exists(os.environ["CURATION_POLICY"]):
        return os.path.abspath(os.environ["CURATION_POLICY"])
    current = Path(start_path or __file__).resolve()
    for parent in [current] + list(current.parents):
        candidate = parent / "assets" / "curation_policy.json"
        if candidate.is_file():
            return str(candidate)
        candidate_direct = parent / "curation_policy.json"
        if candidate_direct.is_file():
            return str(candidate_direct)
    return None


def load_policy(policy_path=None, caller_file=None):
    """Load curation_policy.json (Single Source of Truth)."""
    resolved_path = find_curation_policy(policy_path, caller_file)
    if not resolved_path:
        raise FileNotFoundError(
            "Could not locate assets/curation_policy.json (the Single Source of Truth). "
            "Policy parameters must not be hardcoded or duplicated across components."
        )
    with open(resolved_path, "r", encoding="utf-8") as f:
        return json.load(f)


def __getattr__(name):
    if name == "DEFAULT_POLICY":
        return load_policy()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def get_script_fingerprint(entry_file=None):
    """Compute a deterministic SHA-256 fingerprint across the entrypoint and library modules."""
    try:
        hasher = hashlib.sha256()
        pkg_dir = Path(__file__).resolve().parent
        if entry_file and os.path.exists(entry_file):
            with open(entry_file, "rb") as f:
                hasher.update(f.read())
        for py_file in sorted(pkg_dir.glob("*.py")):
            with open(py_file, "rb") as f:
                hasher.update(f.read())
        return hasher.hexdigest()[:16]
    except Exception:
        return "unknown"
