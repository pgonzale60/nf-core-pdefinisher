from __future__ import annotations

import hashlib
from typing import Dict, Any

__version__ = "1.0.0"

DEFAULT_VALIDATION_POLICY: Dict[str, Any] = {
    "version": __version__,
    "min_edit_flank_bp": 500,
    "min_edit_depth": 5.0,
    "min_spanning_reads": 2,
    "flagrant_hom_vaf_threshold": 0.80,
    "flagrant_min_support": 5,
    "flagrant_min_sv_len": 50,
    "flagrant_dropout_min_window": 10000,
    "flagrant_dropout_max_cov_ratio": 0.10,
}


def get_script_fingerprint() -> str:
    """Computes SHA256 digest of package components."""
    try:
        from pathlib import Path
        pkg_dir = Path(__file__).resolve().parent
        hasher = hashlib.sha256()
        for p in sorted(pkg_dir.glob("*.py")):
            hasher.update(p.read_bytes())
        return hasher.hexdigest()[:16]
    except Exception:
        return "unknown"
