"""Auto-synced Single Source of Truth (SSOT) parameters from conf/cbm_params.json."""

import json
from pathlib import Path

_FALLBACK_SSOT = {'aux_evalue_cutoff': 0.05,
 'cbm_catalytic_proximity_decay_bp': 75.0,
 'cleavage_mh_protection_radius': 6,
 'core_evalue_cutoff': 0.05,
 'core_grid_models': ['zoops', 'anr', 'oops'],
 'core_grid_radii': [30, 50, 100, 150],
 'core_grid_width_ranges': [{'args': ['-minw', '6', '-maxw', '15'], 'label': 'w6_15'},
                            {'args': ['-minw', '12', '-maxw', '30'], 'label': 'w12_30'},
                            {'args': ['-minw', '15', '-maxw', '50'],
                             'label': 'w15_50'}],
 'distal_bg_exclusion_radius': 60,
 'fimo_aux_threshold': 0.001,
 'fimo_core_threshold': 0.0001,
 'itr_min_consecutive_repeats': 4,
 'kmer_size': 31,
 'max_anr_site_density_credit': 2.0,
 'max_germ_flank_len': 300,
 'max_soma_flank_len': 300,
 'min_cluster_size_hifi': 1,
 'min_cluster_size_illumina': 3,
 'min_coherence': 0.85,
 'min_germ_reads': 3,
 'min_telo_reads': 5,
 'satellite_max_period': 100,
 'satellite_min_copies': 3,
 'satellite_min_identity': 0.82,
 'satellite_min_period': 10,
 'somatic_dedup_cbm_exclusion_bp': 20,
 'stage1_min_identity': 0.98,
 'stage3_min_identity': 0.85,
 'telomere_repeat': 'TTAGGC',
 'telomeric_echo_1mm_cov_cutoff': 0.7,
 'telomeric_echo_exact_kmer_cutoff': 0.5,
 'terminal_bait_len': 60}


def load_ssot_params() -> dict:
    """Load canonical parameters from conf/cbm_params.json if accessible, else synced SSOT."""
    cur = Path(__file__).resolve().parent
    for _ in range(8):
        cand = cur / "conf" / "cbm_params.json"
        if cand.is_file():
            return json.loads(cand.read_text(encoding="utf-8"))
        if cur.parent == cur:
            break
        cur = cur.parent
    return dict(_FALLBACK_SSOT)


SSOT = load_ssot_params()
