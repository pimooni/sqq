from __future__ import annotations

"""Configuration defaults and YAML/JSON loading."""

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    yaml = None


DEFAULT_MODE = "50"
MODE_PRESETS: dict[str, dict[str, Any]] = {
    "00": {
        "label": "rigorous",
        "worker_fraction": 0.25,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
    },
    "50": {
        "label": "standard",
        "worker_fraction": 0.50,
        "bond_mode": "auto",
        "ring_sizes": [5, 6],
    },
    "99": {
        "label": "performance",
        "worker_fraction": 0.90,
        "bond_mode": "oo",
        "ring_sizes": [5, 6],
    },
}


# Defaults are intentionally explicit so a run can be reproduced from
# run_config.yaml without relying on hidden command-line assumptions.
DEFAULT_CONFIG: dict[str, Any] = {
    "mode": DEFAULT_MODE,
    "input": {
        "pattern": "*.gro",
        "recursive": False,
        "first_file_time_ps": 0.0,
        "frame_time_step_ps": 100.0,
        "xtc_stride": 1,
    },
    "water": {
        "resnames": ["SOL", "TIP", "WAT", "HOH"],
        "oxygen_names": ["OW", "O", "OH2"],
        "hydrogen_names": ["HW1", "HW2", "H1", "H2", "HW", "HT1", "HT2"],
    },
    "guest": {
        "resnames": ["CH4", "CO2", "MET", "ETH"],
        "center_atoms": {"CH4": ["C"], "CO2": ["C"]},
        "center_mode": "center_atom",
    },
    "graph": {
        "bond_mode": "auto",
        "oo_cutoff_nm": 0.35,
        "hbond_distance_nm": 0.35,
        "hbond_angle_deg": 30.0,
        "pair_file": None,
        "pair_id": "resid",
    },
    "pbc": {
        "box_mode": "orthorhombic",
    },
    "ring": {
        "sizes": [5, 6],
        "report_sizes": "auto",
        "chordless": True,
        "definition": "chordless",
    },
    "quasi_cage": {
        "search_policy": "bounded",
        "enabled": True,
        "base_sizes": "auto",
        "side_sizes": "auto",
        "max_combinations_per_base": 50000,
        "max_layers": 1,
        "max_rings_per_layer": 6,
        "max_layer_states_per_seed": 200,
        "max_candidates_per_edge": 4,
        "max_layer_candidates": 24,
    },
    "cage": {
        "report_types": "auto",
        "max_faces": 20,
        "enabled": True,
        "search_mode": "grow",
        "seed_mode": "ring",
        "max_states_per_seed": 20000,
        "max_total_states": 5000000,
        "max_boundary_candidates": 8,
        "fast_closure": True,
        "fast_closure_max_states": 20000,
        "scientific_validation": False,
        "max_face_planarity_rms_nm": 0.06,
        "max_face_edge_cv": 0.35,
        "min_cage_volume_nm3": 1.0e-6,
        "occupancy_mode": "polyhedron",
        "occupancy_radius_nm": 0.5,
    },
    "hydrate_cluster": {
        "enabled": False,
        "min_cage": 2,
        "detail": False,
    },
    "hydrate_order": {
        "mcg1_enabled": True,
        "mcg3_enabled": False,
        "mcg_guest_resnames": ["CH4", "MET"],
        "mcg_guest_cutoff_nm": 0.90,
        "mcg_water_cutoff_nm": 0.60,
        "mcg_cone_half_angle_deg": 45.0,
        "mcg_min_waters": 5,
        "dhop35_enabled": True,
        "dhop30_enabled": False,
        "dhop_neighbor_cutoff_nm": 0.35,
        "dhop_planar_counts": [11, 12],
        "dhop_min_qualified_neighbors": 3,
    },
    "order": {
        "f3f4_enabled": True,
        "q_enabled": True,
        "q_degree": [6, 12],
        "q_neighbor_mode": "graph",
        "q_cutoff_nm": 0.35,
        "q_n_neighbor": None,
        "focus_waters": [],
    },
    "ice": {
        "enabled": True,
        "method": "chill",
        "min_six_rings": 2,
        "require_four_coord_neighbors": True,
    },
    "output": {
        "write_info": True,
        "write_gro": True,
        "write_ring_gro": True,
        "write_half_cage_gro": True,
        "write_quasi_cage_gro": True,
        "write_cage_gro": True,
        "write_ice_gro": True,
        "write_tsv": False,
        "write_order_tsv": False,
        "write_vmd": False,
        "write_xlsx_summary": True,
        "write_empty_files": False,
        "structure_layout": "grouped",
        "gro_atom_mode": "full_water",
        "center_resname": "CNT",
    },
    "parallel": {
        "backend": "process",
        "workers": "auto",
        "math_threads": 1,
    },
    "debug": {
        "use_networkx_checks": False,
    },
}


def normalize_mode(value: Any) -> str:
    """Normalize and validate a two-digit analysis mode."""
    text = str(value).strip()
    if text.isdigit():
        text = text.zfill(2)
    if text not in MODE_PRESETS:
        choices = ", ".join(MODE_PRESETS)
        raise ValueError(f"mode must be one of: {choices}")
    return text


def mode_label(mode: Any) -> str:
    """Return the human-readable mode label."""
    return str(MODE_PRESETS[normalize_mode(mode)]["label"])


def mode_worker_fraction(mode: Any) -> float:
    """Return the automatic worker fraction for a mode."""
    return float(MODE_PRESETS[normalize_mode(mode)]["worker_fraction"])


def apply_mode_preset(config: dict[str, Any], mode: Any) -> dict[str, Any]:
    """Apply the scientific and worker-policy base settings for one mode."""
    normalized = normalize_mode(mode)
    preset = MODE_PRESETS[normalized]
    config["mode"] = normalized
    config["graph"]["bond_mode"] = preset["bond_mode"]
    config["ring"]["sizes"] = list(preset["ring_sizes"])
    config["ring"]["report_sizes"] = "auto"
    config["quasi_cage"]["base_sizes"] = "auto"
    config["quasi_cage"]["side_sizes"] = "auto"
    return config


def load_config(path: Path | None, mode: Any = None) -> dict[str, Any]:
    """Load a mode preset, then merge user configuration over it."""
    user_config: dict[str, Any] = {}
    if path is not None:
        with path.open("r", encoding="utf-8-sig") as handle:
            text = handle.read()
        if yaml is not None:
            user_config = yaml.safe_load(text) or {}
        else:
            # Source-tree smoke tests can run without PyYAML; installed SQQ uses YAML.
            try:
                user_config = json.loads(text) if text.strip() else {}
            except json.JSONDecodeError as exc:
                raise RuntimeError("Reading YAML config files requires PyYAML. Install with `pip install -e .`.") from exc
        if not isinstance(user_config, dict):
            raise ValueError(f"Config file must contain a YAML mapping: {path}")

    selected_mode = normalize_mode(mode if mode is not None else user_config.get("mode", DEFAULT_MODE))
    config = apply_mode_preset(deepcopy(DEFAULT_CONFIG), selected_mode)
    merge_config(config, user_config)
    config["mode"] = selected_mode
    return config


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge user configuration into defaults."""
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge_config(base[key], value)
        else:
            base[key] = value
    return base


def write_default_config(path: Path) -> None:
    """Write the default configuration template."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        dump_config(DEFAULT_CONFIG, handle)


def dump_config(config: dict[str, Any], handle) -> None:
    """Write YAML when available, otherwise a JSON-compatible fallback."""
    if yaml is not None:
        yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
    else:
        json.dump(config, handle, ensure_ascii=False, indent=2)
        handle.write("\n")






