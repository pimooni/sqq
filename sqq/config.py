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


# Defaults are intentionally explicit so a run can be reproduced from
# run_config.yaml without relying on hidden command-line assumptions.
DEFAULT_CONFIG: dict[str, Any] = {
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
        "primitive": True,
        "chordless": True,
    },
    "cup": {
        "mode": "general",
        "enabled": True,
        "base_sizes": "auto",
        "side_sizes": "auto",
        "max_combinations_per_base": 50000,
    },
    "cage": {
        "ring_sizes": [5, 6],
        "target_types": ["512", "51262", "51263", "51264"],
        "output_other": False,
        "other_max_faces": 20,
        "enabled": True,
        "search_mode": "grow",
        "seed_mode": "ring",
        "max_states_per_seed": 2000,
        "max_total_states": 250000,
        "occupancy_mode": "polyhedron",
        "occupancy_radius_nm": 0.5,
    },
    "order": {
        "f3f4_enabled": True,
        "focus_waters": [],
    },
    "ice": {
        "enabled": True,
        "method": "chill",
        "min_six_rings": 2,
        "require_four_coord_neighbors": True,
    },
    "output": {
        "write_gro": True,
        "write_tsv": False,
        "write_vmd": False,
        "write_xlsx_summary": True,
        "write_empty_files": False,
        "gro_atom_mode": "full_water",
        "center_resname": "CNT",
    },
    "parallel": {
        "n_jobs": "auto",
    },
    "debug": {
        "use_networkx_checks": False,
    },
}


def load_config(path: Path | None) -> dict[str, Any]:
    """Load a config file and merge it over built-in defaults."""
    config = deepcopy(DEFAULT_CONFIG)
    if path is None:
        return config
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
    return merge_config(config, user_config)


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







