from __future__ import annotations

"""Configuration defaults and YAML/JSON loading."""

from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    yaml = None


DEFAULT_MODE = "50"
DEFAULT_ORDER_PARAMETERS = ("f3", "f4")
ALL_ORDER_PARAMETERS = (
    "f3",
    "f4",
    "q6",
    "q12",
    "mcg1",
    "mcg3",
    "dhop35",
    "dhop30",
)
ORDER_PARAMETER_ALIASES = {
    "mcg-1": "mcg1",
    "mcg_1": "mcg1",
    "mcg-3": "mcg3",
    "mcg_3": "mcg3",
    "dhop-35": "dhop35",
    "dhop_35": "dhop35",
    "dhop-30": "dhop30",
    "dhop_30": "dhop30",
}
DEFAULT_DISABLED_OUTPUTS: tuple[str, ...] = ()
ALL_DISABLED_OUTPUTS = (
    "info",
    "membership-tsv",
    "order-tsv",
    "vmd",
    "gro",
    "xlsx",
    "summary-detail",
)
OUTPUT_TYPE_ORDER = (
    "info",
    "membership-tsv",
    "order-tsv",
    "vmd",
    "gro",
    "ring-gro",
    "half-gro",
    "quasi-gro",
    "cage-gro",
    "ice-gro",
    "xlsx",
    "summary-detail",
)
OUTPUT_TYPE_ALIASES = {
    "info-md": "info",
    "half-cage-gro": "half-gro",
    "quasi-cage-gro": "quasi-gro",
    "detail-csv": "summary-detail",
}
GRO_OUTPUT_TYPES = {
    "ring-gro",
    "half-gro",
    "quasi-gro",
    "cage-gro",
    "ice-gro",
}
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


# Explicit defaults keep run_config.yaml reproducible.
DEFAULT_CONFIG: dict[str, Any] = {
    "mode": DEFAULT_MODE,
    "input": {
        "pattern": "*.gro",
        "recursive": False,
        "first_file_time_ps": 0.0,
        "frame_time_step_ps": 100.0,
        "xtc_stride": 1,
        "xyz_scale": 0.1,
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
        "mcg_guest_resnames": ["CH4", "MET"],
        "mcg_guest_cutoff_nm": 0.90,
        "mcg_water_cutoff_nm": 0.60,
        "mcg_cone_half_angle_deg": 45.0,
        "mcg_min_waters": 5,
        "dhop_neighbor_cutoff_nm": 0.35,
        "dhop_planar_counts": [11, 12],
        "dhop_min_qualified_neighbors": 3,
    },
    "order": {
        "parameters": ["f3", "f4"],
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
        "disabled_outputs": [],
        "write_tsv": False,
        "write_order_tsv": False,
        "write_vmd": False,
        "summary_detail_dir": "summary_detail",
        "cage_isomer_rows": "nonzero",
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


def normalize_order_parameters(value: Any = None) -> tuple[str, ...]:
    """Normalize the unified order-parameter selection into stable output order."""
    if value is None or value == "":
        raw_items: list[Any] = list(DEFAULT_ORDER_PARAMETERS)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "order.parameters / --order-parameter must be a comma-separated list."
            ) from exc
    if not raw_items:
        return ()

    cleaned = [
        ORDER_PARAMETER_ALIASES.get(str(item).strip().lower(), str(item).strip().lower())
        for item in raw_items
    ]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in order.parameters / --order-parameter.")
        return ALL_ORDER_PARAMETERS if cleaned[0] == "all" else ()

    supported_fixed = {"f3", "f4", "mcg1", "mcg3", "dhop35", "dhop30"}
    normalized: set[str] = set()
    for name in cleaned:
        if name in supported_fixed:
            normalized.add(name)
            continue
        match = re.fullmatch(r"q(\d+)", name)
        if match:
            normalized.add(f"q{int(match.group(1))}")
            continue
        raise ValueError(
            f"Unsupported order parameter '{name}'. Use f3, f4, qN, mcg1, mcg3, "
            "dhop35, dhop30, all, or none."
        )
    return tuple(sorted(normalized, key=order_parameter_sort_key))


def order_parameter_sort_key(name: str) -> tuple[int, int]:
    """Return the canonical display order for one normalized parameter name."""
    fixed_order = {
        "f3": (0, 0),
        "f4": (1, 0),
        "mcg1": (3, 0),
        "mcg3": (4, 0),
        "dhop35": (5, 0),
        "dhop30": (6, 0),
    }
    if name in fixed_order:
        return fixed_order[name]
    match = re.fullmatch(r"q(\d+)", name)
    if match:
        return 2, int(match.group(1))
    return 99, 0


def q_degrees_from_order_parameters(value: Any) -> tuple[int, ...]:
    """Extract the selected Q_l degrees from a unified parameter selection."""
    return tuple(
        int(name[1:])
        for name in normalize_order_parameters(value)
        if name.startswith("q")
    )


def order_parameter_display(value: Any) -> str:
    """Render a normalized selection for terminal and report metadata."""
    parameters = normalize_order_parameters(value)
    return ", ".join(parameters) if parameters else "none"


def normalize_disabled_outputs(value: Any = None) -> tuple[str, ...]:
    """Normalize the unified output-suppression list."""
    if value is None or value == "":
        raw_items: list[Any] = list(DEFAULT_DISABLED_OUTPUTS)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "output.disabled_outputs / --no-output must be a comma-separated list."
            ) from exc
    if not raw_items:
        return ()

    cleaned = [
        OUTPUT_TYPE_ALIASES.get(str(item).strip().lower(), str(item).strip().lower())
        for item in raw_items
    ]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in output.disabled_outputs / --no-output.")
        return ALL_DISABLED_OUTPUTS if cleaned[0] == "all" else ()

    supported = set(OUTPUT_TYPE_ORDER)
    unknown = sorted(set(cleaned) - supported)
    if unknown:
        raise ValueError(
            f"Unsupported output type(s) {unknown}. Use info, membership-tsv, "
            "order-tsv, vmd, gro, ring-gro, half-gro, quasi-gro, cage-gro, "
            "ice-gro, xlsx, summary-detail, all, or none."
        )
    normalized = set(cleaned)
    if "gro" in normalized:
        normalized.difference_update(GRO_OUTPUT_TYPES)
    return tuple(name for name in OUTPUT_TYPE_ORDER if name in normalized)


def disabled_output_display(value: Any) -> str:
    """Render normalized disabled outputs for terminal and metadata."""
    outputs = normalize_disabled_outputs(value)
    return ", ".join(outputs) if outputs else "none"


def output_enabled(config: dict[str, Any], output_type: str) -> bool:
    """Return whether one output category is enabled by the unified policy."""
    normalized_type = OUTPUT_TYPE_ALIASES.get(
        str(output_type).strip().lower(),
        str(output_type).strip().lower(),
    )
    if normalized_type not in set(OUTPUT_TYPE_ORDER):
        raise ValueError(f"Unsupported output type: {output_type}")
    disabled = set(
        normalize_disabled_outputs(
            config.get("output", {}).get("disabled_outputs", DEFAULT_DISABLED_OUTPUTS)
        )
    )
    if normalized_type in GRO_OUTPUT_TYPES and "gro" in disabled:
        return False
    return normalized_type not in disabled


def migrate_legacy_disabled_outputs(
    user_config: dict[str, Any],
) -> tuple[str, ...] | None:
    """Translate pre-0.2.7 output booleans when the unified list is absent."""
    output = user_config.get("output", {})
    if not isinstance(output, dict) or "disabled_outputs" in output:
        return None
    legacy_keys = {
        "write_info",
        "write_gro",
        "write_ring_gro",
        "write_half_cage_gro",
        "write_quasi_cage_gro",
        "write_cage_gro",
        "write_ice_gro",
        "write_xlsx_summary",
        "write_summary_detail_csv",
    }
    if not legacy_keys & set(output):
        return None
    disabled: list[str] = []
    for key, name in (
        ("write_info", "info"),
        ("write_gro", "gro"),
        ("write_ring_gro", "ring-gro"),
        ("write_half_cage_gro", "half-gro"),
        ("write_quasi_cage_gro", "quasi-gro"),
        ("write_cage_gro", "cage-gro"),
        ("write_ice_gro", "ice-gro"),
        ("write_xlsx_summary", "xlsx"),
        ("write_summary_detail_csv", "summary-detail"),
    ):
        if key in output and not legacy_enabled(output[key]):
            disabled.append(name)
    return normalize_disabled_outputs(disabled)


def strip_legacy_selection_keys(user_config: dict[str, Any]) -> None:
    """Remove migrated selector booleans so run_config has one source of truth."""
    order = user_config.get("order", {})
    if isinstance(order, dict):
        for key in ("f3f4_enabled", "q_enabled", "q_degree"):
            order.pop(key, None)
    hydrate_order = user_config.get("hydrate_order", {})
    if isinstance(hydrate_order, dict):
        for key in (
            "mcg1_enabled",
            "mcg3_enabled",
            "dhop35_enabled",
            "dhop30_enabled",
        ):
            hydrate_order.pop(key, None)
    output = user_config.get("output", {})
    if isinstance(output, dict):
        for key in (
            "write_info",
            "write_gro",
            "write_ring_gro",
            "write_half_cage_gro",
            "write_quasi_cage_gro",
            "write_cage_gro",
            "write_ice_gro",
            "write_xlsx_summary",
            "write_summary_detail_csv",
        ):
            output.pop(key, None)


def migrate_legacy_order_parameters(user_config: dict[str, Any]) -> tuple[str, ...] | None:
    """Translate explicit pre-0.2.7 enable flags when no unified list is present."""
    order = user_config.get("order", {})
    hydrate_order = user_config.get("hydrate_order", {})
    if not isinstance(order, dict) or not isinstance(hydrate_order, dict):
        return None
    if "parameters" in order:
        return None
    legacy_order_keys = {"f3f4_enabled", "q_enabled", "q_degree"}
    legacy_hydrate_keys = {
        "mcg1_enabled",
        "mcg3_enabled",
        "dhop35_enabled",
        "dhop30_enabled",
    }
    if not (legacy_order_keys & set(order) or legacy_hydrate_keys & set(hydrate_order)):
        return None

    selected: list[str] = []
    if legacy_enabled(order.get("f3f4_enabled", True)):
        selected.extend(("f3", "f4"))
    if legacy_enabled(order.get("q_enabled", True)):
        degrees = order.get("q_degree", [6, 12])
        if isinstance(degrees, str):
            degrees = [item.strip() for item in degrees.split(",") if item.strip()]
        selected.extend(f"q{int(degree)}" for degree in degrees)
    for key, name, default in (
        ("mcg1_enabled", "mcg1", True),
        ("mcg3_enabled", "mcg3", False),
        ("dhop35_enabled", "dhop35", True),
        ("dhop30_enabled", "dhop30", False),
    ):
        if legacy_enabled(hydrate_order.get(key, default)):
            selected.append(name)
    return normalize_order_parameters(selected or ["none"])


def legacy_enabled(value: Any) -> bool:
    """Interpret old YAML booleans without treating the string 'false' as true."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"false", "off", "no", "0", "none", ""}:
        return False
    if text in {"true", "on", "yes", "1"}:
        return True
    raise ValueError(
        f"Legacy enable/disable value must be boolean-like; got {value!r}."
    )


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
            # JSON fallback supports source-tree smoke tests without PyYAML.
            try:
                user_config = json.loads(text) if text.strip() else {}
            except json.JSONDecodeError as exc:
                raise RuntimeError("Reading YAML config files requires PyYAML. Install with `pip install -e .`.") from exc
        if not isinstance(user_config, dict):
            raise ValueError(f"Config file must contain a YAML mapping: {path}")

    migrated_parameters = migrate_legacy_order_parameters(user_config)
    if migrated_parameters is not None:
        user_config.setdefault("order", {})["parameters"] = list(migrated_parameters)
    migrated_outputs = migrate_legacy_disabled_outputs(user_config)
    if migrated_outputs is not None:
        user_config.setdefault("output", {})["disabled_outputs"] = list(
            migrated_outputs
        )
    strip_legacy_selection_keys(user_config)

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






