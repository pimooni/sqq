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
CPP_MODE = "cpp"
CPP_DEFAULT_OUTPUT_TYPES = ("info", "cage-gro", "summary-csv")
CPP_ALL_OUTPUT_TYPES = ("info", "cage-gro", "summary-csv", "summary-xlsx")
CPP_OUTPUT_TYPES = frozenset({"info", "gro", "cage-gro", "summary-csv", "summary-xlsx"})
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
DEFAULT_OUTPUT_TYPES = (
    "info",
    "gro",
    "summary-xlsx",
)
ALL_OUTPUT_TYPES = (
    "info",
    "membership-tsv",
    "order-tsv",
    "vmd",
    "gro",
    "cluster-gro",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-detail",
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
    "cluster-gro",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-detail",
)
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
        "find_cluster": True,
    },
    "09": {
        "label": "rigorous-performance",
        "worker_fraction": 0.90,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
        "find_cluster": True,
    },
    "50": {
        "label": "standard",
        "worker_fraction": 0.50,
        "bond_mode": "auto",
        "ring_sizes": [5, 6],
        "find_cluster": False,
    },
    "99": {
        "label": "performance",
        "worker_fraction": 0.90,
        "bond_mode": "oo",
        "ring_sizes": [5, 6],
        "find_cluster": False,
    },
    CPP_MODE: {
        "label": "sqq-cpp",
        "worker_fraction": 0.90,
        "bond_mode": "auto",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
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
        "types": list(DEFAULT_OUTPUT_TYPES),
        "summary_csv_dir": "summary_csv",
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
    """Normalize and validate an analysis mode."""
    text = str(value).strip().lower()
    if text.isdigit():
        text = text.zfill(2)
    if text not in MODE_PRESETS:
        choices = ", ".join(MODE_PRESETS)
        raise ValueError(f"mode must be one of: {choices}")
    return text


def mode_label(mode: Any) -> str:
    """Return the human-readable mode label."""
    return str(MODE_PRESETS[normalize_mode(mode)]["label"])


def is_cpp_mode(mode: Any) -> bool:
    """Return whether the native C++ backend was selected."""
    return normalize_mode(mode) == CPP_MODE


def mode_display(mode: Any) -> str:
    """Return the unified public mode/engine label."""
    normalized = normalize_mode(mode)
    return "sqq-cpp" if normalized == CPP_MODE else f"{normalized} (sqq-py)"


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


def normalize_cpp_order_parameters(value: Any = None) -> tuple[str, ...]:
    """Normalize the F3/F4-only selector used by mode cpp."""
    if value is None or value == "":
        raw_items: list[Any] = list(DEFAULT_ORDER_PARAMETERS)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "mode cpp --order-parameter must be f3, f4, f3,f4, all, or none."
            ) from exc
    if not raw_items:
        return ()
    cleaned = [str(item).strip().lower() for item in raw_items]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in mode cpp.")
        return DEFAULT_ORDER_PARAMETERS if cleaned[0] == "all" else ()
    normalized = normalize_order_parameters(cleaned)
    unsupported = [name for name in normalized if name not in {"f3", "f4"}]
    if unsupported:
        names = ", ".join(unsupported)
        raise ValueError(
            f"order parameter(s) {names} are not supported in mode cpp; use f3 and/or f4."
        )
    return normalized


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


def normalize_output_types(value: Any = None) -> tuple[str, ...]:
    """Normalize the positive output allowlist."""
    if value is None:
        raw_items: list[Any] = list(DEFAULT_OUTPUT_TYPES)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "output.types / --output-type must be a comma-separated list."
            ) from exc
    if not raw_items:
        return ()

    cleaned = [str(item).strip().lower() for item in raw_items]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in output.types / --output-type.")
        return ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()

    supported = set(OUTPUT_TYPE_ORDER)
    unknown = sorted(set(cleaned) - supported)
    if unknown:
        raise ValueError(
            f"Unsupported output type(s) {unknown}. Use info, membership-tsv, "
            "order-tsv, vmd, gro, ring-gro, half-gro, quasi-gro, cage-gro, "
            "ice-gro, cluster-gro, summary-xlsx, summary-csv, summary-detail-csv, cluster-detail, all, or none."
        )
    normalized = set(cleaned)
    if "gro" in normalized:
        normalized.difference_update(GRO_OUTPUT_TYPES)
    return tuple(name for name in OUTPUT_TYPE_ORDER if name in normalized)


def output_type_display(value: Any, *, cpp_mode: bool = False) -> str:
    """Render normalized output types for terminal and metadata."""
    outputs = (
        normalize_cpp_output_types(value)
        if cpp_mode
        else normalize_output_types(value)
    )
    return ", ".join(outputs) if outputs else "none"


def output_enabled(config: dict[str, Any], output_type: str) -> bool:
    """Return whether one output category is selected."""
    normalized_type = str(output_type).strip().lower()
    if normalized_type not in set(OUTPUT_TYPE_ORDER):
        raise ValueError(f"Unsupported output type: {output_type}")
    enabled = set(
        normalize_output_types(
            config.get("output", {}).get("types", DEFAULT_OUTPUT_TYPES)
        )
    )
    if normalized_type in GRO_OUTPUT_TYPES and "gro" in enabled:
        return True
    return normalized_type in enabled


def strip_legacy_selection_keys(user_config: dict[str, Any]) -> None:
    """Remove migrated order-selector booleans."""
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
    config["hydrate_cluster"]["enabled"] = bool(preset["find_cluster"])
    if normalized == CPP_MODE:
        config["quasi_cage"]["enabled"] = False
        config["ice"]["enabled"] = False
        config["order"]["parameters"] = list(DEFAULT_ORDER_PARAMETERS)
        config["output"]["types"] = list(CPP_DEFAULT_OUTPUT_TYPES)
        config["cage"]["fast_closure"] = False
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



def normalize_cpp_output_types(value: Any = None) -> tuple[str, ...]:
    """Normalize the compact output allowlist used by mode cpp."""
    if value is None:
        raw_items: list[Any] = list(CPP_DEFAULT_OUTPUT_TYPES)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "mode cpp --output-type must contain info, gro, cage-gro, summary-csv, or summary-xlsx."
            ) from exc
    if not raw_items:
        return ()
    cleaned = [str(item).strip().lower() for item in raw_items]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in mode cpp.")
        return CPP_ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()
    unsupported = sorted(set(cleaned) - CPP_OUTPUT_TYPES)
    if unsupported:
        names = ", ".join(unsupported)
        raise ValueError(
            f"output type(s) {names} are not supported in mode cpp; "
            "use info, gro, cage-gro, summary-csv, summary-xlsx, all, or none."
        )
    normalized = set(cleaned)
    if "gro" in normalized:
        normalized.remove("gro")
        normalized.add("cage-gro")
    return tuple(name for name in CPP_ALL_OUTPUT_TYPES if name in normalized)


def validate_cpp_cli(args: Any, config: dict[str, Any]) -> None:
    """Reject unsupported settings and finalize the mode-cpp subset."""
    if not is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        return
    errors: list[str] = []
    explicit_unsupported: set[str] = set()
    unsupported_args = (
        ("ring_size", "--ring-size"),
        ("quasi_size", "--quasi-size"),
        ("quasi_base_size", "--quasi-base-size"),
        ("quasi_side_size", "--quasi-side-size"),
        ("quasi_max_layer", "--quasi-max-layer"),
        ("quasi_search_policy", "--quasi-search-policy"),
        ("no_q", "--no-q"),
        ("q_degree", "--q-degree"),
        ("q_neighbor_mode", "--q-neighbor-mode"),
        ("q_cutoff", "--q-cutoff"),
        ("q_n_neighbor", "--q-n-neighbor"),
        ("mcg3", "--mcg3"),
        ("dhop30", "--dhop30"),
        ("cage_fast_closure", "--cage-fast-closure"),
        ("find_cluster", "--find-cluster"),
        ("cluster_min_cage", "--cluster-min-cage"),
    )
    for attribute, option in unsupported_args:
        value = getattr(args, attribute, None)
        if value not in (None, False):
            errors.append(f"{option} is not supported in mode cpp")
            explicit_unsupported.add(attribute)

    ring = config.setdefault("ring", {})
    explicit_ring_definition = getattr(args, "ring_definition", None) not in (None, "chordless")
    if explicit_ring_definition:
        errors.append("--ring-definition shortest_path is not supported in mode cpp")
    if not explicit_ring_definition and str(ring.get("definition", "chordless")).strip().lower() != "chordless":
        errors.append("ring.definition must be chordless in mode cpp")
    if not bool(ring.get("chordless", True)):
        errors.append("ring.chordless=false is not supported in mode cpp")
    if "ring_size" not in explicit_unsupported and ring.get("report_sizes", "auto") not in (None, "", "auto"):
        errors.append("public ring reporting is not supported in mode cpp")
    try:
        ring["sizes"] = _normalize_cpp_ring_sizes(ring.get("sizes", (4, 5, 6)))
    except (TypeError, ValueError) as exc:
        errors.append(str(exc))

    quasi = config.setdefault("quasi_cage", {})
    if getattr(args, "size", None):
        quasi["base_sizes"] = "auto"
        quasi["side_sizes"] = "auto"
    explicit_quasi = any(name.startswith("quasi_") for name in explicit_unsupported)
    if not explicit_quasi:
        for key, default in DEFAULT_CONFIG["quasi_cage"].items():
            if key != "enabled" and quasi.get(key, default) != default:
                errors.append(f"quasi_cage.{key} is not supported in mode cpp")
    quasi["enabled"] = False

    cage = config.setdefault("cage", {})
    if str(cage.get("search_mode", "grow")).strip().lower() != "grow":
        errors.append("cage.search_mode must be grow in mode cpp")
    if str(cage.get("seed_mode", "ring")).strip().lower() != "ring":
        errors.append("cage.seed_mode must be ring in mode cpp")
    default_fast_limit = DEFAULT_CONFIG["cage"]["fast_closure_max_states"]
    if cage.get("fast_closure_max_states", default_fast_limit) != default_fast_limit:
        errors.append("cage.fast_closure_max_states is not supported in mode cpp")
    cage["fast_closure"] = False

    cluster = config.setdefault("hydrate_cluster", {})
    if "find_cluster" not in explicit_unsupported and _cpp_requested_on(cluster.get("enabled", False)):
        errors.append("hydrate_cluster.enabled is not supported in mode cpp")
    default_min_cage = DEFAULT_CONFIG["hydrate_cluster"]["min_cage"]
    if (
        "cluster_min_cage" not in explicit_unsupported
        and cluster.get("min_cage", default_min_cage) != default_min_cage
    ):
        errors.append("hydrate_cluster.min_cage is not supported in mode cpp")
    cluster["enabled"] = False

    ice = config.setdefault("ice", {})
    for key, default in DEFAULT_CONFIG["ice"].items():
        if key != "enabled" and ice.get(key, default) != default:
            errors.append(f"ice.{key} is not supported in mode cpp")
    ice["enabled"] = False

    order = config.setdefault("order", {})
    explicit_legacy_order = bool(
        {"no_q", "q_degree", "mcg3", "dhop30"} & explicit_unsupported
    )
    order_source = getattr(args, "order_parameter", None)
    if order_source is None and explicit_legacy_order:
        order_source = DEFAULT_ORDER_PARAMETERS
    elif order_source is None:
        order_source = order.get("parameters", DEFAULT_ORDER_PARAMETERS)
    try:
        order["parameters"] = list(normalize_cpp_order_parameters(order_source))
    except ValueError as exc:
        errors.append(str(exc))
    q_option_for_key = {
        "q_neighbor_mode": "q_neighbor_mode",
        "q_cutoff_nm": "q_cutoff",
        "q_n_neighbor": "q_n_neighbor",
    }
    for key, option_name in q_option_for_key.items():
        default = DEFAULT_CONFIG["order"][key]
        if option_name not in explicit_unsupported and order.get(key, default) != default:
            errors.append(f"order.{key} is not supported in mode cpp")
    for key, default in DEFAULT_CONFIG["hydrate_order"].items():
        if config.get("hydrate_order", {}).get(key, default) != default:
            errors.append(f"hydrate_order.{key} is not supported in mode cpp")

    output = config.setdefault("output", {})
    output_source = getattr(args, "output_type", None)
    if output_source is None:
        output_source = output.get("types", CPP_DEFAULT_OUTPUT_TYPES)
    try:
        output["types"] = list(normalize_cpp_output_types(output_source))
    except ValueError as exc:
        errors.append(str(exc))
    default_detail_dir = DEFAULT_CONFIG["output"]["summary_detail_dir"]
    if output.get("summary_detail_dir", default_detail_dir) != default_detail_dir:
        errors.append("output.summary_detail_dir is not supported in mode cpp")

    parallel_backend = str(config.get("parallel", {}).get("backend", "process"))
    if parallel_backend.strip().lower() == "thread":
        errors.append("--parallel-backend thread is not supported in mode cpp")
    if errors:
        raise ValueError("; ".join(dict.fromkeys(errors)))


def _normalize_cpp_ring_sizes(value: Any) -> list[int]:
    if isinstance(value, str):
        raw = [item.strip() for item in value.split(",") if item.strip()]
    else:
        raw = list(value)
    sizes = sorted({int(item) for item in raw})
    if not sizes:
        raise ValueError("mode cpp requires at least one ring size from 4, 5, and 6")
    unsupported = [size for size in sizes if size not in {4, 5, 6}]
    if unsupported:
        names = ", ".join(str(size) for size in unsupported)
        raise ValueError(f"ring size(s) {names} are not supported in mode cpp")
    return sizes


def _cpp_requested_on(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off", "none"}
    return bool(value)


