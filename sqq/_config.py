from __future__ import annotations

"""Compatibility implementation for SQQ configuration handling."""

from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any
import warnings

try:
    import yaml
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    yaml = None


if yaml is not None:
    class _UniqueKeySafeLoader(yaml.SafeLoader):
        """Safe YAML loader that rejects duplicate mapping keys."""

        def construct_mapping(self, node, deep: bool = False):
            self.flatten_mapping(node)
            mapping: dict[Any, Any] = {}
            for key_node, value_node in node.value:
                key = self.construct_object(key_node, deep=deep)
                try:
                    duplicate = key in mapping
                except TypeError as exc:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping",
                        node.start_mark,
                        "found an unhashable mapping key",
                        key_node.start_mark,
                    ) from exc
                if duplicate:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping",
                        node.start_mark,
                        f"found duplicate key {key!r}",
                        key_node.start_mark,
                    )
                mapping[key] = self.construct_object(value_node, deep=deep)
            return mapping
else:  # pragma: no cover - used only by the JSON source-tree fallback.
    _UniqueKeySafeLoader = None


DEFAULT_MODE = "py"
CONFIG_SCHEMA_VERSION = "0.5.1"
DEFAULT_ORDER_PARAMETERS = ("f3", "f4")
CPP_MODE = "cpp"
CPP_MODES = frozenset({"99", CPP_MODE})
CPP_DEFAULT_OUTPUT_TYPES = (
    "info",
    "sqq-render",
    "summary-csv",
    "summary-detail-csv",
)
CPP_ALL_OUTPUT_TYPES = (
    "info",
    "gro",
    "cage-gro",
    "f3-gro",
    "f4-gro",
    "sqq-render",
    "summary-csv",
    "summary-xlsx",
    "summary-detail-csv",
)
CPP_OUTPUT_TYPES = frozenset(CPP_ALL_OUTPUT_TYPES)
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
    "sqq-render",
    "summary-xlsx",
)
ALL_OUTPUT_TYPES = (
    "info",
    "membership-tsv",
    "order-tsv",
    "f3-gro",
    "f4-gro",
    "gro",
    "sqq-render",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-gro",
    "cluster-detail",
)
OUTPUT_TYPE_ORDER = (
    "info",
    "membership-tsv",
    "order-tsv",
    "f3-gro",
    "f4-gro",
    "gro",
    "sqq-render",
    "ring-gro",
    "half-gro",
    "quasi-gro",
    "cage-gro",
    "ice-gro",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-gro",
    "cluster-detail",
)
GRO_OUTPUT_TYPES = {
    "ring-gro",
    "half-gro",
    "quasi-gro",
    "cage-gro",
    "ice-gro",
}

# YAML uses singular collection names. Runtime code keeps the established
# internal names so older callers and the analysis pipeline remain compatible.
_YAML_KEY_ALIASES: dict[tuple[str, ...], dict[str, str]] = {
    ("input", "lammps"): {
        "unit": "units",
    },
    ("graph",): {
        "mode": "bond_mode",
    },
    ("water",): {
        "resname": "resnames",
        "oxygen_name": "oxygen_names",
        "hydrogen_name": "hydrogen_names",
    },
    ("guest",): {
        "resname": "resnames",
        "center_atom": "center_atoms",
    },
    ("additive",): {
        "resname": "resnames",
    },
    ("environment",): {
        "resname": "resnames",
    },
    ("ring",): {
        "size": "sizes",
        "report_size": "report_sizes",
    },
    ("quasi_cage",): {
        "base_size": "base_sizes",
        "side_size": "side_sizes",
        "max_combination_per_base": "max_combinations_per_base",
        "max_layer": "max_layers",
        "max_ring_per_layer": "max_rings_per_layer",
        "max_layer_state_per_seed": "max_layer_states_per_seed",
        "max_candidate_per_edge": "max_candidates_per_edge",
        "max_layer_candidate": "max_layer_candidates",
    },
    ("cage",): {
        "report_type": "report_types",
        "max_face": "max_faces",
        "max_state_per_seed": "max_states_per_seed",
        "max_total_state": "max_total_states",
        "max_boundary_candidate": "max_boundary_candidates",
        "fast_closure_max_state": "fast_closure_max_states",
    },
    ("hydrate_order",): {
        "mcg_guest_resname": "mcg_guest_resnames",
        "mcg_min_water": "mcg_min_waters",
        "dhop_planar_count": "dhop_planar_counts",
        "dhop_min_qualified_neighbor": "dhop_min_qualified_neighbors",
    },
    ("order",): {
        "enabled": "parameters",
        "focus_water": "focus_waters",
    },
    ("ice",): {
        "min_six_ring": "min_six_rings",
        "require_four_coord_neighbor": "require_four_coord_neighbors",
    },
    ("output",): {
        "type": "types",
        "cage_isomer_row": "cage_isomer_rows",
        "write_empty_file": "write_empty_files",
        "context_role": "context_roles",
    },
    ("parallel",): {
        "worker": "workers",
        "math_thread": "math_threads",
    },
    ("track",): {
        "min_shared_water": "min_shared_waters",
    },
    ("debug",): {
        "use_networkx_check": "use_networkx_checks",
    },
}
MODE_PRESETS: dict[str, dict[str, Any]] = {
    "00": {
        "label": "rigorous",
        "worker_fraction": 1.0,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
        "find_cluster": True,
        "output_types": [
            "info",
            "sqq-render",
            "summary-xlsx",
        ],
    },
    "py": {
        "label": "standard",
        "worker_count": 1,
        "bond_mode": "auto",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": list(DEFAULT_OUTPUT_TYPES),
    },
    "99": {
        "label": "sqq-cpp-performance",
        "worker_fraction": 1.0,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": [
            "info",
            "sqq-render",
            "summary-csv",
            "summary-detail-csv",
        ],
    },
    CPP_MODE: {
        "label": "sqq-cpp",
        "worker_count": 1,
        "bond_mode": "auto",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": list(CPP_DEFAULT_OUTPUT_TYPES),
    },
}


# Explicit defaults keep the generated sqq_config.yaml reproducible.
DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": CONFIG_SCHEMA_VERSION,
    "mode": DEFAULT_MODE,
    "run": {
        "strict": False,
    },
    "input": {
        "pattern": "*.gro",
        "recursive": False,
        "first_file_time_ps": 0.0,
        "frame_time_step_ps": 100.0,
        "delta_time_ps": None,
        "xyz_scale": 0.1,
        "lammps": {
            "units": "real",
            "timestep": 1.0,
            "atom_style": "full",
            "coordinate_convention": "auto",
            "type_map": {},
        },
    },
    "component": {
        "auto_classify": True,
        "unknown_role": "other",
        "unknown_action": "warn",
        "role_map": {},
    },
    "additive": {
        "resnames": [],
    },
    "environment": {
        "resnames": [],
    },
    "water": {
        "resnames": ["SOL", "TIP", "WAT", "HOH"],
        "oxygen_names": ["OW", "O", "OH2"],
        "hydrogen_names": ["HW1", "HW2", "H1", "H2", "HW", "HT1", "HT2"],
    },
    "guest": {
        "resnames": ["CH4", "CO2", "MET", "ETH"],
        "center_atoms": {"CH4": ["C"], "CO2": ["C"], "MET": ["C"]},
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
        "sizes": [4, 5, 6],
        "report_sizes": "auto",
        "chordless": True,
        "definition": "chordless",
    },
    "half_cage": {
        "enabled": "auto",
    },
    "quasi_cage": {
        "search_policy": "bounded",
        "enabled": "auto",
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
        "max_states_per_seed": 0,
        "max_total_states": 0,
        "max_boundary_candidates": 8,
        "fast_closure": False,
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
        "summary_csv_dir": "summary",
        "cage_isomer_rows": "nonzero",
        "write_empty_files": False,
        "structure_layout": "grouped",
        "gro_atom_mode": "cage_oxygen_guest",
        "context_roles": [],
        "center_resname": "CNT",
    },
    "render": {
        "atom_scope": "full",
    },
    "parallel": {
        "backend": "process",
        "workers": "auto",
        "math_threads": 1,
    },
    "track": {
        "target": "all",
        "source": None,
        "min_jaccard": 0.50,
        "min_shared_fraction": 0.60,
        "min_shared_waters": 3,
        "max_center_distance_nm": None,
        "gap_frame": 0,
        "guest_tiebreak": True,
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
        raise ValueError(f"engine must be one of: {choices}")
    return text


def mode_label(mode: Any) -> str:
    """Return the human-readable mode label."""
    return str(MODE_PRESETS[normalize_mode(mode)]["label"])


def is_cpp_mode(mode: Any) -> bool:
    """Return whether the native C++ backend was selected."""
    return normalize_mode(mode) in CPP_MODES


def engine_display(selector: Any) -> str:
    """Return the effective public engine name for a selector."""
    normalized = normalize_mode(selector)
    return "sqq-cpp" if normalized in CPP_MODES else "sqq-py"


def mode_display(mode: Any) -> str:
    """Return the unified public mode/engine label."""
    normalized = normalize_mode(mode)
    if normalized in CPP_MODES:
        return "sqq-cpp" if normalized == CPP_MODE else f"{normalized} (sqq-cpp)"
    return f"{normalized} (sqq-py)"


def mode_worker_count(mode: Any) -> int | None:
    """Return a fixed default worker count, if the mode defines one."""
    value = MODE_PRESETS[normalize_mode(mode)].get("worker_count")
    return None if value is None else int(value)


def mode_worker_fraction(mode: Any) -> float:
    """Return the automatic worker fraction for a fraction-based mode."""
    preset = MODE_PRESETS[normalize_mode(mode)]
    if "worker_fraction" not in preset:
        raise ValueError(f"engine {normalize_mode(mode)} uses a fixed worker count")
    return float(preset["worker_fraction"])


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
    """Normalize the F3/F4-only selector used by engine cpp."""
    if value is None or value == "":
        raw_items: list[Any] = list(DEFAULT_ORDER_PARAMETERS)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "engine cpp --order-parameter must be f3, f4, f3,f4, all, or none."
            ) from exc
    if not raw_items:
        return ()
    cleaned = [str(item).strip().lower() for item in raw_items]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in engine cpp.")
        return DEFAULT_ORDER_PARAMETERS if cleaned[0] == "all" else ()
    normalized = normalize_order_parameters(cleaned)
    unsupported = [name for name in normalized if name not in {"f3", "f4"}]
    if unsupported:
        names = ", ".join(unsupported)
        raise ValueError(
            f"order parameter(s) {names} are not supported in engine cpp; use f3 and/or f4."
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
    """Normalize the positive SQQ-Py output allowlist."""
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
    exclusive = set(cleaned) & {"all", "none"}
    if exclusive:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in output.types / --output-type.")
        return ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()
    if "default" in cleaned:
        cleaned = [
            *DEFAULT_OUTPUT_TYPES,
            *(item for item in cleaned if item != "default"),
        ]

    supported = set(OUTPUT_TYPE_ORDER)
    unknown = sorted(set(cleaned) - supported)
    removed = {"sqq-cage-gro", "vmd"}.intersection(unknown)
    if removed:
        names = ", ".join(sorted(removed))
        raise ValueError(
            f"Output type(s) {names} were removed; use sqq-render instead."
        )
    if unknown:
        raise ValueError(
            f"Unsupported output type(s) {unknown}. Use default, info, membership-tsv, "
            "order-tsv, f3-gro, f4-gro, sqq-render, gro, ring-gro, half-gro, quasi-gro, cage-gro, "
            "ice-gro, cluster-gro, summary-xlsx, summary-csv, summary-detail-csv, "
            "cluster-detail, all, or none."
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


def migrate_yaml_keys(config: dict[str, Any]) -> dict[str, Any]:
    """Convert canonical singular YAML keys to stable runtime keys in place."""
    config.pop("schema_version", None)
    if "engine" in config:
        if "mode" in config:
            raise ValueError(
                "Config must not contain both top-level 'engine' and legacy 'mode'; "
                "use 'engine' only."
            )
        config["mode"] = config.pop("engine")
    elif "mode" in config:
        warnings.warn(
            "Top-level config key 'mode' is deprecated; use 'engine'.",
            UserWarning,
            stacklevel=2,
        )

    if "order_parameter" in config:
        if "order" in config:
            raise ValueError(
                "Config must not contain both top-level 'order_parameter' and "
                "legacy 'order'; use 'order_parameter' only."
            )
        config["order"] = config.pop("order_parameter")
    elif "order" in config:
        warnings.warn(
            "Top-level config key 'order' is deprecated; use 'order_parameter'.",
            UserWarning,
            stacklevel=2,
        )
    order = config.get("order")
    if isinstance(order, dict) and "parameter" in order:
        if "enabled" in order:
            raise ValueError(
                "Config must not contain both order_parameter.enabled and legacy "
                "order.parameter; use order_parameter.enabled only."
            )
        warnings.warn(
            "Config key order.parameter is deprecated; use order_parameter.enabled.",
            UserWarning,
            stacklevel=2,
        )
        order["enabled"] = order.pop("parameter")
    graph = config.get("graph")
    if isinstance(graph, dict) and "bond_mode" in graph and "mode" not in graph:
        warnings.warn(
            "Config key graph.bond_mode is deprecated; use graph.mode.",
            UserWarning,
            stacklevel=2,
        )
    if isinstance(graph, dict) and "pairs_file" in graph:
        if "pair_file" in graph:
            raise ValueError(
                "Config must not contain both graph.pair_file and legacy "
                "graph.pairs_file; use graph.pair_file only."
            )
        graph["pair_file"] = graph.pop("pairs_file")

    for path, aliases in _YAML_KEY_ALIASES.items():
        section = _nested_mapping(config, path)
        if section is None:
            continue
        for canonical, internal in aliases.items():
            if canonical in section:
                if internal in section and internal != canonical:
                    dotted = ".".join((*path, canonical))
                    legacy = ".".join((*path, internal))
                    raise ValueError(
                        f"Config must not contain both {dotted} and legacy {legacy}; "
                        f"use {dotted} only."
                    )
                section[internal] = section.pop(canonical)
    return config


def canonical_config(config: dict[str, Any]) -> dict[str, Any]:
    """Return a YAML-facing copy with canonical singular field names."""
    data = deepcopy(config)
    engine = data.pop("engine", None)
    internal_mode = data.pop("mode", None)
    if internal_mode is not None:
        engine = internal_mode
    if engine is None:
        engine = DEFAULT_MODE

    graph = data.get("graph")
    if isinstance(graph, dict) and "pairs_file" in graph:
        if "pair_file" not in graph:
            graph["pair_file"] = graph["pairs_file"]
        graph.pop("pairs_file", None)

    for path, aliases in _YAML_KEY_ALIASES.items():
        section = _nested_mapping(data, path)
        if section is None:
            continue
        reverse = {internal: canonical for canonical, internal in aliases.items()}
        converted: dict[str, Any] = {}
        for key, value in section.items():
            converted[reverse.get(key, key)] = value
        section.clear()
        section.update(converted)
    data.pop("schema_version", None)
    data = {
        ("order_parameter" if key == "order" else key): value
        for key, value in data.items()
    }
    return {
        "schema_version": CONFIG_SCHEMA_VERSION,
        "engine": normalize_mode(engine),
        **data,
    }


def _nested_mapping(
    config: dict[str, Any], path: tuple[str, ...]
) -> dict[str, Any] | None:
    section: Any = config
    for key in path:
        if not isinstance(section, dict):
            return None
        section = section.get(key)
    return section if isinstance(section, dict) else None


def normalize_engine_capabilities(
    config: dict[str, Any],
    *,
    user_config: dict[str, Any] | None = None,
    emit_warnings: bool = True,
) -> dict[str, Any]:
    """Resolve auto switches and safely disable unsupported C++ analyses."""
    cpp = is_cpp_mode(config.get("mode", DEFAULT_MODE))
    adjustments = list(config.get("adjustments", ()))

    half = config.setdefault("half_cage", {})
    half_value = _resolve_auto_toggle(half.get("enabled", "auto"), not cpp)
    quasi = config.setdefault("quasi_cage", {})
    quasi_value = _resolve_auto_toggle(quasi.get("enabled", "auto"), not cpp)
    cage = config.setdefault("cage", {})
    explicit_cage = _explicit_section(user_config, "cage")
    default_fast_limit = DEFAULT_CONFIG["cage"]["fast_closure_max_states"]
    if (
        _explicit_toggle_on(explicit_cage.get("fast_closure"))
        or explicit_cage.get("fast_closure_max_states", default_fast_limit)
        != default_fast_limit
    ):
        _record_capability_adjustment(
            adjustments,
            "legacy cage fast closure disabled; exact sparse GROW is used",
            emit_warnings,
        )
    cage["fast_closure"] = False
    cage["fast_closure_max_states"] = default_fast_limit

    if cpp:
        explicit_half = _explicit_section(user_config, "half_cage")
        explicit_quasi = _explicit_section(user_config, "quasi_cage")
        if _explicit_toggle_on(explicit_half.get("enabled")):
            _record_capability_adjustment(
                adjustments,
                "half_cage disabled by SQQ-CPP",
                emit_warnings,
            )
        quasi_customized = any(
            key != "enabled"
            and value != DEFAULT_CONFIG["quasi_cage"].get(key)
            for key, value in explicit_quasi.items()
        )
        if _explicit_toggle_on(explicit_quasi.get("enabled")) or quasi_customized:
            _record_capability_adjustment(
                adjustments,
                "quasi_cage disabled by SQQ-CPP; quasi settings ignored",
                emit_warnings,
            )
        _remove_cpp_yaml_outputs(config, adjustments, emit_warnings)
        half_value = False
        quasi_value = False
        for key, value in DEFAULT_CONFIG["quasi_cage"].items():
            if key != "enabled":
                quasi[key] = deepcopy(value)

    half["enabled"] = half_value
    quasi["enabled"] = quasi_value
    if adjustments:
        config["adjustments"] = list(dict.fromkeys(str(item) for item in adjustments))
    return config


def _remove_cpp_yaml_outputs(
    config: dict[str, Any], adjustments: list[Any], emit_warnings: bool
) -> None:
    output = config.setdefault("output", {})
    value = output.get("types", ())
    if isinstance(value, str):
        items = [item.strip().lower() for item in value.split(",") if item.strip()]
    else:
        try:
            items = [str(item).strip().lower() for item in value if str(item).strip()]
        except TypeError:
            return
    unsupported = [name for name in ("half-gro", "quasi-gro") if name in items]
    if not unsupported:
        return
    output["types"] = [name for name in items if name not in set(unsupported)]
    _record_capability_adjustment(
        adjustments,
        "SQQ-CPP removed unsupported YAML output type(s): " + ", ".join(unsupported),
        emit_warnings,
    )

def _resolve_auto_toggle(value: Any, auto_value: bool) -> bool:
    if isinstance(value, str) and value.strip().lower() == "auto":
        return auto_value
    return legacy_enabled(value)


def _explicit_section(
    config: dict[str, Any] | None, section_name: str
) -> dict[str, Any]:
    if not isinstance(config, dict):
        return {}
    section = config.get(section_name, {})
    return section if isinstance(section, dict) else {}


def _explicit_toggle_on(value: Any) -> bool:
    if value is None or (isinstance(value, str) and value.strip().lower() == "auto"):
        return False
    return legacy_enabled(value)


def _record_capability_adjustment(
    adjustments: list[Any], message: str, emit_warning: bool
) -> None:
    if message not in adjustments:
        adjustments.append(message)
    if emit_warning:
        warnings.warn(message, UserWarning, stacklevel=3)

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
    config["output"]["types"] = list(preset["output_types"])
    if is_cpp_mode(normalized):
        config["half_cage"]["enabled"] = False
        config["quasi_cage"]["enabled"] = False
        config["ice"]["enabled"] = False
        config["order"]["parameters"] = list(DEFAULT_ORDER_PARAMETERS)
        config["cage"]["fast_closure"] = False
    return config


_RESOLVED_RUN_KEYS = {
    "sqq_version",
    "release_date",
    "engine_selector",
    "engine",
    "config_output",
    "status",
    "error",
    "graph_mode_requested",
    "graph_mode_effective",
    "graph_mode_display",
    "graph_mode_by_group",
    "order_parameters",
    "find_half",
    "find_quasi",
    "find_cluster",
    "output_types",
    "sampling_interval",
    "native_frame_interval_ps",
    "delta_time_ps",
    "raw_frame_step",
    "selected_frames",
    "source_frames_total",
    "frames_total",
    "frames_ok",
    "frames_failed",
    "failures",
    "worker_request",
    "worker_policy",
    "workers_resolved",
    "parallel_backend",
    "math_threads_per_worker",
    "summary_write",
    "topology_group_count",
    "topology_group_limit",
    "topology_group_limit_exceeded",
    "topology_group_labels_enabled",
    "info_only_fallback_required",
    "topology_grouping",
    "topology_groups",
    "topology_source_mapping",
    "topology_group",
    "topology_fingerprint",
    "requested_output_types",
    "output_policy",
    "warnings",
}
_CONFIG_EXTRA_KEYS: dict[tuple[str, ...], set[str]] = {
    (): {"adjustments"},
    ("run",): _RESOLVED_RUN_KEYS,
    ("input",): {"format", "topology", "sampling"},
    ("input", "lammps"): {
        "resolved_type_map",
        "type_map_source",
        "rebuilt_molecules",
    },
    ("graph",): {"effective_bond_mode", "effective_bond_mode_by_group"},
}
_CONFIG_FLEXIBLE_PATHS = {
    ("component", "role_map"),
    ("guest", "center_atoms"),
    ("input", "sampling"),
    ("input", "lammps", "type_map"),
    ("input", "lammps", "resolved_type_map"),
    ("graph", "effective_bond_mode_by_group"),
    ("run", "graph_mode_by_group"),
    ("run", "summary_write"),
}


def validate_user_config_keys(
    config: dict[str, Any],
    schema: dict[str, Any] | None = None,
    path: tuple[str, ...] = (),
) -> None:
    """Reject unknown YAML keys after canonical/legacy key migration."""
    if path in _CONFIG_FLEXIBLE_PATHS:
        return
    expected = DEFAULT_CONFIG if schema is None else schema
    allowed_extra = _CONFIG_EXTRA_KEYS.get(path, set())
    for key, value in config.items():
        if key not in expected:
            if key in allowed_extra:
                continue
            dotted = ".".join((*path, str(key)))
            raise ValueError(f"Unknown configuration key: {dotted}")
        expected_value = expected[key]
        if isinstance(expected_value, dict):
            if not isinstance(value, dict):
                dotted = ".".join((*path, str(key)))
                raise ValueError(f"Configuration section must be a mapping: {dotted}")
            validate_user_config_keys(value, expected_value, (*path, str(key)))

def load_config(path: Path | None, mode: Any = None) -> dict[str, Any]:
    """Load an engine preset, then merge user configuration over it."""
    user_config: dict[str, Any] = {}
    if path is not None:
        with path.open("r", encoding="utf-8-sig") as handle:
            text = handle.read()
        if yaml is not None:
            user_config = yaml.load(text, Loader=_UniqueKeySafeLoader) or {}
        else:
            # JSON fallback supports source-tree smoke tests without PyYAML.
            try:
                user_config = json.loads(text) if text.strip() else {}
            except json.JSONDecodeError as exc:
                raise RuntimeError("Reading YAML config files requires PyYAML. Install with `pip install -e .`.") from exc
        if not isinstance(user_config, dict):
            raise ValueError(f"Config file must contain a YAML mapping: {path}")

    migrate_yaml_keys(user_config)
    if path is not None:
        graph = user_config.get("graph", {})
        pair_file = graph.get("pair_file") if isinstance(graph, dict) else None
        if pair_file not in (None, ""):
            pair_path = Path(str(pair_file)).expanduser()
            if not pair_path.is_absolute():
                pair_path = path.resolve().parent / pair_path
            graph["pair_file"] = str(pair_path.resolve())
    migrated_parameters = migrate_legacy_order_parameters(user_config)
    if migrated_parameters is not None:
        user_config.setdefault("order", {})["parameters"] = list(migrated_parameters)
    strip_legacy_selection_keys(user_config)
    validate_user_config_keys(user_config)

    selected_mode = normalize_mode(mode if mode is not None else user_config.get("mode", DEFAULT_MODE))
    config = apply_mode_preset(deepcopy(DEFAULT_CONFIG), selected_mode)
    merge_config(config, user_config)
    config["mode"] = selected_mode
    config["schema_version"] = CONFIG_SCHEMA_VERSION
    render = config.setdefault("render", {})
    atom_scope = str(render.get("atom_scope", "full")).strip().lower()
    if atom_scope not in {"full", "compact"}:
        raise ValueError("render.atom_scope must be full or compact.")
    render["atom_scope"] = atom_scope
    normalize_engine_capabilities(config, user_config=user_config)
    return config


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge user configuration into defaults."""
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge_config(base[key], value)
        else:
            base[key] = value
    return base


_DEFAULT_CONFIG_SECTION_COMMENTS: dict[str, str] = {
    "run": "Run behavior",
    "input": "Input discovery, time sampling, and LAMMPS reader settings",
    "component": "Automatic component classification",
    "additive": "Additive residue names",
    "environment": "Environment or wall residue names",
    "water": "Water selection",
    "guest": "Guest selection",
    "graph": "Water-network construction",
    "pbc": "Periodic-boundary settings",
    "ring": "Ring search and reporting",
    "half_cage": "Standard half-cage search (SQQ-Py)",
    "quasi_cage": "Layered quasi-cage search (SQQ-Py)",
    "cage": "Complete-cage recognition",
    "hydrate_cluster": "Hydrate phase and cluster recognition (SQQ-Py)",
    "hydrate_order": "MCG and DHOP definitions",
    "order_parameter": "Order-parameter calculation",
    "ice": "Ice-like water classification (SQQ-Py)",
    "output": "Output selection and layout",
    "render": "VMD render topology and trajectory",
    "parallel": "Worker and math-thread policy",
    "track": "Cross-frame persistent cage tracking",
    "debug": "Developer diagnostics",
}

_DEFAULT_CONFIG_INLINE_COMMENTS: dict[tuple[str, ...], str] = {
    ("schema_version",): "managed by SQQ",
    ("engine",): "choices: py, cpp; compatibility presets: 00, 99",
    ("run", "strict"): "choices: true, false",
    ("input", "pattern"): "glob used for directory input",
    ("input", "recursive"): "choices: true, false",
    ("input", "first_file_time_ps"): "ps; independent GRO series",
    ("input", "frame_time_step_ps"): "ps; independent GRO series",
    ("input", "delta_time_ps"): "ps; null analyzes every stored frame",
    ("input", "xyz_scale"): "coordinate-to-nm scale",
    ("input", "lammps", "unit"): "choices: real, metal, nano",
    ("input", "lammps", "timestep"): "LAMMPS integration timestep",
    ("input", "lammps", "atom_style"): "choices: full, molecular, bond, angle",
    ("input", "lammps", "coordinate_convention"): "choices: auto, x, xs, xu, xsu, unscaled, scaled, unwrapped, scaled_unwrapped",
    ("component", "unknown_action"): "choices: warn, ignore, error",
    ("guest", "center_mode"): "choices: center_atom, centroid, auto",
    ("graph", "mode"): "choices: auto, hbond, oo, pairs",
    ("graph", "pair_file"): "required when graph.mode is pairs",
    ("graph", "pair_id"): "choices: resid, oxygen_index, atomid",
    ("pbc", "box_mode"): "currently orthorhombic only",
    ("ring", "size"): "supported sizes, normally [4, 5, 6]",
    ("ring", "report_size"): "auto or a subset of ring.size",
    ("ring", "definition"): "choices: chordless, shortest_path",
    ("half_cage", "enabled"): "choices: auto, true, false",
    ("quasi_cage", "enabled"): "choices: auto, true, false",
    ("quasi_cage", "search_policy"): "choices: bounded, exact",
    ("quasi_cage", "max_layer"): "maximum layered growth depth; minimum 1",
    ("cage", "report_type"): "auto, all, I, II, H, HS-I, TS-I, I2II, or cage labels",
    ("cage", "max_state_per_seed"): "0 = unlimited; positive values are hard diagnostic guards",
    ("cage", "max_total_state"): "0 = unlimited; positive values are hard diagnostic guards",
    ("cage", "max_boundary_candidate"): "legacy compatibility key; exact search never truncates candidates",
    ("cage", "fast_closure"): "legacy compatibility key; normalized to false",
    ("cage", "fast_closure_max_state"): "legacy compatibility key; ignored by exact search",
    ("cage", "scientific_validation"): "choices: true, false",
    ("hydrate_cluster", "enabled"): "choices: true, false",
    ("order_parameter", "enabled"): "choices: f3, f4, qN, mcg1, mcg3, dhop35, dhop30, all, none",
    ("order_parameter", "q_neighbor_mode"): "choices: graph, cutoff, nearest, lammps",
    ("output", "type"): "default may be combined; all/none are exclusive",
    ("output", "structure_layout"): "choices: grouped, flat",
    ("render", "atom_scope"): "choices: full, compact",
    ("parallel", "backend"): "choices: process, thread, serial",
    ("parallel", "worker"): "auto, integer count, fraction, or percentage",
    ("parallel", "math_thread"): "numeric-library threads per worker",
    ("track", "target"): "all, cage type, hydrate phase, persistent tID, or comma-separated targets",
    ("track", "source"): "Analyze result directory; null uses the current directory",
    ("track", "min_jaccard"): "minimum member-water Jaccard similarity in [0, 1]",
    ("track", "min_shared_fraction"): "minimum smaller-cage shared-water fraction in [0, 1]",
    ("track", "min_shared_water"): "minimum number of shared water molecules",
    ("track", "max_center_distance_nm"): "optional orthorhombic-PBC center-distance guard",
    ("track", "gap_frame"): "maximum explicitly recorded missing selected frames; 0 disables bridging",
    ("track", "guest_tiebreak"): "use guest continuity only to break otherwise comparable matches",
}


def default_config_template() -> str:
    """Return the commented, round-trippable configuration template."""
    payload = canonical_config(DEFAULT_CONFIG)
    if yaml is None:
        return json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    raw = yaml.safe_dump(payload, allow_unicode=True, sort_keys=False)
    output: list[str] = []
    path_by_depth: dict[int, str] = {}
    key_pattern = re.compile(r"^(?P<indent> *)(?P<key>[A-Za-z0-9_]+):(?P<value>.*)$")
    for line in raw.splitlines():
        match = key_pattern.match(line)
        if match is None:
            output.append(line)
            continue
        depth = len(match.group("indent")) // 2
        key = match.group("key")
        path_by_depth[depth] = key
        for stale_depth in tuple(path_by_depth):
            if stale_depth > depth:
                del path_by_depth[stale_depth]
        path = tuple(path_by_depth[index] for index in range(depth + 1))
        if depth == 0 and key in _DEFAULT_CONFIG_SECTION_COMMENTS:
            if output and output[-1] != "":
                output.append("")
            output.append(f"# {_DEFAULT_CONFIG_SECTION_COMMENTS[key]}")
        comment = _DEFAULT_CONFIG_INLINE_COMMENTS.get(path)
        output.append(f"{line}  # {comment}" if comment else line)
    return "\n".join(output) + "\n"


def write_default_config(path: Path) -> None:
    """Write the commented default configuration without overwriting a file."""
    if path.exists():
        raise FileExistsError(f"Configuration file already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(default_config_template(), encoding="utf-8", newline="\n")

def dump_config(config: dict[str, Any], handle) -> None:
    """Write YAML when available, otherwise a JSON-compatible fallback."""
    payload = canonical_config(config)
    if yaml is not None:
        yaml.safe_dump(payload, handle, allow_unicode=True, sort_keys=False)
    else:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")



def normalize_cpp_output_types(value: Any = None) -> tuple[str, ...]:
    """Normalize the compact output allowlist used by engine cpp."""
    if value is None:
        raw_items: list[Any] = list(CPP_DEFAULT_OUTPUT_TYPES)
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "SQQ-CPP --output-type contains an invalid output name."
            ) from exc
    if not raw_items:
        return ()
    cleaned = [str(item).strip().lower() for item in raw_items]
    exclusive = set(cleaned) & {"all", "none"}
    if exclusive:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in engine cpp.")
        return CPP_ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()
    if "default" in cleaned:
        cleaned = [
            *CPP_DEFAULT_OUTPUT_TYPES,
            *(item for item in cleaned if item != "default"),
        ]
    unsupported = sorted(set(cleaned) - CPP_OUTPUT_TYPES)
    removed = {"sqq-cage-gro", "vmd"}.intersection(unsupported)
    if removed:
        names = ", ".join(sorted(removed))
        raise ValueError(
            f"Output type(s) {names} were removed; use sqq-render instead."
        )
    if unsupported:
        names = ", ".join(unsupported)
        raise ValueError(
            f"output type(s) {names} are not supported by SQQ-CPP; "
            "use default, info, gro, cage-gro, f3-gro, f4-gro, sqq-render, summary-csv, "
            "summary-xlsx, summary-detail-csv, all, or none."
        )
    normalized = set(cleaned)
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
        ("cluster_min_cage", "--cluster-min-cage"),
    )
    for attribute, option in unsupported_args:
        value = getattr(args, attribute, None)
        if value not in (None, False):
            errors.append(f"{option} is not supported in engine cpp")
            explicit_unsupported.add(attribute)

    find_half = getattr(args, "find_half", None)
    if find_half not in (None, False, "off"):
        errors.append("--find-half on is not supported by SQQ-CPP")
    find_quasi = getattr(args, "find_quasi", None)
    if find_quasi not in (None, False, "off"):
        errors.append("--find-quasi on is not supported by SQQ-CPP")

    find_cluster = getattr(args, "find_cluster", None)
    if find_cluster not in (None, False, "off"):
        errors.append("--find-cluster on is not supported by SQQ-CPP")
        explicit_unsupported.add("find_cluster")

    ring = config.setdefault("ring", {})
    explicit_ring_definition = getattr(args, "ring_definition", None) not in (None, "chordless")
    if explicit_ring_definition:
        errors.append("--ring-definition shortest_path is not supported in engine cpp")
    if not explicit_ring_definition and str(ring.get("definition", "chordless")).strip().lower() != "chordless":
        errors.append("ring.definition must be chordless in engine cpp")
    if not bool(ring.get("chordless", True)):
        errors.append("ring.chordless=false is not supported in engine cpp")
    if "ring_size" not in explicit_unsupported and ring.get("report_sizes", "auto") not in (None, "", "auto"):
        errors.append("public ring reporting is not supported in engine cpp")
    try:
        ring["sizes"] = _normalize_cpp_ring_sizes(ring.get("sizes", (4, 5, 6)))
    except (TypeError, ValueError) as exc:
        errors.append(str(exc))

    half = config.setdefault("half_cage", {})
    if _cpp_requested_on(half.get("enabled", False)):
        errors.append("half_cage.enabled is not supported in engine cpp")
    half["enabled"] = False

    quasi = config.setdefault("quasi_cage", {})
    if getattr(args, "size", None):
        quasi["base_sizes"] = "auto"
        quasi["side_sizes"] = "auto"
    explicit_quasi = any(name.startswith("quasi_") for name in explicit_unsupported)
    if not explicit_quasi:
        for key, default in DEFAULT_CONFIG["quasi_cage"].items():
            if key != "enabled" and quasi.get(key, default) != default:
                errors.append(f"quasi_cage.{key} is not supported in engine cpp")
    if _cpp_requested_on(quasi.get("enabled", False)):
        errors.append("quasi_cage.enabled is not supported in engine cpp")
    quasi["enabled"] = False

    cage = config.setdefault("cage", {})
    if str(cage.get("search_mode", "grow")).strip().lower() != "grow":
        errors.append("cage.search_mode must be grow in engine cpp")
    if str(cage.get("seed_mode", "ring")).strip().lower() != "ring":
        errors.append("cage.seed_mode must be ring in engine cpp")
    default_fast_limit = DEFAULT_CONFIG["cage"]["fast_closure_max_states"]
    if cage.get("fast_closure_max_states", default_fast_limit) != default_fast_limit:
        errors.append("cage.fast_closure_max_states is not supported in engine cpp")
    cage["fast_closure"] = False

    cluster = config.setdefault("hydrate_cluster", {})
    if "find_cluster" not in explicit_unsupported and _cpp_requested_on(cluster.get("enabled", False)):
        errors.append("hydrate_cluster.enabled is not supported in engine cpp")
    default_min_cage = DEFAULT_CONFIG["hydrate_cluster"]["min_cage"]
    if (
        "cluster_min_cage" not in explicit_unsupported
        and cluster.get("min_cage", default_min_cage) != default_min_cage
    ):
        errors.append("hydrate_cluster.min_cage is not supported in engine cpp")
    cluster["enabled"] = False

    ice = config.setdefault("ice", {})
    for key, default in DEFAULT_CONFIG["ice"].items():
        if key != "enabled" and ice.get(key, default) != default:
            errors.append(f"ice.{key} is not supported in engine cpp")
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
            errors.append(f"order.{key} is not supported in engine cpp")
    for key, default in DEFAULT_CONFIG["hydrate_order"].items():
        if config.get("hydrate_order", {}).get(key, default) != default:
            errors.append(f"hydrate_order.{key} is not supported in engine cpp")

    output = config.setdefault("output", {})
    output_source = getattr(args, "output_type", None)
    if output_source is None:
        output_source = output.get("types", CPP_DEFAULT_OUTPUT_TYPES)
    try:
        output["types"] = list(normalize_cpp_output_types(output_source))
    except ValueError as exc:
        errors.append(str(exc))

    parallel_backend = str(config.get("parallel", {}).get("backend", "process"))
    if parallel_backend.strip().lower() == "thread":
        errors.append("--parallel-backend thread is not supported in engine cpp")
    if errors:
        raise ValueError("; ".join(dict.fromkeys(errors)))


def _normalize_cpp_ring_sizes(value: Any) -> list[int]:
    if isinstance(value, str):
        raw = [item.strip() for item in value.split(",") if item.strip()]
    else:
        raw = list(value)
    sizes = sorted({int(item) for item in raw})
    if not sizes:
        raise ValueError("engine cpp requires at least one ring size from 4, 5, and 6")
    unsupported = [size for size in sizes if size not in {4, 5, 6}]
    if unsupported:
        names = ", ".join(str(size) for size in unsupported)
        raise ValueError(f"ring size(s) {names} are not supported in engine cpp")
    return sizes


def _cpp_requested_on(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off", "none"}
    return bool(value)


