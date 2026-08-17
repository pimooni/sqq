"""Configuration value normalization and validation."""

from __future__ import annotations

import re
from typing import Any

from .defaults import (
    ALL_ORDER_PARAMETERS,
    ALL_OUTPUT_TYPES,
    CPP_ALL_OUTPUT_TYPES,
    CPP_DEFAULT_OUTPUT_TYPES,
    CPP_OUTPUT_TYPES,
    DEFAULT_CONFIG,
    DEFAULT_MODE,
    DEFAULT_ORDER_PARAMETERS,
    DEFAULT_OUTPUT_TYPES,
    GRO_OUTPUT_TYPES,
    ORDER_PARAMETER_ALIASES,
    ORDER_PARAMETER_CHOICES,
    OUTPUT_TYPE_ORDER,
)
from .presets import is_cpp_mode


def order_parameter_sort_key(name: str) -> tuple[int, int]:
    fixed = {
        "f3": (0, 0),
        "f4": (1, 0),
        "mcg1": (3, 0),
        "mcg3": (4, 0),
        "dhop35": (5, 0),
        "dhop30": (6, 0),
    }
    if name in fixed:
        return fixed[name]
    match = re.fullmatch(r"q(\d+)", name)
    return (2, int(match.group(1))) if match else (99, 0)


def normalize_order_parameters(value: Any = None) -> tuple[str, ...]:
    if value is None or value == "":
        raw: list[Any] = list(DEFAULT_ORDER_PARAMETERS)
    elif isinstance(value, str):
        raw = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "order_parameter.enabled / --order-parameter must be a comma-separated list."
            ) from exc
    if not raw:
        return ()
    cleaned = [
        ORDER_PARAMETER_ALIASES.get(str(item).strip().lower(), str(item).strip().lower())
        for item in raw
    ]
    keywords = set(cleaned) & {"all", "none"}
    if keywords:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in order_parameter.enabled / --order-parameter.")
        return ALL_ORDER_PARAMETERS if cleaned[0] == "all" else ()
    supported = {"f3", "f4", "mcg1", "mcg3", "dhop35", "dhop30"}
    normalized: set[str] = set()
    for name in cleaned:
        if name in supported:
            normalized.add(name)
        elif (match := re.fullmatch(r"q(\d+)", name)):
            normalized.add(f"q{int(match.group(1))}")
        else:
            raise ValueError(
                f"Unsupported order parameter '{name}'. Use "
                f"{', '.join(ORDER_PARAMETER_CHOICES)}."
            )
    return tuple(sorted(normalized, key=order_parameter_sort_key))


def normalize_cpp_order_parameters(value: Any = None) -> tuple[str, ...]:
    if value is None or value == "":
        raw: list[Any] = list(DEFAULT_ORDER_PARAMETERS)
    elif isinstance(value, str):
        raw = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw = [item for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError(
                "engine cpp --order-parameter must be f3, f4, f3,f4, all, or none."
            ) from exc
    if not raw:
        return ()
    cleaned = [str(item).strip().lower() for item in raw]
    if set(cleaned) & {"all", "none"}:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in engine cpp.")
        return DEFAULT_ORDER_PARAMETERS if cleaned[0] == "all" else ()
    normalized = normalize_order_parameters(cleaned)
    unsupported = [name for name in normalized if name not in {"f3", "f4"}]
    if unsupported:
        raise ValueError(
            f"order parameter(s) {', '.join(unsupported)} are not supported in "
            "engine cpp; use f3 and/or f4."
        )
    return normalized


def q_degrees_from_order_parameters(value: Any) -> tuple[int, ...]:
    return tuple(
        int(name[1:])
        for name in normalize_order_parameters(value)
        if name.startswith("q")
    )


def order_parameter_display(value: Any) -> str:
    values = normalize_order_parameters(value)
    return ", ".join(values) if values else "none"


def _output_items(value: Any, defaults: tuple[str, ...]) -> list[str]:
    if value is None:
        return list(defaults)
    if isinstance(value, str):
        return [item.strip().lower() for item in value.split(",") if item.strip()]
    try:
        return [str(item).strip().lower() for item in value if str(item).strip()]
    except TypeError as exc:
        raise ValueError("output.type / --output-type must be a comma-separated list.") from exc


def normalize_output_types(value: Any = None) -> tuple[str, ...]:
    cleaned = _output_items(value, DEFAULT_OUTPUT_TYPES)
    if not cleaned:
        return ()
    if set(cleaned) & {"all", "none"}:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in output.type / --output-type.")
        return ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()
    if "default" in cleaned:
        cleaned = [*DEFAULT_OUTPUT_TYPES, *(item for item in cleaned if item != "default")]
    unknown = sorted(set(cleaned) - set(OUTPUT_TYPE_ORDER))
    removed = {"sqq-cage-gro", "vmd"}.intersection(unknown)
    if removed:
        raise ValueError(
            f"Output type(s) {', '.join(sorted(removed))} were removed; use sqq-render instead."
        )
    if unknown:
        raise ValueError(f"Unsupported output type(s) {unknown}.")
    normalized = set(cleaned)
    if "gro" in normalized:
        normalized.difference_update(GRO_OUTPUT_TYPES)
    return tuple(name for name in OUTPUT_TYPE_ORDER if name in normalized)


def normalize_cpp_output_types(value: Any = None) -> tuple[str, ...]:
    cleaned = _output_items(value, CPP_DEFAULT_OUTPUT_TYPES)
    if not cleaned:
        return ()
    if set(cleaned) & {"all", "none"}:
        if len(cleaned) != 1:
            raise ValueError("Use 'all' or 'none' alone in engine cpp.")
        return CPP_ALL_OUTPUT_TYPES if cleaned[0] == "all" else ()
    if "default" in cleaned:
        cleaned = [*CPP_DEFAULT_OUTPUT_TYPES, *(item for item in cleaned if item != "default")]
    unsupported = sorted(set(cleaned) - CPP_OUTPUT_TYPES)
    removed = {"sqq-cage-gro", "vmd"}.intersection(unsupported)
    if removed:
        raise ValueError(
            f"Output type(s) {', '.join(sorted(removed))} were removed; use sqq-render instead."
        )
    if unsupported:
        raise ValueError(
            f"output type(s) {', '.join(unsupported)} are not supported by SQQ-CPP."
        )
    return tuple(name for name in CPP_ALL_OUTPUT_TYPES if name in set(cleaned))


def output_type_display(value: Any, *, cpp_mode: bool = False) -> str:
    values = normalize_cpp_output_types(value) if cpp_mode else normalize_output_types(value)
    return ", ".join(values) if values else "none"


def output_enabled(config: dict[str, Any], output_type: str) -> bool:
    name = str(output_type).strip().lower()
    if name not in set(OUTPUT_TYPE_ORDER):
        raise ValueError(f"Unsupported output type: {output_type}")
    enabled = set(normalize_output_types(config.get("output", {}).get("types", DEFAULT_OUTPUT_TYPES)))
    return name in enabled or (name in GRO_OUTPUT_TYPES and "gro" in enabled)


def normalize_parallel_backend(value: Any) -> str:
    """Normalize the configured frame-execution backend."""
    backend = str(value or "process").strip().lower()
    if backend not in {"process", "thread", "serial"}:
        raise ValueError(
            "parallel.backend / --parallel-backend must be process, thread, or serial."
        )
    return backend


_RESOLVED_RUN_KEYS = {
    "sqq_version", "release_date", "engine_selector", "engine", "config_output",
    "status", "error", "graph_mode_requested", "graph_mode_effective",
    "graph_mode_reason", "graph_mode_display", "graph_mode_by_group",
    "graph_mode_reason_by_group", "resolution_adjustments",
    "order_parameters", "find_half", "find_quasi", "find_cluster", "output_types",
    "sampling_interval", "native_frame_interval_ps", "delta_time_ps", "raw_frame_step",
    "selected_frames", "source_frames_total", "frames_total", "frames_ok", "frames_failed",
    "failures", "worker_request", "worker_policy", "workers_resolved", "parallel_backend",
    "math_threads_per_worker", "summary_write", "topology_group_count", "topology_group_limit",
    "topology_group_limit_exceeded", "topology_group_labels_enabled",
    "info_only_fallback_required", "topology_grouping", "topology_groups",
    "topology_source_mapping", "topology_group", "topology_fingerprint",
    "requested_output_types", "output_policy", "warnings", "profile",
}
_CONFIG_EXTRA_KEYS: dict[tuple[str, ...], set[str]] = {
    (): {"adjustments", "resolution_report"},
    ("run",): _RESOLVED_RUN_KEYS,
    ("input",): {"format", "topology", "sampling"},
    ("input", "lammps"): {"resolved_type_map", "type_map_source", "rebuilt_molecules"},
    ("graph",): {
        "effective_bond_mode", "effective_bond_mode_reason",
        "effective_bond_mode_by_group", "effective_bond_mode_reason_by_group",
    },
}
_CONFIG_FLEXIBLE_PATHS = {
    ("component", "role_map"), ("guest", "center_atoms"), ("input", "sampling"),
    ("input", "lammps", "type_map"), ("input", "lammps", "resolved_type_map"),
    ("graph", "effective_bond_mode_by_group"),
    ("graph", "effective_bond_mode_reason_by_group"),
    ("run", "graph_mode_by_group"), ("run", "graph_mode_reason_by_group"),
    ("run", "summary_write"),
}


def validate_user_config_keys(
    config: dict[str, Any],
    schema: dict[str, Any] | None = None,
    path: tuple[str, ...] = (),
) -> None:
    if path in _CONFIG_FLEXIBLE_PATHS:
        return
    expected = DEFAULT_CONFIG if schema is None else schema
    for key, value in config.items():
        if key not in expected:
            if key in _CONFIG_EXTRA_KEYS.get(path, set()):
                continue
            raise ValueError(f"Unknown configuration key: {'.'.join((*path, str(key)))}")
        expected_value = expected[key]
        if isinstance(expected_value, dict):
            if not isinstance(value, dict):
                raise ValueError(
                    f"Configuration section must be a mapping: {'.'.join((*path, str(key)))}"
                )
            validate_user_config_keys(value, expected_value, (*path, str(key)))


def _normalize_cpp_ring_sizes(value: Any) -> list[int]:
    raw = [item.strip() for item in value.split(",") if item.strip()] if isinstance(value, str) else list(value)
    sizes = sorted({int(item) for item in raw})
    if not sizes:
        raise ValueError("engine cpp requires at least one ring size from 4, 5, and 6")
    unsupported = [size for size in sizes if size not in {4, 5, 6}]
    if unsupported:
        raise ValueError(
            f"ring size(s) {', '.join(str(size) for size in unsupported)} are not supported in engine cpp"
        )
    return sizes


def _requested_on(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off", "none"}
    return bool(value)


def validate_cpp_cli(args: Any, config: dict[str, Any]) -> None:
    """Finalize and validate the supported C++ subset.

    Only current public CLI fields are inspected. Removed options are handled by
    the CLI migration guard and are deliberately absent here.
    """
    errors: list[str] = []
    if is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        if getattr(args, "find_half", None) not in (None, False, "off"):
            errors.append("--find-half on is not supported by SQQ-CPP")
        if getattr(args, "find_quasi", None) not in (None, False, "off"):
            errors.append("--find-quasi on is not supported by SQQ-CPP")
        if getattr(args, "find_cluster", None) not in (None, False, "off"):
            errors.append("--find-cluster on is not supported by SQQ-CPP")
    try:
        validate_cpp_config(config)
    except ValueError as exc:
        errors.append(str(exc))
    if errors:
        raise ValueError("; ".join(dict.fromkeys(errors)))


def validate_cpp_config(config: dict[str, Any]) -> None:
    """Validate the normalized C++ capability subset for CLI and API use."""
    if not is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        return
    errors: list[str] = []

    ring = config.setdefault("ring", {})
    if str(ring.get("definition", "chordless")).strip().lower() != "chordless":
        errors.append("ring.definition must be chordless in engine cpp")
    if not bool(ring.get("chordless", True)):
        errors.append("ring.chordless=false is not supported in engine cpp")
    if ring.get("report_sizes", "auto") not in (None, "", "auto"):
        errors.append("public ring reporting is not supported in engine cpp")
    try:
        ring["sizes"] = _normalize_cpp_ring_sizes(ring.get("sizes", (4, 5, 6)))
    except (TypeError, ValueError) as exc:
        errors.append(str(exc))

    half = config.setdefault("half_cage", {})
    if _requested_on(half.get("enabled", False)):
        errors.append("half_cage.enabled is not supported in engine cpp")
    half["enabled"] = False
    quasi = config.setdefault("quasi_cage", {})
    if _requested_on(quasi.get("enabled", False)):
        errors.append("quasi_cage.enabled is not supported in engine cpp")
    quasi["enabled"] = False
    cluster = config.setdefault("hydrate_cluster", {})
    if _requested_on(cluster.get("enabled", False)):
        errors.append("hydrate_cluster.enabled is not supported in engine cpp")
    cluster["enabled"] = False
    config.setdefault("ice", {})["enabled"] = False

    cage = config.setdefault("cage", {})
    if str(cage.get("search_mode", "grow")).strip().lower() != "grow":
        errors.append("cage.search_mode must be grow in engine cpp")
    if str(cage.get("seed_mode", "ring")).strip().lower() != "ring":
        errors.append("cage.seed_mode must be ring in engine cpp")

    order = config.setdefault("order", {})
    source = order.get("parameters", DEFAULT_ORDER_PARAMETERS)
    try:
        order["parameters"] = list(normalize_cpp_order_parameters(source))
    except ValueError as exc:
        errors.append(str(exc))
    for key in ("q_neighbor_mode", "q_cutoff_nm", "q_n_neighbor"):
        default = DEFAULT_CONFIG["order"][key]
        if order.get(key, default) != default:
            errors.append(f"order.{key} is not supported in engine cpp")
    for key, default in DEFAULT_CONFIG["hydrate_order"].items():
        if config.get("hydrate_order", {}).get(key, default) != default:
            errors.append(f"hydrate_order.{key} is not supported in engine cpp")

    output = config.setdefault("output", {})
    source = output.get("types", CPP_DEFAULT_OUTPUT_TYPES)
    try:
        output["types"] = list(normalize_cpp_output_types(source))
    except ValueError as exc:
        errors.append(str(exc))
    if str(config.get("parallel", {}).get("backend", "process")).strip().lower() == "thread":
        errors.append("parallel.backend=thread is not supported in engine cpp")
    if errors:
        raise ValueError("; ".join(dict.fromkeys(errors)))


__all__ = [
    "normalize_cpp_order_parameters", "normalize_cpp_output_types",
    "normalize_order_parameters", "normalize_output_types", "normalize_parallel_backend",
    "order_parameter_display",
    "order_parameter_sort_key", "output_enabled", "output_type_display",
    "q_degrees_from_order_parameters", "validate_cpp_cli", "validate_cpp_config", "validate_user_config_keys",
]
