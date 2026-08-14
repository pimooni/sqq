"""Backward-compatible YAML migration."""

from __future__ import annotations

from typing import Any
import warnings

from .validation import normalize_order_parameters


YAML_KEY_ALIASES: dict[tuple[str, ...], dict[str, str]] = {
    ("input", "lammps"): {"unit": "units"},
    ("graph",): {"mode": "bond_mode"},
    ("water",): {
        "resname": "resnames", "oxygen_name": "oxygen_names",
        "hydrogen_name": "hydrogen_names",
    },
    ("guest",): {"resname": "resnames", "center_atom": "center_atoms"},
    ("additive",): {"resname": "resnames"},
    ("environment",): {"resname": "resnames"},
    ("ring",): {"size": "sizes", "report_size": "report_sizes"},
    ("quasi_cage",): {
        "base_size": "base_sizes", "side_size": "side_sizes",
        "max_combination_per_base": "max_combinations_per_base",
        "max_layer": "max_layers", "max_ring_per_layer": "max_rings_per_layer",
        "max_layer_state_per_seed": "max_layer_states_per_seed",
        "max_candidate_per_edge": "max_candidates_per_edge",
        "max_layer_candidate": "max_layer_candidates",
    },
    ("cage",): {
        "report_type": "report_types", "max_face": "max_faces",
        "max_state_per_seed": "max_states_per_seed",
        "max_total_state": "max_total_states",
        "max_boundary_candidate": "max_boundary_candidates",
    },
    ("hydrate_order",): {
        "mcg_guest_resname": "mcg_guest_resnames", "mcg_min_water": "mcg_min_waters",
        "dhop_planar_count": "dhop_planar_counts",
        "dhop_min_qualified_neighbor": "dhop_min_qualified_neighbors",
    },
    ("order",): {"enabled": "parameters", "focus_water": "focus_waters"},
    ("ice",): {
        "min_six_ring": "min_six_rings",
        "require_four_coord_neighbor": "require_four_coord_neighbors",
    },
    ("output",): {
        "type": "types", "cage_isomer_row": "cage_isomer_rows",
        "write_empty_file": "write_empty_files", "context_role": "context_roles",
    },
    ("parallel",): {"worker": "workers", "math_thread": "math_threads"},
    ("track",): {"min_shared_water": "min_shared_waters"},
    ("debug",): {"use_networkx_check": "use_networkx_checks"},
}


def nested_mapping(config: dict[str, Any], path: tuple[str, ...]) -> dict[str, Any] | None:
    section: Any = config
    for key in path:
        if not isinstance(section, dict):
            return None
        section = section.get(key)
    return section if isinstance(section, dict) else None


def legacy_enabled(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"false", "off", "no", "0", "none", ""}:
        return False
    if text in {"true", "on", "yes", "1"}:
        return True
    raise ValueError(f"Legacy enable/disable value must be boolean-like; got {value!r}.")


def strip_legacy_selection_keys(user_config: dict[str, Any]) -> None:
    order = user_config.get("order", {})
    if isinstance(order, dict):
        for key in ("f3f4_enabled", "q_enabled", "q_degree"):
            order.pop(key, None)
    hydrate = user_config.get("hydrate_order", {})
    if isinstance(hydrate, dict):
        for key in ("mcg1_enabled", "mcg3_enabled", "dhop35_enabled", "dhop30_enabled"):
            hydrate.pop(key, None)


def migrate_legacy_order_parameters(user_config: dict[str, Any]) -> tuple[str, ...] | None:
    order = user_config.get("order", {})
    hydrate = user_config.get("hydrate_order", {})
    if not isinstance(order, dict) or not isinstance(hydrate, dict) or "parameters" in order:
        return None
    order_keys = {"f3f4_enabled", "q_enabled", "q_degree"}
    hydrate_keys = {"mcg1_enabled", "mcg3_enabled", "dhop35_enabled", "dhop30_enabled"}
    if not (order_keys & set(order) or hydrate_keys & set(hydrate)):
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
        ("mcg1_enabled", "mcg1", True), ("mcg3_enabled", "mcg3", False),
        ("dhop35_enabled", "dhop35", True), ("dhop30_enabled", "dhop30", False),
    ):
        if legacy_enabled(hydrate.get(key, default)):
            selected.append(name)
    return normalize_order_parameters(selected or ["none"])


def _remove_fast_closure(config: dict[str, Any]) -> None:
    cage = config.get("cage")
    if not isinstance(cage, dict):
        return
    removed = False
    for key in ("fast_closure", "fast_closure_max_state", "fast_closure_max_states"):
        removed = cage.pop(key, None) is not None or removed
    if removed:
        message = (
            "Legacy cage.fast_closure settings ignored; exact cage search used."
        )
        adjustments = config.setdefault("adjustments", [])
        if not isinstance(adjustments, list):
            adjustments = list(adjustments) if adjustments else []
            config["adjustments"] = adjustments
        if message not in adjustments:
            adjustments.append(message)
        warnings.warn(
            message,
            UserWarning,
            stacklevel=3,
        )


def migrate_yaml_keys(config: dict[str, Any]) -> dict[str, Any]:
    """Convert canonical and supported legacy YAML keys to runtime keys."""
    config.pop("schema_version", None)
    if "engine" in config:
        if "mode" in config:
            raise ValueError("Config must not contain both top-level 'engine' and legacy 'mode'.")
        config["mode"] = config.pop("engine")
    elif "mode" in config:
        warnings.warn("Top-level config key 'mode' is deprecated; use 'engine'.", UserWarning, stacklevel=2)
    if "order_parameter" in config:
        if "order" in config:
            raise ValueError("Config must not contain both 'order_parameter' and legacy 'order'.")
        config["order"] = config.pop("order_parameter")
    elif "order" in config:
        warnings.warn("Top-level config key 'order' is deprecated; use 'order_parameter'.", UserWarning, stacklevel=2)
    order = config.get("order")
    if isinstance(order, dict) and "parameter" in order:
        if "enabled" in order:
            raise ValueError("Use order_parameter.enabled only; legacy order.parameter is duplicated.")
        warnings.warn("Config key order.parameter is deprecated; use order_parameter.enabled.", UserWarning, stacklevel=2)
        order["enabled"] = order.pop("parameter")
    graph = config.get("graph")
    if isinstance(graph, dict) and "bond_mode" in graph and "mode" not in graph:
        warnings.warn("Config key graph.bond_mode is deprecated; use graph.mode.", UserWarning, stacklevel=2)
    if isinstance(graph, dict) and "pairs_file" in graph:
        if "pair_file" in graph:
            raise ValueError("Use graph.pair_file only; legacy graph.pairs_file is duplicated.")
        graph["pair_file"] = graph.pop("pairs_file")
    _remove_fast_closure(config)
    for path, aliases in YAML_KEY_ALIASES.items():
        section = nested_mapping(config, path)
        if section is None:
            continue
        for canonical, internal in aliases.items():
            if canonical not in section:
                continue
            if internal in section and internal != canonical:
                raise ValueError(
                    f"Config must not contain both {'.'.join((*path, canonical))} and "
                    f"legacy {'.'.join((*path, internal))}."
                )
            section[internal] = section.pop(canonical)
    return config


__all__ = [
    "legacy_enabled", "migrate_legacy_order_parameters", "migrate_yaml_keys",
    "strip_legacy_selection_keys",
]
