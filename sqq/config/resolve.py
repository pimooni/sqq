"""Final configuration and graph-mode resolution helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .capabilities import normalize_engine_capabilities
from .presets import (
    apply_mode_preset,
    engine_display,
    is_cpp_mode,
    mode_display,
    mode_label,
    mode_worker_count,
    mode_worker_fraction,
    normalize_mode,
    profile_name,
)
from .validation import (
    normalize_cpp_order_parameters,
    normalize_cpp_output_types,
    normalize_order_parameters,
    normalize_output_types,
    order_parameter_display,
    output_enabled,
    output_type_display,
    q_degrees_from_order_parameters,
    validate_cpp_cli,
)


@dataclass(frozen=True)
class GraphModeResolution:
    requested: str
    effective: str
    reason: str


def resolve_graph_mode(
    requested: str,
    waters: Iterable[Any],
    pair_file: str | Path | None = None,
) -> GraphModeResolution:
    """Resolve graph mode without silently accepting incomplete topology."""
    mode = str(requested).strip().lower()
    if mode not in {"auto", "hbond", "oo", "pairs"}:
        raise ValueError(f"Unsupported bond_mode: {requested}")
    if mode == "pairs":
        if pair_file is None:
            raise ValueError("bond_mode=pairs requires graph.pair_file or --pair.")
        return GraphModeResolution(mode, mode, "explicit pair map")
    if mode != "auto":
        return GraphModeResolution(mode, mode, "explicit graph mode")

    selected = tuple(waters)
    if not selected:
        raise ValueError("graph.mode=auto cannot be resolved because no water molecules were selected.")
    counts = tuple(len(getattr(water, "hydrogens", ())) for water in selected)
    if all(count >= 2 for count in counts):
        return GraphModeResolution("auto", "hbond", "complete water hydrogen topology")
    if all(count == 0 for count in counts):
        return GraphModeResolution("auto", "oo", "oxygen-only water topology")
    raise ValueError(
        "graph.mode=auto found mixed or incomplete water hydrogen topology; "
        "use a consistent topology or choose graph.mode explicitly."
    )


__all__ = [
    "GraphModeResolution", "apply_mode_preset", "engine_display", "is_cpp_mode",
    "mode_display", "mode_label", "mode_worker_count", "mode_worker_fraction",
    "normalize_cpp_order_parameters", "normalize_cpp_output_types",
    "normalize_engine_capabilities", "normalize_mode", "normalize_order_parameters",
    "normalize_output_types", "order_parameter_display", "output_enabled",
    "output_type_display", "profile_name", "q_degrees_from_order_parameters",
    "resolve_graph_mode", "validate_cpp_cli",
]
