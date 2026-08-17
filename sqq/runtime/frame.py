"""Worker-safe analysis dispatch for one already loaded frame."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from ..config import (
    is_cpp_mode,
    resolve_cage_report_types,
    resolve_graph_mode,
    resolve_size_list,
)
from ..core.selection import select_guests, select_waters
from ..core.sqq_cpp import analyze_frame_cpp
from ..core.sqq_py import analyze_frame as analyze_frame_py
from ..io.pairs import read_pair_edges
from ..models import Frame, FrameResult, Guest, Water


StageCallback = Callable[[str], None]


def analyze_frame(
    frame: Frame,
    resolved_config: Mapping[str, Any],
    stage_callback: StageCallback | None = None,
) -> FrameResult:
    """Analyze one frame without output or reporting side effects."""
    config = dict(resolved_config)
    config["graph"] = dict(resolved_config["graph"])
    _report_stage(stage_callback, "resolving settings")
    ring_sizes = resolve_size_list(config["ring"]["sizes"], fallback=[], key="ring.sizes")
    ring_report_sizes = resolve_size_list(
        config["ring"].get("report_sizes", "auto"),
        fallback=ring_sizes,
        key="ring.report_sizes",
    )
    quasi_base_sizes = resolve_size_list(
        config["quasi_cage"].get("base_sizes", "auto"),
        fallback=ring_sizes,
        key="quasi_cage.base_sizes",
    )
    quasi_side_sizes = resolve_size_list(
        config["quasi_cage"].get("side_sizes", "auto"),
        fallback=ring_sizes,
        key="quasi_cage.side_sizes",
    )
    cage_report_types = resolve_cage_report_types(
        config["cage"].get("report_types", []),
        ring_sizes,
        int(config["cage"].get("max_faces", 20)),
    )

    _report_stage(stage_callback, "selecting molecules")
    waters = select_waters(
        frame.atoms,
        resnames=set(config["water"]["resnames"]),
        oxygen_names=set(config["water"]["oxygen_names"]),
        hydrogen_names=set(config["water"]["hydrogen_names"]),
    )
    guests = select_guests(
        frame.atoms,
        resnames=set(config["guest"]["resnames"]),
        center_atoms=config["guest"].get("center_atoms", {}),
        center_mode=str(config["guest"].get("center_mode", "center_atom")),
    )
    graph = config["graph"]
    graph_resolution = resolve_graph_mode(
        str(graph.get("bond_mode", "auto")),
        waters,
        graph.get("pair_file"),
    )
    graph["effective_bond_mode"] = graph_resolution.effective
    graph["effective_bond_mode_reason"] = graph_resolution.reason
    pair_edges = _read_pair_edges(frame, waters, guests, config)

    if is_cpp_mode(config.get("mode")):
        _report_stage(stage_callback, "building water graph")
        _report_stage(stage_callback, "searching rings")
        _report_stage(stage_callback, "searching cage")
        result = analyze_frame_cpp(
            frame,
            waters,
            guests,
            config,
            pair_edges=pair_edges,
            cage_report_types=cage_report_types,
            ring_report_sizes=tuple(ring_report_sizes),
        )
        _report_stage(stage_callback, "computing order parameters")
        return result

    return analyze_frame_py(
        frame,
        waters,
        guests,
        config,
        pair_edges=pair_edges,
        cage_report_types=cage_report_types,
        ring_report_sizes=tuple(ring_report_sizes),
        ring_sizes=ring_sizes,
        quasi_base_sizes=quasi_base_sizes,
        quasi_side_sizes=quasi_side_sizes,
        stage_callback=stage_callback,
    )


def _read_pair_edges(
    frame: Frame,
    waters: list[Water],
    guests: list[Guest],
    config: Mapping[str, Any],
) -> list[tuple[int, int]] | None:
    """Normalize an optional pair file before dispatching either engine."""
    graph = config["graph"]
    graph_mode = str(
        graph.get("effective_bond_mode", graph.get("bond_mode", "auto"))
    ).strip().lower()
    if graph_mode != "pairs":
        return None
    pair_path = graph.get("pair_file")
    if not pair_path:
        raise ValueError(
            "bond_mode=pairs requires --pair or graph.pair_file in sqq_config.yaml."
        )
    return read_pair_edges(
        pair_path,
        frame.atoms,
        waters,
        str(graph.get("pair_id", "resid")),
        guests,
    )


def _report_stage(callback: StageCallback | None, stage: str) -> None:
    if callback is not None:
        callback(stage)


__all__ = [
    "StageCallback",
    "analyze_frame",
]
