from __future__ import annotations

"""Python adapter for the optional C++17 analysis core."""

from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from .. import __version__
from ..config import normalize_cpp_order_parameters
from ..models import Cage, F3F4Result, Frame, FrameResult, GraphResult, Guest, Ring, Water, WaterOrder
from .graph import read_pair_edges


def analyze_frame_cpp(
    frame: Frame,
    waters: list[Water],
    guests: list[Guest],
    config: dict[str, Any],
    *,
    cage_report_types: tuple[str, ...] | None = None,
    ring_report_sizes: tuple[int, ...] = (),
) -> FrameResult:
    """Run one normalized frame through ``sqq.core._sqq_cpp``."""
    native = _load_native_module()
    graph_config = config["graph"]
    ring_config = config["ring"]
    cage_config = config["cage"]
    selected_order = set(
        normalize_cpp_order_parameters(config.get("order", {}).get("parameters"))
    )

    pair_edges: list[tuple[int, int]] = []
    if str(graph_config.get("bond_mode", "auto")) == "pairs":
        pair_path = graph_config.get("pair_file")
        if not pair_path:
            raise ValueError("mode cpp with bond_mode=pairs requires --pairs or graph.pair_file.")
        pair_edges = read_pair_edges(
            Path(pair_path),
            frame.atoms,
            waters,
            str(graph_config.get("pair_id", "resid")),
        )

    native_frame = {
        "positions": [tuple(float(value) for value in atom.xyz) for atom in frame.atoms],
        "waters": [
            {
                "oxygen": int(water.oxygen),
                "hydrogens": [int(index) for index in water.hydrogens],
            }
            for water in waters
        ],
        "guests": [
            {
                "resid": int(guest.resid),
                "resname": str(guest.resname),
                "atoms": [int(index) for index in guest.atoms],
                "center_atom": None if guest.center_atom is None else int(guest.center_atom),
            }
            for guest in guests
        ],
        "box": _native_box(frame.box),
        "pair_edges": pair_edges,
    }
    native_options = {
        "bond_mode": str(graph_config.get("bond_mode", "auto")),
        "oo_cutoff_nm": float(graph_config.get("oo_cutoff_nm", 0.35)),
        "hbond_distance_nm": float(graph_config.get("hbond_distance_nm", 0.35)),
        "hbond_angle_deg": float(graph_config.get("hbond_angle_deg", 30.0)),
        "ring_sizes": [int(size) for size in ring_config.get("sizes", (4, 5, 6))],
        "chordless": True,
        "cage_enabled": bool(cage_config.get("enabled", True)),
        "max_faces": int(cage_config.get("max_faces", 20)),
        "cage_target_types": [],
        "cage_report_types": [],
        "max_states_per_seed": int(cage_config.get("max_states_per_seed", 20000)),
        "max_total_states": int(cage_config.get("max_total_states", 5000000)),
        "max_boundary_candidates": int(cage_config.get("max_boundary_candidates", 8)),
        "scientific_validation": bool(cage_config.get("scientific_validation", False)),
        "max_face_planarity_rms_nm": float(cage_config.get("max_face_planarity_rms_nm", 0.06)),
        "max_face_edge_cv": float(cage_config.get("max_face_edge_cv", 0.35)),
        "min_cage_volume_nm3": float(cage_config.get("min_cage_volume_nm3", 1.0e-6)),
        "occupancy_mode": str(cage_config.get("occupancy_mode", "polyhedron")),
        "occupancy_radius_nm": float(cage_config.get("occupancy_radius_nm", 0.5)),
        "compute_f3": "f3" in selected_order,
        "compute_f4": "f4" in selected_order,
    }

    raw = native.analyze(native_frame, native_options)
    if not isinstance(raw, dict):
        raise RuntimeError("SQQ-CPP returned an invalid non-mapping result.")
    graph = _graph_result(waters, raw)
    rings, ring_by_native_index = _ring_results(raw.get("rings", ()))
    all_cages = _cage_results(raw.get("cages", ()), ring_by_native_index, guests)
    cages = _reported_cages(all_cages, cage_report_types)
    f3f4 = _f3f4_result(frame, waters, raw.get("f3f4"), selected_order, config)
    warnings = [str(item) for item in raw.get("warnings", ())]
    if "f4" in selected_order and f3f4 is not None and f3f4.f4_valid == 0:
        message = "F4 is unavailable because usable water-hydrogen coordinates are missing."
        if message not in warnings:
            warnings.append(message)

    return FrameResult(
        frame=frame,
        waters=waters,
        guests=guests,
        graph=graph,
        rings=rings,
        ring_report_sizes=tuple(ring_report_sizes),
        cages=cages,
        all_cages=all_cages,
        cage_report_types=cage_report_types,
        f3f4=f3f4,
        warnings=warnings,
    )


def _load_native_module():
    try:
        from . import _sqq_cpp
    except ImportError as exc:  # pragma: no cover - platform/package dependent.
        raise RuntimeError(
            "SQQ-CPP native extension is unavailable. Install a wheel built for this "
            "Python version and platform; mode cpp never falls back to sqq-py."
        ) from exc
    try:
        native_version = str(_sqq_cpp.core_version())
    except (AttributeError, TypeError) as exc:
        raise RuntimeError("SQQ-CPP native extension does not expose a valid core version.") from exc
    if native_version != __version__:
        raise RuntimeError(
            f"SQQ-CPP version mismatch: Python package {__version__}, native core "
            f"{native_version}. Reinstall SQQ for this Python environment."
        )
    return _sqq_cpp


def _native_box(box: np.ndarray | None) -> tuple[float, float, float] | None:
    if box is None or len(box) < 3:
        return None
    values = tuple(float(value) for value in np.asarray(box, dtype=float)[:3])
    return values if all(value > 0.0 for value in values) else None


def _graph_result(waters: list[Water], raw: dict[str, Any]) -> GraphResult:
    edges = sorted(
        {
            tuple(sorted((int(pair[0]), int(pair[1]))))
            for pair in raw.get("edges", ())
            if int(pair[0]) != int(pair[1])
        }
    )
    adjacency: dict[int, set[int]] = {int(water.oxygen): set() for water in waters}
    for left, right in edges:
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)
    mode = str(raw.get("effective_bond_mode", raw.get("graph_mode", "")))
    return GraphResult(mode=mode, edges=edges, adjacency=adjacency)


def _ring_results(
    records: Any,
) -> tuple[dict[int, list[Ring]], dict[int, Ring]]:
    normalized = sorted(
        [
            (
                int(record.get("index", index)),
                int(record["size"]),
                tuple(int(node) for node in record["nodes"]),
            )
            for index, record in enumerate(records or ())
        ],
        key=lambda item: (item[1], item[2], item[0]),
    )
    counts: dict[int, int] = defaultdict(int)
    by_size: dict[int, list[Ring]] = defaultdict(list)
    by_native_index: dict[int, Ring] = {}
    for native_index, size, nodes in normalized:
        counts[size] += 1
        ring = Ring(object_id=f"ring{size}_{counts[size]:05d}", nodes=nodes)
        by_size[size].append(ring)
        by_native_index[native_index] = ring
    return dict(by_size), by_native_index


def _cage_results(
    records: Any,
    ring_by_native_index: dict[int, Ring],
    guests: list[Guest],
) -> list[Cage]:
    normalized = sorted(
        records or (),
        key=lambda record: (
            str(record.get("cage_type", "")),
            tuple(int(node) for node in record.get("waters", ())),
            tuple(int(index) for index in record.get("ring_indices", ())),
        ),
    )
    type_counts: dict[str, int] = defaultdict(int)
    cages: list[Cage] = []
    for record in normalized:
        cage_type = str(record["cage_type"])
        type_counts[cage_type] += 1
        ring_ids = tuple(
            sorted(
                ring_by_native_index[int(index)].object_id
                for index in record.get("ring_indices", ())
            )
        )
        guest_ids = _guest_ids(record, guests)
        cages.append(
            Cage(
                object_id=f"{cage_type}_{type_counts[cage_type]:05d}",
                cage_type=cage_type,
                rings=ring_ids,
                waters=tuple(sorted(int(node) for node in record.get("waters", ()))),
                center=np.asarray(record.get("center", (0.0, 0.0, 0.0)), dtype=float),
                guest_ids=guest_ids,
                isomer=None if record.get("isomer") is None else str(record["isomer"]),
            )
        )
    return cages


def _guest_ids(record: dict[str, Any], guests: list[Guest]) -> tuple[str, ...]:
    indexes = record.get("guest_indices")
    if indexes is not None:
        return tuple(
            f"{guests[int(index)].resname}{guests[int(index)].resid}"
            for index in indexes
            if 0 <= int(index) < len(guests)
        )
    return tuple(str(item) for item in record.get("guest_ids", ()))


def _reported_cages(
    cages: list[Cage],
    report_types: tuple[str, ...] | None,
) -> list[Cage]:
    if report_types is None:
        return list(cages)
    allowed = set(report_types)
    return [cage for cage in cages if cage.cage_type in allowed]


def _f3f4_result(
    frame: Frame,
    waters: list[Water],
    raw: Any,
    selected_order: set[str],
    config: dict[str, Any],
) -> F3F4Result | None:
    if not selected_order:
        return None
    data = raw if isinstance(raw, dict) else {}
    native_rows = list(data.get("per_water", ()))
    by_index = {int(row.get("water_index", index)): row for index, row in enumerate(native_rows)}
    focus_resids = {int(item) for item in config.get("order", {}).get("focus_waters", ())}
    per_water: list[WaterOrder] = []
    f3_values: list[float] = []
    f4_values: list[float] = []
    f3_focus: list[float] = []
    f4_focus: list[float] = []
    for index, water in enumerate(waters):
        row = by_index.get(index, {})
        f3 = _optional_float(row.get("f3")) if "f3" in selected_order else None
        f4 = _optional_float(row.get("f4")) if "f4" in selected_order else None
        if f3 is not None:
            f3_values.append(f3)
            if water.resid in focus_resids:
                f3_focus.append(f3)
        if f4 is not None:
            f4_values.append(f4)
            if water.resid in focus_resids:
                f4_focus.append(f4)
        atom = frame.atoms[water.oxygen]
        per_water.append(
            WaterOrder(
                oxygen=water.oxygen,
                resid=water.resid,
                atomid=atom.atomid,
                xyz=atom.xyz,
                f3=f3,
                f4=f4,
            )
        )
    return F3F4Result(
        per_water=tuple(per_water),
        f3_mean=_mean_or_none(f3_values),
        f4_mean=_mean_or_none(f4_values),
        f3_valid=len(f3_values),
        f4_valid=len(f4_values),
        focus_resids=tuple(sorted(focus_resids)),
        f3_focus_mean=_mean_or_none(f3_focus) if focus_resids else None,
        f4_focus_mean=_mean_or_none(f4_focus) if focus_resids else None,
        f3_focus_valid=len(f3_focus),
        f4_focus_valid=len(f4_focus),
    )


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)


def _mean_or_none(values: list[float]) -> float | None:
    return None if not values else float(sum(values) / len(values))
