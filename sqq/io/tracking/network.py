"""Additive cage-type transition tables and an optional compact plot."""

from __future__ import annotations

from collections import Counter, defaultdict
from math import cos, pi, sin
import os
from pathlib import Path
from typing import Sequence
from uuid import uuid4
import warnings

from ...models.tracking import CageObservation, Row, TargetSelection, TrackingResult

__all__ = [
    "cage_transition_edge_rows",
    "cage_transition_node_rows",
    "write_cage_transition_plot",
]


def cage_transition_edge_rows(
    data: TrackingResult | TargetSelection,
) -> list[Row]:
    """Aggregate directed type changes over directly adjacent observations."""
    counts: Counter[tuple[str, str]] = Counter()
    tracks_by_edge: dict[tuple[str, str], set[str]] = defaultdict(set)
    exposure: Counter[str] = Counter()
    for track in data.tracks:
        ordered = sorted(track.observations, key=_observation_key)
        for left, right in zip(ordered, ordered[1:]):
            if right.frame_index != left.frame_index + 1:
                continue
            if left.time_ps is not None and right.time_ps is not None:
                exposure[left.cage_type] += max(0.0, right.time_ps - left.time_ps)
            if left.cage_type == right.cage_type:
                continue
            edge = (left.cage_type, right.cage_type)
            counts[edge] += 1
            tracks_by_edge[edge].add(track.track_id)
    outgoing = Counter()
    for (source, _destination), count in counts.items():
        outgoing[source] += count
    target = _target_label(data)
    return [
        {
            "target": target,
            "source_cage_type": source,
            "destination_cage_type": destination,
            "transition_count": count,
            "independent_track_count": len(tracks_by_edge[(source, destination)]),
            "source_state_exposure_ps": exposure.get(source, 0.0),
            "conditional_probability": count / outgoing[source],
            "rate_per_ns": (
                None
                if exposure.get(source, 0.0) <= 0.0
                else count / (exposure[source] / 1000.0)
            ),
        }
        for (source, destination), count in sorted(counts.items())
    ]


def cage_transition_node_rows(
    data: TrackingResult | TargetSelection,
    *,
    edges: Sequence[Row] | None = None,
) -> list[Row]:
    """Summarize cage-type exposure and incoming/outgoing transition degree.

    ``edges`` accepts the already computed ``cage_transition_edge_rows(data)``.
    """
    observations: dict[str, list[CageObservation]] = defaultdict(list)
    tracks_by_type: dict[str, set[str]] = defaultdict(set)
    occupancy: dict[str, Counter[str]] = defaultdict(Counter)
    phases: dict[str, Counter[str]] = defaultdict(Counter)
    exposure: Counter[str] = Counter()
    for track in data.tracks:
        ordered = sorted(track.observations, key=_observation_key)
        for item in ordered:
            observations[item.cage_type].append(item)
            tracks_by_type[item.cage_type].add(track.track_id)
            occupancy[item.cage_type][_occupancy_class(item.guest_ids)] += 1
            labels = item.phase_labels or ("unassigned",)
            for phase in labels:
                phases[item.cage_type][phase] += 1
        for left, right in zip(ordered, ordered[1:]):
            if (
                right.frame_index == left.frame_index + 1
                and left.time_ps is not None
                and right.time_ps is not None
            ):
                exposure[left.cage_type] += max(0.0, right.time_ps - left.time_ps)

    edges = cage_transition_edge_rows(data) if edges is None else list(edges)
    incoming_neighbors: dict[str, set[str]] = defaultdict(set)
    outgoing_neighbors: dict[str, set[str]] = defaultdict(set)
    incoming_count: Counter[str] = Counter()
    outgoing_count: Counter[str] = Counter()
    for edge in edges:
        source = str(edge["source_cage_type"])
        destination = str(edge["destination_cage_type"])
        count = int(edge["transition_count"])
        outgoing_neighbors[source].add(destination)
        incoming_neighbors[destination].add(source)
        outgoing_count[source] += count
        incoming_count[destination] += count
    cage_types = sorted(
        set(observations) | set(incoming_neighbors) | set(outgoing_neighbors)
    )
    target = _target_label(data)
    return [
        {
            "target": target,
            "cage_type": cage_type,
            "track_count": len(tracks_by_type[cage_type]),
            "observation_count": len(observations[cage_type]),
            "residence_time_ps": exposure.get(cage_type, 0.0),
            "occupancy_composition": _counter_text(occupancy[cage_type]),
            "phase_composition": _counter_text(phases[cage_type]),
            "incoming_degree": len(incoming_neighbors[cage_type]),
            "outgoing_degree": len(outgoing_neighbors[cage_type]),
            "incoming_transition_count": incoming_count[cage_type],
            "outgoing_transition_count": outgoing_count[cage_type],
        }
        for cage_type in cage_types
    ]


def write_cage_transition_plot(
    data: TrackingResult | TargetSelection,
    path: str | Path,
    *,
    edges: Sequence[Row] | None = None,
    nodes: Sequence[Row] | None = None,
) -> Path | None:
    """Write a small directed network when matplotlib is available."""
    edges = cage_transition_edge_rows(data) if edges is None else list(edges)
    nodes = (
        cage_transition_node_rows(data, edges=edges) if nodes is None else list(nodes)
    )
    if not edges or not nodes:
        return None
    try:
        # The object-oriented figure renders through the Agg canvas and never
        # touches an interactive GUI backend, which may be unavailable on a
        # headless machine even when Matplotlib itself is installed.
        from matplotlib.figure import Figure
    except (ImportError, RuntimeError):
        return None
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.stem}.{uuid4().hex}.tmp.png")
    labels = [str(row["cage_type"]) for row in nodes]
    count = len(labels)
    positions = {
        label: (
            cos(2.0 * pi * index / count),
            sin(2.0 * pi * index / count),
        )
        for index, label in enumerate(labels)
    }
    try:
        figure = Figure(figsize=(max(5.0, count * 0.8), 5.0))
        axes = figure.add_subplot(111)
        axes.set_axis_off()
        for label, (x, y) in positions.items():
            axes.scatter([x], [y], s=850, color="#dbeafe", edgecolor="#1d4ed8", zorder=3)
            axes.text(x, y, label, ha="center", va="center", fontsize=8, zorder=4)
        maximum = max(int(row["transition_count"]) for row in edges)
        for row in edges:
            source = str(row["source_cage_type"])
            destination = str(row["destination_cage_type"])
            x1, y1 = positions[source]
            x2, y2 = positions[destination]
            width = 0.8 + 2.2 * int(row["transition_count"]) / maximum
            axes.annotate(
                "",
                xy=(x2, y2),
                xytext=(x1, y1),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "#475569",
                    "linewidth": width,
                    "shrinkA": 24,
                    "shrinkB": 24,
                    "connectionstyle": "arc3,rad=0.12",
                },
                zorder=2,
            )
        axes.set_xlim(-1.35, 1.35)
        axes.set_ylim(-1.35, 1.35)
        figure.tight_layout()
        figure.savefig(temporary, dpi=180, bbox_inches="tight")
        os.replace(temporary, target)
    except Exception as exc:
        # The CSV representation is authoritative; the plot is best effort,
        # but a broken plotting backend must not fail silently.
        warnings.warn(
            f"SQQ could not write {target.name} ({type(exc).__name__}: {exc}); "
            "the cage-transition CSV tables remain complete.",
            UserWarning,
            stacklevel=2,
        )
        return None
    finally:
        temporary.unlink(missing_ok=True)
    return target


def _observation_key(item: CageObservation) -> tuple[int, str]:
    return item.frame_index, item.local_cage_id


def _occupancy_class(guest_ids: Sequence[str]) -> str:
    count = len(guest_ids)
    return "empty" if count == 0 else "single" if count == 1 else "multiple"


def _counter_text(values: Counter[str]) -> str:
    return ";".join(f"{key}:{values[key]}" for key in sorted(values))


def _target_label(data: TrackingResult | TargetSelection) -> str:
    return "all" if isinstance(data, TrackingResult) else data.target.value
