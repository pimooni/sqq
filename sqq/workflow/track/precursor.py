"""Persistent-cage precursor reconstruction for the Track workflow.

Precursor history follows the water set of a persistent cage backwards through
the frames before its birth. The per-frame facts needed for that (water-oxygen
identities, coordinates, water-water edges, ring/half-cage/quasi-cage/cage
water sets, and the render atom mapping) form one JSON-safe *precursor record*.
The first raw-Track pass spools these records while it analyzes, so the prefix
never has to be re-analyzed; when no spool exists the prefix is re-analyzed
through the shared runner and reduced to the same records.
"""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import csv
from dataclasses import replace
import os
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np

from ...config import DEFAULT_MODE, is_cpp_mode
from ...core.common.pbc import minimum_image
from ...io.render import TRACK_MEMBERSHIP_NAME, TRACK_RENDER_DIRECTORY
from ...io.tracking import (
    add_precursor_membership,
    append_track_info_section,
    target_directory_name,
)
from ...io.tracking.precursor_record import (
    PRECURSOR_RECORD_FORMAT,
    PRECURSOR_RECORD_VERSION,
    PrecursorRecord,
    precursor_record_from_result,
    precursor_record_path,
    read_precursor_record,
    write_precursor_record,
)
from ...models.tracking import TargetSelection
from ...runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ...runtime.session import AnalysisEvent, AnalysisRunner, AnalysisSink
from ...ui.progress import RunProgressDisplay


_PRECURSOR_STATE_FIELDS = (
    "status", "reason", "track_id", "frame_index", "frame", "time_ps",
    "state", "target_water_count", "present_water_count", "bond_count",
    "component_count", "largest_component", "ring_count", "half_cage_count",
    "quasi_cage_count", "cage_count",
)
_WATER_HISTORY_FIELDS = (
    "status", "reason", "track_id", "frame_index", "frame", "time_ps",
    "state", "water_atomid", "present", "atom_index", "resid", "x_nm",
    "y_nm", "z_nm", "target_degree",
)
PrecursorData = tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[int, tuple[tuple[int, ...], tuple[float, float, float]]],
]

# --- outputs ------------------------------------------------------------------------


def write_precursor_outputs(
    selections: Sequence[TargetSelection],
    written: Mapping[str, Path],
    histories: Mapping[str, PrecursorData] | None,
    *,
    source_mode: bool,
) -> None:
    for selection in selections:
        if selection.target.kind != "track":
            continue
        directory = written.get(target_directory_name(selection.target))
        if directory is None:
            continue
        target_dir = Path(directory)
        if histories is None:
            reason = (
                "unavailable from imported Track state because per-frame graph, "
                "ring, half-cage, and quasi-cage objects were not saved"
                if source_mode
                else "per-frame analysis data are unavailable"
            )
            precursor = [_unavailable_row(selection.target.value, reason)]
            water = [_unavailable_row(selection.target.value, reason)]
            _append_precursor_status(target_dir / "track_info.md", "unavailable", reason)
        else:
            precursor, water, render_frames = histories.get(
                selection.target.value, ([], [], {})
            )
            status = "available" if precursor else "unavailable"
            reason = (
                "" if precursor else "the selected cage is left-censored at the first frame"
            )
            if not precursor:
                precursor = [_unavailable_row(selection.target.value, reason)]
                water = [_unavailable_row(selection.target.value, reason)]
            _append_precursor_status(target_dir / "track_info.md", status, reason)
            if render_frames:
                add_precursor_membership(
                    target_dir / TRACK_RENDER_DIRECTORY / TRACK_MEMBERSHIP_NAME,
                    selection.target.value,
                    render_frames,
                )
        _write_csv(target_dir / "precursor_state.csv", precursor, _PRECURSOR_STATE_FIELDS)
        _write_csv(target_dir / "water_history.csv", water, _WATER_HISTORY_FIELDS)


# --- histories ----------------------------------------------------------------------


def raw_precursor_histories(
    selections: Sequence[TargetSelection],
    plan: RunPlan,
    config: Mapping[str, Any],
) -> dict[str, PrecursorData]:
    """Classify the pre-birth prefix of every persistent-ID target.

    Spooled records from the first pass are used when they exist for every
    required frame; otherwise the prefix is re-analyzed through the shared
    runner and reduced to the same records, so both paths produce identical
    tables.
    """
    targets: dict[str, tuple[int, frozenset[int]]] = {}
    output: dict[str, PrecursorData] = {}
    for selection in selections:
        if selection.target.kind != "track":
            continue
        if len(selection.tracks) != 1:
            raise ValueError(
                f"Persistent target {selection.target.value} must select exactly one cage."
            )
        track = selection.tracks[0]
        birth = track.first
        if birth.frame_index <= 0:
            output[track.track_id] = ([], [], {})
            continue
        targets[track.track_id] = (
            int(birth.frame_index),
            frozenset(int(value) for value in birth.water_atomids),
        )
        output[track.track_id] = ([], [], {})
    if not targets:
        return output

    last_required = max(item[0] for item in targets.values())
    spool = plan.context.precursor_spool_dir
    records = _spooled_records(spool, last_required) if spool is not None else None
    if records is None:
        _reanalyze_and_classify(
            plan,
            config,
            last_required,
            targets,
            output,
        )
        return output
    for frame_index, record in enumerate(records):
        _classify_into(record, frame_index, targets, output)
    return output


def _spooled_records(
    spool: Path,
    last_required: int,
) -> Iterator[PrecursorRecord] | None:
    """Return a lazy ordered reader when the complete required prefix exists."""
    if any(
        not precursor_record_path(spool, index).is_file()
        for index in range(last_required + 1)
    ):
        return None

    def records() -> Iterator[PrecursorRecord]:
        for frame_index in range(last_required + 1):
            record = read_precursor_record(spool, frame_index)
            if record is None:
                raise RuntimeError(
                    "A validated precursor record disappeared while it was being read: "
                    f"frame {frame_index}."
                )
            yield record

    return records()


def _classify_into(
    record: PrecursorRecord,
    frame_index: int,
    targets: Mapping[str, tuple[int, frozenset[int]]],
    output: dict[str, PrecursorData],
) -> None:
    for track_id, (birth_index, water_ids) in targets.items():
        if frame_index > birth_index:
            continue
        state, water_rows = _classify_precursor_record(record, frame_index, track_id, water_ids)
        states, waters, render_frames = output[track_id]
        states.append(state)
        waters.extend(water_rows)
        if frame_index < birth_index:
            render_frames[frame_index] = _precursor_render_frame(record, water_ids)


class _PrecursorRecordSink(AnalysisSink):
    """Classify each fallback frame immediately without retaining its record."""

    def __init__(
        self,
        atom_scope: str,
        targets: Mapping[str, tuple[int, frozenset[int]]],
        output: dict[str, PrecursorData],
    ) -> None:
        self.atom_scope = atom_scope
        self.targets = targets
        self.output = output
        self.consumed = 0

    def start(self, plan: RunPlan) -> None:
        self.consumed = 0

    def consume(self, task: FrameTask, outcome: TaskOutcome) -> None:
        if not outcome.ok or outcome.result is None:
            raise RuntimeError(
                f"Precursor tracing failed for frame {task.display_name!r}: "
                f"{outcome.error_message or 'analysis result unavailable'}"
            )
        frame_index = int(task.frame_index)
        record = precursor_record_from_result(
            outcome.result,
            frame_index,
            atom_scope=self.atom_scope,
        )
        _classify_into(record, frame_index, self.targets, self.output)
        self.consumed += 1

    def finish(self, plan: RunPlan, outcomes: Sequence[TaskOutcome]) -> None:
        if self.consumed != len(plan.tasks):
            raise RuntimeError(
                "Precursor tracing analyzed "
                f"{self.consumed} frames; expected {len(plan.tasks)}."
            )


class _PrecursorProgressSink:
    """Translate runner events into the precursor progress panel."""

    def __init__(self, display: RunProgressDisplay) -> None:
        self.display = display

    def __call__(self, event: AnalysisEvent) -> None:
        task = event.task
        if event.kind == "task-start" and task is not None:
            self.display.start_frame(task.frame_index, task.display_name)
        elif event.kind == "stage" and event.stage is not None:
            self.display.update_stage(event.stage)
        elif event.kind == "task-complete":
            self.display.complete_frame(event.status == "ok")
        elif event.kind == "task-cancelled":
            self.display.complete_frame(False)


def _reanalyze_and_classify(
    plan: RunPlan,
    config: Mapping[str, Any],
    last_required: int,
    targets: Mapping[str, tuple[int, frozenset[int]]],
    output: dict[str, PrecursorData],
) -> None:
    """Re-run and immediately classify a prefix whose spool is incomplete."""
    prefix_tasks = tuple(
        sorted(
            (task for task in plan.tasks if int(task.frame_index) <= last_required),
            key=lambda task: int(task.frame_index),
        )
    )
    if [int(task.frame_index) for task in prefix_tasks] != list(range(last_required + 1)):
        raise RuntimeError(
            "The first-pass plan does not contain every prefix frame required for "
            "precursor tracing."
        )
    precursor_config = deepcopy(dict(plan.context.config))
    precursor_config.setdefault("output", {})["types"] = []
    precursor_plan = replace(
        plan,
        tasks=prefix_tasks,
        context=replace(
            plan.context,
            config=precursor_config,
            strict=True,
            retain_results=False,
            stream_results=True,
            tracking_snapshots=False,
            precursor_spool_dir=None,
            fragment_dir=None,
            group_fragment_dirs={},
        ),
        policy=replace(plan.policy, strict=True),
    )
    sink = _PrecursorRecordSink(
        str(config.get("render", {}).get("atom_scope", "full")),
        targets,
        output,
    )
    progress = RunProgressDisplay(
        total=len(prefix_tasks),
        total_started_at=perf_counter(),
        include_cluster_stage=bool(
            config.get("hydrate_cluster", {}).get("enabled", False)
        ),
        cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        include_patch_stage=bool(
            config.get("half_cage", {}).get("enabled", False)
            or config.get("quasi_cage", {}).get("enabled", False)
        ),
        title="Precursor Pass",
    )
    try:
        outcomes = AnalysisRunner(
            precursor_plan,
            event_sink=_PrecursorProgressSink(progress),
            sinks=(sink,),
        ).run()
    finally:
        progress.close()
    if len(outcomes) != len(prefix_tasks) or any(not item.ok for item in outcomes):
        raise RuntimeError(
            "Precursor tracing returned "
            f"{len(outcomes)} frames; expected {len(prefix_tasks)}."
        )


# --- record classification -----------------------------------------------------------


def _record_oxygen_positions(record: PrecursorRecord) -> dict[int, int]:
    """Map water identity (atom index + 1) to its position in the record arrays."""
    return {int(index) + 1: position for position, index in enumerate(record["oxygen_index"])}


def _precursor_render_frame(
    record: PrecursorRecord,
    target_water_ids: frozenset[int],
) -> tuple[tuple[int, ...], tuple[float, float, float]]:
    positions = _record_oxygen_positions(record)
    ordered_ids = sorted(target_water_ids)
    if any(water_id not in positions for water_id in ordered_ids):
        raise ValueError(
            "A precursor target water is absent from the trajectory topology."
        )
    render_index = record["render_index"]
    atom_indexes: list[int] = []
    for water_id in ordered_ids:
        value = render_index[positions[water_id]]
        if value is None:
            raise ValueError(
                "A precursor target water is absent from the SQQ render topology."
            )
        atom_indexes.append(int(value))
    coordinates = record["oxygen_xyz"]
    box = None if record.get("box") is None else np.asarray(record["box"], dtype=float)
    anchor = np.asarray(coordinates[positions[ordered_ids[0]]], dtype=float)
    points = [anchor]
    for water_id in ordered_ids[1:]:
        delta = minimum_image(
            np.asarray(coordinates[positions[water_id]], dtype=float) - anchor, box
        )
        points.append(anchor + delta)
    center_nm = np.mean(np.asarray(points, dtype=float), axis=0)
    if box is not None:
        center_nm = np.mod(center_nm, box)
    center_angstrom = tuple(float(value) * 10.0 for value in center_nm)
    return tuple(atom_indexes), (center_angstrom[0], center_angstrom[1], center_angstrom[2])


def _classify_precursor_record(
    record: PrecursorRecord,
    frame_index: int,
    track_id: str,
    target_atomids: frozenset[int],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    positions = _record_oxygen_positions(record)
    oxygen_index = record["oxygen_index"]
    index_by_id = {
        atomid: int(oxygen_index[positions[atomid]])
        for atomid in target_atomids
        if atomid in positions
    }
    target_indexes = frozenset(index_by_id.values())
    edges = {
        (int(left), int(right))
        for left, right in record["edges"]
        if left in target_indexes and right in target_indexes
    }
    components = _connected_components(target_indexes, edges)
    rings = [nodes for nodes in record["rings"] if set(nodes).issubset(target_indexes)]
    half = [waters for waters in record["half_cages"] if set(waters).issubset(target_indexes)]
    quasi = [waters for waters in record["quasi_cages"] if set(waters).issubset(target_indexes)]
    cages = [
        waters
        for waters in record["cages"]
        if frozenset(int(index) + 1 for index in waters) == target_atomids
    ]
    if cages:
        state_name = "cage"
    elif quasi:
        state_name = "quasi"
    elif half:
        state_name = "half"
    elif rings:
        state_name = "ring"
    elif edges:
        state_name = "connected"
    else:
        state_name = "dispersed"
    common = {
        "status": "available", "reason": "", "track_id": track_id,
        "frame_index": frame_index, "frame": record["frame_name"],
        "time_ps": record["time_ps"], "state": state_name,
    }
    state_row: dict[str, object] = {
        **common,
        "target_water_count": len(target_atomids),
        "present_water_count": len(target_indexes),
        "bond_count": len(edges),
        "component_count": len(components),
        "largest_component": max((len(item) for item in components), default=0),
        "ring_count": len(rings),
        "half_cage_count": len(half),
        "quasi_cage_count": len(quasi),
        "cage_count": len(cages),
    }
    degree: dict[int, int] = defaultdict(int)
    for left, right in edges:
        degree[left] += 1
        degree[right] += 1
    water_rows: list[dict[str, object]] = []
    for atomid in sorted(target_atomids):
        position = positions.get(atomid)
        if position is None:
            water_rows.append(
                {
                    **common,
                    "water_atomid": atomid,
                    "present": False,
                    "atom_index": "",
                    "resid": "",
                    "x_nm": "",
                    "y_nm": "",
                    "z_nm": "",
                    "target_degree": 0,
                }
            )
            continue
        index = int(oxygen_index[position])
        xyz = record["oxygen_xyz"][position]
        water_rows.append(
            {
                **common,
                "water_atomid": atomid,
                "present": True,
                "atom_index": index,
                "resid": int(record["oxygen_resid"][position]),
                "x_nm": float(xyz[0]),
                "y_nm": float(xyz[1]),
                "z_nm": float(xyz[2]),
                "target_degree": degree[index],
            }
        )
    return state_row, water_rows


def _connected_components(
    nodes: Iterable[int],
    edges: Iterable[tuple[int, int]],
) -> list[frozenset[int]]:
    adjacency: dict[int, set[int]] = {int(node): set() for node in nodes}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(adjacency)
    output: list[frozenset[int]] = []
    while remaining:
        seed = min(remaining)
        stack = [seed]
        members: set[int] = set()
        while stack:
            node = stack.pop()
            if node in members:
                continue
            members.add(node)
            stack.extend(adjacency[node] - members)
        remaining.difference_update(members)
        output.append(frozenset(members))
    return output


def _unavailable_row(track_id: str, reason: str) -> dict[str, object]:
    return {"status": "unavailable", "reason": reason, "track_id": track_id}


def _append_precursor_status(path: Path, status: str, reason: str) -> None:
    suffix = ["## Precursor History", "", f"- status: {status}"]
    if reason:
        suffix.append(f"- reason: {reason}")
    append_track_info_section(path, "\n".join(suffix))


def _write_csv(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(fields), extrasaction="ignore"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


__all__ = [
    "PRECURSOR_RECORD_FORMAT",
    "PRECURSOR_RECORD_VERSION",
    "PrecursorData",
    "PrecursorRecord",
    "precursor_record_from_result",
    "precursor_record_path",
    "raw_precursor_histories",
    "read_precursor_record",
    "write_precursor_outputs",
    "write_precursor_record",
]
