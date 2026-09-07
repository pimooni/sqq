"""Persistent-cage precursor reconstruction for the Track workflow."""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import csv
from dataclasses import replace
import os
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from ...config import DEFAULT_MODE, is_cpp_mode
from ...core.common.geometry import pbc_aware_centroid
from ...io.render import TRACK_MEMBERSHIP_NAME, TRACK_RENDER_DIRECTORY
from ...io.render.frame import visualization_atoms
from ...io.tracking import (
    add_precursor_membership,
    append_track_info_section,
    target_directory_name,
)
from ...models import FrameResult
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


class _PrecursorSink(AnalysisSink):
    """Classify the target water set in every prefix frame, in frame order."""

    def __init__(
        self,
        targets: Mapping[str, tuple[int, frozenset[int]]],
        output: dict[str, PrecursorData],
        *,
        atom_scope: str,
    ) -> None:
        self.targets = dict(targets)
        self.output = output
        self.atom_scope = atom_scope
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
        analyzed = outcome.result
        for track_id, (birth_index, water_ids) in self.targets.items():
            if frame_index > birth_index:
                continue
            state, water_rows = _classify_precursor_frame(
                analyzed, frame_index, track_id, water_ids
            )
            states, waters, render_frames = self.output[track_id]
            states.append(state)
            waters.extend(water_rows)
            if frame_index < birth_index:
                render_frames[frame_index] = _precursor_render_frame(
                    analyzed, water_ids, atom_scope=self.atom_scope
                )
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


def raw_precursor_histories(
    selections: Sequence[TargetSelection],
    plan: RunPlan,
    config: Mapping[str, Any],
) -> dict[str, PrecursorData]:
    """Reanalyze the required pre-birth prefix without retaining full results."""
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
            fragment_dir=None,
            group_fragment_dirs={},
        ),
        policy=replace(plan.policy, strict=True),
    )
    sink = _PrecursorSink(
        targets,
        output,
        atom_scope=str(config.get("render", {}).get("atom_scope", "full")),
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
    return output


def _precursor_render_frame(
    result: FrameResult,
    target_water_ids: frozenset[int],
    *,
    atom_scope: str,
) -> tuple[tuple[int, ...], tuple[float, float, float]]:
    atom_by_id = {int(atom.index) + 1: atom for atom in result.frame.atoms}
    source_indexes = tuple(
        int(atom_by_id[water_id].index)
        for water_id in sorted(target_water_ids)
        if water_id in atom_by_id
    )
    if len(source_indexes) != len(target_water_ids):
        raise ValueError(
            "A precursor target water is absent from the trajectory topology."
        )
    render_index = {
        int(atom.index): index
        for index, atom in enumerate(
            visualization_atoms(result, atom_scope=atom_scope)
        )
    }
    try:
        atom_indexes = tuple(render_index[index] for index in source_indexes)
    except KeyError as exc:
        raise ValueError(
            "A precursor target water is absent from the SQQ render topology."
        ) from exc
    center_nm = np.asarray(
        pbc_aware_centroid(result.frame, list(source_indexes)), dtype=float
    )
    if result.frame.box is not None:
        box = np.asarray(result.frame.box, dtype=float).reshape(-1)
        if (
            len(box) >= 3
            and np.all(np.isfinite(box[:3]))
            and np.all(box[:3] > 0.0)
        ):
            center_nm = np.mod(center_nm, box[:3])
    center_angstrom = tuple(float(value) * 10.0 for value in center_nm)
    return atom_indexes, (
        center_angstrom[0], center_angstrom[1], center_angstrom[2]
    )


def _classify_precursor_frame(
    result: FrameResult,
    frame_index: int,
    track_id: str,
    target_atomids: frozenset[int],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    atom_by_id = {int(atom.index) + 1: atom for atom in result.frame.atoms}
    index_by_id = {
        atomid: int(atom_by_id[atomid].index)
        for atomid in target_atomids
        if atomid in atom_by_id
    }
    target_indexes = frozenset(index_by_id.values())
    edges = {
        tuple(sorted((int(left), int(right))))
        for left, right in result.graph.edges
        if left in target_indexes and right in target_indexes
    }
    components = _connected_components(target_indexes, edges)
    rings = [
        ring for values in result.rings.values() for ring in values
        if set(ring.nodes).issubset(target_indexes)
    ]
    half = [
        patch for patch in result.half_cages
        if set(patch.waters).issubset(target_indexes)
    ]
    quasi = [
        patch for patch in result.quasi_cages
        if set(patch.waters).issubset(target_indexes)
    ]
    cages = [
        cage for cage in (result.all_cages or result.cages)
        if _cage_atomids(result, cage.waters) == target_atomids
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
        "frame_index": frame_index, "frame": result.frame.name,
        "time_ps": result.frame.time_ps, "state": state_name,
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
    degree = defaultdict(int)
    for left, right in edges:
        degree[left] += 1
        degree[right] += 1
    water_rows: list[dict[str, object]] = []
    for atomid in sorted(target_atomids):
        atom = atom_by_id.get(atomid)
        water_rows.append(
            {
                **common,
                "water_atomid": atomid,
                "present": atom is not None,
                "atom_index": "" if atom is None else int(atom.index),
                "resid": "" if atom is None else int(atom.resid),
                "x_nm": "" if atom is None else float(atom.xyz[0]),
                "y_nm": "" if atom is None else float(atom.xyz[1]),
                "z_nm": "" if atom is None else float(atom.xyz[2]),
                "target_degree": 0 if atom is None else degree[int(atom.index)],
            }
        )
    return state_row, water_rows


def _cage_atomids(result: FrameResult, waters: Iterable[int]) -> frozenset[int]:
    atomid_by_index = {
        int(atom.index): int(atom.index) + 1 for atom in result.frame.atoms
    }
    return frozenset(
        atomid_by_index[int(index)]
        for index in waters
        if int(index) in atomid_by_index
    )


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


__all__ = ["PrecursorData", "raw_precursor_histories", "write_precursor_outputs"]
