from __future__ import annotations

"""Build one typed execution plan for every Analyze input form."""

from collections import Counter
from collections.abc import Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

from ..config import DEFAULT_MODE
from ..config.resolve import resolve_graph_mode
from ..core.selection import select_waters
from ..io.gro_grouping import GroGroupingResult, scan_and_group_gro_inputs
from ..io.lammps import LAMMPS_TRAJECTORY_SUFFIXES
from ..io.trajectory import read_frames, read_gro, trajectory_frame_selection
from ..runtime.contracts import (
    ExecutionPolicy,
    FrameTask,
    InputKind,
    RunContext,
    RunPlan,
)
from ..runtime.parallel.policy import (
    normalize_parallel_backend,
    process_in_flight_limit,
    resolve_workers,
)


_TRAJECTORY_SUFFIXES = {".xtc", ".trr"}


def build_run_plan(
    paths: Sequence[str | Path],
    config: Mapping[str, Any],
    output_root: str | Path,
    *,
    topology: str | Path | None = None,
    grouping: GroGroupingResult | None = None,
) -> RunPlan:
    """Resolve inputs, topology groups, frames, graph modes, and execution policy.

    The builder is intentionally side-effect free apart from reading input metadata;
    it never creates output directories or fragment workspaces.
    """
    sources = tuple(Path(item) for item in paths)
    if not sources:
        raise ValueError("Analyze requires at least one input source.")
    root = Path(output_root)
    top = Path(topology) if topology is not None else None
    resolved = deepcopy(dict(config))
    strict = bool(resolved.get("run", {}).get("strict", False))

    if _is_independent_gro_batch(sources):
        return _build_gro_batch_plan(
            sources,
            resolved,
            root,
            top,
            grouping=grouping,
            strict=strict,
        )
    if len(sources) == 1:
        return _build_single_source_plan(sources[0], resolved, root, top, strict)
    return _build_independent_file_plan(sources, resolved, root, top, strict)


def _build_gro_batch_plan(
    sources: tuple[Path, ...],
    config: dict[str, Any],
    output_root: Path,
    topology: Path | None,
    *,
    grouping: GroGroupingResult | None,
    strict: bool,
) -> RunPlan:
    groups = grouping or scan_and_group_gro_inputs(sources, strict=strict)
    execution_config = deepcopy(config)
    if groups.info_only_fallback_required:
        execution_config.setdefault("output", {})["types"] = ["info"]

    global_names = _unique_output_names(groups.assignments)
    group_configs: dict[int, dict[str, Any]] = {}
    group_roots: dict[int, Path] = {}
    task_indexes: dict[int, tuple[int, ...]] = {}
    requested_modes: dict[int, str] = {}
    effective_modes: dict[int, str] = {}
    tasks: list[FrameTask] = []

    for group in groups.groups:
        group_root = (
            output_root
            if groups.info_only_fallback_required or groups.group_count == 1
            else output_root / f"result_{group.label}"
        )
        group_config = deepcopy(execution_config)
        requested, effective = _resolve_graph_for_source(
            group.paths[0], None, group_config, raw_frame_index=None
        )
        group_configs[group.group_index] = group_config
        group_roots[group.group_index] = group_root
        requested_modes[group.group_index] = requested
        effective_modes[group.group_index] = effective

        names = (
            global_names
            if groups.info_only_fallback_required
            else _unique_output_names(group.inputs)
        )
        indexes: list[int] = []
        for local_index, assignment in enumerate(group.inputs):
            indexes.append(assignment.source_index)
            tasks.append(
                FrameTask(
                    task_index=assignment.source_index,
                    frame_index=local_index,
                    source=assignment.path,
                    group_key=group.group_index,
                    output_name=names[assignment.source_index],
                    output_root=group_root,
                    separated_output=True,
                )
            )
        task_indexes[group.group_index] = tuple(indexes)

    tasks.sort(key=lambda item: item.task_index)
    policy = _execution_policy(execution_config, len(tasks), strict)
    context = RunContext(
        config=execution_config,
        output_root=output_root,
        strict=strict,
        topology=topology,
        input_kind=InputKind.GRO_BATCH,
        group_configs=group_configs,
        group_output_roots=group_roots,
    )
    roots = tuple(dict.fromkeys(group_roots.values())) or (output_root,)
    return RunPlan(
        input_kind=InputKind.GRO_BATCH,
        tasks=tuple(tasks),
        context=context,
        policy=policy,
        topology_groups=task_indexes,
        sampling={
            "selected_frames": len(tasks),
            "total_frames": len(sources),
            "failed_sources": len(groups.failures),
            **groups.limit_metadata(),
        },
        requested_graph_modes=requested_modes,
        effective_graph_modes=effective_modes,
        output_roots=roots,
    )


def _build_single_source_plan(
    source: Path,
    config: dict[str, Any],
    output_root: Path,
    topology: Path | None,
    strict: bool,
) -> RunPlan:
    suffix = source.suffix.lower()
    selection = None
    if suffix == ".gro" or suffix in _TRAJECTORY_SUFFIXES or suffix in LAMMPS_TRAJECTORY_SUFFIXES:
        selection = trajectory_frame_selection(
            source,
            topology,
            delta_time_ps=config.get("input", {}).get("delta_time_ps"),
            lammps_config=config.get("input", {}).get("lammps", {}),
        )
    elif config.get("input", {}).get("delta_time_ps") is not None:
        raise ValueError(
            "--delta-time requires one XTC, TRR, LAMMPS trajectory, or stacked GRO input."
        )

    if suffix in LAMMPS_TRAJECTORY_SUFFIXES:
        input_kind = InputKind.LAMMPS
    elif suffix in _TRAJECTORY_SUFFIXES:
        input_kind = InputKind.TRAJECTORY
    elif suffix == ".gro" and selection is not None and selection.total_frames > 1:
        input_kind = InputKind.STACKED_GRO
    else:
        input_kind = InputKind.GRO

    raw_indexes = selection.raw_indexes if selection is not None else (0,)
    requested, effective = _resolve_graph_for_source(
        source,
        topology,
        config,
        raw_frame_index=raw_indexes[0] if raw_indexes else None,
    )
    separated = input_kind is not InputKind.GRO or len(raw_indexes) > 1
    tasks = tuple(
        FrameTask(
            task_index=task_index,
            frame_index=task_index,
            source=source,
            raw_frame_index=(raw_index if selection is not None else None),
            output_root=output_root,
            separated_output=separated,
        )
        for task_index, raw_index in enumerate(raw_indexes)
    )
    sampling = _selection_metadata(selection, len(tasks))
    if sampling:
        config.setdefault("input", {})["sampling"] = dict(sampling)
    else:
        config.setdefault("input", {}).pop("sampling", None)
    policy = _execution_policy(config, len(tasks), strict)
    context = RunContext(
        config=config,
        output_root=output_root,
        strict=strict,
        topology=topology,
        # Stacked GRO files use the native multi-frame reader and need no
        # separate topology.  Only MDAnalysis/LAMMPS trajectories are opened
        # once as worker-owned trajectory readers.
        trajectory=(
            source
            if input_kind in {InputKind.TRAJECTORY, InputKind.LAMMPS}
            else None
        ),
        input_kind=input_kind,
    )
    return RunPlan(
        input_kind=input_kind,
        tasks=tasks,
        context=context,
        policy=policy,
        sampling=sampling,
        requested_graph_modes={"run": requested},
        effective_graph_modes={"run": effective},
        output_roots=(output_root,),
    )


def _build_independent_file_plan(
    sources: tuple[Path, ...],
    config: dict[str, Any],
    output_root: Path,
    topology: Path | None,
    strict: bool,
) -> RunPlan:
    if config.get("input", {}).get("delta_time_ps") is not None:
        raise ValueError(
            "--delta-time applies to one trajectory or one stacked GRO, not multiple files."
        )
    requested, effective = _resolve_graph_for_source(
        sources[0], topology, config, raw_frame_index=None
    )
    names = _unique_path_names(sources)
    tasks = tuple(
        FrameTask(
            task_index=index,
            frame_index=index,
            source=source,
            output_name=names[index],
            output_root=output_root,
            separated_output=True,
        )
        for index, source in enumerate(sources)
    )
    policy = _execution_policy(config, len(tasks), strict)
    return RunPlan(
        input_kind=InputKind.GRO_BATCH,
        tasks=tasks,
        context=RunContext(
            config=config,
            output_root=output_root,
            strict=strict,
            topology=topology,
            input_kind=InputKind.GRO_BATCH,
        ),
        policy=policy,
        sampling={"selected_frames": len(tasks), "total_frames": len(tasks)},
        requested_graph_modes={"run": requested},
        effective_graph_modes={"run": effective},
        output_roots=(output_root,),
    )


def _resolve_graph_for_source(
    source: Path,
    topology: Path | None,
    config: dict[str, Any],
    *,
    raw_frame_index: int | None,
) -> tuple[str, str]:
    graph = config.setdefault("graph", {})
    requested = str(graph.get("bond_mode", "auto"))
    pair_file = graph.get("pair_file")
    if requested == "auto":
        if source.suffix.lower() == ".gro" and raw_frame_index in (None, 0):
            try:
                frame = read_gro(source)
            except ValueError:
                frame = next(
                    iter(
                        read_frames(
                            [source],
                            topology=topology,
                            frame_indexes=[0],
                            lammps_config=config.get("input", {}).get("lammps", {}),
                        )
                    )
                )
        else:
            indexes = None if raw_frame_index is None else [raw_frame_index]
            frame = next(
                iter(
                    read_frames(
                        [source],
                        topology=topology,
                        xyz_scale=float(config.get("input", {}).get("xyz_scale", 0.1)),
                        frame_indexes=indexes,
                        lammps_config=config.get("input", {}).get("lammps", {}),
                    )
                )
            )
        water = config.get("water", {})
        waters = select_waters(
            frame.atoms,
            resnames=set(water.get("resnames", ())),
            oxygen_names=set(water.get("oxygen_names", ())),
            hydrogen_names=set(water.get("hydrogen_names", ())),
        )
        resolution = resolve_graph_mode("auto", waters, pair_file)
    else:
        resolution = resolve_graph_mode(requested, [], pair_file)
    graph["effective_bond_mode"] = resolution.effective
    graph["effective_bond_mode_reason"] = resolution.reason
    return requested, resolution.effective


def _execution_policy(
    config: Mapping[str, Any],
    task_count: int,
    strict: bool,
) -> ExecutionPolicy:
    parallel = config.get("parallel", {})
    backend = normalize_parallel_backend(parallel.get("backend", "process"))
    workers = resolve_workers(
        parallel.get("workers"),
        max(1, task_count),
        mode=config.get("mode", DEFAULT_MODE),
        backend=backend,
    )
    if task_count <= 1 or workers <= 1:
        backend = "serial"
        workers = 1
    return ExecutionPolicy(
        backend=backend,
        workers=workers,
        strict=strict,
        in_flight_limit=process_in_flight_limit(workers),
        math_threads=max(1, int(parallel.get("math_threads", 1))),
    )


def _selection_metadata(selection: Any, task_count: int) -> dict[str, Any]:
    if selection is None:
        return {}
    return {
        "native_frame_interval_ps": selection.native_interval_ps,
        "delta_time_ps": selection.delta_time_ps,
        "raw_frame_step": selection.raw_frame_step,
        "selected_frames": task_count,
        "total_frames": selection.total_frames,
    }


def _is_independent_gro_batch(sources: Sequence[Path]) -> bool:
    return len(sources) > 1 and all(path.suffix.lower() == ".gro" for path in sources)


def _unique_output_names(assignments: Sequence[Any]) -> dict[int, str]:
    ordered = sorted(assignments, key=lambda item: item.source_index)
    counts = Counter(item.path.stem.casefold() for item in ordered)
    seen: Counter[str] = Counter()
    used: set[str] = set()
    names: dict[int, str] = {}
    for item in ordered:
        key = item.path.stem.casefold()
        seen[key] += 1
        preferred = item.path.stem if counts[key] == 1 else f"{item.path.stem}_{seen[key]:03d}"
        candidate = preferred
        serial = 2
        while candidate.casefold() in used:
            candidate = f"{preferred}_{serial:03d}"
            serial += 1
        used.add(candidate.casefold())
        names[int(item.source_index)] = candidate
    return names


def _unique_path_names(paths: Sequence[Path]) -> dict[int, str]:
    class _PathRecord:
        def __init__(self, source_index: int, path: Path) -> None:
            self.source_index = source_index
            self.path = path

    return _unique_output_names(
        tuple(_PathRecord(index, path) for index, path in enumerate(paths))
    )


__all__ = ["build_run_plan"]
