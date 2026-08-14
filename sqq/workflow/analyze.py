from __future__ import annotations

"""Analyze workflow built on the shared typed runtime."""

from argparse import Namespace
from copy import deepcopy
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping, Sequence
import warnings

from ..config import (
    apply_cli_overrides,
    is_cpp_mode,
    load_config,
    normalize_analysis_scopes,
    output_enabled,
    refresh_resolution_report,
    validate_cpp_cli,
)
from ..io.gro_grouping import (
    GroGroupingResult,
    GroTopologyGroup,
    gro_topology_descriptor,
    scan_and_group_gro_inputs,
)
from ..io.lammps import inspect_lammps_topology_mapping
from ..io.output_cleanup import cleanup_previous_analyze_outputs
from ..io.tracking import write_tracking_result, write_tracking_tables
from ..io.render import (
    RenderSession,
    RenderSpec,
)
from ..io.reporting import (
    failed_row,
    write_run_config,
    write_summary,
)
from ..io.trajectory import expand_inputs, read_gro
from ..runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ..core.tracking import TrackingAccumulator, snapshot_from_frame_result
from ..models.tracking import TrackingConfig, TrackingResult
from ..runtime.frame import analyze_frame
from ..runtime.output_lock import output_lock
from ..ui.diagnostics import RunDiagnostics, capture_run_warnings
from ..ui.final_results import print_final_results, refresh_terminal
from ..ui.progress import (
    ParallelRunProgressDisplay,
    RunProgressDisplay,
    print_output_write_status,
)
from ..ui.run_header import (
    build_run_info,
    input_format_label,
    print_run_banner,
    print_run_header,
)
from ..ui.run_statistics import completed_run_statistics
from .analyze_plan import build_run_plan
from .session import AnalysisEvent, AnalysisRunner, AnalysisSink


def analyze(args: Namespace) -> None:
    """Run Analyze while exclusively owning its output root."""
    print_run_banner(getattr(args, "engine", None))
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    with output_lock(outdir):
        diagnostics = RunDiagnostics()
        with capture_run_warnings(diagnostics):
            try:
                _analyze_locked(args, diagnostics)
            finally:
                diagnostics.emit()


def _analyze_locked(args: Namespace, diagnostics: RunDiagnostics) -> None:
    run_started_at = datetime.now().astimezone()
    started_at = perf_counter()
    config = load_config(
        Path(args.config) if args.config else None,
        mode=getattr(args, "engine", None),
    )
    apply_cli_overrides(config, args)
    validate_cpp_cli(args, config)
    normalize_analysis_scopes(config)
    refresh_resolution_report(config)

    input_path = Path(args.input)
    paths = expand_inputs(
        input_path,
        pattern=str(config["input"]["pattern"]),
        recursive=bool(config["input"]["recursive"]),
    )
    topology = Path(args.topology) if args.topology else None
    config["input"]["format"] = input_format_label(paths)
    config["input"]["topology"] = str(topology.resolve()) if topology else None
    _resolve_lammps_mapping(config, args, topology)

    strict = bool(config.get("run", {}).get("strict", False))
    grouping = (
        scan_and_group_gro_inputs(paths, strict=strict)
        if _is_multi_gro(paths)
        else None
    )
    if grouping is not None:
        _validate_shared_gro_topology(topology, grouping)
    requested_output_types = list(config.get("output", {}).get("types", ()))
    grouping_warnings = _grouping_warnings(grouping)

    plan = build_run_plan(
        paths,
        config,
        Path(args.output),
        topology=topology,
        grouping=grouping,
    )
    plan, execution_config, group_configs = _resolve_plan_metadata(plan, grouping)
    cleanup_previous_analyze_outputs(Path(args.output), execution_config)
    for root in plan.output_roots:
        root.mkdir(parents=True, exist_ok=True)

    render_sessions, plan = _prepare_render_sessions(plan, execution_config, group_configs)
    tracking_sink, plan = _prepare_tracking_sink(
        plan,
        execution_config,
        group_configs,
        enabled=not bool(
            grouping is not None and grouping.info_only_fallback_required
        ),
    )
    workers = plan.policy.workers
    backend = plan.policy.backend
    extra_fields = None
    if grouping is not None:
        extra_fields = [
            ("Topology groups", grouping.group_count),
            ("Grouping policy", _grouping_policy(grouping)),
        ]
    print_run_header(
        args,
        execution_config,
        input_path,
        Path(args.output),
        paths,
        topology,
        workers,
        backend,
        run_started_at,
        extra_configuration_fields=extra_fields,
    )

    initial_info = build_run_info(
        args,
        execution_config,
        input_path,
        Path(args.output),
        paths,
        topology,
        workers,
        backend,
        0.0,
        run_started_at,
        run_started_at,
        [],
    )
    if grouping is not None:
        _add_grouping_metadata(
            initial_info,
            grouping,
            plan,
            group_configs,
            requested_output_types,
            grouping_warnings,
        )
    initial_info.update(status="running", error="")
    write_run_config(Path(args.output), execution_config, initial_info)
    _write_initial_group_configs(
        args,
        input_path,
        topology,
        grouping,
        plan,
        group_configs,
        requested_output_types,
        grouping_warnings,
        workers,
        backend,
        run_started_at,
    )

    progress = _ProgressBridge(plan, execution_config, group_configs, started_at)
    try:
        outcomes = AnalysisRunner(
            plan,
            event_sink=progress,
            sinks=((tracking_sink,) if tracking_sink is not None else ()),
        ).run()
        progress.close()
        rows_by_index = _rows_by_input_index(outcomes, grouping)
        analysis_completed_at = perf_counter()
        print_output_write_status()
        tracking_results = (
            tracking_sink.results if tracking_sink is not None else {}
        )
        for key, session in render_sessions.items():
            session.finalize(tracking=tracking_results.get(key))
        _write_analyze_tracking_outputs(plan, tracking_results)
    except Exception as exc:
        progress.close()
        for session in render_sessions.values():
            session.abort()
        _write_failed_run_configs(
            args,
            execution_config,
            input_path,
            paths,
            Path(args.output),
            topology,
            workers,
            backend,
            run_started_at,
            started_at,
            grouping,
            plan,
            group_configs,
            requested_output_types,
            grouping_warnings,
            locals().get("rows_by_index", {}),
            exc,
        )
        raise

    rows = [rows_by_index[index] for index in sorted(rows_by_index)]
    finished_at = datetime.now().astimezone()
    root_info = build_run_info(
        args,
        execution_config,
        input_path,
        Path(args.output),
        paths,
        topology,
        workers,
        backend,
        perf_counter() - started_at,
        run_started_at,
        finished_at,
        rows,
    )
    if grouping is not None:
        _add_grouping_metadata(
            root_info,
            grouping,
            plan,
            group_configs,
            requested_output_types,
            grouping_warnings,
        )
    root_info.update(status="completed", error="")
    try:
        _write_completed_reports(
            args,
            root_info,
            execution_config,
            input_path,
            topology,
            grouping,
            plan,
            group_configs,
            requested_output_types,
            grouping_warnings,
            rows_by_index,
            workers,
            backend,
            run_started_at,
            finished_at,
            perf_counter() - started_at,
        )
    except Exception as exc:
        root_info.update(status="failed", error=str(exc))
        write_run_config(Path(args.output), execution_config, root_info)
        raise

    final_finished_at = datetime.now().astimezone()
    total_seconds = perf_counter() - started_at
    analysis_seconds = max(0.0, analysis_completed_at - started_at)
    write_seconds = max(0.0, total_seconds - analysis_seconds)
    root_info.update(
        elapsed_seconds=round(total_seconds, 3),
        analysis_seconds=round(analysis_seconds, 3),
        write_seconds=round(write_seconds, 3),
        finish_time=final_finished_at.strftime("%H:%M:%S"),
        finished_at=final_finished_at.isoformat(timespec="seconds"),
    )
    final_config = _final_root_config(grouping, execution_config, group_configs)
    write_run_config(Path(args.output), final_config, root_info)
    statistics = completed_run_statistics(
        root_info,
        final_config,
        result_path=Path(args.output),
        requested_frames=int(plan.sampling.get("selected_frames", len(paths))),
        analysis_seconds=analysis_seconds,
        write_seconds=write_seconds,
        total_seconds=total_seconds,
        track=_has_tracking_state(plan),
    )
    statistics["diagnostic_messages"] = diagnostics.consume()
    refresh_terminal()
    print_final_results(root_info, final_config, statistics)


class _ProgressBridge:
    """Translate workflow events into the established terminal display."""

    def __init__(
        self,
        plan: RunPlan,
        config: Mapping[str, Any],
        group_configs: Mapping[int, Mapping[str, Any]],
        started_at: float,
    ) -> None:
        configs = tuple(group_configs.values()) or (config,)
        include_cluster = any(
            bool(item.get("hydrate_cluster", {}).get("enabled", False))
            for item in configs
        )
        include_patch = any(
            bool(item.get("half_cage", {}).get("enabled", False))
            or bool(item.get("quasi_cage", {}).get("enabled", False))
            for item in configs
        )
        cpp = all(is_cpp_mode(item.get("mode")) for item in configs)
        if plan.policy.backend != "serial" and plan.policy.workers > 1:
            self._display: RunProgressDisplay | ParallelRunProgressDisplay = (
                ParallelRunProgressDisplay(
                    len(plan.tasks),
                    plan.policy.workers,
                    started_at,
                    include_cluster,
                    cpp_mode=cpp,
                    include_patch_stage=include_patch,
                )
            )
            self._parallel = True
        else:
            self._display = RunProgressDisplay(
                len(plan.tasks),
                started_at,
                include_cluster,
                cpp_mode=cpp,
                include_patch_stage=include_patch,
            )
            self._parallel = False
        self._closed = False

    def __call__(self, event: AnalysisEvent) -> None:
        task = event.task
        if task is None:
            return
        if event.kind == "task-start":
            if self._parallel:
                self._display.start_file(task.task_index, task.display_name)  # type: ignore[union-attr]
            else:
                self._display.start_frame(task.task_index, task.display_name)  # type: ignore[union-attr]
        elif event.kind == "stage" and event.stage:
            if self._parallel:
                self._display.update_stage(task.task_index, event.stage)  # type: ignore[union-attr]
            else:
                self._display.update_stage(event.stage)  # type: ignore[union-attr]
        elif event.kind in {"task-complete", "task-cancelled"}:
            success = event.kind == "task-complete" and event.status == "ok"
            if self._parallel:
                self._display.complete_file(task.task_index, success)  # type: ignore[union-attr]
            else:
                self._display.complete_frame(success)  # type: ignore[union-attr]

    def close(self) -> None:
        if not self._closed:
            self._display.close()
            self._closed = True


def _resolve_lammps_mapping(
    config: dict[str, Any],
    args: Namespace,
    topology: Path | None,
) -> None:
    if not str(config["input"]["format"]).startswith("lammps-"):
        return
    lammps = config["input"]["lammps"]
    explicit = bool(lammps.get("type_map"))
    resolved = {}
    rebuilt = False
    if topology is not None:
        resolved, rebuilt = inspect_lammps_topology_mapping(topology, lammps)
    lammps["resolved_type_map"] = {
        type_id: (
            {"ignore": True}
            if entry.ignore
            else {"resname": entry.resname, "atomname": entry.atomname}
        )
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    }
    mapping = ", ".join(
        f"{type_id}=ignore"
        if entry.ignore
        else f"{type_id}={entry.resname}/{entry.atomname}"
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    )
    if explicit:
        source = str(Path(args.config).resolve()) if args.config else "<configuration>"
    else:
        source = "auto (DATA topology)"
        if rebuilt:
            source += "; molecule IDs rebuilt from Bonds"
    lammps["type_map_source"] = f"{source}: {mapping}" if mapping else source


def _resolve_plan_metadata(
    plan: RunPlan,
    grouping: GroGroupingResult | None,
) -> tuple[RunPlan, dict[str, Any], dict[int, dict[str, Any]]]:
    config = deepcopy(dict(plan.context.config))
    groups = {
        int(key): deepcopy(dict(value))
        for key, value in plan.context.group_configs.items()
    }
    if grouping is not None and groups:
        modes: dict[str, str] = {}
        reasons: dict[str, str] = {}
        for group in grouping.groups:
            label = "result" if grouping.group_count == 1 else f"result_{group.label}"
            graph = groups[group.group_index]["graph"]
            modes[label] = str(graph["effective_bond_mode"])
            reason = str(graph.get("effective_bond_mode_reason", ""))
            if reason:
                reasons[label] = reason
        config["graph"]["effective_bond_mode_by_group"] = modes
        config["graph"]["effective_bond_mode_reason_by_group"] = reasons
        unique_modes = set(modes.values())
        if len(unique_modes) == 1:
            config["graph"]["effective_bond_mode"] = next(iter(unique_modes))
            unique_reasons = set(reasons.values())
            if len(unique_reasons) == 1:
                config["graph"]["effective_bond_mode_reason"] = next(iter(unique_reasons))
        else:
            config["graph"].pop("effective_bond_mode", None)
            config["graph"].pop("effective_bond_mode_reason", None)
    context = replace(plan.context, config=config, group_configs=groups)
    return replace(plan, context=context), config, groups


def _prepare_render_sessions(
    plan: RunPlan,
    config: dict[str, Any],
    group_configs: Mapping[int, dict[str, Any]],
) -> tuple[dict[int | str, RenderSession], RunPlan]:
    sessions: dict[int | str, RenderSession] = {}
    group_fragment_dirs: dict[int, Path] = {}
    if group_configs:
        for key, group_config in group_configs.items():
            root = Path(plan.context.group_output_roots[key])
            RenderSession.cleanup_output(root)
            if not output_enabled(group_config, "sqq-render"):
                continue
            graph = group_config["graph"]
            session = RenderSession.create(
                root,
                RenderSpec(
                    component_roles=group_config,
                    requested_graph_mode=graph.get("bond_mode"),
                    effective_graph_mode=graph.get("effective_bond_mode"),
                    atom_scope=group_config.get("render", {}).get("atom_scope", "full"),
                ),
            )
            sessions[key] = session
            group_fragment_dirs[key] = session.fragment_dir
        context = replace(plan.context, group_fragment_dirs=group_fragment_dirs)
    else:
        root = plan.context.output_root
        RenderSession.cleanup_output(root)
        if output_enabled(config, "sqq-render"):
            graph = config["graph"]
            session = RenderSession.create(
                root,
                RenderSpec(
                    component_roles=config,
                    requested_graph_mode=graph.get("bond_mode"),
                    effective_graph_mode=graph.get("effective_bond_mode"),
                    atom_scope=config.get("render", {}).get("atom_scope", "full"),
                ),
            )
            sessions["run"] = session
            context = replace(plan.context, fragment_dir=session.fragment_dir)
        else:
            context = plan.context
    return sessions, replace(plan, context=context)


_TRACKING_CONFIG_FIELDS = frozenset(
    {
        "min_jaccard",
        "min_shared_fraction",
        "min_shared_waters",
        "max_center_distance_nm",
        "gap_frame",
        "guest_tiebreak",
    }
)


class _AnalyzeTrackingSink(AnalysisSink):
    """Build one bounded-memory persistent-ID stream per topology group."""

    def __init__(
        self,
        configs: Mapping[int | str, Mapping[str, Any]],
    ) -> None:
        self._configs = dict(configs)
        self._accumulators: dict[int | str, TrackingAccumulator] = {}
        self._expected: dict[int | str, int] = {}
        self._failed: set[int | str] = set()
        self.results: dict[int | str, TrackingResult] = {}

    def start(self, plan: RunPlan) -> None:
        self._accumulators = {
            key: TrackingAccumulator(_tracking_config(config))
            for key, config in self._configs.items()
        }
        self._expected = {key: 0 for key in self._configs}
        self._failed.clear()
        self.results.clear()

    def consume(self, task: FrameTask, outcome: TaskOutcome) -> None:
        key: int | str = task.group_key if task.group_key is not None else "run"
        if key not in self._accumulators:
            raise RuntimeError(f"No tracking accumulator for topology group {key!r}.")
        expected = self._expected[key]
        if int(task.frame_index) != expected:
            raise RuntimeError(
                f"Topology group {key!r} tracking frame order is not contiguous: "
                f"expected {expected}, got {task.frame_index}."
            )
        self._expected[key] = expected + 1
        if key in self._failed:
            return
        if not outcome.ok or outcome.result is None:
            self._failed.add(key)
            return
        self._accumulators[key].add(
            snapshot_from_frame_result(outcome.result, int(task.frame_index))
        )

    def finish(
        self,
        plan: RunPlan,
        outcomes: Sequence[TaskOutcome],
    ) -> None:
        for key, accumulator in self._accumulators.items():
            if key in self._failed:
                warnings.warn(
                    "Persistent Track state was not written for topology group "
                    f"{key!r} because one or more frames failed analysis.",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            self.results[key] = accumulator.result()


def _tracking_config(config: Mapping[str, Any]) -> TrackingConfig:
    values = config.get("track", {})
    if not isinstance(values, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    return TrackingConfig.from_mapping(
        {
            name: values[name]
            for name in _TRACKING_CONFIG_FIELDS
            if name in values
        }
    )


def _prepare_tracking_sink(
    plan: RunPlan,
    config: Mapping[str, Any],
    group_configs: Mapping[int, Mapping[str, Any]],
    *,
    enabled: bool,
) -> tuple[_AnalyzeTrackingSink | None, RunPlan]:
    if not enabled:
        return None, plan
    configs: dict[int | str, Mapping[str, Any]] = (
        dict(group_configs) if group_configs else {"run": config}
    )
    sink = _AnalyzeTrackingSink(configs)
    context = replace(
        plan.context,
        retain_results=False,
        stream_results=True,
    )
    return sink, replace(plan, context=context)


def _write_analyze_tracking_outputs(
    plan: RunPlan,
    results: Mapping[int | str, TrackingResult],
) -> None:
    for key, result in results.items():
        root = (
            Path(plan.context.group_output_roots[key])
            if key != "run"
            else Path(plan.context.output_root)
        )
        track_root = root / "track"
        write_tracking_tables(result, track_root)
        write_tracking_result(result, track_root / "track_state.json")


def _rows_by_input_index(
    outcomes: tuple[TaskOutcome, ...],
    grouping: GroGroupingResult | None,
) -> dict[int, dict[str, Any]]:
    rows = {outcome.task_index: dict(outcome.row) for outcome in outcomes}
    if grouping is not None:
        for failure in grouping.failures:
            rows[failure.source_index] = failed_row(
                failure.path.stem,
                str(failure.path),
                failure.error,
            )
    return rows


def _write_completed_reports(
    args: Namespace,
    root_info: dict[str, Any],
    config: dict[str, Any],
    input_path: Path,
    topology: Path | None,
    grouping: GroGroupingResult | None,
    plan: RunPlan,
    group_configs: Mapping[int, dict[str, Any]],
    requested_output_types: list[str],
    grouping_warnings: list[str],
    rows_by_index: Mapping[int, dict[str, Any]],
    workers: int,
    backend: str,
    started_at: datetime,
    finished_at: datetime,
    elapsed: float,
) -> None:
    root = plan.context.output_root
    all_rows = [rows_by_index[index] for index in sorted(rows_by_index)]
    if grouping is None:
        root_info["summary_write"] = write_summary(
            all_rows,
            root,
            config,
            write_xlsx=output_enabled(config, "summary-xlsx"),
            run_info=root_info,
        )
        write_run_config(root, config, root_info)
        return
    if grouping.info_only_fallback_required or not grouping.groups:
        write_run_config(root, config, root_info)
        return
    if grouping.group_count == 1:
        group = grouping.groups[0]
        group_config = group_configs[group.group_index]
        root_info["topology_group"] = "single"
        root_info["topology_fingerprint"] = group.fingerprint
        root_info["summary_write"] = write_summary(
            all_rows,
            root,
            group_config,
            write_xlsx=output_enabled(group_config, "summary-xlsx"),
            run_info=root_info,
        )
        write_run_config(root, group_config, root_info)
        return

    metrics: dict[str, Any] = {}
    for group in grouping.groups:
        group_rows = [
            rows_by_index[index]
            for index in group.source_indices
            if index in rows_by_index
        ]
        group_info = _group_run_info(
            args,
            input_path,
            topology,
            group,
            group_configs[group.group_index],
            Path(plan.context.group_output_roots[group.group_index]),
            workers,
            backend,
            started_at,
            finished_at,
            elapsed,
            group_rows,
            grouping,
            requested_output_types,
            grouping_warnings,
        )
        group_info.update(status="completed", error="")
        metric = write_summary(
            group_rows,
            Path(plan.context.group_output_roots[group.group_index]),
            group_configs[group.group_index],
            write_xlsx=output_enabled(group_configs[group.group_index], "summary-xlsx"),
            run_info=group_info,
        )
        group_info["summary_write"] = metric
        write_run_config(
            Path(plan.context.group_output_roots[group.group_index]),
            group_configs[group.group_index],
            group_info,
        )
        metrics[str(group.label)] = metric
    root_info["summary_write"] = {"groups": metrics}
    write_run_config(root, config, root_info)


def _write_initial_group_configs(
    args: Namespace,
    input_path: Path,
    topology: Path | None,
    grouping: GroGroupingResult | None,
    plan: RunPlan,
    group_configs: Mapping[int, dict[str, Any]],
    requested_output_types: list[str],
    grouping_warnings: list[str],
    workers: int,
    backend: str,
    started_at: datetime,
) -> None:
    if grouping is None or grouping.group_count <= 1 or grouping.info_only_fallback_required:
        return
    for group in grouping.groups:
        info = _group_run_info(
            args,
            input_path,
            topology,
            group,
            group_configs[group.group_index],
            Path(plan.context.group_output_roots[group.group_index]),
            workers,
            backend,
            started_at,
            started_at,
            0.0,
            [],
            grouping,
            requested_output_types,
            grouping_warnings,
        )
        info.update(status="running", error="")
        write_run_config(
            Path(plan.context.group_output_roots[group.group_index]),
            group_configs[group.group_index],
            info,
        )


def _write_failed_run_configs(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    paths: list[Path],
    root: Path,
    topology: Path | None,
    workers: int,
    backend: str,
    started_at_wall: datetime,
    started_at: float,
    grouping: GroGroupingResult | None,
    plan: RunPlan,
    group_configs: Mapping[int, dict[str, Any]],
    requested_output_types: list[str],
    grouping_warnings: list[str],
    rows_by_index: Mapping[int, dict[str, Any]],
    error: Exception,
) -> None:
    failed_at = datetime.now().astimezone()
    rows = [rows_by_index[index] for index in sorted(rows_by_index)]
    info = build_run_info(
        args,
        config,
        input_path,
        root,
        paths,
        topology,
        workers,
        backend,
        perf_counter() - started_at,
        started_at_wall,
        failed_at,
        rows,
    )
    if grouping is not None:
        _add_grouping_metadata(
            info,
            grouping,
            plan,
            group_configs,
            requested_output_types,
            grouping_warnings,
        )
    info.update(status="failed", error=str(error))
    try:
        write_run_config(root, config, info)
    except Exception:
        pass
    if grouping is None or grouping.group_count <= 1 or grouping.info_only_fallback_required:
        return
    for group in grouping.groups:
        group_rows = [rows_by_index[i] for i in group.source_indices if i in rows_by_index]
        group_info = _group_run_info(
            args,
            input_path,
            topology,
            group,
            group_configs[group.group_index],
            Path(plan.context.group_output_roots[group.group_index]),
            workers,
            backend,
            started_at_wall,
            failed_at,
            perf_counter() - started_at,
            group_rows,
            grouping,
            requested_output_types,
            grouping_warnings,
        )
        group_info.update(status="failed", error=str(error))
        try:
            write_run_config(
                Path(plan.context.group_output_roots[group.group_index]),
                group_configs[group.group_index],
                group_info,
            )
        except Exception:
            pass


def _group_run_info(
    args: Namespace,
    input_path: Path,
    topology: Path | None,
    group: GroTopologyGroup,
    config: dict[str, Any],
    outdir: Path,
    workers: int,
    backend: str,
    started_at: datetime,
    finished_at: datetime,
    elapsed: float,
    rows: list[dict[str, Any]],
    grouping: GroGroupingResult,
    requested_output_types: list[str],
    grouping_warnings: list[str],
) -> dict[str, Any]:
    info = build_run_info(
        args,
        config,
        input_path,
        outdir,
        list(group.paths),
        topology,
        workers,
        backend,
        elapsed,
        started_at,
        finished_at,
        rows,
    )
    info.update(grouping.limit_metadata())
    info.update(
        topology_grouping=_grouping_policy(grouping),
        topology_group=group.label,
        topology_fingerprint=group.fingerprint,
        topology_source_mapping=[
            item
            for item in grouping.source_mapping()
            if item.get("group_index") == group.group_index
        ],
        requested_output_types=list(requested_output_types),
        output_policy="mode/configured outputs",
        warnings=list(grouping_warnings),
    )
    return info


def _add_grouping_metadata(
    info: dict[str, Any],
    grouping: GroGroupingResult,
    plan: RunPlan,
    group_configs: Mapping[int, dict[str, Any]],
    requested_output_types: list[str],
    grouping_warnings: list[str],
) -> None:
    info.update(grouping.limit_metadata())
    info["topology_grouping"] = _grouping_policy(grouping)
    info["topology_groups"] = [
        {
            "group_index": group.group_index,
            "group_label": group.label,
            "fingerprint": group.fingerprint,
            "source_count": len(group.inputs),
            "output_dir": str(
                Path(plan.context.group_output_roots[group.group_index]).resolve()
            ),
            "effective_graph_mode": group_configs[group.group_index]["graph"].get(
                "effective_bond_mode",
                group_configs[group.group_index]["graph"]["bond_mode"],
            ),
            "effective_graph_mode_reason": group_configs[group.group_index]["graph"].get(
                "effective_bond_mode_reason", ""
            ),
            "sources": [str(path.resolve()) for path in group.paths],
        }
        for group in grouping.groups
    ]
    info["topology_source_mapping"] = list(grouping.source_mapping())
    info["requested_output_types"] = list(requested_output_types)
    info["output_policy"] = (
        "forced info-only"
        if grouping.info_only_fallback_required
        else "mode/configured outputs"
    )
    info["warnings"] = list(grouping_warnings)


def _grouping_warnings(grouping: GroGroupingResult | None) -> list[str]:
    if grouping is None or not grouping.info_only_fallback_required:
        return []
    message = (
        f"Detected {grouping.group_count} independent GRO topology groups, above "
        f"the supported A-Z limit ({grouping.group_limit}). All readable files "
        "will be analyzed with info-only output; summaries, GRO output, and VMD "
        "output are disabled for this run."
    )
    warnings.warn(message, UserWarning, stacklevel=3)
    return [message]


def _grouping_policy(grouping: GroGroupingResult) -> str:
    if grouping.info_only_fallback_required:
        return "info-only (>26 topology groups)"
    if grouping.group_count <= 1:
        return "single result root"
    return "separate result_A ... result_Z roots"


def _validate_shared_gro_topology(
    topology: Path | None,
    grouping: GroGroupingResult,
) -> None:
    if topology is None:
        return
    if topology.suffix.lower() != ".gro":
        raise ValueError(
            "Multiple independent GRO inputs require a GRO file for --top/-t "
            "so SQQ can validate every source topology."
        )
    descriptor = gro_topology_descriptor(read_gro(topology))
    mismatches = [
        str(item.path)
        for item in grouping.assignments
        if item.descriptor != descriptor
    ]
    if mismatches:
        raise ValueError(
            f"Shared topology {topology} does not match these GRO source(s):\n  - "
            + "\n  - ".join(mismatches)
        )


def _final_root_config(
    grouping: GroGroupingResult | None,
    config: dict[str, Any],
    group_configs: Mapping[int, dict[str, Any]],
) -> dict[str, Any]:
    if grouping is not None and grouping.group_count == 1 and grouping.groups:
        return group_configs[grouping.groups[0].group_index]
    return config


def _has_tracking_state(plan: RunPlan) -> bool:
    if plan.topology_groups:
        for key, task_indexes in plan.topology_groups.items():
            if len(task_indexes) < 2:
                continue
            root = Path(plan.context.group_output_roots[key])
            if (root / "track" / "track_state.json").is_file():
                return True
        return False
    return len(plan.tasks) > 1 and any(
        (root / "track" / "track_state.json").is_file()
        for root in plan.output_roots
    )


def _is_multi_gro(paths: list[Path]) -> bool:
    return len(paths) > 1 and all(path.suffix.lower() == ".gro" for path in paths)


__all__ = [
    "AnalysisEvent",
    "AnalysisRunner",
    "AnalysisSink",
    "analyze",
    "analyze_frame",
    "build_run_plan",
]
