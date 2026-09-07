"""Cross-frame cage tracking workflow for ``sqq track``."""

from __future__ import annotations

from argparse import Namespace
from copy import deepcopy
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from typing import Any, Mapping, Sequence
import warnings

import numpy as np

from ... import __release_date__, __version__
from ...config import (
    DEFAULT_MODE,
    apply_cli_overrides,
    engine_display,
    is_cpp_mode,
    load_config,
    normalize_engine_capabilities,
    normalize_analysis_scopes,
    normalize_mode,
    refresh_resolution_report,
    validate_cpp_cli,
)
from ...presentation.citation import completed_citation_evidence
from ...core.sqq_cpp import require_native
from ...core.tracking import (
    TrackingAccumulator,
    parse_targets,
    select_targets,
)
from ...io.render import (
    RenderSession,
    RenderSpec,
    RenderBundle,
    TRACK_MEMBERSHIP_NAME,
    TRACK_RENDER_DIRECTORY,
    discover_sqq_cage_bundle,
    publish_target_render_bundle,
    validate_tracking_source_bundle,
)
from ...io.reporting import write_run_config
from ...io.tracking import (
    append_track_info_section,
    discover_track_state,
    read_tracking_result,
    target_directory_name,
    write_track_outputs,
)
from ...io.input.trajectory import effective_frame_time_ps, expand_inputs
from ...models.tracking import (
    TargetSelection,
    TargetSpec,
    TrackingConfig,
    TrackingResult,
)
from ...runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ...runtime.output_lock import reserve_output_directory
from ...ui.diagnostics import RunDiagnostics, capture_run_warnings
from ...ui.formatting import format_time_zone
from ...ui.progress import (
    RunProgressDisplay,
    print_output_directory_notice,
)
from ...ui.run_header import input_format_label, print_run_banner, print_track_header
from ...ui import completed_run_statistics, print_final_results, refresh_terminal
from ...runtime.session import AnalysisEvent, AnalysisRunner, AnalysisSink
from ..analyze.plan import build_run_plan
from ..analyze.tracking import tracking_snapshot
from .precursor import (
    raw_precursor_histories as _raw_precursor_histories,
    write_precursor_outputs as _write_precursor_outputs,
)
from .request import (
    explicit_tracking_fields as _explicit_tracking_fields,
    prepare_lammps_metadata as _prepare_lammps_metadata,
    prepare_target_capabilities as _prepare_target_capabilities,
    resolve_track_request as _resolve_track_request,
    source_engine_selector as _source_engine_selector,
    tracking_config as _tracking_config,
    update_track_config as _update_track_config,
    validate_requested_tracks as _validate_requested_tracks,
)


__all__ = ["track"]


class _TrackingAnalysisSink(AnalysisSink):
    """Link successful frame results in authoritative plan order."""

    def __init__(self, accumulator: TrackingAccumulator) -> None:
        self.accumulator = accumulator
        self.consumed = 0

    def start(self, plan: RunPlan) -> None:
        self.consumed = 0

    def consume(self, task: FrameTask, outcome: TaskOutcome) -> None:
        if not outcome.ok:
            raise RuntimeError(
                "Tracking cannot skip failed frame "
                f"{task.display_name!r}: {outcome.error_message or 'analysis failed'}"
            )
        if outcome.snapshot_error:
            raise RuntimeError(
                f"Tracking frame {task.display_name!r} could not be reduced to "
                f"a persistent cage snapshot: {outcome.snapshot_error}"
            )
        snapshot = tracking_snapshot(task, outcome)
        if snapshot is None:
            raise RuntimeError(
                f"Tracking frame {task.display_name!r} did not retain its analysis result."
            )
        self.accumulator.add(snapshot)
        self.consumed += 1

    def finish(
        self,
        plan: RunPlan,
        outcomes: Sequence[TaskOutcome],
    ) -> None:
        if self.consumed != len(plan.tasks):
            raise RuntimeError(
                "Tracking linked "
                f"{self.consumed} frames; expected {len(plan.tasks)}."
            )


class _TrackProgressSink:
    """Translate shared runner events into the Track progress panel."""

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


def _track_executed_features(
    config: Mapping[str, Any],
    *,
    frames: int,
    raw_mode: bool,
    has_guest_observations: bool,
) -> dict[str, bool]:
    """Describe work that this Track invocation actually executed."""
    ran = int(frames) > 0
    temporal_ran = int(frames) >= 2
    features = {
        "water_network": False,
        "ring_topology": False,
        "cage_topology": False,
        "half_cage": False,
        "quasi_cage": False,
        "cage_isomer": False,
        "cage_occupancy": ran and has_guest_observations,
        "hydrate_phase_domain": False,
        "vmd_rendering": ran,
        "cage_tracking": ran,
        "cage_transition": temporal_ran,
        "cage_lifetime": temporal_ran,
        "guest_residence": temporal_ran and has_guest_observations,
    }
    if not ran or not raw_mode:
        return features
    features.update(
        {
            "water_network": True,
            "ring_topology": True,
            "cage_topology": True,
            "half_cage": _section_enabled(config, "half_cage"),
            "quasi_cage": _section_enabled(config, "quasi_cage"),
            "cage_isomer": True,
            "cage_occupancy": has_guest_observations,
            "hydrate_phase_domain": _section_enabled(
                config, "hydrate_cluster"
            ),
        }
    )
    return features


def _section_enabled(config: Mapping[str, Any], name: str) -> bool:
    section = config.get(name, {})
    return isinstance(section, Mapping) and bool(section.get("enabled", False))


def track(args: Namespace) -> None:
    """Run Track while exclusively owning the output root."""
    print_run_banner(getattr(args, "engine", None))
    started_wall = datetime.now().astimezone()
    started = perf_counter()
    diagnostics = RunDiagnostics()
    try:
        with capture_run_warnings(diagnostics):
            _preflight_track_paths(args)
            prepared = _prepare_track(args)
            with reserve_output_directory(Path(args.output)) as output_selection:
                output_selection.apply_to_config(prepared.config)
                args.output = str(output_selection.resolved)
                args.output_requested = str(output_selection.requested)
                args.output_auto_renamed = output_selection.auto_renamed
                print_output_directory_notice(
                    output_selection.requested,
                    output_selection.resolved,
                    auto_renamed=output_selection.auto_renamed,
                )
                _track_locked(
                    args,
                    diagnostics,
                    prepared,
                    started_wall=started_wall,
                    started=started,
                )
    except BaseException:
        diagnostics.consume()
        raise
    finally:
        diagnostics.emit()


def _preflight_track_paths(args: Namespace) -> None:
    """Reject missing explicit paths before creating the output root."""
    for attribute, label in (
        ("config", "Configuration file"),
        ("input", "Input"),
        ("topology", "Topology file"),
        ("pair", "Pair file"),
        ("source", "Analyze result"),
    ):
        raw = getattr(args, attribute, None)
        if not raw:
            continue
        path = Path(raw)
        if attribute in {"input", "source"}:
            exists = path.exists()
        else:
            exists = path.is_file()
        if not exists:
            raise FileNotFoundError(f"{label} does not exist: {path}")
    output = Path(args.output)
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(f"Output path is not a directory: {output}")


@dataclass(frozen=True, slots=True)
class _TrackPreparation:
    config: dict[str, Any]
    targets: tuple[TargetSpec, ...]
    source_state_path: Path | None
    requested_tracking: TrackingConfig


def _prepare_track(args: Namespace) -> _TrackPreparation:
    """Resolve Track configuration and sources without creating outputs."""
    config = load_config(
        Path(args.config) if getattr(args, "config", None) else None,
        mode=getattr(args, "engine", None),
    )
    apply_cli_overrides(config, args)
    target_value, source_value = _resolve_track_request(config, args)
    setattr(args, "source", source_value)
    if getattr(args, "input", None) and source_value:
        raise ValueError("Use either --input or --source, not both.")
    targets = parse_targets(target_value)
    _update_track_config(config, args, targets)
    source_state_path: Path | None = None
    if not getattr(args, "input", None):
        source_root = Path(source_value or Path.cwd())
        source_state_path = discover_track_state(source_root)
        source_engine = _source_engine_selector(source_state_path)
        if source_engine is not None:
            _check_source_engine_request(
                getattr(args, "engine", None), source_engine
            )
            config["mode"] = source_engine
            normalize_engine_capabilities(config, emit_warnings=False)
    _prepare_target_capabilities(config, targets)
    validate_cpp_cli(args, config)
    if getattr(args, "input", None) and is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        # Raw Track analyzes frames; fail once here instead of per frame.
        require_native()
    normalize_analysis_scopes(config)
    refresh_resolution_report(config)
    requested_tracking = _tracking_config(config)
    return _TrackPreparation(
        config=config,
        targets=tuple(targets),
        source_state_path=source_state_path,
        requested_tracking=requested_tracking,
    )


def _check_source_engine_request(requested: object, source_engine: str) -> None:
    """Never let an explicit ``-e`` silently pretend to change imported science."""
    if requested in (None, ""):
        return
    requested_mode = normalize_mode(str(requested))
    if requested_mode == source_engine:
        return
    if is_cpp_mode(requested_mode) != is_cpp_mode(source_engine):
        raise ValueError(
            f"Explicit engine {requested_mode!r} conflicts with the "
            f"{engine_display(source_engine)} Analyze source. --source imports "
            "completed analysis and cannot rerun it with another engine; omit -e "
            "or rerun Track from the raw trajectory."
        )
    warnings.warn(
        f"Explicit engine {requested_mode!r} differs from the source engine "
        f"{source_engine!r}; the imported Analyze engine is used.",
        UserWarning,
        stacklevel=2,
    )


def _track_locked(
    args: Namespace,
    diagnostics: RunDiagnostics,
    prepared: _TrackPreparation,
    *,
    started_wall: datetime,
    started: float,
) -> None:
    output = Path(args.output)
    config = prepared.config
    targets = prepared.targets
    source_state_path = prepared.source_state_path
    requested_tracking = prepared.requested_tracking
    published_render_scripts: list[str] = []

    print_track_header(args, config, targets, started_wall)

    try:
        with TemporaryDirectory(
            prefix=".sqq-track-",
            dir=output,
            ignore_cleanup_errors=True,
        ) as temporary:
            if getattr(args, "input", None):
                result, source_bundle, analysis_plan, metadata = _track_input(
                    args,
                    config,
                    targets,
                    started,
                    requested_tracking,
                    source_root=Path(temporary) / "source",
                )
                source = None
            else:
                source = Path(getattr(args, "source", None) or Path.cwd())
                if source_state_path is None:
                    raise RuntimeError("Track source state was not resolved.")
                state_path = source_state_path
                result = read_tracking_result(state_path)
                explicit_fields = _explicit_tracking_fields(args)
                mismatches = [
                    name
                    for name in explicit_fields
                    if getattr(requested_tracking, name) != getattr(result.config, name)
                ]
                if mismatches:
                    raise ValueError(
                        "Explicit tracking setting(s) differ from the imported "
                        "track_state.json: "
                        + ", ".join(sorted(mismatches))
                        + ". Matching settings cannot be reapplied in --source mode; "
                        "rerun Track from the raw trajectory."
                    )
                config.setdefault("track", {}).update(result.config.to_dict())
                source_bundle = discover_sqq_cage_bundle(
                    source, state_path=state_path
                )
                analysis_plan = None
                metadata = {
                    "source": str(source.resolve()),
                    "state": str(state_path.resolve()),
                    "input": None,
                    "topology": None,
                    "input_format": "track-state",
                    "selected_frames": len(result.frames),
                    "source_frames_total": len(result.frames),
                    "native_frame_interval_ps": _native_interval(result),
                    "raw_frame_step": 1,
                }

            _validate_requested_tracks(result, targets)
            selections = select_targets(
                result, (target.raw for target in targets)
            )
            guest_ids = {
                guest_id
                for selection in selections
                for observation in selection.observations
                for guest_id in observation.guest_ids
            }
            has_guest_observations = bool(guest_ids)
            precursor_data = (
                _raw_precursor_histories(selections, analysis_plan, config)
                if analysis_plan is not None
                else None
            )
            analysis_completed_at = perf_counter()
            validate_tracking_source_bundle(
                source_bundle,
                result=result,
            )
            written = write_track_outputs(
                result,
                output,
                targets=[target.raw for target in targets],
            )
            for selection in selections:
                directory = written[target_directory_name(selection.target)]
                bundle = publish_target_render_bundle(
                    selection,
                    directory,
                    source_bundle,
                )
                if bundle.complete and bundle.script_path is not None:
                    published_render_scripts.append(str(bundle.script_path.resolve()))

        _write_precursor_outputs(
            selections,
            written,
            precursor_data,
            source_mode=not bool(getattr(args, "input", None)),
        )
        elapsed = perf_counter() - started
        run_info = _run_info(
            config,
            args,
            metadata,
            targets,
            result,
            started_wall,
            datetime.now().astimezone(),
            elapsed,
            status="completed",
        )
        raw_mode = bool(getattr(args, "input", None))
        track_features = _track_executed_features(
            config,
            frames=len(result.frames),
            raw_mode=raw_mode,
            has_guest_observations=has_guest_observations,
        )
        citation_evidence = completed_citation_evidence(
            config,
            successful_frames=len(result.frames),
            completed_outputs=("sqq-render",),
            track=True,
            occupancy_evaluated=has_guest_observations,
            guest_residence_evaluated=has_guest_observations,
        )
        citation_evidence["executed_features"] = track_features
        if not raw_mode:
            citation_evidence["executed_order_parameters"] = ()
        run_info.update(citation_evidence)
        run_info["guest_molecules"] = len(guest_ids)
        write_run_config(output, config, run_info)
        _write_root_track_info(output, config, metadata, targets, result, elapsed)

        final_finished_at = datetime.now().astimezone()
        total_seconds = perf_counter() - started
        analysis_seconds = max(0.0, analysis_completed_at - started)
        write_seconds = max(0.0, total_seconds - analysis_seconds)
        run_info["elapsed_seconds"] = round(total_seconds, 3)
        run_info["analysis_seconds"] = round(analysis_seconds, 3)
        run_info["write_seconds"] = round(write_seconds, 3)
        run_info["finish_time"] = final_finished_at.strftime("%H:%M:%S")
        run_info["finished_at"] = final_finished_at.isoformat(timespec="seconds")
        write_run_config(output, config, run_info)
        statistics = completed_run_statistics(
            run_info,
            config,
            result_path=output,
            requested_frames=int(metadata.get("selected_frames", len(result.frames))),
            analysis_seconds=analysis_seconds,
            write_seconds=write_seconds,
            total_seconds=total_seconds,
            track=True,
        )
        completed_outputs = tuple(statistics.get("completed_outputs", ()))
        if "sqq-render" not in completed_outputs:
            statistics["completed_outputs"] = completed_outputs + ("sqq-render",)
        statistics["render_script_paths"] = tuple(published_render_scripts)
        statistics["diagnostic_messages"] = diagnostics.consume()
        refresh_terminal()
        print_final_results(run_info, config, statistics)
    except Exception as exc:
        run_info = _run_info(
            config,
            args,
            {
                "source": str(getattr(args, "source", "") or ""),
                "input": str(getattr(args, "input", "") or ""),
                "topology": str(getattr(args, "topology", "") or ""),
                "input_format": "",
                "selected_frames": 0,
                "source_frames_total": 0,
                "native_frame_interval_ps": None,
                "raw_frame_step": 1,
            },
            targets,
            None,
            started_wall,
            datetime.now().astimezone(),
            perf_counter() - started,
            status="failed",
            error=str(exc),
        )
        write_run_config(output, config, run_info)
        raise

def _track_input(
    args: Namespace,
    config: dict[str, Any],
    targets: Sequence[TargetSpec],
    started: float,
    tracking_config: TrackingConfig,
    *,
    source_root: Path,
) -> tuple[
    TrackingResult,
    RenderBundle,
    RunPlan,
    dict[str, Any],
]:
    input_path = Path(args.input)
    pattern = getattr(args, "pattern", None) or config["input"]["pattern"]
    recursive = bool(
        getattr(args, "recursive", False) or config["input"]["recursive"]
    )
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)
    if len(paths) != 1:
        raise ValueError(
            "sqq track accepts one trajectory or one stacked GRO file; "
            f"{len(paths)} files were matched."
        )
    trajectory = paths[0]
    topology = Path(args.topology) if getattr(args, "topology", None) else None
    config["input"]["format"] = input_format_label(paths)
    config["input"]["topology"] = (
        str(topology.resolve()) if topology is not None else None
    )
    _prepare_lammps_metadata(config, trajectory, topology, args)
    config["parallel"]["workers"] = 1
    config["parallel"]["backend"] = "serial"
    config.setdefault("adjustments", [])
    adjustment = "sqq track analyzes selected frames serially before linking"
    if adjustment not in config["adjustments"]:
        config["adjustments"].append(adjustment)

    # Track owns a temporary Analyze workspace.  Only the render fragments are
    # needed there; target-specific reports are published later by io.tracking.
    execution_config = deepcopy(config)
    execution_config.setdefault("output", {})["types"] = ["sqq-render"]
    plan = build_run_plan(
        paths,
        execution_config,
        source_root,
        topology=topology,
    )
    selected_frames = len(plan.tasks)
    if selected_frames < 1:
        raise ValueError("sqq track did not select any trajectory frames.")
    _require_nondecreasing_frame_times(plan, trajectory)

    planned_config = plan.context.config
    planned_graph = planned_config.get("graph", {})
    config.setdefault("graph", {}).update(
        {
            key: planned_graph[key]
            for key in (
                "effective_bond_mode",
                "effective_bond_mode_reason",
            )
            if key in planned_graph
        }
    )
    config.setdefault("input", {})["sampling"] = dict(plan.sampling)

    effective_modes = tuple(dict.fromkeys(plan.effective_graph_modes.values()))
    render_session = RenderSession.create(
        source_root,
        RenderSpec(
            atom_scope=str(config.get("render", {}).get("atom_scope", "full")),
            component_roles=planned_config,
            requested_graph_mode=str(config["graph"]["bond_mode"]),
            effective_graph_mode=(
                effective_modes[0] if len(effective_modes) == 1 else None
            ),
        ),
    )
    plan = replace(
        plan,
        context=replace(
            plan.context,
            strict=True,
            retain_results=False,
            stream_results=True,
            tracking_snapshots=True,
            fragment_dir=render_session.fragment_dir,
        ),
        policy=replace(
            plan.policy,
            backend="serial",
            workers=1,
            strict=True,
            in_flight_limit=1,
        ),
    )

    accumulator = TrackingAccumulator(tracking_config)
    tracking_sink = _TrackingAnalysisSink(accumulator)
    progress = RunProgressDisplay(
        total=selected_frames,
        total_started_at=started,
        include_cluster_stage=bool(
            config.get("hydrate_cluster", {}).get("enabled", False)
        ),
        cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        include_patch_stage=bool(
            config.get("half_cage", {}).get("enabled", False)
            or config.get("quasi_cage", {}).get("enabled", False)
        ),
    )
    try:
        outcomes = AnalysisRunner(
            plan,
            event_sink=_TrackProgressSink(progress),
            sinks=(tracking_sink,),
        ).run()
        if len(outcomes) != selected_frames or any(not item.ok for item in outcomes):
            raise RuntimeError(
                "Tracking cannot skip failed or missing trajectory frames."
            )
        result = accumulator.result()
        bundle = render_session.finalize(tracking=result)
    except Exception:
        render_session.abort()
        raise
    finally:
        progress.close()

    if not bundle.complete:
        raise RuntimeError("No SQQ render bundle was produced for tracking.")
    metadata = {
        "source": None,
        "input": str(trajectory.resolve()),
        "topology": str(topology.resolve()) if topology is not None else None,
        "input_format": config["input"]["format"],
        "selected_frames": selected_frames,
        "source_frames_total": int(plan.sampling.get("total_frames", selected_frames)),
        "native_frame_interval_ps": plan.sampling.get("native_frame_interval_ps"),
        "raw_frame_step": int(plan.sampling.get("raw_frame_step", 1)),
    }
    return result, bundle, plan, metadata


def _require_nondecreasing_frame_times(plan: RunPlan, trajectory: Path) -> None:
    """Fail before analysis when stored frame times cannot form a time axis.

    The tracking accumulator rejects a snapshot whose time precedes the
    previous one; checking the planned times here reports the offending frames
    immediately instead of after the prefix has been analyzed.
    """
    stored_times = plan.frame_times_ps
    if not stored_times:
        return
    if len(stored_times) != len(plan.tasks):
        raise RuntimeError(
            "Tracking frame-time metadata does not match the selected frame count."
        )
    input_config = plan.context.config.get("input", {})
    if not isinstance(input_config, Mapping):
        raise ValueError("input must be a mapping in sqq_config.yaml.")
    effective_times: list[float] = []
    synthesized: list[bool] = []
    for task, stored in zip(plan.tasks, stored_times):
        try:
            value = effective_frame_time_ps(
                stored,
                raw_frame_index=task.raw_frame_index,
                frame_index=task.frame_index,
                input_config=input_config,
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Tracking frame {task.display_name!r} has an invalid "
                f"effective time ({exc})."
            ) from exc
        effective_times.append(value)
        synthesized.append(stored is None)
    for position in range(1, len(effective_times)):
        previous, current = effective_times[position - 1], effective_times[position]
        if current < previous:
            previous_task = plan.tasks[position - 1]
            current_task = plan.tasks[position]
            previous_source = " (synthesized)" if synthesized[position - 1] else ""
            current_source = " (synthesized)" if synthesized[position] else ""
            raise ValueError(
                "Tracking requires nondecreasing frame times, but selected frame "
                f"{position} ({current_task.display_name}, t={current:g} ps"
                f"{current_source}) "
                f"precedes frame {position - 1} ({previous_task.display_name}, "
                f"t={previous:g} ps{previous_source}) in {trajectory}. "
                "Reorder the frames or fix "
                "their time metadata before tracking."
            )


def _native_interval(result: TrackingResult) -> float | None:
    values = [frame.time_ps for frame in result.frames]
    if len(values) < 2 or any(value is None for value in values):
        return None
    differences = np.diff(np.asarray(values, dtype=float))
    if len(differences) == 0 or np.any(differences <= 0):
        return None
    return (
        float(differences[0])
        if np.allclose(differences, differences[0], rtol=1.0e-6, atol=1.0e-6)
        else None
    )


def _run_info(
    config: Mapping[str, Any],
    args: Namespace,
    metadata: Mapping[str, Any],
    targets: Sequence[TargetSpec],
    result: TrackingResult | None,
    started: datetime,
    finished: datetime,
    elapsed: float,
    *,
    status: str,
    error: str = "",
) -> dict[str, Any]:
    frames = len(result.frames) if result is not None else 0
    graph_mode = str(config.get("graph", {}).get("bond_mode", ""))
    input_value = metadata.get("input") or metadata.get("source")
    output_metadata = config.get("output", {})
    if not isinstance(output_metadata, Mapping):
        output_metadata = {}
    requested_output_path = output_metadata.get("requested_path") or str(args.output)
    resolved_output_path = output_metadata.get("resolved_path") or str(
        Path(args.output).resolve()
    )
    output_auto_renamed = bool(output_metadata.get("auto_renamed", False))
    return {
        "sqq_version": __version__,
        "release_date": __release_date__,
        "status": status,
        "error": error,
        "command": "track",
        "started_at": started.isoformat(timespec="seconds"),
        "finished_at": finished.isoformat(timespec="seconds"),
        "date": started.strftime("%Y-%m-%d"),
        "start_time": started.strftime("%H:%M:%S"),
        "finish_time": finished.strftime("%H:%M:%S"),
        "time_zone": format_time_zone(started),
        "working_dir": str(Path.cwd()),
        "elapsed_seconds": round(float(elapsed), 3),
        "input": input_value,
        "topology": metadata.get("topology"),
        "input_format": metadata.get("input_format"),
        "source": metadata.get("source"),
        "matched_files": 1 if input_value else 0,
        "first_file": input_value,
        "output": str(Path(resolved_output_path).resolve()),
        "output_requested_path": str(requested_output_path),
        "output_resolved_path": str(Path(resolved_output_path).resolve()),
        "output_auto_renamed": output_auto_renamed,
        "engine_selector": config.get("mode", DEFAULT_MODE),
        "sqq_engine": engine_display(config.get("mode", DEFAULT_MODE)),
        "graph_mode": graph_mode,
        "effective_graph_modes": "",
        "graph_mode_display": graph_mode,
        "delta_time_ps": config.get("input", {}).get("delta_time_ps"),
        "native_frame_interval_ps": metadata.get("native_frame_interval_ps"),
        "raw_frame_step": metadata.get("raw_frame_step", 1),
        "selected_frames": metadata.get("selected_frames", frames),
        "source_frames_total": metadata.get("source_frames_total", frames),
        "frames_total": frames,
        "frames_ok": frames if status == "completed" else 0,
        "frames_failed": 0 if status == "completed" else 1,
        "failures": [] if not error else [{"frame": "", "source": "", "error": error}],
        "workers": 1,
        "parallel_backend": "serial",
        "math_threads": int(config.get("parallel", {}).get("math_threads", 1)),
        "target": ",".join(target.raw for target in targets),
        "track_count": len(result.tracks) if result is not None else 0,
        "gap_frame": result.config.gap_frame if result is not None else None,
        "config_file": getattr(args, "config", None) or "<built-in defaults>",
        "config_output": str((Path(args.output) / "sqq_config_resolved.yaml").resolve()),
    }


def _markdown_table_cell(value: object) -> str:
    return str(value).replace("|", r"\|")


def _write_root_track_info(
    output: Path,
    config: Mapping[str, Any],
    metadata: Mapping[str, Any],
    targets: Sequence[TargetSpec],
    result: TrackingResult,
    elapsed: float,
) -> None:
    root = output / "track"
    root.mkdir(parents=True, exist_ok=True)
    values = (
        ("SQQ version", __version__),
        ("release date", __release_date__),
        ("SQQ engine", engine_display(config.get("mode", DEFAULT_MODE))),
        ("input", metadata.get("input") or metadata.get("source") or ""),
        ("topology", metadata.get("topology") or ""),
        ("targets", ", ".join(target.raw for target in targets)),
        ("frames", len(result.frames)),
        ("tracks", len(result.tracks)),
        ("events", len(result.events)),
        ("gap frame", result.config.gap_frame),
        ("elapsed (s)", f"{elapsed:.3f}"),
    )
    lines = [
        "## Run metadata",
        "",
        "| item | value |",
        "| --- | --- |",
        *(
            f"| {_markdown_table_cell(item)} | {_markdown_table_cell(value)} |"
            for item, value in values
        ),
    ]
    append_track_info_section(root / "track_info.md", "\n".join(lines))
