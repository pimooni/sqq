from __future__ import annotations

"""Cross-frame cage tracking workflow for ``sqq track``."""

from argparse import Namespace
from collections import defaultdict
from copy import deepcopy
import csv
from dataclasses import replace
from datetime import datetime
import os
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml

from .. import __release_date__, __version__
from ..config import (
    DEFAULT_MODE,
    apply_cli_overrides,
    engine_display,
    is_cpp_mode,
    load_config,
    normalize_engine_capabilities,
    normalize_analysis_scopes,
    refresh_resolution_report,
    normalize_mode,
    validate_cpp_cli,
)
from ..core.tracking import (
    TrackingAccumulator,
    parse_targets,
    select_targets,
    snapshot_from_frame_result,
)
from ..core.geometry import pbc_aware_centroid
from ..io.lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    inspect_lammps_topology_mapping,
)
from ..io.render import (
    RenderSession,
    RenderSpec,
    RenderBundle,
    TRACK_MEMBERSHIP_NAME,
    TRACK_RENDER_DIRECTORY,
    discover_sqq_cage_bundle,
    publish_target_render_bundle,
    validate_tracking_source_bundle,
)
from ..io.render.frame import visualization_atoms
from ..io.reporting import write_run_config
from ..io.tracking import (
    add_precursor_membership,
    discover_track_state,
    read_tracking_result,
    target_directory_name,
    write_track_outputs,
)
from ..io.trajectory import expand_inputs, read_frames, trajectory_frame_selection
from ..models import FrameResult
from ..models.tracking import (
    TargetSelection,
    TargetSpec,
    TrackingConfig,
    TrackingResult,
)
from ..runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ..runtime.frame import analyze_frame
from ..runtime.output_lock import output_lock
from ..ui.diagnostics import RunDiagnostics, capture_run_warnings
from ..ui.formatting import format_time_zone
from ..ui.progress import RunProgressDisplay
from ..ui.run_header import input_format_label, print_run_banner
from ..ui import completed_run_statistics, print_final_results, refresh_terminal
from .analyze_plan import build_run_plan
from .session import AnalysisEvent, AnalysisRunner, AnalysisSink


__all__ = ["track"]

_TRACK_ID = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)
_TRACKING_FIELDS = {
    "min_jaccard",
    "min_shared_fraction",
    "min_shared_waters",
    "max_center_distance_nm",
    "gap_frame",
    "guest_tiebreak",
}
_TRACK_RUNTIME_FIELDS = {"target", "source"}
_RAW_TRACK_FIELD_ALIASES = {"min_shared_water": "min_shared_waters"}
_PRECURSOR_STATE_FIELDS = (
    "status",
    "reason",
    "track_id",
    "frame_index",
    "frame",
    "time_ps",
    "state",
    "target_water_count",
    "present_water_count",
    "bond_count",
    "component_count",
    "largest_component",
    "ring_count",
    "half_cage_count",
    "quasi_cage_count",
    "cage_count",
)
_WATER_HISTORY_FIELDS = (
    "status",
    "reason",
    "track_id",
    "frame_index",
    "frame",
    "time_ps",
    "state",
    "water_atomid",
    "present",
    "atom_index",
    "resid",
    "x_nm",
    "y_nm",
    "z_nm",
    "target_degree",
)
_PrecursorData = tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[int, tuple[tuple[int, ...], tuple[float, float, float]]],
]


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
        if outcome.result is None:
            raise RuntimeError(
                f"Tracking frame {task.display_name!r} did not retain its analysis result."
            )
        self.accumulator.add(
            snapshot_from_frame_result(outcome.result, int(task.frame_index))
        )
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
) -> dict[str, bool]:
    """Describe work that this Track invocation actually executed."""
    ran = int(frames) > 0
    features = {
        "water_network": False,
        "ring_topology": False,
        "cage_topology": False,
        "half_cage": False,
        "quasi_cage": False,
        "cage_isomer": False,
        "cage_occupancy": False,
        "hydrate_phase_domain": False,
        "vmd_rendering": ran,
        "cage_tracking": ran,
        "cage_lifetime": ran,
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
            "cage_occupancy": True,
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
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    with output_lock(output):
        diagnostics = RunDiagnostics()
        with capture_run_warnings(diagnostics):
            try:
                _track_locked(args, diagnostics)
            finally:
                diagnostics.emit()


def _track_locked(args: Namespace, diagnostics: RunDiagnostics) -> None:
    started_wall = datetime.now().astimezone()
    started = perf_counter()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
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
            config["mode"] = source_engine
            normalize_engine_capabilities(config, emit_warnings=False)
    _prepare_target_capabilities(config, targets)
    validate_cpp_cli(args, config)
    normalize_analysis_scopes(config)
    refresh_resolution_report(config)
    requested_tracking = _tracking_config(config)

    try:
        with TemporaryDirectory(
            prefix=".sqq-track-",
            dir=output,
            ignore_cleanup_errors=True,
        ) as temporary:
            if getattr(args, "input", None):
                result, source_bundle, retained_results, metadata = _track_input(
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
                retained_results = None
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
            precursor_data = (
                _raw_precursor_histories(selections, args, config)
                if getattr(args, "input", None)
                else None
            )
            analysis_completed_at = perf_counter()
            validate_tracking_source_bundle(
                source_bundle,
                frame_count=len(result.frames),
            )
            written = write_track_outputs(
                result,
                output,
                targets=[target.raw for target in targets],
            )
            for selection in selections:
                directory = written[target_directory_name(selection.target)]
                publish_target_render_bundle(
                    selection,
                    directory,
                    source_bundle,
                )

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
        raw_mode = bool(getattr(args, "input", None))
        statistics["executed_features"] = _track_executed_features(
            config,
            frames=len(result.frames),
            raw_mode=raw_mode,
        )
        if not raw_mode:
            statistics["executed_order_parameters"] = ()
        completed_outputs = tuple(statistics.get("completed_outputs", ()))
        if "sqq-render" not in completed_outputs:
            statistics["completed_outputs"] = completed_outputs + ("sqq-render",)
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
    list[FrameResult] | None,
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
    return result, bundle, None, metadata


def _prepare_target_capabilities(
    config: dict[str, Any],
    targets: Sequence[TargetSpec],
) -> None:
    if not any(target.kind == "phase" for target in targets):
        return
    if is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        raise ValueError(
            "Phase targets require SQQ-Py hydrate-cluster classification; "
            "use -e py or -e 00."
        )
    config.setdefault("hydrate_cluster", {})["enabled"] = True


def _update_track_config(
    config: dict[str, Any],
    args: Namespace,
    targets: Sequence[TargetSpec],
) -> None:
    section = config.setdefault("track", {})
    if not isinstance(section, dict):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    section["target"] = ",".join(target.raw for target in targets)
    section["source"] = (
        str(Path(args.source).resolve())
        if getattr(args, "source", None)
        else None
    )
    for name in _TRACKING_FIELDS:
        value = getattr(args, name, None)
        if value is not None:
            section[name] = value


def _resolve_track_request(
    config: Mapping[str, Any],
    args: Namespace,
) -> tuple[str | Iterable[str], str | Path | None]:
    """Resolve CLI-over-YAML target/source values before mutating config."""
    section = config.get("track", {})
    if section is None:
        section = {}
    if not isinstance(section, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    cli_target = getattr(args, "target", None)
    target = (
        section.get("target", "all")
        if cli_target is None or cli_target == ""
        else cli_target
    )
    cli_source = getattr(args, "source", None)
    configured_source = section.get("source")
    source = cli_source if cli_source not in (None, "") else configured_source
    if source in (None, ""):
        source = None
    return target, source  # type: ignore[return-value]


def _source_engine_selector(state_path: Path) -> str | None:
    """Read the Analyze engine recorded beside an imported Track state."""
    candidates = (
        state_path.parent.parent / "sqq_config_resolved.yaml",
        state_path.parent / "sqq_config_resolved.yaml",
    )
    config_path = next((path for path in candidates if path.is_file()), None)
    if config_path is None:
        return None
    try:
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ValueError(
            f"Cannot read the Analyze engine from {config_path}: {exc}"
        ) from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"Invalid resolved SQQ configuration: {config_path}")
    raw = payload.get("engine", payload.get("mode"))
    if raw in (None, ""):
        run = payload.get("run", {})
        if isinstance(run, Mapping):
            raw = run.get("engine_selector")
    if raw in (None, ""):
        return None
    return normalize_mode(raw)


def _tracking_config(config: Mapping[str, Any]) -> TrackingConfig:
    values = config.get("track", {})
    if not isinstance(values, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    unknown = sorted(set(values).difference(_TRACKING_FIELDS | _TRACK_RUNTIME_FIELDS))
    if unknown:
        raise ValueError(
            "Unsupported track configuration field(s): "
            + ", ".join(unknown)
            + ". Remove fields that have no implemented effect."
        )
    return TrackingConfig.from_mapping(
        {name: values[name] for name in _TRACKING_FIELDS if name in values}
    )


def _explicit_tracking_fields(args: Namespace) -> set[str]:
    """Return fields explicitly supplied by CLI or the user's YAML file."""
    fields = {
        name
        for name in _TRACKING_FIELDS
        if getattr(args, name, None) is not None
    }
    config_path = getattr(args, "config", None)
    if not config_path:
        return fields
    try:
        payload = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ValueError(f"Cannot inspect Track settings in {config_path}: {exc}") from exc
    if payload is None:
        return fields
    if not isinstance(payload, Mapping):
        raise ValueError("The configuration root must be a mapping.")
    section = payload.get("track", {})
    if section is None:
        return fields
    if not isinstance(section, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    for raw_name in section:
        normalized = _RAW_TRACK_FIELD_ALIASES.get(str(raw_name), str(raw_name))
        if normalized in _TRACKING_FIELDS:
            fields.add(normalized)
    return fields


def _prepare_lammps_metadata(
    config: dict[str, Any],
    trajectory: Path,
    topology: Path | None,
    args: Namespace,
) -> None:
    if trajectory.suffix.lower() not in LAMMPS_TRAJECTORY_SUFFIXES:
        return
    lammps = config["input"]["lammps"]
    explicit = bool(lammps.get("type_map"))
    resolved: Mapping[str, Any] = {}
    rebuilt = False
    if topology is not None:
        resolved, rebuilt = inspect_lammps_topology_mapping(topology, lammps)
    mapping = ", ".join(
        (
            f"{type_id}=ignore"
            if entry.ignore
            else f"{type_id}={entry.resname}/{entry.atomname}"
        )
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    )
    lammps["resolved_type_map"] = {
        type_id: (
            {"ignore": True}
            if entry.ignore
            else {"resname": entry.resname, "atomname": entry.atomname}
        )
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    }
    if explicit:
        source = (
            str(Path(args.config).resolve())
            if getattr(args, "config", None)
            else "<configuration>"
        )
    else:
        source = "auto (DATA topology)"
        if rebuilt:
            source += "; molecule IDs rebuilt from Bonds"
    lammps["type_map_source"] = f"{source}: {mapping}" if mapping else source


def _validate_requested_tracks(
    result: TrackingResult,
    targets: Sequence[TargetSpec],
) -> None:
    known = {track.track_id for track in result.tracks}
    missing = [
        target.value
        for target in targets
        if target.kind == "track" and target.value not in known
    ]
    if missing:
        preview = ", ".join(sorted(known, key=_track_sort_key)[:12]) or "none"
        raise ValueError(
            "Unknown persistent cage ID(s): "
            + ", ".join(missing)
            + f". Available IDs begin with: {preview}."
        )


def _track_sort_key(value: str) -> int:
    match = _TRACK_ID.fullmatch(value)
    return int(match.group(1)) if match is not None else 2**63 - 1


def _write_precursor_outputs(
    selections: Sequence[TargetSelection],
    written: Mapping[str, Path],
    histories: Mapping[str, _PrecursorData] | None,
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
                    target_dir
                    / TRACK_RENDER_DIRECTORY
                    / TRACK_MEMBERSHIP_NAME,
                    selection.target.value,
                    render_frames,
                )
        _write_csv(target_dir / "precursor_state.csv", precursor, _PRECURSOR_STATE_FIELDS)
        _write_csv(target_dir / "water_history.csv", water, _WATER_HISTORY_FIELDS)


def _raw_precursor_histories(
    selections: Sequence[TargetSelection],
    args: Namespace,
    config: dict[str, Any],
) -> dict[str, _PrecursorData]:
    """Reanalyze only the required prefix for persistent-ID precursors."""
    targets: dict[str, tuple[int, frozenset[int]]] = {}
    output: dict[str, _PrecursorData] = {}
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

    input_path = Path(args.input)
    pattern = getattr(args, "pattern", None) or config["input"]["pattern"]
    recursive = bool(
        getattr(args, "recursive", False) or config["input"]["recursive"]
    )
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)
    if len(paths) != 1:
        raise ValueError("Persistent-ID precursor tracing requires one raw trajectory.")
    trajectory = paths[0]
    topology = Path(args.topology) if getattr(args, "topology", None) else None
    frame_selection = trajectory_frame_selection(
        trajectory,
        topology,
        delta_time_ps=config["input"].get("delta_time_ps"),
        lammps_config=config["input"].get("lammps", {}),
    )
    last_required = max(item[0] for item in targets.values())
    raw_indexes = frame_selection.raw_indexes[: last_required + 1]
    progress = RunProgressDisplay(
        total=len(raw_indexes),
        total_started_at=perf_counter(),
        include_cluster_stage=bool(
            config.get("hydrate_cluster", {}).get("enabled", False)
        ),
        cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        include_patch_stage=bool(
            config.get("half_cage", {}).get("enabled", False)
            or config.get("quasi_cage", {}).get("enabled", False)
        ),
    )
    analyzed_count = 0
    try:
        frames = read_frames(
            paths,
            topology=topology,
            xyz_scale=float(config["input"].get("xyz_scale", 0.1)),
            frame_indexes=raw_indexes,
            lammps_config=config["input"].get("lammps", {}),
        )
        for frame_index, frame in enumerate(frames):
            if frame.time_ps is None:
                raw_index = raw_indexes[frame_index]
                frame.time_ps = (
                    float(config["input"]["first_file_time_ps"])
                    + raw_index * float(config["input"]["frame_time_step_ps"])
                )
            callback = progress.start_frame(frame_index, frame.name)
            try:
                analyzed = analyze_frame(
                    frame,
                    config,
                    stage_callback=callback,
                )
                for track_id, (birth_index, water_ids) in targets.items():
                    if frame_index > birth_index:
                        continue
                    state, water_rows = _classify_precursor_frame(
                        analyzed, frame_index, track_id, water_ids
                    )
                    states, waters, render_frames = output[track_id]
                    states.append(state)
                    waters.extend(water_rows)
                    if frame_index < birth_index:
                        render_frames[frame_index] = _precursor_render_frame(
                            analyzed,
                            water_ids,
                            atom_scope=config.get("render", {}).get(
                                "atom_scope", "full"
                            ),
                        )
                analyzed_count += 1
            except Exception as exc:
                progress.complete_frame(False)
                raise RuntimeError(
                    f"Precursor tracing failed for frame {frame.name!r}: {exc}"
                ) from exc
            progress.complete_frame(True)
    finally:
        progress.close()
    if analyzed_count != len(raw_indexes):
        raise RuntimeError(
            "Trajectory reader returned "
            f"{analyzed_count} precursor frames; expected {len(raw_indexes)}."
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
        center_angstrom[0],
        center_angstrom[1],
        center_angstrom[2],
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
        ring
        for values in result.rings.values()
        for ring in values
        if set(ring.nodes).issubset(target_indexes)
    ]
    half = [
        patch
        for patch in result.half_cages
        if set(patch.waters).issubset(target_indexes)
    ]
    quasi = [
        patch
        for patch in result.quasi_cages
        if set(patch.waters).issubset(target_indexes)
    ]
    cages = [
        cage
        for cage in (result.all_cages or result.cages)
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
        "status": "available",
        "reason": "",
        "track_id": track_id,
        "frame_index": frame_index,
        "frame": result.frame.name,
        "time_ps": result.frame.time_ps,
        "state": state_name,
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
    suffix.append("")
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(suffix))


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
        "output": str(Path(args.output).resolve()),
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
        "# SQQ Track Run",
        "",
        "| item | value |",
        "| --- | --- |",
        *(
            f"| {_markdown_table_cell(item)} | {_markdown_table_cell(value)} |"
            for item, value in values
        ),
        "",
    ]
    (root / "track_info.md").write_text(
        "\n".join(lines), encoding="utf-8", newline="\n"
    )
