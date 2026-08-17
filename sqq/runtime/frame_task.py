"""Acquire, analyze, and transactionally publish one frame task."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

from ..config import DEFAULT_MODE, is_cpp_mode, output_enabled
from ..io.gro_writer import (
    write_cage_gro_files,
    write_half_cage_gro_files,
    write_hydrate_cluster_gro_files,
    write_ice_gro_file,
    write_quasi_cage_gro_files,
    write_ring_gro_files,
    write_water_order_gro_file,
)
from ..io.render.frame import write_sqq_cage_fragment
from ..io.render.models import FRAGMENT_DIRECTORY
from ..io.reporting import failed_row, result_row, write_frame_info, write_membership, write_order_parameter
from ..io.trajectory import read_frames
from ..models import Frame, FrameResult
from .contracts import FrameTask, RunContext, TaskOutcome, TaskStatus
from .frame import StageCallback, analyze_frame


def execute_frame_task(
    task: FrameTask,
    run_context: RunContext,
    *,
    stage_callback: StageCallback | None = None,
) -> TaskOutcome:
    """Execute one task with strict/non-strict error and output semantics."""
    config = run_context.config_for(task)
    output_root = run_context.output_root_for(task)
    fragment_dir = run_context.fragment_dir_for(task)
    frame: Frame | None = None
    try:
        _report_stage(stage_callback, "reading frame")
        frame = _obtain_frame(task, run_context, config)
        if task.output_name:
            frame.name = task.output_name
        if frame.time_ps is None:
            input_config = config.get("input", {})
            time_index = (
                task.raw_frame_index
                if task.raw_frame_index is not None
                else task.frame_index
            )
            frame.time_ps = float(input_config.get("first_file_time_ps", 0.0)) + (
                time_index * float(input_config.get("frame_time_step_ps", 1.0))
            )

        result = analyze_frame(frame, config, stage_callback=stage_callback)
        _publish_frame_result(
            task,
            run_context,
            config,
            result,
            output_root,
            fragment_dir,
            stage_callback,
        )
        return TaskOutcome(
            task_index=task.task_index,
            frame_index=task.frame_index,
            status=TaskStatus.OK,
            row=result_row(result),
            result=(
                result
                if run_context.retain_results or run_context.stream_results
                else None
            ),
        )
    except Exception as exc:
        if run_context.strict:
            raise
        frame_name = frame.name if frame is not None else task.display_name
        source = str(
            (frame.source if frame is not None else None)
            or task.source
            or run_context.trajectory
            or ""
        )
        row = failed_row(frame_name, source, str(exc))
        return TaskOutcome(
            task_index=task.task_index,
            frame_index=task.frame_index,
            status=TaskStatus.FAILED,
            row=row,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _obtain_frame(
    task: FrameTask,
    context: RunContext,
    config: Mapping[str, Any],
) -> Frame:
    if task.frame is not None:
        return task.frame
    source = task.source or context.trajectory
    if source is None:
        raise ValueError(f"{task.display_name} has no frame or input source.")
    indexes = None if task.raw_frame_index is None else [task.raw_frame_index]
    lammps_config = config.get("input", {}).get("lammps")
    frames = iter(
        read_frames(
            [Path(source)],
            topology=context.topology,
            xyz_scale=float(config.get("input", {}).get("xyz_scale", 0.1)),
            frame_indexes=indexes,
            lammps_config=lammps_config,
        )
    )
    try:
        return next(frames)
    except StopIteration as exc:
        raise ValueError(f"No frame was read for {task.display_name}.") from exc
    finally:
        close = getattr(frames, "close", None)
        if callable(close):
            close()


def _publish_frame_result(
    task: FrameTask,
    context: RunContext,
    config: Mapping[str, Any],
    result: FrameResult,
    output_root: Path,
    fragment_dir: Path | None,
    stage_callback: StageCallback | None,
) -> None:
    if task.separated_output:
        report_dir = output_root / "info"
        frame_dir = output_root / "gro" / result.frame.name
    else:
        report_dir = output_root / result.frame.name
        frame_dir = report_dir

    _report_stage(stage_callback, "writing outputs")
    output_root.mkdir(parents=True, exist_ok=True)
    staging_root = Path(tempfile.mkdtemp(prefix=".sqq-frame-stage-", dir=output_root))
    staged_frame_dir = staging_root / "frame"
    staged_report_dir = staging_root / "info" if report_dir != frame_dir else staged_frame_dir
    fragment = None
    try:
        write_frame_outputs(result, staged_frame_dir, config, report_dir=staged_report_dir)
        if output_enabled(config, "sqq-render"):
            fragment = write_sqq_cage_fragment(
                result,
                fragment_dir or output_root / FRAGMENT_DIRECTORY,
                task.frame_index,
                requested_graph_mode=config["graph"]["bond_mode"],
                atom_scope=config.get("render", {}).get("atom_scope", "full"),
                component_config=config,
            )
        commit_staged_frame_outputs(
            result.frame.name,
            staged_report_dir,
            staged_frame_dir,
            report_dir,
            frame_dir,
            staging_root,
        )
    except Exception:
        if fragment is not None:
            fragment.gro_path.unlink(missing_ok=True)
            fragment.manifest_path.unlink(missing_ok=True)
        raise
    finally:
        shutil.rmtree(staging_root, ignore_errors=True)
    _report_stage(stage_callback, "done")


def write_frame_outputs(
    result: FrameResult,
    frame_dir: Path,
    config: Mapping[str, Any],
    *,
    report_dir: Path | None = None,
) -> None:
    """Write configured per-frame files into a staging directory."""
    output = config.get("output", {})
    report_dir = report_dir or frame_dir
    order_parameters = config.get("order", {}).get("parameters", ["f3", "f4"])
    cpp_mode = is_cpp_mode(config.get("mode"))
    if output_enabled(config, "info"):
        write_frame_info(
            result,
            report_dir,
            ring_sizes=list(result.ring_report_sizes),
            requested_bond_mode=config["graph"]["bond_mode"],
            order_parameters=order_parameters,
            analysis_mode=config.get("mode", DEFAULT_MODE),
            input_metadata=frame_input_metadata(config),
        )
    else:
        _remove_optional_info_output(result, report_dir)

    if not cpp_mode and output_enabled(config, "membership-tsv"):
        write_membership(result, report_dir)
    else:
        _remove_optional_tsv_outputs(result, report_dir, remove_membership=True, remove_order=False)
    if not cpp_mode and output_enabled(config, "order-tsv"):
        write_order_parameter(result, report_dir, order_parameters=order_parameters)
    else:
        _remove_optional_tsv_outputs(result, report_dir, remove_membership=False, remove_order=True)
    _remove_optional_vmd_output(result, report_dir)

    layout = str(output.get("structure_layout", "grouped"))
    write_empty = bool(output.get("write_empty_files", False))
    center_resname = str(output.get("center_resname", "CNT"))
    for parameter in ("f3", "f4"):
        _remove_water_order_gro_output(result, frame_dir, parameter, layout)
        if output_enabled(config, f"{parameter}-gro"):
            write_water_order_gro_file(
                result, frame_dir, parameter, write_empty=write_empty, layout=layout
            )
    _remove_generated_gro_outputs(result, frame_dir, "ring-gro", layout)
    if not cpp_mode and output_enabled(config, "ring-gro"):
        write_ring_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            sizes=set(result.ring_report_sizes),
            center_resname=center_resname,
        )
    _remove_generated_gro_outputs(result, frame_dir, "half-gro", layout)
    if not cpp_mode and output_enabled(config, "half-gro"):
        write_half_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            center_resname=center_resname,
        )
    _remove_generated_gro_outputs(result, frame_dir, "quasi-gro", layout)
    if not cpp_mode and output_enabled(config, "quasi-gro"):
        write_quasi_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            center_resname=center_resname,
        )
    _remove_generated_gro_outputs(result, frame_dir, "cage-gro", layout)
    if output_enabled(config, "cage-gro"):
        write_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            include_centers=not cpp_mode,
            center_resname=center_resname,
        )
    _remove_generated_gro_outputs(result, frame_dir, "cluster-gro", layout)
    if not cpp_mode and result.hydrate_cluster_enabled and output_enabled(config, "cluster-gro"):
        write_hydrate_cluster_gro_files(
            result, frame_dir, write_empty=write_empty, layout=layout
        )
    _remove_generated_gro_outputs(result, frame_dir, "ice-gro", layout)
    if not cpp_mode and output_enabled(config, "ice-gro"):
        write_ice_gro_file(result, frame_dir, write_empty=write_empty, layout=layout)
    _remove_frame_directory_if_empty(frame_dir)
    if report_dir != frame_dir:
        _remove_frame_directory_if_empty(report_dir)
        _remove_frame_directory_if_empty(report_dir.parent)
        _remove_frame_directory_if_empty(frame_dir.parent)


def commit_staged_frame_outputs(
    frame_name: str,
    staged_report_dir: Path,
    staged_frame_dir: Path,
    report_dir: Path,
    frame_dir: Path,
    staging_root: Path,
) -> None:
    """Commit one complete frame output set and restore it on failure."""
    root_pairs: list[tuple[Path, Path]] = []
    for source, target in ((staged_report_dir, report_dir), (staged_frame_dir, frame_dir)):
        if not any(existing == target for _, existing in root_pairs):
            root_pairs.append((source, target))

    backup_root = staging_root / "backup"
    backups: list[tuple[Path, Path]] = []
    installed: list[Path] = []
    try:
        for root_index, (_, target_root) in enumerate(root_pairs):
            if not target_root.is_dir():
                continue
            for existing in sorted(
                (path for path in target_root.rglob("*") if path.is_file()),
                key=lambda path: path.as_posix(),
            ):
                if not is_generated_frame_output(existing, frame_name):
                    continue
                backup = backup_root / str(root_index) / existing.relative_to(target_root)
                backup.parent.mkdir(parents=True, exist_ok=True)
                os.replace(existing, backup)
                backups.append((backup, existing))

        for source_root, target_root in root_pairs:
            if not source_root.is_dir():
                continue
            for source in sorted(
                (path for path in source_root.rglob("*") if path.is_file()),
                key=lambda path: path.as_posix(),
            ):
                target = target_root / source.relative_to(source_root)
                target.parent.mkdir(parents=True, exist_ok=True)
                os.replace(source, target)
                installed.append(target)
    except Exception:
        for target in reversed(installed):
            target.unlink(missing_ok=True)
        for backup, target in reversed(backups):
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(backup, target)
        raise
    finally:
        _cleanup_empty_output_directories(report_dir, frame_dir)


def is_generated_frame_output(path: Path, frame_name: str) -> bool:
    name = path.name.casefold()
    return name.startswith(f"{frame_name}_".casefold()) and name.endswith(
        (".md", ".tsv", ".gro", ".vmd.tcl")
    )


def frame_input_metadata(config: Mapping[str, Any]) -> dict[str, Any]:
    input_config = config.get("input", {})
    sampling = input_config.get("sampling", {})
    selected = sampling.get("selected_frames", 0)
    total = sampling.get("total_frames", 0)
    delta = input_config.get("delta_time_ps")
    native = sampling.get("native_frame_interval_ps")
    metadata: dict[str, Any] = {
        "input_format": input_config.get("format", ""),
        "topology": input_config.get("topology"),
        "sampling_interval": (
            f"all [native {native:g} ps; analyzing {selected} of {total} frames]"
            if delta is None and native is not None
            else f"{delta:g} ps [native {native:g} ps; analyzing {selected} of {total} frames]"
            if delta is not None and native is not None
            else ""
        ),
        "native_frame_interval_ps": native,
        "delta_time_ps": delta,
        "raw_frame_step": sampling.get("raw_frame_step"),
        "selected_frames": f"{selected} / {total}",
        "find_half": _on_off(config.get("half_cage", {}).get("enabled", False)),
        "find_quasi": _on_off(config.get("quasi_cage", {}).get("enabled", False)),
    }
    if str(metadata["input_format"]).startswith("lammps-"):
        lammps = input_config.get("lammps", {})
        metadata.update(
            {
                "lammps_units": lammps.get("units", "real"),
                "lammps_timestep": lammps.get("timestep", 1.0),
                "lammps_atom_style": lammps.get("atom_style", "full"),
                "lammps_type_map_source": lammps.get("type_map_source", "<configuration>"),
            }
        )
    return metadata


def _remove_optional_info_output(result: FrameResult, frame_dir: Path) -> None:
    (frame_dir / f"{result.frame.name}_info.md").unlink(missing_ok=True)


def _remove_optional_tsv_outputs(
    result: FrameResult,
    frame_dir: Path,
    *,
    remove_membership: bool,
    remove_order: bool,
) -> None:
    suffixes: list[str] = []
    if remove_membership:
        suffixes.append("membership")
    if remove_order:
        suffixes.extend(("f3f4", "order_parameter"))
    for suffix in suffixes:
        (frame_dir / f"{result.frame.name}_{suffix}.tsv").unlink(missing_ok=True)


def _remove_optional_vmd_output(result: FrameResult, frame_dir: Path) -> None:
    (frame_dir / f"{result.frame.name}_view.vmd.tcl").unlink(missing_ok=True)


def _remove_water_order_gro_output(
    result: FrameResult,
    frame_dir: Path,
    parameter: str,
    layout: str,
) -> None:
    filename = f"{result.frame.name}_{parameter}.gro"
    path = frame_dir / filename if layout == "flat" else frame_dir / "order" / filename
    path.unlink(missing_ok=True)
    if layout == "grouped":
        _remove_frame_directory_if_empty(path.parent)


def _remove_generated_gro_outputs(
    result: FrameResult,
    frame_dir: Path,
    output_type: str,
    layout: str,
) -> None:
    grouped = {
        "ring-gro": "ring",
        "half-gro": "half_cage",
        "quasi-gro": "quasi_cage",
        "cage-gro": "cage",
        "ice-gro": "ice",
        "cluster-gro": "hydrate_cluster",
    }
    patterns = {
        "ring-gro": f"{result.frame.name}_ring_*.gro",
        "half-gro": f"{result.frame.name}_hc_*.gro",
        "quasi-gro": f"{result.frame.name}_qc_*.gro",
        "cage-gro": f"{result.frame.name}_cage_*.gro",
        "ice-gro": f"{result.frame.name}_ice*.gro",
        "cluster-gro": f"{result.frame.name}_cluster_*.gro",
    }
    if layout not in {"grouped", "flat"}:
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    root = frame_dir / grouped[output_type]
    if root.exists():
        parents: set[Path] = set()
        for path in root.rglob("*.gro"):
            if not path.name.startswith(f"{result.frame.name}_"):
                continue
            parent = path.parent
            while parent == root or root in parent.parents:
                parents.add(parent)
                if parent == root:
                    break
                parent = parent.parent
            path.unlink(missing_ok=True)
        for directory in sorted(parents, key=lambda item: len(item.parts), reverse=True):
            _remove_frame_directory_if_empty(directory)
    for path in frame_dir.glob(patterns[output_type]):
        path.unlink(missing_ok=True)


def _cleanup_empty_output_directories(report_dir: Path, frame_dir: Path) -> None:
    for root in {report_dir, frame_dir}:
        if root.is_dir():
            for directory in sorted(
                (path for path in root.rglob("*") if path.is_dir()),
                key=lambda path: len(path.parts),
                reverse=True,
            ):
                _remove_frame_directory_if_empty(directory)
        _remove_frame_directory_if_empty(root)
    if frame_dir.parent.name == "gro":
        _remove_frame_directory_if_empty(frame_dir.parent)


def _remove_frame_directory_if_empty(path: Path) -> None:
    try:
        path.rmdir()
    except OSError:
        pass


def _on_off(value: Any) -> str:
    if isinstance(value, str):
        return "on" if value.strip().lower() in {"on", "true", "yes", "1"} else "off"
    return "on" if bool(value) else "off"


def _report_stage(callback: Callable[[str], None] | None, stage: str) -> None:
    if callback is not None:
        callback(stage)


__all__ = [
    "commit_staged_frame_outputs",
    "execute_frame_task",
    "frame_input_metadata",
    "is_generated_frame_output",
    "write_frame_outputs",
]
