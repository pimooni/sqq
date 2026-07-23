from __future__ import annotations

"""Topology-aware orchestration for multiple independent GRO snapshots."""

from argparse import Namespace
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import get_context
from pathlib import Path
from queue import Empty
import sys
from time import perf_counter
from typing import Any

from .config import is_cpp_mode, output_enabled
from .core.graph import resolve_bond_mode
from .core.selection import select_waters
from .io.gro_grouping import GroGroupingResult, gro_topology_descriptor, scan_and_group_gro_inputs
from .io.summary import (
    failed_row,
    remove_summary_csvs,
    remove_summary_detail_csvs,
    write_run_config,
    write_summary,
)
from .io.trajectory import read_frames, read_gro
from .io.vmd import (
    cleanup_sqq_cage_bundle,
    finalize_sqq_cage_bundle,
    prepare_sqq_cage_fragments,
)
from .parallel import initialize_file_worker, limited_math_threads, process_file_task


@dataclass(frozen=True)
class StandaloneFileTask:
    """One file task with global progress and group-local frame identity."""

    global_index: int
    local_index: int
    path: Path
    group_key: int
    outdir: Path
    output_name: str


def is_multi_gro_batch(paths: list[Path]) -> bool:
    """Return whether inputs are multiple independent GRO snapshots."""
    return len(paths) > 1 and all(path.suffix.lower() == ".gro" for path in paths)


def analyze_multi_gro_batch(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    paths: list[Path],
    outdir: Path,
    topology: Path | None,
    run_started_at: datetime,
    started_at: float,
) -> None:
    """Analyze heterogeneous GRO inputs in topology-compatible result roots."""
    from . import pipeline as pipeline_api

    strict = bool(args.strict)
    grouping = scan_and_group_gro_inputs(paths, strict=strict)
    validate_shared_gro_topology(topology, grouping)

    requested_output_types = list(config.get("output", {}).get("types", []))
    execution_config = deepcopy(config)
    warnings: list[str] = []
    if grouping.info_only_fallback_required:
        warning = (
            f"Detected {grouping.group_count} independent GRO topology groups, above "
            f"the supported A-Z limit ({grouping.group_limit}). All readable files "
            "will be analyzed with info-only output; summaries, GRO output, and VMD "
            "output are disabled for this run."
        )
        warnings.append(warning)
        print(f"Warning: {warning}", file=sys.stderr)
        execution_config["output"]["types"] = ["info"]

    cleanup_previous_multi_gro_outputs(outdir, execution_config)

    group_configs: dict[int, dict[str, Any]] = {}
    group_outdirs: dict[int, Path] = {}
    group_names: dict[int, dict[int, str]] = {}
    global_names = unique_assignment_output_names(grouping.assignments)
    for group in grouping.groups:
        group_outdir = (
            outdir
            if grouping.info_only_fallback_required or grouping.group_count == 1
            else outdir / f"result_{group.label}"
        )
        group_outdir.mkdir(parents=True, exist_ok=True)
        group_config = deepcopy(execution_config)
        group_config["graph"]["effective_bond_mode"] = resolve_group_graph_mode(
            group.paths[0], group_config
        )
        group_configs[group.group_index] = group_config
        group_outdirs[group.group_index] = group_outdir
        group_names[group.group_index] = (
            global_names
            if grouping.info_only_fallback_required
            else unique_assignment_output_names(group.inputs)
        )

    tasks = [
        StandaloneFileTask(
            global_index=assignment.source_index,
            local_index=local_index,
            path=assignment.path,
            group_key=group.group_index,
            outdir=group_outdirs[group.group_index],
            output_name=group_names[group.group_index][assignment.source_index],
        )
        for group in grouping.groups
        for local_index, assignment in enumerate(group.inputs)
    ]
    tasks.sort(key=lambda task: task.global_index)

    parallel_backend = pipeline_api.normalize_parallel_backend(
        execution_config.get("parallel", {}).get("backend", "process")
    )
    requested_workers = pipeline_api.resolve_workers(
        execution_config["parallel"].get("workers"),
        max(1, len(tasks)),
        mode=execution_config.get("mode", "50"),
        backend=parallel_backend,
    )
    workers = requested_workers if len(tasks) > 1 and parallel_backend != "serial" else 1
    active_backend = parallel_backend if workers > 1 else "serial"

    pipeline_api.print_run_header(
        args,
        execution_config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        active_backend,
        run_started_at,
    )
    pipeline_api.print_terminal_field("Topology groups", grouping.group_count)
    pipeline_api.print_terminal_field("Grouping policy", topology_grouping_policy(grouping))
    print("")

    initial_info = pipeline_api.build_run_info(
        args,
        execution_config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        active_backend,
        0.0,
        run_started_at,
        run_started_at,
        [],
    )
    add_topology_grouping_metadata(
        initial_info,
        grouping,
        group_outdirs,
        group_configs,
        requested_output_types,
        warnings,
    )
    initial_info["status"] = "running"
    initial_info["error"] = ""
    write_run_config(outdir, execution_config, initial_info)

    if grouping.group_count > 1 and not grouping.info_only_fallback_required:
        for group in grouping.groups:
            group_info = build_group_run_info(
                pipeline_api,
                args,
                input_path,
                topology,
                group,
                group_configs[group.group_index],
                group_outdirs[group.group_index],
                workers,
                active_backend,
                run_started_at,
                run_started_at,
                0.0,
                [],
                grouping,
                requested_output_types,
                warnings,
            )
            group_info["status"] = "running"
            group_info["error"] = ""
            write_run_config(
                group_outdirs[group.group_index],
                group_configs[group.group_index],
                group_info,
            )

    bundle_groups: list[int] = []
    if not grouping.info_only_fallback_required:
        for group in grouping.groups:
            group_config = group_configs[group.group_index]
            group_outdir = group_outdirs[group.group_index]
            bundle_gro = output_enabled(group_config, "sqq-cage-gro")
            bundle_script = output_enabled(group_config, "sqq-render")
            cleanup_sqq_cage_bundle(group_outdir)
            if bundle_gro or bundle_script:
                prepare_sqq_cage_fragments(group_outdir)
                bundle_groups.append(group.group_index)
    else:
        cleanup_sqq_cage_bundle(outdir)

    rows_by_index: dict[int, dict[str, Any]] = {
        failure.source_index: failed_row(
            failure.path.stem, str(failure.path), failure.error
        )
        for failure in grouping.failures
    }
    try:
        rows_by_index.update(
            analyze_grouped_gro_tasks(
                tasks,
                execution_config,
                group_configs,
                outdir,
                workers=workers,
                backend=active_backend,
                strict=strict,
                total_started_at=started_at,
            )
        )
        for group_index in bundle_groups:
            group_config = group_configs[group_index]
            finalize_sqq_cage_bundle(
                group_outdirs[group_index],
                write_gro=output_enabled(group_config, "sqq-cage-gro"),
                write_script=output_enabled(group_config, "sqq-render"),
            )
    except Exception as exc:
        for group_outdir in set(group_outdirs.values()):
            cleanup_sqq_cage_bundle(group_outdir)
        write_failed_manifests(
            pipeline_api,
            args,
            execution_config,
            input_path,
            paths,
            outdir,
            topology,
            workers,
            active_backend,
            run_started_at,
            started_at,
            rows_by_index,
            grouping,
            group_outdirs,
            group_configs,
            requested_output_types,
            warnings,
            exc,
        )
        raise

    finished_at = datetime.now().astimezone()
    elapsed_seconds = perf_counter() - started_at
    all_rows = [rows_by_index[index] for index in sorted(rows_by_index)]
    root_info = pipeline_api.build_run_info(
        args,
        execution_config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        active_backend,
        elapsed_seconds,
        run_started_at,
        finished_at,
        all_rows,
    )
    add_topology_grouping_metadata(
        root_info,
        grouping,
        group_outdirs,
        group_configs,
        requested_output_types,
        warnings,
    )
    root_info["status"] = "completed"
    root_info["error"] = ""

    try:
        if grouping.info_only_fallback_required or not grouping.groups:
            write_run_config(outdir, execution_config, root_info)
        elif grouping.group_count == 1:
            group = grouping.groups[0]
            group_config = group_configs[group.group_index]
            root_info["topology_group"] = "single"
            root_info["topology_fingerprint"] = group.fingerprint
            root_info["summary_write"] = write_summary(
                all_rows,
                outdir,
                group_config,
                write_xlsx=output_enabled(group_config, "summary-xlsx"),
                run_info=root_info,
            )
            write_run_config(outdir, group_config, root_info)
        else:
            summary_metrics: dict[str, Any] = {}
            for group in grouping.groups:
                group_rows = [
                    rows_by_index[index]
                    for index in group.source_indices
                    if index in rows_by_index
                ]
                group_info = build_group_run_info(
                    pipeline_api,
                    args,
                    input_path,
                    topology,
                    group,
                    group_configs[group.group_index],
                    group_outdirs[group.group_index],
                    workers,
                    active_backend,
                    run_started_at,
                    finished_at,
                    elapsed_seconds,
                    group_rows,
                    grouping,
                    requested_output_types,
                    warnings,
                )
                group_info["status"] = "completed"
                group_info["error"] = ""
                metrics = write_summary(
                    group_rows,
                    group_outdirs[group.group_index],
                    group_configs[group.group_index],
                    write_xlsx=output_enabled(group_configs[group.group_index], "summary-xlsx"),
                    run_info=group_info,
                )
                group_info["summary_write"] = metrics
                write_run_config(
                    group_outdirs[group.group_index],
                    group_configs[group.group_index],
                    group_info,
                )
                summary_metrics[str(group.label)] = metrics
            root_info["summary_write"] = {"groups": summary_metrics}
            write_run_config(outdir, execution_config, root_info)
    except Exception as exc:
        write_failed_manifests(
            pipeline_api,
            args,
            execution_config,
            input_path,
            paths,
            outdir,
            topology,
            workers,
            active_backend,
            run_started_at,
            started_at,
            rows_by_index,
            grouping,
            group_outdirs,
            group_configs,
            requested_output_types,
            warnings,
            exc,
        )
        raise
    pipeline_api.print_run_summary(root_info)
    print(f"Wrote SQQ results: {outdir}")


def analyze_grouped_gro_tasks(
    tasks: list[StandaloneFileTask],
    base_config: dict[str, Any],
    group_configs: dict[int, dict[str, Any]],
    outdir: Path,
    *,
    workers: int,
    backend: str,
    strict: bool,
    total_started_at: float,
) -> dict[int, dict[str, Any]]:
    """Analyze every topology group through one shared scheduling pool."""
    if not tasks:
        return {}
    if workers <= 1 or backend == "serial":
        return analyze_grouped_gro_serial(tasks, group_configs, strict, total_started_at)
    if backend == "thread":
        return analyze_grouped_gro_threaded(
            tasks, group_configs, workers, strict, total_started_at
        )
    if backend != "process":
        raise ValueError("Parallel analysis requires backend=process or backend=thread.")
    return analyze_grouped_gro_processes(
        tasks,
        base_config,
        group_configs,
        outdir,
        workers,
        strict,
        total_started_at,
    )


def analyze_grouped_gro_serial(
    tasks: list[StandaloneFileTask],
    group_configs: dict[int, dict[str, Any]],
    strict: bool,
    total_started_at: float,
) -> dict[int, dict[str, Any]]:
    """Analyze grouped GRO tasks serially while retaining local frame indexes."""
    from . import pipeline as pipeline_api

    rows: dict[int, dict[str, Any]] = {}
    progress = pipeline_api.RunProgressDisplay(
        total=len(tasks),
        total_started_at=total_started_at,
        include_cluster_stage=any(
            bool(item.get("hydrate_cluster", {}).get("enabled", False))
            for item in group_configs.values()
        ),
        cpp_mode=any(is_cpp_mode(item.get("mode")) for item in group_configs.values()),
    )
    try:
        for display_index, task in enumerate(tasks):
            callback = progress.start_frame(display_index, task.output_name)
            _, row = process_thread_task(
                task,
                group_configs[task.group_key],
                strict,
                callback=callback,
            )
            rows[task.global_index] = row
            progress.complete_frame(row.get("status") == "ok")
    finally:
        progress.close()
    return rows


def analyze_grouped_gro_threaded(
    tasks: list[StandaloneFileTask],
    group_configs: dict[int, dict[str, Any]],
    workers: int,
    strict: bool,
    total_started_at: float,
) -> dict[int, dict[str, Any]]:
    """Analyze grouped GRO tasks in one shared thread pool."""
    from . import pipeline as pipeline_api

    rows: dict[int, dict[str, Any]] = {}
    progress = pipeline_api.ParallelRunProgressDisplay(
        total=len(tasks),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=any(
            bool(item.get("hydrate_cluster", {}).get("enabled", False))
            for item in group_configs.values()
        ),
        cpp_mode=any(is_cpp_mode(item.get("mode")) for item in group_configs.values()),
    )
    display_by_global = {task.global_index: index for index, task in enumerate(tasks)}
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            iterator = iter(tasks)
            futures: dict[Any, StandaloneFileTask] = {}

            def fill_queue() -> None:
                while len(futures) < pipeline_api.process_in_flight_limit(workers):
                    try:
                        task = next(iterator)
                    except StopIteration:
                        return
                    display_index = display_by_global[task.global_index]
                    future = executor.submit(
                        process_thread_task,
                        task,
                        group_configs[task.group_key],
                        strict,
                        None,
                        progress,
                        display_index,
                    )
                    futures[future] = task

            fill_queue()
            while futures:
                done, _ = wait(set(futures), return_when=FIRST_COMPLETED)
                for future in done:
                    task = futures.pop(future)
                    display_index = display_by_global[task.global_index]
                    try:
                        global_index, row = future.result()
                    except Exception:
                        progress.complete_file(display_index, False)
                        for queued in futures:
                            queued.cancel()
                        raise
                    rows[global_index] = row
                    progress.complete_file(display_index, row.get("status") == "ok")
                fill_queue()
    finally:
        progress.close()
    return rows


def analyze_grouped_gro_processes(
    tasks: list[StandaloneFileTask],
    base_config: dict[str, Any],
    group_configs: dict[int, dict[str, Any]],
    outdir: Path,
    workers: int,
    strict: bool,
    total_started_at: float,
) -> dict[int, dict[str, Any]]:
    """Analyze grouped GRO tasks in one shared spawned-process pool."""
    from . import pipeline as pipeline_api

    rows: dict[int, dict[str, Any]] = {}
    progress = pipeline_api.ParallelRunProgressDisplay(
        total=len(tasks),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=any(
            bool(item.get("hydrate_cluster", {}).get("enabled", False))
            for item in group_configs.values()
        ),
        cpp_mode=any(is_cpp_mode(item.get("mode")) for item in group_configs.values()),
    )
    display_by_global = {task.global_index: index for index, task in enumerate(tasks)}
    context = get_context("spawn")
    stage_queue = context.Queue()
    math_threads = int(base_config.get("parallel", {}).get("math_threads", 1))
    try:
        with limited_math_threads(math_threads):
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=context,
                initializer=initialize_file_worker,
                initargs=(base_config, str(outdir), strict, stage_queue, group_configs),
            ) as executor:
                iterator = iter(tasks)
                futures: dict[Any, StandaloneFileTask] = {}

                def fill_queue() -> None:
                    while len(futures) < pipeline_api.process_in_flight_limit(workers):
                        try:
                            task = next(iterator)
                        except StopIteration:
                            return
                        display_index = display_by_global[task.global_index]
                        future = executor.submit(
                            process_file_task,
                            display_index,
                            str(task.path),
                            task.group_key,
                            str(task.outdir),
                            task.local_index,
                            task.output_name,
                            True,
                        )
                        futures[future] = task

                fill_queue()
                while futures:
                    drain_group_stage_events(stage_queue, progress)
                    done, _ = wait(
                        set(futures), timeout=0.1, return_when=FIRST_COMPLETED
                    )
                    for future in done:
                        task = futures.pop(future)
                        display_index = display_by_global[task.global_index]
                        try:
                            _, row = future.result()
                        except Exception:
                            progress.complete_file(display_index, False)
                            for queued in futures:
                                queued.cancel()
                            raise
                        rows[task.global_index] = row
                        progress.complete_file(display_index, row.get("status") == "ok")
                    fill_queue()
                drain_group_stage_events(stage_queue, progress)
    finally:
        progress.close()
        stage_queue.close()
        stage_queue.join_thread()
    return rows


def drain_group_stage_events(stage_queue: Any, progress: Any) -> None:
    """Apply worker events emitted with compact display indexes."""
    while True:
        try:
            kind, frame_index, value, timestamp = stage_queue.get_nowait()
        except Empty:
            return
        if kind == "start":
            progress.start_file(frame_index, value, started_at=timestamp)
        elif kind == "stage":
            progress.update_stage(frame_index, value, started_at=timestamp)
        elif kind == "complete":
            progress.complete_file(frame_index, value == "ok")


def process_thread_task(
    task: StandaloneFileTask,
    config: dict[str, Any],
    strict: bool,
    callback: Any = None,
    progress: Any = None,
    display_index: int | None = None,
) -> tuple[int, dict[str, Any]]:
    """Read and analyze one grouped GRO task in the current process."""
    from .pipeline import process_frame

    if callback is None and progress is not None and display_index is not None:
        callback = progress.start_file(display_index, task.output_name)
    try:
        frame = next(
            iter(
                read_frames(
                    [task.path],
                    xyz_scale=float(config.get("input", {}).get("xyz_scale", 0.1)),
                )
            )
        )
        frame.name = task.output_name
        row = process_frame(
            task.local_index,
            frame,
            config,
            task.outdir,
            strict=strict,
            stage_callback=callback,
            separated_output=True,
        )
    except Exception as exc:
        if strict:
            raise
        row = failed_row(task.output_name, str(task.path), str(exc))
    return task.global_index, row


def cleanup_previous_multi_gro_outputs(
    outdir: Path,
    config: dict[str, Any],
) -> None:
    """Remove known stale SQQ files before selecting a new group layout."""
    roots = [
        outdir / f"result_{label}"
        for label in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        if (outdir / f"result_{label}").is_dir()
    ]
    for root in [*roots, outdir]:
        cleanup_generated_output_root(root, config)
    for root in roots:
        remove_empty_tree(root)


def cleanup_generated_output_root(root: Path, config: dict[str, Any]) -> None:
    """Remove known generated files while preserving unrelated user files."""
    if not root.exists():
        return
    cleanup_sqq_cage_bundle(root)
    for name in ("summary.xlsx", "summary.md", "run_config.yaml"):
        (root / name).unlink(missing_ok=True)
    remove_summary_csvs(root, config)
    legacy_config = deepcopy(config)
    legacy_config.setdefault("output", {})["summary_csv_dir"] = "summary_csv"
    remove_summary_csvs(root, legacy_config)
    remove_summary_detail_csvs(root, config)

    report_suffixes = (
        "_info.md",
        "_membership.tsv",
        "_order_parameter.tsv",
        "_f3f4.tsv",
        "_view.vmd.tcl",
    )
    gro_markers = ("_ring_", "_hc_", "_qc_", "_cage_", "_ice", "_cluster_")
    for candidate in root.rglob("*"):
        if not candidate.is_file():
            continue
        name = candidate.name.casefold()
        if name.endswith(report_suffixes):
            candidate.unlink(missing_ok=True)
        elif name.endswith(".gro") and any(marker in name for marker in gro_markers):
            candidate.unlink(missing_ok=True)
    remove_empty_tree(root, keep_root=True)


def remove_empty_tree(root: Path, *, keep_root: bool = False) -> None:
    """Remove empty descendants and optionally the empty root itself."""
    if not root.exists():
        return
    directories = sorted(
        (path for path in root.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for directory in directories:
        try:
            directory.rmdir()
        except OSError:
            pass
    if not keep_root:
        try:
            root.rmdir()
        except OSError:
            pass


def write_failed_manifests(
    pipeline_api: Any,
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    paths: list[Path],
    outdir: Path,
    topology: Path | None,
    workers: int,
    active_backend: str,
    run_started_at: datetime,
    started_at: float,
    rows_by_index: dict[int, dict[str, Any]],
    grouping: GroGroupingResult,
    group_outdirs: dict[int, Path],
    group_configs: dict[int, dict[str, Any]],
    requested_output_types: list[str],
    warnings: list[str],
    error: Exception,
) -> None:
    """Best-effort finalization of root and group manifests after failure."""
    failed_at = datetime.now().astimezone()
    elapsed = perf_counter() - started_at
    rows = [rows_by_index[index] for index in sorted(rows_by_index)]
    failed_info = pipeline_api.build_run_info(
        args,
        config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        active_backend,
        elapsed,
        run_started_at,
        failed_at,
        rows,
    )
    add_topology_grouping_metadata(
        failed_info,
        grouping,
        group_outdirs,
        group_configs,
        requested_output_types,
        warnings,
    )
    failed_info["status"] = "failed"
    failed_info["error"] = str(error)
    try:
        write_run_config(outdir, config, failed_info)
    except Exception:
        pass

    if grouping.group_count <= 1 or grouping.info_only_fallback_required:
        return
    for group in grouping.groups:
        group_rows = [
            rows_by_index[index]
            for index in group.source_indices
            if index in rows_by_index
        ]
        group_info = build_group_run_info(
            pipeline_api,
            args,
            input_path,
            topology,
            group,
            group_configs[group.group_index],
            group_outdirs[group.group_index],
            workers,
            active_backend,
            run_started_at,
            failed_at,
            elapsed,
            group_rows,
            grouping,
            requested_output_types,
            warnings,
        )
        group_info["status"] = "failed"
        group_info["error"] = str(error)
        try:
            write_run_config(
                group_outdirs[group.group_index],
                group_configs[group.group_index],
                group_info,
            )
        except Exception:
            pass


def validate_shared_gro_topology(
    topology: Path | None,
    grouping: GroGroupingResult,
) -> None:
    """Require one shared GRO topology to match every readable source."""
    if topology is None:
        return
    if topology.suffix.lower() != ".gro":
        raise ValueError(
            "Multiple independent GRO inputs require a GRO file for --top/-t "
            "so SQQ can validate every source topology."
        )
    descriptor = gro_topology_descriptor(read_gro(topology))
    mismatches = [
        str(assignment.path)
        for assignment in grouping.assignments
        if assignment.descriptor != descriptor
    ]
    if mismatches:
        joined = "\n  - ".join(mismatches)
        raise ValueError(
            f"Shared topology {topology} does not match these GRO source(s):\n  - {joined}"
        )


def resolve_group_graph_mode(path: Path, config: dict[str, Any]) -> str:
    """Resolve auto graph mode once for a topology-compatible GRO group."""
    requested = str(config["graph"]["bond_mode"])
    if requested != "auto":
        return resolve_bond_mode(requested, [], config["graph"].get("pair_file"))
    frame = read_gro(path)
    waters = select_waters(
        frame.atoms,
        resnames=set(config["water"]["resnames"]),
        oxygen_names=set(config["water"]["oxygen_names"]),
        hydrogen_names=set(config["water"]["hydrogen_names"]),
    )
    return resolve_bond_mode("auto", waters, config["graph"].get("pair_file"))


def unique_assignment_output_names(assignments: Any) -> dict[int, str]:
    """Return collision-free deterministic names for one output root."""
    ordered = sorted(assignments, key=lambda item: item.source_index)
    counts = Counter(item.path.stem.casefold() for item in ordered)
    seen: Counter[str] = Counter()
    used: set[str] = set()
    names: dict[int, str] = {}
    for item in ordered:
        key = item.path.stem.casefold()
        seen[key] += 1
        preferred = (
            item.path.stem
            if counts[key] == 1
            else f"{item.path.stem}_{seen[key]:03d}"
        )
        candidate = preferred
        disambiguator = 2
        while candidate.casefold() in used:
            candidate = f"{preferred}_{disambiguator:03d}"
            disambiguator += 1
        used.add(candidate.casefold())
        names[item.source_index] = candidate
    return names


def topology_grouping_policy(grouping: GroGroupingResult) -> str:
    """Return the public label for the selected grouping policy."""
    if grouping.info_only_fallback_required:
        return "info-only (>26 topology groups)"
    if grouping.group_count <= 1:
        return "single result root"
    return "separate result_A ... result_Z roots"


def topology_group_records(
    grouping: GroGroupingResult,
    group_outdirs: dict[int, Path],
    group_configs: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build YAML-safe group records for the root run manifest."""
    return [
        {
            "group_index": group.group_index,
            "group_label": group.label,
            "fingerprint": group.fingerprint,
            "source_count": len(group.inputs),
            "output_dir": str(group_outdirs[group.group_index].resolve()),
            "effective_graph_mode": group_configs[group.group_index]["graph"].get(
                "effective_bond_mode",
                group_configs[group.group_index]["graph"]["bond_mode"],
            ),
            "sources": [str(path.resolve()) for path in group.paths],
        }
        for group in grouping.groups
    ]


def add_topology_grouping_metadata(
    run_info: dict[str, Any],
    grouping: GroGroupingResult,
    group_outdirs: dict[int, Path],
    group_configs: dict[int, dict[str, Any]],
    requested_output_types: list[str],
    warnings: list[str],
) -> None:
    """Attach the complete multi-GRO grouping contract to run metadata."""
    run_info.update(grouping.limit_metadata())
    run_info["topology_grouping"] = topology_grouping_policy(grouping)
    run_info["topology_groups"] = topology_group_records(
        grouping, group_outdirs, group_configs
    )
    run_info["topology_source_mapping"] = list(grouping.source_mapping())
    run_info["requested_output_types"] = list(requested_output_types)
    run_info["output_policy"] = (
        "forced info-only"
        if grouping.info_only_fallback_required
        else "mode/configured outputs"
    )
    run_info["warnings"] = list(warnings)


def build_group_run_info(
    pipeline_api: Any,
    args: Namespace,
    input_path: Path,
    topology: Path | None,
    group: Any,
    config: dict[str, Any],
    outdir: Path,
    workers: int,
    active_backend: str,
    started_at: datetime,
    finished_at: datetime,
    elapsed_seconds: float,
    rows: list[dict[str, Any]],
    grouping: GroGroupingResult,
    requested_output_types: list[str],
    warnings: list[str],
) -> dict[str, Any]:
    """Build one topology group's run metadata and source mapping."""
    info = pipeline_api.build_run_info(
        args,
        config,
        input_path,
        outdir,
        list(group.paths),
        topology,
        workers,
        active_backend,
        elapsed_seconds,
        started_at,
        finished_at,
        rows,
    )
    info.update(grouping.limit_metadata())
    info["topology_grouping"] = topology_grouping_policy(grouping)
    info["topology_group"] = group.label
    info["topology_fingerprint"] = group.fingerprint
    info["topology_source_mapping"] = [
        record
        for record in grouping.source_mapping()
        if record.get("group_index") == group.group_index
    ]
    info["requested_output_types"] = list(requested_output_types)
    info["output_policy"] = "mode/configured outputs"
    info["warnings"] = list(warnings)
    return info
