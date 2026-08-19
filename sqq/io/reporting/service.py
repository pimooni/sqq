"""Coordinate report selection, staging, atomic publication, and rollback."""

from __future__ import annotations

from copy import deepcopy
import os
from pathlib import Path
import shutil
import tempfile
from time import perf_counter
from typing import Any

import pandas as pd

from ... import __release_date__, __version__
from ...config import (
    DEFAULT_MODE,
    dump_config,
    is_cpp_mode,
    load_config,
    order_parameter_display,
    output_enabled,
    output_type_display,
)
from .csv_writer import (
    remove_summary_csvs,
    remove_summary_detail_csvs,
    write_summary_csvs,
    write_summary_detail_csvs,
)
from .excel_writer import ensure_excel_table_size, format_summary_workbook
from .models import (
    LEGACY_SUMMARY_DETAIL_DIRECTORY,
    SUMMARY_COLUMNS,
    SUMMARY_DETAIL_TABLE_NAMES,
    SUMMARY_MAIN_TABLE_NAMES,
    ReportTable,
    ReportTables,
    table_metric,
)
from .tables import stable_extra_columns, summary_detail_tables, summary_output_tables
from .transaction import commit_output_bundle, temporary_output_path


def write_summary(
    rows: list[dict[str, Any]],
    outdir: Path,
    config: dict[str, Any],
    write_xlsx: bool = True,
    run_info: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Publish config, XLSX, and CSV summaries as one output transaction."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    directory_names = summary_directory_names(root, config)
    staging_root = Path(
        tempfile.mkdtemp(prefix=".sqq-summary-stage-", dir=root)
    )
    try:
        metrics = _write_summary_staged(
            rows,
            staging_root,
            config,
            write_xlsx=write_xlsx,
            run_info=run_info,
        )
        completed_info = deepcopy(run_info or {})
        completed_info["summary_write"] = metrics
        write_run_config(staging_root, config, completed_info)
        pending = [
            (source, root / source.relative_to(staging_root))
            for source in staging_root.rglob("*")
            if source.is_file()
        ]
        removals = summary_generated_paths(root, directory_names)
        commit_output_bundle(pending, removals)
    finally:
        shutil.rmtree(staging_root, ignore_errors=True)
    remove_empty_summary_directories(root, directory_names)
    return metrics


def clear_previous_summary_outputs(outdir: Path, config: dict[str, Any]) -> None:
    """Remove one previous SQQ summary generation before a new run starts."""
    root = Path(outdir)
    if not root.exists():
        return
    directory_names = summary_directory_names(root, config)
    commit_output_bundle([], summary_generated_paths(root, directory_names))
    remove_empty_summary_directories(root, directory_names)


def summary_directory_names(outdir: Path, config: dict[str, Any]) -> tuple[str, ...]:
    """Return current, previous, default, and legacy summary directories."""
    names = {
        summary_directory_name(config),
        "summary",
        "summary_csv",
        LEGACY_SUMMARY_DETAIL_DIRECTORY,
    }
    previous_path = Path(outdir) / "sqq_config_resolved.yaml"
    if previous_path.is_file():
        try:
            previous = load_config(previous_path)
        except Exception:
            previous = None
        if isinstance(previous, dict):
            names.add(summary_directory_name(previous))
    return tuple(sorted(name for name in names if safe_summary_directory_name(name)))


def summary_directory_name(config: dict[str, Any]) -> str:
    """Return the configured relative summary directory."""
    return (
        str(config.get("output", {}).get("summary_csv_dir", "summary")).strip()
        or "summary"
    )


def safe_summary_directory_name(name: str) -> bool:
    """Restrict cleanup to relative paths contained by the result root."""
    path = Path(name)
    return not path.is_absolute() and ".." not in path.parts


def summary_generated_paths(
    outdir: Path,
    directory_names: tuple[str, ...],
) -> list[Path]:
    """Enumerate only known SQQ summary files."""
    paths = [
        outdir / "summary.xlsx",
        outdir / "summary.md",
        outdir / "sqq_config_resolved.yaml",
        outdir / "run_config.yaml",
    ]
    csv_names = tuple(
        dict.fromkeys(
            (*SUMMARY_MAIN_TABLE_NAMES, *SUMMARY_DETAIL_TABLE_NAMES)
        )
    )
    for directory_name in directory_names:
        directory = outdir / directory_name
        paths.extend(directory / f"{name}.csv" for name in csv_names)
    return list(dict.fromkeys(paths))


def remove_empty_summary_directories(
    outdir: Path,
    directory_names: tuple[str, ...],
) -> None:
    """Remove empty known summary directories without touching user files."""
    for directory_name in directory_names:
        directory = outdir / directory_name
        if directory == outdir:
            continue
        try:
            directory.rmdir()
        except OSError:
            pass


def _write_summary_staged(
    rows: list[dict[str, Any]],
    outdir: Path,
    config: dict[str, Any],
    write_xlsx: bool = True,
    run_info: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write global summaries and return deterministic output timing metadata."""
    started = perf_counter()
    outdir.mkdir(parents=True, exist_ok=True)
    metrics: dict[str, Any] = {
        "dataframe_seconds": 0.0,
        "config_initial_seconds": 0.0,
        "detail_table_build_seconds": 0.0,
        "detail_csv": {"enabled": False, "total_seconds": 0.0, "tables": []},
        "summary_csv": {"enabled": False, "total_seconds": 0.0, "tables": []},
        "xlsx": {"enabled": False, "total_seconds": 0.0, "sheets": []},
    }

    dataframe_started = perf_counter()
    columns = list(SUMMARY_COLUMNS)
    extra_columns = stable_extra_columns(rows, columns)
    public_rows = [
        {key: value for key, value in row.items() if not str(key).startswith("_")}
        for row in rows
    ]
    data = pd.DataFrame(public_rows, columns=columns + extra_columns)
    metrics["dataframe_seconds"] = round(perf_counter() - dataframe_started, 6)

    summary_md = outdir / "summary.md"
    if summary_md.exists():
        summary_md.unlink()

    config_started = perf_counter()
    write_run_config(outdir, config, run_info or {})
    metrics["config_initial_seconds"] = round(perf_counter() - config_started, 6)

    detail_index = pd.DataFrame()
    detail_report_tables: tuple[ReportTable, ...] = ()
    if (
        output_enabled(config, "summary-detail-csv")
        or output_enabled(config, "cluster-detail")
    ):
        detail_build_started = perf_counter()
        detail_tables = summary_detail_tables(data, config, raw_rows=rows)
        detail_report_tables = tuple(
            ReportTable(name, table, detail=True)
            for name, table in detail_tables.items()
        )
        metrics["detail_table_build_seconds"] = round(perf_counter() - detail_build_started, 6)
        detail_index, detail_metrics = write_summary_detail_csvs(
            outdir,
            detail_report_tables,
            config,
            return_metrics=True,
        )
        metrics["detail_csv"] = detail_metrics
    else:
        remove_summary_detail_csvs(outdir, config)

    write_summary_csv = output_enabled(config, "summary-csv")
    write_summary_xlsx = bool(write_xlsx) and output_enabled(config, "summary-xlsx")
    summary_tables: tuple[ReportTable, ...] = ()
    if write_summary_csv or write_summary_xlsx:
        summary_tables = summary_output_tables(
            data,
            config,
            run_info or {},
            detail_index,
        )
    report_tables = ReportTables(main=summary_tables, detail=detail_report_tables)

    if write_summary_csv:
        metrics["summary_csv"] = write_summary_csvs(
            outdir, report_tables.main, config, return_metrics=True
        )
    else:
        remove_summary_csvs(outdir, config)

    summary_xlsx = outdir / "summary.xlsx"
    if write_summary_xlsx:
        xlsx_started = perf_counter()
        xlsx_metrics: dict[str, Any] = {
            "enabled": True,
            "table_write_seconds": 0.0,
            "format_seconds": 0.0,
            "save_seconds": 0.0,
            "total_seconds": 0.0,
            "bytes": 0,
            "sheets": [],
        }
        temp_path = temporary_output_path(summary_xlsx)
        writer: pd.ExcelWriter | None = None
        try:
            writer = pd.ExcelWriter(temp_path, engine="openpyxl")
            tables = report_tables.main

            for report_table in tables:
                sheet_name = report_table.name
                table = report_table.data
                include_header = report_table.include_header
                ensure_excel_table_size(table, sheet_name, include_header=include_header)
                table_started = perf_counter()
                table.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                    header=include_header,
                )
                elapsed = perf_counter() - table_started
                xlsx_metrics["table_write_seconds"] += elapsed
                xlsx_metrics["sheets"].append(
                    table_metric(
                        sheet_name,
                        table,
                        write_seconds=elapsed,
                    )
                )

            format_started = perf_counter()
            format_metrics = format_summary_workbook(writer.book)
            xlsx_metrics["format_seconds"] = perf_counter() - format_started
            format_by_sheet = {item["sheet"]: item for item in format_metrics}
            for item in xlsx_metrics["sheets"]:
                item.update(format_by_sheet.get(item["sheet"], {}))

            save_started = perf_counter()
            writer.close()
            writer = None
            xlsx_metrics["save_seconds"] = perf_counter() - save_started
            os.replace(temp_path, summary_xlsx)
            xlsx_metrics["bytes"] = summary_xlsx.stat().st_size
        except Exception:
            if writer is not None:
                try:
                    writer.close()
                except Exception:
                    pass
            raise
        finally:
            temp_path.unlink(missing_ok=True)
        xlsx_metrics["table_write_seconds"] = round(xlsx_metrics["table_write_seconds"], 6)
        xlsx_metrics["format_seconds"] = round(xlsx_metrics["format_seconds"], 6)
        xlsx_metrics["save_seconds"] = round(xlsx_metrics["save_seconds"], 6)
        xlsx_metrics["total_seconds"] = round(perf_counter() - xlsx_started, 6)
        metrics["xlsx"] = xlsx_metrics
    else:
        summary_xlsx.unlink(missing_ok=True)

    metrics["total_seconds"] = round(perf_counter() - started, 6)
    return metrics


def write_run_config(
    outdir: Path,
    config: dict[str, Any],
    run_info: dict[str, Any],
) -> dict[str, Any]:
    """Atomically write mandatory run metadata and return the dumped mapping."""
    outdir.mkdir(parents=True, exist_ok=True)
    config_for_dump = config_with_run_metadata(config, run_info)
    target = outdir / "sqq_config_resolved.yaml"
    legacy_targets = (outdir / "run_config.yaml",)
    temp_path = temporary_output_path(target)
    try:
        with temp_path.open("w", encoding="utf-8", newline="\n") as handle:
            dump_config(config_for_dump, handle)
        os.replace(temp_path, target)
        for legacy_target in legacy_targets:
            legacy_target.unlink(missing_ok=True)
    finally:
        temp_path.unlink(missing_ok=True)
    return config_for_dump


def config_with_run_metadata(config: dict[str, Any], run_info: dict[str, Any]) -> dict[str, Any]:
    """Return a config copy that preserves raw settings and records resolved run metadata."""
    enriched = deepcopy(config)
    if not run_info:
        return enriched
    run = enriched.setdefault("run", {})
    parallel = config.get("parallel", {})
    run.update({
        "sqq_version": run_info.get("sqq_version", __version__),
        "release_date": run_info.get("release_date", __release_date__),
        "engine_selector": config.get("mode", DEFAULT_MODE),
        "engine": "sqq-cpp" if is_cpp_mode(config.get("mode", DEFAULT_MODE)) else "sqq-py",
        "profile": run_info.get("profile", config.get("run", {}).get("profile", "")),
        "resolution_adjustments": run_info.get("resolution_adjustments", ()),
        "config_output": run_info.get(
            "config_output",
            "sqq_config_resolved.yaml",
        ),
        "status": run_info.get("status", ""),
        "error": run_info.get("error", ""),
        "graph_mode_requested": run_info.get("graph_mode", config.get("graph", {}).get("bond_mode", "")),
        "graph_mode_effective": run_info.get("effective_graph_modes", ""),
        "graph_mode_reason": run_info.get("graph_mode_reason", ""),
        "graph_mode_display": run_info.get("graph_mode_display", ""),
        "graph_mode_by_group": run_info.get("graph_mode_by_group", {}),
        "graph_mode_reason_by_group": run_info.get("graph_mode_reason_by_group", {}),
        "order_parameters": run_info.get(
            "order_parameters",
            order_parameter_display(config.get("order", {}).get("parameters")),
        ),
        "find_half": run_info.get("find_half", "on" if config.get("half_cage", {}).get("enabled", False) else "off"),
        "find_quasi": run_info.get("find_quasi", "on" if config.get("quasi_cage", {}).get("enabled", False) else "off"),
        "find_cluster": run_info.get(
            "find_cluster",
            "on" if config.get("hydrate_cluster", {}).get("enabled", False) else "off",
        ),
        "output_types": run_info.get(
            "output_types",
            output_type_display(
                config.get("output", {}).get("types"),
                cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
            ),
        ),
        "output_requested_path": run_info.get("output_requested_path", ""),
        "output_resolved_path": run_info.get("output_resolved_path", ""),
        "output_auto_renamed": bool(run_info.get("output_auto_renamed", False)),
        "sampling_interval": run_info.get("sampling_interval", ""),
        "native_frame_interval_ps": run_info.get("native_frame_interval_ps"),
        "delta_time_ps": run_info.get("delta_time_ps"),
        "raw_frame_step": run_info.get("raw_frame_step", 1),
        "selected_frames": run_info.get("selected_frames", ""),
        "source_frames_total": run_info.get("source_frames_total", ""),
        "frames_total": run_info.get("frames_total", ""),
        "frames_ok": run_info.get("frames_ok", ""),
        "frames_failed": run_info.get("frames_failed", ""),
        "failures": run_info.get("failures", []),
        "worker_request": parallel.get("workers", "auto"),
        "worker_policy": run_info.get("worker_policy", ""),
        "workers_resolved": run_info.get("workers", ""),
        "parallel_backend": run_info.get("parallel_backend", "serial"),
        "math_threads_per_worker": run_info.get("math_threads", 1),
        "summary_write": run_info.get("summary_write", {}),
    })
    if run["graph_mode_by_group"] and not run["graph_mode_display"]:
        run.pop("graph_mode_display")
    for key in (
        "topology_group_count",
        "topology_group_limit",
        "topology_group_limit_exceeded",
        "topology_group_labels_enabled",
        "info_only_fallback_required",
        "topology_grouping",
        "topology_groups",
        "topology_source_mapping",
        "topology_group",
        "topology_fingerprint",
        "requested_output_types",
        "output_policy",
        "warnings",
    ):
        if key in run_info:
            run[key] = deepcopy(run_info[key])
    return enriched
