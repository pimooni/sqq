from __future__ import annotations

"""Atomic writers for main and detail CSV report tables."""

import os
from pathlib import Path
import shutil
import tempfile
from time import perf_counter
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from .models import (
    LEGACY_SUMMARY_DETAIL_DIRECTORY,
    SUMMARY_DETAIL_TABLE_NAMES,
    SUMMARY_MAIN_TABLE_NAMES,
    ReportTable,
    table_metric,
)
from .transaction import commit_output_bundle, temporary_output_path


def write_summary_csvs(
    outdir: Path,
    tables: Sequence[ReportTable | tuple[str, pd.DataFrame, bool]],
    config: dict[str, Any],
    *,
    return_metrics: bool = False,
) -> dict[str, Any] | None:
    """Atomically write one main-summary CSV per workbook-equivalent table."""
    started = perf_counter()
    dir_name = (
        str(config.get("output", {}).get("summary_csv_dir", "summary")).strip()
        or "summary"
    )
    summary_dir = outdir / dir_name
    summary_dir.mkdir(parents=True, exist_ok=True)
    metrics: dict[str, Any] = {
        "enabled": True,
        "total_seconds": 0.0,
        "tables": [],
    }
    pending: list[tuple[Path, Path]] = []
    written_names: set[str] = set()
    try:
        for item in tables:
            report_table = (
                item
                if isinstance(item, ReportTable)
                else ReportTable(item[0], item[1], include_header=bool(item[2]))
            )
            name = report_table.name
            table = report_table.data
            target = summary_dir / f"{name}.csv"
            temp_path = temporary_output_path(target)
            pending.append((temp_path, target))
            write_started = perf_counter()
            table.to_csv(
                temp_path,
                index=False,
                header=report_table.include_header,
                encoding="utf-8-sig",
            )
            elapsed = perf_counter() - write_started
            written_names.add(name)
            metric = table_metric(name, table, write_seconds=elapsed)
            metric["bytes"] = temp_path.stat().st_size
            metrics["tables"].append(metric)

        stale_paths = [
            summary_dir / f"{name}.csv"
            for name in SUMMARY_MAIN_TABLE_NAMES
            if name not in written_names
        ]
        commit_output_bundle(pending, stale_paths)
    finally:
        for temp_path, _ in pending:
            temp_path.unlink(missing_ok=True)

    metrics["total_seconds"] = round(perf_counter() - started, 6)
    if return_metrics:
        return metrics
    return None


def remove_summary_csvs(outdir: Path, config: dict[str, Any]) -> None:
    """Remove known main-summary CSVs while preserving unrelated files."""
    dir_name = (
        str(config.get("output", {}).get("summary_csv_dir", "summary")).strip()
        or "summary"
    )
    summary_dir = outdir / dir_name
    if not summary_dir.exists():
        return
    commit_output_bundle(
        [],
        [summary_dir / f"{name}.csv" for name in SUMMARY_MAIN_TABLE_NAMES],
    )
    try:
        summary_dir.rmdir()
    except OSError:
        pass


def write_summary_detail_csvs(
    outdir: Path,
    tables: Sequence[ReportTable] | Mapping[str, pd.DataFrame],
    config: dict[str, Any],
    *,
    return_metrics: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, dict[str, Any]]:
    """Atomically write detail CSVs and optionally return their timing metadata."""
    started = perf_counter()
    detail_dir_name = (
        str(config.get("output", {}).get("summary_csv_dir", "summary")).strip()
        or "summary"
    )
    detail_dir = outdir / detail_dir_name
    detail_dir.mkdir(parents=True, exist_ok=True)
    _remove_legacy_summary_detail_csvs(outdir, detail_dir)
    rows: list[dict[str, Any]] = []
    metrics: dict[str, Any] = {"enabled": True, "total_seconds": 0.0, "tables": []}
    pending: list[tuple[Path, Path]] = []
    written_names: set[str] = set()
    try:
        tables_by_name = (
            {name: ReportTable(name, table, detail=True) for name, table in tables.items()}
            if isinstance(tables, Mapping)
            else {report_table.name: report_table for report_table in tables}
        )
        for name in SUMMARY_DETAIL_TABLE_NAMES:
            report_table = tables_by_name.get(name)
            if report_table is None:
                continue
            table = report_table.data
            path = detail_dir / f"{name}.csv"
            temp_path = temporary_output_path(path)
            pending.append((temp_path, path))
            write_started = perf_counter()
            table.to_csv(temp_path, index=False, encoding="utf-8-sig")
            elapsed = perf_counter() - write_started
            written_names.add(name)
            relative_dir = detail_dir_name.rstrip("/\\")
            relative_file = f"{relative_dir}/{name}.csv".replace("\\", "/")
            rows.append(
                {
                    "table": name,
                    "file": relative_file,
                    "rows": len(table),
                    "columns": len(table.columns),
                }
            )
            metric = table_metric(name, table, write_seconds=elapsed)
            metric["bytes"] = temp_path.stat().st_size
            metrics["tables"].append(metric)
        stale_paths = [
            detail_dir / f"{name}.csv"
            for name in SUMMARY_DETAIL_TABLE_NAMES
            if name not in written_names
        ]
        commit_output_bundle(pending, stale_paths)
    finally:
        for temp_path, _ in pending:
            temp_path.unlink(missing_ok=True)
    metrics["total_seconds"] = round(perf_counter() - started, 6)
    detail_index = pd.DataFrame(rows, columns=["table", "file", "rows", "columns"])
    if return_metrics:
        return detail_index, metrics
    return detail_index


def remove_summary_detail_csvs(outdir: Path, config: dict[str, Any]) -> None:
    """Remove known stale detail CSVs when summary-detail-csv is disabled."""
    detail_dir_name = (
        str(config.get("output", {}).get("summary_csv_dir", "summary")).strip()
        or "summary"
    )
    detail_dir = outdir / detail_dir_name
    _remove_legacy_summary_detail_csvs(outdir, detail_dir)
    if not detail_dir.exists():
        return
    commit_output_bundle(
        [],
        [detail_dir / f"{name}.csv" for name in SUMMARY_DETAIL_TABLE_NAMES],
    )
    try:
        detail_dir.rmdir()
    except OSError:
        pass


def _remove_legacy_summary_detail_csvs(outdir: Path, current_dir: Path) -> None:
    """Remove known SQQ detail files from the retired detail directory."""
    legacy_dir = outdir / LEGACY_SUMMARY_DETAIL_DIRECTORY
    if legacy_dir == current_dir or not legacy_dir.exists():
        return
    commit_output_bundle(
        [],
        [legacy_dir / f"{name}.csv" for name in SUMMARY_DETAIL_TABLE_NAMES],
    )
    try:
        legacy_dir.rmdir()
    except OSError:
        pass
