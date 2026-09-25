"""Build top-level dashboard, workbook, and Markdown summary tables."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .... import __version__
from ....config import (
    DEFAULT_MODE,
    is_cpp_mode,
    normalize_order_parameters,
    order_parameter_display,
    output_type_display,
    q_degrees_from_order_parameters,
)
from ....presentation import SQQ_AUTHOR, SQQ_TITLE
from ....presentation.citation import (
    build_citation_recommendation,
    completed_citation_evidence,
)
from ....presentation.graph_mode import graph_mode_display
from ..models import ReportTable
from .details import (
    cage_isomer_summary_table,
    cage_occupancy_summary_table,
    cage_summary_table,
    connection_sheet_name,
    connection_summary_table,
    hydrate_cluster_summary_table,
    molecule_summary_table,
    order_parameter_summary_table,
    patch_summary_table,
    summary_simple_table,
)
from .formatting import (
    configured_ring_report_sizes,
    dashboard_cage_targets,
    excel_scalar,
    frames_failed_count,
    frames_ok_count,
    has_selected_guests,
    markdown_table,
    min_mean_max_column,
    sum_numeric_column,
)


def summary_output_tables(
    data: pd.DataFrame,
    config: dict[str, Any],
    run_info: dict[str, Any],
    detail_index: pd.DataFrame,
) -> tuple[ReportTable, ...]:
    """Build the shared table list for summary XLSX and CSV outputs."""
    tables: list[ReportTable] = [
        ReportTable(
            "summary",
            summary_dashboard_table(data, run_info, config),
            include_header=False,
        ),
    ]
    tables.extend(
        ReportTable(sheet_name, table)
        for sheet_name, table in summary_sheet_tables(data, config).items()
    )
    if not detail_index.empty:
        tables.append(ReportTable("detail_index", detail_index))
    return tuple(tables)


def summary_dashboard_table(data: pd.DataFrame, run_info: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    """Build a compact human-facing dashboard for the summary sheet."""
    cpp_mode = is_cpp_mode(config.get("mode", DEFAULT_MODE))
    title = SQQ_TITLE
    author = SQQ_AUTHOR
    matched_files = run_info.get("matched_files", "")
    try:
        matched_count = int(matched_files)
    except (TypeError, ValueError):
        matched_count = 0
    graph_mode_value = run_info.get("graph_mode_display") or graph_mode_display(
        run_info.get("graph_mode", config.get("graph", {}).get("bond_mode", "auto")),
        data.get("connection_mode", []),
    )
    selected_order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters", ["f3", "f4"])
    )
    selected_order_set = set(selected_order_parameters)
    q_degrees = q_degrees_from_order_parameters(selected_order_parameters)

    rows: list[list[Any]] = [
        [title, ""],
        [author, ""],
        ["", ""],
        ["Basic Information", ""],
        ["Date", run_info.get("date", "")],
        ["Start time", run_info.get("start_time", "")],
        ["Finish time", run_info.get("finish_time", "")],
        ["Time zone", run_info.get("time_zone", "")],
        ["Duration (s)", run_info.get("elapsed_seconds", "")],
        ["Working directory", run_info.get("working_dir", "")],
        ["Input", run_info.get("input", "")],
        ["Matched files", matched_files],
        ["Input format", run_info.get("input_format", "")],
    ]
    if matched_count > 1:
        rows.extend([
            ["First file", run_info.get("first_file", "")],
            ["Last file", run_info.get("last_file", "")],
        ])
    else:
        rows.append(["Current file", run_info.get("first_file", "")])
    rows.extend([
        ["Output directory", run_info.get("output_dir", "")],
        ["summary.xlsx", run_info.get("summary_xlsx", "")],
        ["summary_csv", run_info.get("summary_csv", "")],
        ["summary_detail_csv", run_info.get("summary_detail_csv", "")],
        ["sqq_config_resolved.yaml", run_info.get("config_output", "")],
        ["", ""],
        ["Configuration", ""],
        ["SQQ version", run_info.get("sqq_version", __version__)],
        ["SQQ engine", "sqq-cpp" if cpp_mode else "sqq-py"],
        ["Engine selector", run_info.get("engine_selector", config.get("mode", DEFAULT_MODE))],
        ["Profile", run_info.get("profile", config.get("run", {}).get("profile", ""))],
        ["Config file", run_info.get("config_file", "<built-in defaults>")],
        ["Topology", run_info.get("topology", "<none>")],
        ["Graph mode", graph_mode_value],
        ["Graph mode reason", run_info.get("graph_mode_reason", "")],
        ["Search sizes", excel_scalar(config.get("ring", {}).get("sizes", ""))],
    ])
    for adjustment in run_info.get("resolution_adjustments", ()):
        if not isinstance(adjustment, dict):
            continue
        parameter = adjustment.get("parameter", "parameter")
        effective = adjustment.get("effective")
        reason = adjustment.get("reason", "automatic adjustment")
        value = str(parameter)
        if effective is not None:
            value += f" -> {excel_scalar(effective)}"
        value += f" [{reason}]"
        rows.append(["Adjustment", value])
    if not cpp_mode:
        rows.append(["Ring report sizes", excel_scalar(configured_ring_report_sizes(config))])
    rows.append(["Ring definition", config.get("ring", {}).get("definition", "chordless")])
    if run_info.get("sampling_interval"):
        native = run_info.get("native_frame_interval_ps")
        delta = run_info.get("delta_time_ps")
        rows.extend([
            ["Native frame interval", "unknown" if native is None else f"{float(native):g} ps"],
            ["Delta time", "all" if delta is None else f"{float(delta):g} ps"],
            ["Raw frame step", run_info.get("raw_frame_step", 1)],
            ["Selected frames", f"{run_info.get('selected_frames', 0)} / {run_info.get('source_frames_total', 0)}"],
        ])
    if str(run_info.get("input_format", "")).startswith("lammps-"):
        rows.extend([
            ["LAMMPS units", run_info.get("lammps_units", "")],
            ["LAMMPS timestep", run_info.get("lammps_timestep", "")],
            ["LAMMPS atom style", run_info.get("lammps_atom_style", "")],
            ["LAMMPS type map", run_info.get("lammps_type_map_source", "")],
        ])
    if not cpp_mode:
        rows.extend([
            ["Find half", "on" if config.get("half_cage", {}).get("enabled", False) else "off"],
            ["Find quasi", "on" if config.get("quasi_cage", {}).get("enabled", False) else "off"],
        ])
        if config.get("quasi_cage", {}).get("enabled", False):
            rows.extend([
                ["Quasi-cage sizes", f"{excel_scalar(config.get('quasi_cage', {}).get('base_sizes', 'auto'))} / {excel_scalar(config.get('quasi_cage', {}).get('side_sizes', 'auto'))}"],
                ["Quasi max layer", config.get("quasi_cage", {}).get("max_layers", "")],
                ["Quasi search policy", config.get("quasi_cage", {}).get("search_policy", "bounded")],
            ])
    rows.extend([
        ["Cage report types", dashboard_cage_targets(config)],
        ["Maximum cage face", config.get("cage", {}).get("max_faces", 20)],
    ])
    if not cpp_mode:
        rows.extend([
            ["Find cluster", "on" if config.get("hydrate_cluster", {}).get("enabled", False) else "off"],
            ["Cluster min cage", config.get("hydrate_cluster", {}).get("min_cage", 2)],
        ])
    rows.append(["Order parameters", order_parameter_display(selected_order_parameters)])
    if not cpp_mode and selected_order_set & {"mcg1", "mcg3"}:
        rows.append([
            "MCG guest / water cutoff (nm)",
            f"{config.get('hydrate_order', {}).get('mcg_guest_cutoff_nm', 0.90)} / "
            f"{config.get('hydrate_order', {}).get('mcg_water_cutoff_nm', 0.60)}",
        ])
    if not cpp_mode and selected_order_set & {"dhop35", "dhop30"}:
        rows.append([
            "DHOP O-O cutoff (nm)",
            config.get("hydrate_order", {}).get("dhop_neighbor_cutoff_nm", 0.35),
        ])
    if not cpp_mode and q_degrees:
        rows.extend([
            ["Q_l degree", excel_scalar(q_degrees)],
            ["Q_l neighbor mode", config.get("order", {}).get("q_neighbor_mode", "graph")],
            ["Q_l cutoff (nm)", config.get("order", {}).get("q_cutoff_nm", 0.35)],
            ["Q_l n neighbor", config.get("order", {}).get("q_n_neighbor", "NULL")],
        ])
    rows.extend([
        [
            "Output types",
            run_info.get("output_types")
            or output_type_display(
                config.get("output", {}).get("types"),
                cpp_mode=cpp_mode,
            ),
        ],
        ["Output layout", run_info.get("output_layout", "")],
        ["Worker policy", run_info.get("worker_policy", "")],
        ["Parallel backend", run_info.get("parallel_backend", "serial")],
        ["Math threads per worker", run_info.get("math_threads", 1)],
        ["Workers", run_info.get("workers", "")],
        ["", ""],
        ["Analysis Results (min / mean / max)", ""],
        ["Frames total / ok / failed", f"{len(data)} / {frames_ok_count(data)} / {frames_failed_count(data)}"],
        ["Water molecules", min_mean_max_column(data, "n_waters")],
        ["Guest molecules", min_mean_max_column(data, "n_guests")],
        ["Connections", min_mean_max_column(data, "connection_count")],
    ])
    if not cpp_mode:
        for size in configured_ring_report_sizes(config):
            rows.append([f"Ring{size}", min_mean_max_column(data, f"ring{size}")])
        rows.extend([
            ["Half cage", min_mean_max_column(data, "half_cage_total")],
            ["Quasi cage", min_mean_max_column(data, "quasi_cage_total")],
        ])
    rows.append(["Cage total", min_mean_max_column(data, "cage_total")])
    if cpp_mode and not has_selected_guests(data):
        rows.append(["Cage occupancy", "not evaluated"])
    else:
        rows.extend([
            ["Empty cage", min_mean_max_column(data, "cage_empty")],
            ["Occupied cage", min_mean_max_column(data, "cage_occupied")],
        ])
    if not cpp_mode and config.get("hydrate_cluster", {}).get("enabled", False):
        rows.extend([
            ["Hydrate cluster", min_mean_max_column(data, "hydrate_cluster_count")],
            ["Isolated cage", min_mean_max_column(data, "isolated_cage_count")],
        ])
    if not cpp_mode:
        rows.append(["Ice-like waters", min_mean_max_column(data, "ice_like_waters")])

    completed_outputs = run_info.get("completed_outputs", run_info.get("output_types", ()))
    if "executed_features" in run_info:
        citation_evidence = {
            key: run_info[key]
            for key in (
                "successful_frames",
                "executed_features",
                "executed_order_parameters",
                "completed_outputs",
                "occupancy_evaluated",
                "guest_residence_evaluated",
            )
            if key in run_info
        }
    else:
        citation_evidence = completed_citation_evidence(
            config,
            successful_frames=frames_ok_count(data),
            completed_outputs=completed_outputs,
            track=bool(run_info.get("track", False)),
            occupancy_evaluated=has_selected_guests(data),
        )
    citation_statistics: dict[str, Any] = {
        **citation_evidence,
        "failed_frames": frames_failed_count(data),
        "status": run_info.get("status", "completed"),
        "guest_molecules": sum_numeric_column(data, "n_guests"),
    }
    citation = build_citation_recommendation(run_info, config, citation_statistics)
    rows.extend(
        [
            ["", ""],
            ["Citation Recommendation", ""],
            ["Recommended text", citation.sentence],
            ["Publication", citation.publication.removeprefix("Publication: ")],
            ["DOI", citation.doi.removeprefix("DOI        : ")],
            ["GitHub", citation.github.removeprefix("GitHub     : ")],
        ]
    )
    return pd.DataFrame(rows)


def summary_markdown(data: pd.DataFrame) -> str:
    """Render the global summary as readable grouped markdown tables."""
    lines = ["# SQQ summary", ""]
    for title, table in summary_markdown_tables(data):
        if table.empty:
            continue
        lines.extend(["", f"## {title}", "", markdown_table(table).rstrip()])
    return "\n".join(lines).strip() + "\n"


def summary_markdown_tables(data: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    """Build human-facing markdown summary tables."""
    frame_columns = [
        "frame",
        "time_ps",
        "source",
        "status",
        "n_atoms",
        "n_waters",
        "n_guests",
        "connection_mode",
        "connection_count",
        "hbond_count",
        "oo_connection_count",
        "pair_connection_count",
    ]
    tables: list[tuple[str, pd.DataFrame]] = [
        ("Failures", failure_summary_table(data)),
        ("Frames", summary_simple_table(data, frame_columns)),
        ("Molecules", molecule_summary_table(data)),
        ("Rings", summary_simple_table(data, ["frame", "time_ps", "ring4", "ring5", "ring6", "ring7", "free_ring4", "free_ring5", "free_ring6", "free_ring7"])),
        ("Half Cage", patch_summary_table(data, "half_cage")),
        ("Quasi Cage", patch_summary_table(data, "quasi_cage")),
    ]
    tables.extend(
        [
            ("Cages", cage_summary_table(data)),
            ("Cage Occupancy", cage_occupancy_summary_table(data, markdown_style=True)),
            ("Cage Isomers", cage_isomer_summary_table(data, include_zero_rows=False)),
            ("Order Parameters", order_parameter_summary_table(data)),
            ("Ice", summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"])),
        ]
    )
    return tables


def summary_sheet_tables(data: pd.DataFrame, config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """Build lightweight main-summary tables using the configured scopes."""
    if is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        tables: dict[str, pd.DataFrame] = {
            "failures": failure_summary_table(data),
            "cage": cage_summary_table(data),
            "order_parameter": order_parameter_summary_table(
                data,
                config.get("order", {}).get("parameters", ["f3", "f4"]),
                include_focus=bool(config.get("order", {}).get("focus_waters", [])),
            ),
        }
        return {name: table for name, table in tables.items() if not table.empty}

    ring_sizes = configured_ring_report_sizes(config)
    ring_columns = ["frame", "time_ps"]
    for size in ring_sizes:
        ring_columns.extend([f"ring{size}", f"free_ring{size}"])
    tables: dict[str, pd.DataFrame] = {
        "failures": failure_summary_table(data),
        connection_sheet_name(data): connection_summary_table(data),
        "ring": summary_simple_table(data, ring_columns),
        "half_cage": patch_summary_table(data, "half_cage"),
        "quasi_cage": patch_summary_table(data, "quasi_cage"),
        "cage": cage_summary_table(data),
        "hydrate_cluster": hydrate_cluster_summary_table(data),
        "order_parameter": order_parameter_summary_table(
            data,
            config.get("order", {}).get("parameters", ["f3", "f4"]),
            include_focus=bool(config.get("order", {}).get("focus_waters", [])),
        ),
        "ice": summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"]),
    }
    return {name: table for name, table in tables.items() if not table.empty}


def failure_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Return one diagnostic row per failed input frame."""
    if "status" not in data.columns:
        return pd.DataFrame()
    columns = [
        column
        for column in ("frame", "time_ps", "source", "status", "error")
        if column in data.columns
    ]
    failed = data.loc[
        data["status"].astype(str).str.lower().eq("failed"),
        columns,
    ]
    return failed.reset_index(drop=True)

__all__ = (
    "summary_output_tables",
    "summary_dashboard_table",
    "summary_markdown",
    "summary_markdown_tables",
    "summary_sheet_tables",
    "failure_summary_table",
)
