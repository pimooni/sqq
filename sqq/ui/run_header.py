"""Run banner, configuration header, and normalized run metadata."""

from __future__ import annotations

from argparse import Namespace
from datetime import datetime
from pathlib import Path
from typing import Any

from .. import __version__
from ..banner import banner_for_engine
from ..config import (
    DEFAULT_MODE,
    is_cpp_mode,
    normalize_order_parameters,
    normalize_output_types,
    order_parameter_display,
    output_enabled,
    output_type_display,
    parse_on_off,
    profile_name,
    q_degrees_from_order_parameters,
)
from ..display import graph_mode_display, ordered_unique_graph_modes
from ..io.reporting.tables import dashboard_cage_targets
from ..io.trajectory import TrajectorySelection
from ..runtime.parallel.policy import worker_policy_text
from .formatting import format_time_zone, terminal_field_line, write_terminal_block


BOND_MODE_DISPLAY_NAMES = {
    "auto": "auto",
    "hbond": "hydrogen bond",
    "oo": "O-O connectivity",
    "pairs": "user-defined pairs",
}


def print_run_banner(engine: str | None = None) -> None:
    """Print the SQQ banner before preparation can emit diagnostics."""
    write_terminal_block([banner_for_engine(engine)])


def sampling_metadata(selection: TrajectorySelection) -> dict[str, Any]:
    """Return compact YAML-safe metadata for a resolved frame selection."""
    return {
        "native_frame_interval_ps": selection.native_interval_ps,
        "delta_time_ps": selection.delta_time_ps,
        "raw_frame_step": selection.raw_frame_step,
        "selected_frames": selection.selected_frames,
        "total_frames": selection.total_frames,
    }


def sampling_interval_display(config: dict[str, Any]) -> str:
    """Render the one-line terminal summary for resolved sampling."""
    sampling = config.get("input", {}).get("sampling", {})
    if not sampling:
        return ""
    requested = sampling.get("delta_time_ps")
    native = sampling.get("native_frame_interval_ps")
    selected = sampling.get("selected_frames", 0)
    total = sampling.get("total_frames", 0)
    interval_text = "all" if requested is None else f"{format_ps(requested)} ps"
    native_text = "unknown" if native is None else f"{format_ps(native)} ps"
    return f"{interval_text} [native {native_text}; {selected} of {total} frames]"


def sampling_selected_frames_text(config: dict[str, Any]) -> str:
    """Render selected and available frame counts for reports."""
    sampling = config.get("input", {}).get("sampling", {})
    if not sampling:
        return ""
    return f"{sampling.get('selected_frames', 0)} / {sampling.get('total_frames', 0)}"


def format_ps(value: Any) -> str:
    """Format a finite ps value without unnecessary decimal zeros."""
    try:
        return f"{float(value):g}"
    except (TypeError, ValueError):
        return str(value)


def row_effective_graph_modes(rows: list[dict[str, Any]]) -> list[str]:
    """Collect effective graph modes from successful summary rows."""
    modes: list[str] = []
    for row in rows:
        if str(row.get("status", "ok")).lower() == "failed":
            continue
        mode = str(row.get("connection_mode", "")).strip()
        if mode:
            modes.append(mode)
    return modes


def build_run_info(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    outdir: Path,
    paths: list[Path],
    topology: Path | None,
    workers: int,
    parallel_backend: str,
    elapsed_seconds: float,
    started_at_wall: datetime,
    finished_at_wall: datetime,
    rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Collect run metadata for the terminal, config, and summaries."""
    requested_graph_mode = config["graph"]["bond_mode"]
    effective_graph_modes = row_effective_graph_modes(rows or [])
    graph_config = config.get("graph", {})
    configured_graph_mode = str(graph_config.get("effective_bond_mode", "")).strip()
    raw_group_modes = graph_config.get("effective_bond_mode_by_group", {})
    raw_group_reasons = graph_config.get("effective_bond_mode_reason_by_group", {})
    configured_group_modes = (
        {
            str(label): str(mode).strip()
            for label, mode in raw_group_modes.items()
            if str(mode).strip()
        }
        if isinstance(raw_group_modes, dict)
        else {}
    )
    configured_group_reasons = (
        {
            str(label): str(reason).strip()
            for label, reason in raw_group_reasons.items()
            if str(reason).strip()
        }
        if isinstance(raw_group_reasons, dict)
        else {}
    )
    configured_graph_reason = str(
        graph_config.get("effective_bond_mode_reason", "")
    ).strip()
    allowed_modes = set(configured_group_modes.values())
    if configured_graph_mode:
        allowed_modes.add(configured_graph_mode)
    if effective_graph_modes and allowed_modes:
        if any(mode not in allowed_modes for mode in effective_graph_modes):
            raise RuntimeError("Per-frame graph mode disagrees with preflight resolution.")
    elif not effective_graph_modes:
        effective_graph_modes = (
            list(configured_group_modes.values())
            if configured_group_modes
            else [configured_graph_mode]
        )
    unique_effective_modes = ordered_unique_graph_modes(effective_graph_modes)
    graph_mode_by_group = {
        label: graph_mode_display(requested_graph_mode, [mode])
        for label, mode in configured_group_modes.items()
    }
    graph_mode_reason_by_group = {
        label: configured_group_reasons.get(label, "")
        for label in configured_group_modes
        if configured_group_reasons.get(label, "")
    }
    graph_display = (
        graph_mode_display(requested_graph_mode, unique_effective_modes)
        if len(unique_effective_modes) == 1 or requested_graph_mode != "auto"
        else ""
    )
    selected_order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters")
    )
    q_degrees = q_degrees_from_order_parameters(selected_order_parameters)
    output_types = normalize_output_types(config.get("output", {}).get("types"))
    result_rows = rows or []
    failures = [
        {
            "frame": str(row.get("frame", "")),
            "source": str(row.get("source", "")),
            "error": str(row.get("error", "")),
        }
        for row in result_rows
        if str(row.get("status", "")).lower() == "failed"
    ]
    input_format = input_format_label(paths)
    info: dict[str, Any] = {
        "working_dir": str(Path.cwd()),
        "input": str(input_path),
        "input_format": input_format,
        "output_dir": str(outdir.resolve()),
        "date": started_at_wall.strftime("%Y-%m-%d"),
        "start_time": started_at_wall.strftime("%H:%M:%S"),
        "finish_time": finished_at_wall.strftime("%H:%M:%S"),
        "started_at": started_at_wall.isoformat(timespec="seconds"),
        "finished_at": finished_at_wall.isoformat(timespec="seconds"),
        "time_zone": format_time_zone(started_at_wall),
        "config_file": args.config or "<built-in defaults>",
        "sqq_version": __version__,
        "engine_selector": config.get("mode", DEFAULT_MODE),
        "sqq_engine": "sqq-cpp"
        if is_cpp_mode(config.get("mode", DEFAULT_MODE))
        else "sqq-py",
        "profile": profile_name(config.get("mode", DEFAULT_MODE)),
        "resolution_adjustments": list(
            config.get("resolution_report", {}).get("adjustments", ())
        ),
        "worker_policy": worker_policy_text(config),
        "topology": str(topology) if topology else "<none>",
        "matched_files": len(paths),
        "delta_time_ps": config["input"].get("delta_time_ps"),
        "sampling_interval": sampling_interval_display(config),
        "native_frame_interval_ps": config["input"].get("sampling", {}).get("native_frame_interval_ps"),
        "raw_frame_step": config["input"].get("sampling", {}).get("raw_frame_step", 1),
        "selected_frames": config["input"].get("sampling", {}).get("selected_frames", len(result_rows)),
        "source_frames_total": config["input"].get("sampling", {}).get("total_frames", len(result_rows)),
        "frames_total": len(result_rows),
        "frames_ok": sum(str(row.get("status", "")).lower() == "ok" for row in result_rows),
        "frames_failed": len(failures),
        "failures": failures,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "graph_mode": requested_graph_mode,
        "effective_graph_modes": ", ".join(unique_effective_modes),
        "graph_mode_display": graph_display,
        "graph_mode_reason": config.get("graph", {}).get("effective_bond_mode_reason", configured_graph_reason),
        "graph_mode_by_group": graph_mode_by_group,
        "graph_mode_reason_by_group": graph_mode_reason_by_group,
        "search_sizes": config["ring"]["sizes"],
        "ring_report_sizes": config["ring"]["report_sizes"],
        "find_half": on_off_text(config.get("half_cage", {}).get("enabled", False)),
        "find_quasi": on_off_text(config.get("quasi_cage", {}).get("enabled", False)),
        "quasi_cage_base_sizes": config["quasi_cage"].get("base_sizes", "auto"),
        "quasi_cage_side_sizes": config["quasi_cage"].get("side_sizes", "auto"),
        "cage_report_types": config["cage"].get("report_types", []),
        "max_cage_face": config["cage"].get("max_faces", 20),
        "cage_scientific_validation": on_off_text(config["cage"].get("scientific_validation", False)),
        "find_cluster": on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)),
        "cluster_min_cage": config.get("hydrate_cluster", {}).get("min_cage", 2),
        "order_parameters": order_parameter_display(selected_order_parameters),
        "hydrate_order": hydrate_order_config_text(config),
        "mcg3": on_off_text("mcg3" in selected_order_parameters),
        "dhop30": on_off_text("dhop30" in selected_order_parameters),
        "dhop_neighbor_cutoff_nm": config.get("hydrate_order", {}).get("dhop_neighbor_cutoff_nm", 0.35),
        "q_enabled": bool(q_degrees),
        "q_degree": list(q_degrees),
        "q_neighbor_mode": config["order"].get("q_neighbor_mode", "graph"),
        "q_cutoff_nm": config["order"].get("q_cutoff_nm", 0.35),
        "q_n_neighbor": config["order"].get("q_n_neighbor", None),
        "output_types": output_type_display(
            output_types, cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE))
        ),
        "output_layout": config["output"].get("structure_layout", "grouped"),
        "workers": workers,
        "parallel_backend": parallel_backend,
        "math_threads": int(config.get("parallel", {}).get("math_threads", 1)),
        "summary_xlsx": str((outdir / "summary.xlsx").resolve())
        if output_enabled(config, "summary-xlsx")
        else "<disabled>",
        "summary_csv": _summary_csv_path(outdir, config)
        if output_enabled(config, "summary-csv")
        else "<disabled>",
        "summary_detail_csv": _summary_csv_path(outdir, config)
        if output_enabled(config, "summary-detail-csv")
        or output_enabled(config, "cluster-detail")
        else "<disabled>",
        "config_output": str((outdir / "sqq_config_resolved.yaml").resolve()),
    }
    if input_format.startswith("lammps-"):
        lammps = config["input"].get("lammps", {})
        info.update(
            {
                "lammps_units": lammps.get("units", "real"),
                "lammps_timestep": lammps.get("timestep", 1.0),
                "lammps_atom_style": lammps.get("atom_style", "full"),
                "lammps_type_map_source": lammps.get("type_map_source", "<configuration>"),
            }
        )
    if paths:
        info["first_file"] = str(paths[0].resolve())
        info["last_file"] = str(paths[-1].resolve())
    return info


def _summary_csv_path(outdir: Path, config: dict[str, Any]) -> str:
    directory = str(config.get("output", {}).get("summary_csv_dir", "summary"))
    return str((outdir / directory).resolve())


def print_run_header(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    outdir: Path,
    paths: list[Path],
    topology: Path | None,
    workers: int,
    parallel_backend: str,
    started_at_wall: datetime,
    *,
    extra_configuration_fields: list[tuple[str, Any]] | None = None,
) -> None:
    """Print one atomic block of static run information."""
    lines = ["Basic Information"]

    def add_field(label: str, value: Any) -> None:
        lines.append(terminal_field_line(label, value))

    add_field("date", started_at_wall.strftime("%Y-%m-%d"))
    add_field("start_time", started_at_wall.strftime("%H:%M:%S"))
    add_field("time_zone", format_time_zone(started_at_wall))
    add_field("working_dir", Path.cwd())
    add_field("input", input_path)
    add_field("input_format", input_format_label(paths))
    add_field("matched_files", len(paths))
    add_field("output", outdir)
    lines.extend(["", "Configuration"])
    add_field("SQQ version", __version__)
    add_field(
        "SQQ engine",
        "sqq-cpp" if is_cpp_mode(config.get("mode", DEFAULT_MODE)) else "sqq-py",
    )
    add_field("Engine selector", config.get("mode", DEFAULT_MODE))
    add_field("Profile", profile_name(config.get("mode", DEFAULT_MODE)))
    add_field("Config file", args.config or "<built-in defaults>")
    add_field("Topology", topology or "<none>")
    if config["input"].get("sampling"):
        add_field("Sampling Interval", sampling_interval_display(config))
    if input_format_label(paths).startswith("lammps-"):
        lammps = config["input"].get("lammps", {})
        add_field("LAMMPS units", lammps.get("units", "real"))
        add_field("LAMMPS timestep", lammps.get("timestep", 1.0))
        add_field("LAMMPS atom style", lammps.get("atom_style", "full"))
        add_field("LAMMPS type map", lammps.get("type_map_source", "<configuration>"))
    requested_graph_mode = config["graph"]["bond_mode"]
    effective_graph_mode = config["graph"].get("effective_bond_mode", "")
    effective_graph_reason = config["graph"].get("effective_bond_mode_reason", "")
    group_graph_modes = config["graph"].get("effective_bond_mode_by_group", {})
    group_graph_reasons = config["graph"].get("effective_bond_mode_reason_by_group", {})
    if requested_graph_mode != "auto" or effective_graph_mode:
        add_field(
            "Graph mode",
            graph_mode_display(requested_graph_mode, [effective_graph_mode]),
        )
        if effective_graph_reason:
            add_field("Graph mode reason", effective_graph_reason)
    elif isinstance(group_graph_modes, dict) and group_graph_modes:
        for label, mode in group_graph_modes.items():
            add_field(
                f"Graph mode ({label})",
                graph_mode_display(requested_graph_mode, [mode]),
            )
            reason = (
                group_graph_reasons.get(label, "")
                if isinstance(group_graph_reasons, dict)
                else ""
            )
            if reason:
                add_field(f"Graph reason ({label})", reason)
    else:
        graph_mode_display(requested_graph_mode, [])
    add_field("Search sizes", config["ring"]["sizes"])
    add_field("Ring definition", config["ring"].get("definition", "chordless"))
    if not is_cpp_mode(config.get("mode")):
        add_field("Ring report sizes", config["ring"]["report_sizes"])
        add_field("Find half", on_off_text(config.get("half_cage", {}).get("enabled", False)))
        add_field("Find quasi", on_off_text(config.get("quasi_cage", {}).get("enabled", False)))
        if config.get("quasi_cage", {}).get("enabled", False):
            add_field(
                "Quasi-cage sizes",
                f"{config['quasi_cage'].get('base_sizes', 'auto')} / "
                f"{config['quasi_cage'].get('side_sizes', 'auto')}",
            )
            add_field("Quasi max layer", config["quasi_cage"].get("max_layers", ""))
            add_field(
                "Quasi search policy",
                config["quasi_cage"].get("search_policy", "bounded"),
            )
    add_field("Cage report types", dashboard_cage_targets(config))
    add_field("Maximum cage face", config["cage"].get("max_faces", 20))
    add_field(
        "Scientific validation",
        on_off_text(config["cage"].get("scientific_validation", False)),
    )
    if not is_cpp_mode(config.get("mode")):
        add_field(
            "Find cluster",
            on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)),
        )
        add_field("Cluster min cage", config.get("hydrate_cluster", {}).get("min_cage", 2))
    add_field("Order parameters", order_parameter_config_text(config))
    if q_degrees_from_order_parameters(config.get("order", {}).get("parameters")):
        add_field("Q_l settings", q_config_text(config))
    add_field(
        "Output types",
        output_type_display(
            config.get("output", {}).get("types"),
            cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        ),
    )
    add_field("Output layout", config["output"].get("structure_layout", "grouped"))
    add_field("Worker policy", worker_policy_text(config))
    add_field("Parallel backend", parallel_backend)
    add_field("Math threads per worker", config.get("parallel", {}).get("math_threads", 1))
    add_field("Workers", workers)
    for adjustment in config.get("resolution_report", {}).get("adjustments", ()):
        if not isinstance(adjustment, dict):
            continue
        parameter = adjustment.get("parameter", "parameter")
        effective = adjustment.get("effective")
        reason = adjustment.get("reason", "automatic adjustment")
        value = f"{parameter}"
        if effective is not None:
            value += f" -> {effective}"
        value += f" [{reason}]"
        add_field("Adjustment", value)
    for label, value in extra_configuration_fields or []:
        add_field(label, value)
    lines.append("")
    write_terminal_block(lines)


def print_run_summary(run_info: dict[str, Any]) -> None:
    """Print final effective run metadata as one terminal block."""
    lines = ["Run Summary"]

    def add_field(label: str, value: Any) -> None:
        lines.append(terminal_field_line(label, value))

    add_field("Finish time", run_info.get("finish_time", ""))
    add_field("Duration (s)", run_info.get("elapsed_seconds", ""))
    add_field("SQQ version", run_info.get("sqq_version", __version__))
    add_field("SQQ engine", run_info.get("sqq_engine", ""))
    add_field("Engine selector", run_info.get("engine_selector", ""))
    add_field("Profile", run_info.get("profile", ""))
    graph_mode_by_group = run_info.get("graph_mode_by_group", {})
    if isinstance(graph_mode_by_group, dict) and len(graph_mode_by_group) > 1:
        for label, mode in graph_mode_by_group.items():
            add_field(f"Graph mode ({label})", mode)
    else:
        add_field(
            "Graph mode",
            run_info.get("graph_mode_display", run_info.get("graph_mode", "")),
        )
        if run_info.get("graph_mode_reason"):
            add_field("Graph mode reason", run_info.get("graph_mode_reason", ""))
    add_field("Order parameters", run_info.get("order_parameters", ""))
    if run_info.get("sqq_engine") != "sqq-cpp":
        add_field("Find cluster", run_info.get("find_cluster", "off"))
    add_field("Output types", run_info.get("output_types", "none"))
    add_field("Worker policy", run_info.get("worker_policy", ""))
    add_field("Parallel backend", run_info.get("parallel_backend", "serial"))
    add_field("Workers", run_info.get("workers", ""))
    summary_write = run_info.get("summary_write", {})
    if isinstance(summary_write, dict) and "total_seconds" in summary_write:
        add_field("Summary write (s)", summary_write.get("total_seconds", ""))
    lines.append("")
    write_terminal_block(lines)


def frame_input_metadata(config: dict[str, Any]) -> dict[str, Any]:
    """Return normalized input provenance for one frame report."""
    input_config = config.get("input", {})
    metadata = {
        "input_format": input_config.get("format", ""),
        "topology": input_config.get("topology"),
        "sampling_interval": sampling_interval_display(config),
        "native_frame_interval_ps": input_config.get("sampling", {}).get("native_frame_interval_ps"),
        "delta_time_ps": input_config.get("delta_time_ps"),
        "raw_frame_step": input_config.get("sampling", {}).get("raw_frame_step"),
        "selected_frames": sampling_selected_frames_text(config),
        "find_half": on_off_text(config.get("half_cage", {}).get("enabled", False)),
        "find_quasi": on_off_text(config.get("quasi_cage", {}).get("enabled", False)),
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


def input_format_label(paths: list[Path]) -> str:
    """Return one compact source-format label for run metadata."""
    mapping = {
        ".gro": "gromacs-gro",
        ".xyz": "xyz",
        ".xtc": "gromacs-xtc",
        ".trr": "gromacs-trr",
        ".dump": "lammps-dump",
        ".lammpstrj": "lammps-dump",
        ".dcd": "lammps-dcd",
    }
    labels: list[str] = []
    for path in paths:
        label = mapping.get(path.suffix.lower(), path.suffix.lower().lstrip("."))
        if label not in labels:
            labels.append(label)
    if not labels:
        return "unknown"
    if len(labels) == 1:
        return labels[0]
    return "mixed (" + ", ".join(labels) + ")"


def bond_mode_display_name(value: Any) -> str:
    """Return a readable label without changing config identifiers."""
    mode = str(value)
    return BOND_MODE_DISPLAY_NAMES.get(mode, mode)


def hydrate_order_config_text(config: dict[str, Any]) -> str:
    """Render the selected MCG/DHOP subset for compatibility metadata."""
    parameters = normalize_order_parameters(config.get("order", {}).get("parameters"))
    active = [
        name for name in parameters if name in {"mcg1", "mcg3", "dhop35", "dhop30"}
    ]
    return ", ".join(active) if active else "disabled"


def order_parameter_config_text(config: dict[str, Any]) -> str:
    """Render the unified order-parameter selection."""
    return order_parameter_display(config.get("order", {}).get("parameters"))


def q_config_text(config: dict[str, Any]) -> str:
    """Render Steinhardt Q_l settings for the run header."""
    order = config.get("order", {})
    degrees = q_degrees_from_order_parameters(order.get("parameters"))
    if not degrees:
        return "disabled"
    n_neighbors = order.get("q_n_neighbor", None)
    n_text = "NULL" if n_neighbors in (None, "", "null", "NULL") else str(n_neighbors)
    degree = ",".join(str(item) for item in degrees)
    return (
        f"degree={degree}; mode={order.get('q_neighbor_mode', 'graph')}; "
        f"cutoff={order.get('q_cutoff_nm', 0.35)} nm; n={n_text}"
    )


def on_off_text(value: Any) -> str:
    """Render on/off settings using the CLI vocabulary."""
    return "on" if parse_on_off(value, "on/off setting") else "off"


__all__ = [
    "BOND_MODE_DISPLAY_NAMES",
    "bond_mode_display_name",
    "build_run_info",
    "format_ps",
    "frame_input_metadata",
    "hydrate_order_config_text",
    "input_format_label",
    "on_off_text",
    "order_parameter_config_text",
    "print_run_banner",
    "print_run_header",
    "print_run_summary",
    "q_config_text",
    "row_effective_graph_modes",
    "sampling_interval_display",
    "sampling_metadata",
    "sampling_selected_frames_text",
]
