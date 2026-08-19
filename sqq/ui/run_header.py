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
from .formatting import (
    format_started,
    format_time_zone,
    terminal_field_line,
    write_terminal_block,
)


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
    if requested_graph_mode == "auto":
        if not unique_effective_modes:
            raise RuntimeError(
                "Graph mode auto was not resolved before run information was built."
            )
        invalid_modes = [
            mode for mode in unique_effective_modes if mode not in {"hbond", "oo"}
        ]
        if invalid_modes:
            raise RuntimeError(
                "Graph mode auto resolved to an unsupported effective mode: "
                + ", ".join(invalid_modes)
            )
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
    successful_rows = [
        row for row in result_rows if str(row.get("status", "")).lower() == "ok"
    ]
    guest_molecules = sum(
        int(row.get("n_guests", 0) or 0) for row in successful_rows
    )
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
        "guest_molecules": guest_molecules,
        "occupancy_evaluated": guest_molecules > 0,
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


def compact_sqq_display(selector: Any, *, version: str = __version__) -> str:
    """Return one combined engine, selector, and version label."""
    mode = str(selector or DEFAULT_MODE).strip().lower()
    engine = "sqq-cpp" if is_cpp_mode(mode) else "sqq-py"
    if mode in {"00", "99"}:
        engine = f"{engine}-{mode}"
    return f"{engine} ({version})"


def compact_graph_reason(reason: Any) -> str:
    """Reduce graph preflight explanations to one stable terminal phrase."""
    text = " ".join(str(reason or "").split())
    normalized = text.casefold()
    if normalized == "complete water hydrogen topology":
        return "all waters have >=2 H"
    if normalized == "oxygen-only water topology":
        return "waters only have O"
    if normalized == "explicit graph mode":
        return "explicit"
    if normalized == "explicit pair map":
        return "pair map"
    return text


def compact_graph_display(mode: Any, reason: Any = "") -> str:
    """Join requested/effective graph mode and its concise reason."""
    display = str(mode or "").strip()
    concise = compact_graph_reason(reason)
    return f"{display} [{concise}]" if display and concise else display


def compact_group_graph_display(
    requested: Any,
    modes: dict[Any, Any],
    reasons: dict[Any, Any] | None = None,
) -> str:
    """Group topology labels that resolved to the same graph mode/reason."""
    grouped: dict[tuple[str, str], list[str]] = {}
    reason_map = reasons or {}
    for label, mode in modes.items():
        reason = compact_graph_reason(reason_map.get(label, ""))
        key = (str(mode), reason)
        grouped.setdefault(key, []).append(str(label))
    parts: list[str] = []
    for (mode, reason), labels in grouped.items():
        base_display = (
            mode
            if "->" in mode
            else graph_mode_display(requested, [mode])
        )
        display = compact_graph_display(
            base_display, reason
        )
        parts.append(f"{'/'.join(labels)}: {display}")
    return "; ".join(parts)


def compact_input_display(path: Any, input_format: Any) -> str:
    """Join a resolved input path and human-facing format label."""
    raw = str(path or "").strip()
    try:
        resolved = str(Path(raw).expanduser().resolve()) if raw else ""
    except (OSError, RuntimeError):
        resolved = raw
    format_text = human_input_format(input_format)
    return f"{resolved} [{format_text}]" if format_text else resolved


def human_input_format(value: Any) -> str:
    """Return a compact display form for normalized input-format identifiers."""
    text = str(value or "").strip()
    mapping = {
        "gromacs-gro": "GROMACS GRO",
        "gromacs-xtc": "GROMACS XTC",
        "gromacs-trr": "GROMACS TRR",
        "lammps-dump": "LAMMPS trajectory",
        "lammps-dcd": "LAMMPS DCD",
        "track-state": "SQQ Analyze result",
        "xyz": "XYZ",
    }
    return mapping.get(text.casefold(), text)


def compact_ring_scope(config: dict[str, Any]) -> str:
    """Combine ring search and report sizes on one terminal row."""
    ring = config.get("ring", {})
    search = _size_scope_text(ring.get("sizes", ()))
    report_value = ring.get("report_sizes", ring.get("sizes", ()))
    if report_value in (None, "", "auto"):
        report_value = ring.get("sizes", ())
    report = _size_scope_text(report_value)
    if search == report:
        return f"{search} [search and report]"
    return f"search {search}; report {report}"


def compact_additional_search(config: dict[str, Any]) -> str:
    """Render HALF, QUASI, and CLUSTER as one compact feature line."""
    half = bool(config.get("half_cage", {}).get("enabled", False))
    quasi = bool(config.get("quasi_cage", {}).get("enabled", False))
    cluster = bool(config.get("hydrate_cluster", {}).get("enabled", False))
    half_text = f"HALF {'on' if half else 'off'}"
    quasi_text = f"QUASI {'on' if quasi else 'off'}"
    if quasi:
        quasi_text += f" (L{config.get('quasi_cage', {}).get('max_layers', 1)})"
    cluster_text = f"CLUSTER {'on' if cluster else 'off'}"
    if cluster:
        cluster_text += (
            f" (>={config.get('hydrate_cluster', {}).get('min_cage', 2)} cages)"
        )
    return "; ".join((half_text, quasi_text, cluster_text))


def _size_scope_text(value: Any) -> str:
    if isinstance(value, str):
        parts = [item.strip() for item in value.replace("/", ",").split(",") if item.strip()]
    elif isinstance(value, (list, tuple, set)):
        parts = [str(item) for item in value]
    elif value in (None, ""):
        parts = []
    else:
        parts = [str(value)]
    return "/".join(parts) if parts else "none"


def _explicit_cage_report(config: dict[str, Any]) -> str:
    raw = config.get("cage", {}).get("report_types", "auto")
    if raw in (None, "", (), [], "auto", "all"):
        return ""
    value = dashboard_cage_targets(config)
    return "" if value.startswith("all detected") else value


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
    """Print the compact static Analyze header used above the live panel."""
    lines = ["Basic Information"]

    def add_field(label: str, value: Any) -> None:
        lines.append(terminal_field_line(label, value))

    input_format = input_format_label(paths)
    add_field("Started", format_started(started_at_wall))
    add_field("Input", compact_input_display(input_path, input_format))
    add_field("Output", outdir.resolve())
    lines.extend(["", "Configuration"])
    add_field("SQQ", compact_sqq_display(config.get("mode", DEFAULT_MODE)))
    add_field("Config file", args.config or "<built-in defaults>")
    add_field("Topology", topology or "<none>")
    if config["input"].get("sampling"):
        add_field("Sampling Interval", sampling_interval_display(config))
    if input_format.startswith("lammps-"):
        lammps = config["input"].get("lammps", {})
        add_field(
            "LAMMPS",
            (
                f"units {lammps.get('units', 'real')}; "
                f"timestep {lammps.get('timestep', 1.0)}; "
                f"style {lammps.get('atom_style', 'full')}; "
                f"type map {lammps.get('type_map_source', '<configuration>')}"
            ),
        )
    requested_graph_mode = config["graph"]["bond_mode"]
    effective_graph_mode = config["graph"].get("effective_bond_mode", "")
    effective_graph_reason = config["graph"].get("effective_bond_mode_reason", "")
    group_graph_modes = config["graph"].get("effective_bond_mode_by_group", {})
    group_graph_reasons = config["graph"].get("effective_bond_mode_reason_by_group", {})
    if requested_graph_mode != "auto" or effective_graph_mode:
        add_field("Graph", compact_graph_display(
            graph_mode_display(requested_graph_mode, [effective_graph_mode]),
            effective_graph_reason,
        ))
    elif isinstance(group_graph_modes, dict) and group_graph_modes:
        add_field(
            "Graph",
            compact_group_graph_display(
                requested_graph_mode,
                group_graph_modes,
                group_graph_reasons if isinstance(group_graph_reasons, dict) else {},
            ),
        )
    else:
        raise RuntimeError(
            "Graph mode auto was not resolved before the run header was printed."
        )
    add_field("Ring sizes", compact_ring_scope(config))
    add_field("Ring definition", config["ring"].get("definition", "chordless"))
    add_field("Additional search", compact_additional_search(config))
    explicit_cages = _explicit_cage_report(config)
    if explicit_cages:
        add_field("Cage report", explicit_cages)
    add_field("Maximum cage face", config["cage"].get("max_faces", 20))
    add_field(
        "Cage validation",
        on_off_text(config["cage"].get("scientific_validation", False)),
    )
    add_field("Order parameters", order_parameter_config_text(config))
    if q_degrees_from_order_parameters(config.get("order", {}).get("parameters")):
        add_field("Q_l settings", q_config_text(config))
    add_field(
        "Output",
        output_type_display(
            config.get("output", {}).get("types"),
            cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        ),
    )
    adjustment_values: list[str] = []
    for adjustment in config.get("resolution_report", {}).get("adjustments", ()):
        if not isinstance(adjustment, dict):
            continue
        parameter = adjustment.get("parameter", "parameter")
        effective = adjustment.get("effective")
        value = f"{parameter}"
        if effective is not None:
            value += f" -> {effective}"
        adjustment_values.append(value)
    if adjustment_values:
        add_field("Adjustments", "; ".join(adjustment_values))
    for label, value in extra_configuration_fields or []:
        add_field(label, value)
    lines.append("")
    write_terminal_block(lines)


def print_track_header(
    args: Namespace,
    config: dict[str, Any],
    targets: list[Any] | tuple[Any, ...],
    started_at_wall: datetime,
) -> None:
    """Print the compact static Track header before raw/source processing."""
    raw_input = getattr(args, "input", None)
    raw_source = getattr(args, "source", None)
    input_value = raw_input or raw_source or Path.cwd()
    if raw_input:
        input_format = input_format_label([Path(raw_input)])
    else:
        input_format = "track-state"
    target_text = ",".join(str(getattr(item, "raw", item)) for item in targets)
    lines = [
        "Basic Information",
        terminal_field_line("Started", format_started(started_at_wall)),
        terminal_field_line("Input", compact_input_display(input_value, input_format)),
        terminal_field_line("Output", Path(args.output).resolve()),
        "",
        "Configuration",
        terminal_field_line("SQQ", compact_sqq_display(config.get("mode", DEFAULT_MODE))),
        terminal_field_line("Config file", getattr(args, "config", None) or "<built-in defaults>"),
        terminal_field_line("Topology", getattr(args, "topology", None) or "<none>"),
        terminal_field_line("Target", target_text or "all"),
    ]
    if raw_input and config.get("input", {}).get("delta_time_ps") is not None:
        lines.append(terminal_field_line(
            "Sampling interval",
            f"{format_ps(config['input']['delta_time_ps'])} ps",
        ))
    lines.append("")
    write_terminal_block(lines)


def print_run_summary(run_info: dict[str, Any]) -> None:
    """Print a compact compatibility summary for legacy API callers."""
    lines = ["Run Summary"]

    def add_field(label: str, value: Any) -> None:
        lines.append(terminal_field_line(label, value))

    add_field(
        "SQQ",
        compact_sqq_display(
            run_info.get("engine_selector", DEFAULT_MODE),
            version=str(run_info.get("sqq_version", __version__)),
        ),
    )
    graph_mode_by_group = run_info.get("graph_mode_by_group", {})
    if isinstance(graph_mode_by_group, dict) and len(graph_mode_by_group) > 1:
        add_field(
            "Graph",
            compact_group_graph_display(
                run_info.get("graph_mode", "auto"),
                graph_mode_by_group,
                run_info.get("graph_mode_reason_by_group", {}),
            ),
        )
    else:
        add_field(
            "Graph",
            compact_graph_display(
                run_info.get("graph_mode_display", run_info.get("graph_mode", "")),
                run_info.get("graph_mode_reason", ""),
            ),
        )
    add_field("Order parameters", run_info.get("order_parameters", ""))
    add_field("Output", run_info.get("output_types", "none"))
    add_field(
        "Frames",
        (
            f"analyzed {run_info.get('frames_total', 0)}; "
            f"ok {run_info.get('frames_ok', 0)}; "
            f"failed {run_info.get('frames_failed', 0)}"
        ),
    )
    add_field("Time", f"total {run_info.get('elapsed_seconds', 0)} s")
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
    "compact_additional_search",
    "compact_graph_display",
    "compact_group_graph_display",
    "compact_graph_reason",
    "compact_input_display",
    "compact_ring_scope",
    "compact_sqq_display",
    "format_ps",
    "frame_input_metadata",
    "hydrate_order_config_text",
    "input_format_label",
    "on_off_text",
    "order_parameter_config_text",
    "print_run_banner",
    "print_run_header",
    "print_run_summary",
    "print_track_header",
    "q_config_text",
    "row_effective_graph_modes",
    "sampling_interval_display",
    "sampling_metadata",
    "sampling_selected_frames_text",
]
