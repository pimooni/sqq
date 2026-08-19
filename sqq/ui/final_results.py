"""Build the final SQQ terminal results screen.

This module deliberately has no dependency on the Analyze workflow, pandas, or
the output writers.  Callers pass the resolved run metadata, resolved
configuration, and a small mapping of final statistics.  Values in
``statistics`` take precedence over their legacy ``run_info`` equivalents.
"""

from __future__ import annotations

import math
import re
import sys
from collections.abc import Iterable, Mapping
from datetime import datetime
from pathlib import Path
from typing import Any, TextIO

from .. import __version__
from ..citation import build_citation_recommendation
from ..io.render.inspect import inspect_render_script, render_launch_commands
from .formatting import format_started
from .run_header import (
    compact_additional_search,
    compact_graph_display,
    compact_group_graph_display,
    compact_input_display,
    compact_output_display,
    compact_ring_scope,
    compact_sqq_display,
)

_ANSI_BOLD = "\x1b[1m"
_ANSI_RESET = "\x1b[0m"
_LABEL_WIDTH = 24
_MISSING = object()
_TRUE_TEXT = {"1", "true", "yes", "on", "enabled", "complete", "completed"}
_FALSE_TEXT = {"", "0", "false", "no", "off", "disabled", "none", "null", "auto"}
_TOKEN_SPLIT = re.compile(r"[,;/\s]+")


def render_final_results(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
    *,
    ansi: bool = True,
) -> str:
    """Return the complete final terminal screen without a trailing newline.

    Canonical statistic keys are ``requested_frames``, ``analyzed_frames``,
    ``successful_frames``, ``failed_frames``, ``total_seconds``,
    ``analysis_seconds``, ``write_seconds``, ``status``, and ``result_path``.
    Current ``run_info`` spellings remain accepted as fallbacks.  Optional
    ``executed_features``, ``executed_order_parameters``, and
    ``completed_outputs`` entries make citation generation authoritative.
    """
    _require_mapping("run_info", run_info)
    _require_mapping("config", config)
    _require_mapping("statistics", statistics)

    totals = _result_totals(run_info, statistics)
    track_run = _is_track_run(run_info)
    lines: list[str] = [_bold("Basic Information", ansi)]
    _append_basic_information(lines, run_info, config)

    lines.extend(["", _bold("Configuration", ansi)])
    if track_run:
        _append_track_configuration(lines, run_info, config)
    else:
        _append_configuration(lines, run_info, config)

    result_heading = "Tracking Results" if track_run else "Analysis Results"
    lines.extend(["", _bold(result_heading, ansi)])
    _add_field(
        lines,
        "Frames",
        (
            f"requested {totals['requested']}; analyzed {totals['analyzed']}; "
            f"ok {totals['successful']}; failed {totals['failed']}"
        ),
    )
    time_value = (
        f"total {_format_seconds(totals['total_seconds'])}; "
        f"analysis {_format_seconds(totals['analysis_seconds'])}; "
        f"output {_format_seconds(totals['write_seconds'])}"
    )
    if totals["successful"] > 1 and totals["mean_seconds"] is not None:
        time_value += f"; mean {_format_seconds(totals['mean_seconds'])}/frame"
    _add_field(lines, "Time", time_value)

    launch_commands = _validated_render_launch_commands(statistics)
    if launch_commands is not None:
        lines.extend(["", _bold("VMD Rendering", ansi)])
        _add_field(lines, "VMD Tk Console", launch_commands[0])
        _add_field(lines, "Terminal", launch_commands[1])
    if track_run:
        _add_present_field(lines, "Tracks", _lookup(run_info, "track_count"))
    if totals["status"] not in {"completed", "ok", "successful"} or totals["failed"]:
        _add_field(lines, "Status", totals["status"])

    diagnostics = _diagnostic_messages(statistics)
    if diagnostics:
        preview = diagnostics[0]
        if len(diagnostics) > 1:
            preview += f"; +{len(diagnostics) - 1} more (see result files)"
        _add_field(lines, "Warnings", preview)

    citation = build_citation_recommendation(run_info, config, statistics)
    lines.extend(["", _bold("Citation Recommendation", ansi)])
    lines.append(f"  {citation.sentence}")
    lines.append(f"  {citation.publication}")
    lines.append(f"  {citation.github}")
    return "\n".join(lines)


def _diagnostic_messages(statistics: Mapping[str, Any]) -> tuple[str, ...]:
    """Return normalized, de-duplicated diagnostics for the final screen."""
    raw = statistics.get("diagnostic_messages", ())
    if raw is None:
        return ()
    if isinstance(raw, str):
        values: Iterable[Any] = (raw,)
    elif isinstance(raw, Iterable):
        values = raw
    else:
        values = (raw,)
    messages: list[str] = []
    seen: set[str] = set()
    for value in values:
        message = " ".join(str(value).split())
        if message and message not in seen:
            seen.add(message)
            messages.append(message)
    return tuple(messages)


def _validated_render_launch_commands(
    statistics: Mapping[str, Any],
) -> tuple[str, str] | None:
    """Return commands only for one complete, manifest-backed render package."""
    raw = statistics.get(
        "render_script_paths",
        statistics.get("render_script_path", ()),
    )
    if isinstance(raw, (str, Path)):
        candidates = (Path(raw),)
    elif isinstance(raw, Iterable) and not isinstance(raw, (bytes, Mapping)):
        candidates = tuple(Path(item) for item in raw if str(item).strip())
    else:
        return None
    if len(candidates) != 1:
        return None
    script = candidates[0]
    try:
        if not script.is_file() or script.stat().st_size == 0:
            return None
    except OSError:
        return None
    inspection = inspect_render_script(script)
    if inspection is None or not inspection.complete or inspection.source != "manifest":
        return None
    required_roles = {
        reference.role for reference in inspection.references if reference.required
    }
    if not {"topology", "trajectory", "membership"}.issubset(required_roles):
        return None
    return render_launch_commands(inspection)


def print_final_results(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
    *,
    stream: TextIO | None = None,
    ansi: bool | None = None,
) -> None:
    """Write one final terminal screen and flush it.

    ``ansi=None`` enables bold text only for a TTY.  Pass ``ansi=False`` for
    snapshots, redirected output, and terminals without ANSI support.
    """
    target = sys.stdout if stream is None else stream
    use_ansi = _stream_supports_ansi(target) if ansi is None else bool(ansi)
    target.write(render_final_results(run_info, config, statistics, ansi=use_ansi))
    target.write("\n")
    target.flush()


def refresh_terminal(*, stream: TextIO | None = None) -> bool:
    """Clear an interactive terminal before drawing the final screen."""
    target = sys.stdout if stream is None else stream
    if not _stream_supports_ansi(target):
        return False
    target.write("\x1b[2J\x1b[H")
    target.flush()
    return True


def build_citation_sentence(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> str:
    """Return the shared feature-derived citation sentence."""
    return build_citation_recommendation(run_info, config, statistics).sentence


def _started_display(run_info: Mapping[str, Any]) -> Any:
    value = _lookup(run_info, "started_at")
    if value is not _MISSING and value is not None:
        try:
            return format_started(datetime.fromisoformat(str(value)))
        except ValueError:
            pass
    date = _lookup(run_info, "date")
    time = _lookup(run_info, "start_time")
    zone = _lookup(run_info, "time_zone")
    parts = [
        str(item).strip()
        for item in (date, time, zone)
        if item is not _MISSING and item is not None and str(item).strip()
    ]
    return ", ".join(parts) if parts else _MISSING


def _explicit_cage_report_value(
    run_info: Mapping[str, Any], config: Mapping[str, Any]
) -> str:
    value = _first(
        _lookup(run_info, "cage_report_types"),
        _config_value(config, "cage", "report_type"),
        _config_value(config, "cage", "report_types"),
    )
    if value is _MISSING or value is None:
        return ""
    if isinstance(value, str):
        text = value.strip()
        if not text or text.casefold() in {"auto", "all"} or text.casefold().startswith("all detected"):
            return ""
        return text
    if isinstance(value, Iterable) and not isinstance(value, (bytes, Mapping)):
        values = [str(item).strip() for item in value if str(item).strip()]
        return ", ".join(values)
    return str(value)


def _missing_as_empty(value: Any) -> Any:
    return "" if value is _MISSING else value


def _append_basic_information(
    lines: list[str],
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> None:
    _add_present_field(lines, "Started", _started_display(run_info))
    input_value = _first(_lookup(run_info, "input"), _lookup(run_info, "source"))
    if input_value is not _MISSING:
        _add_field(
            lines,
            "Input",
            compact_input_display(
                input_value,
                _missing_as_empty(_lookup(run_info, "input_format")),
            ),
        )
    output_value = _first(
        _lookup(run_info, "output_resolved_path", "output_dir", "output"),
        _lookup(run_info, "result_path"),
    )
    output_config = config.get("output", {})
    metadata = dict(output_config) if isinstance(output_config, Mapping) else {}
    resolved_path = _lookup(run_info, "output_resolved_path")
    if resolved_path is not _MISSING:
        metadata["resolved_path"] = resolved_path
    auto_renamed = _lookup(run_info, "output_auto_renamed")
    if auto_renamed is not _MISSING:
        metadata["auto_renamed"] = _as_bool(auto_renamed, False)
    if output_value is not _MISSING:
        _add_field(lines, "Output", compact_output_display(output_value, metadata))


def _append_configuration(
    lines: list[str],
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> None:
    selector = _first(_lookup(run_info, "engine_selector"), _lookup(config, "mode"), "py")
    _add_field(
        lines,
        "SQQ",
        compact_sqq_display(
            selector,
            version=str(_first(_lookup(run_info, "sqq_version"), __version__)),
        ),
    )
    _add_present_field(lines, "Config file", _lookup(run_info, "config_file"))
    _add_present_field(lines, "Sampling Interval", _lookup(run_info, "sampling_interval"))
    if str(_missing_as_empty(_lookup(run_info, "input_format"))).startswith("lammps-"):
        _add_field(
            lines,
            "LAMMPS",
            (
                f"units {_missing_as_empty(_lookup(run_info, 'lammps_units'))}; "
                f"timestep {_missing_as_empty(_lookup(run_info, 'lammps_timestep'))}; "
                f"style {_missing_as_empty(_lookup(run_info, 'lammps_atom_style'))}; "
                f"type map {_missing_as_empty(_lookup(run_info, 'lammps_type_map_source'))}"
            ),
        )

    group_modes = _lookup(run_info, "graph_mode_by_group")
    if isinstance(group_modes, Mapping) and len(group_modes) > 1:
        group_reasons = _lookup(run_info, "graph_mode_reason_by_group")
        _add_field(
            lines,
            "Graph",
            compact_group_graph_display(
                _first(_lookup(run_info, "graph_mode"), "auto"),
                dict(group_modes),
                dict(group_reasons) if isinstance(group_reasons, Mapping) else {},
            ),
        )
    else:
        _add_present_field(
            lines,
            "Graph",
            compact_graph_display(
                _graph_mode(run_info, config),
                _missing_as_empty(_lookup(run_info, "graph_mode_reason")),
            ),
        )

    _add_field(lines, "Ring", compact_ring_scope(dict(config)))
    _add_field(lines, "Additional search", compact_additional_search(dict(config)))
    cage_report = _explicit_cage_report_value(run_info, config)
    if cage_report:
        _add_field(lines, "Cage report", cage_report)
    _add_present_field(
        lines,
        "Maximum cage face",
        _first(_lookup(run_info, "max_cage_face"), _config_value(config, "cage", "max_face"), _config_value(config, "cage", "max_faces")),
    )
    cage_validation = _first(
        _lookup(run_info, "cage_scientific_validation"),
        _config_value(config, "cage", "scientific_validation"),
    )
    _add_present_field(
        lines,
        "Cage validation",
        _on_off(_as_bool(cage_validation, False)),
    )

    _add_field(lines, "Order parameters", _display_tokens(_selected_order_parameters(run_info, config)))
    _add_field(lines, "Output", _display_tokens(_selected_outputs(run_info, config)))
    topology_groups = _lookup(run_info, "topology_groups")
    if (
        topology_groups is not _MISSING
        and not isinstance(topology_groups, (str, bytes, Mapping))
    ):
        try:
            topology_groups = len(topology_groups)
        except TypeError:
            pass
    if topology_groups is not _MISSING:
        grouping = str(
            _first(
                _lookup(run_info, "topology_grouping", "grouping_policy"),
                "automatic",
            )
        )
        _add_present_field(lines, "Topology groups", f"{topology_groups} [{grouping}]")


def _append_track_configuration(
    lines: list[str],
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> None:
    """Render only settings that have an implemented effect in Track."""
    track_config = config.get("track", {})
    if not isinstance(track_config, Mapping):
        track_config = {}

    selector = _first(_lookup(run_info, "engine_selector"), _lookup(config, "mode"), "py")
    _add_field(
        lines,
        "SQQ",
        compact_sqq_display(
            selector,
            version=str(_first(_lookup(run_info, "sqq_version"), __version__)),
        ),
    )
    _add_present_field(lines, "Config file", _lookup(run_info, "config_file"))

    _add_field(
        lines,
        "Target",
        _first(_lookup(run_info, "target"), track_config.get("target"), "all"),
    )
    jaccard = _first(track_config.get("min_jaccard"), 0.5)
    shared_fraction = _first(track_config.get("min_shared_fraction"), 0.6)
    shared_water = _first(
        track_config.get("min_shared_water"),
        track_config.get("min_shared_waters"),
        3,
    )
    gap = _first(track_config.get("gap_frame"), 0)
    _add_field(
        lines,
        "Track matching",
        (
            f"Jaccard >={jaccard}; shared >={shared_fraction} / "
            f"{shared_water} waters; gap {gap} frames"
        ),
    )
    maximum_distance = _first(
        track_config.get("max_center_distance_nm"),
        "none",
    )
    _add_field(
        lines,
        "Track options",
        (
            f"center distance {maximum_distance} nm; guest tie-break "
            f"{_on_off(_as_bool(track_config.get('guest_tiebreak'), True))}"
        ),
    )
    _add_field(lines, "Output", "Track CSV, sqq-render")


def _is_track_run(run_info: Mapping[str, Any]) -> bool:
    command = _lookup(run_info, "command")
    return command is not _MISSING and str(command).strip().lower() == "track"


def _result_totals(
    run_info: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> dict[str, Any]:
    successful = _nonnegative_int(
        "successful_frames",
        _first(
            _lookup(statistics, "successful_frames", "successful", "frames_ok"),
            _lookup(run_info, "successful_frames", "frames_ok"),
            0,
        ),
    )
    failed = _nonnegative_int(
        "failed_frames",
        _first(
            _lookup(statistics, "failed_frames", "failed", "frames_failed"),
            _lookup(run_info, "failed_frames", "frames_failed"),
            0,
        ),
    )
    analyzed = _nonnegative_int(
        "analyzed_frames",
        _first(
            _lookup(statistics, "analyzed_frames", "analyzed", "frames_total"),
            _lookup(run_info, "analyzed_frames", "frames_total"),
            successful + failed,
        ),
    )
    requested = _nonnegative_int(
        "requested_frames",
        _first(
            _lookup(statistics, "requested_frames", "requested", "frames_requested"),
            _lookup(run_info, "requested_frames", "frames_requested", "selected_frames"),
            analyzed,
        ),
    )

    write_seconds = _nonnegative_float_or_none(
        "write_seconds",
        _first(
            _lookup(statistics, "write_seconds", "output_write_seconds", "summary_write_seconds"),
            _nested_lookup(statistics, "timing", "write_seconds"),
            _nested_lookup(statistics, "summary_write", "total_seconds"),
            _lookup(run_info, "write_seconds", "output_write_seconds"),
            _nested_lookup(run_info, "summary_write", "total_seconds"),
            0.0,
        ),
    )
    total_seconds = _nonnegative_float_or_none(
        "total_seconds",
        _first(
            _lookup(statistics, "total_seconds", "total_elapsed_seconds", "elapsed_seconds"),
            _nested_lookup(statistics, "timing", "total_seconds"),
            _lookup(run_info, "total_seconds", "total_elapsed_seconds", "elapsed_seconds"),
        ),
    )
    analysis_value = _first(
        _lookup(statistics, "analysis_seconds", "analysis_time_seconds"),
        _nested_lookup(statistics, "timing", "analysis_seconds"),
        _lookup(run_info, "analysis_seconds", "analysis_time_seconds"),
    )
    if analysis_value is _MISSING:
        analysis_seconds = (
            max(0.0, total_seconds - write_seconds)
            if total_seconds is not None and write_seconds is not None
            else None
        )
    else:
        analysis_seconds = _nonnegative_float_or_none("analysis_seconds", analysis_value)
    if total_seconds is None and analysis_seconds is not None and write_seconds is not None:
        total_seconds = analysis_seconds + write_seconds

    mean_value = _first(
        _lookup(statistics, "mean_seconds", "mean_seconds_per_successful_frame"),
        _nested_lookup(statistics, "timing", "mean_seconds"),
    )
    if mean_value is _MISSING:
        mean_seconds = (
            analysis_seconds / successful
            if analysis_seconds is not None and successful > 0
            else None
        )
    else:
        mean_seconds = _nonnegative_float_or_none("mean_seconds", mean_value)

    raw_status = str(
        _first(_lookup(statistics, "status", "run_status"), _lookup(run_info, "status"), "unknown")
    ).strip().lower()
    status = raw_status or "unknown"
    if status in {"complete", "completed", "success", "successful", "ok"} and failed:
        status = "completed with failures"
    elif status in {"complete", "success", "successful", "ok"}:
        status = "completed"

    result_path = _first(
        _lookup(statistics, "result_path", "output_path", "output_dir"),
        _lookup(run_info, "result_path", "output_path", "output_dir"),
        "N/A",
    )
    if isinstance(result_path, Path):
        result_path = str(result_path)

    return {
        "requested": requested,
        "analyzed": analyzed,
        "successful": successful,
        "failed": failed,
        "total_seconds": total_seconds,
        "analysis_seconds": analysis_seconds,
        "write_seconds": write_seconds,
        "mean_seconds": mean_seconds,
        "status": status,
        "result_path": result_path,
    }


def _selected_order_parameters(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> tuple[str, ...]:
    return _tokens(
        _first(
            _lookup(run_info, "order_parameters"),
            _config_value(config, "order_parameter", "enabled"),
            _config_value(config, "order", "parameter"),
            _config_value(config, "order", "parameters"),
        )
    )


def _selected_outputs(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> tuple[str, ...]:
    return _tokens(
        _first(
            _lookup(run_info, "output_types"),
            _config_value(config, "output", "type"),
            _config_value(config, "output", "types"),
        )
    )


def _graph_mode(run_info: Mapping[str, Any], config: Mapping[str, Any]) -> Any:
    display = _lookup(run_info, "graph_mode_display")
    if display is not _MISSING and str(display).strip():
        return display
    requested = _first(
        _lookup(run_info, "graph_mode"),
        _config_value(config, "graph", "mode"),
        _config_value(config, "graph", "bond_mode"),
    )
    effective = _first(
        _lookup(run_info, "effective_graph_modes", "effective_graph_mode"),
        _config_value(config, "graph", "effective_mode"),
        _config_value(config, "graph", "effective_bond_mode"),
    )
    if requested is _MISSING:
        return effective
    if str(requested).strip().lower() == "auto" and effective is not _MISSING and str(effective).strip():
        return f"auto -> {effective}"
    return requested


def _tokens(value: Any) -> tuple[str, ...]:
    if value is _MISSING or value is None:
        return ()
    raw: list[Any]
    if isinstance(value, str):
        raw = [item for item in _TOKEN_SPLIT.split(value.strip()) if item]
    elif isinstance(value, Mapping):
        raw = [key for key, enabled in value.items() if _as_bool(enabled, False)]
    elif isinstance(value, Iterable):
        raw = list(value)
    else:
        raw = [value]
    result: list[str] = []
    seen: set[str] = set()
    for item in raw:
        text = str(item).strip()
        normalized = text.lower()
        if not text or normalized in {"none", "null", "off", "false", "[]"} or normalized in seen:
            continue
        seen.add(normalized)
        result.append(text)
    return tuple(result)


def _display_tokens(values: tuple[str, ...]) -> str:
    return ", ".join(values) if values else "none"


def _format_seconds(value: float | None) -> str:
    return "N/A" if value is None else f"{value:.3f} s"


def _on_off(value: bool) -> str:
    return "on" if value else "off"


def _as_bool(value: Any, default: bool) -> bool:
    if value is _MISSING or value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    text = str(value).strip().lower()
    if text in _TRUE_TEXT:
        return True
    if text in _FALSE_TEXT:
        return False
    return default


def _nonnegative_int(name: str, value: Any) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if number < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return number


def _nonnegative_float_or_none(name: str, value: Any) -> float | None:
    if value is _MISSING or value is None or str(value).strip() == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite nonnegative number") from exc
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{name} must be a finite nonnegative number")
    return number


def _lookup(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return _MISSING


def _nested_lookup(mapping: Mapping[str, Any], *keys: str) -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, Mapping) or key not in current:
            return _MISSING
        current = current[key]
    return current


def _config_value(config: Mapping[str, Any], section: str, key: str) -> Any:
    return _nested_lookup(config, section, key)


def _first(*values: Any) -> Any:
    for value in values:
        if value is not _MISSING and value is not None:
            return value
    return _MISSING


def _add_present_field(lines: list[str], label: str, value: Any) -> None:
    if value is _MISSING or value is None or str(value).strip() == "":
        return
    _add_field(lines, label, value)


def _add_field(lines: list[str], label: str, value: Any) -> None:
    lines.append(f"  {label:<{_LABEL_WIDTH}}: {_terminal_text(value)}")


def _terminal_text(value: Any) -> str:
    if value is _MISSING or value is None:
        return "N/A"
    if isinstance(value, bool):
        text = str(value).lower()
    elif isinstance(value, (list, tuple, set)):
        text = ", ".join(str(item) for item in value)
    else:
        text = str(value)
    return "".join(character if character >= " " and character != "\x7f" else " " for character in text)


def _bold(text: str, ansi: bool) -> str:
    return f"{_ANSI_BOLD}{text}{_ANSI_RESET}" if ansi else text


def _stream_supports_ansi(stream: TextIO) -> bool:
    isatty = getattr(stream, "isatty", None)
    try:
        return bool(isatty()) if callable(isatty) else False
    except OSError:
        return False


def _require_mapping(name: str, value: Any) -> None:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
