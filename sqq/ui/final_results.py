from __future__ import annotations

"""Build the final SQQ terminal results screen.

This module deliberately has no dependency on the Analyze workflow, pandas, or
the output writers.  Callers pass the resolved run metadata, resolved
configuration, and a small mapping of final statistics.  Values in
``statistics`` take precedence over their legacy ``run_info`` equivalents.
"""

import math
import re
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, TextIO

from .. import __version__


PUBLICATION_LINE = (
    "Publication: J. Pang & Q. Sun, SQQ: Python Joint Toolkit for "
    "Water-Shell Topology Analysis, in submission."
)
GITHUB_LINE = "Github     : https://github.com/pimooni/sqq"

_ANSI_BOLD = "\x1b[1m"
_ANSI_RESET = "\x1b[0m"
_LABEL_WIDTH = 32
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
    _append_basic_information(lines, run_info)

    lines.extend(["", _bold("Configuration", ansi)])
    if track_run:
        _append_track_configuration(lines, run_info, config)
    else:
        _append_configuration(lines, run_info, config)

    diagnostics = _diagnostic_messages(statistics)
    if diagnostics:
        lines.extend(["", _bold("Diagnostics", ansi)])
        _add_field(lines, "Warnings", len(diagnostics))
        lines.extend(f"  Warning: {message}" for message in diagnostics)

    result_heading = "Tracking Results" if track_run else "Analysis Results"
    lines.extend(["", _bold(result_heading, ansi)])
    _add_field(lines, "Requested frames", totals["requested"])
    _add_field(lines, "Analyzed frames", totals["analyzed"])
    _add_field(lines, "Successful frames", totals["successful"])
    _add_field(lines, "Failed frames", totals["failed"])
    _add_field(lines, "Total elapsed time", _format_seconds(totals["total_seconds"]))
    _add_field(lines, "Analysis time", _format_seconds(totals["analysis_seconds"]))
    _add_field(lines, "Output-writing time", _format_seconds(totals["write_seconds"]))
    _add_field(
        lines,
        "Mean time / successful frame",
        _format_seconds(totals["mean_seconds"]),
    )
    _add_field(lines, "Run status", totals["status"])
    _add_field(lines, "Result path", totals["result_path"])

    sentence = build_citation_sentence(run_info, config, statistics)
    lines.extend(
        [
            "",
            _bold("Citation Recommendation", ansi),
            _bold(sentence, ansi),
            PUBLICATION_LINE,
            GITHUB_LINE,
        ]
    )
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
    """Build a deterministic recommendation from work that actually ran.

    An enabled search is described as an analysis even when it found zero
    matching objects.  The sentence intentionally never claims that a
    structure was identified.  ``executed_features`` may be a mapping of
    feature names to booleans or an authoritative iterable of feature names.
    """
    _require_mapping("run_info", run_info)
    _require_mapping("config", config)
    _require_mapping("statistics", statistics)

    totals = _result_totals(run_info, statistics)
    completed_scientific_work = totals["successful"] > 0
    features: list[str] = []

    if _feature_ran(
        "water_network",
        run_info,
        config,
        statistics,
        default=completed_scientific_work and _configured_enabled(config, ("graph", "enabled"), True),
    ):
        features.append("water-network analysis")
    if _feature_ran(
        "ring_topology",
        run_info,
        config,
        statistics,
        aliases=("ring", "rings"),
        default=completed_scientific_work and _configured_enabled(config, ("ring", "enabled"), True),
    ):
        features.append("ring-topology analysis")

    cage_ran = _feature_ran(
        "cage_topology",
        run_info,
        config,
        statistics,
        aliases=("cage", "cages"),
        default=completed_scientific_work and _configured_enabled(config, ("cage", "enabled"), True),
    )
    if cage_ran:
        features.append("cage-topology analysis")

    if _feature_ran(
        "half_cage",
        run_info,
        config,
        statistics,
        default=completed_scientific_work
        and _effective_enabled(run_info, config, "find_half", "half_cage"),
    ):
        features.append("half-cage analysis")
    if _feature_ran(
        "quasi_cage",
        run_info,
        config,
        statistics,
        default=completed_scientific_work
        and _effective_enabled(run_info, config, "find_quasi", "quasi_cage"),
    ):
        features.append("quasi-cage analysis")
    if _feature_ran(
        "cage_isomer",
        run_info,
        config,
        statistics,
        aliases=("isomer", "cage_isomers"),
        default=cage_ran,
    ):
        features.append("cage-isomer analysis")
    if _feature_ran(
        "cage_occupancy",
        run_info,
        config,
        statistics,
        aliases=("occupancy",),
        default=completed_scientific_work
        and _explicit_or_inferred_occupancy(run_info, statistics),
    ):
        features.append("cage-occupancy analysis")

    parameters = _executed_order_parameters(run_info, config, statistics)
    if completed_scientific_work and parameters:
        features.append(f"{_slash_join(parameters)} order-parameter analysis")

    if _feature_ran(
        "hydrate_phase_domain",
        run_info,
        config,
        statistics,
        aliases=("hydrate_cluster", "cluster", "phase_domain"),
        default=completed_scientific_work
        and _effective_enabled(run_info, config, "find_cluster", "hydrate_cluster"),
    ):
        features.append("hydrate phase/domain analysis")

    if _feature_ran(
        "vmd_rendering",
        run_info,
        config,
        statistics,
        aliases=("sqq_render", "render", "vmd"),
        default=_completed_render_output(run_info, config, statistics, totals["status"]),
    ):
        features.append("VMD rendering")
    if _feature_ran(
        "cage_tracking",
        run_info,
        config,
        statistics,
        aliases=("tracking", "track"),
        default=completed_scientific_work and _future_feature_enabled(config, "track"),
    ):
        features.append("cage-tracking analysis")
    if _feature_ran(
        "cage_lifetime",
        run_info,
        config,
        statistics,
        aliases=("lifetime", "lifetimes"),
        default=completed_scientific_work and _future_feature_enabled(config, "lifetime"),
    ):
        features.append("cage-lifetime analysis")

    if not features:
        return "Please cite SQQ when using results from this run."
    return f"Please cite SQQ when using results from this run, including {_natural_join(features)}."


def _append_basic_information(lines: list[str], run_info: Mapping[str, Any]) -> None:
    fields = (
        ("Date", _lookup(run_info, "date")),
        ("Start time", _lookup(run_info, "start_time")),
        ("Finish time", _lookup(run_info, "finish_time")),
        ("Time zone", _lookup(run_info, "time_zone")),
        ("Working directory", _lookup(run_info, "working_dir")),
        ("Input", _lookup(run_info, "input")),
        ("Input format", _lookup(run_info, "input_format")),
        ("Matched files", _lookup(run_info, "matched_files")),
    )
    for label, value in fields:
        _add_present_field(lines, label, value)

    matched = _optional_int(_lookup(run_info, "matched_files"))
    first_file = _lookup(run_info, "first_file")
    last_file = _lookup(run_info, "last_file")
    if matched is not None and matched > 1:
        _add_present_field(lines, "First file", first_file)
        _add_present_field(lines, "Last file", last_file)
    elif first_file is not _MISSING:
        _add_present_field(lines, "Current file", first_file)


def _append_configuration(
    lines: list[str],
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> None:
    engine = _engine_name(run_info, config)
    _add_field(lines, "SQQ version", _first(_lookup(run_info, "sqq_version"), __version__))
    _add_field(lines, "SQQ engine", engine)
    _add_present_field(lines, "Config file", _lookup(run_info, "config_file"))
    _add_present_field(lines, "Topology", _lookup(run_info, "topology"))
    _add_present_field(lines, "Sampling interval", _lookup(run_info, "sampling_interval"))

    group_modes = _lookup(run_info, "graph_mode_by_group")
    if isinstance(group_modes, Mapping) and len(group_modes) > 1:
        for label, mode in group_modes.items():
            _add_field(lines, f"Graph mode ({label})", mode)
    else:
        _add_present_field(lines, "Graph mode", _graph_mode(run_info, config))

    _add_present_field(
        lines,
        "Search sizes",
        _first(_config_value(config, "ring", "size"), _config_value(config, "ring", "sizes")),
    )
    _add_present_field(
        lines,
        "Ring definition",
        _first(_config_value(config, "ring", "definition"), _config_value(config, "ring", "chordless")),
    )
    if engine != "sqq-cpp":
        _add_present_field(
            lines,
            "Ring report sizes",
            _first(
                _lookup(run_info, "ring_report_sizes"),
                _config_value(config, "ring", "report_size"),
                _config_value(config, "ring", "report_sizes"),
            ),
        )
        _add_field(
            lines,
            "Find half",
            _on_off(_effective_enabled(run_info, config, "find_half", "half_cage")),
        )
        _add_field(
            lines,
            "Find quasi",
            _on_off(_effective_enabled(run_info, config, "find_quasi", "quasi_cage")),
        )

    _add_present_field(
        lines,
        "Cage report types",
        _first(
            _lookup(run_info, "cage_report_types"),
            _config_value(config, "cage", "report_type"),
            _config_value(config, "cage", "report_types"),
        ),
    )
    _add_present_field(
        lines,
        "Maximum cage face",
        _first(_lookup(run_info, "max_cage_face"), _config_value(config, "cage", "max_face"), _config_value(config, "cage", "max_faces")),
    )
    if engine != "sqq-cpp":
        _add_field(
            lines,
            "Find cluster",
            _on_off(_effective_enabled(run_info, config, "find_cluster", "hydrate_cluster")),
        )

    _add_field(lines, "Order parameters", _display_tokens(_selected_order_parameters(run_info, config)))
    _add_field(lines, "Output types", _display_tokens(_selected_outputs(run_info, config)))
    _add_present_field(
        lines,
        "Output layout",
        _first(_lookup(run_info, "output_layout"), _config_value(config, "output", "structure_layout")),
    )
    _add_present_field(lines, "Worker policy", _lookup(run_info, "worker_policy"))
    _add_present_field(lines, "Parallel backend", _lookup(run_info, "parallel_backend"))
    _add_present_field(lines, "Math threads per worker", _lookup(run_info, "math_threads"))
    _add_present_field(
        lines,
        "Workers",
        _first(_lookup(run_info, "workers"), _config_value(config, "parallel", "worker"), _config_value(config, "parallel", "workers")),
    )
    topology_groups = _lookup(run_info, "topology_groups")
    if (
        topology_groups is not _MISSING
        and not isinstance(topology_groups, (str, bytes, Mapping))
    ):
        try:
            topology_groups = len(topology_groups)
        except TypeError:
            pass
    _add_present_field(lines, "Topology groups", topology_groups)
    _add_present_field(lines, "Grouping policy", _lookup(run_info, "grouping_policy"))


def _append_track_configuration(
    lines: list[str],
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
) -> None:
    """Render only settings that have an implemented effect in Track."""
    track_config = config.get("track", {})
    if not isinstance(track_config, Mapping):
        track_config = {}

    _add_field(lines, "SQQ version", _first(_lookup(run_info, "sqq_version"), __version__))
    _add_field(lines, "SQQ engine", _engine_name(run_info, config))
    _add_present_field(lines, "Config file", _lookup(run_info, "config_file"))

    source = _lookup(run_info, "source")
    if source is not _MISSING and source is not None and str(source).strip():
        _add_field(lines, "Source", source)
    else:
        _add_present_field(lines, "Input", _lookup(run_info, "input"))

    _add_field(
        lines,
        "Target",
        _first(_lookup(run_info, "target"), track_config.get("target"), "all"),
    )
    _add_field(
        lines,
        "Minimum Jaccard",
        _first(track_config.get("min_jaccard"), 0.5),
    )
    _add_field(
        lines,
        "Minimum shared fraction",
        _first(track_config.get("min_shared_fraction"), 0.6),
    )
    _add_field(
        lines,
        "Minimum shared water",
        _first(
            track_config.get("min_shared_water"),
            track_config.get("min_shared_waters"),
            3,
        ),
    )
    maximum_distance = _first(
        track_config.get("max_center_distance_nm"),
        "none",
    )
    _add_field(lines, "Maximum center distance (nm)", maximum_distance)
    _add_field(
        lines,
        "Gap frames",
        _first(track_config.get("gap_frame"), 0),
    )
    _add_field(
        lines,
        "Guest tie-break",
        _on_off(_as_bool(track_config.get("guest_tiebreak"), True)),
    )
    _add_field(lines, "Output types", "Track CSV, sqq-render")


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


def _feature_ran(
    name: str,
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
    *,
    aliases: tuple[str, ...] = (),
    default: bool,
) -> bool:
    del config  # Kept in the signature so every feature resolver has one contract.
    names = tuple(_normalize_feature_name(item) for item in (name, *aliases))
    executed = _lookup(statistics, "executed_features")
    if isinstance(executed, Mapping):
        normalized = {
            _normalize_feature_name(str(key)): value for key, value in executed.items()
        }
        for candidate in names:
            if candidate in normalized:
                return _as_bool(normalized[candidate], False)
        return False
    elif executed is not _MISSING and not isinstance(executed, (str, bytes)):
        try:
            normalized_names = {_normalize_feature_name(str(item)) for item in executed}
        except TypeError as exc:
            raise TypeError("statistics['executed_features'] must be a mapping or iterable") from exc
        return any(candidate in normalized_names for candidate in names)

    direct_names = tuple(
        item
        for candidate in names
        for item in (candidate, f"{candidate}_analyzed", f"{candidate}_executed")
    )
    direct = _lookup(statistics, *direct_names)
    if direct is not _MISSING:
        return _as_bool(direct, False)
    run_value = _lookup(run_info, *direct_names)
    if run_value is not _MISSING:
        return _as_bool(run_value, False)
    return bool(default)


def _executed_order_parameters(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> tuple[str, ...]:
    explicit = _lookup(statistics, "executed_order_parameters", "order_parameters_executed")
    if explicit is not _MISSING:
        return _tokens(explicit)
    feature_set = _lookup(statistics, "executed_features")
    if isinstance(feature_set, Mapping):
        normalized = {
            _normalize_feature_name(str(key))
            for key, enabled in feature_set.items()
            if _as_bool(enabled, False)
        }
        if not normalized & {"order_parameter", "order_parameters"}:
            return ()
    elif feature_set is not _MISSING:
        normalized = {_normalize_feature_name(str(item)) for item in feature_set}
        if not normalized & {"order_parameter", "order_parameters"}:
            return ()
    return _selected_order_parameters(run_info, config)


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


def _completed_render_output(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
    status: str,
) -> bool:
    completed = _first(
        _lookup(statistics, "completed_outputs", "written_outputs"),
        _lookup(run_info, "completed_outputs", "written_outputs"),
    )
    if completed is not _MISSING:
        outputs = {_normalize_feature_name(item) for item in _tokens(completed)}
        return bool(outputs & {"sqq_render", "vmd", "vmd_render", "vmd_rendering"})
    if status not in {"completed", "completed with failures", "partial"}:
        return False
    outputs = {_normalize_feature_name(item) for item in _selected_outputs(run_info, config)}
    return "sqq_render" in outputs


def _explicit_or_inferred_occupancy(
    run_info: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> bool:
    explicit = _first(
        _lookup(statistics, "occupancy_evaluated", "cage_occupancy_evaluated"),
        _lookup(run_info, "occupancy_evaluated", "cage_occupancy_evaluated"),
    )
    if explicit is not _MISSING:
        return _as_bool(explicit, False)
    guests = _first(
        _lookup(statistics, "guest_molecules", "guest_count", "n_guests", "has_selected_guests"),
        _lookup(run_info, "guest_molecules", "guest_count", "n_guests", "has_selected_guests"),
    )
    if guests is _MISSING:
        return False
    if isinstance(guests, bool):
        return guests
    try:
        return float(guests) > 0
    except (TypeError, ValueError):
        return _as_bool(guests, False)


def _effective_enabled(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    run_key: str,
    config_section: str,
) -> bool:
    run_value = _lookup(run_info, run_key)
    if run_value is not _MISSING:
        return _as_bool(run_value, False)
    return _as_bool(_config_value(config, config_section, "enabled"), False)


def _future_feature_enabled(config: Mapping[str, Any], feature: str) -> bool:
    direct = _config_value(config, feature, "enabled")
    if direct is not _MISSING:
        return _as_bool(direct, False)
    track = _config_value(config, "track", feature)
    if isinstance(track, Mapping):
        return _as_bool(_lookup(track, "enabled"), False)
    return _as_bool(track, False)


def _configured_enabled(
    config: Mapping[str, Any],
    path: tuple[str, str],
    default: bool,
) -> bool:
    value = _config_value(config, *path)
    return default if value is _MISSING else _as_bool(value, default)


def _engine_name(run_info: Mapping[str, Any], config: Mapping[str, Any]) -> str:
    value = _first(
        _lookup(run_info, "sqq_engine", "engine"),
        _lookup(config, "engine", "mode"),
        "sqq-py",
    )
    text = str(value).strip().lower()
    if text in {"99", "cpp", "sqq-cpp"} or "cpp" in text:
        return "sqq-cpp"
    return "sqq-py"


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


def _slash_join(values: tuple[str, ...]) -> str:
    return "/".join(_order_parameter_label(value) for value in values)


def _order_parameter_label(value: str) -> str:
    text = value.strip()
    lowered = text.lower().replace("_", "")
    if lowered.startswith("q") and lowered[1:].isdigit():
        return f"Q{lowered[1:]}"
    if lowered.startswith("mcg"):
        return lowered.upper()
    if lowered.startswith("dhop"):
        return lowered.upper()
    if lowered in {"f3", "f4"}:
        return lowered.upper()
    return text


def _natural_join(values: list[str]) -> str:
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return f"{values[0]} and {values[1]}"
    return f"{', '.join(values[:-1])}, and {values[-1]}"


def _normalize_feature_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


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


def _optional_int(value: Any) -> int | None:
    if value is _MISSING or value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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
