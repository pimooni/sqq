from __future__ import annotations

"""Feature-derived citation recommendations shared by all output surfaces."""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import re
from typing import Any


PUBLICATION_LINE = (
    "Publication: J. Pang & Q. Sun, SQQ: Python Joint Toolkit for "
    "Water-Shell Topology Analysis, in submission."
)
GITHUB_LINE = "GitHub     : https://github.com/pimooni/sqq"

_MISSING = object()
_TOKEN_SPLIT = re.compile(r"[,;/\s]+")
_TRUE_TEXT = {"1", "true", "yes", "on", "enabled", "complete", "completed"}
_FALSE_TEXT = {"", "0", "false", "no", "off", "disabled", "none", "null", "auto"}


@dataclass(frozen=True, slots=True)
class CitationRecommendation:
    """One recommendation rendered identically in terminal and reports."""

    sentence: str
    publication: str = PUBLICATION_LINE
    github: str = GITHUB_LINE


def build_citation_recommendation(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> CitationRecommendation:
    """Build a recommendation from work that completed successfully.

    ``statistics['executed_features']`` is authoritative when supplied.  An
    enabled analysis that completed successfully may be cited even when it
    found zero matching structures; this function never claims detection.
    """
    _require_mapping("run_info", run_info)
    _require_mapping("config", config)
    _require_mapping("statistics", statistics)

    completed = _successful_frames(run_info, statistics) > 0
    features: list[str] = []

    if _feature_ran(
        "water_network",
        run_info,
        statistics,
        default=completed and _configured_enabled(config, "graph", True),
    ):
        features.append("water-network analysis")
    if _feature_ran(
        "ring_topology",
        run_info,
        statistics,
        aliases=("ring", "rings"),
        default=completed and _configured_enabled(config, "ring", True),
    ):
        features.append("ring-topology analysis")

    cage_ran = _feature_ran(
        "cage_topology",
        run_info,
        statistics,
        aliases=("cage", "cages"),
        default=completed and _configured_enabled(config, "cage", True),
    )
    if cage_ran:
        features.append("cage-topology analysis")
    if _feature_ran(
        "half_cage",
        run_info,
        statistics,
        default=completed and _effective_enabled(run_info, config, "find_half", "half_cage"),
    ):
        features.append("half-cage analysis")
    if _feature_ran(
        "quasi_cage",
        run_info,
        statistics,
        default=completed and _effective_enabled(run_info, config, "find_quasi", "quasi_cage"),
    ):
        features.append("quasi-cage analysis")
    if _feature_ran(
        "cage_isomer",
        run_info,
        statistics,
        aliases=("isomer", "cage_isomers"),
        default=cage_ran,
    ):
        features.append("cage-isomer analysis")
    if _feature_ran(
        "cage_occupancy",
        run_info,
        statistics,
        aliases=("occupancy",),
        default=completed and _occupancy_evaluated(run_info, statistics),
    ):
        features.append("cage-occupancy analysis")

    parameters = _executed_order_parameters(run_info, config, statistics)
    if completed and parameters:
        features.append(f"{_slash_join(parameters)} order-parameter analysis")
    if _feature_ran(
        "hydrate_phase_domain",
        run_info,
        statistics,
        aliases=("hydrate_cluster", "cluster", "phase_domain"),
        default=completed
        and _effective_enabled(run_info, config, "find_cluster", "hydrate_cluster"),
    ):
        features.append("hydrate phase/domain analysis")
    if _feature_ran(
        "vmd_rendering",
        run_info,
        statistics,
        aliases=("sqq_render", "render", "vmd"),
        default=_completed_render_output(run_info, config, statistics),
    ):
        features.append("VMD rendering")
    if _feature_ran(
        "cage_tracking",
        run_info,
        statistics,
        aliases=("tracking", "track"),
        default=completed and _track_enabled(config),
    ):
        features.append("cage-tracking analysis")
    if _feature_ran(
        "cage_lifetime",
        run_info,
        statistics,
        aliases=("lifetime", "lifetimes"),
        default=completed and _lifetime_enabled(config),
    ):
        features.append("cage-lifetime analysis")

    if features:
        sentence = (
            "Please cite SQQ when using results from this run, including "
            f"{_natural_join(features)}."
        )
    else:
        sentence = "Please cite SQQ when using results from this run."
    return CitationRecommendation(sentence=sentence)


def build_citation_sentence(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> str:
    """Compatibility helper returning only the recommended sentence."""
    return build_citation_recommendation(run_info, config, statistics).sentence


def completed_citation_evidence(
    config: Mapping[str, Any],
    *,
    successful_frames: int,
    completed_outputs: Iterable[Any] = (),
    track: bool = False,
) -> dict[str, Any]:
    """Return authoritative feature evidence for a completed workflow."""
    successful = max(0, int(successful_frames))
    ran = successful > 0
    outputs = _tokens(completed_outputs)
    order_parameters = _tokens(
        _first(
            _nested(config, "order_parameter", "enabled"),
            _nested(config, "order", "parameters"),
            _nested(config, "order", "parameter"),
        )
    )
    return {
        "successful_frames": successful,
        "executed_features": {
            "water_network": ran,
            "ring_topology": ran,
            "cage_topology": ran,
            "half_cage": ran and _as_bool(_nested(config, "half_cage", "enabled"), False),
            "quasi_cage": ran and _as_bool(_nested(config, "quasi_cage", "enabled"), False),
            "cage_isomer": ran,
            "cage_occupancy": ran,
            "hydrate_phase_domain": ran
            and _as_bool(_nested(config, "hydrate_cluster", "enabled"), False),
            "vmd_rendering": "sqq-render" in outputs,
            "cage_tracking": ran and bool(track),
            "cage_lifetime": ran and bool(track),
        },
        "executed_order_parameters": order_parameters if ran else (),
        "completed_outputs": outputs,
    }


def _feature_ran(
    name: str,
    run_info: Mapping[str, Any],
    statistics: Mapping[str, Any],
    *,
    aliases: tuple[str, ...] = (),
    default: bool,
) -> bool:
    names = tuple(_normalize_feature_name(item) for item in (name, *aliases))
    executed = statistics.get("executed_features", _MISSING)
    if isinstance(executed, Mapping):
        normalized = {
            _normalize_feature_name(key): value for key, value in executed.items()
        }
        return any(_as_bool(normalized[item], False) for item in names if item in normalized)
    if executed is not _MISSING and not isinstance(executed, (str, bytes)):
        normalized = {_normalize_feature_name(item) for item in executed}
        return any(item in normalized for item in names)

    direct_names = tuple(
        value
        for item in names
        for value in (item, f"{item}_analyzed", f"{item}_executed")
    )
    for source in (statistics, run_info):
        for key in direct_names:
            if key in source:
                return _as_bool(source[key], False)
    return bool(default)


def _successful_frames(
    run_info: Mapping[str, Any], statistics: Mapping[str, Any]
) -> int:
    for source in (statistics, run_info):
        for key in ("successful_frames", "frames_ok", "completed_files"):
            if key in source:
                try:
                    return max(0, int(source[key]))
                except (TypeError, ValueError):
                    return 0
    return 0


def _executed_order_parameters(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> tuple[str, ...]:
    for key in ("executed_order_parameters", "order_parameters_executed"):
        if key in statistics:
            return _tokens(statistics[key])
    executed = statistics.get("executed_features", _MISSING)
    if isinstance(executed, Mapping):
        names = {
            _normalize_feature_name(key)
            for key, enabled in executed.items()
            if _as_bool(enabled, False)
        }
        if not names & {"order_parameter", "order_parameters"}:
            return ()
    elif executed is not _MISSING:
        names = {_normalize_feature_name(item) for item in executed}
        if not names & {"order_parameter", "order_parameters"}:
            return ()
    return _tokens(
        _first(
            run_info.get("order_parameters", _MISSING),
            _nested(config, "order_parameter", "enabled"),
            _nested(config, "order", "parameters"),
            _nested(config, "order", "parameter"),
        )
    )


def _completed_render_output(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    statistics: Mapping[str, Any],
) -> bool:
    completed = _first(
        statistics.get("completed_outputs", _MISSING),
        statistics.get("written_outputs", _MISSING),
        run_info.get("completed_outputs", _MISSING),
        run_info.get("written_outputs", _MISSING),
    )
    if completed is not _MISSING:
        names = {_normalize_feature_name(item) for item in _tokens(completed)}
        return bool(names & {"sqq_render", "vmd", "vmd_render", "vmd_rendering"})
    status = str(_first(statistics.get("status", _MISSING), run_info.get("status", ""))).lower()
    if status not in {"completed", "completed with failures", "partial", "ok"}:
        return False
    outputs = _tokens(
        _first(
            run_info.get("output_types", _MISSING),
            _nested(config, "output", "types"),
            _nested(config, "output", "type"),
        )
    )
    return "sqq_render" in {_normalize_feature_name(item) for item in outputs}


def _occupancy_evaluated(
    run_info: Mapping[str, Any], statistics: Mapping[str, Any]
) -> bool:
    explicit = _first(
        statistics.get("occupancy_evaluated", _MISSING),
        statistics.get("cage_occupancy_evaluated", _MISSING),
        run_info.get("occupancy_evaluated", _MISSING),
        run_info.get("cage_occupancy_evaluated", _MISSING),
    )
    if explicit is not _MISSING:
        return _as_bool(explicit, False)
    guests = _first(
        statistics.get("guest_molecules", _MISSING),
        statistics.get("n_guests", _MISSING),
        statistics.get("has_selected_guests", _MISSING),
        run_info.get("n_guests", _MISSING),
    )
    try:
        return guests is not _MISSING and float(guests) > 0
    except (TypeError, ValueError):
        return _as_bool(guests, False)


def _effective_enabled(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    run_key: str,
    section: str,
) -> bool:
    if run_key in run_info:
        return _as_bool(run_info[run_key], False)
    return _as_bool(_nested(config, section, "enabled"), False)


def _configured_enabled(config: Mapping[str, Any], section: str, default: bool) -> bool:
    value = _nested(config, section, "enabled")
    return default if value is _MISSING else _as_bool(value, default)


def _track_enabled(config: Mapping[str, Any]) -> bool:
    return _as_bool(_nested(config, "track", "enabled"), False)


def _lifetime_enabled(config: Mapping[str, Any]) -> bool:
    value = _nested(config, "track", "lifetime")
    if isinstance(value, Mapping):
        value = value.get("enabled", _MISSING)
    return _as_bool(value, False)


def _tokens(value: Any) -> tuple[str, ...]:
    if value is _MISSING or value is None:
        return ()
    if isinstance(value, str):
        raw: Iterable[Any] = (item for item in _TOKEN_SPLIT.split(value.strip()) if item)
    elif isinstance(value, Mapping):
        raw = (key for key, enabled in value.items() if _as_bool(enabled, False))
    elif isinstance(value, Iterable):
        raw = value
    else:
        raw = (value,)
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


def _slash_join(values: tuple[str, ...]) -> str:
    return "/".join(_order_parameter_label(value) for value in values)


def _order_parameter_label(value: str) -> str:
    text = value.strip()
    normalized = text.lower().replace("_", "")
    if normalized.startswith("q") and normalized[1:].isdigit():
        return f"Q{normalized[1:]}"
    if normalized.startswith(("mcg", "dhop")) or normalized in {"f3", "f4"}:
        return normalized.upper()
    return text


def _natural_join(values: list[str]) -> str:
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return f"{values[0]} and {values[1]}"
    return f"{', '.join(values[:-1])}, and {values[-1]}"


def _normalize_feature_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def _nested(mapping: Mapping[str, Any], *keys: str) -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, Mapping) or key not in current:
            return _MISSING
        current = current[key]
    return current


def _first(*values: Any) -> Any:
    for value in values:
        if value is not _MISSING and value is not None:
            return value
    return _MISSING


def _as_bool(value: Any, default: bool) -> bool:
    if value is _MISSING or value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in _TRUE_TEXT:
        return True
    if text in _FALSE_TEXT:
        return False
    return default


def _require_mapping(name: str, value: Any) -> None:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
