"""Engine capability adjustment."""

from __future__ import annotations

from copy import deepcopy
from typing import Any
import warnings

from .defaults import DEFAULT_CONFIG, DEFAULT_MODE
from .migrate import legacy_enabled
from .presets import is_cpp_mode
from .validation import (
    normalize_cpp_order_parameters,
    normalize_cpp_output_types,
    validate_cpp_cli,
)


def _auto_toggle(value: Any, automatic: bool) -> bool:
    if isinstance(value, str) and value.strip().lower() == "auto":
        return automatic
    return legacy_enabled(value)


def _explicit(config: dict[str, Any] | None, section: str) -> dict[str, Any]:
    value = config.get(section, {}) if isinstance(config, dict) else {}
    return value if isinstance(value, dict) else {}


def _explicitly_enabled(value: Any) -> bool:
    if value is None or (isinstance(value, str) and value.strip().lower() == "auto"):
        return False
    return legacy_enabled(value)


def _section_customized(
    explicit: dict[str, Any], defaults: dict[str, Any]
) -> bool:
    for key, value in explicit.items():
        default = defaults.get(key)
        if key == "enabled":
            if isinstance(value, str) and value.strip().lower() == "auto":
                continue
            try:
                if legacy_enabled(value) == legacy_enabled(default):
                    continue
            except ValueError:
                pass
        if value != default:
            return True
    return False


def _record(adjustments: list[Any], message: str, emit: bool) -> None:
    if message not in adjustments:
        adjustments.append(message)
    if emit:
        warnings.warn(message, UserWarning, stacklevel=3)


def _remove_cpp_outputs(config: dict[str, Any], adjustments: list[Any], emit: bool) -> None:
    output = config.setdefault("output", {})
    value = output.get("types", ())
    if isinstance(value, str):
        items = [item.strip().lower() for item in value.split(",") if item.strip()]
    else:
        try:
            items = [str(item).strip().lower() for item in value if str(item).strip()]
        except TypeError:
            return
    unsupported = [name for name in ("half-gro", "quasi-gro") if name in items]
    if unsupported:
        output["types"] = [name for name in items if name not in set(unsupported)]
        _record(
            adjustments,
            "SQQ-CPP removed unsupported YAML output type(s): " + ", ".join(unsupported),
            emit,
        )


def normalize_engine_capabilities(
    config: dict[str, Any],
    *,
    user_config: dict[str, Any] | None = None,
    emit_warnings: bool = True,
) -> dict[str, Any]:
    """Resolve auto switches and safely disable unsupported C++ analyses."""
    cpp = is_cpp_mode(config.get("mode", DEFAULT_MODE))
    adjustments = list(config.get("adjustments", ()))
    half = config.setdefault("half_cage", {})
    quasi = config.setdefault("quasi_cage", {})
    half_value = _auto_toggle(half.get("enabled", "auto"), not cpp)
    quasi_value = _auto_toggle(quasi.get("enabled", "auto"), not cpp)
    if cpp:
        explicit_half = _explicit(user_config, "half_cage")
        explicit_quasi = _explicit(user_config, "quasi_cage")
        if _explicitly_enabled(explicit_half.get("enabled")):
            _record(adjustments, "half_cage disabled by SQQ-CPP", emit_warnings)
        customized = any(
            key != "enabled" and value != DEFAULT_CONFIG["quasi_cage"].get(key)
            for key, value in explicit_quasi.items()
        )
        if _explicitly_enabled(explicit_quasi.get("enabled")) or customized:
            _record(
                adjustments,
                "quasi_cage disabled by SQQ-CPP; quasi settings ignored",
                emit_warnings,
            )
        _remove_cpp_outputs(config, adjustments, emit_warnings)
        half_value = quasi_value = False
        for key, value in DEFAULT_CONFIG["quasi_cage"].items():
            if key != "enabled":
                quasi[key] = deepcopy(value)
        cluster = config.setdefault("hydrate_cluster", {})
        explicit_cluster = _explicit(user_config, "hydrate_cluster")
        if explicit_cluster and _section_customized(
            explicit_cluster, DEFAULT_CONFIG["hydrate_cluster"]
        ):
            _record(
                adjustments,
                "hydrate_cluster disabled by SQQ-CPP; cluster settings ignored",
                emit_warnings,
            )
        cluster.clear()
        cluster.update(deepcopy(DEFAULT_CONFIG["hydrate_cluster"]))
        cluster["enabled"] = False

        ice = config.setdefault("ice", {})
        explicit_ice = _explicit(user_config, "ice")
        if explicit_ice and _section_customized(
            explicit_ice, DEFAULT_CONFIG["ice"]
        ):
            _record(
                adjustments,
                "ice analysis disabled by SQQ-CPP; ice settings ignored",
                emit_warnings,
            )
        ice.clear()
        ice.update(deepcopy(DEFAULT_CONFIG["ice"]))
        ice["enabled"] = False
    half["enabled"] = half_value
    quasi["enabled"] = quasi_value
    if adjustments:
        config["adjustments"] = list(dict.fromkeys(str(item) for item in adjustments))
    return config


__all__ = [
    "normalize_cpp_order_parameters", "normalize_cpp_output_types",
    "normalize_engine_capabilities", "validate_cpp_cli",
]
