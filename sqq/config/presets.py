"""Public engine selectors and compatibility presets."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from .defaults import (
    CPP_MODES,
    DEFAULT_ORDER_PARAMETERS,
    MODE_PRESETS,
)


def normalize_mode(value: Any) -> str:
    """Normalize and validate a public engine selector."""
    text = str(value).strip().lower()
    if text.isdigit():
        text = text.zfill(2)
    if text not in MODE_PRESETS:
        raise ValueError(f"engine must be one of: {', '.join(MODE_PRESETS)}")
    return text


def mode_label(mode: Any) -> str:
    return str(MODE_PRESETS[normalize_mode(mode)]["label"])


def is_cpp_mode(mode: Any) -> bool:
    return normalize_mode(mode) in CPP_MODES


def engine_display(selector: Any) -> str:
    return "sqq-cpp" if is_cpp_mode(selector) else "sqq-py"


def profile_name(selector: Any) -> str:
    return str(MODE_PRESETS[normalize_mode(selector)]["profile"])


def mode_display(mode: Any) -> str:
    normalized = normalize_mode(mode)
    if normalized in CPP_MODES:
        return "sqq-cpp" if normalized == "cpp" else f"{normalized} (sqq-cpp)"
    return f"{normalized} (sqq-py)"


def mode_worker_count(mode: Any) -> int | None:
    value = MODE_PRESETS[normalize_mode(mode)].get("worker_count")
    return None if value is None else int(value)


def mode_worker_fraction(mode: Any) -> float:
    normalized = normalize_mode(mode)
    preset = MODE_PRESETS[normalized]
    if "worker_fraction" not in preset:
        raise ValueError(f"engine {normalized} uses a fixed worker count")
    return float(preset["worker_fraction"])


def apply_mode_preset(config: dict[str, Any], mode: Any) -> dict[str, Any]:
    """Apply the compatibility preset selected by ``-e``/``--engine``."""
    normalized = normalize_mode(mode)
    preset = MODE_PRESETS[normalized]
    config["mode"] = normalized
    config["graph"]["bond_mode"] = preset["bond_mode"]
    config["ring"]["sizes"] = list(preset["ring_sizes"])
    config["ring"]["report_sizes"] = "auto"
    config["quasi_cage"]["base_sizes"] = "auto"
    config["quasi_cage"]["side_sizes"] = "auto"
    config["hydrate_cluster"]["enabled"] = bool(preset["find_cluster"])
    config["output"]["types"] = deepcopy(preset["output_types"])
    if is_cpp_mode(normalized):
        config["half_cage"]["enabled"] = False
        config["quasi_cage"]["enabled"] = False
        config["ice"]["enabled"] = False
        config["order"]["parameters"] = list(DEFAULT_ORDER_PARAMETERS)
    return config


__all__ = [
    "apply_mode_preset",
    "engine_display",
    "is_cpp_mode",
    "mode_display",
    "mode_label",
    "mode_worker_count",
    "mode_worker_fraction",
    "normalize_mode",
    "profile_name",
]
