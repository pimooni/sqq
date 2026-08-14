"""Configuration file loading and deterministic resolution."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

from .capabilities import normalize_engine_capabilities
from .defaults import CONFIG_SCHEMA_VERSION, DEFAULT_CONFIG, DEFAULT_MODE, MODE_PRESETS
from .migrate import (
    migrate_legacy_order_parameters,
    migrate_yaml_keys,
    strip_legacy_selection_keys,
)
from .overrides import merge_config
from .presets import apply_mode_preset, engine_display, normalize_mode, profile_name
from .validation import validate_user_config_keys


if yaml is not None:
    class UniqueKeySafeLoader(yaml.SafeLoader):
        """Safe YAML loader that rejects duplicate mapping keys."""

        def construct_mapping(self, node, deep: bool = False):
            self.flatten_mapping(node)
            mapping: dict[Any, Any] = {}
            for key_node, value_node in node.value:
                key = self.construct_object(key_node, deep=deep)
                try:
                    duplicate = key in mapping
                except TypeError as exc:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping", node.start_mark,
                        "found an unhashable mapping key", key_node.start_mark,
                    ) from exc
                if duplicate:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping", node.start_mark,
                        f"found duplicate key {key!r}", key_node.start_mark,
                    )
                mapping[key] = self.construct_object(value_node, deep=deep)
            return mapping
else:  # pragma: no cover
    UniqueKeySafeLoader = None


def _read_user_config(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    text = path.read_text(encoding="utf-8-sig")
    if yaml is not None:
        value = yaml.load(text, Loader=UniqueKeySafeLoader) or {}
    else:
        try:
            value = json.loads(text) if text.strip() else {}
        except json.JSONDecodeError as exc:
            raise RuntimeError("Reading YAML config files requires PyYAML.") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Config file must contain a YAML mapping: {path}")
    return value


def _resolve_pair_file(user_config: dict[str, Any], path: Path | None) -> None:
    if path is None:
        return
    graph = user_config.get("graph", {})
    pair_file = graph.get("pair_file") if isinstance(graph, dict) else None
    if pair_file in (None, ""):
        return
    pair_path = Path(str(pair_file)).expanduser()
    if not pair_path.is_absolute():
        pair_path = path.resolve().parent / pair_path
    graph["pair_file"] = str(pair_path.resolve())


def _preset_adjustments(selector: str, config: dict[str, Any]) -> list[dict[str, Any]]:
    if selector not in {"00", "99"}:
        return []
    preset = MODE_PRESETS[selector]
    candidates = (
        ("graph.mode", config["graph"].get("bond_mode"), preset["bond_mode"]),
        ("ring.size", config["ring"].get("sizes"), preset["ring_sizes"]),
        (
            "hydrate_cluster.enabled",
            config["hydrate_cluster"].get("enabled"),
            preset["find_cluster"],
        ),
        ("output.type", config["output"].get("types"), preset["output_types"]),
    )
    result = [
        {
            "parameter": parameter,
            "requested": selector,
            "effective": effective,
            "reason": f"forced by compatibility profile {selector}",
        }
        for parameter, effective, expected in candidates
        if effective == expected
    ]
    worker = config.get("parallel", {}).get("workers", "auto")
    if worker in (None, "", "auto"):
        result.append(
            {
                "parameter": "parallel.worker",
                "requested": selector,
                "effective": "100%",
                "reason": f"forced by compatibility profile {selector}",
            }
        )
    return result


def refresh_resolution_report(config: dict[str, Any]) -> None:
    """Refresh auditable adjustments after overrides and normalization."""
    selector = normalize_mode(config.get("mode", DEFAULT_MODE))
    adjustments = _preset_adjustments(selector, config)
    for message in config.get("adjustments", ()):
        reason = str(message)
        legacy_fast_closure = reason.startswith("Legacy cage.fast_closure")
        adjustments.append(
            {
                "parameter": (
                    "cage.fast_closure"
                    if legacy_fast_closure
                    else "engine.capability"
                ),
                "requested": "legacy setting" if legacy_fast_closure else selector,
                "effective": "ignored" if legacy_fast_closure else None,
                "reason": reason,
            }
        )
    config["resolution_report"] = {
        "requested_selector": selector,
        "engine": engine_display(selector),
        "profile": profile_name(selector),
        "adjustments": adjustments,
    }


def load_config(path: Path | None, mode: Any = None) -> dict[str, Any]:
    """Resolve defaults, preset, migrated YAML, and engine capabilities."""
    return resolve_config(path, mode=mode)


def resolve_config(
    source: str | Path | Mapping[str, Any] | None = None,
    mode: Any = None,
) -> dict[str, Any]:
    """Resolve a configuration path, partial mapping, or defaults."""
    if isinstance(source, Mapping):
        user_config = deepcopy(dict(source))
        source_path = None
    elif source is None:
        user_config = {}
        source_path = None
    else:
        source_path = Path(source).expanduser()
        user_config = _read_user_config(source_path)
    return _resolve_user_config(user_config, source_path, mode)


def _resolve_user_config(
    user_config: dict[str, Any],
    source_path: Path | None,
    mode: Any,
) -> dict[str, Any]:
    """Apply the deterministic resolver to one isolated user mapping."""
    selector_source = (
        mode
        if mode is not None
        else user_config.get("engine", user_config.get("mode", DEFAULT_MODE))
    )
    selector = normalize_mode(selector_source)
    config = apply_mode_preset(deepcopy(DEFAULT_CONFIG), selector)

    migrate_yaml_keys(user_config)
    _resolve_pair_file(user_config, source_path)
    migrated = migrate_legacy_order_parameters(user_config)
    if migrated is not None:
        user_config.setdefault("order", {})["parameters"] = list(migrated)
    strip_legacy_selection_keys(user_config)
    validate_user_config_keys(user_config)

    merge_config(config, user_config)
    config["mode"] = selector
    config["schema_version"] = CONFIG_SCHEMA_VERSION
    render = config.setdefault("render", {})
    scope = str(render.get("atom_scope", "full")).strip().lower()
    if scope not in {"full", "compact"}:
        raise ValueError("render.atom_scope must be full or compact.")
    render["atom_scope"] = scope
    normalize_engine_capabilities(config, user_config=user_config)

    run = config.setdefault("run", {})
    run["engine_selector"] = selector
    run["engine"] = engine_display(selector)
    run["profile"] = profile_name(selector)
    refresh_resolution_report(config)
    return config


__all__ = ["load_config", "refresh_resolution_report", "resolve_config"]
