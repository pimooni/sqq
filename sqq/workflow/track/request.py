"""Configuration and target-request preparation for the Track workflow."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import yaml

from ...config import DEFAULT_MODE, is_cpp_mode, normalize_mode
from ...io.input.lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    inspect_lammps_topology_mapping,
)
from ...models.tracking import TargetSpec, TrackingConfig, TrackingResult


TRACKING_FIELDS = {
    "min_jaccard",
    "min_shared_fraction",
    "min_shared_waters",
    "max_center_distance_nm",
    "gap_frame",
    "max_gap_ps",
    "guest_tiebreak",
    "ambiguity_score_margin",
    "near_threshold_tolerance",
}
TRACK_RUNTIME_FIELDS = {"target", "source"}
RAW_TRACK_FIELD_ALIASES = {"min_shared_water": "min_shared_waters"}
_TRACK_ID = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)


def prepare_target_capabilities(
    config: dict[str, Any],
    targets: Sequence[TargetSpec],
) -> None:
    if not any(target.kind == "phase" for target in targets):
        return
    if is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        raise ValueError(
            "Phase targets require SQQ-Py hydrate-cluster classification; "
            "use -e py or -e 00."
        )
    config.setdefault("hydrate_cluster", {})["enabled"] = True


def update_track_config(
    config: dict[str, Any],
    args: Namespace,
    targets: Sequence[TargetSpec],
) -> None:
    section = config.setdefault("track", {})
    if not isinstance(section, dict):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    section["target"] = ",".join(target.raw for target in targets)
    section["source"] = (
        str(Path(args.source).resolve())
        if getattr(args, "source", None)
        else None
    )
    for name in TRACKING_FIELDS:
        value = getattr(args, name, None)
        if value is not None:
            section[name] = value


def resolve_track_request(
    config: Mapping[str, Any],
    args: Namespace,
) -> tuple[str | Iterable[str], str | Path | None]:
    """Resolve CLI-over-YAML target/source values before mutating config."""
    section = config.get("track", {})
    if section is None:
        section = {}
    if not isinstance(section, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    cli_target = getattr(args, "target", None)
    target = (
        section.get("target", "all")
        if cli_target is None or cli_target == ""
        else cli_target
    )
    cli_source = getattr(args, "source", None)
    configured_source = section.get("source")
    source = cli_source if cli_source not in (None, "") else configured_source
    if source in (None, ""):
        source = None
    return target, source  # type: ignore[return-value]


def source_engine_selector(state_path: Path) -> str | None:
    """Read the Analyze engine recorded beside an imported Track state."""
    candidates = (
        state_path.parent.parent / "sqq_config_resolved.yaml",
        state_path.parent / "sqq_config_resolved.yaml",
    )
    config_path = next((path for path in candidates if path.is_file()), None)
    if config_path is None:
        return None
    try:
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ValueError(
            f"Cannot read the Analyze engine from {config_path}: {exc}"
        ) from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"Invalid resolved SQQ configuration: {config_path}")
    raw = payload.get("engine", payload.get("mode"))
    if raw in (None, ""):
        run = payload.get("run", {})
        if isinstance(run, Mapping):
            raw = run.get("engine_selector")
    if raw in (None, ""):
        return None
    return normalize_mode(raw)


def tracking_config(config: Mapping[str, Any]) -> TrackingConfig:
    values = config.get("track", {})
    if not isinstance(values, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    unknown = sorted(set(values).difference(TRACKING_FIELDS | TRACK_RUNTIME_FIELDS))
    if unknown:
        raise ValueError(
            "Unsupported track configuration field(s): "
            + ", ".join(unknown)
            + ". Remove fields that have no implemented effect."
        )
    return TrackingConfig.from_mapping(
        {name: values[name] for name in TRACKING_FIELDS if name in values}
    )


def explicit_tracking_fields(args: Namespace) -> set[str]:
    """Return fields explicitly supplied by CLI or the user's YAML file."""
    fields = {
        name
        for name in TRACKING_FIELDS
        if getattr(args, name, None) is not None
    }
    config_path = getattr(args, "config", None)
    if not config_path:
        return fields
    try:
        payload = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ValueError(f"Cannot inspect Track settings in {config_path}: {exc}") from exc
    if payload is None:
        return fields
    if not isinstance(payload, Mapping):
        raise ValueError("The configuration root must be a mapping.")
    section = payload.get("track", {})
    if section is None:
        return fields
    if not isinstance(section, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    for raw_name in section:
        normalized = RAW_TRACK_FIELD_ALIASES.get(str(raw_name), str(raw_name))
        if normalized in TRACKING_FIELDS:
            fields.add(normalized)
    return fields


def prepare_lammps_metadata(
    config: dict[str, Any],
    trajectory: Path,
    topology: Path | None,
    args: Namespace,
) -> None:
    if trajectory.suffix.lower() not in LAMMPS_TRAJECTORY_SUFFIXES:
        return
    lammps = config["input"]["lammps"]
    explicit = bool(lammps.get("type_map"))
    resolved: Mapping[str, Any] = {}
    rebuilt = False
    if topology is not None:
        resolved, rebuilt = inspect_lammps_topology_mapping(topology, lammps)
    mapping = ", ".join(
        (
            f"{type_id}=ignore"
            if entry.ignore
            else f"{type_id}={entry.resname}/{entry.atomname}"
        )
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    )
    lammps["resolved_type_map"] = {
        type_id: (
            {"ignore": True}
            if entry.ignore
            else {"resname": entry.resname, "atomname": entry.atomname}
        )
        for type_id, entry in sorted(resolved.items(), key=lambda item: int(item[0]))
    }
    if explicit:
        source = (
            str(Path(args.config).resolve())
            if getattr(args, "config", None)
            else "<configuration>"
        )
    else:
        source = "auto (DATA topology)"
        if rebuilt:
            source += "; molecule IDs rebuilt from Bonds"
    lammps["type_map_source"] = f"{source}: {mapping}" if mapping else source


def validate_requested_tracks(
    result: TrackingResult,
    targets: Sequence[TargetSpec],
) -> None:
    known = {track.track_id for track in result.tracks}
    missing = [
        target.value
        for target in targets
        if target.kind == "track" and target.value not in known
    ]
    if missing:
        preview = ", ".join(sorted(known, key=track_sort_key)[:12]) or "none"
        raise ValueError(
            "Unknown persistent cage ID(s): "
            + ", ".join(missing)
            + f". Available IDs begin with: {preview}."
        )


def track_sort_key(value: str) -> int:
    match = _TRACK_ID.fullmatch(value)
    return int(match.group(1)) if match is not None else 2**63 - 1


__all__ = [
    "TRACKING_FIELDS",
    "explicit_tracking_fields",
    "prepare_lammps_metadata",
    "prepare_target_capabilities",
    "resolve_track_request",
    "source_engine_selector",
    "tracking_config",
    "update_track_config",
    "validate_requested_tracks",
]
