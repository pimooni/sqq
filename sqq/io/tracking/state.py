"""Serialization and discovery for persistent cage-tracking state."""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
import csv
import json
import os
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence
from uuid import uuid4

import numpy as np

from ...core.tracking import (
    cage_lineage_rows,
    event_rows,
    guest_event_rows,
    guest_residence_lifetime_rows,
    guest_residence_rows,
    lifetime_distribution_rows,
    lifetime_rows,
    lifetime_survival_rows,
    observation_rows,
    occupancy_state_lifetime_rows,
    occupancy_transition_rows,
    population_rows,
    select_targets,
    tracking_quality_rows,
)
from ...models.tracking import (
    CageObservation,
    CageTrack,
    FrameStamp,
    Row,
    TargetSelection,
    TargetSpec,
    TrackEvent,
    TrackingConfig,
    TrackingResult,
)
from .network import (
    cage_transition_edge_rows,
    cage_transition_node_rows,
    write_cage_transition_plot,
)


TRACK_DIRECTORY_NAME = "track"
TRACK_STATE_NAME = "track_state.json"
TRACK_INFO_NAME = "track_info.md"

_STATE_FORMAT = "SQQ track state"
_STATE_VERSION = 4
_SUPPORTED_STATE_VERSIONS = {1, 2, 3, 4}
_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
_TRACK_ID = re.compile(r"^t[1-9][0-9]*$", re.IGNORECASE)
_CAGE_TERM = re.compile(r"([456])\^([0-9]+)")
_NUMERIC_SUFFIX = re.compile(r"(\d+)$")
_EVENT_KINDS = {
    "birth",
    "death",
    "type_change",
    "phase_change",
    "type_change_unresolved",
    "phase_change_unresolved",
    "gap",
    "split",
    "merge",
    "split_candidate",
    "merge_candidate",
    "split_confirmed",
    "merge_confirmed",
    "guest_enter",
    "guest_exit",
    "guest_exchange",
    "occupancy_change",
    "empty_to_occupied",
    "occupied_to_empty",
    "single_to_multiple",
    "multiple_to_single",
    "unresolved_across_gap",
}

_OBSERVATION_FIELDS = (
    "target", "track_id", "frame_index", "frame", "time_ps", "cage_id",
    "cage_type", "phase", "water_atomids", "water_count", "guest_ids",
    "guest_count", "center_x_nm", "center_y_nm", "center_z_nm",
    "match_jaccard", "match_shared_fraction", "match_center_distance_nm",
    "match_topology_similarity", "match_shared_waters", "match_candidate_count",
    "match_score", "match_runner_up_score", "match_score_margin",
    "match_jaccard_margin", "match_shared_fraction_margin",
    "match_shared_waters_margin", "match_center_distance_margin_nm",
    "match_near_threshold", "match_status", "match_diagnostic_source",
    "gap_frames", "gap_time_ps",
)
_TRACK_FIELDS = (
    "target", "track_id", "first_frame_index", "last_frame_index",
    "first_time_ps", "last_time_ps", "end_time_ps", "observed_frames",
    "span_frames", "gap_frames", "gap_time_ps", "observed_span_ps",
    "lifetime_ps", "lifetime_lower_ps", "lifetime_upper_ps",
    "occupancy_time_ps", "duration_status", "left_censored", "right_censored",
    "matching_observations", "match_score_min", "match_score_mean",
    "match_score_margin_min", "match_score_margin_mean", "match_jaccard_min",
    "match_jaccard_mean", "ambiguous_observations", "gap_bridge_observations",
    "near_threshold_observations", "initial_cage_type", "final_cage_type",
    "cage_types", "phases",
)
_LIFETIME_DISTRIBUTION_FIELDS = (
    "target", "lifetime_ps", "span_frames", "track_count", "fraction",
    "cumulative_fraction", "survival_fraction", "uncensored_count",
    "left_censored_count", "right_censored_count",
)
_EVENT_FIELDS = (
    "target", "event_id", "event", "frame_index", "frame", "time_ps",
    "source_track_ids", "destination_track_ids", "source_cage_types",
    "destination_cage_types", "source_phases", "destination_phases",
    "gap_frames", "gap_time_ps", "censored", "guest_ids_entered",
    "guest_ids_exited", "source_occupancy", "destination_occupancy",
    "guest_composition_before", "guest_composition_after", "evidence_status",
    "water_conservation", "center_distance_nm", "persistence_frames",
)
_POPULATION_FIELDS = (
    "target", "frame_index", "frame", "time_ps", "group", "label", "cage_count",
)
_GUEST_RESIDENCE_FIELDS = (
    "target", "track_id", "guest_id", "episode", "start_frame_index",
    "end_frame_index", "start_time_ps", "end_time_ps", "observed_frames",
    "residence_time_ps", "left_censored", "right_censored",
)
_TRACKING_QUALITY_FIELDS = (
    "target", "track_id", "frame_index", "frame", "time_ps", "cage_id",
    "cage_type", "candidate_count", "match_score", "runner_up_score",
    "score_margin", "match_jaccard", "match_shared_fraction",
    "match_shared_waters", "match_center_distance_nm",
    "match_topology_similarity", "jaccard_margin", "shared_fraction_margin",
    "shared_waters_margin", "center_distance_margin_nm", "near_threshold",
    "gap_frames", "gap_time_ps", "match_status", "diagnostic_source",
)
_GUEST_EVENT_FIELDS = (
    "target", "event_id", "event", "frame_index", "frame", "time_ps",
    "track_id", "cage_type", "phase", "guest_ids_entered", "guest_ids_exited",
    "guest_composition_before", "guest_composition_after", "source_occupancy",
    "destination_occupancy", "gap_frames", "gap_time_ps", "censored",
)
_GUEST_RESIDENCE_LIFETIME_FIELDS = (
    "target", "track_id", "guest_id", "residence_index", "start_frame_index",
    "end_frame_index", "start_time_ps", "end_time_ps", "observed_frames",
    "observed_span_ps", "residence_lifetime_lower_ps",
    "residence_lifetime_upper_ps", "residence_time_ps", "cage_types", "phases",
    "duration_status", "left_censored", "right_censored",
    "left_gap_unresolved", "right_gap_unresolved", "gap_censored",
)
_OCCUPANCY_STATE_LIFETIME_FIELDS = (
    "target", "track_id", "occupancy_state", "guest_composition",
    "residence_index", "start_frame_index", "end_frame_index", "start_time_ps",
    "end_time_ps", "observed_frames", "observed_span_ps",
    "residence_lifetime_lower_ps", "residence_lifetime_upper_ps",
    "residence_time_ps", "cage_types", "phases", "duration_status",
    "left_censored", "right_censored",
    "left_gap_unresolved", "right_gap_unresolved", "gap_censored",
)
_OCCUPANCY_TRANSITION_FIELDS = (
    "target", "source_state", "destination_state", "transition_count",
    "independent_track_count", "source_residence_ps", "conditional_probability",
    "rate_per_ns",
)
_LIFETIME_SURVIVAL_FIELDS = (
    "target", "time_basis", "time_ps", "number_at_risk", "deaths",
    "right_censored", "left_censored_excluded", "time_unavailable_excluded",
    "survival_probability", "greenwood_variance", "greenwood_standard_error",
    "confidence_lower", "confidence_upper",
)
_CAGE_LINEAGE_FIELDS = (
    "target", "event_id", "event", "frame_index", "frame", "time_ps",
    "parent_track_ids", "child_track_ids", "source_cage_types",
    "destination_cage_types", "water_conservation", "center_distance_nm",
    "persistence_frames", "confirmation_status", "censored",
)
_CAGE_TRANSITION_NODE_FIELDS = (
    "target", "cage_type", "track_count", "observation_count",
    "residence_time_ps", "occupancy_composition", "phase_composition",
    "incoming_degree", "outgoing_degree", "incoming_transition_count",
    "outgoing_transition_count",
)
_CAGE_TRANSITION_EDGE_FIELDS = (
    "target", "source_cage_type", "destination_cage_type", "transition_count",
    "independent_track_count", "source_state_exposure_ps",
    "conditional_probability", "rate_per_ns",
)

# Relative output path -> CSV columns, in publication order.
_TABLE_FIELDS: dict[str, tuple[str, ...]] = {
    "cage_observation.csv": _OBSERVATION_FIELDS,
    "cage_track.csv": _TRACK_FIELDS,
    "lifetime_distribution.csv": _LIFETIME_DISTRIBUTION_FIELDS,
    "cage_event.csv": _EVENT_FIELDS,
    "cage_population.csv": _POPULATION_FIELDS,
    "guest_residence.csv": _GUEST_RESIDENCE_FIELDS,
    "statistics/tracking_quality.csv": _TRACKING_QUALITY_FIELDS,
    "statistics/guest_event.csv": _GUEST_EVENT_FIELDS,
    "statistics/guest_residence_lifetime.csv": _GUEST_RESIDENCE_LIFETIME_FIELDS,
    "statistics/occupancy_state_lifetime.csv": _OCCUPANCY_STATE_LIFETIME_FIELDS,
    "statistics/occupancy_transition.csv": _OCCUPANCY_TRANSITION_FIELDS,
    "statistics/lifetime_survival.csv": _LIFETIME_SURVIVAL_FIELDS,
    "statistics/cage_lineage.csv": _CAGE_LINEAGE_FIELDS,
    "network/cage_transition_nodes.csv": _CAGE_TRANSITION_NODE_FIELDS,
    "network/cage_transition_edges.csv": _CAGE_TRANSITION_EDGE_FIELDS,
}
_STATISTICS_TABLE_NAMES = tuple(
    name.split("/", 1)[1] for name in _TABLE_FIELDS if name.startswith("statistics/")
)

__all__ = [
    "TRACK_DIRECTORY_NAME",
    "TRACK_INFO_NAME",
    "TRACK_STATE_NAME",
    "TrackTables",
    "append_track_info_section",
    "build_track_tables",
    "deserialize_tracking_result",
    "discover_track_state",
    "add_precursor_membership",
    "read_tracking_result",
    "rewrite_membership_track_ids",
    "serialize_tracking_result",
    "target_directory_name",
    "write_target_selection",
    "write_track_info",
    "write_track_outputs",
    "write_tracking_result",
    "write_tracking_tables",
]


def append_track_info_section(path: str | Path, text: str) -> Path:
    """Crash-safely append one Markdown section for a single owned run.

    SQQ Track holds the output-root lock while this read/replace operation is
    used.  The helper prevents partial files; it is not a substitute for
    caller-level exclusion between unrelated writers.
    """
    target = Path(path)
    existing = (
        target.read_text(encoding="utf-8").rstrip()
        if target.is_file()
        else "# SQQ Track"
    )
    _atomic_write_text(target, existing + "\n\n" + text.rstrip() + "\n")
    return target


def serialize_tracking_result(result: TrackingResult) -> dict[str, object]:
    """Return the current JSON-safe tracking state."""
    if not isinstance(result, TrackingResult):
        raise TypeError("result must be a TrackingResult.")
    return {
        "format": _STATE_FORMAT,
        "version": _STATE_VERSION,
        "tracking_config": result.config.to_dict(),
        "source_provenance": dict(result.source_provenance),
        "frames": [_frame_to_dict(frame) for frame in result.frames],
        "tracks": [_track_to_dict(track) for track in result.tracks],
        "events": [_event_to_dict(event) for event in result.events],
    }


def deserialize_tracking_result(payload: Mapping[str, object]) -> TrackingResult:
    """Read state and migrate archived gap-change event semantics."""
    if not isinstance(payload, Mapping):
        raise TypeError("Tracking state must be a mapping.")
    if payload.get("format") != _STATE_FORMAT:
        raise ValueError("Unsupported tracking-state format.")
    version = payload.get("version")
    if version not in _SUPPORTED_STATE_VERSIONS:
        raise ValueError(f"Unsupported tracking-state version: {version!r}.")
    raw_frames = _mapping_sequence(payload.get("frames"), "frames")
    raw_tracks = _mapping_sequence(payload.get("tracks"), "tracks")
    raw_events = _mapping_sequence(payload.get("events"), "events")
    config = (
        TrackingConfig()
        if version == 1
        else TrackingConfig.from_mapping(
            _required_mapping(payload.get("tracking_config"), "tracking_config")
        )
    )
    frames = tuple(_frame_from_dict(value) for value in raw_frames)
    tracks = tuple(_track_from_dict(value) for value in raw_tracks)
    events = tuple(_event_from_dict(value) for value in raw_events)
    if version < 4:
        events = _migrate_gap_change_events(events, tracks)
        tracks = _migrate_partial_diagnostics(tracks, config, state_version=version)
    provenance = dict(
        _required_mapping(
            payload.get("source_provenance", {}),
            "source_provenance",
        )
    )
    _validate_result_identity(frames, tracks, events, config)
    return TrackingResult(
        frames=frames, tracks=tracks, events=events, config=config,
        source_provenance=provenance,
    )


def write_tracking_result(result: TrackingResult, path: str | Path) -> Path:
    """Atomically write ``track_state.json`` or an explicit JSON path."""
    target = Path(path)
    if target.exists() and target.is_dir():
        target = target / TRACK_STATE_NAME
    elif not target.suffix:
        target = target / TRACK_STATE_NAME
    text = json.dumps(
        serialize_tracking_result(result),
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    )
    _atomic_write_text(target, text + "\n")
    return target


def read_tracking_result(source: str | Path | None = None) -> TrackingResult:
    path = discover_track_state(source)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid tracking-state JSON: {path}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"Tracking-state root must be a mapping: {path}")
    return deserialize_tracking_result(payload)


def discover_track_state(source: str | Path | None = None) -> Path:
    """Find one Track state, without pretending an old Analyze result has one."""
    root = Path.cwd() if source is None else Path(source)
    if root.is_file():
        if root.name != TRACK_STATE_NAME:
            raise ValueError(f"Expected {TRACK_STATE_NAME}, got: {root}")
        return root.resolve()
    candidates = [
        root / TRACK_DIRECTORY_NAME / TRACK_STATE_NAME,
        root / TRACK_STATE_NAME,
    ]
    if root.name == TRACK_DIRECTORY_NAME:
        candidates.insert(0, root / TRACK_STATE_NAME)
    existing = _unique_existing_files(candidates)
    if len(existing) == 1:
        return existing[0]
    if len(existing) > 1:
        raise ValueError(
            "Multiple tracking-state files were found: "
            + ", ".join(str(path) for path in existing)
        )
    grouped = sorted(
        (
            path.resolve()
            for path in root.glob(
                f"result_*/{TRACK_DIRECTORY_NAME}/{TRACK_STATE_NAME}"
            )
            if path.is_file()
        ),
        key=lambda path: str(path).casefold(),
    )
    if len(grouped) == 1:
        return grouped[0]
    if len(grouped) > 1:
        raise ValueError(
            "The source contains multiple systems; select one result_A/result_B directory."
        )
    raise FileNotFoundError(
        f"Cannot find {TRACK_DIRECTORY_NAME}/{TRACK_STATE_NAME} in {root}. "
        "This source was not saved with Track state; rerun Analyze with a "
        "current SQQ release."
    )


def target_directory_name(target: TargetSpec) -> str:
    if target.kind == "all":
        name = "all"
    elif target.kind == "cage_type":
        name = "type_" + _compact_cage_type(target.value)
    elif target.kind == "phase":
        name = "phase_" + target.value
    elif target.kind == "track":
        name = "cage_" + target.value
    else:
        raise ValueError(f"Unsupported target kind: {target.kind!r}.")
    return _safe_component(name)


TrackTables = dict[str, list[Row]]


def build_track_tables(data: TrackingResult | TargetSelection) -> TrackTables:
    """Build every Track table once, keyed by its relative output path.

    Derived tables (lifetime distribution, survival, transition nodes, the
    Markdown summary) reuse the shared per-track and edge rows instead of
    recomputing them for every consumer.
    """
    lifetimes = lifetime_rows(data)
    edges = cage_transition_edge_rows(data)
    return {
        "cage_observation.csv": observation_rows(data),
        "cage_track.csv": lifetimes,
        "lifetime_distribution.csv": lifetime_distribution_rows(
            data, samples=lifetimes
        ),
        "cage_event.csv": event_rows(data),
        "cage_population.csv": population_rows(data),
        "guest_residence.csv": guest_residence_rows(data),
        "statistics/tracking_quality.csv": tracking_quality_rows(data),
        "statistics/guest_event.csv": guest_event_rows(data),
        "statistics/guest_residence_lifetime.csv": guest_residence_lifetime_rows(data),
        "statistics/occupancy_state_lifetime.csv": occupancy_state_lifetime_rows(data),
        "statistics/occupancy_transition.csv": occupancy_transition_rows(data),
        "statistics/lifetime_survival.csv": lifetime_survival_rows(
            data, samples=lifetimes
        ),
        "statistics/cage_lineage.csv": cage_lineage_rows(data),
        "network/cage_transition_nodes.csv": cage_transition_node_rows(
            data, edges=edges
        ),
        "network/cage_transition_edges.csv": edges,
    }


def write_tracking_tables(
    data: TrackingResult | TargetSelection,
    directory: str | Path,
    *,
    tables: TrackTables | None = None,
) -> dict[str, Path]:
    """Write legacy root tables plus additive normalized statistics tables."""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    tables = build_track_tables(data) if tables is None else tables
    written: dict[str, Path] = {}
    for relative_name, fields in _TABLE_FIELDS.items():
        path = root / relative_name
        _atomic_write_csv(path, tables[relative_name], fields)
        written[relative_name] = path
    plot = write_cage_transition_plot(
        data,
        root / "network" / "cage_transition_network.png",
        edges=tables["network/cage_transition_edges.csv"],
        nodes=tables["network/cage_transition_nodes.csv"],
    )
    if plot is not None:
        written["network/cage_transition_network.png"] = plot
    return written


def write_track_info(
    data: TrackingResult | TargetSelection,
    directory: str | Path,
    *,
    tables: TrackTables | None = None,
) -> Path:
    """Atomically write the Markdown summary for one Track result or target."""
    path = Path(directory) / TRACK_INFO_NAME
    tables = build_track_tables(data) if tables is None else tables
    _atomic_write_text(path, _track_info_text(data, tables))
    return path


def write_target_selection(
    selection: TargetSelection,
    track_root: str | Path,
) -> Path:
    """Write normalized tables and metadata for one Track target."""
    root = Path(track_root)
    directory = _safe_child(root, target_directory_name(selection.target))
    directory.mkdir(parents=True, exist_ok=True)
    tables = build_track_tables(selection)
    write_tracking_tables(selection, directory, tables=tables)
    write_track_info(selection, directory, tables=tables)
    return directory


def write_track_outputs(
    result: TrackingResult,
    outdir: str | Path,
    *,
    targets: str | Iterable[str] = "all",
) -> dict[str, Path]:
    """Write Track state, normalized tables, and target metadata."""
    root = Path(outdir) / TRACK_DIRECTORY_NAME
    root.mkdir(parents=True, exist_ok=True)
    tables = build_track_tables(result)
    written = write_tracking_tables(result, root, tables=tables)
    written[TRACK_STATE_NAME] = write_tracking_result(result, root / TRACK_STATE_NAME)
    written[TRACK_INFO_NAME] = write_track_info(result, root, tables=tables)
    for selection in select_targets(result, targets):
        directory = write_target_selection(selection, root)
        written[target_directory_name(selection.target)] = directory
    return written


def rewrite_membership_track_ids(
    path: str | Path,
    data: TrackingResult | TargetSelection,
) -> Path:
    """Atomically map cage C/M records in a copied TSV to persistent IDs."""
    target = Path(path)
    observations = data.observations
    mapping: dict[tuple[int, str, str], str] = {}
    for item in observations:
        key = (
            int(item.frame_index),
            _compact_object_id(item.local_cage_id),
            str(item.cage_type),
        )
        previous = mapping.setdefault(key, item.track_id)
        if previous != item.track_id:
            raise ValueError(
                f"Frame {item.frame_index} cage {key[2]}:{key[1]} maps to both "
                f"{previous} and {item.track_id}."
            )
    expected_persistent = {
        (int(item.frame_index), str(item.track_id), str(item.cage_type))
        for item in observations
    }
    temporary = _temporary_path(target)
    try:
        with target.open("r", encoding="utf-8", newline="") as source_handle:
            reader = csv.DictReader(source_handle, delimiter="\t")
            required = {"record", "render_frame", "family", "cage_id", "cage_type"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Invalid SQQ membership TSV: {target}")
            fields = list(reader.fieldnames)
            frame_by_render: dict[int, int] = {}
            frame_stamps = tuple(data.frames)
            seen_persistent: set[tuple[int, str, str]] = set()
            with temporary.open("w", encoding="ascii", newline="") as output:
                writer = csv.DictWriter(
                    output,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                for row in reader:
                    record = row.get("record")
                    render_index = int(row["render_frame"])
                    if record == "F":
                        if render_index < 0 or render_index >= len(frame_stamps):
                            raise ValueError(
                                f"Membership render frame {render_index} has no "
                                "tracking frame."
                            )
                        if render_index in frame_by_render:
                            raise ValueError(
                                f"Membership render frame {render_index} is repeated."
                            )
                        stamp = frame_stamps[render_index]
                        frame_by_render[render_index] = stamp.frame_index
                        raw_time = row.get("time_ps", "-")
                        if stamp.time_ps is not None and raw_time not in {None, "", "-"}:
                            if abs(float(raw_time) - stamp.time_ps) > max(
                                1.0e-6, abs(stamp.time_ps) * 1.0e-9
                            ):
                                raise ValueError(
                                    f"Membership time for render frame {render_index} "
                                    "does not match tracking state."
                                )
                    elif (
                        (record == "C" and row.get("family") == "cage")
                        or (
                            record == "M"
                            and row.get("family") in {"cage", "guest"}
                        )
                    ):
                        if render_index not in frame_by_render:
                            raise ValueError(
                                f"Cage membership precedes frame record {render_index}."
                            )
                        frame_index = frame_by_render[render_index]
                        cage_id = str(row.get("cage_id", ""))
                        cage_type = str(row.get("cage_type", ""))
                        counts_as_cage = row.get("family") == "cage"
                        key = (frame_index, cage_id, cage_type)
                        track_id = mapping.get(key)
                        if track_id is None and _TRACK_ID.fullmatch(cage_id):
                            if (frame_index, cage_id, cage_type) not in expected_persistent:
                                if isinstance(data, TargetSelection):
                                    continue
                                raise ValueError(
                                    f"Persistent membership {cage_id} does not match "
                                    f"tracking state in frame {frame_index}."
                                )
                            if counts_as_cage:
                                seen_persistent.add((frame_index, cage_id, cage_type))
                        elif track_id is None:
                            if isinstance(data, TargetSelection):
                                continue
                            raise ValueError(
                                f"No persistent ID for cage {cage_type}:{cage_id} "
                                f"in frame {frame_index}."
                            )
                        else:
                            row["cage_id"] = track_id
                            if counts_as_cage:
                                seen_persistent.add((frame_index, track_id, cage_type))
                    writer.writerow(row)
                output.flush()
                os.fsync(output.fileno())
        if len(frame_by_render) != len(frame_stamps):
            raise ValueError(
                "Membership TSV and tracking state have different frame records."
            )
        missing_cages = expected_persistent.difference(seen_persistent)
        if missing_cages:
            preview = ", ".join(
                f"frame {frame}:{track_id}:{cage_type}"
                for frame, track_id, cage_type in sorted(missing_cages)[:5]
            )
            raise ValueError(
                "Tracking state contains cages absent from membership TSV: "
                + preview
            )
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def add_precursor_membership(
    path: str | Path,
    track_id: str,
    frames: Mapping[int, tuple[Sequence[int], Sequence[float]]],
) -> Path:
    """Add pre-birth target-water memberships to one Track render TSV."""
    target = Path(path)
    identifier = str(track_id).strip().lower()
    if not _TRACK_ID.fullmatch(identifier):
        raise ValueError(f"Invalid persistent cage ID for precursor render: {track_id!r}.")
    normalized: dict[int, tuple[tuple[int, ...], tuple[float, float, float]]] = {}
    for raw_frame, value in frames.items():
        frame = int(raw_frame)
        if frame < 0:
            raise ValueError("Precursor render frame indexes must be nonnegative.")
        atom_indexes = tuple(sorted({int(item) for item in value[0]}))
        if not atom_indexes or atom_indexes[0] < 0:
            raise ValueError(
                f"Precursor render frame {frame} must contain nonnegative atom indexes."
            )
        center = tuple(float(item) for item in value[1])
        if len(center) != 3 or any(not np.isfinite(item) for item in center):
            raise ValueError(f"Precursor render frame {frame} has an invalid center.")
        normalized[frame] = (atom_indexes, (center[0], center[1], center[2]))
    if not normalized:
        return target

    with target.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Invalid SQQ membership TSV: {target}")
        fields = list(reader.fieldnames)
        required = {
            "record",
            "render_frame",
            "family",
            "cage_id",
            "cage_type",
            "atom_indices",
            "center_x_angstrom",
            "center_y_angstrom",
            "center_z_angstrom",
        }
        if not required.issubset(fields):
            raise ValueError(f"Invalid SQQ membership TSV: {target}")
        known_frames: set[int] = set()
        for row in reader:
            frame = int(row["render_frame"])
            if row.get("record") == "F":
                known_frames.add(frame)
            elif (
                row.get("record") in {"C", "M"}
                and row.get("family") == "cage"
                and str(row.get("cage_id", "")).lower() == identifier
                and frame in normalized
            ):
                raise ValueError(
                    f"Precursor render frame {frame} already contains cage {identifier}."
                )
    missing = sorted(set(normalized).difference(known_frames))
    if missing:
        raise ValueError(
            "Precursor render references missing membership frame(s): "
            + ", ".join(map(str, missing[:10]))
        )

    temporary = _temporary_path(target)
    inserted: set[int] = set()
    try:
        with target.open("r", encoding="utf-8", newline="") as source_handle:
            reader = csv.DictReader(source_handle, delimiter="\t")
            with temporary.open("w", encoding="ascii", newline="") as output:
                writer = csv.DictWriter(
                    output,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                for row in reader:
                    writer.writerow(row)
                    if row.get("record") != "F":
                        continue
                    frame = int(row["render_frame"])
                    if frame not in normalized:
                        continue
                    atom_indexes, center = normalized[frame]
                    base = {field: "-" for field in fields}
                    base.update(
                        {
                            "render_frame": str(frame),
                            "family": "cage",
                            "cage_id": identifier,
                            "cage_type": "precursor",
                        }
                    )
                    center_row = dict(base)
                    center_row.update(
                        {
                            "record": "C",
                            "atom_indices": "-",
                            "center_x_angstrom": format(center[0], ".17g"),
                            "center_y_angstrom": format(center[1], ".17g"),
                            "center_z_angstrom": format(center[2], ".17g"),
                        }
                    )
                    writer.writerow(center_row)
                    member_row = dict(base)
                    member_row.update(
                        {
                            "record": "M",
                            "atom_indices": ",".join(map(str, atom_indexes)),
                        }
                    )
                    writer.writerow(member_row)
                    inserted.add(frame)
                output.flush()
                os.fsync(output.fileno())
        if inserted != set(normalized):
            raise ValueError("Not all precursor render frames were written.")
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def _frame_to_dict(frame: FrameStamp) -> dict[str, object]:
    return {
        "frame_index": frame.frame_index,
        "frame_name": frame.frame_name,
        "time_ps": frame.time_ps,
        "source": frame.source,
    }


def _frame_from_dict(value: Mapping[str, object]) -> FrameStamp:
    return FrameStamp(
        frame_index=int(_required(value, "frame_index")),
        frame_name=str(_required(value, "frame_name")),
        time_ps=_optional_float(value.get("time_ps")),
        source=str(value.get("source", "")),
    )


def _observation_to_dict(item: CageObservation) -> dict[str, object]:
    return {
        "track_id": item.track_id,
        "frame_index": item.frame_index,
        "frame_name": item.frame_name,
        "time_ps": item.time_ps,
        "local_cage_id": item.local_cage_id,
        "cage_type": item.cage_type,
        "phase": item.phase,
        "phase_labels": list(item.phase_labels),
        "water_atomids": list(item.water_atomids),
        "center": list(item.center),
        "topology": list(item.topology),
        "guest_ids": list(item.guest_ids),
        "match_jaccard": item.match_jaccard,
        "match_shared_fraction": item.match_shared_fraction,
        "match_center_distance_nm": item.match_center_distance_nm,
        "match_topology_similarity": item.match_topology_similarity,
        "match_shared_waters": getattr(item, "match_shared_waters", None),
        "match_candidate_count": getattr(item, "match_candidate_count", 0),
        "match_score": getattr(item, "match_score", None),
        "match_runner_up_score": getattr(item, "match_runner_up_score", None),
        "match_score_margin": getattr(item, "match_score_margin", None),
        "match_jaccard_margin": getattr(item, "match_jaccard_margin", None),
        "match_shared_fraction_margin": getattr(
            item, "match_shared_fraction_margin", None
        ),
        "match_shared_waters_margin": getattr(
            item, "match_shared_waters_margin", None
        ),
        "match_center_distance_margin_nm": getattr(
            item, "match_center_distance_margin_nm", None
        ),
        "match_near_threshold": getattr(item, "match_near_threshold", False),
        "match_status": getattr(item, "match_status", "unavailable"),
        "match_diagnostic_source": getattr(item, "match_diagnostic_source", ""),
        "gap_frames": item.gap_frames,
        "gap_time_ps": getattr(item, "gap_time_ps", None),
    }


def _observation_from_dict(value: Mapping[str, object]) -> CageObservation:
    center = tuple(float(item) for item in _sequence(value.get("center"), "center"))
    if len(center) != 3:
        raise ValueError("Cage observation center must contain three coordinates.")
    phase = str(value.get("phase", "unassigned"))
    raw_phases = value.get("phase_labels", [phase])
    gap = int(value.get("gap_frames", 0))
    if gap < 0:
        raise ValueError("Cage observation gap_frames must be nonnegative.")
    return CageObservation(
        track_id=str(_required(value, "track_id")),
        frame_index=int(_required(value, "frame_index")),
        frame_name=str(_required(value, "frame_name")),
        time_ps=_optional_float(value.get("time_ps")),
        local_cage_id=str(_required(value, "local_cage_id")),
        cage_type=str(_required(value, "cage_type")),
        phase=phase,
        phase_labels=tuple(str(item) for item in _sequence(raw_phases, "phase_labels")),
        water_atomids=tuple(
            int(item) for item in _sequence(value.get("water_atomids"), "water_atomids")
        ),
        center=(center[0], center[1], center[2]),
        topology=tuple(
            int(item) for item in _sequence(value.get("topology"), "topology")
        ),
        guest_ids=tuple(
            str(item) for item in _sequence(value.get("guest_ids"), "guest_ids")
        ),
        match_jaccard=_optional_float(value.get("match_jaccard")),
        match_shared_fraction=_optional_float(value.get("match_shared_fraction")),
        match_shared_waters=_optional_int(value.get("match_shared_waters")),
        match_center_distance_nm=_optional_float(
            value.get("match_center_distance_nm")
        ),
        match_topology_similarity=_optional_float(
            value.get("match_topology_similarity")
        ),
        match_candidate_count=int(value.get("match_candidate_count", 0)),
        match_score=_optional_float(value.get("match_score")),
        match_runner_up_score=_optional_float(value.get("match_runner_up_score")),
        match_score_margin=_optional_float(value.get("match_score_margin")),
        match_jaccard_margin=_optional_float(value.get("match_jaccard_margin")),
        match_shared_fraction_margin=_optional_float(
            value.get("match_shared_fraction_margin")
        ),
        match_shared_waters_margin=_optional_int(
            value.get("match_shared_waters_margin")
        ),
        match_center_distance_margin_nm=_optional_float(
            value.get("match_center_distance_margin_nm")
        ),
        match_near_threshold=_strict_bool(
            value.get("match_near_threshold", False), "match_near_threshold"
        ),
        match_status=str(value.get("match_status", "unavailable")),  # type: ignore[arg-type]
        match_diagnostic_source=str(value.get("match_diagnostic_source") or ""),
        gap_frames=gap,
        gap_time_ps=_optional_float(value.get("gap_time_ps")),
    )


def _track_to_dict(track: CageTrack) -> dict[str, object]:
    return {
        "track_id": track.track_id,
        "left_censored": track.left_censored,
        "right_censored": track.right_censored,
        "observations": [_observation_to_dict(item) for item in track.observations],
    }


def _track_from_dict(value: Mapping[str, object]) -> CageTrack:
    observations = tuple(
        _observation_from_dict(item)
        for item in _mapping_sequence(value.get("observations"), "observations")
    )
    if not observations:
        raise ValueError("A cage track must contain at least one observation.")
    track_id = str(_required(value, "track_id"))
    if any(item.track_id != track_id for item in observations):
        raise ValueError(f"Track {track_id} contains another track ID.")
    return CageTrack(
        track_id=track_id,
        observations=observations,
        left_censored=_strict_bool(value.get("left_censored"), "left_censored"),
        right_censored=_strict_bool(value.get("right_censored"), "right_censored"),
    )


def _event_to_dict(event: TrackEvent) -> dict[str, object]:
    return {
        "event_id": event.event_id,
        "kind": event.kind,
        "frame_index": event.frame_index,
        "frame_name": event.frame_name,
        "time_ps": event.time_ps,
        "source_track_ids": list(event.source_track_ids),
        "destination_track_ids": list(event.destination_track_ids),
        "source_cage_types": list(event.source_cage_types),
        "destination_cage_types": list(event.destination_cage_types),
        "source_phases": list(event.source_phases),
        "destination_phases": list(event.destination_phases),
        "gap_frames": event.gap_frames,
        "gap_time_ps": getattr(event, "gap_time_ps", None),
        "censored": event.censored,
        "guest_ids_entered": list(getattr(event, "guest_ids_entered", ())),
        "guest_ids_exited": list(getattr(event, "guest_ids_exited", ())),
        "source_occupancy": getattr(event, "source_occupancy", None),
        "destination_occupancy": getattr(event, "destination_occupancy", None),
        "guest_composition_before": list(
            getattr(event, "guest_composition_before", ())
        ),
        "guest_composition_after": list(
            getattr(event, "guest_composition_after", ())
        ),
        "evidence_status": getattr(event, "evidence_status", None),
        "water_conservation": getattr(event, "water_conservation", None),
        "center_distance_nm": getattr(event, "center_distance_nm", None),
        "persistence_frames": getattr(event, "persistence_frames", None),
    }


def _event_from_dict(value: Mapping[str, object]) -> TrackEvent:
    kind = str(_required(value, "kind"))
    if kind not in _EVENT_KINDS:
        raise ValueError(f"Unsupported tracking event: {kind!r}.")
    gap = int(value.get("gap_frames", 0))
    if gap < 0:
        raise ValueError("Track event gap_frames must be nonnegative.")
    return TrackEvent(
        event_id=str(_required(value, "event_id")),
        kind=kind,  # type: ignore[arg-type]
        frame_index=int(_required(value, "frame_index")),
        frame_name=str(_required(value, "frame_name")),
        time_ps=_optional_float(value.get("time_ps")),
        source_track_ids=tuple(
            str(item)
            for item in _sequence(value.get("source_track_ids", ()), "source_track_ids")
        ),
        destination_track_ids=tuple(
            str(item)
            for item in _sequence(
                value.get("destination_track_ids", ()), "destination_track_ids"
            )
        ),
        source_cage_types=tuple(
            str(item)
            for item in _sequence(
                value.get("source_cage_types", ()), "source_cage_types"
            )
        ),
        destination_cage_types=tuple(
            str(item)
            for item in _sequence(
                value.get("destination_cage_types", ()),
                "destination_cage_types",
            )
        ),
        source_phases=tuple(
            str(item)
            for item in _sequence(value.get("source_phases", ()), "source_phases")
        ),
        destination_phases=tuple(
            str(item)
            for item in _sequence(
                value.get("destination_phases", ()), "destination_phases"
            )
        ),
        gap_frames=gap,
        gap_time_ps=_optional_float(value.get("gap_time_ps")),
        censored=_strict_bool(value.get("censored"), "censored"),
        guest_ids_entered=tuple(
            str(item)
            for item in _sequence(
                value.get("guest_ids_entered", ()), "guest_ids_entered"
            )
        ),
        guest_ids_exited=tuple(
            str(item)
            for item in _sequence(
                value.get("guest_ids_exited", ()), "guest_ids_exited"
            )
        ),
        source_occupancy=str(value.get("source_occupancy") or ""),
        destination_occupancy=str(value.get("destination_occupancy") or ""),
        guest_composition_before=tuple(
            str(item)
            for item in _sequence(
                value.get("guest_composition_before", ()),
                "guest_composition_before",
            )
        ),
        guest_composition_after=tuple(
            str(item)
            for item in _sequence(
                value.get("guest_composition_after", ()),
                "guest_composition_after",
            )
        ),
        evidence_status=str(value.get("evidence_status") or ""),
        water_conservation=_optional_float(value.get("water_conservation")),
        center_distance_nm=_optional_float(value.get("center_distance_nm")),
        persistence_frames=int(value.get("persistence_frames") or 0),
    )


def _migrate_gap_change_events(
    events: Sequence[TrackEvent],
    tracks: Sequence[CageTrack],
) -> tuple[TrackEvent, ...]:
    """Make pre-v4 type/phase changes across missing cages explicitly unresolved."""
    gap_by_track_frame = {
        (track.track_id, observation.frame_index): observation
        for track in tracks
        for observation in track.observations
        if observation.gap_frames > 0
    }
    migrated: list[TrackEvent] = []
    for event in events:
        if event.kind not in {"type_change", "phase_change"}:
            migrated.append(event)
            continue
        track_ids = event.track_ids
        observation = (
            gap_by_track_frame.get((track_ids[0], event.frame_index))
            if len(track_ids) == 1
            else None
        )
        if observation is None:
            migrated.append(event)
            continue
        migrated.append(
            replace(
                event,
                kind=(
                    "type_change_unresolved"
                    if event.kind == "type_change"
                    else "phase_change_unresolved"
                ),
                gap_frames=observation.gap_frames,
                gap_time_ps=observation.gap_time_ps,
                censored=True,
                evidence_status="gap_unresolved",
            )
        )
    return tuple(migrated)


def _migrate_partial_diagnostics(
    tracks: Sequence[CageTrack],
    config: TrackingConfig,
    *,
    state_version: int,
) -> tuple[CageTrack, ...]:
    """Recover threshold diagnostics when archived evidence permits.

    Schema 2/3 retain the required evidence and thresholds.  Schema 1 does
    not, so its raw evidence remains available without derived diagnostics.
    """
    if state_version <= 1:
        return tuple(tracks)
    has_evidence = any(
        item.match_jaccard is not None
        for track in tracks
        for item in track.observations
    )
    if not has_evidence:
        return tuple(tracks)
    tolerance = config.near_threshold_tolerance
    migrated: list[CageTrack] = []
    for track in tracks:
        observations: list[CageObservation] = []
        for item in track.observations:
            if item.match_status != "unavailable":
                observations.append(item)
                continue
            if item.match_jaccard is None or item.match_shared_fraction is None:
                observations.append(replace(item, match_status="new"))
                continue
            jaccard_margin = item.match_jaccard - config.min_jaccard
            fraction_margin = item.match_shared_fraction - config.min_shared_fraction
            waters_margin = (
                None
                if item.match_shared_waters is None
                else item.match_shared_waters - config.min_shared_waters
            )
            distance_margin = (
                None
                if config.max_center_distance_nm is None
                or item.match_center_distance_nm is None
                else config.max_center_distance_nm - item.match_center_distance_nm
            )
            near_threshold = (
                jaccard_margin <= tolerance
                or fraction_margin <= tolerance
                or (
                    waters_margin is not None
                    and waters_margin / max(1, config.min_shared_waters) <= tolerance
                )
                or (
                    distance_margin is not None
                    and config.max_center_distance_nm not in (None, 0.0)
                    and distance_margin / float(config.max_center_distance_nm)
                    <= tolerance
                )
            )
            if item.gap_frames:
                status = "gap_bridge"
            elif near_threshold:
                status = "ambiguous"
            else:
                status = "secure"
            observations.append(
                replace(
                    item,
                    match_jaccard_margin=jaccard_margin,
                    match_shared_fraction_margin=fraction_margin,
                    match_shared_waters_margin=waters_margin,
                    match_center_distance_margin_nm=distance_margin,
                    match_near_threshold=near_threshold,
                    match_status=status,
                    match_diagnostic_source="migrated_partial",
                )
            )
        migrated.append(replace(track, observations=tuple(observations)))
    return tuple(migrated)


def _validate_result_identity(
    frames: Sequence[FrameStamp],
    tracks: Sequence[CageTrack],
    events: Sequence[TrackEvent],
    config: TrackingConfig,
) -> None:
    frame_indexes = [frame.frame_index for frame in frames]
    if frame_indexes != sorted(set(frame_indexes)):
        raise ValueError("Tracking-state frame indexes must be unique and increasing.")
    if any(
        right != left + 1
        for left, right in zip(frame_indexes, frame_indexes[1:])
    ):
        raise ValueError(
            "Tracking-state frame indexes must be consecutive; empty analyzed "
            "frames must be represented explicitly."
        )
    frame_set = set(frame_indexes)
    frame_positions = {value: position for position, value in enumerate(frame_indexes)}
    frame_by_index = {frame.frame_index: frame for frame in frames}
    numeric_times = [frame.time_ps for frame in frames if frame.time_ps is not None]
    if any(right < left for left, right in zip(numeric_times, numeric_times[1:])):
        raise ValueError("Tracking-state frame times must be nondecreasing.")
    track_ids = [track.track_id for track in tracks]
    if len(track_ids) != len(set(track_ids)):
        raise ValueError("Tracking-state track IDs must be unique.")
    known_tracks = set(track_ids)
    for track in tracks:
        indexes = [item.frame_index for item in track.observations]
        if indexes != sorted(set(indexes)):
            raise ValueError(
                f"Track {track.track_id} frame indexes must be unique and increasing."
            )
        if not set(indexes).issubset(frame_set):
            raise ValueError(f"Track {track.track_id} refers to an unknown frame.")
        if any(item.gap_frames > config.gap_frame for item in track.observations):
            raise ValueError(
                f"Track {track.track_id} records a gap larger than configured."
            )
        for observation_index, item in enumerate(track.observations):
            stamp = frame_by_index[item.frame_index]
            if item.frame_name != stamp.frame_name or not _same_optional_time(
                item.time_ps,
                stamp.time_ps,
            ):
                raise ValueError(
                    f"Track {track.track_id} observation metadata does not match "
                    f"frame {item.frame_index}."
                )
            if item.match_status not in {
                "unavailable",
                "new",
                "secure",
                "ambiguous",
                "gap_bridge",
            }:
                raise ValueError(
                    f"Track {track.track_id} has invalid match status "
                    f"{item.match_status!r}."
                )
            expected_gap = 0
            if observation_index:
                previous = track.observations[observation_index - 1]
                expected_gap = (
                    frame_positions[item.frame_index]
                    - frame_positions[previous.frame_index]
                    - 1
                )
            if item.gap_frames != expected_gap:
                raise ValueError(
                    f"Track {track.track_id} observation in frame {item.frame_index} "
                    f"records gap_frames={item.gap_frames}, expected {expected_gap}."
                )
            if (
                item.gap_frames
                and config.max_gap_ps is not None
                and item.gap_time_ps is not None
                and item.gap_time_ps > config.max_gap_ps
            ):
                raise ValueError(
                    f"Track {track.track_id} records a physical gap larger than configured."
                )
    event_ids = [event.event_id for event in events]
    if len(event_ids) != len(set(event_ids)):
        raise ValueError("Tracking-state event IDs must be unique.")
    for event in events:
        unknown = set(event.track_ids) - known_tracks
        if unknown:
            raise ValueError(
                f"Event {event.event_id} refers to unknown track IDs: "
                + ", ".join(sorted(unknown))
            )
        if event.frame_index not in frame_set:
            raise ValueError(f"Event {event.event_id} refers to an unknown frame.")
        stamp = frame_by_index[event.frame_index]
        if event.frame_name != stamp.frame_name or not _same_optional_time(
            event.time_ps,
            stamp.time_ps,
        ):
            raise ValueError(
                f"Event {event.event_id} metadata does not match frame "
                f"{event.frame_index}."
            )
        gap_event = event.kind in {
            "gap",
            "type_change_unresolved",
            "phase_change_unresolved",
            "unresolved_across_gap",
        }
        if gap_event and not 1 <= event.gap_frames <= config.gap_frame:
            raise ValueError(
                f"Gap event {event.event_id} has an invalid gap_frames value."
            )
        if not gap_event and event.gap_frames != 0:
            raise ValueError(
                f"Non-gap event {event.event_id} must record gap_frames=0."
            )
        if event.persistence_frames < 0:
            raise ValueError(
                f"Event {event.event_id} has a negative persistence frame count."
            )


def _same_optional_time(left: float | None, right: float | None) -> bool:
    if left is None or right is None:
        return left is right
    return abs(left - right) <= max(1.0e-9, abs(left) * 1.0e-12)


def _track_info_text(
    data: TrackingResult | TargetSelection,
    tables: TrackTables,
) -> str:
    lifetimes = tables["cage_track.csv"]
    numeric_lifetimes = [
        float(row["lifetime_lower_ps"])
        for row in lifetimes
        if row.get("lifetime_lower_ps") is not None
    ]
    quality = tables["statistics/tracking_quality.csv"]
    quality_counts = Counter(str(row.get("match_status", "unavailable")) for row in quality)
    match_scores = [
        float(row["match_score"])
        for row in quality
        if row.get("match_score") is not None
    ]
    match_margins = [
        float(row["score_margin"])
        for row in quality
        if row.get("score_margin") is not None
    ]
    guest_residence = tables["statistics/guest_residence_lifetime.csv"]
    guest_lifetimes = _row_numeric_values(
        guest_residence,
        "residence_lifetime_lower_ps",
        "observed_span_ps",
    )
    occupancy_residence = tables["statistics/occupancy_state_lifetime.csv"]
    occupancy_lifetimes = _row_numeric_values(
        occupancy_residence,
        "residence_lifetime_lower_ps",
        "observed_span_ps",
    )
    event_counts = Counter(event.kind for event in data.events)
    transition_edges = sorted(
        tables["network/cage_transition_edges.csv"],
        key=lambda row: (
            -int(row.get("transition_count", 0)),
            str(row.get("source_cage_type", "")),
            str(row.get("destination_cage_type", "")),
        ),
    )
    cage_types = sorted({item.cage_type for item in data.observations})
    phases = sorted(
        {phase for item in data.observations for phase in item.phase_labels},
        key=str.casefold,
    )
    if isinstance(data, TargetSelection):
        target_value = data.target.value
        target_kind = data.target.kind
    else:
        target_value = "all"
        target_kind = "all"
    values: list[tuple[str, object]] = [
        ("target", target_value),
        ("target kind", target_kind),
        ("frames", len(data.frames)),
        ("tracks", len(data.tracks)),
        ("observations", len(data.observations)),
        ("configured gap frames", data.config.gap_frame),
        (
            "configured maximum gap (ps)",
            getattr(data.config, "max_gap_ps", None) or "none",
        ),
        ("cage types", ", ".join(cage_types) or "none"),
        ("phases", ", ".join(phases) or "none"),
    ]
    if numeric_lifetimes:
        values.extend(
            (
                ("cage lifetime lower bound min (ps)", f"{min(numeric_lifetimes):.9g}"),
                (
                    "cage lifetime lower bound mean (ps)",
                    f"{sum(numeric_lifetimes) / len(numeric_lifetimes):.9g}",
                ),
                ("cage lifetime lower bound max (ps)", f"{max(numeric_lifetimes):.9g}"),
            )
        )
    values.extend(
        (
            ("secure matches", quality_counts.get("secure", 0)),
            ("ambiguous matches", quality_counts.get("ambiguous", 0)),
            ("gap-bridge matches", quality_counts.get("gap_bridge", 0)),
        )
    )
    if match_scores:
        values.extend(
            (
                ("match score minimum", f"{min(match_scores):.9g}"),
                ("match score mean", f"{sum(match_scores) / len(match_scores):.9g}"),
            )
        )
    if match_margins:
        values.extend(
            (
                ("match score margin minimum", f"{min(match_margins):.9g}"),
                (
                    "match score margin mean",
                    f"{sum(match_margins) / len(match_margins):.9g}",
                ),
            )
        )
    duration_counts = Counter(str(row.get("duration_status")) for row in lifetimes)
    for status in sorted(duration_counts):
        values.append((f"{status} cage lifetimes", duration_counts[status]))
    if guest_residence:
        values.append(("guest residence count", len(guest_residence)))
    if guest_lifetimes:
        values.extend(
            (
                (
                    "guest residence lower bound mean (ps)",
                    f"{sum(guest_lifetimes) / len(guest_lifetimes):.9g}",
                ),
                (
                    "guest residence lower bound max (ps)",
                    f"{max(guest_lifetimes):.9g}",
                ),
            )
        )
    if occupancy_residence:
        values.append(("occupancy-state residence count", len(occupancy_residence)))
    if occupancy_lifetimes:
        values.extend(
            (
                (
                    "occupancy-state lower bound mean (ps)",
                    f"{sum(occupancy_lifetimes) / len(occupancy_lifetimes):.9g}",
                ),
                (
                    "occupancy-state lower bound max (ps)",
                    f"{max(occupancy_lifetimes):.9g}",
                ),
            )
        )
    for kind in (
        "birth",
        "death",
        "type_change",
        "phase_change",
        "type_change_unresolved",
        "phase_change_unresolved",
        "gap",
        "split",
        "merge",
        "split_candidate",
        "merge_candidate",
        "split_confirmed",
        "merge_confirmed",
        "guest_enter",
        "guest_exit",
        "guest_exchange",
        "occupancy_change",
        "empty_to_occupied",
        "occupied_to_empty",
        "single_to_multiple",
        "multiple_to_single",
        "unresolved_across_gap",
    ):
        values.append((f"{kind} events", event_counts.get(kind, 0)))
    lines = ["# SQQ Track", "", "| item | value |", "| --- | --- |"]
    lines.extend(f"| {_markdown(item)} | {_markdown(value)} |" for item, value in values)
    lines.extend(
        (
            "",
            "## Statistics files",
            "",
            *(f"- `statistics/{filename}`" for filename in _STATISTICS_TABLE_NAMES),
        )
    )
    if transition_edges:
        lines.extend(
            (
                "",
                "## Leading cage-type transitions",
                "",
                "| source | destination | count | tracks | rate (ns^-1) |",
                "| --- | --- | ---: | ---: | ---: |",
            )
        )
        for row in transition_edges[:10]:
            rate = row.get("rate_per_ns")
            rate_text = "" if rate is None else f"{float(rate):.9g}"
            lines.append(
                "| "
                + " | ".join(
                    (
                        _markdown(row.get("source_cage_type", "")),
                        _markdown(row.get("destination_cage_type", "")),
                        _markdown(row.get("transition_count", 0)),
                        _markdown(row.get("independent_track_count", 0)),
                        rate_text,
                    )
                )
                + " |"
            )
    lines.extend(
        (
            "",
            "## Cage-type transition network",
            "",
            "- `network/cage_transition_nodes.csv`",
            "- `network/cage_transition_edges.csv`",
            "- `network/cage_transition_network.png` (only when plotting is available and transitions exist)",
        )
    )
    lines.append("")
    return "\n".join(lines)


def _row_numeric_values(
    rows: Iterable[Mapping[str, object]],
    *fields: str,
) -> list[float]:
    output: list[float] = []
    for row in rows:
        value = next((row.get(field) for field in fields if row.get(field) is not None), None)
        if value is not None:
            output.append(float(value))
    return output


def _compact_cage_type(value: str) -> str:
    text = str(value).strip()
    matches = _CAGE_TERM.findall(text)
    if matches:
        compact = "".join(symbol + exponent for symbol, exponent in matches)
        residue = _CAGE_TERM.sub("", text).replace("-", "").replace("_", "")
        if residue:
            raise ValueError(f"Cannot compact cage type for a path: {value!r}.")
    else:
        compact = text.replace("^", "").replace("-", "").replace("_", "")
    return _safe_component(compact)


def _compact_object_id(value: object) -> str:
    text = str(value)
    match = _NUMERIC_SUFFIX.search(text)
    return str(int(match.group(1))) if match is not None else text


def _safe_component(value: str) -> str:
    if not _SAFE_COMPONENT.fullmatch(value) or value in {".", ".."}:
        raise ValueError(f"Unsafe output path component: {value!r}.")
    return value


def _safe_child(root: Path, component: str) -> Path:
    safe = _safe_component(component)
    root_resolved = root.resolve()
    child = (root_resolved / safe).resolve()
    if child.parent != root_resolved:
        raise ValueError(f"Output path escapes its tracking directory: {child}")
    return child


def _atomic_write_text(path: Path, text: str, *, encoding: str = "utf-8") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_write_csv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(fields), extrasaction="ignore"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: _csv_value(row.get(field)) for field in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _temporary_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{uuid4().hex}.tmp")


def _unique_existing_files(candidates: Iterable[Path]) -> list[Path]:
    output: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        if not candidate.is_file():
            continue
        resolved = candidate.resolve()
        identity = os.path.normcase(str(resolved))
        if identity not in seen:
            seen.add(identity)
            output.append(resolved)
    return output


def _mapping_sequence(value: object, label: str) -> tuple[Mapping[str, object], ...]:
    items = _sequence(value, label)
    if any(not isinstance(item, Mapping) for item in items):
        raise ValueError(f"Tracking-state {label} must contain mappings.")
    return tuple(item for item in items if isinstance(item, Mapping))


def _sequence(value: object, label: str) -> tuple[object, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Tracking-state {label} must be a sequence.")
    return tuple(value)


def _required_mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Tracking-state {label} must be a mapping.")
    return value


def _required(value: Mapping[str, object], key: str) -> object:
    if key not in value or value[key] is None:
        raise ValueError(f"Tracking-state field {key!r} is required.")
    return value[key]


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    result = float(value)
    if not (-float("inf") < result < float("inf")):
        raise ValueError("Tracking-state numeric values must be finite.")
    return result


def _optional_int(value: object) -> int | None:
    if value is None or value == "":
        return None
    return int(value)


def _strict_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"Tracking-state {label} must be true or false.")
    return value


def _csv_value(value: object) -> object:
    return "" if value is None else value


def _markdown(value: object) -> str:
    return str(value).replace("|", r"\|").replace("\r", " ").replace("\n", " ")
