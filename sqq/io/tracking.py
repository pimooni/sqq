from __future__ import annotations

"""State, normalized tables, and target bundles for cage tracking."""

from collections import Counter
import csv
import json
import os
from pathlib import Path
import re
import shutil
from typing import Iterable, Mapping, Sequence
from uuid import uuid4

import numpy as np

from .render import (
    SQQ_CAGE_GRO_NAME,
    SQQ_CAGE_MEMBERSHIP_NAME,
    SQQ_CAGE_XTC_NAME,
    SQQ_RENDER_DIRECTORY,
    SQQ_RENDER_SCRIPT_NAME,
    SqqCageBundle,
    vmd_script_text,
)
from ..core.tracking import (
    CageObservation,
    CageTrack,
    FrameStamp,
    TargetSelection,
    TargetSpec,
    TrackEvent,
    TrackingConfig,
    TrackingResult,
    event_rows,
    guest_residence_rows,
    lifetime_distribution_rows,
    lifetime_rows,
    observation_rows,
    population_rows,
    select_targets,
)


TRACK_DIRECTORY_NAME = "track"
TRACK_STATE_NAME = "track_state.json"
TRACK_RENDER_DIRECTORY = SQQ_RENDER_DIRECTORY
TRACK_GRO_NAME = "sqq_track.gro"
TRACK_XTC_NAME = "sqq_track.xtc"
TRACK_MEMBERSHIP_NAME = "sqq_track.membership.tsv"
TRACK_TCL_NAME = "sqq_track.vmd.tcl"

_STATE_FORMAT = "SQQ track state"
_STATE_VERSION = 2
_SUPPORTED_STATE_VERSIONS = {1, 2}
_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
_TRACK_ID = re.compile(r"^t[1-9][0-9]*$", re.IGNORECASE)
_CAGE_TERM = re.compile(r"([456])\^([0-9]+)")
_NUMERIC_SUFFIX = re.compile(r"(\d+)$")
_EVENT_KINDS = {
    "birth",
    "death",
    "type_change",
    "phase_change",
    "gap",
    "split",
    "merge",
}

_OBSERVATION_FIELDS = (
    "target",
    "track_id",
    "frame_index",
    "frame",
    "time_ps",
    "cage_id",
    "cage_type",
    "phase",
    "water_atomids",
    "water_count",
    "guest_ids",
    "guest_count",
    "center_x_nm",
    "center_y_nm",
    "center_z_nm",
    "match_jaccard",
    "match_shared_fraction",
    "match_center_distance_nm",
    "match_topology_similarity",
    "gap_frames",
)
_TRACK_FIELDS = (
    "target",
    "track_id",
    "first_frame_index",
    "last_frame_index",
    "first_time_ps",
    "last_time_ps",
    "end_time_ps",
    "observed_frames",
    "span_frames",
    "gap_frames",
    "observed_span_ps",
    "lifetime_ps",
    "duration_status",
    "left_censored",
    "right_censored",
    "initial_cage_type",
    "final_cage_type",
    "cage_types",
    "phases",
)
_LIFETIME_DISTRIBUTION_FIELDS = (
    "target",
    "lifetime_ps",
    "span_frames",
    "track_count",
    "fraction",
    "cumulative_fraction",
    "survival_fraction",
    "uncensored_count",
    "left_censored_count",
    "right_censored_count",
)
_EVENT_FIELDS = (
    "target",
    "event_id",
    "event",
    "frame_index",
    "frame",
    "time_ps",
    "source_track_ids",
    "destination_track_ids",
    "source_cage_types",
    "destination_cage_types",
    "source_phases",
    "destination_phases",
    "gap_frames",
    "censored",
)
_POPULATION_FIELDS = (
    "target",
    "frame_index",
    "frame",
    "time_ps",
    "group",
    "label",
    "cage_count",
)
_GUEST_RESIDENCE_FIELDS = (
    "target",
    "track_id",
    "guest_id",
    "episode",
    "start_frame_index",
    "end_frame_index",
    "start_time_ps",
    "end_time_ps",
    "observed_frames",
    "residence_time_ps",
    "left_censored",
    "right_censored",
)

_TABLE_SPECS = (
    ("cage_observation.csv", observation_rows, _OBSERVATION_FIELDS),
    ("cage_track.csv", lifetime_rows, _TRACK_FIELDS),
    (
        "lifetime_distribution.csv",
        lifetime_distribution_rows,
        _LIFETIME_DISTRIBUTION_FIELDS,
    ),
    ("cage_event.csv", event_rows, _EVENT_FIELDS),
    ("cage_population.csv", population_rows, _POPULATION_FIELDS),
    ("guest_residence.csv", guest_residence_rows, _GUEST_RESIDENCE_FIELDS),
)

__all__ = [
    "TRACK_DIRECTORY_NAME",
    "TRACK_GRO_NAME",
    "TRACK_MEMBERSHIP_NAME",
    "TRACK_RENDER_DIRECTORY",
    "TRACK_STATE_NAME",
    "TRACK_TCL_NAME",
    "TRACK_XTC_NAME",
    "deserialize_tracking_result",
    "discover_sqq_cage_bundle",
    "discover_sqq_cage_gro",
    "discover_track_state",
    "read_tracking_result",
    "rewrite_membership_track_ids",
    "serialize_tracking_result",
    "target_directory_name",
    "write_target_selection",
    "write_track_outputs",
    "write_tracking_result",
    "write_tracking_tables",
]


def serialize_tracking_result(result: TrackingResult) -> dict[str, object]:
    """Return the v2 JSON-safe tracking state."""
    if not isinstance(result, TrackingResult):
        raise TypeError("result must be a TrackingResult.")
    return {
        "format": _STATE_FORMAT,
        "version": _STATE_VERSION,
        "tracking_config": result.config.to_dict(),
        "frames": [_frame_to_dict(frame) for frame in result.frames],
        "tracks": [_track_to_dict(track) for track in result.tracks],
        "events": [_event_to_dict(event) for event in result.events],
    }


def deserialize_tracking_result(payload: Mapping[str, object]) -> TrackingResult:
    """Read v2 state and migrate the archived v1 representation."""
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
    _validate_result_identity(frames, tracks, events, config)
    return TrackingResult(frames=frames, tracks=tracks, events=events, config=config)


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
        "This source was not saved with Track state; rerun Analyze with SQQ 0.5.1."
    )


def discover_sqq_cage_bundle(
    source: str | Path | None = None,
    *,
    state_path: str | Path | None = None,
) -> SqqCageBundle:
    """Find the current underscore-named compact visualization bundle."""
    roots: list[Path] = []
    if source is not None:
        root = Path(source)
        if root.is_file():
            root = root.parent
        roots.extend((root / SQQ_RENDER_DIRECTORY, root))
    if state_path is not None:
        state = Path(state_path)
        roots.extend(
            (
                state.parent / SQQ_RENDER_DIRECTORY,
                state.parent.parent / SQQ_RENDER_DIRECTORY,
            )
        )
    if source is None and state_path is None:
        roots.extend(
            (
                Path.cwd() / SQQ_RENDER_DIRECTORY,
                Path.cwd() / TRACK_DIRECTORY_NAME / SQQ_RENDER_DIRECTORY,
            )
        )
    bundles: list[SqqCageBundle] = []
    seen: set[str] = set()
    for candidate in roots:
        render_dir = candidate.resolve()
        identity = os.path.normcase(str(render_dir))
        if identity in seen:
            continue
        seen.add(identity)
        gro = render_dir / SQQ_CAGE_GRO_NAME
        xtc = render_dir / SQQ_CAGE_XTC_NAME
        membership = render_dir / SQQ_CAGE_MEMBERSHIP_NAME
        script = render_dir / SQQ_RENDER_SCRIPT_NAME
        if all(path.is_file() for path in (gro, xtc, membership, script)):
            bundles.append(
                SqqCageBundle(
                    gro_path=gro,
                    script_path=script,
                    frame_count=_membership_frame_count(membership),
                    xtc_path=xtc,
                    membership_path=membership,
                    render_dir=render_dir,
                )
            )
    if len(bundles) == 1:
        return bundles[0]
    if len(bundles) > 1:
        raise ValueError(
            "Multiple SQQ render bundles were found: "
            + ", ".join(str(bundle.render_dir) for bundle in bundles)
        )
    raise FileNotFoundError(
        f"Cannot find {SQQ_RENDER_DIRECTORY}/{{{SQQ_CAGE_GRO_NAME}, "
        f"{SQQ_CAGE_XTC_NAME}, {SQQ_CAGE_MEMBERSHIP_NAME}, "
        f"{SQQ_RENDER_SCRIPT_NAME}}}."
    )


def discover_sqq_cage_gro(
    source: str | Path | None = None,
    *,
    state_path: str | Path | None = None,
) -> Path:
    bundle = discover_sqq_cage_bundle(source, state_path=state_path)
    if bundle.gro_path is None:
        raise FileNotFoundError("The SQQ render bundle has no topology GRO.")
    return bundle.gro_path


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


def write_tracking_tables(
    data: TrackingResult | TargetSelection,
    directory: str | Path,
) -> dict[str, Path]:
    """Write six normalized tables, including a distinct lifetime distribution."""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for filename, row_builder, fields in _TABLE_SPECS:
        path = root / filename
        _atomic_write_csv(path, row_builder(data), fields)
        written[filename] = path
    return written


def write_target_selection(
    selection: TargetSelection,
    track_root: str | Path,
    *,
    source_bundle: SqqCageBundle,
) -> Path:
    """Write one target and convert copied membership IDs to persistent IDs."""
    root = Path(track_root)
    directory = _safe_child(root, target_directory_name(selection.target))
    directory.mkdir(parents=True, exist_ok=True)
    write_tracking_tables(selection, directory)
    _atomic_write_text(directory / "track_info.md", _track_info_text(selection))

    gro_source, xtc_source, membership_source = _required_render_paths(source_bundle)
    render_dir = directory / TRACK_RENDER_DIRECTORY
    render_dir.mkdir(parents=True, exist_ok=True)
    _atomic_link_or_copy(gro_source, render_dir / TRACK_GRO_NAME)
    _atomic_link_or_copy(xtc_source, render_dir / TRACK_XTC_NAME)
    membership_target = render_dir / TRACK_MEMBERSHIP_NAME
    _atomic_copy(membership_source, membership_target)
    rewrite_membership_track_ids(membership_target, selection)
    _atomic_write_text(
        render_dir / TRACK_TCL_NAME,
        _target_vmd_script(selection),
        encoding="ascii",
    )
    (directory / TRACK_GRO_NAME).unlink(missing_ok=True)
    (directory / TRACK_TCL_NAME).unlink(missing_ok=True)
    return directory


def write_track_outputs(
    result: TrackingResult,
    outdir: str | Path,
    *,
    targets: str | Iterable[str] = "all",
    source: str | Path | None = None,
    source_bundle: SqqCageBundle | None = None,
    source_gro: str | Path | None = None,
) -> dict[str, Path]:
    """Write run state/tables and one self-contained bundle per target."""
    if source_bundle is None:
        if source_gro is not None:
            raise ValueError(
                "Track visualization requires the complete render bundle, "
                "not a standalone GRO."
            )
        source_bundle = discover_sqq_cage_bundle(
            source if source is not None else outdir
        )
    _required_render_paths(source_bundle)
    if source_bundle.frame_count != len(result.frames):
        raise ValueError(
            "Tracking state and render bundle frame counts differ: "
            f"{len(result.frames)} versus {source_bundle.frame_count}."
        )

    root = Path(outdir) / TRACK_DIRECTORY_NAME
    root.mkdir(parents=True, exist_ok=True)
    written = write_tracking_tables(result, root)
    written[TRACK_STATE_NAME] = write_tracking_result(result, root / TRACK_STATE_NAME)
    for selection in select_targets(result, targets):
        directory = write_target_selection(
            selection,
            root,
            source_bundle=source_bundle,
        )
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
        "gap_frames": item.gap_frames,
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
        match_center_distance_nm=_optional_float(
            value.get("match_center_distance_nm")
        ),
        match_topology_similarity=_optional_float(
            value.get("match_topology_similarity")
        ),
        gap_frames=gap,
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
        "censored": event.censored,
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
        censored=_strict_bool(value.get("censored"), "censored"),
    )


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
        if event.kind == "gap" and not 1 <= event.gap_frames <= config.gap_frame:
            raise ValueError(
                f"Gap event {event.event_id} has an invalid gap_frames value."
            )
        if event.kind != "gap" and event.gap_frames != 0:
            raise ValueError(
                f"Non-gap event {event.event_id} must record gap_frames=0."
            )


def _track_info_text(selection: TargetSelection) -> str:
    lifetimes = lifetime_rows(selection)
    numeric_lifetimes = [
        float(row["lifetime_ps"])
        for row in lifetimes
        if row.get("lifetime_ps") is not None
    ]
    event_counts = Counter(event.kind for event in selection.events)
    cage_types = sorted({item.cage_type for item in selection.observations})
    phases = sorted(
        {phase for item in selection.observations for phase in item.phase_labels},
        key=str.casefold,
    )
    values: list[tuple[str, object]] = [
        ("target", selection.target.value),
        ("target kind", selection.target.kind),
        ("frames", len(selection.frames)),
        ("tracks", len(selection.tracks)),
        ("observations", len(selection.observations)),
        ("configured gap frames", selection.config.gap_frame),
        ("cage types", ", ".join(cage_types) or "none"),
        ("phases", ", ".join(phases) or "none"),
    ]
    if numeric_lifetimes:
        values.extend(
            (
                ("lifetime min (ps)", f"{min(numeric_lifetimes):.9g}"),
                (
                    "lifetime mean (ps)",
                    f"{sum(numeric_lifetimes) / len(numeric_lifetimes):.9g}",
                ),
                ("lifetime max (ps)", f"{max(numeric_lifetimes):.9g}"),
            )
        )
    for kind in (
        "birth",
        "death",
        "type_change",
        "phase_change",
        "gap",
        "split",
        "merge",
    ):
        values.append((f"{kind} events", event_counts.get(kind, 0)))
    lines = ["# SQQ Track", "", "| item | value |", "| --- | --- |"]
    lines.extend(f"| {_markdown(item)} | {_markdown(value)} |" for item, value in values)
    lines.append("")
    return "\n".join(lines)


def _target_vmd_script(selection: TargetSelection) -> str:
    target_name = target_directory_name(selection.target)
    script = vmd_script_text(
        gro_filename=TRACK_GRO_NAME,
        xtc_filename=TRACK_XTC_NAME,
        membership_filename=TRACK_MEMBERSHIP_NAME,
        molecule_name=f"SQQ track {target_name}",
        render_kind="track",
    ).rstrip()
    return (
        script
        + "\n\n# Default lifecycle view for this tracking target.\n"
        + "\n".join(_selection_commands(selection))
        + "\n"
    )


def _selection_commands(selection: TargetSelection) -> list[str]:
    if selection.target.kind == "all":
        return ["sqq show cage all"]
    object_ids = sorted(
        {track.track_id for track in selection.tracks},
        key=lambda value: int(value[1:]),
    )
    if not object_ids:
        return [
            "set ::SQQ::active_families {}",
            "set ::SQQ::custom_show_active 1",
            "::SQQ::render_current",
            f'puts "SQQ track target {target_directory_name(selection.target)} matched no cages."',
        ]
    return [
        "sqq show cage " + " ".join(object_ids[start : start + 100])
        for start in range(0, len(object_ids), 100)
    ]


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


def _required_render_paths(bundle: SqqCageBundle) -> tuple[Path, Path, Path]:
    paths = (bundle.gro_path, bundle.xtc_path, bundle.membership_path)
    labels = ("topology GRO", "XTC trajectory", "membership TSV")
    missing = [
        label
        for label, path in zip(labels, paths)
        if path is None or not Path(path).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete SQQ render bundle; missing " + ", ".join(missing) + "."
        )
    return tuple(Path(path).resolve() for path in paths)  # type: ignore[return-value]


def _membership_frame_count(path: Path) -> int:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "record" not in reader.fieldnames:
            raise ValueError(f"Invalid SQQ membership TSV: {path}")
        return sum(1 for row in reader if row.get("record") == "F")


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


def _atomic_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and source.samefile(target):
        return
    temporary = _temporary_path(target)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_link_or_copy(source: Path, target: Path) -> None:
    source = Path(source).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and source.samefile(target):
        return
    temporary = _temporary_path(target)
    try:
        try:
            os.link(source, temporary)
        except OSError:
            temporary.unlink(missing_ok=True)
            shutil.copyfile(source, temporary)
        os.replace(temporary, target)
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


def _strict_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"Tracking-state {label} must be true or false.")
    return value


def _csv_value(value: object) -> object:
    return "" if value is None else value


def _markdown(value: object) -> str:
    return str(value).replace("|", r"\|").replace("\r", " ").replace("\n", " ")
