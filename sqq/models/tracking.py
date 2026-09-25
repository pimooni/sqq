"""Persistent cage-tracking data contracts."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, fields
from functools import cached_property
import re
from typing import Iterable, Literal, Mapping

import numpy as np


EventKind = Literal[
    "birth",
    "death",
    "type_change",
    "phase_change",
    "type_change_unresolved",
    "phase_change_unresolved",
    "gap",
    # ``split``/``merge`` remain readable for pre-0.5.6 track-state files.
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
]
MatchStatus = Literal["unavailable", "new", "secure", "ambiguous", "gap_bridge"]
TargetKind = Literal["all", "cage_type", "phase", "track"]
Row = dict[str, object]

__all__ = [
    "CageObservation",
    "CageTrack",
    "EventKind",
    "FrameStamp",
    "MatchStatus",
    "Row",
    "TargetKind",
    "TargetSelection",
    "TargetSpec",
    "TrackCageSnapshot",
    "TrackEvent",
    "TrackFrameSnapshot",
    "TrackingConfig",
    "TrackingResult",
]

_TRACK_PATTERN = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)
_PHASE_ALIASES = {
    "i": "sI",
    "si": "sI",
    "ii": "sII",
    "sii": "sII",
    "h": "sH",
    "sh": "sH",
    "boundary": "boundary",
    "unclassified": "unclassified",
    "ambiguous": "ambiguous",
    "isolated": "isolated",
    "unassigned": "unassigned",
    "mixed": "mixed",
}
_PHASE_ORDER = {
    "sI": 0,
    "sII": 1,
    "sH": 2,
    "boundary": 3,
    "ambiguous": 4,
    "unclassified": 5,
    "isolated": 6,
    "unassigned": 7,
}


def _pickle_compatible_slots(cls):
    """Keep slotted records readable from pre-0.5.7 dataclass pickles.

    Unslotted frozen dataclasses were pickled as name-to-value dictionaries.
    Python's generated pickle hook for a frozen slotted dataclass instead
    expects a positional sequence and silently assigns dictionary keys as
    values.  Store names going forward and accept both representations.
    """
    names = tuple(item.name for item in fields(cls))

    def __getstate__(self):
        return {name: getattr(self, name) for name in names}

    def __setstate__(self, state):
        if isinstance(state, Mapping):
            missing = [name for name in names if name not in state]
            if missing:
                raise ValueError(
                    "Pickled tracking record is missing field(s): "
                    + ", ".join(missing)
                )
            values = tuple(state[name] for name in names)
        elif isinstance(state, (tuple, list)):
            if len(state) != len(names):
                raise ValueError(
                    f"Pickled tracking record has {len(state)} values; "
                    f"expected {len(names)}."
                )
            values = tuple(state)
        else:
            raise TypeError("Pickled tracking record state must be a mapping or sequence.")
        for name, value in zip(names, values):
            object.__setattr__(self, name, value)

    cls.__getstate__ = __getstate__
    cls.__setstate__ = __setstate__
    return cls


@dataclass(frozen=True)
class TrackingConfig:
    """Thresholds used to link cages between selected frames."""

    min_jaccard: float = 0.50
    min_shared_fraction: float = 0.60
    min_shared_waters: int = 3
    max_center_distance_nm: float | None = None
    gap_frame: int = 0
    max_gap_ps: float | None = None
    guest_tiebreak: bool = True
    ambiguity_score_margin: float = 1.0
    near_threshold_tolerance: float = 0.05

    def __post_init__(self) -> None:
        # YAML ``true``/``false`` must not be silently promoted to 1.0/0.0.
        for name in (
            "min_jaccard",
            "min_shared_fraction",
            "max_center_distance_nm",
            "max_gap_ps",
            "ambiguity_score_margin",
            "near_threshold_tolerance",
        ):
            if isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a number, not true/false.")
        for name in ("min_jaccard", "min_shared_fraction"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be finite and between 0 and 1.")
            object.__setattr__(self, name, value)
        if isinstance(self.min_shared_waters, bool):
            raise ValueError("min_shared_waters must be a positive integer.")
        shared_waters = int(self.min_shared_waters)
        if shared_waters != self.min_shared_waters or shared_waters < 1:
            raise ValueError("min_shared_waters must be a positive integer.")
        object.__setattr__(self, "min_shared_waters", shared_waters)
        if self.max_center_distance_nm is not None:
            maximum = float(self.max_center_distance_nm)
            if not np.isfinite(maximum) or maximum <= 0.0:
                raise ValueError(
                    "max_center_distance_nm must be finite and positive when provided."
                )
            object.__setattr__(self, "max_center_distance_nm", maximum)
        if isinstance(self.gap_frame, bool):
            raise ValueError("gap_frame must be a nonnegative integer.")
        gap = int(self.gap_frame)
        if gap != self.gap_frame or gap < 0:
            raise ValueError("gap_frame must be a nonnegative integer.")
        object.__setattr__(self, "gap_frame", gap)
        if self.max_gap_ps is not None:
            maximum_gap = float(self.max_gap_ps)
            if not np.isfinite(maximum_gap) or maximum_gap <= 0.0:
                raise ValueError("max_gap_ps must be finite and positive when provided.")
            object.__setattr__(self, "max_gap_ps", maximum_gap)
        if not isinstance(self.guest_tiebreak, bool):
            raise ValueError("guest_tiebreak must be true or false.")
        ambiguity_margin = float(self.ambiguity_score_margin)
        if not np.isfinite(ambiguity_margin) or ambiguity_margin < 0.0:
            raise ValueError("ambiguity_score_margin must be finite and nonnegative.")
        object.__setattr__(self, "ambiguity_score_margin", ambiguity_margin)
        tolerance = float(self.near_threshold_tolerance)
        if not np.isfinite(tolerance) or not 0.0 <= tolerance <= 1.0:
            raise ValueError(
                "near_threshold_tolerance must be finite and between 0 and 1."
            )
        object.__setattr__(self, "near_threshold_tolerance", tolerance)

    @classmethod
    def from_mapping(cls, values: Mapping[str, object] | None) -> "TrackingConfig":
        """Build settings while rejecting keys that would otherwise be inert."""
        if values is None:
            return cls()
        if not isinstance(values, Mapping):
            raise TypeError("tracking configuration must be a mapping.")
        supported = {
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
        unknown = sorted(set(values).difference(supported))
        if unknown:
            raise ValueError(
                "Unsupported tracking configuration field(s): " + ", ".join(unknown)
            )
        for name in (
            "min_jaccard",
            "min_shared_fraction",
            "max_center_distance_nm",
            "max_gap_ps",
            "ambiguity_score_margin",
            "near_threshold_tolerance",
        ):
            if isinstance(values.get(name), bool):
                raise ValueError(f"{name} must be a number, not true/false.")
        maximum = values.get("max_center_distance_nm")
        maximum_gap = values.get("max_gap_ps")
        return cls(
            min_jaccard=float(values.get("min_jaccard", 0.50)),
            min_shared_fraction=float(values.get("min_shared_fraction", 0.60)),
            min_shared_waters=values.get("min_shared_waters", 3),  # type: ignore[arg-type]
            max_center_distance_nm=(
                None if maximum in (None, "") else float(maximum)
            ),
            gap_frame=values.get("gap_frame", 0),  # type: ignore[arg-type]
            max_gap_ps=(
                None if maximum_gap in (None, "") else float(maximum_gap)
            ),
            guest_tiebreak=_strict_mapping_bool(
                values.get("guest_tiebreak", True), "guest_tiebreak"
            ),
            ambiguity_score_margin=float(values.get("ambiguity_score_margin", 1.0)),
            near_threshold_tolerance=float(
                values.get("near_threshold_tolerance", 0.05)
            ),
        )

    def to_dict(self) -> Row:
        return {
            "min_jaccard": self.min_jaccard,
            "min_shared_fraction": self.min_shared_fraction,
            "min_shared_waters": self.min_shared_waters,
            "max_center_distance_nm": self.max_center_distance_nm,
            "gap_frame": self.gap_frame,
            "max_gap_ps": self.max_gap_ps,
            "guest_tiebreak": self.guest_tiebreak,
            "ambiguity_score_margin": self.ambiguity_score_margin,
            "near_threshold_tolerance": self.near_threshold_tolerance,
        }


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class FrameStamp:
    frame_index: int
    frame_name: str
    time_ps: float | None
    source: str


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class TrackCageSnapshot:
    """JSON-safe representation of one complete cage."""

    local_cage_id: str
    cage_type: str
    phase_labels: tuple[str, ...]
    water_atomids: tuple[int, ...]
    center: tuple[float, float, float]
    topology: tuple[int, ...]
    guest_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        local_id = str(self.local_cage_id).strip()
        cage_type = str(self.cage_type).strip()
        if not local_id or not cage_type:
            raise ValueError("Snapshot cage ID and type must be nonempty.")
        center = tuple(float(value) for value in self.center)
        if len(center) != 3 or any(not np.isfinite(value) for value in center):
            raise ValueError(f"Snapshot cage {local_id} has an invalid center.")
        waters = tuple(sorted(int(value) for value in self.water_atomids))
        if len(waters) != len(set(waters)):
            raise ValueError(f"Snapshot cage {local_id} repeats water atom IDs.")
        phases = _normalized_phases(self.phase_labels)
        topology = tuple(sorted(int(value) for value in self.topology))
        guests = tuple(sorted(set(str(value) for value in self.guest_ids)))
        object.__setattr__(self, "local_cage_id", local_id)
        object.__setattr__(self, "cage_type", cage_type)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "water_atomids", waters)
        object.__setattr__(self, "phase_labels", phases)
        object.__setattr__(self, "topology", topology)
        object.__setattr__(self, "guest_ids", guests)

    def to_dict(self) -> Row:
        return {
            "local_cage_id": self.local_cage_id,
            "cage_type": self.cage_type,
            "phase_labels": list(self.phase_labels),
            "water_atomids": list(self.water_atomids),
            "center": list(self.center),
            "topology": list(self.topology),
            "guest_ids": list(self.guest_ids),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "TrackCageSnapshot":
        center = tuple(float(item) for item in _as_sequence(value.get("center"), "center"))
        if len(center) != 3:
            raise ValueError("TrackCageSnapshot.center must contain three coordinates.")
        return cls(
            local_cage_id=str(value.get("local_cage_id", "")),
            cage_type=str(value.get("cage_type", "")),
            phase_labels=tuple(
                str(item)
                for item in _as_sequence(value.get("phase_labels", ()), "phase_labels")
            ),
            water_atomids=tuple(
                int(item)
                for item in _as_sequence(value.get("water_atomids", ()), "water_atomids")
            ),
            center=(center[0], center[1], center[2]),
            topology=tuple(
                int(item)
                for item in _as_sequence(value.get("topology", ()), "topology")
            ),
            guest_ids=tuple(
                str(item)
                for item in _as_sequence(value.get("guest_ids", ()), "guest_ids")
            ),
        )


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class TrackFrameSnapshot:
    """Frame metadata and cages consumed by the serial tracker."""

    frame_index: int
    frame_name: str
    time_ps: float | None
    source: str
    box: tuple[float, ...] | None
    cages: tuple[TrackCageSnapshot, ...]

    def __post_init__(self) -> None:
        index = int(self.frame_index)
        if index < 0:
            raise ValueError("Snapshot frame_index must be nonnegative.")
        time = None if self.time_ps is None else float(self.time_ps)
        if time is not None and not np.isfinite(time):
            raise ValueError("Snapshot time must be finite when provided.")
        box = None if self.box is None else tuple(float(value) for value in self.box)
        if box is not None and (
            len(box) < 3 or any(not np.isfinite(value) for value in box)
        ):
            raise ValueError("Snapshot box must contain finite box values.")
        cages = tuple(self.cages)
        local_ids = [cage.local_cage_id for cage in cages]
        if len(local_ids) != len(set(local_ids)):
            raise ValueError(f"Snapshot cage IDs are not unique in frame {self.frame_name}.")
        object.__setattr__(self, "frame_index", index)
        object.__setattr__(self, "time_ps", time)
        object.__setattr__(self, "box", box)
        object.__setattr__(self, "cages", cages)

    def to_dict(self) -> Row:
        return {
            "frame_index": self.frame_index,
            "frame_name": self.frame_name,
            "time_ps": self.time_ps,
            "source": self.source,
            "box": None if self.box is None else list(self.box),
            "cages": [cage.to_dict() for cage in self.cages],
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "TrackFrameSnapshot":
        raw_box = value.get("box")
        raw_cages = _as_sequence(value.get("cages", ()), "cages")
        if any(not isinstance(item, Mapping) for item in raw_cages):
            raise ValueError("Each snapshot cage must be a mapping.")
        return cls(
            frame_index=int(value.get("frame_index", -1)),
            frame_name=str(value.get("frame_name", "")),
            time_ps=None if value.get("time_ps") is None else float(value["time_ps"]),
            source=str(value.get("source", "")),
            box=(
                None
                if raw_box is None
                else tuple(float(item) for item in _as_sequence(raw_box, "box"))
            ),
            cages=tuple(
                TrackCageSnapshot.from_dict(item)
                for item in raw_cages
                if isinstance(item, Mapping)
            ),
        )


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class CageObservation:
    track_id: str
    frame_index: int
    frame_name: str
    time_ps: float | None
    local_cage_id: str
    cage_type: str
    phase: str
    phase_labels: tuple[str, ...]
    water_atomids: tuple[int, ...]
    center: tuple[float, float, float]
    topology: tuple[int, ...]
    guest_ids: tuple[str, ...]
    match_jaccard: float | None = None
    match_shared_fraction: float | None = None
    match_shared_waters: int | None = None
    match_center_distance_nm: float | None = None
    match_topology_similarity: float | None = None
    match_candidate_count: int = 0
    match_score: float | None = None
    match_runner_up_score: float | None = None
    match_score_margin: float | None = None
    match_jaccard_margin: float | None = None
    match_shared_fraction_margin: float | None = None
    match_shared_waters_margin: int | None = None
    match_center_distance_margin_nm: float | None = None
    match_near_threshold: bool = False
    match_status: MatchStatus = "unavailable"
    # ``computed`` for live diagnostics, ``migrated_partial`` when thresholds
    # and partial evidence were both archived, and empty for births or legacy
    # observations whose threshold-dependent diagnosis is not recoverable.
    match_diagnostic_source: str = ""
    gap_frames: int = 0
    gap_time_ps: float | None = None


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class CageTrack:
    track_id: str
    observations: Sequence[CageObservation]
    left_censored: bool
    right_censored: bool

    @property
    def first(self) -> CageObservation:
        return self.observations[0]

    @property
    def last(self) -> CageObservation:
        return self.observations[-1]

    @property
    def gap_frames(self) -> int:
        return sum(item.gap_frames for item in self.observations)


@_pickle_compatible_slots
@dataclass(frozen=True, slots=True)
class TrackEvent:
    event_id: str
    kind: EventKind
    frame_index: int
    frame_name: str
    time_ps: float | None
    source_track_ids: tuple[str, ...] = ()
    destination_track_ids: tuple[str, ...] = ()
    source_cage_types: tuple[str, ...] = ()
    destination_cage_types: tuple[str, ...] = ()
    source_phases: tuple[str, ...] = ()
    destination_phases: tuple[str, ...] = ()
    gap_frames: int = 0
    gap_time_ps: float | None = None
    censored: bool = False
    guest_ids_entered: tuple[str, ...] = ()
    guest_ids_exited: tuple[str, ...] = ()
    source_occupancy: str = ""
    destination_occupancy: str = ""
    guest_composition_before: tuple[str, ...] = ()
    guest_composition_after: tuple[str, ...] = ()
    evidence_status: str = ""
    water_conservation: float | None = None
    center_distance_nm: float | None = None
    persistence_frames: int = 0

    @property
    def track_ids(self) -> tuple[str, ...]:
        return _unique(self.source_track_ids + self.destination_track_ids)


@dataclass(frozen=True)
class TrackingResult:
    frames: Sequence[FrameStamp]
    tracks: Sequence[CageTrack]
    events: Sequence[TrackEvent]
    config: TrackingConfig = field(default_factory=TrackingConfig)
    source_provenance: Mapping[str, object] = field(default_factory=dict)
    _observation_view: Sequence[CageObservation] | None = field(
        default=None, repr=False, compare=False
    )

    # The frozen instance never changes, so the sorted view is computed once
    # instead of on every table builder that consumes it.
    @cached_property
    def observations(self) -> Sequence[CageObservation]:
        if self._observation_view is not None:
            return self._observation_view
        rows = (item for track in self.tracks for item in track.observations)
        return tuple(sorted(rows, key=_observation_sort_key))


@dataclass(frozen=True)
class TargetSpec:
    raw: str
    kind: TargetKind
    value: str

    @property
    def key(self) -> str:
        return "all" if self.kind == "all" else f"{self.kind}_{self.value}"


@dataclass(frozen=True)
class TargetSelection:
    target: TargetSpec
    frames: Sequence[FrameStamp]
    tracks: Sequence[CageTrack]
    events: Sequence[TrackEvent]
    config: TrackingConfig = field(default_factory=TrackingConfig)
    _observation_view: Sequence[CageObservation] | None = field(
        default=None, repr=False, compare=False
    )

    @cached_property
    def observations(self) -> Sequence[CageObservation]:
        if self._observation_view is not None:
            return self._observation_view
        rows = (item for track in self.tracks for item in track.observations)
        return tuple(sorted(rows, key=_observation_sort_key))


def _canonical_phase(value: str) -> str:
    text = str(value).strip()
    return _PHASE_ALIASES.get(text.casefold(), text or "unassigned")


def _normalized_phases(values: Iterable[str]) -> tuple[str, ...]:
    phases = {_canonical_phase(value) for value in values if str(value).strip()}
    if not phases:
        phases.add("unassigned")
    return tuple(sorted(phases, key=lambda value: (_PHASE_ORDER.get(value, 100), value)))


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values))


def _track_number(track_id: str) -> int:
    match = _TRACK_PATTERN.fullmatch(track_id)
    if match is None:
        raise ValueError(f"Invalid persistent track ID: {track_id!r}.")
    return int(match.group(1))


def _observation_sort_key(item: CageObservation) -> tuple[int, int, str]:
    return (
        item.frame_index,
        _track_number(item.track_id),
        item.local_cage_id,
    )


def _as_sequence(value: object, label: str) -> tuple[object, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Tracking snapshot {label} must be a sequence.")
    return tuple(value)


def _strict_mapping_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false.")
    return value
