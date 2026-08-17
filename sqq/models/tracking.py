"""Persistent cage-tracking data contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Iterable, Literal, Mapping

import numpy as np


EventKind = Literal[
    "birth",
    "death",
    "type_change",
    "phase_change",
    "gap",
    "split",
    "merge",
]
TargetKind = Literal["all", "cage_type", "phase", "track"]
Row = dict[str, object]

__all__ = [
    "CageObservation",
    "CageTrack",
    "EventKind",
    "FrameStamp",
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


@dataclass(frozen=True)
class TrackingConfig:
    """Thresholds used to link cages between selected frames."""

    min_jaccard: float = 0.50
    min_shared_fraction: float = 0.60
    min_shared_waters: int = 3
    max_center_distance_nm: float | None = None
    gap_frame: int = 0
    guest_tiebreak: bool = True

    def __post_init__(self) -> None:
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
        if not isinstance(self.guest_tiebreak, bool):
            raise ValueError("guest_tiebreak must be true or false.")

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
            "guest_tiebreak",
        }
        unknown = sorted(set(values).difference(supported))
        if unknown:
            raise ValueError(
                "Unsupported tracking configuration field(s): " + ", ".join(unknown)
            )
        maximum = values.get("max_center_distance_nm")
        return cls(
            min_jaccard=float(values.get("min_jaccard", 0.50)),
            min_shared_fraction=float(values.get("min_shared_fraction", 0.60)),
            min_shared_waters=values.get("min_shared_waters", 3),  # type: ignore[arg-type]
            max_center_distance_nm=(
                None if maximum in (None, "") else float(maximum)
            ),
            gap_frame=values.get("gap_frame", 0),  # type: ignore[arg-type]
            guest_tiebreak=_strict_mapping_bool(
                values.get("guest_tiebreak", True), "guest_tiebreak"
            ),
        )

    def to_dict(self) -> Row:
        return {
            "min_jaccard": self.min_jaccard,
            "min_shared_fraction": self.min_shared_fraction,
            "min_shared_waters": self.min_shared_waters,
            "max_center_distance_nm": self.max_center_distance_nm,
            "gap_frame": self.gap_frame,
            "guest_tiebreak": self.guest_tiebreak,
        }


@dataclass(frozen=True)
class FrameStamp:
    frame_index: int
    frame_name: str
    time_ps: float | None
    source: str


@dataclass(frozen=True)
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


@dataclass(frozen=True)
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


@dataclass(frozen=True)
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
    match_center_distance_nm: float | None = None
    match_topology_similarity: float | None = None
    gap_frames: int = 0


@dataclass(frozen=True)
class CageTrack:
    track_id: str
    observations: tuple[CageObservation, ...]
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


@dataclass(frozen=True)
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
    censored: bool = False

    @property
    def track_ids(self) -> tuple[str, ...]:
        return _unique(self.source_track_ids + self.destination_track_ids)


@dataclass(frozen=True)
class TrackingResult:
    frames: tuple[FrameStamp, ...]
    tracks: tuple[CageTrack, ...]
    events: tuple[TrackEvent, ...]
    config: TrackingConfig = field(default_factory=TrackingConfig)

    @property
    def observations(self) -> tuple[CageObservation, ...]:
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
    frames: tuple[FrameStamp, ...]
    tracks: tuple[CageTrack, ...]
    events: tuple[TrackEvent, ...]
    config: TrackingConfig = field(default_factory=TrackingConfig)

    @property
    def observations(self) -> tuple[CageObservation, ...]:
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
    return item.frame_index, _track_number(item.track_id), item.local_cage_id


def _as_sequence(value: object, label: str) -> tuple[object, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Tracking snapshot {label} must be a sequence.")
    return tuple(value)


def _strict_mapping_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false.")
    return value
