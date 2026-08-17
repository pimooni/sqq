"""Deterministic, streaming cage tracking across analyzed frames."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from math import inf
import re
from typing import Iterable, Mapping, Sequence

import numpy as np

from .cage import canonical_cage_type
from .pbc import distance
from ..models import Cage, FrameResult, Ring
from ..models.tracking import (
    CageObservation,
    CageTrack,
    EventKind,
    FrameStamp,
    Row,
    TargetSelection,
    TargetSpec,
    TrackCageSnapshot,
    TrackEvent,
    TrackFrameSnapshot,
    TrackingConfig,
    TrackingResult,
)


__all__ = [
    "CageObservation",
    "CageTrack",
    "EventKind",
    "FrameStamp",
    "TargetSelection",
    "TargetSpec",
    "TrackCageSnapshot",
    "TrackEvent",
    "TrackFrameSnapshot",
    "TrackingAccumulator",
    "TrackingConfig",
    "TrackingResult",
    "event_rows",
    "guest_residence_rows",
    "lifetime_distribution_rows",
    "lifetime_rows",
    "observation_rows",
    "parse_targets",
    "population_rows",
    "select_targets",
    "snapshot_from_frame_result",
    "track_cages",
    "track_snapshots",
]

_TRACK_PATTERN = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)
_COMPACT_GENERIC_CAGE_PATTERN = re.compile(r"^4([0-9]+)5([0-9]+)6([0-9]+)$")
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
class _CageState:
    local_cage_id: str
    cage_type: str
    phase: str
    phase_labels: tuple[str, ...]
    water_atomids: frozenset[int]
    center: tuple[float, float, float]
    topology: tuple[int, ...]
    guest_ids: tuple[str, ...]


@dataclass(frozen=True)
class _TrackedState:
    state: _CageState
    track_id: str
    last_position: int


@dataclass(frozen=True)
class _Candidate:
    previous_index: int
    current_index: int
    jaccard: float
    shared_fraction: float
    center_distance_nm: float | None
    topology_similarity: float
    guest_similarity: float
    gap_frames: int


class TrackingAccumulator:
    """Incrementally track snapshots without retaining the snapshot sequence."""

    def __init__(self, config: TrackingConfig | None = None) -> None:
        self.config = config or TrackingConfig()
        self._frames: list[FrameStamp] = []
        self._track_observations: dict[str, list[CageObservation]] = {}
        self._events: list[TrackEvent] = []
        self._active: list[_TrackedState] = []
        self._expired: set[str] = set()
        self._next_track_number = 1
        self._finished = False

    def add(self, snapshot: TrackFrameSnapshot) -> None:
        if self._finished:
            raise RuntimeError("Cannot add a snapshot after tracking was finalized.")
        self._validate_next(snapshot)
        stamp = _snapshot_stamp(snapshot)
        position = len(self._frames)
        self._frames.append(stamp)
        current_states = _snapshot_states(snapshot)
        if position == 0:
            self._add_first_frame(stamp, current_states)
            return

        box = None if snapshot.box is None else np.asarray(snapshot.box, dtype=float)
        candidates = _match_candidates(
            self._active,
            current_states,
            box,
            self.config,
            position,
        )
        assignments = _maximum_weight_assign(
            candidates, self._active, current_states, self.config
        )
        evidence_by_current = {
            candidate.current_index: candidate for candidate in assignments
        }
        previous_by_current = {
            candidate.current_index: candidate.previous_index
            for candidate in assignments
        }
        current_track_ids: dict[int, str] = {}

        for current_index, state in enumerate(current_states):
            previous_index = previous_by_current.get(current_index)
            if previous_index is None:
                track_id = self._new_track_id()
                self._track_observations[track_id] = []
                self._events.append(
                    _event(
                        "birth",
                        stamp,
                        destination=(track_id,),
                        destination_types=(state.cage_type,),
                        destination_phases=state.phase_labels,
                    )
                )
            else:
                track_id = self._active[previous_index].track_id
            current_track_ids[current_index] = track_id
            evidence = evidence_by_current.get(current_index)
            self._track_observations[track_id].append(
                _observation(track_id, stamp, state, evidence=evidence)
            )
            if evidence is not None and evidence.gap_frames:
                previous = self._active[evidence.previous_index]
                self._events.append(
                    _event(
                        "gap",
                        stamp,
                        source=(track_id,),
                        destination=(track_id,),
                        source_types=(previous.state.cage_type,),
                        destination_types=(state.cage_type,),
                        source_phases=previous.state.phase_labels,
                        destination_phases=state.phase_labels,
                        gap_frames=evidence.gap_frames,
                    )
                )

        self._add_change_events(
            stamp, current_states, current_track_ids, assignments
        )
        self._add_split_merge_events(
            stamp, current_states, current_track_ids, candidates
        )

        matched_previous = {candidate.previous_index for candidate in assignments}
        dormant: list[_TrackedState] = []
        for previous_index, previous in enumerate(self._active):
            if previous_index in matched_previous:
                continue
            if position - previous.last_position > self.config.gap_frame:
                self._expire(previous)
            else:
                dormant.append(previous)

        current_active = [
            _TrackedState(
                state=state,
                track_id=current_track_ids[index],
                last_position=position,
            )
            for index, state in enumerate(current_states)
        ]
        self._active = current_active + sorted(
            dormant, key=lambda item: _track_number(item.track_id)
        )

    def result(self) -> TrackingResult:
        self._finished = True
        active_ids = {item.track_id for item in self._active}
        tracks = tuple(
            CageTrack(
                track_id=track_id,
                observations=tuple(self._track_observations[track_id]),
                left_censored=bool(
                    self._frames
                    and self._track_observations[track_id][0].frame_index
                    == self._frames[0].frame_index
                ),
                right_censored=track_id in active_ids,
            )
            for track_id in sorted(self._track_observations, key=_track_number)
        )
        events = tuple(
            replace(event, event_id=f"e{index}")
            for index, event in enumerate(self._events, start=1)
        )
        return TrackingResult(
            frames=tuple(self._frames),
            tracks=tracks,
            events=events,
            config=self.config,
        )

    def _validate_next(self, snapshot: TrackFrameSnapshot) -> None:
        if not self._frames:
            return
        previous = self._frames[-1]
        if snapshot.frame_index <= previous.frame_index:
            raise ValueError("Tracking snapshots must have strictly increasing frame indexes.")
        if snapshot.frame_index != previous.frame_index + 1:
            raise ValueError(
                "Tracking snapshot frame_index must be consecutive; represent an "
                "analyzed frame with no cages as an empty snapshot so gap_frame is "
                "recorded explicitly."
            )
        if (
            previous.time_ps is not None
            and snapshot.time_ps is not None
            and snapshot.time_ps < previous.time_ps
        ):
            raise ValueError("Tracking snapshots must be sorted by nondecreasing time.")

    def _new_track_id(self) -> str:
        track_id = f"t{self._next_track_number}"
        self._next_track_number += 1
        return track_id

    def _add_first_frame(
        self, stamp: FrameStamp, states: Sequence[_CageState]
    ) -> None:
        for state in states:
            track_id = self._new_track_id()
            self._track_observations[track_id] = [
                _observation(track_id, stamp, state)
            ]
            self._active.append(
                _TrackedState(state=state, track_id=track_id, last_position=0)
            )
            self._events.append(
                _event(
                    "birth",
                    stamp,
                    destination=(track_id,),
                    destination_types=(state.cage_type,),
                    destination_phases=state.phase_labels,
                    censored=True,
                )
            )

    def _expire(self, previous: _TrackedState) -> None:
        if previous.track_id in self._expired:
            return
        first_absent_position = previous.last_position + 1
        stamp = self._frames[first_absent_position]
        self._events.append(
            _event(
                "death",
                stamp,
                source=(previous.track_id,),
                source_types=(previous.state.cage_type,),
                source_phases=previous.state.phase_labels,
            )
        )
        self._expired.add(previous.track_id)

    def _add_change_events(
        self,
        stamp: FrameStamp,
        current_states: Sequence[_CageState],
        current_track_ids: Mapping[int, str],
        assignments: Sequence[_Candidate],
    ) -> None:
        for candidate in assignments:
            previous = self._active[candidate.previous_index]
            current = current_states[candidate.current_index]
            track_id = current_track_ids[candidate.current_index]
            if previous.state.cage_type != current.cage_type:
                self._events.append(
                    _event(
                        "type_change",
                        stamp,
                        source=(track_id,),
                        destination=(track_id,),
                        source_types=(previous.state.cage_type,),
                        destination_types=(current.cage_type,),
                        source_phases=previous.state.phase_labels,
                        destination_phases=current.phase_labels,
                    )
                )
            if previous.state.phase_labels != current.phase_labels:
                self._events.append(
                    _event(
                        "phase_change",
                        stamp,
                        source=(track_id,),
                        destination=(track_id,),
                        source_types=(previous.state.cage_type,),
                        destination_types=(current.cage_type,),
                        source_phases=previous.state.phase_labels,
                        destination_phases=current.phase_labels,
                    )
                )

    def _add_split_merge_events(
        self,
        stamp: FrameStamp,
        current_states: Sequence[_CageState],
        current_track_ids: Mapping[int, str],
        candidates: Sequence[_Candidate],
    ) -> None:
        previous_links: dict[int, set[int]] = defaultdict(set)
        current_links: dict[int, set[int]] = defaultdict(set)
        for candidate in candidates:
            previous_links[candidate.previous_index].add(candidate.current_index)
            current_links[candidate.current_index].add(candidate.previous_index)
        for previous_index in sorted(previous_links):
            current_indexes = sorted(previous_links[previous_index])
            destinations = _unique(
                current_track_ids[index] for index in current_indexes
            )
            if len(destinations) > 1:
                previous = self._active[previous_index]
                self._events.append(
                    _event(
                        "split",
                        stamp,
                        source=(previous.track_id,),
                        destination=destinations,
                        source_types=(previous.state.cage_type,),
                        destination_types=tuple(
                            current_states[index].cage_type
                            for index in current_indexes
                        ),
                        source_phases=previous.state.phase_labels,
                        destination_phases=_unique(
                            phase
                            for index in current_indexes
                            for phase in current_states[index].phase_labels
                        ),
                    )
                )
        for current_index in sorted(current_links):
            previous_indexes = sorted(current_links[current_index])
            sources = _unique(
                self._active[index].track_id for index in previous_indexes
            )
            if len(sources) > 1:
                current = current_states[current_index]
                self._events.append(
                    _event(
                        "merge",
                        stamp,
                        source=sources,
                        destination=(current_track_ids[current_index],),
                        source_types=tuple(
                            self._active[index].state.cage_type
                            for index in previous_indexes
                        ),
                        destination_types=(current.cage_type,),
                        source_phases=_unique(
                            phase
                            for index in previous_indexes
                            for phase in self._active[index].state.phase_labels
                        ),
                        destination_phases=current.phase_labels,
                    )
                )


def track_cages(
    frame_results: Iterable[FrameResult],
    config: TrackingConfig | None = None,
) -> TrackingResult:
    """Reduce and track frame results without retaining the result sequence."""
    return track_snapshots(
        (
            snapshot_from_frame_result(result, frame_index)
            for frame_index, result in enumerate(frame_results)
        ),
        config,
    )


def track_snapshots(
    snapshots: Iterable[TrackFrameSnapshot],
    config: TrackingConfig | None = None,
) -> TrackingResult:
    """Assign persistent ``tN`` IDs while consuming snapshots once."""
    accumulator = TrackingAccumulator(config)
    for snapshot in snapshots:
        accumulator.add(snapshot)
    return accumulator.result()


def parse_targets(value: str | Iterable[str]) -> tuple[TargetSpec, ...]:
    """Parse comma-separated all, cage-type, phase, and ``tN`` targets."""
    raw_parts = value.split(",") if isinstance(value, str) else tuple(value)
    parts = [
        piece.strip()
        for raw in raw_parts
        for piece in str(raw).split(",")
        if piece.strip()
    ]
    if not parts:
        raise ValueError("At least one tracking target is required.")
    targets: list[TargetSpec] = []
    seen: set[tuple[str, str]] = set()
    for raw in parts:
        folded = raw.casefold()
        if folded == "all":
            target = TargetSpec(raw=raw, kind="all", value="all")
        elif match := _TRACK_PATTERN.fullmatch(raw):
            target = TargetSpec(raw=raw, kind="track", value=f"t{int(match.group(1))}")
        elif folded in _PHASE_ALIASES:
            target = TargetSpec(raw=raw, kind="phase", value=_PHASE_ALIASES[folded])
        else:
            compact = _COMPACT_GENERIC_CAGE_PATTERN.fullmatch(raw)
            cage_label = (
                f"4^{compact.group(1)}-5^{compact.group(2)}-6^{compact.group(3)}"
                if compact is not None
                else raw
            )
            try:
                cage_type = canonical_cage_type(cage_label)
            except ValueError as exc:
                raise ValueError(f"Unsupported tracking target: {raw!r}.") from exc
            target = TargetSpec(raw=raw, kind="cage_type", value=cage_type)
        identity = (target.kind, target.value)
        if identity not in seen:
            seen.add(identity)
            targets.append(target)
    return tuple(targets)


def select_targets(
    result: TrackingResult,
    targets: str | Iterable[str],
) -> tuple[TargetSelection, ...]:
    """Select full lifecycles for every requested target."""
    selections: list[TargetSelection] = []
    for target in parse_targets(targets):
        selected = tuple(
            track for track in result.tracks if _track_matches(track, target)
        )
        track_ids = {track.track_id for track in selected}
        events = tuple(
            event for event in result.events if track_ids.intersection(event.track_ids)
        )
        selections.append(
            TargetSelection(
                target=target,
                frames=result.frames,
                tracks=selected,
                events=events,
                config=result.config,
            )
        )
    return tuple(selections)


def observation_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    target = _target_label(data)
    return [
        {
            "target": target,
            "track_id": item.track_id,
            "frame_index": item.frame_index,
            "frame": item.frame_name,
            "time_ps": item.time_ps,
            "cage_id": item.local_cage_id,
            "cage_type": item.cage_type,
            "phase": item.phase,
            "water_atomids": ",".join(map(str, item.water_atomids)),
            "water_count": len(item.water_atomids),
            "guest_ids": ",".join(item.guest_ids),
            "guest_count": len(item.guest_ids),
            "center_x_nm": item.center[0],
            "center_y_nm": item.center[1],
            "center_z_nm": item.center[2],
            "match_jaccard": item.match_jaccard,
            "match_shared_fraction": item.match_shared_fraction,
            "match_center_distance_nm": item.match_center_distance_nm,
            "match_topology_similarity": item.match_topology_similarity,
            "gap_frames": item.gap_frames,
        }
        for item in data.observations
    ]


def lifetime_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    """Return one sample row per persistent cage."""
    target = _target_label(data)
    positions = {frame.frame_index: index for index, frame in enumerate(data.frames)}
    death_by_track = {
        track_id: event
        for event in data.events
        if event.kind == "death"
        for track_id in event.source_track_ids
    }
    rows: list[Row] = []
    for track in data.tracks:
        first, last = track.first, track.last
        observed_span_ps = (
            None
            if first.time_ps is None or last.time_ps is None
            else float(last.time_ps - first.time_ps)
        )
        death = death_by_track.get(track.track_id)
        end_time_ps = (
            death.time_ps
            if death is not None and death.time_ps is not None
            else last.time_ps
        )
        lifetime_ps = (
            None
            if first.time_ps is None or end_time_ps is None
            else float(end_time_ps - first.time_ps)
        )
        first_position = positions.get(first.frame_index)
        last_position = positions.get(last.frame_index)
        span_frames = (
            None
            if first_position is None or last_position is None
            else last_position - first_position + 1
        )
        rows.append(
            {
                "target": target,
                "track_id": track.track_id,
                "first_frame_index": first.frame_index,
                "last_frame_index": last.frame_index,
                "first_time_ps": first.time_ps,
                "last_time_ps": last.time_ps,
                "end_time_ps": end_time_ps,
                "observed_frames": len(track.observations),
                "span_frames": span_frames,
                "gap_frames": track.gap_frames,
                "observed_span_ps": observed_span_ps,
                "lifetime_ps": lifetime_ps,
                "duration_status": _duration_status(track, lifetime_ps),
                "left_censored": track.left_censored,
                "right_censored": track.right_censored,
                "initial_cage_type": first.cage_type,
                "final_cage_type": last.cage_type,
                "cage_types": ",".join(
                    _unique(item.cage_type for item in track.observations)
                ),
                "phases": ",".join(
                    _unique(
                        phase
                        for item in track.observations
                        for phase in item.phase_labels
                    )
                ),
            }
        )
    return rows


def lifetime_distribution_rows(
    data: TrackingResult | TargetSelection,
) -> list[Row]:
    """Aggregate exact lifetime samples into a standalone discrete distribution."""
    samples = lifetime_rows(data)
    if not samples:
        return []
    grouped: dict[tuple[int | None, float | None], list[Row]] = defaultdict(list)
    for sample in samples:
        span = sample.get("span_frames")
        lifetime = sample.get("lifetime_ps")
        key = (
            None if span is None else int(span),
            None if lifetime is None else float(lifetime),
        )
        grouped[key].append(sample)
    ordered = sorted(
        grouped,
        key=lambda item: (
            item[1] is None,
            inf if item[1] is None else item[1],
            inf if item[0] is None else item[0],
        ),
    )
    total = len(samples)
    cumulative = 0
    rows: list[Row] = []
    for span_frames, lifetime_ps in ordered:
        values = grouped[(span_frames, lifetime_ps)]
        count = len(values)
        cumulative += count
        rows.append(
            {
                "target": _target_label(data),
                "lifetime_ps": lifetime_ps,
                "span_frames": span_frames,
                "track_count": count,
                "fraction": count / total,
                "cumulative_fraction": cumulative / total,
                "survival_fraction": sum(
                    len(grouped[key]) for key in ordered
                    if _lifetime_key_ge(key, (span_frames, lifetime_ps))
                )
                / total,
                "uncensored_count": sum(
                    not bool(item["left_censored"])
                    and not bool(item["right_censored"])
                    for item in values
                ),
                "left_censored_count": sum(
                    bool(item["left_censored"]) for item in values
                ),
                "right_censored_count": sum(
                    bool(item["right_censored"]) for item in values
                ),
            }
        )
    return rows


def population_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    target = _target_label(data)
    by_frame: dict[int, list[CageObservation]] = defaultdict(list)
    for item in data.observations:
        by_frame[item.frame_index].append(item)
    rows: list[Row] = []
    for stamp in data.frames:
        observations = by_frame.get(stamp.frame_index, [])
        rows.append(_population_row(target, stamp, "total", "all", len(observations)))
        for cage_type, count in sorted(Counter(item.cage_type for item in observations).items()):
            rows.append(_population_row(target, stamp, "cage_type", cage_type, count))
        phase_counts = Counter(
            phase for item in observations for phase in item.phase_labels
        )
        for phase in sorted(phase_counts, key=_phase_sort_key):
            rows.append(
                _population_row(target, stamp, "phase", phase, phase_counts[phase])
            )
    return rows


def event_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    target = _target_label(data)
    return [
        {
            "target": target,
            "event_id": item.event_id,
            "event": item.kind,
            "frame_index": item.frame_index,
            "frame": item.frame_name,
            "time_ps": item.time_ps,
            "source_track_ids": ",".join(item.source_track_ids),
            "destination_track_ids": ",".join(item.destination_track_ids),
            "source_cage_types": ",".join(item.source_cage_types),
            "destination_cage_types": ",".join(item.destination_cage_types),
            "source_phases": ",".join(item.source_phases),
            "destination_phases": ",".join(item.destination_phases),
            "gap_frames": item.gap_frames,
            "censored": item.censored,
        }
        for item in data.events
    ]


def guest_residence_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    """Return contiguous guest episodes; a cage gap always splits an episode."""
    target = _target_label(data)
    frame_positions = {
        stamp.frame_index: position for position, stamp in enumerate(data.frames)
    }
    rows: list[Row] = []
    for track in data.tracks:
        by_guest: dict[str, list[CageObservation]] = defaultdict(list)
        for item in track.observations:
            for guest_id in item.guest_ids:
                by_guest[guest_id].append(item)
        for guest_id in sorted(by_guest):
            for episode_index, episode in enumerate(
                _contiguous_episodes(by_guest[guest_id], frame_positions), start=1
            ):
                first, last = episode[0], episode[-1]
                residence = (
                    None
                    if first.time_ps is None or last.time_ps is None
                    else float(last.time_ps - first.time_ps)
                )
                rows.append(
                    {
                        "target": target,
                        "track_id": track.track_id,
                        "guest_id": guest_id,
                        "episode": episode_index,
                        "start_frame_index": first.frame_index,
                        "end_frame_index": last.frame_index,
                        "start_time_ps": first.time_ps,
                        "end_time_ps": last.time_ps,
                        "observed_frames": len(episode),
                        "residence_time_ps": residence,
                        "left_censored": (
                            track.left_censored
                            and first.frame_index == track.first.frame_index
                        ),
                        "right_censored": (
                            track.right_censored
                            and last.frame_index == track.last.frame_index
                        ),
                    }
                )
    return rows


def snapshot_from_frame_result(
    result: FrameResult,
    frame_index: int,
) -> TrackFrameSnapshot:
    """Reduce one 0.4.3 ``FrameResult`` to stable identity metadata."""
    states = _frame_states(result)
    raw_box = result.frame.box
    box = (
        None
        if raw_box is None
        else tuple(float(value) for value in np.asarray(raw_box, dtype=float).ravel())
    )
    return TrackFrameSnapshot(
        frame_index=int(frame_index),
        frame_name=result.frame.name,
        time_ps=None if result.frame.time_ps is None else float(result.frame.time_ps),
        source="" if result.frame.source is None else str(result.frame.source),
        box=box,
        cages=tuple(
            TrackCageSnapshot(
                local_cage_id=state.local_cage_id,
                cage_type=state.cage_type,
                phase_labels=state.phase_labels,
                water_atomids=tuple(sorted(state.water_atomids)),
                center=state.center,
                topology=state.topology,
                guest_ids=state.guest_ids,
            )
            for state in states
        ),
    )


def _snapshot_stamp(snapshot: TrackFrameSnapshot) -> FrameStamp:
    return FrameStamp(
        frame_index=snapshot.frame_index,
        frame_name=snapshot.frame_name,
        time_ps=snapshot.time_ps,
        source=snapshot.source,
    )


def _snapshot_states(snapshot: TrackFrameSnapshot) -> list[_CageState]:
    return sorted(
        (
            _CageState(
                local_cage_id=cage.local_cage_id,
                cage_type=cage.cage_type,
                phase=_phase_label(cage.phase_labels),
                phase_labels=cage.phase_labels,
                water_atomids=frozenset(cage.water_atomids),
                center=cage.center,
                topology=cage.topology,
                guest_ids=cage.guest_ids,
            )
            for cage in snapshot.cages
        ),
        key=_state_sort_key,
    )


def _frame_states(result: FrameResult) -> list[_CageState]:
    atom_indexes = [atom.index for atom in result.frame.atoms]
    if len(atom_indexes) != len(set(atom_indexes)):
        raise ValueError(f"Atom indexes are not unique in frame {result.frame.name}.")
    oxygen_indexes = [water.oxygen for water in result.waters]
    if len(oxygen_indexes) != len(set(oxygen_indexes)):
        raise ValueError(
            f"Water oxygen indexes are not unique in frame {result.frame.name}."
        )
    atom_by_index = {int(atom.index): atom for atom in result.frame.atoms}
    oxygen_to_identity: dict[int, int] = {}
    for water in result.waters:
        oxygen_index = int(water.oxygen)
        if oxygen_index not in atom_by_index:
            raise ValueError(
                f"Water oxygen index {water.oxygen} is absent from frame "
                f"{result.frame.name}."
            )
        oxygen_to_identity[oxygen_index] = oxygen_index + 1
    if len(set(oxygen_to_identity.values())) != len(oxygen_to_identity):
        raise ValueError(
            f"Water topology identities are not unique in frame {result.frame.name}."
        )
    ring_by_id = {
        ring.object_id: ring for rings in result.rings.values() for ring in rings
    }
    phase_by_cage = _phase_labels_by_cage(result)
    cages = result.all_cages or result.cages
    states = [
        _cage_state(cage, oxygen_to_identity, ring_by_id, phase_by_cage)
        for cage in cages
    ]
    local_ids = [state.local_cage_id for state in states]
    if len(local_ids) != len(set(local_ids)):
        raise ValueError(f"Cage IDs are not unique in frame {result.frame.name}.")
    return sorted(states, key=_state_sort_key)


def _cage_state(
    cage: Cage,
    oxygen_to_identity: Mapping[int, int],
    ring_by_id: Mapping[str, Ring],
    phase_by_cage: Mapping[str, tuple[str, ...]],
) -> _CageState:
    missing = sorted(set(cage.waters).difference(oxygen_to_identity))
    if missing:
        raise ValueError(
            f"Cage {cage.object_id} contains unknown water oxygen indexes: {missing}."
        )
    water_atomids = frozenset(oxygen_to_identity[index] for index in cage.waters)
    if len(water_atomids) != len(cage.waters):
        raise ValueError(f"Cage {cage.object_id} repeats a water topology identity.")
    missing_rings = sorted(set(cage.rings).difference(ring_by_id))
    if missing_rings:
        raise ValueError(
            f"Cage {cage.object_id} references unknown rings: {missing_rings}."
        )
    center = np.asarray(cage.center, dtype=float)
    if center.shape != (3,) or np.any(~np.isfinite(center)):
        raise ValueError(f"Cage {cage.object_id} has an invalid center.")
    phases = phase_by_cage.get(cage.object_id, ("unassigned",))
    return _CageState(
        local_cage_id=cage.object_id,
        cage_type=cage.cage_type,
        phase=_phase_label(phases),
        phase_labels=phases,
        water_atomids=water_atomids,
        center=tuple(float(value) for value in center),
        topology=tuple(sorted(ring_by_id[ring_id].size for ring_id in cage.rings)),
        guest_ids=tuple(sorted(set(cage.guest_ids))),
    )


def _phase_labels_by_cage(result: FrameResult) -> dict[str, tuple[str, ...]]:
    domain_phases: dict[str, set[str]] = defaultdict(set)
    fallback: dict[str, set[str]] = defaultdict(set)
    for domain in result.hydrate_domains:
        phase = _canonical_phase(domain.hydrate_type)
        for cage_id in domain.cage_ids:
            domain_phases[cage_id].add(phase)
    for cluster in result.hydrate_clusters:
        cluster_phase = _canonical_phase(cluster.hydrate_type)
        if cluster_phase in {"sI", "sII", "sH"}:
            for cage_id in cluster.classified_cage_ids:
                fallback[cage_id].add(cluster_phase)
        for cage_id in cluster.boundary_cage_ids:
            fallback[cage_id].add("boundary")
        for cage_id in cluster.ambiguous_cage_ids:
            fallback[cage_id].add("ambiguous")
        for cage_id in cluster.unclassified_cage_ids:
            fallback[cage_id].add("unclassified")
    for cage_id in result.isolated_cage_ids:
        fallback[cage_id].add("isolated")
    labels: dict[str, tuple[str, ...]] = {}
    for cage_id in set(domain_phases) | set(fallback):
        labels[cage_id] = _normalized_phases(
            domain_phases.get(cage_id) or fallback[cage_id]
        )
    return labels


def _match_candidates(
    previous: Sequence[_TrackedState],
    current: Sequence[_CageState],
    box: np.ndarray | None,
    config: TrackingConfig,
    current_position: int,
) -> list[_Candidate]:
    current_by_water: dict[int, list[int]] = defaultdict(list)
    for current_index, state in enumerate(current):
        for water_id in state.water_atomids:
            current_by_water[water_id].append(current_index)
    candidates: list[_Candidate] = []
    for previous_index, old in enumerate(previous):
        gap_frames = current_position - old.last_position - 1
        if gap_frames > config.gap_frame:
            continue
        shared_counts: Counter[int] = Counter()
        for water_id in old.state.water_atomids:
            shared_counts.update(current_by_water.get(water_id, ()))
        for current_index in sorted(shared_counts):
            shared = shared_counts[current_index]
            if shared < config.min_shared_waters:
                continue
            new = current[current_index]
            union = len(old.state.water_atomids | new.water_atomids)
            smaller = min(len(old.state.water_atomids), len(new.water_atomids))
            jaccard = shared / union if union else 0.0
            shared_fraction = shared / smaller if smaller else 0.0
            if (
                jaccard < config.min_jaccard
                or shared_fraction < config.min_shared_fraction
            ):
                continue
            center_distance = distance(
                np.asarray(old.state.center), np.asarray(new.center), box
            )
            if (
                config.max_center_distance_nm is not None
                and center_distance > config.max_center_distance_nm
            ):
                continue
            candidates.append(
                _Candidate(
                    previous_index=previous_index,
                    current_index=current_index,
                    jaccard=jaccard,
                    shared_fraction=shared_fraction,
                    center_distance_nm=center_distance,
                    topology_similarity=_multiset_similarity(
                        old.state.topology, new.topology
                    ),
                    guest_similarity=_set_similarity(
                        old.state.guest_ids, new.guest_ids
                    ),
                    gap_frames=gap_frames,
                )
            )
    return candidates


def _maximum_weight_assign(
    candidates: Sequence[_Candidate],
    previous: Sequence[_TrackedState],
    current: Sequence[_CageState],
    config: TrackingConfig,
) -> list[_Candidate]:
    if not candidates or not previous or not current:
        return []
    candidate_by_pair = {
        (candidate.previous_index, candidate.current_index): candidate
        for candidate in candidates
    }
    row_count = len(previous)
    real_column_count = len(current)
    column_count = real_column_count + row_count
    match_bonus = float((max(row_count, real_column_count) + 1) * 1000)
    blocked_cost = match_bonus * 2.0
    costs = [[0.0] * column_count for _ in range(row_count)]
    for row in range(row_count):
        for column in range(real_column_count):
            candidate = candidate_by_pair.get((row, column))
            if candidate is None:
                costs[row][column] = blocked_cost
                continue
            center_score = (
                0.0
                if candidate.center_distance_nm is None
                else 1.0 / (1.0 + candidate.center_distance_nm)
            )
            quality = (
                100.0 * candidate.jaccard
                + 10.0 * candidate.shared_fraction
                + center_score
                + candidate.topology_similarity
                + (
                    0.001 * candidate.guest_similarity
                    if config.guest_tiebreak
                    else 0.0
                )
                - 0.01 * candidate.gap_frames
            )
            costs[row][column] = -(match_bonus + quality)
    columns = _hungarian_columns(costs)
    return sorted(
        (
            candidate_by_pair[(row, column)]
            for row, column in enumerate(columns)
            if column < real_column_count and (row, column) in candidate_by_pair
        ),
        key=lambda item: item.current_index,
    )


def _hungarian_columns(costs: Sequence[Sequence[float]]) -> list[int]:
    row_count = len(costs)
    column_count = len(costs[0]) if row_count else 0
    if row_count > column_count:
        raise ValueError("Hungarian assignment requires at least as many columns as rows.")
    u = [0.0] * (row_count + 1)
    v = [0.0] * (column_count + 1)
    matched_row = [0] * (column_count + 1)
    previous_column = [0] * (column_count + 1)
    tolerance = 1.0e-12
    for row in range(1, row_count + 1):
        matched_row[0] = row
        minimum = [inf] * (column_count + 1)
        used = [False] * (column_count + 1)
        column = 0
        while True:
            used[column] = True
            active_row = matched_row[column]
            delta = inf
            next_column = 0
            for candidate_column in range(1, column_count + 1):
                if used[candidate_column]:
                    continue
                reduced = (
                    costs[active_row - 1][candidate_column - 1]
                    - u[active_row]
                    - v[candidate_column]
                )
                if reduced < minimum[candidate_column] - tolerance:
                    minimum[candidate_column] = reduced
                    previous_column[candidate_column] = column
                if minimum[candidate_column] < delta - tolerance:
                    delta = minimum[candidate_column]
                    next_column = candidate_column
            for candidate_column in range(column_count + 1):
                if used[candidate_column]:
                    u[matched_row[candidate_column]] += delta
                    v[candidate_column] -= delta
                else:
                    minimum[candidate_column] -= delta
            column = next_column
            if matched_row[column] == 0:
                break
        while True:
            next_column = previous_column[column]
            matched_row[column] = matched_row[next_column]
            column = next_column
            if column == 0:
                break
    assignment = [-1] * row_count
    for column in range(1, column_count + 1):
        if matched_row[column] != 0:
            assignment[matched_row[column] - 1] = column - 1
    return assignment


def _observation(
    track_id: str,
    stamp: FrameStamp,
    state: _CageState,
    evidence: _Candidate | None = None,
) -> CageObservation:
    return CageObservation(
        track_id=track_id,
        frame_index=stamp.frame_index,
        frame_name=stamp.frame_name,
        time_ps=stamp.time_ps,
        local_cage_id=state.local_cage_id,
        cage_type=state.cage_type,
        phase=state.phase,
        phase_labels=state.phase_labels,
        water_atomids=tuple(sorted(state.water_atomids)),
        center=state.center,
        topology=state.topology,
        guest_ids=state.guest_ids,
        match_jaccard=None if evidence is None else evidence.jaccard,
        match_shared_fraction=None if evidence is None else evidence.shared_fraction,
        match_center_distance_nm=(
            None if evidence is None else evidence.center_distance_nm
        ),
        match_topology_similarity=(
            None if evidence is None else evidence.topology_similarity
        ),
        gap_frames=0 if evidence is None else evidence.gap_frames,
    )


def _event(
    kind: EventKind,
    stamp: FrameStamp,
    source: Iterable[str] = (),
    destination: Iterable[str] = (),
    source_types: Iterable[str] = (),
    destination_types: Iterable[str] = (),
    source_phases: Iterable[str] = (),
    destination_phases: Iterable[str] = (),
    gap_frames: int = 0,
    censored: bool = False,
) -> TrackEvent:
    return TrackEvent(
        event_id="",
        kind=kind,
        frame_index=stamp.frame_index,
        frame_name=stamp.frame_name,
        time_ps=stamp.time_ps,
        source_track_ids=_unique(source),
        destination_track_ids=_unique(destination),
        source_cage_types=tuple(source_types),
        destination_cage_types=tuple(destination_types),
        source_phases=_unique(source_phases),
        destination_phases=_unique(destination_phases),
        gap_frames=int(gap_frames),
        censored=bool(censored),
    )


def _track_matches(track: CageTrack, target: TargetSpec) -> bool:
    if target.kind == "all":
        return True
    if target.kind == "track":
        return track.track_id == target.value
    if target.kind == "cage_type":
        return any(item.cage_type == target.value for item in track.observations)
    if target.value == "mixed":
        return any(len(item.phase_labels) > 1 for item in track.observations)
    return any(target.value in item.phase_labels for item in track.observations)


def _population_row(
    target: str,
    stamp: FrameStamp,
    group: str,
    label: str,
    count: int,
) -> Row:
    return {
        "target": target,
        "frame_index": stamp.frame_index,
        "frame": stamp.frame_name,
        "time_ps": stamp.time_ps,
        "group": group,
        "label": label,
        "cage_count": count,
    }


def _contiguous_episodes(
    observations: Sequence[CageObservation],
    frame_positions: Mapping[int, int],
) -> list[list[CageObservation]]:
    ordered = sorted(observations, key=_observation_sort_key)
    episodes: list[list[CageObservation]] = []
    for item in ordered:
        previous_position = (
            None if not episodes else frame_positions.get(episodes[-1][-1].frame_index)
        )
        current_position = frame_positions.get(item.frame_index)
        if (
            not episodes
            or previous_position is None
            or current_position != previous_position + 1
        ):
            episodes.append([item])
        else:
            episodes[-1].append(item)
    return episodes


def _target_label(data: TrackingResult | TargetSelection) -> str:
    return "all" if isinstance(data, TrackingResult) else data.target.value


def _state_sort_key(state: _CageState) -> tuple[object, ...]:
    return (
        min(state.water_atomids, default=-1),
        tuple(sorted(state.water_atomids)),
        state.cage_type,
        state.local_cage_id,
    )


def _observation_sort_key(item: CageObservation) -> tuple[int, int, str]:
    return item.frame_index, _track_number(item.track_id), item.local_cage_id


def _track_number(track_id: str) -> int:
    match = _TRACK_PATTERN.fullmatch(track_id)
    if match is None:
        raise ValueError(f"Invalid persistent track ID: {track_id!r}.")
    return int(match.group(1))


def _canonical_phase(value: str) -> str:
    text = str(value).strip()
    return _PHASE_ALIASES.get(text.casefold(), text or "unassigned")


def _normalized_phases(values: Iterable[str]) -> tuple[str, ...]:
    phases = {_canonical_phase(value) for value in values if str(value).strip()}
    if not phases:
        phases.add("unassigned")
    return tuple(sorted(phases, key=_phase_sort_key))


def _phase_label(values: Sequence[str]) -> str:
    unique = _normalized_phases(values)
    return unique[0] if len(unique) == 1 else "mixed(" + "+".join(unique) + ")"


def _phase_sort_key(value: str) -> tuple[int, str]:
    return _PHASE_ORDER.get(value, 100), value


def _multiset_similarity(left: Sequence[int], right: Sequence[int]) -> float:
    if not left and not right:
        return 0.0
    left_counts = Counter(left)
    right_counts = Counter(right)
    intersection = sum((left_counts & right_counts).values())
    union = sum((left_counts | right_counts).values())
    return intersection / union if union else 0.0


def _set_similarity(left: Sequence[str], right: Sequence[str]) -> float:
    left_set, right_set = set(left), set(right)
    if not left_set and not right_set:
        return 0.0
    return len(left_set & right_set) / len(left_set | right_set)


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values))



def _lifetime_key_ge(
    left: tuple[int | None, float | None],
    right: tuple[int | None, float | None],
) -> bool:
    left_frames, left_ps = left
    right_frames, right_ps = right
    if left_ps is not None and right_ps is not None:
        return left_ps >= right_ps
    if left_frames is not None and right_frames is not None:
        return left_frames >= right_frames
    return left_ps is None and right_ps is not None


def _duration_status(track: CageTrack, lifetime_ps: float | None) -> str:
    if lifetime_ps is None:
        return "time_unavailable"
    if track.left_censored and track.right_censored:
        return "left_and_right_censored_lower_bound"
    if track.left_censored:
        return "left_censored_lower_bound"
    if track.right_censored:
        return "right_censored_lower_bound"
    return "complete"
