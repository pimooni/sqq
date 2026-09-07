"""Deterministic streaming cage tracking and frame-to-frame matching."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from math import inf
import re
from typing import Iterable, Mapping, Sequence

import numpy as np

from ...models import FrameResult
from ...models.tracking import (
    CageObservation,
    CageTrack,
    EventKind,
    FrameStamp,
    TrackEvent,
    TrackFrameSnapshot,
    TrackingConfig,
    TrackingResult,
)
from ..common.pbc import distance
from .snapshot import (
    _CageState,
    _snapshot_stamp,
    _snapshot_states,
    snapshot_from_frame_result,
)

_TRACK_PATTERN = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)
_ASSIGNMENT_TOLERANCE = 1.0e-12

__all__ = ["TrackingAccumulator", "track_cages", "track_snapshots"]


@dataclass(frozen=True)
class _TrackedState:
    state: _CageState
    track_id: str
    last_position: int
    last_time_ps: float | None

@dataclass(frozen=True)
class _Candidate:
    previous_index: int
    current_index: int
    jaccard: float
    shared_fraction: float
    shared_waters: int
    center_distance_nm: float | None
    topology_similarity: float
    guest_similarity: float
    gap_frames: int
    gap_time_ps: float | None

@dataclass(frozen=True)
class _PendingLineage:
    kind: str
    source_track_ids: tuple[str, ...]
    destination_track_ids: tuple[str, ...]
    source_cage_types: tuple[str, ...]
    destination_cage_types: tuple[str, ...]
    source_phases: tuple[str, ...]
    destination_phases: tuple[str, ...]
    water_conservation: float
    center_distance_nm: float | None
    evidence_status: str
    source_states: tuple[_CageState, ...]
    destination_states: tuple[_CageState, ...]

class TrackingAccumulator:
    """Incrementally track snapshots without retaining the snapshot sequence."""

    def __init__(self, config: TrackingConfig | None = None) -> None:
        self.config = config or TrackingConfig()
        self._frames: list[FrameStamp] = []
        self._track_observations: dict[str, list[CageObservation]] = {}
        self._events: list[TrackEvent] = []
        self._active: list[_TrackedState] = []
        self._expired: set[str] = set()
        self._pending_lineage: list[_PendingLineage] = []
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
            stamp.time_ps,
        )
        assignments = _maximum_weight_assign(
            candidates, self._active, current_states, self.config
        )
        evidence_by_current = {
            candidate.current_index: candidate for candidate in assignments
        }
        candidates_by_current: dict[int, list[_Candidate]] = defaultdict(list)
        candidates_by_previous: dict[int, list[_Candidate]] = defaultdict(list)
        for candidate in candidates:
            candidates_by_current[candidate.current_index].append(candidate)
            candidates_by_previous[candidate.previous_index].append(candidate)
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
            diagnostic_candidates = candidates_by_current.get(current_index, ())
            if evidence is not None:
                diagnostic_candidates = tuple(
                    {
                        (item.previous_index, item.current_index): item
                        for item in (
                            *diagnostic_candidates,
                            *candidates_by_previous.get(evidence.previous_index, ()),
                        )
                    }.values()
                )
            self._track_observations[track_id].append(
                _observation(
                    track_id,
                    stamp,
                    state,
                    evidence=evidence,
                    competitors=diagnostic_candidates,
                    config=self.config,
                )
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
                        gap_time_ps=evidence.gap_time_ps,
                    )
                )

        self._confirm_pending_lineage(stamp, current_states, current_track_ids, box)
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
                last_time_ps=stamp.time_ps,
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
        # A dormant cage expires only after its gap allowance, so its death is
        # appended later than events of the frames in between although it is
        # stamped at the first absent frame.  The stable sort restores time
        # order while preserving the within-frame emission order.
        events = tuple(
            replace(event, event_id=f"e{index}")
            for index, event in enumerate(
                sorted(self._events, key=lambda item: item.frame_index), start=1
            )
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
        if self.config.max_gap_ps is not None and (
            previous.time_ps is None or snapshot.time_ps is None
        ):
            raise ValueError(
                "track.max_gap_ps requires physical time for every selected frame."
            )
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
                _TrackedState(
                    state=state,
                    track_id=track_id,
                    last_position=0,
                    last_time_ps=stamp.time_ps,
                )
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
            boundary = {
                "gap_frames": candidate.gap_frames,
                "gap_time_ps": candidate.gap_time_ps,
                "censored": bool(candidate.gap_frames),
                "evidence_status": "gap_unresolved" if candidate.gap_frames else "",
            }
            if previous.state.cage_type != current.cage_type:
                self._events.append(
                    _event(
                        "type_change_unresolved" if candidate.gap_frames else "type_change",
                        stamp,
                        source=(track_id,),
                        destination=(track_id,),
                        source_types=(previous.state.cage_type,),
                        destination_types=(current.cage_type,),
                        source_phases=previous.state.phase_labels,
                        destination_phases=current.phase_labels,
                        **boundary,
                    )
                )
            if previous.state.phase_labels != current.phase_labels:
                self._events.append(
                    _event(
                        "phase_change_unresolved" if candidate.gap_frames else "phase_change",
                        stamp,
                        source=(track_id,),
                        destination=(track_id,),
                        source_types=(previous.state.cage_type,),
                        destination_types=(current.cage_type,),
                        source_phases=previous.state.phase_labels,
                        destination_phases=current.phase_labels,
                        **boundary,
                    )
                )
            self._add_guest_events(stamp, previous, current, track_id, candidate)

    def _add_guest_events(
        self,
        stamp: FrameStamp,
        previous: _TrackedState,
        current: _CageState,
        track_id: str,
        candidate: _Candidate,
    ) -> None:
        before = tuple(sorted(previous.state.guest_ids))
        after = tuple(sorted(current.guest_ids))
        if before == after:
            return
        entered = tuple(sorted(set(after).difference(before)))
        exited = tuple(sorted(set(before).difference(after)))
        common = {
            "source": (track_id,),
            "destination": (track_id,),
            "source_types": (previous.state.cage_type,),
            "destination_types": (current.cage_type,),
            "source_phases": previous.state.phase_labels,
            "destination_phases": current.phase_labels,
            "gap_frames": candidate.gap_frames,
            "gap_time_ps": candidate.gap_time_ps,
            "guest_ids_entered": entered,
            "guest_ids_exited": exited,
            "source_occupancy": _occupancy_class(before),
            "destination_occupancy": _occupancy_class(after),
            "guest_composition_before": before,
            "guest_composition_after": after,
        }
        if candidate.gap_frames:
            self._events.append(
                _event(
                    "unresolved_across_gap",
                    stamp,
                    censored=True,
                    evidence_status="gap_unresolved",
                    **common,
                )
            )
            return
        if entered and exited:
            guest_kind: EventKind = "guest_exchange"
        elif entered:
            guest_kind = "guest_enter"
        else:
            guest_kind = "guest_exit"
        self._events.append(_event(guest_kind, stamp, **common))
        before_class = _occupancy_class(before)
        after_class = _occupancy_class(after)
        if before_class != after_class:
            self._events.append(_event("occupancy_change", stamp, **common))
            detailed_kind = _occupancy_transition_kind(before_class, after_class)
            if detailed_kind is not None:
                self._events.append(_event(detailed_kind, stamp, **common))

    def _confirm_pending_lineage(
        self,
        stamp: FrameStamp,
        current_states: Sequence[_CageState],
        current_track_ids: Mapping[int, str],
        box: np.ndarray | None,
    ) -> None:
        pending, self._pending_lineage = self._pending_lineage, []
        states_by_track = {
            current_track_ids[index]: state
            for index, state in enumerate(current_states)
        }
        for item in pending:
            if any(track_id not in states_by_track for track_id in item.destination_track_ids):
                continue
            if item.kind == "merge" and any(
                track_id in states_by_track
                for track_id in item.source_track_ids
                if track_id not in item.destination_track_ids
            ):
                continue
            destinations = tuple(states_by_track[track_id] for track_id in item.destination_track_ids)
            # A persistent ID alone is not lineage evidence: recheck the actual
            # next-frame water shell, topology and (optionally) center distance.
            continuation = [
                _lineage_continuation(before, after, box, self.config)
                for before, after in zip(item.destination_states, destinations)
            ]
            if not all(continuation):
                continue
            if item.kind == "split":
                conservation = _overlap_fraction(
                    item.source_states[0].water_atomids,
                    set().union(*(state.water_atomids for state in destinations)),
                )
                original = _independent_contribution_sets(
                    item.source_states[0].water_atomids, item.destination_states
                )
                retained = _independent_contribution_sets(
                    item.source_states[0].water_atomids, destinations
                )
                contributions = tuple(len(before & after) for before, after in zip(original, retained))
                kind: EventKind = "split_confirmed"
            else:
                conservation = _overlap_fraction(
                    destinations[0].water_atomids,
                    set().union(*(state.water_atomids for state in item.source_states)),
                )
                contributions = _independent_contributions(
                    destinations[0].water_atomids, item.source_states
                )
                kind = "merge_confirmed"
            if (
                conservation < self.config.min_shared_fraction
                or min(contributions, default=0) < self.config.min_shared_waters
            ):
                continue
            self._events.append(
                _event(
                    kind,
                    stamp,
                    source=item.source_track_ids,
                    destination=item.destination_track_ids,
                    source_types=item.source_cage_types,
                    destination_types=item.destination_cage_types,
                    source_phases=item.source_phases,
                    destination_phases=item.destination_phases,
                    evidence_status="confirmed_water_persistence",
                    water_conservation=min(item.water_conservation, conservation),
                    center_distance_nm=item.center_distance_nm,
                    persistence_frames=1,
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
        previous_ids = {
            state.track_id for state in self._active
            if state.last_position == len(self._frames) - 2
        }
        current_ids = set(current_track_ids.values())
        for candidate in candidates:
            if candidate.gap_frames:
                # Missing observations cannot establish when a lineage event occurred.
                continue
            previous_links[candidate.previous_index].add(candidate.current_index)
            current_links[candidate.current_index].add(candidate.previous_index)
        for previous_index in sorted(previous_links):
            current_indexes = sorted(previous_links[previous_index])
            destinations = _unique(
                current_track_ids[index] for index in current_indexes
            )
            if len(destinations) > 1:
                previous = self._active[previous_index]
                relevant = [
                    candidate
                    for candidate in candidates
                    if candidate.previous_index == previous_index
                    and candidate.current_index in current_indexes
                ]
                destination_water = set().union(
                    *(current_states[index].water_atomids for index in current_indexes)
                )
                conservation = _overlap_fraction(
                    previous.state.water_atomids, destination_water
                )
                center_distance_nm = _maximum_candidate_distance(relevant)
                evidence_status = _lineage_evidence_status(
                    conservation, relevant, self.config
                )
                if all(track_id in previous_ids for track_id in destinations):
                    evidence_status = "candidate_stable_neighbors"
                elif min(_independent_contributions(
                    previous.state.water_atomids,
                    tuple(current_states[index] for index in current_indexes),
                ), default=0) < self.config.min_shared_waters:
                    evidence_status = "candidate_shared_water_only"
                self._events.append(
                    _event(
                        "split_candidate",
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
                        evidence_status=evidence_status,
                        water_conservation=conservation,
                        center_distance_nm=center_distance_nm,
                    )
                )
                if evidence_status == "candidate_supported":
                    self._pending_lineage.append(
                        _PendingLineage(
                            kind="split",
                            source_track_ids=(previous.track_id,),
                            destination_track_ids=destinations,
                            source_cage_types=(previous.state.cage_type,),
                            destination_cage_types=tuple(
                                current_states[index].cage_type
                                for index in current_indexes
                            ),
                            source_phases=previous.state.phase_labels,
                            destination_phases=_unique(
                                phase
                                for index in current_indexes
                                for phase in current_states[index].phase_labels
                            ),
                            water_conservation=conservation,
                            center_distance_nm=center_distance_nm,
                            evidence_status=evidence_status,
                            source_states=(previous.state,),
                            destination_states=tuple(current_states[index] for index in current_indexes),
                        )
                    )
        for current_index in sorted(current_links):
            previous_indexes = sorted(current_links[current_index])
            sources = _unique(
                self._active[index].track_id for index in previous_indexes
            )
            if len(sources) > 1:
                current = current_states[current_index]
                relevant = [
                    candidate
                    for candidate in candidates
                    if candidate.current_index == current_index
                    and candidate.previous_index in previous_indexes
                ]
                source_water = set().union(
                    *(self._active[index].state.water_atomids for index in previous_indexes)
                )
                conservation = _overlap_fraction(current.water_atomids, source_water)
                center_distance_nm = _maximum_candidate_distance(relevant)
                evidence_status = _lineage_evidence_status(
                    conservation, relevant, self.config
                )
                if all(track_id in current_ids for track_id in sources):
                    evidence_status = "candidate_stable_neighbors"
                elif min(_independent_contributions(
                    current.water_atomids,
                    tuple(self._active[index].state for index in previous_indexes),
                ), default=0) < self.config.min_shared_waters:
                    evidence_status = "candidate_shared_water_only"
                self._events.append(
                    _event(
                        "merge_candidate",
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
                        evidence_status=evidence_status,
                        water_conservation=conservation,
                        center_distance_nm=center_distance_nm,
                    )
                )
                if evidence_status == "candidate_supported":
                    self._pending_lineage.append(
                        _PendingLineage(
                            kind="merge",
                            source_track_ids=sources,
                            destination_track_ids=(current_track_ids[current_index],),
                            source_cage_types=tuple(
                                self._active[index].state.cage_type
                                for index in previous_indexes
                            ),
                            destination_cage_types=(current.cage_type,),
                            source_phases=_unique(
                                phase
                                for index in previous_indexes
                                for phase in self._active[index].state.phase_labels
                            ),
                            destination_phases=current.phase_labels,
                            water_conservation=conservation,
                            center_distance_nm=center_distance_nm,
                            evidence_status=evidence_status,
                            source_states=tuple(self._active[index].state for index in previous_indexes),
                            destination_states=(current,),
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

def _match_candidates(
    previous: Sequence[_TrackedState],
    current: Sequence[_CageState],
    box: np.ndarray | None,
    config: TrackingConfig,
    current_position: int,
    current_time_ps: float | None,
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
        elapsed_ps = (
            None
            if current_time_ps is None or old.last_time_ps is None
            else float(current_time_ps - old.last_time_ps)
        )
        gap_time_ps = elapsed_ps if gap_frames else 0.0 if elapsed_ps is not None else None
        if (
            gap_frames
            and config.max_gap_ps is not None
            and elapsed_ps is not None
            and elapsed_ps > config.max_gap_ps
        ):
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
                    shared_waters=shared,
                    center_distance_nm=center_distance,
                    topology_similarity=_multiset_similarity(
                        old.state.topology, new.topology
                    ),
                    guest_similarity=_set_similarity(
                        old.state.guest_ids, new.guest_ids
                    ),
                    gap_frames=gap_frames,
                    gap_time_ps=gap_time_ps,
                )
            )
    return candidates

def _maximum_weight_assign(
    candidates: Sequence[_Candidate],
    previous: Sequence[_TrackedState],
    current: Sequence[_CageState],
    config: TrackingConfig,
) -> list[_Candidate]:
    """Return the global maximum-weight one-to-one continuation set.

    Connected candidate components are independent.  Zero-cost dummy columns
    preserve unmatched choices, and sorted rows and columns preserve ties.
    """
    if not candidates or not previous or not current:
        return []
    candidate_by_pair = {
        (candidate.previous_index, candidate.current_index): candidate
        for candidate in candidates
    }
    # One frame-wide bonus makes match count dominate match quality.
    match_bonus = float((max(len(previous), len(current)) + 1) * 1000)
    blocked_cost = match_bonus * 2.0
    assigned: list[_Candidate] = []
    for rows, columns in _candidate_components(candidates):
        if len(rows) == 1 and len(columns) == 1:
            candidate = candidate_by_pair[(rows[0], columns[0])]
            candidate_cost = -(
                match_bonus + _candidate_quality(candidate, config)
            )
            # Retain the dummy-column choice when a long gap favors no match.
            if candidate_cost <= _ASSIGNMENT_TOLERANCE:
                assigned.append(candidate)
            continue
        real_column_count = len(columns)
        costs = [[0.0] * (real_column_count + len(rows)) for _ in rows]
        for local_row, row in enumerate(rows):
            for local_column, column in enumerate(columns):
                candidate = candidate_by_pair.get((row, column))
                if candidate is None:
                    costs[local_row][local_column] = blocked_cost
                    continue
                quality = _candidate_quality(candidate, config)
                costs[local_row][local_column] = -(match_bonus + quality)
        for local_row, local_column in enumerate(_hungarian_columns(costs)):
            if local_column >= real_column_count:
                continue
            pair = (rows[local_row], columns[local_column])
            if pair in candidate_by_pair:
                assigned.append(candidate_by_pair[pair])
    return sorted(assigned, key=lambda item: item.current_index)


def _candidate_components(
    candidates: Sequence[_Candidate],
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Group candidate pairs into connected previous/current components."""
    parent: dict[tuple[str, int], tuple[str, int]] = {}

    def find(node: tuple[str, int]) -> tuple[str, int]:
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != root:
            parent[node], node = root, parent[node]
        return root

    for candidate in candidates:
        row_node = ("r", candidate.previous_index)
        column_node = ("c", candidate.current_index)
        parent.setdefault(row_node, row_node)
        parent.setdefault(column_node, column_node)
        row_root, column_root = find(row_node), find(column_node)
        if row_root != column_root:
            parent[column_root] = row_root
    rows_by_root: dict[tuple[str, int], set[int]] = defaultdict(set)
    columns_by_root: dict[tuple[str, int], set[int]] = defaultdict(set)
    for kind, index in parent:
        root = find((kind, index))
        (rows_by_root if kind == "r" else columns_by_root)[root].add(index)
    return sorted(
        (tuple(sorted(rows_by_root[root])), tuple(sorted(columns_by_root[root])))
        for root in rows_by_root
    )

def _candidate_quality(candidate: _Candidate, config: TrackingConfig) -> float:
    """Return the unchanged pre-0.5.6 assignment quality as a diagnostic value."""
    center_score = (
        0.0
        if candidate.center_distance_nm is None
        else 1.0 / (1.0 + candidate.center_distance_nm)
    )
    return (
        100.0 * candidate.jaccard
        + 10.0 * candidate.shared_fraction
        + center_score
        + candidate.topology_similarity
        + (0.001 * candidate.guest_similarity if config.guest_tiebreak else 0.0)
        - 0.01 * candidate.gap_frames
    )

def _match_diagnostic(
    evidence: _Candidate | None,
    competitors: Sequence[_Candidate],
    config: TrackingConfig | None,
) -> dict[str, object]:
    if evidence is None or config is None:
        return {
            "candidate_count": 0,
            "score": None,
            "runner_up_score": None,
            "score_margin": None,
            "jaccard_margin": None,
            "shared_fraction_margin": None,
            "shared_waters_margin": None,
            "center_distance_margin_nm": None,
            "near_threshold": False,
            "status": "new",
        }
    score = _candidate_quality(evidence, config)
    other_scores = sorted(
        (
            _candidate_quality(item, config)
            for item in competitors
            if item is not evidence
        ),
        reverse=True,
    )
    runner_up = other_scores[0] if other_scores else None
    score_margin = None if runner_up is None else score - runner_up
    jaccard_margin = evidence.jaccard - config.min_jaccard
    fraction_margin = evidence.shared_fraction - config.min_shared_fraction
    waters_margin = evidence.shared_waters - config.min_shared_waters
    distance_margin = (
        None
        if config.max_center_distance_nm is None
        or evidence.center_distance_nm is None
        else config.max_center_distance_nm - evidence.center_distance_nm
    )
    tolerance = config.near_threshold_tolerance
    normalized_waters_margin = waters_margin / max(1, config.min_shared_waters)
    normalized_distance_margin = (
        None
        if distance_margin is None or config.max_center_distance_nm in (None, 0.0)
        else distance_margin / float(config.max_center_distance_nm)
    )
    near_threshold = (
        jaccard_margin <= tolerance
        or fraction_margin <= tolerance
        or normalized_waters_margin <= tolerance
        or (
            normalized_distance_margin is not None
            and normalized_distance_margin <= tolerance
        )
    )
    if evidence.gap_frames:
        status = "gap_bridge"
    elif near_threshold or (
        score_margin is not None and score_margin <= config.ambiguity_score_margin
    ):
        status = "ambiguous"
    else:
        status = "secure"
    return {
        "candidate_count": len(competitors),
        "score": score,
        "runner_up_score": runner_up,
        "score_margin": score_margin,
        "jaccard_margin": jaccard_margin,
        "shared_fraction_margin": fraction_margin,
        "shared_waters_margin": waters_margin,
        "center_distance_margin_nm": distance_margin,
        "near_threshold": near_threshold,
        "status": status,
    }

def _hungarian_columns(costs: Sequence[Sequence[float]]) -> list[int]:
    row_count = len(costs)
    column_count = len(costs[0]) if row_count else 0
    if row_count > column_count:
        raise ValueError("Hungarian assignment requires at least as many columns as rows.")
    u = [0.0] * (row_count + 1)
    v = [0.0] * (column_count + 1)
    matched_row = [0] * (column_count + 1)
    previous_column = [0] * (column_count + 1)
    tolerance = _ASSIGNMENT_TOLERANCE
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
    competitors: Sequence[_Candidate] = (),
    config: TrackingConfig | None = None,
) -> CageObservation:
    diagnostic = _match_diagnostic(evidence, competitors, config)
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
        match_shared_waters=None if evidence is None else evidence.shared_waters,
        match_center_distance_nm=(
            None if evidence is None else evidence.center_distance_nm
        ),
        match_topology_similarity=(
            None if evidence is None else evidence.topology_similarity
        ),
        match_candidate_count=diagnostic["candidate_count"],
        match_score=diagnostic["score"],
        match_runner_up_score=diagnostic["runner_up_score"],
        match_score_margin=diagnostic["score_margin"],
        match_jaccard_margin=diagnostic["jaccard_margin"],
        match_shared_fraction_margin=diagnostic["shared_fraction_margin"],
        match_shared_waters_margin=diagnostic["shared_waters_margin"],
        match_center_distance_margin_nm=diagnostic["center_distance_margin_nm"],
        match_near_threshold=diagnostic["near_threshold"],
        match_status=diagnostic["status"],
        match_diagnostic_source="" if evidence is None else "computed",
        gap_frames=0 if evidence is None else evidence.gap_frames,
        gap_time_ps=None if evidence is None else evidence.gap_time_ps,
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
    gap_time_ps: float | None = None,
    censored: bool = False,
    guest_ids_entered: Iterable[str] = (),
    guest_ids_exited: Iterable[str] = (),
    source_occupancy: str = "",
    destination_occupancy: str = "",
    guest_composition_before: Iterable[str] = (),
    guest_composition_after: Iterable[str] = (),
    evidence_status: str = "",
    water_conservation: float | None = None,
    center_distance_nm: float | None = None,
    persistence_frames: int = 0,
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
        gap_time_ps=gap_time_ps,
        censored=bool(censored),
        guest_ids_entered=_unique(guest_ids_entered),
        guest_ids_exited=_unique(guest_ids_exited),
        source_occupancy=str(source_occupancy),
        destination_occupancy=str(destination_occupancy),
        guest_composition_before=_unique(guest_composition_before),
        guest_composition_after=_unique(guest_composition_after),
        evidence_status=str(evidence_status),
        water_conservation=water_conservation,
        center_distance_nm=center_distance_nm,
        persistence_frames=int(persistence_frames),
    )

def _track_number(track_id: str) -> int:
    match = _TRACK_PATTERN.fullmatch(track_id)
    if match is None:
        raise ValueError(f"Invalid persistent track ID: {track_id!r}.")
    return int(match.group(1))

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

def _overlap_fraction(left: Iterable[int], right: Iterable[int]) -> float:
    left_set, right_set = set(left), set(right)
    return len(left_set & right_set) / len(left_set) if left_set else 0.0

def _maximum_candidate_distance(candidates: Sequence[_Candidate]) -> float | None:
    values = [
        item.center_distance_nm
        for item in candidates
        if item.center_distance_nm is not None
    ]
    return max(values) if values else None

def _lineage_evidence_status(
    conservation: float,
    candidates: Sequence[_Candidate],
    config: TrackingConfig,
) -> str:
    if conservation < config.min_shared_fraction:
        return "candidate_low_conservation"
    if not candidates or any(item.topology_similarity <= 0.0 for item in candidates):
        return "candidate_low_topology"
    return "candidate_supported"

def _independent_contribution_sets(
    parent_waters: Iterable[int], branches: Sequence[_CageState]
) -> tuple[frozenset[int], ...]:
    """Exclude waters shared by branches from evidence of a lineage change.

    Shared-face adjacency by itself is insufficient: each branch must supply
    at least ``min_shared_waters`` distinct source/target waters of its own.
    This conservative membership criterion does not infer a physical mechanism.
    """
    parent_set = set(parent_waters)
    counts: Counter[int] = Counter()
    for branch in branches:
        counts.update(parent_set.intersection(branch.water_atomids))
    return tuple(
        frozenset(water for water in parent_set.intersection(branch.water_atomids) if counts[water] == 1)
        for branch in branches
    )

def _independent_contributions(
    parent_waters: Iterable[int], branches: Sequence[_CageState]
) -> tuple[int, ...]:
    return tuple(len(waters) for waters in _independent_contribution_sets(parent_waters, branches))

def _lineage_continuation(
    before: _CageState,
    after: _CageState,
    box: np.ndarray | None,
    config: TrackingConfig,
) -> bool:
    """Revalidate one destination shell using the unchanged matching thresholds."""
    shared = len(before.water_atomids & after.water_atomids)
    union = len(before.water_atomids | after.water_atomids)
    smaller = min(len(before.water_atomids), len(after.water_atomids))
    if (
        shared < config.min_shared_waters
        or not union
        or not smaller
        or shared / union < config.min_jaccard
        or shared / smaller < config.min_shared_fraction
        or _multiset_similarity(before.topology, after.topology) <= 0.0
    ):
        return False
    return config.max_center_distance_nm is None or distance(
        np.asarray(before.center), np.asarray(after.center), box
    ) <= config.max_center_distance_nm

def _occupancy_class(guest_ids: Sequence[str]) -> str:
    count = len(guest_ids)
    return "empty" if count == 0 else "single" if count == 1 else "multiple"

def _occupancy_transition_kind(
    source: str,
    destination: str,
) -> EventKind | None:
    mapping: dict[tuple[str, str], EventKind] = {
        ("empty", "single"): "empty_to_occupied",
        ("empty", "multiple"): "empty_to_occupied",
        ("single", "empty"): "occupied_to_empty",
        ("multiple", "empty"): "occupied_to_empty",
        ("single", "multiple"): "single_to_multiple",
        ("multiple", "single"): "multiple_to_single",
    }
    return mapping.get((source, destination))

def _unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values))
