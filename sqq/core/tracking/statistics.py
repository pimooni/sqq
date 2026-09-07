"""Tabular statistics derived from cage tracking results."""

from __future__ import annotations

from collections import Counter, defaultdict
from math import inf
import re
from typing import Iterable, Mapping, Sequence

from ...models.tracking import (
    CageObservation,
    CageTrack,
    FrameStamp,
    Row,
    TargetSelection,
    TrackingResult,
)

_TRACK_PATTERN = re.compile(r"^t0*([1-9][0-9]*)$", re.IGNORECASE)
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

__all__ = [
    "cage_lineage_rows",
    "event_rows",
    "guest_event_rows",
    "guest_residence_lifetime_rows",
    "guest_residence_rows",
    "lifetime_distribution_rows",
    "lifetime_rows",
    "lifetime_survival_rows",
    "observation_rows",
    "occupancy_state_lifetime_rows",
    "occupancy_transition_rows",
    "population_rows",
    "tracking_quality_rows",
]

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
            "match_shared_waters": item.match_shared_waters,
            "match_center_distance_nm": item.match_center_distance_nm,
            "match_topology_similarity": item.match_topology_similarity,
            "match_candidate_count": item.match_candidate_count,
            "match_score": item.match_score,
            "match_runner_up_score": item.match_runner_up_score,
            "match_score_margin": item.match_score_margin,
            "match_jaccard_margin": item.match_jaccard_margin,
            "match_shared_fraction_margin": item.match_shared_fraction_margin,
            "match_shared_waters_margin": item.match_shared_waters_margin,
            "match_center_distance_margin_nm": item.match_center_distance_margin_nm,
            "match_near_threshold": item.match_near_threshold,
            "match_status": item.match_status,
            "match_diagnostic_source": item.match_diagnostic_source,
            "gap_frames": item.gap_frames,
            "gap_time_ps": item.gap_time_ps,
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
        matched = [
            item
            for item in track.observations
            if item.match_status in {"secure", "ambiguous", "gap_bridge"}
        ]
        match_scores = [
            float(item.match_score) for item in matched if item.match_score is not None
        ]
        score_margins = [
            float(item.match_score_margin)
            for item in matched
            if item.match_score_margin is not None
        ]
        jaccards = [
            float(item.match_jaccard)
            for item in matched
            if item.match_jaccard is not None
        ]
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
        lifetime_lower_ps = observed_span_ps
        lifetime_upper_ps = _interval_upper_bound(
            data.frames,
            first_position,
            last_position,
            left_censored=track.left_censored,
            right_censored=track.right_censored,
        )
        occupancy_time_ps = _sampled_occupancy_time(
            data.frames,
            (positions[item.frame_index] for item in track.observations),
        )
        gap_time_ps = (
            0.0
            if track.gap_frames == 0
            else _sum_optional(
                item.gap_time_ps
                for item in track.observations
                if item.gap_frames > 0
            )
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
                "gap_time_ps": gap_time_ps,
                "observed_span_ps": observed_span_ps,
                "lifetime_ps": lifetime_ps,
                "lifetime_lower_ps": lifetime_lower_ps,
                "lifetime_upper_ps": lifetime_upper_ps,
                "occupancy_time_ps": occupancy_time_ps,
                "duration_status": _duration_status(
                    track,
                    lifetime_lower_ps,
                    lifetime_upper_ps,
                ),
                "left_censored": track.left_censored,
                "right_censored": track.right_censored,
                "matching_observations": len(matched),
                "match_score_min": min(match_scores) if match_scores else None,
                "match_score_mean": (
                    sum(match_scores) / len(match_scores) if match_scores else None
                ),
                "match_score_margin_min": (
                    min(score_margins) if score_margins else None
                ),
                "match_score_margin_mean": (
                    sum(score_margins) / len(score_margins)
                    if score_margins
                    else None
                ),
                "match_jaccard_min": min(jaccards) if jaccards else None,
                "match_jaccard_mean": (
                    sum(jaccards) / len(jaccards) if jaccards else None
                ),
                "ambiguous_observations": sum(
                    item.match_status == "ambiguous" for item in matched
                ),
                "gap_bridge_observations": sum(
                    item.match_status == "gap_bridge" for item in matched
                ),
                "near_threshold_observations": sum(
                    item.match_near_threshold for item in matched
                ),
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
    *,
    samples: Sequence[Row] | None = None,
) -> list[Row]:
    """Aggregate exact lifetime samples into a standalone discrete distribution.

    ``samples`` accepts the already computed ``lifetime_rows(data)`` so one
    table set does not rebuild the per-track rows several times.
    """
    samples = lifetime_rows(data) if samples is None else list(samples)
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
            "gap_time_ps": item.gap_time_ps,
            "censored": item.censored,
            "guest_ids_entered": ",".join(item.guest_ids_entered),
            "guest_ids_exited": ",".join(item.guest_ids_exited),
            "source_occupancy": item.source_occupancy,
            "destination_occupancy": item.destination_occupancy,
            "guest_composition_before": ",".join(item.guest_composition_before),
            "guest_composition_after": ",".join(item.guest_composition_after),
            "evidence_status": item.evidence_status,
            "water_conservation": item.water_conservation,
            "center_distance_nm": item.center_distance_nm,
            "persistence_frames": item.persistence_frames,
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

def tracking_quality_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    """Return one diagnostics row for every accepted cross-frame match."""
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
            "candidate_count": item.match_candidate_count,
            "match_score": item.match_score,
            "runner_up_score": item.match_runner_up_score,
            "score_margin": item.match_score_margin,
            "match_jaccard": item.match_jaccard,
            "match_shared_fraction": item.match_shared_fraction,
            "match_shared_waters": item.match_shared_waters,
            "match_center_distance_nm": item.match_center_distance_nm,
            "match_topology_similarity": item.match_topology_similarity,
            "jaccard_margin": item.match_jaccard_margin,
            "shared_fraction_margin": item.match_shared_fraction_margin,
            "shared_waters_margin": item.match_shared_waters_margin,
            "center_distance_margin_nm": item.match_center_distance_margin_nm,
            "near_threshold": item.match_near_threshold,
            "gap_frames": item.gap_frames,
            "gap_time_ps": item.gap_time_ps,
            "match_status": item.match_status,
            "diagnostic_source": item.match_diagnostic_source,
        }
        for item in data.observations
        if item.match_status in {"secure", "ambiguous", "gap_bridge"}
    ]

def guest_event_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    """Return guest identity entry, exit, exchange, and unresolved-gap events."""
    kinds = {
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
    return [
        {
            "target": _target_label(data),
            "event_id": item.event_id,
            "event": item.kind,
            "frame_index": item.frame_index,
            "frame": item.frame_name,
            "time_ps": item.time_ps,
            "track_id": next(iter(item.track_ids), ""),
            "cage_type": next(iter(item.destination_cage_types or item.source_cage_types), ""),
            "phase": next(iter(item.destination_phases or item.source_phases), ""),
            "guest_ids_entered": ",".join(item.guest_ids_entered),
            "guest_ids_exited": ",".join(item.guest_ids_exited),
            "guest_composition_before": ",".join(item.guest_composition_before),
            "guest_composition_after": ",".join(item.guest_composition_after),
            "source_occupancy": item.source_occupancy,
            "destination_occupancy": item.destination_occupancy,
            "gap_frames": item.gap_frames,
            "gap_time_ps": item.gap_time_ps,
            "censored": item.censored,
        }
        for item in data.events
        if item.kind in kinds
    ]

def guest_residence_lifetime_rows(
    data: TrackingResult | TargetSelection,
) -> list[Row]:
    """Return visible guest episodes without treating cage gaps as guest exits."""
    target = _target_label(data)
    positions = {frame.frame_index: index for index, frame in enumerate(data.frames)}
    rows: list[Row] = []
    for track in data.tracks:
        track_positions = {positions[item.frame_index] for item in track.observations}
        first_track_position = positions[track.first.frame_index]
        last_track_position = positions[track.last.frame_index]
        by_guest: dict[str, list[CageObservation]] = defaultdict(list)
        for item in track.observations:
            for guest_id in item.guest_ids:
                by_guest[guest_id].append(item)
        for guest_id in sorted(by_guest):
            residences = _contiguous_episodes(by_guest[guest_id], positions)
            for residence_index, residence in enumerate(residences, start=1):
                first, last = residence[0], residence[-1]
                first_position = positions[first.frame_index]
                last_position = positions[last.frame_index]
                left_censored = bool(
                    track.left_censored and first.frame_index == track.first.frame_index
                )
                right_censored = bool(
                    track.right_censored and last.frame_index == track.last.frame_index
                )
                left_gap_unresolved = (
                    first_position > first_track_position
                    and first_position - 1 not in track_positions
                )
                right_gap_unresolved = (
                    last_position < last_track_position
                    and last_position + 1 not in track_positions
                )
                gap_censored = left_gap_unresolved or right_gap_unresolved
                lower = _time_difference(first.time_ps, last.time_ps)
                upper = _interval_upper_bound(
                    data.frames,
                    first_position,
                    last_position,
                    left_censored=left_censored or left_gap_unresolved,
                    right_censored=right_censored or right_gap_unresolved,
                )
                rows.append(
                    {
                        "target": target,
                        "track_id": track.track_id,
                        "guest_id": guest_id,
                        "residence_index": residence_index,
                        "start_frame_index": first.frame_index,
                        "end_frame_index": last.frame_index,
                        "start_time_ps": first.time_ps,
                        "end_time_ps": last.time_ps,
                        "observed_frames": len(residence),
                        "observed_span_ps": lower,
                        "residence_lifetime_lower_ps": lower,
                        "residence_lifetime_upper_ps": upper,
                        "residence_time_ps": _sampled_occupancy_time(
                            data.frames,
                            (positions[item.frame_index] for item in residence),
                        ),
                        "cage_types": ",".join(
                            _unique(item.cage_type for item in residence)
                        ),
                        "phases": ",".join(
                            _unique(
                                phase for item in residence for phase in item.phase_labels
                            )
                        ),
                        "duration_status": _residence_duration_status(
                            lower, upper, left_censored, right_censored, gap_censored
                        ),
                        "left_censored": left_censored,
                        "right_censored": right_censored,
                        "left_gap_unresolved": left_gap_unresolved,
                        "right_gap_unresolved": right_gap_unresolved,
                        "gap_censored": gap_censored,
                    }
                )
    return rows

def occupancy_state_lifetime_rows(
    data: TrackingResult | TargetSelection,
) -> list[Row]:
    """Return visible exact-guest states with unresolved cage-gap boundaries."""
    target = _target_label(data)
    positions = {frame.frame_index: index for index, frame in enumerate(data.frames)}
    rows: list[Row] = []
    for track in data.tracks:
        track_positions = {positions[item.frame_index] for item in track.observations}
        first_track_position = positions[track.first.frame_index]
        last_track_position = positions[track.last.frame_index]
        residences: list[list[CageObservation]] = []
        for item in track.observations:
            previous = residences[-1][-1] if residences else None
            if (
                previous is None
                or positions[item.frame_index] != positions[previous.frame_index] + 1
                or tuple(sorted(item.guest_ids)) != tuple(sorted(previous.guest_ids))
            ):
                residences.append([item])
            else:
                residences[-1].append(item)
        for residence_index, residence in enumerate(residences, start=1):
            first, last = residence[0], residence[-1]
            first_position = positions[first.frame_index]
            last_position = positions[last.frame_index]
            left_censored = bool(
                track.left_censored and first.frame_index == track.first.frame_index
            )
            right_censored = bool(
                track.right_censored and last.frame_index == track.last.frame_index
            )
            left_gap_unresolved = (
                first_position > first_track_position
                and first_position - 1 not in track_positions
            )
            right_gap_unresolved = (
                last_position < last_track_position
                and last_position + 1 not in track_positions
            )
            gap_censored = left_gap_unresolved or right_gap_unresolved
            lower = _time_difference(first.time_ps, last.time_ps)
            upper = _interval_upper_bound(
                data.frames,
                first_position,
                last_position,
                left_censored=left_censored or left_gap_unresolved,
                right_censored=right_censored or right_gap_unresolved,
            )
            composition = tuple(sorted(first.guest_ids))
            rows.append(
                {
                    "target": target,
                    "track_id": track.track_id,
                    "occupancy_state": _occupancy_class(composition),
                    "guest_composition": ",".join(composition),
                    "residence_index": residence_index,
                    "start_frame_index": first.frame_index,
                    "end_frame_index": last.frame_index,
                    "start_time_ps": first.time_ps,
                    "end_time_ps": last.time_ps,
                    "observed_frames": len(residence),
                    "observed_span_ps": lower,
                    "residence_lifetime_lower_ps": lower,
                    "residence_lifetime_upper_ps": upper,
                    "residence_time_ps": _sampled_occupancy_time(
                        data.frames,
                        (positions[item.frame_index] for item in residence),
                    ),
                    "cage_types": ",".join(
                        _unique(item.cage_type for item in residence)
                    ),
                    "phases": ",".join(
                        _unique(
                            phase for item in residence for phase in item.phase_labels
                        )
                    ),
                    "duration_status": _residence_duration_status(
                        lower, upper, left_censored, right_censored, gap_censored
                    ),
                    "left_censored": left_censored,
                    "right_censored": right_censored,
                    "left_gap_unresolved": left_gap_unresolved,
                    "right_gap_unresolved": right_gap_unresolved,
                    "gap_censored": gap_censored,
                }
            )
    return rows

def occupancy_transition_rows(
    data: TrackingResult | TargetSelection,
) -> list[Row]:
    """Aggregate directed exact-composition occupancy transitions."""
    counts: Counter[tuple[str, str]] = Counter()
    tracks_by_edge: dict[tuple[str, str], set[str]] = defaultdict(set)
    exposure: Counter[str] = Counter()
    for track in data.tracks:
        ordered = sorted(track.observations, key=_observation_sort_key)
        for left, right in zip(ordered, ordered[1:]):
            if right.frame_index != left.frame_index + 1:
                continue
            source = _occupancy_state_label(left.guest_ids)
            destination = _occupancy_state_label(right.guest_ids)
            delta = _time_difference(left.time_ps, right.time_ps)
            if delta is not None:
                exposure[source] += delta
            if source == destination:
                continue
            edge = (source, destination)
            counts[edge] += 1
            tracks_by_edge[edge].add(track.track_id)
    outgoing = Counter()
    for (source, _destination), count in counts.items():
        outgoing[source] += count
    return [
        {
            "target": _target_label(data),
            "source_state": source,
            "destination_state": destination,
            "transition_count": count,
            "independent_track_count": len(tracks_by_edge[(source, destination)]),
            "source_residence_ps": exposure.get(source, 0.0),
            "conditional_probability": count / outgoing[source],
            "rate_per_ns": (
                None
                if exposure.get(source, 0.0) <= 0.0
                else count / (exposure[source] / 1000.0)
            ),
        }
        for (source, destination), count in sorted(counts.items())
    ]

_SURVIVAL_TIME_BASES: tuple[tuple[str, str], ...] = (
    # Death placed at the first absent selected frame (interval upper bound;
    # the compatibility ``lifetime_ps`` convention).
    ("first_absence", "lifetime_ps"),
    # Death placed at the last observed selected frame (interval lower bound).
    ("last_observed", "lifetime_lower_ps"),
)


def lifetime_survival_rows(
    data: TrackingResult | TargetSelection,
    *,
    samples: Sequence[Row] | None = None,
) -> list[Row]:
    """Return bounding Kaplan-Meier tables on the selected-frame grid.

    Death is placed at either the first absent or last observed frame.
    Right-censored tracks remain; left-censored and untimed tracks are counted
    but excluded.  ``samples`` may contain precomputed lifetime rows.
    """
    samples = lifetime_rows(data) if samples is None else list(samples)
    if not samples:
        return []
    target = _target_label(data)
    rows: list[Row] = []
    for basis, event_field in _SURVIVAL_TIME_BASES:
        rows.extend(_survival_table(target, samples, basis, event_field))
    return rows


def _survival_table(
    target: str,
    samples: Sequence[Row],
    basis: str,
    event_field: str,
) -> list[Row]:
    left_excluded = sum(bool(row["left_censored"]) for row in samples)
    deaths_at: Counter[float] = Counter()
    censored_at: Counter[float] = Counter()
    time_unavailable = 0
    for row in samples:
        if bool(row["left_censored"]):
            continue
        if bool(row["right_censored"]):
            censor_time = row.get("lifetime_lower_ps")
            if censor_time is None:
                time_unavailable += 1
                continue
            censored_at[float(censor_time)] += 1
            continue
        event_time = row.get(event_field)
        if event_time is None:
            time_unavailable += 1
            continue
        deaths_at[float(event_time)] += 1
    included = sum(deaths_at.values()) + sum(censored_at.values())
    times = sorted(set(deaths_at) | set(censored_at))
    if not times:
        return [
            {
                "target": target,
                "time_basis": basis,
                "time_ps": None,
                "number_at_risk": 0,
                "deaths": 0,
                "right_censored": 0,
                "left_censored_excluded": left_excluded,
                "time_unavailable_excluded": time_unavailable,
                "survival_probability": None,
                "greenwood_variance": None,
                "greenwood_standard_error": None,
                "confidence_lower": None,
                "confidence_upper": None,
            }
        ]
    at_risk = included
    survival = 1.0
    greenwood_sum = 0.0
    rows: list[Row] = []
    for index, time_ps in enumerate(times):
        deaths = deaths_at[time_ps]
        right_censored = censored_at[time_ps]
        if deaths and at_risk:
            survival *= 1.0 - deaths / at_risk
            if at_risk > deaths:
                greenwood_sum += deaths / (at_risk * (at_risk - deaths))
        variance = survival * survival * greenwood_sum
        standard_error = variance ** 0.5
        rows.append(
            {
                "target": target,
                "time_basis": basis,
                "time_ps": time_ps,
                "number_at_risk": at_risk,
                "deaths": deaths,
                "right_censored": right_censored,
                "left_censored_excluded": left_excluded if index == 0 else 0,
                "time_unavailable_excluded": time_unavailable if index == 0 else 0,
                "survival_probability": survival,
                "greenwood_variance": variance,
                "greenwood_standard_error": standard_error,
                "confidence_lower": max(0.0, survival - 1.96 * standard_error),
                "confidence_upper": min(1.0, survival + 1.96 * standard_error),
            }
        )
        at_risk -= deaths + right_censored
    return rows

def cage_lineage_rows(data: TrackingResult | TargetSelection) -> list[Row]:
    """Return candidate and one-frame-confirmed cage lineage events."""
    kinds = {
        "split",
        "merge",
        "split_candidate",
        "merge_candidate",
        "split_confirmed",
        "merge_confirmed",
    }
    return [
        {
            "target": _target_label(data),
            "event_id": item.event_id,
            "event": item.kind,
            "frame_index": item.frame_index,
            "frame": item.frame_name,
            "time_ps": item.time_ps,
            "parent_track_ids": ",".join(item.source_track_ids),
            "child_track_ids": ",".join(item.destination_track_ids),
            "source_cage_types": ",".join(item.source_cage_types),
            "destination_cage_types": ",".join(item.destination_cage_types),
            "water_conservation": item.water_conservation,
            "center_distance_nm": item.center_distance_nm,
            "persistence_frames": item.persistence_frames,
            "confirmation_status": (
                "confirmed" if item.kind.endswith("_confirmed") else "candidate"
            ),
            "censored": item.censored,
        }
        for item in data.events
        if item.kind in kinds
    ]

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

def _observation_sort_key(item: CageObservation) -> tuple[int, int, str]:
    return item.frame_index, _track_number(item.track_id), item.local_cage_id

def _track_number(track_id: str) -> int:
    match = _TRACK_PATTERN.fullmatch(track_id)
    if match is None:
        raise ValueError(f"Invalid persistent track ID: {track_id!r}.")
    return int(match.group(1))

def _phase_sort_key(value: str) -> tuple[int, str]:
    return _PHASE_ORDER.get(value, 100), value

def _occupancy_class(guest_ids: Sequence[str]) -> str:
    count = len(guest_ids)
    return "empty" if count == 0 else "single" if count == 1 else "multiple"

def _occupancy_state_label(guest_ids: Sequence[str]) -> str:
    composition = tuple(sorted(guest_ids))
    state = _occupancy_class(composition)
    return state if not composition else f"{state}:" + "+".join(composition)

def _time_difference(start: float | None, end: float | None) -> float | None:
    return None if start is None or end is None else float(end - start)

def _sum_optional(values: Iterable[float | None]) -> float | None:
    present = [float(value) for value in values if value is not None]
    return sum(present) if present else None

def _interval_upper_bound(
    frames: Sequence[FrameStamp],
    first_position: int | None,
    last_position: int | None,
    *,
    left_censored: bool,
    right_censored: bool,
) -> float | None:
    if (
        first_position is None
        or last_position is None
        or left_censored
        or right_censored
        or first_position <= 0
        or last_position + 1 >= len(frames)
    ):
        return None
    previous_absent = frames[first_position - 1].time_ps
    next_absent = frames[last_position + 1].time_ps
    return _time_difference(previous_absent, next_absent)

def _sampled_occupancy_time(
    frames: Sequence[FrameStamp],
    positions: Iterable[int],
) -> float | None:
    selected = sorted(set(int(position) for position in positions))
    if not selected:
        return 0.0
    times = [frame.time_ps for frame in frames]
    if any(value is None for value in times):
        return None
    numeric = [float(value) for value in times if value is not None]
    if len(numeric) == 1:
        return 0.0
    widths: list[float] = []
    for position in selected:
        if position == 0:
            width = numeric[1] - numeric[0]
        elif position == len(numeric) - 1:
            width = numeric[-1] - numeric[-2]
        else:
            width = (numeric[position + 1] - numeric[position - 1]) / 2.0
        widths.append(max(0.0, width))
    return float(sum(widths))

def _interval_duration_status(
    lower: float | None,
    upper: float | None,
    left_censored: bool,
    right_censored: bool,
) -> str:
    if lower is None:
        return "time_unavailable"
    if left_censored and right_censored:
        return "left_and_right_censored_lower_bound"
    if left_censored:
        return "left_censored_lower_bound"
    if right_censored:
        return "right_censored_lower_bound"
    if upper is None:
        return "interval_bound_unavailable"
    if abs(upper - lower) <= 1.0e-12:
        return "exact"
    return "interval_censored"

def _residence_duration_status(
    lower: float | None,
    upper: float | None,
    left_censored: bool,
    right_censored: bool,
    gap_censored: bool,
) -> str:
    # A missing cage does not establish absence of its guest or guest state.
    # Retain the visible span as a lower bound, but do not invent a finite end.
    if lower is not None and gap_censored:
        return "gap_censored_lower_bound"
    return _interval_duration_status(lower, upper, left_censored, right_censored)

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

def _duration_status(
    track: CageTrack,
    lower: float | None,
    upper: float | None,
) -> str:
    return _interval_duration_status(
        lower,
        upper,
        track.left_censored,
        track.right_censored,
    )
