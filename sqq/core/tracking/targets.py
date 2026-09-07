"""Tracking-target parsing and lifecycle selection."""

from __future__ import annotations

import re
from typing import Iterable

from ...models.cage_type import canonical_cage_type
from ...models.tracking import CageTrack, TargetSelection, TargetSpec, TrackingResult

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

__all__ = ["parse_targets", "select_targets"]

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
