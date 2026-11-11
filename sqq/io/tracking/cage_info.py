"""Shared, bounded-memory per-frame cage lists for Analyze and Track."""

from __future__ import annotations

from itertools import groupby
from dataclasses import replace
from pathlib import Path
import re
from typing import Sequence
import warnings

from ...core.tracking.statistics import iter_observation_rows
from ...models.tracking import TargetSelection, TrackingResult


def write_cage_info_frames(
    result: TrackingResult,
    destinations: Sequence[tuple[TargetSelection | None, Path]],
) -> None:
    """Write each frame once per selection, preserving full-system neighbors."""
    from .state import _OBSERVATION_FIELDS, _atomic_write_csv

    if not destinations:
        return
    names: dict[str, str] = {}
    for frame in result.frames:
        name = frame.frame_name
        if (not name or name in {".", ".."} or
                re.search(r'[<>:"/\\|?*\x00-\x1f]', name) or name.endswith((".", " "))):
            raise ValueError(f"Unsafe cage-info frame name: {name!r}.")
        if name.casefold() in names:
            raise ValueError(f"Duplicate cage-info frame name: {name!r}.")
        names[name.casefold()] = name
    selections = [
        (None if selection is None or selection.target.kind == "all" else
         frozenset(track.track_id for track in selection.tracks), Path(root))
        for selection, root in destinations
    ]
    fields = tuple(field for field in _OBSERVATION_FIELDS if field != "target") + (
        "connectivity", "neighbor_count",
    )
    grouped = iter(groupby(result.observations, key=lambda item: item.frame_index))
    pending = next(grouped, None)
    unavailable = False
    for frame in result.frames:
        pairs = []
        if pending is not None and pending[0] == frame.frame_index:
            items = tuple(pending[1])
            frame_data = replace(result, _observation_view=items)
            pairs = list(zip(items, iter_observation_rows(frame_data)))
            pending = next(grouped, None)
        if pending is not None and pending[0] < frame.frame_index:
            raise ValueError("Cage observations are not in frame order.")
        by_local = {item.local_cage_id: item for item, _row in pairs}
        if len(by_local) != len(pairs):
            raise ValueError("A cage-info frame contains duplicate local cage IDs.")
        for item, row in pairs:
            neighbors = item.neighbor_cage_ids
            if neighbors is None:
                row.update(connectivity="NA", neighbor_count="NA")
                unavailable = True
                continue
            if item.local_cage_id in neighbors or len(set(neighbors)) != len(neighbors):
                raise ValueError("Invalid cage-info connectivity: self-link or duplicate neighbor.")
            for neighbor in neighbors:
                other = by_local.get(neighbor)
                if other is None or other.neighbor_cage_ids is None or item.local_cage_id not in other.neighbor_cage_ids:
                    raise ValueError("Cage-info connectivity must be symmetric within the same frame.")
            identifiers = sorted(
                (by_local[neighbor].track_id for neighbor in neighbors),
                key=lambda identifier: int(identifier[1:]),
            )
            row.update(connectivity=",".join(identifiers), neighbor_count=len(identifiers))
        for track_ids, root in selections:
            path = root / "gro" / frame.frame_name / f"cage_info_{frame.frame_name}.csv"
            rows = (row for item, row in pairs if track_ids is None or item.track_id in track_ids)
            _atomic_write_csv(path, rows, fields)
    if pending is not None:
        raise ValueError("Cage observations refer to an unknown frame.")
    if unavailable:
        warnings.warn(
            "Shared-face connectivity is unavailable in legacy Track state; "
            "cage_info uses NA, not zero neighbors. Reanalyze raw input to obtain connectivity.",
            UserWarning, stacklevel=2,
        )
