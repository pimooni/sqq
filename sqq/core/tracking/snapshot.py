"""Frame-to-snapshot conversion and normalized per-cage state."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np

from ...models import Cage, FrameResult, Ring
from ...models.tracking import FrameStamp, TrackCageSnapshot, TrackFrameSnapshot

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

__all__ = ["snapshot_from_frame_result"]


@dataclass(frozen=True, slots=True)
class _CageState:
    local_cage_id: str
    cage_type: str
    phase: str
    phase_labels: tuple[str, ...]
    water_atomids: frozenset[int]
    center: tuple[float, float, float]
    topology: tuple[int, ...]
    guest_ids: tuple[str, ...]

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

def _state_sort_key(state: _CageState) -> tuple[object, ...]:
    return (
        min(state.water_atomids, default=-1),
        tuple(sorted(state.water_atomids)),
        state.cage_type,
        state.local_cage_id,
    )

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
