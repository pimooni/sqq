"""F3/F4 orientational order parameters."""

from __future__ import annotations

from itertools import combinations
from math import atan2, cos, radians
from typing import Any

import numpy as np

from ...models import Atom, F3F4Result, Frame, GraphResult, Water, WaterOrder
from ..pbc import minimum_image
from .steinhardt import (
    build_graph_neighbor_vector_cache,
    build_q_candidate_cache,
    normalize_q_degree,
    normalize_q_neighbor_mode,
    q_values_from_candidates,
    resolve_q_neighbor_count,
)


TETRAHEDRAL_COS2 = cos(radians(109.47)) ** 2


def compute_f3f4(frame: Frame, waters: list[Water], graph: GraphResult, focus_resids: set[int] | None = None) -> F3F4Result:
    """Compute per-water F3/F4 with Q_l disabled for backward-compatible callers."""
    return compute_order_parameters(
        frame,
        waters,
        graph,
        f3f4_enabled=True,
        q_enabled=False,
        focus_resids=focus_resids,
    )


def compute_order_parameters(
    frame: Frame,
    waters: list[Water],
    graph: GraphResult,
    *,
    f3f4_enabled: bool | None = None,
    f3_enabled: bool | None = None,
    f4_enabled: bool | None = None,
    q_enabled: bool = True,
    q_neighbor_mode: str = "graph",
    q_cutoff_nm: float = 0.35,
    q_n_neighbor: int | None = None,
    q_degree: Any = None,
    focus_resids: set[int] | None = None,
) -> F3F4Result:
    """Compute per-water F3/F4 and requested Steinhardt Q_l values."""
    legacy_f3f4_enabled = True if f3f4_enabled is None else bool(f3f4_enabled)
    resolved_f3_enabled = legacy_f3f4_enabled if f3_enabled is None else bool(f3_enabled)
    resolved_f4_enabled = legacy_f3f4_enabled if f4_enabled is None else bool(f4_enabled)
    water_by_oxygen = {water.oxygen: water for water in waters}
    focus_resids = set(focus_resids or set())
    rows: list[WaterOrder] = []
    f3_values: list[float] = []
    f4_values: list[float] = []
    f3_focus_values: list[float] = []
    f4_focus_values: list[float] = []
    q_degree_resolved = normalize_q_degree(q_degree) if q_enabled else ()
    q_value_lists: dict[int, list[float]] = {degree: [] for degree in q_degree_resolved}
    q_focus_value_lists: dict[int, list[float]] = {degree: [] for degree in q_degree_resolved}
    q_mode = normalize_q_neighbor_mode(q_neighbor_mode)
    q_fixed_neighbor = resolve_q_neighbor_count(q_mode, q_n_neighbor)
    graph_vectors = (
        build_graph_neighbor_vector_cache(frame, waters, graph)
        if resolved_f3_enabled or resolved_f4_enabled or (q_degree_resolved and q_mode == "graph")
        else {}
    )
    hydrogen_vectors = (
        build_hydrogen_vector_cache(frame, waters)
        if resolved_f4_enabled
        else {}
    )
    q_candidates = (
        build_q_candidate_cache(
            frame,
            waters,
            graph,
            mode=q_mode,
            cutoff_nm=q_cutoff_nm,
            graph_vectors=graph_vectors,
        )
        if q_degree_resolved
        else {}
    )

    for water in waters:
        neighbors = [idx for idx in sorted(graph.adjacency.get(water.oxygen, set())) if idx in water_by_oxygen]
        f3 = f3_for_water(
            frame.atoms,
            water.oxygen,
            neighbors,
            frame.box,
            neighbor_vectors=graph_vectors.get(water.oxygen),
        ) if resolved_f3_enabled else None
        f4 = f4_for_water(
            frame.atoms,
            water,
            [water_by_oxygen[idx] for idx in neighbors],
            frame.box,
            neighbor_vectors=graph_vectors.get(water.oxygen),
            hydrogen_vectors=hydrogen_vectors,
        ) if resolved_f4_enabled else None
        q_values, q_neighbors = q_values_from_candidates(
            q_candidates.get(water.oxygen, ()),
            n_neighbor=q_fixed_neighbor,
            degree=q_degree_resolved,
        ) if q_degree_resolved else ({}, 0)
        if f3 is not None:
            f3_values.append(f3)
            if water.resid in focus_resids:
                f3_focus_values.append(f3)
        if f4 is not None:
            f4_values.append(f4)
            if water.resid in focus_resids:
                f4_focus_values.append(f4)
        for degree, value in q_values.items():
            if value is None:
                continue
            q_value_lists.setdefault(degree, []).append(value)
            if water.resid in focus_resids:
                q_focus_value_lists.setdefault(degree, []).append(value)
        atom = frame.atoms[water.oxygen]
        rows.append(
            WaterOrder(
                oxygen=water.oxygen,
                resid=water.resid,
                atomid=atom.atomid,
                xyz=atom.xyz,
                f3=f3,
                f4=f4,
                q_values=dict(q_values),
                q_neighbors=q_neighbors,
            )
        )

    q_means = {degree: mean_or_none(q_value_lists.get(degree, [])) for degree in q_degree_resolved}
    q_valid_counts = {degree: len(q_value_lists.get(degree, [])) for degree in q_degree_resolved}
    q_focus_means = {
        degree: mean_or_none(q_focus_value_lists.get(degree, [])) if focus_resids else None
        for degree in q_degree_resolved
    }
    q_focus_valid_counts = {degree: len(q_focus_value_lists.get(degree, [])) for degree in q_degree_resolved}
    return F3F4Result(
        per_water=tuple(rows),
        f3_mean=mean_or_none(f3_values),
        f4_mean=mean_or_none(f4_values),
        f3_valid=len(f3_values),
        f4_valid=len(f4_values),
        focus_resids=tuple(sorted(focus_resids)),
        f3_focus_mean=mean_or_none(f3_focus_values) if focus_resids else None,
        f4_focus_mean=mean_or_none(f4_focus_values) if focus_resids else None,
        f3_focus_valid=len(f3_focus_values),
        f4_focus_valid=len(f4_focus_values),
        q_degree=q_degree_resolved,
        q_means=q_means,
        q_valid_counts=q_valid_counts,
        q_focus_means=q_focus_means,
        q_focus_valid_counts=q_focus_valid_counts,
        q_neighbor_mode=q_mode,
        q_cutoff_nm=q_cutoff_nm if q_degree_resolved else None,
        q_n_neighbor=q_fixed_neighbor if q_degree_resolved else None,
    )


def f3_for_water(
    atoms: list[Atom],
    oxygen: int,
    neighbors: list[int],
    box: np.ndarray | None,
    neighbor_vectors: dict[int, np.ndarray] | None = None,
) -> float | None:
    """F3: average squared deviation over neighbor-neighbor angles."""
    if len(neighbors) < 2:
        return None
    center = atoms[oxygen].xyz
    terms = []
    for a, b in combinations(neighbors, 2):
        va = neighbor_vectors[a] if neighbor_vectors is not None else minimum_image(atoms[a].xyz - center, box)
        vb = neighbor_vectors[b] if neighbor_vectors is not None else minimum_image(atoms[b].xyz - center, box)
        norm = float(np.linalg.norm(va) * np.linalg.norm(vb))
        if norm <= 1e-12:
            continue
        r = float(np.dot(va, vb) / norm)
        terms.append((r * abs(r) + TETRAHEDRAL_COS2) ** 2)
    return mean_or_none(terms)


def f4_for_water(
    atoms: list[Atom],
    water: Water,
    neighbors: list[Water],
    box: np.ndarray | None,
    *,
    neighbor_vectors: dict[int, np.ndarray] | None = None,
    hydrogen_vectors: dict[int, dict[int, np.ndarray]] | None = None,
) -> float | None:
    """F4 with reusable local O-H and O-O vectors when supplied."""
    if len(water.hydrogens) < 2:
        return None
    terms = []
    for neighbor in neighbors:
        if len(neighbor.hydrogens) < 2:
            continue
        h1, h2 = farthest_hydrogen_pair(atoms, water, neighbor, box)
        oo_vector = None if neighbor_vectors is None else neighbor_vectors.get(neighbor.oxygen)
        first_h_vector = None if hydrogen_vectors is None else hydrogen_vectors.get(water.oxygen, {}).get(h1)
        second_h_vector = None if hydrogen_vectors is None else hydrogen_vectors.get(neighbor.oxygen, {}).get(h2)
        if oo_vector is not None and first_h_vector is not None and second_h_vector is not None:
            angle = dihedral_from_local_vectors(first_h_vector, oo_vector, second_h_vector)
        else:
            angle = dihedral_h_o_o_h(atoms, h1, water.oxygen, neighbor.oxygen, h2, box)
        terms.append(cos(3.0 * angle))
    return mean_or_none(terms)


def build_hydrogen_vector_cache(frame: Frame, waters: list[Water]) -> dict[int, dict[int, np.ndarray]]:
    """Cache each selected water's local O-H vectors for F4 dihedrals."""
    return {
        water.oxygen: {
            hydrogen: minimum_image(frame.atoms[hydrogen].xyz - frame.atoms[water.oxygen].xyz, frame.box)
            for hydrogen in water.hydrogens[:2]
        }
        for water in waters
        if water.hydrogens
    }


def dihedral_from_local_vectors(
    first_h_vector: np.ndarray,
    oo_vector: np.ndarray,
    second_h_vector: np.ndarray,
) -> float:
    """Evaluate H-O-O-H after reusing the three historical minimum-image vectors."""
    origin = np.zeros(3, dtype=float)
    return dihedral(
        np.asarray(first_h_vector, dtype=float),
        origin,
        np.asarray(oo_vector, dtype=float),
        np.asarray(oo_vector, dtype=float) + np.asarray(second_h_vector, dtype=float),
    )


def farthest_hydrogen_pair(atoms: list[Atom], water_a: Water, water_b: Water, box: np.ndarray | None) -> tuple[int, int]:
    """Match the shell script by choosing the largest inter-water H-H distance."""
    best_pair = (water_a.hydrogens[0], water_b.hydrogens[0])
    best_dist2 = -1.0
    for h_a in water_a.hydrogens[:2]:
        for h_b in water_b.hydrogens[:2]:
            delta = minimum_image(atoms[h_b].xyz - atoms[h_a].xyz, box)
            dist2 = float(np.dot(delta, delta))
            if dist2 > best_dist2:
                best_dist2 = dist2
                best_pair = (h_a, h_b)
    return best_pair


def dihedral_h_o_o_h(
    atoms: list[Atom],
    h1: int,
    o1: int,
    o2: int,
    h2: int,
    box: np.ndarray | None,
) -> float:
    """Return the locally unwrapped H-O-O-H dihedral angle in radians."""
    p2 = atoms[o1].xyz
    p1 = p2 + minimum_image(atoms[h1].xyz - atoms[o1].xyz, box)
    p3 = p2 + minimum_image(atoms[o2].xyz - atoms[o1].xyz, box)
    p4 = p3 + minimum_image(atoms[h2].xyz - atoms[o2].xyz, box)
    return dihedral(p1, p2, p3, p4)


def dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Compute a signed four-point dihedral angle in radians."""
    b0 = p1 - p2
    b1 = p3 - p2
    b2 = p4 - p3
    b1_norm = float(np.linalg.norm(b1))
    if b1_norm <= 1e-12:
        return 0.0
    b1 = b1 / b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return atan2(y, x)


def mean_or_none(values: list[float]) -> float | None:
    """Return a float mean while preserving empty values as None."""
    if not values:
        return None
    return float(np.mean(values))


__all__ = [
    "compute_f3f4",
    "compute_order_parameters",
    "f3_for_water",
    "f4_for_water",
]
