from __future__ import annotations

"""F3/F4 AOP order parameters from the user's reference script."""

from functools import lru_cache
from itertools import combinations
from math import atan2, cos, factorial, pi, radians, sin, sqrt
from typing import Any

import numpy as np

from ..models import Atom, F3F4Result, Frame, GraphResult, Water, WaterOrder
from .pbc import minimum_image


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
    f3f4_enabled: bool = True,
    q_enabled: bool = True,
    q_neighbor_mode: str = "graph",
    q_cutoff_nm: float = 0.35,
    q_n_neighbor: int | None = None,
    q_degree: Any = None,
    focus_resids: set[int] | None = None,
) -> F3F4Result:
    """Compute per-water F3/F4 and requested Steinhardt Q_l values."""
    water_by_oxygen = {water.oxygen: water for water in waters}
    focus_resids = set(focus_resids or set())
    rows: list[WaterOrder] = []
    f3_values: list[float] = []
    f4_values: list[float] = []
    f3_focus_values: list[float] = []
    f4_focus_values: list[float] = []
    q_degree_resolved = normalize_q_degree(q_degree)
    if not q_enabled:
        q_degree_resolved = ()
    q_value_lists: dict[int, list[float]] = {degree: [] for degree in q_degree_resolved}
    q_focus_value_lists: dict[int, list[float]] = {degree: [] for degree in q_degree_resolved}
    q_mode = normalize_q_neighbor_mode(q_neighbor_mode)
    q_fixed_neighbor = resolve_q_neighbor_count(q_mode, q_n_neighbor)
    graph_vectors = (
        build_graph_neighbor_vector_cache(frame, waters, graph)
        if f3f4_enabled or (q_degree_resolved and q_mode == "graph")
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
        ) if f3f4_enabled else None
        f4 = f4_for_water(frame.atoms, water, [water_by_oxygen[idx] for idx in neighbors], frame.box) if f3f4_enabled else None
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


def f4_for_water(atoms: list[Atom], water: Water, neighbors: list[Water], box: np.ndarray | None) -> float | None:
    """F4: average cos(3*dihedral) after selecting the farthest H-H pair."""
    if len(water.hydrogens) < 2:
        return None
    terms = []
    for neighbor in neighbors:
        if len(neighbor.hydrogens) < 2:
            continue
        h1, h2 = farthest_hydrogen_pair(atoms, water, neighbor, box)
        angle = dihedral_h_o_o_h(atoms, h1, water.oxygen, neighbor.oxygen, h2, box)
        terms.append(cos(3.0 * angle))
    return mean_or_none(terms)


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


def normalize_q_degree(value: Any = None) -> tuple[int, ...]:
    """Resolve requested Steinhardt Q_l degree list; default reports Q6 and Q12."""
    if value in (None, "", "auto"):
        raw_items: list[Any] = [6, 12]
    elif isinstance(value, str):
        raw_items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            raw_items = list(value)
        except TypeError as exc:
            raise ValueError("order.q_degree must be a comma-separated list or a list of integers.") from exc
    degrees = sorted({int(item) for item in raw_items})
    if not degrees:
        raise ValueError("order.q_degree must contain at least one degree.")
    if any(degree < 0 for degree in degrees):
        raise ValueError("order.q_degree must contain non-negative integers.")
    return tuple(degrees)


def normalize_q_neighbor_mode(mode: str) -> str:
    """Normalize Q_l neighbor source names."""
    value = str(mode or "graph").strip().lower().replace("-", "_")
    aliases = {"water": "graph", "oo": "cutoff", "cutoff_all": "cutoff"}
    value = aliases.get(value, value)
    if value not in {"graph", "cutoff", "nearest", "lammps"}:
        raise ValueError("order.q_neighbor_mode must be graph, cutoff, nearest, or lammps.")
    return value


def resolve_q_neighbor_count(mode: str, value: int | str | None) -> int | None:
    """Resolve fixed-neighbor Q_l settings; lammps mode defaults to nnn=12."""
    if value in (None, "", "auto", "none", "None", "NULL", "null"):
        return 12 if mode == "lammps" else None
    try:
        count = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("order.q_n_neighbor must be an integer or null.") from exc
    if count < 1:
        raise ValueError("order.q_n_neighbor must be at least 1 when set.")
    return count


def q_values_for_water(
    frame: Frame,
    waters: list[Water],
    graph: GraphResult,
    water: Water,
    *,
    mode: str,
    cutoff_nm: float,
    n_neighbor: int | None,
    degree: tuple[int, ...],
) -> tuple[dict[int, float | None], int]:
    """Compute requested Q_l values for one water oxygen from one neighbor list."""
    candidates = q_candidate_vectors(frame, waters, graph, water, mode=mode, cutoff_nm=cutoff_nm)
    return q_values_from_candidates(candidates, n_neighbor=n_neighbor, degree=degree)


def q_values_from_candidates(
    candidates: list[tuple[float, np.ndarray]] | tuple[tuple[float, np.ndarray], ...],
    *,
    n_neighbor: int | None,
    degree: tuple[int, ...],
) -> tuple[dict[int, float | None], int]:
    """Compute Q_l from one already ordered per-water candidate list."""
    if n_neighbor is not None:
        if len(candidates) < n_neighbor:
            # LAMMPS sets Q_l to zero when fewer than nnn neighbors are available.
            return {item: 0.0 for item in degree}, len(candidates)
        candidates = candidates[:n_neighbor]
    if not candidates:
        return {item: None for item in degree}, 0
    vectors = [vector for _, vector in candidates]
    return q_values_from_vectors(vectors, degree), len(vectors)


def build_q_candidate_cache(
    frame: Frame,
    waters: list[Water],
    graph: GraphResult,
    *,
    mode: str,
    cutoff_nm: float,
    graph_vectors: dict[int, dict[int, np.ndarray]] | None = None,
) -> dict[int, list[tuple[float, np.ndarray]]]:
    """Build every ordered Q_l neighbor list once for the current frame."""
    ordered_oxygens = [water.oxygen for water in waters]
    oxygen_rank = {oxygen: index for index, oxygen in enumerate(ordered_oxygens)}
    candidates: dict[int, list[tuple[float, np.ndarray]]] = {
        oxygen: [] for oxygen in ordered_oxygens
    }

    if mode == "graph":
        resolved_graph_vectors = graph_vectors or build_graph_neighbor_vector_cache(frame, waters, graph)
        for oxygen in ordered_oxygens:
            neighbors = graph.adjacency.get(oxygen, set())
            ordered_neighbors = sorted(
                (other for other in neighbors if other in oxygen_rank),
                key=oxygen_rank.__getitem__,
            )
            for other_oxygen in ordered_neighbors:
                vector = resolved_graph_vectors[oxygen][other_oxygen]
                pair = q_vector_candidate(vector)
                if pair is not None:
                    candidates[oxygen].append(pair)
    else:
        for left_index, left in enumerate(waters):
            left_center = frame.atoms[left.oxygen].xyz
            for right in waters[left_index + 1 :]:
                vector = minimum_image(frame.atoms[right.oxygen].xyz - left_center, frame.box)
                distance = float(np.linalg.norm(vector))
                if distance <= 1e-12 or distance > cutoff_nm:
                    continue
                candidates[left.oxygen].append((distance, vector))
                candidates[right.oxygen].append((distance, -vector))

    for oxygen in candidates:
        candidates[oxygen].sort(key=lambda item: item[0])
    return candidates


def build_graph_neighbor_vector_cache(
    frame: Frame,
    waters: list[Water],
    graph: GraphResult,
) -> dict[int, dict[int, np.ndarray]]:
    """Compute each undirected graph-edge vector once and expose both directions."""
    oxygen_set = {water.oxygen for water in waters}
    vectors: dict[int, dict[int, np.ndarray]] = {oxygen: {} for oxygen in oxygen_set}
    for left in sorted(oxygen_set):
        for right in sorted(graph.adjacency.get(left, set())):
            if right not in oxygen_set or right <= left:
                continue
            vector = minimum_image(frame.atoms[right].xyz - frame.atoms[left].xyz, frame.box)
            vectors[left][right] = vector
            vectors[right][left] = -vector
    return vectors


def q_vector_candidate(vector: np.ndarray) -> tuple[float, np.ndarray] | None:
    """Return a reusable nonzero Q_l vector and its distance."""
    distance = float(np.linalg.norm(vector))
    if distance <= 1e-12:
        return None
    return distance, vector


def q_candidate_vectors(
    frame: Frame,
    waters: list[Water],
    graph: GraphResult,
    water: Water,
    *,
    mode: str,
    cutoff_nm: float,
) -> list[tuple[float, np.ndarray]]:
    """Return distance-sorted oxygen-neighbor vectors for Q_l."""
    center = frame.atoms[water.oxygen].xyz
    candidates: list[tuple[float, np.ndarray]] = []
    if mode == "graph":
        neighbor_oxygens = sorted(graph.adjacency.get(water.oxygen, set()))
        allowed = set(neighbor_oxygens)
        source = [other for other in waters if other.oxygen in allowed]
    else:
        source = [other for other in waters if other.oxygen != water.oxygen]
    for other in source:
        if other.oxygen == water.oxygen:
            continue
        vector = minimum_image(frame.atoms[other.oxygen].xyz - center, frame.box)
        distance = float(np.linalg.norm(vector))
        if distance <= 1e-12:
            continue
        if mode in {"cutoff", "nearest", "lammps"} and distance > cutoff_nm:
            continue
        candidates.append((distance, vector))
    return sorted(candidates, key=lambda item: item[0])


def q_l_from_vectors(vectors: list[np.ndarray], degree: int = 6) -> float:
    """LAMMPS/Steinhardt Q_l from unweighted bond vectors."""
    return q_values_from_vectors(vectors, (degree,))[degree]


def q_values_from_vectors(vectors: list[np.ndarray], degrees: tuple[int, ...]) -> dict[int, float]:
    """Compute multiple Q_l degrees while reusing vector angles and constants."""
    if not vectors:
        return {degree: 0.0 for degree in degrees}
    totals = {degree: [0j for _ in range(2 * degree + 1)] for degree in degrees}
    for vector in vectors:
        norm = float(np.linalg.norm(vector))
        if norm <= 1e-12:
            continue
        x, y, z = (float(value) / norm for value in vector)
        theta_cos = max(-1.0, min(1.0, z))
        phi = atan2(y, x)
        for degree in degrees:
            degree_totals = totals[degree]
            for order in range(0, degree + 1):
                positive = spherical_harmonic_from_angles(degree, order, theta_cos, phi)
                degree_totals[degree + order] += positive
                if order:
                    degree_totals[degree - order] += ((-1) ** order) * positive.conjugate()
    divisor = len(vectors)
    values: dict[int, float] = {}
    for degree in degrees:
        total_norm = sum(abs(value / divisor) ** 2 for value in totals[degree])
        values[degree] = float(sqrt(4.0 * pi / (2 * degree + 1) * total_norm))
    return values


def spherical_harmonic(degree: int, order: int, vector: np.ndarray) -> complex:
    """Return Y_lm(theta, phi) for the direction of a Cartesian vector."""
    norm = float(np.linalg.norm(vector))
    if norm <= 1e-12:
        return 0j
    x, y, z = (float(value) / norm for value in vector)
    theta_cos = max(-1.0, min(1.0, z))
    phi = atan2(y, x)
    if order < 0:
        positive = spherical_harmonic_from_angles(degree, -order, theta_cos, phi)
        return ((-1) ** order) * positive.conjugate()
    return spherical_harmonic_from_angles(degree, order, theta_cos, phi)


def spherical_harmonic_from_angles(degree: int, order: int, theta_cos: float, phi: float) -> complex:
    """Return positive-order Y_lm from one normalized direction."""
    legendre = associated_legendre(degree, order, theta_cos)
    norm_factor = spherical_harmonic_normalization(degree, order)
    phase = complex(cos(order * phi), sin(order * phi))
    return norm_factor * legendre * phase


@lru_cache(maxsize=None)
def spherical_harmonic_normalization(degree: int, order: int) -> float:
    """Cache the degree/order normalization shared by every water and frame."""
    return sqrt((2 * degree + 1) / (4 * pi) * factorial(degree - order) / factorial(degree + order))


def associated_legendre(degree: int, order: int, x: float) -> float:
    """Associated Legendre polynomial with the Condon-Shortley phase."""
    if order < 0 or order > degree:
        return 0.0
    pmm = 1.0
    if order > 0:
        somx2 = sqrt(max(0.0, (1.0 - x) * (1.0 + x)))
        fact = 1.0
        for _ in range(1, order + 1):
            pmm *= -fact * somx2
            fact += 2.0
    if degree == order:
        return pmm
    pmmp1 = x * (2 * order + 1) * pmm
    if degree == order + 1:
        return pmmp1
    pll = 0.0
    p_l_minus_two = pmm
    p_l_minus_one = pmmp1
    for ell in range(order + 2, degree + 1):
        pll = ((2 * ell - 1) * x * p_l_minus_one - (ell + order - 1) * p_l_minus_two) / (ell - order)
        p_l_minus_two = p_l_minus_one
        p_l_minus_one = pll
    return pll
