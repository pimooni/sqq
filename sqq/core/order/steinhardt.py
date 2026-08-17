"""Steinhardt bond-orientational order parameters."""

from __future__ import annotations

from collections.abc import Iterable
from functools import lru_cache
from math import atan2, cos, factorial, pi, sin, sqrt
from typing import Any

import numpy as np

from ...models import Frame, GraphResult, Water
from ..pbc import minimum_image
from ..spatial import self_cutoff_pairs


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


def fixed_neighbor_shortfall_warning(
    neighbor_counts: Iterable[int],
    *,
    mode: str,
    n_neighbor: int | None,
    cutoff_nm: float | None,
    threshold_fraction: float = 0.5,
) -> str | None:
    """Summarize widespread fixed-neighbor shortfalls without changing Q_l."""
    if mode not in {"nearest", "lammps"} or n_neighbor is None:
        return None
    counts = tuple(int(value) for value in neighbor_counts)
    if not counts:
        return None
    shortfall = sum(value < n_neighbor for value in counts)
    if shortfall == 0 or shortfall / len(counts) < threshold_fraction:
        return None
    cutoff_text = "the configured cutoff" if cutoff_nm is None else f"{cutoff_nm:g} nm"
    return (
        f"Q_l {mode} mode found fewer than {n_neighbor} neighbors within "
        f"{cutoff_text} for {shortfall}/{len(counts)} waters; their Q_l values "
        "follow the configured fixed-neighbor rule and are reported as 0."
    )


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
        coordinates = np.asarray(
            [frame.atoms[water.oxygen].xyz for water in waters],
            dtype=float,
        )
        for left_index, right_index in self_cutoff_pairs(
            coordinates,
            frame.box,
            cutoff_nm,
        ):
            left = waters[left_index]
            right = waters[right_index]
            vector = minimum_image(
                coordinates[right_index] - coordinates[left_index],
                frame.box,
            )
            distance = float(np.linalg.norm(vector))
            if distance <= 1e-12:
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
    valid_vector_count = 0
    for vector in vectors:
        norm = float(np.linalg.norm(vector))
        if norm <= 1e-12:
            continue
        valid_vector_count += 1
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
    if valid_vector_count == 0:
        return {degree: 0.0 for degree in degrees}
    divisor = valid_vector_count
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


__all__ = [
    "build_graph_neighbor_vector_cache",
    "build_q_candidate_cache",
    "fixed_neighbor_shortfall_warning",
    "normalize_q_degree",
    "normalize_q_neighbor_mode",
    "q_l_from_vectors",
    "q_values_for_water",
    "q_values_from_candidates",
    "q_values_from_vectors",
    "resolve_q_neighbor_count",
]
