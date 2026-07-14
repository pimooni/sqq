from __future__ import annotations

"""DHOP30/DHOP35 planar-water hydrate order parameters."""

from collections import defaultdict
from typing import Any

import numpy as np

from .pbc import minimum_image
from .spatial import self_cutoff_pairs
from ..models import ClusterOrderValue, Frame, Water


def compute_dhop_order(
    frame: Frame,
    waters: list[Water],
    config: dict[str, Any],
) -> tuple[ClusterOrderValue, ClusterOrderValue | None]:
    """Compute DHOP35/DHOP30 with batched plane-normal comparisons per O-O bond."""
    dhop35_enabled = bool(config.get("dhop35_enabled", True))
    dhop30_enabled = bool(config.get("dhop30_enabled", False))
    if not dhop35_enabled and not dhop30_enabled:
        return ClusterOrderValue(None, member_type="water"), None
    coords = np.asarray([frame.atoms[water.oxygen].xyz for water in waters], dtype=float)
    cutoff = float(config.get("dhop_neighbor_cutoff_nm", 0.35))
    pairs = self_cutoff_pairs(coords, frame.box, cutoff)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for first, second in pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)

    counts35 = np.zeros(len(waters), dtype=int) if dhop35_enabled else None
    counts30 = np.zeros(len(waters), dtype=int) if dhop30_enabled else None
    cosine35 = float(np.cos(np.deg2rad(35.0)))
    cosine30 = float(np.cos(np.deg2rad(30.0)))
    for center_first, center_second in pairs:
        axis = minimum_image(coords[center_second] - coords[center_first], frame.box)
        first_outer = sorted(adjacency[center_first] - {center_second})
        second_outer = sorted(adjacency[center_second] - {center_first})
        if not first_outer or not second_outer:
            continue
        first_vectors = minimum_image(coords[first_outer] - coords[center_first], frame.box)
        second_vectors = minimum_image(coords[second_outer] - coords[center_second], frame.box)
        first_normals = np.cross(first_vectors, axis)
        second_normals = np.cross(-axis, second_vectors)
        events35, events30 = planar_pair_event_counts(
            first_normals,
            second_normals,
            cosine35 if counts35 is not None else None,
            cosine30 if counts30 is not None else None,
        )
        if counts35 is not None:
            counts35[center_first] += events35
            counts35[center_second] += events35
        if counts30 is not None:
            counts30[center_first] += events30
            counts30[center_second] += events30

    planar_counts = {int(value) for value in config.get("dhop_planar_counts", (11, 12))}
    min_neighbors = int(config.get("dhop_min_qualified_neighbors", 3))
    if not planar_counts or min(planar_counts) < 0:
        raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.")
    if min_neighbors < 1:
        raise ValueError("hydrate_order.dhop_min_qualified_neighbors must be at least 1.")
    dhop35 = _cluster_value(counts35, planar_counts, min_neighbors, adjacency, waters) if counts35 is not None else ClusterOrderValue(None, member_type="water")
    dhop30 = _cluster_value(counts30, planar_counts, min_neighbors, adjacency, waters) if counts30 is not None else None
    return dhop35, dhop30


def planar_pair_event_counts(
    first_normals: np.ndarray,
    second_normals: np.ndarray,
    cosine35: float | None,
    cosine30: float | None,
) -> tuple[int, int]:
    """Count qualifying normal pairs while retaining scalar behavior at thresholds."""
    first_norms = np.linalg.norm(first_normals, axis=1)
    second_norms = np.linalg.norm(second_normals, axis=1)
    first_valid = first_norms > 1.0e-14
    second_valid = second_norms > 1.0e-14
    if not np.any(first_valid) or not np.any(second_valid):
        return 0, 0
    left = first_normals[first_valid]
    right = second_normals[second_valid]
    denominators = np.outer(first_norms[first_valid], second_norms[second_valid])
    cosines = (left @ right.T) / denominators
    np.clip(cosines, -1.0, 1.0, out=cosines)
    thresholds = [value for value in (cosine35, cosine30) if value is not None]
    if thresholds:
        nearest = np.minimum.reduce([np.abs(cosines - value) for value in thresholds])
        for row, column in zip(*np.where(nearest <= 1.0e-12), strict=True):
            scalar = float(np.dot(left[row], right[column]) / denominators[row, column])
            cosines[row, column] = min(1.0, max(-1.0, scalar))
    events35 = int(np.count_nonzero(cosines >= cosine35)) if cosine35 is not None else 0
    events30 = int(np.count_nonzero(cosines >= cosine30)) if cosine30 is not None else 0
    return events35, events30


def _cluster_value(
    planar_events: np.ndarray,
    planar_counts: set[int],
    min_neighbors: int,
    adjacency: dict[int, set[int]],
    waters: list[Water],
) -> ClusterOrderValue:
    qualified = {index for index, count in enumerate(planar_events) if int(count) in planar_counts}
    seeds = {
        index
        for index in qualified
        if len(adjacency.get(index, set()) & qualified) >= min_neighbors
    }
    tagged = set(seeds)
    for index in seeds:
        tagged.update(adjacency.get(index, set()))
    components = _components(tagged, adjacency)
    largest = min(components, key=lambda item: (-len(item), tuple(item))) if components else ()
    members = tuple(waters[index].oxygen for index in largest)
    return ClusterOrderValue(
        largest_cluster_size=len(largest),
        members=members,
        eligible_count=len(tagged),
        member_type="water",
    )


def _components(nodes: set[int], adjacency: dict[int, set[int]]) -> list[tuple[int, ...]]:
    remaining = set(nodes)
    components: list[tuple[int, ...]] = []
    while remaining:
        root = min(remaining)
        stack = [root]
        remaining.remove(root)
        component: list[int] = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in sorted(adjacency.get(node, set()) & remaining, reverse=True):
                remaining.remove(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components
