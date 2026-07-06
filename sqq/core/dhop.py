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
    """Compute DHOP35 and, when requested, DHOP30 using a dedicated O--O graph."""
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

    counts35 = np.zeros(len(waters), dtype=int)
    counts30 = np.zeros(len(waters), dtype=int) if dhop30_enabled else None
    cosine35 = float(np.cos(np.deg2rad(35.0)))
    cosine30 = float(np.cos(np.deg2rad(30.0)))
    for center_first, center_second in pairs:
        axis = minimum_image(coords[center_second] - coords[center_first], frame.box)
        for first_outer in sorted(adjacency[center_first] - {center_second}):
            first_vector = minimum_image(coords[first_outer] - coords[center_first], frame.box)
            first_normal = np.cross(first_vector, axis)
            first_norm = float(np.linalg.norm(first_normal))
            if first_norm <= 1.0e-14:
                continue
            for second_outer in sorted(adjacency[center_second] - {center_first}):
                second_vector = minimum_image(coords[second_outer] - coords[center_second], frame.box)
                second_normal = np.cross(-axis, second_vector)
                second_norm = float(np.linalg.norm(second_normal))
                if second_norm <= 1.0e-14:
                    continue
                cosine = float(np.dot(first_normal, second_normal) / (first_norm * second_norm))
                cosine = min(1.0, max(-1.0, cosine))
                if cosine >= cosine35:
                    counts35[center_first] += 1
                    counts35[center_second] += 1
                if counts30 is not None and cosine >= cosine30:
                    counts30[center_first] += 1
                    counts30[center_second] += 1

    planar_counts = {int(value) for value in config.get("dhop_planar_counts", (11, 12))}
    min_neighbors = int(config.get("dhop_min_qualified_neighbors", 3))
    if not planar_counts or min(planar_counts) < 0:
        raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.")
    if min_neighbors < 1:
        raise ValueError("hydrate_order.dhop_min_qualified_neighbors must be at least 1.")
    dhop35 = _cluster_value(counts35, planar_counts, min_neighbors, adjacency, waters) if dhop35_enabled else ClusterOrderValue(None, member_type="water")
    dhop30 = _cluster_value(counts30, planar_counts, min_neighbors, adjacency, waters) if counts30 is not None else None
    return dhop35, dhop30


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
