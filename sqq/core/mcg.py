from __future__ import annotations

"""Mutually coordinated guest (MCG) hydrate order parameters."""

from collections import defaultdict
from typing import Any

import numpy as np

from .pbc import minimum_image
from .spatial import cross_cutoff_pairs, self_cutoff_pairs
from ..models import ClusterOrderValue, Frame, Guest, Water


def compute_mcg_order(
    frame: Frame,
    waters: list[Water],
    guests: list[Guest],
    config: dict[str, Any],
) -> tuple[ClusterOrderValue, ClusterOrderValue | None]:
    """Compute MCG-1 and, when requested, MCG-3 from one shared MCG graph."""
    mcg1_enabled = bool(config.get("mcg1_enabled", True))
    mcg3_enabled = bool(config.get("mcg3_enabled", False))
    selected_names = {str(value).upper() for value in config.get("mcg_guest_resnames", ("CH4", "MET"))}
    selected = [(index, guest) for index, guest in enumerate(guests) if guest.resname.upper() in selected_names]
    if not selected:
        unavailable = ClusterOrderValue(None, member_type="guest")
        return unavailable, ClusterOrderValue(None, member_type="guest") if mcg3_enabled else None
    if not mcg1_enabled and not mcg3_enabled:
        return ClusterOrderValue(None, member_type="guest"), None

    guest_coords = np.asarray([_guest_center(frame, guest) for _, guest in selected], dtype=float)
    water_coords = np.asarray([frame.atoms[water.oxygen].xyz for water in waters], dtype=float)
    guest_cutoff = float(config.get("mcg_guest_cutoff_nm", 0.90))
    water_cutoff = float(config.get("mcg_water_cutoff_nm", 0.60))
    half_angle = float(config.get("mcg_cone_half_angle_deg", 45.0))
    min_waters = int(config.get("mcg_min_waters", 5))
    if not 0 < half_angle < 90:
        raise ValueError("hydrate_order.mcg_cone_half_angle_deg must be between 0 and 90.")
    if min_waters < 1:
        raise ValueError("hydrate_order.mcg_min_waters must be at least 1.")

    nearby_waters: dict[int, set[int]] = defaultdict(set)
    for guest_index, water_index in cross_cutoff_pairs(guest_coords, water_coords, frame.box, water_cutoff):
        nearby_waters[guest_index].add(water_index)
    cosine_limit = float(np.cos(np.deg2rad(half_angle)))
    qualifying_edges: list[tuple[int, int]] = []
    for first, second in self_cutoff_pairs(guest_coords, frame.box, guest_cutoff):
        shared = nearby_waters.get(first, set()) & nearby_waters.get(second, set())
        if len(shared) < min_waters:
            continue
        axis = minimum_image(guest_coords[second] - guest_coords[first], frame.box)
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm <= 0:
            continue
        coordinated = 0
        for water_index in shared:
            from_first = minimum_image(water_coords[water_index] - guest_coords[first], frame.box)
            from_second = minimum_image(water_coords[water_index] - guest_coords[second], frame.box)
            norm_first = float(np.linalg.norm(from_first))
            norm_second = float(np.linalg.norm(from_second))
            if norm_first <= 0 or norm_second <= 0:
                continue
            first_cosine = float(np.dot(from_first, axis) / (norm_first * axis_norm))
            second_cosine = float(np.dot(from_second, -axis) / (norm_second * axis_norm))
            if first_cosine >= cosine_limit and second_cosine >= cosine_limit:
                coordinated += 1
        if coordinated >= min_waters:
            qualifying_edges.append((first, second))

    mcg1 = _cluster_value(qualifying_edges, selected, min_degree=1) if mcg1_enabled else ClusterOrderValue(None, member_type="guest")
    mcg3 = _cluster_value(qualifying_edges, selected, min_degree=3) if mcg3_enabled else None
    return mcg1, mcg3


def _guest_center(frame: Frame, guest: Guest) -> np.ndarray:
    if guest.center_atom is not None:
        return np.asarray(frame.atoms[guest.center_atom].xyz, dtype=float)
    if not guest.atoms:
        raise ValueError(f"Guest resid {guest.resid} has no atoms.")
    anchor = np.asarray(frame.atoms[guest.atoms[0]].xyz, dtype=float)
    points = [anchor]
    for atom_index in guest.atoms[1:]:
        delta = minimum_image(np.asarray(frame.atoms[atom_index].xyz, dtype=float) - anchor, frame.box)
        points.append(anchor + delta)
    return np.mean(np.asarray(points), axis=0)


def _cluster_value(
    edges: list[tuple[int, int]],
    selected: list[tuple[int, Guest]],
    min_degree: int,
) -> ClusterOrderValue:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    eligible = {node for node, neighbors in adjacency.items() if len(neighbors) >= min_degree}
    components = _components(eligible, adjacency)
    largest = min(components, key=lambda item: (-len(item), tuple(item))) if components else ()
    original_members = tuple(selected[index][0] for index in largest)
    return ClusterOrderValue(
        largest_cluster_size=len(largest),
        members=original_members,
        eligible_count=len(eligible),
        member_type="guest",
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
            neighbors = sorted(adjacency.get(node, set()) & remaining, reverse=True)
            for neighbor in neighbors:
                remaining.remove(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components
