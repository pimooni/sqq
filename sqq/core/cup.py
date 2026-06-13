from __future__ import annotations

"""General cup search from base rings and side-ring closure."""

from collections import defaultdict
from itertools import product

import numpy as np

from ..core.pbc import minimum_image
from ..models import Cup, Frame, Ring


def find_cups(
    frame: Frame,
    rings: dict[int, list[Ring]],
    enabled: bool = False,
    base_sizes: list[int] | None = None,
    side_sizes: list[int] | None = None,
    max_combinations_per_base: int = 50000,
) -> list[Cup]:
    """Find cups where each base edge grows one side ring and the side wall closes."""
    if not enabled:
        return []

    base_allowed = set(base_sizes or [5, 6])
    side_allowed = set(side_sizes or [5, 6])
    all_rings = [ring for group in rings.values() for ring in group]
    ring_by_id = {ring.object_id: ring for ring in all_rings}
    edge_to_rings = build_edge_to_rings(all_rings)
    cups: list[Cup] = []
    seen: set[frozenset[str]] = set()
    type_counts: dict[str, int] = defaultdict(int)

    for base in all_rings:
        if base.size not in base_allowed:
            continue
        # A cup is seeded by a base ring. Each base edge must have one side
        # ring sharing exactly that base edge with the base.
        base_edges = ordered_edges(base.nodes)
        candidate_lists: list[list[Ring]] = []
        for edge in base_edges:
            candidates = [
                ring
                for ring in edge_to_rings.get(edge, [])
                if ring.object_id != base.object_id
                and ring.size in side_allowed
                and (base.edges & ring.edges) == {edge}
            ]
            if not candidates:
                candidate_lists = []
                break
            candidate_lists.append(candidates)
        if not candidate_lists:
            continue

        # A malformed local network can create a large Cartesian product; skip
        # those rare seeds and let the summary remain conservative.
        combination_count = 1
        for candidates in candidate_lists:
            combination_count *= len(candidates)
            if combination_count > max_combinations_per_base:
                break
        if combination_count > max_combinations_per_base:
            continue

        for side_tuple in product(*candidate_lists):
            side_rings = list(side_tuple)
            side_ids = [ring.object_id for ring in side_rings]
            if len(set(side_ids)) != base.size:
                continue
            cup_key = frozenset([base.object_id, *side_ids])
            if cup_key in seen:
                continue
            # Adjacent side rings must share a non-base edge, forming the cup wall.
            if not side_wall_closed(base, side_rings):
                continue
            unique_nodes = tuple(sorted({node for ring in [base, *side_rings] for node in ring.nodes}))
            # Exact water count catches overlapped or shifted side-ring choices.
            expected_nodes = sum(ring.size for ring in side_rings) - 2 * base.size
            if len(unique_nodes) != expected_nodes:
                continue
            sequence = canonical_sequence([ring.size for ring in side_rings])
            cup_type = f"cup{base.size}_{sequence}"
            type_counts[cup_type] += 1
            ring_ids = tuple([base.object_id, *side_ids])
            center = cup_center(frame, [ring_by_id[ring_id] for ring_id in ring_ids])
            cups.append(
                Cup(
                    object_id=f"{cup_type}_{type_counts[cup_type]:05d}",
                    cup_type=cup_type,
                    rings=ring_ids,
                    waters=unique_nodes,
                    center=center,
                )
            )
            seen.add(cup_key)
    return cups


def build_edge_to_rings(rings: list[Ring]) -> dict[tuple[int, int], list[Ring]]:
    """Map each graph edge to the rings that use it."""
    edge_to_rings: dict[tuple[int, int], list[Ring]] = defaultdict(list)
    for ring in rings:
        for edge in ring.edges:
            edge_to_rings[edge].append(ring)
    return dict(edge_to_rings)


def ordered_edges(nodes: tuple[int, ...]) -> list[tuple[int, int]]:
    """Return ring edges in cyclic order with sorted edge endpoints."""
    edges = []
    for idx, a in enumerate(nodes):
        b = nodes[(idx + 1) % len(nodes)]
        edges.append((a, b) if a < b else (b, a))
    return edges


def side_wall_closed(base: Ring, side_rings: list[Ring]) -> bool:
    """Check that neighboring side rings share exactly one wall edge."""
    base_edges = base.edges
    n = len(side_rings)
    for idx in range(n):
        current_ring = side_rings[idx]
        next_ring = side_rings[(idx + 1) % n]
        shared_nonbase = (current_ring.edges & next_ring.edges) - base_edges
        if len(shared_nonbase) != 1:
            return False
    return True


def canonical_sequence(values: list[int]) -> str:
    """Canonicalize side-ring sizes while preserving topological isomers."""
    n = len(values)
    rotations = [values[i:] + values[:i] for i in range(n)]
    rev = list(reversed(values))
    rotations.extend(rev[i:] + rev[:i] for i in range(n))
    return "".join(str(value) for value in min(rotations))


def cup_center(frame: Frame, rings: list[Ring]) -> np.ndarray:
    """Use the centroid of unique, locally unwrapped oxygen nodes."""
    nodes = sorted({node for ring in rings for node in ring.nodes})
    edges = sorted({edge for ring in rings for edge in ring.edges})
    unwrapped = unwrap_connected_nodes(frame, nodes, edges)
    return np.mean([unwrapped[node] for node in nodes], axis=0)


def unwrap_connected_nodes(frame: Frame, nodes: list[int], edges: list[tuple[int, int]]) -> dict[int, np.ndarray]:
    """Unwrap a connected face patch so its centroid is not split by PBC."""
    adjacency: dict[int, set[int]] = {node: set() for node in nodes}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    start = nodes[0]
    unwrapped = {start: frame.atoms[start].xyz}
    stack = [start]
    while stack:
        current = stack.pop()
        for nb in adjacency.get(current, set()):
            if nb in unwrapped:
                continue
            delta = minimum_image(frame.atoms[nb].xyz - frame.atoms[current].xyz, frame.box)
            unwrapped[nb] = unwrapped[current] + delta
            stack.append(nb)
    for node in nodes:
        unwrapped.setdefault(node, frame.atoms[node].xyz)
    return unwrapped
