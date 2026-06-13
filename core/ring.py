from __future__ import annotations

"""Ring search on the water-network graph."""

from collections import defaultdict

from ..models import Ring


def find_rings(adjacency: dict[int, set[int]], sizes: list[int], chordless: bool = True) -> dict[int, list[Ring]]:
    """Find primitive candidate rings with a bounded, non-recursive DFS."""
    if not sizes:
        return {}
    min_size = min(sizes)
    max_size = max(sizes)
    allowed = set(sizes)
    found: set[tuple[int, ...]] = set()
    nodes = sorted(adjacency)

    for start in nodes:
        # Store partial paths explicitly instead of relying on Python recursion.
        stack: list[list[int]] = [[start]]
        while stack:
            path = stack.pop()
            current = path[-1]
            if len(path) > max_size:
                continue
            for nb in sorted(adjacency.get(current, ())):
                if nb == start and len(path) in allowed:
                    if not chordless or is_chordless(path, adjacency):
                        found.add(canonical_cycle(path))
                    continue
                if len(path) >= max_size:
                    continue
                if nb in path:
                    continue
                # This ordering rule prevents rediscovering the same cycle from
                # every node while preserving all canonical cycles.
                if nb < start:
                    continue
                stack.append(path + [nb])

    by_size: dict[int, list[Ring]] = defaultdict(list)
    counts: dict[int, int] = defaultdict(int)
    for nodes_tuple in sorted(found, key=lambda item: (len(item), item)):
        # Stable object ids make downstream membership and GRO output readable.
        size = len(nodes_tuple)
        counts[size] += 1
        by_size[size].append(Ring(object_id=f"ring{size}_{counts[size]:05d}", nodes=nodes_tuple))
    return dict(by_size)


def canonical_cycle(path: list[int]) -> tuple[int, ...]:
    """Return the lexicographically smallest rotation over both directions."""
    n = len(path)
    rotations = [tuple(path[i:] + path[:i]) for i in range(n)]
    rev = list(reversed(path))
    rotations.extend(tuple(rev[i:] + rev[:i]) for i in range(n))
    return min(rotations)


def is_chordless(path: list[int], adjacency: dict[int, set[int]]) -> bool:
    """Reject cycles that contain an internal graph edge between non-neighbors."""
    n = len(path)
    for i, a in enumerate(path):
        for j in range(i + 1, n):
            if (j - i) % n in (1, n - 1):
                continue
            if path[j] in adjacency.get(a, set()):
                return False
    return True


def occupied_ring_ids_by_water(rings: dict[int, list[Ring]]) -> dict[int, set[str]]:
    """Index ring ids by oxygen node for later ownership filtering."""
    by_water: dict[int, set[str]] = defaultdict(set)
    for group in rings.values():
        for ring in group:
            for node in ring.nodes:
                by_water[node].add(ring.object_id)
    return dict(by_water)
