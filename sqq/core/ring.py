from __future__ import annotations

"""Ring search on the water-network graph."""

from collections import defaultdict

from ..models import Ring


def find_rings(
    adjacency: dict[int, set[int]],
    sizes: list[int],
    chordless: bool = True,
    definition: str = "chordless",
) -> dict[int, list[Ring]]:
    """Find bounded cycles with incremental chord pruning and symmetry breaking."""
    if not sizes:
        return {}
    ring_definition = str(definition or "chordless").strip().lower()
    if ring_definition not in {"chordless", "shortest_path"}:
        raise ValueError("ring.definition must be chordless or shortest_path.")
    require_chordless = bool(chordless) or ring_definition == "shortest_path"
    max_size = max(sizes)
    allowed = set(sizes)
    found: set[tuple[int, ...]] = set()
    shortest_path_cache: dict[tuple[int, int], dict[int, int]] = {}
    nodes = sorted(adjacency)
    ordered_neighbors = {node: tuple(sorted(adjacency.get(node, ()))) for node in nodes}
    bit_for_node = {node: 1 << index for index, node in enumerate(nodes)}

    for start in nodes:
        start_bit = bit_for_node[start]
        stack: list[tuple[tuple[int, ...], int]] = [((start,), start_bit)]
        while stack:
            path, visited = stack.pop()
            current = path[-1]
            if len(path) >= max_size:
                continue
            for neighbor in reversed(ordered_neighbors.get(current, ())):
                if neighbor <= start or visited & bit_for_node.get(neighbor, 0):
                    continue
                if require_chordless:
                    earlier_neighbors = adjacency.get(neighbor, set()).intersection(path[:-1])
                    closes_cycle = len(path) >= 2 and start in earlier_neighbors
                    if earlier_neighbors - ({start} if closes_cycle else set()):
                        continue
                    new_path = (*path, neighbor)
                    if closes_cycle:
                        if len(new_path) in allowed and new_path[1] < new_path[-1]:
                            if ring_definition != "shortest_path" or is_shortest_path_ring(
                                new_path,
                                adjacency,
                                distance_cache=shortest_path_cache,
                            ):
                                found.add(canonical_cycle(list(new_path)))
                        continue
                    stack.append((new_path, visited | bit_for_node[neighbor]))
                    continue

                new_path = (*path, neighbor)
                if start in adjacency.get(neighbor, set()) and len(new_path) in allowed and new_path[1] < new_path[-1]:
                    found.add(canonical_cycle(list(new_path)))
                stack.append((new_path, visited | bit_for_node[neighbor]))

    by_size: dict[int, list[Ring]] = defaultdict(list)
    counts: dict[int, int] = defaultdict(int)
    for nodes_tuple in sorted(found, key=lambda item: (len(item), item)):
        size = len(nodes_tuple)
        counts[size] += 1
        by_size[size].append(Ring(object_id=f"ring{size}_{counts[size]:05d}", nodes=nodes_tuple))
    return dict(by_size)


def is_shortest_path_ring(
    path: tuple[int, ...] | list[int],
    adjacency: dict[int, set[int]],
    distance_cache: dict[tuple[int, int], dict[int, int]] | None = None,
) -> bool:
    """Apply Franzblau's all-pairs shortest-path criterion to one cycle."""
    cycle = tuple(path)
    size = len(cycle)
    max_required = size // 2
    for source_index, source in enumerate(cycle):
        cache_key = (source, max_required)
        distances = None if distance_cache is None else distance_cache.get(cache_key)
        if distances is None:
            distances = bounded_graph_distances(source, adjacency, max_required)
            if distance_cache is not None:
                distance_cache[cache_key] = distances
        for target_index in range(source_index + 1, size):
            along = target_index - source_index
            cycle_distance = min(along, size - along)
            if distances.get(cycle[target_index]) != cycle_distance:
                return False
    return True


def bounded_graph_distances(start: int, adjacency: dict[int, set[int]], limit: int) -> dict[int, int]:
    """Return BFS distances no farther than the largest distance needed by a ring."""
    distances = {start: 0}
    frontier = [start]
    for depth in range(1, limit + 1):
        next_frontier: list[int] = []
        for node in frontier:
            for neighbor in adjacency.get(node, set()):
                if neighbor in distances:
                    continue
                distances[neighbor] = depth
                next_frontier.append(neighbor)
        frontier = next_frontier
        if not frontier:
            break
    return distances

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
