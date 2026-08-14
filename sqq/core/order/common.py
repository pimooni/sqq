from __future__ import annotations

from collections.abc import Mapping, Set


def connected_components(
    nodes: Set[int],
    adjacency: Mapping[int, Set[int]],
) -> list[tuple[int, ...]]:
    """Return deterministic connected components induced by ``nodes``."""
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
