"""Deterministic connected-component helpers shared by scientific modules."""

from __future__ import annotations

from collections.abc import Hashable, Mapping, Set
from typing import TypeVar


NodeT = TypeVar("NodeT", bound=Hashable)


def connected_components(
    nodes: Set[NodeT],
    adjacency: Mapping[NodeT, Set[NodeT]],
) -> list[tuple[NodeT, ...]]:
    """Return deterministic components of the graph induced by ``nodes``."""
    remaining = set(nodes)
    components: list[tuple[NodeT, ...]] = []
    while remaining:
        root = min(remaining)
        stack = [root]
        remaining.remove(root)
        component: list[NodeT] = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in sorted(adjacency.get(node, set()) & remaining, reverse=True):
                remaining.remove(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


__all__ = ["connected_components"]
