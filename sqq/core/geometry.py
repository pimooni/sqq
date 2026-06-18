from __future__ import annotations

"""Shared geometric helpers for local topology objects."""

import numpy as np

from ..models import Frame
from .pbc import minimum_image


def unwrap_connected_nodes(frame: Frame, nodes: list[int], edges: list[tuple[int, int]]) -> dict[int, np.ndarray]:
    """Unwrap a connected oxygen-node patch so it is not split by PBC."""
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
