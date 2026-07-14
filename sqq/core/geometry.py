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

    missing = sorted(set(nodes) - set(unwrapped))
    if missing:
        raise ValueError(
            "Cannot unwrap a disconnected topology object; "
            f"unreachable oxygen nodes: {missing[:10]}"
        )
    return unwrapped


def pbc_aware_centroid(frame: Frame, atom_indices: tuple[int, ...] | list[int]) -> np.ndarray:
    """Return a molecular centroid after unwrapping every atom near one anchor."""
    if not atom_indices:
        raise ValueError("Cannot compute a centroid for an empty atom selection.")
    anchor = np.asarray(frame.atoms[atom_indices[0]].xyz, dtype=float)
    points = [anchor]
    for atom_index in atom_indices[1:]:
        delta = minimum_image(
            np.asarray(frame.atoms[atom_index].xyz, dtype=float) - anchor,
            frame.box,
        )
        points.append(anchor + delta)
    return np.mean(np.asarray(points, dtype=float), axis=0)
