from __future__ import annotations

"""Ring, patch, cage, and graph data contracts."""

from dataclasses import dataclass
from functools import cached_property

import numpy as np


@dataclass(frozen=True)
class Ring:
    object_id: str
    nodes: tuple[int, ...]

    @cached_property
    def size(self) -> int:
        return len(self.nodes)

    @cached_property
    def edges(self) -> frozenset[tuple[int, int]]:
        pairs = []
        for index, left in enumerate(self.nodes):
            right = self.nodes[(index + 1) % len(self.nodes)]
            pairs.append((left, right) if left < right else (right, left))
        return frozenset(pairs)


@dataclass(frozen=True)
class CagePatch:
    object_id: str
    patch_type: str
    kind: str
    rings: tuple[str, ...]
    waters: tuple[int, ...]
    center: np.ndarray
    layers: tuple[str, ...] = ()


@dataclass(frozen=True)
class Cage:
    object_id: str
    cage_type: str
    rings: tuple[str, ...]
    waters: tuple[int, ...]
    center: np.ndarray
    guest_ids: tuple[str, ...] = ()
    isomer: str | None = None

    @property
    def occupied(self) -> bool:
        return bool(self.guest_ids)


@dataclass
class GraphResult:
    mode: str
    edges: list[tuple[int, int]]
    adjacency: dict[int, set[int]]


__all__ = ["Ring", "CagePatch", "Cage", "GraphResult"]
