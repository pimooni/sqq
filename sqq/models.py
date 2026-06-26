from __future__ import annotations

"""Shared data models for one SQQ analysis run."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


# Atom, Frame, Water, and Guest describe parsed trajectory input.
@dataclass(frozen=True)
class Atom:
    index: int
    resid: int
    resname: str
    atomname: str
    atomid: int
    xyz: np.ndarray


@dataclass
class Frame:
    name: str
    atoms: list[Atom]
    box: np.ndarray | None = None
    time_ps: float | None = None
    source: Path | None = None


@dataclass(frozen=True)
class Water:
    resid: int
    resname: str
    oxygen: int
    hydrogens: tuple[int, ...]
    atoms: tuple[int, ...]


@dataclass(frozen=True)
class Guest:
    resid: int
    resname: str
    atoms: tuple[int, ...]
    center_atom: int | None = None


# Ring, CagePatch, and Cage describe detected topology objects by oxygen-node ids.
@dataclass(frozen=True)
class Ring:
    object_id: str
    nodes: tuple[int, ...]

    @property
    def size(self) -> int:
        return len(self.nodes)

    @property
    def edges(self) -> frozenset[tuple[int, int]]:
        """Return undirected graph edges around the ring."""
        pairs = []
        for i, a in enumerate(self.nodes):
            b = self.nodes[(i + 1) % len(self.nodes)]
            pairs.append((a, b) if a < b else (b, a))
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
        """Whether at least one guest is assigned to the cage."""
        return bool(self.guest_ids)


@dataclass(frozen=True)
class WaterOrder:
    oxygen: int
    resid: int
    atomid: int
    xyz: np.ndarray
    f3: float | None
    f4: float | None
    q_values: dict[int, float | None] = field(default_factory=dict)
    q_neighbors: int = 0


@dataclass(frozen=True)
class F3F4Result:
    per_water: tuple[WaterOrder, ...]
    f3_mean: float | None
    f4_mean: float | None
    f3_valid: int
    f4_valid: int
    focus_resids: tuple[int, ...] = ()
    f3_focus_mean: float | None = None
    f4_focus_mean: float | None = None
    f3_focus_valid: int = 0
    f4_focus_valid: int = 0
    q_degree: tuple[int, ...] = ()
    q_means: dict[int, float | None] = field(default_factory=dict)
    q_valid_counts: dict[int, int] = field(default_factory=dict)
    q_focus_means: dict[int, float | None] = field(default_factory=dict)
    q_focus_valid_counts: dict[int, int] = field(default_factory=dict)
    q_neighbor_mode: str = ""
    q_cutoff_nm: float | None = None
    q_n_neighbor: int | None = None


# FrameResult carries both raw selections and derived topology for exporters.
@dataclass
class GraphResult:
    mode: str
    edges: list[tuple[int, int]]
    adjacency: dict[int, set[int]]


@dataclass
class FrameResult:
    frame: Frame
    waters: list[Water]
    guests: list[Guest]
    graph: GraphResult
    rings: dict[int, list[Ring]]
    ring_report_sizes: tuple[int, ...] = ()
    half_cages: list[CagePatch] = field(default_factory=list)
    quasi_cages: list[CagePatch] = field(default_factory=list)
    cages: list[Cage] = field(default_factory=list)
    all_cages: list[Cage] = field(default_factory=list)
    cage_report_types: tuple[str, ...] | None = None
    f3f4: F3F4Result | None = None
    ice_like_waters: tuple[int, ...] = ()
    ice_i_waters: tuple[int, ...] = ()
    interfacial_ice_waters: tuple[int, ...] = ()
    warnings: list[str] = field(default_factory=list)



