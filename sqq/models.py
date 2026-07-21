from __future__ import annotations

"""Shared data models for one SQQ analysis run."""

from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np


# Parsed trajectory input.
@dataclass(frozen=True)
class Atom:
    index: int
    resid: int
    resname: str
    atomname: str
    atomid: int
    xyz: np.ndarray
    velocity: np.ndarray | None = None
    molecule_id: int | None = None


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


# Detected topology objects use oxygen-node ids.
@dataclass(frozen=True)
class Ring:
    object_id: str
    nodes: tuple[int, ...]

    @cached_property
    def size(self) -> int:
        return len(self.nodes)

    @cached_property
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
class HydrateMotif:
    object_id: str
    hydrate_type: str
    anchor_cage_ids: tuple[str, ...]
    cage_ids: tuple[str, ...]
    shared_face_ids: tuple[str, ...] = ()
    cluster_id: str = ""
    domain_id: str = ""
    status: str = "partial"
    completeness: float = 0.0
    consistency: float = 0.0
    confidence: float = 0.0
    classification_method: str = "SQQ local face topology"

    @property
    def cage_count(self) -> int:
        """Number of cages in this overlapping local topology motif."""
        return len(self.cage_ids)


@dataclass(frozen=True)
class HydrateDomain:
    object_id: str
    cluster_id: str
    hydrate_type: str
    cage_ids: tuple[str, ...]
    motif_ids: tuple[str, ...]
    waters: tuple[int, ...]
    guest_ids: tuple[str, ...]
    boundary_cage_ids: tuple[str, ...] = ()
    confidence: float = 0.0
    status: str = "growth"
    seed_count: int = 0
    seed_cage_ids: tuple[str, ...] = ()
    classified_fraction: float = 0.0

    @property
    def cage_count(self) -> int:
        """Number of uniquely phase-classified cages in this domain."""
        return len(self.cage_ids)

    @property
    def motif_count(self) -> int:
        """Number of phase-core motifs assigned to this domain."""
        return len(self.motif_ids)

    @property
    def water_count(self) -> int:
        """Number of unique water oxygens used by the domain cages."""
        return len(self.waters)

    @property
    def guest_count(self) -> int:
        """Number of unique guest molecules assigned to the domain cages."""
        return len(self.guest_ids)


@dataclass(frozen=True)
class HydrateCluster:
    object_id: str
    cage_ids: tuple[str, ...]
    cage_types: tuple[str, ...]
    waters: tuple[int, ...]
    guest_ids: tuple[str, ...]
    shared_faces: tuple[tuple[str, str, str], ...] = ()
    hydrate_type: str = "unclassified"
    motif_ids: tuple[str, ...] = ()
    domain_ids: tuple[str, ...] = ()
    classified_cage_ids: tuple[str, ...] = ()
    unclassified_cage_ids: tuple[str, ...] = ()
    ambiguous_cage_ids: tuple[str, ...] = ()
    boundary_cage_ids: tuple[str, ...] = ()
    hydrate_type_counts: tuple[tuple[str, int], ...] = ()

    @property
    def cage_count(self) -> int:
        """Number of cages in this connected cage cluster."""
        return len(self.cage_ids)

    @property
    def water_count(self) -> int:
        """Number of unique water oxygens used by the cluster cages."""
        return len(self.waters)

    @property
    def guest_count(self) -> int:
        """Number of unique guest molecules assigned to the cluster cages."""
        return len(self.guest_ids)

    @property
    def motif_count(self) -> int:
        """Number of phase-core motifs in this cluster."""
        return len(self.motif_ids)

    @property
    def domain_count(self) -> int:
        """Number of phase domains in this cluster."""
        return len(self.domain_ids)

    @property
    def boundary_cage_count(self) -> int:
        """Number of non-phase cages in the external boundary layer."""
        return len(self.boundary_cage_ids)


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


@dataclass(frozen=True)
class ClusterOrderValue:
    """Largest connected component for one hydrate nucleation order parameter."""

    largest_cluster_size: int | None
    members: tuple[int, ...] = ()
    eligible_count: int = 0
    member_type: str = ""


@dataclass(frozen=True)
class HydrateOrderResult:
    """MCG and DHOP values reported for one frame."""

    mcg1: ClusterOrderValue
    dhop35: ClusterOrderValue
    mcg3: ClusterOrderValue | None = None
    dhop30: ClusterOrderValue | None = None


# Frame analysis result for exporters.
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
    hydrate_cluster_enabled: bool = False
    hydrate_cluster_detail: bool = False
    hydrate_clusters: list[HydrateCluster] = field(default_factory=list)
    hydrate_motifs: list[HydrateMotif] = field(default_factory=list)
    hydrate_domains: list[HydrateDomain] = field(default_factory=list)
    isolated_cage_ids: tuple[str, ...] = ()
    f3f4: F3F4Result | None = None
    hydrate_order: HydrateOrderResult | None = None
    ice_like_waters: tuple[int, ...] = ()
    ice_i_waters: tuple[int, ...] = ()
    interfacial_ice_waters: tuple[int, ...] = ()
    warnings: list[str] = field(default_factory=list)
