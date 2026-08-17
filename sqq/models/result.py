"""Complete per-frame analysis result contract."""

from __future__ import annotations

from dataclasses import dataclass, field

from .order import F3F4Result, HydrateOrderResult
from .phase import HydrateCluster, HydrateDomain, HydrateMotif
from .structure import Frame, Guest, Water
from .topology import Cage, CagePatch, GraphResult, Ring


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


__all__ = ["FrameResult"]
