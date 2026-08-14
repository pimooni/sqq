from __future__ import annotations

"""Hydrate motif, domain, and cluster data contracts."""

from dataclasses import dataclass


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
        return len(self.cage_ids)

    @property
    def motif_count(self) -> int:
        return len(self.motif_ids)

    @property
    def water_count(self) -> int:
        return len(self.waters)

    @property
    def guest_count(self) -> int:
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
        return len(self.cage_ids)

    @property
    def water_count(self) -> int:
        return len(self.waters)

    @property
    def guest_count(self) -> int:
        return len(self.guest_ids)

    @property
    def motif_count(self) -> int:
        return len(self.motif_ids)

    @property
    def domain_count(self) -> int:
        return len(self.domain_ids)

    @property
    def boundary_cage_count(self) -> int:
        return len(self.boundary_cage_ids)


__all__ = ["HydrateMotif", "HydrateDomain", "HydrateCluster"]
