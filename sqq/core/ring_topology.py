from __future__ import annotations

"""Frame-local indexes and geometry shared by ring-face analyses."""

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from ..models import Frame, Ring
from .geometry import unwrap_connected_nodes


Edge = tuple[int, int]


@dataclass(frozen=True)
class RingFaceQuality:
    """PBC-aware geometry diagnostics for one ordered ring face."""

    planarity_rms_nm: float
    edge_length_cv: float
    projected_area_nm2: float


@dataclass
class RingTopologyIndex:
    """Reusable topology and geometry for all rings in one frame."""

    ring_by_id: dict[str, Ring]
    ring_centers: dict[str, np.ndarray]
    edge_to_ring_ids: dict[Edge, tuple[str, ...]]
    ring_adjacency: dict[str, frozenset[str]]
    ring_normals: dict[str, np.ndarray] = field(default_factory=dict)
    face_quality: dict[str, RingFaceQuality] = field(default_factory=dict)
    distance_cache: dict[tuple[str, str], float] = field(default_factory=dict)

    def edge_to_rings(self) -> dict[Edge, tuple[Ring, ...]]:
        """Return the edge incidence map in the object form used by quasi search."""
        return {
            edge: tuple(self.ring_by_id[ring_id] for ring_id in ring_ids)
            for edge, ring_ids in self.edge_to_ring_ids.items()
        }

    def face_geometries(self) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Return cached ring centers and least-squares normals when available."""
        return {
            ring_id: (self.ring_centers[ring_id], normal)
            for ring_id, normal in self.ring_normals.items()
        }


def build_ring_topology_index(
    frame: Frame,
    rings: dict[int, list[Ring]] | list[Ring],
    *,
    compute_face_quality: bool = False,
    compute_face_normals: bool = False,
) -> RingTopologyIndex:
    """Build one deterministic ring index for a frame."""
    all_rings = (
        [ring for group in rings.values() for ring in group]
        if isinstance(rings, dict)
        else list(rings)
    )
    all_rings.sort(key=lambda ring: ring.object_id)
    ring_by_id = {ring.object_id: ring for ring in all_rings}
    ring_centers: dict[str, np.ndarray] = {}
    ring_normals: dict[str, np.ndarray] = {}
    quality: dict[str, RingFaceQuality] = {}
    raw_edge_to_ids: dict[Edge, list[str]] = defaultdict(list)

    for ring in all_rings:
        coords = ring_unwrapped_coordinates(frame, ring)
        center = np.mean(coords, axis=0)
        ring_centers[ring.object_id] = center
        normal: np.ndarray | None = None
        if (compute_face_normals or compute_face_quality) and len(coords) >= 3:
            try:
                _, _, axes = np.linalg.svd(coords - center, full_matrices=False)
                candidate = np.asarray(axes[-1], dtype=float)
                norm = float(np.linalg.norm(candidate))
                if norm > 1.0e-12:
                    normal = candidate / norm
            except np.linalg.LinAlgError:
                if compute_face_quality:
                    raise
        if compute_face_normals and normal is not None:
            ring_normals[ring.object_id] = normal
        for edge in ring.edges:
            raw_edge_to_ids[edge].append(ring.object_id)
        if compute_face_quality:
            quality[ring.object_id] = measure_ring_face_quality(
                coords,
                center=center,
                normal=normal,
            )

    edge_to_ring_ids = {
        edge: tuple(sorted(ring_ids))
        for edge, ring_ids in raw_edge_to_ids.items()
    }
    adjacency: dict[str, set[str]] = {ring_id: set() for ring_id in ring_by_id}
    for ring_ids in edge_to_ring_ids.values():
        for left_pos, left in enumerate(ring_ids):
            for right in ring_ids[left_pos + 1 :]:
                adjacency[left].add(right)
                adjacency[right].add(left)

    return RingTopologyIndex(
        ring_by_id=ring_by_id,
        ring_centers=ring_centers,
        edge_to_ring_ids=edge_to_ring_ids,
        ring_adjacency={key: frozenset(sorted(value)) for key, value in adjacency.items()},
        ring_normals=ring_normals,
        face_quality=quality,
    )


def ring_unwrapped_coordinates(frame: Frame, ring: Ring) -> np.ndarray:
    """Return ordered ring coordinates after local PBC unwrapping."""
    unwrapped = unwrap_connected_nodes(frame, list(ring.nodes), list(ring.edges))
    return np.asarray([unwrapped[node] for node in ring.nodes], dtype=float)


def measure_ring_face_quality(
    coords: np.ndarray,
    *,
    center: np.ndarray | None = None,
    normal: np.ndarray | None = None,
) -> RingFaceQuality:
    """Measure best-plane deviation, edge variation, and projected area."""
    resolved_center = np.mean(coords, axis=0) if center is None else center
    centered = coords - resolved_center
    if normal is None:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        normal = vh[-1]
    deviations = centered @ normal
    planarity_rms = float(np.sqrt(np.mean(deviations * deviations)))

    rolled = np.roll(coords, -1, axis=0)
    edge_lengths = np.linalg.norm(rolled - coords, axis=1)
    edge_mean = float(np.mean(edge_lengths))
    edge_cv = float(np.std(edge_lengths) / edge_mean) if edge_mean > 1e-12 else float("inf")

    area_vector = np.zeros(3, dtype=float)
    for index in range(len(centered)):
        area_vector += np.cross(centered[index], centered[(index + 1) % len(centered)])
    projected_area = 0.5 * abs(float(np.dot(area_vector, normal)))
    return RingFaceQuality(planarity_rms, edge_cv, projected_area)
