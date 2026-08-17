"""Closed cage search by ring-face growth and polyhedron validation."""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from dataclasses import replace
from itertools import combinations

import numpy as np

from ..models import Cage, CagePatch, Frame, Guest, Ring, guest_id
from .geometry import pbc_aware_centroid, unwrap_connected_nodes
from .pbc import minimum_image
from .spatial import PointSpatialIndex
from .ring_topology import (
    RingFaceQuality,
    RingTopologyIndex,
    build_ring_topology_index,
    measure_ring_face_quality,
    ring_unwrapped_coordinates,
)


KNOWN_CAGE_TYPES = ["512", "51262", "51263", "51264", "51268", "435663"]
FACE_COUNT_SIZES = (4, 5, 6)
FACE_COUNT_INDEX = {size: index for index, size in enumerate(FACE_COUNT_SIZES)}

TARGET_FACE_COUNTS = {
    "512": {4: 0, 5: 12, 6: 0},
    "51262": {4: 0, 5: 12, 6: 2},
    "51263": {4: 0, 5: 12, 6: 3},
    "51264": {4: 0, 5: 12, 6: 4},
    "51268": {4: 0, 5: 12, 6: 8},
    "435663": {4: 3, 5: 6, 6: 3},
}

CAGE_REPORT_GROUPS = {
    "I": ("512", "51262"),
    "II": ("512", "51264"),
    "H": ("512", "51268", "435663"),
    "HS-I": ("512", "51262", "51263"),
    "TS-I": ("512", "51262", "51263"),
    "I2II": ("51263",),
}


class CageSearchLimitError(RuntimeError):
    """Exact cage search stopped by a diagnostic state guard."""


def find_cages(
    frame: Frame,
    rings: dict[int, list[Ring]],
    patches: list[CagePatch],
    guests: list[Guest],
    enabled: bool = False,
    ring_sizes: list[int] | None = None,
    max_faces: int = 20,
    search_mode: str = "grow",
    seed_mode: str = "patch",
    max_states_per_seed: int = 0,
    max_total_states: int = 0,
    max_boundary_candidates: int = 8,
    occupancy_radius_nm: float = 0.5,
    occupancy_mode: str = "polyhedron",
    scientific_validation: bool = False,
    max_face_planarity_rms_nm: float = 0.06,
    max_face_edge_cv: float = 0.35,
    min_cage_volume_nm3: float = 1.0e-6,
    topology_index: RingTopologyIndex | None = None,
) -> list[Cage]:
    """Find every valid closed cage in the configured face-size scope."""
    if not enabled:
        return []

    # Cage closure uses the established 4/5/6 face model.
    allowed_sizes = set(ring_sizes or [5, 6]) & {4, 5, 6}
    if not allowed_sizes:
        return []
    targets = build_cage_targets(allowed_sizes, max_faces)
    active_sizes = {size for counts in targets.values() for size, count in counts.items() if count > 0}
    all_rings = [ring for group in rings.values() for ring in group if ring.size in allowed_sizes and ring.size in active_sizes]
    topology = topology_index or build_ring_topology_index(
        frame,
        all_rings,
        compute_face_quality=scientific_validation,
        compute_ring_centers=False,
        compute_adjacency=False,
    )
    active_ring_ids = {ring.object_id for ring in all_rings}
    ring_by_id = {
        ring_id: ring
        for ring_id, ring in topology.ring_by_id.items()
        if ring_id in active_ring_ids
    }
    # Grow cages from shared boundary edges.
    edge_to_ring_ids = {
        edge: {ring_id for ring_id in ring_ids if ring_id in active_ring_ids}
        for edge, ring_ids in topology.edge_to_ring_ids.items()
        if any(ring_id in active_ring_ids for ring_id in ring_ids)
    }
    found: dict[frozenset[str], Cage] = {}
    seen_water_keys: set[tuple[str, tuple[int, ...]]] = set()
    type_counts: dict[str, int] = defaultdict(int)
    guest_centers = tuple(guest_center(frame, guest) for guest in guests)
    guest_center_index = (
        PointSpatialIndex(
            np.asarray(guest_centers, dtype=float),
            frame.box,
            cell_size=max(float(occupancy_radius_nm) + 0.25, 0.25),
        )
        if guest_centers
        else None
    )

    def add_candidate(face_ids: frozenset[str], cage_type: str) -> None:
        if face_ids in found:
            return
        face_rings = [ring_by_id[ring_id] for ring_id in sorted(face_ids)]
        # A matching face count still needs a closed shell.
        if not is_closed_polyhedron(face_rings):
            return
        unwrapped = cage_unwrapped_nodes(frame, face_rings)
        if scientific_validation:
            qualities = {}
            for ring in face_rings:
                quality = topology.face_quality.get(ring.object_id)
                if quality is None:
                    quality = measure_ring_face_quality(ring_unwrapped_coordinates(frame, ring))
                    topology.face_quality[ring.object_id] = quality
                qualities[ring.object_id] = quality
            geometry = scientific_polyhedron_geometry(
                face_rings,
                unwrapped,
                qualities,
                max_face_planarity_rms_nm=max_face_planarity_rms_nm,
                max_face_edge_cv=max_face_edge_cv,
                min_cage_volume_nm3=min_cage_volume_nm3,
            )
            if geometry is None:
                return
            center, _ = geometry
        else:
            center = np.mean([unwrapped[node] for node in sorted(unwrapped)], axis=0)
        waters = tuple(sorted({node for ring in face_rings for node in ring.nodes}))
        water_key = (cage_type, waters)
        if water_key in seen_water_keys:
            return
        seen_water_keys.add(water_key)
        type_counts[cage_type] += 1
        guest_ids = assigned_guests(
            frame,
            guests,
            center,
            face_rings,
            unwrapped,
            occupancy_radius_nm,
            occupancy_mode,
            guest_centers=guest_centers,
            guest_center_index=guest_center_index,
        )
        found[face_ids] = Cage(
            object_id=f"{cage_type}_{type_counts[cage_type]:05d}",
            cage_type=cage_type,
            rings=tuple(sorted(face_ids)),
            waters=waters,
            center=center,
            guest_ids=guest_ids,
            isomer=cage_isomer_label(face_rings),
        )

    mode = search_mode.lower().strip()
    if mode in {"expand", "hybrid"}:
        # Legacy names map to grow.
        mode = "grow"

    if mode == "grow":
        seed_face_sets = grow_seed_face_sets(patches, ring_by_id, seed_mode)
        for face_ids, cage_type in grow_cage_candidates(
            ring_by_id,
            edge_to_ring_ids,
            targets,
            seed_face_sets,
            max_states_per_seed=max_states_per_seed,
            max_total_states=max_total_states,
            max_boundary_candidates=max_boundary_candidates,
            status=None,
        ):
            add_candidate(face_ids, cage_type)
    elif mode in {"pair", "patch_pair"}:
        # Compatibility/debug closure of two open patches.
        patch_face_sets = [frozenset(patch.rings) for patch in patches]
        for face_ids, cage_type in fast_pair_cage_candidates(patch_face_sets, ring_by_id, targets):
            add_candidate(face_ids, cage_type)
    else:
        raise ValueError(f"Unsupported cage.search_mode: {search_mode}")

    ordered = sorted(
        found.values(),
        key=lambda cage: (cage.cage_type, cage.waters, cage.rings),
    )
    canonical_counts: dict[str, int] = defaultdict(int)
    result: list[Cage] = []
    for cage in ordered:
        canonical_counts[cage.cage_type] += 1
        result.append(
            replace(
                cage,
                object_id=f"{cage.cage_type}_{canonical_counts[cage.cage_type]:05d}",
            )
        )
    return result


def build_cage_targets(
    allowed_sizes: set[int],
    max_faces: int,
) -> dict[str, dict[int, int]]:
    """Build all Euler-compatible target compositions in the search scope."""
    targets: dict[str, dict[int, int]] = {}
    for counts in generated_other_face_counts(allowed_sizes, max_faces=max_faces):
        targets[cage_type_for_counts(counts)] = counts
    return targets


def generated_other_face_counts(allowed_sizes: set[int], max_faces: int) -> list[dict[int, int]]:
    """Generate Euler-compatible 4/5/6 cage compositions."""
    sizes = allowed_sizes & {4, 5, 6}
    if not sizes:
        return []
    targets: list[dict[int, int]] = []
    max_faces = max(1, int(max_faces))
    max_four = 6 if 4 in sizes else 0
    for n4 in range(max_four + 1):
        n5 = 12 - 2 * n4
        if n5 < 0:
            continue
        if n5 and 5 not in sizes:
            continue
        base_faces = n4 + n5
        if base_faces > max_faces:
            continue
        max_six = max_faces - base_faces if 6 in sizes else 0
        for n6 in range(max_six + 1):
            counts = {4: n4, 5: n5, 6: n6}
            if sum(counts.values()) == 0:
                continue
            targets.append(counts)
    return targets


def parse_cage_face_label(label: str) -> dict[int, int] | None:
    """Parse a named cage, numeric 1-10-2, or generic 4^1-5^10-6^2 label."""
    text = label.strip()
    if not text:
        return None
    if text in TARGET_FACE_COUNTS:
        return dict(TARGET_FACE_COUNTS[text])
    if text.count("-") == 2 and all(part.strip().isdigit() for part in text.split("-")):
        n4, n5, n6 = (int(part) for part in text.split("-"))
        return {4: n4, 5: n5, 6: n6}
    counts: dict[int, int] = {}
    for token in text.replace("_", "-").split("-"):
        token = token.strip()
        if not token:
            continue
        if "^" not in token:
            return None
        size_text, count_text = token.split("^", 1)
        if not size_text.isdigit() or not count_text.isdigit():
            return None
        counts[int(size_text)] = int(count_text)
    return counts or None


def cage_type_for_counts(counts: dict[int, int]) -> str:
    """Use a compact named label when available, otherwise a generic label."""
    for name in KNOWN_CAGE_TYPES:
        if counts_match(counts, TARGET_FACE_COUNTS[name]):
            return name
    return canonical_cage_face_label(counts)


def canonical_cage_type(label: str) -> str:
    """Normalize one report label to the cage type emitted by the search."""
    text = str(label).strip()
    counts = parse_cage_face_label(text)
    if counts is None:
        raise ValueError(
            f"Unsupported cage type '{label}'. Use a named type such as 51268 "
            "or a face-count label such as 4^1-5^10-6^2."
        )
    return cage_type_for_counts(counts)


def canonical_cage_face_label(counts: dict[int, int]) -> str:
    """Return an ASCII cage type label safe for summary columns and filenames."""
    return "-".join(f"{size}^{counts.get(size, 0)}" for size in sorted(counts) if counts.get(size, 0) > 0)


def grow_cage_candidates(
    ring_by_id: dict[str, Ring],
    edge_to_ring_ids: dict[tuple[int, int], set[str]],
    targets: dict[str, dict[int, int]],
    seed_face_sets: list[frozenset[str]],
    max_states_per_seed: int,
    max_total_states: int,
    max_boundary_candidates: int,
    status: dict[str, bool] | None = None,
):
    """Yield exact closed shells using sparse local edge incidences."""
    ring_ids = sorted(ring_by_id)
    rank = {ring_id: index for index, ring_id in enumerate(ring_ids)}
    ring_sizes = {ring_id: ring.size for ring_id, ring in ring_by_id.items()}
    ring_edges = {
        ring_id: tuple(sorted(ring.edges))
        for ring_id, ring in ring_by_id.items()
    }
    sparse_edge_to_rings = {
        edge: tuple(sorted((ring_id for ring_id in ring_ids if ring_id in rank), key=rank.__getitem__))
        for edge, ring_ids in edge_to_ring_ids.items()
    }
    target_names = tuple(targets)
    target_face_counts = tuple(
        tuple(targets[name].get(size, 0) for size in FACE_COUNT_SIZES)
        for name in target_names
    )
    target_total_incidences = tuple(
        sum(size * count for size, count in zip(FACE_COUNT_SIZES, counts, strict=True))
        for counts in target_face_counts
    )
    all_target_mask = (1 << len(target_names)) - 1
    target_count_masks = build_target_count_masks(target_names, targets)
    seed_index_by_anchor = build_seed_index_by_anchor(seed_face_sets)
    total_states = 0

    for seed_index, seed_face_ids in enumerate(seed_face_sets):
        if not seed_face_ids or any(ring_id not in rank for ring_id in seed_face_ids):
            continue
        seed_rank = rank[next(iter(seed_face_ids))] if len(seed_face_ids) == 1 else -1
        edge_counts: dict[tuple[int, int], int] = {}
        open_edges: set[tuple[int, int]] = set()
        valid_seed = True
        for ring_id in sorted(seed_face_ids, key=rank.__getitem__):
            for edge in ring_edges[ring_id]:
                previous = edge_counts.get(edge, 0)
                if previous >= 2:
                    valid_seed = False
                    break
                if previous == 0:
                    edge_counts[edge] = 1
                    open_edges.add(edge)
                else:
                    edge_counts[edge] = 2
                    open_edges.remove(edge)
            if not valid_seed:
                break
        if not valid_seed:
            continue

        seed_counts = face_count_tuple(seed_face_ids, ring_sizes)
        seed_incidence = sum(
            size * count for size, count in zip(FACE_COUNT_SIZES, seed_counts, strict=True)
        )
        compatible_targets = compatible_target_mask_tuple(
            seed_counts,
            all_target_mask,
            target_count_masks,
        )
        compatible_targets = prune_target_mask_by_edge_budget(
            compatible_targets,
            seed_incidence,
            len(open_edges),
            target_total_incidences,
        )
        if not compatible_targets:
            continue

        faces = sorted(seed_face_ids, key=rank.__getitem__)
        face_set = set(faces)
        seen_states: set[tuple[str, ...]] = set()
        local_states = 0

        def visit(
            counts: tuple[int, int, int],
            used_incidence: int,
            compatible: int,
        ):
            nonlocal local_states, total_states
            state_key = tuple(faces)
            if state_key in seen_states:
                return
            if len(seed_face_ids) > 1 and contains_earlier_seed(
                frozenset(face_set), seed_index, seed_index_by_anchor
            ):
                return
            if max_states_per_seed > 0 and local_states >= max_states_per_seed:
                if status is not None:
                    status["hit_seed_limit"] = True
                raise CageSearchLimitError(
                    "Exact cage search reached "
                    f"cage.max_state_per_seed={max_states_per_seed} at seed "
                    f"{min(seed_face_ids)}; no partial frame result was accepted."
                )
            if max_total_states > 0 and total_states >= max_total_states:
                if status is not None:
                    status["hit_total_limit"] = True
                raise CageSearchLimitError(
                    "Exact cage search reached "
                    f"cage.max_total_state={max_total_states}; "
                    "no partial frame result was accepted."
                )
            seen_states.add(state_key)
            local_states += 1
            total_states += 1

            if not open_edges:
                cage_type = target_type_for_count_tuple(
                    counts, target_names, target_face_counts
                )
                if cage_type is not None:
                    yield frozenset(faces), cage_type
                return

            next_ids = ordered_boundary_candidates_sparse(
                faces,
                face_set,
                open_edges,
                edge_counts,
                counts,
                compatible,
                target_count_masks,
                sparse_edge_to_rings,
                ring_sizes,
                ring_edges,
                rank,
                seed_rank,
                max_boundary_candidates,
            )
            for next_id in next_ids:
                ring_size = ring_sizes[next_id]
                count_index = FACE_COUNT_INDEX[ring_size]
                next_count = counts[count_index] + 1
                next_counts = counts[:count_index] + (next_count,) + counts[count_index + 1 :]
                next_compatible = compatible & target_count_masks.get(
                    (ring_size, next_count), 0
                )
                if not next_compatible:
                    continue

                if any(edge_counts.get(edge, 0) >= 2 for edge in ring_edges[next_id]):
                    continue
                changes: list[tuple[tuple[int, int], int]] = []
                for edge in ring_edges[next_id]:
                    previous = edge_counts.get(edge, 0)
                    changes.append((edge, previous))
                    if previous == 0:
                        edge_counts[edge] = 1
                        open_edges.add(edge)
                    else:
                        edge_counts[edge] = 2
                        open_edges.remove(edge)
                next_compatible = prune_target_mask_by_edge_budget(
                    next_compatible,
                    used_incidence + ring_size,
                    len(open_edges),
                    target_total_incidences,
                )
                if next_compatible:
                    position = bisect_left(faces, next_id)
                    faces.insert(position, next_id)
                    face_set.add(next_id)
                    try:
                        yield from visit(
                            next_counts,
                            used_incidence + ring_size,
                            next_compatible,
                        )
                    finally:
                        face_set.remove(next_id)
                        faces.pop(position)

                for edge, previous in reversed(changes):
                    if previous == 0:
                        del edge_counts[edge]
                        open_edges.discard(edge)
                    else:
                        edge_counts[edge] = 1
                        open_edges.add(edge)

        yield from visit(seed_counts, seed_incidence, compatible_targets)


def face_count_tuple(face_ids: frozenset[str], ring_sizes: dict[str, int]) -> tuple[int, int, int]:
    """Count 4/5/6 faces without allocating a per-state dictionary."""
    counts = [0, 0, 0]
    for ring_id in face_ids:
        counts[FACE_COUNT_INDEX[ring_sizes[ring_id]]] += 1
    return tuple(counts)


def compatible_target_mask_tuple(
    counts: tuple[int, int, int],
    all_target_mask: int,
    target_count_masks: dict[tuple[int, int], int],
) -> int:
    """Return target bits compatible with fixed 4/5/6 face counts."""
    mask = all_target_mask
    for size, count in zip(FACE_COUNT_SIZES, counts, strict=True):
        mask &= target_count_masks.get((size, count), 0)
        if not mask:
            break
    return mask


def target_type_for_count_tuple(
    counts: tuple[int, int, int],
    target_names: tuple[str, ...],
    target_face_counts: tuple[tuple[int, int, int], ...],
) -> str | None:
    """Resolve a closed fixed-size face-count state to a configured target."""
    for name, target_counts in zip(target_names, target_face_counts, strict=True):
        if counts == target_counts:
            return name
    return None


def build_target_count_masks(
    target_names: tuple[str, ...],
    targets: dict[str, dict[int, int]],
) -> dict[tuple[int, int], int]:
    """Precompute target bits that permit each face-size count."""
    sizes = sorted({size for target in targets.values() for size in target})
    masks: dict[tuple[int, int], int] = {}
    for size in sizes:
        maximum = max((target.get(size, 0) for target in targets.values()), default=0)
        for count in range(maximum + 1):
            mask = 0
            for index, name in enumerate(target_names):
                if count <= targets[name].get(size, 0):
                    mask |= 1 << index
            masks[(size, count)] = mask
    return masks


def prune_target_mask_by_edge_budget(
    target_mask: int,
    used_incidence: int,
    open_edge_count: int,
    target_total_incidences: tuple[int, ...],
) -> int:
    """Keep targets whose remaining incidences can close every open edge."""
    kept = 0
    remaining = target_mask
    while remaining:
        bit = remaining & -remaining
        remaining ^= bit
        remaining_incidence = target_total_incidences[bit.bit_length() - 1] - used_incidence
        if remaining_incidence < open_edge_count:
            continue
        if (remaining_incidence - open_edge_count) % 2:
            continue
        kept |= bit
    return kept


def ordered_boundary_candidates_sparse(
    _face_ids: list[str],
    face_set: set[str],
    open_edges: set[tuple[int, int]],
    edge_counts: dict[tuple[int, int], int],
    counts: tuple[int, int, int],
    compatible_targets: int,
    target_count_masks: dict[tuple[int, int], int],
    edge_to_ring_ids: dict[tuple[int, int], tuple[str, ...]],
    ring_sizes: dict[str, int],
    ring_edges: dict[str, tuple[tuple[int, int], ...]],
    rank: dict[str, int],
    seed_rank: int,
    _candidate_batch_width: int,
) -> list[str]:
    """Choose the minimum-candidate open edge without dropping branches."""
    best: list[str] | None = None
    for edge in sorted(open_edges):
        candidates: list[str] = []
        for ring_id in edge_to_ring_ids.get(edge, ()):
            if ring_id in face_set or rank[ring_id] < seed_rank:
                continue
            if any(edge_counts.get(item, 0) >= 2 for item in ring_edges[ring_id]):
                continue
            ring_size = ring_sizes[ring_id]
            next_count = counts[FACE_COUNT_INDEX[ring_size]] + 1
            if not compatible_targets & target_count_masks.get((ring_size, next_count), 0):
                continue
            candidates.append(ring_id)
        if not candidates:
            return []
        if best is None or len(candidates) < len(best):
            best = candidates
    if not best:
        return []
    return sorted(best, key=rank.__getitem__)


def grow_seed_face_sets(patches: list[CagePatch], ring_by_id: dict[str, Ring], seed_mode: str) -> list[frozenset[str]]:
    """Build cage-grow seeds from open cage patches, or rings for comparison."""
    mode = seed_mode.lower().strip()
    if mode not in {"patch", "ring", "auto"}:
        raise ValueError(f"Unsupported cage.seed_mode: {seed_mode}")

    patch_seeds: list[frozenset[str]] = []
    seen: set[frozenset[str]] = set()
    for patch in patches:
        face_ids = frozenset(patch.rings)
        if face_ids in seen:
            continue
        if face_ids and all(ring_id in ring_by_id for ring_id in face_ids):
            patch_seeds.append(face_ids)
            seen.add(face_ids)

    if mode == "auto":
        mode = "patch" if patch_seeds else "ring"
    if mode == "patch":
        return sorted(patch_seeds, key=lambda item: (min(item), len(item), sorted(item)))
    return [frozenset([ring_id]) for ring_id in sorted(ring_by_id)]


def build_seed_index_by_anchor(seed_face_sets: list[frozenset[str]]) -> dict[str, list[tuple[int, frozenset[str]]]]:
    """Index patch seeds by their smallest ring id for fast subset checks."""
    by_anchor: dict[str, list[tuple[int, frozenset[str]]]] = defaultdict(list)
    for index, seed in enumerate(seed_face_sets):
        if len(seed) <= 1:
            continue
        by_anchor[min(seed)].append((index, seed))
    return dict(by_anchor)


def contains_earlier_seed(
    face_ids: frozenset[str],
    seed_index: int,
    seed_index_by_anchor: dict[str, list[tuple[int, frozenset[str]]]],
) -> bool:
    """Prune duplicate patch-seed growth once an earlier seed is contained."""
    for ring_id in face_ids:
        for earlier_index, earlier in seed_index_by_anchor.get(ring_id, []):
            if earlier_index >= seed_index:
                continue
            if earlier <= face_ids:
                return True
    return False


def fast_pair_cage_candidates(
    patch_face_sets: list[frozenset[str]],
    ring_by_id: dict[str, Ring],
    targets: dict[str, dict[int, int]],
):
    """Yield cage candidates from patch pairs with matching open boundaries."""
    buckets: dict[frozenset[tuple[int, int]], list[tuple[int, frozenset[str], dict[int, int]]]] = defaultdict(list)
    for idx, face_ids in enumerate(patch_face_sets):
        edge_counts = candidate_edge_counts(face_ids, ring_by_id)
        if any(count > 2 for count in edge_counts.values()):
            continue
        boundary = frozenset(edge for edge, count in edge_counts.items() if count == 1)
        if boundary:
            buckets[boundary].append((idx, face_ids, face_counts(face_ids, ring_by_id)))

    emitted: set[frozenset[str]] = set()
    for infos in buckets.values():
        if len(infos) < 2:
            continue
        for left_pos, (idx_a, faces_a, counts_a) in enumerate(infos):
            for idx_b, faces_b, counts_b in infos[left_pos + 1 :]:
                if idx_b <= idx_a or faces_a & faces_b:
                    continue
                merged_counts = merge_face_counts(counts_a, counts_b)
                cage_type = next((name for name, target in targets.items() if counts_match(merged_counts, target)), None)
                if cage_type is None:
                    continue
                candidate = faces_a | faces_b
                if candidate in emitted:
                    continue
                emitted.add(candidate)
                yield candidate, cage_type


def fast_patch_closure_candidates(
    patch_face_sets: list[frozenset[str]],
    ring_by_id: dict[str, Ring],
    targets: dict[str, dict[int, int]],
    *,
    max_patches: int = 4,
    max_states: int = 20000,
    status: dict[str, bool] | None = None,
):
    """Yield closed cages assembled from two to four indexed half-cage patches."""
    unique = sorted(
        {
            face_ids
            for face_ids in patch_face_sets
            if face_ids and all(ring_id in ring_by_id for ring_id in face_ids)
        },
        key=lambda item: (len(item), tuple(sorted(item))),
    )
    if len(unique) < 2 or max_patches < 2 or max_states <= 0:
        return

    boundaries = []
    for face_ids in unique:
        edge_counts = candidate_edge_counts(face_ids, ring_by_id)
        boundaries.append({edge for edge, count in edge_counts.items() if count == 1})

    neighbors: dict[int, set[int]] = {index: set() for index in range(len(unique))}
    for left, right in combinations(range(len(unique)), 2):
        if unique[left] & unique[right] or boundaries[left] & boundaries[right]:
            neighbors[left].add(right)
            neighbors[right].add(left)

    max_counts = max_target_face_counts(targets)
    emitted: set[frozenset[str]] = set()
    seen_combinations: set[tuple[int, ...]] = set()
    state_count = 0

    for start in range(len(unique)):
        stack = [((start,), unique[start], {item for item in neighbors[start] if item > start})]
        while stack:
            combo, face_ids, frontier = stack.pop()
            if combo in seen_combinations:
                continue
            seen_combinations.add(combo)
            state_count += 1
            if state_count > max_states:
                if status is not None:
                    status["hit_state_limit"] = True
                return

            if len(combo) >= 2:
                cage_type = target_type_for_faces(face_ids, ring_by_id, targets)
                if cage_type is not None and face_ids not in emitted:
                    face_rings = [ring_by_id[ring_id] for ring_id in sorted(face_ids)]
                    if is_closed_polyhedron(face_rings):
                        emitted.add(face_ids)
                        yield face_ids, cage_type
            if len(combo) >= max_patches:
                continue

            for next_index in sorted(frontier, reverse=True):
                next_combo = tuple(sorted((*combo, next_index)))
                next_faces = face_ids | unique[next_index]
                if next_combo in seen_combinations:
                    continue
                if not counts_fit(face_counts(next_faces, ring_by_id), max_counts):
                    continue
                if not no_edge_overuse(next_faces, ring_by_id):
                    continue
                next_frontier = set(frontier)
                next_frontier.update(neighbors[next_index])
                next_frontier.difference_update(next_combo)
                next_frontier = {item for item in next_frontier if item > start}
                stack.append((next_combo, next_faces, next_frontier))


def build_edge_to_ring_ids(rings: list[Ring]) -> dict[tuple[int, int], set[str]]:
    """Map each oxygen-network edge to all ring faces that use it."""
    edge_to_ring_ids: dict[tuple[int, int], set[str]] = defaultdict(set)
    for ring in rings:
        for edge in ring.edges:
            edge_to_ring_ids[edge].add(ring.object_id)
    return dict(edge_to_ring_ids)


def target_type_for_faces(face_ids: frozenset[str], ring_by_id: dict[str, Ring], targets: dict[str, dict[int, int]]) -> str | None:
    """Return the target cage type matching a face set, if any."""
    counts = face_counts(face_ids, ring_by_id)
    return target_type_for_counts(counts, targets)


def target_type_for_counts(counts: dict[int, int], targets: dict[str, dict[int, int]]) -> str | None:
    """Return the target cage type matching a face-count map, if any."""
    for cage_type, target_counts in targets.items():
        if counts_match(counts, target_counts):
            return cage_type
    return None


def max_target_face_counts(targets: dict[str, dict[int, int]]) -> dict[int, int]:
    """Merge target cage limits so multiple cage types can share one grow pass."""
    merged: dict[int, int] = {}
    for target_counts in targets.values():
        for size, count in target_counts.items():
            merged[size] = max(merged.get(size, 0), count)
    return merged


def face_counts(face_ids: frozenset[str], ring_by_id: dict[str, Ring]) -> dict[int, int]:
    """Count ring-face sizes in a candidate shell."""
    counts: dict[int, int] = {}
    for ring_id in face_ids:
        size = ring_by_id[ring_id].size
        counts[size] = counts.get(size, 0) + 1
    return counts


def counts_fit(counts: dict[int, int], target_counts: dict[int, int]) -> bool:
    """Check whether a partial shell can still grow into the target type."""
    sizes = set(counts) | set(target_counts)
    return all(counts.get(size, 0) <= target_counts.get(size, 0) for size in sizes)


def counts_match(counts: dict[int, int], target_counts: dict[int, int]) -> bool:
    """Compare face counts while treating missing sizes as zero."""
    sizes = set(counts) | set(target_counts)
    return all(counts.get(size, 0) == target_counts.get(size, 0) for size in sizes)


def merge_face_counts(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    """Merge two face-count maps."""
    sizes = set(left) | set(right)
    return {size: left.get(size, 0) + right.get(size, 0) for size in sizes}


def boundary_candidate_ring_ids(
    face_ids: frozenset[str],
    ring_by_id: dict[str, Ring],
    edge_to_ring_ids: dict[tuple[int, int], set[str]],
) -> set[str]:
    """Find rings touching current boundary edges."""
    edge_counts = candidate_edge_counts(face_ids, ring_by_id)
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    candidates: set[str] = set()
    for edge in boundary_edges:
        candidates.update(edge_to_ring_ids.get(edge, set()))
    return candidates - set(face_ids)


def candidate_edge_counts(face_ids: frozenset[str], ring_by_id: dict[str, Ring]) -> dict[tuple[int, int], int]:
    """Count how many candidate faces use each edge."""
    edge_counts: dict[tuple[int, int], int] = defaultdict(int)
    for ring_id in face_ids:
        for edge in ring_by_id[ring_id].edges:
            edge_counts[edge] += 1
    return dict(edge_counts)


def no_edge_overuse(face_ids: frozenset[str], ring_by_id: dict[str, Ring]) -> bool:
    """Reject partial shells where any edge is already used by more than two faces."""
    return all(count <= 2 for count in candidate_edge_counts(face_ids, ring_by_id).values())


def is_closed_polyhedron(rings: list[Ring]) -> bool:
    """Require a connected, trivalent, manifold spherical shell."""
    edge_counts: dict[tuple[int, int], int] = defaultdict(int)
    nodes: set[int] = set()
    neighbors: dict[int, set[int]] = defaultdict(set)
    for ring in rings:
        nodes.update(ring.nodes)
        for edge in ring.edges:
            edge_counts[edge] += 1
            left, right = edge
            neighbors[left].add(right)
            neighbors[right].add(left)
    if not edge_counts or any(count != 2 for count in edge_counts.values()):
        return False
    if len(nodes) - len(edge_counts) + len(rings) != 2:
        return False
    if any(len(neighbors[node]) != 3 for node in nodes):
        return False
    return face_adjacency_connected(rings) and vertex_links_are_manifold(rings)



def scientific_polyhedron_geometry(
    rings: list[Ring],
    unwrapped_nodes: dict[int, np.ndarray],
    face_quality: dict[str, RingFaceQuality],
    *,
    max_face_planarity_rms_nm: float,
    max_face_edge_cv: float,
    min_cage_volume_nm3: float,
) -> tuple[np.ndarray, float] | None:
    """Apply opt-in face-quality and volume checks."""
    for ring in rings:
        quality = face_quality.get(ring.object_id)
        if quality is None:
            return None
        if quality.projected_area_nm2 <= 1.0e-12:
            return None
        if quality.planarity_rms_nm > max_face_planarity_rms_nm:
            return None
        if quality.edge_length_cv > max_face_edge_cv:
            return None

    reference = np.mean([unwrapped_nodes[node] for node in sorted(unwrapped_nodes)], axis=0)
    triangles = triangulate_cage_faces(rings, unwrapped_nodes, reference)
    geometry = polyhedron_volume_centroid(triangles, reference)
    if geometry is None:
        return None
    center, volume = geometry
    if volume < min_cage_volume_nm3:
        return None
    return center, volume


def face_adjacency_connected(rings: list[Ring]) -> bool:
    """Require all shell faces to belong to one edge-connected component."""
    if not rings:
        return False
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_index, ring in enumerate(rings):
        for edge in ring.edges:
            edge_faces[edge].append(face_index)
    adjacency = {index: set() for index in range(len(rings))}
    for face_indexes in edge_faces.values():
        if len(face_indexes) != 2:
            return False
        left, right = face_indexes
        adjacency[left].add(right)
        adjacency[right].add(left)

    visited = set()
    stack = [0]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        stack.extend(adjacency[current] - visited)
    return len(visited) == len(rings)


def vertex_links_are_manifold(rings: list[Ring]) -> bool:
    """Require the incident faces around every shell vertex to form one cycle."""
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    vertex_faces: dict[int, set[int]] = defaultdict(set)
    for face_index, ring in enumerate(rings):
        for node in ring.nodes:
            vertex_faces[node].add(face_index)
        for edge in ring.edges:
            edge_faces[edge].append(face_index)

    for node, incident in vertex_faces.items():
        link = {face_index: set() for face_index in incident}
        for edge, face_indexes in edge_faces.items():
            if node not in edge:
                continue
            if len(face_indexes) != 2:
                return False
            left, right = face_indexes
            if left in link and right in link:
                link[left].add(right)
                link[right].add(left)
        if any(len(neighbors) != 2 for neighbors in link.values()):
            return False
        visited = set()
        stack = [next(iter(link))]
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            stack.extend(link[current] - visited)
        if visited != set(link):
            return False
    return True


def polyhedron_volume_centroid(
    triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    reference: np.ndarray,
) -> tuple[np.ndarray, float] | None:
    """Return volume centroid and absolute volume from an oriented triangle shell."""
    signed_volume = 0.0
    weighted_centroid = np.zeros(3, dtype=float)
    for a, b, c in triangles:
        ar = a - reference
        br = b - reference
        cr = c - reference
        tetra_volume = float(np.dot(ar, np.cross(br, cr))) / 6.0
        signed_volume += tetra_volume
        weighted_centroid += tetra_volume * (ar + br + cr) / 4.0
    if abs(signed_volume) <= 1.0e-12:
        return None
    center = reference + weighted_centroid / signed_volume
    return center, abs(signed_volume)



def cage_unwrapped_nodes(frame: Frame, rings: list[Ring]) -> dict[int, np.ndarray]:
    """Return locally unwrapped oxygen coordinates for a cage shell."""
    nodes = sorted({node for ring in rings for node in ring.nodes})
    edges = sorted({edge for ring in rings for edge in ring.edges})
    return unwrap_connected_nodes(frame, nodes, edges)


def cage_center(frame: Frame, rings: list[Ring]) -> np.ndarray:
    """Compute a locally unwrapped oxygen centroid for a cage."""
    nodes = sorted({node for ring in rings for node in ring.nodes})
    edges = sorted({edge for ring in rings for edge in ring.edges})
    unwrapped = unwrap_connected_nodes(frame, nodes, edges)
    return np.mean([unwrapped[node] for node in nodes], axis=0)


def assigned_guests(
    frame: Frame,
    guests: list[Guest],
    center: np.ndarray,
    rings: list[Ring],
    unwrapped_nodes: dict[int, np.ndarray],
    radius_nm: float,
    occupancy_mode: str,
    *,
    guest_centers: tuple[np.ndarray, ...] | None = None,
    guest_center_index: PointSpatialIndex | None = None,
) -> tuple[str, ...]:
    """Assign guests with one spatial candidate query and batched solid angles."""
    mode = occupancy_mode.lower().strip()
    if mode not in {"polyhedron", "center", "auto"}:
        raise ValueError(f"Unsupported cage.occupancy_mode: {occupancy_mode}")
    if not guests:
        return ()

    shell_radius = max(float(np.linalg.norm(pos - center)) for pos in unwrapped_nodes.values()) if unwrapped_nodes else radius_nm
    triangles = triangulate_cage_faces(rings, unwrapped_nodes, center) if mode in {"polyhedron", "auto"} else []
    resolved_centers = guest_centers or tuple(guest_center(frame, guest) for guest in guests)
    candidate_radius = float(radius_nm) if mode == "center" else shell_radius + 0.25
    if mode == "auto":
        candidate_radius = max(candidate_radius, float(radius_nm))
    candidate_indices = (
        guest_center_index.query(center, candidate_radius)
        if guest_center_index is not None
        else tuple(range(len(guests)))
    )
    if not candidate_indices:
        return ()

    raw_centers = np.asarray([resolved_centers[index] for index in candidate_indices], dtype=float)
    local_deltas = minimum_image(raw_centers - center, frame.box)
    local_distances = np.linalg.norm(local_deltas, axis=1)
    local_centers = center + local_deltas
    center_hits = local_distances <= float(radius_nm)
    poly_hits = np.zeros(len(candidate_indices), dtype=bool)
    if triangles:
        poly_mask = local_distances <= shell_radius + 0.25
        if np.any(poly_mask):
            poly_hits[poly_mask] = points_in_polyhedron(local_centers[poly_mask], triangles)

    assigned: list[str] = []
    for index, center_hit, poly_hit in zip(candidate_indices, center_hits, poly_hits, strict=True):
        if (
            (mode == "polyhedron" and bool(poly_hit))
            or (mode == "center" and bool(center_hit))
            or (mode == "auto" and (bool(poly_hit) or bool(center_hit)))
        ):
            assigned.append(guest_id(guests[index]))
    return tuple(assigned)


def triangulate_cage_faces(
    rings: list[Ring],
    unwrapped_nodes: dict[int, np.ndarray],
    center: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Triangulate each ring face and orient triangles away from the cage center."""
    triangles = []
    for ring in rings:
        vertices = [unwrapped_nodes[node] for node in ring.nodes]
        if len(vertices) < 3:
            continue
        anchor = vertices[0]
        for idx in range(1, len(vertices) - 1):
            a = anchor
            b = vertices[idx]
            c = vertices[idx + 1]
            normal = np.cross(b - a, c - a)
            tri_center = (a + b + c) / 3.0
            if float(np.dot(normal, tri_center - center)) < 0.0:
                b, c = c, b
            triangles.append((a, b, c))
    return triangles


def point_in_polyhedron(point: np.ndarray, triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> bool:
    """Use the oriented solid-angle sum to test whether a point is inside."""
    solid_angle = 0.0
    for a, b, c in triangles:
        solid_angle += triangle_solid_angle(a - point, b - point, c - point)
    return abs(solid_angle) > 2.0 * np.pi


def points_in_polyhedron(points: np.ndarray, triangles: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> np.ndarray:
    """Vectorize solid-angle membership across candidate points.

    Triangle accumulation stays in source order. Candidates numerically close to
    the historical 2*pi decision boundary use the scalar implementation.
    """
    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    total = np.zeros(len(values), dtype=float)
    for a, b, c in triangles:
        av = a - values
        bv = b - values
        cv = c - values
        la = np.linalg.norm(av, axis=1)
        lb = np.linalg.norm(bv, axis=1)
        lc = np.linalg.norm(cv, axis=1)
        valid = np.minimum(np.minimum(la, lb), lc) > 1.0e-12
        if not np.any(valid):
            continue
        numerator = np.einsum("ij,ij->i", av, np.cross(bv, cv))
        denominator = (
            la * lb * lc
            + np.einsum("ij,ij->i", av, bv) * lc
            + np.einsum("ij,ij->i", bv, cv) * la
            + np.einsum("ij,ij->i", cv, av) * lb
        )
        angles = np.zeros(len(values), dtype=float)
        angles[valid] = 2.0 * np.arctan2(numerator[valid], denominator[valid])
        total += angles
    hits = np.abs(total) > 2.0 * np.pi
    uncertain = np.abs(np.abs(total) - 2.0 * np.pi) <= 1.0e-9
    for index in np.flatnonzero(uncertain):
        hits[index] = point_in_polyhedron(values[index], triangles)
    return hits


def triangle_solid_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Return the signed solid angle of a triangle at the origin."""
    la = float(np.linalg.norm(a))
    lb = float(np.linalg.norm(b))
    lc = float(np.linalg.norm(c))
    if min(la, lb, lc) <= 1e-12:
        return 0.0
    numerator = float(np.dot(a, np.cross(b, c)))
    denominator = la * lb * lc + float(np.dot(a, b)) * lc + float(np.dot(b, c)) * la + float(np.dot(c, a)) * lb
    return 2.0 * float(np.arctan2(numerator, denominator))


def guest_center(frame: Frame, guest: Guest) -> np.ndarray:
    """Use the configured center atom when available, otherwise geometry center."""
    if guest.center_atom is not None:
        return frame.atoms[guest.center_atom].xyz
    return pbc_aware_centroid(frame, guest.atoms)


def cage_isomer_label(rings: list[Ring]) -> str | None:
    """Describe the adjacency pattern among 6-ring cage faces."""
    hex_rings = [ring for ring in rings if ring.size == 6]
    if not hex_rings:
        return None
    adjacency: dict[int, set[int]] = {idx: set() for idx in range(len(hex_rings))}
    for left, right in combinations(range(len(hex_rings)), 2):
        if hex_rings[left].edges & hex_rings[right].edges:
            adjacency[left].add(right)
            adjacency[right].add(left)
    return hex_adjacency_label(adjacency)


def hex_adjacency_label(adjacency: dict[int, set[int]]) -> str:
    """Return a readable cage-isomer label without opaque iso numbers."""
    n_hex = len(adjacency)
    edge_count = sum(len(neighbors) for neighbors in adjacency.values()) // 2
    degrees = sorted((len(neighbors) for neighbors in adjacency.values()), reverse=True)
    if n_hex == 1:
        return "6single"
    if edge_count == 0:
        return f"{n_hex}x6sep"
    if n_hex == 2:
        return "6adj"
    if n_hex == 3:
        if edge_count == 1:
            return "6pair+single"
        if edge_count == 2:
            return "6chain3"
        if edge_count == 3:
            return "6tri3"
    if n_hex == 4:
        if edge_count == 1:
            return "6pair+2single"
        if edge_count == 2:
            return "2x6pair" if degrees == [1, 1, 1, 1] else "6chain3+single"
        if edge_count == 3:
            if degrees == [3, 1, 1, 1]:
                return "6star3"
            if degrees == [2, 2, 1, 1]:
                return "6chain4"
            if degrees == [2, 2, 2, 0]:
                return "6tri3+single"
        if edge_count == 4:
            return "6cycle4" if degrees == [2, 2, 2, 2] else "6tri3+tail"
        if edge_count == 5:
            return "6K4-e"
        if edge_count == 6:
            return "6K4"
    degree_text = "".join(str(degree) for degree in degrees)
    return f"6n{n_hex}e{edge_count}d{degree_text}"


