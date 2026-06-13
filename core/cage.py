from __future__ import annotations

"""Closed cage search by ring-face growth and polyhedron validation."""

from collections import defaultdict
from itertools import combinations

import numpy as np

from ..models import Cage, Cup, Frame, Guest, Ring
from .pbc import distance, minimum_image


KNOWN_CAGE_TYPES = ["512", "51262", "51263", "51264"]

TARGET_FACE_COUNTS = {
    "512": {4: 0, 5: 12, 6: 0},
    "51262": {4: 0, 5: 12, 6: 2},
    "51263": {4: 0, 5: 12, 6: 3},
    "51264": {4: 0, 5: 12, 6: 4},
}


def find_cages(
    frame: Frame,
    rings: dict[int, list[Ring]],
    cups: list[Cup],
    guests: list[Guest],
    enabled: bool = False,
    target_types: list[str] | None = None,
    ring_sizes: list[int] | None = None,
    output_other: bool = False,
    other_max_faces: int = 20,
    search_mode: str = "grow",
    seed_mode: str = "cup",
    max_states_per_seed: int = 2000,
    max_total_states: int = 250000,
    occupancy_radius_nm: float = 0.5,
    occupancy_mode: str = "polyhedron",
) -> list[Cage]:
    """Find closed hydrate cages with target face counts."""
    if not enabled:
        return []

    allowed_sizes = set(ring_sizes or [5, 6])
    unsupported_sizes = allowed_sizes - {4, 5, 6}
    if unsupported_sizes:
        raise ValueError(f"cage.ring_sizes supports only 4, 5, and 6; got {sorted(unsupported_sizes)}")
    targets = build_cage_targets(target_types, allowed_sizes, output_other, other_max_faces)
    active_sizes = {size for counts in targets.values() for size, count in counts.items() if count > 0}
    all_rings = [ring for group in rings.values() for ring in group if ring.size in allowed_sizes and ring.size in active_sizes]
    ring_by_id = {ring.object_id: ring for ring in all_rings}
    edge_to_ring_ids = build_edge_to_ring_ids(all_rings)
    found: dict[frozenset[str], Cage] = {}
    seen_water_keys: set[tuple[str, tuple[int, ...]]] = set()
    type_counts: dict[str, int] = defaultdict(int)

    def add_candidate(face_ids: frozenset[str], cage_type: str) -> None:
        if face_ids in found:
            return
        face_rings = [ring_by_id[ring_id] for ring_id in sorted(face_ids)]
        if not is_closed_polyhedron(face_rings):
            return
        waters = tuple(sorted({node for ring in face_rings for node in ring.nodes}))
        water_key = (cage_type, waters)
        if water_key in seen_water_keys:
            return
        seen_water_keys.add(water_key)
        type_counts[cage_type] += 1
        unwrapped = cage_unwrapped_nodes(frame, face_rings)
        center = np.mean([unwrapped[node] for node in sorted(unwrapped)], axis=0)
        guest_ids = assigned_guests(
            frame,
            guests,
            center,
            face_rings,
            unwrapped,
            occupancy_radius_nm,
            occupancy_mode,
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
        # Older drafts used these names. In SQQ they are aliases for grow.
        mode = "grow"

    if mode == "grow":
        seed_face_sets = grow_seed_face_sets(cups, ring_by_id, seed_mode)
        for face_ids, cage_type in grow_cage_candidates(
            ring_by_id,
            edge_to_ring_ids,
            targets,
            seed_face_sets,
            max_states_per_seed=max_states_per_seed,
            max_total_states=max_total_states,
        ):
            add_candidate(face_ids, cage_type)
    elif mode in {"pair", "cup_pair"}:
        # Compatibility/debug path only: two cup boundaries close into a shell.
        cup_face_sets = [frozenset(cup.rings) for cup in cups]
        for face_ids, cage_type in fast_pair_cage_candidates(cup_face_sets, ring_by_id, targets):
            add_candidate(face_ids, cage_type)
    else:
        raise ValueError(f"Unsupported cage.search_mode: {search_mode}")

    return sorted(found.values(), key=lambda cage: (cage.cage_type, cage.object_id))


def build_cage_targets(
    target_types: list[str] | None,
    allowed_sizes: set[int],
    output_other: bool,
    other_max_faces: int,
) -> dict[str, dict[int, int]]:
    """Build target face-count maps for standard and optional other cages."""
    targets: dict[str, dict[int, int]] = {}
    requested = target_types or KNOWN_CAGE_TYPES
    for raw_name in requested:
        name = str(raw_name)
        if name in TARGET_FACE_COUNTS:
            counts = TARGET_FACE_COUNTS[name]
            if face_counts_use_allowed_sizes(counts, allowed_sizes):
                targets[name] = dict(counts)
            continue
        counts = parse_cage_face_label(name)
        if counts and face_counts_use_allowed_sizes(counts, allowed_sizes):
            targets[canonical_cage_face_label(counts)] = counts

    if output_other:
        for counts in generated_other_face_counts(allowed_sizes, max_faces=other_max_faces):
            if any(counts_match(counts, existing) for existing in targets.values()):
                continue
            targets[canonical_cage_face_label(counts)] = counts
    return targets


def face_counts_use_allowed_sizes(counts: dict[int, int], allowed_sizes: set[int]) -> bool:
    """Reject target types that require ring sizes disabled for cage search."""
    return all(size in allowed_sizes or count == 0 for size, count in counts.items())


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
    """Parse HTR-style 1-10-2 or generic 4^1-5^10-6^2 labels."""
    text = label.strip()
    if not text:
        return None
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
):
    """Yield closed shells grown from cup or ring face seeds."""
    ring_ids = sorted(ring_by_id)
    rank = {ring_id: idx for idx, ring_id in enumerate(ring_ids)}
    seed_index_by_anchor = build_seed_index_by_anchor(seed_face_sets)
    seen_states: set[tuple[str, frozenset[str]]] = set()
    total_states = 0

    for seed_index, seed_face_ids in enumerate(seed_face_sets):
        # Single-ring seeds use rank pruning to avoid rediscovering a shell
        # from every face. Cup seeds already represent larger directed patches;
        # allowing lower-ranked added faces avoids missing cages whose complete
        # shell has an earlier ring outside the selected cup seed.
        seed_rank = rank[next(iter(seed_face_ids))] if len(seed_face_ids) == 1 else -1
        seed_edge_counts = candidate_edge_counts(seed_face_ids, ring_by_id)
        if any(count > 2 for count in seed_edge_counts.values()):
            continue
        seed_counts = face_counts(seed_face_ids, ring_by_id)
        for cage_type, target_counts in targets.items():
            if total_states >= max_total_states:
                return
            if not counts_fit(seed_counts, target_counts):
                continue
            stack = [(seed_face_ids, seed_edge_counts, seed_counts)]
            local_states = 0

            while stack and local_states < max_states_per_seed and total_states < max_total_states:
                face_ids, edge_counts, counts = stack.pop()
                state_key = (cage_type, face_ids)
                if state_key in seen_states:
                    continue
                if len(seed_face_ids) > 1 and contains_earlier_seed(face_ids, seed_index, seed_index_by_anchor):
                    continue
                seen_states.add(state_key)
                local_states += 1
                total_states += 1

                boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
                if not boundary_edges:
                    if counts_match(counts, target_counts):
                        yield face_ids, cage_type
                    continue
                if counts_match(counts, target_counts):
                    continue

                next_ids = ordered_boundary_candidates(
                    face_ids,
                    edge_counts,
                    counts,
                    target_counts,
                    edge_to_ring_ids,
                    ring_by_id,
                    rank,
                    seed_rank,
                )
                for next_id in reversed(next_ids):
                    ring = ring_by_id[next_id]
                    next_edge_counts = add_ring_edges(edge_counts, ring)
                    if next_edge_counts is None:
                        continue
                    next_counts = dict(counts)
                    next_counts[ring.size] = next_counts.get(ring.size, 0) + 1
                    if not counts_fit(next_counts, target_counts):
                        continue
                    stack.append((frozenset([*face_ids, next_id]), next_edge_counts, next_counts))


def grow_seed_face_sets(cups: list[Cup], ring_by_id: dict[str, Ring], seed_mode: str) -> list[frozenset[str]]:
    """Build cage-grow seeds from complete cups by default, or rings for comparison."""
    mode = seed_mode.lower().strip()
    if mode not in {"cup", "ring", "auto"}:
        raise ValueError(f"Unsupported cage.seed_mode: {seed_mode}")

    cup_seeds: list[frozenset[str]] = []
    seen: set[frozenset[str]] = set()
    for cup in cups:
        face_ids = frozenset(cup.rings)
        if face_ids in seen:
            continue
        if face_ids and all(ring_id in ring_by_id for ring_id in face_ids):
            cup_seeds.append(face_ids)
            seen.add(face_ids)

    if mode == "auto":
        mode = "cup" if cup_seeds else "ring"
    if mode == "cup":
        return sorted(cup_seeds, key=lambda item: (min(item), len(item), sorted(item)))
    return [frozenset([ring_id]) for ring_id in sorted(ring_by_id)]


def build_seed_index_by_anchor(seed_face_sets: list[frozenset[str]]) -> dict[str, list[tuple[int, frozenset[str]]]]:
    """Index cup seeds by their smallest ring id for fast subset checks."""
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
    """Prune duplicate cup-seed growth once an earlier seed is contained."""
    for ring_id in face_ids:
        for earlier_index, earlier in seed_index_by_anchor.get(ring_id, []):
            if earlier_index >= seed_index:
                continue
            if earlier <= face_ids:
                return True
    return False


def ordered_boundary_candidates(
    face_ids: frozenset[str],
    edge_counts: dict[tuple[int, int], int],
    counts: dict[int, int],
    target_counts: dict[int, int],
    edge_to_ring_ids: dict[tuple[int, int], set[str]],
    ring_by_id: dict[str, Ring],
    rank: dict[str, int],
    seed_rank: int,
) -> list[str]:
    """Choose the most constrained boundary edge and return addable rings."""
    best: list[str] | None = None
    for edge, edge_count in edge_counts.items():
        if edge_count != 1:
            continue
        candidates = []
        for ring_id in edge_to_ring_ids.get(edge, set()):
            if ring_id in face_ids or rank[ring_id] < seed_rank:
                continue
            ring = ring_by_id[ring_id]
            if counts.get(ring.size, 0) + 1 > target_counts.get(ring.size, 0):
                continue
            if can_add_ring(edge_counts, ring):
                candidates.append(ring_id)
        if not candidates:
            return []
        candidates.sort(key=rank.__getitem__)
        if best is None or len(candidates) < len(best):
            best = candidates
    return best or []


def can_add_ring(edge_counts: dict[tuple[int, int], int], ring: Ring) -> bool:
    """Check whether a face can be added without any edge exceeding degree two."""
    return all(edge_counts.get(edge, 0) < 2 for edge in ring.edges)


def add_ring_edges(edge_counts: dict[tuple[int, int], int], ring: Ring) -> dict[tuple[int, int], int] | None:
    """Return updated edge counts, or None if the face overuses an edge."""
    updated = dict(edge_counts)
    for edge in ring.edges:
        count = updated.get(edge, 0) + 1
        if count > 2:
            return None
        updated[edge] = count
    return updated


def fast_pair_cage_candidates(
    cup_face_sets: list[frozenset[str]],
    ring_by_id: dict[str, Ring],
    targets: dict[str, dict[int, int]],
):
    """Yield cage candidates from cup pairs with matching open boundaries."""
    buckets: dict[frozenset[tuple[int, int]], list[tuple[int, frozenset[str], dict[int, int]]]] = defaultdict(list)
    for idx, face_ids in enumerate(cup_face_sets):
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
    for cage_type, target_counts in targets.items():
        if counts_match(counts, target_counts):
            return cage_type
    return None


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
    """Validate a closed shell with edge counts and Euler characteristic."""
    edge_counts: dict[tuple[int, int], int] = defaultdict(int)
    nodes = set()
    for ring in rings:
        nodes.update(ring.nodes)
        for edge in ring.edges:
            edge_counts[edge] += 1
    if not edge_counts or any(count != 2 for count in edge_counts.values()):
        return False
    return len(nodes) - len(edge_counts) + len(rings) == 2




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


def unwrap_connected_nodes(frame: Frame, nodes: list[int], edges: list[tuple[int, int]]) -> dict[int, np.ndarray]:
    """Unwrap all cage oxygen nodes through the face-edge graph."""
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


def assigned_guests(
    frame: Frame,
    guests: list[Guest],
    center: np.ndarray,
    rings: list[Ring],
    unwrapped_nodes: dict[int, np.ndarray],
    radius_nm: float,
    occupancy_mode: str,
) -> tuple[str, ...]:
    """Assign guest residues by point-in-polyhedron or center-distance mode."""
    mode = occupancy_mode.lower().strip()
    if mode not in {"polyhedron", "center", "auto"}:
        raise ValueError(f"Unsupported cage.occupancy_mode: {occupancy_mode}")

    assigned = []
    shell_radius = max(float(np.linalg.norm(pos - center)) for pos in unwrapped_nodes.values()) if unwrapped_nodes else radius_nm
    triangles = triangulate_cage_faces(rings, unwrapped_nodes, center) if mode in {"polyhedron", "auto"} else []
    for guest in guests:
        raw_center = guest_center(frame, guest)
        local_center = center + minimum_image(raw_center - center, frame.box)
        center_hit = distance(raw_center, center, frame.box) <= radius_nm
        poly_hit = False
        if triangles and float(np.linalg.norm(local_center - center)) <= shell_radius + 0.25:
            poly_hit = point_in_polyhedron(local_center, triangles)
        if (mode == "polyhedron" and poly_hit) or (mode == "center" and center_hit) or (mode == "auto" and (poly_hit or center_hit)):
            assigned.append(f"{guest.resname}{guest.resid}")
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


def triangle_solid_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Return the signed solid angle of a triangle at the origin."""
    la = float(np.linalg.norm(a))
    lb = float(np.linalg.norm(b))
    lc = float(np.linalg.norm(c))
    if min(la, lb, lc) <= 1e-12:
        return 4.0 * np.pi
    numerator = float(np.dot(a, np.cross(b, c)))
    denominator = la * lb * lc + float(np.dot(a, b)) * lc + float(np.dot(b, c)) * la + float(np.dot(c, a)) * lb
    return 2.0 * float(np.arctan2(numerator, denominator))


def guest_center(frame: Frame, guest: Guest) -> np.ndarray:
    """Use the configured center atom when available, otherwise geometry center."""
    if guest.center_atom is not None:
        return frame.atoms[guest.center_atom].xyz
    return np.mean([frame.atoms[idx].xyz for idx in guest.atoms], axis=0)


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


