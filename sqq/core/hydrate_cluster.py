from __future__ import annotations

"""Fast face-graph classification of hydrate clusters and phase domains."""

from collections import Counter, defaultdict, deque
from dataclasses import dataclass, replace
from itertools import combinations
import re

import numpy as np

from ..models import Cage, Frame, HydrateCluster, HydrateDomain, HydrateMotif, Ring
from .geometry import unwrap_connected_nodes
from .pbc import minimum_image


PHASE_TYPES = ("sI", "sII", "sH")
PHASE_ORDER = {name: index for index, name in enumerate(PHASE_TYPES)}
RING_ID_PATTERN = re.compile(r"^ring(\d+)_")

# Ideal first-shell fingerprints keyed by neighbor type and face size.
PHASE_TEMPLATES: dict[str, dict[str, Counter[tuple[str, int]]]] = {
    "sI": {
        "512": Counter({("51262", 5): 12}),
        "51262": Counter({("512", 5): 4, ("51262", 5): 8, ("51262", 6): 2}),
    },
    "sII": {
        "512": Counter({("512", 5): 6, ("51264", 5): 6}),
        "51264": Counter({("512", 5): 12, ("51264", 6): 4}),
    },
    "sH": {
        "512": Counter({("512", 5): 4, ("435663", 5): 4, ("51268", 5): 4}),
        "435663": Counter({("435663", 4): 3, ("512", 5): 6, ("51268", 6): 3}),
        "51268": Counter({("512", 5): 12, ("435663", 6): 6, ("51268", 6): 2}),
    },
}
STRICT_COUNT_TOLERANCE = 1
EXPANSION_COUNT_TOLERANCE = 1
EXPANSION_MIN_PHASE_CONTACTS = 2
SPATIAL_CORE_MIN_COVERAGE = 0.50
SPATIAL_CORE_MIN_PURITY = 0.50
SPATIAL_CORE_MIN_SCORE = 0.55
SPATIAL_CORE_MIN_MEAN_SCORE = 0.60
SPATIAL_CORE_MIN_SIZE = 3
PHASE_CORE_EDGE = {
    "sI": ("51262", "51262", 6),
    "sII": ("51264", "51264", 6),
    "sH": ("435663", "51268", 6),
}
PHASE_CORE_ANCHOR_TYPES = {
    "sI": {"51262"},
    "sII": {"51264"},
    "sH": {"435663", "51268"},
}


@dataclass(frozen=True)
class PhaseSeed:
    """Internal phase evidence; seeds are deliberately not exported."""

    hydrate_type: str
    anchor_indexes: tuple[int, ...]
    cage_indexes: tuple[int, ...]
    shared_face_ids: tuple[str, ...]
    cluster_id: str
    source: str = "strict"


@dataclass(frozen=True)
class DomainSpec:
    """Internal domain component before stable output ids are assigned."""

    cluster_id: str
    hydrate_type: str
    cage_indexes: tuple[int, ...]
    seed_cage_indexes: tuple[int, ...]
    seed_count: int
    boundary_indexes: tuple[int, ...] = ()


def find_hydrate_clusters(
    cages: list[Cage],
    min_cage: int = 2,
    ring_sizes: dict[str, int] | None = None,
) -> tuple[list[HydrateCluster], tuple[str, ...]]:
    """Compatibility wrapper returning clusters and sub-threshold cage ids."""
    clusters, _, _, isolated = analyze_hydrate_clusters(
        cages,
        min_cage=min_cage,
        ring_sizes=ring_sizes,
    )
    return clusters, isolated


def analyze_hydrate_clusters(
    cages: list[Cage],
    min_cage: int = 2,
    ring_sizes: dict[str, int] | None = None,
    frame: Frame | None = None,
    rings_by_id: dict[str, Ring] | None = None,
    face_geometries: dict[str, tuple[np.ndarray, np.ndarray]] | None = None,
) -> tuple[list[HydrateCluster], list[HydrateMotif], list[HydrateDomain], tuple[str, ...]]:
    """Build clusters, identify strict and spatial evidence, and expand phase domains.

    Motifs are internal implementation details. The second return value remains
    for API compatibility and is always empty.
    """
    if min_cage < 1:
        raise ValueError("hydrate_cluster.min_cage / --cluster-min-cage must be at least 1.")
    if not cages:
        return [], [], [], ()

    if face_geometries is None and frame is not None and rings_by_id is not None:
        face_geometries = build_ring_face_geometries(frame, rings_by_id)
    adjacency, shared_faces, face_sizes = build_cage_graph(
        cages,
        ring_sizes,
        face_geometries=face_geometries,
        box=None if frame is None else frame.box,
    )
    signatures, evidence = build_local_fingerprints(
        cages,
        adjacency,
        shared_faces,
        face_sizes,
    )

    clusters: list[HydrateCluster] = []
    cluster_components: dict[str, tuple[int, ...]] = {}
    isolated: list[str] = []
    for component in connected_components(adjacency):
        ordered = tuple(sorted(component))
        if len(ordered) < min_cage:
            isolated.extend(cages[index].object_id for index in ordered)
            continue
        cluster = build_cluster(len(clusters) + 1, ordered, cages, shared_faces)
        clusters.append(cluster)
        cluster_components[cluster.object_id] = ordered

    all_seeds: list[PhaseSeed] = []
    accepted_by_cluster: dict[str, dict[str, set[int]]] = {}
    seed_members_by_cluster: dict[str, dict[str, set[int]]] = {}
    phase_edges_by_cluster: dict[str, dict[str, set[tuple[int, int]]]] = {}

    for cluster in clusters:
        component = cluster_components[cluster.object_id]
        seeds = find_strict_phase_seeds(
            cluster.object_id,
            component,
            cages,
            adjacency,
            shared_faces,
            face_sizes,
            signatures,
            evidence,
        )
        seeds.extend(
            find_spatial_phase_cores(
                cluster.object_id,
                component,
                cages,
                adjacency,
                shared_faces,
                face_sizes,
                signatures,
            )
        )
        all_seeds.extend(seeds)
        seeds_by_phase = {
            phase: [seed for seed in seeds if seed.hydrate_type == phase]
            for phase in PHASE_TYPES
        }
        seed_members = {
            phase: {
                index
                for seed in seeds_by_phase[phase]
                for index in seed.cage_indexes
            }
            for phase in PHASE_TYPES
        }
        phase_edges = {
            phase: {
                edge_key(left, right)
                for seed in seeds_by_phase[phase]
                for left, right in combinations(seed.cage_indexes, 2)
                if right in adjacency[left]
                and (
                    seed.source == "strict"
                    or phase_edge_compatible(
                        phase,
                        left,
                        right,
                        cages,
                        shared_faces,
                        face_sizes,
                        set(),
                    )
                )
            }
            for phase in PHASE_TYPES
        }
        accepted = {
            phase: expand_phase_from_seeds(
                phase,
                set(component),
                seed_members[phase],
                cages,
                adjacency,
                shared_faces,
                face_sizes,
                signatures,
                phase_edges[phase],
            )
            for phase in PHASE_TYPES
        }
        accepted_by_cluster[cluster.object_id] = accepted
        seed_members_by_cluster[cluster.object_id] = seed_members
        phase_edges_by_cluster[cluster.object_id] = phase_edges

    # Evaluate each phase independently before resolving overlaps.
    claims: dict[int, set[str]] = defaultdict(set)
    for cluster in clusters:
        for phase, indexes in accepted_by_cluster[cluster.object_id].items():
            for index in indexes:
                claims[index].add(phase)

    specs = build_domain_specs(
        clusters,
        cluster_components,
        cages,
        adjacency,
        shared_faces,
        face_sizes,
        claims,
        all_seeds,
        seed_members_by_cluster,
        phase_edges_by_cluster,
    )

    domain_members = {index for spec in specs for index in spec.cage_indexes}
    boundary_indexes, ambiguous_indexes = resolve_boundary_cages(
        clusters,
        cluster_components,
        claims,
        adjacency,
        domain_members,
    )
    specs = attach_domain_boundaries(
        specs,
        boundary_indexes,
        adjacency,
    )
    domains = materialize_domains(specs, cages)
    clusters = enrich_clusters(
        clusters,
        cluster_components,
        domains,
        cages,
        boundary_indexes,
        ambiguous_indexes,
    )
    return clusters, [], domains, tuple(isolated)


def build_cage_graph(
    cages: list[Cage],
    ring_sizes: dict[str, int] | None = None,
    face_geometries: dict[str, tuple[np.ndarray, np.ndarray]] | None = None,
    box: np.ndarray | None = None,
) -> tuple[dict[int, set[int]], dict[tuple[int, int], set[str]], dict[str, int]]:
    """Build a cage graph with at most one physical pair per shared face."""
    adjacency: dict[int, set[int]] = {index: set() for index in range(len(cages))}
    shared_faces: dict[tuple[int, int], set[str]] = defaultdict(set)
    ring_to_cages: dict[str, list[int]] = defaultdict(list)
    for index, cage in enumerate(cages):
        for ring_id in cage.rings:
            ring_to_cages[ring_id].append(index)

    for ring_id, indexes in ring_to_cages.items():
        pair = resolve_shared_face_pair(
            ring_id,
            sorted(set(indexes)),
            cages,
            face_geometries,
            box,
        )
        if pair is None:
            continue
        left, right = pair
        adjacency[left].add(right)
        adjacency[right].add(left)
        shared_faces[pair].add(ring_id)

    resolved_sizes = dict(ring_sizes or {})
    for ring_id in ring_to_cages:
        if ring_id in resolved_sizes:
            continue
        match = RING_ID_PATTERN.match(ring_id)
        if match:
            resolved_sizes[ring_id] = int(match.group(1))
    return adjacency, dict(shared_faces), resolved_sizes


def build_ring_face_geometries(
    frame: Frame,
    rings_by_id: dict[str, Ring],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return locally unwrapped ring centers and least-squares plane normals."""
    geometries: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ring_id, ring in rings_by_id.items():
        unwrapped = unwrap_connected_nodes(frame, list(ring.nodes), list(ring.edges))
        points = np.asarray([unwrapped[node] for node in ring.nodes], dtype=float)
        if len(points) < 3:
            continue
        center = np.mean(points, axis=0)
        centered = points - center
        try:
            _, _, axes = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        normal = np.asarray(axes[-1], dtype=float)
        norm = float(np.linalg.norm(normal))
        if norm <= 1.0e-12:
            continue
        geometries[ring_id] = (center, normal / norm)
    return geometries


def resolve_shared_face_pair(
    ring_id: str,
    indexes: list[int],
    cages: list[Cage],
    face_geometries: dict[str, tuple[np.ndarray, np.ndarray]] | None,
    box: np.ndarray | None,
) -> tuple[int, int] | None:
    """Select one cage on each side of a ring plane.

    A physical face can separate at most two cages. When geometry is absent,
    only an already unambiguous two-cage face is accepted.
    """
    if len(indexes) < 2:
        return None
    geometry = None if face_geometries is None else face_geometries.get(ring_id)
    if geometry is None:
        return edge_key(indexes[0], indexes[1]) if len(indexes) == 2 else None

    center, normal = geometry
    by_side: dict[int, list[tuple[tuple[float, float, float, int], int]]] = {
        -1: [],
        1: [],
    }
    for index in indexes:
        delta = minimum_image(np.asarray(cages[index].center) - center, box)
        distance = float(np.linalg.norm(delta))
        if distance <= 1.0e-12:
            continue
        signed = float(np.dot(delta, normal))
        if abs(signed) <= 1.0e-5:
            continue
        side = 1 if signed > 0.0 else -1
        alignment = abs(signed) / distance
        lateral = float(np.sqrt(max(distance * distance - signed * signed, 0.0)))
        score = (-alignment, lateral, -abs(signed), index)
        by_side[side].append((score, index))
    if not by_side[-1] or not by_side[1]:
        return None
    left = min(by_side[-1])[1]
    right = min(by_side[1])[1]
    return edge_key(left, right)


def build_local_fingerprints(
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
) -> tuple[
    dict[int, Counter[tuple[str, int]]],
    dict[int, dict[tuple[str, int], tuple[tuple[int, str], ...]]],
]:
    """Compute every first-shell labelled fingerprint in one graph pass."""
    signatures: dict[int, Counter[tuple[str, int]]] = {
        index: Counter() for index in range(len(cages))
    }
    raw_evidence: dict[int, dict[tuple[str, int], list[tuple[int, str]]]] = {
        index: defaultdict(list) for index in range(len(cages))
    }
    for (left, right), ring_ids in shared_faces.items():
        for ring_id in ring_ids:
            face_size = face_sizes.get(ring_id)
            if face_size is None:
                continue
            left_key = (cages[right].cage_type, face_size)
            right_key = (cages[left].cage_type, face_size)
            signatures[left][left_key] += 1
            signatures[right][right_key] += 1
            raw_evidence[left][left_key].append((right, ring_id))
            raw_evidence[right][right_key].append((left, ring_id))
    evidence = {
        index: {
            key: tuple(sorted(items))
            for key, items in by_key.items()
        }
        for index, by_key in raw_evidence.items()
    }
    return signatures, evidence


def find_strict_phase_seeds(
    cluster_id: str,
    component: tuple[int, ...],
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
    signatures: dict[int, Counter[tuple[str, int]]],
    evidence: dict[int, dict[tuple[str, int], tuple[tuple[int, str], ...]]],
) -> list[PhaseSeed]:
    """Find strict local phase fingerprints plus composite sH evidence."""
    found: dict[tuple[str, tuple[int, ...], tuple[int, ...]], PhaseSeed] = {}
    for center in component:
        cage_type = cages[center].cage_type
        for phase in PHASE_TYPES:
            expected = PHASE_TEMPLATES[phase].get(cage_type)
            if expected is None or not strict_signature_match(signatures[center], expected):
                continue
            members = {
                center,
                *(
                    neighbor
                    for key, items in evidence[center].items()
                    if key in expected
                    for neighbor, _ in items
                ),
            }
            faces = {
                ring_id
                for key, items in evidence[center].items()
                if key in expected
                for _, ring_id in items
            }
            seed = PhaseSeed(
                hydrate_type=phase,
                anchor_indexes=(center,),
                cage_indexes=tuple(sorted(members)),
                shared_face_ids=tuple(sorted(faces)),
                cluster_id=cluster_id,
            )
            found[(phase, seed.anchor_indexes, seed.cage_indexes)] = seed

    for seed in find_s_h_seeds(
        cluster_id,
        component,
        cages,
        adjacency,
        shared_faces,
        face_sizes,
    ):
        found[(seed.hydrate_type, seed.anchor_indexes, seed.cage_indexes)] = seed
    return [found[key] for key in sorted(found, key=lambda item: (PHASE_ORDER[item[0]], item[1], item[2]))]


def strict_signature_match(
    observed: Counter[tuple[str, int]],
    expected: Counter[tuple[str, int]],
) -> bool:
    """Return whether a fingerprint is a complete phase seed signature."""
    if not all(observed[key] > 0 for key in expected):
        return False
    if any(key not in expected for key in observed):
        return False
    return all(
        abs(observed[key] - expected_count) <= STRICT_COUNT_TOLERANCE
        for key, expected_count in expected.items()
    )


def phase_signature_support(
    observed: Counter[tuple[str, int]],
    expected: Counter[tuple[str, int]],
) -> tuple[float, float, float, int]:
    """Score how much of a local phase fingerprint is present and phase-pure."""
    matched = sum(
        min(observed.get(key, 0), expected_count)
        for key, expected_count in expected.items()
    )
    expected_total = sum(expected.values())
    observed_total = sum(observed.values())
    coverage = matched / expected_total if expected_total else 0.0
    purity = matched / observed_total if observed_total else 0.0
    score = (
        2.0 * coverage * purity / (coverage + purity)
        if coverage > 0.0 and purity > 0.0
        else 0.0
    )
    return score, coverage, purity, matched


def find_spatial_phase_cores(
    cluster_id: str,
    component: tuple[int, ...],
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
    signatures: dict[int, Counter[tuple[str, int]]],
) -> list[PhaseSeed]:
    """Find coherent per-frame phase cores without requiring one perfect cage."""
    output: list[PhaseSeed] = []
    component_set = set(component)
    for phase in PHASE_TYPES:
        scores: dict[int, float] = {}
        for index in component:
            expected = PHASE_TEMPLATES[phase].get(cages[index].cage_type)
            if expected is None:
                continue
            score, coverage, purity, _ = phase_signature_support(
                signatures[index],
                expected,
            )
            if (
                coverage >= SPATIAL_CORE_MIN_COVERAGE
                and purity >= SPATIAL_CORE_MIN_PURITY
                and score >= SPATIAL_CORE_MIN_SCORE
            ):
                scores[index] = score

        candidate_set = set(scores)
        phase_adjacency = {
            index: {
                neighbor
                for neighbor in adjacency[index].intersection(
                    candidate_set,
                    component_set,
                )
                if phase_edge_compatible(
                    phase,
                    index,
                    neighbor,
                    cages,
                    shared_faces,
                    face_sizes,
                    set(),
                )
            }
            for index in candidate_set
        }
        core = graph_two_core(phase_adjacency)
        if not core:
            continue
        for indexes in induced_components(phase_adjacency, core):
            members = tuple(sorted(indexes))
            if len(members) < SPATIAL_CORE_MIN_SIZE:
                continue
            mean_score = sum(scores[index] for index in members) / len(members)
            if mean_score < SPATIAL_CORE_MIN_MEAN_SCORE:
                continue
            if not has_phase_core_edge(
                phase,
                members,
                cages,
                shared_faces,
                face_sizes,
            ):
                continue
            anchors = tuple(
                index
                for index in members
                if cages[index].cage_type in PHASE_CORE_ANCHOR_TYPES[phase]
            )
            output.append(
                PhaseSeed(
                    hydrate_type=phase,
                    anchor_indexes=anchors,
                    cage_indexes=anchors,
                    shared_face_ids=internal_shared_faces(anchors, shared_faces),
                    cluster_id=cluster_id,
                    source="spatial_core",
                )
            )
    return output


def graph_two_core(adjacency: dict[int, set[int]]) -> set[int]:
    """Return the maximal subgraph whose nodes have at least two neighbors."""
    remaining = set(adjacency)
    degrees = {
        index: len(neighbors.intersection(remaining))
        for index, neighbors in adjacency.items()
    }
    queue = deque(sorted(index for index, degree in degrees.items() if degree < 2))
    while queue:
        index = queue.popleft()
        if index not in remaining:
            continue
        remaining.remove(index)
        for neighbor in adjacency[index].intersection(remaining):
            degrees[neighbor] -= 1
            if degrees[neighbor] == 1:
                queue.append(neighbor)
    return remaining


def has_phase_core_edge(
    phase: str,
    indexes: tuple[int, ...],
    cages: list[Cage],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
) -> bool:
    """Require the phase-defining large-cage connection in a spatial core."""
    left_type, right_type, required_size = PHASE_CORE_EDGE[phase]
    for left, right in combinations(indexes, 2):
        pair_types = (cages[left].cage_type, cages[right].cage_type)
        if not (
            pair_types == (left_type, right_type)
            or pair_types == (right_type, left_type)
        ):
            continue
        if pair_shares_face_size(
            left,
            right,
            required_size,
            shared_faces,
            face_sizes,
        ):
            return True
    return False


def partial_signature_match(
    observed: Counter[tuple[str, int]],
    expected: Counter[tuple[str, int]],
) -> bool:
    """Accept an incomplete internal fingerprint and ignore external contacts."""
    internal = Counter({key: count for key, count in observed.items() if key in expected})
    if not internal:
        return False
    return all(
        count <= expected[key] + EXPANSION_COUNT_TOLERANCE
        for key, count in internal.items()
    )


def find_s_h_seeds(
    cluster_id: str,
    component: tuple[int, ...],
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
) -> list[PhaseSeed]:
    """Find supplemental high-confidence sH composite seeds."""
    component_set = set(component)
    candidate_pairs: set[tuple[int, int]] = set()
    for small in component:
        if cages[small].cage_type != "512":
            continue
        large_neighbors = sorted(
            neighbor
            for neighbor in adjacency[small].intersection(component_set)
            if cages[neighbor].cage_type == "51268"
            and pair_shares_face_size(small, neighbor, 5, shared_faces, face_sizes)
        )
        candidate_pairs.update(edge_key(left, right) for left, right in combinations(large_neighbors, 2))

    output: list[PhaseSeed] = []
    for left, right in sorted(candidate_pairs):
        if right in adjacency[left]:
            continue
        common_small = {
            index
            for index in adjacency[left].intersection(adjacency[right], component_set)
            if cages[index].cage_type == "512"
            and pair_shares_face_size(left, index, 5, shared_faces, face_sizes)
            and pair_shares_face_size(right, index, 5, shared_faces, face_sizes)
        }
        if len(common_small) != 6:
            continue
        medium = {
            index
            for index in component
            if cages[index].cage_type == "435663"
            and adjacency[index].intersection({left, right})
            and adjacency[index].intersection(common_small)
        }
        if len(medium) != 6 or not has_s_h_bridge(left, right, medium, adjacency):
            continue
        members = tuple(sorted({left, right, *common_small, *medium}))
        output.append(
            PhaseSeed(
                hydrate_type="sH",
                anchor_indexes=(left, right),
                cage_indexes=members,
                shared_face_ids=internal_shared_faces(members, shared_faces),
                cluster_id=cluster_id,
            )
        )
    return output


def expand_phase_from_seeds(
    phase: str,
    component: set[int],
    seed_members: set[int],
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
    signatures: dict[int, Counter[tuple[str, int]]],
    seed_edges: set[tuple[int, int]],
) -> set[int]:
    """Expand one phase from strict seeds and spatial-core anchors."""
    if not seed_members:
        return set()
    accepted = set(seed_members)

    compatible_nodes = {
        index
        for index in component.difference(accepted)
        if (expected := PHASE_TEMPLATES[phase].get(cages[index].cage_type)) is not None
        and partial_signature_match(signatures[index], expected)
    }
    contact_counts: Counter[int] = Counter()
    queue = deque(sorted(accepted))
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current].intersection(compatible_nodes):
            if neighbor in accepted:
                continue
            if not phase_edge_compatible(
                phase,
                current,
                neighbor,
                cages,
                shared_faces,
                face_sizes,
                seed_edges,
            ):
                continue
            contact_counts[neighbor] += 1
            if contact_counts[neighbor] < EXPANSION_MIN_PHASE_CONTACTS:
                continue
            accepted.add(neighbor)
            compatible_nodes.remove(neighbor)
            queue.append(neighbor)
    return accepted


def phase_edge_compatible(
    phase: str,
    left: int,
    right: int,
    cages: list[Cage],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
    seed_edges: set[tuple[int, int]],
) -> bool:
    """Return whether both endpoints allow a labelled phase edge."""
    pair = edge_key(left, right)
    if pair in seed_edges:
        return True
    left_expected = PHASE_TEMPLATES[phase].get(cages[left].cage_type)
    right_expected = PHASE_TEMPLATES[phase].get(cages[right].cage_type)
    if left_expected is None or right_expected is None:
        return False
    return any(
        (cages[right].cage_type, face_sizes.get(ring_id)) in left_expected
        and (cages[left].cage_type, face_sizes.get(ring_id)) in right_expected
        for ring_id in shared_faces.get(pair, ())
    )


def build_domain_specs(
    clusters: list[HydrateCluster],
    cluster_components: dict[str, tuple[int, ...]],
    cages: list[Cage],
    adjacency: dict[int, set[int]],
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
    claims: dict[int, set[str]],
    seeds: list[PhaseSeed],
    seed_members_by_cluster: dict[str, dict[str, set[int]]],
    phase_edges_by_cluster: dict[str, dict[str, set[tuple[int, int]]]],
) -> list[DomainSpec]:
    """Create exclusive phase components while retaining internal seed evidence."""
    specs: list[DomainSpec] = []
    for cluster in clusters:
        component = cluster_components[cluster.object_id]
        for phase in PHASE_TYPES:
            exclusive = {
                index for index in component if claims.get(index, set()) == {phase}
            }
            if not exclusive:
                continue
            phase_adjacency = {
                index: {
                    neighbor
                    for neighbor in adjacency[index].intersection(exclusive)
                    if phase_edge_compatible(
                        phase,
                        index,
                        neighbor,
                        cages,
                        shared_faces,
                        face_sizes,
                        phase_edges_by_cluster[cluster.object_id][phase],
                    )
                }
                for index in exclusive
            }
            phase_seeds = [
                seed
                for seed in seeds
                if seed.cluster_id == cluster.object_id and seed.hydrate_type == phase
            ]
            for domain_component in induced_components(phase_adjacency, exclusive):
                domain_set = set(domain_component)
                contained_seeds = [
                    seed
                    for seed in phase_seeds
                    if set(seed.anchor_indexes).issubset(domain_set)
                ]
                if not contained_seeds:
                    continue
                seed_cages = domain_set.intersection(
                    seed_members_by_cluster[cluster.object_id][phase]
                )
                specs.append(
                    DomainSpec(
                        cluster_id=cluster.object_id,
                        hydrate_type=phase,
                        cage_indexes=tuple(sorted(domain_set)),
                        seed_cage_indexes=tuple(sorted(seed_cages)),
                        seed_count=len(contained_seeds),
                    )
                )
    return sorted(
        specs,
        key=lambda item: (
            item.cluster_id,
            PHASE_ORDER[item.hydrate_type],
            item.cage_indexes,
        ),
    )


def resolve_boundary_cages(
    clusters: list[HydrateCluster],
    cluster_components: dict[str, tuple[int, ...]],
    claims: dict[int, set[str]],
    adjacency: dict[int, set[int]],
    domain_members: set[int],
) -> tuple[set[int], set[int]]:
    """Return the external first cage layer and unresolved phase claims."""
    boundary: set[int] = set()
    ambiguous: set[int] = set()
    for cluster in clusters:
        component = cluster_components[cluster.object_id]
        component_set = set(component)
        for index in component:
            if index in domain_members:
                continue
            touches_domain = any(
                neighbor in domain_members
                for neighbor in adjacency[index].intersection(component_set)
            )
            if touches_domain:
                boundary.add(index)
            elif len(claims.get(index, set())) > 1:
                ambiguous.add(index)
    return boundary, ambiguous


def attach_domain_boundaries(
    specs: list[DomainSpec],
    boundary_indexes: set[int],
    adjacency: dict[int, set[int]],
) -> list[DomainSpec]:
    """Attach unique external boundary contacts to every domain."""
    output: list[DomainSpec] = []
    for spec in specs:
        members = set(spec.cage_indexes)
        external_boundary = boundary_indexes.difference(members)
        contacts = {
            neighbor
            for index in members
            for neighbor in adjacency[index].intersection(external_boundary)
        }
        output.append(replace(spec, boundary_indexes=tuple(sorted(contacts))))
    return output


def materialize_domains(specs: list[DomainSpec], cages: list[Cage]) -> list[HydrateDomain]:
    """Convert exclusive domain components into public models."""
    domains: list[HydrateDomain] = []
    for number, spec in enumerate(specs, start=1):
        waters = tuple(
            sorted({water for index in spec.cage_indexes for water in cages[index].waters})
        )
        guests = unique_in_order(
            guest for index in spec.cage_indexes for guest in cages[index].guest_ids
        )
        expanded_count = len(spec.cage_indexes) - len(spec.seed_cage_indexes)
        classified_fraction = (
            len(spec.cage_indexes)
            / (len(spec.cage_indexes) + len(spec.boundary_indexes))
            if spec.cage_indexes or spec.boundary_indexes
            else 0.0
        )
        domains.append(
            HydrateDomain(
                object_id=f"domain_{number:05d}",
                cluster_id=spec.cluster_id,
                hydrate_type=spec.hydrate_type,
                cage_ids=tuple(cages[index].object_id for index in spec.cage_indexes),
                motif_ids=(),
                waters=waters,
                guest_ids=guests,
                boundary_cage_ids=tuple(
                    cages[index].object_id for index in spec.boundary_indexes
                ),
                confidence=round(classified_fraction, 6),
                status="complete" if expanded_count == 0 else "expanded",
                seed_count=spec.seed_count,
                seed_cage_ids=tuple(
                    cages[index].object_id for index in spec.seed_cage_indexes
                ),
                classified_fraction=round(classified_fraction, 6),
            )
        )
    return domains


def enrich_clusters(
    clusters: list[HydrateCluster],
    cluster_components: dict[str, tuple[int, ...]],
    domains: list[HydrateDomain],
    cages: list[Cage],
    boundary_indexes: set[int],
    ambiguous_indexes: set[int],
) -> list[HydrateCluster]:
    """Partition each cluster into phase, boundary, ambiguous, and residual cages."""
    cage_index_by_id = {cage.object_id: index for index, cage in enumerate(cages)}
    output: list[HydrateCluster] = []
    for cluster in clusters:
        component = cluster_components[cluster.object_id]
        cluster_domains = [
            domain for domain in domains if domain.cluster_id == cluster.object_id
        ]
        classified_indexes = {
            cage_index_by_id[cage_id]
            for domain in cluster_domains
            for cage_id in domain.cage_ids
        }
        boundary = set(component).intersection(boundary_indexes)
        ambiguous = set(component).intersection(ambiguous_indexes).difference(boundary)
        unclassified = set(component).difference(
            classified_indexes,
            boundary,
            ambiguous,
        )
        domain_types = {domain.hydrate_type for domain in cluster_domains}
        if not domain_types:
            hydrate_type = "unclassified"
        elif len(domain_types) == 1:
            hydrate_type = next(iter(domain_types))
        else:
            hydrate_type = "mixed"
        type_counts = Counter(
            domain.hydrate_type
            for domain in cluster_domains
            for _ in domain.cage_ids
        )
        output.append(
            replace(
                cluster,
                hydrate_type=hydrate_type,
                motif_ids=(),
                domain_ids=tuple(domain.object_id for domain in cluster_domains),
                classified_cage_ids=tuple(
                    cages[index].object_id for index in component if index in classified_indexes
                ),
                unclassified_cage_ids=tuple(
                    cages[index].object_id for index in component if index in unclassified
                ),
                ambiguous_cage_ids=tuple(
                    cages[index].object_id for index in component if index in ambiguous
                ),
                boundary_cage_ids=tuple(
                    cages[index].object_id for index in component if index in boundary
                ),
                hydrate_type_counts=tuple(
                    (phase, type_counts[phase])
                    for phase in PHASE_TYPES
                    if type_counts[phase]
                ),
            )
        )
    return output


def pair_shares_face_size(
    left: int,
    right: int,
    size: int,
    shared_faces: dict[tuple[int, int], set[str]],
    face_sizes: dict[str, int],
) -> bool:
    """Return whether a cage pair shares a face of the requested size."""
    return any(
        face_sizes.get(ring_id) == size
        for ring_id in shared_faces.get(edge_key(left, right), ())
    )


def has_s_h_bridge(
    left: int,
    right: int,
    medium: set[int],
    adjacency: dict[int, set[int]],
) -> bool:
    """Require adjacent medium cages to bridge the separated sH anchors."""
    for first, second in combinations(sorted(medium), 2):
        if second not in adjacency[first]:
            continue
        if first in adjacency[left] and second in adjacency[right]:
            return True
        if second in adjacency[left] and first in adjacency[right]:
            return True
    return False


def internal_shared_faces(
    indexes: tuple[int, ...],
    shared_faces: dict[tuple[int, int], set[str]],
) -> tuple[str, ...]:
    """Collect graph-edge evidence internal to one seed."""
    result: set[str] = set()
    index_set = set(indexes)
    for (left, right), ring_ids in shared_faces.items():
        if left in index_set and right in index_set:
            result.update(ring_ids)
    return tuple(sorted(result))


def connected_components(adjacency: dict[int, set[int]]) -> list[tuple[int, ...]]:
    """Return deterministic connected components from an undirected graph."""
    return list(induced_components(adjacency, set(adjacency)))


def induced_components(
    adjacency: dict[int, set[int]],
    nodes: set[int],
) -> tuple[tuple[int, ...], ...]:
    """Return components of the graph induced by a selected node set."""
    seen: set[int] = set()
    components: list[tuple[int, ...]] = []
    for start in sorted(nodes):
        if start in seen:
            continue
        stack = [start]
        component: set[int] = set()
        seen.add(start)
        while stack:
            current = stack.pop()
            component.add(current)
            for neighbor in sorted(
                adjacency.get(current, set()).intersection(nodes),
                reverse=True,
            ):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(components)


def build_cluster(
    cluster_index: int,
    cage_indexes: tuple[int, ...],
    cages: list[Cage],
    shared_faces: dict[tuple[int, int], set[str]],
) -> HydrateCluster:
    """Create one cluster model from a connected cage component."""
    cage_ids = tuple(cages[index].object_id for index in cage_indexes)
    cage_types = tuple(cages[index].cage_type for index in cage_indexes)
    waters = tuple(
        sorted({water for index in cage_indexes for water in cages[index].waters})
    )
    guest_ids = unique_in_order(
        guest for index in cage_indexes for guest in cages[index].guest_ids
    )
    index_set = set(cage_indexes)
    face_edges: list[tuple[str, str, str]] = []
    for (left, right), ring_ids in shared_faces.items():
        if left not in index_set or right not in index_set:
            continue
        for ring_id in sorted(ring_ids):
            face_edges.append(
                (cages[left].object_id, cages[right].object_id, ring_id)
            )
    return HydrateCluster(
        object_id=f"cluster_{cluster_index:05d}",
        cage_ids=cage_ids,
        cage_types=cage_types,
        waters=waters,
        guest_ids=guest_ids,
        shared_faces=tuple(face_edges),
    )


def edge_key(left: int, right: int) -> tuple[int, int]:
    """Return the canonical key for one undirected cage edge."""
    return (left, right) if left < right else (right, left)


def unique_in_order(values) -> tuple[str, ...]:
    """Deduplicate values without losing first-seen order."""
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return tuple(result)
