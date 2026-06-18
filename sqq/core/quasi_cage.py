from __future__ import annotations

"""Half-cage and quasi-cage search from layered ring-face patches."""

from collections import Counter, defaultdict
import numpy as np

from ..models import CagePatch, Frame, Ring
from .geometry import unwrap_connected_nodes
from .pbc import distance


SUPERSCRIPT_DIGITS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
SUBSCRIPT_DIGITS = str.maketrans("0123456789-", "₀₁₂₃₄₅₆₇₈₉₋")


def find_cage_patches(
    frame: Frame,
    rings: dict[int, list[Ring]],
    enabled: bool = False,
    base_sizes: list[int] | None = None,
    side_sizes: list[int] | None = None,
    max_combinations_per_base: int = 50000,
    max_layers: int = 1,
    max_rings_per_layer: int = 6,
    max_layer_states_per_seed: int = 200,
    max_candidates_per_edge: int = 4,
    max_layer_candidates: int = 24,
) -> tuple[list[CagePatch], list[CagePatch]]:
    """Find standard half-cages and nonstandard open quasi-cage patches."""
    if not enabled:
        return [], []

    base_allowed = set(base_sizes or [5, 6])
    side_allowed = set(side_sizes or [5, 6])
    all_rings = [ring for group in rings.values() for ring in group]
    ring_by_id = {ring.object_id: ring for ring in all_rings}
    ring_centers = build_ring_centers(frame, all_rings)
    # Shared-edge lookup is the hard topology filter; centers only rank candidates.
    edge_to_rings = build_edge_to_rings(all_rings)
    half_cages: list[CagePatch] = []
    quasi_cages: list[CagePatch] = []
    seen_half: set[frozenset[str]] = set()
    seen_quasi: set[frozenset[str]] = set()
    type_counts: dict[str, int] = defaultdict(int)

    for base in all_rings:
        if base.size not in base_allowed:
            continue
        # L1 candidates must share the corresponding base-ring edge.
        candidate_lists = side_ring_candidate_lists(
            base,
            edge_to_rings,
            side_allowed,
            ring_centers,
            frame.box,
            max_candidates_per_edge,
        )
        if not candidate_lists:
            continue

        for side_rings in iter_closed_side_walls(base, candidate_lists, max_combinations_per_base):
            side_ids = [ring.object_id for ring in side_rings]
            base_patch_rings = [base, *side_rings]
            if not has_expected_nodes(base, side_rings):
                continue

            l1_sequence = canonical_sequence([ring.size for ring in side_rings])
            # Classify the closed L1 side wall first; standard forms become half_cage.
            add_classified_patch(
                frame,
                ring_by_id,
                base.size,
                [base.object_id, *side_ids],
                (l1_sequence,),
                half_cages,
                quasi_cages,
                seen_half,
                seen_quasi,
                type_counts,
            )
            if half_cage_type(base.size, (l1_sequence, "6")) is not None:
                l2_rings = next_layer_candidates(
                    base_patch_rings,
                    side_rings,
                    edge_to_rings,
                    side_allowed,
                    ring_centers,
                    frame.box,
                    lower_rings=[base],
                    max_candidates=max_layer_candidates,
                )
                l2_sixes = [ring for ring in l2_rings if ring.size == 6]
            else:
                l2_sixes = []
            for l2_ring in l2_sixes:
                # Only the standard 6r + 5^6 + 6^1 form may bypass the
                # max_layers=1 strict-L1 quasi-cage default.
                add_classified_patch(
                    frame,
                    ring_by_id,
                    base.size,
                    [base.object_id, *side_ids, l2_ring.object_id],
                    (l1_sequence, "6"),
                    half_cages,
                    quasi_cages,
                    seen_half,
                    seen_quasi,
                    type_counts,
                )

            # L2/L3 grow from exposed frontier edges and may remain dangling.
            add_layered_quasi_cages(
                frame,
                ring_by_id,
                edge_to_rings,
                side_allowed,
                ring_centers,
                base.size,
                [base.object_id, *side_ids],
                side_ids,
                (l1_sequence,),
                max_layers=max_layers,
                max_rings_per_layer=max_rings_per_layer,
                max_layer_states_per_seed=max_layer_states_per_seed,
                max_layer_candidates=max_layer_candidates,
                half_cages=half_cages,
                quasi_cages=quasi_cages,
                seen_quasi=seen_quasi,
                seen_half=seen_half,
                type_counts=type_counts,
            )

    half_cages = remove_subset_patches(half_cages)
    half_sets = {frozenset(patch.rings) for patch in half_cages}
    quasi_cages = [patch for patch in quasi_cages if frozenset(patch.rings) not in half_sets]
    return sorted(half_cages, key=lambda item: (item.patch_type, item.object_id)), sorted(quasi_cages, key=lambda item: (item.patch_type, item.object_id))


def build_ring_centers(frame: Frame, rings: list[Ring]) -> dict[str, np.ndarray]:
    """Compute one locally unwrapped oxygen centroid per ring."""
    return {ring.object_id: ring_center(frame, ring) for ring in rings}


def ring_center(frame: Frame, ring: Ring) -> np.ndarray:
    """Compute a ring centroid for spatial candidate pruning."""
    unwrapped = unwrap_connected_nodes(frame, list(ring.nodes), list(ring.edges))
    return np.mean([unwrapped[node] for node in ring.nodes], axis=0)


def nearest_rings(
    candidates: list[Ring],
    references: list[Ring],
    ring_centers: dict[str, np.ndarray],
    box: np.ndarray | None,
    limit: int,
) -> list[Ring]:
    """Keep the nearest topological candidates by ring-center distance."""
    if limit <= 0 or len(candidates) <= limit:
        return candidates

    def score(ring: Ring) -> tuple[float, str]:
        center = ring_centers[ring.object_id]
        nearest = min(distance(center, ring_centers[reference.object_id], box) for reference in references)
        return nearest, ring.object_id

    return sorted(candidates, key=score)[:limit]


def build_edge_to_rings(rings: list[Ring]) -> dict[tuple[int, int], list[Ring]]:
    """Map each graph edge to the rings that use it."""
    edge_to_rings: dict[tuple[int, int], list[Ring]] = defaultdict(list)
    for ring in rings:
        for edge in ring.edges:
            edge_to_rings[edge].append(ring)
    return dict(edge_to_rings)


def side_ring_candidate_lists(
    base: Ring,
    edge_to_rings: dict[tuple[int, int], list[Ring]],
    allowed_sizes: set[int],
    ring_centers: dict[str, np.ndarray],
    box: np.ndarray | None,
    max_candidates_per_edge: int,
) -> list[list[Ring]]:
    """Find L1 candidates by base-edge lookup, then order them geometrically."""
    candidate_lists: list[list[Ring]] = []
    for edge in ordered_edges(base.nodes):
        candidates = [
            ring
            for ring in shared_edge_ring_candidates(
                [edge],
                edge_to_rings,
                exclude_ids={base.object_id},
                allowed_sizes=allowed_sizes,
            )
            if (base.edges & ring.edges) == {edge}
        ]
        candidates = nearest_rings(candidates, [base], ring_centers, box, max_candidates_per_edge)
        if not candidates:
            return []
        candidate_lists.append(candidates)
    return candidate_lists


def shared_edge_ring_candidates(
    edges: set[tuple[int, int]] | list[tuple[int, int]],
    edge_to_rings: dict[tuple[int, int], list[Ring]],
    exclude_ids: set[str],
    allowed_sizes: set[int],
    blocked_edges: set[tuple[int, int]] | None = None,
) -> list[Ring]:
    """Collect ring candidates by edge-to-ring reverse lookup."""
    blocked_edges = blocked_edges or set()
    found: dict[str, Ring] = {}
    for edge in edges:
        for ring in edge_to_rings.get(edge, []):
            if ring.object_id in exclude_ids or ring.size not in allowed_sizes:
                continue
            if blocked_edges and ring.edges & blocked_edges:
                continue
            found[ring.object_id] = ring
    return [found[ring_id] for ring_id in sorted(found)]


def ordered_edges(nodes: tuple[int, ...]) -> list[tuple[int, int]]:
    """Return ring edges in cyclic order with sorted edge endpoints."""
    edges = []
    for idx, a in enumerate(nodes):
        b = nodes[(idx + 1) % len(nodes)]
        edges.append((a, b) if a < b else (b, a))
    return edges


def side_wall_closed(base: Ring, side_rings: list[Ring]) -> bool:
    """Check that neighboring side rings share exactly one non-base wall edge."""
    base_edges = base.edges
    n = len(side_rings)
    for idx in range(n):
        current_ring = side_rings[idx]
        next_ring = side_rings[(idx + 1) % n]
        shared_nonbase = (current_ring.edges & next_ring.edges) - base_edges
        if len(shared_nonbase) != 1:
            return False
    return True


def iter_closed_side_walls(base: Ring, candidate_lists: list[list[Ring]], max_states: int):
    """Yield side-ring walls while pruning non-adjacent choices early."""
    if not candidate_lists or max_states <= 0:
        return
    base_edges = base.edges
    emitted = 0
    stack: list[tuple[int, list[Ring], set[str]]] = []
    for ring in reversed(candidate_lists[0]):
        stack.append((1, [ring], {ring.object_id}))

    while stack and emitted < max_states:
        edge_index, selected, used = stack.pop()
        if edge_index == len(candidate_lists):
            if side_rings_touch(selected[-1], selected[0], base_edges):
                emitted += 1
                yield list(selected)
            continue
        previous = selected[-1]
        for ring in reversed(candidate_lists[edge_index]):
            if ring.object_id in used:
                continue
            if not side_rings_touch(previous, ring, base_edges):
                continue
            stack.append((edge_index + 1, [*selected, ring], {*used, ring.object_id}))


def side_rings_touch(left: Ring, right: Ring, base_edges: frozenset[tuple[int, int]]) -> bool:
    """Whether two neighboring side rings share one non-base wall edge."""
    return len((left.edges & right.edges) - base_edges) == 1


def has_expected_nodes(base: Ring, side_rings: list[Ring]) -> bool:
    """Reject overlapped or shifted side-ring choices."""
    unique_nodes = {node for ring in [base, *side_rings] for node in ring.nodes}
    expected_nodes = sum(ring.size for ring in side_rings) - 2 * base.size
    return len(unique_nodes) == expected_nodes


def next_layer_candidates(
    patch_rings: list[Ring],
    frontier_rings: list[Ring],
    edge_to_rings: dict[tuple[int, int], list[Ring]],
    allowed_sizes: set[int],
    ring_centers: dict[str, np.ndarray],
    box: np.ndarray | None,
    lower_rings: list[Ring] | None = None,
    max_candidates: int = 24,
) -> list[Ring]:
    """Collect rings attached to frontier boundary edges by reverse lookup."""
    patch_ids = {ring.object_id for ring in patch_rings}
    lower_edges = {edge for ring in lower_rings or [] for edge in ring.edges}
    growth_edges = growth_edges_for_frontier(patch_rings, frontier_rings)
    candidates = shared_edge_ring_candidates(
        growth_edges,
        edge_to_rings,
        exclude_ids=patch_ids,
        allowed_sizes=allowed_sizes,
        blocked_edges=lower_edges,
    )
    return nearest_rings(candidates, frontier_rings, ring_centers, box, max_candidates)


def growth_edges_for_frontier(patch_rings: list[Ring], frontier_rings: list[Ring]) -> set[tuple[int, int]]:
    """Return exposed frontier edges that can grow the next quasi-cage layer."""
    edge_counts: Counter[tuple[int, int]] = Counter()
    for ring in patch_rings:
        edge_counts.update(ring.edges)
    boundary_edges = {edge for edge, count in edge_counts.items() if count == 1}
    frontier_edges = {edge for ring in frontier_rings for edge in ring.edges}
    return boundary_edges & frontier_edges


def add_layered_quasi_cages(
    frame: Frame,
    ring_by_id: dict[str, Ring],
    edge_to_rings: dict[tuple[int, int], list[Ring]],
    allowed_sizes: set[int],
    ring_centers: dict[str, np.ndarray],
    base_size: int,
    seed_ids: list[str],
    frontier_ids: list[str],
    layer_sequences: tuple[str, ...],
    max_layers: int,
    max_rings_per_layer: int,
    max_layer_states_per_seed: int,
    max_layer_candidates: int,
    half_cages: list[CagePatch],
    quasi_cages: list[CagePatch],
    seen_quasi: set[frozenset[str]],
    seen_half: set[frozenset[str]],
    type_counts: dict[str, int],
) -> None:
    """Grow and report L2/L3 quasi-cage patches from one L0+L1 seed."""
    if max_layers < 2:
        return
    states: list[tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]] = [
        (tuple(sorted(set(seed_ids))), tuple(sorted(set(frontier_ids))), layer_sequences)
    ]
    for _layer_index in range(2, max_layers + 1):
        next_states: list[tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]] = []
        for patch_ids, current_frontier_ids, sequences in states:
            patch_rings = [ring_by_id[ring_id] for ring_id in patch_ids]
            frontier_rings = [ring_by_id[ring_id] for ring_id in current_frontier_ids]
            lower_rings = [ring_by_id[ring_id] for ring_id in set(patch_ids) - set(current_frontier_ids)]
            candidates = next_layer_candidates(
                patch_rings,
                frontier_rings,
                edge_to_rings,
                allowed_sizes,
                ring_centers,
                frame.box,
                lower_rings=lower_rings,
                max_candidates=max_layer_candidates,
            )
            remaining_budget = max_layer_states_per_seed - len(next_states)
            for layer_rings in connected_layer_subsets(candidates, max_rings_per_layer, remaining_budget):
                layer_ids = tuple(sorted(ring.object_id for ring in layer_rings))
                all_ids = tuple(sorted(set(patch_ids) | set(layer_ids)))
                if len(all_ids) == len(patch_ids):
                    continue
                key = frozenset(all_ids)
                if key in seen_quasi or key in seen_half:
                    continue
                layer_seq = layer_sequence(layer_rings)
                new_sequences = (*sequences, layer_seq)
                add_classified_patch(
                    frame,
                    ring_by_id,
                    base_size,
                    list(all_ids),
                    new_sequences,
                    half_cages,
                    quasi_cages,
                    seen_half,
                    seen_quasi,
                    type_counts,
                )
                next_states.append((all_ids, layer_ids, new_sequences))
                if len(next_states) >= max_layer_states_per_seed:
                    break
            if len(next_states) >= max_layer_states_per_seed:
                break
        states = unique_layer_states(next_states)
        if not states:
            return


def connected_layer_subsets(candidates: list[Ring], max_size: int, max_states: int) -> list[list[Ring]]:
    """Return bounded connected growth units for one dangling layer."""
    if max_size <= 0 or max_states <= 0:
        return []
    ring_by_id = {ring.object_id: ring for ring in candidates}
    adjacency = layer_adjacency(candidates)
    units: list[list[Ring]] = []
    seen: set[frozenset[str]] = set()

    def add_unit(ring_ids: set[str] | frozenset[str]) -> None:
        if len(units) >= max_states:
            return
        key = frozenset(ring_ids)
        if not key or key in seen:
            return
        seen.add(key)
        units.append([ring_by_id[ring_id] for ring_id in sorted(key)])

    for component in connected_components(candidates, adjacency):
        if len(units) >= max_states:
            break
        if len(component) <= max_size:
            add_unit(component)
            continue

        # Large frontiers are represented by every dangling ring plus local
        # connected neighborhoods, rather than every possible subset.
        for ring_id in sorted(component):
            add_unit({ring_id})
        for ring_id in sorted(component):
            add_unit(limited_neighborhood(ring_id, adjacency, max_size))
    return units


def connected_components(candidates: list[Ring], adjacency: dict[str, set[str]]) -> list[set[str]]:
    """Return connected components among layer candidates."""
    remaining = {ring.object_id for ring in candidates}
    components: list[set[str]] = []
    while remaining:
        start = min(remaining)
        stack = [start]
        component: set[str] = set()
        while stack:
            ring_id = stack.pop()
            if ring_id in component:
                continue
            component.add(ring_id)
            stack.extend(sorted(adjacency.get(ring_id, set()) & remaining, reverse=True))
        remaining -= component
        components.append(component)
    return components


def limited_neighborhood(start: str, adjacency: dict[str, set[str]], max_size: int) -> set[str]:
    """Grow a deterministic connected neighborhood around one candidate ring."""
    selected = {start}
    queue = [start]
    while queue and len(selected) < max_size:
        current = queue.pop(0)
        for neighbor in sorted(adjacency.get(current, set())):
            if neighbor in selected:
                continue
            selected.add(neighbor)
            queue.append(neighbor)
            if len(selected) >= max_size:
                break
    return selected


def layer_adjacency(rings: list[Ring]) -> dict[str, set[str]]:
    """Connect layer candidates that share a ring edge."""
    adjacency = {ring.object_id: set() for ring in rings}
    for idx, left in enumerate(rings):
        for right in rings[idx + 1 :]:
            if left.edges & right.edges:
                adjacency[left.object_id].add(right.object_id)
                adjacency[right.object_id].add(left.object_id)
    return adjacency


def unique_layer_states(states: list[tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]]) -> list[tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]]:
    """Keep one deterministic state for each full patch ring set."""
    seen: set[tuple[str, ...]] = set()
    unique: list[tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]] = []
    for state in sorted(states, key=lambda item: (len(item[0]), item[0], item[1])):
        patch_ids = state[0]
        if patch_ids in seen:
            continue
        seen.add(patch_ids)
        unique.append(state)
    return unique


def canonical_sequence(values: list[int]) -> str:
    """Canonicalize cyclic layer sizes while preserving ring-size isomers."""
    n = len(values)
    rotations = [values[i:] + values[:i] for i in range(n)]
    rev = list(reversed(values))
    rotations.extend(rev[i:] + rev[:i] for i in range(n))
    return "".join(str(value) for value in min(rotations))


def layer_sequence(rings: list[Ring]) -> str:
    """Return a deterministic sequence for a non-cyclic outer layer."""
    return "".join(str(ring.size) for ring in sorted(rings, key=lambda ring: (ring.size, ring.object_id)))


def quasi_label(base_size: int, layer_sequences: tuple[str, ...]) -> str:
    """Build a quasi-cage label from L1/L2/L3 layer sequences."""
    return f"qc_{base_size}r_" + "_".join(layer_label(sequence, True) for sequence in layer_sequences)


def add_classified_patch(
    frame: Frame,
    ring_by_id: dict[str, Ring],
    base_size: int,
    ring_ids: list[str],
    layer_sequences: tuple[str, ...],
    half_cages: list[CagePatch],
    quasi_cages: list[CagePatch],
    seen_half: set[frozenset[str]],
    seen_quasi: set[frozenset[str]],
    type_counts: dict[str, int],
) -> bool:
    """Classify one layered patch as a standard half-cage or quasi-cage."""
    half_type = half_cage_type(base_size, layer_sequences)
    if half_type is not None:
        layers = (f"{base_size}r", *[layer_label(sequence, False) for sequence in layer_sequences])
        patch = build_patch(frame, ring_by_id, "half_cage", half_type, ring_ids, layers)
        return add_patch(patch, half_cages, seen_half, type_counts)

    label = quasi_label(base_size, layer_sequences)
    layers = (f"{base_size}r", *[layer_label(sequence, True) for sequence in layer_sequences])
    patch = build_patch(frame, ring_by_id, "quasi_cage", label, ring_ids, layers)
    return add_patch(patch, quasi_cages, seen_quasi, type_counts)


def half_cage_type(base_size: int, layer_sequences: tuple[str, ...]) -> str | None:
    """Return the standard half-cage type for a layered patch, if any."""
    if base_size == 5 and layer_sequences == ("55555",):
        return "hc_5r_5⁵"
    if base_size == 6 and layer_sequences == ("555555",):
        return "hc_6r_5⁶"
    if base_size == 6 and layer_sequences == ("555555", "6"):
        return "hc_6r_5⁶_6¹"
    return None


def layer_label(sequence: str, include_isomer: bool) -> str:
    """Render one layer as composition plus optional subscript isomer."""
    counts = Counter(sequence)
    composition = "".join(f"{size}{superscript_number(counts[size])}" for size in sorted(counts, key=int))
    return f"{composition}{subscript_text(sequence)}" if include_isomer else composition


def superscript_number(value: int) -> str:
    """Return an integer as superscript Arabic numerals."""
    return str(value).translate(SUPERSCRIPT_DIGITS)


def subscript_text(value: str) -> str:
    """Return a digit string as subscript Arabic numerals."""
    return str(value).translate(SUBSCRIPT_DIGITS)


def build_patch(
    frame: Frame,
    ring_by_id: dict[str, Ring],
    kind: str,
    patch_type: str,
    ring_ids: list[str],
    layers: tuple[str, ...],
) -> CagePatch:
    """Build a CagePatch and compute its center from unwrapped oxygen nodes."""
    unique_ring_ids = tuple(sorted(set(ring_ids)))
    face_rings = [ring_by_id[ring_id] for ring_id in unique_ring_ids]
    waters = tuple(sorted({node for ring in face_rings for node in ring.nodes}))
    center = patch_center(frame, face_rings)
    return CagePatch(
        object_id=patch_type,
        patch_type=patch_type,
        kind=kind,
        rings=unique_ring_ids,
        waters=waters,
        center=center,
        layers=layers,
    )


def add_patch(
    patch: CagePatch,
    output: list[CagePatch],
    seen: set[frozenset[str]],
    type_counts: dict[str, int],
) -> bool:
    """Append a patch once and assign a stable object id."""
    key = frozenset(patch.rings)
    if key in seen:
        return False
    type_counts[patch.patch_type] += 1
    output.append(
        CagePatch(
            object_id=f"{patch.patch_type}_{type_counts[patch.patch_type]:05d}",
            patch_type=patch.patch_type,
            kind=patch.kind,
            rings=patch.rings,
            waters=patch.waters,
            center=patch.center,
            layers=patch.layers,
        )
    )
    seen.add(key)
    return True


def remove_subset_patches(patches: list[CagePatch]) -> list[CagePatch]:
    """Keep maximal patches within one class to avoid nested double counts."""
    ring_sets = [set(patch.rings) for patch in patches]
    maximal: list[CagePatch] = []
    for index, patch in enumerate(patches):
        patch_rings = ring_sets[index]
        if any(patch_rings < other_rings for other_index, other_rings in enumerate(ring_sets) if other_index != index):
            continue
        maximal.append(patch)
    return maximal


def patch_center(frame: Frame, rings: list[Ring]) -> np.ndarray:
    """Use the centroid of unique, locally unwrapped oxygen nodes."""
    nodes = sorted({node for ring in rings for node in ring.nodes})
    edges = sorted({edge for ring in rings for edge in ring.edges})
    unwrapped = unwrap_connected_nodes(frame, nodes, edges)
    return np.mean([unwrapped[node] for node in nodes], axis=0)
