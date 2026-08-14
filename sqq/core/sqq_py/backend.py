from __future__ import annotations

"""SQQ-Py frame-analysis backend."""

from typing import Any, Callable

from ...config import (
    normalize_order_parameters,
    output_enabled,
    q_degrees_from_order_parameters,
)
from ...models import Cage, CagePatch, Frame, FrameResult, Guest, HydrateOrderResult, Water
from ..cage import find_cages
from ..graph import build_water_graph
from ..phase import analyze_hydrate_clusters
from ..ice import classify_ice_waters
from ..order import compute_dhop_order, compute_mcg_order, compute_order_parameters
from ..half_quasi import find_cage_patches
from ..ring import find_rings
from ..ring_topology import build_ring_topology_index


def analyze_frame(
    frame: Frame,
    waters: list[Water],
    guests: list[Guest],
    config: dict[str, Any],
    *,
    pair_edges: list[tuple[int, int]] | None = None,
    cage_report_types: tuple[str, ...] | None,
    ring_report_sizes: tuple[int, ...],
    ring_sizes: list[int],
    quasi_base_sizes: list[int],
    quasi_side_sizes: list[int],
    stage_callback: Callable[[str], None] | None = None,
) -> FrameResult:
    """Analyze a selected frame with the SQQ-Py scientific kernel.

    ``config`` must already be normalized, and ``waters``/``guests`` must be
    selected by the compatibility dispatcher.  Keeping those responsibilities
    outside the backend gives SQQ-Py and SQQ-CPP the same dispatch boundary.
    """
    _report_stage(stage_callback, "building water graph")
    graph_config = config["graph"]
    graph_mode = str(
        graph_config.get("effective_bond_mode", graph_config["bond_mode"])
    )
    graph = build_water_graph(
        frame.atoms,
        waters,
        frame.box,
        bond_mode=graph_mode,
        oo_cutoff_nm=float(graph_config["oo_cutoff_nm"]),
        hbond_distance_nm=float(graph_config["hbond_distance_nm"]),
        hbond_angle_deg=float(graph_config["hbond_angle_deg"]),
        pair_edges=pair_edges,
    )
    _report_stage(stage_callback, "searching rings")
    rings = find_rings(
        graph.adjacency,
        sizes=ring_sizes,
        chordless=bool(config["ring"]["chordless"]),
        definition=str(config["ring"].get("definition", "chordless")),
    )
    scientific_validation = bool(config["cage"].get("scientific_validation", False))
    hydrate_cluster_enabled = bool(config.get("hydrate_cluster", {}).get("enabled", False))
    find_half = bool(config.get("half_cage", {}).get("enabled", False))
    find_quasi = bool(config["quasi_cage"].get("enabled", False))
    ring_topology = build_ring_topology_index(
        frame,
        rings,
        compute_face_quality=scientific_validation,
        compute_face_normals=hydrate_cluster_enabled,
        compute_ring_centers=find_half or find_quasi or hydrate_cluster_enabled,
        compute_adjacency=False,
    )
    warnings: list[str] = []
    if find_half or find_quasi:
        _report_stage(stage_callback, "searching half/quasi cage")
    half_cages, quasi_cages = find_cage_patches(
        frame,
        rings,
        find_half=find_half,
        find_quasi=find_quasi,
        base_sizes=quasi_base_sizes,
        side_sizes=quasi_side_sizes,
        max_combinations_per_base=int(config["quasi_cage"].get("max_combinations_per_base", 50000)),
        max_layers=int(config["quasi_cage"].get("max_layers", 1)),
        max_rings_per_layer=int(config["quasi_cage"].get("max_rings_per_layer", config["quasi_cage"].get("max_outer_layer_rings", 6))),
        max_layer_states_per_seed=int(config["quasi_cage"].get("max_layer_states_per_seed", 200)),
        max_candidates_per_edge=int(config["quasi_cage"].get("max_candidates_per_edge", 4)),
        max_layer_candidates=int(config["quasi_cage"].get("max_layer_candidates", 24)),
        topology_index=ring_topology,
        search_policy=str(config["quasi_cage"].get("search_policy", "bounded")),
        warnings=warnings,
    )
    cage_seed_patches = [*half_cages, *quasi_cages]
    _report_stage(stage_callback, "searching cage")
    cage_ring_sizes = [size for size in ring_sizes if size in {4, 5, 6}]
    all_cages = find_cages(
        frame,
        rings,
        cage_seed_patches,
        guests,
        enabled=bool(config["cage"].get("enabled", False)),
        ring_sizes=cage_ring_sizes,
        max_faces=int(config["cage"].get("max_faces", 20)),
        search_mode=str(config["cage"].get("search_mode", "grow")),
        seed_mode=str(config["cage"].get("seed_mode", "ring")),
        max_states_per_seed=int(config["cage"].get("max_states_per_seed", 0)),
        max_total_states=int(config["cage"].get("max_total_states", 0)),
        max_boundary_candidates=int(config["cage"].get("max_boundary_candidates", 8)),
        occupancy_radius_nm=float(config["cage"].get("occupancy_radius_nm", 0.5)),
        occupancy_mode=str(config["cage"].get("occupancy_mode", "polyhedron")),
        scientific_validation=scientific_validation,
        max_face_planarity_rms_nm=float(config["cage"].get("max_face_planarity_rms_nm", 0.06)),
        max_face_edge_cv=float(config["cage"].get("max_face_edge_cv", 0.35)),
        min_cage_volume_nm3=float(config["cage"].get("min_cage_volume_nm3", 1.0e-6)),
        topology_index=ring_topology,
    )
    cages = select_reported_cages(all_cages, cage_report_types)
    hydrate_cluster_detail = output_enabled(config, "cluster-detail")
    if hydrate_cluster_enabled:
        _report_stage(stage_callback, "classifying hydrate cluster")
        rings_by_id = ring_topology.ring_by_id
        ring_sizes_by_id = {ring_id: ring.size for ring_id, ring in rings_by_id.items()}
        hydrate_clusters, hydrate_motifs, hydrate_domains, isolated_cage_ids = analyze_hydrate_clusters(
            all_cages,
            min_cage=int(config.get("hydrate_cluster", {}).get("min_cage", 2)),
            ring_sizes=ring_sizes_by_id,
            frame=frame,
            rings_by_id=rings_by_id,
            face_geometries=ring_topology.face_geometries(),
        )
    else:
        hydrate_clusters, hydrate_motifs, hydrate_domains, isolated_cage_ids = [], [], [], ()
    _report_stage(stage_callback, "filtering free patches")
    quasi_cages = filter_free_patches(quasi_cages, all_cages)
    half_cages = filter_free_patches(half_cages, all_cages, higher_priority_patches=quasi_cages)
    focus_resids = {int(item) for item in config["order"].get("focus_waters", [])}
    _report_stage(stage_callback, "computing order parameters")
    order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters", ["f3", "f4"])
    )
    selected_order_parameters = set(order_parameters)
    q_degrees = q_degrees_from_order_parameters(order_parameters)
    if selected_order_parameters & {"f3", "f4"} or q_degrees:
        f3f4 = compute_order_parameters(
            frame,
            waters,
            graph,
            f3_enabled="f3" in selected_order_parameters,
            f4_enabled="f4" in selected_order_parameters,
            q_enabled=bool(q_degrees),
            q_neighbor_mode=str(config["order"].get("q_neighbor_mode", "graph")),
            q_cutoff_nm=float(config["order"].get("q_cutoff_nm", 0.35)),
            q_n_neighbor=config["order"].get("q_n_neighbor", None),
            q_degree=q_degrees,
            focus_resids=focus_resids,
        )
    else:
        f3f4 = None

    hydrate_parameters = selected_order_parameters & {
        "mcg1",
        "mcg3",
        "dhop35",
        "dhop30",
    }
    if hydrate_parameters:
        hydrate_order_config = {
            **config.get("hydrate_order", {}),
            "mcg1_enabled": "mcg1" in hydrate_parameters,
            "mcg3_enabled": "mcg3" in hydrate_parameters,
            "dhop35_enabled": "dhop35" in hydrate_parameters,
            "dhop30_enabled": "dhop30" in hydrate_parameters,
        }
        mcg1, mcg3 = compute_mcg_order(frame, waters, guests, hydrate_order_config)
        dhop35, dhop30 = compute_dhop_order(frame, waters, hydrate_order_config)
        hydrate_order = HydrateOrderResult(
            mcg1=mcg1,
            dhop35=dhop35,
            mcg3=mcg3,
            dhop30=dhop30,
        )
    else:
        hydrate_order = None
    _report_stage(stage_callback, "classifying ice")
    ice_classes = classify_ice_waters(
        graph,
        waters,
        rings,
        enabled=bool(config["ice"].get("enabled", False)),
        min_six_rings=int(config["ice"].get("min_six_rings", 2)),
        require_four_coord_neighbors=bool(config["ice"].get("require_four_coord_neighbors", True)),
    )
    if (find_half or find_quasi) and not half_cages and not quasi_cages:
        warnings.append("No enabled half_cage or quasi_cage was found with the current patch criteria.")
    return FrameResult(
        frame=frame,
        waters=waters,
        guests=guests,
        graph=graph,
        rings=rings,
        ring_report_sizes=tuple(ring_report_sizes),
        half_cages=half_cages,
        quasi_cages=quasi_cages,
        cages=cages,
        all_cages=all_cages,
        cage_report_types=cage_report_types,
        hydrate_cluster_enabled=hydrate_cluster_enabled,
        hydrate_cluster_detail=hydrate_cluster_detail,
        hydrate_clusters=hydrate_clusters,
        hydrate_motifs=hydrate_motifs,
        hydrate_domains=hydrate_domains,
        isolated_cage_ids=isolated_cage_ids,
        f3f4=f3f4,
        hydrate_order=hydrate_order,
        ice_like_waters=ice_classes.ice_like,
        ice_i_waters=ice_classes.ice_i,
        interfacial_ice_waters=ice_classes.interfacial,
        warnings=warnings,
    )


def select_reported_cages(
    cages: list[Cage],
    report_types: tuple[str, ...] | None,
) -> list[Cage]:
    """Filter detected cages for reports without changing topology filtering."""
    if report_types is None:
        return list(cages)
    allowed = set(report_types)
    return [cage for cage in cages if cage.cage_type in allowed]


def filter_free_patches(
    patches: list[CagePatch],
    cages: list[Cage],
    higher_priority_patches: list[CagePatch] | None = None,
) -> list[CagePatch]:
    """Remove consumed patches using ring-to-owner inverted indexes."""
    cage_ring_sets = [frozenset(cage.rings) for cage in cages]
    higher_priority_ring_sets = [frozenset(patch.rings) for patch in higher_priority_patches or []]
    cage_index = subset_owner_index(cage_ring_sets)
    higher_index = subset_owner_index(higher_priority_ring_sets)
    free_patches = []
    for patch in patches:
        patch_rings = frozenset(patch.rings)
        if is_subset_of_indexed_owner(patch_rings, cage_ring_sets, cage_index, strict=False):
            continue
        if is_subset_of_indexed_owner(patch_rings, higher_priority_ring_sets, higher_index, strict=True):
            continue
        free_patches.append(patch)
    return free_patches


def subset_owner_index(ring_sets: list[frozenset[str]]) -> dict[str, set[int]]:
    """Index candidate supersets by every ring they contain."""
    owners: dict[str, set[int]] = {}
    for index, ring_ids in enumerate(ring_sets):
        for ring_id in ring_ids:
            owners.setdefault(ring_id, set()).add(index)
    return owners


def is_subset_of_indexed_owner(
    ring_ids: frozenset[str],
    owners: list[frozenset[str]],
    index: dict[str, set[int]],
    *,
    strict: bool,
) -> bool:
    """Test subset ownership after narrowing candidates by the rarest ring."""
    if not ring_ids:
        return any((not strict) or bool(owner) for owner in owners)
    if any(ring_id not in index for ring_id in ring_ids):
        return False
    anchor = min(ring_ids, key=lambda ring_id: len(index[ring_id]))
    return any(
        ring_ids < owners[owner_index] if strict else ring_ids <= owners[owner_index]
        for owner_index in index[anchor]
    )


def _report_stage(callback: Callable[[str], None] | None, stage: str) -> None:
    if callback is not None:
        callback(stage)


analyze_frame_py = analyze_frame


__all__ = [
    "analyze_frame",
    "analyze_frame_py",
    "filter_free_patches",
    "is_subset_of_indexed_owner",
    "select_reported_cages",
    "subset_owner_index",
]
