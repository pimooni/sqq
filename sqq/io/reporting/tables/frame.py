"""Convert one scientific frame result into report records."""

from __future__ import annotations

from collections import Counter
from typing import Any

from ....models import FrameResult
from ....presentation.occupancy import (
    guest_composition_label,
    guest_lookup as build_guest_lookup,
)
from .formatting import (
    atom_resname_counts,
    cage_display_label,
    cluster_guest_composition,
    cluster_size_distribution,
    format_cage_type_counts,
    guest_resname_order,
    ordered_cage_types,
    patch_breakdown,
    patch_composition_label,
    patch_counts,
    present_cage_types,
    source_label,
)

QUASI_ISOMER_DETAIL_KEY = "_quasi_cage_isomer_detail"


def result_row(
    result: FrameResult,
    *,
    include_cluster_details: bool = True,
) -> dict[str, Any]:
    """Flatten a FrameResult into one summary-table row."""
    cage_counts: dict[str, int] = {}
    cage_detail_counts: dict[str, int] = {}
    cage_isomers: dict[str, dict[str, int]] = {}
    guests_by_id = build_guest_lookup(result.guests)
    guest_order = guest_resname_order(result)
    molecule_counts = atom_resname_counts(result)

    for cage in result.cages:
        cage_counts[cage.cage_type] = cage_counts.get(cage.cage_type, 0) + 1
        occupancy_key = f"cage_{cage.cage_type}_{'occupied' if cage.occupied else 'empty'}"
        cage_detail_counts[occupancy_key] = cage_detail_counts.get(occupancy_key, 0) + 1
        isomer = cage.isomer or "plain"
        type_isomers = cage_isomers.setdefault(cage.cage_type, {})
        type_isomers[isomer] = type_isomers.get(isomer, 0) + 1
        if cage.occupied:
            composition = guest_composition_label(cage, guests_by_id, guest_order)
            if composition:
                guest_key = f"cage_{cage.cage_type}_{composition}"
                cage_detail_counts[guest_key] = cage_detail_counts.get(guest_key, 0) + 1
            if len(cage.guest_ids) > 1:
                multi_key = f"cage_{cage.cage_type}_multi"
                cage_detail_counts[multi_key] = cage_detail_counts.get(multi_key, 0) + 1

    empty = sum(1 for cage in result.cages if not cage.occupied)
    occupied = sum(1 for cage in result.cages if cage.occupied)
    used_ring_ids = {ring_id for patch in [*result.half_cages, *result.quasi_cages] for ring_id in patch.rings}
    filtering_cages = result.all_cages or result.cages
    used_ring_ids.update(ring_id for cage in filtering_cages for ring_id in cage.rings)
    half_cage_counts = patch_counts(result.half_cages)
    quasi_cage_counts = patch_counts(result.quasi_cages)
    quasi_composition_counts: dict[str, int] = {}
    for patch_type, count in quasi_cage_counts.items():
        composition = patch_composition_label(patch_type)
        quasi_composition_counts[composition] = quasi_composition_counts.get(composition, 0) + count
    cluster_details = (
        hydrate_cluster_detail_records(result, guests_by_id, guest_order)
        if include_cluster_details and result.hydrate_cluster_detail
        else []
    )
    domain_details = (
        hydrate_domain_records(result, guests_by_id, guest_order)
        if include_cluster_details and result.hydrate_cluster_detail
        else []
    )
    largest_cluster = max(result.hydrate_clusters, key=lambda cluster: cluster.cage_count, default=None)
    cluster_type_counts = Counter(cluster.hydrate_type for cluster in result.hydrate_clusters)
    domain_type_counts = Counter(domain.hydrate_type for domain in result.hydrate_domains)
    classified_cage_count = sum(
        len(cluster.classified_cage_ids) for cluster in result.hydrate_clusters
    )
    boundary_cage_count = sum(
        len(cluster.boundary_cage_ids) for cluster in result.hydrate_clusters
    )
    ambiguous_cage_count = sum(
        len(cluster.ambiguous_cage_ids) for cluster in result.hydrate_clusters
    )
    unclassified_cage_count = sum(
        len(cluster.unclassified_cage_ids) for cluster in result.hydrate_clusters
    )

    def free_count(size: int) -> int:
        # Free rings are rings not consumed by open patches or cages.
        return sum(1 for ring in result.rings.get(size, []) if ring.object_id not in used_ring_ids)

    f3f4 = result.f3f4
    connection_counts = graph_connection_counts(result)
    row: dict[str, Any] = {
        "frame": result.frame.name,
        "time_ps": result.frame.time_ps,
        "source": source_label(result.frame.source),
        "status": "ok",
        "error": "",
        "n_atoms": len(result.frame.atoms),
        "n_waters": len(result.waters),
        "n_guests": len(result.guests),
        "bond_mode": result.graph.mode,
        "n_edges": len(result.graph.edges),
        **connection_counts,
        "ring4": len(result.rings.get(4, [])),
        "ring5": len(result.rings.get(5, [])),
        "ring6": len(result.rings.get(6, [])),
        "ring7": len(result.rings.get(7, [])),
        "free_ring4": free_count(4),
        "free_ring5": free_count(5),
        "free_ring6": free_count(6),
        "free_ring7": free_count(7),
        "half_cage_total": len(result.half_cages),
        "half_cage_breakdown": patch_breakdown(half_cage_counts),
        "quasi_cage_total": len(result.quasi_cages),
        "quasi_cage_breakdown": patch_breakdown(quasi_cage_counts),
        "cage_report_types": "all" if result.cage_report_types is None else ";".join(result.cage_report_types),
        "cage_512": cage_counts.get("512", 0),
        "cage_51262": cage_counts.get("51262", 0),
        "cage_51263": cage_counts.get("51263", 0),
        "cage_51264": cage_counts.get("51264", 0),
        "cage_51268": cage_counts.get("51268", 0),
        "cage_435663": cage_counts.get("435663", 0),
        "cage_total": len(result.cages),
        "cage_empty": empty,
        "cage_occupied": occupied,
        "hydrate_cluster_enabled": "on" if result.hydrate_cluster_enabled else "off",
        "hydrate_cluster_detail_enabled": "on" if result.hydrate_cluster_detail else "off",
        "hydrate_cluster_count": len(result.hydrate_clusters) if result.hydrate_cluster_enabled else 0,
        "sI_cluster_count": cluster_type_counts.get("sI", 0),
        "sII_cluster_count": cluster_type_counts.get("sII", 0),
        "sH_cluster_count": cluster_type_counts.get("sH", 0),
        "mixed_cluster_count": cluster_type_counts.get("mixed", 0),
        "unclassified_cluster_count": cluster_type_counts.get("unclassified", 0),
        "hydrate_domain_count": len(result.hydrate_domains),
        "sI_domain_count": domain_type_counts.get("sI", 0),
        "sII_domain_count": domain_type_counts.get("sII", 0),
        "sH_domain_count": domain_type_counts.get("sH", 0),
        "classified_cage_count": classified_cage_count,
        "boundary_cage_count": boundary_cage_count,
        "ambiguous_cage_count": ambiguous_cage_count,
        "unclassified_cage_count": unclassified_cage_count,
        "isolated_cage_count": len(result.isolated_cage_ids) if result.hydrate_cluster_enabled else 0,
        "largest_cluster_cage_count": 0 if largest_cluster is None else largest_cluster.cage_count,
        "largest_cluster_water_count": 0 if largest_cluster is None else largest_cluster.water_count,
        "cluster_size_distribution": cluster_size_distribution(result.hydrate_clusters),
        "MCG1_largest_cluster": None if result.hydrate_order is None else result.hydrate_order.mcg1.largest_cluster_size,
        "DHOP35_largest_cluster": None if result.hydrate_order is None else result.hydrate_order.dhop35.largest_cluster_size,
        "MCG3_largest_cluster": None if result.hydrate_order is None or result.hydrate_order.mcg3 is None else result.hydrate_order.mcg3.largest_cluster_size,
        "DHOP30_largest_cluster": None if result.hydrate_order is None or result.hydrate_order.dhop30 is None else result.hydrate_order.dhop30.largest_cluster_size,
        "MCG3_enabled": result.hydrate_order is not None and result.hydrate_order.mcg3 is not None,
        "DHOP30_enabled": result.hydrate_order is not None and result.hydrate_order.dhop30 is not None,
        "hydrate_cluster_detail": cluster_details,
        "hydrate_domain_detail": domain_details,
        "F3_mean": None if f3f4 is None else f3f4.f3_mean,
        "F4_mean": None if f3f4 is None else f3f4.f4_mean,
        "F3_count": None if f3f4 is None else f3f4.f3_valid,
        "F4_count": None if f3f4 is None else f3f4.f4_valid,
        "F3_valid_waters": None if f3f4 is None else f3f4.f3_valid,
        "F4_valid_waters": None if f3f4 is None else f3f4.f4_valid,
        "F3_focus_mean": None if f3f4 is None else f3f4.f3_focus_mean,
        "F4_focus_mean": None if f3f4 is None else f3f4.f4_focus_mean,
        "F3_focus_count": None if f3f4 is None else f3f4.f3_focus_valid,
        "F4_focus_count": None if f3f4 is None else f3f4.f4_focus_valid,
        "F3_focus_valid_waters": None if f3f4 is None else f3f4.f3_focus_valid,
        "F4_focus_valid_waters": None if f3f4 is None else f3f4.f4_focus_valid,
        "ice_like_waters": len(result.ice_like_waters),
        "ice_i_waters": len(result.ice_i_waters),
        "interfacial_ice_waters": len(result.interfacial_ice_waters),
        # Keep exact quasi isomers out of the wide summary table.
        QUASI_ISOMER_DETAIL_KEY: tuple(sorted(quasi_cage_counts.items())),
    }
    if f3f4 is not None:
        for degree in f3f4.q_degree:
            prefix = f"q{degree}"
            row[f"{prefix}_mean"] = f3f4.q_means.get(degree)
            row[f"{prefix}_count"] = f3f4.q_valid_counts.get(degree, 0)
            row[f"{prefix}_valid_waters"] = f3f4.q_valid_counts.get(degree, 0)
            row[f"{prefix}_focus_mean"] = f3f4.q_focus_means.get(degree)
            row[f"{prefix}_focus_count"] = f3f4.q_focus_valid_counts.get(degree, 0)
            row[f"{prefix}_focus_valid_waters"] = f3f4.q_focus_valid_counts.get(degree, 0)
    else:
        for degree in (6, 12):
            prefix = f"q{degree}"
            row[f"{prefix}_mean"] = None
            row[f"{prefix}_count"] = None
            row[f"{prefix}_valid_waters"] = None
            row[f"{prefix}_focus_mean"] = None
            row[f"{prefix}_focus_count"] = None
            row[f"{prefix}_focus_valid_waters"] = None
    for resname, count in molecule_counts.items():
        row[f"mol_{resname}"] = count
    row["mol_TOTAL"] = len(result.frame.atoms)
    row["guest_order"] = ";".join(guest_order)
    for patch_type, count in half_cage_counts.items():
        row[f"half_cage_{patch_type}"] = count
    for composition, count in quasi_composition_counts.items():
        row[f"quasi_cage_{composition}"] = count

    cage_types = ordered_cage_types(cage_counts)
    for cage_type in cage_types:
        prefix = f"cage_{cage_type}"
        row[prefix] = cage_counts.get(cage_type, row.get(prefix, 0))
        row[f"{prefix}_empty"] = cage_detail_counts.get(f"{prefix}_empty", 0)
        row[f"{prefix}_occupied"] = cage_detail_counts.get(f"{prefix}_occupied", 0)
        row[f"{prefix}_multi"] = cage_detail_counts.get(f"{prefix}_multi", 0)
        if cage_type in cage_isomers:
            parts = [f"{key}:{cage_isomers[cage_type][key]}" for key in sorted(cage_isomers[cage_type])]
            row[f"{prefix}_isomers"] = "; ".join(parts)
            for isomer, count in cage_isomers[cage_type].items():
                row[f"{prefix}_isomer_{isomer}"] = count

    for key in sorted(cage_detail_counts):
        if key not in row:
            row[key] = cage_detail_counts[key]
    return row


def graph_connection_counts(result: FrameResult) -> dict[str, Any]:
    """Return graph counts and a diagnostic-only coordination distribution."""
    mode = result.graph.mode
    edge_count = len(result.graph.edges)
    degrees = [len(result.graph.adjacency.get(water.oxygen, set())) for water in result.waters]
    water_count = len(degrees)
    bins = {degree: sum(value == degree for value in degrees) for degree in range(5)}
    over_four = sum(value > 4 for value in degrees)

    def fraction(count: int) -> float:
        return 0.0 if water_count == 0 else count / water_count

    values: dict[str, Any] = {
        "connection_mode": mode,
        "connection_count": edge_count,
        "hbond_count": edge_count if mode == "hbond" else None,
        "oo_connection_count": edge_count if mode == "oo" else None,
        "pair_connection_count": edge_count if mode == "pairs" else None,
        "mean_coordination": 0.0 if water_count == 0 else sum(degrees) / water_count,
        "coordination_gt4": over_four,
        "coordination_gt4_fraction": fraction(over_four),
        "degree_le2_fraction": fraction(sum(bins[degree] for degree in range(3))),
        "degree4_fraction": fraction(bins[4]),
        "over4_fraction": fraction(over_four),
    }
    for degree in range(5):
        values[f"coordination_{degree}"] = bins[degree]
        values[f"coordination_{degree}_fraction"] = fraction(bins[degree])
    return values


def hydrate_cluster_detail_records(
    result: FrameResult,
    guests_by_id: dict[str, Any],
    guest_order: list[str],
) -> list[dict[str, Any]]:
    """Return one plotting-friendly detail record per hydrate cluster."""
    if not result.hydrate_cluster_enabled:
        return []
    cage_by_id = {cage.object_id: cage for cage in (result.all_cages or result.cages)}
    records: list[dict[str, Any]] = []
    for cluster in result.hydrate_clusters:
        cluster_cages = [
            cage_by_id[cage_id]
            for cage_id in cluster.cage_ids
            if cage_id in cage_by_id
        ]
        boundary_cages = [
            cage_by_id[cage_id]
            for cage_id in cluster.boundary_cage_ids
            if cage_id in cage_by_id
        ]
        type_counts = Counter(cage.cage_type for cage in cluster_cages)
        boundary_counts = Counter(cage.cage_type for cage in boundary_cages)
        ordered_types = present_cage_types(type_counts)
        classified_fraction = (
            0.0
            if not cluster.cage_count
            else len(cluster.classified_cage_ids) / cluster.cage_count
        )
        records.append(
            {
                "cluster_id": cluster.object_id,
                "hydrate_type": cluster.hydrate_type,
                "cage_count": cluster.cage_count,
                "water_count": cluster.water_count,
                "guest_count": cluster.guest_count,
                "empty_cage_count": sum(
                    1 for cage in cluster_cages if not cage.occupied
                ),
                "occupied_cage_count": sum(
                    1 for cage in cluster_cages if cage.occupied
                ),
                "classified_cage_count": len(cluster.classified_cage_ids),
                "boundary_cage_count": len(cluster.boundary_cage_ids),
                "ambiguous_cage_count": len(cluster.ambiguous_cage_ids),
                "unclassified_cage_count": len(cluster.unclassified_cage_ids),
                "classified_cage_fraction": classified_fraction,
                "domain_count": cluster.domain_count,
                "cage_type_counts": {
                    cage_type: type_counts[cage_type] for cage_type in ordered_types
                },
                "cage_composition": ";".join(
                    f"{cage_display_label(cage_type)}:{type_counts[cage_type]}"
                    for cage_type in ordered_types
                ),
                "boundary_composition": format_cage_type_counts(boundary_counts),
                "guest_composition": cluster_guest_composition(
                    cluster.guest_ids,
                    guests_by_id,
                    guest_order,
                ),
                "domain_ids": ";".join(cluster.domain_ids),
                "cage_ids": ";".join(cluster.cage_ids),
                "classified_cage_ids": ";".join(cluster.classified_cage_ids),
                "boundary_cage_ids": ";".join(cluster.boundary_cage_ids),
                "ambiguous_cage_ids": ";".join(cluster.ambiguous_cage_ids),
                "unclassified_cage_ids": ";".join(cluster.unclassified_cage_ids),
                "shared_face_count": len(cluster.shared_faces),
            }
        )
    return records


def hydrate_domain_records(
    result: FrameResult,
    guests_by_id: dict[str, Any],
    guest_order: list[str],
) -> list[dict[str, Any]]:
    """Return one plotting-friendly record per hydrate domain."""
    if not result.hydrate_cluster_enabled:
        return []
    cage_by_id = {cage.object_id: cage for cage in (result.all_cages or result.cages)}
    records: list[dict[str, Any]] = []
    for domain in result.hydrate_domains:
        domain_cages = [cage_by_id[cage_id] for cage_id in domain.cage_ids if cage_id in cage_by_id]
        type_counts = Counter(cage.cage_type for cage in domain_cages)
        ordered_types = present_cage_types(type_counts)
        seed_cage_ids = set(domain.seed_cage_ids)
        records.append(
            {
                "domain_id": domain.object_id,
                "cluster_id": domain.cluster_id,
                "hydrate_type": domain.hydrate_type,
                "status": domain.status,
                "cage_count": domain.cage_count,
                "seed_count": domain.seed_count,
                "seed_cage_count": len(seed_cage_ids),
                "expanded_cage_count": domain.cage_count - len(seed_cage_ids),
                "classified_fraction": domain.classified_fraction,
                "water_count": domain.water_count,
                "guest_count": domain.guest_count,
                "external_boundary_contact_count": len(domain.boundary_cage_ids),
                "cage_composition": ";".join(f"{cage_display_label(cage_type)}:{type_counts[cage_type]}" for cage_type in ordered_types),
                "guest_composition": cluster_guest_composition(domain.guest_ids, guests_by_id, guest_order),
                "cage_ids": ";".join(domain.cage_ids),
                "seed_cage_ids": ";".join(domain.seed_cage_ids),
                "external_boundary_contact_ids": ";".join(domain.boundary_cage_ids),
            }
        )
    return records


def hydrate_motif_records(result: FrameResult) -> list[dict[str, Any]]:
    """Return one record per overlapping local topology motif."""
    if not result.hydrate_cluster_enabled:
        return []
    cage_by_id = {cage.object_id: cage for cage in (result.all_cages or result.cages)}
    records: list[dict[str, Any]] = []
    for motif in result.hydrate_motifs:
        member_cages = [cage_by_id[cage_id] for cage_id in motif.cage_ids if cage_id in cage_by_id]
        type_counts = Counter(cage.cage_type for cage in member_cages)
        ordered_types = present_cage_types(type_counts)
        core_cages = [cage_by_id[cage_id] for cage_id in motif.anchor_cage_ids if cage_id in cage_by_id]
        core_type_counts = Counter(cage.cage_type for cage in core_cages)
        ordered_core_types = present_cage_types(core_type_counts)
        anchor_types = [cage_display_label(cage_by_id[cage_id].cage_type) for cage_id in motif.anchor_cage_ids if cage_id in cage_by_id]
        records.append(
            {
                "motif_id": motif.object_id,
                "cluster_id": motif.cluster_id,
                "domain_id": motif.domain_id,
                "hydrate_type": motif.hydrate_type,
                "status": motif.status,
                "completeness": motif.completeness,
                "consistency": motif.consistency,
                "confidence": motif.confidence,
                "anchor_cage_types": ";".join(anchor_types),
                "anchor_cage_ids": ";".join(motif.anchor_cage_ids),
                "member_cage_count": motif.cage_count,
                "cage_count": motif.cage_count,
                "support_cage_count": motif.cage_count - len(motif.anchor_cage_ids),
                "cage_composition": ";".join(f"{cage_display_label(cage_type)}:{type_counts[cage_type]}" for cage_type in ordered_types),
                "core_cage_count": len(motif.anchor_cage_ids),
                "core_cage_composition": ";".join(f"{cage_display_label(cage_type)}:{core_type_counts[cage_type]}" for cage_type in ordered_core_types),
                "core_cage_ids": ";".join(motif.anchor_cage_ids),
                "motif_cage_ids": ";".join(motif.cage_ids),
                "member_cage_ids": ";".join(motif.cage_ids),
                "shared_face_count": len(motif.shared_face_ids),
                "internal_shared_face_count": len(motif.shared_face_ids),
                "shared_face_ids": ";".join(motif.shared_face_ids),
                "internal_shared_face_ids": ";".join(motif.shared_face_ids),
                "classification_method": motif.classification_method,
            }
        )
    return records

__all__ = (
    "result_row",
    "graph_connection_counts",
    "hydrate_cluster_detail_records",
    "hydrate_domain_records",
    "hydrate_motif_records",
)
