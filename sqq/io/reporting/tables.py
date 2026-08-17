"""Convert scientific frame results into shared report tables."""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import re
from typing import Any
import unicodedata

import pandas as pd

from ... import __version__
from ...banner import SQQ_AUTHOR, SQQ_TITLE
from ...citation import build_citation_recommendation, completed_citation_evidence
from ...config import (
    DEFAULT_MODE,
    is_cpp_mode,
    normalize_order_parameters,
    order_parameter_display,
    output_enabled,
    output_type_display,
    q_degrees_from_order_parameters,
)
from ...core.cage import KNOWN_CAGE_TYPES, parse_cage_face_label
from ...display import graph_mode_display
from ...models import FrameResult
from ..occupancy import (
    guest_composition_label,
    guest_lookup as build_guest_lookup,
    guest_resname_order as guest_resname_order_from_guests,
)
from .models import ReportTable

SUBSCRIPT_DIGIT_DELETE = dict.fromkeys(range(0x2080, 0x208A))
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


def cluster_guest_composition(guest_ids: tuple[str, ...], guests_by_id: dict[str, Any], guest_order: list[str]) -> str:
    """Summarize all guest residue names inside one hydrate cluster."""
    names = [guests_by_id[item].resname for item in guest_ids if item in guests_by_id]
    if not names:
        return ""
    counts = Counter(names)
    order_index = {name: index for index, name in enumerate(guest_order)}
    ordered_names = sorted(counts, key=lambda name: (order_index.get(name, 10_000), name))
    return "+".join(name if counts[name] == 1 else f"{name}x{counts[name]}" for name in ordered_names)


def cluster_size_distribution(clusters) -> str:
    """Summarize cluster sizes as cage_count:number_of_clusters."""
    if not clusters:
        return ""
    counts = Counter(cluster.cage_count for cluster in clusters)
    return ";".join(f"{size}:{counts[size]}" for size in sorted(counts))


def patch_counts(patches) -> dict[str, int]:
    """Count open cage patches by patch_type."""
    counts: dict[str, int] = {}
    for patch in patches:
        counts[patch.patch_type] = counts.get(patch.patch_type, 0) + 1
    return counts


def patch_breakdown(counts: dict[str, int]) -> str:
    """Render a compact patch count list for broad summary rows."""
    return "; ".join(f"{key}:{counts[key]}" for key in sorted(counts))


def atom_resname_counts(result: FrameResult) -> dict[str, int]:
    """Count atoms by residue name while preserving source-file order."""
    counts: dict[str, int] = {}
    for atom in result.frame.atoms:
        counts[atom.resname] = counts.get(atom.resname, 0) + 1
    return counts


def superscript_number(value: int) -> str:
    """Render small integer counts with Unicode superscript digits."""
    superscripts = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    return str(value).translate(superscripts)


def ordered_cage_types(types) -> list[str]:
    """Return canonical hydrate cages first, followed by other cage labels."""
    values = set(types.keys() if isinstance(types, dict) else types)
    ordered = list(KNOWN_CAGE_TYPES)
    extras = sorted((cage_type for cage_type in values if cage_type not in KNOWN_CAGE_TYPES), key=cage_sort_key)
    return ordered + extras


def present_cage_types(types) -> list[str]:
    """Return canonical cage ordering restricted to values that are present."""
    values = set(types.keys() if isinstance(types, dict) else types)
    return [cage_type for cage_type in ordered_cage_types(values) if cage_type in values]


def cage_sort_key(cage_type: str) -> tuple[int, tuple[tuple[int, int], ...], str]:
    """Sort other cages by face count and then by label."""
    counts = parse_cage_face_label(cage_type)
    if counts is None:
        return (999, (), cage_type)
    return (sum(counts.values()), tuple(sorted(counts.items())), cage_type)


def cage_display_label(cage_type: str) -> str:
    """Render compact cage labels for human-facing markdown tables."""
    known = {
        "512": f"5{superscript_number(12)}",
        "51262": f"5{superscript_number(12)}6{superscript_number(2)}",
        "51263": f"5{superscript_number(12)}6{superscript_number(3)}",
        "51264": f"5{superscript_number(12)}6{superscript_number(4)}",
        "51268": f"5{superscript_number(12)}6{superscript_number(8)}",
        "435663": f"4{superscript_number(3)}5{superscript_number(6)}6{superscript_number(3)}",
    }
    if cage_type in known:
        return known[cage_type]
    counts = parse_cage_face_label(cage_type)
    if counts is None:
        return cage_type
    return "".join(f"{size}{superscript_number(count)}" for size, count in sorted(counts.items()) if count > 0)


def format_cage_type_counts(counts: Counter) -> str:
    """Format nonzero cage-type counts in canonical display order."""
    return ", ".join(
        f"{cage_display_label(cage_type)}={int(counts[cage_type])}"
        for cage_type in present_cage_types(counts)
        if int(counts[cage_type]) > 0
    )


def guest_composition_sort_key(label: str, guest_order: list[str]) -> tuple[Any, ...]:
    """Sort exact guest compositions by occupancy size and source guest order."""
    order_index = {name: index for index, name in enumerate(guest_order)}
    components: list[tuple[int, str, int]] = []
    total_guests = 0
    for part in label.split("+"):
        name, separator, count_text = part.rpartition("x")
        if separator and count_text.isdigit():
            count = int(count_text)
        else:
            name = part
            count = 1
        total_guests += count
        components.append((order_index.get(name, 10_000), name, count))
    return total_guests, len(components), tuple(components), label


def guest_resname_order(result: FrameResult) -> list[str]:
    """Return guest residue names by their first atom position in the frame."""
    return guest_resname_order_from_guests(result.guests)


def patch_display_label(patch_type: str) -> str:
    """Remove internal HC/QC prefixes from human-facing patch labels."""
    for prefix in ("hc_", "qc_"):
        if patch_type.startswith(prefix):
            return patch_type.removeprefix(prefix)
    return patch_type


def source_label(source: Path | None) -> str:
    """Return an absolute source path for human-facing reports."""
    if source is None:
        return ""
    return str(Path(source).resolve())


def summary_output_tables(
    data: pd.DataFrame,
    config: dict[str, Any],
    run_info: dict[str, Any],
    detail_index: pd.DataFrame,
) -> tuple[ReportTable, ...]:
    """Build the shared table list for summary XLSX and CSV outputs."""
    tables: list[ReportTable] = [
        ReportTable(
            "summary",
            summary_dashboard_table(data, run_info, config),
            include_header=False,
        ),
    ]
    tables.extend(
        ReportTable(sheet_name, table)
        for sheet_name, table in summary_sheet_tables(data, config).items()
    )
    if not detail_index.empty:
        tables.append(ReportTable("detail_index", detail_index))
    return tuple(tables)


def summary_dashboard_table(data: pd.DataFrame, run_info: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    """Build a compact human-facing dashboard for the summary sheet."""
    cpp_mode = is_cpp_mode(config.get("mode", DEFAULT_MODE))
    title = SQQ_TITLE
    author = SQQ_AUTHOR
    matched_files = run_info.get("matched_files", "")
    try:
        matched_count = int(matched_files)
    except (TypeError, ValueError):
        matched_count = 0
    graph_mode_value = run_info.get("graph_mode_display") or graph_mode_display(
        run_info.get("graph_mode", config.get("graph", {}).get("bond_mode", "auto")),
        data.get("connection_mode", []),
    )
    selected_order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters", ["f3", "f4"])
    )
    selected_order_set = set(selected_order_parameters)
    q_degrees = q_degrees_from_order_parameters(selected_order_parameters)

    rows: list[list[Any]] = [
        [title, ""],
        [author, ""],
        ["", ""],
        ["Basic Information", ""],
        ["Date", run_info.get("date", "")],
        ["Start time", run_info.get("start_time", "")],
        ["Finish time", run_info.get("finish_time", "")],
        ["Time zone", run_info.get("time_zone", "")],
        ["Duration (s)", run_info.get("elapsed_seconds", "")],
        ["Working directory", run_info.get("working_dir", "")],
        ["Input", run_info.get("input", "")],
        ["Matched files", matched_files],
        ["Input format", run_info.get("input_format", "")],
    ]
    if matched_count > 1:
        rows.extend([
            ["First file", run_info.get("first_file", "")],
            ["Last file", run_info.get("last_file", "")],
        ])
    else:
        rows.append(["Current file", run_info.get("first_file", "")])
    rows.extend([
        ["Output directory", run_info.get("output_dir", "")],
        ["summary.xlsx", run_info.get("summary_xlsx", "")],
        ["summary_csv", run_info.get("summary_csv", "")],
        ["summary_detail_csv", run_info.get("summary_detail_csv", "")],
        ["sqq_config_resolved.yaml", run_info.get("config_output", "")],
        ["", ""],
        ["Configuration", ""],
        ["SQQ version", run_info.get("sqq_version", __version__)],
        ["SQQ engine", "sqq-cpp" if cpp_mode else "sqq-py"],
        ["Engine selector", run_info.get("engine_selector", config.get("mode", DEFAULT_MODE))],
        ["Profile", run_info.get("profile", config.get("run", {}).get("profile", ""))],
        ["Config file", run_info.get("config_file", "<built-in defaults>")],
        ["Topology", run_info.get("topology", "<none>")],
        ["Graph mode", graph_mode_value],
        ["Graph mode reason", run_info.get("graph_mode_reason", "")],
        ["Search sizes", excel_scalar(config.get("ring", {}).get("sizes", ""))],
    ])
    for adjustment in run_info.get("resolution_adjustments", ()):
        if not isinstance(adjustment, dict):
            continue
        parameter = adjustment.get("parameter", "parameter")
        effective = adjustment.get("effective")
        reason = adjustment.get("reason", "automatic adjustment")
        value = str(parameter)
        if effective is not None:
            value += f" -> {excel_scalar(effective)}"
        value += f" [{reason}]"
        rows.append(["Adjustment", value])
    if not cpp_mode:
        rows.append(["Ring report sizes", excel_scalar(configured_ring_report_sizes(config))])
    rows.append(["Ring definition", config.get("ring", {}).get("definition", "chordless")])
    if run_info.get("sampling_interval"):
        native = run_info.get("native_frame_interval_ps")
        delta = run_info.get("delta_time_ps")
        rows.extend([
            ["Native frame interval", "unknown" if native is None else f"{float(native):g} ps"],
            ["Delta time", "all" if delta is None else f"{float(delta):g} ps"],
            ["Raw frame step", run_info.get("raw_frame_step", 1)],
            ["Selected frames", f"{run_info.get('selected_frames', 0)} / {run_info.get('source_frames_total', 0)}"],
        ])
    if str(run_info.get("input_format", "")).startswith("lammps-"):
        rows.extend([
            ["LAMMPS units", run_info.get("lammps_units", "")],
            ["LAMMPS timestep", run_info.get("lammps_timestep", "")],
            ["LAMMPS atom style", run_info.get("lammps_atom_style", "")],
            ["LAMMPS type map", run_info.get("lammps_type_map_source", "")],
        ])
    if not cpp_mode:
        rows.extend([
            ["Find half", "on" if config.get("half_cage", {}).get("enabled", False) else "off"],
            ["Find quasi", "on" if config.get("quasi_cage", {}).get("enabled", False) else "off"],
        ])
        if config.get("quasi_cage", {}).get("enabled", False):
            rows.extend([
                ["Quasi-cage sizes", f"{excel_scalar(config.get('quasi_cage', {}).get('base_sizes', 'auto'))} / {excel_scalar(config.get('quasi_cage', {}).get('side_sizes', 'auto'))}"],
                ["Quasi max layer", config.get("quasi_cage", {}).get("max_layers", "")],
                ["Quasi search policy", config.get("quasi_cage", {}).get("search_policy", "bounded")],
            ])
    rows.extend([
        ["Cage report types", dashboard_cage_targets(config)],
        ["Maximum cage face", config.get("cage", {}).get("max_faces", 20)],
    ])
    if not cpp_mode:
        rows.extend([
            ["Find cluster", "on" if config.get("hydrate_cluster", {}).get("enabled", False) else "off"],
            ["Cluster min cage", config.get("hydrate_cluster", {}).get("min_cage", 2)],
        ])
    rows.append(["Order parameters", order_parameter_display(selected_order_parameters)])
    if not cpp_mode and selected_order_set & {"mcg1", "mcg3"}:
        rows.append([
            "MCG guest / water cutoff (nm)",
            f"{config.get('hydrate_order', {}).get('mcg_guest_cutoff_nm', 0.90)} / "
            f"{config.get('hydrate_order', {}).get('mcg_water_cutoff_nm', 0.60)}",
        ])
    if not cpp_mode and selected_order_set & {"dhop35", "dhop30"}:
        rows.append([
            "DHOP O-O cutoff (nm)",
            config.get("hydrate_order", {}).get("dhop_neighbor_cutoff_nm", 0.35),
        ])
    if not cpp_mode and q_degrees:
        rows.extend([
            ["Q_l degree", excel_scalar(q_degrees)],
            ["Q_l neighbor mode", config.get("order", {}).get("q_neighbor_mode", "graph")],
            ["Q_l cutoff (nm)", config.get("order", {}).get("q_cutoff_nm", 0.35)],
            ["Q_l n neighbor", config.get("order", {}).get("q_n_neighbor", "NULL")],
        ])
    rows.extend([
        [
            "Output types",
            run_info.get("output_types")
            or output_type_display(
                config.get("output", {}).get("types"),
                cpp_mode=cpp_mode,
            ),
        ],
        ["Output layout", run_info.get("output_layout", "")],
        ["Worker policy", run_info.get("worker_policy", "")],
        ["Parallel backend", run_info.get("parallel_backend", "serial")],
        ["Math threads per worker", run_info.get("math_threads", 1)],
        ["Workers", run_info.get("workers", "")],
        ["", ""],
        ["Analysis Results (min / mean / max)", ""],
        ["Frames total / ok / failed", f"{len(data)} / {frames_ok_count(data)} / {frames_failed_count(data)}"],
        ["Water molecules", min_mean_max_column(data, "n_waters")],
        ["Guest molecules", min_mean_max_column(data, "n_guests")],
        ["Connections", min_mean_max_column(data, "connection_count")],
    ])
    if not cpp_mode:
        for size in configured_ring_report_sizes(config):
            rows.append([f"Ring{size}", min_mean_max_column(data, f"ring{size}")])
        rows.extend([
            ["Half cage", min_mean_max_column(data, "half_cage_total")],
            ["Quasi cage", min_mean_max_column(data, "quasi_cage_total")],
        ])
    rows.append(["Cage total", min_mean_max_column(data, "cage_total")])
    if cpp_mode and not has_selected_guests(data):
        rows.append(["Cage occupancy", "not evaluated"])
    else:
        rows.extend([
            ["Empty cage", min_mean_max_column(data, "cage_empty")],
            ["Occupied cage", min_mean_max_column(data, "cage_occupied")],
        ])
    if not cpp_mode and config.get("hydrate_cluster", {}).get("enabled", False):
        rows.extend([
            ["Hydrate cluster", min_mean_max_column(data, "hydrate_cluster_count")],
            ["Isolated cage", min_mean_max_column(data, "isolated_cage_count")],
        ])
    if not cpp_mode:
        rows.append(["Ice-like waters", min_mean_max_column(data, "ice_like_waters")])

    completed_outputs = run_info.get("completed_outputs", run_info.get("output_types", ()))
    citation_statistics: dict[str, Any] = {
        **completed_citation_evidence(
            config,
            successful_frames=frames_ok_count(data),
            completed_outputs=completed_outputs,
            track=bool(run_info.get("track", False)),
        ),
        "failed_frames": frames_failed_count(data),
        "status": run_info.get("status", "completed"),
        "guest_molecules": sum_numeric_column(data, "n_guests"),
    }
    for key in ("executed_features", "executed_order_parameters", "occupancy_evaluated"):
        if key in run_info:
            citation_statistics[key] = run_info[key]
    citation = build_citation_recommendation(run_info, config, citation_statistics)
    rows.extend(
        [
            ["", ""],
            ["Citation Recommendation", ""],
            ["Recommended text", citation.sentence],
            ["Publication", citation.publication.removeprefix("Publication: ")],
            ["GitHub", citation.github.removeprefix("GitHub     : ")],
        ]
    )
    return pd.DataFrame(rows)


def frames_ok_count(data: pd.DataFrame) -> int:
    """Count successfully analyzed frames."""
    return int((data.get("status") == "ok").sum()) if "status" in data else len(data)


def frames_failed_count(data: pd.DataFrame) -> int:
    """Count failed frames."""
    return int((data.get("status") == "failed").sum()) if "status" in data else 0


def has_selected_guests(data: pd.DataFrame) -> bool:
    """Return whether any analyzed frame contains a selected guest molecule."""
    if "n_guests" not in data.columns:
        return False
    return bool((pd.to_numeric(data["n_guests"], errors="coerce").fillna(0) > 0).any())


def first_data_value(data: pd.DataFrame, column: str, fallback: Any = "") -> Any:
    """Return the first non-empty value in a summary column."""
    if column not in data:
        return fallback
    for value in data[column]:
        if pd.notna(value) and value != "":
            return value
    return fallback


def sum_numeric_column(data: pd.DataFrame, column: str) -> int:
    """Sum a numeric summary column while ignoring blanks."""
    if column not in data:
        return 0
    return int(pd.to_numeric(data[column], errors="coerce").fillna(0).sum())


def min_mean_max_column(data: pd.DataFrame, column: str) -> str:
    """Render per-frame min / mean / max statistics for a numeric summary column."""
    if column not in data:
        return "0 / 0.0 / 0"
    values = pd.to_numeric(data[column], errors="coerce").dropna()
    if values.empty:
        return "0 / 0.0 / 0"
    return " / ".join([
        format_stat_value(values.min()),
        format_stat_value(values.mean(), force_decimal=True),
        format_stat_value(values.max()),
    ])


def format_stat_value(value: Any, *, force_decimal: bool = False) -> str:
    """Format dashboard min/mean/max values compactly but readably."""
    numeric = float(value)
    if force_decimal and numeric.is_integer():
        return f"{numeric:.1f}"
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.3f}".rstrip("0").rstrip(".")


def molecule_totals(data: pd.DataFrame) -> dict[str, int]:
    """Collect total atom counts by residue name for the dashboard."""
    totals: dict[str, int] = {}
    for column in data.columns:
        if not column.startswith("mol_") or column == "mol_TOTAL":
            continue
        totals[column.removeprefix("mol_")] = sum_numeric_column(data, column)
    if "mol_TOTAL" in data:
        totals["TOTAL"] = sum_numeric_column(data, "mol_TOTAL")
    return totals


def configured_ring_report_sizes(config: dict[str, Any]) -> list[int]:
    """Return normalized ring report sizes for dashboards and data sheets."""
    search_sizes = config.get("ring", {}).get("sizes", [5, 6])
    value = config.get("ring", {}).get("report_sizes", "auto")
    if value in (None, "", "auto"):
        value = search_sizes
    if isinstance(value, str):
        return sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    return sorted({int(item) for item in value})


def dashboard_cage_targets(config: dict[str, Any]) -> str:
    """Render exact cage report types with human-facing superscripts."""
    targets = config.get("cage", {}).get("report_types", [])
    if isinstance(targets, str):
        if targets.strip().lower() in {"auto", "all"}:
            return "all detected cages (follows --size)"
        raw_targets = [item.strip() for item in targets.split(",") if item.strip()]
    else:
        raw_targets = [str(item) for item in targets or []]
    return ", ".join(cage_display_label(target) for target in raw_targets)


def stable_extra_columns(rows: list[dict[str, Any]], base_columns: list[str]) -> list[str]:
    """Collect non-core columns in first-seen order."""
    seen = set(base_columns)
    extras: list[str] = []
    for row in rows:
        for key in row:
            if str(key).startswith("_"):
                continue
            if key in seen:
                continue
            seen.add(key)
            extras.append(key)
    return extras


def summary_markdown(data: pd.DataFrame) -> str:
    """Render the global summary as readable grouped markdown tables."""
    lines = ["# SQQ summary", ""]
    for title, table in summary_markdown_tables(data):
        if table.empty:
            continue
        lines.extend(["", f"## {title}", "", markdown_table(table).rstrip()])
    return "\n".join(lines).strip() + "\n"


def summary_markdown_tables(data: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    """Build human-facing markdown summary tables."""
    frame_columns = [
        "frame",
        "time_ps",
        "source",
        "status",
        "n_atoms",
        "n_waters",
        "n_guests",
        "connection_mode",
        "connection_count",
        "hbond_count",
        "oo_connection_count",
        "pair_connection_count",
    ]
    tables: list[tuple[str, pd.DataFrame]] = [
        ("Failures", failure_summary_table(data)),
        ("Frames", summary_simple_table(data, frame_columns)),
        ("Molecules", molecule_summary_table(data)),
        ("Rings", summary_simple_table(data, ["frame", "time_ps", "ring4", "ring5", "ring6", "ring7", "free_ring4", "free_ring5", "free_ring6", "free_ring7"])),
        ("Half Cage", patch_summary_table(data, "half_cage")),
        ("Quasi Cage", patch_summary_table(data, "quasi_cage")),
    ]
    tables.extend(
        [
            ("Cages", cage_summary_table(data)),
            ("Cage Occupancy", cage_occupancy_summary_table(data, markdown_style=True)),
            ("Cage Isomers", cage_isomer_summary_table(data, include_zero_rows=False)),
            ("Order Parameters", order_parameter_summary_table(data)),
            ("Ice", summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"])),
        ]
    )
    return tables


def summary_sheet_tables(data: pd.DataFrame, config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """Build lightweight main-summary tables using the configured scopes."""
    if is_cpp_mode(config.get("mode", DEFAULT_MODE)):
        tables: dict[str, pd.DataFrame] = {
            "failures": failure_summary_table(data),
            "cage": cage_summary_table(data),
            "order_parameter": order_parameter_summary_table(
                data,
                config.get("order", {}).get("parameters", ["f3", "f4"]),
                include_focus=bool(config.get("order", {}).get("focus_waters", [])),
            ),
        }
        return {name: table for name, table in tables.items() if not table.empty}

    ring_sizes = configured_ring_report_sizes(config)
    ring_columns = ["frame", "time_ps"]
    for size in ring_sizes:
        ring_columns.extend([f"ring{size}", f"free_ring{size}"])
    tables: dict[str, pd.DataFrame] = {
        "failures": failure_summary_table(data),
        connection_sheet_name(data): connection_summary_table(data),
        "ring": summary_simple_table(data, ring_columns),
        "half_cage": patch_summary_table(data, "half_cage"),
        "quasi_cage": patch_summary_table(data, "quasi_cage"),
        "cage": cage_summary_table(data),
        "hydrate_cluster": hydrate_cluster_summary_table(data),
        "order_parameter": order_parameter_summary_table(
            data,
            config.get("order", {}).get("parameters", ["f3", "f4"]),
            include_focus=bool(config.get("order", {}).get("focus_waters", [])),
        ),
        "ice": summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"]),
    }
    return {name: table for name, table in tables.items() if not table.empty}


def failure_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Return one diagnostic row per failed input frame."""
    if "status" not in data.columns:
        return pd.DataFrame()
    columns = [
        column
        for column in ("frame", "time_ps", "source", "status", "error")
        if column in data.columns
    ]
    failed = data.loc[
        data["status"].astype(str).str.lower().eq("failed"),
        columns,
    ]
    return failed.reset_index(drop=True)


def summary_detail_tables(
    data: pd.DataFrame,
    config: dict[str, Any],
    *,
    raw_rows: list[dict[str, Any]] | None = None,
) -> dict[str, pd.DataFrame]:
    """Build detailed CSV tables selected independently from the main summary."""
    include_zero_isomers = (
        str(config.get("output", {}).get("cage_isomer_rows", "nonzero")).lower()
        == "all"
    )
    tables: dict[str, pd.DataFrame] = {}
    if output_enabled(config, "summary-detail-csv"):
        tables.update(
            {
                "cage_occupancy": cage_occupancy_summary_table(
                    data,
                    markdown_style=False,
                ),
                "cage_isomer": cage_isomer_summary_table(
                    data,
                    include_zero_rows=include_zero_isomers,
                ),
            }
        )
        if not is_cpp_mode(config.get("mode", DEFAULT_MODE)):
            tables["quasi_cage_isomer"] = quasi_cage_isomer_summary_table(
                data,
                raw_rows=raw_rows,
            )
    cluster_detail_enabled = output_enabled(config, "cluster-detail")
    if cluster_detail_enabled:
        tables["hydrate_domain"] = hydrate_domain_table(data)
        tables["hydrate_cluster_detail"] = hydrate_cluster_detail_table(data)
    keep_empty: set[str] = set()
    if cluster_detail_enabled and hydrate_cluster_is_enabled(data):
        keep_empty.add("hydrate_domain")
        keep_empty.add("hydrate_cluster_detail")
    return {
        name: table
        for name, table in tables.items()
        if not table.empty or name in keep_empty
    }


def hydrate_cluster_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Build the per-frame hydrate_cluster summary sheet."""
    if not hydrate_cluster_is_enabled(data):
        return pd.DataFrame()
    columns = [
        "frame",
        "time_ps",
        "hydrate_cluster_count",
        "sI_cluster_count",
        "sII_cluster_count",
        "sH_cluster_count",
        "mixed_cluster_count",
        "unclassified_cluster_count",
        "hydrate_domain_count",
        "sI_domain_count",
        "sII_domain_count",
        "sH_domain_count",
        "classified_cage_count",
        "boundary_cage_count",
        "ambiguous_cage_count",
        "unclassified_cage_count",
        "isolated_cage_count",
        "largest_cluster_cage_count",
        "largest_cluster_water_count",
        "cluster_size_distribution",
    ]
    return summary_simple_table(data, columns)


def hydrate_cluster_is_enabled(data: pd.DataFrame) -> bool:
    """Return whether any frame requested hydrate_cluster reporting."""
    if "hydrate_cluster_enabled" not in data.columns:
        return False
    return bool(data["hydrate_cluster_enabled"].astype(str).str.lower().eq("on").any())


def hydrate_cluster_detail_table(data: pd.DataFrame) -> pd.DataFrame:
    """Expand stored per-frame cluster details into one row per cluster."""
    return expanded_hydrate_table(data, "hydrate_cluster_detail")


def hydrate_domain_table(data: pd.DataFrame) -> pd.DataFrame:
    """Expand stored per-frame domain details into one row per domain."""
    return expanded_hydrate_table(data, "hydrate_domain_detail")


def hydrate_motif_table(data: pd.DataFrame) -> pd.DataFrame:
    """Expand stored per-frame motif evidence into one row per motif."""
    return expanded_hydrate_table(data, "hydrate_motif_detail")


def expanded_hydrate_table(data: pd.DataFrame, column: str) -> pd.DataFrame:
    """Expand one stored list of hydrate records into a flat workbook table."""
    schemas = {
        "hydrate_cluster_detail": [
            "cluster_id",
            "hydrate_type",
            "cage_count",
            "water_count",
            "guest_count",
            "empty_cage_count",
            "occupied_cage_count",
            "classified_cage_count",
            "boundary_cage_count",
            "ambiguous_cage_count",
            "unclassified_cage_count",
            "classified_cage_fraction",
            "domain_count",
            "cage_type_counts",
            "cage_composition",
            "boundary_composition",
            "guest_composition",
            "domain_ids",
            "cage_ids",
            "classified_cage_ids",
            "boundary_cage_ids",
            "ambiguous_cage_ids",
            "unclassified_cage_ids",
            "shared_face_count",
        ],
        "hydrate_domain_detail": [
            "domain_id",
            "cluster_id",
            "hydrate_type",
            "status",
            "cage_count",
            "seed_count",
            "seed_cage_count",
            "expanded_cage_count",
            "classified_fraction",
            "water_count",
            "guest_count",
            "external_boundary_contact_count",
            "cage_composition",
            "guest_composition",
            "cage_ids",
            "seed_cage_ids",
            "external_boundary_contact_ids",
        ],
        "hydrate_motif_detail": [
            "motif_id",
            "cluster_id",
            "domain_id",
            "hydrate_type",
            "status",
            "completeness",
            "consistency",
            "confidence",
            "cage_count",
            "core_cage_count",
            "support_cage_count",
            "cage_composition",
            "core_cage_composition",
            "core_cage_ids",
            "motif_cage_ids",
            "internal_shared_face_count",
            "internal_shared_face_ids",
            "classification_method",
        ],
    }
    columns = ["frame", "time_ps", *schemas.get(column, [])]
    if column not in data.columns:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for _, record in data.iterrows():
        details = record.get(column, [])
        if not isinstance(details, list):
            continue
        for item in details:
            row = {
                "frame": record.get("frame", ""),
                "time_ps": record.get("time_ps", ""),
            }
            row.update(item)
            rows.append(row)
    return pd.DataFrame(rows).reindex(columns=columns)


def connection_sheet_name(data: pd.DataFrame) -> str:
    """Name the connection sheet after the active graph mode."""
    modes = [str(value) for value in data.get("connection_mode", pd.Series(dtype=str)).dropna() if str(value)]
    mode = modes[0] if modes else "connection"
    if mode == "hbond":
        return "hbond"
    if mode == "oo":
        return "oo_connection"
    if mode == "pairs":
        return "pair_connection"
    return "connection"


def connection_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Build a per-frame connection and coordination diagnostic table."""
    columns = [
        "frame",
        "time_ps",
        "connection_mode",
        "connection_count",
        "mean_coordination",
        "coordination_0",
        "coordination_1",
        "coordination_2",
        "coordination_3",
        "coordination_4",
        "coordination_gt4",
        "coordination_0_fraction",
        "coordination_1_fraction",
        "coordination_2_fraction",
        "coordination_3_fraction",
        "coordination_4_fraction",
        "coordination_gt4_fraction",
        "degree_le2_fraction",
        "degree4_fraction",
        "over4_fraction",
    ]
    modes = [str(value) for value in data.get("connection_mode", pd.Series(dtype=str)).dropna() if str(value)]
    mode = modes[0] if modes else ""
    mode_column = {
        "hbond": "hbond_count",
        "oo": "oo_connection_count",
        "pairs": "pair_connection_count",
    }.get(mode)
    if mode_column:
        columns.insert(4, mode_column)
    return summary_simple_table(data, columns)


def summary_simple_table(data: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Return a table with only columns present in the run data."""
    existing = [column for column in columns if column in data.columns]
    return data.loc[:, existing].copy() if existing else pd.DataFrame()


def order_parameter_summary_table(
    data: pd.DataFrame,
    parameters: Any | None = None,
    *,
    include_focus: bool | None = None,
) -> pd.DataFrame:
    """Build the per-frame F3/F4/Q_l and hydrate order-parameter table."""
    selected = (
        normalize_order_parameters(parameters)
        if parameters is not None
        else infer_order_parameters_from_data(data)
    )
    if not selected:
        return pd.DataFrame()
    if include_focus is None:
        focus_count_columns = [
            column
            for column in data.columns
            if str(column).endswith("_focus_count")
        ]
        include_focus = any(
            bool((pd.to_numeric(data[column], errors="coerce").fillna(0) > 0).any())
            for column in focus_count_columns
        )
    columns = ["frame", "time_ps"]
    for name in selected:
        if name == "f3":
            columns.extend(["F3_mean", "F3_count"])
        elif name == "f4":
            columns.extend(["F4_mean", "F4_count"])
        elif name.startswith("q"):
            columns.extend([f"{name}_mean", f"{name}_count"])
        elif name == "mcg1":
            columns.append("MCG1_largest_cluster")
        elif name == "mcg3":
            columns.append("MCG3_largest_cluster")
        elif name == "dhop35":
            columns.append("DHOP35_largest_cluster")
        elif name == "dhop30":
            columns.append("DHOP30_largest_cluster")
    if include_focus:
        for name in selected:
            if name == "f3":
                columns.extend(["F3_focus_mean", "F3_focus_count"])
            elif name == "f4":
                columns.extend(["F4_focus_mean", "F4_focus_count"])
            elif name.startswith("q"):
                columns.extend([f"{name}_focus_mean", f"{name}_focus_count"])
    table = data.reindex(columns=columns).copy()
    for column in (
        "MCG1_largest_cluster",
        "MCG3_largest_cluster",
        "DHOP35_largest_cluster",
        "DHOP30_largest_cluster",
    ):
        if column in table:
            table[column] = table[column].where(table[column].notna(), "N/A")
    return table.rename(
        columns={
            "MCG1_largest_cluster": "MCG-1",
            "DHOP35_largest_cluster": "DHOP35",
            "MCG3_largest_cluster": "MCG-3",
            "DHOP30_largest_cluster": "DHOP30",
        }
    )


def infer_order_parameters_from_data(data: pd.DataFrame) -> tuple[str, ...]:
    """Infer legacy table selection when a caller does not provide config."""
    inferred: list[str] = []
    for name, count_column in (("f3", "F3_count"), ("f4", "F4_count")):
        if count_column in data and not data[count_column].replace("", pd.NA).isna().all():
            inferred.append(name)
    inferred.extend(f"q{degree}" for degree in q_degree_from_data(data))
    for name, column in (
        ("mcg1", "MCG1_largest_cluster"),
        ("mcg3", "MCG3_largest_cluster"),
        ("dhop35", "DHOP35_largest_cluster"),
        ("dhop30", "DHOP30_largest_cluster"),
    ):
        if column in data and not data[column].replace("", pd.NA).isna().all():
            inferred.append(name)
    return normalize_order_parameters(inferred or ["none"])


def q_degree_from_data(data: pd.DataFrame) -> list[int]:
    """Infer reported Q_l degree values from summary columns."""
    degree_values: set[int] = set()
    for column in data.columns:
        match = re.fullmatch(r"q(\d+)_(?:mean|count)", str(column))
        if match and not data[column].replace("", pd.NA).isna().all():
            degree_values.add(int(match.group(1)))
    return sorted(degree_values)


def molecule_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Build the global molecule-count table using source-file residue order."""
    molecule_columns = [column for column in data.columns if column.startswith("mol_")]
    if not molecule_columns:
        return pd.DataFrame()
    output = pd.DataFrame({"frame": data["frame"], "time_ps": data["time_ps"]})
    for column in molecule_columns:
        output[column.removeprefix("mol_")] = [count_cell(value) for value in data[column]]
    return output


def patch_summary_table(data: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Build a plotting-friendly half_cage or quasi_cage table."""
    columns = [column for column in data.columns if column.startswith(f"{prefix}_") and column not in {f"{prefix}_total", f"{prefix}_breakdown"}]
    if not columns:
        return pd.DataFrame()
    sorted_columns = sorted(columns)
    labels = [patch_summary_label(column.removeprefix(f"{prefix}_"), prefix) for column in sorted_columns]
    output_labels = sorted(set(labels))
    rows: list[dict[str, Any]] = []
    for _, record in data.iterrows():
        row: dict[str, Any] = {
            "frame": record.get("frame", ""),
            "time_ps": record.get("time_ps", ""),
        }
        total = 0
        for column, label in zip(sorted_columns, labels):
            count = count_cell(record.get(column, 0))
            row[label] = row.get(label, 0) + count
            total += count
        row["total"] = total
        rows.append(row)
    return pd.DataFrame(rows).reindex(columns=["frame", "time_ps", *output_labels, "total"])


def patch_summary_label(patch_type: str, prefix: str) -> str:
    """Return the workbook label for an open-patch type."""
    if prefix != "quasi_cage":
        return patch_type
    return patch_composition_label(patch_type)


def patch_composition_label(patch_type: str) -> str:
    """Return a patch label without internal prefix or ring-sequence isomer marks."""
    return patch_display_label(patch_type).translate(SUBSCRIPT_DIGIT_DELETE)


def summary_cage_types_from_data(data: pd.DataFrame) -> list[str]:
    """Collect exact cage report types, including requested zero-count types."""
    requested: set[str] = set()
    report_all = False
    for value in data.get("cage_report_types", pd.Series(dtype=str)).dropna():
        marker = str(value).strip()
        if not marker:
            continue
        if marker.lower() == "all":
            report_all = True
            continue
        requested.update(item for item in marker.split(";") if item)
    if requested and not report_all:
        return [item for item in ordered_cage_types(requested) if item in requested]

    detected: set[str] = set()
    for column in data.columns:
        if not column.startswith("cage_"):
            continue
        label = column.removeprefix("cage_")
        if label in {"empty", "occupied", "total", "report_types"} or "_" in label:
            continue
        values = pd.to_numeric(data[column], errors="coerce").fillna(0)
        if bool((values > 0).any()):
            detected.add(label)
    return [item for item in ordered_cage_types(detected) if item in detected]


def cage_summary_table(data: pd.DataFrame) -> pd.DataFrame:
    """Build the global cage-count table with superscript cage headers."""
    cage_types = summary_cage_types_from_data(data)
    output = pd.DataFrame({"frame": data["frame"], "time_ps": data["time_ps"]})
    for cage_type in cage_types:
        output[cage_display_label(cage_type)] = [count_cell(value) for value in data.get(f"cage_{cage_type}", pd.Series([0] * len(data)))]
    output["total"] = output[[cage_display_label(cage_type) for cage_type in cage_types]].sum(axis=1)
    return output


def cage_occupancy_summary_table(data: pd.DataFrame, markdown_style: bool) -> pd.DataFrame:
    """Build global cage occupancy rows, optionally with tree markers."""
    cage_types = summary_cage_types_from_data(data)
    guest_labels = global_guest_labels(data, cage_types)
    rows: list[dict[str, Any]] = []
    for _, record in data.iterrows():
        child_labels = row_guest_labels(record, guest_labels, cage_types)
        if sum(cage_count(record, cage_type, "multi") for cage_type in cage_types) > 0:
            child_labels.append("multi")
        for label in ["empty", "occupied", *child_labels]:
            counts = [cage_count(record, cage_type, label) for cage_type in cage_types]
            branch = ""
            display_label = label
            if label not in {"empty", "occupied"} and markdown_style:
                branch = "└" if label == child_labels[-1] else "├"
                display_label = f"{branch} {label}"
            row: dict[str, Any] = {
                "frame": record.get("frame", ""),
                "time_ps": record.get("time_ps", ""),
                "occupancy": display_label,
            }
            if not markdown_style and label not in {"empty", "occupied"}:
                row["level"] = "detail"
            elif not markdown_style:
                row["level"] = "class"
            for cage_type, count in zip(cage_types, counts):
                row[cage_display_label(cage_type)] = f"{branch} {count}" if branch else count
            row["total"] = f"{branch} {sum(counts)}" if branch else sum(counts)
            rows.append(row)
        total_counts = [cage_count(record, cage_type, "empty") + cage_count(record, cage_type, "occupied") for cage_type in cage_types]
        row = {"frame": record.get("frame", ""), "time_ps": record.get("time_ps", ""), "occupancy": "total"}
        if not markdown_style:
            row["level"] = "total"
        for cage_type, count in zip(cage_types, total_counts):
            row[cage_display_label(cage_type)] = count
        row["total"] = sum(total_counts)
        rows.append(row)
    return pd.DataFrame(rows)


def global_guest_labels(data: pd.DataFrame, cage_types: list[str]) -> list[str]:
    """Collect guest residue labels from occupancy columns in first-seen order."""
    labels: list[str] = []
    for value in data.get("guest_order", pd.Series(dtype=str)):
        if pd.isna(value) or value == "":
            continue
        for label in str(value).split(";"):
            if label and label not in labels:
                labels.append(label)
    for column in data.columns:
        label = cage_occupancy_label_from_column(column, cage_types)
        if label and label not in {"empty", "occupied", "multi"} and label not in labels:
            labels.append(label)
    return labels


def row_guest_labels(record: pd.Series, labels: list[str], cage_types: list[str]) -> list[str]:
    """Return guest labels that occur in this frame, preserving global order."""
    result: list[str] = []
    for label in labels:
        if any(cage_count(record, cage_type, label) for cage_type in cage_types):
            result.append(label)
    return result


def cage_occupancy_label_from_column(column: str, cage_types: list[str]) -> str | None:
    """Extract MET from cage_512_MET style columns."""
    if not column.startswith("cage_"):
        return None
    label = None
    for cage_type in sorted(cage_types, key=len, reverse=True):
        prefix = f"cage_{cage_type}_"
        if column.startswith(prefix):
            label = column.removeprefix(prefix)
            break
    if label is None:
        return None
    if label.startswith("isomer") or label == "isomers":
        return None
    return label


def cage_count(record: pd.Series, cage_type: str, label: str) -> int:
    """Read one cage occupancy count from a summary row."""
    return count_cell(record.get(f"cage_{cage_type}_{label}", 0))


def cage_isomer_summary_table(data: pd.DataFrame, include_zero_rows: bool) -> pd.DataFrame:
    """Build global cage-isomer rows with separate columns."""
    cage_types = summary_cage_types_from_data(data)
    labels = global_cage_isomer_labels(data, cage_types)
    if not labels:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for _, record in data.iterrows():
        for label in labels:
            row: dict[str, Any] = {
                "frame": record.get("frame", ""),
                "time_ps": record.get("time_ps", ""),
                "isomer": label,
            }
            total = 0
            for cage_type in cage_types:
                column = f"cage_{cage_type}_isomer_{label}"
                if column not in data.columns:
                    row[cage_display_label(cage_type)] = "-"
                    continue
                count = count_cell(record.get(column, 0))
                row[cage_display_label(cage_type)] = count
                total += count
            if total == 0 and not include_zero_rows:
                continue
            row["total"] = total
            rows.append(row)
        total_row: dict[str, Any] = {"frame": record.get("frame", ""), "time_ps": record.get("time_ps", ""), "isomer": "total"}
        total = 0
        for cage_type in cage_types:
            count = count_cell(record.get(f"cage_{cage_type}", 0))
            total_row[cage_display_label(cage_type)] = count
            total += count
        total_row["total"] = total
        rows.append(total_row)
    return pd.DataFrame(rows)


def global_cage_isomer_labels(data: pd.DataFrame, cage_types: list[str]) -> list[str]:
    """Collect cage-isomer labels in canonical cage-type order."""
    labels: list[str] = []
    for cage_type in cage_types:
        prefix = f"cage_{cage_type}_isomer_"
        for column in data.columns:
            if column.startswith(prefix):
                label = column.removeprefix(prefix)
                if label not in labels:
                    labels.append(label)
    return labels


def quasi_cage_isomer_summary_table(
    data: pd.DataFrame,
    *,
    raw_rows: list[dict[str, Any]] | None = None,
) -> pd.DataFrame:
    """Build long-form quasi-cage isomer rows for CSV detail output."""
    output_columns = ["frame", "time_ps", "quasi_cage_type", "isomer", "count"]
    rows: list[dict[str, Any]] = []
    if raw_rows is not None:
        for record in raw_rows:
            details = record.get(QUASI_ISOMER_DETAIL_KEY, ())
            for isomer, count in details:
                count = count_cell(count)
                if count:
                    rows.append(
                        {
                            "frame": record.get("frame", ""),
                            "time_ps": record.get("time_ps", ""),
                            "quasi_cage_type": patch_composition_label(isomer),
                            "isomer": patch_display_label(isomer),
                            "count": count,
                        }
                    )
        return pd.DataFrame(rows, columns=output_columns)

    # Support legacy wide-row callers.
    columns = [
        column
        for column in data.columns
        if column.startswith("quasi_cage_") and column not in {"quasi_cage_total", "quasi_cage_breakdown"}
    ]
    for _, record in data.iterrows():
        for column in sorted(columns):
            count = count_cell(record.get(column, 0))
            if count == 0:
                continue
            isomer = column.removeprefix("quasi_cage_")
            rows.append(
                {
                    "frame": record.get("frame", ""),
                    "time_ps": record.get("time_ps", ""),
                    "quasi_cage_type": patch_composition_label(isomer),
                    "isomer": patch_display_label(isomer),
                    "count": count,
                }
            )
    return pd.DataFrame(rows, columns=output_columns)


def count_cell(value: Any) -> int:
    """Normalize missing numeric count cells to zero."""
    if value is None or value == "":
        return 0
    try:
        if pd.isna(value):
            return 0
    except TypeError:
        pass
    return int(value)


def markdown_table(data: pd.DataFrame) -> str:
    """Render a pandas DataFrame as a simple GitHub-style table."""
    headers = [str(col) for col in data.columns]
    body = [[format_summary_cell(value) for value in row] for row in data.itertuples(index=False, name=None)]
    widths = [len(header) for header in headers]
    for row in body:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    lines = [
        "| " + " | ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)) + " |",
        "| " + " | ".join("-" * widths[idx] for idx in range(len(headers))) + " |",
    ]
    for row in body:
        lines.append("| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + " |")
    return "\n".join(lines) + "\n"


def format_summary_cell(value: Any) -> str:
    """Format markdown cells while keeping count columns readable."""
    if isinstance(value, (list, tuple, dict)):
        return repr(value)
    if pd.isna(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def excel_scalar(value: Any) -> Any:
    """Convert containers into readable scalar values for XLSX cells."""
    if isinstance(value, (list, tuple, set)):
        return ",".join(str(item) for item in value)
    if isinstance(value, dict):
        return repr(value)
    return value


def flatten_config(config: dict[str, Any], prefix: str = "") -> list[dict[str, str]]:
    """Flatten nested config keys for the main summary config table."""
    rows = []
    for key, value in config.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            rows.extend(flatten_config(value, name))
        else:
            rows.append({"parameter": name, "value": repr(value)})
    return rows
