from __future__ import annotations

"""Per-frame Markdown and auxiliary machine-readable outputs."""

from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import re
from typing import Any
import unicodedata

from ... import __release_date__, __version__
from ...config import (
    DEFAULT_MODE,
    is_cpp_mode,
    normalize_order_parameters,
    output_enabled,
    q_degrees_from_order_parameters,
)
from ...core.cage import parse_cage_face_label
from ...display import graph_mode_display
from ...models import CagePatch, FrameResult
from ..occupancy import guest_composition_label, guest_lookup as build_guest_lookup
from .tables import (
    atom_resname_counts,
    cage_display_label,
    format_summary_cell,
    guest_composition_sort_key,
    guest_resname_order,
    format_cage_type_counts,
    hydrate_cluster_detail_records,
    hydrate_domain_records,
    hydrate_motif_records,
    ordered_cage_types,
    patch_breakdown,
    patch_counts,
    patch_display_label,
    present_cage_types,
    result_row,
    source_label,
)
from .models import SUMMARY_COLUMNS

import pandas as pd

TREE_MIDDLE = "\u251c"
TREE_LAST = "\u2514"
TREE_PIPE = "\u2502"
SUBSCRIPT_DIGIT_DELETE = dict.fromkeys(range(0x2080, 0x208A))


def hydrate_cluster_info_section(
    result: FrameResult,
    row: dict[str, Any] | None = None,
) -> list[str]:
    """Render the compact, mutually exclusive cluster hierarchy."""
    if not result.hydrate_cluster_enabled:
        return []

    cage_by_id = {
        cage.object_id: cage for cage in (result.all_cages or result.cages)
    }
    domains_by_cluster: dict[str, list[Any]] = defaultdict(list)
    for domain in result.hydrate_domains:
        domains_by_cluster[domain.cluster_id].append(domain)

    rows: list[list[Any]] = []
    displayed_cluster_ids: set[str] = set()
    for cluster in sorted(result.hydrate_clusters, key=lambda item: item.object_id):
        cluster_ids = tuple(dict.fromkeys(cluster.cage_ids))
        cluster_id_set = set(cluster_ids)
        displayed_cluster_ids.update(cluster_id_set)
        remaining_ids = set(cluster_ids)
        rows.append(
            [cluster.object_id, cluster.hydrate_type, len(cluster_id_set)]
        )

        children: list[tuple[str, str, tuple[str, ...]]] = []
        domains = sorted(
            domains_by_cluster.get(cluster.object_id, []),
            key=lambda item: item.object_id,
        )
        for domain in domains:
            domain_ids = tuple(
                cage_id
                for cage_id in dict.fromkeys(domain.cage_ids)
                if cage_id in remaining_ids
            )
            if not domain_ids:
                continue
            remaining_ids.difference_update(domain_ids)
            children.append(
                (domain.object_id, domain.hydrate_type, domain_ids)
            )

        boundary_ids = tuple(
            cage_id
            for cage_id in dict.fromkeys(cluster.boundary_cage_ids)
            if cage_id in remaining_ids
        )
        if boundary_ids:
            remaining_ids.difference_update(boundary_ids)
            children.append(("boundary", "boundary", boundary_ids))

        unclassified_ids = tuple(
            cage_id for cage_id in cluster_ids if cage_id in remaining_ids
        )
        if unclassified_ids:
            children.append(
                ("unclassified", "unclassified", unclassified_ids)
            )

        for child_index, (name, hydrate_type, cage_ids) in enumerate(children):
            child_branch = (
                TREE_LAST if child_index == len(children) - 1 else TREE_MIDDLE
            )
            rows.append(
                [
                    f"{child_branch} {name}",
                    hydrate_type,
                    f"{child_branch} {len(cage_ids)}",
                ]
            )
            type_counts = Counter(
                cage_by_id[cage_id].cage_type
                for cage_id in cage_ids
                if cage_id in cage_by_id
            )
            cage_types = present_cage_types(type_counts)
            for type_index, cage_type in enumerate(cage_types):
                type_branch = (
                    TREE_LAST if type_index == len(cage_types) - 1 else TREE_MIDDLE
                )
                rows.append(
                    [
                        f"  {type_branch} {cage_display_label(cage_type)}",
                        "",
                        f"  {type_branch} {type_counts[cage_type]}",
                    ]
                )

    isolated_count = len(
        set(result.isolated_cage_ids).difference(displayed_cluster_ids)
    )
    if isolated_count:
        rows.append(["isolated", "isolated", isolated_count])

    if not rows:
        return ["", "## Hydrate Cluster", "", "no hydrate cluster"]
    return section_table(
        "Hydrate Cluster",
        ["item", "type", "cage_qty"],
        rows,
    )


def hydrate_cluster_detail_section(records: list[dict[str, Any]]) -> list[str]:
    """Render one vertical table per cluster, including cage composition."""
    if not records:
        return section_table(
            "Hydrate Cluster Detail",
            ["item", "value"],
            [["cluster_count", 0]],
        )
    lines = ["", "## Hydrate Cluster Detail", ""]
    metric_keys = (
        "hydrate_type",
        "cage_count",
        "classified_cage_count",
        "boundary_cage_count",
        "ambiguous_cage_count",
        "unclassified_cage_count",
        "classified_cage_fraction",
        "domain_count",
        "water_count",
        "guest_count",
        "empty_cage_count",
        "occupied_cage_count",
        "boundary_composition",
        "guest_composition",
    )
    for record in records:
        rows = [[key, record.get(key, "")] for key in metric_keys]
        type_counts = record.get("cage_type_counts", {})
        if type_counts:
            rows.append(["cage composition", ""])
            ordered_types = present_cage_types(type_counts)
            for index, cage_type in enumerate(ordered_types):
                branch = TREE_LAST if index == len(ordered_types) - 1 else TREE_MIDDLE
                rows.append(
                    [f"{branch} {cage_display_label(cage_type)}", type_counts[cage_type]]
                )
        lines.extend(object_vertical_table(str(record["cluster_id"]), rows))
    return lines


def hydrate_cluster_hierarchy_section(result: FrameResult) -> list[str]:
    """Render mutually exclusive cluster categories."""
    cage_by_id = {cage.object_id: cage for cage in (result.all_cages or result.cages)}
    domains_by_cluster: dict[str, list[Any]] = defaultdict(list)
    for domain in result.hydrate_domains:
        domains_by_cluster[domain.cluster_id].append(domain)
    if not result.hydrate_clusters:
        return ["", "## Hydrate Cluster Hierarchy", "", "no hydrate cluster"]

    lines = ["", "## Hydrate Cluster Hierarchy", ""]
    for cluster_index, cluster in enumerate(result.hydrate_clusters):
        if cluster_index:
            lines.append("")
        cluster_counts = Counter(
            cage_by_id[cage_id].cage_type
            for cage_id in cluster.cage_ids
            if cage_id in cage_by_id
        )
        cage_types = present_cage_types(cluster_counts)
        headers = [
            "item",
            "type",
            "cage_qty",
            *(cage_display_label(cage_type) for cage_type in cage_types),
        ]
        rows: list[list[Any]] = [
            hierarchy_table_row(
                cluster.object_id,
                cluster.hydrate_type,
                cluster.cage_ids,
                cage_types,
                cage_by_id,
            )
        ]
        children: list[tuple[str, str, tuple[str, ...]]] = [
            (domain.object_id, domain.hydrate_type, domain.cage_ids)
            for domain in domains_by_cluster.get(cluster.object_id, [])
        ]
        if cluster.boundary_cage_ids:
            children.append(("boundary", "boundary", cluster.boundary_cage_ids))
        if cluster.ambiguous_cage_ids:
            children.append(("ambiguous", "ambiguous", cluster.ambiguous_cage_ids))
        if cluster.unclassified_cage_ids:
            children.append(
                ("unclassified", "unclassified", cluster.unclassified_cage_ids)
            )
        for child_index, (name, hydrate_type, cage_ids) in enumerate(children):
            branch = TREE_LAST if child_index == len(children) - 1 else TREE_MIDDLE
            rows.append(
                hierarchy_table_row(
                    f"{branch} {name}",
                    hydrate_type,
                    cage_ids,
                    cage_types,
                    cage_by_id,
                )
            )
        lines.append(markdown_rows(headers, rows).rstrip())
    return lines


def hierarchy_table_row(
    item: str,
    hydrate_type: str,
    cage_ids: tuple[str, ...],
    cage_types: list[str],
    cage_by_id: dict[str, Any],
) -> list[Any]:
    """Return one hierarchy row; absent cage types are rendered as dashes."""
    counts = Counter(
        cage_by_id[cage_id].cage_type
        for cage_id in cage_ids
        if cage_id in cage_by_id
    )
    return [
        item,
        hydrate_type,
        len(cage_ids),
        *(counts[cage_type] if counts[cage_type] else "-" for cage_type in cage_types),
    ]


def hierarchy_label(
    object_id: str,
    hydrate_type: str,
    cage_ids: tuple[str, ...],
    cage_by_id: dict[str, Any],
) -> str:
    """Return one readable tree label with nonzero cage-type counts."""
    counts = Counter(cage_by_id[cage_id].cage_type for cage_id in cage_ids if cage_id in cage_by_id)
    parts = [hydrate_type] if hydrate_type else []
    parts.append(f"cages={len(cage_ids)}")
    parts.extend(f"{cage_display_label(cage_type)}={counts[cage_type]}" for cage_type in present_cage_types(counts))
    return f"{object_id} [{', '.join(parts)}]"


def domain_hierarchy_label(domain: Any, cage_by_id: dict[str, Any]) -> str:
    """Return a domain label with internal seed and expansion counts."""
    seed_ids = set(domain.seed_cage_ids)
    counts = Counter(cage_by_id[cage_id].cage_type for cage_id in domain.cage_ids if cage_id in cage_by_id)
    parts = [
        domain.hydrate_type,
        domain.status,
        f"cages={domain.cage_count}",
        f"seeds={domain.seed_count}",
        f"seed_cages={len(seed_ids)}",
        f"expanded={domain.cage_count - len(seed_ids)}",
    ]
    parts.extend(f"{cage_display_label(cage_type)}={counts[cage_type]}" for cage_type in present_cage_types(counts))
    return f"{domain.object_id} [{', '.join(parts)}]"


def motif_hierarchy_label(motif: Any, cage_by_id: dict[str, Any]) -> str:
    """Return a compact label for one overlapping local topology motif."""
    counts = Counter(cage_by_id[cage_id].cage_type for cage_id in motif.cage_ids if cage_id in cage_by_id)
    core_count = len(motif.anchor_cage_ids)
    parts = [
        motif.hydrate_type,
        motif.status,
        f"cages={motif.cage_count}",
        f"completeness={motif.completeness:.2f}",
        f"core={core_count}",
        f"support={motif.cage_count - core_count}",
    ]
    parts.extend(f"{cage_display_label(cage_type)}={counts[cage_type]}" for cage_type in present_cage_types(counts))
    return f"{motif.object_id} [{', '.join(parts)}]"


def object_vertical_table(object_id: str, rows: list[list[Any]]) -> list[str]:
    """Render one level-three object heading followed by a narrow table."""
    return [f"### {object_id}", "", markdown_rows(["item", object_id], rows).rstrip(), ""]


def hydrate_domain_info_section(
    result: FrameResult,
    guests_by_id: dict[str, Any],
    guest_order: list[str],
) -> list[str]:
    """Render one vertical detail table per hydrate domain."""
    records = hydrate_domain_records(result, guests_by_id, guest_order)
    if not records:
        return section_table("Hydrate Domain", ["item", "value"], [["domain_count", 0]])
    lines = ["", "## Hydrate Domain", ""]
    keys = (
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
    )
    for record in records:
        rows = [[key, record.get(key, "")] for key in keys]
        lines.extend(object_vertical_table(str(record["domain_id"]), rows))
    return lines


def hydrate_motif_info_section(result: FrameResult) -> list[str]:
    """Render one vertical evidence table per hydrate motif."""
    records = hydrate_motif_records(result)
    if not records:
        return section_table("Hydrate Motif", ["item", "value"], [["motif_count", 0]])
    lines = ["", "## Hydrate Motif", ""]
    keys = (
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
        "internal_shared_face_count",
        "classification_method",
    )
    for record in records:
        rows = [[key, record.get(key, "")] for key in keys]
        lines.extend(object_vertical_table(str(record["motif_id"]), rows))
    return lines


def hydrate_boundary_info_section(result: FrameResult) -> list[str]:
    """Render exclusive boundary totals and composition."""
    clusters = [
        cluster for cluster in result.hydrate_clusters if cluster.boundary_cage_ids
    ]
    if not clusters:
        return section_table(
            "Hydrate Boundary",
            ["item", "value"],
            [["boundary_cage_count", 0]],
        )
    lines = ["", "## Hydrate Boundary", ""]
    cage_by_id = {cage.object_id: cage for cage in (result.all_cages or result.cages)}
    for cluster in clusters:
        counts = Counter(
            cage_by_id[cage_id].cage_type
            for cage_id in cluster.boundary_cage_ids
            if cage_id in cage_by_id
        )
        rows = [
            ["boundary_cage_count", len(cluster.boundary_cage_ids)],
            ["boundary_composition", format_cage_type_counts(counts)],
        ]
        lines.extend(object_vertical_table(cluster.object_id, rows))
    return lines


def append_tree_value_rows(rows: list[list[Any]], label: str, item_label: str, values: list[str]) -> None:
    """Append a vertical, consistently branched list to a two-column table."""
    rows.append([label, ""])
    for index, value in enumerate(values):
        branch = TREE_LAST if index == len(values) - 1 else TREE_MIDDLE
        rows.append([f"{branch} {item_label}", value])


def split_record_ids(value: Any) -> list[str]:
    """Split a semicolon-delimited workbook field for vertical info output."""
    return [item for item in str(value).split(";") if item]


def failed_row(frame_name: str, source: str, error: str) -> dict[str, Any]:
    """Create a summary row for a skipped or failed frame."""
    row = {column: "" for column in SUMMARY_COLUMNS}
    row.update({"frame": frame_name, "source": source, "status": "failed", "error": error})
    return row


def write_frame_info(
    result: FrameResult,
    frame_dir: Path,
    ring_sizes: list[int] | None = None,
    requested_bond_mode: Any | None = None,
    order_parameters: Any | None = None,
    analysis_mode: Any = DEFAULT_MODE,
    input_metadata: dict[str, Any] | None = None,
) -> None:
    """Write the per-frame Markdown report with inspection-oriented tables."""
    frame_dir.mkdir(parents=True, exist_ok=True)
    row = result_row(result, include_cluster_details=False)
    selected_order_parameters = selected_order_parameters_for_result(
        result,
        order_parameters,
    )
    selected_order_set = set(selected_order_parameters)
    cage_values = {cage.cage_type for cage in result.cages}
    cage_types = [cage_type for cage_type in ordered_cage_types(cage_values) if cage_type in cage_values]
    default_ring_sizes = result.ring_report_sizes or tuple(result.rings)
    enabled_ring_sizes = sorted(set(ring_sizes if ring_sizes is not None else default_ring_sizes))
    cpp_mode = is_cpp_mode(analysis_mode)
    metadata = input_metadata or {}
    metadata_rows = [
        [key, metadata[key]]
        for key in (
            "input_format", "topology", "sampling_interval",
            "native_frame_interval_ps", "delta_time_ps", "raw_frame_step",
            "selected_frames", "find_half", "find_quasi", "lammps_units",
            "lammps_timestep", "lammps_atom_style", "lammps_type_map_source",
            "graph_mode_reason",
        )
        if metadata.get(key) not in (None, "", "<none>")
    ]
    frame_information_rows = [
        ["sqq version", __version__],
        ["sqq engine", "sqq-cpp" if cpp_mode else "sqq-py"],
        ["date & time", report_datetime_label()],
        ["source", source_label(result.frame.source)],
        *metadata_rows,
        ["frame", result.frame.name],
        ["time_ps", result.frame.time_ps],
        ["graph_mode", graph_mode_display(requested_bond_mode or row["connection_mode"], [row["connection_mode"]])],
        ["bond_mode", row["connection_mode"]],
        ["ring_sizes", ", ".join(str(size) for size in enabled_ring_sizes)],
    ]
    if not cpp_mode:
        frame_information_rows.append(["find_cluster", "on" if result.hydrate_cluster_enabled else "off"])
    frame_information_rows.extend([
        ["status", "ok"],
        ["n_atoms", len(result.frame.atoms)],
        ["n_waters", len(result.waters)],
        ["n_guests", len(result.guests)],
    ])
    lines = [
        f"# SQQ Frame Report: {result.frame.name}",
        "",
        *section_table("Frame Information", ["item", "value"], frame_information_rows),
    ]
    lines.extend(section_table("Molecules", ["resname", "molecules", "atoms"], molecule_count_rows(result)))
    lines.extend(connection_info_section(result, row))

    if not cpp_mode:
        ring_rows = [
            [size, row.get(f"ring{size}", 0), row.get(f"free_ring{size}", 0)]
            for size in enabled_ring_sizes
        ]
        ring_rows.append(
            [
                "total",
                sum(int(item[1]) for item in ring_rows),
                sum(int(item[2]) for item in ring_rows),
            ]
        )
        lines.extend(section_table("Ring", ["ring size", "total", "free"], ring_rows))
        if str(metadata.get("find_half", "on")).lower() == "on":
            lines.extend(patch_info_section("Half Cage", result.half_cages))
        if str(metadata.get("find_quasi", "on")).lower() == "on":
            lines.extend(patch_info_section("Quasi Cage", result.quasi_cages))
            lines.extend(patch_isomer_description_section("Quasi Cage Isomer Description", result.quasi_cages))

    lines.extend(cage_info_section(result, cage_types))
    lines.extend(cage_isomer_description_section(result, cage_types))
    lines.extend(cage_occupancy_section(result, cage_types, evaluated=not cpp_mode or bool(result.guests)))
    if not cpp_mode:
        lines.extend(hydrate_cluster_info_section(result))
    has_focus = result.f3f4 is not None and bool(result.f3f4.focus_resids)
    order_headers = ["metric", "count", "mean"]
    if has_focus:
        order_headers.extend(["focus_count", "focus_mean"])
    order_rows: list[list[Any]] = []
    for name, label, prefix in (
        ("f3", "F3", "F3"),
        ("f4", "F4", "F4"),
    ):
        if name not in selected_order_set:
            continue
        metric_row = [label, row.get(f"{prefix}_count"), row.get(f"{prefix}_mean")]
        if has_focus:
            metric_row.extend(
                [
                    row.get(f"{prefix}_focus_count"),
                    row.get(f"{prefix}_focus_mean"),
                ]
            )
        order_rows.append(metric_row)
    for degree in (() if cpp_mode else q_degrees_from_order_parameters(selected_order_parameters)):
        prefix = f"q{degree}"
        metric_row = [
            f"Q{degree}",
            row.get(f"{prefix}_count"),
            row.get(f"{prefix}_mean"),
        ]
        if has_focus:
            metric_row.extend(
                [
                    row.get(f"{prefix}_focus_count"),
                    row.get(f"{prefix}_focus_mean"),
                ]
            )
        order_rows.append(metric_row)
    lines.extend(section_table("Order Parameters", order_headers, order_rows))

    if not cpp_mode and result.hydrate_order is not None:
        hydrate_order_rows: list[list[Any]] = []
        if "mcg1" in selected_order_set:
            hydrate_order_rows.append(
                ["MCG-1", order_size_label(result.hydrate_order.mcg1.largest_cluster_size), result.hydrate_order.mcg1.member_type]
            )
        if "mcg3" in selected_order_set and result.hydrate_order.mcg3 is not None:
            hydrate_order_rows.append(["MCG-3", order_size_label(result.hydrate_order.mcg3.largest_cluster_size), result.hydrate_order.mcg3.member_type])
        if "dhop35" in selected_order_set:
            hydrate_order_rows.append(
                ["DHOP35", order_size_label(result.hydrate_order.dhop35.largest_cluster_size), result.hydrate_order.dhop35.member_type]
            )
        if "dhop30" in selected_order_set and result.hydrate_order.dhop30 is not None:
            hydrate_order_rows.append(["DHOP30", order_size_label(result.hydrate_order.dhop30.largest_cluster_size), result.hydrate_order.dhop30.member_type])
        lines.extend(
            section_table(
                "Hydrate Nucleation Order Parameters",
                ["parameter", "largest cluster", "member type"],
                hydrate_order_rows,
            )
        )

    if not cpp_mode:
        ice_rows = [
            ["ice_like_waters", row["ice_like_waters"]],
            ["ice_i_waters", row["ice_i_waters"]],
            ["interfacial_ice_waters", row["interfacial_ice_waters"]],
        ]
        lines.extend(section_table("Ice", ["structure", "water molecules"], ice_rows))
    if result.warnings:
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- {warning}" for warning in result.warnings)
    (frame_dir / f"{result.frame.name}_info.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def selected_order_parameters_for_result(
    result: FrameResult,
    value: Any | None,
) -> tuple[str, ...]:
    """Resolve explicit output selection or infer legacy direct-call behavior."""
    if value is not None:
        return normalize_order_parameters(value)
    inferred: list[str] = []
    if result.f3f4 is not None:
        inferred.extend(("f3", "f4"))
        inferred.extend(f"q{degree}" for degree in result.f3f4.q_degree)
    if result.hydrate_order is not None:
        inferred.extend(("mcg1", "dhop35"))
        if result.hydrate_order.mcg3 is not None:
            inferred.append("mcg3")
        if result.hydrate_order.dhop30 is not None:
            inferred.append("dhop30")
    return normalize_order_parameters(inferred or ["none"])


def order_size_label(value: int | None) -> int | str:
    """Render unavailable order parameters explicitly without conflating them with zero."""
    return "N/A" if value is None else value


def section_table(title: str, headers: list[str], rows: list[list[Any]]) -> list[str]:
    """Render one small markdown section."""
    if not rows:
        return []
    lines = ["", f"## {title}", "", markdown_rows(headers, rows).rstrip()]
    return lines


def markdown_rows(headers: list[str], rows: list[list[Any]]) -> str:
    """Render a source-aligned Markdown table using Unicode display widths."""
    text_rows = [[format_summary_cell(value) for value in row] for row in rows]
    header_text = [str(header) for header in headers]
    widths = [max(3, display_width(header)) for header in header_text]
    for row in text_rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], display_width(value))
    lines = [
        "| " + " | ".join(pad_display(header, widths[idx]) for idx, header in enumerate(header_text)) + " |",
        "| " + " | ".join("-" * widths[idx] for idx in range(len(headers))) + " |",
    ]
    for row in text_rows:
        padded = [*row, *([""] * (len(headers) - len(row)))]
        lines.append("| " + " | ".join(pad_display(padded[idx], widths[idx]) for idx in range(len(headers))) + " |")
    return "\n".join(lines) + "\n"


def display_width(value: Any) -> int:
    """Return a practical monospace width for Markdown source alignment."""
    width = 0
    for char in str(value):
        if unicodedata.combining(char):
            continue
        width += 2 if unicodedata.east_asian_width(char) in {"F", "W"} else 1
    return width


def pad_display(value: Any, width: int) -> str:
    """Right-pad text to a requested Unicode display width."""
    text = str(value)
    return text + " " * max(0, width - display_width(text))


def molecule_count_rows(result: FrameResult) -> list[list[Any]]:
    """Count molecules and atoms by residue name using source identity."""
    atom_counts = atom_resname_counts(result)
    molecule_counts: dict[str, int] = {}
    atoms = result.frame.atoms
    explicit_ids = [atom.molecule_id for atom in atoms]
    if any(value is not None for value in explicit_ids):
        if any(value is None for value in explicit_ids):
            raise ValueError("Atom records mix explicit and implicit molecule identities.")
        for resname, _ in dict.fromkeys(
            (atom.resname, int(atom.molecule_id)) for atom in atoms
        ):
            molecule_counts[resname] = molecule_counts.get(resname, 0) + 1
    else:
        previous_residue: tuple[int, str] | None = None
        for atom in atoms:
            residue = (atom.resid, atom.resname)
            if residue != previous_residue:
                molecule_counts[atom.resname] = molecule_counts.get(atom.resname, 0) + 1
                previous_residue = residue
    rows = [[resname, molecule_counts.get(resname, 0), atom_counts[resname]] for resname in atom_counts]
    rows.append(["total", sum(molecule_counts.values()), len(result.frame.atoms)])
    return rows


def connection_info_section(result: FrameResult, row: dict[str, Any]) -> list[str]:
    """Render network coordination diagnostics without modifying graph edges."""
    mode = result.graph.mode
    title, count_label = {
        "hbond": ("Hydrogen-Bond Coordination", "hydrogen bonds"),
        "oo": ("O-O Connectivity Coordination", "O-O connections"),
        "pairs": ("Pair Connectivity Coordination", "user-defined pairs"),
    }.get(mode, ("Network Coordination", "connections"))
    rows: list[list[Any]] = [
        ["water molecules", len(result.waters), ""],
        [count_label, row["connection_count"], ""],
        ["mean coordination", row["mean_coordination"], ""],
    ]
    for degree in range(5):
        rows.append([
            f"degree {degree}",
            row[f"coordination_{degree}"],
            row[f"coordination_{degree}_fraction"],
        ])
    rows.append(["degree >4", row["coordination_gt4"], row["coordination_gt4_fraction"]])
    rows.append(["degree <=2", "", row["degree_le2_fraction"]])
    return section_table(title, ["item", "count/value", "fraction"], rows)


def patch_info_section(title: str, patches: list[CagePatch]) -> list[str]:
    """Render open patches as composition totals with nested isomers."""
    counts = patch_counts(patches)
    type_header = f"{title.lower().replace(' ', '-')} type"
    if not counts:
        return section_table(title, [type_header, "count"], [["total", 0]])

    grouped: dict[str, dict[str, int]] = {}
    for patch_type, count in counts.items():
        label = patch_display_label(patch_type)
        parent = label.translate(SUBSCRIPT_DIGIT_DELETE)
        children = grouped.setdefault(parent, {})
        children[label] = children.get(label, 0) + count

    rows: list[list[Any]] = []
    for parent in sorted(grouped):
        children = grouped[parent]
        rows.append([parent, sum(children.values())])
        if any(child != parent for child in children):
            labels = sorted(children)
            for index, label in enumerate(labels):
                branch = TREE_LAST if index == len(labels) - 1 else TREE_MIDDLE
                rows.append([f"{branch} {label}", f"{branch} {children[label]}"])
    rows.append(["total", sum(counts.values())])
    return section_table(title, [type_header, "count"], rows)


def patch_isomer_description_section(title: str, patches: list[CagePatch]) -> list[str]:
    """Explain each reported patch isomer as a layered ring sequence."""
    if not patches:
        return []
    counts: dict[str, int] = {}
    descriptions: dict[str, str] = {}
    for patch in patches:
        label = patch_display_label(patch.patch_type)
        counts[label] = counts.get(label, 0) + 1
        descriptions[label] = describe_patch_isomer(patch)
    rows = [[label, counts[label], descriptions[label]] for label in sorted(counts)]
    return section_table(title, ["isomer", "count", "description"], rows)


def describe_patch_isomer(patch: CagePatch) -> str:
    """Describe base/L1/L2/L3 ring layers for one quasi-cage label."""
    if not patch.layers:
        return "Layer information is not available."
    parts = [f"base ring: {patch.layers[0].removesuffix('r')}"]
    for index, layer in enumerate(patch.layers[1:], start=1):
        sequence = subscript_digit_text(layer)
        composition = layer.translate(SUBSCRIPT_DIGIT_DELETE)
        if sequence:
            layer_name = "closed side-ring sequence" if index == 1 else "outer-layer ring sequence"
            parts.append(f"L{index}: {layer_name} {sequence} ({composition})")
        else:
            parts.append(f"L{index}: composition {composition}")
    return "; ".join(parts) + "."


def subscript_digit_text(text: str) -> str:
    """Return normal digits from Unicode subscript digits in a label."""
    digits: list[str] = []
    for char in text:
        code = ord(char)
        if 0x2080 <= code <= 0x2089:
            digits.append(str(code - 0x2080))
    return "".join(digits)


def cage_occupancy_section(
    result: FrameResult,
    cage_types: list[str],
    *,
    evaluated: bool = True,
) -> list[str]:
    """Show one cage type per row with dynamic guest-composition columns."""
    if not evaluated:
        return section_table("Cage Occupancy", ["status"], [["not evaluated (no selected guests)"]])
    guests_by_id = build_guest_lookup(result.guests)
    guest_order = guest_resname_order(result)
    counts: dict[str, dict[str, int]] = {
        cage_type: {"empty": 0, "occupied": 0} for cage_type in cage_types
    }
    compositions: set[str] = set()
    for cage in result.cages:
        cage_counts = counts.setdefault(cage.cage_type, {"empty": 0, "occupied": 0})
        cage_counts["occupied" if cage.occupied else "empty"] += 1
        if not cage.occupied:
            continue
        composition = guest_composition_label(cage, guests_by_id, guest_order) or "unknown"
        compositions.add(composition)
        cage_counts[composition] = cage_counts.get(composition, 0) + 1

    guest_labels = [label for label in guest_order if label in compositions]
    extra_guest_labels = sorted(
        compositions.difference(guest_labels),
        key=lambda label: guest_composition_sort_key(label, guest_order),
    )
    child_labels = [*guest_labels, *extra_guest_labels]
    child_markers = [TREE_LAST if index == len(child_labels) - 1 else TREE_MIDDLE for index in range(len(child_labels))]
    headers = [
        "cage type",
        "total",
        "empty",
        "occupied",
        *[f"{marker} {label}" for marker, label in zip(child_markers, child_labels)],
    ]
    table_rows: list[list[Any]] = []
    for cage_type in cage_types:
        cage_counts = counts.get(cage_type, {})
        empty = cage_counts.get("empty", 0)
        occupied = cage_counts.get("occupied", 0)
        table_rows.append(
            [
                cage_display_label(cage_type),
                empty + occupied,
                empty,
                occupied,
                *[f"{marker} {cage_counts.get(label, 0)}" for marker, label in zip(child_markers, child_labels)],
            ]
        )

    total_empty = sum(counts.get(cage_type, {}).get("empty", 0) for cage_type in cage_types)
    total_occupied = sum(counts.get(cage_type, {}).get("occupied", 0) for cage_type in cage_types)
    table_rows.append(
        [
            "total",
            total_empty + total_occupied,
            total_empty,
            total_occupied,
            *[
                f"{marker} {sum(counts.get(cage_type, {}).get(label, 0) for cage_type in cage_types)}"
                for marker, label in zip(child_markers, child_labels)
            ],
        ]
    )
    return section_table("Cage Occupancy", headers, table_rows)


def cage_info_section(result: FrameResult, cage_types: list[str]) -> list[str]:
    """Show cage totals with nested structural isomers in one section."""
    isomers: dict[str, dict[str, int]] = {cage_type: {} for cage_type in cage_types}
    for cage in result.cages:
        label = cage.isomer or "plain"
        type_isomers = isomers.setdefault(cage.cage_type, {})
        type_isomers[label] = type_isomers.get(label, 0) + 1

    rows: list[list[Any]] = []
    for cage_type in cage_types:
        cage_label = cage_display_label(cage_type)
        type_isomers = isomers.get(cage_type, {})
        rows.append([cage_label, sum(type_isomers.values())])
        if len(type_isomers) <= 1:
            continue
        labels = sorted(type_isomers)
        for index, label in enumerate(labels):
            branch = TREE_LAST if index == len(labels) - 1 else TREE_MIDDLE
            rows.append([f"{branch} {cage_label}_{label}", f"{branch} {type_isomers[label]}"])
    rows.append(["total", len(result.cages)])
    return section_table("Cage", ["cage type", "count"], rows)


def cage_isomer_description_section(result: FrameResult, cage_types: list[str]) -> list[str]:
    """Explain each reported cage isomer as a 6-face adjacency pattern."""
    if not result.cages:
        return []
    counts: dict[tuple[str, str], int] = {}
    for cage in result.cages:
        label = cage.isomer or "plain"
        key = (cage.cage_type, label)
        counts[key] = counts.get(key, 0) + 1

    rows: list[list[Any]] = []
    for cage_type in cage_types:
        cage_label = cage_display_label(cage_type)
        labels = sorted(label for type_name, label in counts if type_name == cage_type)
        for label in labels:
            display = cage_label if label == "plain" else f"{cage_label}_{label}"
            rows.append([display, counts[(cage_type, label)], describe_cage_isomer(cage_type, label)])
    return section_table("Cage Isomer Description", ["isomer", "count", "description"], rows)


def describe_cage_isomer(cage_type: str, label: str) -> str:
    """Describe a cage isomer label in human-facing terms."""
    composition = describe_cage_face_composition(cage_type)
    arrangement = describe_hex_adjacency_label(label)
    return f"{composition}; {arrangement}"


def describe_cage_face_composition(cage_type: str) -> str:
    """Return a compact text description of cage face counts."""
    counts = parse_cage_face_label(cage_type)
    if not counts:
        return f"face composition: {cage_display_label(cage_type)}"
    parts = [
        f"{count} {size}-ring face{'s' if count != 1 else ''}"
        for size, count in sorted(counts.items())
        if count > 0
    ]
    return "face composition: " + ", ".join(parts)


def describe_hex_adjacency_label(label: str) -> str:
    """Return a readable explanation of the 6-ring face adjacency label."""
    if label == "plain":
        return "no 6-ring face arrangement isomer is reported"
    descriptions = {
        "6single": "one 6-ring face; no 6-6 shared edge is possible",
        "6adj": "two 6-ring faces share one edge",
        "6pair+single": "three 6-ring faces contain one adjacent pair and one separated single face",
        "6chain3": "three 6-ring faces form a chain with two 6-6 shared edges",
        "6tri3": "three 6-ring faces are mutually adjacent",
        "6pair+2single": "four 6-ring faces contain one adjacent pair and two separated single faces",
        "2x6pair": "four 6-ring faces form two separated adjacent pairs",
        "6chain3+single": "four 6-ring faces contain one three-face chain and one separated single face",
        "6star3": "four 6-ring faces form a star: one face touches three others",
        "6chain4": "four 6-ring faces form a four-face chain",
        "6tri3+single": "four 6-ring faces contain one mutually adjacent triple and one separated single face",
        "6cycle4": "four 6-ring faces form a four-face cycle",
        "6tri3+tail": "four 6-ring faces contain one mutually adjacent triple with one attached tail face",
        "6K4-e": "four 6-ring faces are almost fully connected, with one 6-6 adjacency missing",
        "6K4": "four 6-ring faces are all mutually adjacent",
    }
    if label in descriptions:
        return descriptions[label]
    separated = re.fullmatch(r"(\d+)x6sep", label)
    if separated:
        return f"{separated.group(1)} 6-ring faces are all separated from each other"
    generic = re.fullmatch(r"6n(\d+)e(\d+)d(\d+)", label)
    if generic:
        n_hex, edge_count, degree_text = generic.groups()
        return f"{n_hex} 6-ring faces with {edge_count} shared 6-6 edges; degree sequence {degree_text}"
    return f"6-ring face arrangement label: {label}"


def report_datetime_label(value: datetime | None = None) -> str:
    """Return local report-generation time for per-frame Markdown output."""
    moment = value or datetime.now().astimezone()
    zone = report_timezone_label(moment)
    return f"{moment:%Y-%m-%d %H:%M:%S} {zone}"


def report_timezone_label(value: datetime) -> str:
    """Return a compact human-facing timezone label."""
    name = value.tzname() or "UTC"
    offset = value.utcoffset()
    total_minutes = int(offset.total_seconds() / 60) if offset is not None else 0
    if total_minutes == 480 and name in {"CST", "China Standard Time", "\u4e2d\u56fd\u6807\u51c6\u65f6\u95f4"}:
        return "Asia/Shanghai"
    if offset is None:
        return name
    sign = "+" if total_minutes >= 0 else "-"
    hours, minutes = divmod(abs(total_minutes), 60)
    offset_text = f"UTC{sign}{hours:02d}:{minutes:02d}"
    return name if name and name != "UTC" else offset_text


def write_membership(result: FrameResult, frame_dir: Path) -> None:
    """Write object-to-water membership for plotting and debugging."""
    rows: list[dict[str, Any]] = []
    water_resid_by_oxygen = {water.oxygen: water.resid for water in result.waters}
    reported_ring_sizes = set(result.ring_report_sizes or tuple(result.rings))
    for size, rings in sorted(result.rings.items()):
        if size not in reported_ring_sizes:
            continue
        for ring in rings:
            rows.append(
                {
                    "object_id": ring.object_id,
                    "object_type": f"ring{size}",
                    "center_atom_name": f"R{size}",
                    "oxygen_indices": ",".join(str(idx) for idx in ring.nodes),
                    "water_resids": ",".join(str(water_resid_by_oxygen[idx]) for idx in ring.nodes),
                    "guest_ids": "",
                    "isomer": "",
                    "unwrap_conflict": "false",
                }
            )
    for patch in [*result.half_cages, *result.quasi_cages]:
        rows.append(
            {
                "object_id": patch.object_id,
                "object_type": patch.patch_type,
                "center_atom_name": "HC" if patch.kind == "half_cage" else "QC",
                "oxygen_indices": ",".join(str(idx) for idx in patch.waters),
                "water_resids": ",".join(str(water_resid_by_oxygen[idx]) for idx in patch.waters),
                "guest_ids": "",
                "isomer": "",
                "unwrap_conflict": "false",
            }
        )
    for cage in result.cages:
        rows.append(
            {
                "object_id": cage.object_id,
                "object_type": cage.cage_type,
                "center_atom_name": cage_center_name(cage.cage_type),
                "oxygen_indices": ",".join(str(idx) for idx in cage.waters),
                "water_resids": ",".join(str(water_resid_by_oxygen[idx]) for idx in cage.waters),
                "guest_ids": ",".join(cage.guest_ids),
                "isomer": cage.isomer or "",
                "unwrap_conflict": "false",
            }
        )
    data = pd.DataFrame(rows)
    data.to_csv(frame_dir / f"{result.frame.name}_membership.tsv", sep="\t", index=False)


def write_order_parameter(
    result: FrameResult,
    frame_dir: Path,
    order_parameters: Any | None = None,
) -> None:
    """Write per-water F3/F4/Q_l values for custom plotting or focus-water checks."""
    path = frame_dir / f"{result.frame.name}_order_parameter.tsv"
    selected = selected_order_parameters_for_result(result, order_parameters)
    selected_set = set(selected)
    q_degrees = q_degrees_from_order_parameters(selected)
    if result.f3f4 is None or not (selected_set & {"f3", "f4"} or q_degrees):
        path.unlink(missing_ok=True)
        return
    rows = []
    for item in result.f3f4.per_water:
        row = {
            "resid": item.resid,
            "atomid": item.atomid,
            "oxygen_index": item.oxygen,
            "x_nm": item.xyz[0],
            "y_nm": item.xyz[1],
            "z_nm": item.xyz[2],
        }
        if "f3" in selected_set:
            row["F3"] = item.f3
        if "f4" in selected_set:
            row["F4"] = item.f4
        for degree in q_degrees:
            row[f"q{degree}"] = item.q_values.get(degree)
        if q_degrees:
            row["q_neighbors"] = item.q_neighbors
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def cage_center_name(cage_type: str) -> str:
    """Return the short CNT atom name used for a cage center."""
    return {"512": "G512", "51262": "G62", "51263": "G63", "51264": "G64", "51268": "G68", "435663": "G436"}.get(cage_type, "CAGE")[:5]


def write_vmd_script(result: FrameResult, frame_dir: Path) -> None:
    """Write a small VMD helper script with default colors."""
    path = frame_dir / f"{result.frame.name}_view.vmd.tcl"
    lines = [
        "# SQQ VMD helper",
        "# Source this after loading SQQ GRO files.",
        "color Name R4 gray",
        "color Name R5 purple",
        "color Name R6 tan",
        "color Name R7 black",
        "color Name CP5 cyan",
        "color Name CP6 yellow",
        "color Name G512 blue",
        "color Name G62 green",
        "color Name G63 orange",
        "color Name G64 red",
        "display projection Orthographic",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
