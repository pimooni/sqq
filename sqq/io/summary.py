from __future__ import annotations

"""Markdown, TSV, VMD, and XLSX summary writers."""

from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

from ..banner import SQQ_BANNER
from ..config import dump_config
from ..core.cage import KNOWN_CAGE_TYPES, parse_cage_face_label
from ..models import CagePatch, FrameResult
from .occupancy import guest_composition_label, guest_lookup as build_guest_lookup, guest_resname_order as guest_resname_order_from_guests


# Column order is stable so downstream plotting scripts can rely on it.
SUMMARY_COLUMNS = [
    "frame",
    "time_ps",
    "source",
    "status",
    "error",
    "n_atoms",
    "n_waters",
    "n_guests",
    "bond_mode",
    "n_edges",
    "connection_mode",
    "connection_count",
    "hbond_count",
    "oo_connection_count",
    "pair_connection_count",
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
    "ring4",
    "ring5",
    "ring6",
    "ring7",
    "free_ring4",
    "free_ring5",
    "free_ring6",
    "free_ring7",
    "half_cage_total",
    "half_cage_breakdown",
    "quasi_cage_total",
    "quasi_cage_breakdown",
    "cage_report_types",
    "cage_512",
    "cage_51262",
    "cage_51263",
    "cage_51264",
    "cage_51268",
    "cage_435663",
    "cage_total",
    "cage_empty",
    "cage_occupied",
    "F3_mean",
    "F4_mean",
    "F3_valid_waters",
    "F4_valid_waters",
    "F3_focus_mean",
    "F4_focus_mean",
    "F3_focus_valid_waters",
    "F4_focus_valid_waters",
    "ice_like_waters",
    "ice_i_waters",
    "interfacial_ice_waters",
]

TREE_MIDDLE = "\u251c"
TREE_LAST = "\u2514"
COLUMN_CHILD = "\u21b3"
SUBSCRIPT_DIGIT_DELETE = dict.fromkeys(range(0x2080, 0x208A))


def result_row(result: FrameResult) -> dict[str, Any]:
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
        "F3_mean": None if f3f4 is None else f3f4.f3_mean,
        "F4_mean": None if f3f4 is None else f3f4.f4_mean,
        "F3_valid_waters": None if f3f4 is None else f3f4.f3_valid,
        "F4_valid_waters": None if f3f4 is None else f3f4.f4_valid,
        "F3_focus_mean": None if f3f4 is None else f3f4.f3_focus_mean,
        "F4_focus_mean": None if f3f4 is None else f3f4.f4_focus_mean,
        "F3_focus_valid_waters": None if f3f4 is None else f3f4.f3_focus_valid,
        "F4_focus_valid_waters": None if f3f4 is None else f3f4.f4_focus_valid,
        "ice_like_waters": len(result.ice_like_waters),
        "ice_i_waters": len(result.ice_i_waters),
        "interfacial_ice_waters": len(result.interfacial_ice_waters),
    }
    for resname, count in molecule_counts.items():
        row[f"mol_{resname}"] = count
    row["mol_TOTAL"] = len(result.frame.atoms)
    row["guest_order"] = ";".join(guest_order)
    for patch_type, count in half_cage_counts.items():
        row[f"half_cage_{patch_type}"] = count
    for patch_type, count in quasi_cage_counts.items():
        row[f"quasi_cage_{patch_type}"] = count

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


def patch_counts(patches) -> dict[str, int]:
    """Count open cage patches by patch_type."""
    counts: dict[str, int] = {}
    for patch in patches:
        counts[patch.patch_type] = counts.get(patch.patch_type, 0) + 1
    return counts


def patch_breakdown(counts: dict[str, int]) -> str:
    """Render a compact patch count list for broad summary rows."""
    return "; ".join(f"{key}:{counts[key]}" for key in sorted(counts))


def failed_row(frame_name: str, source: str, error: str) -> dict[str, Any]:
    """Create a summary row for a skipped or failed frame."""
    row = {column: "" for column in SUMMARY_COLUMNS}
    row.update({"frame": frame_name, "source": source, "status": "failed", "error": error})
    return row


def write_frame_info(result: FrameResult, frame_dir: Path, ring_sizes: list[int] | None = None) -> None:
    """Write the per-frame Markdown report with inspection-oriented tables."""
    frame_dir.mkdir(parents=True, exist_ok=True)
    row = result_row(result)
    cage_values = {cage.cage_type for cage in result.cages}
    cage_types = [cage_type for cage_type in ordered_cage_types(cage_values) if cage_type in cage_values]
    default_ring_sizes = result.ring_report_sizes or tuple(result.rings)
    enabled_ring_sizes = sorted(set(ring_sizes if ring_sizes is not None else default_ring_sizes))
    lines = [
        f"# SQQ Frame Report: {result.frame.name}",
        "",
    ]

    lines.extend(
        section_table(
            "Frame Information",
            ["item", "value"],
            [
                ["frame", result.frame.name],
                ["time_ps", result.frame.time_ps],
                ["source", source_label(result.frame.source)],
                ["bond_mode", row["connection_mode"]],
                ["ring_sizes", ", ".join(str(size) for size in enabled_ring_sizes)],
                ["status", "ok"],
                ["n_atoms", len(result.frame.atoms)],
                ["n_waters", len(result.waters)],
                ["n_guests", len(result.guests)],
            ],
        )
    )
    lines.extend(section_table("Molecules", ["resname", "molecules", "atoms"], molecule_count_rows(result)))
    lines.extend(connection_info_section(result, row))

    ring_rows = [[size, row.get(f"free_ring{size}", 0)] for size in enabled_ring_sizes]
    ring_rows.append(["total", sum(int(item[1]) for item in ring_rows)])
    lines.extend(section_table("Ring", ["ring type", "count"], ring_rows))
    lines.extend(patch_info_section("Half Cage", result.half_cages))
    lines.extend(patch_info_section("Quasi Cage", result.quasi_cages))

    cage_rows = [[cage_display_label(cage_type), row.get(f"cage_{cage_type}", 0)] for cage_type in cage_types]
    cage_rows.append(["total", len(result.cages)])
    lines.extend(section_table("Cage", ["cage type", "count"], cage_rows))
    lines.extend(cage_occupancy_section(result, cage_types))
    lines.extend(cage_isomer_section(result, cage_types))

    f3f4_headers = ["metric", "valid_waters", "mean"]
    f3f4_rows = [
        ["F3", row["F3_valid_waters"], row["F3_mean"]],
        ["F4", row["F4_valid_waters"], row["F4_mean"]],
    ]
    if result.f3f4 is not None and result.f3f4.focus_resids:
        f3f4_headers.extend(["focus_valid_waters", "focus_mean"])
        f3f4_rows[0].extend([row["F3_focus_valid_waters"], row["F3_focus_mean"]])
        f3f4_rows[1].extend([row["F4_focus_valid_waters"], row["F4_focus_mean"]])
    lines.extend(section_table("F3/F4", f3f4_headers, f3f4_rows))

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


def section_table(title: str, headers: list[str], rows: list[list[Any]]) -> list[str]:
    """Render one small markdown section."""
    if not rows:
        return []
    lines = ["", f"## {title}", "", markdown_rows(headers, rows).rstrip()]
    return lines


def markdown_rows(headers: list[str], rows: list[list[Any]]) -> str:
    """Render rows as a compact markdown table."""
    text_rows = [[format_summary_cell(value) for value in row] for row in rows]
    widths = [len(str(header)) for header in headers]
    for row in text_rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    lines = [
        "| " + " | ".join(str(header).ljust(widths[idx]) for idx, header in enumerate(headers)) + " |",
        "| " + " | ".join("-" * widths[idx] for idx in range(len(headers))) + " |",
    ]
    for row in text_rows:
        lines.append("| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + " |")
    return "\n".join(lines) + "\n"


def molecule_count_rows(result: FrameResult) -> list[list[Any]]:
    """Count contiguous molecules and atoms by residue name in source order."""
    atom_counts = atom_resname_counts(result)
    molecule_counts: dict[str, int] = {}
    previous_residue: tuple[int, str] | None = None
    for atom in result.frame.atoms:
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


def atom_resname_counts(result: FrameResult) -> dict[str, int]:
    """Count atoms by residue name while preserving source-file order."""
    counts: dict[str, int] = {}
    for atom in result.frame.atoms:
        counts[atom.resname] = counts.get(atom.resname, 0) + 1
    return counts


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


def patch_display_label(patch_type: str) -> str:
    """Remove internal HC/QC prefixes from human-facing patch labels."""
    for prefix in ("hc_", "qc_"):
        if patch_type.startswith(prefix):
            return patch_type.removeprefix(prefix)
    return patch_type


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


def cage_occupancy_section(result: FrameResult, cage_types: list[str]) -> list[str]:
    """Show one cage type per row with dynamic guest-composition columns."""
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
    headers = [
        "cage type",
        "total",
        "empty",
        "occupied",
        *[f"{COLUMN_CHILD} {label}" for label in child_labels],
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
                *[f"{COLUMN_CHILD} {cage_counts.get(label, 0)}" for label in child_labels],
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
                f"{COLUMN_CHILD} {sum(counts.get(cage_type, {}).get(label, 0) for cage_type in cage_types)}"
                for label in child_labels
            ],
        ]
    )
    return section_table("Cage Occupancy", headers, table_rows)


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


def cage_isomer_section(result: FrameResult, cage_types: list[str]) -> list[str]:
    """Show cage composition totals with nested structural isomers."""
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
    if not rows:
        rows.append(["total", 0])
    return section_table("Cage Isomer", ["cage isomer", "count"], rows)


def source_label(source: Path | None) -> str:
    """Return an absolute source path for human-facing reports."""
    if source is None:
        return ""
    return str(Path(source).resolve())


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


def write_f3f4(result: FrameResult, frame_dir: Path) -> None:
    """Write per-water F3/F4 values for custom plotting or focus-water checks."""
    if result.f3f4 is None:
        return
    rows = [
        {
            "resid": item.resid,
            "atomid": item.atomid,
            "oxygen_index": item.oxygen,
            "x_nm": item.xyz[0],
            "y_nm": item.xyz[1],
            "z_nm": item.xyz[2],
            "F3": item.f3,
            "F4": item.f4,
        }
        for item in result.f3f4.per_water
    ]
    pd.DataFrame(rows).to_csv(frame_dir / f"{result.frame.name}_f3f4.tsv", sep="\t", index=False)


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


def write_summary(
    rows: list[dict[str, Any]],
    outdir: Path,
    config: dict[str, Any],
    write_xlsx: bool = True,
    run_info: dict[str, Any] | None = None,
) -> None:
    """Write global XLSX summaries and the final run config."""
    outdir.mkdir(parents=True, exist_ok=True)
    columns = list(SUMMARY_COLUMNS)
    extra_columns = stable_extra_columns(rows, columns)
    data = pd.DataFrame(rows, columns=columns + extra_columns)
    summary_md = outdir / "summary.md"
    if summary_md.exists():
        summary_md.unlink()
    with (outdir / "run_config.yaml").open("w", encoding="utf-8", newline="\n") as handle:
        dump_config(config, handle)
    if write_xlsx:
        with pd.ExcelWriter(outdir / "summary.xlsx", engine="openpyxl") as writer:
            summary_dashboard_table(data, run_info or {}, config).to_excel(writer, sheet_name="summary", index=False, header=False)
            for sheet_name, table in summary_sheet_tables(data, config).items():
                table.to_excel(writer, sheet_name=sheet_name, index=False)
            pd.DataFrame(flatten_config(config)).to_excel(writer, sheet_name="config", index=False)
            format_summary_workbook(writer.book)


def summary_dashboard_table(data: pd.DataFrame, run_info: dict[str, Any], config: dict[str, Any]) -> pd.DataFrame:
    """Build a compact human-facing dashboard for the summary sheet."""
    banner_lines = [line.strip("| ") for line in SQQ_BANNER.splitlines() if line.startswith("|")]
    title = banner_lines[0] if banner_lines else "Shell  Quant  Qualifier"
    author = banner_lines[1] if len(banner_lines) > 1 else "by J. PANG & Q. SUN"
    matched_files = run_info.get("matched_files", "")
    try:
        matched_count = int(matched_files)
    except (TypeError, ValueError):
        matched_count = 0

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
        ["run_config.yaml", run_info.get("run_config", "")],
        ["", ""],
        ["Configuration", ""],
        ["Config", run_info.get("config_file", "<built-in defaults>")],
        ["Topology", run_info.get("topology", "<none>")],
        ["Mode", run_info.get("mode", "")],
        ["Graph mode", first_data_value(data, "connection_mode", run_info.get("graph_mode", ""))],
        ["Search sizes", excel_scalar(config.get("ring", {}).get("sizes", ""))],
        ["Ring report sizes", excel_scalar(configured_ring_report_sizes(config))],
        ["Quasi-cage sizes", f"{excel_scalar(config.get('quasi_cage', {}).get('base_sizes', 'auto'))} / {excel_scalar(config.get('quasi_cage', {}).get('side_sizes', 'auto'))}"],
        ["Quasi max layers", config.get("quasi_cage", {}).get("max_layers", "")],
        ["Cage report types", dashboard_cage_targets(config)],
        ["Maximum cage faces", config.get("cage", {}).get("max_faces", 20)],
        ["Output layout", run_info.get("output_layout", "")],
        ["Worker policy", run_info.get("worker_policy", "")],
        ["Workers", run_info.get("workers", "")],
        ["", ""],
        ["Analysis Results", ""],
        ["Frames total / ok / failed", f"{len(data)} / {frames_ok_count(data)} / {frames_failed_count(data)}"],
        ["Water molecules", sum_numeric_column(data, "n_waters")],
        ["Guest molecules", sum_numeric_column(data, "n_guests")],
        ["Connections", sum_numeric_column(data, "connection_count")],
    ])
    for size in configured_ring_report_sizes(config):
        rows.append([f"Ring{size}", sum_numeric_column(data, f"ring{size}")])
    rows.extend([
        ["Half cage", sum_numeric_column(data, "half_cage_total")],
        ["Quasi cage", sum_numeric_column(data, "quasi_cage_total")],
        ["Cage total", sum_numeric_column(data, "cage_total")],
        ["Empty cage", sum_numeric_column(data, "cage_empty")],
        ["Occupied cage", sum_numeric_column(data, "cage_occupied")],
        ["Ice-like waters", sum_numeric_column(data, "ice_like_waters")],
    ])
    return pd.DataFrame(rows)


def frames_ok_count(data: pd.DataFrame) -> int:
    """Count successfully analyzed frames."""
    return int((data.get("status") == "ok").sum()) if "status" in data else len(data)


def frames_failed_count(data: pd.DataFrame) -> int:
    """Count failed frames."""
    return int((data.get("status") == "failed").sum()) if "status" in data else 0


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
        if targets.strip().lower() == "all":
            return "all detected cages"
        raw_targets = [item.strip() for item in targets.split(",") if item.strip()]
    else:
        raw_targets = [str(item) for item in targets or []]
    return ", ".join(cage_display_label(target) for target in raw_targets)

def format_summary_workbook(workbook) -> None:
    """Apply readable formatting to the generated summary workbook."""
    for worksheet in workbook.worksheets:
        worksheet.sheet_view.showGridLines = False
        if worksheet.title == "summary":
            format_summary_dashboard_sheet(worksheet)
        else:
            format_table_sheet(worksheet)


def format_summary_dashboard_sheet(worksheet) -> None:
    """Style the human-facing dashboard sheet."""
    title_fill = PatternFill("solid", fgColor="F7E7C6")
    author_fill = PatternFill("solid", fgColor="FFF7E6")
    section_fill = PatternFill("solid", fgColor="2563EB")
    label_fill = PatternFill("solid", fgColor="EFF6FF")
    thin_border = Border(
        left=Side(style="thin", color="CBD5E1"),
        right=Side(style="thin", color="CBD5E1"),
        top=Side(style="thin", color="CBD5E1"),
        bottom=Side(style="thin", color="CBD5E1"),
    )
    widths = {"A": 28, "B": 120}
    for column, width in widths.items():
        worksheet.column_dimensions[column].width = width
    worksheet.row_dimensions[1].height = 30
    worksheet.row_dimensions[2].height = 22
    worksheet.freeze_panes = "A4"

    worksheet.merge_cells("A1:B1")
    worksheet.merge_cells("A2:B2")
    for row in (1, 2):
        for cell in worksheet[row]:
            cell.fill = title_fill if row == 1 else author_fill
            cell.font = Font(color="3B2A14", bold=True, size=18 if row == 1 else 11)
            cell.alignment = Alignment(horizontal="center", vertical="center")

    section_labels = {"Basic Information", "Configuration", "Analysis Results"}
    for row in worksheet.iter_rows():
        if row[0].row <= 2:
            continue
        for cell in row:
            cell.alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
            if cell.value not in (None, ""):
                cell.border = thin_border
        if row[0].value in section_labels:
            for cell in row:
                cell.fill = section_fill
                cell.font = Font(color="FFFFFF", bold=True)
                cell.alignment = Alignment(horizontal="left", vertical="center")
        else:
            if row[0].value not in (None, ""):
                row[0].fill = label_fill
                row[0].font = Font(bold=True, color="0F172A")
            if row[1].value not in (None, ""):
                row[1].font = Font(color="111827")
                row[1].alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)


def format_table_sheet(worksheet) -> None:
    """Style plotting-friendly data sheets without changing their data shape."""
    if worksheet.max_row < 1 or worksheet.max_column < 1:
        return
    header_fill = PatternFill("solid", fgColor="1E3A8A")
    header_font = Font(color="FFFFFF", bold=True)
    thin_border = Border(bottom=Side(style="thin", color="CBD5E1"))
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions
    for cell in worksheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin_border
    for column_index in range(1, worksheet.max_column + 1):
        letter = get_column_letter(column_index)
        width = estimated_column_width(worksheet, column_index)
        worksheet.column_dimensions[letter].width = width
    for row in worksheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)


def estimated_column_width(worksheet, column_index: int) -> int:
    """Estimate a bounded Excel column width from visible cell values."""
    max_length = 8
    for row_index in range(1, min(worksheet.max_row, 200) + 1):
        value = worksheet.cell(row=row_index, column=column_index).value
        if value is None:
            continue
        max_length = max(max_length, len(str(value)))
    return min(max(max_length + 2, 10), 48)


def stable_extra_columns(rows: list[dict[str, Any]], base_columns: list[str]) -> list[str]:
    """Collect non-core columns in first-seen order."""
    seen = set(base_columns)
    extras: list[str] = []
    for row in rows:
        for key in row:
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
            ("F3/F4", summary_simple_table(data, ["frame", "time_ps", "F3_mean", "F4_mean", "F3_valid_waters", "F4_valid_waters", "F3_focus_mean", "F4_focus_mean", "F3_focus_valid_waters", "F4_focus_valid_waters"])),
            ("Ice", summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"])),
        ]
    )
    return tables


def summary_sheet_tables(data: pd.DataFrame, config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """Build XLSX sheets using the configured report scopes."""
    ring_sizes = configured_ring_report_sizes(config)
    ring_columns = ["frame", "time_ps"]
    for size in ring_sizes:
        ring_columns.extend([f"ring{size}", f"free_ring{size}"])
    tables: dict[str, pd.DataFrame] = {
        "frame": frame_summary_table(data, ring_sizes),
        connection_sheet_name(data): connection_summary_table(data),
        "ring": summary_simple_table(data, ring_columns),
        "half_cage": patch_summary_table(data, "half_cage"),
        "quasi_cage": patch_summary_table(data, "quasi_cage"),
        "cage": cage_summary_table(data),
        "cage_occupancy": cage_occupancy_summary_table(data, markdown_style=False),
        "cage_isomer": cage_isomer_summary_table(data, include_zero_rows=True),
        "f3f4": summary_simple_table(data, ["frame", "time_ps", "F3_mean", "F4_mean", "F3_valid_waters", "F4_valid_waters", "F3_focus_mean", "F4_focus_mean", "F3_focus_valid_waters", "F4_focus_valid_waters"]),
        "ice": summary_simple_table(data, ["frame", "time_ps", "ice_like_waters", "ice_i_waters", "interfacial_ice_waters"]),
    }
    return {name: table for name, table in tables.items() if not table.empty}


def frame_summary_table(data: pd.DataFrame, ring_sizes: list[int]) -> pd.DataFrame:
    """Build the main per-frame table using the active report scopes."""
    columns = [
        "frame",
        "time_ps",
        "source",
        "status",
        "error",
        "n_atoms",
        "n_waters",
        "n_guests",
    ]
    for size in ring_sizes:
        columns.extend([f"ring{size}", f"free_ring{size}"])
    columns.extend(["half_cage_total", "quasi_cage_total"])
    columns.extend(f"cage_{cage_type}" for cage_type in summary_cage_types_from_data(data))
    columns.extend([
        "cage_total",
        "cage_empty",
        "cage_occupied",
        "F3_mean",
        "F4_mean",
        "F3_valid_waters",
        "F4_valid_waters",
        "ice_like_waters",
        "ice_i_waters",
        "interfacial_ice_waters",
    ])
    return summary_simple_table(data, columns)


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
    labels = [column.removeprefix(f"{prefix}_") for column in sorted_columns]
    rows: list[dict[str, Any]] = []
    for _, record in data.iterrows():
        row: dict[str, Any] = {
            "frame": record.get("frame", ""),
            "time_ps": record.get("time_ps", ""),
        }
        total = 0
        for column, label in zip(sorted_columns, labels):
            count = count_cell(record.get(column, 0))
            row[label] = count
            total += count
        row["total"] = total
        rows.append(row)
    return pd.DataFrame(rows)


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
    """Flatten nested config keys for the XLSX config sheet."""
    rows = []
    for key, value in config.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            rows.extend(flatten_config(value, name))
        else:
            rows.append({"parameter": name, "value": repr(value)})
    return rows









