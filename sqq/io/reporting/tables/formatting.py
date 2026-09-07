"""Shared labels, ordering, scalar formatting, and table utilities."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from ....models import FrameResult
from ....models.cage_type import KNOWN_CAGE_TYPES, parse_cage_face_label
from ....presentation.occupancy import (
    guest_resname_order as guest_resname_order_from_guests,
)

SUBSCRIPT_DIGIT_DELETE = dict.fromkeys(range(0x2080, 0x208A))


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


def patch_summary_label(patch_type: str, prefix: str) -> str:
    """Return the workbook label for an open-patch type."""
    if prefix != "quasi_cage":
        return patch_type
    return patch_composition_label(patch_type)


def patch_composition_label(patch_type: str) -> str:
    """Return a patch label without internal prefix or ring-sequence isomer marks."""
    return patch_display_label(patch_type).translate(SUBSCRIPT_DIGIT_DELETE)


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

__all__ = (
    "cluster_guest_composition",
    "cluster_size_distribution",
    "patch_counts",
    "patch_breakdown",
    "atom_resname_counts",
    "superscript_number",
    "ordered_cage_types",
    "present_cage_types",
    "cage_sort_key",
    "cage_display_label",
    "format_cage_type_counts",
    "guest_composition_sort_key",
    "guest_resname_order",
    "patch_display_label",
    "source_label",
    "frames_ok_count",
    "frames_failed_count",
    "has_selected_guests",
    "first_data_value",
    "sum_numeric_column",
    "min_mean_max_column",
    "format_stat_value",
    "molecule_totals",
    "configured_ring_report_sizes",
    "dashboard_cage_targets",
    "stable_extra_columns",
    "patch_summary_label",
    "patch_composition_label",
    "count_cell",
    "markdown_table",
    "format_summary_cell",
    "excel_scalar",
    "flatten_config",
)
