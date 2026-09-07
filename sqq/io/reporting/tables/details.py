"""Build detailed scientific summary and auxiliary report tables."""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from ....config import DEFAULT_MODE, is_cpp_mode, normalize_order_parameters, output_enabled
from .formatting import (
    cage_display_label,
    count_cell,
    ordered_cage_types,
    patch_composition_label,
    patch_display_label,
    patch_summary_label,
)
from .frame import QUASI_ISOMER_DETAIL_KEY


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

__all__ = (
    "summary_detail_tables",
    "hydrate_cluster_summary_table",
    "hydrate_cluster_is_enabled",
    "hydrate_cluster_detail_table",
    "hydrate_domain_table",
    "hydrate_motif_table",
    "expanded_hydrate_table",
    "connection_sheet_name",
    "connection_summary_table",
    "summary_simple_table",
    "order_parameter_summary_table",
    "infer_order_parameters_from_data",
    "q_degree_from_data",
    "molecule_summary_table",
    "patch_summary_table",
    "summary_cage_types_from_data",
    "cage_summary_table",
    "cage_occupancy_summary_table",
    "global_guest_labels",
    "row_guest_labels",
    "cage_occupancy_label_from_column",
    "cage_count",
    "cage_isomer_summary_table",
    "global_cage_isomer_labels",
    "quasi_cage_isomer_summary_table",
)
