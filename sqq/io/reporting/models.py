"""Data objects and stable schemas used by all report writers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True, slots=True)
class ReportTable:
    """One named table shared by CSV and XLSX writers."""

    name: str
    data: pd.DataFrame
    include_header: bool = True
    detail: bool = False


@dataclass(frozen=True, slots=True)
class ReportTables:
    """All materialized tables for one completed analysis run."""

    main: tuple[ReportTable, ...]
    detail: tuple[ReportTable, ...] = ()

    def by_name(self, name: str) -> ReportTable | None:
        return next((table for table in (*self.main, *self.detail) if table.name == name), None)


# Stable column order is part of the machine-readable output contract.
SUMMARY_COLUMNS = (
    "frame", "time_ps", "source", "status", "error", "n_atoms",
    "n_waters", "n_guests", "bond_mode", "n_edges", "connection_mode",
    "connection_count", "hbond_count", "oo_connection_count",
    "pair_connection_count", "mean_coordination", "coordination_0",
    "coordination_1", "coordination_2", "coordination_3", "coordination_4",
    "coordination_gt4", "coordination_0_fraction", "coordination_1_fraction",
    "coordination_2_fraction", "coordination_3_fraction",
    "coordination_4_fraction", "coordination_gt4_fraction",
    "degree_le2_fraction", "degree4_fraction", "over4_fraction", "ring4",
    "ring5", "ring6", "ring7", "free_ring4", "free_ring5", "free_ring6",
    "free_ring7", "half_cage_total", "half_cage_breakdown",
    "quasi_cage_total", "quasi_cage_breakdown", "cage_report_types",
    "cage_512", "cage_51262", "cage_51263", "cage_51264", "cage_51268",
    "cage_435663", "cage_total", "cage_empty", "cage_occupied",
    "hydrate_cluster_enabled", "hydrate_cluster_detail_enabled",
    "hydrate_cluster_count", "sI_cluster_count", "sII_cluster_count",
    "sH_cluster_count", "mixed_cluster_count", "unclassified_cluster_count",
    "hydrate_domain_count", "sI_domain_count", "sII_domain_count",
    "sH_domain_count", "classified_cage_count", "boundary_cage_count",
    "ambiguous_cage_count", "unclassified_cage_count", "isolated_cage_count",
    "largest_cluster_cage_count", "largest_cluster_water_count",
    "cluster_size_distribution", "MCG1_largest_cluster",
    "DHOP35_largest_cluster", "MCG3_largest_cluster", "DHOP30_largest_cluster",
    "F3_mean", "F4_mean", "q6_mean", "q12_mean", "F3_count", "F4_count",
    "q6_count", "q12_count", "F3_valid_waters", "F4_valid_waters",
    "q6_valid_waters", "q12_valid_waters", "F3_focus_mean", "F4_focus_mean",
    "q6_focus_mean", "q12_focus_mean", "F3_focus_count", "F4_focus_count",
    "q6_focus_count", "q12_focus_count", "F3_focus_valid_waters",
    "F4_focus_valid_waters", "q6_focus_valid_waters", "q12_focus_valid_waters",
    "ice_like_waters", "ice_i_waters", "interfacial_ice_waters",
)

SUMMARY_DETAIL_TABLE_NAMES = (
    "cage_occupancy",
    "cage_isomer",
    "quasi_cage_isomer",
    "hydrate_domain",
    "hydrate_cluster_detail",
)
SUMMARY_MAIN_TABLE_NAMES = (
    "summary", "failures", "connection", "hbond", "oo_connection",
    "pair_connection", "ring", "half_cage", "quasi_cage", "cage",
    "hydrate_cluster", "order_parameter", "ice", "detail_index", "config",
)

LEGACY_SUMMARY_DETAIL_DIRECTORY = "summary_detail"
QUASI_ISOMER_DETAIL_KEY = "_quasi_cage_isomer_detail"
EXCEL_MAX_ROWS = 1_048_576
EXCEL_MAX_COLUMNS = 16_384
FULL_TABLE_FORMAT_MAX_CELLS = 200_000
FULL_TABLE_FORMAT_MAX_COLUMNS = 128
LIGHTWEIGHT_TABLE_COLUMN_WIDTH = 16


def table_metric(
    sheet_name: str,
    table: pd.DataFrame,
    *,
    write_seconds: float = 0.0,
) -> dict[str, Any]:
    return {
        "sheet": sheet_name,
        "rows": int(len(table.index)),
        "columns": int(len(table.columns)),
        "cells": int(len(table.index) * len(table.columns)),
        "write_seconds": float(write_seconds),
    }
