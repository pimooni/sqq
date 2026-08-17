"""Order-parameter result contracts."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class WaterOrder:
    oxygen: int
    resid: int
    atomid: int
    xyz: np.ndarray
    f3: float | None
    f4: float | None
    q_values: dict[int, float | None] = field(default_factory=dict)
    q_neighbors: int = 0


@dataclass(frozen=True)
class F3F4Result:
    per_water: tuple[WaterOrder, ...]
    f3_mean: float | None
    f4_mean: float | None
    f3_valid: int
    f4_valid: int
    focus_resids: tuple[int, ...] = ()
    f3_focus_mean: float | None = None
    f4_focus_mean: float | None = None
    f3_focus_valid: int = 0
    f4_focus_valid: int = 0
    q_degree: tuple[int, ...] = ()
    q_means: dict[int, float | None] = field(default_factory=dict)
    q_valid_counts: dict[int, int] = field(default_factory=dict)
    q_focus_means: dict[int, float | None] = field(default_factory=dict)
    q_focus_valid_counts: dict[int, int] = field(default_factory=dict)
    q_neighbor_mode: str = ""
    q_cutoff_nm: float | None = None
    q_n_neighbor: int | None = None


@dataclass(frozen=True)
class ClusterOrderValue:
    """Largest connected component for one hydrate order parameter."""

    largest_cluster_size: int | None
    members: tuple[int, ...] = ()
    eligible_count: int = 0
    member_type: str = ""


@dataclass(frozen=True)
class HydrateOrderResult:
    """MCG and DHOP values reported for one frame."""

    mcg1: ClusterOrderValue
    dhop35: ClusterOrderValue
    mcg3: ClusterOrderValue | None = None
    dhop30: ClusterOrderValue | None = None


__all__ = ["WaterOrder", "F3F4Result", "ClusterOrderValue", "HydrateOrderResult"]
