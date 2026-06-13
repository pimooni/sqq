from __future__ import annotations

"""F3/F4 AOP order parameters from the user's reference script."""

from itertools import combinations
from math import atan2, cos, radians

import numpy as np

from ..models import Atom, F3F4Result, Frame, GraphResult, Water, WaterOrder
from .pbc import minimum_image


TETRAHEDRAL_COS2 = cos(radians(109.47)) ** 2


def compute_f3f4(frame: Frame, waters: list[Water], graph: GraphResult, focus_resids: set[int] | None = None) -> F3F4Result:
    """Compute per-water F3/F4 plus global and optional focus-water means."""
    water_by_oxygen = {water.oxygen: water for water in waters}
    focus_resids = set(focus_resids or set())
    rows: list[WaterOrder] = []
    f3_values: list[float] = []
    f4_values: list[float] = []
    f3_focus_values: list[float] = []
    f4_focus_values: list[float] = []

    for water in waters:
        neighbors = [idx for idx in sorted(graph.adjacency.get(water.oxygen, set())) if idx in water_by_oxygen]
        f3 = f3_for_water(frame.atoms, water.oxygen, neighbors, frame.box)
        f4 = f4_for_water(frame.atoms, water, [water_by_oxygen[idx] for idx in neighbors], frame.box)
        if f3 is not None:
            f3_values.append(f3)
            if water.resid in focus_resids:
                f3_focus_values.append(f3)
        if f4 is not None:
            f4_values.append(f4)
            if water.resid in focus_resids:
                f4_focus_values.append(f4)
        atom = frame.atoms[water.oxygen]
        rows.append(WaterOrder(oxygen=water.oxygen, resid=water.resid, atomid=atom.atomid, xyz=atom.xyz, f3=f3, f4=f4))

    return F3F4Result(
        per_water=tuple(rows),
        f3_mean=mean_or_none(f3_values),
        f4_mean=mean_or_none(f4_values),
        f3_valid=len(f3_values),
        f4_valid=len(f4_values),
        focus_resids=tuple(sorted(focus_resids)),
        f3_focus_mean=mean_or_none(f3_focus_values) if focus_resids else None,
        f4_focus_mean=mean_or_none(f4_focus_values) if focus_resids else None,
        f3_focus_valid=len(f3_focus_values),
        f4_focus_valid=len(f4_focus_values),
    )


def f3_for_water(atoms: list[Atom], oxygen: int, neighbors: list[int], box: np.ndarray | None) -> float | None:
    """F3: average squared deviation over neighbor-neighbor angles."""
    if len(neighbors) < 2:
        return None
    center = atoms[oxygen].xyz
    terms = []
    for a, b in combinations(neighbors, 2):
        va = minimum_image(atoms[a].xyz - center, box)
        vb = minimum_image(atoms[b].xyz - center, box)
        norm = float(np.linalg.norm(va) * np.linalg.norm(vb))
        if norm <= 1e-12:
            continue
        r = float(np.dot(va, vb) / norm)
        terms.append((r * abs(r) + TETRAHEDRAL_COS2) ** 2)
    return mean_or_none(terms)


def f4_for_water(atoms: list[Atom], water: Water, neighbors: list[Water], box: np.ndarray | None) -> float | None:
    """F4: average cos(3*dihedral) after selecting the farthest H-H pair."""
    if len(water.hydrogens) < 2:
        return None
    terms = []
    for neighbor in neighbors:
        if len(neighbor.hydrogens) < 2:
            continue
        h1, h2 = farthest_hydrogen_pair(atoms, water, neighbor, box)
        angle = dihedral_h_o_o_h(atoms, h1, water.oxygen, neighbor.oxygen, h2, box)
        terms.append(cos(3.0 * angle))
    return mean_or_none(terms)


def farthest_hydrogen_pair(atoms: list[Atom], water_a: Water, water_b: Water, box: np.ndarray | None) -> tuple[int, int]:
    """Match the shell script by choosing the largest inter-water H-H distance."""
    best_pair = (water_a.hydrogens[0], water_b.hydrogens[0])
    best_dist2 = -1.0
    for ha in water_a.hydrogens[:2]:
        for hb in water_b.hydrogens[:2]:
            delta = minimum_image(atoms[hb].xyz - atoms[ha].xyz, box)
            dist2 = float(np.dot(delta, delta))
            if dist2 > best_dist2:
                best_dist2 = dist2
                best_pair = (ha, hb)
    return best_pair


def dihedral_h_o_o_h(
    atoms: list[Atom],
    h1: int,
    o1: int,
    o2: int,
    h2: int,
    box: np.ndarray | None,
) -> float:
    """Return the locally unwrapped H-O-O-H dihedral angle in radians."""
    p2 = atoms[o1].xyz
    p1 = p2 + minimum_image(atoms[h1].xyz - atoms[o1].xyz, box)
    p3 = p2 + minimum_image(atoms[o2].xyz - atoms[o1].xyz, box)
    p4 = p3 + minimum_image(atoms[h2].xyz - atoms[o2].xyz, box)
    return dihedral(p1, p2, p3, p4)


def dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Compute a signed four-point dihedral angle in radians."""
    b0 = p1 - p2
    b1 = p3 - p2
    b2 = p4 - p3
    b1_norm = float(np.linalg.norm(b1))
    if b1_norm <= 1e-12:
        return 0.0
    b1 = b1 / b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return atan2(y, x)


def mean_or_none(values: list[float]) -> float | None:
    """Return a float mean while preserving empty values as None."""
    if not values:
        return None
    return float(np.mean(values))
