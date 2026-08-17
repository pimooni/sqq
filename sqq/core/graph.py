"""Build the water graph used by rings, open cage patches, cages, and order metrics."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations
from math import cos, radians
from pathlib import Path

import numpy as np

from ..models import Atom, GraphResult, Water
from ..config.resolve import resolve_graph_mode
from .pbc import minimum_image


def build_water_graph(
    atoms: list[Atom],
    waters: list[Water],
    box: np.ndarray | None,
    bond_mode: str,
    oo_cutoff_nm: float,
    hbond_distance_nm: float,
    hbond_angle_deg: float,
    pair_edges: Iterable[tuple[int, int]] | None = None,
) -> GraphResult:
    """Construct graph edges between water oxygen nodes."""
    if str(bond_mode).strip().lower() == "pairs":
        if pair_edges is None:
            raise ValueError("bond_mode=pairs requires normalized oxygen pair edges.")
        mode = "pairs"
    else:
        mode = resolve_bond_mode(bond_mode, waters)
    if mode == "pairs":
        edges = _canonical_oxygen_edges(pair_edges, waters)
        adjacency = adjacency_from_edges(waters, edges)
        return GraphResult(mode=mode, edges=edges, adjacency=adjacency)

    cutoff = hbond_distance_nm if mode == "hbond" else oo_cutoff_nm
    edges: list[tuple[int, int]] = []

    for wa, wb in iter_water_pairs(atoms, waters, box, cutoff):
        oa = atoms[wa.oxygen].xyz
        ob = atoms[wb.oxygen].xyz
        oo = minimum_image(ob - oa, box)
        dist = float(np.linalg.norm(oo))
        if dist > cutoff:
            continue
        if mode == "hbond" and not hbond_angle_ok(atoms, wa, wb, oo, dist, box, hbond_angle_deg):
            continue
        a, b = sorted((wa.oxygen, wb.oxygen))
        edges.append((a, b))

    edges.sort()
    adjacency = adjacency_from_edges(waters, edges)
    return GraphResult(mode=mode, edges=edges, adjacency=adjacency)


def resolve_bond_mode(bond_mode: str, waters: list[Water], pair_file: str | Path | None = None) -> str:
    """Return the effective graph mode under the shared resolution rule."""
    return resolve_graph_mode(bond_mode, waters, pair_file).effective


def adjacency_from_edges(waters: list[Water], edges: list[tuple[int, int]]) -> dict[int, set[int]]:
    """Build an adjacency dictionary over all selected water oxygens."""
    adjacency: dict[int, set[int]] = {water.oxygen: set() for water in waters}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def _canonical_oxygen_edges(
    pair_edges: Iterable[tuple[int, int]], waters: list[Water]
) -> list[tuple[int, int]]:
    """Enforce the scientific-kernel contract for pair graph endpoints."""
    oxygen_indices = {int(water.oxygen) for water in waters}
    normalized: set[tuple[int, int]] = set()
    for edge_number, edge in enumerate(pair_edges, start=1):
        try:
            left, right = edge
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"pair_edges edge {edge_number} must contain exactly two oxygen atom indices."
            ) from exc
        left, right = int(left), int(right)
        unknown = [value for value in (left, right) if value not in oxygen_indices]
        if unknown:
            values = ", ".join(str(value) for value in unknown)
            raise ValueError(
                "pair_edges must contain selected-water oxygen atom indices only; "
                f"edge {edge_number} contains {values}."
            )
        if left != right:
            normalized.add(tuple(sorted((left, right))))
    return sorted(normalized)


def iter_water_pairs(atoms: list[Atom], waters: list[Water], box: np.ndarray | None, cutoff: float):
    """Yield cutoff candidates using MDAnalysis, with a deterministic cell-list fallback."""
    accelerated = mdanalysis_water_pairs(atoms, waters, box, cutoff)
    if accelerated is not None:
        yield from accelerated
        return
    yield from cell_list_water_pairs(atoms, waters, box, cutoff)


def mdanalysis_water_pairs(
    atoms: list[Atom],
    waters: list[Water],
    box: np.ndarray | None,
    cutoff: float,
) -> list[tuple[Water, Water]] | None:
    """Use the installed MDAnalysis neighbor search for an orthorhombic frame."""
    try:
        from MDAnalysis.lib.distances import self_capped_distance
    except ImportError:
        return None
    if not waters:
        return []
    coordinates = np.asarray([atoms[water.oxygen].xyz for water in waters], dtype=float)
    dimensions = None
    if box is not None and len(box) >= 3 and np.all(np.asarray(box[:3], dtype=float) > 0):
        lengths = np.asarray(box[:3], dtype=float)
        dimensions = np.asarray([*lengths, 90.0, 90.0, 90.0], dtype=float)
    try:
        pairs = self_capped_distance(
            coordinates,
            max_cutoff=float(cutoff) + 1.0e-7,
            box=dimensions,
            return_distances=False,
        )
    except (RuntimeError, ValueError, TypeError):
        return None
    normalized = sorted({tuple(sorted((int(left), int(right)))) for left, right in np.asarray(pairs) if left != right})
    return [(waters[left], waters[right]) for left, right in normalized]


def cell_list_water_pairs(atoms: list[Atom], waters: list[Water], box: np.ndarray | None, cutoff: float):
    """Yield water pairs using the established orthorhombic cell list."""
    if box is None or len(box) < 3 or np.any(np.asarray(box[:3]) <= 0):
        yield from combinations(waters, 2)
        return

    lengths = np.asarray(box[:3], dtype=float)
    n_cells = np.maximum(np.floor(lengths / cutoff).astype(int), 1)
    cells: dict[tuple[int, int, int], list[Water]] = defaultdict(list)
    for water in waters:
        coord = atoms[water.oxygen].xyz % lengths
        cell = np.floor(coord / lengths * n_cells).astype(int)
        cell = np.minimum(np.maximum(cell, 0), n_cells - 1)
        key = tuple(cell)
        cells[key].append(water)

    visited: set[tuple[int, int]] = set()
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for key, local in cells.items():
        for dx, dy, dz in offsets:
            other_key = ((key[0] + dx) % n_cells[0], (key[1] + dy) % n_cells[1], (key[2] + dz) % n_cells[2])
            for wa in local:
                for wb in cells.get(other_key, []):
                    if wa.oxygen == wb.oxygen:
                        continue
                    pair = tuple(sorted((wa.oxygen, wb.oxygen)))
                    if pair in visited:
                        continue
                    visited.add(pair)
                    yield wa, wb

def hbond_angle_ok(
    atoms: list[Atom],
    wa: Water,
    wb: Water,
    oo_vec: np.ndarray,
    oo_dist: float,
    box: np.ndarray | None,
    angle_deg: float,
) -> bool:
    """Accept a hydrogen bond if either water donates toward the other oxygen."""
    cos_limit = cos(radians(angle_deg))
    return donor_angle_ok(atoms, wa, oo_vec, oo_dist, box, cos_limit) or donor_angle_ok(
        atoms, wb, -oo_vec, oo_dist, box, cos_limit
    )


def donor_angle_ok(
    atoms: list[Atom],
    donor: Water,
    oo_vec: np.ndarray,
    oo_dist: float,
    box: np.ndarray | None,
    cos_limit: float,
) -> bool:
    """Check the donor O-H direction against the donor-acceptor O-O vector."""
    origin = atoms[donor.oxygen].xyz
    for hidx in donor.hydrogens[:2]:
        oh = minimum_image(atoms[hidx].xyz - origin, box)
        oh_norm = float(np.linalg.norm(oh))
        if oh_norm <= 1e-12:
            continue
        if float(np.dot(oo_vec, oh)) >= oo_dist * oh_norm * cos_limit:
            return True
    return False
