from __future__ import annotations

"""Build the water graph used by rings, open cage patches, cages, and order metrics."""

from collections import defaultdict
from itertools import combinations
from math import cos, radians
from pathlib import Path
import re

import numpy as np

from ..models import Atom, GraphResult, Water
from .pbc import minimum_image


def build_water_graph(
    atoms: list[Atom],
    waters: list[Water],
    box: np.ndarray | None,
    bond_mode: str,
    oo_cutoff_nm: float,
    hbond_distance_nm: float,
    hbond_angle_deg: float,
    pair_file: str | Path | None = None,
    pair_id: str = "resid",
) -> GraphResult:
    """Construct graph edges between water oxygen nodes."""
    mode = resolve_bond_mode(bond_mode, waters, pair_file)
    if mode == "pairs":
        edges = read_pair_edges(Path(pair_file), atoms, waters, pair_id)  # type: ignore[arg-type]
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
    """Resolve auto mode from available hydrogen atoms."""
    if bond_mode not in {"auto", "hbond", "oo", "pairs"}:
        raise ValueError(f"Unsupported bond_mode: {bond_mode}")
    if bond_mode == "auto":
        return "hbond" if any(water.hydrogens for water in waters) else "oo"
    if bond_mode == "pairs" and pair_file is None:
        raise ValueError("bond_mode=pairs requires graph.pair_file or --pairs.")
    return bond_mode


def adjacency_from_edges(waters: list[Water], edges: list[tuple[int, int]]) -> dict[int, set[int]]:
    """Build an adjacency dictionary over all selected water oxygens."""
    adjacency: dict[int, set[int]] = {water.oxygen: set() for water in waters}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def read_pair_edges(path: Path, atoms: list[Atom], waters: list[Water], pair_id: str) -> list[tuple[int, int]]:
    """Read a user-provided water-neighbor pair file."""
    id_map = water_id_map(atoms, waters, pair_id)
    edges: set[tuple[int, int]] = set()
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = [part for part in re.split(r"[\s,;]+", line) if part]
        if len(parts) < 2:
            raise ValueError(f"Invalid pair line {lineno} in {path}: {raw_line!r}")
        try:
            left = int(parts[0])
            right = int(parts[1])
        except ValueError as exc:
            raise ValueError(f"Pair ids must be integers at {path}:{lineno}") from exc
        if left not in id_map or right not in id_map:
            raise ValueError(f"Pair id not found at {path}:{lineno}: {left}, {right}")
        a, b = sorted((id_map[left], id_map[right]))
        if a != b:
            edges.add((a, b))
    return sorted(edges)


def water_id_map(atoms: list[Atom], waters: list[Water], pair_id: str) -> dict[int, int]:
    """Map unique external pair ids to internal oxygen node indices."""
    if pair_id == "resid":
        identifiers = ((water.resid, water.oxygen) for water in waters)
    elif pair_id == "oxygen_index":
        identifiers = ((water.oxygen, water.oxygen) for water in waters)
    elif pair_id == "atomid":
        identifiers = ((atoms[water.oxygen].atomid, water.oxygen) for water in waters)
    else:
        raise ValueError("graph.pair_id must be one of: resid, oxygen_index, atomid")

    id_map: dict[int, int] = {}
    for identifier, oxygen in identifiers:
        if identifier in id_map and id_map[identifier] != oxygen:
            raise ValueError(
                f"graph.pair_id={pair_id!r} is not unique among selected waters; "
                "use oxygen_index or atomid."
            )
        id_map[identifier] = oxygen
    return id_map


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
        key = tuple(np.floor(coord / lengths * n_cells).astype(int))
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
