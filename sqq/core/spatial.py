from __future__ import annotations

"""Deterministic cutoff-pair searches for orthorhombic coordinate sets."""

from itertools import product
from math import ceil

import numpy as np

from .pbc import minimum_image


NEIGHBOR_OFFSETS = tuple(product((-1, 0, 1), repeat=3))


class PointSpatialIndex:
    """Reusable deterministic cell index for repeated point-radius queries."""

    def __init__(self, coordinates: np.ndarray, box: np.ndarray | None, cell_size: float):
        self.coordinates = _coordinates(coordinates)
        _validate_cutoff(cell_size)
        self.periodic, self.lengths = _periodic_box(box)
        self.cell_size = float(cell_size)
        if self.periodic:
            assert self.lengths is not None
            shape_array = np.maximum(1, np.floor(self.lengths / self.cell_size).astype(int))
            self.shape: tuple[int, int, int] | None = tuple(int(value) for value in shape_array)
            self.cell_width = self.lengths / shape_array
            self.origin = np.zeros(3, dtype=float)
        else:
            self.shape = None
            self.cell_width = np.full(3, self.cell_size, dtype=float)
            self.origin = np.min(self.coordinates, axis=0) if len(self.coordinates) else np.zeros(3, dtype=float)
        self.cells: dict[tuple[int, int, int], list[int]] = {}
        for index, coordinate in enumerate(self.coordinates):
            self.cells.setdefault(self._key(coordinate), []).append(index)

    def query(self, point: np.ndarray, cutoff: float) -> tuple[int, ...]:
        """Return sorted coordinate indices within an exact PBC-aware cutoff."""
        _validate_cutoff(cutoff)
        if not len(self.coordinates):
            return ()
        key = self._key(point)
        spans = np.maximum(1, np.asarray([ceil(cutoff / width) for width in self.cell_width], dtype=int))
        candidate_keys: set[tuple[int, int, int]] = set()
        for offset in product(*(range(-int(span), int(span) + 1) for span in spans)):
            if self.periodic:
                assert self.shape is not None
                other = tuple((key[axis] + offset[axis]) % self.shape[axis] for axis in range(3))
            else:
                other = tuple(key[axis] + offset[axis] for axis in range(3))
            candidate_keys.add(other)
        candidates = sorted(index for other in candidate_keys for index in self.cells.get(other, ()))
        if not candidates:
            return ()
        deltas = minimum_image(
            self.coordinates[candidates] - np.asarray(point, dtype=float),
            self.lengths if self.periodic else None,
        )
        distances2 = np.einsum("ij,ij->i", deltas, deltas)
        cutoff2 = float(cutoff) ** 2
        return tuple(
            index
            for index, distance2 in zip(candidates, distances2, strict=True)
            if float(distance2) <= cutoff2
        )

    def _key(self, point: np.ndarray) -> tuple[int, int, int]:
        value = np.asarray(point, dtype=float)
        if self.periodic:
            assert self.lengths is not None and self.shape is not None
            scaled = np.floor(np.mod(value, self.lengths) / self.cell_width).astype(int)
            shape = np.asarray(self.shape, dtype=int)
            scaled = np.minimum(np.maximum(scaled, 0), shape - 1)
        else:
            scaled = np.floor((value - self.origin) / self.cell_width).astype(int)
        return tuple(int(item) for item in scaled)


def self_cutoff_pairs(
    coordinates: np.ndarray,
    box: np.ndarray | None,
    cutoff: float,
) -> list[tuple[int, int]]:
    """Return sorted unique ``i < j`` pairs separated by at most ``cutoff``."""
    coords = _coordinates(coordinates)
    _validate_cutoff(cutoff)
    if len(coords) < 2:
        return []
    periodic, lengths = _periodic_box(box)
    accelerated_pairs = _mda_capped_pairs(coords, coords, lengths if periodic else None, cutoff, self_pairs=True)
    if accelerated_pairs is not None:
        return _filter_self_pairs(coords, accelerated_pairs, lengths if periodic else None, cutoff)
    keys, shape = _cell_keys(coords, cutoff, periodic, lengths)
    cells: dict[tuple[int, int, int], list[int]] = {}
    for index, key in enumerate(keys):
        cells.setdefault(key, []).append(index)
    cutoff2 = float(cutoff) ** 2
    pairs: list[tuple[int, int]] = []
    for key in sorted(cells):
        left = cells[key]
        for other_key in _neighbor_keys(key, periodic, shape):
            if other_key not in cells or other_key < key:
                continue
            right = cells[other_key]
            for i in left:
                for j in right:
                    if other_key == key and j <= i:
                        continue
                    delta = minimum_image(coords[j] - coords[i], lengths if periodic else None)
                    if float(np.dot(delta, delta)) <= cutoff2:
                        pairs.append((i, j) if i < j else (j, i))
    return sorted(set(pairs))


def cross_cutoff_pairs(
    left_coordinates: np.ndarray,
    right_coordinates: np.ndarray,
    box: np.ndarray | None,
    cutoff: float,
) -> list[tuple[int, int]]:
    """Return sorted ``(left_index, right_index)`` pairs within ``cutoff``."""
    left = _coordinates(left_coordinates)
    right = _coordinates(right_coordinates)
    _validate_cutoff(cutoff)
    if not len(left) or not len(right):
        return []
    periodic, lengths = _periodic_box(box)
    accelerated_pairs = _mda_capped_pairs(left, right, lengths if periodic else None, cutoff, self_pairs=False)
    if accelerated_pairs is not None:
        return _filter_cross_pairs(left, right, accelerated_pairs, lengths if periodic else None, cutoff)
    combined = np.vstack((left, right))
    keys, shape = _cell_keys(combined, cutoff, periodic, lengths)
    left_keys = keys[: len(left)]
    right_keys = keys[len(left) :]
    right_cells: dict[tuple[int, int, int], list[int]] = {}
    for index, key in enumerate(right_keys):
        right_cells.setdefault(key, []).append(index)
    cutoff2 = float(cutoff) ** 2
    pairs: list[tuple[int, int]] = []
    for i, key in enumerate(left_keys):
        for other_key in _neighbor_keys(key, periodic, shape):
            for j in right_cells.get(other_key, ()):
                delta = minimum_image(right[j] - left[i], lengths if periodic else None)
                if float(np.dot(delta, delta)) <= cutoff2:
                    pairs.append((i, j))
    return sorted(set(pairs))


def _mda_capped_pairs(
    left: np.ndarray,
    right: np.ndarray,
    lengths: np.ndarray | None,
    cutoff: float,
    *,
    self_pairs: bool,
) -> np.ndarray | None:
    """Return MDAnalysis cutoff candidates, or ``None`` for the pure-Python fallback."""
    # Use MDAnalysis only for large candidate searches.
    if len(left) < 32 or len(right) < 32:
        return None
    try:
        from MDAnalysis.lib.distances import capped_distance
    except ImportError:
        return None
    box = None
    if lengths is not None:
        box = np.asarray([*lengths, 90.0, 90.0, 90.0], dtype=np.float32)
    try:
        pairs = capped_distance(
            left,
            right,
            # Expand in float32 so boundary pairs reach the float64 recheck.
            max_cutoff=float(np.nextafter(np.float32(cutoff), np.float32(np.inf))),
            box=box,
            return_distances=False,
        )
    except Exception:
        # Keep the deterministic fallback available.
        return None
    pairs = np.asarray(pairs, dtype=int)
    if pairs.size == 0:
        return np.empty((0, 2), dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        return None
    if self_pairs:
        pairs = pairs[pairs[:, 0] < pairs[:, 1]]
    return pairs


def _filter_self_pairs(
    coordinates: np.ndarray,
    candidates: np.ndarray,
    box: np.ndarray | None,
    cutoff: float,
) -> list[tuple[int, int]]:
    """Apply SQQ's exact inclusive cutoff after accelerated candidate search."""
    if not len(candidates):
        return []
    deltas = minimum_image(coordinates[candidates[:, 1]] - coordinates[candidates[:, 0]], box)
    distances2 = np.einsum("ij,ij->i", deltas, deltas)
    cutoff2 = float(cutoff) ** 2
    return sorted(
        {
            (int(i), int(j))
            for (i, j), distance2 in zip(candidates, distances2, strict=True)
            if float(distance2) <= cutoff2
        }
    )


def _filter_cross_pairs(
    left: np.ndarray,
    right: np.ndarray,
    candidates: np.ndarray,
    box: np.ndarray | None,
    cutoff: float,
) -> list[tuple[int, int]]:
    """Apply SQQ's exact inclusive cross-set cutoff after acceleration."""
    if not len(candidates):
        return []
    deltas = minimum_image(right[candidates[:, 1]] - left[candidates[:, 0]], box)
    distances2 = np.einsum("ij,ij->i", deltas, deltas)
    cutoff2 = float(cutoff) ** 2
    return sorted(
        (int(i), int(j))
        for (i, j), distance2 in zip(candidates, distances2, strict=True)
        if float(distance2) <= cutoff2
    )


def _coordinates(values: np.ndarray) -> np.ndarray:
    coords = np.asarray(values, dtype=float)
    if coords.size == 0:
        return np.empty((0, 3), dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coordinates must have shape (N, 3)")
    return coords


def _validate_cutoff(cutoff: float) -> None:
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be positive and finite")


def _periodic_box(box: np.ndarray | None) -> tuple[bool, np.ndarray | None]:
    if box is None or len(box) < 3:
        return False, None
    lengths = np.asarray(box[:3], dtype=float)
    if np.any(~np.isfinite(lengths)) or np.any(lengths <= 0):
        return False, None
    return True, lengths


def _cell_keys(
    coords: np.ndarray,
    cutoff: float,
    periodic: bool,
    lengths: np.ndarray | None,
) -> tuple[list[tuple[int, int, int]], tuple[int, int, int] | None]:
    if periodic:
        assert lengths is not None
        shape_array = np.maximum(1, np.floor(lengths / cutoff).astype(int))
        wrapped = np.mod(coords, lengths)
        scaled = np.floor(wrapped / (lengths / shape_array)).astype(int)
        shape = tuple(int(value) for value in shape_array)
        return [tuple(int(value) for value in row) for row in scaled], shape
    origin = np.min(coords, axis=0)
    scaled = np.floor((coords - origin) / cutoff).astype(int)
    return [tuple(int(value) for value in row) for row in scaled], None


def _neighbor_keys(
    key: tuple[int, int, int],
    periodic: bool,
    shape: tuple[int, int, int] | None,
) -> list[tuple[int, int, int]]:
    if periodic:
        assert shape is not None
        keys = {
            tuple((key[axis] + offset[axis]) % shape[axis] for axis in range(3))
            for offset in NEIGHBOR_OFFSETS
        }
    else:
        keys = {
            tuple(key[axis] + offset[axis] for axis in range(3))
            for offset in NEIGHBOR_OFFSETS
        }
    return sorted(keys)
