from __future__ import annotations

"""Deterministic cutoff-pair searches for orthorhombic coordinate sets."""

from itertools import product

import numpy as np

from .pbc import minimum_image


NEIGHBOR_OFFSETS = tuple(product((-1, 0, 1), repeat=3))


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
