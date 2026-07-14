from __future__ import annotations

"""Orthorhombic PBC helpers."""

import numpy as np


def minimum_image(delta: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    """Apply the minimum-image convention for an orthorhombic box."""
    result = np.asarray(delta, dtype=float).copy()
    if box is None:
        return result
    if len(box) < 3:
        return result
    ortho = np.asarray(box[:3], dtype=float)
    periodic = np.isfinite(ortho) & (ortho > 0)
    if np.all(periodic):
        # Fast path for a valid orthorhombic box.
        result -= ortho * np.round(result / ortho)
    elif np.any(periodic):
        # Keep valid axes periodic for incomplete boxes.
        result[..., periodic] -= ortho[periodic] * np.round(
            result[..., periodic] / ortho[periodic]
        )
    return result


def distance(a: np.ndarray, b: np.ndarray, box: np.ndarray | None) -> float:
    """Return a minimum-image distance."""
    return float(np.linalg.norm(minimum_image(a - b, box)))


def unwrap_path(points: list[np.ndarray], box: np.ndarray | None) -> np.ndarray:
    """Unwrap a path by keeping each next point near the previous point."""
    if not points:
        return np.empty((0, 3), dtype=float)
    unwrapped = [np.asarray(points[0], dtype=float)]
    for previous, current in zip(points, points[1:], strict=False):
        unwrapped.append(unwrapped[-1] + minimum_image(current - previous, box))
    return np.asarray(unwrapped, dtype=float)
