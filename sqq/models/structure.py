from __future__ import annotations

"""Atomic, molecular, and frame data contracts."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Atom:
    index: int
    resid: int
    resname: str
    atomname: str
    atomid: int
    xyz: np.ndarray
    velocity: np.ndarray | None = None
    molecule_id: int | None = None


@dataclass
class Frame:
    name: str
    atoms: list[Atom]
    box: np.ndarray | None = None
    time_ps: float | None = None
    source: Path | None = None


@dataclass(frozen=True)
class Water:
    resid: int
    resname: str
    oxygen: int
    hydrogens: tuple[int, ...]
    atoms: tuple[int, ...]


@dataclass(frozen=True)
class Guest:
    resid: int
    resname: str
    atoms: tuple[int, ...]
    center_atom: int | None = None


__all__ = ["Atom", "Frame", "Water", "Guest"]
