"""Atomic, molecular, and frame data contracts."""

from __future__ import annotations

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


def guest_id(guest: Guest) -> str:
    """Return an unambiguous topology-stable guest identifier."""
    if not guest.atoms:
        raise ValueError("A guest molecule must contain at least one atom.")
    return f"g{min(int(index) for index in guest.atoms):08d}"


__all__ = ["Atom", "Frame", "Water", "Guest", "guest_id"]
