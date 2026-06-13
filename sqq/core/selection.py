from __future__ import annotations

"""Molecule selection helpers."""

from collections import defaultdict

import numpy as np

from ..models import Atom, Guest, Water


def select_waters(atoms: list[Atom], resnames: set[str], oxygen_names: set[str], hydrogen_names: set[str]) -> list[Water]:
    """Group atoms into water molecules and identify oxygen/hydrogen atoms."""
    groups: dict[tuple[str, int], list[int]] = defaultdict(list)
    for atom in atoms:
        if atom.resname in resnames:
            groups[(atom.resname, atom.resid)].append(atom.index)

    waters: list[Water] = []
    for (resname, resid), indices in sorted(groups.items(), key=lambda item: (item[0][1], item[0][0])):
        # Single-site water models are allowed when no explicit oxygen name exists.
        oxygen = next((idx for idx in indices if atoms[idx].atomname in oxygen_names), None)
        if oxygen is None and len(indices) == 1:
            oxygen = indices[0]
        if oxygen is None:
            continue
        hydrogens = tuple(idx for idx in indices if atoms[idx].atomname in hydrogen_names)
        waters.append(Water(resid=resid, resname=resname, oxygen=oxygen, hydrogens=hydrogens, atoms=tuple(indices)))
    return waters


def select_guests(atoms: list[Atom], resnames: set[str], center_atoms: dict[str, list[str]]) -> list[Guest]:
    """Group guest residues and select their preferred center atom when present."""
    groups: dict[tuple[str, int], list[int]] = defaultdict(list)
    for atom in atoms:
        if atom.resname in resnames:
            groups[(atom.resname, atom.resid)].append(atom.index)

    guests: list[Guest] = []
    for (resname, resid), indices in sorted(groups.items(), key=lambda item: (item[0][1], item[0][0])):
        preferred = set(center_atoms.get(resname, []))
        center = next((idx for idx in indices if atoms[idx].atomname in preferred), None)
        guests.append(Guest(resid=resid, resname=resname, atoms=tuple(indices), center_atom=center))
    return guests


def water_by_oxygen(waters: list[Water]) -> dict[int, Water]:
    """Index water molecules by oxygen atom index."""
    return {water.oxygen: water for water in waters}


def centroid(atoms: list[Atom], indices: tuple[int, ...] | list[int]) -> np.ndarray:
    """Return the geometric center of selected atoms."""
    return np.mean([atoms[idx].xyz for idx in indices], axis=0)
