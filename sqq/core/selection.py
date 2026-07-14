from __future__ import annotations

"""Molecule selection helpers."""

import numpy as np

from ..models import Atom, Guest, Water


def normalize_name(name: str) -> str:
    """Normalize residue and atom names for case-insensitive matching."""
    return name.upper()


def normalize_names(names: set[str]) -> set[str]:
    """Normalize a configured residue or atom-name set."""
    return {normalize_name(name) for name in names}


def selected_residue_groups(
    atoms: list[Atom],
    resnames: set[str],
) -> list[tuple[str, int, list[int]]]:
    """Group contiguous selected residues without merging wrapped residue ids."""
    selected = normalize_names(resnames)
    groups: list[tuple[str, int, list[int]]] = []
    current_key: tuple[str, int] | None = None
    current_indices: list[int] = []
    for atom in atoms:
        key = (normalize_name(atom.resname), atom.resid)
        if key != current_key:
            if current_key is not None and current_indices:
                groups.append((*current_key, current_indices))
            current_key = key
            current_indices = []
        if key[0] in selected:
            current_indices.append(atom.index)
    if current_key is not None and current_indices:
        groups.append((*current_key, current_indices))
    return groups


def select_waters(atoms: list[Atom], resnames: set[str], oxygen_names: set[str], hydrogen_names: set[str]) -> list[Water]:
    """Group atoms into water molecules and identify oxygen/hydrogen atoms."""
    oxygen_keys = normalize_names(oxygen_names)
    hydrogen_keys = normalize_names(hydrogen_names)

    waters: list[Water] = []
    for _, resid, indices in selected_residue_groups(atoms, resnames):
        resname = atoms[indices[0]].resname
        # Single-site water models are allowed when no explicit oxygen name exists.
        oxygen = next((idx for idx in indices if normalize_name(atoms[idx].atomname) in oxygen_keys), None)
        if oxygen is None and len(indices) == 1:
            oxygen = indices[0]
        if oxygen is None:
            continue
        hydrogens = tuple(idx for idx in indices if normalize_name(atoms[idx].atomname) in hydrogen_keys)
        waters.append(Water(resid=resid, resname=resname, oxygen=oxygen, hydrogens=hydrogens, atoms=tuple(indices)))
    return waters


def select_guests(atoms: list[Atom], resnames: set[str], center_atoms: dict[str, list[str]], center_mode: str = "center_atom") -> list[Guest]:
    """Group guest residues and select their preferred center atom when present."""
    mode = str(center_mode or "center_atom").strip().lower()
    if mode not in {"center_atom", "centroid", "auto"}:
        raise ValueError("guest.center_mode must be center_atom, centroid, or auto.")
    center_atom_keys = {
        normalize_name(resname): normalize_names(set(atom_names))
        for resname, atom_names in center_atoms.items()
    }
    guests: list[Guest] = []
    for resname_key, resid, indices in selected_residue_groups(atoms, resnames):
        resname = atoms[indices[0]].resname
        preferred = center_atom_keys.get(resname_key, set())
        center = None if mode == "centroid" else next((idx for idx in indices if normalize_name(atoms[idx].atomname) in preferred), None)
        guests.append(Guest(resid=resid, resname=resname, atoms=tuple(indices), center_atom=center))
    return guests


def water_by_oxygen(waters: list[Water]) -> dict[int, Water]:
    """Index water molecules by oxygen atom index."""
    return {water.oxygen: water for water in waters}


def centroid(atoms: list[Atom], indices: tuple[int, ...] | list[int]) -> np.ndarray:
    """Return the geometric center of selected atoms."""
    return np.mean([atoms[idx].xyz for idx in indices], axis=0)
