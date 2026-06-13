from __future__ import annotations

"""GRO writers for SQQ visualization outputs."""

from pathlib import Path

import numpy as np

from ..core.pbc import unwrap_path
from ..core.selection import water_by_oxygen
from ..models import Atom, Cage, Cup, Frame, FrameResult, Guest, Ring, Water


RING_CENTER_NAMES = {4: "R4", 5: "R5", 6: "R6", 7: "R7"}
CAGE_CENTER_NAMES = {"512": "G512", "51262": "G62", "51263": "G63", "51264": "G64"}


def write_ring_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False) -> None:
    """Write free-ring GRO files after removing rings used by cups/cages."""
    water_lookup = water_by_oxygen(result.waters)
    used_ring_ids = {ring_id for cup in result.cups for ring_id in cup.rings}
    used_ring_ids.update(ring_id for cage in result.cages for ring_id in cage.rings)
    for size, rings in sorted(result.rings.items()):
        rings = [ring for ring in rings if ring.object_id not in used_ring_ids]
        if not rings and not write_empty:
            continue
        atoms = aggregate_ring_atoms(result.frame, rings, water_lookup, result.frame.box)
        if atoms or write_empty:
            write_gro(frame_dir / f"{result.frame.name}_ring{size}.gro", f"{result.frame.name} free ring{size}", atoms, result.frame.box)


def aggregate_ring_atoms(frame: Frame, rings: list[Ring], water_lookup: dict[int, Water], box: np.ndarray | None) -> list[Atom]:
    """Collect full water molecules and one CNT atom per ring."""
    atom_indices: dict[int, Atom] = {}
    centers: list[Atom] = []
    for serial, ring in enumerate(rings, start=1):
        for oxygen in ring.nodes:
            water = water_lookup.get(oxygen)
            if water is None:
                continue
            for atom_idx in water.atoms:
                atom_indices[atom_idx] = frame.atoms[atom_idx]
        points = unwrap_path([frame.atoms[idx].xyz for idx in ring.nodes], box)
        center = np.mean(points, axis=0)
        # Ring center atom names are short so VMD selections stay simple.
        centers.append(
            Atom(
                index=-serial,
                resid=90000 + serial,
                resname="CNT",
                atomname=RING_CENTER_NAMES.get(ring.size, "CNT"),
                atomid=90000 + serial,
                xyz=center,
            )
        )
    return [atom_indices[idx] for idx in sorted(atom_indices)] + centers


def write_cup_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False) -> None:
    """Write one GRO file per cup type."""
    if not result.cups and not write_empty:
        return
    water_lookup = water_by_oxygen(result.waters)
    by_type: dict[str, list[Cup]] = {}
    for cup in result.cups:
        by_type.setdefault(cup.cup_type, []).append(cup)
    for cup_type, cups in sorted(by_type.items()):
        atoms = aggregate_cup_atoms(result.frame, cups, water_lookup)
        if atoms or write_empty:
            write_gro(frame_dir / f"{result.frame.name}_{cup_type}.gro", f"{result.frame.name} {cup_type}", atoms, result.frame.box)


def aggregate_cup_atoms(frame: Frame, cups: list[Cup], water_lookup: dict[int, Water]) -> list[Atom]:
    """Collect full water molecules and one CNT atom per cup."""
    atom_indices: dict[int, Atom] = {}
    centers: list[Atom] = []
    for serial, cup in enumerate(cups, start=1):
        for oxygen in cup.waters:
            water = water_lookup.get(oxygen)
            if water is None:
                continue
            for atom_idx in water.atoms:
                atom_indices[atom_idx] = frame.atoms[atom_idx]
        base_size = cup.cup_type.removeprefix("cup").split("_", 1)[0]
        centers.append(
            Atom(
                index=-100000 - serial,
                resid=91000 + serial,
                resname="CNT",
                atomname=f"CP{base_size}"[:5],
                atomid=91000 + serial,
                xyz=cup.center,
            )
        )
    return [atom_indices[idx] for idx in sorted(atom_indices)] + centers


def write_cage_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False) -> None:
    """Write cage GRO files by type and occupancy state."""
    if not result.cages and not write_empty:
        return
    water_lookup = water_by_oxygen(result.waters)
    guest_lookup = {guest_id(guest): guest for guest in result.guests}
    groups: dict[str, list[Cage]] = {}
    for cage in result.cages:
        groups.setdefault(cage.cage_type, []).append(cage)
        groups.setdefault(f"{cage.cage_type}_{'occupied' if cage.occupied else 'empty'}", []).append(cage)
        if cage.occupied:
            guest_names = sorted({guest_lookup[item].resname for item in cage.guest_ids if item in guest_lookup})
            if len(cage.guest_ids) > 1 or len(guest_names) > 1:
                groups.setdefault(f"{cage.cage_type}_multi", []).append(cage)
            elif guest_names:
                groups.setdefault(f"{cage.cage_type}_{guest_names[0]}", []).append(cage)

    for label, cages in sorted(groups.items()):
        atoms = aggregate_cage_atoms(result.frame, cages, water_lookup, guest_lookup)
        if atoms or write_empty:
            write_gro(frame_dir / f"{result.frame.name}_{label}.gro", f"{result.frame.name} {label}", atoms, result.frame.box)


def aggregate_cage_atoms(
    frame: Frame,
    cages: list[Cage],
    water_lookup: dict[int, Water],
    guest_lookup: dict[str, Guest],
) -> list[Atom]:
    """Collect cage waters, assigned guests, and one CNT atom per cage."""
    atom_indices: dict[int, Atom] = {}
    centers: list[Atom] = []
    for serial, cage in enumerate(cages, start=1):
        for oxygen in cage.waters:
            water = water_lookup.get(oxygen)
            if water is not None:
                for atom_idx in water.atoms:
                    atom_indices[atom_idx] = frame.atoms[atom_idx]
        for gid in cage.guest_ids:
            guest = guest_lookup.get(gid)
            if guest is not None:
                for atom_idx in guest.atoms:
                    atom_indices[atom_idx] = frame.atoms[atom_idx]
        centers.append(
            Atom(
                index=-200000 - serial,
                resid=92000 + serial,
                resname="CNT",
                atomname=CAGE_CENTER_NAMES.get(cage.cage_type, "CAGE")[:5],
                atomid=92000 + serial,
                xyz=cage.center,
            )
        )
    return [atom_indices[idx] for idx in sorted(atom_indices)] + centers


def guest_id(guest: Guest) -> str:
    """Use the same guest identifier as the cage assignment code."""
    return f"{guest.resname}{guest.resid}"


def write_ice_gro_file(result: FrameResult, frame_dir: Path, write_empty: bool = False) -> None:
    """Write ice-like, ice-I-like, and interfacial ice GRO files."""
    groups = {
        "ice": result.ice_like_waters,
        "iceI": result.ice_i_waters,
        "ice_interfacial": result.interfacial_ice_waters,
    }
    if not any(groups.values()) and not write_empty:
        return
    for label, oxygens in groups.items():
        atoms = aggregate_water_atoms(result.frame, result.waters, oxygens)
        if atoms or write_empty:
            write_gro(frame_dir / f"{result.frame.name}_{label}.gro", f"{result.frame.name} {label}", atoms, result.frame.box)


def aggregate_water_atoms(frame: Frame, waters: list[Water], oxygens: tuple[int, ...]) -> list[Atom]:
    """Collect full water molecules for a set of oxygen nodes."""
    water_lookup = water_by_oxygen(waters)
    atom_indices: dict[int, Atom] = {}
    for oxygen in oxygens:
        water = water_lookup.get(oxygen)
        if water is None:
            continue
        for atom_idx in water.atoms:
            atom_indices[atom_idx] = frame.atoms[atom_idx]
    return [atom_indices[idx] for idx in sorted(atom_indices)]


def write_gro(path: Path, title: str, atoms: list[Atom], box: np.ndarray | None) -> None:
    """Write a minimal, VMD-friendly GRO file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(atoms):5d}\n")
        for serial, atom in enumerate(atoms, start=1):
            resid = atom.resid % 100000
            atomid = serial % 100000
            handle.write(
                f"{resid:5d}{atom.resname[:5]:>5}{atom.atomname[:5]:>5}{atomid:5d}"
                f"{atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}\n"
            )
        if box is not None and len(box) >= 3:
            handle.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")
        else:
            handle.write("   0.00000   0.00000   0.00000\n")

