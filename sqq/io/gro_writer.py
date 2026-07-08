from __future__ import annotations

"""GRO writers for SQQ visualization outputs."""

from pathlib import Path
import re
import unicodedata

import numpy as np

from ..core.pbc import unwrap_path
from ..core.selection import water_by_oxygen
from ..models import Atom, Cage, CagePatch, Frame, FrameResult, Guest, Ring, Water
from .occupancy import guest_composition_label, guest_lookup, guest_resname_order


RING_CENTER_NAMES = {4: "R4", 5: "R5", 6: "R6", 7: "R7"}
CAGE_CENTER_NAMES = {"512": "G512", "51262": "G62", "51263": "G63", "51264": "G64", "51268": "G68", "435663": "G436"}
SUPERSCRIPT_DIGITS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
SUPERSCRIPT_TO_ASCII = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
SUBSCRIPT_TO_ASCII = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
SUPERSCRIPT_RUN = re.compile(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+")
SUBSCRIPT_RUN = re.compile(r"[₀₁₂₃₄₅₆₇₈₉]+")


def write_ring_gro_files(
    result: FrameResult,
    frame_dir: Path,
    write_empty: bool = False,
    layout: str = "grouped",
    sizes: set[int] | None = None,
) -> None:
    """Write reported free-ring GRO files after topology-wide filtering."""
    water_lookup = water_by_oxygen(result.waters)
    used_ring_ids = {ring_id for patch in [*result.half_cages, *result.quasi_cages] for ring_id in patch.rings}
    filtering_cages = result.all_cages or result.cages
    used_ring_ids.update(ring_id for cage in filtering_cages for ring_id in cage.rings)
    for size, rings in sorted(result.rings.items()):
        if sizes is not None and size not in sizes:
            continue
        rings = [ring for ring in rings if ring.object_id not in used_ring_ids]
        if not rings and not write_empty:
            continue
        atoms = aggregate_ring_atoms(result.frame, rings, water_lookup, result.frame.box)
        if atoms or write_empty:
            grouped_name = f"{result.frame.name}_ring_{size}.gro"
            flat_name = grouped_name
            path = ring_structure_path(frame_dir, f"ring{size}", grouped_name, flat_name, layout)
            write_gro(path, f"{result.frame.name} free ring{size}", atoms, result.frame.box)


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


def write_half_cage_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False, layout: str = "grouped") -> None:
    """Write one GRO file per half_cage type."""
    write_patch_gro_files(result, result.half_cages, frame_dir, "half_cage", write_empty, layout)


def write_quasi_cage_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False, layout: str = "grouped") -> None:
    """Write one GRO file per quasi_cage type."""
    write_patch_gro_files(result, result.quasi_cages, frame_dir, "quasi_cage", write_empty, layout)


def write_patch_gro_files(
    result: FrameResult,
    patches: list[CagePatch],
    frame_dir: Path,
    category: str,
    write_empty: bool = False,
    layout: str = "grouped",
) -> None:
    """Write one GRO file per open cage-patch type."""
    if not patches and not write_empty:
        return
    water_lookup = water_by_oxygen(result.waters)
    by_type: dict[str, list[CagePatch]] = {}
    for patch in patches:
        by_type.setdefault(patch.patch_type, []).append(patch)
    for patch_type, group in sorted(by_type.items()):
        atoms = aggregate_patch_atoms(result.frame, group, water_lookup)
        if atoms or write_empty:
            file_label = ascii_gro_text(patch_type)
            grouped_name = f"{result.frame.name}_{file_label}.gro"
            flat_name = grouped_name
            path = patch_structure_path(frame_dir, category, file_label, grouped_name, flat_name, layout)
            write_gro(path, f"{result.frame.name} {patch_type}", atoms, result.frame.box)


def aggregate_patch_atoms(frame: Frame, patches: list[CagePatch], water_lookup: dict[int, Water]) -> list[Atom]:
    """Collect full water molecules and one CNT atom per open patch."""
    atom_indices: dict[int, Atom] = {}
    centers: list[Atom] = []
    for serial, patch in enumerate(patches, start=1):
        for oxygen in patch.waters:
            water = water_lookup.get(oxygen)
            if water is None:
                continue
            for atom_idx in water.atoms:
                atom_indices[atom_idx] = frame.atoms[atom_idx]
        atomname = "HC" if patch.kind == "half_cage" else "QC"
        centers.append(
            Atom(
                index=-100000 - serial,
                resid=91000 + serial,
                resname="CNT",
                atomname=atomname,
                atomid=91000 + serial,
                xyz=patch.center,
            )
        )
    return [atom_indices[idx] for idx in sorted(atom_indices)] + centers


def write_cage_gro_files(result: FrameResult, frame_dir: Path, write_empty: bool = False, layout: str = "grouped") -> None:
    """Write cage GRO files by type and occupancy state."""
    if not result.cages and not write_empty:
        return
    water_lookup = water_by_oxygen(result.waters)
    guests_by_id = guest_lookup(result.guests)
    guest_order = guest_resname_order(result.guests)
    groups: dict[tuple[str, str], list[Cage]] = {}
    for cage in result.cages:
        groups.setdefault((cage.cage_type, ""), []).append(cage)
        groups.setdefault((cage.cage_type, "occupied" if cage.occupied else "empty"), []).append(cage)
        if cage.occupied:
            composition = guest_composition_label(cage, guests_by_id, guest_order)
            if composition:
                groups.setdefault((cage.cage_type, composition), []).append(cage)
            if len(cage.guest_ids) > 1:
                groups.setdefault((cage.cage_type, "multi"), []).append(cage)

    for (cage_type, suffix), cages in sorted(groups.items()):
        atoms = aggregate_cage_atoms(result.frame, cages, water_lookup, guests_by_id)
        if atoms or write_empty:
            grouped_label = ascii_gro_text(cage_file_label(cage_type))
            display_label = f"{grouped_label}{'_' + suffix if suffix else ''}"
            grouped_name = f"{result.frame.name}_cage_{display_label}.gro"
            flat_label = f"{cage_type}{'_' + suffix if suffix else ''}"
            flat_name = grouped_name
            path = cage_structure_path(frame_dir, grouped_label, grouped_name, flat_name, layout)
            write_gro(path, f"{result.frame.name} {flat_label}", atoms, result.frame.box)


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


def write_ice_gro_file(result: FrameResult, frame_dir: Path, write_empty: bool = False, layout: str = "grouped") -> None:
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
            grouped_name = f"{result.frame.name}_{label}.gro"
            path = structure_path(frame_dir, "ice", grouped_name, grouped_name, layout)
            write_gro(path, f"{result.frame.name} {label}", atoms, result.frame.box)


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
    """Write a minimal GRO file whose title is safe for ASCII-oriented readers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{ascii_gro_text(title)}\n")
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


def structure_path(frame_dir: Path, category: str, grouped_name: str, flat_name: str, layout: str) -> Path:
    """Return the output path for grouped or flat structure layouts."""
    if layout == "flat":
        return frame_dir / flat_name
    if layout != "grouped":
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    return frame_dir / category / grouped_name


def ring_structure_path(frame_dir: Path, ring_label: str, grouped_name: str, flat_name: str, layout: str) -> Path:
    """Place grouped ring files directly under ring/."""
    if layout == "flat":
        return frame_dir / flat_name
    if layout != "grouped":
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    return frame_dir / "ring" / grouped_name


def patch_structure_path(frame_dir: Path, category: str, patch_label: str, grouped_name: str, flat_name: str, layout: str) -> Path:
    """Place grouped patch files under half_cage/<type>/ or quasi_cage/<type>/."""
    if layout == "flat":
        return frame_dir / flat_name
    if layout != "grouped":
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    return frame_dir / category / patch_label / grouped_name


def cage_structure_path(frame_dir: Path, cage_label: str, grouped_name: str, flat_name: str, layout: str) -> Path:
    """Place grouped cage files under cage/<cage_type>/."""
    if layout == "flat":
        return frame_dir / flat_name
    if layout != "grouped":
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    return frame_dir / "cage" / cage_label / grouped_name


def cage_file_label(cage_type: str) -> str:
    """Return a readable cage filename label with superscript counts."""
    known = {
        "512": "5¹²",
        "51262": "5¹²6²",
        "51263": "5¹²6³",
        "51264": "5¹²6⁴",
    }
    if cage_type in known:
        return known[cage_type]
    counts = parse_generic_cage_label(cage_type)
    if counts:
        return "".join(f"{size}{superscript_number(count)}" for size, count in sorted(counts.items()) if count > 0)
    return cage_type


def composition_label(sequence: str) -> str:
    """Summarize a ring-size sequence with superscript counts."""
    parts = []
    for size in sorted(set(sequence), key=int):
        count = sequence.count(size)
        parts.append(f"{size}{superscript_number(count)}")
    return "".join(parts)


def superscript_number(value: int) -> str:
    """Return an integer as superscript Arabic numerals."""
    return str(value).translate(SUPERSCRIPT_DIGITS)


def ascii_gro_text(value: str) -> str:
    """Convert display-oriented Unicode structure labels to portable ASCII."""
    text = SUPERSCRIPT_RUN.sub(lambda match: "^" + match.group().translate(SUPERSCRIPT_TO_ASCII), str(value))
    text = SUBSCRIPT_RUN.sub(lambda match: "_" + match.group().translate(SUBSCRIPT_TO_ASCII), text)
    normalized = unicodedata.normalize("NFKD", text)
    return "".join(character if ord(character) < 128 else "_" for character in normalized)


def parse_generic_cage_label(label: str) -> dict[int, int]:
    """Parse labels like 4^1-5^10-6^2 into face-count maps."""
    counts: dict[int, int] = {}
    for token in label.split("-"):
        if "^" not in token:
            return {}
        size_text, count_text = token.split("^", 1)
        if not size_text.isdigit() or not count_text.isdigit():
            return {}
        counts[int(size_text)] = int(count_text)
    return counts

