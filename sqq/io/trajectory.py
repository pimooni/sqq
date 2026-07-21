from __future__ import annotations

"""Trajectory and coordinate readers."""

import re
from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np

from ..models import Atom, Frame
from .lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    lammps_trajectory_frame_indices,
    read_lammps,
)


SUPPORTED_SUFFIXES = {".gro", ".xyz", ".xtc", ".trr"} | set(
    LAMMPS_TRAJECTORY_SUFFIXES
)
BOX_TOLERANCE = 1.0e-8


def expand_inputs(input_path: Path, pattern: str, recursive: bool) -> list[Path]:
    """Expand a single input file, a directory pattern, or a direct glob."""
    if input_path.is_dir():
        iterator = input_path.rglob(pattern) if recursive else input_path.glob(pattern)
        paths = [path for path in iterator if path.is_file() and path.suffix.lower() in SUPPORTED_SUFFIXES]
    elif has_glob_magic(str(input_path)):
        parent = input_path.parent if str(input_path.parent) not in {"", "."} else Path(".")
        iterator = parent.glob(input_path.name)
        paths = [path for path in iterator if path.is_file() and path.suffix.lower() in SUPPORTED_SUFFIXES]
    else:
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")
        if not input_path.is_file():
            raise ValueError(f"Input path is not a file: {input_path}")
        if input_path.suffix.lower() not in SUPPORTED_SUFFIXES:
            raise ValueError(f"Unsupported input format: {input_path}")
        paths = [input_path]
    paths = sorted(paths, key=natural_key)
    if not paths:
        raise FileNotFoundError(f"No input files matched: {input_path} / {pattern}")
    return paths


def has_glob_magic(text: str) -> bool:
    """Return whether an input path contains glob wildcards."""
    return any(char in text for char in "*?[")


def natural_key(path: Path) -> list[object]:
    """Sort names like 1.gro, 2.gro, 10.gro in numeric order."""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", path.name)]


def read_frames(
    paths: list[Path],
    topology: Path | None = None,
    xyz_scale: float = 0.1,
    trajectory_stride: int = 1,
    lammps_config: Mapping[str, object] | None = None,
) -> Iterable[Frame]:
    """Yield frames from supported coordinate and trajectory formats."""
    stride = validated_stride(trajectory_stride)
    for path in paths:
        suffix = path.suffix.lower()
        if suffix == ".gro":
            yield read_gro(path)
        elif suffix == ".xyz":
            yield read_xyz(path, scale=xyz_scale)
        elif suffix in {".xtc", ".trr"}:
            yield from read_mdanalysis(path, topology, stride=stride)
        elif suffix in LAMMPS_TRAJECTORY_SUFFIXES:
            values = dict(lammps_config or {})
            values["stride"] = stride
            yield from read_lammps(path, topology, values)
        else:
            raise ValueError(f"Unsupported input format: {path}")


def read_gro(path: Path) -> Frame:
    """Read one GROMACS GRO frame."""
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file: {path}")
    title = lines[0].strip()
    try:
        natoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid atom count in GRO file: {path}") from exc
    if natoms < 0:
        raise ValueError(f"Negative atom count in GRO file: {path}")
    atom_end = 2 + natoms
    if len(lines) <= atom_end:
        remaining = max(0, len(lines) - 2)
        raise ValueError(
            f"Truncated GRO file {path}: declared {natoms} atoms and requires "
            f"a separate box line, but found only {remaining} lines after the header."
        )

    atoms = [_parse_gro_atom(i, line) for i, line in enumerate(lines[2:atom_end])]
    box = parse_gro_box(lines[atom_end], path)
    if any(line.strip() for line in lines[atom_end + 1 :]):
        raise ValueError(
            f"Invalid GRO file: {path} contains extra non-empty records after "
            "the box line; SQQ accepts one GRO frame per file."
        )
    return Frame(name=path.stem, atoms=atoms, box=box, time_ps=_parse_title_time_ps(title), source=path)


def parse_gro_box(line: str, path: Path | None = None) -> np.ndarray | None:
    """Parse an orthorhombic 3/9-value GRO box and reject triclinic tilt."""
    parts = line.split()
    if len(parts) not in {3, 9}:
        source = f" in {path}" if path is not None else ""
        raise ValueError(f"GRO box must contain 3 or 9 values{source}; got {len(parts)}.")
    try:
        values = np.asarray([float(value) for value in parts], dtype=float)
    except ValueError as exc:
        source = f" in {path}" if path is not None else ""
        raise ValueError(f"Invalid GRO box values{source}: {line!r}") from exc
    if np.any(~np.isfinite(values)):
        source = f" in {path}" if path is not None else ""
        raise ValueError(f"Non-finite GRO box values{source}.")
    if len(values) == 9 and np.any(np.abs(values[3:]) > BOX_TOLERANCE):
        source = f" in {path}" if path is not None else ""
        raise ValueError(
            "Triclinic GRO boxes are not supported"
            f"{source}; convert the frame to an orthorhombic representation first."
        )
    lengths = values[:3]
    if np.all(np.abs(lengths) <= BOX_TOLERANCE):
        return None
    if np.any(lengths <= 0):
        source = f" in {path}" if path is not None else ""
        raise ValueError(f"GRO box lengths must be positive or all zero{source}.")
    return lengths


def _parse_gro_atom(index: int, line: str) -> Atom:
    record = line.split(";", 1)[0].rstrip()
    fixed_width = True
    try:
        # Accept standard fixed-width GRO records.
        resid = int(record[0:5])
        resname = record[5:10].strip()
        atomname = record[10:15].strip()
        atomid = int(record[15:20])
        xyz = np.asarray([float(record[20:28]), float(record[28:36]), float(record[36:44])], dtype=float)
    except ValueError:
        fixed_width = False
        # Also accept whitespace-separated generated records.
        parts = line.split()
        if len(parts) < 6:
            raise
        head = parts[0]
        match = re.fullmatch(r"([+-]?\d+)(.+)", head)
        if match is None:
            raise ValueError(f"Invalid GRO residue token: {head!r}")
        resid = int(match.group(1))
        resname = match.group(2)
        atomname = parts[1]
        atomid = int(parts[2])
        xyz = np.asarray([float(parts[3]), float(parts[4]), float(parts[5])], dtype=float)
    velocity = None
    if fixed_width and len(line) >= 68 and line[44:68].strip():
        try:
            velocity = np.asarray(
                [float(line[44:52]), float(line[52:60]), float(line[60:68])],
                dtype=float,
            )
        except ValueError as exc:
            raise ValueError("Invalid GRO atom velocity fields.") from exc
    elif not fixed_width and len(parts) >= 9:
        try:
            velocity = np.asarray([float(parts[6]), float(parts[7]), float(parts[8])], dtype=float)
        except ValueError as exc:
            raise ValueError("Invalid GRO atom velocity fields.") from exc
    if not resname or not atomname:
        raise ValueError("GRO atom records require non-empty residue and atom names.")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("GRO atom coordinates must be finite.")
    if velocity is not None and not np.all(np.isfinite(velocity)):
        raise ValueError("GRO atom velocities must be finite.")
    return Atom(index=index, resid=resid, resname=resname, atomname=atomname, atomid=atomid, xyz=xyz, velocity=velocity)


def _parse_title_time_ps(title: str) -> float | None:
    """Extract `t=` from a GRO title when present."""
    match = re.search(r"\bt\s*=\s*([-+0-9.eE]+)", title)
    if not match:
        return None
    try:
        value = float(match.group(1))
        return value if np.isfinite(value) else None
    except ValueError:
        return None


def read_xyz(path: Path, scale: float = 0.1) -> Frame:
    """Read an XYZ file and convert coordinates to nm by default."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("XYZ coordinate scale must be positive and finite.")
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {path}")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom count in {path}: {lines[0]!r}") from exc
    if natoms < 0:
        raise ValueError(f"Invalid XYZ atom count in {path}: {natoms}.")
    expected_lines = natoms + 2
    if len(lines) < expected_lines:
        raise ValueError(
            f"Invalid XYZ file: {path} declares {natoms} atoms but contains "
            f"only {max(0, len(lines) - 2)} atom records."
        )
    if any(line.strip() for line in lines[expected_lines:]):
        raise ValueError(
            f"Invalid XYZ file: {path} contains extra non-empty records after "
            f"its declared {natoms} atoms; SQQ accepts one XYZ frame per file."
        )
    atoms = []
    for index, line in enumerate(lines[2:expected_lines]):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom line in {path}: {line!r}")
        try:
            xyz = np.asarray([float(parts[1]), float(parts[2]), float(parts[3])], dtype=float) * scale
        except ValueError as exc:
            raise ValueError(
                f"Invalid XYZ coordinates in {path} atom record {index + 1}: {line!r}"
            ) from exc
        if not np.all(np.isfinite(xyz)):
            raise ValueError(
                f"Invalid XYZ coordinates in {path} atom record {index + 1}: coordinates must be finite."
            )
        atoms.append(Atom(index=index, resid=index + 1, resname=parts[0], atomname=parts[0], atomid=index + 1, xyz=xyz))
    return Frame(name=path.stem, atoms=atoms, source=path)


def open_mdanalysis_universe(path: Path, topology: Path | None):
    """Open one trajectory with its topology using the optional MDAnalysis runtime."""
    if topology is None:
        raise ValueError("XTC/TRR input requires --top, for example --top topol.gro.")
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError("Reading XTC/TRR requires MDAnalysis.") from exc
    return mda.Universe(str(topology), str(path))


def close_mdanalysis_universe(universe) -> None:
    """Close a Universe trajectory reader when the backend exposes close()."""
    trajectory = getattr(universe, "trajectory", None)
    close = getattr(trajectory, "close", None)
    if callable(close):
        close()


def trajectory_frame_indices(
    path: Path,
    topology: Path | None,
    stride: int = 1,
    lammps_config: Mapping[str, object] | None = None,
) -> list[int]:
    """Return raw trajectory indexes selected by the configured stride."""
    step = validated_stride(stride)
    if path.suffix.lower() in LAMMPS_TRAJECTORY_SUFFIXES:
        values = dict(lammps_config or {})
        values["stride"] = step
        return lammps_trajectory_frame_indices(path, topology, values)
    universe = open_mdanalysis_universe(path, topology)
    try:
        return list(range(0, len(universe.trajectory), step))
    finally:
        close_mdanalysis_universe(universe)


def trajectory_atom_metadata(universe) -> tuple[tuple[int, int, str, str, int], ...]:
    """Return immutable atom metadata shared by every trajectory frame."""
    return tuple(
        (index, int(atom.resid), str(atom.resname), str(atom.name), int(atom.id))
        for index, atom in enumerate(universe.atoms)
    )


def frame_from_mdanalysis_universe(
    universe,
    path: Path,
    raw_frame_index: int,
    atom_metadata: tuple[tuple[int, int, str, str, int], ...] | None = None,
) -> Frame:
    """Materialize one selected MDAnalysis frame as the SQQ data model."""
    ts = universe.trajectory[int(raw_frame_index)]
    positions_nm = np.asarray(universe.atoms.positions, dtype=float) / 10.0
    if np.any(~np.isfinite(positions_nm)):
        raise ValueError(
            f"Non-finite trajectory coordinates in {path} frame {raw_frame_index}."
        )
    metadata = atom_metadata if atom_metadata is not None else trajectory_atom_metadata(universe)
    if len(metadata) != len(positions_nm):
        raise ValueError(f"Trajectory atom metadata does not match coordinates in {path}.")
    atoms = [
        Atom(
            index=index,
            resid=resid,
            resname=resname,
            atomname=atomname,
            atomid=atomid,
            xyz=np.asarray(xyz, dtype=float),
        )
        for (index, resid, resname, atomname, atomid), xyz in zip(metadata, positions_nm, strict=True)
    ]
    box = None
    if ts.dimensions is not None and len(ts.dimensions) >= 3:
        dimensions = np.asarray(ts.dimensions, dtype=float)
        lengths = dimensions[:3] / 10.0
        if np.all(np.abs(lengths) <= BOX_TOLERANCE):
            box = None
        elif np.any(~np.isfinite(lengths)) or np.any(lengths <= 0):
            raise ValueError(
                f"Invalid trajectory box lengths in {path} frame {raw_frame_index}: "
                f"{lengths.tolist()}."
            )
        elif len(dimensions) >= 6:
            angles = dimensions[3:6]
            if np.any(~np.isfinite(angles)) or not np.allclose(
                angles,
                90.0,
                atol=1.0e-5,
                rtol=0.0,
            ):
                raise ValueError(
                    f"Triclinic trajectory boxes are not supported: {path} "
                    f"frame {raw_frame_index} has angles {angles.tolist()}."
                )
            box = lengths
        else:
            box = lengths
    try:
        time_ps = float(ts.time)
    except (TypeError, ValueError):
        time_ps = None
    if time_ps is not None and not np.isfinite(time_ps):
        time_ps = None
    return Frame(
        name=f"{path.stem}_frame{ts.frame:06d}",
        atoms=atoms,
        box=box,
        time_ps=time_ps,
        source=path,
    )


def read_mdanalysis(path: Path, topology: Path | None, stride: int = 1) -> Iterable[Frame]:
    """Read XTC/TRR through MDAnalysis using a separate topology file."""
    universe = open_mdanalysis_universe(path, topology)
    try:
        atom_metadata = trajectory_atom_metadata(universe)
        for raw_frame_index in trajectory_frame_indices_from_length(len(universe.trajectory), stride):
            yield frame_from_mdanalysis_universe(
                universe,
                path,
                raw_frame_index,
                atom_metadata=atom_metadata,
            )
    finally:
        close_mdanalysis_universe(universe)


def trajectory_frame_indices_from_length(length: int, stride: int = 1) -> range:
    """Return a lazy raw-index range for an already-open trajectory."""
    return range(0, int(length), validated_stride(stride))


def validated_stride(value: int) -> int:
    """Return a positive trajectory stride."""
    try:
        stride = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("input.trajectory_stride must be a positive integer.") from exc
    if stride < 1:
        raise ValueError("input.trajectory_stride must be a positive integer.")
    return stride

