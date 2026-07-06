from __future__ import annotations

"""Trajectory and coordinate readers."""

import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from ..models import Atom, Frame


SUPPORTED_SUFFIXES = {".gro", ".xyz", ".xtc", ".trr"}


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


def read_frames(paths: list[Path], topology: Path | None = None, xyz_scale: float = 0.1, xtc_stride: int = 1) -> Iterable[Frame]:
    """Yield frames from supported coordinate and trajectory formats."""
    for path in paths:
        suffix = path.suffix.lower()
        if suffix == ".gro":
            yield read_gro(path)
        elif suffix == ".xyz":
            yield read_xyz(path, scale=xyz_scale)
        elif suffix in {".xtc", ".trr"}:
            yield from read_mdanalysis(path, topology, stride=xtc_stride)
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

    atoms = [_parse_gro_atom(i, line) for i, line in enumerate(lines[2 : 2 + natoms])]
    box = None
    if len(lines) > 2 + natoms:
        parts = lines[2 + natoms].split()
        if len(parts) >= 3:
            box = np.asarray([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
    return Frame(name=path.stem, atoms=atoms, box=box, time_ps=_parse_title_time_ps(title), source=path)


def _parse_gro_atom(index: int, line: str) -> Atom:
    try:
        # Standard GRO uses fixed-width columns for residue, atom, and xyz.
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomid = int(line[15:20])
        xyz = np.asarray([float(line[20:28]), float(line[28:36]), float(line[36:44])], dtype=float)
    except ValueError:
        # Some generated GRO-like files are whitespace-separated; accept them.
        parts = line.split()
        if len(parts) < 6:
            raise
        head = parts[0]
        digits = "".join(ch for ch in head if ch.isdigit())
        letters = "".join(ch for ch in head if not ch.isdigit())
        resid = int(digits or 0)
        resname = letters
        atomname = parts[1]
        atomid = int(parts[2])
        xyz = np.asarray([float(parts[3]), float(parts[4]), float(parts[5])], dtype=float)
    return Atom(index=index, resid=resid, resname=resname, atomname=atomname, atomid=atomid, xyz=xyz)


def _parse_title_time_ps(title: str) -> float | None:
    """Extract `t=` from a GRO title when present."""
    match = re.search(r"\bt\s*=\s*([-+0-9.eE]+)", title)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def read_xyz(path: Path, scale: float = 0.1) -> Frame:
    """Read an XYZ file and convert coordinates to nm by default."""
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {path}")
    natoms = int(lines[0].strip())
    atoms = []
    for index, line in enumerate(lines[2 : 2 + natoms]):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom line in {path}: {line!r}")
        xyz = np.asarray([float(parts[1]), float(parts[2]), float(parts[3])], dtype=float) * scale
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


def trajectory_frame_indices(path: Path, topology: Path | None, stride: int = 1) -> list[int]:
    """Return raw trajectory indexes selected by the configured stride."""
    universe = open_mdanalysis_universe(path, topology)
    try:
        step = max(1, int(stride))
        return list(range(0, len(universe.trajectory), step))
    finally:
        close_mdanalysis_universe(universe)


def frame_from_mdanalysis_universe(universe, path: Path, raw_frame_index: int) -> Frame:
    """Materialize one selected MDAnalysis frame as the SQQ data model."""
    ts = universe.trajectory[int(raw_frame_index)]
    positions_nm = universe.atoms.positions / 10.0
    atoms = [
        Atom(
            index=index,
            resid=int(atom.resid),
            resname=str(atom.resname),
            atomname=str(atom.name),
            atomid=int(atom.id),
            xyz=np.asarray(xyz, dtype=float),
        )
        for index, (atom, xyz) in enumerate(zip(universe.atoms, positions_nm, strict=True))
    ]
    box = None
    if ts.dimensions is not None and len(ts.dimensions) >= 3:
        box = np.asarray(ts.dimensions[:3], dtype=float) / 10.0
    return Frame(
        name=f"{path.stem}_frame{ts.frame:06d}",
        atoms=atoms,
        box=box,
        time_ps=float(ts.time),
        source=path,
    )


def read_mdanalysis(path: Path, topology: Path | None, stride: int = 1) -> Iterable[Frame]:
    """Read XTC/TRR through MDAnalysis using a separate topology file."""
    universe = open_mdanalysis_universe(path, topology)
    try:
        for raw_frame_index in trajectory_frame_indices_from_length(len(universe.trajectory), stride):
            yield frame_from_mdanalysis_universe(universe, path, raw_frame_index)
    finally:
        close_mdanalysis_universe(universe)


def trajectory_frame_indices_from_length(length: int, stride: int = 1) -> range:
    """Return a lazy raw-index range for an already-open trajectory."""
    return range(0, int(length), max(1, int(stride)))

