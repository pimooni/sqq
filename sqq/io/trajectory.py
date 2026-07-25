from __future__ import annotations

"""Trajectory and coordinate readers."""

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from ..models import Atom, Frame
from .lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    lammps_trajectory_times_ps,
    read_lammps,
)


SUPPORTED_SUFFIXES = {".gro", ".xyz", ".xtc", ".trr"} | set(
    LAMMPS_TRAJECTORY_SUFFIXES
)
BOX_TOLERANCE = 1.0e-8
TIME_TOLERANCE = 1.0e-6


@dataclass(frozen=True)
class TrajectorySelection:
    """Resolved time-based frame selection for one trajectory-like input."""

    raw_indexes: tuple[int, ...]
    total_frames: int
    native_interval_ps: float | None
    delta_time_ps: float | None
    raw_frame_step: int

    @property
    def selected_frames(self) -> int:
        return len(self.raw_indexes)


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


def natural_key(path: Path) -> tuple[tuple[tuple[int, int | str], ...], str]:
    """Naturally sort filenames with a stable full-path tie breaker."""
    parts = tuple(
        (0, int(part)) if part.isdigit() else (1, part.lower())
        for part in re.split(r"(\d+)", path.name)
    )
    return parts, path.as_posix().casefold()


def read_frames(
    paths: list[Path],
    topology: Path | None = None,
    xyz_scale: float = 0.1,
    frame_indexes: Sequence[int] | None = None,
    lammps_config: Mapping[str, object] | None = None,
) -> Iterable[Frame]:
    """Yield selected frames from supported coordinate and trajectory formats."""
    if frame_indexes is not None and len(paths) != 1:
        raise ValueError("Frame indexes can only select frames from one input file.")
    for path in paths:
        suffix = path.suffix.lower()
        if suffix == ".gro":
            yield from read_gro_frames(path, frame_indexes)
        elif suffix == ".xyz":
            if frame_indexes not in (None, (), [0], (0,)):
                raise ValueError("XYZ input contains one frame and only accepts frame index 0.")
            yield read_xyz(path, scale=xyz_scale)
        elif suffix in {".xtc", ".trr"}:
            yield from read_mdanalysis(path, topology, raw_indexes=frame_indexes)
        elif suffix in LAMMPS_TRAJECTORY_SUFFIXES:
            yield from read_lammps(path, topology, lammps_config, raw_indexes=frame_indexes)
        else:
            raise ValueError(f"Unsupported input format: {path}")


def read_gro(path: Path) -> Frame:
    """Read one GROMACS GRO frame and reject a stacked trajectory."""
    iterator = _iter_gro_blocks(path)
    try:
        frame = next(iterator)
    except StopIteration as exc:
        raise ValueError(f"Invalid GRO file: {path}") from exc
    try:
        next(iterator)
    except StopIteration:
        frame.name = path.stem
        return frame
    raise ValueError(
        f"GRO file contains multiple frames: {path}; use it as an analysis input "
        "instead of a single-frame topology file."
    )


def read_gro_frames(
    path: Path,
    raw_indexes: Sequence[int] | None = None,
) -> Iterable[Frame]:
    """Yield one or more topology-consistent frames from a stacked GRO file."""
    selected = None if raw_indexes is None else _validated_raw_indexes(raw_indexes)
    selected_set = None if selected is None else set(selected)
    iterator = _iter_gro_blocks(path)
    try:
        first = next(iterator)
    except StopIteration as exc:
        raise ValueError(f"Invalid GRO file: {path}") from exc
    signature = _gro_topology_signature(first)
    try:
        second = next(iterator)
    except StopIteration:
        if selected_set is None or 0 in selected_set:
            first.name = path.stem
            yield first
        if selected is not None and any(index != 0 for index in selected):
            raise ValueError(f"GRO frame index exceeds the single frame in {path}.")
        return

    frames = ((0, first), (1, second))
    last_index = -1
    for raw_index, frame in frames:
        _validate_gro_topology(frame, signature, path, raw_index)
        frame.name = f"{path.stem}_frame{raw_index:06d}"
        if selected_set is None or raw_index in selected_set:
            yield frame
        last_index = raw_index
    for raw_index, frame in enumerate(iterator, start=2):
        _validate_gro_topology(frame, signature, path, raw_index)
        frame.name = f"{path.stem}_frame{raw_index:06d}"
        if selected_set is None or raw_index in selected_set:
            yield frame
        last_index = raw_index
    if selected is not None and selected and selected[-1] > last_index:
        raise ValueError(
            f"GRO frame index {selected[-1]} exceeds the {last_index + 1} frames in {path}."
        )


def _iter_gro_blocks(path: Path) -> Iterable[Frame]:
    """Parse repeated title/count/atoms/box GRO blocks without loading the file."""
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        frame_index = 0
        while True:
            title_line = handle.readline()
            while title_line != "" and not title_line.strip():
                title_line = handle.readline()
            if title_line == "":
                return
            title = title_line.rstrip("\r\n")
            count_line = handle.readline()
            if count_line == "":
                if not title.strip():
                    return
                raise ValueError(
                    f"Truncated GRO file {path}: frame {frame_index} lacks an atom count."
                )
            try:
                natoms = int(count_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"Invalid GRO file {path}: atom count in frame {frame_index} is "
                    f"{count_line.strip()!r}."
                ) from exc
            if natoms < 0:
                raise ValueError(
                    f"Negative atom count in GRO file {path} frame {frame_index}."
                )
            atoms: list[Atom] = []
            for atom_index in range(natoms):
                line = handle.readline()
                if line == "":
                    raise ValueError(
                        f"Truncated GRO file {path}: frame {frame_index} declares "
                        f"{natoms} atoms but ends at atom {atom_index}."
                    )
                try:
                    atoms.append(_parse_gro_atom(atom_index, line.rstrip("\r\n")))
                except (TypeError, ValueError, IndexError) as exc:
                    raise ValueError(
                        f"Invalid GRO atom record in {path} frame {frame_index}, "
                        f"atom {atom_index + 1}: {exc}"
                    ) from exc
            box_line = handle.readline()
            if box_line == "":
                raise ValueError(
                    f"Truncated GRO file {path}: frame {frame_index} lacks a box line."
                )
            box = parse_gro_box(box_line, path)
            yield Frame(
                name=f"{path.stem}_frame{frame_index:06d}",
                atoms=atoms,
                box=box,
                time_ps=_parse_title_time_ps(title),
                source=path,
            )
            frame_index += 1


def _gro_topology_signature(frame: Frame) -> tuple[tuple[int, str, str, int], ...]:
    """Return the ordered atom identity required to stay fixed across GRO frames."""
    return tuple(
        (int(atom.resid), str(atom.resname), str(atom.atomname), int(atom.atomid))
        for atom in frame.atoms
    )


def _validate_gro_topology(
    frame: Frame,
    expected: tuple[tuple[int, str, str, int], ...],
    path: Path,
    raw_index: int,
) -> None:
    """Reject atom-count or ordered-identity changes within one stacked GRO."""
    observed = _gro_topology_signature(frame)
    if observed != expected:
        if len(observed) != len(expected):
            detail = f"atom count changed from {len(expected)} to {len(observed)}"
        else:
            mismatch = next(
                index
                for index, (left, right) in enumerate(zip(expected, observed, strict=True))
                if left != right
            )
            detail = f"ordered atom identity changed at atom {mismatch + 1}"
        raise ValueError(
            f"Stacked GRO topology mismatch in {path} frame {raw_index}: {detail}."
        )


def _validated_raw_indexes(values: Sequence[int]) -> tuple[int, ...]:
    """Normalize sorted unique nonnegative raw frame indexes."""
    indexes = tuple(int(value) for value in values)
    if any(index < 0 for index in indexes):
        raise ValueError("Raw frame indexes must be nonnegative.")
    if indexes != tuple(sorted(set(indexes))):
        raise ValueError("Raw frame indexes must be sorted and unique.")
    return indexes


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
        parts = record.split()
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
    if fixed_width and len(record) >= 68 and record[44:68].strip():
        try:
            velocity = np.asarray(
                [float(record[44:52]), float(record[52:60]), float(record[60:68])],
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
    match = re.search(r"\b(?:time_ps|t)\s*=\s*([-+0-9.eE]+)", title, re.IGNORECASE)
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


def trajectory_frame_selection(
    path: Path,
    topology: Path | None,
    delta_time_ps: float | None = None,
    lammps_config: Mapping[str, object] | None = None,
) -> TrajectorySelection:
    """Resolve selected raw indexes from physical frame times."""
    suffix = path.suffix.lower()
    if suffix == ".gro":
        times = gro_frame_times_ps(path)
    elif suffix in LAMMPS_TRAJECTORY_SUFFIXES:
        times = lammps_trajectory_times_ps(path, topology, lammps_config)
    elif suffix in {".xtc", ".trr"}:
        universe = open_mdanalysis_universe(path, topology)
        try:
            times = []
            for raw_index, ts in enumerate(universe.trajectory):
                try:
                    value = float(ts.time)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Trajectory frame {raw_index} in {path} lacks a valid time."
                    ) from exc
                times.append(value if np.isfinite(value) else None)
        finally:
            close_mdanalysis_universe(universe)
    else:
        raise ValueError(f"Time-based frame selection is not supported for {path}.")
    return select_time_frames(times, delta_time_ps, path)


def gro_frame_times_ps(path: Path) -> list[float | None]:
    """Read and validate stacked-GRO topology while collecting frame times."""
    times: list[float | None] = []
    signature: tuple[tuple[int, str, str, int], ...] | None = None
    for raw_index, frame in enumerate(_iter_gro_blocks(path)):
        if signature is None:
            signature = _gro_topology_signature(frame)
        else:
            _validate_gro_topology(frame, signature, path, raw_index)
        times.append(frame.time_ps)
    if not times:
        raise ValueError(f"Invalid GRO file: {path}")
    return times


def select_time_frames(
    times_ps: Sequence[float | None],
    delta_time_ps: float | None,
    source: Path | str,
) -> TrajectorySelection:
    """Select a regular physical-time interval without rounding to nearby frames."""
    total = len(times_ps)
    if total < 1:
        raise ValueError(f"Trajectory contains no frames: {source}")
    requested = _optional_positive_time(delta_time_ps)
    complete_times = all(value is not None and np.isfinite(float(value)) for value in times_ps)
    if requested is not None and not complete_times:
        raise ValueError(
            f"--delta-time requires valid time metadata on every frame in {source}."
        )

    native: float | None = None
    regular = False
    if complete_times and total > 1:
        values = np.asarray(times_ps, dtype=float)
        differences = np.diff(values)
        if np.any(differences <= 0.0):
            if requested is not None:
                raise ValueError(
                    f"--delta-time requires strictly increasing frame times in {source}."
                )
        else:
            native = float(differences[0])
            tolerance = max(TIME_TOLERANCE, abs(native) * TIME_TOLERANCE)
            regular = bool(
                np.allclose(
                    differences,
                    native,
                    rtol=TIME_TOLERANCE,
                    atol=tolerance,
                )
            )
            if not regular:
                native = None
                if requested is not None:
                    raise ValueError(
                        f"--delta-time requires a regular native frame interval in {source}."
                    )

    step = 1
    if requested is not None and total > 1:
        if native is None or not regular:
            raise ValueError(
                f"Cannot resolve the native frame interval required by --delta-time in {source}."
            )
        tolerance = max(TIME_TOLERANCE, abs(native) * TIME_TOLERANCE)
        if requested + tolerance < native:
            raise ValueError(
                f"--delta-time {requested:g} ps is shorter than the native "
                f"interval {native:g} ps in {source}."
            )
        ratio = requested / native
        nearest = int(round(ratio))
        if nearest < 1 or not math_isclose(ratio, nearest):
            raise ValueError(
                f"--delta-time {requested:g} ps must be an integer multiple of "
                f"the native interval {native:g} ps in {source}."
            )
        step = nearest
    return TrajectorySelection(
        raw_indexes=tuple(range(0, total, step)),
        total_frames=total,
        native_interval_ps=native,
        delta_time_ps=requested,
        raw_frame_step=step,
    )


def _optional_positive_time(value: Any) -> float | None:
    """Normalize an optional positive physical-time value in ps."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("input.delta_time_ps / -dt / --delta-time must be positive.") from exc
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError("input.delta_time_ps / -dt / --delta-time must be positive.")
    return number


def math_isclose(left: float, right: float) -> bool:
    """Compare a sampling ratio with an integer using one strict tolerance."""
    return bool(np.isclose(left, right, rtol=TIME_TOLERANCE, atol=TIME_TOLERANCE))


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


def read_mdanalysis(
    path: Path,
    topology: Path | None,
    raw_indexes: Sequence[int] | None = None,
) -> Iterable[Frame]:
    """Read selected XTC/TRR frames through MDAnalysis."""
    universe = open_mdanalysis_universe(path, topology)
    try:
        atom_metadata = trajectory_atom_metadata(universe)
        indexes = (
            tuple(range(len(universe.trajectory)))
            if raw_indexes is None
            else _validated_raw_indexes(raw_indexes)
        )
        if indexes and indexes[-1] >= len(universe.trajectory):
            raise ValueError(
                f"Trajectory frame index {indexes[-1]} exceeds the "
                f"{len(universe.trajectory)} frames in {path}."
            )
        for raw_frame_index in indexes:
            yield frame_from_mdanalysis_universe(
                universe,
                path,
                raw_frame_index,
                atom_metadata=atom_metadata,
            )
    finally:
        close_mdanalysis_universe(universe)
