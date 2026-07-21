from __future__ import annotations

"""Strict orthorhombic LAMMPS trajectory input for SQQ."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..models import Atom, Frame


LAMMPS_DUMP_SUFFIXES = frozenset({".dump", ".lammpstrj"})
LAMMPS_DCD_SUFFIXES = frozenset({".dcd"})
LAMMPS_TRAJECTORY_SUFFIXES = LAMMPS_DUMP_SUFFIXES | LAMMPS_DCD_SUFFIXES
BOX_TOLERANCE = 1.0e-8

# Native LAMMPS length/time to nm/ps, followed by MDAnalysis unit names.
_UNIT_FACTORS = {
    "real": (0.1, 0.001, "Angstrom", "fs"),
    "metal": (0.1, 1.0, "Angstrom", "ps"),
    "nano": (1.0, 1000.0, "nm", "ns"),
}
_ATOM_STYLES = {
    "full": "id resid type charge x y z",
    "molecular": "id resid type x y z",
    "bond": "id resid type x y z",
    "angle": "id resid type x y z",
}
_COORDINATE_ALIASES = {
    "auto": "auto",
    "x": "unscaled",
    "xyz": "unscaled",
    "x/y/z": "unscaled",
    "x y z": "unscaled",
    "unscaled": "unscaled",
    "xs": "scaled",
    "xsyszs": "scaled",
    "xs/ys/zs": "scaled",
    "xs ys zs": "scaled",
    "scaled": "scaled",
    "xu": "unwrapped",
    "xuyuzu": "unwrapped",
    "xu/yu/zu": "unwrapped",
    "xu yu zu": "unwrapped",
    "unwrapped": "unwrapped",
    "xsu": "scaled_unwrapped",
    "xsuysuzsu": "scaled_unwrapped",
    "xsu/ysu/zsu": "scaled_unwrapped",
    "xsu ysu zsu": "scaled_unwrapped",
    "scaled_unwrapped": "scaled_unwrapped",
    "scaled-unwrapped": "scaled_unwrapped",
}
_COORDINATE_COLUMNS = {
    "unscaled": ("x", "y", "z"),
    "scaled": ("xs", "ys", "zs"),
    "unwrapped": ("xu", "yu", "zu"),
    "scaled_unwrapped": ("xsu", "ysu", "zsu"),
}


@dataclass(frozen=True)
class LammpsTypeMapEntry:
    """SQQ names assigned to one numeric LAMMPS atom type."""

    resname: str
    atomname: str
    ignore: bool = False


@dataclass(frozen=True)
class LammpsInputConfig:
    """Normalized first-release LAMMPS reader configuration."""

    units: str
    timestep: float
    atom_style: str
    coordinate_convention: str
    type_map: Mapping[str, LammpsTypeMapEntry]
    stride: int

    @property
    def length_to_nm(self) -> float:
        return _UNIT_FACTORS[self.units][0]

    @property
    def time_to_ps(self) -> float:
        return _UNIT_FACTORS[self.units][1]

    @property
    def mda_length_unit(self) -> str:
        return _UNIT_FACTORS[self.units][2]

    @property
    def mda_time_unit(self) -> str:
        return _UNIT_FACTORS[self.units][3]


@dataclass(frozen=True)
class LammpsDumpFrameInfo:
    """Validated dump-frame header metadata."""

    step: int
    atom_count: int
    box_native: tuple[float, float, float]
    coordinate_convention: str


def normalize_lammps_config(
    config: Mapping[str, Any] | LammpsInputConfig | None,
) -> LammpsInputConfig:
    """Return a strict configuration accepted by all public reader helpers."""
    if isinstance(config, LammpsInputConfig):
        return config
    values = dict(config or {})
    units = str(values.get("units", "real")).strip().lower()
    if units not in _UNIT_FACTORS:
        allowed = ", ".join(_UNIT_FACTORS)
        raise ValueError(
            f"LAMMPS units must be one of {allowed}; got {units!r}. "
            "units lj and unknown unit styles are not supported."
        )
    return LammpsInputConfig(
        units=units,
        timestep=_positive_float(values.get("timestep", 1.0), "LAMMPS timestep"),
        atom_style=normalize_atom_style(values.get("atom_style", "full")),
        coordinate_convention=normalize_coordinate_convention(
            values.get("coordinate_convention", "auto")
        ),
        type_map=normalize_type_map(values.get("type_map", {})),
        stride=_positive_int(values.get("stride", 1), "LAMMPS stride"),
    )


def normalize_atom_style(value: Any) -> str:
    """Convert a named/custom style to MDAnalysis DATA column names."""
    text = str(value).strip().lower()
    if text in _ATOM_STYLES:
        return _ATOM_STYLES[text]
    tokens = text.replace("molecule-id", "resid").replace("molecule", "resid").split()
    tokens = ["resid" if item in {"mol", "molid", "molecule_id"} else item for item in tokens]
    tokens = ["charge" if item == "q" else item for item in tokens]
    required = {"id", "resid", "type", "x", "y", "z"}
    missing = sorted(required.difference(tokens))
    if missing:
        raise ValueError(
            "LAMMPS atom_style must expose atom id, molecule id, type, and x/y/z; "
            f"missing {', '.join(missing)}. First-release SQQ does not support "
            "atomic/charge styles without molecule IDs."
        )
    if len(tokens) != len(set(tokens)):
        raise ValueError("LAMMPS atom_style contains duplicate columns.")
    return " ".join(tokens)


def normalize_coordinate_convention(value: Any) -> str:
    """Normalize x/xs/xu/xsu aliases to MDAnalysis names."""
    text = str(value).strip().lower()
    try:
        return _COORDINATE_ALIASES[text]
    except KeyError as exc:
        raise ValueError(
            "LAMMPS coordinate_convention must be auto, x, xs, xu, xsu, "
            "unscaled, scaled, unwrapped, or scaled_unwrapped."
        ) from exc


def normalize_type_map(value: Any) -> dict[str, LammpsTypeMapEntry]:
    """Normalize numeric atom-type mappings and explicit ignore entries."""
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError("LAMMPS type_map must be a mapping keyed by numeric atom type.")
    result: dict[str, LammpsTypeMapEntry] = {}
    for raw_type, raw_entry in value.items():
        type_id = _numeric_type_key(raw_type)
        if type_id in result:
            raise ValueError(f"Duplicate LAMMPS type_map entry for atom type {type_id}.")
        result[type_id] = _normalize_type_entry(type_id, raw_entry)
    return result


def _normalize_type_entry(type_id: str, value: Any) -> LammpsTypeMapEntry:
    if value is None or (
        isinstance(value, str) and value.strip().lower() in {"ignore", "ignored", "skip"}
    ):
        return _ignored_type_entry(type_id)
    if isinstance(value, Mapping):
        if bool(value.get("ignore", False)):
            return _ignored_type_entry(type_id)
        resname = value.get("resname")
        atomname = value.get("atomname", value.get("name"))
    elif isinstance(value, str):
        parts = [part.strip() for part in value.replace("/", ":").split(":")]
        if len(parts) != 2:
            raise ValueError(
                f"LAMMPS type_map[{type_id}] string must be 'RESNAME:ATOMNAME' or 'ignore'."
            )
        resname, atomname = parts
    elif (
        isinstance(value, Sequence)
        and not isinstance(value, (bytes, bytearray))
        and len(value) == 2
    ):
        resname, atomname = value
    else:
        raise ValueError(
            f"LAMMPS type_map[{type_id}] must provide resname/atomname or explicit ignore."
        )
    return LammpsTypeMapEntry(
        resname=_valid_name(resname, f"LAMMPS type_map[{type_id}].resname"),
        atomname=_valid_name(atomname, f"LAMMPS type_map[{type_id}].atomname"),
    )


def _ignored_type_entry(type_id: str) -> LammpsTypeMapEntry:
    # Reserved names keep the atom in Frame/GRO but outside default selections.
    return LammpsTypeMapEntry(resname="IGN", atomname=f"T{type_id}", ignore=True)


def validate_lammps_data_box(path: Path, config: LammpsInputConfig) -> np.ndarray:
    """Validate an orthorhombic DATA box and return its lengths in nm."""
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"LAMMPS DATA topology does not exist: {path}")
    bounds: dict[str, tuple[float, float]] = {}
    tilt: tuple[float, float, float] | None = None
    general_triclinic = False
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 4 and tuple(parts[-2:]) in {
                ("xlo", "xhi"),
                ("ylo", "yhi"),
                ("zlo", "zhi"),
            }:
                axis = parts[-2][0]
                bounds[axis] = (
                    _finite_float(parts[0], f"LAMMPS DATA {axis}lo"),
                    _finite_float(parts[1], f"LAMMPS DATA {axis}hi"),
                )
            elif len(parts) >= 6 and parts[-3:] == ["xy", "xz", "yz"]:
                tilt = tuple(
                    _finite_float(parts[index], "LAMMPS DATA tilt") for index in range(3)
                )
            elif any(item in {"avec", "bvec", "cvec", "abc"} for item in parts):
                general_triclinic = True
    if general_triclinic:
        raise ValueError(f"General triclinic LAMMPS DATA boxes are not supported: {path}")
    if tilt is not None and any(abs(item) > BOX_TOLERANCE for item in tilt):
        raise ValueError(f"Tilted LAMMPS DATA boxes are not supported: {path} has tilt {tilt}.")
    missing = [axis for axis in "xyz" if axis not in bounds]
    if missing:
        raise ValueError(f"LAMMPS DATA topology lacks {'/'.join(missing)} box bounds: {path}")
    lengths = np.asarray(
        [bounds[axis][1] - bounds[axis][0] for axis in "xyz"], dtype=float
    )
    _validate_box_lengths(lengths, f"LAMMPS DATA topology {path}")
    return lengths * config.length_to_nm


def inspect_lammps_dump(
    path: Path,
    config: LammpsInputConfig,
    topology_atom_ids: Sequence[int] | None = None,
) -> tuple[LammpsDumpFrameInfo, ...]:
    """Validate every dump box, boundary, coordinate schema, and atom-ID set."""
    expected_ids = (
        None
        if topology_atom_ids is None
        else tuple(sorted(int(item) for item in topology_atom_ids))
    )
    frames: list[LammpsDumpFrameInfo] = []
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        frame_index = 0
        while True:
            header = handle.readline()
            if not header:
                break
            if not header.strip():
                continue
            _expect_header(header, "ITEM: TIMESTEP", path, frame_index)
            step_line = _required_line(handle, path, frame_index, "timestep")
            try:
                step = int(step_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"Invalid LAMMPS timestep in {path} frame {frame_index}: {step_line!r}"
                ) from exc

            _expect_header(
                _required_line(handle, path, frame_index, "atom-count header"),
                "ITEM: NUMBER OF ATOMS",
                path,
                frame_index,
            )
            count_line = _required_line(handle, path, frame_index, "atom count")
            try:
                atom_count = int(count_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"Invalid LAMMPS atom count in {path} frame {frame_index}: {count_line!r}"
                ) from exc
            if atom_count < 1:
                raise ValueError(
                    f"LAMMPS dump atom count must be positive in {path} frame {frame_index}."
                )
            if expected_ids is not None and atom_count != len(expected_ids):
                raise ValueError(
                    f"LAMMPS dump/topology atom-count mismatch in {path} frame {frame_index}: "
                    f"{atom_count} != {len(expected_ids)}."
                )

            box_header = _required_line(handle, path, frame_index, "box header").split()
            if box_header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
                raise ValueError(f"Expected ITEM: BOX BOUNDS in {path} frame {frame_index}.")
            box_tokens = box_header[3:]
            if any(
                item in {"xy", "xz", "yz", "abc", "avec", "bvec", "cvec"}
                for item in box_tokens
            ):
                raise ValueError(
                    f"Triclinic/tilted LAMMPS dump boxes are not supported: "
                    f"{path} frame {frame_index}."
                )
            if box_tokens != ["pp", "pp", "pp"]:
                actual = " ".join(box_tokens) or "no boundary flags"
                raise ValueError(
                    "LAMMPS dump requires fully periodic boundaries 'pp pp pp': "
                    f"{path} frame {frame_index} has {actual}."
                )
            lengths: list[float] = []
            for axis in "xyz":
                fields = _required_line(handle, path, frame_index, f"{axis} bounds").split()
                if len(fields) != 2:
                    raise ValueError(
                        f"Orthorhombic LAMMPS {axis} bounds need two values in "
                        f"{path} frame {frame_index}."
                    )
                lower = _finite_float(fields[0], f"LAMMPS {axis}lo")
                upper = _finite_float(fields[1], f"LAMMPS {axis}hi")
                lengths.append(upper - lower)
            _validate_box_lengths(
                np.asarray(lengths), f"LAMMPS dump {path} frame {frame_index}"
            )

            atom_header = _required_line(handle, path, frame_index, "atom header").split()
            if atom_header[:2] != ["ITEM:", "ATOMS"]:
                raise ValueError(f"Expected ITEM: ATOMS in {path} frame {frame_index}.")
            columns = atom_header[2:]
            if "id" not in columns:
                raise ValueError(
                    f"LAMMPS dump must contain an atom id column: {path} frame {frame_index}."
                )
            convention = _dump_coordinate_convention(
                columns, config.coordinate_convention, path, frame_index
            )
            id_column = columns.index("id")
            ids: list[int] = []
            for atom_row in range(atom_count):
                fields = _required_line(
                    handle, path, frame_index, f"atom row {atom_row + 1}"
                ).split()
                if len(fields) < len(columns):
                    raise ValueError(
                        f"Truncated LAMMPS atom row in {path} frame {frame_index}."
                    )
                try:
                    atom_id = int(fields[id_column])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid LAMMPS atom id in {path} frame {frame_index}."
                    ) from exc
                if atom_id < 1:
                    raise ValueError(
                        f"LAMMPS atom IDs must be positive in {path} frame {frame_index}."
                    )
                ids.append(atom_id)
            ordered_ids = tuple(sorted(ids))
            if len(set(ordered_ids)) != atom_count:
                raise ValueError(f"Duplicate atom IDs in {path} frame {frame_index}.")
            if expected_ids is not None and ordered_ids != expected_ids:
                raise ValueError(
                    f"LAMMPS dump atom IDs do not match DATA topology in "
                    f"{path} frame {frame_index}."
                )
            frames.append(
                LammpsDumpFrameInfo(
                    step=step,
                    atom_count=atom_count,
                    box_native=tuple(lengths),
                    coordinate_convention=convention,
                )
            )
            frame_index += 1
    if not frames:
        raise ValueError(f"LAMMPS dump contains no frames: {path}")
    conventions = {item.coordinate_convention for item in frames}
    if len(conventions) != 1:
        raise ValueError(
            f"LAMMPS coordinate columns change between frames in {path}: "
            f"{sorted(conventions)}"
        )
    return tuple(frames)


def open_lammps_universe(
    path: Path,
    topology: Path | None,
    config: Mapping[str, Any] | LammpsInputConfig | None,
):
    """Open a LAMMPS DATA + dump/DCD trajectory through MDAnalysis."""
    settings = normalize_lammps_config(config)
    topology_path = _required_topology(topology)
    validate_lammps_data_box(topology_path, settings)
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError("Reading LAMMPS trajectories requires MDAnalysis.") from exc

    common = {"topology_format": "DATA", "atom_style": settings.atom_style}
    suffix = path.suffix.lower()
    if suffix in LAMMPS_DUMP_SUFFIXES:
        universe = mda.Universe(
            str(topology_path),
            str(path),
            format="LAMMPSDUMP",
            lammps_coordinate_convention=settings.coordinate_convention,
            **common,
        )
    elif suffix in LAMMPS_DCD_SUFFIXES:
        universe = mda.Universe(
            str(topology_path),
            str(path),
            format="LAMMPS",
            lengthunit=settings.mda_length_unit,
            timeunit=settings.mda_time_unit,
            **common,
        )
    else:
        raise ValueError(f"Unsupported LAMMPS trajectory format: {path}")
    try:
        lammps_atom_metadata(universe, settings)
    except Exception:
        close_lammps_universe(universe)
        raise
    return universe


def close_lammps_universe(universe) -> None:
    """Close the MDAnalysis reader when it owns an open trajectory handle."""
    trajectory = getattr(universe, "trajectory", None)
    close = getattr(trajectory, "close", None)
    if callable(close):
        close()


def lammps_atom_metadata(
    universe,
    config: Mapping[str, Any] | LammpsInputConfig,
) -> tuple[tuple[int, int, str, str, int], ...]:
    """Build immutable SQQ metadata in sorted LAMMPS atom-ID order."""
    settings = normalize_lammps_config(config)
    if not settings.type_map:
        raise ValueError(
            "LAMMPS input requires a non-empty type_map; every numeric atom type "
            "must map to resname/atomname or explicit ignore."
        )
    raw: list[tuple[int, int, str]] = []
    for atom in universe.atoms:
        raw.append((int(atom.id), int(atom.resid), _numeric_type_key(atom.type)))
    raw.sort(key=lambda item: item[0])
    atom_ids = [item[0] for item in raw]
    if any(item < 1 for item in atom_ids) or len(atom_ids) != len(set(atom_ids)):
        raise ValueError("LAMMPS DATA atom IDs must be positive and unique.")
    missing = sorted(
        {type_id for _, _, type_id in raw if type_id not in settings.type_map}, key=int
    )
    if missing:
        raise ValueError(f"LAMMPS type_map is missing atom type(s): {', '.join(missing)}.")
    molecule_resnames: dict[int, set[str]] = {}
    for _, resid, type_id in raw:
        entry = settings.type_map[type_id]
        if entry.ignore:
            continue
        if resid < 1:
            raise ValueError("Mapped LAMMPS atoms require positive molecule IDs.")
        molecule_resnames.setdefault(resid, set()).add(entry.resname)
    conflicts = sorted(str(resid) for resid, names in molecule_resnames.items() if len(names) > 1)
    if conflicts:
        raise ValueError("LAMMPS molecule IDs map to multiple residue names: " + ", ".join(conflicts))

    return tuple(
        (
            index,
            resid,
            settings.type_map[type_id].resname,
            settings.type_map[type_id].atomname,
            atom_id,
        )
        for index, (atom_id, resid, type_id) in enumerate(raw)
    )


def frame_from_lammps_universe(
    universe,
    path: Path,
    raw_frame_index: int,
    config: Mapping[str, Any] | LammpsInputConfig,
    atom_metadata: tuple[tuple[int, int, str, str, int], ...] | None = None,
) -> Frame:
    """Materialize one LAMMPS frame in stable atom order, nm, and ps."""
    settings = normalize_lammps_config(config)
    ts = universe.trajectory[int(raw_frame_index)]
    suffix = path.suffix.lower()
    if suffix in LAMMPS_DUMP_SUFFIXES:
        position_factor = settings.length_to_nm
    elif suffix in LAMMPS_DCD_SUFFIXES:
        # The configured LAMMPS DCD reader converts positions to Angstrom.
        position_factor = 0.1
    else:
        raise ValueError(f"Unsupported LAMMPS trajectory format: {path}")
    atom_ids = np.asarray(universe.atoms.ids, dtype=int)
    if np.any(atom_ids < 1) or len(atom_ids) != len(set(atom_ids.tolist())):
        raise ValueError(
            f"LAMMPS trajectory atom IDs must be positive and unique: {path}."
        )
    atom_order = np.argsort(atom_ids, kind="stable")
    positions_nm = np.asarray(universe.atoms.positions, dtype=float)[atom_order] * position_factor
    if np.any(~np.isfinite(positions_nm)):
        raise ValueError(
            f"Non-finite LAMMPS coordinates in {path} frame {raw_frame_index}."
        )
    metadata = atom_metadata or lammps_atom_metadata(universe, settings)
    if len(metadata) != len(positions_nm):
        raise ValueError(f"LAMMPS topology does not match coordinates in {path}.")
    atoms = [
        Atom(
            index=index,
            resid=resid,
            resname=resname,
            atomname=atomname,
            atomid=atomid,
            xyz=np.asarray(xyz, dtype=float),
            molecule_id=resid,
        )
        for (index, resid, resname, atomname, atomid), xyz in zip(
            metadata, positions_nm, strict=True
        )
    ]

    dimensions = getattr(ts, "dimensions", None)
    if dimensions is None or len(dimensions) < 6:
        raise ValueError(
            f"LAMMPS trajectory frame lacks a periodic cell: {path} frame {raw_frame_index}."
        )
    dimensions = np.asarray(dimensions, dtype=float)
    box_factor = settings.length_to_nm if suffix in LAMMPS_DUMP_SUFFIXES else 0.1
    box = dimensions[:3] * box_factor
    _validate_box_lengths(box, f"LAMMPS trajectory {path} frame {raw_frame_index}")
    angles = dimensions[3:6]
    if np.any(~np.isfinite(angles)) or not np.allclose(
        angles, 90.0, atol=1.0e-5, rtol=0.0
    ):
        raise ValueError(
            f"Triclinic LAMMPS boxes are not supported: {path} frame "
            f"{raw_frame_index} has angles {angles.tolist()}."
        )

    if suffix in LAMMPS_DUMP_SUFFIXES:
        try:
            step = int(ts.data["step"])
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"LAMMPS dump frame lacks a valid timestep: {path} frame {raw_frame_index}."
            ) from exc
        time_ps = step * settings.timestep * settings.time_to_ps
    else:
        try:
            time_ps = float(ts.time)
        except (TypeError, ValueError):
            time_ps = raw_frame_index * settings.timestep * settings.time_to_ps
    if not np.isfinite(time_ps):
        raise ValueError(f"Non-finite LAMMPS frame time in {path} frame {raw_frame_index}.")
    return Frame(
        name=f"{path.stem}_frame{int(ts.frame):06d}",
        atoms=atoms,
        box=np.asarray(box, dtype=float),
        time_ps=float(time_ps),
        source=path,
    )


def lammps_trajectory_frame_indices(
    path: Path,
    topology: Path | None,
    config: Mapping[str, Any] | LammpsInputConfig | None,
) -> list[int]:
    """Return raw frame indexes selected by the normalized stride."""
    settings = normalize_lammps_config(config)
    universe = open_lammps_universe(path, topology, settings)
    try:
        metadata = lammps_atom_metadata(universe, settings)
        if path.suffix.lower() in LAMMPS_DUMP_SUFFIXES:
            frames = inspect_lammps_dump(
                path, settings, [item[4] for item in metadata]
            )
            if len(frames) != len(universe.trajectory):
                raise ValueError(
                    f"LAMMPS dump frame-count mismatch between validation and MDAnalysis: {path}."
                )
        return list(range(0, len(universe.trajectory), settings.stride))
    finally:
        close_lammps_universe(universe)


def read_lammps(
    path: Path,
    topology: Path | None,
    config: Mapping[str, Any] | LammpsInputConfig | None,
) -> Iterable[Frame]:
    """Yield validated LAMMPS frames through the common SQQ model."""
    settings = normalize_lammps_config(config)
    universe = open_lammps_universe(path, topology, settings)
    try:
        metadata = lammps_atom_metadata(universe, settings)
        if path.suffix.lower() in LAMMPS_DUMP_SUFFIXES:
            frames = inspect_lammps_dump(
                path, settings, [item[4] for item in metadata]
            )
            if len(frames) != len(universe.trajectory):
                raise ValueError(
                    f"LAMMPS dump frame-count mismatch between validation and MDAnalysis: {path}."
                )
        for raw_index in range(0, len(universe.trajectory), settings.stride):
            yield frame_from_lammps_universe(
                universe,
                path,
                raw_index,
                settings,
                atom_metadata=metadata,
            )
    finally:
        close_lammps_universe(universe)


def _dump_coordinate_convention(
    columns: Sequence[str], requested: str, path: Path, frame_index: int
) -> str:
    available = [
        convention
        for convention, required in _COORDINATE_COLUMNS.items()
        if all(item in columns for item in required)
    ]
    if requested == "auto":
        if not available:
            raise ValueError(
                f"No x/xs/xu/xsu coordinate triplet in {path} frame {frame_index}."
            )
        return available[0]
    if requested not in available:
        required = "/".join(_COORDINATE_COLUMNS[requested])
        raise ValueError(
            f"Requested LAMMPS coordinates {required} are absent in "
            f"{path} frame {frame_index}."
        )
    return requested


def _required_topology(topology: Path | None) -> Path:
    if topology is None:
        raise ValueError(
            "LAMMPS trajectory input requires a DATA topology via --top system.data."
        )
    return Path(topology)


def _expect_header(line: str, expected: str, path: Path, frame_index: int) -> None:
    if line.strip() != expected:
        raise ValueError(
            f"Expected {expected!r} in {path} frame {frame_index}; got {line.strip()!r}."
        )


def _required_line(handle, path: Path, frame_index: int, label: str) -> str:
    line = handle.readline()
    if not line:
        raise ValueError(
            f"Truncated LAMMPS dump {path} in frame {frame_index} ({label})."
        )
    return line


def _validate_box_lengths(lengths: np.ndarray, label: str) -> None:
    if len(lengths) != 3 or np.any(~np.isfinite(lengths)) or np.any(lengths <= 0):
        raise ValueError(
            f"{label} requires three positive finite orthorhombic box lengths; "
            f"got {lengths.tolist()}."
        )


def _numeric_type_key(value: Any) -> str:
    text = str(value).strip()
    try:
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"LAMMPS atom type must be numeric; got {value!r}.") from exc
    if number < 1 or text not in {str(number), f"+{number}"}:
        raise ValueError(f"LAMMPS atom type must be a positive integer; got {value!r}.")
    return str(number)


def _valid_name(value: Any, label: str) -> str:
    if value is None:
        raise ValueError(f"{label} is required.")
    text = str(value).strip()
    if not text or any(item.isspace() for item in text):
        raise ValueError(f"{label} must be a non-empty name without whitespace.")
    return text


def _positive_float(value: Any, label: str) -> float:
    number = _finite_float(value, label)
    if number <= 0:
        raise ValueError(f"{label} must be positive.")
    return number


def _finite_float(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number.") from exc
    if not np.isfinite(number):
        raise ValueError(f"{label} must be a finite number.")
    return number


def _positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer.")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer.") from exc
    if str(value).strip() not in {str(number), f"+{number}"} or number < 1:
        raise ValueError(f"{label} must be a positive integer.")
    return number
