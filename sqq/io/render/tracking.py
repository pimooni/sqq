"""Discovery and publication of render packages for Track targets."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from typing import Mapping, Sequence
import warnings
from uuid import uuid4

from ...models.tracking import TargetSelection, TrackingResult
from ..tracking import TRACK_DIRECTORY_NAME, rewrite_membership_for_targets
from .inspect import inspect_render_script
from .models import (
    SQQ_CAGE_GRO_NAME,
    SQQ_CAGE_MEMBERSHIP_NAME,
    SQQ_CAGE_XTC_NAME,
    SQQ_RENDER_DIRECTORY,
    SQQ_RENDER_SCRIPT_NAME,
    RenderBundle,
    RenderNames,
)
from .tcl import vmd_script_text


TRACK_RENDER_NAMES = RenderNames(
    directory=SQQ_RENDER_DIRECTORY,
    topology="sqq_track.gro",
    trajectory="sqq_track.xtc",
    membership="sqq_track.membership.tsv",
    script="sqq_track.vmd.tcl",
)
TRACK_RENDER_DIRECTORY = TRACK_RENDER_NAMES.directory
TRACK_GRO_NAME = TRACK_RENDER_NAMES.topology
TRACK_XTC_NAME = TRACK_RENDER_NAMES.trajectory
TRACK_MEMBERSHIP_NAME = TRACK_RENDER_NAMES.membership
TRACK_TCL_NAME = TRACK_RENDER_NAMES.script

_SOURCE_PROVENANCE_FORMAT = "SQQ render source provenance"
_SOURCE_PROVENANCE_VERSION = 1


def discover_sqq_cage_bundle(
    source: str | Path | None = None,
    *,
    state_path: str | Path | None = None,
) -> RenderBundle:
    """Find one complete Analyze render package associated with Track input."""
    roots: list[Path] = []
    if source is not None:
        root = Path(source)
        if root.is_file():
            root = root.parent
        roots.extend((root / SQQ_RENDER_DIRECTORY, root))
    if state_path is not None:
        state = Path(state_path)
        roots.extend(
            (
                state.parent / SQQ_RENDER_DIRECTORY,
                state.parent.parent / SQQ_RENDER_DIRECTORY,
            )
        )
    if source is None and state_path is None:
        roots.extend(
            (
                Path.cwd() / SQQ_RENDER_DIRECTORY,
                Path.cwd() / TRACK_DIRECTORY_NAME / SQQ_RENDER_DIRECTORY,
            )
        )

    bundles: list[RenderBundle] = []
    seen: set[str] = set()
    for candidate in roots:
        render_dir = candidate.resolve()
        identity = os.path.normcase(str(render_dir))
        if identity in seen:
            continue
        seen.add(identity)
        gro = render_dir / SQQ_CAGE_GRO_NAME
        xtc = render_dir / SQQ_CAGE_XTC_NAME
        membership = render_dir / SQQ_CAGE_MEMBERSHIP_NAME
        script = render_dir / SQQ_RENDER_SCRIPT_NAME
        if all(path.is_file() for path in (gro, xtc, membership, script)):
            bundles.append(
                RenderBundle(
                    gro_path=gro,
                    script_path=script,
                    frame_count=_bundle_frame_count(gro, membership),
                    xtc_path=xtc,
                    membership_path=membership,
                    render_dir=render_dir,
                )
            )
    if len(bundles) == 1:
        return bundles[0]
    if len(bundles) > 1:
        raise ValueError(
            "Multiple SQQ render bundles were found: "
            + ", ".join(str(bundle.render_dir) for bundle in bundles)
        )
    if _looks_like_track_result(source, state_path):
        raise FileNotFoundError(
            "The --source directory is a Track result: it contains target "
            f"{TRACK_TCL_NAME} packages but no Analyze "
            f"{SQQ_RENDER_DIRECTORY}/{SQQ_CAGE_GRO_NAME} package. Track results "
            "cannot be imported again; use the Analyze result directory that "
            "produced them."
        )
    raise FileNotFoundError(
        f"Cannot find {SQQ_RENDER_DIRECTORY}/{{{SQQ_CAGE_GRO_NAME}, "
        f"{SQQ_CAGE_XTC_NAME}, {SQQ_CAGE_MEMBERSHIP_NAME}, "
        f"{SQQ_RENDER_SCRIPT_NAME}}}."
    )


def _looks_like_track_result(
    source: str | Path | None,
    state_path: str | Path | None,
) -> bool:
    """Detect a Track output root, whose render packages are per-target."""
    track_roots: list[Path] = []
    if state_path is not None:
        track_roots.append(Path(state_path).parent)
    if source is not None:
        root = Path(source)
        if root.is_file():
            root = root.parent
        track_roots.extend((root / TRACK_DIRECTORY_NAME, root))
    for track_root in track_roots:
        if not track_root.is_dir():
            continue
        if any(
            (child / SQQ_RENDER_DIRECTORY / TRACK_TCL_NAME).is_file()
            for child in track_root.iterdir()
            if child.is_dir()
        ):
            return True
    return False


def discover_sqq_cage_gro(
    source: str | Path | None = None,
    *,
    state_path: str | Path | None = None,
) -> Path:
    """Return the topology GRO from one discovered Analyze render package."""
    bundle = discover_sqq_cage_bundle(source, state_path=state_path)
    if bundle.gro_path is None:
        raise FileNotFoundError("The SQQ render bundle has no topology GRO.")
    return bundle.gro_path


def validate_tracking_source_bundle(
    bundle: RenderBundle,
    *,
    frame_count: int | None = None,
    result: TrackingResult | None = None,
) -> tuple[Path, Path, Path]:
    """Validate a render package against its Track state and saved provenance."""
    expected_frames = len(result.frames) if result is not None else frame_count
    if expected_frames is None:
        raise ValueError("Tracking source validation requires a frame count or state.")
    if bundle.frame_count != int(expected_frames):
        raise ValueError(
            "Tracking state and render bundle frame counts differ: "
            f"{expected_frames} versus {bundle.frame_count}."
        )
    paths, current = _source_bundle_facts(bundle, result=result)
    if result is not None and result.source_provenance:
        _validate_source_provenance(result.source_provenance, current)
    return paths


def build_tracking_source_provenance(
    bundle: RenderBundle,
    result: TrackingResult,
) -> dict[str, object]:
    """Build the reproducible identity record stored in ``track_state.json``."""
    _, provenance = _source_bundle_facts(bundle, result=result)
    return provenance


def publish_target_render_bundle(
    selection: TargetSelection,
    target_directory: str | Path,
    source_bundle: RenderBundle,
) -> RenderBundle:
    """Publish one self-contained four-file render package for a Track target."""
    return publish_target_render_bundles(
        ((selection, target_directory),), source_bundle
    )[0]


def publish_target_render_bundles(
    targets: Sequence[tuple[TargetSelection, str | Path]],
    source_bundle: RenderBundle,
) -> list[RenderBundle]:
    """Publish the render packages of several targets from one source read.

    The topology GRO and XTC are hard-linked (or copied) per target; the
    membership TSV is read once and rewritten for every target in the same
    pass; each target then receives its own default-view Tcl script.
    """
    if not targets:
        return []
    gro_source, xtc_source, membership_source = _required_render_paths(source_bundle)
    layouts: list[tuple[TargetSelection, Path, Path, Path, Path, Path, Path]] = []
    for selection, target_directory in targets:
        target_root = Path(target_directory)
        render_dir, gro, xtc, membership, script = TRACK_RENDER_NAMES.paths(target_root)
        render_dir.mkdir(parents=True, exist_ok=True)
        _atomic_link_or_copy(gro_source, gro)
        _atomic_link_or_copy(xtc_source, xtc)
        layouts.append((selection, target_root, render_dir, gro, xtc, membership, script))
    rewrite_membership_for_targets(
        membership_source,
        [(membership, selection) for selection, _root, _dir, _gro, _xtc, membership, _s in layouts],
    )
    bundles: list[RenderBundle] = []
    for selection, target_root, render_dir, gro, xtc, membership, script in layouts:
        _atomic_write_text(
            script,
            _target_vmd_script(selection, target_root.name),
            encoding="ascii",
        )
        # Remove obsolete root-level files from pre-package Track layouts.
        (target_root / TRACK_GRO_NAME).unlink(missing_ok=True)
        (target_root / TRACK_TCL_NAME).unlink(missing_ok=True)
        bundles.append(
            RenderBundle(
                gro_path=gro,
                script_path=script,
                frame_count=source_bundle.frame_count,
                xtc_path=xtc,
                membership_path=membership,
                render_dir=render_dir,
            )
        )
    return bundles


def _target_vmd_script(selection: TargetSelection, target_name: str) -> str:
    script = vmd_script_text(
        gro_filename=TRACK_GRO_NAME,
        xtc_filename=TRACK_XTC_NAME,
        membership_filename=TRACK_MEMBERSHIP_NAME,
        molecule_name=f"SQQ track {target_name}",
        render_kind="track",
    ).rstrip()
    commands = "\n".join(
        "    " + line for line in _selection_commands(selection, target_name)
    )
    return (
        script
        + "\n\n# Default lifecycle view for this tracking target.\n"
        + "::SQQ::when_ready {\n"
        + commands
        + "\n}\n"
    )


def _selection_commands(
    selection: TargetSelection,
    target_name: str,
) -> list[str]:
    if selection.target.kind == "all":
        return ["sqq show cage all"]
    object_ids = sorted(
        {track.track_id for track in selection.tracks},
        key=lambda value: int(value[1:]),
    )
    if not object_ids:
        return [
            "set ::SQQ::active_families {}",
            "set ::SQQ::custom_show_active 1",
            "::SQQ::render_current",
            f'puts "SQQ track target {target_name} matched no cages."',
        ]
    return [
        "set ::SQQ::default_track_targets {" + " ".join(object_ids) + "}",
        "::SQQ::set_show [linsert $::SQQ::default_track_targets 0 cage]",
        "unset ::SQQ::default_track_targets",
    ]


def _required_render_paths(bundle: RenderBundle) -> tuple[Path, Path, Path]:
    paths = (bundle.gro_path, bundle.xtc_path, bundle.membership_path)
    labels = ("topology GRO", "XTC trajectory", "membership TSV")
    missing = [
        label
        for label, path in zip(labels, paths)
        if path is None or not Path(path).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete SQQ render bundle; missing " + ", ".join(missing) + "."
        )
    resolved = tuple(Path(path).resolve() for path in paths)
    empty = [
        label
        for label, path in zip(labels, resolved)
        if path.stat().st_size <= 0
    ]
    if empty:
        raise ValueError(
            "Incomplete SQQ render bundle; empty " + ", ".join(empty) + "."
        )
    return resolved  # type: ignore[return-value]


def _source_bundle_facts(
    bundle: RenderBundle,
    *,
    result: TrackingResult | None,
) -> tuple[tuple[Path, Path, Path], dict[str, object]]:
    gro, xtc, membership = _required_render_paths(bundle)
    script = _required_script_path(bundle)
    _validate_script_references(script, gro, xtc, membership)
    # A bundle produced by this process carries digests, atom count, and
    # topology identity from the writer; a discovered package is re-read and
    # re-hashed in full because nothing about it can be trusted.
    fresh = (
        bundle.file_digests is not None
        and bundle.atom_count is not None
        and bundle.topology_identity is not None
        and set(bundle.file_digests) >= {"topology", "trajectory", "membership", "script"}
    )
    if fresh:
        atom_count = int(bundle.atom_count)  # type: ignore[arg-type]
        topology_identity = str(bundle.topology_identity)
        xtc_frames = bundle.frame_count
        topology_digest: str | None = None
    else:
        atom_count, topology_identity, topology_digest = _gro_identity(gro)
        xtc_atoms, xtc_frames = _xtc_shape(xtc)
        if xtc_atoms != atom_count:
            raise ValueError(
                "SQQ render topology and trajectory atom counts differ: "
                f"{atom_count} versus {xtc_atoms}."
            )
    membership_facts = _membership_facts(
        membership,
        atom_count=atom_count,
        result=result,
    )
    membership_frames = int(membership_facts["frame_count"])
    if membership_frames != xtc_frames or bundle.frame_count != membership_frames:
        raise ValueError(
            "SQQ render trajectory, membership, and bundle frame counts differ: "
            f"{xtc_frames}, {membership_frames}, and {bundle.frame_count}."
        )
    if fresh:
        files = {
            role: dict(bundle.file_digests[role])  # type: ignore[index]
            for role in ("topology", "trajectory", "membership", "script")
        }
    else:
        files = {
            "topology": _file_identity(gro, digest=topology_digest),
            "trajectory": _file_identity(xtc),
            "membership": _file_identity(
                membership, digest=str(membership_facts["sha256"])
            ),
            "script": _file_identity(script),
        }
    provenance: dict[str, object] = {
        "format": _SOURCE_PROVENANCE_FORMAT,
        "version": _SOURCE_PROVENANCE_VERSION,
        "atom_count": atom_count,
        "frame_count": membership_frames,
        "files": files,
        "topology_identity_sha256": topology_identity,
        "component_signature_sha256": membership_facts[
            "component_signature_sha256"
        ],
        "frame_mapping_sha256": membership_facts["frame_mapping_sha256"],
        "cage_mapping_sha256": membership_facts["cage_mapping_sha256"],
    }
    if result is not None:
        provenance["tracking_config_sha256"] = _json_digest(
            result.config.to_dict()
        )
    return (gro, xtc, membership), provenance


def _required_script_path(bundle: RenderBundle) -> Path:
    if bundle.script_path is None or not Path(bundle.script_path).is_file():
        raise FileNotFoundError("Incomplete SQQ render bundle; missing VMD Tcl script.")
    script = Path(bundle.script_path).resolve()
    if script.stat().st_size <= 0:
        raise ValueError("Incomplete SQQ render bundle; empty VMD Tcl script.")
    return script


def _validate_script_references(
    script: Path,
    gro: Path,
    xtc: Path,
    membership: Path,
) -> None:
    inspection = inspect_render_script(script)
    if inspection is None or not inspection.complete:
        raise ValueError(f"Invalid or incomplete SQQ VMD script: {script}")
    references = {
        reference.role: (script.parent / reference.path).resolve()
        for reference in inspection.references
    }
    expected = {
        "topology": gro,
        "trajectory": xtc,
        "membership": membership,
    }
    for role, path in expected.items():
        if references.get(role) != path:
            raise ValueError(
                f"SQQ VMD script {role} reference does not match its render bundle."
            )


def _gro_identity(path: Path) -> tuple[int, str, str]:
    """Return atom count, atom-identity digest, and the file SHA-256 in one read."""
    try:
        with path.open("rb") as handle:
            reader = _HashingLineReader(handle)
            if not reader.readline():
                raise ValueError("missing title")
            atom_count = int(reader.readline().strip())
            if atom_count < 1:
                raise ValueError("invalid atom count")
            identity = hashlib.sha256()
            for _ in range(atom_count):
                line = reader.readline().rstrip("\r\n")
                if len(line) < 20:
                    raise ValueError("truncated atom record")
                identity.update(line[:20].encode("utf-8"))
                identity.update(b"\n")
            if not reader.readline():
                raise ValueError("missing box")
            reader.consume_rest()
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"Invalid or truncated SQQ topology GRO: {path}") from exc
    return atom_count, identity.hexdigest(), reader.digest.hexdigest()


class _HashingLineReader:
    """Iterate decoded lines of a binary file while hashing the raw bytes."""

    def __init__(self, handle: object, encoding: str = "utf-8") -> None:
        self._handle = handle
        self._encoding = encoding
        self.digest = hashlib.sha256()

    def readline(self) -> str:
        raw = self._handle.readline()  # type: ignore[attr-defined]
        self.digest.update(raw)
        return raw.decode(self._encoding)

    def __iter__(self):
        for raw in self._handle:  # type: ignore[attr-defined]
            self.digest.update(raw)
            yield raw.decode(self._encoding)

    def consume_rest(self) -> None:
        for block in iter(lambda: self._handle.read(1024 * 1024), b""):  # type: ignore[attr-defined]
            self.digest.update(block)


def _xtc_shape(path: Path) -> tuple[int, int]:
    # The low-level XDR handle computes frame offsets in memory.  The
    # high-level XTCReader would persist an ``.<name>_offsets.npz`` cache
    # beside the trajectory, which must never appear inside a published
    # render package or a read-only ``--source`` directory.
    try:
        from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    except ImportError as exc:
        raise RuntimeError("Validating an SQQ XTC requires MDAnalysis.") from exc
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with XTCFile(str(path)) as handle:
                return int(handle.n_atoms), int(len(handle))
    except Exception as exc:
        raise ValueError(f"Invalid SQQ XTC trajectory: {path}") from exc


def _membership_facts(
    path: Path,
    *,
    atom_count: int,
    result: TrackingResult | None,
) -> dict[str, object]:
    frames: list[tuple[int, int, float | None]] = []
    cages: set[tuple[int, str, str]] = set()
    guests: dict[str, tuple[str, tuple[int, ...]]] = {}
    components: dict[tuple[str, str], set[int]] = {}
    current_sources: dict[int, int] = {}
    with path.open("rb") as binary_handle:
        hashing = _HashingLineReader(binary_handle)
        reader = csv.DictReader(hashing, delimiter="\t")
        required = {
            "record",
            "render_frame",
            "source_frame",
            "time_ps",
            "family",
            "cage_id",
            "cage_type",
            "atom_indices",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid SQQ membership TSV: {path}")
        for row in reader:
            try:
                render_index = int(row["render_frame"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid render-frame index in {path}.") from exc
            record = str(row.get("record", ""))
            if record == "F":
                if render_index != len(frames) or render_index in current_sources:
                    raise ValueError(
                        "SQQ membership frame records must be unique and consecutive."
                    )
                try:
                    source_frame = int(row["source_frame"])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Invalid source-frame index in {path}.") from exc
                time_ps = _membership_time(row.get("time_ps"), path)
                frames.append((render_index, source_frame, time_ps))
                current_sources[render_index] = source_frame
                continue
            if render_index not in current_sources:
                raise ValueError(
                    f"Membership record precedes frame {render_index} in {path}."
                )
            indexes = _membership_atom_indexes(
                row.get("atom_indices"),
                atom_count=atom_count,
                path=path,
            )
            if record == "C" and row.get("family") == "cage":
                key = (
                    current_sources[render_index],
                    str(row.get("cage_id", "")),
                    str(row.get("cage_type", "")),
                )
                if key in cages:
                    raise ValueError(f"Duplicate cage membership record in {path}.")
                cages.add(key)
            elif record == "G" and row.get("family") == "guest":
                guest_id = str(row.get("cage_id", ""))
                identity = (str(row.get("cage_type", "")), indexes)
                previous = guests.setdefault(guest_id, identity)
                if previous != identity:
                    raise ValueError(
                        f"Guest {guest_id} changes atom identity in {path}."
                    )
            elif record == "P" and row.get("family") == "component":
                key = (str(row.get("cage_id", "")), str(row.get("cage_type", "")))
                component = components.setdefault(key, set())
                if component.intersection(indexes):
                    raise ValueError(f"Duplicate component atom index in {path}.")
                component.update(indexes)

    if result is not None:
        expected_frames = [
            (index, frame.frame_index, frame.time_ps)
            for index, frame in enumerate(result.frames)
        ]
        if len(frames) != len(expected_frames):
            raise ValueError(
                "Membership TSV and tracking state have different frame records."
            )
        for actual, expected in zip(frames, expected_frames):
            if actual[:2] != expected[:2] or not _same_optional_time(
                actual[2], expected[2]
            ):
                raise ValueError(
                    f"Membership frame/time mapping differs at render frame {actual[0]}."
                )
        expected_cages = {
            (item.frame_index, item.track_id, item.cage_type)
            for item in result.observations
        }
        if cages != expected_cages:
            missing = len(expected_cages.difference(cages))
            extra = len(cages.difference(expected_cages))
            raise ValueError(
                "Membership cage identities do not match tracking state "
                f"({missing} missing, {extra} unexpected)."
            )

    component_rows = [
        (role, name, sorted(indexes))
        for (role, name), indexes in sorted(components.items())
    ]
    return {
        "frame_count": len(frames),
        "component_signature_sha256": _json_digest(component_rows),
        "frame_mapping_sha256": _json_digest(frames),
        "cage_mapping_sha256": _json_digest(sorted(cages)),
        "sha256": hashing.digest.hexdigest(),
    }


def _membership_time(value: object, path: Path) -> float | None:
    if value in {None, "", "-"}:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid physical time in {path}.") from exc
    if not (float("-inf") < number < float("inf")):
        raise ValueError(f"Invalid physical time in {path}.")
    return number


def _membership_atom_indexes(
    value: object,
    *,
    atom_count: int,
    path: Path,
) -> tuple[int, ...]:
    if value in {None, "", "-"}:
        return ()
    try:
        indexes = tuple(int(item) for item in str(value).split(","))
    except ValueError as exc:
        raise ValueError(f"Invalid atom indexes in {path}.") from exc
    if len(indexes) != len(set(indexes)) or any(
        index < 0 or index >= atom_count for index in indexes
    ):
        raise ValueError(f"Out-of-range or duplicate atom indexes in {path}.")
    return indexes


def _file_identity(path: Path, *, digest: str | None = None) -> dict[str, object]:
    """Return the provenance record of one file, hashing it unless already done."""
    if digest is None:
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(block)
        digest = hasher.hexdigest()
    return {
        "name": path.name,
        "size": path.stat().st_size,
        "sha256": digest,
    }


def _json_digest(value: object) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _same_optional_time(left: float | None, right: float | None) -> bool:
    if left is None or right is None:
        return left is right
    return abs(left - right) <= max(1.0e-6, abs(right) * 1.0e-9)


def _validate_source_provenance(
    expected: Mapping[str, object],
    current: Mapping[str, object],
) -> None:
    if expected.get("format") != _SOURCE_PROVENANCE_FORMAT:
        raise ValueError("Unsupported tracking source-provenance format.")
    if expected.get("version") != _SOURCE_PROVENANCE_VERSION:
        raise ValueError("Unsupported tracking source-provenance version.")
    for key in (
        "atom_count",
        "frame_count",
        "topology_identity_sha256",
        "component_signature_sha256",
        "frame_mapping_sha256",
        "cage_mapping_sha256",
        "tracking_config_sha256",
    ):
        if expected.get(key) != current.get(key):
            raise ValueError(
                f"SQQ render source does not match track_state.json provenance: {key}."
            )
    expected_files = expected.get("files")
    current_files = current.get("files")
    if not isinstance(expected_files, Mapping) or not isinstance(
        current_files, Mapping
    ):
        raise ValueError("Invalid tracking source-provenance file records.")
    for role in ("topology", "trajectory", "membership", "script"):
        if expected_files.get(role) != current_files.get(role):
            raise ValueError(
                "SQQ render source does not match track_state.json provenance: "
                f"{role} file."
            )


_GRO_TITLE_FRAMES = re.compile(r"\bframes=(\d+)\b")


def _bundle_frame_count(gro: Path, membership: Path) -> int:
    """Read the frame count from the GRO title; scan the TSV only when absent.

    ``_write_topology_gro`` records ``frames=N`` in the title line, so a
    discovered package normally costs one line read instead of a full
    membership pass. Validation later cross-checks the TSV, XTC, and state.
    """
    try:
        with Path(gro).open("r", encoding="utf-8", errors="replace") as handle:
            title = handle.readline()
    except OSError:
        title = ""
    match = _GRO_TITLE_FRAMES.search(title)
    if match is not None:
        return int(match.group(1))
    return _membership_frame_count(membership)


def _membership_frame_count(path: Path) -> int:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "record" not in reader.fieldnames:
            raise ValueError(f"Invalid SQQ membership TSV: {path}")
        return sum(1 for row in reader if row.get("record") == "F")


def _atomic_write_text(path: Path, text: str, *, encoding: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding=encoding, newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and source.samefile(target):
        return
    temporary = _temporary_path(target)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_link_or_copy(source: Path, target: Path) -> None:
    source = Path(source).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and source.samefile(target):
        return
    temporary = _temporary_path(target)
    try:
        try:
            os.link(source, temporary)
        except OSError:
            temporary.unlink(missing_ok=True)
            shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def _temporary_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{uuid4().hex}.tmp")


__all__ = [
    "TRACK_GRO_NAME",
    "TRACK_MEMBERSHIP_NAME",
    "TRACK_RENDER_DIRECTORY",
    "TRACK_RENDER_NAMES",
    "TRACK_TCL_NAME",
    "TRACK_XTC_NAME",
    "build_tracking_source_provenance",
    "discover_sqq_cage_bundle",
    "discover_sqq_cage_gro",
    "publish_target_render_bundle",
    "publish_target_render_bundles",
    "validate_tracking_source_bundle",
]
