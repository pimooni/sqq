"""Atomic assembly and publication of complete SQQ render packages."""

from __future__ import annotations

import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import warnings
from tempfile import mkdtemp
from time import sleep
from typing import TYPE_CHECKING, Any, Iterable
from uuid import uuid4

import numpy as np

from .frame import (
    _ascii_annotation,
    _membership_token,
    _render_role,
    normalize_render_atom_scope,
    validate_render_fragment_lines,
)
from .models import (
    ANNOTATION_PREFIX,
    ATOM_PREFIX_WIDTH,
    COMPONENT_INDEX_CHUNK,
    FRAGMENT_DIRECTORY,
    FRAGMENT_DIRECTORY_GLOB,
    LEGACY_SQQ_CAGE_GRO_NAME,
    LEGACY_SQQ_CAGE_MEMBERSHIP_NAME,
    LEGACY_SQQ_CAGE_SCRIPT_NAME,
    LEGACY_SQQ_CAGE_XTC_NAME,
    LEGACY_SQQ_RENDER_SCRIPT_NAME,
    SQQ_CAGE_GRO_NAME,
    SQQ_CAGE_MEMBERSHIP_NAME,
    SQQ_CAGE_XTC_NAME,
    SQQ_RENDER_DIRECTORY,
    SQQ_RENDER_SCRIPT_NAME,
    RenderBundle,
    RenderFragment,
    RenderNames,
)
from .tcl import vmd_script_text

if TYPE_CHECKING:
    from ...models.tracking import TargetSelection, TrackingResult


def _prepare_fragment_workspace(outdir: Path) -> Path:
    """Create an isolated run-level fragment workspace."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    return Path(mkdtemp(prefix=f"{FRAGMENT_DIRECTORY}-", dir=root))


def _finalize_bundle(
    outdir: Path,
    fragments: Iterable[RenderFragment | Path] | None = None,
    *,
    fragment_dir: Path | None = None,
    tracking: TrackingResult | TargetSelection | None = None,
    cleanup: bool = True,
    names: RenderNames | None = None,
    render_kind: str = "analyze",
    molecule_name: str = "SQQ cages",
) -> RenderBundle:
    """Build and atomically publish one complete four-file render package."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    workspace = (
        Path(fragment_dir)
        if fragment_dir is not None
        else root / FRAGMENT_DIRECTORY
    )
    output_names = names or RenderNames()
    render_dir, gro_path, xtc_path, membership_path, script_path = (
        output_names.paths(root)
    )
    stage_dir = root / f".sqq-render-stage-{uuid4().hex}"
    stage_gro = stage_dir / output_names.topology
    stage_xtc = stage_dir / output_names.trajectory
    stage_membership = stage_dir / output_names.membership
    stage_script = stage_dir / output_names.script
    for legacy_name in (
        LEGACY_SQQ_CAGE_GRO_NAME,
        LEGACY_SQQ_CAGE_XTC_NAME,
        LEGACY_SQQ_CAGE_MEMBERSHIP_NAME,
        LEGACY_SQQ_CAGE_SCRIPT_NAME,
        LEGACY_SQQ_RENDER_SCRIPT_NAME,
    ):
        _best_effort_unlink(
            render_dir / legacy_name,
            "legacy render output",
        )
    manifests = _fragment_manifests(workspace, fragments)
    try:
        if not manifests:
            _remove_legacy_root_render_outputs(root)
            _resilient_rmtree(
                render_dir,
                warn=True,
                description="SQQ render directory",
            )
            return RenderBundle(None, None, 0, None, None, None)

        records = [_read_fragment_manifest(path) for path in manifests]
        records.sort(key=lambda item: item["frame_index"])
        _validate_fragment_records(records)
        stage_dir.mkdir(parents=True, exist_ok=False)
        # Digests are taken while the data are written so provenance never has
        # to re-read the published package.
        topology_digest, topology_identity, membership_digest = _write_render_data(
            stage_gro, stage_xtc, stage_membership, records
        )
        if tracking is not None:
            from ..tracking import rewrite_membership_track_ids

            membership_digest = rewrite_membership_track_ids(
                stage_membership, tracking, return_digest=True
            )[1]
        script_text = vmd_script_text(
            gro_filename=output_names.topology,
            xtc_filename=output_names.trajectory,
            membership_filename=output_names.membership,
            molecule_name=molecule_name,
            render_kind=render_kind,
        )
        _atomic_write_text(stage_script, script_text, encoding="ascii")
        file_digests = {
            "topology": _digest_record(stage_gro, output_names.topology, topology_digest),
            "trajectory": _digest_record(
                stage_xtc, output_names.trajectory, _sha256_file(stage_xtc)
            ),
            "membership": _digest_record(
                stage_membership, output_names.membership, membership_digest
            ),
            "script": _digest_record(
                stage_script,
                output_names.script,
                hashlib.sha256(script_text.encode("ascii")).hexdigest(),
            ),
        }

        _publish_render_directory(stage_dir, render_dir)
        _remove_legacy_root_render_outputs(root)
        return RenderBundle(
            gro_path=gro_path,
            script_path=script_path,
            frame_count=len(records),
            xtc_path=xtc_path,
            membership_path=membership_path,
            render_dir=render_dir,
            file_digests=file_digests,
            atom_count=int(records[0]["atom_count"]),
            topology_identity=topology_identity,
        )
    except Exception:
        for path, description in (
            (stage_gro, "partial topology GRO"),
            (stage_xtc, "partial XTC"),
            (stage_membership, "partial membership TSV"),
            (stage_script, "partial VMD script"),
        ):
            _best_effort_unlink(path, description)
        _resilient_rmtree(
            stage_dir,
            warn=False,
            description="partial render staging directory",
        )
        raise
    finally:
        if cleanup:
            _cleanup_fragment_workspace(workspace)

def _cleanup_fragment_workspace(fragment_dir: Path) -> bool:
    """Best-effort removal that cannot invalidate completed outputs."""
    return _resilient_rmtree(Path(fragment_dir), warn=True)


def _cleanup_fragment_workspaces(outdir: Path) -> None:
    """Remove abandoned legacy and run-isolated fragment workspaces."""
    root = Path(outdir)
    candidates = [root / FRAGMENT_DIRECTORY]
    if root.exists():
        candidates.extend(sorted(root.glob(FRAGMENT_DIRECTORY_GLOB)))
    for path in candidates:
        _cleanup_fragment_workspace(path)


def _remove_legacy_root_render_outputs(root: Path) -> None:
    for path in (
        root / SQQ_CAGE_GRO_NAME,
        root / SQQ_CAGE_XTC_NAME,
        root / SQQ_CAGE_MEMBERSHIP_NAME,
        root / SQQ_RENDER_SCRIPT_NAME,
        root / LEGACY_SQQ_CAGE_GRO_NAME,
        root / LEGACY_SQQ_CAGE_XTC_NAME,
        root / LEGACY_SQQ_CAGE_MEMBERSHIP_NAME,
        root / LEGACY_SQQ_CAGE_SCRIPT_NAME,
        root / LEGACY_SQQ_RENDER_SCRIPT_NAME,
    ):
        _best_effort_unlink(path, "legacy root-level render output")


def _remove_visible_render_outputs(root: Path) -> None:
    _remove_legacy_root_render_outputs(root)
    _resilient_rmtree(
        root / SQQ_RENDER_DIRECTORY,
        warn=True,
        description="SQQ render directory",
    )


def _cleanup_bundle(
    outdir: Path,
    *,
    fragment_dir: Path | None = None,
) -> None:
    """Remove visible bundle outputs and any abandoned fragments."""
    root = Path(outdir)
    _remove_visible_render_outputs(root)
    if fragment_dir is None:
        _cleanup_fragment_workspaces(root)
    else:
        _cleanup_fragment_workspace(fragment_dir)


def _publish_render_directory(stage: Path, target: Path) -> None:
    """Publish a complete package while preserving the previous one on failure."""
    target.parent.mkdir(parents=True, exist_ok=True)
    backup = target.parent / f".{target.name}.backup-{uuid4().hex}"
    moved_previous = False
    try:
        if target.exists():
            _replace_path(target, backup)
            moved_previous = True
        _replace_path(stage, target)
    except Exception:
        if moved_previous and backup.exists() and not target.exists():
            _replace_path(backup, target)
        raise
    finally:
        if backup.exists():
            _resilient_rmtree(
                backup,
                warn=True,
                description="previous SQQ render directory",
            )


def _replace_path(source: Path, target: Path, attempts: int = 5) -> None:
    """Retry transient directory-renaming failures on shared/Windows filesystems."""
    transient = {errno.EACCES, errno.EPERM, errno.EBUSY}
    for attempt in range(max(1, int(attempts))):
        try:
            os.replace(source, target)
            return
        except OSError as exc:
            if exc.errno not in transient or attempt + 1 >= attempts:
                raise
            sleep(0.05 * (2**attempt))


def _best_effort_unlink(path: Path, description: str) -> bool:
    try:
        path.unlink(missing_ok=True)
        return True
    except OSError as exc:
        warnings.warn(
            f"SQQ could not remove {description} {path}: {exc}",
            UserWarning,
            stacklevel=2,
        )
        return False


def _resilient_rmtree(
    path: Path,
    *,
    attempts: int = 5,
    initial_delay: float = 0.05,
    warn: bool,
    description: str = "temporary fragment directory",
) -> bool:
    transient = {
        errno.ENOTEMPTY,
        errno.EBUSY,
        errno.EACCES,
        errno.EPERM,
    }
    last_error: OSError | None = None
    for attempt in range(max(1, int(attempts))):
        try:
            shutil.rmtree(path)
            if not path.exists():
                return True
        except FileNotFoundError:
            return True
        except OSError as exc:
            last_error = exc
            if exc.errno not in transient:
                break
        if attempt + 1 < attempts:
            sleep(initial_delay * (2**attempt))
    if warn:
        detail = f": {last_error}" if last_error is not None else ""
        warnings.warn(
            f"SQQ could not remove {description} "
            f"{path}{detail}. Finalized outputs, if present, remain valid.",
            UserWarning,
            stacklevel=2,
        )
    return False

def _fragment_manifests(
    fragment_dir: Path,
    fragments: Iterable[RenderFragment | Path] | None,
) -> list[Path]:
    if fragments is None:
        return sorted(fragment_dir.glob("frame_*.json")) if fragment_dir.exists() else []
    paths: list[Path] = []
    for fragment in fragments:
        path = fragment.manifest_path if isinstance(fragment, RenderFragment) else Path(fragment)
        paths.append(path)
    return paths


def _read_fragment_manifest(path: Path) -> dict[str, Any]:
    try:
        record = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid SQQ cage fragment manifest: {path}") from exc
    required = {
        "status",
        "frame_index",
        "frame_name",
        "atom_count",
        "atom_signature",
        "effective_graph_mode",
        "gro_file",
    }
    missing = sorted(required.difference(record))
    if missing:
        raise ValueError(
            f"SQQ cage fragment manifest {path} is missing: {', '.join(missing)}"
        )
    if record["status"] != "ok":
        raise ValueError(f"SQQ cage fragment is not successful: {path}")
    record["manifest_path"] = path
    record["gro_path"] = path.parent / str(record["gro_file"])
    return record


def _validate_fragment_records(records: list[dict[str, Any]]) -> None:
    indexes = [int(record["frame_index"]) for record in records]
    if len(set(indexes)) != len(indexes):
        raise ValueError("Duplicate frame indexes in SQQ cage fragments.")
    reference = records[0]
    reference_scope = normalize_render_atom_scope(
        reference.get("atom_scope", "compact")
    )
    for record in records:
        scope = normalize_render_atom_scope(record.get("atom_scope", "compact"))
        if scope != reference_scope:
            raise ValueError(
                "sqq_cage.gro requires one render atom scope across frames; "
                f"{reference['frame_name']} uses {reference_scope} but "
                f"{record['frame_name']} uses {scope}."
            )
        if int(record["atom_count"]) != int(reference["atom_count"]):
            raise ValueError(
                "sqq_cage.gro requires a compatible atom topology across frames; "
                f"{reference['frame_name']} has {reference['atom_count']} atoms but "
                f"{record['frame_name']} has {record['atom_count']}."
            )
        if record["atom_signature"] != reference["atom_signature"]:
            raise ValueError(
                "sqq_cage.gro requires identical atom identity and order across "
                f"frames; {record['frame_name']} does not match "
                f"{reference['frame_name']}."
            )
        gro_path = Path(record["gro_path"])
        if not gro_path.is_file():
            raise ValueError(f"Missing SQQ cage fragment: {gro_path}")
        # Content validation happens in the single writing pass, where each
        # fragment is read exactly once.


def _write_render_data(
    gro_path: Path,
    xtc_path: Path,
    membership_path: Path,
    records: list[dict[str, Any]],
) -> tuple[str, str, str]:
    """Write topology, membership, and trajectory from one pass over the fragments.

    Every fragment GRO is parsed exactly once: the first fragment provides the
    topology, and each fragment's lines yield both its membership rows and its
    XTC coordinates. Returns the SHA-256 of the topology GRO, the atom-identity
    digest, and the SHA-256 of the membership TSV (all computed while writing).
    """
    frame_count = len(records)
    topology_digest = topology_identity = ""
    with _MembershipTsvWriter(membership_path) as membership, _XtcWriter(
        xtc_path, int(records[0]["atom_count"])
    ) as trajectory:
        for render_index, record in enumerate(records):
            lines = _fragment_lines(record)
            validate_render_fragment_lines(
                lines, Path(record["gro_path"]), int(record["atom_count"])
            )
            if render_index == 0:
                topology_digest, topology_identity = _write_topology_gro(
                    gro_path, record, frame_count, lines=lines
                )
            groups = _fragment_membership_groups(record, lines=lines)
            positions, box = _fragment_coordinates_and_box(record, lines=lines)
            membership.write_frame(render_index, record, groups)
            trajectory.write_frame(render_index, record, positions, box)
    return topology_digest, topology_identity, membership.digest


def _digest_record(path: Path, name: str, sha256_hex: str) -> dict[str, object]:
    return {"name": name, "size": path.stat().st_size, "sha256": sha256_hex}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class _HashingTextWriter:
    """Text-file proxy that hashes the encoded bytes of everything written."""

    def __init__(self, handle: Any, encoding: str) -> None:
        self._handle = handle
        self._encoding = encoding
        self.digest = hashlib.sha256()

    def write(self, text: str) -> int:
        self.digest.update(text.encode(self._encoding))
        return self._handle.write(text)

    def flush(self) -> None:
        self._handle.flush()

    def fileno(self) -> int:
        return self._handle.fileno()


def _fragment_lines(record: dict[str, Any]) -> list[str]:
    path = Path(record["gro_path"])
    lines = path.read_text(encoding="ascii").splitlines()
    atom_count = int(record["atom_count"])
    if len(lines) != atom_count + 3:
        raise ValueError(f"Invalid SQQ cage fragment: {path}")
    return lines


def _write_topology_gro(
    path: Path,
    record: dict[str, Any],
    frame_count: int,
    *,
    lines: list[str] | None = None,
) -> tuple[str, str]:
    """Write the topology GRO; return its SHA-256 and the atom-identity digest.

    The identity digest hashes the 20-character residue/atom prefix of every
    atom record exactly as ``io.render.tracking._gro_identity`` does when it
    re-reads a package, so both paths produce the same provenance value.
    """
    lines = _fragment_lines(record) if lines is None else lines
    atom_count = int(record["atom_count"])
    output = [
        f"SQQ cage topology frames={int(frame_count)}",
        f"{atom_count:5d}",
    ]
    atom_lines = [line[:ATOM_PREFIX_WIDTH] for line in lines[2 : 2 + atom_count]]
    output.extend(atom_lines)
    output.append(lines[2 + atom_count])
    text = "\n".join(output) + "\n"
    _atomic_write_text(path, text, encoding="ascii")
    identity = hashlib.sha256()
    for line in atom_lines:
        identity.update(line[:20].encode("utf-8"))
        identity.update(b"\n")
    return hashlib.sha256(text.encode("ascii")).hexdigest(), identity.hexdigest()


def _fragment_coordinates_and_box(
    record: dict[str, Any],
    *,
    lines: list[str] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    lines = _fragment_lines(record) if lines is None else lines
    atom_count = int(record["atom_count"])
    positions = np.empty((atom_count, 3), dtype=np.float32)
    for atom_index, line in enumerate(lines[2 : 2 + atom_count]):
        try:
            positions[atom_index] = (
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid GRO coordinates in {record['gro_path']} at atom "
                f"{atom_index}."
            ) from exc
    try:
        box = np.asarray(
            [float(value) for value in lines[2 + atom_count].split()[:3]],
            dtype=np.float32,
        )
    except ValueError as exc:
        raise ValueError(f"Invalid GRO box in {record['gro_path']}.") from exc
    if box.shape != (3,) or np.any(~np.isfinite(box)):
        raise ValueError(f"Invalid GRO box in {record['gro_path']}.")
    return positions, box


class _XtcWriter:
    """Frame-by-frame XTC writer with atomic publication."""

    def __init__(self, path: Path, atom_count: int) -> None:
        self.path = path
        self.atom_count = int(atom_count)
        self.temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}.xtc")
        self._writer = None
        self._timestep = None

    def __enter__(self) -> "_XtcWriter":
        try:
            import MDAnalysis as mda
            from MDAnalysis.coordinates.XTC import XTCWriter
        except ImportError as exc:
            raise RuntimeError(
                "Writing sqq_cage.xtc requires MDAnalysis."
            ) from exc
        universe = mda.Universe.empty(self.atom_count, trajectory=True)
        self._universe = universe
        self._timestep = universe.trajectory.ts
        self._writer = XTCWriter(
            str(self.temporary),
            self.atom_count,
            convert_units=True,
            precision=3,
        )
        return self

    def write_frame(
        self,
        render_index: int,
        record: dict[str, Any],
        positions_nm: np.ndarray,
        box_nm: np.ndarray,
    ) -> None:
        timestep = self._timestep
        timestep.positions = positions_nm * 10.0
        if np.any(box_nm > 0.0):
            timestep.dimensions = np.asarray(
                [
                    box_nm[0] * 10.0,
                    box_nm[1] * 10.0,
                    box_nm[2] * 10.0,
                    90.0,
                    90.0,
                    90.0,
                ],
                dtype=np.float32,
            )
        else:
            timestep.dimensions = None
        timestep.frame = render_index
        raw_time = record.get("time_ps")
        timestep.time = (
            float(render_index) if raw_time is None else float(raw_time)
        )
        timestep.data["step"] = int(record["frame_index"])
        self._writer.write(self._universe.atoms)

    def __exit__(self, exc_type, exc, traceback) -> None:
        try:
            if self._writer is not None:
                self._writer.close()
            if exc_type is None:
                os.replace(self.temporary, self.path)
        finally:
            self.temporary.unlink(missing_ok=True)


class _MembershipTsvWriter:
    """Frame-by-frame membership TSV writer that hashes what it writes."""

    _HEADER = (
        "record\trender_frame\tsource_frame\ttime_ps\tgraph_mode\tfamily\t"
        "cage_id\tcage_type\tphase\tdomain\tcluster\tatom_indices\t"
        "center_x_angstrom\tcenter_y_angstrom\tcenter_z_angstrom\n"
    )

    def __init__(self, path: Path) -> None:
        self.path = path
        self.temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
        self._raw = None
        self._output: _HashingTextWriter | None = None
        self._canonical_components: (
            tuple[tuple[str, str, tuple[int, ...]], ...] | None
        ) = None
        self.digest = ""

    def __enter__(self) -> "_MembershipTsvWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._raw = self.temporary.open("w", encoding="ascii", newline="\n")
        self._output = _HashingTextWriter(self._raw, "ascii")
        self._output.write(self._HEADER)
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        try:
            if self._raw is not None:
                self._raw.close()
            if exc_type is None:
                assert self._output is not None
                self.digest = self._output.digest.hexdigest()
                os.replace(self.temporary, self.path)
        finally:
            self.temporary.unlink(missing_ok=True)

    def write_frame(
        self,
        render_index: int,
        record: dict[str, Any],
        groups: dict[tuple[str, str], list[int]],
    ) -> None:
        output = self._output
        assert output is not None
        time_value = record.get("time_ps")
        time_text = "-" if time_value is None else f"{float(time_value):.9g}"
        graph_text = _tsv_field(record.get("graph_mode_display", "unknown"))
        output.write(
            "\t".join(
                (
                    "F",
                    str(render_index),
                    str(int(record["frame_index"])),
                    time_text,
                    graph_text,
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                )
            )
            + "\n"
        )
        for center in sorted(
            record.get("cage_centers", ()),
            key=lambda item: (
                str(item.get("cage_type", "")),
                str(item.get("cage_id", "")),
            ),
        ):
            required = {
                "cage_id",
                "cage_type",
                "phase",
                "domain",
                "cluster",
                "center_angstrom",
            }
            missing = sorted(required.difference(center))
            if missing:
                raise ValueError(
                    "Invalid SQQ cage-center metadata; missing "
                    + ", ".join(missing)
                    + "."
                )
            xyz = np.asarray(center["center_angstrom"], dtype=float)
            if xyz.shape != (3,) or np.any(~np.isfinite(xyz)):
                raise ValueError(
                    "Invalid SQQ cage center in fragment metadata."
                )
            output.write(
                "\t".join(
                    (
                        "C",
                        str(render_index),
                        "-",
                        "-",
                        "-",
                        "cage",
                        _membership_token(center["cage_id"]),
                        _membership_token(center["cage_type"]),
                        _membership_token(center["phase"]),
                        _membership_token(center["domain"]),
                        _membership_token(center["cluster"]),
                        "-",
                        *(f"{float(value):.17g}" for value in xyz),
                    )
                )
                + "\n"
            )
        for guest in sorted(
            record.get("guest_groups", ()),
            key=lambda item: str(item.get("guest_id", "")),
        ):
            required = {"guest_id", "resname", "atom_indices"}
            missing = sorted(required.difference(guest))
            if missing:
                raise ValueError(
                    "Invalid SQQ guest-group metadata; missing "
                    + ", ".join(missing)
                    + "."
                )
            atom_indexes = sorted(
                {int(value) for value in guest["atom_indices"]}
            )
            if (
                not atom_indexes
                or atom_indexes[0] < 0
                or atom_indexes[-1] >= int(record["atom_count"])
            ):
                raise ValueError(
                    "Invalid SQQ guest atom indexes in fragment metadata."
                )
            output.write(
                "\t".join(
                    (
                        "G",
                        str(render_index),
                        "-",
                        "-",
                        "-",
                        "guest",
                        _membership_token(guest["guest_id"]),
                        _membership_token(guest["resname"]),
                        "-",
                        "-",
                        "-",
                        ",".join(str(value) for value in atom_indexes),
                        "-",
                        "-",
                        "-",
                    )
                )
                + "\n"
            )
        components: list[tuple[str, str, tuple[int, ...]]] = []
        for component in sorted(
            record.get("component_groups", ()),
            key=lambda item: (
                str(item.get("role", "")),
                str(item.get("resname", "")),
            ),
        ):
            required = {"role", "resname", "atom_indices"}
            missing = sorted(required.difference(component))
            if missing:
                raise ValueError(
                    "Invalid SQQ component-group metadata; missing "
                    + ", ".join(missing)
                    + "."
                )
            role = _render_role(component["role"])
            atom_indexes = sorted(
                {int(value) for value in component["atom_indices"]}
            )
            if (
                not atom_indexes
                or atom_indexes[0] < 0
                or atom_indexes[-1] >= int(record["atom_count"])
            ):
                raise ValueError(
                    "Invalid SQQ component atom indexes in fragment metadata."
                )
            components.append(
                (
                    role,
                    _membership_token(component["resname"]),
                    tuple(atom_indexes),
                )
            )
        component_signature = tuple(components)
        if self._canonical_components is None:
            self._canonical_components = component_signature
        elif component_signature != self._canonical_components:
            raise ValueError(
                "SQQ render component topology changes between frames; "
                f"frame 0 and frame {render_index} do not match."
            )
        if render_index == 0:
            for role, resname, atom_indexes_tuple in components:
                atom_indexes = list(atom_indexes_tuple)
                for start in range(0, len(atom_indexes), COMPONENT_INDEX_CHUNK):
                    chunk = atom_indexes[start : start + COMPONENT_INDEX_CHUNK]
                    output.write(
                        "\t".join(
                            (
                                "P",
                                "0",
                                "-",
                                "-",
                                "-",
                                "component",
                                role,
                                resname,
                                "-",
                                "-",
                                "-",
                                ",".join(str(value) for value in chunk),
                                "-",
                                "-",
                                "-",
                            )
                        )
                        + "\n"
                    )
        for (family, membership), atom_indexes in sorted(groups.items()):
            cage_id, cage_type, phase, domain_id, cluster_id = (
                membership.split(":")
            )
            atoms = ",".join(str(value) for value in sorted(set(atom_indexes)))
            output.write(
                "\t".join(
                    (
                        "M",
                        str(render_index),
                        "-",
                        "-",
                        "-",
                        family,
                        cage_id,
                        cage_type,
                        phase,
                        domain_id,
                        cluster_id,
                        atoms,
                        "-",
                        "-",
                        "-",
                    )
                )
                + "\n"
            )


def _fragment_membership_groups(
    record: dict[str, Any],
    *,
    lines: list[str] | None = None,
) -> dict[tuple[str, str], list[int]]:
    lines = _fragment_lines(record) if lines is None else lines
    atom_count = int(record["atom_count"])
    groups: dict[tuple[str, str], list[int]] = {}
    for atom_index, line in enumerate(lines[2 : 2 + atom_count]):
        if ANNOTATION_PREFIX not in line:
            continue
        payload = line.split(ANNOTATION_PREFIX, 1)[1]
        if " g=" not in payload:
            raise ValueError(f"Invalid SQQ annotation in {record['gro_path']}.")
        cage_payload, guest_payload = payload.split(" g=", 1)
        for family, family_payload in (
            ("cage", cage_payload),
            ("guest", guest_payload),
        ):
            if family_payload in {"", "-"}:
                continue
            for membership in family_payload.split(","):
                if len(membership.split(":")) != 5:
                    raise ValueError(
                        f"Invalid SQQ membership in {record['gro_path']}."
                    )
                groups.setdefault((family, membership), []).append(atom_index)
    return groups


def _tsv_field(value: Any) -> str:
    text = _ascii_annotation(str(value))
    if any(character in text for character in "\t\r\n"):
        raise ValueError(f"Invalid SQQ TSV field: {value!r}")
    return text

def _atomic_write_text(path: Path, text: str, *, encoding: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_text(text, encoding=encoding, newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
