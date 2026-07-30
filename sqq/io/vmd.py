from __future__ import annotations

"""Run-level annotated GRO and VMD rendering outputs."""

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import json
import os
from pathlib import Path
import re
import shutil
import socket
import sys
from tempfile import mkdtemp
from time import sleep
from typing import Any, BinaryIO, Iterable, Iterator
from uuid import uuid4

import numpy as np

from ..display import graph_mode_display
from ..models import Atom, Frame, FrameResult
from .gro_grouping import gro_topology_fingerprint
from .gro_writer import ascii_gro_text
from .occupancy import guest_id


FRAGMENT_DIRECTORY = ".sqq-cage-fragments"
FRAGMENT_DIRECTORY_GLOB = f"{FRAGMENT_DIRECTORY}-*"
OUTPUT_LOCK_NAME = ".sqq.lock"
SQQ_RENDER_DIRECTORY = "sqq_render"
SQQ_CAGE_GRO_NAME = "sqq-cage.gro"
SQQ_CAGE_XTC_NAME = "sqq-cage.xtc"
SQQ_CAGE_MEMBERSHIP_NAME = "sqq-cage.membership.tsv"
SQQ_RENDER_SCRIPT_NAME = "sqq-cage.vmd.tcl"
LEGACY_SQQ_RENDER_SCRIPT_NAME = "sqq-render.vmd.tcl"
ANNOTATION_PREFIX = "; SQQ1 m="
ATOM_PREFIX_WIDTH = 44
EMPTY_VELOCITY_WIDTH = 24
ANNOTATION_COLUMN = ATOM_PREFIX_WIDTH + EMPTY_VELOCITY_WIDTH + 1
MEMBERSHIP_CLASSES = {
    "sI": "I",
    "sII": "II",
    "sH": "H",
    "boundary": "B",
    "unclassified": "U",
    "ambiguous": "A",
    "isolated": "X",
}


@dataclass(frozen=True)
class SqqCageFragment:
    """One complete annotated GRO frame and its compact manifest."""

    frame_index: int
    gro_path: Path
    manifest_path: Path
    atom_count: int
    atom_signature: str
    effective_graph_mode: str


@dataclass(frozen=True)
class SqqCageBundle:
    """Visible files produced by run-level bundle finalization."""

    gro_path: Path | None
    script_path: Path | None
    frame_count: int
    xtc_path: Path | None = None
    membership_path: Path | None = None
    render_dir: Path | None = None


@dataclass
class SqqOutputLock:
    """One process-held lock for an SQQ output root."""

    path: Path
    handle: BinaryIO
    token: str

    def release(self) -> None:
        if self.handle.closed:
            return
        try:
            self.handle.seek(0)
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(self.handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()


@dataclass(frozen=True)
class CageMembership:
    cage_id: str
    cage_type: str
    class_code: str
    domain_id: str
    cluster_id: str

    def encode(self) -> str:
        return ":".join(
            _membership_token(value)
            for value in (
                self.cage_id,
                self.cage_type,
                self.class_code,
                self.domain_id,
                self.cluster_id,
            )
        )


@contextmanager
def sqq_output_lock(outdir: Path) -> Iterator[SqqOutputLock]:
    """Prevent concurrent SQQ runs from sharing one output root."""
    lock = _acquire_output_lock(Path(outdir))
    try:
        yield lock
    finally:
        lock.release()


def prepare_sqq_cage_fragments(outdir: Path) -> Path:
    """Create an isolated run-level fragment workspace."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    return Path(mkdtemp(prefix=f"{FRAGMENT_DIRECTORY}-", dir=root))


def write_sqq_cage_fragment(
    result: FrameResult,
    fragment_dir: Path,
    frame_index: int,
    requested_graph_mode: str | None = None,
) -> SqqCageFragment:
    """Atomically write one complete annotated GRO block for later merging."""
    index = int(frame_index)
    if index < 0:
        raise ValueError("SQQ cage fragment frame_index must be non-negative.")
    root = Path(fragment_dir)
    root.mkdir(parents=True, exist_ok=True)
    stem = f"frame_{index:09d}"
    gro_path = root / f"{stem}.gro"
    manifest_path = root / f"{stem}.json"

    memberships = water_cage_memberships(result)
    output_atoms = visualization_atoms(result)
    guest_memberships = guest_cage_memberships(result)
    graph_display = _frame_graph_display(result, requested_graph_mode)
    block = annotated_gro_block(
        result,
        memberships,
        graph_display,
        atoms=output_atoms,
        guest_memberships=guest_memberships,
    )
    signature = atom_signature(output_atoms)
    manifest: dict[str, Any] = {
        "format": "SQQ cage fragment",
        "version": 1,
        "status": "ok",
        "frame_index": index,
        "frame_name": ascii_gro_text(result.frame.name),
        "time_ps": result.frame.time_ps,
        "atom_count": len(output_atoms),
        "atom_signature": signature,
        "effective_graph_mode": str(result.graph.mode),
        "requested_graph_mode": (
            None if requested_graph_mode is None else str(requested_graph_mode)
        ),
        "graph_mode_display": graph_display,
        "gro_file": gro_path.name,
    }
    _atomic_write_text(gro_path, block, encoding="ascii")
    _atomic_write_text(
        manifest_path,
        json.dumps(manifest, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
        + "\n",
        encoding="ascii",
    )
    return SqqCageFragment(
        frame_index=index,
        gro_path=gro_path,
        manifest_path=manifest_path,
        atom_count=len(output_atoms),
        atom_signature=signature,
        effective_graph_mode=str(result.graph.mode),
    )


def finalize_sqq_cage_bundle(
    outdir: Path,
    fragments: Iterable[SqqCageFragment | Path] | None = None,
    *,
    fragment_dir: Path | None = None,
    write_gro: bool = True,
    write_script: bool = True,
    cleanup: bool = True,
) -> SqqCageBundle:
    """Build the compact run-level VMD visualization bundle."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    workspace = (
        Path(fragment_dir)
        if fragment_dir is not None
        else root / FRAGMENT_DIRECTORY
    )
    render_dir = root / SQQ_RENDER_DIRECTORY
    gro_path = render_dir / SQQ_CAGE_GRO_NAME
    xtc_path = render_dir / SQQ_CAGE_XTC_NAME
    membership_path = render_dir / SQQ_CAGE_MEMBERSHIP_NAME
    script_path = render_dir / SQQ_RENDER_SCRIPT_NAME
    manifests = _fragment_manifests(workspace, fragments)
    try:
        if not manifests:
            _remove_visible_render_outputs(root)
            return SqqCageBundle(None, None, 0, None, None, None)

        records = [_read_fragment_manifest(path) for path in manifests]
        records.sort(key=lambda item: item["frame_index"])
        _validate_fragment_records(records)
        if write_script and not write_gro:
            raise ValueError(
                "sqq-render requires sqq-cage-gro; enable sqq-cage-gro first."
            )

        if write_script:
            render_dir.mkdir(parents=True, exist_ok=True)
            _write_compact_render_data(
                gro_path,
                xtc_path,
                membership_path,
                records,
            )
            _atomic_write_text(
                script_path,
                vmd_script_text(
                    xtc_filename=SQQ_CAGE_XTC_NAME,
                    membership_filename=SQQ_CAGE_MEMBERSHIP_NAME,
                ),
                encoding="ascii",
            )
        elif write_gro:
            render_dir.mkdir(parents=True, exist_ok=True)
            _write_topology_gro(gro_path, records[0], len(records))
            for path in (xtc_path, membership_path, script_path):
                path.unlink(missing_ok=True)
        else:
            _resilient_rmtree(
                render_dir,
                warn=True,
                description="SQQ render directory",
            )

        _remove_legacy_root_render_outputs(root)
        return SqqCageBundle(
            gro_path=gro_path if write_gro else None,
            script_path=script_path if write_script else None,
            frame_count=len(records),
            xtc_path=xtc_path if write_script else None,
            membership_path=membership_path if write_script else None,
            render_dir=render_dir if write_gro or write_script else None,
        )
    except Exception:
        for path, description in (
            (gro_path, "partial topology GRO"),
            (xtc_path, "partial XTC"),
            (membership_path, "partial membership TSV"),
            (script_path, "partial VMD script"),
        ):
            _best_effort_unlink(path, description)
        raise
    finally:
        if cleanup:
            cleanup_sqq_cage_fragment_workspace(workspace)

def cleanup_sqq_cage_fragment_workspace(fragment_dir: Path) -> bool:
    """Best-effort removal that cannot invalidate completed outputs."""
    return _resilient_rmtree(Path(fragment_dir), warn=True)


def cleanup_sqq_cage_fragments(outdir: Path) -> None:
    """Remove abandoned legacy and run-isolated fragment workspaces."""
    root = Path(outdir)
    candidates = [root / FRAGMENT_DIRECTORY]
    if root.exists():
        candidates.extend(sorted(root.glob(FRAGMENT_DIRECTORY_GLOB)))
    for path in candidates:
        cleanup_sqq_cage_fragment_workspace(path)


def _remove_legacy_root_render_outputs(root: Path) -> None:
    for path in (
        root / SQQ_CAGE_GRO_NAME,
        root / SQQ_CAGE_XTC_NAME,
        root / SQQ_CAGE_MEMBERSHIP_NAME,
        root / SQQ_RENDER_SCRIPT_NAME,
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


def cleanup_sqq_cage_bundle(
    outdir: Path,
    *,
    fragment_dir: Path | None = None,
) -> None:
    """Remove visible bundle outputs and any abandoned fragments."""
    root = Path(outdir)
    _remove_visible_render_outputs(root)
    if fragment_dir is None:
        cleanup_sqq_cage_fragments(root)
    else:
        cleanup_sqq_cage_fragment_workspace(fragment_dir)

def _acquire_output_lock(root: Path) -> SqqOutputLock:
    root.mkdir(parents=True, exist_ok=True)
    path = root / OUTPUT_LOCK_NAME
    handle = path.open("a+b")
    try:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        handle.close()
        owner = _output_lock_owner(path)
        detail = f" ({owner})" if owner else ""
        raise RuntimeError(
            f"SQQ output directory is already in use: {root}{detail}. "
            "Wait for the active run or choose another --output directory."
        ) from exc

    token = uuid4().hex
    metadata = {
        "format": "SQQ output lock",
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started_at": datetime.now(timezone.utc).isoformat(),
        "token": token,
    }
    handle.seek(0)
    handle.truncate()
    handle.write(
        (json.dumps(metadata, ensure_ascii=True, sort_keys=True) + "\n").encode(
            "ascii"
        )
    )
    handle.flush()
    try:
        os.fsync(handle.fileno())
    except OSError:
        pass
    handle.seek(0)
    return SqqOutputLock(path=path, handle=handle, token=token)


def _output_lock_owner(path: Path) -> str:
    try:
        metadata = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return ""
    values = []
    if metadata.get("pid") is not None:
        values.append(f"PID {metadata['pid']}")
    if metadata.get("host"):
        values.append(f"host {metadata['host']}")
    return ", ".join(values)


def _best_effort_unlink(path: Path, description: str) -> bool:
    try:
        path.unlink(missing_ok=True)
        return True
    except OSError as exc:
        print(
            f"Warning: SQQ could not remove {description} {path}: {exc}",
            file=sys.stderr,
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
        print(
            f"Warning: SQQ could not remove {description} "
            f"{path}{detail}. Finalized outputs, if present, remain valid.",
            file=sys.stderr,
        )
    return False


def _cage_membership_records(
    result: FrameResult,
) -> tuple[list[Any], dict[str, CageMembership]]:
    """Build one validated membership record for every output cage."""
    cages = list(result.all_cages or result.cages)
    cage_by_id: dict[str, Any] = {}
    for cage in cages:
        if cage.object_id in cage_by_id:
            raise ValueError(f"Duplicate cage id in SQQ cage output: {cage.object_id}")
        cage_by_id[cage.object_id] = cage
    classification = _cage_classification(result, cage_by_id)
    records: dict[str, CageMembership] = {}
    for cage in cages:
        class_code, domain_id, cluster_id = classification.get(
            cage.object_id, ("-", "-", "-")
        )
        records[cage.object_id] = CageMembership(
            cage_id=_compact_object_id(cage.object_id),
            cage_type=cage.cage_type,
            class_code=class_code,
            domain_id=_compact_object_id(domain_id),
            cluster_id=_compact_object_id(cluster_id),
        )
    return cages, records


def water_cage_memberships(result: FrameResult) -> dict[int, tuple[CageMembership, ...]]:
    """Map water-oxygen atom indexes to every cage membership."""
    cages, records = _cage_membership_records(result)
    memberships: dict[int, list[CageMembership]] = {}
    for cage in cages:
        item = records[cage.object_id]
        for oxygen in cage.waters:
            oxygen_index = int(oxygen)
            if oxygen_index < 0 or oxygen_index >= len(result.frame.atoms):
                raise ValueError(
                    f"Cage {cage.object_id} references invalid oxygen index "
                    f"{oxygen_index}."
                )
            memberships.setdefault(oxygen_index, []).append(item)
    return {index: tuple(items) for index, items in memberships.items()}


def guest_cage_memberships(result: FrameResult) -> dict[int, tuple[CageMembership, ...]]:
    """Map every atom of an assigned guest to all containing cages."""
    cages, records = _cage_membership_records(result)
    guests: dict[str, Any] = {}
    for guest in result.guests:
        identifier = guest_id(guest)
        if identifier in guests:
            raise ValueError(f"Duplicate guest id in SQQ cage output: {identifier}")
        guests[identifier] = guest

    by_guest: dict[str, list[CageMembership]] = {}
    seen: dict[str, set[str]] = {}
    for cage in cages:
        item = records[cage.object_id]
        for identifier in cage.guest_ids:
            if identifier not in guests:
                raise ValueError(
                    f"Cage {cage.object_id} references unknown guest id: {identifier}"
                )
            if cage.object_id in seen.setdefault(identifier, set()):
                continue
            seen[identifier].add(cage.object_id)
            by_guest.setdefault(identifier, []).append(item)

    memberships: dict[int, tuple[CageMembership, ...]] = {}
    for identifier, items in by_guest.items():
        guest = guests[identifier]
        for atom_index in guest.atoms:
            index = int(atom_index)
            if index < 0 or index >= len(result.frame.atoms):
                raise ValueError(
                    f"Guest {identifier} references invalid atom index {index}."
                )
            memberships[index] = tuple(items)
    return memberships


def visualization_atoms(result: FrameResult) -> list[Atom]:
    """Return stable VMD atoms: every water oxygen and every complete guest."""
    included = {int(water.oxygen) for water in result.waters}
    included.update(
        int(atom_index)
        for guest in result.guests
        for atom_index in guest.atoms
    )
    atoms = [atom for atom in result.frame.atoms if int(atom.index) in included]
    if len(atoms) != len(included):
        missing = sorted(included.difference(int(atom.index) for atom in atoms))
        raise ValueError(f"SQQ visualization references missing atom indexes: {missing}.")
    return atoms


def annotated_gro_block(
    result: FrameResult,
    memberships: dict[int, tuple[CageMembership, ...]],
    graph_display: str,
    *,
    atoms: Iterable[Atom] | None = None,
    guest_memberships: dict[int, tuple[CageMembership, ...]] | None = None,
) -> str:
    """Return one complete ASCII GRO block with SQQ annotations."""
    title_parts = ["SQQ cage", f"frame={ascii_gro_text(result.frame.name)}"]
    if result.frame.time_ps is not None:
        time_value = float(result.frame.time_ps)
        if not np.isfinite(time_value):
            raise ValueError("SQQ cage GRO time_ps must be finite when provided.")
        title_parts.append(f"time_ps={time_value:.9g}")
    title_parts.append("graph=" + ascii_gro_text(graph_display))
    output_atoms = list(result.frame.atoms if atoms is None else atoms)
    lines = [" ".join(title_parts), f"{len(output_atoms):5d}"]
    guest_memberships = guest_memberships or {}
    for atom in output_atoms:
        encoded = ",".join(item.encode() for item in memberships.get(int(atom.index), ()))
        guest_encoded = ",".join(
            item.encode() for item in guest_memberships.get(int(atom.index), ())
        )
        lines.append(_annotated_atom_line(atom, encoded or "-", guest_encoded or "-"))
    lines.append(_box_line(result.frame.box))
    return "\n".join(lines) + "\n"


def atom_signature(atoms: Iterable[Atom]) -> str:
    """Hash topology-compatible atom and molecule ordering for VMD frames."""
    frame = Frame(name="", atoms=list(atoms))
    return gro_topology_fingerprint(frame)


def _tcl_braced_literal(value: str, *, label: str) -> str:
    """Return a restricted, substitution-free Tcl literal."""
    text = str(value)
    if not text:
        raise ValueError(f"{label} must not be empty.")
    if any(character in text for character in "{}\\\r\n\0"):
        raise ValueError(f"{label} contains unsupported Tcl characters.")
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"{label} must contain ASCII text only.") from exc
    return "{" + text + "}"


def vmd_script_text(
    gro_filename: str = SQQ_CAGE_GRO_NAME,
    xtc_filename: str = SQQ_CAGE_XTC_NAME,
    membership_filename: str = SQQ_CAGE_MEMBERSHIP_NAME,
    molecule_name: str = "SQQ cages",
) -> str:
    """Return the shared, parameterized ASCII Tcl renderer."""
    replacements = {
        "__SQQ_GRO_FILENAME__": (gro_filename, "VMD GRO filename"),
        "__SQQ_XTC_FILENAME__": (xtc_filename, "VMD XTC filename"),
        "__SQQ_MEMBERSHIP_FILENAME__": (
            membership_filename,
            "VMD membership filename",
        ),
        "__SQQ_MOLECULE_NAME__": (molecule_name, "VMD molecule name"),
    }
    script = _VMD_SCRIPT
    for placeholder, (filename, label) in replacements.items():
        value = str(filename)
        if value in {".", ".."} or "/" in value or "\\" in value:
            raise ValueError(f"{label} must name a file in the script directory.")
        script = script.replace(
            placeholder,
            _tcl_braced_literal(value, label=label),
        )
    if "__SQQ_" in script:
        raise AssertionError("Unresolved SQQ Tcl template placeholder.")
    return script

def _cage_classification(
    result: FrameResult,
    cage_by_id: dict[str, Any],
) -> dict[str, tuple[str, str, str]]:
    if not result.hydrate_cluster_enabled:
        return {}

    cluster_for: dict[str, str] = {}
    category_for: dict[str, str] = {}
    for cluster in result.hydrate_clusters:
        cluster_id = str(cluster.object_id)
        _require_known_cages(cluster.cage_ids, cage_by_id, f"cluster {cluster_id}")
        for cage_id in cluster.cage_ids:
            _claim_unique(cluster_for, cage_id, cluster_id, "cluster")
        for label, cage_ids in (
            ("B", cluster.boundary_cage_ids),
            ("A", cluster.ambiguous_cage_ids),
            ("U", cluster.unclassified_cage_ids),
        ):
            _require_known_cages(cage_ids, cage_by_id, f"cluster {cluster_id}")
            for cage_id in cage_ids:
                _claim_unique(category_for, cage_id, label, "cluster category")
                assigned_cluster = cluster_for.get(cage_id)
                if assigned_cluster != cluster_id:
                    raise ValueError(
                        f"Cage {cage_id} category belongs to {cluster_id}, but its "
                        f"cluster membership is {assigned_cluster or 'missing'}."
                    )

    output: dict[str, tuple[str, str, str]] = {}
    for domain in result.hydrate_domains:
        domain_id = str(domain.object_id)
        cluster_id = str(domain.cluster_id)
        class_code = MEMBERSHIP_CLASSES.get(str(domain.hydrate_type))
        if class_code not in {"I", "II", "H"}:
            raise ValueError(
                f"Unsupported hydrate domain type for {domain_id}: "
                f"{domain.hydrate_type}"
            )
        _require_known_cages(domain.cage_ids, cage_by_id, f"domain {domain_id}")
        for cage_id in domain.cage_ids:
            assigned_cluster = cluster_for.get(cage_id)
            if assigned_cluster != cluster_id:
                raise ValueError(
                    f"Cage {cage_id} in {domain_id} references cluster {cluster_id}, "
                    f"but its cluster membership is {assigned_cluster or 'missing'}."
                )
            if cage_id in category_for:
                raise ValueError(
                    f"Cage {cage_id} is both phase-classified and "
                    f"{category_for[cage_id]}-classified."
                )
            value = (class_code, domain_id, cluster_id)
            _claim_unique(output, cage_id, value, "hydrate domain")

    for cage_id, class_code in category_for.items():
        if cage_id in output:
            raise ValueError(f"Cage {cage_id} has conflicting hydrate classifications.")
        output[cage_id] = (class_code, "-", cluster_for[cage_id])

    for cage_id, cluster_id in cluster_for.items():
        output.setdefault(cage_id, ("U", "-", cluster_id))

    _require_known_cages(result.isolated_cage_ids, cage_by_id, "isolated cage list")
    for cage_id in result.isolated_cage_ids:
        if cage_id in cluster_for or cage_id in output:
            raise ValueError(
                f"Isolated cage {cage_id} also belongs to a hydrate cluster."
            )
        output[cage_id] = ("X", "-", "-")
    return output


def _claim_unique(
    mapping: dict[str, Any],
    key: str,
    value: Any,
    label: str,
) -> None:
    previous = mapping.setdefault(key, value)
    if previous != value:
        raise ValueError(
            f"Cage {key} has conflicting {label} assignments: "
            f"{previous} and {value}."
        )


def _require_known_cages(
    cage_ids: Iterable[str],
    cage_by_id: dict[str, Any],
    owner: str,
) -> None:
    missing = [str(cage_id) for cage_id in cage_ids if cage_id not in cage_by_id]
    if missing:
        raise ValueError(
            f"{owner} references unknown cage ids: " + ", ".join(missing[:10])
        )
def _compact_object_id(value: Any) -> str:
    """Use the numeric suffix of a stable SQQ object id in GRO annotations."""
    text = _ascii_annotation(str(value) if value not in {None, ""} else "-")
    if text == "-":
        return text
    match = re.search(r"(\d+)$", text)
    return str(int(match.group(1))) if match else text




def _annotated_atom_line(
    atom: Atom,
    encoded_memberships: str,
    encoded_guest_memberships: str = "-",
) -> str:
    xyz = np.asarray(atom.xyz, dtype=float)
    if xyz.shape != (3,) or np.any(~np.isfinite(xyz)):
        raise ValueError(f"Invalid GRO coordinates for atom index {atom.index}.")
    coordinates = "".join(_gro_coordinate(value) for value in xyz)
    prefix = (
        f"{int(atom.resid) % 100000:5d}"
        f"{ascii_gro_text(atom.resname)[:5]:>5}"
        f"{ascii_gro_text(atom.atomname)[:5]:>5}"
        f"{int(atom.atomid) % 100000:5d}"
        f"{coordinates}"
    )
    if len(prefix) != ATOM_PREFIX_WIDTH:
        raise ValueError(f"Invalid fixed-width GRO atom record for index {atom.index}.")
    annotation = _ascii_annotation(encoded_memberships)
    if atom.velocity is None:
        velocity_text = " " * EMPTY_VELOCITY_WIDTH
    else:
        velocity = np.asarray(atom.velocity, dtype=float)
        if velocity.shape != (3,) or np.any(~np.isfinite(velocity)):
            raise ValueError(f"Invalid GRO velocities for atom index {atom.index}.")
        velocity_text = "".join(_gro_velocity(value) for value in velocity)
    guest_annotation = _ascii_annotation(encoded_guest_memberships)
    line = prefix + velocity_text + ANNOTATION_PREFIX + annotation + " g=" + guest_annotation
    if line.index(";") + 1 != ANNOTATION_COLUMN:
        raise AssertionError("SQQ GRO annotation column is not 69.")
    return line


def _gro_coordinate(value: float) -> str:
    field = f"{float(value):8.3f}"
    if len(field) != 8:
        raise ValueError(
            f"Coordinate {value!r} does not fit the GRO 8.3 fixed-width field."
        )
    return field


def _gro_velocity(value: float) -> str:
    field = f"{float(value):8.4f}"
    if len(field) != 8:
        raise ValueError(
            f"Velocity {value!r} does not fit the GRO 8.4 fixed-width field."
        )
    return field


def _box_line(box: np.ndarray | None) -> str:
    if box is None:
        return "   0.00000   0.00000   0.00000"
    values = np.asarray(box, dtype=float).reshape(-1)
    if len(values) < 3 or np.any(~np.isfinite(values[:3])):
        raise ValueError("SQQ cage GRO box requires three finite lengths.")
    return f"{values[0]:10.5f}{values[1]:10.5f}{values[2]:10.5f}"


def _membership_token(value: Any) -> str:
    text = _ascii_annotation(str(value) if value not in {None, ""} else "-")
    if any(character in text for character in ",:\t\r\n "):
        raise ValueError(f"Invalid SQQ membership token: {value!r}")
    return text


def _ascii_annotation(value: str) -> str:
    text = str(value)
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"SQQ GRO annotation is not ASCII: {value!r}") from exc
    return text


def _frame_graph_display(
    result: FrameResult,
    requested_graph_mode: str | None,
) -> str:
    if requested_graph_mode is None:
        return str(result.graph.mode)
    return graph_mode_display(requested_graph_mode, [result.graph.mode])


def _fragment_manifests(
    fragment_dir: Path,
    fragments: Iterable[SqqCageFragment | Path] | None,
) -> list[Path]:
    if fragments is None:
        return sorted(fragment_dir.glob("frame_*.json")) if fragment_dir.exists() else []
    paths: list[Path] = []
    for fragment in fragments:
        path = fragment.manifest_path if isinstance(fragment, SqqCageFragment) else Path(fragment)
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
    for record in records:
        if int(record["atom_count"]) != int(reference["atom_count"]):
            raise ValueError(
                "sqq-cage.gro requires a compatible atom topology across frames; "
                f"{reference['frame_name']} has {reference['atom_count']} atoms but "
                f"{record['frame_name']} has {record['atom_count']}."
            )
        if record["atom_signature"] != reference["atom_signature"]:
            raise ValueError(
                "sqq-cage.gro requires identical atom identity and order across "
                f"frames; {record['frame_name']} does not match "
                f"{reference['frame_name']}."
            )
        gro_path = Path(record["gro_path"])
        if not gro_path.is_file():
            raise ValueError(f"Missing SQQ cage fragment: {gro_path}")
        _validate_gro_fragment(gro_path, int(record["atom_count"]))


def _validate_gro_fragment(path: Path, atom_count: int) -> None:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Invalid SQQ cage GRO fragment: {path}") from exc
    if len(lines) != atom_count + 3:
        raise ValueError(
            f"SQQ cage fragment {path} is incomplete: expected "
            f"{atom_count + 3} records, got {len(lines)}."
        )
    try:
        declared = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid atom count in SQQ cage fragment: {path}") from exc
    if declared != atom_count:
        raise ValueError(
            f"SQQ cage fragment {path} declares {declared} atoms, expected "
            f"{atom_count}."
        )
    for line in lines[2 : 2 + atom_count]:
        if len(line) < ANNOTATION_COLUMN or line[ANNOTATION_COLUMN - 1] != ";":
            raise ValueError(
                f"SQQ annotation is not in column {ANNOTATION_COLUMN}: {path}"
            )


def _merge_gro_fragments(path: Path, records: list[dict[str, Any]]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as output:
            for record in records:
                data = Path(record["gro_path"]).read_bytes()
                if not data.endswith(b"\n"):
                    raise ValueError(
                        f"SQQ cage fragment does not end with a newline: "
                        f"{record['gro_path']}"
                    )
                output.write(data)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_compact_render_data(
    gro_path: Path,
    xtc_path: Path,
    membership_path: Path,
    records: list[dict[str, Any]],
) -> None:
    _write_topology_gro(gro_path, records[0], len(records))
    _write_membership_tsv(membership_path, records)
    _write_xtc(xtc_path, records)


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
) -> None:
    lines = _fragment_lines(record)
    atom_count = int(record["atom_count"])
    output = [
        f"SQQ cage topology frames={int(frame_count)}",
        f"{atom_count:5d}",
    ]
    output.extend(line[:ATOM_PREFIX_WIDTH] for line in lines[2 : 2 + atom_count])
    output.append(lines[2 + atom_count])
    _atomic_write_text(path, "\n".join(output) + "\n", encoding="ascii")


def _fragment_coordinates_and_box(
    record: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    lines = _fragment_lines(record)
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


def _write_xtc(path: Path, records: list[dict[str, Any]]) -> None:
    try:
        import MDAnalysis as mda
        from MDAnalysis.coordinates.XTC import XTCWriter
    except ImportError as exc:
        raise RuntimeError(
            "Writing sqq-cage.xtc requires MDAnalysis."
        ) from exc

    atom_count = int(records[0]["atom_count"])
    universe = mda.Universe.empty(atom_count, trajectory=True)
    timestep = universe.trajectory.ts
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}.xtc")
    try:
        with XTCWriter(
            str(temporary),
            atom_count,
            convert_units=True,
            precision=3,
        ) as writer:
            for render_index, record in enumerate(records):
                positions_nm, box_nm = _fragment_coordinates_and_box(record)
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
                writer.write(universe.atoms)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_membership_tsv(
    path: Path,
    records: list[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    header = (
        "record\trender_frame\tsource_frame\ttime_ps\tgraph_mode\tfamily\t"
        "cage_id\tcage_type\tphase\tdomain\tcluster\tatom_indices\n"
    )
    try:
        with temporary.open("w", encoding="ascii", newline="\n") as output:
            output.write(header)
            for render_index, record in enumerate(records):
                time_value = record.get("time_ps")
                time_text = "-" if time_value is None else f"{float(time_value):.9g}"
                graph_text = _tsv_field(record.get("graph_mode_display", "unknown"))
                output.write(
                    f"F\t{render_index}\t{int(record['frame_index'])}\t"
                    f"{time_text}\t{graph_text}\t-\t-\t-\t-\t-\t-\t-\n"
                )
                groups = _fragment_membership_groups(record)
                for (family, membership), atom_indexes in sorted(groups.items()):
                    cage_id, cage_type, phase, domain_id, cluster_id = (
                        membership.split(":")
                    )
                    atoms = ",".join(str(value) for value in sorted(set(atom_indexes)))
                    output.write(
                        f"M\t{render_index}\t-\t-\t-\t{family}\t{cage_id}\t"
                        f"{cage_type}\t{phase}\t{domain_id}\t{cluster_id}\t{atoms}\n"
                    )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _fragment_membership_groups(
    record: dict[str, Any],
) -> dict[tuple[str, str], list[int]]:
    lines = _fragment_lines(record)
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


_VMD_SCRIPT = r'''# SQQ annotated cage and guest renderer for VMD.
namespace eval ::SQQ {
    catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
    catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_changed}
    if {[info exists frame_after_id] && $frame_after_id ne ""} {
        catch {after cancel $frame_after_id}
    }
    variable molid -1
    variable gro_path ""
    variable xtc_path ""
    variable membership_path ""
    variable active_families {cage}
    variable custom_show_active 0
    variable active_targets
    variable representation_names {}
    variable frame_after_id ""
    variable displayed_graph_mode "__unset__"
    variable label_visible 0
    variable pick_mode off
    variable selected_target {}
    variable graphics_ids {}
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable color_overrides
    variable known_objects
    variable object_aliases
    variable cage_types
    foreach name {group_keys group_atoms graph_mode color_overrides known_objects object_aliases cage_types active_targets} {
        catch {array unset $name}
        array set $name {}
    }
    set active_targets(cage) [list [list cage *]]
}

proc ::SQQ::add_member {frame source key atom_index} {
    variable group_keys
    variable group_atoms
    variable known_objects
    set known_objects($source,$key) 1
    set group_key "$frame,$source"
    set atom_key "$frame,$source,$key"
    if {![info exists group_atoms($atom_key)]} {
        lappend group_keys($group_key) $key
        set group_atoms($atom_key) {}
    }
    lappend group_atoms($atom_key) $atom_index
}

proc ::SQQ::register_object {source key} {
    variable known_objects
    set known_objects($source,$key) 1
}

proc ::SQQ::register_cage_object {frame family cage_type object_id} {
    variable object_aliases
    variable cage_types
    set id_source "${family}-id"
    ::SQQ::register_object $family $cage_type
    ::SQQ::register_object $id_source $object_id
    set cage_types($frame,$object_id) $cage_type
    set compact [string map [list "^" "" "-" "" "_" ""] $cage_type]
    if {$compact ne $cage_type} {
        set object_aliases($family,$compact) $cage_type
        if {![string match "${cage_type}_*" $object_id]} { return }
        set suffix [string range $object_id [string length $cage_type] end]
        set object_aliases($id_source,${compact}${suffix}) $object_id
    }
}

proc ::SQQ::register_cage {frame cage_type object_id} {
    ::SQQ::register_cage_object $frame cage $cage_type $object_id
}

proc ::SQQ::register_guest_cage {frame cage_type object_id} {
    ::SQQ::register_cage_object $frame guest $cage_type $object_id
}

proc ::SQQ::deduplicate_frame_memberships {frame} {
    variable group_keys
    variable group_atoms
    foreach source {cage-id guest-id phase cluster domain} {
        set group_key "$frame,$source"
        if {![info exists group_keys($group_key)]} { continue }
        foreach key $group_keys($group_key) {
            set atom_key "$frame,$source,$key"
            set group_atoms($atom_key) [lsort -integer -unique $group_atoms($atom_key)]
        }
    }
}

proc ::SQQ::numbered_id {prefix value} {
    if {[regexp {^[0-9]+$} $value]} {
        scan $value %d number
        return [format "%s_%05d" $prefix $number]
    }
    return [string tolower $value]
}

proc ::SQQ::cage_object_id {cage_type cage_id} {
    if {[regexp {^[0-9]+$} $cage_id]} {
        scan $cage_id %d number
        return [format "%s_%05d" $cage_type $number]
    }
    return "${cage_type}_$cage_id"
}

proc ::SQQ::cage_type_from_object_id {frame object_id} {
    variable cage_types
    if {[info exists cage_types($frame,$object_id)]} {
        return $cage_types($frame,$object_id)
    }
    if {[regexp {^(.+)_([0-9]+)$} $object_id -> cage_type number]} {
        return $cage_type
    }
    return ""
}

proc ::SQQ::read_memberships {frame atom_index family payload} {
    if {$payload eq "" || $payload eq "-"} { return }
    foreach membership [split $payload ,] {
        set fields [split $membership :]
        if {[llength $fields] != 5} {
            error "Invalid SQQ $family membership in frame $frame: $membership"
        }
        lassign $fields cage_id cage_type phase domain_id cluster_id
        if {$cage_type eq "-"} { continue }
        set object_id [::SQQ::cage_object_id $cage_type $cage_id]
        if {$family eq "cage"} {
            ::SQQ::register_cage $frame $cage_type $object_id
            ::SQQ::add_member $frame cage-id $object_id $atom_index
            if {$phase ne "-"} { ::SQQ::add_member $frame phase $phase $atom_index }
            if {$cluster_id ne "-"} {
                ::SQQ::add_member $frame cluster [::SQQ::numbered_id cluster $cluster_id] $atom_index
            }
            if {$domain_id ne "-"} {
                ::SQQ::add_member $frame domain [::SQQ::numbered_id domain $domain_id] $atom_index
            }
        } elseif {$family eq "guest"} {
            ::SQQ::register_guest_cage $frame $cage_type $object_id
            ::SQQ::add_member $frame guest-id $object_id $atom_index
        } else {
            error "Invalid SQQ membership family '$family'"
        }
    }
}

proc ::SQQ::read_membership_tsv {path} {
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable known_objects
    variable object_aliases
    variable cage_types
    foreach name {group_keys group_atoms graph_mode known_objects object_aliases cage_types} {
        array unset $name
        array set $name {}
    }
    set handle [open $path r]
    fconfigure $handle -encoding ascii -translation auto
    if {[gets $handle header] < 0 || $header ne "record\trender_frame\tsource_frame\ttime_ps\tgraph_mode\tfamily\tcage_id\tcage_type\tphase\tdomain\tcluster\tatom_indices"} {
        close $handle
        error "Invalid SQQ membership TSV header: $path"
    }
    array set seen_frames {}
    set frame_count 0
    set line_number 1
    while {[gets $handle line] >= 0} {
        incr line_number
        if {$line eq ""} { continue }
        set fields [split $line "\t"]
        if {[llength $fields] != 12} {
            close $handle
            error "Invalid SQQ membership TSV row $line_number"
        }
        lassign $fields record frame source_frame time_ps graph_value family cage_id cage_type phase domain_id cluster_id atom_indexes
        if {![string is integer -strict $frame] || $frame < 0} {
            close $handle
            error "Invalid SQQ render frame at TSV row $line_number"
        }
        if {$record eq "F"} {
            if {[info exists seen_frames($frame)]} {
                close $handle
                error "Duplicate SQQ frame metadata for frame $frame"
            }
            set seen_frames($frame) 1
            set graph_mode($frame) $graph_value
            if {$frame + 1 > $frame_count} { set frame_count [expr {$frame + 1}] }
        } elseif {$record eq "M"} {
            if {$family ni {cage guest}} {
                close $handle
                error "Invalid SQQ membership family at TSV row $line_number"
            }
            set payload "$cage_id:$cage_type:$phase:$domain_id:$cluster_id"
            if {$atom_indexes eq "-" || $atom_indexes eq ""} {
                close $handle
                error "Missing SQQ atom indexes at TSV row $line_number"
            }
            foreach atom_index [split $atom_indexes ,] {
                if {![string is integer -strict $atom_index] || $atom_index < 0} {
                    close $handle
                    error "Invalid SQQ atom index at TSV row $line_number"
                }
                if {[catch {::SQQ::read_memberships $frame $atom_index $family $payload} message]} {
                    close $handle
                    error $message
                }
            }
        } else {
            close $handle
            error "Invalid SQQ membership record at TSV row $line_number: $record"
        }
    }
    close $handle
    for {set frame 0} {$frame < $frame_count} {incr frame} {
        if {![info exists seen_frames($frame)]} {
            error "Missing SQQ frame metadata for frame $frame"
        }
        ::SQQ::deduplicate_frame_memberships $frame
    }
    return $frame_count
}
proc ::SQQ::key_rank {family key} {
    if {$family in {cage guest}} {
        set order {512 51262 51263 51264 435663 51268}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    if {$family eq "phase"} {
        set order {I II H B A U X}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    return 0
}

proc ::SQQ::ordered_keys {family keys} {
    set decorated {}
    foreach key [lsort -unique $keys] {
        lappend decorated [list [::SQQ::key_rank $family $key] $key]
    }
    set output {}
    foreach item [lsort -integer -index 0 $decorated] { lappend output [lindex $item 1] }
    return $output
}

proc ::SQQ::standard_cage_rank {cage_type} {
    set order {512 51262 51263 51264 435663 51268}
    set position [lsearch -exact $order $cage_type]
    if {$position < 0} { return 0 }
    return [expr {$position + 1}]
}

proc ::SQQ::generic_cage_rank {cage_type} {
    set total_faces 0
    set ring_kinds 0
    set fields [regexp -all -inline {([0-9]+)\^([0-9]+)} $cage_type]
    for {set index 0} {$index < [llength $fields]} {incr index 3} {
        set count [lindex $fields [expr {$index + 2}]]
        if {[string is integer -strict $count]} {
            set total_faces [expr {$total_faces + $count}]
            incr ring_kinds
        }
    }
    return [list $total_faces $ring_kinds]
}

proc ::SQQ::object_render_key {frame object_id color_priority color_id explicit} {
    set cage_type [::SQQ::cage_type_from_object_id $frame $object_id]
    set exact [expr {$explicit || $color_priority == 3}]
    set standard_rank [::SQQ::standard_cage_rank $cage_type]
    if {$standard_rank > 0} {
        set standard 1
        set primary $standard_rank
        set secondary 0
    } else {
        set standard 0
        lassign [::SQQ::generic_cage_rank $cage_type] primary secondary
    }
    set exact_id [expr {$exact ? $object_id : ""}]
    return [list $exact $standard $primary $secondary $cage_type $exact_id $color_id]
}

proc ::SQQ::compare_object_render_keys {left right} {
    foreach index {0 1 2 3} {
        set left_value [lindex $left $index]
        set right_value [lindex $right $index]
        if {$left_value != $right_value} {
            return [expr {$left_value < $right_value ? -1 : 1}]
        }
    }
    foreach index {4 5} {
        set comparison [string compare [lindex $left $index] [lindex $right $index]]
        if {$comparison != 0} { return $comparison }
    }
    set left_color [lindex $left 6]
    set right_color [lindex $right 6]
    return [expr {$left_color < $right_color ? -1 : ($left_color > $right_color)}]
}

proc ::SQQ::cage_radius_tier {render_key} {
    if {[lindex $render_key 0]} { return 7 }
    if {[lindex $render_key 1]} { return [lindex $render_key 2] }
    return 0
}

proc ::SQQ::cage_layer_radius {tier tiers} {
    set count [llength $tiers]
    if {$count <= 1} { return 0.125 }
    set index [lsearch -exact $tiers $tier]
    return [format "%.3f" [expr {0.125 + 0.005 * $index / double($count - 1)}]]
}

proc ::SQQ::stable_color {key} {
    set palette {0 1 7 3 11 10 4 9 5 6}
    set hash 0
    foreach character [split $key ""] {
        scan $character %c code
        set hash [expr {(($hash * 33) + $code) & 0x7fffffff}]
    }
    return [lindex $palette [expr {$hash % [llength $palette]}]]
}

proc ::SQQ::color_id {family key} {
    if {$family in {cage guest}} {
        switch -- $key {
            512 { return 7 }
            51262 { return 0 }
            51263 { return 1 }
            51264 { return 3 }
            51268 { return 11 }
            435663 { return 10 }
            default { return 2 }
        }
    }
    if {$family eq "phase"} {
        switch -- $key {
            I { return 1 }
            II { return 0 }
            H { return 7 }
            B { return 3 }
            A { return 11 }
            U { return 2 }
            X { return 4 }
            default { return 2 }
        }
    }
    return [::SQQ::stable_color $key]
}

proc ::SQQ::source_family {source} {
    switch -- $source {
        cage - cage-id { return cage }
        guest - guest-id { return guest }
        phase - cluster - domain { return $source }
    }
    error "Unknown SQQ source '$source'"
}

proc ::SQQ::default_color_id {frame source key} {
    set family [::SQQ::source_family $source]
    if {$source in {cage-id guest-id}} {
        set cage_type [::SQQ::cage_type_from_object_id $frame $key]
        if {$cage_type ne ""} { return [::SQQ::color_id $family $cage_type] }
        return 2
    }
    return [::SQQ::color_id $family $key]
}

proc ::SQQ::effective_color {frame source key} {
    variable color_overrides
    set exact "$source,$key"
    if {[info exists color_overrides($exact)]} {
        if {$color_overrides($exact) eq "default"} {
            return [list [::SQQ::default_color_id $frame $source $key] 3]
        }
        return [list $color_overrides($exact) 3]
    }
    set family [::SQQ::source_family $source]
    if {$source in {cage-id guest-id}} {
        set cage_type [::SQQ::cage_type_from_object_id $frame $key]
        if {$cage_type ne ""} {
            set type_key "$family,$cage_type"
            if {[info exists color_overrides($type_key)]} {
                if {$color_overrides($type_key) eq "default"} {
                    return [list [::SQQ::default_color_id $frame $source $key] 2]
                }
                return [list $color_overrides($type_key) 2]
            }
        }
    }
    set category_key "$family,*"
    if {[info exists color_overrides($category_key)]} {
        return [list $color_overrides($category_key) 1]
    }
    return [list [::SQQ::default_color_id $frame $source $key] 0]
}

proc ::SQQ::phase_key {value} {
    set aliases [dict create si I i I sii II ii II sh H h H boundary B b B ambiguous A a A unclassified U u U isolated X x X]
    set key [string tolower $value]
    if {[dict exists $aliases $key]} { return [dict get $aliases $key] }
    return ""
}

proc ::SQQ::normalize_family {value} {
    set family [string tolower [string trim $value]]
    if {$family ni {cage guest phase cluster domain}} {
        error "SQQ family must be cage, guest, phase, cluster, or domain"
    }
    return $family
}

proc ::SQQ::is_family_token {value} {
    return [expr {[string tolower [string trim $value]] in {cage guest phase cluster domain}}]
}

proc ::SQQ::parse_target {family value} {
    variable known_objects
    variable object_aliases
    set family [::SQQ::normalize_family $family]
    set token [string trim $value]
    if {[string equal -nocase $token all]} { return [list $family *] }
    if {$family in {cage guest}} {
        set id_source "${family}-id"
        foreach source [list $id_source $family] {
            if {[info exists known_objects($source,$token)]} { return [list $source $token] }
            if {[info exists object_aliases($source,$token)]} {
                return [list $source $object_aliases($source,$token)]
            }
        }
        if {[regexp {^(.+)_([0-9]+)$} $token -> cage_type number]} {
            return [list $id_source [::SQQ::cage_object_id $cage_type $number]]
        }
        if {[regexp {^[0-9]+$} $token]} { return [list $family $token] }
    } elseif {$family eq "phase"} {
        set key [::SQQ::phase_key $token]
        if {$key ne ""} { return [list phase $key] }
    } elseif {$family eq "cluster"} {
        if {[regexp -nocase {^cluster_([0-9]+)$} $token -> number]} {
            return [list cluster [::SQQ::numbered_id cluster $number]]
        }
        if {[regexp {^[0-9]+$} $token]} {
            return [list cluster [::SQQ::numbered_id cluster $token]]
        }
    } elseif {$family eq "domain"} {
        if {[regexp -nocase {^domain_([0-9]+)$} $token -> number]} {
            return [list domain [::SQQ::numbered_id domain $number]]
        }
        if {[regexp {^[0-9]+$} $token]} {
            return [list domain [::SQQ::numbered_id domain $token]]
        }
    }
    error "Unknown SQQ $family target '$value'"
}

proc ::SQQ::require_known_target {family target} {
    variable known_objects
    lassign $target source key
    if {$key eq "*"} { return }
    if {![info exists known_objects($source,$key)]} {
        error "SQQ $family target '[::SQQ::target_label $family $target]' does not exist in the loaded trajectory"
    }
}

proc ::SQQ::phase_label {key} {
    switch -- $key {
        I { return sI }
        II { return sII }
        H { return sH }
        B { return boundary }
        A { return ambiguous }
        U { return unclassified }
        X { return isolated }
    }
    return $key
}

proc ::SQQ::target_label {family target} {
    lassign $target source key
    if {$key eq "*"} { return all }
    if {$family eq "phase"} { return [::SQQ::phase_label $key] }
    return $key
}

proc ::SQQ::parse_targets {family values} {
    if {[llength $values] == 0} { error "At least one SQQ target is required" }
    set targets {}
    foreach value $values {
        set target [::SQQ::parse_target $family $value]
        ::SQQ::require_known_target $family $target
        if {$target ni $targets} { lappend targets $target }
    }
    if {[llength $targets] > 1} {
        foreach target $targets {
            if {[lindex $target 1] eq "*"} { error "The all target must be used alone" }
        }
    }
    return $targets
}

proc ::SQQ::parse_show_groups {values} {
    if {[llength $values] < 2} {
        error "Usage: sqq show <family> <target...> ?<family> <target...> ...?"
    }
    set groups {}
    set family ""
    set raw_targets {}
    foreach value $values {
        if {[::SQQ::is_family_token $value]} {
            if {$family ne ""} {
                if {[llength $raw_targets] == 0} {
                    error "SQQ show family '$family' requires at least one target"
                }
                lappend groups [list $family [::SQQ::parse_targets $family $raw_targets]]
            }
            set family [::SQQ::normalize_family $value]
            set raw_targets {}
        } else {
            if {$family eq ""} {
                error "SQQ show must begin with cage, guest, phase, cluster, or domain"
            }
            lappend raw_targets $value
        }
    }
    if {$family eq "" || [llength $raw_targets] == 0} {
        if {$family eq ""} {
            error "SQQ show must begin with cage, guest, phase, cluster, or domain"
        }
        error "SQQ show family '$family' requires at least one target"
    }
    lappend groups [list $family [::SQQ::parse_targets $family $raw_targets]]
    return $groups
}

proc ::SQQ::merge_show_group {family targets} {
    variable active_families
    variable active_targets
    if {$family ni $active_families} {
        lappend active_families $family
        set active_targets($family) {}
    }
    set includes_all 0
    foreach target $targets {
        if {[lindex $target 1] eq "*"} {
            set includes_all 1
            break
        }
    }
    set already_all 0
    foreach target $active_targets($family) {
        if {[lindex $target 1] eq "*"} {
            set already_all 1
            break
        }
    }
    if {$includes_all} {
        set active_targets($family) $targets
    } elseif {!$already_all} {
        foreach target $targets {
            if {$target ni $active_targets($family)} {
                lappend active_targets($family) $target
            }
        }
    }
}

proc ::SQQ::set_show {values args} {
    variable active_families
    variable active_targets
    variable custom_show_active
    if {[llength $args] > 0} {
        set combined [list $values]
        foreach argument $args {
            foreach value $argument { lappend combined $value }
        }
        set values $combined
    }
    set groups [::SQQ::parse_show_groups $values]
    if {!$custom_show_active} {
        set active_families {}
        array unset active_targets
        array set active_targets {}
        set custom_show_active 1
    }
    foreach group $groups {
        lassign $group family targets
        ::SQQ::merge_show_group $family $targets
    }
    ::SQQ::render_current
}

proc ::SQQ::reset_show {{announce 0}} {
    variable active_families
    variable active_targets
    variable color_overrides
    variable custom_show_active
    variable label_visible
    variable pick_mode
    variable selected_target
    set active_families {cage}
    array unset active_targets
    array set active_targets {}
    set active_targets(cage) [list [list cage *]]
    array unset color_overrides
    array set color_overrides {}
    set custom_show_active 0
    set label_visible 0
    set pick_mode off
    set selected_target {}
    ::SQQ::render_current
    if {$announce} { puts "SQQ clear: restored source-time cage view" }
}

proc ::SQQ::set_label {values} {
    variable label_visible
    set count [llength $values]
    if {$count == 0} {
        set label_visible [expr {!$label_visible}]
    } elseif {$count == 1} {
        set value [string tolower [string trim [lindex $values 0]]]
        if {$value ni {on off}} { error "Usage: sqq show label ?on|off?" }
        set label_visible [expr {$value eq "on"}]
    } else {
        error "Usage: sqq show label ?on|off?"
    }
    ::SQQ::render_current
    puts "SQQ label: [expr {$label_visible ? "on" : "off"}]"
}

proc ::SQQ::set_pick_mode {value} {
    variable pick_mode
    variable selected_target
    set mode [string tolower [string trim $value]]
    if {$mode ni {center guest off}} {
        error "Usage: sqq pick center|guest|off"
    }
    set pick_mode $mode
    set selected_target {}
    ::SQQ::render_current
    if {$mode eq "off"} {
        puts "SQQ pick: off; restored opaque objects"
    } else {
        puts "SQQ pick: $mode; objects are transparent until one is selected"
    }
}

proc ::SQQ::base_material {} {
    variable pick_mode
    return [expr {$pick_mode eq "off" ? "Opaque" : "Transparent"}]
}

proc ::SQQ::color_value {value} {
    if {[string equal -nocase $value default]} { return default }
    set names [colorinfo colors]
    if {[string is integer -strict $value]} {
        scan $value %d color_id
        if {$color_id < 0 || $color_id >= [llength $names]} {
            error "VMD ColorID must be between 0 and [expr {[llength $names] - 1}]"
        }
        return $color_id
    }
    set color_id [lsearch -nocase -exact $names $value]
    if {$color_id < 0} { error "Unknown VMD color '$value'" }
    return $color_id
}

proc ::SQQ::clear_family_colors {family} {
    variable color_overrides
    foreach name [array names color_overrides] {
        set source [lindex [split $name ,] 0]
        if {[::SQQ::source_family $source] eq $family} { unset color_overrides($name) }
    }
}

proc ::SQQ::set_colors {family values color} {
    variable color_overrides
    set family [::SQQ::normalize_family $family]
    set targets [::SQQ::parse_targets $family $values]
    set value [::SQQ::color_value $color]
    foreach target $targets {
        lassign $target source key
        if {$key eq "*"} {
            ::SQQ::clear_family_colors $family
            if {$value ne "default"} { set color_overrides($family,*) $value }
        } else {
            set color_overrides($source,$key) $value
        }
    }
    ::SQQ::render_current
    set labels {}
    foreach target $targets { lappend labels [::SQQ::target_label $family $target] }
    puts "SQQ color $family: [join $labels { }] -> $color"
}

proc ::SQQ::track_representation {rep_index} {
    variable molid
    variable representation_names
    if {[catch {mol repname $molid $rep_index} name]} { return }
    if {$name ni $representation_names} { lappend representation_names $name }
}

proc ::SQQ::adopt_initial_representations {} {
    variable molid
    set count [molinfo $molid get numreps]
    for {set rep 0} {$rep < $count} {incr rep} { ::SQQ::track_representation $rep }
}

proc ::SQQ::clear_representations {} {
    variable molid
    variable representation_names
    if {$molid < 0} { return }
    set indexes {}
    foreach name $representation_names {
        if {[catch {mol repindex $molid $name} rep]} { continue }
        if {[string is integer -strict $rep] && $rep >= 0} { lappend indexes $rep }
    }
    foreach rep [lsort -integer -decreasing -unique $indexes] { mol delrep $rep $molid }
    set representation_names {}
}

proc ::SQQ::expanded_targets {frame family targets} {
    variable group_keys
    set expanded {}
    foreach target $targets {
        lassign $target source key
        if {$family in {cage guest} && $source eq $family} {
            set id_source "${family}-id"
            set group_key "$frame,$id_source"
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_id [lsort -unique $group_keys($group_key)] {
                if {$key eq "*" || [::SQQ::cage_type_from_object_id $frame $object_id] eq $key} {
                    set item [list $id_source $object_id]
                    if {$item ni $expanded} { lappend expanded $item }
                }
            }
        } elseif {$key eq "*"} {
            set group_key "$frame,$source"
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_key [::SQQ::ordered_keys $family $group_keys($group_key)] {
                set item [list $source $object_key]
                if {$item ni $expanded} { lappend expanded $item }
            }
        } else {
            set item [list $source $key]
            if {$item ni $expanded} { lappend expanded $item }
        }
    }
    return $expanded
}

proc ::SQQ::compare_render_keys {left right} {
    lassign $left left_priority left_color
    lassign $right right_priority right_color
    if {$left_priority != $right_priority} {
        return [expr {$left_priority < $right_priority ? -1 : 1}]
    }
    return [expr {$left_color < $right_color ? -1 : ($left_color > $right_color)}]
}

proc ::SQQ::add_dynamic_bonds_representation {indexes color_id radius {material Opaque}} {
    variable molid
    if {[llength $indexes] == 0} { return 0 }
    mol representation DynamicBonds 3.5 $radius 12.0
    mol color ColorID $color_id
    mol selection "index [join $indexes { }]"
    mol material $material
    mol addrep $molid
    ::SQQ::track_representation [expr {[molinfo $molid get numreps] - 1}]
    return 1
}

proc ::SQQ::add_guest_representation {indexes color_id {material Opaque}} {
    variable molid
    if {[llength $indexes] == 0} { return 0 }
    mol representation CPK 1.0 0.3 12.0 12.0
    mol color ColorID $color_id
    mol selection "index [join $indexes { }]"
    mol material $material
    mol addrep $molid
    ::SQQ::track_representation [expr {[molinfo $molid get numreps] - 1}]
    return 1
}

proc ::SQQ::clear_graphics {} {
    variable graphics_ids
    variable molid
    if {$molid >= 0} {
        foreach graphics_id $graphics_ids {
            catch {graphics $molid delete $graphics_id}
        }
    }
    set graphics_ids {}
}

proc ::SQQ::object_center {indexes frame} {
    variable molid
    if {[llength $indexes] == 0} { return {} }
    if {[catch {set selection [atomselect $molid "index [join $indexes { }]" frame $frame]}]} {
        return {}
    }
    if {[catch {set coordinates [$selection get {x y z}]}]} {
        catch {$selection delete}
        return {}
    }
    catch {$selection delete}
    if {[llength $coordinates] == 0} { return {} }
    set x 0.0
    set y 0.0
    set z 0.0
    foreach coordinate $coordinates {
        lassign $coordinate cx cy cz
        set x [expr {$x + $cx}]
        set y [expr {$y + $cy}]
        set z [expr {$z + $cz}]
    }
    set count [expr {double([llength $coordinates])}]
    return [list [expr {$x / $count}] [expr {$y / $count}] [expr {$z / $count}]]
}

proc ::SQQ::render_labels {frame} {
    variable graphics_ids
    variable group_atoms
    variable group_keys
    variable label_visible
    variable molid
    variable pick_mode
    if {!$label_visible && $pick_mode ne "center"} { return }
    set group_key "$frame,cage-id"
    if {![info exists group_keys($group_key)]} { return }
    foreach object_id [lsort -dictionary -unique $group_keys($group_key)] {
        set atom_key "$frame,cage-id,$object_id"
        if {![info exists group_atoms($atom_key)]} { continue }
        set center [::SQQ::object_center $group_atoms($atom_key) $frame]
        if {[llength $center] != 3} { continue }
        if {$pick_mode eq "center"} {
            catch {graphics $molid color yellow}
            if {![catch {graphics $molid sphere $center radius 0.35 resolution 10} graphics_id]} {
                lappend graphics_ids $graphics_id
            }
        }
        catch {graphics $molid color black}
        if {![catch {graphics $molid text $center $object_id size 1.0 thickness 1.0} graphics_id]} {
            lappend graphics_ids $graphics_id
        }
    }
}

proc ::SQQ::find_pick_target {frame atom_index mode} {
    variable group_atoms
    variable group_keys
    set source [expr {$mode eq "guest" ? "guest-id" : "cage-id"}]
    set group_key "$frame,$source"
    if {![info exists group_keys($group_key)]} { return {} }
    set candidates {}
    foreach object_id $group_keys($group_key) {
        set atom_key "$frame,$source,$object_id"
        if {[info exists group_atoms($atom_key)] &&
            [lsearch -integer -exact $group_atoms($atom_key) $atom_index] >= 0} {
            lappend candidates $object_id
        }
    }
    set candidates [lsort -dictionary -unique $candidates]
    if {[llength $candidates] == 0} { return {} }
    if {[llength $candidates] > 1} {
        puts "SQQ pick: atom $atom_index matches [join $candidates { }]; selected [lindex $candidates 0]"
    }
    return [list cage-id [lindex $candidates 0]]
}

proc ::SQQ::pick_atom_changed {name1 name2 operation} {
    variable molid
    variable pick_mode
    variable selected_target
    if {$pick_mode eq "off" || $molid < 0} { return }
    if {![info exists ::vmd_pick_atom] ||
        ![string is integer -strict $::vmd_pick_atom] || $::vmd_pick_atom < 0} {
        return
    }
    if {[info exists ::vmd_pick_mol] && $::vmd_pick_mol != $molid} { return }
    set frame [molinfo $molid get frame]
    set target [::SQQ::find_pick_target $frame $::vmd_pick_atom $pick_mode]
    if {[llength $target] != 2} {
        puts "SQQ pick: atom $::vmd_pick_atom has no $pick_mode target in frame $frame"
        return
    }
    set selected_target $target
    ::SQQ::render_current
    puts "SQQ selected: [lindex $selected_target 1] (frame $frame)"
}

proc ::SQQ::render_selected {frame} {
    variable group_atoms
    variable pick_mode
    variable selected_target
    if {$pick_mode eq "off" || [llength $selected_target] != 2} { return }
    lassign $selected_target source key
    set atom_key "$frame,cage-id,$key"
    if {![info exists group_atoms($atom_key)]} { return }
    lassign [::SQQ::effective_color $frame cage-id $key] color_id priority
    set indexes [lsort -integer -unique $group_atoms($atom_key)]
    ::SQQ::add_dynamic_bonds_representation $indexes $color_id 0.250 Opaque
    if {$pick_mode eq "guest"} {
        set guest_key "$frame,guest-id,$key"
        if {[info exists group_atoms($guest_key)]} {
            set guest_indexes [lsort -integer -unique $group_atoms($guest_key)]
            ::SQQ::add_guest_representation $guest_indexes $color_id Opaque
        }
    }
}

proc ::SQQ::announce_graph_mode {frame} {
    variable graph_mode
    variable displayed_graph_mode
    set value [expr {[info exists graph_mode($frame)] ? $graph_mode($frame) : "unknown"}]
    if {$value ne $displayed_graph_mode} {
        puts "SQQ graph: $value"
        set displayed_graph_mode $value
    }
}

proc ::SQQ::ordered_active_families {} {
    variable active_families
    set ordered {}
    foreach family {phase cluster domain cage guest} {
        if {$family in $active_families} { lappend ordered $family }
    }
    return $ordered
}

proc ::SQQ::render_family {frame family targets} {
    variable group_atoms
    set representation_count 0
    set material [::SQQ::base_material]
    if {$family in {cage guest}} {
        array set explicit_ids {}
        foreach target $targets {
            lassign $target source key
            if {$source eq "${family}-id"} { set explicit_ids($key) 1 }
        }
        array set layer_atoms {}
        set layer_keys {}
        foreach item [::SQQ::expanded_targets $frame $family $targets] {
            lassign $item source key
            set atom_key "$frame,$source,$key"
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $frame $source $key] color_id color_priority
            set explicit [info exists explicit_ids($key)]
            set layer_key [::SQQ::object_render_key $frame $key $color_priority $color_id $explicit]
            if {![info exists layer_atoms($layer_key)]} {
                lappend layer_keys $layer_key
                set layer_atoms($layer_key) {}
            }
            foreach atom_index $group_atoms($atom_key) { lappend layer_atoms($layer_key) $atom_index }
        }
        set layer_keys [lsort -command ::SQQ::compare_object_render_keys $layer_keys]
        if {$family eq "cage"} {
            set radius_tiers {}
            foreach layer_key $layer_keys { lappend radius_tiers [::SQQ::cage_radius_tier $layer_key] }
            set radius_tiers [lsort -integer -unique $radius_tiers]
            foreach layer_key $layer_keys {
                set indexes [lsort -integer -unique $layer_atoms($layer_key)]
                set radius [::SQQ::cage_layer_radius [::SQQ::cage_radius_tier $layer_key] $radius_tiers]
                incr representation_count [::SQQ::add_dynamic_bonds_representation $indexes [lindex $layer_key 6] $radius $material]
            }
        } else {
            foreach layer_key $layer_keys {
                set indexes [lsort -integer -unique $layer_atoms($layer_key)]
                incr representation_count [::SQQ::add_guest_representation $indexes [lindex $layer_key 6] $material]
            }
        }
    } else {
        array set color_atoms {}
        set render_keys {}
        foreach item [::SQQ::expanded_targets $frame $family $targets] {
            lassign $item source key
            set atom_key "$frame,$source,$key"
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $frame $source $key] color_id priority
            set render_key "$priority,$color_id"
            if {![info exists color_atoms($render_key)]} { lappend render_keys [list $priority $color_id] }
            foreach atom_index $group_atoms($atom_key) { lappend color_atoms($render_key) $atom_index }
        }
        foreach render_key [lsort -command ::SQQ::compare_render_keys $render_keys] {
            lassign $render_key priority color_id
            set indexes [lsort -integer -unique $color_atoms($priority,$color_id)]
            incr representation_count [::SQQ::add_dynamic_bonds_representation $indexes $color_id 0.125 $material]
        }
    }
    return $representation_count
}

proc ::SQQ::render_current {} {
    ::SQQ::cancel_pending_render
    variable molid
    variable active_families
    variable active_targets
    variable pick_mode
    if {$molid < 0 || $molid ni [molinfo list]} { return }
    set frame [molinfo $molid get frame]
    ::SQQ::announce_graph_mode $frame
    ::SQQ::clear_representations
    ::SQQ::clear_graphics
    foreach family [::SQQ::ordered_active_families] {
        set targets $active_targets($family)
        set representation_count [::SQQ::render_family $frame $family $targets]
        set labels {}
        foreach target $targets { lappend labels [::SQQ::target_label $family $target] }
        if {$representation_count == 0} {
            puts "SQQ show $family: no memberships for [join $labels { }] in frame $frame"
        } else {
            puts "SQQ show $family: [join $labels { }] (frame $frame)"
        }
    }
    if {$pick_mode eq "guest" && "guest" ni $active_families} {
        ::SQQ::render_family $frame guest [list [list guest *]]
    }
    ::SQQ::render_selected $frame
    ::SQQ::render_labels $frame
    display update
}

proc ::SQQ::cancel_pending_render {} {
    variable frame_after_id
    if {$frame_after_id ne ""} {
        catch {after cancel $frame_after_id}
        set frame_after_id ""
    }
}

proc ::SQQ::render_pending {} {
    variable frame_after_id
    set frame_after_id ""
    ::SQQ::render_current
}

proc ::SQQ::frame_changed {name1 name2 operation} {
    variable molid
    variable frame_after_id
    variable selected_target
    if {$name2 ne "$molid"} { return }
    set selected_target {}
    if {$frame_after_id ne ""} { catch {after cancel $frame_after_id} }
    set frame_after_id [after idle [list ::SQQ::render_pending]]
}

proc ::SQQ::startup_help {} {
    puts "SQQ VMD Renderer"
    puts ""
    puts "Default view : cage all (opaque)"
    puts "Show mode    : additive"
    puts ""
    puts "Commands:"
    puts "  sqq show <family> <target...> ?<family> <target...> ...?"
    puts "  sqq show label ?on|off?"
    puts "  sqq color <family> <target...> <color>"
    puts "  sqq pick center|guest|off"
    puts "  sqq clear"
    puts "  sqq -h"
    puts ""
    puts "Examples:"
    puts "  sqq show cage 512"
    puts "  sqq show cage 512 guest 512"
    puts "  sqq pick guest"
}

proc ::SQQ::help {} {
    puts "SQQ VMD commands"
    puts ""
    puts "Usage:"
    puts "  sqq show <family> <target...> ?<family> <target...> ...?"
    puts "  sqq show label ?on|off?       (lable is accepted)"
    puts "  sqq color <family> <target...> <VMD-color|ColorID|default>"
    puts "  sqq pick center|guest|off"
    puts "  sqq clear"
    puts "  sqq help | sqq -h | sqq --help"
    puts ""
    puts "Families:"
    puts "  cage      Cage topology or exact cage ID"
    puts "  guest     Guests assigned to a cage topology or exact cage ID"
    puts "  phase     sI, sII, sH, boundary, ambiguous, unclassified, isolated"
    puts "  cluster   Exact cluster ID"
    puts "  domain    Exact domain ID"
    puts ""
    puts "Interaction:"
    puts "  Labels are off by default; sqq show label toggles their state."
    puts "  Pick mode makes objects transparent; the selected cage is opaque."
    puts "  sqq clear restores the source-time cage-all opaque view."
    puts ""
    puts "Examples:"
    puts "  sqq show cage all"
    puts "  sqq show cage 512 51264"
    puts "  sqq show cage 512 guest 512"
    puts "  sqq show label"
    puts "  sqq pick center"
    puts "  sqq color cage 512 green"
    puts "  sqq clear"
}

proc sqq {{command help} args} {
    switch -- [string tolower $command] {
        show {
            if {[llength $args] >= 1 &&
                [string tolower [lindex $args 0]] in {label lable}} {
                ::SQQ::set_label [lrange $args 1 end]
            } else {
                if {[llength $args] < 2} {
                    error "Usage: sqq show <family> <target...> ?<family> <target...> ...?"
                }
                ::SQQ::set_show $args
            }
        }
        color {
            if {[llength $args] < 3} {
                error "Usage: sqq color <family> <target> ?target ...? <VMD-color|ColorID|default>"
            }
            ::SQQ::set_colors [lindex $args 0] [lrange $args 1 end-1] [lindex $args end]
        }
        pick {
            if {[llength $args] != 1} { error "Usage: sqq pick center|guest|off" }
            ::SQQ::set_pick_mode [lindex $args 0]
        }
        clear {
            if {[llength $args] != 0} { error "Usage: sqq clear" }
            ::SQQ::reset_show 1
        }
        help - -h - --help {
            if {[llength $args] != 0} { error "Usage: sqq help" }
            ::SQQ::help
        }
        default { error "Unknown SQQ command '$command'; use show, color, pick, clear, or help" }
    }
}

set script_dir [file dirname [file normalize [info script]]]
set ::SQQ::gro_path [file join $script_dir __SQQ_GRO_FILENAME__]
set ::SQQ::xtc_path [file join $script_dir __SQQ_XTC_FILENAME__]
set ::SQQ::membership_path [file join $script_dir __SQQ_MEMBERSHIP_FILENAME__]
foreach {label path} [list GRO $::SQQ::gro_path XTC $::SQQ::xtc_path membership $::SQQ::membership_path] {
    if {![file isfile $path]} { error "SQQ $label file not found: $path" }
}
set parsed_frames [::SQQ::read_membership_tsv $::SQQ::membership_path]
if {$parsed_frames == 0} { error "SQQ membership TSV contains no frames: $::SQQ::membership_path" }
set ::SQQ::molid [mol new $::SQQ::gro_path type gro waitfor all]
if {$parsed_frames > 1} {
    mol addfile $::SQQ::xtc_path type xtc first 1 waitfor all molid $::SQQ::molid
}
set loaded_frames [molinfo $::SQQ::molid get numframes]
molinfo $::SQQ::molid set frame 0
mol rename $::SQQ::molid __SQQ_MOLECULE_NAME__
if {$loaded_frames != $parsed_frames} {
    error "SQQ frame count mismatch: parsed $parsed_frames, VMD loaded $loaded_frames"
}
::SQQ::adopt_initial_representations
display projection Orthographic
catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
trace add variable ::vmd_frame write ::SQQ::frame_changed
catch {trace remove variable ::vmd_pick_atom write ::SQQ::pick_atom_changed}
trace add variable ::vmd_pick_atom write ::SQQ::pick_atom_changed
::SQQ::startup_help
::SQQ::reset_show
catch {color Display Background white}
'''
