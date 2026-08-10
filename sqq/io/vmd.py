from __future__ import annotations

"""Run-level annotated GRO and VMD rendering outputs."""

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
from importlib import resources
import json
import os
from pathlib import Path
import re
import shutil
import socket
import warnings
from tempfile import mkdtemp
from time import sleep
from typing import Any, BinaryIO, Iterable, Iterator, Mapping
from uuid import uuid4

import numpy as np

from ..display import graph_mode_display
from ..models import Atom, Frame, FrameResult
from ..vmd_command import (
    RenderFileReference,
    render_manifest_block,
    tcl_help_body,
)
from .gro_grouping import gro_topology_fingerprint
from .gro_writer import ascii_gro_text
from .occupancy import guest_id


FRAGMENT_DIRECTORY = ".sqq-cage-fragments"
FRAGMENT_DIRECTORY_GLOB = f"{FRAGMENT_DIRECTORY}-*"
OUTPUT_LOCK_NAME = ".sqq.lock"
SQQ_RENDER_DIRECTORY = "sqq_render"
SQQ_CAGE_GRO_NAME = "sqq_cage.gro"
SQQ_CAGE_XTC_NAME = "sqq_cage.xtc"
SQQ_CAGE_MEMBERSHIP_NAME = "sqq_cage.membership.tsv"
SQQ_RENDER_SCRIPT_NAME = "sqq_cage.vmd.tcl"
LEGACY_SQQ_CAGE_GRO_NAME = "sqq-cage.gro"
LEGACY_SQQ_CAGE_XTC_NAME = "sqq-cage.xtc"
LEGACY_SQQ_CAGE_MEMBERSHIP_NAME = "sqq-cage.membership.tsv"
LEGACY_SQQ_CAGE_SCRIPT_NAME = "sqq-cage.vmd.tcl"
LEGACY_SQQ_RENDER_SCRIPT_NAME = "sqq-render.vmd.tcl"
ANNOTATION_PREFIX = "; SQQ1 m="
ATOM_PREFIX_WIDTH = 44
EMPTY_VELOCITY_WIDTH = 24
ANNOTATION_COLUMN = ATOM_PREFIX_WIDTH + EMPTY_VELOCITY_WIDTH + 1
COMPONENT_INDEX_CHUNK = 8192
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
    atom_scope: str = "full",
    component_config: Mapping[str, Any] | None = None,
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

    scope = normalize_render_atom_scope(atom_scope)
    memberships = water_cage_memberships(result)
    output_atoms = visualization_atoms(result, atom_scope=scope)
    guest_memberships = guest_cage_memberships(result)
    cage_centers = _cage_center_manifest_records(result)
    guest_groups = _guest_group_manifest_records(result, output_atoms)
    component_groups = _component_group_manifest_records(
        result,
        output_atoms,
        component_config,
    )
    graph_display = _frame_graph_display(result, requested_graph_mode)
    block = annotated_gro_block(
        result,
        memberships,
        graph_display,
        atoms=output_atoms,
        guest_memberships=guest_memberships,
    )
    signature = atom_signature(output_atoms)
    from ..core.tracking import snapshot_from_frame_result

    manifest: dict[str, Any] = {
        "format": "SQQ cage fragment",
        "version": 1,
        "status": "ok",
        "frame_index": index,
        "frame_name": ascii_gro_text(result.frame.name),
        "time_ps": result.frame.time_ps,
        "atom_count": len(output_atoms),
        "atom_signature": signature,
        "atom_scope": scope,
        "effective_graph_mode": str(result.graph.mode),
        "requested_graph_mode": (
            None if requested_graph_mode is None else str(requested_graph_mode)
        ),
        "graph_mode_display": graph_display,
        "cage_centers": cage_centers,
        "guest_groups": guest_groups,
        "component_groups": component_groups,
        "tracking_snapshot": snapshot_from_frame_result(result, index).to_dict(),
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
    write_tracking: bool = True,
    tracking_config: Mapping[str, Any] | None = None,
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
            _remove_visible_render_outputs(root)
            return SqqCageBundle(None, None, 0, None, None, None)

        records = [_read_fragment_manifest(path) for path in manifests]
        records.sort(key=lambda item: item["frame_index"])
        _validate_fragment_records(records)
        tracking = (
            _tracking_from_fragment_records(records, tracking_config)
            if write_tracking
            else None
        )
        if write_script and not write_gro:
            raise ValueError(
                "sqq-render requires its topology GRO."
            )

        if write_script:
            render_dir.mkdir(parents=True, exist_ok=True)
            _write_render_data(
                gro_path,
                xtc_path,
                membership_path,
                records,
            )
            if tracking is not None:
                from .tracking import rewrite_membership_track_ids

                rewrite_membership_track_ids(membership_path, tracking)
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

        if tracking is not None:
            from .tracking import write_tracking_result, write_tracking_tables

            tracking_root = root / "track"
            write_tracking_tables(tracking, tracking_root)
            write_tracking_result(tracking, tracking_root / "track_state.json")

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
        if write_tracking:
            _resilient_rmtree(
                root / "track",
                warn=True,
                description="partial tracking output directory",
            )
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


def cleanup_sqq_cage_bundle(
    outdir: Path,
    *,
    fragment_dir: Path | None = None,
) -> None:
    """Remove visible bundle outputs and any abandoned fragments."""
    root = Path(outdir)
    _remove_visible_render_outputs(root)
    _resilient_rmtree(
        root / "track",
        warn=True,
        description="tracking output directory",
    )
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


def _cage_center_manifest_records(result: FrameResult) -> list[dict[str, Any]]:
    """Return wrapped PBC-aware cage centers in VMD's angstrom unit."""
    cages, memberships = _cage_membership_records(result)
    box = None
    if result.frame.box is not None:
        values = np.asarray(result.frame.box, dtype=float).reshape(-1)
        if len(values) >= 3 and np.all(np.isfinite(values[:3])) and np.all(
            values[:3] > 0.0
        ):
            box = values[:3]
    output: list[dict[str, Any]] = []
    for cage in cages:
        center_nm = np.asarray(cage.center, dtype=float)
        if center_nm.shape != (3,) or np.any(~np.isfinite(center_nm)):
            raise ValueError(
                f"Cage {cage.object_id} has an invalid center for SQQ rendering."
            )
        if box is not None:
            center_nm = np.mod(center_nm, box)
        membership = memberships[cage.object_id]
        output.append(
            {
                "cage_id": membership.cage_id,
                "cage_type": membership.cage_type,
                "phase": membership.class_code,
                "domain": membership.domain_id,
                "cluster": membership.cluster_id,
                "center_angstrom": [
                    float(value) * 10.0 for value in center_nm
                ],
            }
        )
    return output


def _guest_group_manifest_records(
    result: FrameResult,
    output_atoms: list[Atom],
) -> list[dict[str, Any]]:
    """Map every complete guest molecule to render-topology atom indexes."""
    output_index = {
        int(atom.index): index for index, atom in enumerate(output_atoms)
    }
    records: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    claimed_atoms: dict[int, str] = {}
    for guest in result.guests:
        identifier = _membership_token(guest_id(guest))
        if identifier in seen_ids:
            raise ValueError(
                f"Duplicate guest id in SQQ cage output: {identifier}"
            )
        seen_ids.add(identifier)
        indexes: list[int] = []
        for source_index in guest.atoms:
            atom_index = int(source_index)
            if atom_index not in output_index:
                raise ValueError(
                    f"Guest {identifier} atom {atom_index} is absent from the "
                    "SQQ visualization topology."
                )
            render_index = output_index[atom_index]
            previous = claimed_atoms.setdefault(render_index, identifier)
            if previous != identifier:
                raise ValueError(
                    f"SQQ visualization atom {render_index} belongs to both "
                    f"guests {previous} and {identifier}."
                )
            indexes.append(render_index)
        if not indexes:
            raise ValueError(
                f"Guest {identifier} has no atoms for SQQ visualization."
            )
        records.append(
            {
                "guest_id": identifier,
                "resname": _membership_token(guest.resname),
                "atom_indices": sorted(set(indexes)),
            }
        )
    return records


def _component_group_manifest_records(
    result: FrameResult,
    output_atoms: list[Atom],
    config: Mapping[str, Any] | None,
) -> list[dict[str, Any]]:
    """Group render atoms by scientific role and residue name."""
    values = config if isinstance(config, Mapping) else {}
    component = values.get("component", {})
    component = component if isinstance(component, Mapping) else {}
    role_map = _render_role_map(component.get("role_map", {}))
    additive = _render_resnames(values.get("additive", {}))
    environment = _render_resnames(values.get("environment", {}))
    unknown_role = _render_role(component.get("unknown_role", "other"))

    water_atoms = {
        int(atom_index)
        for water in result.waters
        for atom_index in water.atoms
    }
    guest_atoms = {
        int(atom_index)
        for guest in result.guests
        for atom_index in guest.atoms
    }
    groups: dict[tuple[str, str], list[int]] = {}
    for render_index, atom in enumerate(output_atoms):
        source_index = int(atom.index)
        resname = str(atom.resname)
        normalized_resname = resname.strip().upper()
        if source_index in water_atoms:
            role = "water"
        elif source_index in guest_atoms:
            role = "guest"
        elif normalized_resname in role_map:
            role = role_map[normalized_resname]
        elif normalized_resname in additive:
            role = "additive"
        elif normalized_resname in environment:
            role = "environment"
        else:
            role = unknown_role
        groups.setdefault((role, _membership_token(resname)), []).append(
            render_index
        )
    return [
        {
            "role": role,
            "resname": resname,
            "atom_indices": indexes,
        }
        for (role, resname), indexes in sorted(groups.items())
    ]


def _render_resnames(section: Any) -> set[str]:
    if not isinstance(section, Mapping):
        return set()
    raw = section.get("resnames", section.get("resname", ()))
    if raw in (None, ""):
        return set()
    if isinstance(raw, str):
        items = [item.strip() for item in raw.split(",") if item.strip()]
    else:
        try:
            items = [str(item).strip() for item in raw if str(item).strip()]
        except TypeError as exc:
            raise ValueError("Component residue names must be a list.") from exc
    return {item.upper() for item in items}


def _render_role_map(value: Any) -> dict[str, str]:
    """Accept residue-to-role and role-to-residue-list mappings."""
    if value in (None, ""):
        return {}
    if not isinstance(value, Mapping):
        raise ValueError("component.role_map must be a mapping.")
    output: dict[str, str] = {}
    for key, raw in value.items():
        key_text = str(key).strip()
        if isinstance(raw, str):
            role = _render_role(raw)
            resnames = (key_text,)
        else:
            role = _render_role(key_text)
            try:
                resnames = tuple(str(item).strip() for item in raw)
            except TypeError as exc:
                raise ValueError(
                    "component.role_map values must be roles or residue-name lists."
                ) from exc
        for resname in resnames:
            if not resname:
                continue
            normalized = resname.upper()
            previous = output.setdefault(normalized, role)
            if previous != role:
                raise ValueError(
                    f"component.role_map assigns {resname!r} to both "
                    f"{previous} and {role}."
                )
    return output


def _render_role(value: Any) -> str:
    role = str(value).strip().lower()
    supported = {"water", "guest", "additive", "environment", "other"}
    if role not in supported:
        raise ValueError(
            "Render component role must be water, guest, additive, "
            "environment, or other."
        )
    return role


def normalize_render_atom_scope(value: Any) -> str:
    """Normalize the atom selection used by the render topology."""
    scope = str(value).strip().lower()
    if scope not in {"full", "compact"}:
        raise ValueError("render.atom_scope must be full or compact.")
    return scope


def visualization_atoms(
    result: FrameResult,
    atom_scope: str = "compact",
) -> list[Atom]:
    """Return stable render atoms for a full or compact topology."""
    scope = normalize_render_atom_scope(atom_scope)
    if scope == "full":
        return list(result.frame.atoms)
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
    render_kind: str = "analyze",
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
    script = _VMD_SCRIPT.replace(
        "__SQQ_RENDER_MANIFEST__",
        render_manifest_block(
            kind=render_kind,
            files=(
                RenderFileReference("topology", str(gro_filename)),
                RenderFileReference("trajectory", str(xtc_filename)),
                RenderFileReference("membership", str(membership_filename)),
            ),
        ),
    ).replace("__SQQ_HELP_BODY__", tcl_help_body())
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


def _tracking_from_fragment_records(
    records: list[dict[str, Any]],
    config: Mapping[str, Any] | None,
):
    """Build persistent cage identities from complete Analyze fragments."""
    raw = [record.get("tracking_snapshot") for record in records]
    if all(value is None for value in raw):
        return None
    if any(value is None for value in raw):
        raise ValueError("Some SQQ cage fragments lack tracking snapshots.")

    indexes = [int(record["frame_index"]) for record in records]
    if indexes != list(range(len(records))):
        warnings.warn(
            "Persistent cage IDs were not written because one or more selected "
            "frames failed before render finalization.",
            UserWarning,
            stacklevel=2,
        )
        return None

    from ..core.tracking import (
        TrackFrameSnapshot,
        TrackingConfig,
        track_snapshots,
    )

    supported = {
        "min_jaccard",
        "min_shared_fraction",
        "min_shared_waters",
        "max_center_distance_nm",
        "gap_frame",
        "guest_tiebreak",
    }
    values = {
        key: value
        for key, value in dict(config or {}).items()
        if key in supported
    }
    settings = TrackingConfig.from_mapping(values)
    snapshots = tuple(
        TrackFrameSnapshot.from_dict(value)
        for value in raw
        if isinstance(value, Mapping)
    )
    if len(snapshots) != len(records):
        raise ValueError("Invalid tracking snapshot in SQQ cage fragments.")
    return track_snapshots(snapshots, settings)


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


def _write_render_data(
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
            "Writing sqq_cage.xtc requires MDAnalysis."
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
        "cage_id\tcage_type\tphase\tdomain\tcluster\tatom_indices\t"
        "center_x_angstrom\tcenter_y_angstrom\tcenter_z_angstrom\n"
    )
    try:
        with temporary.open("w", encoding="ascii", newline="\n") as output:
            output.write(header)
            for render_index, record in enumerate(records):
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
                    for start in range(0, len(atom_indexes), COMPONENT_INDEX_CHUNK):
                        chunk = atom_indexes[start : start + COMPONENT_INDEX_CHUNK]
                        output.write(
                            "\t".join(
                                (
                                    "P",
                                    str(render_index),
                                    "-",
                                    "-",
                                    "-",
                                    "component",
                                    role,
                                    _membership_token(component["resname"]),
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
                groups = _fragment_membership_groups(record)
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


_VMD_SCRIPT_RESOURCE = "sqq_cage.tcl"


def _load_vmd_script_template() -> str:
    """Load the shared Tcl renderer from the packaged resource."""
    resource = resources.files(__package__).joinpath("render").joinpath(
        _VMD_SCRIPT_RESOURCE
    )
    return resource.read_text(encoding="ascii")


_VMD_SCRIPT = _load_vmd_script_template()
