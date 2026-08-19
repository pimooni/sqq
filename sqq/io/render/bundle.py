"""Atomic assembly and publication of complete SQQ render packages."""

from __future__ import annotations

import errno
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
    validate_render_fragment,
)
from .models import (
    ANNOTATION_COLUMN,
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
        _write_render_data(
            stage_gro,
            stage_xtc,
            stage_membership,
            records,
        )
        if tracking is not None:
            from ..tracking import rewrite_membership_track_ids

            rewrite_membership_track_ids(stage_membership, tracking)
        _atomic_write_text(
            stage_script,
            vmd_script_text(
                gro_filename=output_names.topology,
                xtc_filename=output_names.trajectory,
                membership_filename=output_names.membership,
                molecule_name=molecule_name,
                render_kind=render_kind,
            ),
            encoding="ascii",
        )

        _publish_render_directory(stage_dir, render_dir)
        _remove_legacy_root_render_outputs(root)
        return RenderBundle(
            gro_path=gro_path,
            script_path=script_path,
            frame_count=len(records),
            xtc_path=xtc_path,
            membership_path=membership_path,
            render_dir=render_dir,
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
        validate_render_fragment(gro_path, int(record["atom_count"]))


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
            canonical_components: (
                tuple[tuple[str, str, tuple[int, ...]], ...] | None
            ) = None
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
                if canonical_components is None:
                    canonical_components = component_signature
                elif component_signature != canonical_components:
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
