from __future__ import annotations

"""Per-frame render selection, membership, and fragment generation."""

import json
import os
from pathlib import Path
import re
from typing import Any, Iterable, Mapping

import numpy as np

from ...display import graph_mode_display
from ...models import Atom, Frame, FrameResult
from ..gro_grouping import gro_topology_fingerprint
from ..gro_writer import ascii_gro_text
from ..occupancy import guest_id
from .models import (
    ANNOTATION_COLUMN,
    ANNOTATION_PREFIX,
    ATOM_PREFIX_WIDTH,
    EMPTY_VELOCITY_WIDTH,
    MEMBERSHIP_CLASSES,
    CageMembership,
    RenderFragment,
)


def write_sqq_cage_fragment(
    result: FrameResult,
    fragment_dir: Path,
    frame_index: int,
    requested_graph_mode: str | None = None,
    atom_scope: str = "full",
    component_config: Mapping[str, Any] | None = None,
) -> RenderFragment:
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
        "gro_file": gro_path.name,
    }
    _atomic_write_text(gro_path, block, encoding="ascii")
    _atomic_write_text(
        manifest_path,
        json.dumps(manifest, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
        + "\n",
        encoding="ascii",
    )
    return RenderFragment(
        frame_index=index,
        gro_path=gro_path,
        manifest_path=manifest_path,
        atom_count=len(output_atoms),
        atom_signature=signature,
        effective_graph_mode=str(result.graph.mode),
    )

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


def validate_render_fragment(path: Path, atom_count: int) -> None:
    """Validate one complete worker-written annotated GRO fragment."""
    try:
        lines = Path(path).read_text(encoding="ascii").splitlines()
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

def membership_token(value: Any) -> str:
    """Return one safe public membership token."""
    return _membership_token(value)


def _atomic_write_text(path: Path, text: str, *, encoding: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_text(text, encoding=encoding, newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


__all__ = [
    "annotated_gro_block",
    "atom_signature",
    "guest_cage_memberships",
    "membership_token",
    "normalize_render_atom_scope",
    "visualization_atoms",
    "validate_render_fragment",
    "water_cage_memberships",
    "write_sqq_cage_fragment",
]
