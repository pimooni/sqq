from __future__ import annotations

"""Run-level annotated GRO and VMD rendering outputs."""

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Iterable

import numpy as np

from ..display import graph_mode_display
from ..models import Atom, Frame, FrameResult
from .gro_grouping import gro_topology_fingerprint
from .gro_writer import ascii_gro_text


FRAGMENT_DIRECTORY = ".sqq-cage-fragments"
SQQ_CAGE_GRO_NAME = "sqq-cage.gro"
SQQ_RENDER_SCRIPT_NAME = "sqq-render.vmd.tcl"
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


def prepare_sqq_cage_fragments(outdir: Path) -> Path:
    """Start a clean run-level fragment workspace."""
    root = Path(outdir)
    fragment_dir = root / FRAGMENT_DIRECTORY
    if fragment_dir.exists():
        shutil.rmtree(fragment_dir)
    fragment_dir.mkdir(parents=True, exist_ok=True)
    return fragment_dir


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
    graph_display = _frame_graph_display(result, requested_graph_mode)
    block = annotated_gro_block(result, memberships, graph_display)
    signature = atom_signature(result.frame.atoms)
    manifest: dict[str, Any] = {
        "format": "SQQ cage fragment",
        "version": 1,
        "status": "ok",
        "frame_index": index,
        "frame_name": ascii_gro_text(result.frame.name),
        "time_ps": result.frame.time_ps,
        "atom_count": len(result.frame.atoms),
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
        atom_count=len(result.frame.atoms),
        atom_signature=signature,
        effective_graph_mode=str(result.graph.mode),
    )


def finalize_sqq_cage_bundle(
    outdir: Path,
    fragments: Iterable[SqqCageFragment | Path] | None = None,
    *,
    write_gro: bool = True,
    write_script: bool = True,
    cleanup: bool = True,
) -> SqqCageBundle:
    """Merge sorted fragments and write the run-level VMD helper."""
    root = Path(outdir)
    root.mkdir(parents=True, exist_ok=True)
    fragment_dir = root / FRAGMENT_DIRECTORY
    gro_path = root / SQQ_CAGE_GRO_NAME
    script_path = root / SQQ_RENDER_SCRIPT_NAME
    manifests = _fragment_manifests(fragment_dir, fragments)
    try:
        if not manifests:
            gro_path.unlink(missing_ok=True)
            script_path.unlink(missing_ok=True)
            return SqqCageBundle(None, None, 0)

        records = [_read_fragment_manifest(path) for path in manifests]
        records.sort(key=lambda item: item["frame_index"])
        _validate_fragment_records(records)

        if write_gro:
            _merge_gro_fragments(gro_path, records)
        else:
            gro_path.unlink(missing_ok=True)

        if write_script:
            if not write_gro and not gro_path.exists():
                raise ValueError(
                    "sqq-render requires sqq-cage.gro; enable sqq-cage-gro first."
                )
            _atomic_write_text(script_path, vmd_script_text(), encoding="ascii")
        else:
            script_path.unlink(missing_ok=True)

        return SqqCageBundle(
            gro_path if write_gro else None,
            script_path if write_script else None,
            len(records),
        )
    except Exception:
        gro_path.unlink(missing_ok=True)
        script_path.unlink(missing_ok=True)
        raise
    finally:
        if cleanup and fragment_dir.exists():
            shutil.rmtree(fragment_dir)


def cleanup_sqq_cage_bundle(outdir: Path) -> None:
    """Remove visible bundle outputs and any abandoned fragments."""
    root = Path(outdir)
    (root / SQQ_CAGE_GRO_NAME).unlink(missing_ok=True)
    (root / SQQ_RENDER_SCRIPT_NAME).unlink(missing_ok=True)
    fragment_dir = root / FRAGMENT_DIRECTORY
    if fragment_dir.exists():
        shutil.rmtree(fragment_dir)


def water_cage_memberships(result: FrameResult) -> dict[int, tuple[CageMembership, ...]]:
    """Map water-oxygen atom indexes to every cage membership."""
    cages = list(result.all_cages or result.cages)
    cage_by_id: dict[str, Any] = {}
    for cage in cages:
        if cage.object_id in cage_by_id:
            raise ValueError(f"Duplicate cage id in SQQ cage output: {cage.object_id}")
        cage_by_id[cage.object_id] = cage

    classification = _cage_classification(result, cage_by_id)
    memberships: dict[int, list[CageMembership]] = {}
    for cage in cages:
        class_code, domain_id, cluster_id = classification.get(
            cage.object_id, ("-", "-", "-")
        )
        item = CageMembership(
            cage_id=_compact_object_id(cage.object_id),
            cage_type=cage.cage_type,
            class_code=class_code,
            domain_id=_compact_object_id(domain_id),
            cluster_id=_compact_object_id(cluster_id),
        )
        for oxygen in cage.waters:
            oxygen_index = int(oxygen)
            if oxygen_index < 0 or oxygen_index >= len(result.frame.atoms):
                raise ValueError(
                    f"Cage {cage.object_id} references invalid oxygen index "
                    f"{oxygen_index}."
                )
            memberships.setdefault(oxygen_index, []).append(item)
    return {index: tuple(items) for index, items in memberships.items()}


def annotated_gro_block(
    result: FrameResult,
    memberships: dict[int, tuple[CageMembership, ...]],
    graph_display: str,
) -> str:
    """Return one complete ASCII GRO block with SQQ annotations."""
    title_parts = ["SQQ cage", f"frame={ascii_gro_text(result.frame.name)}"]
    if result.frame.time_ps is not None:
        time_value = float(result.frame.time_ps)
        if not np.isfinite(time_value):
            raise ValueError("SQQ cage GRO time_ps must be finite when provided.")
        title_parts.append(f"time_ps={time_value:.9g}")
    title_parts.append("graph=" + ascii_gro_text(graph_display))
    lines = [" ".join(title_parts), f"{len(result.frame.atoms):5d}"]
    for atom in result.frame.atoms:
        encoded = ",".join(item.encode() for item in memberships.get(int(atom.index), ()))
        lines.append(_annotated_atom_line(atom, encoded or "-"))
    lines.append(_box_line(result.frame.box))
    return "\n".join(lines) + "\n"


def atom_signature(atoms: Iterable[Atom]) -> str:
    """Hash topology-compatible atom and molecule ordering for VMD frames."""
    frame = Frame(name="", atoms=list(atoms))
    return gro_topology_fingerprint(frame)


def vmd_script_text() -> str:
    """Return the self-contained ASCII Tcl renderer."""
    return _VMD_SCRIPT


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




def _annotated_atom_line(atom: Atom, encoded_memberships: str) -> str:
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
    line = prefix + velocity_text + ANNOTATION_PREFIX + annotation
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


def _atomic_write_text(path: Path, text: str, *, encoding: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_text(text, encoding=encoding, newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


_VMD_SCRIPT = r'''# SQQ annotated cage renderer for VMD.
namespace eval ::SQQ {
    catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
    if {[info exists frame_after_id] && $frame_after_id ne ""} {
        catch {after cancel $frame_after_id}
    }
    variable molid -1
    variable gro_path ""
    variable current_family cage
    variable current_targets {{cage cage *}}
    variable representation_names {}
    variable frame_after_id ""
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable color_overrides
    variable known_objects
    variable object_aliases
    variable cage_types
    foreach name {group_keys group_atoms graph_mode color_overrides known_objects object_aliases cage_types} {
        catch {array unset $name}
        array set $name {}
    }
}

proc ::SQQ::add_member {frame mode key atom_index} {
    variable group_keys
    variable group_atoms
    variable known_objects
    set known_objects($mode,$key) 1
    set group_key "$frame,$mode"
    set atom_key "$frame,$mode,$key"
    if {![info exists group_atoms($atom_key)]} {
        lappend group_keys($group_key) $key
        set group_atoms($atom_key) {}
    }
    lappend group_atoms($atom_key) $atom_index
}

proc ::SQQ::register_object {mode key} {
    variable known_objects
    set known_objects($mode,$key) 1
}

proc ::SQQ::register_cage {cage_type object_id} {
    variable object_aliases
    variable cage_types
    ::SQQ::register_object cage $cage_type
    ::SQQ::register_object cage-id $object_id
    set cage_types($object_id) $cage_type
    set compact [string map [list "^" "" "-" "" "_" ""] $cage_type]
    if {$compact ne $cage_type} {
        set object_aliases(cage,$compact) $cage_type
        set suffix [string range $object_id [string length $cage_type] end]
        set object_aliases(cage-id,${compact}${suffix}) $object_id
    }
}

proc ::SQQ::deduplicate_frame_memberships {frame} {
    variable group_keys
    variable group_atoms
    foreach mode {phase cluster domain} {
        set group_key "$frame,$mode"
        if {![info exists group_keys($group_key)]} { continue }
        foreach key $group_keys($group_key) {
            set atom_key "$frame,$mode,$key"
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

proc ::SQQ::cage_type_from_object_id {object_id} {
    variable cage_types
    if {[info exists cage_types($object_id)]} {
        return $cage_types($object_id)
    }
    if {[regexp {^(.+)_([0-9]+)$} $object_id -> cage_type number]} {
        return $cage_type
    }
    return ""
}

proc ::SQQ::read_annotations {path} {
    variable group_keys
    variable group_atoms
    variable graph_mode
    variable known_objects
    variable object_aliases
    variable cage_types
    array unset group_keys
    array unset group_atoms
    array unset graph_mode
    array unset known_objects
    array unset object_aliases
    array unset cage_types
    set handle [open $path r]
    fconfigure $handle -encoding ascii -translation auto
    set frame 0
    while {[gets $handle title] >= 0} {
        if {[regexp {graph=(.*)$} $title -> graph_value]} {
            set graph_mode($frame) [string trim $graph_value]
        } else {
            set graph_mode($frame) "unknown"
        }
        if {[gets $handle count_line] < 0} {
            close $handle
            error "Truncated sqq-cage.gro after frame title"
        }
        set atom_count [string trim $count_line]
        if {![string is integer -strict $atom_count] || $atom_count < 0} {
            close $handle
            error "Invalid atom count in sqq-cage.gro: $count_line"
        }
        for {set atom_index 0} {$atom_index < $atom_count} {incr atom_index} {
            if {[gets $handle line] < 0} {
                close $handle
                error "Truncated atom records in sqq-cage.gro frame $frame"
            }
            if {![regexp {; SQQ1 m=(.*)$} $line -> payload] || $payload eq "-"} {
                continue
            }
            foreach membership [split $payload ,] {
                set fields [split $membership :]
                if {[llength $fields] != 5} {
                    close $handle
                    error "Invalid SQQ membership in frame $frame: $membership"
                }
                lassign $fields cage_id cage_type phase domain_id cluster_id
                if {$cage_type ne "-"} {
                    set object_id [::SQQ::cage_object_id $cage_type $cage_id]
                    ::SQQ::register_cage $cage_type $object_id
                    ::SQQ::add_member $frame cage-id $object_id $atom_index
                }
                if {$phase ne "-"} {
                    ::SQQ::add_member $frame phase $phase $atom_index
                }
                if {$cluster_id ne "-"} {
                    set key [::SQQ::numbered_id cluster $cluster_id]
                    ::SQQ::add_member $frame cluster $key $atom_index
                }
                if {$domain_id ne "-"} {
                    set key [::SQQ::numbered_id domain $domain_id]
                    ::SQQ::add_member $frame domain $key $atom_index
                }
            }
        }
        if {[gets $handle box_line] < 0} {
            close $handle
            error "Missing box record in sqq-cage.gro frame $frame"
        }
        ::SQQ::deduplicate_frame_memberships $frame
        incr frame
    }
    close $handle
    return $frame
}

proc ::SQQ::key_rank {mode key} {
    if {$mode eq "cage"} {
        set order {512 51262 51263 51264 435663 51268}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    if {$mode eq "phase"} {
        set order {I II H B A U X}
        set position [lsearch -exact $order $key]
        if {$position >= 0} { return $position }
        return 100
    }
    return 0
}

proc ::SQQ::ordered_keys {mode keys} {
    set decorated {}
    foreach key [lsort -unique $keys] {
        lappend decorated [list [::SQQ::key_rank $mode $key] $key]
    }
    set output {}
    foreach item [lsort -integer -index 0 $decorated] {
        lappend output [lindex $item 1]
    }
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

proc ::SQQ::cage_render_key {object_id color_priority color_id explicit} {
    set cage_type [::SQQ::cage_type_from_object_id $object_id]
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

proc ::SQQ::compare_cage_render_keys {left right} {
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
    set radius [expr {0.125 + 0.005 * $index / double($count - 1)}]
    return [format "%.3f" $radius]
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

proc ::SQQ::color_id {mode key} {
    if {$mode eq "cage"} {
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
    if {$mode eq "phase"} {
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

proc ::SQQ::default_color_id {source key} {
    if {$source eq "cage-id"} {
        set cage_type [::SQQ::cage_type_from_object_id $key]
        if {$cage_type ne ""} { return [::SQQ::color_id cage $cage_type] }
        return 2
    }
    return [::SQQ::color_id $source $key]
}

proc ::SQQ::effective_color {source key} {
    variable color_overrides
    set exact "$source,$key"
    if {[info exists color_overrides($exact)]} {
        if {$color_overrides($exact) eq "default"} {
            return [list [::SQQ::default_color_id $source $key] 3]
        }
        return [list $color_overrides($exact) 3]
    }
    set family $source
    if {$source eq "cage-id"} {
        set family cage
        set cage_type [::SQQ::cage_type_from_object_id $key]
        if {$cage_type ne ""} {
            set type_key "cage,$cage_type"
            if {[info exists color_overrides($type_key)]} {
                if {$color_overrides($type_key) eq "default"} {
                    return [list [::SQQ::default_color_id $source $key] 2]
                }
                return [list $color_overrides($type_key) 2]
            }
        }
    }
    set category_key "$family,*"
    if {[info exists color_overrides($category_key)]} {
        return [list $color_overrides($category_key) 1]
    }
    return [list [::SQQ::default_color_id $source $key] 0]
}

proc ::SQQ::effective_color_id {source key} {
    return [lindex [::SQQ::effective_color $source $key] 0]
}

proc ::SQQ::phase_key {value} {
    set aliases [dict create \
        si I i I sii II ii II sh H h H \
        boundary B b B ambiguous A a A \
        unclassified U u U isolated X x X]
    set key [string tolower $value]
    if {[dict exists $aliases $key]} {
        return [dict get $aliases $key]
    }
    return ""
}

proc ::SQQ::parse_target {value} {
    variable known_objects
    variable object_aliases
    set token [string trim $value]
    set lower [string tolower $token]
    if {$lower eq "all"} {
        return [list cage cage *]
    }
    if {$lower in {cage phase cluster domain}} {
        return [list $lower $lower *]
    }
    set phase [::SQQ::phase_key $token]
    if {$phase ne ""} {
        return [list phase phase $phase]
    }
    foreach source {cage-id cage} {
        if {[info exists known_objects($source,$token)]} {
            return [list cage $source $token]
        }
        if {[info exists object_aliases($source,$token)]} {
            return [list cage $source $object_aliases($source,$token)]
        }
    }
    if {[regexp -nocase {^cluster_([0-9]+)$} $token -> number]} {
        return [list cluster cluster [::SQQ::numbered_id cluster $number]]
    }
    if {[regexp -nocase {^domain_([0-9]+)$} $token -> number]} {
        return [list domain domain [::SQQ::numbered_id domain $number]]
    }
    if {[regexp {^([0-9]+)_([0-9]+)$} $token -> cage_type number]} {
        return [list cage cage-id [::SQQ::cage_object_id $cage_type $number]]
    }
    if {[regexp {^[0-9]+$} $token]} {
        return [list cage cage $token]
    }
    error "Unknown SQQ object '$value'"
}

proc ::SQQ::require_known_target {target} {
    variable known_objects
    lassign $target family source key
    if {$key eq "*" || $source eq "phase"} { return }
    if {![info exists known_objects($source,$key)]} {
        set label [::SQQ::target_label $target]
        error "SQQ object '$label' does not exist in the loaded trajectory"
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

proc ::SQQ::target_label {target} {
    lassign $target family source key
    if {$family eq "cage" && $key eq "*"} { return all }
    if {$key eq "*"} { return $family }
    if {$source eq "phase"} { return [::SQQ::phase_label $key] }
    return $key
}

proc ::SQQ::set_show {values} {
    variable current_family
    variable current_targets
    if {[llength $values] == 0} {
        error "Usage: sqq show <object> ?object ...?"
    }
    set family ""
    set targets {}
    foreach value $values {
        if {[string equal -nocase [string trim $value] cage]} {
            error "'cage' is not a show target; use 'sqq show all'"
        }
        set target [::SQQ::parse_target $value]
        ::SQQ::require_known_target $target
        set target_family [lindex $target 0]
        if {$family eq ""} {
            set family $target_family
        } elseif {$target_family ne $family} {
            error "sqq show cannot mix $family and $target_family objects"
        }
        if {$target ni $targets} { lappend targets $target }
    }
    if {[llength $targets] > 1} {
        foreach target $targets {
            if {[lindex $target 2] eq "*"} {
                error "A category name must be used alone with sqq show"
            }
        }
    }
    set current_family $family
    set current_targets $targets
    ::SQQ::render_current
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
    if {$color_id < 0} {
        error "Unknown VMD color '$value'"
    }
    return $color_id
}

proc ::SQQ::clear_family_colors {family} {
    variable color_overrides
    foreach name [array names color_overrides] {
        set source [lindex [split $name ,] 0]
        if {$source eq $family || ($family eq "cage" && $source eq "cage-id")} {
            unset color_overrides($name)
        }
    }
}

proc ::SQQ::set_color {object color} {
    variable color_overrides
    set target [::SQQ::parse_target $object]
    lassign $target family source key
    set value [::SQQ::color_value $color]
    ::SQQ::require_known_target $target
    if {$key eq "*"} {
        ::SQQ::clear_family_colors $family
        if {$value ne "default"} { set color_overrides($family,*) $value }
    } else {
        set color_overrides($source,$key) $value
    }
    ::SQQ::render_current
    puts "SQQ color: [::SQQ::target_label $target] -> $color"
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
    for {set rep 0} {$rep < $count} {incr rep} {
        ::SQQ::track_representation $rep
    }
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
    foreach rep [lsort -integer -decreasing -unique $indexes] {
        mol delrep $rep $molid
    }
    set representation_names {}
}

proc ::SQQ::expanded_targets {frame targets} {
    variable group_keys
    set expanded {}
    foreach target $targets {
        lassign $target family source key
        if {$family eq "cage" && $source ne "cage-id"} {
            set group_key "$frame,cage-id"
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_id [lsort -unique $group_keys($group_key)] {
                if {$key eq "*" || [string match "${key}_*" $object_id]} {
                    set item [list cage-id $object_id]
                    if {$item ni $expanded} { lappend expanded $item }
                }
            }
        } elseif {$key eq "*"} {
            set group_key "$frame,$source"
            if {![info exists group_keys($group_key)]} { continue }
            foreach object_key [::SQQ::ordered_keys $source $group_keys($group_key)] {
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

proc ::SQQ::add_dynamic_bonds_representation {indexes color_id radius} {
    variable molid
    if {[llength $indexes] == 0} { return 0 }
    mol representation DynamicBonds 3.5 $radius 12.0
    mol color ColorID $color_id
    mol selection "index [join $indexes { }]"
    mol material Opaque
    mol addrep $molid
    set rep_index [expr {[molinfo $molid get numreps] - 1}]
    ::SQQ::track_representation $rep_index
    return 1
}

proc ::SQQ::render_current {} {
    ::SQQ::cancel_pending_render
    variable molid
    variable current_family
    variable current_targets
    variable group_atoms
    variable graph_mode
    if {$molid < 0 || $molid ni [molinfo list]} { return }
    set frame [molinfo $molid get frame]
    if {[info exists graph_mode($frame)]} {
        puts "SQQ graph: $graph_mode($frame)"
    }
    ::SQQ::clear_representations
    set representation_count 0
    if {$current_family eq "cage"} {
        array set explicit_ids {}
        foreach target $current_targets {
            lassign $target family source key
            if {$source eq "cage-id"} { set explicit_ids($key) 1 }
        }
        array set layer_atoms {}
        set layer_keys {}
        foreach item [::SQQ::expanded_targets $frame $current_targets] {
            lassign $item source key
            set atom_key "$frame,$source,$key"
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $source $key] color_id color_priority
            set explicit [info exists explicit_ids($key)]
            set layer_key [::SQQ::cage_render_key $key $color_priority $color_id $explicit]
            if {![info exists layer_atoms($layer_key)]} {
                lappend layer_keys $layer_key
                set layer_atoms($layer_key) {}
            }
            foreach atom_index $group_atoms($atom_key) {
                lappend layer_atoms($layer_key) $atom_index
            }
        }
        set layer_keys [lsort -command ::SQQ::compare_cage_render_keys $layer_keys]
        set radius_tiers {}
        foreach layer_key $layer_keys {
            lappend radius_tiers [::SQQ::cage_radius_tier $layer_key]
        }
        set radius_tiers [lsort -integer -unique $radius_tiers]
        foreach layer_key $layer_keys {
            set color_id [lindex $layer_key 6]
            set indexes [lsort -integer -unique $layer_atoms($layer_key)]
            set tier [::SQQ::cage_radius_tier $layer_key]
            set radius [::SQQ::cage_layer_radius $tier $radius_tiers]
            incr representation_count [::SQQ::add_dynamic_bonds_representation $indexes $color_id $radius]
        }
    } else {
        array set color_atoms {}
        set render_keys {}
        foreach item [::SQQ::expanded_targets $frame $current_targets] {
            lassign $item source key
            set atom_key "$frame,$source,$key"
            if {![info exists group_atoms($atom_key)]} { continue }
            lassign [::SQQ::effective_color $source $key] color_id priority
            set render_key "$priority,$color_id"
            if {![info exists color_atoms($render_key)]} {
                lappend render_keys [list $priority $color_id]
            }
            foreach atom_index $group_atoms($atom_key) {
                lappend color_atoms($render_key) $atom_index
            }
        }
        foreach render_key [lsort -command ::SQQ::compare_render_keys $render_keys] {
            lassign $render_key priority color_id
            set indexes [lsort -integer -unique $color_atoms($priority,$color_id)]
            incr representation_count [::SQQ::add_dynamic_bonds_representation $indexes $color_id 0.125]
        }
    }
    display update
    set labels {}
    foreach target $current_targets { lappend labels [::SQQ::target_label $target] }
    if {$representation_count == 0} {
        puts "SQQ show: no $current_family memberships for [join $labels { }] in frame $frame"
    } else {
        puts "SQQ show: [join $labels { }] (frame $frame)"
    }
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
    if {$name2 ne "$molid"} { return }
    if {$frame_after_id ne ""} { catch {after cancel $frame_after_id} }
    set frame_after_id [after idle [list ::SQQ::render_pending]]
}

proc ::SQQ::split_frames {path directory} {
    file mkdir $directory
    set input [open $path r]
    fconfigure $input -encoding ascii -translation auto
    set files {}
    set frame 0
    while {[gets $input title] >= 0} {
        if {[gets $input count_line] < 0} {
            close $input
            error "Truncated sqq-cage.gro after frame title"
        }
        set atom_count [string trim $count_line]
        if {![string is integer -strict $atom_count] || $atom_count < 0} {
            close $input
            error "Invalid atom count in sqq-cage.gro: $count_line"
        }
        set frame_path [file join $directory [format "frame_%09d.gro" $frame]]
        set output [open $frame_path w]
        fconfigure $output -encoding ascii -translation lf
        puts $output $title
        puts $output $count_line
        for {set atom_index 0} {$atom_index < $atom_count} {incr atom_index} {
            if {[gets $input line] < 0} {
                close $output
                close $input
                error "Truncated atom records in sqq-cage.gro frame $frame"
            }
            puts $output $line
        }
        if {[gets $input box_line] < 0} {
            close $output
            close $input
            error "Missing box record in sqq-cage.gro frame $frame"
        }
        puts $output $box_line
        close $output
        lappend files $frame_path
        incr frame
    }
    close $input
    return $files
}

proc ::SQQ::help {} {
    puts "SQQ commands:"
    puts "  sqq show all"
    puts "  sqq show <category>"
    puts "    category: phase, cluster, domain"
    puts "  sqq show <object> ?object ...?"
    puts "  sqq color <object|category> <VMD-color|ColorID|default>"
    puts "  sqq help"
    puts "Examples:"
    puts "  sqq show all"
    puts "  sqq show 51262"
    puts "  sqq show 512 51262"
    puts "  sqq show 51262_00053"
    puts "  sqq show sI boundary"
    puts "  sqq color 51262 blue"
}

proc sqq {{command help} args} {
    switch -- [string tolower $command] {
        show { ::SQQ::set_show $args }
        color {
            if {[llength $args] != 2} {
                error "Usage: sqq color <object|category> <VMD-color|ColorID|default>"
            }
            ::SQQ::set_color [lindex $args 0] [lindex $args 1]
        }
        help {
            if {[llength $args] != 0} { error "Usage: sqq help" }
            ::SQQ::help
        }
        default { error "Unknown SQQ command '$command'; use show, color, or help" }
    }
}

set ::SQQ::gro_path [file join [file dirname [file normalize [info script]]] "sqq-cage.gro"]
if {![file isfile $::SQQ::gro_path]} {
    error "SQQ GRO file not found: $::SQQ::gro_path"
}
set parsed_frames [::SQQ::read_annotations $::SQQ::gro_path]
set temp_dir [file join [file dirname $::SQQ::gro_path] ".sqq-vmd-[pid]-[clock clicks]"]
set frame_files [::SQQ::split_frames $::SQQ::gro_path $temp_dir]
if {[llength $frame_files] == 0} {
    error "SQQ GRO file contains no frames: $::SQQ::gro_path"
}
set ::SQQ::molid [mol new [lindex $frame_files 0] type gro waitfor all]
foreach frame_path [lrange $frame_files 1 end] {
    mol addfile $frame_path type gro waitfor all molid $::SQQ::molid
}
foreach frame_path $frame_files { file delete -force $frame_path }
file delete -force $temp_dir
mol rename $::SQQ::molid "SQQ cages"
set loaded_frames [molinfo $::SQQ::molid get numframes]
if {$loaded_frames != $parsed_frames} {
    error "SQQ frame count mismatch: parsed $parsed_frames, VMD loaded $loaded_frames"
}
::SQQ::adopt_initial_representations
display projection Orthographic
catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
trace add variable ::vmd_frame write ::SQQ::frame_changed
sqq show all
'''
