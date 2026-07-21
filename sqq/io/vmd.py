from __future__ import annotations

"""Run-level annotated GRO and VMD rendering outputs."""

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Iterable

import numpy as np

from ..display import graph_mode_display
from ..models import Atom, FrameResult
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
    """Hash the atom identity and order represented by an output GRO frame."""
    digest = sha256()
    for atom in atoms:
        record = (
            int(atom.index),
            int(atom.resid),
            ascii_gro_text(atom.resname)[:5],
            ascii_gro_text(atom.atomname)[:5],
            int(atom.atomid),
        )
        digest.update(json.dumps(record, ensure_ascii=True, separators=(",", ":")).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


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
    variable molid -1
    variable gro_path ""
    variable current_mode cage
    variable group_keys
    variable group_atoms
    variable graph_mode
    array set group_keys {}
    array set group_atoms {}
    array set graph_mode {}
}

proc ::SQQ::add_member {frame mode key atom_index} {
    variable group_keys
    variable group_atoms
    set group_key "$frame,$mode"
    set atom_key "$frame,$mode,$key"
    if {![info exists group_atoms($atom_key)]} {
        lappend group_keys($group_key) $key
        set group_atoms($atom_key) {}
    }
    lappend group_atoms($atom_key) $atom_index
}

proc ::SQQ::read_annotations {path} {
    variable group_keys
    variable group_atoms
    variable graph_mode
    array unset group_keys
    array unset group_atoms
    array unset graph_mode
    set handle [open $path r]
    fconfigure $handle -encoding ascii -translation lf
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
                    ::SQQ::add_member $frame cage $cage_type $atom_index
                }
                if {$phase ne "-"} {
                    ::SQQ::add_member $frame phase $phase $atom_index
                }
                if {$cluster_id ne "-"} {
                    ::SQQ::add_member $frame cluster $cluster_id $atom_index
                }
                if {$domain_id ne "-"} {
                    ::SQQ::add_member $frame domain $domain_id $atom_index
                }
            }
        }
        if {[gets $handle box_line] < 0} {
            close $handle
            error "Missing box record in sqq-cage.gro frame $frame"
        }
        incr frame
    }
    close $handle
    return $frame
}

proc ::SQQ::key_rank {mode key} {
    if {$mode eq "cage"} {
        set order {512 51262 51263 51264 51268 435663}
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

proc ::SQQ::clear_representations {} {
    variable molid
    if {$molid < 0} { return }
    set count [molinfo $molid get numreps]
    for {set rep [expr {$count - 1}]} {$rep >= 0} {incr rep -1} {
        mol delrep $rep $molid
    }
}

proc ::SQQ::render {mode} {
    variable molid
    variable current_mode
    variable group_keys
    variable group_atoms
    variable graph_mode
    set mode [string tolower $mode]
    if {$mode eq "help"} {
        puts "SQQ commands: sqq cage | phase | cluster | domain | help"
        return
    }
    if {$mode ni {cage phase cluster domain}} {
        error "Unknown SQQ render mode '$mode'; use: cage, phase, cluster, domain, help"
    }
    set current_mode $mode
    set frame [molinfo $molid get frame]
    if {[info exists graph_mode($frame)]} {
        puts "SQQ graph: $graph_mode($frame)"
    }
    set group_key "$frame,$mode"
    ::SQQ::clear_representations
    if {![info exists group_keys($group_key)]} {
        puts "SQQ render: no $mode memberships in frame $frame"
        display update
        return
    }
    foreach key [::SQQ::ordered_keys $mode $group_keys($group_key)] {
        set atom_key "$frame,$mode,$key"
        set indexes [lsort -integer -unique $group_atoms($atom_key)]
        if {[llength $indexes] == 0} { continue }
        mol representation DynamicBonds 3.5 0.12 12.0
        mol color ColorID [::SQQ::color_id $mode $key]
        mol selection "index [join $indexes { }]"
        mol material Opaque
        mol addrep $molid
    }
    display update
    puts "SQQ render: $mode (frame $frame)"
}

proc ::SQQ::frame_changed {name1 name2 operation} {
    variable molid
    variable current_mode
    if {$name2 ne "$molid"} { return }
    after idle [list ::SQQ::render $current_mode]
}

proc ::SQQ::split_frames {path directory} {
    file mkdir $directory
    set input [open $path r]
    fconfigure $input -encoding ascii -translation lf
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

proc sqq {{mode help} args} {
    if {[llength $args] != 0} {
        error "Usage: sqq cage|phase|cluster|domain|help"
    }
    ::SQQ::render $mode
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
display projection Orthographic
catch {trace remove variable ::vmd_frame write ::SQQ::frame_changed}
trace add variable ::vmd_frame write ::SQQ::frame_changed
sqq cage
'''
