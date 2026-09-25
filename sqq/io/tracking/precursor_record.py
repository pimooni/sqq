"""Run-private per-frame records for persistent-cage precursor tracing.

A record holds the facts that precursor reconstruction needs from one analyzed
frame: water-oxygen identities, residues, coordinates, the render atom mapping,
water-water edges, and the water sets of rings, half-cages, quasi-cages, and
cages. Raw Track spools records only until all requested persistent IDs first
appear, so their pre-birth prefixes are never analyzed twice.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from ...models import FrameResult
from ..render.frame import visualization_atoms


PRECURSOR_RECORD_FORMAT = "SQQ precursor record"
PRECURSOR_RECORD_VERSION = 1
PRECURSOR_SPOOL_COMPLETE_NAME = ".complete"
PrecursorRecord = dict[str, Any]


def precursor_record_from_result(
    result: FrameResult,
    frame_index: int,
    *,
    atom_scope: str,
) -> PrecursorRecord:
    """Reduce one frame result to the facts needed for precursor tracing."""
    oxygen_indexes = sorted(int(water.oxygen) for water in result.waters)
    atoms = result.frame.atoms
    render_position = {
        int(atom.index): position
        for position, atom in enumerate(
            visualization_atoms(result, atom_scope=atom_scope)
        )
    }
    box = None
    if result.frame.box is not None:
        values = np.asarray(result.frame.box, dtype=float).reshape(-1)
        if (
            len(values) >= 3
            and np.all(np.isfinite(values[:3]))
            and np.all(values[:3] > 0.0)
        ):
            box = [float(value) for value in values[:3]]
    return {
        "format": PRECURSOR_RECORD_FORMAT,
        "version": PRECURSOR_RECORD_VERSION,
        "frame_index": int(frame_index),
        "frame_name": str(result.frame.name),
        "time_ps": (
            None if result.frame.time_ps is None else float(result.frame.time_ps)
        ),
        "box": box,
        "oxygen_index": oxygen_indexes,
        "oxygen_resid": [int(atoms[index].resid) for index in oxygen_indexes],
        "oxygen_xyz": [
            [float(value) for value in atoms[index].xyz] for index in oxygen_indexes
        ],
        "render_index": [render_position.get(index) for index in oxygen_indexes],
        "edges": sorted(
            [int(min(left, right)), int(max(left, right))]
            for left, right in result.graph.edges
        ),
        "rings": [
            sorted(int(node) for node in ring.nodes)
            for rings in result.rings.values()
            for ring in rings
        ],
        "half_cages": [
            sorted(int(index) for index in patch.waters) for patch in result.half_cages
        ],
        "quasi_cages": [
            sorted(int(index) for index in patch.waters) for patch in result.quasi_cages
        ],
        "cages": [
            sorted(int(index) for index in cage.waters)
            for cage in (result.all_cages or result.cages)
        ],
    }


def precursor_record_path(directory: str | Path, frame_index: int) -> Path:
    return Path(directory) / f"frame_{int(frame_index):09d}.json"


def precursor_spool_is_complete(directory: str | Path) -> bool:
    """Return whether ordered Track no longer needs later precursor records."""
    return (Path(directory) / PRECURSOR_SPOOL_COMPLETE_NAME).is_file()


def mark_precursor_spool_complete(directory: str | Path) -> Path:
    """Atomically stop later frame tasks from extending a precursor spool."""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    marker = root / PRECURSOR_SPOOL_COMPLETE_NAME
    temporary = marker.with_name(f".{marker.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text("complete\n", encoding="ascii")
        os.replace(temporary, marker)
    finally:
        temporary.unlink(missing_ok=True)
    return marker


def write_precursor_record(directory: str | Path, record: PrecursorRecord) -> Path:
    """Atomically write one precursor record into the run-private spool."""
    path = precursor_record_path(directory, int(record["frame_index"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(
                record, ensure_ascii=True, separators=(",", ":"), allow_nan=False
            ),
            encoding="ascii",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def read_precursor_record(
    directory: str | Path,
    frame_index: int,
) -> PrecursorRecord | None:
    """Return the spooled record of one frame, or ``None`` when it is absent."""
    path = precursor_record_path(directory, frame_index)
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid SQQ precursor record: {path}") from exc
    if (
        not isinstance(record, dict)
        or record.get("format") != PRECURSOR_RECORD_FORMAT
        or record.get("version") != PRECURSOR_RECORD_VERSION
        or int(record.get("frame_index", -1)) != int(frame_index)
    ):
        raise ValueError(f"Invalid SQQ precursor record: {path}")
    return record


__all__ = [
    "PRECURSOR_RECORD_FORMAT",
    "PRECURSOR_SPOOL_COMPLETE_NAME",
    "PRECURSOR_RECORD_VERSION",
    "PrecursorRecord",
    "mark_precursor_spool_complete",
    "precursor_record_from_result",
    "precursor_record_path",
    "precursor_spool_is_complete",
    "read_precursor_record",
    "write_precursor_record",
]
