from __future__ import annotations

"""Discovery and publication of render packages for Track targets."""

import csv
import os
from pathlib import Path
import shutil
from uuid import uuid4

from ...models.tracking import TargetSelection
from ..tracking import TRACK_DIRECTORY_NAME, rewrite_membership_track_ids
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
                    frame_count=_membership_frame_count(membership),
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
    raise FileNotFoundError(
        f"Cannot find {SQQ_RENDER_DIRECTORY}/{{{SQQ_CAGE_GRO_NAME}, "
        f"{SQQ_CAGE_XTC_NAME}, {SQQ_CAGE_MEMBERSHIP_NAME}, "
        f"{SQQ_RENDER_SCRIPT_NAME}}}."
    )


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
    frame_count: int,
) -> tuple[Path, Path, Path]:
    """Validate source render files and their correspondence to Track state."""
    paths = _required_render_paths(bundle)
    if bundle.frame_count != int(frame_count):
        raise ValueError(
            "Tracking state and render bundle frame counts differ: "
            f"{frame_count} versus {bundle.frame_count}."
        )
    return paths


def publish_target_render_bundle(
    selection: TargetSelection,
    target_directory: str | Path,
    source_bundle: RenderBundle,
) -> RenderBundle:
    """Publish one self-contained four-file render package for a Track target."""
    gro_source, xtc_source, membership_source = _required_render_paths(source_bundle)
    target_root = Path(target_directory)
    render_dir, gro, xtc, membership, script = TRACK_RENDER_NAMES.paths(target_root)
    render_dir.mkdir(parents=True, exist_ok=True)
    _atomic_link_or_copy(gro_source, gro)
    _atomic_link_or_copy(xtc_source, xtc)
    _atomic_copy(membership_source, membership)
    rewrite_membership_track_ids(membership, selection)
    _atomic_write_text(
        script,
        _target_vmd_script(selection, target_root.name),
        encoding="ascii",
    )

    # Remove obsolete root-level files from pre-package Track layouts.
    (target_root / TRACK_GRO_NAME).unlink(missing_ok=True)
    (target_root / TRACK_TCL_NAME).unlink(missing_ok=True)
    return RenderBundle(
        gro_path=gro,
        script_path=script,
        frame_count=source_bundle.frame_count,
        xtc_path=xtc,
        membership_path=membership,
        render_dir=render_dir,
    )


def _target_vmd_script(selection: TargetSelection, target_name: str) -> str:
    script = vmd_script_text(
        gro_filename=TRACK_GRO_NAME,
        xtc_filename=TRACK_XTC_NAME,
        membership_filename=TRACK_MEMBERSHIP_NAME,
        molecule_name=f"SQQ track {target_name}",
        render_kind="track",
    ).rstrip()
    return (
        script
        + "\n\n# Default lifecycle view for this tracking target.\n"
        + "\n".join(_selection_commands(selection, target_name))
        + "\n"
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
        "sqq show cage " + " ".join(object_ids[start : start + 100])
        for start in range(0, len(object_ids), 100)
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
    return tuple(Path(path).resolve() for path in paths)  # type: ignore[return-value]


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
    "discover_sqq_cage_bundle",
    "discover_sqq_cage_gro",
    "publish_target_render_bundle",
    "validate_tracking_source_bundle",
]
