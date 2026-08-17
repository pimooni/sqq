"""Data contracts for SQQ render packages."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping


FRAGMENT_DIRECTORY = ".sqq-cage-fragments"
FRAGMENT_DIRECTORY_GLOB = f"{FRAGMENT_DIRECTORY}-*"
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
class RenderNames:
    """Names of one complete four-file render package."""

    directory: str = SQQ_RENDER_DIRECTORY
    topology: str = SQQ_CAGE_GRO_NAME
    trajectory: str = SQQ_CAGE_XTC_NAME
    membership: str = SQQ_CAGE_MEMBERSHIP_NAME
    script: str = SQQ_RENDER_SCRIPT_NAME

    def __post_init__(self) -> None:
        for label, value in (
            ("render directory", self.directory),
            ("topology filename", self.topology),
            ("trajectory filename", self.trajectory),
            ("membership filename", self.membership),
            ("script filename", self.script),
        ):
            text = str(value)
            if not text or text in {".", ".."} or Path(text).name != text:
                raise ValueError(f"Invalid SQQ {label}: {value!r}.")

    def paths(self, root: Path) -> tuple[Path, Path, Path, Path, Path]:
        render_dir = Path(root) / self.directory
        return (
            render_dir,
            render_dir / self.topology,
            render_dir / self.trajectory,
            render_dir / self.membership,
            render_dir / self.script,
        )


@dataclass(frozen=True)
class RenderSpec:
    """Requested render behavior shared by Analyze and Track."""

    kind: str = "analyze"
    atom_scope: str = "full"
    names: RenderNames = field(default_factory=RenderNames)
    component_roles: Mapping[str, Any] = field(default_factory=dict)
    requested_graph_mode: str | None = None
    effective_graph_mode: str | None = None
    target: str | None = None
    molecule_name: str = "SQQ cages"

    def __post_init__(self) -> None:
        kind = str(self.kind).strip().lower()
        if kind not in {"analyze", "track"}:
            raise ValueError("Render kind must be analyze or track.")
        scope = str(self.atom_scope).strip().lower()
        if scope not in {"full", "compact"}:
            raise ValueError("Render atom scope must be full or compact.")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "atom_scope", scope)
        object.__setattr__(
            self,
            "component_roles",
            MappingProxyType(dict(self.component_roles)),
        )


@dataclass(frozen=True)
class RenderFragment:
    """One complete annotated GRO frame and its compact manifest."""

    frame_index: int
    gro_path: Path
    manifest_path: Path
    atom_count: int
    atom_signature: str
    effective_graph_mode: str


@dataclass(frozen=True)
class RenderBundle:
    """Files produced by render-package finalization."""

    gro_path: Path | None
    script_path: Path | None
    frame_count: int
    xtc_path: Path | None = None
    membership_path: Path | None = None
    render_dir: Path | None = None

    @property
    def complete(self) -> bool:
        return self.frame_count > 0 and all(
            path is not None
            for path in (
                self.gro_path,
                self.xtc_path,
                self.membership_path,
                self.script_path,
            )
        )


@dataclass(frozen=True)
class RenderFileReference:
    role: str
    path: str
    required: bool = True


@dataclass(frozen=True)
class RenderPackageInspection:
    script_path: Path
    kind: str
    references: tuple[RenderFileReference, ...]
    source: str
    errors: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    empty: tuple[str, ...] = ()
    not_files: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        return not (self.errors or self.missing or self.empty or self.not_files)


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


def _membership_token(value: Any) -> str:
    text = str(value) if value not in {None, ""} else "-"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"SQQ membership token is not ASCII: {value!r}") from exc
    if any(character in text for character in ",:\t\r\n "):
        raise ValueError(f"Invalid SQQ membership token: {value!r}")
    return text


__all__ = [
    "ANNOTATION_COLUMN",
    "ANNOTATION_PREFIX",
    "ATOM_PREFIX_WIDTH",
    "COMPONENT_INDEX_CHUNK",
    "EMPTY_VELOCITY_WIDTH",
    "FRAGMENT_DIRECTORY",
    "FRAGMENT_DIRECTORY_GLOB",
    "LEGACY_SQQ_CAGE_GRO_NAME",
    "LEGACY_SQQ_CAGE_MEMBERSHIP_NAME",
    "LEGACY_SQQ_CAGE_SCRIPT_NAME",
    "LEGACY_SQQ_CAGE_XTC_NAME",
    "LEGACY_SQQ_RENDER_SCRIPT_NAME",
    "MEMBERSHIP_CLASSES",
    "SQQ_CAGE_GRO_NAME",
    "SQQ_CAGE_MEMBERSHIP_NAME",
    "SQQ_CAGE_XTC_NAME",
    "SQQ_RENDER_DIRECTORY",
    "SQQ_RENDER_SCRIPT_NAME",
    "CageMembership",
    "RenderBundle",
    "RenderFileReference",
    "RenderFragment",
    "RenderNames",
    "RenderPackageInspection",
    "RenderSpec",
]
