from __future__ import annotations

"""Single render lifecycle shared by Analyze and Track."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from ...models import FrameResult
from .bundle import (
    _cleanup_bundle,
    _cleanup_fragment_workspace,
    _finalize_bundle,
    _prepare_fragment_workspace,
)
from .frame import write_sqq_cage_fragment
from .models import RenderBundle, RenderFragment, RenderSpec

if TYPE_CHECKING:
    from ...models.tracking import TargetSelection, TrackingResult


@dataclass
class RenderSession:
    """Stage frame fragments and publish exactly one complete render package."""

    root: Path
    spec: RenderSpec
    fragment_dir: Path
    _fragments: list[RenderFragment] = field(default_factory=list, repr=False)
    _state: str = field(default="active", repr=False)

    @classmethod
    def create(
        cls,
        outdir: Path,
        spec: RenderSpec | None = None,
    ) -> "RenderSession":
        root = Path(outdir)
        render_spec = spec or RenderSpec()
        return cls(
            root=root,
            spec=render_spec,
            fragment_dir=_prepare_fragment_workspace(root),
        )

    @classmethod
    def cleanup_output(
        cls,
        outdir: Path,
        *,
        fragment_dir: Path | None = None,
    ) -> None:
        """Remove render files and abandoned workspaces owned by one output root."""
        _cleanup_bundle(Path(outdir), fragment_dir=fragment_dir)

    @staticmethod
    def cleanup_workspace(fragment_dir: Path) -> bool:
        """Best-effort cleanup used by runtime failure handling."""
        return _cleanup_fragment_workspace(Path(fragment_dir))

    def stage_frame(self, result: FrameResult, frame_index: int) -> RenderFragment:
        """Write and register one worker-safe frame fragment."""
        self._require_active()
        if (
            self.spec.effective_graph_mode is not None
            and str(result.graph.mode) != str(self.spec.effective_graph_mode)
        ):
            raise ValueError(
                "Render frame graph mode does not match RenderSpec: "
                f"{result.graph.mode} != {self.spec.effective_graph_mode}."
            )
        fragment = write_sqq_cage_fragment(
            result,
            self.fragment_dir,
            frame_index,
            requested_graph_mode=self.spec.requested_graph_mode,
            atom_scope=self.spec.atom_scope,
            component_config=self.spec.component_roles,
        )
        self._fragments.append(fragment)
        return fragment

    def finalize(
        self,
        tracking: TrackingResult | TargetSelection | None = None,
    ) -> RenderBundle:
        """Atomically publish the four-file package and close the session."""
        self._require_active()
        try:
            bundle = _finalize_bundle(
                self.root,
                # Worker processes register fragments through the shared
                # workspace rather than this process-local list.
                fragments=tuple(self._fragments) or None,
                fragment_dir=self.fragment_dir,
                tracking=tracking,
                names=self.spec.names,
                render_kind=self.spec.kind,
                molecule_name=self.spec.molecule_name,
            )
        except Exception:
            self._state = "aborted"
            raise
        self._state = "finalized"
        return bundle

    def abort(self) -> None:
        """Discard unpublished fragments without touching a previous package."""
        if self._state != "active":
            return
        _cleanup_fragment_workspace(self.fragment_dir)
        self._state = "aborted"

    @property
    def state(self) -> str:
        return self._state

    def _require_active(self) -> None:
        if self._state != "active":
            raise RuntimeError(f"Render session is already {self._state}.")


__all__ = ["RenderSession"]
