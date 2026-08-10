"""Worker-safe entry point for one analyzed frame."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from ..models import Frame


def process_frame(
    frame_index: int,
    frame: Frame,
    config: dict[str, Any],
    outdir: Path,
    strict: bool,
    stage_callback: Callable[[str], None] | None = None,
    *,
    separated_output: bool = False,
    fragment_dir: Path | None = None,
) -> dict[str, Any]:
    """Process one frame through the 0.5.x compatibility implementation."""
    from ..pipeline import process_frame as implementation

    return implementation(
        frame_index,
        frame,
        config,
        outdir,
        strict,
        stage_callback,
        separated_output=separated_output,
        fragment_dir=fragment_dir,
    )


__all__ = ["process_frame"]
