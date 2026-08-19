"""Build authoritative final-screen statistics from a completed workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from ..citation import completed_citation_evidence
from ..config import normalize_output_types


def completed_run_statistics(
    run_info: Mapping[str, Any],
    config: Mapping[str, Any],
    *,
    result_path: str | Path,
    requested_frames: int,
    analysis_seconds: float,
    write_seconds: float,
    total_seconds: float,
    track: bool = False,
) -> dict[str, Any]:
    """Return final counters, timings, outputs, and executed feature flags."""
    successful = int(run_info.get("frames_ok", 0) or 0)
    failed = int(run_info.get("frames_failed", 0) or 0)
    output_types = normalize_output_types(
        config.get("output", {}).get("types")  # type: ignore[union-attr]
    )
    citation_evidence = completed_citation_evidence(
        config,
        successful_frames=successful,
        completed_outputs=output_types,
        track=track,
        occupancy_evaluated=bool(
            run_info.get(
                "occupancy_evaluated",
                run_info.get("guest_molecules", 0),
            )
        ),
    )
    if "executed_features" in run_info:
        for key in (
            "successful_frames",
            "executed_features",
            "executed_order_parameters",
            "completed_outputs",
            "occupancy_evaluated",
            "guest_residence_evaluated",
        ):
            if key in run_info:
                citation_evidence[key] = run_info[key]
    return {
        "requested_frames": max(0, int(requested_frames)),
        "analyzed_frames": successful + failed,
        "successful_frames": successful,
        "failed_frames": failed,
        "total_seconds": max(0.0, float(total_seconds)),
        "analysis_seconds": max(0.0, float(analysis_seconds)),
        "write_seconds": max(0.0, float(write_seconds)),
        "status": str(run_info.get("status", "completed")),
        "result_path": str(Path(result_path).resolve()),
        **citation_evidence,
    }


__all__ = ["completed_run_statistics"]
