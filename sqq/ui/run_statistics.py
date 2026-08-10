"""Build authoritative final-screen statistics from a completed workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from ..config import normalize_order_parameters, normalize_output_types


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
    order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters")  # type: ignore[union-attr]
    )
    features = _executed_features(config, successful, output_types, track=track)
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
        "executed_features": features,
        "executed_order_parameters": order_parameters if successful else (),
        "completed_outputs": output_types,
    }


def _executed_features(
    config: Mapping[str, Any],
    successful: int,
    outputs: tuple[str, ...],
    *,
    track: bool,
) -> dict[str, bool]:
    ran = successful > 0
    half = _section_enabled(config, "half_cage")
    quasi = _section_enabled(config, "quasi_cage")
    cluster = _section_enabled(config, "hydrate_cluster")
    return {
        "water_network": ran,
        "ring_topology": ran,
        "cage_topology": ran,
        "half_cage": ran and half,
        "quasi_cage": ran and quasi,
        "cage_isomer": ran,
        "cage_occupancy": ran,
        "hydrate_phase_domain": ran and cluster,
        "vmd_rendering": "sqq-render" in outputs,
        "cage_tracking": ran and track,
        "cage_lifetime": ran and track,
    }


def _section_enabled(config: Mapping[str, Any], name: str) -> bool:
    section = config.get(name, {})
    return isinstance(section, Mapping) and bool(section.get("enabled", False))


__all__ = ["completed_run_statistics"]
