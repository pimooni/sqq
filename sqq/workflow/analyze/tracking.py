"""Persistent-ID tracking sink used by the Analyze workflow."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any
import warnings

from ...core.tracking import TrackingAccumulator, snapshot_from_frame_result
from ...io.tracking import (
    build_track_tables,
    write_track_info,
    write_tracking_result,
    write_tracking_tables,
)
from ...models.tracking import TrackFrameSnapshot, TrackingConfig, TrackingResult
from ...runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ...runtime.session import AnalysisSink


_TRACKING_CONFIG_FIELDS = frozenset(
    {
        "min_jaccard",
        "min_shared_fraction",
        "min_shared_waters",
        "max_center_distance_nm",
        "gap_frame",
        "max_gap_ps",
        "guest_tiebreak",
        "ambiguity_score_margin",
        "near_threshold_tolerance",
    }
)


class AnalyzeTrackingSink(AnalysisSink):
    """Build one bounded-memory persistent-ID stream per topology group."""

    def __init__(
        self,
        configs: Mapping[int | str, Mapping[str, Any]],
    ) -> None:
        self._configs = dict(configs)
        self._accumulators: dict[int | str, TrackingAccumulator] = {}
        self._expected: dict[int | str, int] = {}
        self._failed: dict[int | str, str] = {}
        self.results: dict[int | str, TrackingResult] = {}

    def start(self, plan: RunPlan) -> None:
        self._accumulators = {
            key: TrackingAccumulator(_tracking_config(config))
            for key, config in self._configs.items()
        }
        self._expected = {key: 0 for key in self._configs}
        self._failed.clear()
        self.results.clear()

    def consume(self, task: FrameTask, outcome: TaskOutcome) -> None:
        key: int | str = task.group_key if task.group_key is not None else "run"
        if key not in self._accumulators:
            raise RuntimeError(f"No tracking accumulator for topology group {key!r}.")
        expected = self._expected[key]
        if int(task.frame_index) != expected:
            raise RuntimeError(
                f"Topology group {key!r} tracking frame order is not contiguous: "
                f"expected {expected}, got {task.frame_index}."
            )
        self._expected[key] = expected + 1
        if key in self._failed:
            return
        if outcome.snapshot_error:
            self._failed[key] = (
                f"frame {task.display_name!r} could not be reduced to a "
                f"persistent cage snapshot ({outcome.snapshot_error})"
            )
            return
        snapshot = tracking_snapshot(task, outcome)
        if snapshot is None:
            self._failed[key] = "one or more frames failed analysis"
            return
        try:
            self._accumulators[key].add(snapshot)
        except ValueError as exc:
            # Preserve frame outputs when only persistent linking is invalid.
            self._failed[key] = (
                f"frame {task.display_name!r} could not be linked to the "
                f"persistent cage stream ({exc})"
            )

    def finish(
        self,
        plan: RunPlan,
        outcomes: Sequence[TaskOutcome],
    ) -> None:
        for key, accumulator in self._accumulators.items():
            reason = self._failed.get(key)
            if reason is not None:
                warnings.warn(
                    "Persistent Track state was not written for topology group "
                    f"{key!r} because {reason}.",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            self.results[key] = accumulator.result()


def tracking_snapshot(task: FrameTask, outcome: TaskOutcome) -> TrackFrameSnapshot | None:
    """Return the frame snapshot of a successful outcome, or ``None``.

    Workers attach the compact snapshot when the run requested it; a caller
    that still streams complete results is reduced here instead.
    """
    if not outcome.ok:
        return None
    if outcome.snapshot is not None:
        return outcome.snapshot
    if outcome.result is None:
        return None
    return snapshot_from_frame_result(outcome.result, int(task.frame_index))


def _tracking_config(config: Mapping[str, Any]) -> TrackingConfig:
    values = config.get("track", {})
    if not isinstance(values, Mapping):
        raise ValueError("track must be a mapping in sqq_config.yaml.")
    return TrackingConfig.from_mapping(
        {
            name: values[name]
            for name in _TRACKING_CONFIG_FIELDS
            if name in values
        }
    )


def prepare_tracking_sink(
    plan: RunPlan,
    config: Mapping[str, Any],
    group_configs: Mapping[int, Mapping[str, Any]],
    *,
    enabled: bool,
) -> tuple[AnalyzeTrackingSink | None, RunPlan]:
    if not enabled:
        return None, plan
    configs: dict[int | str, Mapping[str, Any]] = (
        dict(group_configs) if group_configs else {"run": config}
    )
    sink = AnalyzeTrackingSink(configs)
    context = replace(
        plan.context,
        retain_results=False,
        stream_results=True,
        tracking_snapshots=True,
    )
    return sink, replace(plan, context=context)


def write_analyze_tracking_outputs(
    plan: RunPlan,
    results: Mapping[int | str, TrackingResult],
) -> None:
    for key, result in results.items():
        root = (
            Path(plan.context.group_output_roots[key])
            if key != "run"
            else Path(plan.context.output_root)
        )
        track_root = root / "track"
        tables = build_track_tables(result)
        write_tracking_tables(result, track_root, tables=tables)
        write_tracking_result(result, track_root / "track_state.json")
        write_track_info(result, track_root, tables=tables)


__all__ = [
    "AnalyzeTrackingSink",
    "prepare_tracking_sink",
    "tracking_snapshot",
    "write_analyze_tracking_outputs",
]
