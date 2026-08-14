from __future__ import annotations

"""Spawn-safe worker context with one private trajectory reader per process."""

import atexit
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from ...io.lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    LammpsInputConfig,
    close_lammps_universe,
    frame_from_lammps_universe,
    lammps_atom_metadata,
    normalize_lammps_config,
    open_lammps_universe,
)
from ...io.trajectory import (
    close_mdanalysis_universe,
    frame_from_mdanalysis_universe,
    open_mdanalysis_universe,
    trajectory_atom_metadata,
)
from ..contracts import FrameTask, RunContext, TaskOutcome
from ..frame_task import execute_frame_task
from .events import QueueStageEmitter


@dataclass(slots=True)
class WorkerContext:
    run_context: RunContext
    stage_queue: Any = None
    universe: Any = None
    trajectory_metadata: tuple[tuple[int, int, str, str, int], ...] | None = None
    lammps_config: LammpsInputConfig | None = None

    def open(self) -> None:
        trajectory = self.run_context.trajectory
        topology = self.run_context.topology
        if trajectory is None:
            return
        if topology is None:
            raise ValueError("Trajectory workers require a topology path.")
        trajectory = Path(trajectory)
        if trajectory.suffix.lower() in LAMMPS_TRAJECTORY_SUFFIXES:
            raw = self.run_context.config.get("input", {}).get("lammps")
            self.lammps_config = normalize_lammps_config(raw)
            self.universe = open_lammps_universe(
                trajectory,
                Path(topology),
                self.lammps_config,
            )
            self.trajectory_metadata = lammps_atom_metadata(
                self.universe,
                self.lammps_config,
            )
        else:
            self.universe = open_mdanalysis_universe(trajectory, Path(topology))
            self.trajectory_metadata = trajectory_atom_metadata(self.universe)

    def close(self) -> None:
        if self.universe is None:
            return
        if self.lammps_config is None:
            close_mdanalysis_universe(self.universe)
        else:
            close_lammps_universe(self.universe)
        self.universe = None
        self.trajectory_metadata = None
        self.lammps_config = None

    def load_frame(self, task: FrameTask):
        if task.frame is not None:
            return task.frame
        if self.universe is None:
            return None
        if task.raw_frame_index is None:
            raise ValueError(f"{task.display_name} requires a raw trajectory frame index.")
        trajectory = Path(self.run_context.trajectory or task.source or "")
        if self.lammps_config is not None:
            return frame_from_lammps_universe(
                self.universe,
                trajectory,
                task.raw_frame_index,
                self.lammps_config,
                atom_metadata=self.trajectory_metadata,
            )
        return frame_from_mdanalysis_universe(
            self.universe,
            trajectory,
            task.raw_frame_index,
            atom_metadata=self.trajectory_metadata,
        )


_PROCESS_CONTEXT: WorkerContext | None = None

def initialize_worker_context(run_context: RunContext, stage_queue: Any = None) -> None:
    """Top-level initializer suitable for the multiprocessing spawn method."""
    global _PROCESS_CONTEXT
    close_worker_context()
    if stage_queue is not None:
        try:
            # Workers only produce advisory progress.  Never make process
            # shutdown wait for a saturated queue feeder to flush.
            stage_queue.cancel_join_thread()
        except (AttributeError, OSError, ValueError):
            pass
    context = WorkerContext(run_context=run_context, stage_queue=stage_queue)
    context.open()
    _PROCESS_CONTEXT = context
    atexit.register(close_worker_context)


def current_worker_context() -> WorkerContext:
    if _PROCESS_CONTEXT is None:
        raise RuntimeError("SQQ worker context was not initialized.")
    return _PROCESS_CONTEXT


def execute_worker_task(task: FrameTask) -> TaskOutcome:
    context = current_worker_context()
    frame = context.load_frame(task)
    loaded_task = replace(task, frame=frame) if frame is not None else task
    emitter = (
        QueueStageEmitter(context.stage_queue, task.task_index)
        if context.stage_queue is not None
        else None
    )
    return execute_frame_task(
        loaded_task,
        context.run_context,
        stage_callback=emitter,
    )


def close_worker_context() -> None:
    global _PROCESS_CONTEXT
    if _PROCESS_CONTEXT is not None:
        _PROCESS_CONTEXT.close()
        _PROCESS_CONTEXT = None


__all__ = [
    "WorkerContext",
    "close_worker_context",
    "current_worker_context",
    "execute_worker_task",
    "initialize_worker_context",
]
