"""Typed contracts shared by Analyze, Track, and parallel workers."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from ..models import Frame, FrameResult


class InputKind(str, Enum):
    GRO = "gro"
    GRO_BATCH = "gro-batch"
    STACKED_GRO = "stacked-gro"
    TRAJECTORY = "trajectory"
    LAMMPS = "lammps"


class TaskStatus(str, Enum):
    OK = "ok"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    backend: str = "serial"
    workers: int = 1
    strict: bool = False
    in_flight_limit: int = 1
    math_threads: int = 1


@dataclass(frozen=True, slots=True)
class FrameTask:
    """One ordered unit of frame acquisition, analysis, and publication."""

    task_index: int
    frame_index: int
    source: Path | None = None
    raw_frame_index: int | None = None
    frame: Frame | None = None
    group_key: str | int | None = None
    output_name: str | None = None
    output_root: Path | None = None
    separated_output: bool = False

    @property
    def display_name(self) -> str:
        if self.output_name:
            return self.output_name
        if self.source is not None:
            if self.raw_frame_index is not None:
                return f"{self.source.stem}_frame{self.raw_frame_index:06d}"
            return self.source.name
        if self.frame is not None:
            return self.frame.name
        return f"frame_{self.frame_index:06d}"


@dataclass(frozen=True, slots=True)
class RunContext:
    """Resolved static inputs required to execute frame tasks."""

    config: Mapping[str, Any]
    output_root: Path
    strict: bool = False
    fragment_dir: Path | None = None
    topology: Path | None = None
    trajectory: Path | None = None
    input_kind: InputKind = InputKind.GRO
    retain_results: bool = False
    stream_results: bool = False
    group_configs: Mapping[str | int, Mapping[str, Any]] = field(default_factory=dict)
    group_output_roots: Mapping[str | int, Path] = field(default_factory=dict)
    group_fragment_dirs: Mapping[str | int, Path] = field(default_factory=dict)

    def config_for(self, task: FrameTask) -> Mapping[str, Any]:
        if task.group_key is None:
            return self.config
        try:
            return self.group_configs[task.group_key]
        except KeyError as exc:
            raise RuntimeError(
                f"No resolved configuration for topology group {task.group_key!r}."
            ) from exc

    def output_root_for(self, task: FrameTask) -> Path:
        if task.output_root is not None:
            return task.output_root
        if task.group_key is not None and task.group_key in self.group_output_roots:
            return self.group_output_roots[task.group_key]
        return self.output_root

    def fragment_dir_for(self, task: FrameTask) -> Path | None:
        if task.group_key is not None and task.group_key in self.group_fragment_dirs:
            return self.group_fragment_dirs[task.group_key]
        return self.fragment_dir


@dataclass(frozen=True, slots=True)
class TaskOutcome:
    """Authoritative result of one task, independent of progress events."""

    task_index: int
    frame_index: int
    status: TaskStatus
    row: Mapping[str, Any]
    result: FrameResult | None = None
    error_type: str | None = None
    error_message: str | None = None

    @property
    def ok(self) -> bool:
        return self.status is TaskStatus.OK


@dataclass(frozen=True, slots=True)
class RunPlan:
    input_kind: InputKind
    tasks: tuple[FrameTask, ...]
    context: RunContext
    policy: ExecutionPolicy
    topology_groups: Mapping[str | int, tuple[int, ...]] = field(default_factory=dict)
    sampling: Mapping[str, Any] = field(default_factory=dict)
    requested_graph_modes: Mapping[str | int, str] = field(default_factory=dict)
    effective_graph_modes: Mapping[str | int, str] = field(default_factory=dict)
    output_roots: tuple[Path, ...] = ()


__all__ = [
    "ExecutionPolicy",
    "FrameTask",
    "InputKind",
    "RunContext",
    "RunPlan",
    "TaskOutcome",
    "TaskStatus",
]
