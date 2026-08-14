from __future__ import annotations

"""Workflow-level lifecycle facade over the shared runtime runner."""

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Protocol

from ..runtime.analysis_runner import AnalysisRunner as RuntimeAnalysisRunner
from ..runtime.contracts import FrameTask, RunPlan, TaskOutcome
from ..runtime.parallel.events import StageEvent, StageEventKind
from ..runtime.parallel.executor import TaskHandler


@dataclass(frozen=True, slots=True)
class AnalysisEvent:
    """One workflow event independent of terminal presentation."""

    kind: str
    task: FrameTask | None = None
    stage: str | None = None
    outcome: TaskOutcome | None = None
    status: str | None = None


class AnalysisSink(Protocol):
    """Consumer for reporting, rendering, or tracking state."""

    def start(self, plan: RunPlan) -> None: ...

    def consume(self, task: FrameTask, outcome: TaskOutcome) -> None: ...

    def finish(self, plan: RunPlan, outcomes: Sequence[TaskOutcome]) -> None: ...


class AnalysisRunner:
    """Publish workflow lifecycle around the single runtime executor."""

    def __init__(
        self,
        plan: RunPlan,
        *,
        event_sink: Callable[[AnalysisEvent], None] | None = None,
        sinks: Iterable[AnalysisSink] = (),
        task_handler: TaskHandler | None = None,
    ) -> None:
        self.plan = plan
        self.event_sink = event_sink
        self.sinks = tuple(sinks)
        self.task_handler = task_handler
        self._tasks = {task.task_index: task for task in plan.tasks}

    def run(self) -> tuple[TaskOutcome, ...]:
        """Execute the plan and return outcomes in plan order."""
        self._emit(AnalysisEvent("run-start"))
        for sink in self.sinks:
            sink.start(self.plan)

        def consume(outcome: TaskOutcome) -> None:
            task = self._tasks[outcome.task_index]
            for sink in self.sinks:
                sink.consume(task, outcome)

        try:
            outcomes = RuntimeAnalysisRunner((consume,)).run(
                self.plan,
                task_handler=self.task_handler,
                event_sink=self._runtime_event,
            )
            ordered = tuple(outcomes)
            for sink in self.sinks:
                sink.finish(self.plan, ordered)
            self._emit(AnalysisEvent("run-complete"))
            return ordered
        except Exception:
            self._emit(AnalysisEvent("run-failed"))
            raise

    def _runtime_event(self, event: StageEvent) -> None:
        task = self._tasks.get(event.task_index)
        if event.kind is StageEventKind.START:
            translated = AnalysisEvent("task-start", task=task)
        elif event.kind is StageEventKind.STAGE:
            translated = AnalysisEvent("stage", task=task, stage=event.value)
        elif event.kind is StageEventKind.COMPLETE:
            translated = AnalysisEvent(
                "task-complete",
                task=task,
                status=event.value,
            )
        else:
            translated = AnalysisEvent(
                "task-cancelled",
                task=task,
                status=event.value,
            )
        self._emit(translated)

    def _emit(self, event: AnalysisEvent) -> None:
        if self.event_sink is not None:
            self.event_sink(event)


__all__ = ["AnalysisEvent", "AnalysisRunner", "AnalysisSink"]
