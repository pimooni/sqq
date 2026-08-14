from __future__ import annotations

"""Common Analyze/Track runner over the typed runtime contracts."""

from collections.abc import Callable
from dataclasses import dataclass

from .contracts import RunPlan, TaskOutcome
from .parallel.events import EventSink
from .parallel.executor import TaskHandler, execute


OutcomeSink = Callable[[TaskOutcome], None]


@dataclass(slots=True)
class AnalysisRunner:
    """Run one plan, then publish ordered outcomes to optional consumers."""

    outcome_sinks: tuple[OutcomeSink, ...] = ()

    def run(
        self,
        plan: RunPlan,
        *,
        task_handler: TaskHandler | None = None,
        event_sink: EventSink | None = None,
    ) -> list[TaskOutcome]:
        def publish(outcome: TaskOutcome) -> None:
            for sink in self.outcome_sinks:
                sink(outcome)

        outcomes = execute(
            plan,
            task_handler,
            plan.policy.backend,
            plan.policy.workers,
            event_sink,
            publish if plan.context.stream_results else None,
        )
        if not plan.context.stream_results:
            for outcome in outcomes:
                publish(outcome)
        return outcomes


__all__ = ["AnalysisRunner", "OutcomeSink"]
