"""Small progress events that never carry scientific results."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from queue import Empty, Full
from time import perf_counter
from typing import Any, Protocol


class StageEventKind(str, Enum):
    START = "start"
    STAGE = "stage"
    COMPLETE = "complete"
    CANCELLED = "cancelled"


@dataclass(frozen=True, slots=True)
class StageEvent:
    kind: StageEventKind
    task_index: int
    value: str
    timestamp: float


class EventSink(Protocol):
    def __call__(self, event: StageEvent) -> None: ...


def emit_event(
    sink: EventSink | Callable[[StageEvent], None] | None,
    kind: StageEventKind,
    task_index: int,
    value: str,
) -> None:
    if sink is not None:
        sink(StageEvent(kind, int(task_index), str(value), perf_counter()))


class QueueStageEmitter:
    """Non-blocking worker emitter for redundant stage-only UI updates."""

    def __init__(self, queue: Any, task_index: int) -> None:
        self._queue = queue
        self._task_index = int(task_index)
        self.dropped = 0

    def __call__(self, value: str) -> None:
        event = StageEvent(
            StageEventKind.STAGE,
            self._task_index,
            str(value),
            perf_counter(),
        )
        try:
            self._queue.put_nowait(event)
        except (Full, AttributeError, EOFError, OSError, ValueError):
            try:
                self._queue.put(event, block=False)
            except (Full, AttributeError, EOFError, OSError, ValueError):
                self.dropped += 1


def drain_stage_events(
    queue: Any,
    sink: EventSink | Callable[[StageEvent], None] | None,
    *,
    limit: int | None = None,
) -> int:
    """Drain queued stage events without waiting for a producer."""
    drained = 0
    while limit is None or drained < limit:
        try:
            event = queue.get_nowait()
        except (Empty, AttributeError, EOFError, OSError, ValueError):
            break
        if sink is not None:
            try:
                sink(event)
            except Exception:
                # Progress observers are advisory and cannot fail a task.
                pass
        drained += 1
    return drained


__all__ = [
    "EventSink",
    "QueueStageEmitter",
    "StageEvent",
    "StageEventKind",
    "drain_stage_events",
    "emit_event",
]
