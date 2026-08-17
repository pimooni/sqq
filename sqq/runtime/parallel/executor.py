"""Ordered serial/thread/process execution with deterministic failure cleanup."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace
from multiprocessing import get_context
from pathlib import Path
from typing import Any

from ...io.render import RenderSession
from ...io.reporting import failed_row
from ..contracts import FrameTask, RunPlan, TaskOutcome, TaskStatus
from ..frame_task import execute_frame_task
from .events import EventSink, StageEventKind, drain_stage_events, emit_event
from .policy import (
    limited_math_threads,
    normalize_parallel_backend,
    process_in_flight_limit,
    process_worker_cap,
)
from .worker import WorkerContext, execute_worker_task, initialize_worker_context


TaskHandler = Callable[[FrameTask], TaskOutcome]
OutcomeSink = Callable[[TaskOutcome], None]


def execute(
    plan: RunPlan,
    task_handler: TaskHandler | None = None,
    backend: str | None = None,
    workers: int | None = None,
    event_sink: EventSink | None = None,
    outcome_sink: OutcomeSink | None = None,
) -> list[TaskOutcome]:
    """Execute a plan and return outcomes in plan order.

    Progress events are advisory.  Results and exceptions come only from
    return values and Futures, so a saturated UI queue cannot lose science.
    """
    resolved_backend = normalize_parallel_backend(backend or plan.policy.backend)
    resolved_workers = max(1, int(workers or plan.policy.workers))
    if resolved_backend == "process" and (cap := process_worker_cap()) is not None:
        resolved_workers = min(resolved_workers, cap)
    tasks = tuple(plan.tasks)
    _validate_tasks(tasks)
    if not tasks:
        return []
    with limited_math_threads(plan.policy.math_threads):
        if resolved_backend == "serial" or resolved_workers == 1:
            return _execute_serial(
                plan,
                tasks,
                task_handler,
                event_sink,
                outcome_sink,
            )
        return _execute_concurrent(
            plan,
            tasks,
            task_handler,
            resolved_backend,
            resolved_workers,
            event_sink,
            outcome_sink,
        )


def _execute_serial(
    plan: RunPlan,
    tasks: tuple[FrameTask, ...],
    handler: TaskHandler | None,
    event_sink: EventSink | None,
    outcome_sink: OutcomeSink | None,
) -> list[TaskOutcome]:
    outcomes: list[TaskOutcome] = []
    local_context = WorkerContext(plan.context) if handler is None else None
    try:
        if local_context is not None:
            local_context.open()
        for task in tasks:
            _safe_emit(event_sink, StageEventKind.START, task.task_index, task.display_name)
            try:
                if handler is not None:
                    outcome = handler(task)
                else:
                    frame = local_context.load_frame(task) if local_context is not None else None
                    loaded_task = replace(task, frame=frame) if frame is not None else task
                    outcome = execute_frame_task(
                        loaded_task,
                        plan.context,
                        stage_callback=lambda stage, item=task: _safe_emit(
                            event_sink,
                            StageEventKind.STAGE,
                            item.task_index,
                            stage,
                        ),
                    )
                if plan.policy.strict and not outcome.ok:
                    raise RuntimeError(
                        outcome.error_message or f"Task {task.task_index} failed."
                    )
            except Exception as exc:
                if plan.policy.strict:
                    _cleanup_failed_plan(plan)
                    raise
                outcome = _failed_outcome(task, exc)
            if outcome_sink is not None:
                try:
                    outcome_sink(outcome)
                except Exception:
                    _cleanup_failed_plan(plan)
                    raise
            if plan.context.stream_results and not plan.context.retain_results:
                outcome = replace(outcome, result=None)
            outcomes.append(outcome)
            _safe_emit(
                event_sink,
                StageEventKind.COMPLETE,
                task.task_index,
                outcome.status.value,
            )
    finally:
        if local_context is not None:
            local_context.close()
    return outcomes


def _execute_concurrent(
    plan: RunPlan,
    tasks: tuple[FrameTask, ...],
    handler: TaskHandler | None,
    backend: str,
    workers: int,
    event_sink: EventSink | None,
    outcome_sink: OutcomeSink | None,
) -> list[TaskOutcome]:
    limit = max(
        workers,
        int(plan.policy.in_flight_limit or process_in_flight_limit(workers)),
    )
    task_iterator = iter(tasks)
    pending: dict[Future[TaskOutcome], FrameTask] = {}
    outcomes: dict[int, TaskOutcome] = {}
    waiting: dict[int, TaskOutcome] = {}
    task_order = tuple(task.task_index for task in tasks)
    next_publish = 0
    stage_queue: Any = None
    fatal: BaseException | None = None

    if backend == "process":
        mp_context = get_context("spawn")
        stage_queue = mp_context.Queue(maxsize=max(32, workers * 8))
        executor: ProcessPoolExecutor | ThreadPoolExecutor = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=mp_context,
            initializer=initialize_worker_context,
            initargs=(plan.context, stage_queue),
        )
        active_handler = handler or execute_worker_task
    else:
        executor = ThreadPoolExecutor(max_workers=workers)
        active_handler = handler or _LocalTaskHandler(plan)

    try:
        _submit_until_limit(executor, active_handler, task_iterator, pending, limit, event_sink)
        while pending:
            done, _ = wait(tuple(pending), timeout=0.05, return_when=FIRST_COMPLETED)
            if stage_queue is not None:
                drain_stage_events(stage_queue, event_sink)
            if not done:
                continue
            for future in sorted(done, key=lambda item: pending[item].task_index):
                task = pending.pop(future)
                try:
                    outcome = future.result()
                    if plan.policy.strict and not outcome.ok:
                        raise RuntimeError(
                            outcome.error_message or f"Task {task.task_index} failed."
                        )
                except BaseException as exc:
                    if plan.policy.strict or isinstance(exc, BrokenProcessPool):
                        fatal = exc
                        break
                    outcome = _failed_outcome(task, exc)
                if plan.context.stream_results:
                    waiting[task.task_index] = outcome
                else:
                    outcomes[task.task_index] = outcome
                _safe_emit(
                    event_sink,
                    StageEventKind.COMPLETE,
                    task.task_index,
                    outcome.status.value,
                )
            if fatal is not None:
                break
            if plan.context.stream_results:
                try:
                    while (
                        next_publish < len(task_order)
                        and task_order[next_publish] in waiting
                    ):
                        task_index = task_order[next_publish]
                        outcome = waiting.pop(task_index)
                        if outcome_sink is not None:
                            outcome_sink(outcome)
                        if not plan.context.retain_results:
                            outcome = replace(outcome, result=None)
                        outcomes[task_index] = outcome
                        next_publish += 1
                except BaseException as exc:
                    fatal = exc
                    break
            if fatal is not None:
                break
            _submit_until_limit(
                executor,
                active_handler,
                task_iterator,
                pending,
                max(0, limit - len(waiting)),
                event_sink,
            )

        if fatal is not None:
            raise fatal
    except BaseException:
        for future, task in pending.items():
            if future.cancel():
                _safe_emit(
                    event_sink,
                    StageEventKind.CANCELLED,
                    task.task_index,
                    "cancelled",
                )
        _drain_workers_while_exiting(pending, stage_queue, event_sink)
        _cleanup_failed_plan(plan)
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
        if stage_queue is not None:
            drain_stage_events(stage_queue, event_sink)
            try:
                stage_queue.close()
            except (AttributeError, EOFError, OSError, ValueError):
                pass
            try:
                stage_queue.join_thread()
            except (AttributeError, EOFError, OSError, ValueError):
                pass

    if plan.context.stream_results and waiting:
        raise RuntimeError("Concurrent streaming left unpublished frame results.")
    return [outcomes[task.task_index] for task in tasks]


def _submit_until_limit(
    executor: ProcessPoolExecutor | ThreadPoolExecutor,
    handler: TaskHandler,
    tasks: Iterable[FrameTask],
    pending: dict[Future[TaskOutcome], FrameTask],
    limit: int,
    event_sink: EventSink | None,
) -> None:
    iterator = iter(tasks)
    while len(pending) < limit:
        try:
            task = next(iterator)
        except StopIteration:
            return
        future = executor.submit(handler, task)
        pending[future] = task
        _safe_emit(event_sink, StageEventKind.START, task.task_index, task.display_name)


class _LocalTaskHandler:
    """Thread-safe callable binding a read-only run context."""

    def __init__(self, plan: RunPlan) -> None:
        self._context = plan.context

    def __call__(self, task: FrameTask) -> TaskOutcome:
        return execute_frame_task(task, self._context)


def _drain_workers_while_exiting(
    pending: dict[Future[TaskOutcome], FrameTask],
    stage_queue: Any,
    event_sink: EventSink | None,
) -> None:
    running = {future for future in pending if not future.cancelled() and not future.done()}
    while running:
        done, running = wait(running, timeout=0.05, return_when=FIRST_COMPLETED)
        if stage_queue is not None:
            drain_stage_events(stage_queue, event_sink)
        for future in done:
            try:
                future.result()
            except BaseException:
                pass


def _failed_outcome(task: FrameTask, exc: BaseException) -> TaskOutcome:
    source = str(task.source or "")
    return TaskOutcome(
        task_index=task.task_index,
        frame_index=task.frame_index,
        status=TaskStatus.FAILED,
        row=failed_row(task.display_name, source, str(exc)),
        error_type=type(exc).__name__,
        error_message=str(exc),
    )


def _cleanup_failed_plan(plan: RunPlan) -> None:
    candidates: set[Path] = set()
    if plan.context.fragment_dir is not None:
        candidates.add(plan.context.fragment_dir)
    candidates.update(plan.context.group_fragment_dirs.values())
    for path in candidates:
        RenderSession.cleanup_workspace(Path(path))


def _safe_emit(
    sink: EventSink | None,
    kind: StageEventKind,
    task_index: int,
    value: str,
) -> None:
    try:
        emit_event(sink, kind, task_index, value)
    except Exception:
        # A UI observer cannot invalidate scientific work.
        pass


def _validate_tasks(tasks: tuple[FrameTask, ...]) -> None:
    indexes = [task.task_index for task in tasks]
    if len(indexes) != len(set(indexes)):
        raise ValueError("RunPlan task_index values must be unique.")


__all__ = ["OutcomeSink", "TaskHandler", "execute"]
