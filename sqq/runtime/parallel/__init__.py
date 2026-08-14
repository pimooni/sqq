"""Lightweight parallel policy and event contracts."""

from .events import EventSink, StageEvent, StageEventKind
from .policy import (
    classify_worker_request,
    describe_worker_request,
    effective_cpu_count,
    is_auto_worker_request,
    limited_math_threads,
    normalize_parallel_backend,
    physical_cpu_count,
    process_in_flight_limit,
    process_worker_cap,
    resolve_workers,
    worker_policy_text,
)

__all__ = [
    "effective_cpu_count",
    "physical_cpu_count",
    "process_worker_cap",
    "limited_math_threads",
    "process_in_flight_limit",
    "classify_worker_request",
    "describe_worker_request",
    "is_auto_worker_request",
    "normalize_parallel_backend",
    "resolve_workers",
    "worker_policy_text",
    "EventSink",
    "StageEvent",
    "StageEventKind",
]
