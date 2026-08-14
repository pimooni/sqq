from __future__ import annotations

"""CPU, worker-count, and numeric-thread policy."""

import math
import os
import subprocess
import sys
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

from ...config import (
    DEFAULT_MODE,
    mode_worker_count,
    mode_worker_fraction,
    normalize_parallel_backend,
)


def effective_cpu_count() -> int:
    process_count = getattr(os, "process_cpu_count", None)
    if callable(process_count):
        value = process_count()
        if value:
            return max(1, int(value))
    affinity = getattr(os, "sched_getaffinity", None)
    if callable(affinity):
        try:
            return max(1, len(affinity(0)))
        except (OSError, NotImplementedError):
            pass
    return max(1, int(os.cpu_count() or 1))


@lru_cache(maxsize=1)
def physical_cpu_count() -> int:
    logical_available = effective_cpu_count()
    detected = _detect_physical_cpu_count()
    if detected is None:
        return logical_available
    return max(1, min(int(detected), logical_available))


def _detect_physical_cpu_count() -> int | None:
    try:
        import psutil  # type: ignore[import-not-found]

        count = psutil.cpu_count(logical=False)
        if count:
            return int(count)
    except Exception:
        pass
    if os.name == "nt":
        return _detect_command_int(
            [
                "powershell",
                "-NoProfile",
                "-Command",
                "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum",
            ]
        )
    if sys.platform == "darwin":
        return _detect_command_int(["sysctl", "-n", "hw.physicalcpu"])
    return _detect_linux_cpuinfo_physical_cores()


def _detect_command_int(command: list[str]) -> int | None:
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=2, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return None
    for token in completed.stdout.replace("=", " ").split():
        try:
            value = int(token)
        except ValueError:
            continue
        if value > 0:
            return value
    return None


def _detect_linux_cpuinfo_physical_cores() -> int | None:
    cpuinfo = Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return None
    try:
        blocks = cpuinfo.read_text(encoding="utf-8", errors="ignore").strip().split("\n\n")
    except OSError:
        return None
    seen_cores: set[tuple[str, str]] = set()
    physical_ids: set[str] = set()
    cores_per_socket: list[int] = []
    for block in blocks:
        fields: dict[str, str] = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
        physical_id = fields.get("physical id")
        core_id = fields.get("core id")
        if physical_id is not None and core_id is not None:
            seen_cores.add((physical_id, core_id))
        if physical_id is not None:
            physical_ids.add(physical_id)
        if "cpu cores" in fields:
            try:
                cores_per_socket.append(int(fields["cpu cores"]))
            except ValueError:
                pass
    if seen_cores:
        return len(seen_cores)
    if physical_ids and cores_per_socket:
        return len(physical_ids) * max(cores_per_socket)
    if cores_per_socket:
        return max(cores_per_socket)
    return None


def process_worker_cap() -> int | None:
    return 61 if os.name == "nt" else None


@contextmanager
def limited_math_threads(value: int) -> Iterator[None]:
    thread_count = max(1, int(value))
    names = (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    )
    previous = {name: os.environ.get(name) for name in names}
    try:
        for name in names:
            os.environ[name] = str(thread_count)
        yield
    finally:
        for name, old_value in previous.items():
            if old_value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = old_value


def process_in_flight_limit(workers: int) -> int:
    return max(1, int(workers)) * 3


def worker_policy_text(config: dict[str, Any]) -> str:
    value = config.get("parallel", {}).get("workers", "auto")
    reserve_text = "reserve 1 physical core"
    if is_auto_worker_request(value):
        mode = config.get("mode", DEFAULT_MODE)
        fixed_count = mode_worker_count(mode)
        if fixed_count is not None:
            unit = "worker" if fixed_count == 1 else "workers"
            return f"engine default ({fixed_count} {unit})"
        percent = int(round(mode_worker_fraction(mode) * 100))
        return f"auto ({percent}% of physical cores, {reserve_text})"
    try:
        request_text = describe_worker_request(value)
    except ValueError:
        request_text = str(value)
    return f"explicit ({request_text}, {reserve_text})"


def describe_worker_request(value: Any) -> str:
    kind, amount = classify_worker_request(value)
    if kind == "auto":
        return "auto"
    if kind == "fraction":
        return f"{float(amount) * 100.0:g}% of physical cores"
    if kind == "count":
        return f"{int(amount)} workers"
    raise _worker_value_error()


def is_auto_worker_request(value: Any) -> bool:
    return value is None or str(value).strip().lower() in {"", "auto"}


def classify_worker_request(value: Any) -> tuple[str, float | int | None]:
    text = str(value).strip().lower()
    if text in {"", "auto"}:
        return "auto", None
    if text.endswith("%"):
        percent = _parse_worker_number(text[:-1])
        if percent > 100:
            raise _worker_value_error()
        return "fraction", percent / 100.0
    if "." in text:
        fraction = _parse_worker_number(text)
        if fraction > 1:
            raise _worker_value_error()
        return "fraction", fraction
    try:
        count = int(text)
    except (TypeError, ValueError) as exc:
        raise _worker_value_error() from exc
    if count < 1:
        raise _worker_value_error()
    return "count", count


def _parse_worker_number(text: str) -> float:
    try:
        number = float(text)
    except (TypeError, ValueError) as exc:
        raise _worker_value_error() from exc
    if not math.isfinite(number) or number <= 0:
        raise _worker_value_error()
    return number


def _worker_value_error() -> ValueError:
    return ValueError(
        "parallel.workers / --worker must be 'auto', a positive integer worker "
        "count such as 1 or 4, a decimal CPU fraction in (0, 1] such as 0.5 or "
        "1.0, or a percentage in (0%, 100%]."
    )


def resolve_workers(
    value: Any,
    n_paths: int,
    mode: Any = DEFAULT_MODE,
    cpu_total: int | None = None,
    backend: str = "process",
) -> int:
    resolved_backend = normalize_parallel_backend(backend)
    single_or_serial = max(0, int(n_paths)) <= 1 or resolved_backend == "serial"
    if is_auto_worker_request(value):
        if single_or_serial:
            return 1
        physical_total = max(
            1, int(cpu_total if cpu_total is not None else physical_cpu_count())
        )
        fixed_count = mode_worker_count(mode)
        requested = (
            fixed_count
            if fixed_count is not None
            else max(1, int(physical_total * mode_worker_fraction(mode)))
        )
    else:
        classify_worker_request(value)
        if single_or_serial:
            return 1
        physical_total = max(
            1, int(cpu_total if cpu_total is not None else physical_cpu_count())
        )
        kind, amount = classify_worker_request(value)
        requested = (
            max(1, int(physical_total * float(amount)))
            if kind == "fraction"
            else int(amount)
        )
    requested = min(requested, max(1, physical_total - 1), max(1, n_paths))
    if resolved_backend == "process" and (cap := process_worker_cap()) is not None:
        requested = min(requested, cap)
    return max(1, requested)


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
]
