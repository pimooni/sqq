from __future__ import annotations

"""Process-worker helpers for independent coordinate-file analysis."""

import atexit
import os
import subprocess
import sys
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

from .io.lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    LammpsInputConfig,
    close_lammps_universe,
    frame_from_lammps_universe,
    lammps_atom_metadata,
    normalize_lammps_config,
    open_lammps_universe,
)
from .io.summary import failed_row
from .io.trajectory import (
    close_mdanalysis_universe,
    frame_from_mdanalysis_universe,
    open_mdanalysis_universe,
    trajectory_atom_metadata,
    read_frames,
)


StageEvent = tuple[str, int, str, float]

_WORKER_CONFIG: dict[str, Any] | None = None
_WORKER_OUTDIR: Path | None = None
_WORKER_STRICT = False
_WORKER_STAGE_QUEUE: Any = None
_WORKER_TRAJECTORY_PATH: Path | None = None
_WORKER_UNIVERSE: Any = None
_WORKER_TRAJECTORY_METADATA: tuple[tuple[int, int, str, str, int], ...] | None = None
_WORKER_LAMMPS_CONFIG: LammpsInputConfig | None = None


def initialize_file_worker(
    config: dict[str, Any],
    outdir: str,
    strict: bool,
    stage_queue: Any,
) -> None:
    """Install immutable run settings once in each spawned worker."""
    global _WORKER_CONFIG, _WORKER_OUTDIR, _WORKER_STRICT, _WORKER_STAGE_QUEUE
    _WORKER_CONFIG = config
    _WORKER_OUTDIR = Path(outdir)
    _WORKER_STRICT = bool(strict)
    _WORKER_STAGE_QUEUE = stage_queue


def initialize_trajectory_worker(
    config: dict[str, Any],
    outdir: str,
    strict: bool,
    stage_queue: Any,
    trajectory_path: str,
    topology_path: str,
    lammps_config: dict[str, Any] | None = None,
) -> None:
    """Open one private trajectory handle in every spawned worker."""
    initialize_file_worker(config, outdir, strict, stage_queue)
    global _WORKER_TRAJECTORY_PATH, _WORKER_UNIVERSE
    global _WORKER_TRAJECTORY_METADATA, _WORKER_LAMMPS_CONFIG
    _WORKER_TRAJECTORY_PATH = Path(trajectory_path)
    if _WORKER_TRAJECTORY_PATH.suffix.lower() in LAMMPS_TRAJECTORY_SUFFIXES:
        _WORKER_LAMMPS_CONFIG = normalize_lammps_config(lammps_config)
        _WORKER_UNIVERSE = open_lammps_universe(
            _WORKER_TRAJECTORY_PATH, Path(topology_path), _WORKER_LAMMPS_CONFIG
        )
        _WORKER_TRAJECTORY_METADATA = lammps_atom_metadata(
            _WORKER_UNIVERSE, _WORKER_LAMMPS_CONFIG
        )
    else:
        _WORKER_LAMMPS_CONFIG = None
        _WORKER_UNIVERSE = open_mdanalysis_universe(
            _WORKER_TRAJECTORY_PATH, Path(topology_path)
        )
        _WORKER_TRAJECTORY_METADATA = trajectory_atom_metadata(
            _WORKER_UNIVERSE
        )
    atexit.register(close_trajectory_worker)


def process_trajectory_frame_task(frame_index: int, raw_frame_index: int) -> tuple[int, dict[str, Any]]:
    """Seek, analyze, and write one trajectory frame with a worker-local Universe."""
    if _WORKER_CONFIG is None or _WORKER_OUTDIR is None or _WORKER_TRAJECTORY_PATH is None or _WORKER_UNIVERSE is None:
        raise RuntimeError("SQQ trajectory worker was not initialized.")
    from time import perf_counter
    from .pipeline import process_frame

    display_name = f"{_WORKER_TRAJECTORY_PATH.stem}_frame{raw_frame_index:06d}"
    _emit_stage("start", frame_index, display_name, perf_counter())

    def callback(stage: str) -> None:
        _emit_stage("stage", frame_index, stage, perf_counter())

    try:
        if _WORKER_LAMMPS_CONFIG is not None:
            frame = frame_from_lammps_universe(
                _WORKER_UNIVERSE,
                _WORKER_TRAJECTORY_PATH,
                raw_frame_index,
                _WORKER_LAMMPS_CONFIG,
                atom_metadata=_WORKER_TRAJECTORY_METADATA,
            )
        else:
            frame = frame_from_mdanalysis_universe(
                _WORKER_UNIVERSE,
                _WORKER_TRAJECTORY_PATH,
                raw_frame_index,
                atom_metadata=_WORKER_TRAJECTORY_METADATA,
            )
        row = process_frame(
            frame_index,
            frame,
            _WORKER_CONFIG,
            _WORKER_OUTDIR,
            strict=_WORKER_STRICT,
            stage_callback=callback,
        )
    except Exception as exc:
        if _WORKER_STRICT:
            raise
        row = failed_row(display_name, str(_WORKER_TRAJECTORY_PATH), str(exc))
    return frame_index, row


def process_trajectory_batch_task(
    items: tuple[tuple[int, int], ...],
) -> list[tuple[int, dict[str, Any]]]:
    """Analyze one small ordered trajectory batch with a worker-local reader."""
    from time import perf_counter

    results: list[tuple[int, dict[str, Any]]] = []
    for frame_index, raw_frame_index in items:
        result = process_trajectory_frame_task(frame_index, raw_frame_index)
        results.append(result)
        _emit_stage(
            "complete",
            frame_index,
            "ok" if result[1].get("status") == "ok" else "failed",
            perf_counter(),
        )
    return results


def close_trajectory_worker() -> None:
    """Close the private MDAnalysis reader before a worker exits."""
    global _WORKER_UNIVERSE, _WORKER_TRAJECTORY_METADATA, _WORKER_LAMMPS_CONFIG
    if _WORKER_UNIVERSE is not None:
        if _WORKER_LAMMPS_CONFIG is None:
            close_mdanalysis_universe(_WORKER_UNIVERSE)
        else:
            close_lammps_universe(_WORKER_UNIVERSE)
        _WORKER_UNIVERSE = None
    _WORKER_TRAJECTORY_METADATA = None
    _WORKER_LAMMPS_CONFIG = None


def process_file_task(frame_index: int, path_text: str) -> tuple[int, dict[str, Any]]:
    """Read, analyze, and write one independent GRO/XYZ file in a worker."""
    if _WORKER_CONFIG is None or _WORKER_OUTDIR is None:
        raise RuntimeError("SQQ process worker was not initialized.")

    from time import perf_counter

    from .pipeline import process_frame

    path = Path(path_text)
    _emit_stage("start", frame_index, path.name, perf_counter())

    def callback(stage: str) -> None:
        _emit_stage("stage", frame_index, stage, perf_counter())

    try:
        frame = next(
            iter(
                read_frames(
                    [path],
                    xyz_scale=float(
                        _WORKER_CONFIG.get("input", {}).get("xyz_scale", 0.1)
                    ),
                )
            )
        )
        row = process_frame(
            frame_index,
            frame,
            _WORKER_CONFIG,
            _WORKER_OUTDIR,
            strict=_WORKER_STRICT,
            stage_callback=callback,
        )
    except Exception as exc:
        if _WORKER_STRICT:
            raise
        row = failed_row(path.stem, str(path), str(exc))
    return frame_index, row


def _emit_stage(kind: str, frame_index: int, value: str, timestamp: float) -> None:
    """Send one small progress event without coupling workers to the terminal UI."""
    if _WORKER_STAGE_QUEUE is not None:
        _WORKER_STAGE_QUEUE.put((kind, frame_index, value, timestamp))


def effective_cpu_count() -> int:
    """Return CPUs available to this process, respecting scheduler affinity."""
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
    """Return physical CPU cores available to this process when detectable."""
    logical_available = effective_cpu_count()
    detected = _detect_physical_cpu_count()
    if detected is None:
        return logical_available
    return max(1, min(int(detected), logical_available))


def _detect_physical_cpu_count() -> int | None:
    """Best-effort physical-core detection without requiring optional packages."""
    try:
        import psutil  # type: ignore[import-not-found]

        count = psutil.cpu_count(logical=False)
        if count:
            return int(count)
    except Exception:
        pass

    if os.name == "nt":
        return _detect_command_int([
            "powershell",
            "-NoProfile",
            "-Command",
            "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum",
        ])
    if sys.platform == "darwin":
        return _detect_command_int(["sysctl", "-n", "hw.physicalcpu"])
    return _detect_linux_cpuinfo_physical_cores()


def _detect_command_int(command: list[str]) -> int | None:
    """Run one fixed platform command and parse a positive integer result."""
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=2, check=False)
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
    """Parse /proc/cpuinfo physical/core ids on Linux-like systems."""
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
    """Return the documented ProcessPoolExecutor worker cap on Windows."""
    return 61 if os.name == "nt" else None


@contextmanager
def limited_math_threads(value: int) -> Iterator[None]:
    """Give spawned workers one controlled BLAS/OpenMP thread each."""
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
