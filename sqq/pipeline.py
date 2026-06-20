from __future__ import annotations

"""Top-level analysis pipeline for the SQQ command line."""

import os
import sys
from argparse import Namespace
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from threading import Event, Lock, Thread
from time import perf_counter
from typing import Any, Callable

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    class _NullProgress:
        def update(self, n: int = 1) -> None:
            pass

        def set_postfix_str(self, text: str, refresh: bool = True) -> None:
            pass

        def close(self) -> None:
            pass

    def tqdm(iterable=None, total=None, desc=None, **kwargs):
        return iterable if iterable is not None else _NullProgress()

from .banner import SQQ_BANNER
from .config import load_config, mode_label, mode_worker_fraction
from .core.cage import find_cages
from .core.f3f4 import compute_f3f4
from .core.graph import build_water_graph
from .core.ice import classify_ice_waters
from .core.quasi_cage import find_cage_patches
from .core.ring import find_rings
from .core.selection import select_guests, select_waters
from .io.gro_writer import (
    write_cage_gro_files,
    write_half_cage_gro_files,
    write_ice_gro_file,
    write_quasi_cage_gro_files,
    write_ring_gro_files,
)
from .io.summary import (
    dashboard_cage_targets,
    failed_row,
    result_row,
    write_f3f4,
    write_frame_info,
    write_membership,
    write_summary,
    write_vmd_script,
)
from .io.trajectory import expand_inputs, read_frames
from .models import Cage, CagePatch, Frame, FrameResult


PARALLEL_SUFFIXES = {".gro", ".xyz"}


def analyze(args: Namespace) -> None:
    """Run SQQ analysis from parsed command-line arguments."""
    run_started_at = datetime.now().astimezone()
    started_at = perf_counter()
    config = load_config(Path(args.config) if args.config else None, mode=getattr(args, "mode", None))
    apply_cli_overrides(config, args)

    # Directory input follows a one-file-per-frame workflow.
    input_path = Path(args.input)
    pattern = args.pattern or config["input"]["pattern"]
    recursive = bool(args.recursive or config["input"]["recursive"])
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    topology = Path(args.topology) if args.topology else None
    parallelizable = can_parallelize_paths(paths, topology)
    workers = (
        resolve_workers(config["parallel"].get("workers"), len(paths), mode=config.get("mode", "50"))
        if parallelizable
        else 1
    )
    print_run_header(args, config, input_path, outdir, paths, topology, workers, run_started_at)
    if workers > 1 and parallelizable:
        rows = analyze_paths_parallel(
            paths,
            outdir,
            config,
            workers=workers,
            strict=bool(args.strict),
            total_started_at=started_at,
        )
    else:
        rows = analyze_paths_serial(
            paths,
            outdir,
            config,
            topology=topology,
            strict=bool(args.strict),
            total_started_at=started_at,
        )

    elapsed_seconds = perf_counter() - started_at
    run_finished_at = datetime.now().astimezone()
    run_info = build_run_info(
        args,
        config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        elapsed_seconds,
        run_started_at,
        run_finished_at,
    )
    write_summary(rows, outdir, config, write_xlsx=config["output"]["write_xlsx_summary"], run_info=run_info)
    print(f"Wrote SQQ results: {outdir}")


def build_run_info(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    outdir: Path,
    paths: list[Path],
    topology: Path | None,
    workers: int,
    elapsed_seconds: float,
    started_at_wall: datetime,
    finished_at_wall: datetime,
) -> dict[str, Any]:
    """Collect run-level metadata for the summary workbook."""
    info: dict[str, Any] = {
        "working_dir": str(Path.cwd()),
        "input": str(input_path),
        "output_dir": str(outdir.resolve()),
        "date": started_at_wall.strftime("%Y-%m-%d"),
        "start_time": started_at_wall.strftime("%H:%M:%S"),
        "finish_time": finished_at_wall.strftime("%H:%M:%S"),
        "started_at": started_at_wall.isoformat(timespec="seconds"),
        "finished_at": finished_at_wall.isoformat(timespec="seconds"),
        "time_zone": format_time_zone(started_at_wall),
        "config_file": args.config or "<built-in defaults>",
        "mode": f"{config.get('mode', '50')} ({mode_label(config.get('mode', '50'))})",
        "worker_policy": worker_policy_text(config),
        "topology": str(topology) if topology else "<none>",
        "matched_files": len(paths),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "graph_mode": config["graph"]["bond_mode"],
        "ring_sizes": config["ring"]["sizes"],
        "quasi_cage_base_sizes": config["quasi_cage"].get("base_sizes", "auto"),
        "quasi_cage_side_sizes": config["quasi_cage"].get("side_sizes", "auto"),
        "cage_sizes": config["cage"].get("ring_sizes", [5, 6]),
        "output_layout": config["output"].get("structure_layout", "grouped"),
        "workers": workers,
        "summary_xlsx": str((outdir / "summary.xlsx").resolve()),
        "run_config": str((outdir / "run_config.yaml").resolve()),
    }
    if paths:
        info["first_file"] = str(paths[0].resolve())
        info["last_file"] = str(paths[-1].resolve())
    return info


def print_run_header(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    outdir: Path,
    paths: list[Path],
    topology: Path | None,
    workers: int,
    started_at_wall: datetime,
) -> None:
    """Print static run information before the live progress display starts."""
    print(SQQ_BANNER)
    print("Basic Information")
    print_terminal_field("date", started_at_wall.strftime("%Y-%m-%d"))
    print_terminal_field("start_time", started_at_wall.strftime("%H:%M:%S"))
    print_terminal_field("time_zone", format_time_zone(started_at_wall))
    print_terminal_field("working_dir", Path.cwd())
    print_terminal_field("input", input_path)
    print_terminal_field("matched_files", len(paths))
    print_terminal_field("output", outdir)
    print("")
    print("Configuration")
    print_terminal_field("config", args.config or "<built-in defaults>")
    print_terminal_field("topology", topology or "<none>")
    print_terminal_field("mode", f"{config.get('mode', '50')} ({mode_label(config.get('mode', '50'))})")
    print_terminal_field("graph_mode", config["graph"]["bond_mode"])
    print_terminal_field("ring_sizes", config["ring"]["sizes"])
    print_terminal_field("quasi_cage_sizes", f"{config['quasi_cage'].get('base_sizes', 'auto')} / {config['quasi_cage'].get('side_sizes', 'auto')}")
    print_terminal_field("quasi_max_layers", config["quasi_cage"].get("max_layers", ""))
    print_terminal_field("cage_sizes", config["cage"].get("ring_sizes", [5, 6]))
    print_terminal_field("cage_targets", dashboard_cage_targets(config))
    print_terminal_field("other_cages", config["cage"].get("output_other", False))
    print_terminal_field("output_layout", config["output"].get("structure_layout", "grouped"))
    print_terminal_field("worker_policy", worker_policy_text(config))
    print_terminal_field("workers", workers)
    print("")


TERMINAL_LABEL_WIDTH = 22
PROGRESS_BAR_WIDTH = 25


def print_terminal_field(label: str, value: Any) -> None:
    """Print one aligned terminal key-value row."""
    print(f"  {label:<{TERMINAL_LABEL_WIDTH}}: {safe_terminal_text(value)}")


def format_terminal_value(value: Any) -> str:
    """Format terminal values compactly without Python container brackets."""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (list, tuple, set)):
        return ",".join(str(item) for item in value)
    return str(value)


def safe_terminal_text(value: Any) -> str:
    """Avoid UnicodeEncodeError on legacy Windows consoles."""
    text = ascii_superscript_text(format_terminal_value(value))
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    return text.encode(encoding, errors="replace").decode(encoding, errors="replace")


SUPERSCRIPT_DIGITS = {
    "⁰": "0",
    "¹": "1",
    "²": "2",
    "³": "3",
    "⁴": "4",
    "⁵": "5",
    "⁶": "6",
    "⁷": "7",
    "⁸": "8",
    "⁹": "9",
    "⁻": "-",
}


def ascii_superscript_text(text: str) -> str:
    """Convert Unicode superscripts to compact ASCII exponents for terminal output."""
    result: list[str] = []
    in_superscript = False
    for char in text:
        if char in SUPERSCRIPT_DIGITS:
            if not in_superscript:
                result.append("^")
            result.append(SUPERSCRIPT_DIGITS[char])
            in_superscript = True
            continue
        if in_superscript and char.isdigit():
            result.append(" ")
        in_superscript = False
        result.append(char)
    return "".join(result)


TIME_ZONE_ALIASES = {
    ("CST", 480): "China Standard Time",
    ("\u4e2d\u56fd\u6807\u51c6\u65f6\u95f4", 480): "China Standard Time",
}


def format_time_zone(value: datetime) -> str:
    """Format a time-zone name with its signed UTC offset."""
    name = value.tzname() or "UTC"
    offset = value.utcoffset()
    if offset is None:
        return name
    total_minutes = int(offset.total_seconds() / 60)
    name = TIME_ZONE_ALIASES.get((name, total_minutes), name)
    sign = "+" if total_minutes >= 0 else "-"
    hours, minutes = divmod(abs(total_minutes), 60)
    offset_text = f"{sign}{hours}" if minutes == 0 else f"{sign}{hours}:{minutes:02d}"
    return f"{name} ({offset_text})"


def format_seconds(seconds: float) -> str:
    """Format elapsed seconds for the live terminal display."""
    return f"{max(seconds, 0.0):.1f} s"


class RunProgressDisplay:
    """Render per-run progress with current stage, frame, and total timings."""

    def __init__(self, total: int, total_started_at: float) -> None:
        self.total = total
        self.total_started_at = total_started_at
        self.completed = 0
        self.failed = 0
        self.current_index: int | None = None
        self.current_file = "waiting"
        self.stage = "waiting"
        self.frame_started_at = perf_counter()
        self.stage_started_at = perf_counter()
        self._lock = Lock()
        self._stop_event = Event()
        self._rendered_lines = 0
        self._interactive = bool(getattr(sys.stdout, "isatty", lambda: False)())
        self._progress = None if self._interactive else tqdm(total=total, desc="Files", unit="file")
        self._thread: Thread | None = None
        self._render()
        if self._interactive:
            self._thread = Thread(target=self._tick, daemon=True)
            self._thread.start()

    def start_frame(self, frame_index: int, frame_name: str) -> Callable[[str], None]:
        """Select the active frame and return its stage callback."""
        with self._lock:
            now = perf_counter()
            self.current_index = frame_index
            self.current_file = frame_name
            self.stage = "reading frame"
            self.frame_started_at = now
            self.stage_started_at = now
            self._render_locked()
        return self.update_stage

    def update_stage(self, stage: str) -> None:
        """Update the active stage and reset the stage timer when it changes."""
        with self._lock:
            if stage != self.stage:
                self.stage = stage
                self.stage_started_at = perf_counter()
            self._render_locked()

    def complete_frame(self, success: bool) -> None:
        """Mark the current frame as finished and refresh the display."""
        with self._lock:
            self.completed += 1
            if not success:
                self.failed += 1
            if self._progress is not None:
                self._progress.update(1)
            self._render_locked()

    def close(self) -> None:
        """Stop background refresh and close any progress backend."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=0.2)
        with self._lock:
            self._render_locked()
        if self._progress is not None:
            self._progress.close()

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render()

    def _render(self) -> None:
        with self._lock:
            self._render_locked()

    def _render_locked(self) -> None:
        if self._interactive:
            lines = self._panel_lines()
            if self._rendered_lines:
                sys.stdout.write(f"\033[{self._rendered_lines}F")
            for line in lines:
                sys.stdout.write("\r\033[K" + line + "\n")
            sys.stdout.flush()
            self._rendered_lines = len(lines)
            return
        if self._progress is not None and hasattr(self._progress, "set_postfix_str"):
            self._progress.set_postfix_str(self._postfix_text(), refresh=True)

    def _panel_lines(self) -> list[str]:
        return [
            "Analysis Progress",
            f"  {'completed_files':<{TERMINAL_LABEL_WIDTH}}: {self.completed} / {self.total}  [ {self.failed} failed ]",
            f"  {'current_file':<{TERMINAL_LABEL_WIDTH}}: {self._current_file_text()}",
            f"  {'stage':<{TERMINAL_LABEL_WIDTH}}: {self.stage}",
            f"  {'stage / frame / total':<{TERMINAL_LABEL_WIDTH}}: {self._time_text()}",
            "",
            self._files_bar(),
        ]

    def _postfix_text(self) -> str:
        return (
            f"completed_files: {self.completed} / {self.total} [ {self.failed} failed ]; "
            f"current_file: {self._current_file_text()}; "
            f"stage: {self.stage}; "
            f"stage / frame / total: {self._time_text()}"
        )

    def _current_file_text(self) -> str:
        if self.current_index is None:
            return "waiting"
        if self.total <= 1:
            return self.current_file
        return f"{self.current_index + 1} / {self.total}  {self.current_file}"

    def _time_text(self) -> str:
        now = perf_counter()
        stage_elapsed = now - self.stage_started_at
        frame_elapsed = now - self.frame_started_at if self.current_index is not None else 0.0
        total_elapsed = now - self.total_started_at
        return f"{format_seconds(stage_elapsed)} / {format_seconds(frame_elapsed)} / {format_seconds(total_elapsed)}"

    def _files_bar(self) -> str:
        if self.total <= 0:
            fraction = 1.0
        else:
            fraction = min(max(self.completed / self.total, 0.0), 1.0)
        filled = int(round(PROGRESS_BAR_WIDTH * fraction))
        bar = "█" * filled + " " * (PROGRESS_BAR_WIDTH - filled)
        return f"Files: {fraction * 100:3.0f}%|{bar}| {self.completed}/{self.total} completed"


PARALLEL_STAGE_GROUPS = (
    (
        ("reading frame", "reading"),
        ("resolving settings", "settings"),
        ("selecting molecules", "selecting"),
    ),
    (
        ("building water graph", "graph"),
        ("searching rings", "ring"),
        ("searching half/quasi cage", "half/quasi"),
        ("searching cage", "cage"),
    ),
    (
        ("filtering free patches", "filtering"),
        ("computing F3/F4", "F3/F4"),
        ("classifying ice", "ice"),
        ("writing outputs", "output"),
    ),
)
PARALLEL_FILE_PREVIEW_LIMIT = 6
PARALLEL_STAGE_COLUMN_WIDTH = 18
PARALLEL_FILE_COLUMN_WIDTH = 25
PARALLEL_ACTIVE_STAGE_WIDTH = 30


class ParallelRunProgressDisplay:
    """Render aggregate and per-file progress for concurrent frame analysis."""

    def __init__(self, total: int, workers: int, total_started_at: float) -> None:
        self.total = total
        self.workers = workers
        self.total_started_at = total_started_at
        self.completed = 0
        self.failed = 0
        self._active: dict[int, dict[str, Any]] = {}
        self._lock = Lock()
        self._stop_event = Event()
        self._rendered_lines = 0
        self._interactive = bool(getattr(sys.stdout, "isatty", lambda: False)())
        self._progress = None if self._interactive else tqdm(total=total, desc="Files", unit="file")
        self._thread: Thread | None = None
        self._render()
        if self._interactive:
            self._thread = Thread(target=self._tick, daemon=True)
            self._thread.start()

    def start_file(self, frame_index: int, frame_name: str) -> Callable[[str], None]:
        """Register an active file and return its stage callback."""
        with self._lock:
            now = perf_counter()
            self._active[frame_index] = {
                "name": frame_name,
                "stage": "reading frame",
                "file_started_at": now,
                "stage_started_at": now,
            }
            self._render_locked()
        return lambda stage: self.update_stage(frame_index, stage)

    def update_stage(self, frame_index: int, stage: str) -> None:
        """Update one active file without disturbing other worker states."""
        if stage == "done":
            return
        with self._lock:
            state = self._active.get(frame_index)
            if state is None:
                return
            if stage != state["stage"]:
                state["stage"] = stage
                state["stage_started_at"] = perf_counter()
            self._render_locked()

    def complete_file(self, frame_index: int, success: bool) -> None:
        """Move one file from the active set into completed results."""
        with self._lock:
            self._active.pop(frame_index, None)
            self.completed += 1
            if not success:
                self.failed += 1
            if self._progress is not None:
                self._progress.update(1)
            self._render_locked()

    def close(self) -> None:
        """Stop background refresh and close the progress backend."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=0.2)
        with self._lock:
            self._render_locked()
        if self._progress is not None:
            self._progress.close()

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render()

    def _render(self) -> None:
        with self._lock:
            self._render_locked()

    def _render_locked(self) -> None:
        if self._interactive:
            lines = self._panel_lines()
            if self._rendered_lines:
                sys.stdout.write(f"\033[{self._rendered_lines}F")
            for line in lines:
                sys.stdout.write("\r\033[K" + line + "\n")
            sys.stdout.flush()
            self._rendered_lines = len(lines)
            return
        if self._progress is not None and hasattr(self._progress, "set_postfix_str"):
            self._progress.set_postfix_str(self._postfix_text(), refresh=True)

    def _panel_lines(self) -> list[str]:
        stage_lines = self._stage_summary_lines()
        indent = " " * (TERMINAL_LABEL_WIDTH + 4)
        lines = [
            "Analysis Progress",
            f"  {'completed_files':<{TERMINAL_LABEL_WIDTH}}: {self.completed} / {self.total}  [ {self.failed} failed ]",
            f"  {'active_workers':<{TERMINAL_LABEL_WIDTH}}: {len(self._active)} / {self.workers}",
            f"  {'queued_files':<{TERMINAL_LABEL_WIDTH}}: {self._queued_files()}",
            f"  {'stage_summary':<{TERMINAL_LABEL_WIDTH}}: {stage_lines[0]}",
            indent + stage_lines[1],
            indent + stage_lines[2],
            f"  {'total_elapsed':<{TERMINAL_LABEL_WIDTH}}: {format_seconds(perf_counter() - self.total_started_at)}",
            "",
            "  active files",
            f"    {'file':<{PARALLEL_FILE_COLUMN_WIDTH}} {'stage':<{PARALLEL_ACTIVE_STAGE_WIDTH}} stage / file",
        ]
        active_items = sorted(self._active.items())
        preview_slots = min(PARALLEL_FILE_PREVIEW_LIMIT, self.workers)
        now = perf_counter()
        for slot in range(preview_slots):
            if slot < len(active_items):
                lines.append(self._active_file_line(active_items[slot], now))
            else:
                lines.append("")
        if self.workers > preview_slots:
            overflow = max(0, len(active_items) - preview_slots)
            lines.append(f"    ... {overflow} additional active files" if overflow else "")
        lines.extend(["", self._files_bar()])
        return lines

    def _postfix_text(self) -> str:
        stage_text = " | ".join(
            f"{label} {count}"
            for line in self._stage_summary_values()
            for label, count in line
        )
        return (
            f"completed_files: {self.completed} / {self.total} [ {self.failed} failed ]; "
            f"active_workers: {len(self._active)} / {self.workers}; "
            f"queued_files: {self._queued_files()}; stages: {stage_text}; "
            f"total_elapsed: {format_seconds(perf_counter() - self.total_started_at)}"
        )

    def _stage_counts(self) -> Counter[str]:
        return Counter(str(state["stage"]) for state in self._active.values())

    def _stage_summary_values(self) -> list[list[tuple[str, int]]]:
        counts = self._stage_counts()
        return [[(label, counts.get(stage, 0)) for stage, label in group] for group in PARALLEL_STAGE_GROUPS]

    def _stage_summary_lines(self) -> list[str]:
        rows = []
        for values in self._stage_summary_values():
            cells = [f"{label} {count}" for label, count in values]
            rows.append(" | ".join(f"{cell:<{PARALLEL_STAGE_COLUMN_WIDTH}}" for cell in cells).rstrip())
        return rows

    def _active_file_line(self, item: tuple[int, dict[str, Any]], now: float) -> str:
        frame_index, state = item
        file_text = compact_terminal_text(f"{frame_index + 1}/{self.total}  {state['name']}", PARALLEL_FILE_COLUMN_WIDTH)
        stage_text = compact_terminal_text(str(state["stage"]), PARALLEL_ACTIVE_STAGE_WIDTH)
        stage_elapsed = format_seconds(now - float(state["stage_started_at"]))
        file_elapsed = format_seconds(now - float(state["file_started_at"]))
        return (
            f"    {file_text:<{PARALLEL_FILE_COLUMN_WIDTH}} "
            f"{stage_text:<{PARALLEL_ACTIVE_STAGE_WIDTH}} "
            f"{stage_elapsed:>8} / {file_elapsed:>8}"
        )

    def _queued_files(self) -> int:
        return max(0, self.total - self.completed - len(self._active))

    def _files_bar(self) -> str:
        fraction = 1.0 if self.total <= 0 else min(max(self.completed / self.total, 0.0), 1.0)
        filled = int(round(PROGRESS_BAR_WIDTH * fraction))
        bar = "█" * filled + " " * (PROGRESS_BAR_WIDTH - filled)
        return f"Files: {fraction * 100:3.0f}%|{bar}| {self.completed}/{self.total} completed"


def compact_terminal_text(text: str, width: int) -> str:
    """Truncate a live-panel field without changing its column width."""
    if len(text) <= width:
        return text
    if width <= 3:
        return text[:width]
    return text[: width - 3] + "..."


def analyze_paths_serial(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    topology: Path | None,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Analyze frames in input order."""
    rows: list[dict[str, Any]] = []
    frames = read_frames(paths, topology=topology, xtc_stride=int(config["input"].get("xtc_stride", 1)))
    progress = RunProgressDisplay(total=len(paths), total_started_at=total_started_at)
    try:
        for frame_index, frame in enumerate(frames):
            callback = progress.start_frame(frame_index, frame.name)
            row = process_frame(frame_index, frame, config, outdir, strict=strict, stage_callback=callback)
            rows.append(row)
            progress.complete_frame(row.get("status") == "ok")
    finally:
        progress.close()
    return rows


def analyze_paths_parallel(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    workers: int,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Analyze independent coordinate files with live per-worker stages."""
    rows_by_index: dict[int, dict[str, Any]] = {}
    progress = ParallelRunProgressDisplay(total=len(paths), workers=workers, total_started_at=total_started_at)
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(process_single_file_path, frame_index, path, config, outdir, strict, progress): frame_index
                for frame_index, path in enumerate(paths)
            }
            for future in as_completed(futures):
                expected_index = futures[future]
                try:
                    frame_index, row = future.result()
                except Exception:
                    progress.complete_file(expected_index, False)
                    raise
                rows_by_index[frame_index] = row
                progress.complete_file(frame_index, row.get("status") == "ok")
    finally:
        progress.close()
    return [rows_by_index[index] for index in sorted(rows_by_index)]


def process_single_file_path(
    frame_index: int,
    path: Path,
    config: dict[str, Any],
    outdir: Path,
    strict: bool,
    progress: ParallelRunProgressDisplay | None = None,
) -> tuple[int, dict[str, Any]]:
    """Read and analyze one standalone coordinate file."""
    callback = progress.start_file(frame_index, path.name) if progress is not None else None
    frame = next(iter(read_frames([path])))
    return frame_index, process_frame(
        frame_index,
        frame,
        config,
        outdir,
        strict=strict,
        stage_callback=callback,
    )


def process_frame(
    frame_index: int,
    frame: Frame,
    config: dict[str, Any],
    outdir: Path,
    strict: bool,
    stage_callback: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Analyze one frame, write per-frame files, and return a summary row."""
    if frame.time_ps is None:
        frame.time_ps = config["input"]["first_file_time_ps"] + frame_index * config["input"]["frame_time_step_ps"]
    try:
        result = analyze_frame(frame, config, stage_callback=stage_callback)
        frame_dir = outdir / frame.name
        frame_dir.mkdir(parents=True, exist_ok=True)
        report_stage(stage_callback, "writing outputs")
        write_frame_outputs(result, frame_dir, config)
        report_stage(stage_callback, "done")
        return result_row(result)
    except Exception as exc:
        if strict:
            raise
        return failed_row(frame.name, str(frame.source or ""), str(exc))


def report_stage(callback: Callable[[str], None] | None, stage: str) -> None:
    """Update the terminal stage display when a callback is available."""
    if callback is not None:
        callback(stage)


def write_frame_outputs(result: FrameResult, frame_dir: Path, config: dict[str, Any]) -> None:
    """Write all configured per-frame output files."""
    if config["output"].get("write_info", True):
        write_frame_info(
            result,
            frame_dir,
            ring_sizes=resolve_size_list(
                config["ring"].get("sizes", []),
                fallback=sorted(result.rings),
                key="ring.sizes",
            ),
        )
    else:
        remove_optional_info_output(result, frame_dir)
    if config["output"]["write_tsv"]:
        write_membership(result, frame_dir)
        write_f3f4(result, frame_dir)
    else:
        remove_optional_tsv_outputs(result, frame_dir)
    if config["output"]["write_vmd"]:
        write_vmd_script(result, frame_dir)
    else:
        remove_optional_vmd_output(result, frame_dir)
    if config["output"]["write_gro"]:
        layout = str(config["output"].get("structure_layout", "grouped"))
        if config["output"].get("write_ring_gro", True):
            write_ring_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"], layout=layout)
        if config["output"].get("write_half_cage_gro", True):
            write_half_cage_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"], layout=layout)
        if config["output"].get("write_quasi_cage_gro", True):
            write_quasi_cage_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"], layout=layout)
        if config["output"].get("write_cage_gro", True):
            write_cage_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"], layout=layout)
        if config["output"].get("write_ice_gro", True):
            write_ice_gro_file(result, frame_dir, write_empty=config["output"]["write_empty_files"], layout=layout)


def remove_optional_info_output(result: FrameResult, frame_dir: Path) -> None:
    """Remove stale info output when info output is disabled."""
    (frame_dir / f"{result.frame.name}_info.md").unlink(missing_ok=True)


def remove_optional_tsv_outputs(result: FrameResult, frame_dir: Path) -> None:
    """Remove stale optional TSV files when TSV output is disabled."""
    for suffix in ("membership", "f3f4"):
        path = frame_dir / f"{result.frame.name}_{suffix}.tsv"
        path.unlink(missing_ok=True)


def remove_optional_vmd_output(result: FrameResult, frame_dir: Path) -> None:
    """Remove stale optional VMD helper scripts when disabled."""
    (frame_dir / f"{result.frame.name}_view.vmd.tcl").unlink(missing_ok=True)


def worker_policy_text(config: dict[str, Any]) -> str:
    """Describe the active automatic or explicit worker policy."""
    value = config.get("parallel", {}).get("workers", "auto")
    if value in (None, "", "auto"):
        percent = int(round(mode_worker_fraction(config.get("mode", "50")) * 100))
        return f"auto ({percent}% of logical CPUs)"
    return f"explicit ({value})"


def resolve_workers(value: Any, n_paths: int, mode: Any = "50", cpu_total: int | None = None) -> int:
    """Resolve the frame-level worker count and cap it by independent files."""
    if value in (None, "", "auto"):
        logical_cpus = max(1, int(cpu_total if cpu_total is not None else (os.cpu_count() or 1)))
        requested = max(1, int(logical_cpus * mode_worker_fraction(mode)))
    else:
        try:
            requested = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("parallel.workers / --workers must be 'auto' or a positive integer.") from exc
        if requested < 1:
            raise ValueError("parallel.workers / --workers must be 'auto' or a positive integer.")
    return min(requested, max(1, n_paths))


def can_parallelize_paths(paths: list[Path], topology: Path | None) -> bool:
    """Only standalone coordinate files are parallelized in v0.1."""
    if topology is not None:
        return False
    return bool(paths) and all(path.suffix.lower() in PARALLEL_SUFFIXES for path in paths)


def apply_cli_overrides(config: dict[str, Any], args: Namespace) -> None:
    """Apply command-line options after YAML/default configuration."""
    cage_size_overridden = False
    if args.pattern:
        config["input"]["pattern"] = args.pattern
    if args.recursive:
        config["input"]["recursive"] = True
    if getattr(args, "sizes", None):
        config["ring"]["sizes"] = args.sizes
        config["quasi_cage"]["base_sizes"] = args.sizes
        config["quasi_cage"]["side_sizes"] = args.sizes
        config["cage"]["ring_sizes"] = args.sizes
        cage_size_overridden = True
    if getattr(args, "ring_sizes", None):
        config["ring"]["sizes"] = args.ring_sizes
    if getattr(args, "quasi_sizes", None):
        config["quasi_cage"]["base_sizes"] = args.quasi_sizes
        config["quasi_cage"]["side_sizes"] = args.quasi_sizes
    if getattr(args, "quasi_base_sizes", None):
        config["quasi_cage"]["base_sizes"] = args.quasi_base_sizes
    if getattr(args, "quasi_side_sizes", None):
        config["quasi_cage"]["side_sizes"] = args.quasi_side_sizes
    if getattr(args, "quasi_max_layers", None) is not None:
        if args.quasi_max_layers < 1:
            raise ValueError("--quasi-max-layers must be at least 1.")
        config["quasi_cage"]["max_layers"] = args.quasi_max_layers
    if getattr(args, "cage_sizes", None):
        config["cage"]["ring_sizes"] = args.cage_sizes
        cage_size_overridden = True
    if cage_size_overridden and cage_sizes_need_other_outputs(config["cage"]["ring_sizes"]):
        config["cage"]["output_other"] = True
    if getattr(args, "other_cages", False):
        config["cage"]["output_other"] = True
    if getattr(args, "no_other_cages", False):
        config["cage"]["output_other"] = False
    if getattr(args, "other_max_faces", None) is not None:
        config["cage"]["other_max_faces"] = args.other_max_faces
    bond_mode = getattr(args, "bond_mode", None)
    if args.pairs:
        if bond_mode not in (None, "pairs"):
            raise ValueError("--pairs can only be combined with --bond-mode pairs.")
        config["graph"]["pair_file"] = args.pairs
        config["graph"]["bond_mode"] = "pairs"
    elif bond_mode is not None:
        config["graph"]["bond_mode"] = bond_mode
    if config["graph"]["bond_mode"] == "pairs" and not config["graph"].get("pair_file"):
        raise ValueError("--bond-mode pairs requires --pairs PAIRS.txt or graph.pair_file in config.yaml.")
    if args.pair_id:
        config["graph"]["pair_id"] = args.pair_id
    if args.workers is not None:
        config["parallel"]["workers"] = args.workers
    if getattr(args, "output_layout", None):
        config["output"]["structure_layout"] = args.output_layout
    if getattr(args, "no_info", False):
        config["output"]["write_info"] = False
    if args.no_gro:
        config["output"]["write_gro"] = False
    if getattr(args, "no_ring_gro", False):
        config["output"]["write_ring_gro"] = False
    if getattr(args, "no_half_cage_gro", False):
        config["output"]["write_half_cage_gro"] = False
    if getattr(args, "no_quasi_cage_gro", False):
        config["output"]["write_quasi_cage_gro"] = False
    if getattr(args, "no_cage_gro", False):
        config["output"]["write_cage_gro"] = False
    if getattr(args, "no_ice_gro", False):
        config["output"]["write_ice_gro"] = False
    if args.no_xlsx:
        config["output"]["write_xlsx_summary"] = False


def cage_sizes_need_other_outputs(value: Any) -> bool:
    """Enable generated cage targets when CLI cage sizes include 4-ring faces."""
    try:
        return 4 in resolve_size_list(value, fallback=[], key="cage.ring_sizes")
    except ValueError:
        return False


def analyze_frame(
    frame: Frame,
    config: dict[str, Any],
    stage_callback: Callable[[str], None] | None = None,
) -> FrameResult:
    """Analyze one frame and return all topology objects for export."""
    report_stage(stage_callback, "resolving settings")
    ring_sizes = resolve_size_list(config["ring"]["sizes"], fallback=[], key="ring.sizes")
    quasi_base_sizes = resolve_size_list(config["quasi_cage"].get("base_sizes", "auto"), fallback=ring_sizes, key="quasi_cage.base_sizes")
    quasi_side_sizes = resolve_size_list(config["quasi_cage"].get("side_sizes", "auto"), fallback=ring_sizes, key="quasi_cage.side_sizes")
    cage_ring_sizes = resolve_size_list(config["cage"].get("ring_sizes", [5, 6]), fallback=ring_sizes, key="cage.ring_sizes")
    report_stage(stage_callback, "selecting molecules")
    waters = select_waters(
        frame.atoms,
        resnames=set(config["water"]["resnames"]),
        oxygen_names=set(config["water"]["oxygen_names"]),
        hydrogen_names=set(config["water"]["hydrogen_names"]),
    )
    guests = select_guests(
        frame.atoms,
        resnames=set(config["guest"]["resnames"]),
        center_atoms=config["guest"].get("center_atoms", {}),
    )
    # Ring, half-cage, quasi-cage, cage, F3/F4, and ice all consume the same water graph.
    report_stage(stage_callback, "building water graph")
    graph = build_water_graph(
        frame.atoms,
        waters,
        frame.box,
        bond_mode=config["graph"]["bond_mode"],
        oo_cutoff_nm=float(config["graph"]["oo_cutoff_nm"]),
        hbond_distance_nm=float(config["graph"]["hbond_distance_nm"]),
        hbond_angle_deg=float(config["graph"]["hbond_angle_deg"]),
        pair_file=config["graph"].get("pair_file"),
        pair_id=str(config["graph"].get("pair_id", "resid")),
    )
    report_stage(stage_callback, "searching rings")
    rings = find_rings(
        graph.adjacency,
        sizes=ring_sizes,
        chordless=bool(config["ring"]["chordless"]),
    )
    report_stage(stage_callback, "searching half/quasi cage")
    half_cages, quasi_cages = find_cage_patches(
        frame,
        rings,
        enabled=bool(config["quasi_cage"].get("enabled", False)),
        base_sizes=quasi_base_sizes,
        side_sizes=quasi_side_sizes,
        max_combinations_per_base=int(config["quasi_cage"].get("max_combinations_per_base", 50000)),
        max_layers=int(config["quasi_cage"].get("max_layers", 1)),
        max_rings_per_layer=int(config["quasi_cage"].get("max_rings_per_layer", config["quasi_cage"].get("max_outer_layer_rings", 6))),
        max_layer_states_per_seed=int(config["quasi_cage"].get("max_layer_states_per_seed", 200)),
        max_candidates_per_edge=int(config["quasi_cage"].get("max_candidates_per_edge", 4)),
        max_layer_candidates=int(config["quasi_cage"].get("max_layer_candidates", 24)),
    )
    warnings = []
    cage_seed_patches = [*half_cages, *quasi_cages]
    report_stage(stage_callback, "searching cage")
    cages = find_cages(
        frame,
        rings,
        cage_seed_patches,
        guests,
        enabled=bool(config["cage"].get("enabled", False)),
        target_types=list(config["cage"].get("target_types", [])),
        ring_sizes=cage_ring_sizes,
        output_other=bool(config["cage"].get("output_other", False)),
        other_max_faces=int(config["cage"].get("other_max_faces", 20)),
        search_mode=str(config["cage"].get("search_mode", "grow")),
        seed_mode=str(config["cage"].get("seed_mode", "patch")),
        max_states_per_seed=int(config["cage"].get("max_states_per_seed", 20000)),
        max_total_states=int(config["cage"].get("max_total_states", 5000000)),
        max_boundary_candidates=int(config["cage"].get("max_boundary_candidates", 8)),
        occupancy_radius_nm=float(config["cage"].get("occupancy_radius_nm", 0.5)),
        occupancy_mode=str(config["cage"].get("occupancy_mode", "polyhedron")),
        warnings=warnings,
    )
    report_stage(stage_callback, "filtering free patches")
    quasi_cages = filter_free_patches(quasi_cages, cages)
    half_cages = filter_free_patches(half_cages, cages, higher_priority_patches=quasi_cages)
    focus_resids = {int(item) for item in config["order"].get("focus_waters", [])}
    report_stage(stage_callback, "computing F3/F4")
    f3f4 = compute_f3f4(frame, waters, graph, focus_resids=focus_resids) if bool(config["order"].get("f3f4_enabled", True)) else None
    report_stage(stage_callback, "classifying ice")
    ice_classes = classify_ice_waters(
        graph,
        waters,
        rings,
        enabled=bool(config["ice"].get("enabled", False)),
        min_six_rings=int(config["ice"].get("min_six_rings", 2)),
        require_four_coord_neighbors=bool(config["ice"].get("require_four_coord_neighbors", True)),
    )
    if not half_cages and not quasi_cages:
        warnings.append("No half_cage or quasi_cage was found with the current layered patch criteria.")
    return FrameResult(
        frame=frame,
        waters=waters,
        guests=guests,
        graph=graph,
        rings=rings,
        half_cages=half_cages,
        quasi_cages=quasi_cages,
        cages=cages,
        f3f4=f3f4,
        ice_like_waters=ice_classes.ice_like,
        ice_i_waters=ice_classes.ice_i,
        interfacial_ice_waters=ice_classes.interfacial,
        warnings=warnings,
    )


def filter_free_patches(
    patches: list[CagePatch],
    cages: list[Cage],
    higher_priority_patches: list[CagePatch] | None = None,
) -> list[CagePatch]:
    """Remove patches consumed by cages or already reported as a higher-priority class."""
    cage_ring_sets = [set(cage.rings) for cage in cages]
    higher_priority_ring_sets = [set(patch.rings) for patch in higher_priority_patches or []]
    free_patches = []
    for patch in patches:
        patch_rings = set(patch.rings)
        if any(patch_rings <= cage_rings for cage_rings in cage_ring_sets):
            continue
        if any(patch_rings < higher_priority_rings for higher_priority_rings in higher_priority_ring_sets):
            continue
        free_patches.append(patch)
    return free_patches


def resolve_size_list(value: Any, fallback: list[int], key: str) -> list[int]:
    """Resolve ring-size settings, allowing patch sizes to follow ring sizes."""
    if value in (None, "", "auto"):
        if not fallback:
            raise ValueError(f"{key} cannot be auto without a fallback size list.")
        return list(fallback)
    if isinstance(value, str):
        parts = [part.strip() for part in value.split(",") if part.strip()]
        if not parts:
            raise ValueError(f"{key} must contain at least one ring size.")
        return sorted({int(part) for part in parts})
    try:
        sizes = sorted({int(size) for size in value})
    except TypeError as exc:
        raise ValueError(f"{key} must be a list of integers or 'auto'.") from exc
    if not sizes:
        raise ValueError(f"{key} must contain at least one ring size.")
    return sizes
