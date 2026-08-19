"""Serial and parallel terminal progress displays."""

from __future__ import annotations

from collections import Counter
import shutil
import sys
from pathlib import Path
from threading import Event, Lock, Thread
from time import perf_counter
from typing import Any, Callable

from .formatting import TERMINAL_LABEL_WIDTH, format_seconds, terminal_field_line, write_terminal_block


PROGRESS_BAR_WIDTH = 25
PARALLEL_FILE_PREVIEW_LIMIT = 5
PARALLEL_FILE_COLUMN_WIDTH = 25
PARALLEL_ACTIVE_STAGE_WIDTH = 30
PROGRESS_RENDER_INTERVAL_SECONDS = 0.10

STAGE_GROUPS = (
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
        ("classifying hydrate cluster", "cluster"),
    ),
    (
        ("filtering free patches", "filtering"),
        ("computing order parameters", "order"),
        ("classifying ice", "ice"),
        ("writing outputs", "output"),
    ),
)
CPP_STAGE_GROUPS = (
    (
        ("reading frame", "reading"),
        ("resolving settings", "settings"),
        ("selecting molecules", "selecting"),
    ),
    (
        ("building water graph", "graph"),
        ("searching rings", "ring"),
        ("searching cage", "cage"),
    ),
    (
        ("computing order parameters", "order"),
        ("writing outputs", "output"),
    ),
)
STAGE_LABEL_BY_NAME = {
    stage: label for group in STAGE_GROUPS for stage, label in group
}


def configured_stage_groups(
    include_cluster_stage: bool,
    cpp_mode: bool = False,
    include_patch_stage: bool = True,
) -> list[list[tuple[str, str]]]:
    """Return stages applicable to the effective analysis configuration."""
    groups = CPP_STAGE_GROUPS if cpp_mode else STAGE_GROUPS
    return [
        [
            (stage, label)
            for stage, label in group
            if (include_cluster_stage or label != "cluster")
            and (include_patch_stage or label != "half/quasi")
        ]
        for group in groups
    ]


def stage_column_widths(rows: list[list[str]]) -> list[int]:
    """Measure the widest visible cell for each progress-display column."""
    column_count = max((len(row) for row in rows), default=0)
    return [
        max((len(row[index]) for row in rows if index < len(row)), default=0)
        for index in range(column_count)
    ]


ACTIVE_STAGE_ANSI = f"{chr(27)}[1;38;2;0;0;255m"
ANSI_RESET = f"{chr(27)}[0m"


def format_stage_label(label: str, width: int, *, active: bool, bold: bool) -> str:
    """Format a stage cell without counting its ANSI highlight as width."""
    padding = " " * max(width - len(label), 0)
    if active and bold:
        return f"{ACTIVE_STAGE_ANSI}{label}{ANSI_RESET}{padding}"
    return label + padding


def print_output_write_status() -> None:
    """Tell the user that final output publication is still in progress."""
    write_terminal_block(["Writing output files; please wait and do not close SQQ..."])


def print_output_directory_notice(
    requested: str | Path,
    resolved: str | Path,
    *,
    auto_renamed: bool,
) -> None:
    """Report the one intentional output-root rename without extra UI rows."""
    if not auto_renamed:
        return
    write_terminal_block(
        [
            "Output directory is not empty; using "
            f"{Path(resolved)} [auto-renamed] instead of {Path(requested)}."
        ]
    )


def compact_worker_policy(value: Any) -> str:
    """Shorten the resolved worker policy for the live Execution row."""
    text = " ".join(str(value or "auto").split())
    text = text.replace("1 workers", "1 worker")
    normalized = text.casefold()
    if normalized.startswith("auto"):
        return "auto; reserve 1 core"
    if normalized.startswith("explicit"):
        inner = text[text.find("(") + 1 : text.rfind(")")]
        inner = inner.replace("reserve 1 physical core", "reserve 1 core")
        return inner or "explicit"
    if normalized.startswith("engine default"):
        return text.replace("engine default", "engine default", 1)
    return text.replace("reserve 1 physical core", "reserve 1 core")


class RunProgressDisplay:
    """Render one stable progress panel or append-only static checkpoints."""

    def __init__(
        self,
        total: int,
        total_started_at: float,
        include_cluster_stage: bool,
        cpp_mode: bool = False,
        include_patch_stage: bool = True,
        unit: str = "Frames",
        execution: str | None = None,
    ) -> None:
        self.total = total
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(
            include_cluster_stage, cpp_mode, include_patch_stage
        )
        self.completed = 0
        self.failed = 0
        self.unit = str(unit or "Frames")
        self.execution = str(execution).strip() if execution else ""
        self.current_index: int | None = None
        self.current_file = "waiting"
        self.stage = "waiting"
        self.frame_started_at = perf_counter()
        self.stage_started_at = perf_counter()
        self._stream = sys.stdout
        self._lock = Lock()
        self._close_lock = Lock()
        self._stop_event = Event()
        self._closed = False
        self._rendered_lines = 0
        self._interactive = bool(getattr(self._stream, "isatty", lambda: False)())
        self._thread: Thread | None = None
        self._last_render_at = float("-inf")
        self._last_panel_lines: tuple[str, ...] | None = None
        self._last_static_state: tuple[int, int] | None = None
        # A one-task run has no useful aggregate progress information.  Keep
        # accepting stage callbacks for a uniform runner contract, but emit no
        # heading, Execution row, one-item table, bar, or closing blank line.
        self._visible = self.total > 1
        if self._visible:
            self._render(force=True)
        if self._visible and self._interactive:
            self._thread = Thread(target=self._tick, daemon=True)
            self._thread.start()

    def start_frame(self, frame_index: int, frame_name: str) -> Callable[[str], None]:
        with self._lock:
            if self._closed:
                return lambda stage: None
            now = perf_counter()
            self.current_index = frame_index
            self.current_file = frame_name
            self.stage = "reading frame"
            self.frame_started_at = now
            self.stage_started_at = now
            self._render_locked()
        return self.update_stage

    def update_stage(self, stage: str) -> None:
        with self._lock:
            if self._closed:
                return
            if stage != self.stage:
                self.stage = stage
                self.stage_started_at = perf_counter()
            self._render_locked()

    def complete_frame(self, success: bool) -> None:
        with self._lock:
            if self._closed:
                return
            self.completed += 1
            if not success:
                self.failed += 1
            self._render_locked(force=True)

    def close(self) -> None:
        with self._close_lock:
            if self._closed:
                return
            self._stop_event.set()
            if self._thread is not None:
                self._thread.join(timeout=0.2)
            with self._lock:
                if self._closed:
                    return
                if self._visible:
                    self._render_locked(force=True)
                    self._stream.write("\n")
                    self._stream.flush()
                self._closed = True

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render(force=True)

    def _render(self, *, force: bool = False) -> None:
        with self._lock:
            if not self._closed:
                self._render_locked(force=force)

    def _render_locked(self, *, force: bool = False) -> None:
        if self._closed or not self._visible:
            return
        if not self._interactive:
            state = (self.completed, self.failed)
            if state == self._last_static_state:
                return
            if self._last_static_state is not None and self.completed < self.total:
                return
            if self._last_static_state is None:
                self._stream.write("Analysis Progress\n")
                if self.execution:
                    self._stream.write(
                        terminal_field_line("Execution", self.execution) + "\n"
                    )
                if self.total > 1:
                    self._stream.write(
                        terminal_field_line(
                            self.unit,
                            f"0 / {self.total} [0 failed]",
                        )
                        + "\n"
                    )
            self._stream.flush()
            self._last_static_state = state
            return
        now = perf_counter()
        if not force and now - self._last_render_at < PROGRESS_RENDER_INTERVAL_SECONDS:
            return
        self._last_render_at = now
        lines = tuple(self._panel_lines())
        if lines == self._last_panel_lines:
            return
        if self._rendered_lines:
            self._stream.write(f"{chr(27)}[{self._rendered_lines}F")
        for line in lines:
            self._stream.write(chr(13) + f"{chr(27)}[K" + line + chr(10))
        self._stream.flush()
        self._rendered_lines = len(lines)
        self._last_panel_lines = lines

    def _panel_lines(self) -> list[str]:
        stage_lines = self._stage_lines()
        connector_indent = " " * (TERMINAL_LABEL_WIDTH + 2)
        height = shutil.get_terminal_size(fallback=(120, 40)).lines
        if height < 24:
            lines = ["Analysis Progress"]
            if self.execution:
                lines.append(terminal_field_line("Execution", self.execution))
            lines.extend([
                terminal_field_line("Stage", STAGE_LABEL_BY_NAME.get(self.stage, self.stage)),
                terminal_field_line("Stage / total", self._compact_time_text()),
            ])
        else:
            lines = ["Analysis Progress"]
            if self.execution:
                lines.append(terminal_field_line("Execution", self.execution))
            if self.total > 1:
                item_label = "File" if self.unit.casefold() == "files" else "Frame"
                lines.append(terminal_field_line(item_label, self._current_frame_text()))
            lines.extend([
                f"  {'Stage':<{TERMINAL_LABEL_WIDTH}}: {stage_lines[0]}",
                connector_indent + stage_lines[1],
                connector_indent + stage_lines[2],
                terminal_field_line(
                    "Stage / frame / total" if self.total > 1 else "Stage / total",
                    self._time_text() if self.total > 1 else self._compact_time_text(),
                ),
            ])
        if self.total > 1:
            lines.extend(["", self._files_bar()])
        return lines

    def _current_frame_text(self) -> str:
        if self.current_index is None:
            return "waiting"
        return f"{self.current_index + 1} / {self.total}"

    def _time_text(self) -> str:
        now = perf_counter()
        stage_elapsed = now - self.stage_started_at
        frame_elapsed = now - self.frame_started_at if self.current_index is not None else 0.0
        total_elapsed = now - self.total_started_at
        return f"{format_seconds(stage_elapsed)} / {format_seconds(frame_elapsed)} / {format_seconds(total_elapsed)}"

    def _compact_time_text(self) -> str:
        now = perf_counter()
        return (
            f"{format_seconds(now - self.stage_started_at)} / "
            f"{format_seconds(now - self.total_started_at)}"
        )

    def _stage_lines(self) -> list[str]:
        labels_by_row = [[label for _, label in group] for group in self.stage_groups]
        widths = stage_column_widths(labels_by_row)
        return [
            self._stage_flow_line(group, widths, row_index)
            for row_index, group in enumerate(self.stage_groups)
        ]

    def _stage_flow_line(
        self,
        group: list[tuple[str, str]],
        widths: list[int],
        row_index: int,
    ) -> str:
        cells = [
            format_stage_label(
                label,
                widths[index],
                active=stage == self.stage,
                bold=self._interactive,
            )
            for index, (stage, label) in enumerate(group)
        ]
        line = " > ".join(cells).rstrip()
        return "> " + line if row_index > 0 else line

    def _files_bar(self) -> str:
        fraction = 1.0 if self.total <= 0 else min(max(self.completed / self.total, 0.0), 1.0)
        filled = int(round(PROGRESS_BAR_WIDTH * fraction))
        bar = chr(9608) * filled + " " * (PROGRESS_BAR_WIDTH - filled)
        failed = f" [{self.failed} failed]" if self.failed else ""
        return (
            f"{self.unit}: {fraction * 100:3.0f}%|{bar}| "
            f"{self.completed}/{self.total} completed{failed}"
        )


class ParallelRunProgressDisplay:
    """Render aggregate progress without duplicating completion records."""

    def __init__(
        self,
        total: int,
        workers: int,
        total_started_at: float,
        include_cluster_stage: bool,
        cpp_mode: bool = False,
        include_patch_stage: bool = True,
        math_threads: int = 1,
        policy: str = "auto; reserve 1 core",
        unit: str = "Files",
    ) -> None:
        self.total = total
        self.workers = workers
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(
            include_cluster_stage, cpp_mode, include_patch_stage
        )
        self.completed = 0
        self.failed = 0
        self.math_threads = max(1, int(math_threads))
        self.policy = str(policy or "auto; reserve 1 core")
        self.unit = str(unit or "Files")
        self._active: dict[int, dict[str, Any]] = {}
        self._finished: set[int] = set()
        self._stream = sys.stdout
        self._lock = Lock()
        self._close_lock = Lock()
        self._stop_event = Event()
        self._closed = False
        self._rendered_lines = 0
        self._interactive = bool(getattr(self._stream, "isatty", lambda: False)())
        self._thread: Thread | None = None
        self._last_render_at = float("-inf")
        self._last_panel_lines: tuple[str, ...] | None = None
        self._last_static_state: tuple[int, int] | None = None
        self._render(force=True)
        if self._interactive:
            self._thread = Thread(target=self._tick, daemon=True)
            self._thread.start()

    def start_file(
        self,
        frame_index: int,
        frame_name: str,
        started_at: float | None = None,
    ) -> Callable[[str], None]:
        with self._lock:
            if self._closed or frame_index in self._finished:
                return lambda stage: None
            now = perf_counter() if started_at is None else float(started_at)
            self._active[frame_index] = {
                "name": frame_name,
                "stage": "reading frame",
                "file_started_at": now,
                "stage_started_at": now,
            }
            self._render_locked(force=True)
        return lambda stage: self.update_stage(frame_index, stage)

    def update_stage(
        self,
        frame_index: int,
        stage: str,
        started_at: float | None = None,
    ) -> None:
        if stage == "done":
            return
        with self._lock:
            if self._closed:
                return
            state = self._active.get(frame_index)
            if state is None:
                return
            if stage != state["stage"]:
                state["stage"] = stage
                state["stage_started_at"] = (
                    perf_counter() if started_at is None else float(started_at)
                )
            self._render_locked()

    def complete_file(self, frame_index: int, success: bool) -> None:
        with self._lock:
            if self._closed or frame_index in self._finished:
                return
            self._active.pop(frame_index, None)
            self._finished.add(frame_index)
            self.completed += 1
            if not success:
                self.failed += 1
            self._render_locked(force=True)

    def close(self) -> None:
        with self._close_lock:
            if self._closed:
                return
            self._stop_event.set()
            if self._thread is not None:
                self._thread.join(timeout=0.2)
            with self._lock:
                if self._closed:
                    return
                self._render_locked(force=True)
                self._stream.write("\n")
                self._stream.flush()
                self._closed = True

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render(force=True)

    def _render(self, *, force: bool = False) -> None:
        with self._lock:
            if not self._closed:
                self._render_locked(force=force)

    def _render_locked(self, *, force: bool = False) -> None:
        if self._closed:
            return
        if not self._interactive:
            state = (self.completed, self.failed)
            if state == self._last_static_state:
                return
            if self._last_static_state is not None and self.completed < self.total:
                return
            if self._last_static_state is None:
                self._stream.write("Analysis Progress\n")
                self._stream.write(
                    terminal_field_line("Execution", self._execution_text()) + "\n"
                )
            self._stream.flush()
            self._last_static_state = state
            return
        now = perf_counter()
        if not force and now - self._last_render_at < PROGRESS_RENDER_INTERVAL_SECONDS:
            return
        self._last_render_at = now
        lines = tuple(self._panel_lines())
        if lines == self._last_panel_lines:
            return
        if self._rendered_lines:
            self._stream.write(f"{chr(27)}[{self._rendered_lines}F")
        for line in lines:
            self._stream.write(chr(13) + f"{chr(27)}[K" + line + chr(10))
        self._stream.flush()
        self._rendered_lines = len(lines)
        self._last_panel_lines = lines

    def _panel_lines(self) -> list[str]:
        stage_lines = self._stage_summary_lines()
        indent = " " * (TERMINAL_LABEL_WIDTH + 4)
        terminal_height = shutil.get_terminal_size(fallback=(120, 40)).lines
        lines = [
            "Analysis Progress",
            terminal_field_line("Execution", self._execution_text()),
        ]
        if terminal_height < 24:
            active_stages = ", ".join(
                f"{STAGE_LABEL_BY_NAME.get(stage, stage)}:{count}"
                for stage, count in self._stage_counts().items()
                if count
            ) or "waiting"
            lines.append(terminal_field_line("Stages", active_stages))
        else:
            lines.extend([
                f"  {'Stages':<{TERMINAL_LABEL_WIDTH}}: {stage_lines[0]}",
                indent + stage_lines[1],
                indent + stage_lines[2],
            ])
        lines.extend([
            terminal_field_line("Elapsed", format_seconds(perf_counter() - self.total_started_at)),
            "",
            f"    {'active files':<{PARALLEL_FILE_COLUMN_WIDTH}} {'stage':<{PARALLEL_ACTIVE_STAGE_WIDTH}} stage / file",
        ])
        active_items = sorted(self._active.items())
        height_limit = 1 if terminal_height < 24 else (3 if terminal_height < 32 else 5)
        preview_slots = min(PARALLEL_FILE_PREVIEW_LIMIT, height_limit, self.workers)
        now = perf_counter()
        for slot in range(preview_slots):
            lines.append(
                self._active_file_line(active_items[slot], now)
                if slot < len(active_items)
                else ""
            )
        if self.workers > preview_slots:
            overflow = max(0, len(active_items) - preview_slots)
            lines.append(f"    ... {overflow} additional active files" if overflow else "")
        lines.extend(["", self._files_bar()])
        return lines

    def _stage_counts(self) -> Counter[str]:
        return Counter(str(state["stage"]) for state in self._active.values())

    def _stage_summary_values(self) -> list[list[tuple[str, int]]]:
        counts = self._stage_counts()
        return [
            [(label, counts.get(stage, 0)) for stage, label in group]
            for group in self.stage_groups
        ]

    def _stage_summary_cell_rows(self) -> list[list[str]]:
        return [
            [f"{label}:{count}" for label, count in values]
            for values in self._stage_summary_values()
        ]

    def _stage_summary_lines(self) -> list[str]:
        cell_rows = self._stage_summary_cell_rows()
        widths = stage_column_widths(cell_rows)
        return [
            "  ".join(
                f"{cell:<{widths[index]}}" for index, cell in enumerate(row)
            ).rstrip()
            for row in cell_rows
        ]

    def _active_file_line(
        self, item: tuple[int, dict[str, Any]], now: float
    ) -> str:
        frame_index, state = item
        file_text = compact_terminal_text(
            f"{frame_index + 1}/{self.total}  {state['name']}",
            PARALLEL_FILE_COLUMN_WIDTH,
        )
        stage_text = compact_terminal_text(
            str(state["stage"]), PARALLEL_ACTIVE_STAGE_WIDTH
        )
        stage_elapsed = format_seconds(now - float(state["stage_started_at"]))
        file_elapsed = format_seconds(now - float(state["file_started_at"]))
        return (
            f"    {file_text:<{PARALLEL_FILE_COLUMN_WIDTH}} "
            f"{stage_text:<{PARALLEL_ACTIVE_STAGE_WIDTH}} "
            f"{stage_elapsed:>8} / {file_elapsed:>8}"
        )

    def _files_bar(self) -> str:
        fraction = 1.0 if self.total <= 0 else min(max(self.completed / self.total, 0.0), 1.0)
        filled = int(round(PROGRESS_BAR_WIDTH * fraction))
        bar = chr(9608) * filled + " " * (PROGRESS_BAR_WIDTH - filled)
        failed = f" [{self.failed} failed]" if self.failed else ""
        return (
            f"{self.unit}: {fraction * 100:3.0f}%|{bar}| "
            f"{self.completed}/{self.total} completed{failed}"
        )

    def _execution_text(self) -> str:
        process_word = "process" if self.workers == 1 else "processes"
        thread_word = "thread" if self.math_threads == 1 else "threads"
        return (
            f"{self.workers} {process_word} x {self.math_threads} {thread_word} "
            f"[{self.policy}]"
        )


def compact_terminal_text(text: str, width: int) -> str:
    """Truncate a live-panel field without changing its column width."""
    if len(text) <= width:
        return text
    if width <= 3:
        return text[:width]
    return text[: width - 3] + "..."


__all__ = [
    "ACTIVE_STAGE_ANSI",
    "ANSI_RESET",
    "CPP_STAGE_GROUPS",
    "PARALLEL_ACTIVE_STAGE_WIDTH",
    "PARALLEL_FILE_COLUMN_WIDTH",
    "PARALLEL_FILE_PREVIEW_LIMIT",
    "PROGRESS_BAR_WIDTH",
    "PROGRESS_RENDER_INTERVAL_SECONDS",
    "ParallelRunProgressDisplay",
    "RunProgressDisplay",
    "STAGE_GROUPS",
    "STAGE_LABEL_BY_NAME",
    "compact_terminal_text",
    "compact_worker_policy",
    "configured_stage_groups",
    "format_stage_label",
    "print_output_directory_notice",
    "print_output_write_status",
    "stage_column_widths",
]
