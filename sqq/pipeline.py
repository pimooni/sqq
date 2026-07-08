from __future__ import annotations

"""Top-level analysis pipeline for the SQQ command line."""

import os
import sys
from argparse import Namespace
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, as_completed, wait
from datetime import datetime
from multiprocessing import get_context
from pathlib import Path
from queue import Empty
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
from .core.cage import (
    CAGE_REPORT_GROUPS,
    TARGET_FACE_COUNTS,
    canonical_cage_type,
    find_cages,
    parse_cage_face_label,
)
from .core.f3f4 import compute_order_parameters, normalize_q_degree
from .core.dhop import compute_dhop_order
from .core.mcg import compute_mcg_order
from .core.graph import build_water_graph
from .core.hydrate_cluster import analyze_hydrate_clusters
from .core.ice import classify_ice_waters
from .core.quasi_cage import find_cage_patches
from .core.ring import find_rings
from .core.ring_topology import build_ring_topology_index
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
    write_order_parameter,
    write_frame_info,
    write_membership,
    write_summary,
    write_vmd_script,
)
from .io.trajectory import expand_inputs, read_frames, trajectory_frame_indices
from .models import Cage, CagePatch, Frame, FrameResult, HydrateOrderResult
from .parallel import (
    physical_cpu_count,
    initialize_file_worker,
    initialize_trajectory_worker,
    limited_math_threads,
    process_file_task,
    process_trajectory_batch_task,
    process_worker_cap,
)


PARALLEL_SUFFIXES = {".gro", ".xyz"}
BOND_MODE_DISPLAY_NAMES = {
    "auto": "auto",
    "hbond": "hydrogen bond",
    "oo": "O-O connectivity",
    "pairs": "user-defined pairs",
}


def analyze(args: Namespace) -> None:
    """Run SQQ analysis from parsed command-line arguments."""
    run_started_at = datetime.now().astimezone()
    started_at = perf_counter()
    config = load_config(Path(args.config) if args.config else None, mode=getattr(args, "mode", None))
    apply_cli_overrides(config, args)
    normalize_analysis_scopes(config)

    # Directory input follows a one-file-per-frame workflow.
    input_path = Path(args.input)
    pattern = args.pattern or config["input"]["pattern"]
    recursive = bool(args.recursive or config["input"]["recursive"])
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    topology = Path(args.topology) if args.topology else None
    coordinate_parallelizable = can_parallelize_paths(paths, topology)
    trajectory_parallelizable = can_parallelize_trajectory(paths, topology)
    parallel_backend = normalize_parallel_backend(config.get("parallel", {}).get("backend", "process"))
    trajectory_indexes: list[int] = []
    if coordinate_parallelizable:
        validate_unique_output_names(paths)
        work_items = len(paths)
    elif trajectory_parallelizable:
        trajectory_indexes = trajectory_frame_indices(
            paths[0],
            topology,
            stride=int(config["input"].get("xtc_stride", 1)),
        )
        work_items = len(trajectory_indexes)
    else:
        work_items = len(paths)
    parallelizable = coordinate_parallelizable or (
        trajectory_parallelizable and parallel_backend == "process"
    )
    workers = (
        resolve_workers(
            config["parallel"].get("workers"),
            work_items,
            mode=config.get("mode", "50"),
            backend=parallel_backend,
        )
        if parallelizable and parallel_backend != "serial"
        else 1
    )
    active_backend = parallel_backend if workers > 1 and parallelizable else "serial"
    print_run_header(args, config, input_path, outdir, paths, topology, workers, active_backend, run_started_at)
    if workers > 1 and coordinate_parallelizable:
        rows = analyze_paths_parallel(
            paths,
            outdir,
            config,
            workers=workers,
            backend=parallel_backend,
            strict=bool(args.strict),
            total_started_at=started_at,
        )
    elif workers > 1 and trajectory_parallelizable:
        rows = analyze_trajectory_processes(
            paths[0],
            topology,
            trajectory_indexes,
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
            total_frames=work_items if trajectory_parallelizable else None,
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
        active_backend,
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
    parallel_backend: str,
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
        "search_sizes": config["ring"]["sizes"],
        "ring_report_sizes": config["ring"]["report_sizes"],
        "quasi_cage_base_sizes": config["quasi_cage"].get("base_sizes", "auto"),
        "quasi_cage_side_sizes": config["quasi_cage"].get("side_sizes", "auto"),
        "cage_report_types": config["cage"].get("report_types", []),
        "max_cage_face": config["cage"].get("max_faces", 20),
        "cage_fast_closure": on_off_text(config["cage"].get("fast_closure", True)),
        "cage_scientific_validation": on_off_text(config["cage"].get("scientific_validation", False)),
        "hydrate_cluster": on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)),
        "cluster_min_cage": config.get("hydrate_cluster", {}).get("min_cage", 2),
        "cluster_detail": on_off_text(config.get("hydrate_cluster", {}).get("detail", False)),
        "hydrate_order": hydrate_order_config_text(config),
        "mcg3": on_off_text(config.get("hydrate_order", {}).get("mcg3_enabled", False)),
        "dhop30": on_off_text(config.get("hydrate_order", {}).get("dhop30_enabled", False)),
        "dhop_neighbor_cutoff_nm": config.get("hydrate_order", {}).get("dhop_neighbor_cutoff_nm", 0.35),
        "q_enabled": config["order"].get("q_enabled", True),
        "q_degree": config["order"].get("q_degree", [6, 12]),
        "q_neighbor_mode": config["order"].get("q_neighbor_mode", "graph"),
        "q_cutoff_nm": config["order"].get("q_cutoff_nm", 0.35),
        "q_n_neighbor": config["order"].get("q_n_neighbor", None),
        "output_layout": config["output"].get("structure_layout", "grouped"),
        "workers": workers,
        "parallel_backend": parallel_backend,
        "math_threads": int(config.get("parallel", {}).get("math_threads", 1)),
        "summary_xlsx": str((outdir / "summary.xlsx").resolve()),
        "summary_detail": str((outdir / str(config.get("output", {}).get("summary_detail_dir", "summary_detail"))).resolve()),
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
    parallel_backend: str,
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
    print_terminal_field("graph_mode", bond_mode_display_name(config["graph"]["bond_mode"]))
    print_terminal_field("search_sizes", config["ring"]["sizes"])
    print_terminal_field("ring_report_sizes", config["ring"]["report_sizes"])
    print_terminal_field("ring_definition", config["ring"].get("definition", "chordless"))
    print_terminal_field("quasi_cage_sizes", f"{config['quasi_cage'].get('base_sizes', 'auto')} / {config['quasi_cage'].get('side_sizes', 'auto')}")
    print_terminal_field("quasi_max_layer", config["quasi_cage"].get("max_layers", ""))
    print_terminal_field("quasi_search_policy", config["quasi_cage"].get("search_policy", "bounded"))
    print_terminal_field("cage_report_types", dashboard_cage_targets(config))
    print_terminal_field("max_cage_face", config["cage"].get("max_faces", 20))
    print_terminal_field("cage_fast_closure", on_off_text(config["cage"].get("fast_closure", True)))
    print_terminal_field("scientific_validation", on_off_text(config["cage"].get("scientific_validation", False)))
    print_terminal_field("hydrate_cluster", on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)))
    print_terminal_field("cluster_min_cage", config.get("hydrate_cluster", {}).get("min_cage", 2))
    print_terminal_field("cluster_detail", on_off_text(config.get("hydrate_cluster", {}).get("detail", False)))
    print_terminal_field("hydrate_order", hydrate_order_config_text(config))
    print_terminal_field("Q_l", q_config_text(config))
    print_terminal_field("output_layout", config["output"].get("structure_layout", "grouped"))
    print_terminal_field("worker_policy", worker_policy_text(config))
    print_terminal_field("parallel_backend", parallel_backend)
    print_terminal_field("math_threads", config.get("parallel", {}).get("math_threads", 1))
    print_terminal_field("workers", workers)
    print("")


TERMINAL_LABEL_WIDTH = 22
PROGRESS_BAR_WIDTH = 25


def print_terminal_field(label: str, value: Any) -> None:
    """Print one aligned terminal key-value row."""
    print(f"  {label:<{TERMINAL_LABEL_WIDTH}}: {safe_terminal_text(value)}")


def bond_mode_display_name(value: Any) -> str:
    """Return a readable terminal label without changing config identifiers."""
    mode = str(value)
    return BOND_MODE_DISPLAY_NAMES.get(mode, mode)


def hydrate_order_config_text(config: dict[str, Any]) -> str:
    """Render active MCG/DHOP settings for run metadata and the terminal header."""
    order = config.get("hydrate_order", {})
    active = []
    for key, label, default in (
        ("mcg1_enabled", "MCG-1", True),
        ("dhop35_enabled", "DHOP35", True),
        ("mcg3_enabled", "MCG-3", False),
        ("dhop30_enabled", "DHOP30", False),
    ):
        if bool(order.get(key, default)):
            active.append(label)
    return ",".join(active) if active else "disabled"


def q_config_text(config: dict[str, Any]) -> str:
    """Render Steinhardt Q_l settings for the run header."""
    order = config.get("order", {})
    if not bool(order.get("q_enabled", True)):
        return "disabled"
    n_neighbors = order.get("q_n_neighbor", None)
    n_text = "NULL" if n_neighbors in (None, "", "null", "NULL") else str(n_neighbors)
    degree = ",".join(str(item) for item in normalize_q_degree(order.get("q_degree", [6, 12])))
    return f"degree={degree}; mode={order.get('q_neighbor_mode', 'graph')}; cutoff={order.get('q_cutoff_nm', 0.35)} nm; n={n_text}"


def on_off_text(value: Any) -> str:
    """Render on/off settings using the CLI vocabulary."""
    return "on" if parse_on_off(value, "on/off setting") else "off"


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
STAGE_LABEL_BY_NAME = {
    stage: label
    for group in STAGE_GROUPS
    for stage, label in group
}


def configured_stage_groups(include_cluster_stage: bool) -> list[list[tuple[str, str]]]:
    """Return progress stages, hiding hydrate cluster when it is not enabled."""
    return [
        [
            (stage, label)
            for stage, label in group
            if include_cluster_stage or label != "cluster"
        ]
        for group in STAGE_GROUPS
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
    """Format one stage cell, keeping ANSI highlight from affecting visible width."""
    padding = " " * max(width - len(label), 0)
    if active and bold:
        return f"{ACTIVE_STAGE_ANSI}{label}{ANSI_RESET}{padding}"
    return label + padding


class RunProgressDisplay:
    """Render per-run progress with current stage, frame, and total timings."""

    def __init__(self, total: int, total_started_at: float, include_cluster_stage: bool) -> None:
        self.total = total
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(include_cluster_stage)
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
                sys.stdout.write(f"{chr(27)}[{self._rendered_lines}F")
            for line in lines:
                sys.stdout.write(chr(13) + f"{chr(27)}[K" + line + chr(10))
            sys.stdout.flush()
            self._rendered_lines = len(lines)
            return
        if self._progress is not None and hasattr(self._progress, "set_postfix_str"):
            self._progress.set_postfix_str(self._postfix_text(), refresh=True)

    def _panel_lines(self) -> list[str]:
        stage_lines = self._stage_lines()
        connector_indent = " " * (TERMINAL_LABEL_WIDTH + 2)
        return [
            "Analysis Progress",
            f"  {'completed_files':<{TERMINAL_LABEL_WIDTH}}: {self.completed} / {self.total}  [ {self.failed} failed ]",
            f"  {'current_file':<{TERMINAL_LABEL_WIDTH}}: {self._current_file_text()}",
            f"  {'stage':<{TERMINAL_LABEL_WIDTH}}: {stage_lines[0]}",
            connector_indent + stage_lines[1],
            connector_indent + stage_lines[2],
            f"  {'stage / frame / total':<{TERMINAL_LABEL_WIDTH}}: {self._time_text()}",
            "",
            self._files_bar(),
        ]

    def _postfix_text(self) -> str:
        return (
            f"completed_files: {self.completed} / {self.total} [ {self.failed} failed ]; "
            f"current_file: {self._current_file_text()}; "
            f"stage: {STAGE_LABEL_BY_NAME.get(self.stage, self.stage)}; "
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

    def _stage_lines(self) -> list[str]:
        labels_by_row = [[label for _, label in group] for group in self.stage_groups]
        widths = stage_column_widths(labels_by_row)
        return [
            self._stage_flow_line(group, widths, row_index)
            for row_index, group in enumerate(self.stage_groups)
        ]

    def _stage_flow_line(self, group: list[tuple[str, str]], widths: list[int], row_index: int) -> str:
        cells = []
        for index, (stage, label) in enumerate(group):
            cells.append(
                format_stage_label(
                    label,
                    widths[index],
                    active=stage == self.stage,
                    bold=self._interactive,
                )
            )
        line = " > ".join(cells).rstrip()
        if row_index > 0:
            return "> " + line
        return line

    def _files_bar(self) -> str:
        if self.total <= 0:
            fraction = 1.0
        else:
            fraction = min(max(self.completed / self.total, 0.0), 1.0)
        filled = int(round(PROGRESS_BAR_WIDTH * fraction))
        bar = chr(9608) * filled + " " * (PROGRESS_BAR_WIDTH - filled)
        return f"Files: {fraction * 100:3.0f}%|{bar}| {self.completed}/{self.total} completed"


PARALLEL_FILE_PREVIEW_LIMIT = 6
PARALLEL_FILE_COLUMN_WIDTH = 25
PARALLEL_ACTIVE_STAGE_WIDTH = 30


class ParallelRunProgressDisplay:
    """Render aggregate and per-file progress for concurrent frame analysis."""

    def __init__(self, total: int, workers: int, total_started_at: float, include_cluster_stage: bool) -> None:
        self.total = total
        self.workers = workers
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(include_cluster_stage)
        self.completed = 0
        self.failed = 0
        self._active: dict[int, dict[str, Any]] = {}
        self._finished: set[int] = set()
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

    def start_file(self, frame_index: int, frame_name: str, started_at: float | None = None) -> Callable[[str], None]:
        """Register an active file and return its stage callback."""
        with self._lock:
            if frame_index in self._finished:
                return lambda stage: None
            now = perf_counter() if started_at is None else float(started_at)
            self._active[frame_index] = {
                "name": frame_name,
                "stage": "reading frame",
                "file_started_at": now,
                "stage_started_at": now,
            }
            self._render_locked()
        return lambda stage: self.update_stage(frame_index, stage)

    def update_stage(self, frame_index: int, stage: str, started_at: float | None = None) -> None:
        """Update one active file without disturbing other worker states."""
        if stage == "done":
            return
        with self._lock:
            state = self._active.get(frame_index)
            if state is None:
                return
            if stage != state["stage"]:
                state["stage"] = stage
                state["stage_started_at"] = perf_counter() if started_at is None else float(started_at)
            self._render_locked()

    def complete_file(self, frame_index: int, success: bool) -> None:
        """Move one file from the active set into completed results."""
        with self._lock:
            if frame_index in self._finished:
                return
            self._active.pop(frame_index, None)
            self._finished.add(frame_index)
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
                sys.stdout.write(f"{chr(27)}[{self._rendered_lines}F")
            for line in lines:
                sys.stdout.write(chr(13) + f"{chr(27)}[K" + line + chr(10))
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
        stage_text = " / ".join("  ".join(row) for row in self._stage_summary_cell_rows())
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
        return [[(label, counts.get(stage, 0)) for stage, label in group] for group in self.stage_groups]

    def _stage_summary_cell_rows(self) -> list[list[str]]:
        return [
            [f"{label}:{count}" for label, count in values]
            for values in self._stage_summary_values()
        ]

    def _stage_summary_lines(self) -> list[str]:
        cell_rows = self._stage_summary_cell_rows()
        widths = stage_column_widths(cell_rows)
        return [
            "  ".join(f"{cell:<{widths[index]}}" for index, cell in enumerate(row)).rstrip()
            for row in cell_rows
        ]

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
        bar = chr(9608) * filled + " " * (PROGRESS_BAR_WIDTH - filled)
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
    total_frames: int | None = None,
) -> list[dict[str, Any]]:
    """Analyze frames in input order."""
    rows: list[dict[str, Any]] = []
    frames = read_frames(paths, topology=topology, xtc_stride=int(config["input"].get("xtc_stride", 1)))
    progress = RunProgressDisplay(
        total=int(total_frames if total_frames is not None else len(paths)),
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
    )
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
    backend: str,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Analyze independent coordinate files with the selected concurrency backend."""
    resolved_backend = normalize_parallel_backend(backend)
    if resolved_backend == "thread":
        return analyze_paths_threaded(paths, outdir, config, workers, strict, total_started_at)
    if resolved_backend != "process":
        raise ValueError("Parallel analysis requires backend=process or backend=thread.")
    return analyze_paths_processes(paths, outdir, config, workers, strict, total_started_at)


def analyze_paths_threaded(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    workers: int,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Compatibility backend using the legacy shared-memory thread pool."""
    rows_by_index: dict[int, dict[str, Any]] = {}
    progress = ParallelRunProgressDisplay(
        total=len(paths),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
    )
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


def analyze_paths_processes(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    workers: int,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Use spawned processes so CPU-bound topology search can use multiple cores."""
    rows_by_index: dict[int, dict[str, Any]] = {}
    progress = ParallelRunProgressDisplay(
        total=len(paths),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
    )
    context = get_context("spawn")
    stage_queue = context.Queue()
    math_threads = int(config.get("parallel", {}).get("math_threads", 1))
    try:
        with limited_math_threads(math_threads):
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=context,
                initializer=initialize_file_worker,
                initargs=(config, str(outdir), strict, stage_queue),
            ) as executor:
                task_iterator = iter(enumerate(paths))
                futures: dict[Any, int] = {}
                max_in_flight = process_in_flight_limit(workers)

                def fill_queue() -> None:
                    while len(futures) < max_in_flight:
                        try:
                            frame_index, path = next(task_iterator)
                        except StopIteration:
                            return
                        future = executor.submit(process_file_task, frame_index, str(path))
                        futures[future] = frame_index

                fill_queue()
                while futures:
                    drain_process_stage_events(stage_queue, progress)
                    done, _ = wait(set(futures), timeout=0.1, return_when=FIRST_COMPLETED)
                    for future in done:
                        expected_index = futures.pop(future)
                        try:
                            frame_index, row = future.result()
                        except Exception:
                            progress.complete_file(expected_index, False)
                            for queued in futures:
                                queued.cancel()
                            raise
                        rows_by_index[frame_index] = row
                        progress.complete_file(frame_index, row.get("status") == "ok")
                    fill_queue()
                drain_process_stage_events(stage_queue, progress)
    finally:
        progress.close()
        stage_queue.close()
        stage_queue.join_thread()
    return [rows_by_index[index] for index in sorted(rows_by_index)]


def analyze_trajectory_processes(
    trajectory: Path,
    topology: Path | None,
    raw_frame_indexes: list[int],
    outdir: Path,
    config: dict[str, Any],
    workers: int,
    strict: bool,
    total_started_at: float,
) -> list[dict[str, Any]]:
    """Analyze selected frames with one private MDAnalysis Universe per process."""
    if topology is None:
        raise ValueError("XTC/TRR process analysis requires a topology file.")
    rows_by_index: dict[int, dict[str, Any]] = {}
    progress = ParallelRunProgressDisplay(
        total=len(raw_frame_indexes),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
    )
    context = get_context("spawn")
    stage_queue = context.Queue()
    math_threads = int(config.get("parallel", {}).get("math_threads", 1))
    try:
        with limited_math_threads(math_threads):
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=context,
                initializer=initialize_trajectory_worker,
                initargs=(config, str(outdir), strict, stage_queue, str(trajectory), str(topology)),
            ) as executor:
                batch_iterator = iter(trajectory_task_batches(raw_frame_indexes, workers))
                futures: dict[Any, tuple[tuple[int, int], ...]] = {}
                max_in_flight = process_in_flight_limit(workers)

                def fill_queue() -> None:
                    while len(futures) < max_in_flight:
                        try:
                            batch = next(batch_iterator)
                        except StopIteration:
                            return
                        future = executor.submit(process_trajectory_batch_task, batch)
                        futures[future] = batch

                fill_queue()
                while futures:
                    drain_process_stage_events(stage_queue, progress)
                    done, _ = wait(set(futures), timeout=0.1, return_when=FIRST_COMPLETED)
                    for future in done:
                        batch = futures.pop(future)
                        try:
                            results = future.result()
                        except Exception:
                            for frame_index, _ in batch:
                                progress.complete_file(frame_index, False)
                            for queued in futures:
                                queued.cancel()
                            raise
                        for frame_index, row in results:
                            rows_by_index[frame_index] = row
                    fill_queue()
                drain_process_stage_events(stage_queue, progress)
    finally:
        progress.close()
        stage_queue.close()
        stage_queue.join_thread()
    return [rows_by_index[index] for index in sorted(rows_by_index)]

def process_in_flight_limit(workers: int) -> int:
    """Bound queued process tasks without reducing active worker capacity."""
    return max(1, int(workers)) * 3


def trajectory_task_batches(
    raw_frame_indexes: list[int],
    workers: int,
) -> list[tuple[tuple[int, int], ...]]:
    """Group adjacent selected frames into small ordered process tasks."""
    if not raw_frame_indexes:
        return []
    resolved_workers = max(1, int(workers))
    target_batches = resolved_workers * 4
    batch_size = max(1, min(8, (len(raw_frame_indexes) + target_batches - 1) // target_batches))
    indexed = list(enumerate(raw_frame_indexes))
    return [
        tuple(indexed[start : start + batch_size])
        for start in range(0, len(indexed), batch_size)
    ]


def drain_process_stage_events(stage_queue: Any, progress: ParallelRunProgressDisplay) -> None:
    """Apply every queued worker event in the terminal-owning main process."""
    while True:
        try:
            kind, frame_index, value, timestamp = stage_queue.get_nowait()
        except Empty:
            return
        if kind == "start":
            progress.start_file(frame_index, value, started_at=timestamp)
        elif kind == "stage":
            progress.update_stage(frame_index, value, started_at=timestamp)
        elif kind == "complete":
            progress.complete_file(frame_index, value == "ok")

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
    try:
        frame = next(iter(read_frames([path])))
    except Exception as exc:
        if strict:
            raise
        return frame_index, failed_row(path.stem, str(path), str(exc))
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
            ring_sizes=list(result.ring_report_sizes),
        )
    else:
        remove_optional_info_output(result, frame_dir)
    if config["output"]["write_tsv"]:
        write_membership(result, frame_dir)
        write_order_parameter(result, frame_dir)
    else:
        remove_optional_tsv_outputs(result, frame_dir, remove_membership=True, remove_order=not config["output"].get("write_order_tsv", False))
    if config["output"].get("write_order_tsv", False) and not config["output"].get("write_tsv", False):
        write_order_parameter(result, frame_dir)
    if config["output"]["write_vmd"]:
        write_vmd_script(result, frame_dir)
    else:
        remove_optional_vmd_output(result, frame_dir)
    if config["output"]["write_gro"]:
        layout = str(config["output"].get("structure_layout", "grouped"))
        if config["output"].get("write_ring_gro", True):
            write_ring_gro_files(
                result,
                frame_dir,
                write_empty=config["output"]["write_empty_files"],
                layout=layout,
                sizes=set(result.ring_report_sizes),
            )
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


def remove_optional_tsv_outputs(result: FrameResult, frame_dir: Path, *, remove_membership: bool = True, remove_order: bool = True) -> None:
    """Remove stale optional TSV files when TSV output is disabled."""
    suffixes = []
    if remove_membership:
        suffixes.append("membership")
    if remove_order:
        suffixes.extend(["f3f4", "order_parameter"])
    for suffix in suffixes:
        path = frame_dir / f"{result.frame.name}_{suffix}.tsv"
        path.unlink(missing_ok=True)


def remove_optional_vmd_output(result: FrameResult, frame_dir: Path) -> None:
    """Remove stale optional VMD helper scripts when disabled."""
    (frame_dir / f"{result.frame.name}_view.vmd.tcl").unlink(missing_ok=True)


def worker_policy_text(config: dict[str, Any]) -> str:
    """Describe the active automatic or explicit worker policy."""
    value = config.get("parallel", {}).get("workers", "auto")
    reserve_text = "reserve 1 physical core"
    if value in (None, "", "auto"):
        percent = int(round(mode_worker_fraction(config.get("mode", "50")) * 100))
        return f"auto ({percent}% of physical cores, {reserve_text})"
    try:
        request_text = describe_worker_request(value)
    except ValueError:
        request_text = str(value)
    return f"explicit ({request_text}, {reserve_text})"


def describe_worker_request(value: Any) -> str:
    """Render a user-facing worker request for terminal and workbook metadata."""
    text = str(value).strip()
    if text.endswith("%"):
        fraction = parse_worker_fraction(text[:-1], value) / 100.0
        return f"{format_worker_percent(fraction)} of physical cores"
    number = parse_worker_number(text, value)
    if number <= 1:
        return f"{format_worker_percent(number)} of physical cores"
    if not number.is_integer():
        raise worker_value_error()
    return f"{int(number)} workers"


def format_worker_percent(fraction: float) -> str:
    """Format a worker CPU fraction as a compact percentage string."""
    percent = fraction * 100.0
    return f"{percent:g}%"


def parse_worker_fraction(text: str, original: Any) -> float:
    """Parse a positive percentage number from a worker option."""
    number = parse_worker_number(text, original)
    if number <= 0:
        raise worker_value_error()
    return number


def parse_worker_number(text: str, original: Any) -> float:
    """Parse a positive worker option number."""
    try:
        number = float(text)
    except (TypeError, ValueError) as exc:
        raise worker_value_error() from exc
    if number <= 0:
        raise worker_value_error()
    return number


def worker_value_error() -> ValueError:
    """Build the standard worker-option validation error."""
    return ValueError("parallel.workers / --worker must be 'auto', a positive integer worker count, a fraction <= 1, or a percentage such as 50%.")


def normalize_parallel_backend(value: Any) -> str:
    """Normalize the supported serial, process, and compatibility thread backends."""
    backend = str(value or "process").strip().lower()
    if backend not in {"process", "thread", "serial"}:
        raise ValueError("parallel.backend / --parallel-backend must be process, thread, or serial.")
    return backend


def resolve_workers(
    value: Any,
    n_paths: int,
    mode: Any = "50",
    cpu_total: int | None = None,
    backend: str = "process",
) -> int:
    """Resolve workers from physical cores, task count, and platform caps."""
    physical_total = max(1, int(cpu_total if cpu_total is not None else physical_cpu_count()))
    usable_workers = max(1, physical_total - 1)
    if value in (None, "", "auto"):
        requested = max(1, int(physical_total * mode_worker_fraction(mode)))
    else:
        requested = resolve_explicit_worker_request(value, physical_total)
    requested = min(requested, usable_workers, max(1, n_paths))
    if normalize_parallel_backend(backend) == "process":
        platform_cap = process_worker_cap()
        if platform_cap is not None:
            requested = min(requested, platform_cap)
    return max(1, requested)


def resolve_explicit_worker_request(value: Any, physical_total: int) -> int:
    """Resolve -w/--worker: fractions/percentages or integer worker counts."""
    text = str(value).strip()
    if text.endswith("%"):
        percent = parse_worker_fraction(text[:-1], value)
        return max(1, int(physical_total * (percent / 100.0)))
    number = parse_worker_number(text, value)
    if number <= 1:
        return max(1, int(physical_total * number))
    if not number.is_integer():
        raise worker_value_error()
    return int(number)


def validate_unique_output_names(paths: list[Path]) -> None:
    """Reject standalone inputs whose stems would write the same frame directory."""
    grouped: dict[str, list[Path]] = {}
    for path in paths:
        grouped.setdefault(path.stem.casefold(), []).append(path)
    duplicates = [items for items in grouped.values() if len(items) > 1]
    if not duplicates:
        return
    details = "; ".join(", ".join(str(path) for path in items) for items in duplicates)
    raise ValueError(f"Independent input files must have unique stems; output directory collision: {details}")

def can_parallelize_paths(paths: list[Path], topology: Path | None) -> bool:
    """Return whether every input is an independent coordinate file."""
    if topology is not None:
        return False
    return bool(paths) and all(path.suffix.lower() in PARALLEL_SUFFIXES for path in paths)


def can_parallelize_trajectory(paths: list[Path], topology: Path | None) -> bool:
    """Parallelize frames when one indexed XTC/TRR trajectory has a topology."""
    return (
        topology is not None
        and len(paths) == 1
        and paths[0].suffix.lower() in {".xtc", ".trr"}
    )

def apply_cli_overrides(config: dict[str, Any], args: Namespace) -> None:
    """Apply command-line options after YAML/default configuration."""
    if args.pattern:
        config["input"]["pattern"] = args.pattern
    if args.recursive:
        config["input"]["recursive"] = True
    if getattr(args, "size", None):
        config["ring"]["sizes"] = args.size
        config["quasi_cage"]["base_sizes"] = args.size
        config["quasi_cage"]["side_sizes"] = args.size
    if getattr(args, "ring_size", None):
        config["ring"]["report_sizes"] = args.ring_size
    if getattr(args, "quasi_size", None):
        config["quasi_cage"]["base_sizes"] = args.quasi_size
        config["quasi_cage"]["side_sizes"] = args.quasi_size
    if getattr(args, "quasi_base_size", None):
        config["quasi_cage"]["base_sizes"] = args.quasi_base_size
    if getattr(args, "quasi_side_size", None):
        config["quasi_cage"]["side_sizes"] = args.quasi_side_size
    if getattr(args, "quasi_max_layer", None) is not None:
        if args.quasi_max_layer < 1:
            raise ValueError("--quasi-max-layer must be at least 1.")
        config["quasi_cage"]["max_layers"] = args.quasi_max_layer
    if getattr(args, "quasi_search_policy", None):
        config["quasi_cage"]["search_policy"] = args.quasi_search_policy
    if getattr(args, "ring_definition", None):
        config["ring"]["definition"] = args.ring_definition
    if getattr(args, "no_q", False):
        config["order"]["q_enabled"] = False
    if getattr(args, "q_degree", None):
        config["order"]["q_degree"] = args.q_degree
    if getattr(args, "q_neighbor_mode", None):
        config["order"]["q_neighbor_mode"] = args.q_neighbor_mode
    if getattr(args, "q_cutoff", None) is not None:
        if args.q_cutoff <= 0:
            raise ValueError("--q-cutoff must be positive.")
        config["order"]["q_cutoff_nm"] = args.q_cutoff
    if getattr(args, "q_n_neighbor", None) is not None:
        value = str(args.q_n_neighbor).strip()
        if value.lower() in {"null", "none", "auto"}:
            config["order"]["q_n_neighbor"] = None
        else:
            count = int(value)
            if count < 1:
                raise ValueError("--q-n-neighbor must be at least 1 or NULL.")
            config["order"]["q_n_neighbor"] = count
    if getattr(args, "mcg3", None):
        config["hydrate_order"]["mcg3_enabled"] = parse_on_off(args.mcg3, "--mcg3")
    if getattr(args, "dhop30", None):
        config["hydrate_order"]["dhop30_enabled"] = parse_on_off(args.dhop30, "--dhop30")
    if getattr(args, "cage_size", None):
        config["cage"]["report_types"] = args.cage_size
    if getattr(args, "max_cage_face", None) is not None:
        if args.max_cage_face < 1:
            raise ValueError("--max-cage-face must be at least 1.")
        config["cage"]["max_faces"] = args.max_cage_face
    if getattr(args, "cage_fast_closure", None):
        config["cage"]["fast_closure"] = parse_on_off(args.cage_fast_closure, "--cage-fast-closure")
    if getattr(args, "cage_scientific_validation", None):
        config["cage"]["scientific_validation"] = parse_on_off(
            args.cage_scientific_validation,
            "--cage-scientific-validation",
        )
    if getattr(args, "hydrate_cluster", None):
        config["hydrate_cluster"]["enabled"] = parse_on_off(args.hydrate_cluster, "--hydrate-cluster")
    if getattr(args, "cluster_min_cage", None) is not None:
        if args.cluster_min_cage < 1:
            raise ValueError("--cluster-min-cage must be at least 1.")
        config["hydrate_cluster"]["min_cage"] = args.cluster_min_cage
    if getattr(args, "cluster_detail", None):
        config["hydrate_cluster"]["detail"] = parse_on_off(args.cluster_detail, "--cluster-detail")
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
    if getattr(args, "parallel_backend", None):
        config["parallel"]["backend"] = args.parallel_backend
    if getattr(args, "worker", None) is not None:
        config["parallel"]["workers"] = args.worker
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
    if getattr(args, "no_summary_detail", False):
        config["output"]["write_summary_detail_csv"] = False
    if getattr(args, "cage_isomer_rows", None):
        config["output"]["cage_isomer_rows"] = args.cage_isomer_rows
    if getattr(args, "write_order_tsv", False):
        config["output"]["write_order_tsv"] = True


def normalize_analysis_scopes(config: dict[str, Any]) -> None:
    """Normalize search and report scopes before frames are analyzed."""
    search_sizes = resolve_size_list(config["ring"].get("sizes", []), fallback=[], key="ring.sizes")
    unsupported = set(search_sizes) - {4, 5, 6, 7}
    if unsupported:
        raise ValueError(f"ring.sizes / --size supports only 4, 5, 6, and 7; got {sorted(unsupported)}")
    ring_report_sizes = resolve_size_list(
        config["ring"].get("report_sizes", "auto"),
        fallback=search_sizes,
        key="ring.report_sizes",
    )
    if not set(ring_report_sizes) <= set(search_sizes):
        raise ValueError("ring.report_sizes / --ring-size must be a subset of ring.sizes / --size.")

    max_faces = int(config["cage"].get("max_faces", 20))
    if max_faces < 1:
        raise ValueError("cage.max_faces / --max-cage-face must be at least 1.")
    report_types = resolve_cage_report_types(
        config["cage"].get("report_types", []),
        search_sizes,
        max_faces,
    )
    config["ring"]["sizes"] = search_sizes
    config["ring"]["report_sizes"] = ring_report_sizes
    config["cage"]["report_types"] = "all" if report_types is None else list(report_types)
    config["cage"]["max_faces"] = max_faces
    cage = config.setdefault("cage", {})
    cage["fast_closure"] = parse_on_off(cage.get("fast_closure", True), "cage.fast_closure")
    fast_states = int(cage.get("fast_closure_max_states", 20000))
    if fast_states < 1:
        raise ValueError("cage.fast_closure_max_states must be at least 1.")
    cage["fast_closure_max_states"] = fast_states
    cage["scientific_validation"] = parse_on_off(
        cage.get("scientific_validation", False),
        "cage.scientific_validation",
    )
    planarity = float(cage.get("max_face_planarity_rms_nm", 0.06))
    edge_cv = float(cage.get("max_face_edge_cv", 0.35))
    min_volume = float(cage.get("min_cage_volume_nm3", 1.0e-6))
    if planarity < 0:
        raise ValueError("cage.max_face_planarity_rms_nm must be non-negative.")
    if edge_cv < 0:
        raise ValueError("cage.max_face_edge_cv must be non-negative.")
    if min_volume <= 0:
        raise ValueError("cage.min_cage_volume_nm3 must be positive.")
    cage["max_face_planarity_rms_nm"] = planarity
    cage["max_face_edge_cv"] = edge_cv
    cage["min_cage_volume_nm3"] = min_volume
    hydrate_cluster = config.setdefault("hydrate_cluster", {})
    hydrate_cluster["enabled"] = parse_on_off(hydrate_cluster.get("enabled", False), "hydrate_cluster.enabled")
    min_cage = int(hydrate_cluster.get("min_cage", 2))
    if min_cage < 1:
        raise ValueError("hydrate_cluster.min_cage / --cluster-min-cage must be at least 1.")
    hydrate_cluster["min_cage"] = min_cage
    hydrate_cluster["detail"] = parse_on_off(hydrate_cluster.get("detail", False), "hydrate_cluster.detail")
    ring_definition = str(config.get("ring", {}).get("definition", "chordless")).strip().lower()
    if ring_definition not in {"chordless", "shortest_path"}:
        raise ValueError("ring.definition / --ring-definition must be chordless or shortest_path.")
    config["ring"]["definition"] = ring_definition
    quasi_policy = str(config.get("quasi_cage", {}).get("search_policy", "bounded")).strip().lower()
    if quasi_policy not in {"bounded", "exact"}:
        raise ValueError("quasi_cage.search_policy / --quasi-search-policy must be bounded or exact.")
    config["quasi_cage"]["search_policy"] = quasi_policy
    parallel = config.setdefault("parallel", {})
    parallel["backend"] = normalize_parallel_backend(parallel.get("backend", "process"))
    math_threads = int(parallel.get("math_threads", 1))
    if math_threads < 1:
        raise ValueError("parallel.math_threads must be at least 1.")
    parallel["math_threads"] = math_threads
    output = config.setdefault("output", {})
    output["write_summary_detail_csv"] = parse_on_off(
        output.get("write_summary_detail_csv", True),
        "output.write_summary_detail_csv",
    )
    detail_dir = str(output.get("summary_detail_dir", "summary_detail")).strip() or "summary_detail"
    detail_path = Path(detail_dir)
    if detail_path.is_absolute() or ".." in detail_path.parts:
        raise ValueError("output.summary_detail_dir must be a relative directory inside the output folder.")
    output["summary_detail_dir"] = detail_dir
    cage_isomer_rows = str(output.get("cage_isomer_rows", "nonzero")).strip().lower()
    if cage_isomer_rows not in {"nonzero", "all"}:
        raise ValueError("output.cage_isomer_rows / --cage-isomer-rows must be nonzero or all.")
    output["cage_isomer_rows"] = cage_isomer_rows
    guest_center_mode = str(config.get("guest", {}).get("center_mode", "center_atom")).strip().lower()
    if guest_center_mode not in {"center_atom", "centroid", "auto"}:
        raise ValueError("guest.center_mode must be center_atom, centroid, or auto.")
    config["guest"]["center_mode"] = guest_center_mode
    config["order"]["q_degree"] = list(normalize_q_degree(config.get("order", {}).get("q_degree", [6, 12])))
    hydrate_order = config.setdefault("hydrate_order", {})
    for key, default in (
        ("mcg1_enabled", True),
        ("mcg3_enabled", False),
        ("dhop35_enabled", True),
        ("dhop30_enabled", False),
    ):
        hydrate_order[key] = parse_on_off(hydrate_order.get(key, default), f"hydrate_order.{key}")
    positive_values = (
        ("mcg_guest_cutoff_nm", 0.90),
        ("mcg_water_cutoff_nm", 0.60),
        ("dhop_neighbor_cutoff_nm", 0.35),
    )
    for key, default in positive_values:
        value = float(hydrate_order.get(key, default))
        if value <= 0:
            raise ValueError(f"hydrate_order.{key} must be positive.")
        hydrate_order[key] = value
    cone_angle = float(hydrate_order.get("mcg_cone_half_angle_deg", 45.0))
    if not 0 < cone_angle < 90:
        raise ValueError("hydrate_order.mcg_cone_half_angle_deg must be between 0 and 90.")
    hydrate_order["mcg_cone_half_angle_deg"] = cone_angle
    for key, default in (("mcg_min_waters", 5), ("dhop_min_qualified_neighbors", 3)):
        value = int(hydrate_order.get(key, default))
        if value < 1:
            raise ValueError(f"hydrate_order.{key} must be at least 1.")
        hydrate_order[key] = value
    planar_counts = sorted({int(value) for value in hydrate_order.get("dhop_planar_counts", [11, 12])})
    if not planar_counts or planar_counts[0] < 0:
        raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.")
    hydrate_order["dhop_planar_counts"] = planar_counts
    guest_names = [str(value).strip() for value in hydrate_order.get("mcg_guest_resnames", ["CH4", "MET"]) if str(value).strip()]
    if not guest_names:
        raise ValueError("hydrate_order.mcg_guest_resnames must contain at least one residue name.")
    hydrate_order["mcg_guest_resnames"] = guest_names


def parse_on_off(value: Any, key: str) -> bool:
    """Parse YAML booleans and CLI on/off strings consistently."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"on", "true", "yes", "1"}:
        return True
    if text in {"off", "false", "no", "0", "", "none"}:
        return False
    raise ValueError(f"{key} must be on/off or true/false.")


def resolve_cage_report_types(
    value: Any,
    search_sizes: list[int],
    max_faces: int,
) -> tuple[str, ...] | None:
    """Resolve report groups/types; auto/all return every cage in the search scope."""
    if isinstance(value, str):
        items = [item.strip() for item in value.split(",") if item.strip()]
    else:
        try:
            items = [str(item).strip() for item in value if str(item).strip()]
        except TypeError as exc:
            raise ValueError("cage.report_types / --cage-size must be a comma-separated list or 'all'.") from exc
    if not items:
        raise ValueError("cage.report_types / --cage-size must contain at least one cage type.")
    scope_keywords = {item.lower() for item in items} & {"auto", "all"}
    if scope_keywords:
        if len(items) != 1:
            raise ValueError("Use 'auto' or 'all' alone in cage.report_types / --cage-size.")
        return None

    expanded_items: list[str] = []
    for item in items:
        expanded_items.extend(CAGE_REPORT_GROUPS.get(item.upper(), (item,)))

    allowed_sizes = set(search_sizes) & {4, 5, 6}
    resolved: list[str] = []
    for item in expanded_items:
        cage_type = canonical_cage_type(item)
        counts = TARGET_FACE_COUNTS.get(cage_type) or parse_cage_face_label(cage_type)
        if counts is None:
            raise ValueError(f"Unable to resolve cage type: {item}")
        required_sizes = {size for size, count in counts.items() if count > 0}
        if not required_sizes <= allowed_sizes:
            missing = sorted(required_sizes - allowed_sizes)
            raise ValueError(
                f"Cage type {item} requires ring size(s) {missing}, which are absent from --size."
            )
        if sum(counts.values()) > max_faces:
            raise ValueError(
                f"Cage type {item} has {sum(counts.values())} faces, above --max-cage-face={max_faces}."
            )
        if cage_type not in resolved:
            resolved.append(cage_type)
    return tuple(resolved)


def select_reported_cages(cages: list[Cage], report_types: tuple[str, ...] | None) -> list[Cage]:
    """Filter detected cages for reports without changing topology filtering."""
    if report_types is None:
        return list(cages)
    allowed = set(report_types)
    return [cage for cage in cages if cage.cage_type in allowed]

def analyze_frame(
    frame: Frame,
    config: dict[str, Any],
    stage_callback: Callable[[str], None] | None = None,
) -> FrameResult:
    """Analyze one frame and return all topology objects for export."""
    report_stage(stage_callback, "resolving settings")
    normalize_analysis_scopes(config)
    ring_sizes = resolve_size_list(config["ring"]["sizes"], fallback=[], key="ring.sizes")
    ring_report_sizes = resolve_size_list(config["ring"].get("report_sizes", "auto"), fallback=ring_sizes, key="ring.report_sizes")
    quasi_base_sizes = resolve_size_list(config["quasi_cage"].get("base_sizes", "auto"), fallback=ring_sizes, key="quasi_cage.base_sizes")
    quasi_side_sizes = resolve_size_list(config["quasi_cage"].get("side_sizes", "auto"), fallback=ring_sizes, key="quasi_cage.side_sizes")
    cage_ring_sizes = [size for size in ring_sizes if size in {4, 5, 6}]
    cage_report_types = resolve_cage_report_types(
        config["cage"].get("report_types", []),
        ring_sizes,
        int(config["cage"].get("max_faces", 20)),
    )
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
        center_mode=str(config["guest"].get("center_mode", "center_atom")),
    )
    # All structure classifiers consume the same water graph.
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
        definition=str(config["ring"].get("definition", "chordless")),
    )
    scientific_validation = bool(config["cage"].get("scientific_validation", False))
    hydrate_cluster_enabled = bool(config.get("hydrate_cluster", {}).get("enabled", False))
    ring_topology = build_ring_topology_index(
        frame,
        rings,
        compute_face_quality=scientific_validation,
        compute_face_normals=hydrate_cluster_enabled,
    )
    warnings: list[str] = []
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
        topology_index=ring_topology,
        search_policy=str(config["quasi_cage"].get("search_policy", "bounded")),
        warnings=warnings,
    )
    cage_seed_patches = [*half_cages, *quasi_cages]
    report_stage(stage_callback, "searching cage")
    all_cages = find_cages(
        frame,
        rings,
        cage_seed_patches,
        guests,
        enabled=bool(config["cage"].get("enabled", False)),
        ring_sizes=cage_ring_sizes,
        max_faces=int(config["cage"].get("max_faces", 20)),
        search_mode=str(config["cage"].get("search_mode", "grow")),
        seed_mode=str(config["cage"].get("seed_mode", "ring")),
        max_states_per_seed=int(config["cage"].get("max_states_per_seed", 20000)),
        max_total_states=int(config["cage"].get("max_total_states", 5000000)),
        max_boundary_candidates=int(config["cage"].get("max_boundary_candidates", 8)),
        occupancy_radius_nm=float(config["cage"].get("occupancy_radius_nm", 0.5)),
        occupancy_mode=str(config["cage"].get("occupancy_mode", "polyhedron")),
        fast_closure=bool(config["cage"].get("fast_closure", True)),
        fast_closure_max_states=int(config["cage"].get("fast_closure_max_states", 20000)),
        scientific_validation=scientific_validation,
        max_face_planarity_rms_nm=float(config["cage"].get("max_face_planarity_rms_nm", 0.06)),
        max_face_edge_cv=float(config["cage"].get("max_face_edge_cv", 0.35)),
        min_cage_volume_nm3=float(config["cage"].get("min_cage_volume_nm3", 1.0e-6)),
        topology_index=ring_topology,
        warnings=warnings,
    )
    cages = select_reported_cages(all_cages, cage_report_types)
    hydrate_cluster_detail = bool(config.get("hydrate_cluster", {}).get("detail", False))
    if hydrate_cluster_enabled:
        report_stage(stage_callback, "classifying hydrate cluster")
        rings_by_id = ring_topology.ring_by_id
        ring_sizes_by_id = {ring_id: ring.size for ring_id, ring in rings_by_id.items()}
        hydrate_clusters, hydrate_motifs, hydrate_domains, isolated_cage_ids = analyze_hydrate_clusters(
            cages,
            min_cage=int(config.get("hydrate_cluster", {}).get("min_cage", 2)),
            ring_sizes=ring_sizes_by_id,
            frame=frame,
            rings_by_id=rings_by_id,
            face_geometries=ring_topology.face_geometries(),
        )
    else:
        hydrate_clusters, hydrate_motifs, hydrate_domains, isolated_cage_ids = [], [], [], ()
    report_stage(stage_callback, "filtering free patches")
    quasi_cages = filter_free_patches(quasi_cages, all_cages)
    half_cages = filter_free_patches(half_cages, all_cages, higher_priority_patches=quasi_cages)
    focus_resids = {int(item) for item in config["order"].get("focus_waters", [])}
    report_stage(stage_callback, "computing order parameters")
    f3f4 = compute_order_parameters(
        frame,
        waters,
        graph,
        f3f4_enabled=bool(config["order"].get("f3f4_enabled", True)),
        q_enabled=bool(config["order"].get("q_enabled", True)),
        q_neighbor_mode=str(config["order"].get("q_neighbor_mode", "graph")),
        q_cutoff_nm=float(config["order"].get("q_cutoff_nm", 0.35)),
        q_n_neighbor=config["order"].get("q_n_neighbor", None),
        q_degree=config["order"].get("q_degree", [6, 12]),
        focus_resids=focus_resids,
    )
    hydrate_order_config = config.get("hydrate_order", {})
    mcg1, mcg3 = compute_mcg_order(frame, waters, guests, hydrate_order_config)
    dhop35, dhop30 = compute_dhop_order(frame, waters, hydrate_order_config)
    hydrate_order = HydrateOrderResult(
        mcg1=mcg1,
        dhop35=dhop35,
        mcg3=mcg3,
        dhop30=dhop30,
    )
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
        ring_report_sizes=tuple(ring_report_sizes),
        half_cages=half_cages,
        quasi_cages=quasi_cages,
        cages=cages,
        all_cages=all_cages,
        cage_report_types=cage_report_types,
        hydrate_cluster_enabled=hydrate_cluster_enabled,
        hydrate_cluster_detail=hydrate_cluster_detail,
        hydrate_clusters=hydrate_clusters,
        hydrate_motifs=hydrate_motifs,
        hydrate_domains=hydrate_domains,
        isolated_cage_ids=isolated_cage_ids,
        f3f4=f3f4,
        hydrate_order=hydrate_order,
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
    """Remove consumed patches using ring-to-owner inverted indexes."""
    cage_ring_sets = [frozenset(cage.rings) for cage in cages]
    higher_priority_ring_sets = [frozenset(patch.rings) for patch in higher_priority_patches or []]
    cage_index = subset_owner_index(cage_ring_sets)
    higher_index = subset_owner_index(higher_priority_ring_sets)
    free_patches = []
    for patch in patches:
        patch_rings = frozenset(patch.rings)
        if is_subset_of_indexed_owner(patch_rings, cage_ring_sets, cage_index, strict=False):
            continue
        if is_subset_of_indexed_owner(patch_rings, higher_priority_ring_sets, higher_index, strict=True):
            continue
        free_patches.append(patch)
    return free_patches


def subset_owner_index(ring_sets: list[frozenset[str]]) -> dict[str, set[int]]:
    """Index candidate supersets by every ring they contain."""
    owners: dict[str, set[int]] = {}
    for index, ring_ids in enumerate(ring_sets):
        for ring_id in ring_ids:
            owners.setdefault(ring_id, set()).add(index)
    return owners


def is_subset_of_indexed_owner(
    ring_ids: frozenset[str],
    owners: list[frozenset[str]],
    index: dict[str, set[int]],
    *,
    strict: bool,
) -> bool:
    """Test subset ownership after narrowing candidates by the rarest ring."""
    if not ring_ids:
        return any((not strict) or bool(owner) for owner in owners)
    if any(ring_id not in index for ring_id in ring_ids):
        return False
    anchor = min(ring_ids, key=lambda ring_id: len(index[ring_id]))
    return any(
        ring_ids < owners[owner_index] if strict else ring_ids <= owners[owner_index]
        for owner_index in index[anchor]
    )

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
