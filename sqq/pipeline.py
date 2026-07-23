from __future__ import annotations

"""Top-level analysis pipeline for the SQQ command line."""

import math
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

from . import __version__
from .banner import SQQ_BANNER
from .config import (
    DEFAULT_MODE,
    is_cpp_mode,
    load_config,
    mode_display,
    mode_worker_count,
    mode_worker_fraction,
    normalize_cpp_output_types,
    normalize_order_parameters,
    normalize_output_types,
    order_parameter_display,
    output_enabled,
    output_type_display,
    q_degrees_from_order_parameters,
    validate_cpp_cli,
)
from .core.cpp_backend import analyze_frame_cpp
from .core.cage import (
    CAGE_REPORT_GROUPS,
    TARGET_FACE_COUNTS,
    canonical_cage_type,
    find_cages,
    parse_cage_face_label,
)
from .core.f3f4 import (
    compute_order_parameters,
    normalize_q_degree,
    normalize_q_neighbor_mode,
    resolve_q_neighbor_count,
)
from .core.dhop import compute_dhop_order
from .core.mcg import compute_mcg_order
from .core.graph import build_water_graph
from .display import graph_mode_display, ordered_unique_graph_modes
from .core.hydrate_cluster import analyze_hydrate_clusters
from .core.ice import classify_ice_waters
from .core.quasi_cage import find_cage_patches
from .core.ring import find_rings
from .core.ring_topology import build_ring_topology_index
from .core.selection import select_guests, select_waters
from .io.gro_writer import (
    write_cage_gro_files,
    write_half_cage_gro_files,
    write_hydrate_cluster_gro_files,
    write_ice_gro_file,
    write_quasi_cage_gro_files,
    write_ring_gro_files,
)
from .io.lammps import (
    LAMMPS_TRAJECTORY_SUFFIXES,
    inspect_lammps_topology_mapping,
    normalize_lammps_config,
)
from .io.summary import (
    dashboard_cage_targets,
    failed_row,
    result_row,
    write_order_parameter,
    write_frame_info,
    write_membership,
    write_run_config,
    write_summary,
    write_vmd_script,
)
from .io.trajectory import expand_inputs, read_frames, trajectory_frame_indices
from .io.vmd import (
    FRAGMENT_DIRECTORY,
    cleanup_sqq_cage_bundle,
    finalize_sqq_cage_bundle,
    prepare_sqq_cage_fragments,
    write_sqq_cage_fragment,
)
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
    validate_cpp_cli(args, config)
    normalize_analysis_scopes(config)

    # Directory inputs use one file per frame.
    input_path = Path(args.input)
    pattern = args.pattern or config["input"]["pattern"]
    recursive = bool(args.recursive or config["input"]["recursive"])
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    topology = Path(args.topology) if args.topology else None
    config["input"]["format"] = input_format_label(paths)
    config["input"]["topology"] = str(topology.resolve()) if topology else None
    if str(config["input"]["format"]).startswith("lammps-"):
        lammps_config = config["input"]["lammps"]
        explicit_type_map = bool(lammps_config.get("type_map"))
        resolved_type_map = {}
        rebuilt_molecules = False
        if topology is not None:
            resolved_type_map, rebuilt_molecules = inspect_lammps_topology_mapping(
                topology,
                lammps_config,
            )
        mapping_text = ", ".join(
            (
                f"{type_id}=ignore"
                if entry.ignore
                else f"{type_id}={entry.resname}/{entry.atomname}"
            )
            for type_id, entry in sorted(
                resolved_type_map.items(),
                key=lambda item: int(item[0]),
            )
        )
        lammps_config["resolved_type_map"] = {
            type_id: (
                {"ignore": True}
                if entry.ignore
                else {"resname": entry.resname, "atomname": entry.atomname}
            )
            for type_id, entry in sorted(
                resolved_type_map.items(),
                key=lambda item: int(item[0]),
            )
        }
        if explicit_type_map:
            source = str(Path(args.config).resolve()) if args.config else "<configuration>"
        else:
            source = "auto (DATA topology)"
            if rebuilt_molecules:
                source += "; molecule IDs rebuilt from Bonds"
        lammps_config["type_map_source"] = (
            f"{source}: {mapping_text}" if mapping_text else source
        )
    from .gro_batch import analyze_multi_gro_batch, is_multi_gro_batch

    if is_multi_gro_batch(paths):
        analyze_multi_gro_batch(
            args,
            config,
            input_path,
            paths,
            outdir,
            topology,
            run_started_at,
            started_at,
        )
        return

    coordinate_parallelizable = can_parallelize_paths(paths, topology)
    trajectory_parallelizable = can_parallelize_trajectory(paths, topology)
    parallel_backend = normalize_parallel_backend(config.get("parallel", {}).get("backend", "process"))
    trajectory_indexes: list[int] = []
    validate_unique_output_names(paths)
    if coordinate_parallelizable:
        work_items = len(paths)
    elif trajectory_parallelizable:
        trajectory_indexes = trajectory_frame_indices(
            paths[0],
            topology,
            stride=int(config["input"].get("trajectory_stride", 1)),
            lammps_config=config["input"].get("lammps", {}),
        )
        work_items = len(trajectory_indexes)
    else:
        work_items = len(paths)
    parallelizable = coordinate_parallelizable or (
        trajectory_parallelizable and parallel_backend == "process"
    )
    requested_workers = resolve_workers(
        config["parallel"].get("workers"),
        work_items,
        mode=config.get("mode", DEFAULT_MODE),
        backend=parallel_backend,
    )
    workers = (
        requested_workers
        if parallelizable and work_items > 1 and parallel_backend != "serial"
        else 1
    )
    active_backend = parallel_backend if workers > 1 and parallelizable else "serial"
    print_run_header(args, config, input_path, outdir, paths, topology, workers, active_backend, run_started_at)
    initial_run_info = build_run_info(
        args,
        config,
        input_path,
        outdir,
        paths,
        topology,
        workers,
        active_backend,
        0.0,
        run_started_at,
        run_started_at,
        [],
    )
    initial_run_info["status"] = "running"
    initial_run_info["error"] = ""
    write_run_config(outdir, config, initial_run_info)
    bundle_gro = output_enabled(config, "sqq-cage-gro")
    bundle_script = output_enabled(config, "sqq-render")
    cleanup_sqq_cage_bundle(outdir)
    if bundle_gro or bundle_script:
        prepare_sqq_cage_fragments(outdir)

    try:
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
        if bundle_gro or bundle_script:
            finalize_sqq_cage_bundle(
                outdir,
                write_gro=bundle_gro,
                write_script=bundle_script,
            )
    except Exception as exc:
        cleanup_sqq_cage_bundle(outdir)
        failed_at = datetime.now().astimezone()
        failed_run_info = build_run_info(
            args,
            config,
            input_path,
            outdir,
            paths,
            topology,
            workers,
            active_backend,
            perf_counter() - started_at,
            run_started_at,
            failed_at,
            [],
        )
        failed_run_info["status"] = "failed"
        failed_run_info["error"] = str(exc)
        write_run_config(outdir, config, failed_run_info)
        raise

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
        rows,
    )
    run_info["status"] = "completed"
    run_info["error"] = ""
    try:
        run_info["summary_write"] = write_summary(
            rows,
            outdir,
            config,
            write_xlsx=output_enabled(config, "summary-xlsx"),
            run_info=run_info,
        )
        # Rewrite once final write timing is known.
        write_run_config(outdir, config, run_info)
    except Exception as exc:
        run_info["status"] = "failed"
        run_info["error"] = str(exc)
        write_run_config(outdir, config, run_info)
        raise
    print_run_summary(run_info)
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
    rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Collect run-level metadata for terminal, config.yaml, and summaries."""
    requested_graph_mode = config["graph"]["bond_mode"]
    effective_graph_modes = row_effective_graph_modes(rows or [])
    selected_order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters")
    )
    q_degrees = q_degrees_from_order_parameters(selected_order_parameters)
    output_types = normalize_output_types(
        config.get("output", {}).get("types")
    )
    result_rows = rows or []
    failures = [
        {
            "frame": str(row.get("frame", "")),
            "source": str(row.get("source", "")),
            "error": str(row.get("error", "")),
        }
        for row in result_rows
        if str(row.get("status", "")).lower() == "failed"
    ]
    input_format = input_format_label(paths)
    info: dict[str, Any] = {
        "working_dir": str(Path.cwd()),
        "input": str(input_path),
        "input_format": input_format,
        "output_dir": str(outdir.resolve()),
        "date": started_at_wall.strftime("%Y-%m-%d"),
        "start_time": started_at_wall.strftime("%H:%M:%S"),
        "finish_time": finished_at_wall.strftime("%H:%M:%S"),
        "started_at": started_at_wall.isoformat(timespec="seconds"),
        "finished_at": finished_at_wall.isoformat(timespec="seconds"),
        "time_zone": format_time_zone(started_at_wall),
        "config_file": args.config or "<built-in defaults>",
        "sqq_version": __version__,
        "mode": mode_display(config.get("mode", DEFAULT_MODE)),
        "worker_policy": worker_policy_text(config),
        "topology": str(topology) if topology else "<none>",
        "matched_files": len(paths),
        "trajectory_stride": int(config["input"].get("trajectory_stride", 1)),
        "frames_total": len(result_rows),
        "frames_ok": sum(
            str(row.get("status", "")).lower() == "ok"
            for row in result_rows
        ),
        "frames_failed": len(failures),
        "failures": failures,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "graph_mode": requested_graph_mode,
        "effective_graph_modes": ", ".join(ordered_unique_graph_modes(effective_graph_modes)),
        "graph_mode_display": graph_mode_display(requested_graph_mode, effective_graph_modes),
        "search_sizes": config["ring"]["sizes"],
        "ring_report_sizes": config["ring"]["report_sizes"],
        "quasi_cage_base_sizes": config["quasi_cage"].get("base_sizes", "auto"),
        "quasi_cage_side_sizes": config["quasi_cage"].get("side_sizes", "auto"),
        "cage_report_types": config["cage"].get("report_types", []),
        "max_cage_face": config["cage"].get("max_faces", 20),
        "cage_fast_closure": on_off_text(config["cage"].get("fast_closure", True)),
        "cage_scientific_validation": on_off_text(config["cage"].get("scientific_validation", False)),
        "find_cluster": on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)),
        "cluster_min_cage": config.get("hydrate_cluster", {}).get("min_cage", 2),
        "order_parameters": order_parameter_display(selected_order_parameters),
        "hydrate_order": hydrate_order_config_text(config),
        "mcg3": on_off_text("mcg3" in selected_order_parameters),
        "dhop30": on_off_text("dhop30" in selected_order_parameters),
        "dhop_neighbor_cutoff_nm": config.get("hydrate_order", {}).get("dhop_neighbor_cutoff_nm", 0.35),
        "q_enabled": bool(q_degrees),
        "q_degree": list(q_degrees),
        "q_neighbor_mode": config["order"].get("q_neighbor_mode", "graph"),
        "q_cutoff_nm": config["order"].get("q_cutoff_nm", 0.35),
        "q_n_neighbor": config["order"].get("q_n_neighbor", None),
        "output_types": output_type_display(
            output_types,
            cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        ),
        "output_layout": config["output"].get("structure_layout", "grouped"),
        "workers": workers,
        "parallel_backend": parallel_backend,
        "math_threads": int(config.get("parallel", {}).get("math_threads", 1)),
        "summary_xlsx": (
            str((outdir / "summary.xlsx").resolve())
            if output_enabled(config, "summary-xlsx")
            else "<disabled>"
        ),
        "summary_csv": (
            str(
                (
                    outdir
                    / str(
                        config.get("output", {}).get(
                            "summary_csv_dir",
                            "summary",
                        )
                    )
                ).resolve()
            )
            if output_enabled(config, "summary-csv")
            else "<disabled>"
        ),
        "summary_detail_csv": (
            str(
                (
                    outdir
                    / str(
                        config.get("output", {}).get(
                            "summary_detail_dir",
                            "summary_detail",
                        )
                    )
                ).resolve()
            )
            if (
                output_enabled(config, "summary-detail-csv")
                or output_enabled(config, "cluster-detail")
            )
            else "<disabled>"
        ),
        "config_output": str((outdir / "config.yaml").resolve()),
    }
    if input_format.startswith("lammps-"):
        lammps = config["input"].get("lammps", {})
        info.update(
            {
                "lammps_units": lammps.get("units", "real"),
                "lammps_timestep": lammps.get("timestep", 1.0),
                "lammps_atom_style": lammps.get("atom_style", "full"),
                "lammps_type_map_source": lammps.get("type_map_source", "<configuration>"),
            }
        )
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
    print_terminal_field("input_format", input_format_label(paths))
    print_terminal_field("matched_files", len(paths))
    print_terminal_field("output", outdir)
    print("")
    print("Configuration")
    print_terminal_field("SQQ version", __version__)
    print_terminal_field("Mode", mode_display(config.get("mode", DEFAULT_MODE)))
    print_terminal_field("Config file", args.config or "<built-in defaults>")
    print_terminal_field("Topology", topology or "<none>")
    print_terminal_field("Trajectory stride", config["input"].get("trajectory_stride", 1))
    if input_format_label(paths).startswith("lammps-"):
        lammps = config["input"].get("lammps", {})
        print_terminal_field("LAMMPS units", lammps.get("units", "real"))
        print_terminal_field("LAMMPS timestep", lammps.get("timestep", 1.0))
        print_terminal_field("LAMMPS atom style", lammps.get("atom_style", "full"))
        print_terminal_field("LAMMPS type map", lammps.get("type_map_source", "<configuration>"))
    print_terminal_field("Graph mode", graph_mode_display(config["graph"]["bond_mode"]))
    print_terminal_field("Search sizes", config["ring"]["sizes"])
    print_terminal_field("Ring definition", config["ring"].get("definition", "chordless"))
    if not is_cpp_mode(config.get("mode")):
        print_terminal_field("Ring report sizes", config["ring"]["report_sizes"])
        print_terminal_field("Quasi-cage sizes", f"{config['quasi_cage'].get('base_sizes', 'auto')} / {config['quasi_cage'].get('side_sizes', 'auto')}")
        print_terminal_field("Quasi max layer", config["quasi_cage"].get("max_layers", ""))
        print_terminal_field("Quasi search policy", config["quasi_cage"].get("search_policy", "bounded"))
    print_terminal_field("Cage report types", dashboard_cage_targets(config))
    print_terminal_field("Maximum cage face", config["cage"].get("max_faces", 20))
    if not is_cpp_mode(config.get("mode")):
        print_terminal_field("Cage fast closure", on_off_text(config["cage"].get("fast_closure", True)))
    print_terminal_field("Scientific validation", on_off_text(config["cage"].get("scientific_validation", False)))
    if not is_cpp_mode(config.get("mode")):
        print_terminal_field("Find cluster", on_off_text(config.get("hydrate_cluster", {}).get("enabled", False)))
        print_terminal_field("Cluster min cage", config.get("hydrate_cluster", {}).get("min_cage", 2))
    print_terminal_field("Order parameters", order_parameter_config_text(config))
    if q_degrees_from_order_parameters(config.get("order", {}).get("parameters")):
        print_terminal_field("Q_l settings", q_config_text(config))
    print_terminal_field(
        "Output types",
        output_type_display(
            config.get("output", {}).get("types"),
            cpp_mode=is_cpp_mode(config.get("mode", DEFAULT_MODE)),
        ),
    )
    print_terminal_field("Output layout", config["output"].get("structure_layout", "grouped"))
    print_terminal_field("Worker policy", worker_policy_text(config))
    print_terminal_field("Parallel backend", parallel_backend)
    print_terminal_field("Math threads per worker", config.get("parallel", {}).get("math_threads", 1))
    print_terminal_field("Workers", workers)
    print("")


TERMINAL_LABEL_WIDTH = 24
PROGRESS_BAR_WIDTH = 25


def print_terminal_field(label: str, value: Any) -> None:
    """Print one aligned terminal key-value row."""
    print(f"  {label:<{TERMINAL_LABEL_WIDTH}}: {safe_terminal_text(value)}")


def print_run_summary(run_info: dict[str, Any]) -> None:
    """Print final effective run metadata after frame analysis finishes."""
    print("Run Summary")
    print_terminal_field("Finish time", run_info.get("finish_time", ""))
    print_terminal_field("Duration (s)", run_info.get("elapsed_seconds", ""))
    print_terminal_field("SQQ version", run_info.get("sqq_version", __version__))
    print_terminal_field("Mode", run_info.get("mode", ""))
    print_terminal_field("Graph mode", run_info.get("graph_mode_display", run_info.get("graph_mode", "")))
    print_terminal_field("Order parameters", run_info.get("order_parameters", ""))
    if "sqq-cpp" not in str(run_info.get("mode", "")).strip().lower():
        print_terminal_field("Find cluster", run_info.get("find_cluster", "off"))
    print_terminal_field("Output types", run_info.get("output_types", "none"))
    print_terminal_field("Worker policy", run_info.get("worker_policy", ""))
    print_terminal_field("Parallel backend", run_info.get("parallel_backend", "serial"))
    print_terminal_field("Workers", run_info.get("workers", ""))
    summary_write = run_info.get("summary_write", {})
    if isinstance(summary_write, dict) and "total_seconds" in summary_write:
        print_terminal_field("Summary write (s)", summary_write.get("total_seconds", ""))
    print("")


def row_effective_graph_modes(rows: list[dict[str, Any]]) -> list[str]:
    """Collect effective graph modes from successful per-frame summary rows."""
    modes: list[str] = []
    for row in rows:
        if str(row.get("status", "ok")).lower() == "failed":
            continue
        mode = str(row.get("connection_mode", "")).strip()
        if mode:
            modes.append(mode)
    return modes


def frame_input_metadata(config: dict[str, Any]) -> dict[str, Any]:
    """Return normalized input provenance for one frame report."""
    input_config = config.get("input", {})
    metadata = {
        "input_format": input_config.get("format", ""),
        "topology": input_config.get("topology"),
        "trajectory_stride": input_config.get("trajectory_stride", 1),
    }
    if str(metadata["input_format"]).startswith("lammps-"):
        lammps = input_config.get("lammps", {})
        metadata.update(
            {
                "lammps_units": lammps.get("units", "real"),
                "lammps_timestep": lammps.get("timestep", 1.0),
                "lammps_atom_style": lammps.get("atom_style", "full"),
                "lammps_type_map_source": lammps.get("type_map_source", "<configuration>"),
            }
        )
    return metadata


def input_format_label(paths: list[Path]) -> str:
    """Return one compact source-format label for run metadata."""
    mapping = {".gro": "gromacs-gro", ".xyz": "xyz", ".xtc": "gromacs-xtc", ".trr": "gromacs-trr", ".dump": "lammps-dump", ".lammpstrj": "lammps-dump", ".dcd": "lammps-dcd"}
    labels: list[str] = []
    for path in paths:
        label = mapping.get(path.suffix.lower(), path.suffix.lower().lstrip("."))
        if label not in labels:
            labels.append(label)
    if not labels:
        return "unknown"
    if len(labels) == 1:
        return labels[0]
    return "mixed (" + ", ".join(labels) + ")"


def bond_mode_display_name(value: Any) -> str:
    """Return a readable terminal label without changing config identifiers."""
    mode = str(value)
    return BOND_MODE_DISPLAY_NAMES.get(mode, mode)


def hydrate_order_config_text(config: dict[str, Any]) -> str:
    """Render the selected MCG/DHOP subset for compatibility metadata."""
    parameters = normalize_order_parameters(config.get("order", {}).get("parameters"))
    active = [
        name
        for name in parameters
        if name in {"mcg1", "mcg3", "dhop35", "dhop30"}
    ]
    return ", ".join(active) if active else "disabled"


def order_parameter_config_text(config: dict[str, Any]) -> str:
    """Render the unified order-parameter selection."""
    return order_parameter_display(config.get("order", {}).get("parameters"))


def q_config_text(config: dict[str, Any]) -> str:
    """Render Steinhardt Q_l settings for the run header."""
    order = config.get("order", {})
    degrees = q_degrees_from_order_parameters(order.get("parameters"))
    if not degrees:
        return "disabled"
    n_neighbors = order.get("q_n_neighbor", None)
    n_text = "NULL" if n_neighbors in (None, "", "null", "NULL") else str(n_neighbors)
    degree = ",".join(str(item) for item in degrees)
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
    stage: label
    for group in STAGE_GROUPS
    for stage, label in group
}


def configured_stage_groups(
    include_cluster_stage: bool,
    cpp_mode: bool = False,
) -> list[list[tuple[str, str]]]:
    """Return progress stages, hiding hydrate cluster when it is not enabled."""
    groups = CPP_STAGE_GROUPS if cpp_mode else STAGE_GROUPS
    return [
        [
            (stage, label)
            for stage, label in group
            if include_cluster_stage or label != "cluster"
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
    """Format one stage cell, keeping ANSI highlight from affecting visible width."""
    padding = " " * max(width - len(label), 0)
    if active and bold:
        return f"{ACTIVE_STAGE_ANSI}{label}{ANSI_RESET}{padding}"
    return label + padding


class RunProgressDisplay:
    """Render per-run progress with current stage, frame, and total timings."""

    def __init__(
        self,
        total: int,
        total_started_at: float,
        include_cluster_stage: bool,
        cpp_mode: bool = False,
    ) -> None:
        self.total = total
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(include_cluster_stage, cpp_mode)
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
        self._last_render_at = float("-inf")
        self._render(force=True)
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
            self._render_locked(force=True)

    def close(self) -> None:
        """Stop background refresh and close any progress backend."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=0.2)
        with self._lock:
            self._render_locked(force=True)
        if self._progress is not None:
            self._progress.close()

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render(force=True)

    def _render(self, *, force: bool = False) -> None:
        with self._lock:
            self._render_locked(force=force)

    def _render_locked(self, *, force: bool = False) -> None:
        now = perf_counter()
        if not force and now - self._last_render_at < PROGRESS_RENDER_INTERVAL_SECONDS:
            return
        self._last_render_at = now
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
PROGRESS_RENDER_INTERVAL_SECONDS = 0.10


class ParallelRunProgressDisplay:
    """Render aggregate and per-file progress for concurrent frame analysis."""

    def __init__(
        self,
        total: int,
        workers: int,
        total_started_at: float,
        include_cluster_stage: bool,
        cpp_mode: bool = False,
    ) -> None:
        self.total = total
        self.workers = workers
        self.total_started_at = total_started_at
        self.stage_groups = configured_stage_groups(include_cluster_stage, cpp_mode)
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
        self._last_render_at = float("-inf")
        self._render(force=True)
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
            self._render_locked(force=True)
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
            self._render_locked(force=True)

    def close(self) -> None:
        """Stop background refresh and close the progress backend."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=0.2)
        with self._lock:
            self._render_locked(force=True)
        if self._progress is not None:
            self._progress.close()

    def _tick(self) -> None:
        while not self._stop_event.wait(1.0):
            self._render(force=True)

    def _render(self, *, force: bool = False) -> None:
        with self._lock:
            self._render_locked(force=force)

    def _render_locked(self, *, force: bool = False) -> None:
        now = perf_counter()
        if not force and now - self._last_render_at < PROGRESS_RENDER_INTERVAL_SECONDS:
            return
        self._last_render_at = now
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
    progress = RunProgressDisplay(
        total=int(total_frames if total_frames is not None else len(paths)),
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
        cpp_mode=is_cpp_mode(config.get("mode")),
    )
    try:
        standalone_paths = (
            topology is None
            and bool(paths)
            and all(path.suffix.lower() in PARALLEL_SUFFIXES for path in paths)
        )
        if standalone_paths:
            for frame_index, path in enumerate(paths):
                callback = progress.start_frame(frame_index, path.stem)
                try:
                    frame = next(
                        iter(
                            read_frames(
                                [path],
                                xyz_scale=float(
                                    config["input"].get("xyz_scale", 0.1)
                                ),
                            )
                        )
                    )
                except Exception as exc:
                    if strict:
                        raise
                    row = failed_row(path.stem, str(path), str(exc))
                else:
                    row = process_frame(
                        frame_index,
                        frame,
                        config,
                        outdir,
                        strict=strict,
                        stage_callback=callback,
                    )
                rows.append(row)
                progress.complete_frame(row.get("status") == "ok")
            return rows

        frames = iter(
            read_frames(
                paths,
                topology=topology,
                xyz_scale=float(config["input"].get("xyz_scale", 0.1)),
                trajectory_stride=int(config["input"].get("trajectory_stride", 1)),
                lammps_config=config["input"].get("lammps", {}),
            )
        )
        frame_index = 0
        while True:
            try:
                frame = next(frames)
            except StopIteration:
                break
            except Exception as exc:
                if strict:
                    raise
                source = paths[0] if paths else Path("")
                frame_name = f"{source.stem}_frame{frame_index:06d}"
                progress.start_frame(frame_index, frame_name)
                row = failed_row(frame_name, str(source), str(exc))
                rows.append(row)
                progress.complete_frame(False)
                break
            callback = progress.start_frame(frame_index, frame.name)
            row = process_frame(
                frame_index,
                frame,
                config,
                outdir,
                strict=strict,
                stage_callback=callback,
            )
            rows.append(row)
            progress.complete_frame(row.get("status") == "ok")
            frame_index += 1
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
        cpp_mode=is_cpp_mode(config.get("mode")),
    )
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            task_iterator = iter(enumerate(paths))
            futures: dict[Any, int] = {}
            max_in_flight = process_in_flight_limit(workers)

            def fill_queue() -> None:
                while len(futures) < max_in_flight:
                    try:
                        frame_index, path = next(task_iterator)
                    except StopIteration:
                        return
                    future = executor.submit(
                        process_single_file_path,
                        frame_index,
                        path,
                        config,
                        outdir,
                        strict,
                        progress,
                    )
                    futures[future] = frame_index

            fill_queue()
            while futures:
                done, _ = wait(set(futures), return_when=FIRST_COMPLETED)
                for future in done:
                    expected_index = futures.pop(future)
                    try:
                        frame_index, row = future.result()
                    except Exception:
                        progress.complete_file(expected_index, False)
                        for sibling in futures:
                            sibling.cancel()
                        raise
                    rows_by_index[frame_index] = row
                    progress.complete_file(frame_index, row.get("status") == "ok")
                fill_queue()
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
        cpp_mode=is_cpp_mode(config.get("mode")),
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
        raise ValueError("Trajectory process analysis requires --top.")
    rows_by_index: dict[int, dict[str, Any]] = {}
    progress = ParallelRunProgressDisplay(
        total=len(raw_frame_indexes),
        workers=workers,
        total_started_at=total_started_at,
        include_cluster_stage=bool(config.get("hydrate_cluster", {}).get("enabled", False)),
        cpp_mode=is_cpp_mode(config.get("mode")),
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
                initargs=(
                    config,
                    str(outdir),
                    strict,
                    stage_queue,
                    str(trajectory),
                    str(topology),
                    config["input"].get("lammps", {}),
                ),
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
                            progress.complete_file(
                                frame_index,
                                row.get("status") == "ok",
                            )
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
        frame = next(
            iter(
                read_frames(
                    [path],
                    xyz_scale=float(config["input"].get("xyz_scale", 0.1)),
                )
            )
        )
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
    *,
    separated_output: bool = False,
) -> dict[str, Any]:
    """Analyze one frame, write per-frame files, and return a summary row."""
    if frame.time_ps is None:
        frame.time_ps = config["input"]["first_file_time_ps"] + frame_index * config["input"]["frame_time_step_ps"]
    try:
        result = analyze_frame(
            frame,
            config,
            stage_callback=stage_callback,
            normalize_config=False,
        )
        if separated_output:
            report_dir = outdir / "info"
            frame_dir = outdir / "gro" / frame.name
        else:
            report_dir = outdir / frame.name
            frame_dir = report_dir
        report_dir.mkdir(parents=True, exist_ok=True)
        frame_dir.mkdir(parents=True, exist_ok=True)
        report_stage(stage_callback, "writing outputs")
        write_frame_outputs(result, frame_dir, config, report_dir=report_dir)
        if output_enabled(config, "sqq-cage-gro"):
            write_sqq_cage_fragment(
                result,
                outdir / FRAGMENT_DIRECTORY,
                frame_index,
                requested_graph_mode=config["graph"]["bond_mode"],
            )
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


def write_frame_outputs(
    result: FrameResult,
    frame_dir: Path,
    config: dict[str, Any],
    *,
    report_dir: Path | None = None,
) -> None:
    """Write all configured per-frame output files."""
    output = config.get("output", {})
    report_dir = report_dir or frame_dir
    order_parameters = config.get("order", {}).get("parameters", ["f3", "f4"])
    cpp_mode = is_cpp_mode(config.get("mode"))
    if output_enabled(config, "info"):
        write_frame_info(
            result,
            report_dir,
            ring_sizes=list(result.ring_report_sizes),
            requested_bond_mode=config["graph"]["bond_mode"],
            order_parameters=order_parameters,
            analysis_mode=config.get("mode", DEFAULT_MODE),
            input_metadata=frame_input_metadata(config),
        )
    else:
        remove_optional_info_output(result, report_dir)

    if not cpp_mode and output_enabled(config, "membership-tsv"):
        write_membership(result, report_dir)
    else:
        remove_optional_tsv_outputs(
            result,
            report_dir,
            remove_membership=True,
            remove_order=False,
        )
    if not cpp_mode and output_enabled(config, "order-tsv"):
        write_order_parameter(result, report_dir, order_parameters=order_parameters)
    else:
        remove_optional_tsv_outputs(
            result,
            report_dir,
            remove_membership=False,
            remove_order=True,
        )

    if not cpp_mode and output_enabled(config, "vmd"):
        write_vmd_script(result, report_dir)
    else:
        remove_optional_vmd_output(result, report_dir)

    layout = str(output.get("structure_layout", "grouped"))
    write_empty = bool(output.get("write_empty_files", False))
    remove_generated_gro_outputs(result, frame_dir, "ring-gro", layout)
    if not cpp_mode and output_enabled(config, "ring-gro"):
        write_ring_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            sizes=set(result.ring_report_sizes),
        )
    remove_generated_gro_outputs(result, frame_dir, "half-gro", layout)
    if not cpp_mode and output_enabled(config, "half-gro"):
        write_half_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
        )
    remove_generated_gro_outputs(result, frame_dir, "quasi-gro", layout)
    if not cpp_mode and output_enabled(config, "quasi-gro"):
        write_quasi_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
        )
    remove_generated_gro_outputs(result, frame_dir, "cage-gro", layout)
    if output_enabled(config, "cage-gro"):
        write_cage_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
            include_centers=not cpp_mode,
        )
    remove_generated_gro_outputs(result, frame_dir, "cluster-gro", layout)
    if not cpp_mode and result.hydrate_cluster_enabled and output_enabled(config, "cluster-gro"):
        write_hydrate_cluster_gro_files(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
        )
    remove_generated_gro_outputs(result, frame_dir, "ice-gro", layout)
    if not cpp_mode and output_enabled(config, "ice-gro"):
        write_ice_gro_file(
            result,
            frame_dir,
            write_empty=write_empty,
            layout=layout,
        )
    remove_frame_directory_if_empty(frame_dir)
    if report_dir != frame_dir:
        remove_frame_directory_if_empty(report_dir)
        remove_frame_directory_if_empty(report_dir.parent)
        remove_frame_directory_if_empty(frame_dir.parent)


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


def remove_generated_gro_outputs(
    result: FrameResult,
    frame_dir: Path,
    output_type: str,
    layout: str,
) -> None:
    """Remove known SQQ-generated GRO files before rewriting one category."""
    grouped_directories = {
        "ring-gro": "ring",
        "half-gro": "half_cage",
        "quasi-gro": "quasi_cage",
        "cage-gro": "cage",
        "ice-gro": "ice",
        "cluster-gro": "hydrate_cluster",
    }
    flat_patterns = {
        "ring-gro": f"{result.frame.name}_ring_*.gro",
        "half-gro": f"{result.frame.name}_hc_*.gro",
        "quasi-gro": f"{result.frame.name}_qc_*.gro",
        "cage-gro": f"{result.frame.name}_cage_*.gro",
        "ice-gro": f"{result.frame.name}_ice*.gro",
        "cluster-gro": f"{result.frame.name}_cluster_*.gro",
    }
    if layout not in {"grouped", "flat"}:
        raise ValueError("output.structure_layout must be 'grouped' or 'flat'.")
    root = frame_dir / grouped_directories[output_type]
    if root.exists():
        generated_parents: set[Path] = set()
        for path in root.rglob("*.gro"):
            if path.name.startswith(f"{result.frame.name}_"):
                parent = path.parent
                while parent == root or root in parent.parents:
                    generated_parents.add(parent)
                    if parent == root:
                        break
                    parent = parent.parent
                path.unlink(missing_ok=True)
        for directory in sorted(
            generated_parents,
            key=lambda path: len(path.parts),
            reverse=True,
        ):
            try:
                directory.rmdir()
            except OSError:
                pass
    for path in frame_dir.glob(flat_patterns[output_type]):
        path.unlink(missing_ok=True)


def remove_frame_directory_if_empty(frame_dir: Path) -> None:
    """Remove only the frame directory itself, preserving unrelated subdirectories."""
    try:
        frame_dir.rmdir()
    except OSError:
        pass


def worker_policy_text(config: dict[str, Any]) -> str:
    """Describe the active automatic or explicit worker policy."""
    value = config.get("parallel", {}).get("workers", "auto")
    reserve_text = "reserve 1 physical core"
    if is_auto_worker_request(value):
        mode = config.get("mode", DEFAULT_MODE)
        fixed_count = mode_worker_count(mode)
        if fixed_count is not None:
            unit = "worker" if fixed_count == 1 else "workers"
            return f"mode default ({fixed_count} {unit})"
        percent = int(round(mode_worker_fraction(mode) * 100))
        return f"auto ({percent}% of physical cores, {reserve_text})"
    try:
        request_text = describe_worker_request(value)
    except ValueError:
        request_text = str(value)
    return f"explicit ({request_text}, {reserve_text})"


def describe_worker_request(value: Any) -> str:
    """Render a user-facing worker request for terminal and workbook metadata."""
    kind, amount = classify_worker_request(value)
    if kind == "auto":
        return "auto"
    if kind == "fraction":
        return f"{format_worker_percent(float(amount))} of physical cores"
    if kind == "count":
        return f"{int(amount)} workers"
    raise worker_value_error()


def is_auto_worker_request(value: Any) -> bool:
    """Return true for unset or textual auto worker requests."""
    return value is None or str(value).strip().lower() in {"", "auto"}


def classify_worker_request(value: Any) -> tuple[str, float | int | None]:
    """Classify -w/--worker by input form rather than numeric value alone."""
    text = str(value).strip().lower()
    if text in {"", "auto"}:
        return "auto", None
    if text.endswith("%"):
        percent = parse_worker_number(text[:-1], value)
        if percent <= 0 or percent > 100:
            raise worker_value_error()
        return "fraction", percent / 100.0
    if "." in text:
        fraction = parse_worker_number(text, value)
        if fraction <= 0 or fraction > 1:
            raise worker_value_error()
        return "fraction", fraction
    try:
        count = int(text)
    except (TypeError, ValueError) as exc:
        raise worker_value_error() from exc
    if count < 1:
        raise worker_value_error()
    return "count", count


def format_worker_percent(fraction: float) -> str:
    """Format a worker CPU fraction as a compact percentage string."""
    percent = fraction * 100.0
    return f"{percent:g}%"


def parse_worker_number(text: str, original: Any) -> float:
    """Parse a finite positive worker option number."""
    try:
        number = float(text)
    except (TypeError, ValueError) as exc:
        raise worker_value_error() from exc
    if not math.isfinite(number) or number <= 0:
        raise worker_value_error()
    return number


def worker_value_error() -> ValueError:
    """Build the standard worker-option validation error."""
    return ValueError("parallel.workers / --worker must be 'auto', a positive integer worker count such as 1 or 4, a decimal CPU fraction in (0, 1] such as 0.5 or 1.0, or a percentage in (0%, 100%].")


def normalize_parallel_backend(value: Any) -> str:
    """Normalize the supported serial, process, and compatibility thread backends."""
    backend = str(value or "process").strip().lower()
    if backend not in {"process", "thread", "serial"}:
        raise ValueError("parallel.backend / --parallel-backend must be process, thread, or serial.")
    return backend


def resolve_workers(
    value: Any,
    n_paths: int,
    mode: Any = DEFAULT_MODE,
    cpu_total: int | None = None,
    backend: str = "process",
) -> int:
    """Resolve workers from physical cores, task count, and platform caps."""
    resolved_backend = normalize_parallel_backend(backend)
    single_or_serial = max(0, int(n_paths)) <= 1 or resolved_backend == "serial"
    if is_auto_worker_request(value):
        if single_or_serial:
            return 1
        physical_total = max(1, int(cpu_total if cpu_total is not None else physical_cpu_count()))
        fixed_count = mode_worker_count(mode)
        requested = (
            fixed_count
            if fixed_count is not None
            else max(1, int(physical_total * mode_worker_fraction(mode)))
        )
    else:
        # Validate before the single-task short-circuit.
        classify_worker_request(value)
        if single_or_serial:
            return 1
        physical_total = max(1, int(cpu_total if cpu_total is not None else physical_cpu_count()))
        requested = resolve_explicit_worker_request(value, physical_total)
    usable_workers = max(1, physical_total - 1)
    requested = min(requested, usable_workers, max(1, n_paths))
    if resolved_backend == "process":
        platform_cap = process_worker_cap()
        if platform_cap is not None:
            requested = min(requested, platform_cap)
    return max(1, requested)


def resolve_explicit_worker_request(value: Any, physical_total: int) -> int:
    """Resolve -w/--worker using form-based fraction/count semantics."""
    kind, amount = classify_worker_request(value)
    if kind == "fraction":
        return max(1, int(physical_total * float(amount)))
    if kind == "count":
        return int(amount)
    raise worker_value_error()


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
        and paths[0].suffix.lower() in ({".xtc", ".trr"} | set(LAMMPS_TRAJECTORY_SUFFIXES))
    )

def apply_cli_overrides(config: dict[str, Any], args: Namespace) -> None:
    """Apply command-line options after YAML/default configuration."""
    if args.pattern:
        config["input"]["pattern"] = args.pattern
    if args.recursive:
        config["input"]["recursive"] = True
    if getattr(args, "xyz_scale", None) is not None:
        if not math.isfinite(args.xyz_scale) or args.xyz_scale <= 0:
            raise ValueError("--xyz-scale must be positive and finite.")
        config["input"]["xyz_scale"] = args.xyz_scale
    if getattr(args, "trajectory_stride", None) is not None:
        config["input"]["trajectory_stride"] = args.trajectory_stride
    lammps = config["input"].setdefault("lammps", {})
    if getattr(args, "lammps_units", None) is not None:
        lammps["units"] = args.lammps_units
    if getattr(args, "lammps_timestep", None) is not None:
        lammps["timestep"] = args.lammps_timestep
    if getattr(args, "lammps_atom_style", None) is not None:
        lammps["atom_style"] = args.lammps_atom_style
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
    unified_order = getattr(args, "order_parameter", None)
    legacy_order_options = any(
        (
            getattr(args, "no_q", False),
            getattr(args, "q_degree", None),
            getattr(args, "mcg3", None),
            getattr(args, "dhop30", None),
        )
    )
    if unified_order is not None:
        config["order"]["parameters"] = list(normalize_order_parameters(unified_order))
        if legacy_order_options:
            warn_legacy_order_cli(ignored=True)
    elif legacy_order_options:
        warn_legacy_order_cli(ignored=False)
        selected = set(
            normalize_order_parameters(
                config.get("order", {}).get("parameters", ["f3", "f4"])
            )
        )
        if getattr(args, "no_q", False):
            selected = {name for name in selected if not name.startswith("q")}
        elif getattr(args, "q_degree", None):
            selected = {name for name in selected if not name.startswith("q")}
            selected.update(f"q{degree}" for degree in normalize_q_degree(args.q_degree))
        if getattr(args, "mcg3", None):
            if parse_on_off(args.mcg3, "--mcg3"):
                selected.add("mcg3")
            else:
                selected.discard("mcg3")
        if getattr(args, "dhop30", None):
            if parse_on_off(args.dhop30, "--dhop30"):
                selected.add("dhop30")
            else:
                selected.discard("dhop30")
        config["order"]["parameters"] = list(
            normalize_order_parameters(selected)
        )
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
    if getattr(args, "find_cluster", None):
        config["hydrate_cluster"]["enabled"] = parse_on_off(
            args.find_cluster,
            "--find-cluster",
        )
    if getattr(args, "cluster_min_cage", None) is not None:
        if args.cluster_min_cage < 1:
            raise ValueError("--cluster-min-cage must be at least 1.")
        config["hydrate_cluster"]["min_cage"] = args.cluster_min_cage
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
    if getattr(args, "output_type", None) is not None:
        config["output"]["types"] = args.output_type
        if not config.get("hydrate_cluster", {}).get("enabled", False):
            requested = {
                item.strip().lower()
                for item in str(args.output_type).split(",")
            }
            if (
                "all" not in requested
                and {"cluster-gro", "cluster-detail"}.intersection(requested)
            ):
                raise ValueError("Cluster outputs require --find-cluster on.")
    if getattr(args, "cage_isomer_rows", None):
        config["output"]["cage_isomer_rows"] = args.cage_isomer_rows


def warn_legacy_order_cli(*, ignored: bool) -> None:
    """Print one visible compatibility warning for pre-0.2.7 selector options."""
    suffix = " They were ignored because --order-parameter was also supplied." if ignored else ""
    print(
        "Warning: --no-q, -q/--q-degree, --mcg3, and --dhop30 are deprecated; "
        f"use --order-parameter instead.{suffix}",
        file=sys.stderr,
    )


def normalize_analysis_scopes(config: dict[str, Any]) -> None:
    """Normalize search and report scopes before frames are analyzed."""
    input_config = config.setdefault("input", {})
    input_config["recursive"] = parse_on_off(input_config.get("recursive", False), "input.recursive")
    input_config["xyz_scale"] = finite_float(
        input_config.get("xyz_scale", 0.1),
        "input.xyz_scale / --xyz-scale",
        positive=True,
    )
    input_config["first_file_time_ps"] = finite_float(
        input_config.get("first_file_time_ps", 0.0),
        "input.first_file_time_ps",
    )
    input_config["frame_time_step_ps"] = finite_float(
        input_config.get("frame_time_step_ps", 100.0),
        "input.frame_time_step_ps",
    )
    legacy_stride = input_config.pop("xtc_stride", None)
    stride_value = input_config.get("trajectory_stride", legacy_stride if legacy_stride is not None else 1)
    input_config["trajectory_stride"] = positive_integer(
        stride_value,
        "input.trajectory_stride / --trajectory-stride",
    )
    raw_lammps = input_config.get("lammps", {})
    if not isinstance(raw_lammps, dict):
        raise ValueError("input.lammps must be a mapping.")
    lammps_values = dict(raw_lammps)
    lammps_values["stride"] = input_config["trajectory_stride"]
    settings = normalize_lammps_config(lammps_values)
    type_map: dict[str, dict[str, Any]] = {}
    for type_id, entry in settings.type_map.items():
        if entry.ignore:
            type_map[type_id] = {"ignore": True}
        else:
            type_map[type_id] = {"resname": entry.resname, "atomname": entry.atomname}
    input_config["lammps"] = {
        "units": settings.units,
        "timestep": settings.timestep,
        "atom_style": str(raw_lammps.get("atom_style", "full")).strip().lower(),
        "coordinate_convention": settings.coordinate_convention,
        "type_map": type_map,
    }

    graph = config.setdefault("graph", {})
    # Effective mode is internal run state, never a reusable user setting.
    graph.pop("effective_bond_mode", None)
    graph_mode = str(graph.get("bond_mode", "auto")).strip().lower()
    if graph_mode not in {"auto", "hbond", "oo", "pairs"}:
        raise ValueError("graph.bond_mode must be auto, hbond, oo, or pairs.")
    graph["bond_mode"] = graph_mode
    graph["oo_cutoff_nm"] = finite_float(graph.get("oo_cutoff_nm", 0.35), "graph.oo_cutoff_nm", positive=True)
    graph["hbond_distance_nm"] = finite_float(graph.get("hbond_distance_nm", 0.35), "graph.hbond_distance_nm", positive=True)
    graph["hbond_angle_deg"] = finite_float(graph.get("hbond_angle_deg", 30.0), "graph.hbond_angle_deg")
    if not 0 <= graph["hbond_angle_deg"] <= 180:
        raise ValueError("graph.hbond_angle_deg must be between 0 and 180 degrees.")
    pair_id = str(graph.get("pair_id", "resid")).strip().lower()
    if pair_id not in {"resid", "oxygen_index", "atomid"}:
        raise ValueError("graph.pair_id must be resid, oxygen_index, or atomid.")
    graph["pair_id"] = pair_id

    pbc = config.setdefault("pbc", {})
    box_mode = str(pbc.get("box_mode", "orthorhombic")).strip().lower()
    if box_mode != "orthorhombic":
        raise ValueError("pbc.box_mode must be orthorhombic.")
    pbc["box_mode"] = box_mode

    water = config.setdefault("water", {})
    water["resnames"] = string_list(water.get("resnames", []), "water.resnames")
    water["oxygen_names"] = string_list(water.get("oxygen_names", []), "water.oxygen_names", allow_empty=True)
    water["hydrogen_names"] = string_list(water.get("hydrogen_names", []), "water.hydrogen_names", allow_empty=True)
    guest = config.setdefault("guest", {})
    guest["resnames"] = string_list(guest.get("resnames", []), "guest.resnames", allow_empty=True)
    center_atoms = guest.get("center_atoms", {})
    if not isinstance(center_atoms, dict):
        raise ValueError("guest.center_atoms must be a residue-to-atom-list mapping.")
    guest["center_atoms"] = {
        str(resname).strip(): string_list(atom_names, f"guest.center_atoms.{resname}", allow_empty=True)
        for resname, atom_names in center_atoms.items()
        if str(resname).strip()
    }

    ring = config.setdefault("ring", {})
    ring["chordless"] = parse_on_off(ring.get("chordless", True), "ring.chordless")
    search_sizes = resolve_size_list(ring.get("sizes", []), fallback=[], key="ring.sizes")
    unsupported = set(search_sizes) - {4, 5, 6, 7}
    if unsupported:
        raise ValueError(f"ring.sizes / --size supports only 4, 5, 6, and 7; got {sorted(unsupported)}")
    ring_report_sizes = resolve_size_list(
        ring.get("report_sizes", "auto"),
        fallback=search_sizes,
        key="ring.report_sizes",
    )
    if not set(ring_report_sizes) <= set(search_sizes):
        raise ValueError("ring.report_sizes / --ring-size must be a subset of ring.sizes / --size.")

    ring_definition = str(ring.get("definition", "chordless")).strip().lower()
    if ring_definition not in {"chordless", "shortest_path"}:
        raise ValueError("ring.definition / --ring-definition must be chordless or shortest_path.")
    ring["definition"] = ring_definition

    quasi = config.setdefault("quasi_cage", {})
    quasi["enabled"] = parse_on_off(quasi.get("enabled", True), "quasi_cage.enabled")
    raw_quasi_base_sizes = quasi.get("base_sizes", "auto")
    raw_quasi_side_sizes = quasi.get("side_sizes", "auto")
    quasi_base_sizes = resolve_size_list(raw_quasi_base_sizes, fallback=search_sizes, key="quasi_cage.base_sizes")
    quasi_side_sizes = resolve_size_list(raw_quasi_side_sizes, fallback=search_sizes, key="quasi_cage.side_sizes")
    if not set(quasi_base_sizes) <= set(search_sizes):
        raise ValueError("quasi_cage.base_sizes must be a subset of ring.sizes / --size.")
    if not set(quasi_side_sizes) <= set(search_sizes):
        raise ValueError("quasi_cage.side_sizes must be a subset of ring.sizes / --size.")
    quasi["base_sizes"] = "auto" if str(raw_quasi_base_sizes).strip().lower() == "auto" else quasi_base_sizes
    quasi["side_sizes"] = "auto" if str(raw_quasi_side_sizes).strip().lower() == "auto" else quasi_side_sizes
    quasi_policy = str(quasi.get("search_policy", "bounded")).strip().lower()
    if quasi_policy not in {"bounded", "exact"}:
        raise ValueError("quasi_cage.search_policy / --quasi-search-policy must be bounded or exact.")
    quasi["search_policy"] = quasi_policy
    for key, default in (
        ("max_combinations_per_base", 50000),
        ("max_layers", 1),
        ("max_rings_per_layer", 6),
        ("max_layer_states_per_seed", 200),
        ("max_candidates_per_edge", 4),
        ("max_layer_candidates", 24),
    ):
        quasi[key] = positive_integer(quasi.get(key, default), f"quasi_cage.{key}")

    cage = config.setdefault("cage", {})
    max_faces = positive_integer(cage.get("max_faces", 20), "cage.max_faces / --max-cage-face")
    report_types = resolve_cage_report_types(
        cage.get("report_types", []),
        search_sizes,
        max_faces,
    )
    ring["sizes"] = search_sizes
    ring["report_sizes"] = ring_report_sizes
    cage["report_types"] = "all" if report_types is None else list(report_types)
    cage["max_faces"] = max_faces
    cage["enabled"] = parse_on_off(cage.get("enabled", True), "cage.enabled")
    search_mode = str(cage.get("search_mode", "grow")).strip().lower()
    if search_mode in {"expand", "hybrid"}:
        search_mode = "grow"
    if search_mode not in {"grow", "pair", "patch_pair"}:
        raise ValueError("cage.search_mode must be grow, pair, or patch_pair.")
    cage["search_mode"] = search_mode
    seed_mode = str(cage.get("seed_mode", "ring")).strip().lower()
    if seed_mode not in {"ring", "patch"}:
        raise ValueError("cage.seed_mode must be ring or patch.")
    cage["seed_mode"] = seed_mode
    occupancy_mode = str(cage.get("occupancy_mode", "polyhedron")).strip().lower()
    if occupancy_mode not in {"polyhedron", "center", "auto"}:
        raise ValueError("cage.occupancy_mode must be polyhedron, center, or auto.")
    cage["occupancy_mode"] = occupancy_mode
    cage["occupancy_radius_nm"] = finite_float(cage.get("occupancy_radius_nm", 0.5), "cage.occupancy_radius_nm", positive=True)
    for key, default in (
        ("max_states_per_seed", 20000),
        ("max_total_states", 5000000),
        ("max_boundary_candidates", 8),
    ):
        cage[key] = positive_integer(cage.get(key, default), f"cage.{key}")
    cage["fast_closure"] = parse_on_off(cage.get("fast_closure", True), "cage.fast_closure")
    cage["fast_closure_max_states"] = positive_integer(cage.get("fast_closure_max_states", 20000), "cage.fast_closure_max_states")
    cage["scientific_validation"] = parse_on_off(
        cage.get("scientific_validation", False),
        "cage.scientific_validation",
    )
    cage["max_face_planarity_rms_nm"] = finite_float(cage.get("max_face_planarity_rms_nm", 0.06), "cage.max_face_planarity_rms_nm", nonnegative=True)
    cage["max_face_edge_cv"] = finite_float(cage.get("max_face_edge_cv", 0.35), "cage.max_face_edge_cv", nonnegative=True)
    cage["min_cage_volume_nm3"] = finite_float(cage.get("min_cage_volume_nm3", 1.0e-6), "cage.min_cage_volume_nm3", positive=True)
    hydrate_cluster = config.setdefault("hydrate_cluster", {})
    if "detail" in hydrate_cluster:
        raise ValueError(
            "hydrate_cluster.detail is no longer supported; "
            "add cluster-detail to output.types."
        )
    hydrate_cluster["enabled"] = parse_on_off(hydrate_cluster.get("enabled", False), "hydrate_cluster.enabled")
    hydrate_cluster["min_cage"] = positive_integer(hydrate_cluster.get("min_cage", 2), "hydrate_cluster.min_cage / --cluster-min-cage")
    parallel = config.setdefault("parallel", {})
    parallel["backend"] = normalize_parallel_backend(parallel.get("backend", "process"))
    parallel["math_threads"] = positive_integer(parallel.get("math_threads", 1), "parallel.math_threads")
    output = config.setdefault("output", {})
    removed_output_keys = {
        "disabled_outputs",
        "write_tsv",
        "write_order_tsv",
        "write_vmd",
        "write_info",
        "write_gro",
        "write_ring_gro",
        "write_half_cage_gro",
        "write_quasi_cage_gro",
        "write_cage_gro",
        "write_ice_gro",
        "write_xlsx_summary",
        "write_summary_detail_csv",
    }.intersection(output)
    if removed_output_keys:
        names = ", ".join(sorted(removed_output_keys))
        raise ValueError(
            f"Unsupported output configuration key(s): {names}. Use output.types."
        )
    raw_output_types = output.get("types")
    output_normalizer = (
        normalize_cpp_output_types
        if is_cpp_mode(config.get("mode", DEFAULT_MODE))
        else normalize_output_types
    )
    output_types = list(output_normalizer(raw_output_types))
    if not hydrate_cluster["enabled"]:
        output_types = [
            output_type
            for output_type in output_types
            if output_type not in {"cluster-gro", "cluster-detail"}
        ]
    output["types"] = output_types
    output["write_empty_files"] = parse_on_off(
        output.get("write_empty_files", False),
        "output.write_empty_files",
    )
    for key, default in (
        ("summary_csv_dir", "summary"),
        ("summary_detail_dir", "summary_detail"),
    ):
        directory = str(output.get(key, default)).strip() or default
        directory_path = Path(directory)
        if directory_path.is_absolute() or ".." in directory_path.parts:
            raise ValueError(f"output.{key} must be a relative directory inside the output folder.")
        output[key] = directory
    if (
        Path(output["summary_csv_dir"]).as_posix().casefold()
        == Path(output["summary_detail_dir"]).as_posix().casefold()
    ):
        raise ValueError(
            "output.summary_csv_dir and output.summary_detail_dir must be different."
        )
    cage_isomer_rows = str(output.get("cage_isomer_rows", "nonzero")).strip().lower()
    if cage_isomer_rows not in {"nonzero", "all"}:
        raise ValueError("output.cage_isomer_rows / --cage-isomer-rows must be nonzero or all.")
    output["cage_isomer_rows"] = cage_isomer_rows
    structure_layout = str(output.get("structure_layout", "grouped")).strip().lower()
    if structure_layout not in {"grouped", "flat"}:
        raise ValueError("output.structure_layout must be grouped or flat.")
    output["structure_layout"] = structure_layout
    guest_center_mode = str(config.get("guest", {}).get("center_mode", "center_atom")).strip().lower()
    if guest_center_mode not in {"center_atom", "centroid", "auto"}:
        raise ValueError("guest.center_mode must be center_atom, centroid, or auto.")
    config["guest"]["center_mode"] = guest_center_mode
    order = config.setdefault("order", {})
    order["parameters"] = list(
        normalize_order_parameters(order.get("parameters", ["f3", "f4"]))
    )
    q_neighbor_mode = normalize_q_neighbor_mode(str(order.get("q_neighbor_mode", "graph")))
    order["q_neighbor_mode"] = q_neighbor_mode
    order["q_cutoff_nm"] = finite_float(order.get("q_cutoff_nm", 0.35), "order.q_cutoff_nm", positive=True)
    order["q_n_neighbor"] = resolve_q_neighbor_count(q_neighbor_mode, order.get("q_n_neighbor"))
    focus_value = order.get("focus_waters", [])
    if isinstance(focus_value, str):
        focus_items = [item.strip() for item in focus_value.split(",") if item.strip()]
    else:
        try:
            focus_items = list(focus_value)
        except TypeError as exc:
            raise ValueError("order.focus_waters must be a list of residue ids.") from exc
    try:
        order["focus_waters"] = sorted({int(item) for item in focus_items})
    except (TypeError, ValueError) as exc:
        raise ValueError("order.focus_waters must contain integer residue ids.") from exc
    per_water_order = any(
        name in {"f3", "f4"} or name.startswith("q")
        for name in order["parameters"]
    )
    if output_enabled(config, "order-tsv") and not per_water_order:
        print(
            "Warning: output type 'order-tsv' has no per-water F3/F4/Q_l selection; "
            "no order-parameter TSV will be written.",
            file=sys.stderr,
        )
    hydrate_order = config.setdefault("hydrate_order", {})
    positive_values = (
        ("mcg_guest_cutoff_nm", 0.90),
        ("mcg_water_cutoff_nm", 0.60),
        ("dhop_neighbor_cutoff_nm", 0.35),
    )
    for key, default in positive_values:
        hydrate_order[key] = finite_float(hydrate_order.get(key, default), f"hydrate_order.{key}", positive=True)
    cone_angle = finite_float(hydrate_order.get("mcg_cone_half_angle_deg", 45.0), "hydrate_order.mcg_cone_half_angle_deg")
    if not 0 < cone_angle < 90:
        raise ValueError("hydrate_order.mcg_cone_half_angle_deg must be between 0 and 90.")
    hydrate_order["mcg_cone_half_angle_deg"] = cone_angle
    for key, default in (("mcg_min_waters", 5), ("dhop_min_qualified_neighbors", 3)):
        hydrate_order[key] = positive_integer(hydrate_order.get(key, default), f"hydrate_order.{key}")
    planar_counts: set[int] = set()
    try:
        raw_planar_value = hydrate_order.get("dhop_planar_counts", [11, 12])
        raw_planar_counts = (
            raw_planar_value.split(",") if isinstance(raw_planar_value, str) else list(raw_planar_value)
        )
    except TypeError as exc:
        raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.") from exc
    for raw_value in raw_planar_counts:
        try:
            numeric = float(raw_value)
            count = int(numeric)
        except (TypeError, ValueError) as exc:
            raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.") from exc
        if not math.isfinite(numeric) or numeric != count or count < 0:
            raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.")
        planar_counts.add(count)
    if not planar_counts:
        raise ValueError("hydrate_order.dhop_planar_counts must contain non-negative integers.")
    hydrate_order["dhop_planar_counts"] = sorted(planar_counts)
    hydrate_order["mcg_guest_resnames"] = string_list(hydrate_order.get("mcg_guest_resnames", ["CH4", "MET"]), "hydrate_order.mcg_guest_resnames")

    ice = config.setdefault("ice", {})
    ice["enabled"] = parse_on_off(ice.get("enabled", True), "ice.enabled")
    ice_method = str(ice.get("method", "chill")).strip().lower()
    if ice_method != "chill":
        raise ValueError("ice.method must be chill.")
    ice["method"] = ice_method
    ice["min_six_rings"] = positive_integer(ice.get("min_six_rings", 2), "ice.min_six_rings")
    ice["require_four_coord_neighbors"] = parse_on_off(ice.get("require_four_coord_neighbors", True), "ice.require_four_coord_neighbors")


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


def finite_float(
    value: Any,
    key: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> float:
    """Normalize one finite floating-point configuration value."""
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a finite number.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{key} must be finite.")
    if positive and number <= 0:
        raise ValueError(f"{key} must be positive.")
    if nonnegative and number < 0:
        raise ValueError(f"{key} must be non-negative.")
    return number


def positive_integer(value: Any, key: str) -> int:
    """Normalize one strictly positive integer configuration value."""
    if isinstance(value, bool):
        raise ValueError(f"{key} must be a positive integer.")
    try:
        numeric = float(value)
        number = int(numeric)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a positive integer.") from exc
    if not math.isfinite(numeric) or numeric != number or number < 1:
        raise ValueError(f"{key} must be a positive integer.")
    return number


def string_list(value: Any, key: str, *, allow_empty: bool = False) -> list[str]:
    """Normalize a comma-separated string or sequence of names."""
    if value is None:
        items = []
    elif isinstance(value, str):
        items = value.split(",")
    else:
        try:
            items = list(value)
        except TypeError as exc:
            raise ValueError(f"{key} must be a list or comma-separated string.") from exc
    names = [str(item).strip() for item in items if str(item).strip()]
    if not names and not allow_empty:
        raise ValueError(f"{key} must contain at least one name.")
    return names


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
    *,
    normalize_config: bool = True,
) -> FrameResult:
    """Analyze one frame and return all topology objects for export."""
    report_stage(stage_callback, "resolving settings")
    if normalize_config:
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
    if is_cpp_mode(config.get("mode")):
        report_stage(stage_callback, "building water graph")
        report_stage(stage_callback, "searching rings")
        report_stage(stage_callback, "searching cage")
        result = analyze_frame_cpp(
            frame,
            waters,
            guests,
            config,
            cage_report_types=cage_report_types,
            ring_report_sizes=tuple(ring_report_sizes),
        )
        report_stage(stage_callback, "computing order parameters")
        return result
    # All structure classifiers use this graph.
    report_stage(stage_callback, "building water graph")
    graph = build_water_graph(
        frame.atoms,
        waters,
        frame.box,
        bond_mode=config["graph"].get(
            "effective_bond_mode", config["graph"]["bond_mode"]
        ),
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
    hydrate_cluster_detail = output_enabled(config, "cluster-detail")
    if hydrate_cluster_enabled:
        report_stage(stage_callback, "classifying hydrate cluster")
        rings_by_id = ring_topology.ring_by_id
        ring_sizes_by_id = {ring_id: ring.size for ring_id, ring in rings_by_id.items()}
        hydrate_clusters, hydrate_motifs, hydrate_domains, isolated_cage_ids = analyze_hydrate_clusters(
            all_cages,
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
    order_parameters = normalize_order_parameters(
        config.get("order", {}).get("parameters", ["f3", "f4"])
    )
    selected_order_parameters = set(order_parameters)
    q_degrees = q_degrees_from_order_parameters(order_parameters)
    if selected_order_parameters & {"f3", "f4"} or q_degrees:
        f3f4 = compute_order_parameters(
            frame,
            waters,
            graph,
            f3_enabled="f3" in selected_order_parameters,
            f4_enabled="f4" in selected_order_parameters,
            q_enabled=bool(q_degrees),
            q_neighbor_mode=str(config["order"].get("q_neighbor_mode", "graph")),
            q_cutoff_nm=float(config["order"].get("q_cutoff_nm", 0.35)),
            q_n_neighbor=config["order"].get("q_n_neighbor", None),
            q_degree=q_degrees,
            focus_resids=focus_resids,
        )
    else:
        f3f4 = None

    hydrate_parameters = selected_order_parameters & {
        "mcg1",
        "mcg3",
        "dhop35",
        "dhop30",
    }
    if hydrate_parameters:
        hydrate_order_config = {
            **config.get("hydrate_order", {}),
            "mcg1_enabled": "mcg1" in hydrate_parameters,
            "mcg3_enabled": "mcg3" in hydrate_parameters,
            "dhop35_enabled": "dhop35" in hydrate_parameters,
            "dhop30_enabled": "dhop30" in hydrate_parameters,
        }
        mcg1, mcg3 = compute_mcg_order(frame, waters, guests, hydrate_order_config)
        dhop35, dhop30 = compute_dhop_order(frame, waters, hydrate_order_config)
        hydrate_order = HydrateOrderResult(
            mcg1=mcg1,
            dhop35=dhop35,
            mcg3=mcg3,
            dhop30=dhop30,
        )
    else:
        hydrate_order = None
    report_stage(stage_callback, "classifying ice")
    ice_classes = classify_ice_waters(
        graph,
        waters,
        rings,
        enabled=bool(config["ice"].get("enabled", False)),
        min_six_rings=int(config["ice"].get("min_six_rings", 2)),
        require_four_coord_neighbors=bool(config["ice"].get("require_four_coord_neighbors", True)),
    )
    if bool(config["quasi_cage"].get("enabled", False)) and not half_cages and not quasi_cages:
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
