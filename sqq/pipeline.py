from __future__ import annotations

"""Top-level analysis pipeline for the SQQ command line."""

from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Callable

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    def tqdm(iterable, total=None, desc=None):
        return iterable

from .banner import SQQ_BANNER
from .config import load_config
from .core.cage import find_cages
from .core.f3f4 import compute_f3f4
from .core.graph import build_water_graph
from .core.ice import classify_ice_waters
from .core.quasi_cage import find_cage_patches
from .core.ring import find_rings
from .core.selection import select_guests, select_waters
from .io.gro_writer import write_cage_gro_files, write_half_cage_gro_files, write_ice_gro_file, write_quasi_cage_gro_files, write_ring_gro_files
from .io.summary import failed_row, result_row, write_f3f4, write_frame_info, write_membership, write_summary, write_vmd_script
from .io.trajectory import expand_inputs, read_frames
from .models import Cage, CagePatch, Frame, FrameResult


PARALLEL_SUFFIXES = {".gro", ".xyz"}


def analyze(args: Namespace) -> None:
    """Run SQQ analysis from parsed command-line arguments."""
    started_at = perf_counter()
    config = load_config(Path(args.config) if args.config else None)
    apply_cli_overrides(config, args)

    # Directory input follows a one-file-per-frame workflow.
    input_path = Path(args.input)
    pattern = args.pattern or config["input"]["pattern"]
    recursive = bool(args.recursive or config["input"]["recursive"])
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    n_jobs = resolve_n_jobs(config["parallel"].get("n_jobs"), len(paths))
    topology = Path(args.topology) if args.topology else None
    print_run_header(args, config, input_path, outdir, paths, topology, n_jobs)
    if n_jobs > 1 and can_parallelize_paths(paths, topology):
        rows = analyze_paths_parallel(paths, outdir, config, n_jobs=n_jobs, strict=bool(args.strict))
    else:
        rows = analyze_paths_serial(paths, outdir, config, topology=topology, strict=bool(args.strict))

    elapsed_seconds = perf_counter() - started_at
    run_info = build_run_info(args, config, input_path, outdir, paths, topology, n_jobs, elapsed_seconds)
    write_summary(rows, outdir, config, write_xlsx=config["output"]["write_xlsx_summary"], run_info=run_info)
    print(f"Wrote SQQ results: {outdir}")


def build_run_info(
    args: Namespace,
    config: dict[str, Any],
    input_path: Path,
    outdir: Path,
    paths: list[Path],
    topology: Path | None,
    n_jobs: int,
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Collect run-level metadata for the summary workbook."""
    info: dict[str, Any] = {
        "working_dir": str(Path.cwd()),
        "input": str(input_path),
        "output_dir": str(outdir.resolve()),
        "config_file": args.config or "<built-in defaults>",
        "topology": str(topology) if topology else "<none>",
        "matched_files": len(paths),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "graph_mode": config["graph"]["bond_mode"],
        "ring_sizes": config["ring"]["sizes"],
        "quasi_cage_base_sizes": config["quasi_cage"].get("base_sizes", "auto"),
        "quasi_cage_side_sizes": config["quasi_cage"].get("side_sizes", "auto"),
        "cage_sizes": config["cage"].get("ring_sizes", [5, 6]),
        "output_layout": config["output"].get("structure_layout", "grouped"),
        "workers": n_jobs,
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
    n_jobs: int,
) -> None:
    """Print a compact run header before tqdm starts."""
    print(SQQ_BANNER)
    print("Run information")
    print(f"  working_dir        : {Path.cwd()}")
    print(f"  input              : {input_path}")
    print(f"  output             : {outdir}")
    print(f"  config             : {args.config or '<built-in defaults>'}")
    print(f"  topology           : {topology or '<none>'}")
    print(f"  matched_files      : {len(paths)}")
    if len(paths) == 1:
        print(f"  current_file       : {paths[0]}")
    else:
        print(f"  first_file         : {paths[0]}")
        print(f"  last_file          : {paths[-1]}")
    print(f"  graph_mode         : {config['graph']['bond_mode']}")
    print(f"  ring_sizes         : {config['ring']['sizes']}")
    print(f"  quasi_base/side    : {config['quasi_cage'].get('base_sizes', 'auto')} / {config['quasi_cage'].get('side_sizes', 'auto')}")
    print(f"  cage_sizes         : {config['cage'].get('ring_sizes', [5, 6])}")
    print(f"  other_cages        : {config['cage'].get('output_other', False)}")
    print(f"  output_layout      : {config['output'].get('structure_layout', 'grouped')}")
    print(f"  workers            : {n_jobs}")
    print("")


def analyze_paths_serial(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    topology: Path | None,
    strict: bool,
) -> list[dict[str, Any]]:
    """Analyze frames in input order."""
    rows: list[dict[str, Any]] = []
    frames = read_frames(paths, topology=topology, xtc_stride=int(config["input"].get("xtc_stride", 1)))
    progress = tqdm(frames, total=len(paths), desc="SQQ analyze")
    for frame_index, frame in enumerate(progress):
        if hasattr(progress, "set_description_str"):
            progress.set_description_str(f"SQQ {frame.name}")
        rows.append(process_frame(frame_index, frame, config, outdir, strict=strict, stage_callback=progress_stage_callback(progress, frame.name)))
    return rows


def analyze_paths_parallel(
    paths: list[Path],
    outdir: Path,
    config: dict[str, Any],
    n_jobs: int,
    strict: bool,
) -> list[dict[str, Any]]:
    """Analyze independent coordinate files concurrently."""
    rows_by_index: dict[int, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {
            executor.submit(process_single_file_path, frame_index, path, config, outdir, strict): frame_index
            for frame_index, path in enumerate(paths)
        }
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"SQQ analyze ({n_jobs} jobs)"):
            frame_index, row = future.result()
            rows_by_index[frame_index] = row
    return [rows_by_index[index] for index in sorted(rows_by_index)]


def process_single_file_path(
    frame_index: int,
    path: Path,
    config: dict[str, Any],
    outdir: Path,
    strict: bool,
) -> tuple[int, dict[str, Any]]:
    """Read and analyze one standalone coordinate file."""
    frame = next(iter(read_frames([path])))
    return frame_index, process_frame(frame_index, frame, config, outdir, strict=strict)


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


def progress_stage_callback(progress: Any, frame_name: str) -> Callable[[str], None] | None:
    """Build a tqdm-aware callback for the current frame stage."""
    if not hasattr(progress, "set_postfix_str"):
        return None

    def update(stage: str) -> None:
        if hasattr(progress, "set_description_str"):
            progress.set_description_str(f"SQQ {frame_name}")
        progress.set_postfix_str(stage, refresh=True)

    return update


def report_stage(callback: Callable[[str], None] | None, stage: str) -> None:
    """Update the terminal stage display when a callback is available."""
    if callback is not None:
        callback(stage)


def write_frame_outputs(result: FrameResult, frame_dir: Path, config: dict[str, Any]) -> None:
    """Write all configured per-frame output files."""
    if config["output"].get("write_info", True):
        write_frame_info(result, frame_dir)
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


def resolve_n_jobs(value: Any, n_paths: int) -> int:
    """Parse the frame-level worker count."""
    if value in (None, "", "auto"):
        return 1
    try:
        jobs = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("parallel.n_jobs / --n-jobs must be 'auto' or a positive integer.") from exc
    return max(1, min(jobs, max(1, n_paths)))


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
    if args.pairs:
        config["graph"]["bond_mode"] = "pairs"
        config["graph"]["pair_file"] = args.pairs
    if args.pair_id:
        config["graph"]["pair_id"] = args.pair_id
    if args.n_jobs is not None:
        config["parallel"]["n_jobs"] = args.n_jobs
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
