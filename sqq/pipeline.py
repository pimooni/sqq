from __future__ import annotations

"""Top-level analysis pipeline for the SQQ command line."""

from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - exercised in minimal source-tree runs.
    def tqdm(iterable, total=None, desc=None):
        return iterable

from .config import load_config
from .core.cage import find_cages
from .core.cup import find_cups
from .core.f3f4 import compute_f3f4
from .core.graph import build_water_graph
from .core.ice import classify_ice_waters
from .core.ring import find_rings
from .core.selection import select_guests, select_waters
from .io.gro_writer import write_cage_gro_files, write_cup_gro_files, write_ice_gro_file, write_ring_gro_files
from .io.summary import failed_row, result_row, write_f3f4, write_frame_info, write_membership, write_summary, write_vmd_script
from .io.trajectory import expand_inputs, read_frames
from .models import Cage, Cup, Frame, FrameResult


PARALLEL_SUFFIXES = {".gro", ".xyz"}


def analyze(args: Namespace) -> None:
    """Run SQQ analysis from parsed command-line arguments."""
    config = load_config(Path(args.config) if args.config else None)
    apply_cli_overrides(config, args)

    # Directory input follows the HA-style one-file-per-frame workflow.
    input_path = Path(args.input)
    pattern = args.pattern or config["input"]["pattern"]
    recursive = bool(args.recursive or config["input"]["recursive"])
    paths = expand_inputs(input_path, pattern=pattern, recursive=recursive)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    n_jobs = resolve_n_jobs(config["parallel"].get("n_jobs"), len(paths))
    topology = Path(args.topology) if args.topology else None
    if n_jobs > 1 and can_parallelize_paths(paths, topology):
        rows = analyze_paths_parallel(paths, outdir, config, n_jobs=n_jobs, strict=bool(args.strict))
    else:
        rows = analyze_paths_serial(paths, outdir, config, topology=topology, strict=bool(args.strict))

    write_summary(rows, outdir, config, write_xlsx=config["output"]["write_xlsx_summary"])
    print(f"Wrote SQQ results: {outdir}")


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
    for frame_index, frame in enumerate(tqdm(frames, total=len(paths), desc="SQQ analyze")):
        rows.append(process_frame(frame_index, frame, config, outdir, strict=strict))
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


def process_frame(frame_index: int, frame: Frame, config: dict[str, Any], outdir: Path, strict: bool) -> dict[str, Any]:
    """Analyze one frame, write per-frame files, and return a summary row."""
    if frame.time_ps is None:
        frame.time_ps = config["input"]["first_file_time_ps"] + frame_index * config["input"]["frame_time_step_ps"]
    try:
        result = analyze_frame(frame, config)
        frame_dir = outdir / frame.name
        frame_dir.mkdir(parents=True, exist_ok=True)
        write_frame_outputs(result, frame_dir, config)
        return result_row(result)
    except Exception as exc:
        if strict:
            raise
        return failed_row(frame.name, str(frame.source or ""), str(exc))


def write_frame_outputs(result: FrameResult, frame_dir: Path, config: dict[str, Any]) -> None:
    """Write all configured per-frame output files."""
    write_frame_info(result, frame_dir)
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
        write_ring_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"])
        write_cup_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"])
        write_cage_gro_files(result, frame_dir, write_empty=config["output"]["write_empty_files"])
        write_ice_gro_file(result, frame_dir, write_empty=config["output"]["write_empty_files"])


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
    if args.pattern:
        config["input"]["pattern"] = args.pattern
    if args.recursive:
        config["input"]["recursive"] = True
    if args.pairs:
        config["graph"]["bond_mode"] = "pairs"
        config["graph"]["pair_file"] = args.pairs
    if args.pair_id:
        config["graph"]["pair_id"] = args.pair_id
    if args.n_jobs is not None:
        config["parallel"]["n_jobs"] = args.n_jobs
    if args.no_gro:
        config["output"]["write_gro"] = False
    if args.no_xlsx:
        config["output"]["write_xlsx_summary"] = False


def analyze_frame(frame: Frame, config: dict[str, Any]) -> FrameResult:
    """Analyze one frame and return all topology objects for export."""
    ring_sizes = resolve_size_list(config["ring"]["sizes"], fallback=[], key="ring.sizes")
    cup_base_sizes = resolve_size_list(config["cup"].get("base_sizes", "auto"), fallback=ring_sizes, key="cup.base_sizes")
    cup_side_sizes = resolve_size_list(config["cup"].get("side_sizes", "auto"), fallback=ring_sizes, key="cup.side_sizes")
    cage_ring_sizes = resolve_size_list(config["cage"].get("ring_sizes", [5, 6]), fallback=ring_sizes, key="cage.ring_sizes")
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
    # Ring, cup, cage, F3/F4, and ice all consume the same water graph.
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
    rings = find_rings(
        graph.adjacency,
        sizes=ring_sizes,
        chordless=bool(config["ring"]["chordless"]),
    )
    cups = find_cups(
        frame,
        rings,
        enabled=bool(config["cup"].get("enabled", False)),
        base_sizes=cup_base_sizes,
        side_sizes=cup_side_sizes,
        max_combinations_per_base=int(config["cup"].get("max_combinations_per_base", 50000)),
    )
    cages = find_cages(
        frame,
        rings,
        cups,
        guests,
        enabled=bool(config["cage"].get("enabled", False)),
        target_types=list(config["cage"].get("target_types", [])),
        ring_sizes=cage_ring_sizes,
        output_other=bool(config["cage"].get("output_other", False)),
        other_max_faces=int(config["cage"].get("other_max_faces", 20)),
        search_mode=str(config["cage"].get("search_mode", "grow")),
        seed_mode=str(config["cage"].get("seed_mode", "cup")),
        max_states_per_seed=int(config["cage"].get("max_states_per_seed", 2000)),
        max_total_states=int(config["cage"].get("max_total_states", 250000)),
        occupancy_radius_nm=float(config["cage"].get("occupancy_radius_nm", 0.5)),
        occupancy_mode=str(config["cage"].get("occupancy_mode", "polyhedron")),
    )
    cups = filter_free_cups(cups, cages)
    focus_resids = {int(item) for item in config["order"].get("focus_waters", [])}
    f3f4 = compute_f3f4(frame, waters, graph, focus_resids=focus_resids) if bool(config["order"].get("f3f4_enabled", True)) else None
    ice_classes = classify_ice_waters(
        graph,
        waters,
        rings,
        enabled=bool(config["ice"].get("enabled", False)),
        min_six_rings=int(config["ice"].get("min_six_rings", 2)),
        require_four_coord_neighbors=bool(config["ice"].get("require_four_coord_neighbors", True)),
    )
    warnings = []
    if not cups:
        warnings.append("No closed cup was found with the current general cup criteria.")
    return FrameResult(
        frame=frame,
        waters=waters,
        guests=guests,
        graph=graph,
        rings=rings,
        cups=cups,
        cages=cages,
        f3f4=f3f4,
        ice_like_waters=ice_classes.ice_like,
        ice_i_waters=ice_classes.ice_i,
        interfacial_ice_waters=ice_classes.interfacial,
        warnings=warnings,
    )


def filter_free_cups(cups: list[Cup], cages: list[Cage]) -> list[Cup]:
    """Remove cups whose ring faces are fully consumed by a detected cage."""
    cage_ring_sets = [set(cage.rings) for cage in cages]
    if not cage_ring_sets:
        return cups
    free_cups = []
    for cup in cups:
        cup_rings = set(cup.rings)
        if any(cup_rings <= cage_rings for cage_rings in cage_ring_sets):
            continue
        free_cups.append(cup)
    return free_cups


def resolve_size_list(value: Any, fallback: list[int], key: str) -> list[int]:
    """Resolve ring-size settings, allowing cup sizes to follow ring sizes."""
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
