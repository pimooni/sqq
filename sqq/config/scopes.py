from __future__ import annotations

"""Normalize complete analysis scopes independently of CLI workflows."""

import math
from pathlib import Path
from typing import Any
import warnings as python_warnings

from .defaults import DEFAULT_MODE
from .presets import is_cpp_mode
from .validation import (
    normalize_cpp_output_types,
    normalize_order_parameters,
    normalize_output_types,
    normalize_parallel_backend,
    output_enabled,
)

def normalize_analysis_scopes(config: dict[str, Any]) -> None:
    """Normalize search and report scopes before frames are analyzed."""
    input_config = config.setdefault("input", {})
    input_config["recursive"] = parse_on_off(input_config.get("recursive", False), "input.recursive")
    graph_config = config.setdefault("graph", {})
    if str(graph_config.get("bond_mode", "auto")).strip().lower() == "pairs":
        pair_file = graph_config.get("pair_file")
        if pair_file in (None, ""):
            raise ValueError(
                "bond_mode=pairs requires --pair PAIRS.txt or graph.pair_file "
                "in sqq_config.yaml."
            )
        pair_path = Path(str(pair_file)).expanduser()
        if not pair_path.is_absolute():
            pair_path = Path.cwd() / pair_path
        if not pair_path.is_file():
            raise ValueError(f"Pair file does not exist or is not a file: {pair_path}")
        graph_config["pair_file"] = str(pair_path.resolve())
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
    removed_sampling_keys = {
        "trajectory_stride",
        "xtc_stride",
    }.intersection(input_config)
    if removed_sampling_keys:
        names = ", ".join(sorted(removed_sampling_keys))
        raise ValueError(
            f"Unsupported input sampling key(s): {names}. Use input.delta_time_ps."
        )
    raw_delta_time = input_config.get("delta_time_ps")
    input_config["delta_time_ps"] = (
        None
        if raw_delta_time in (None, "")
        else finite_float(
            raw_delta_time,
            "input.delta_time_ps / -dt / --delta-time",
            positive=True,
        )
    )
    raw_lammps = input_config.get("lammps", {})
    if not isinstance(raw_lammps, dict):
        raise ValueError("input.lammps must be a mapping.")
    from ..io.lammps import normalize_lammps_config

    lammps_values = dict(raw_lammps)
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
    graph.pop("effective_bond_mode_reason", None)
    graph.pop("effective_bond_mode_by_group", None)
    graph.pop("effective_bond_mode_reason_by_group", None)
    graph_mode = str(graph.get("bond_mode", "auto")).strip().lower()
    if graph_mode not in {"auto", "hbond", "oo", "pairs"}:
        raise ValueError("graph.bond_mode must be auto, hbond, oo, or pairs.")
    graph["bond_mode"] = graph_mode
    graph["oo_cutoff_nm"] = finite_float(graph.get("oo_cutoff_nm", 0.35), "graph.oo_cutoff_nm", positive=True)
    graph["hbond_distance_nm"] = finite_float(
        graph.get("hbond_distance_nm", 0.35),
        "graph.hbond_distance_nm",
        positive=True,
    )
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

    component = config.setdefault("component", {})
    component["auto_classify"] = parse_on_off(
        component.get("auto_classify", True),
        "component.auto_classify",
    )
    unknown_role = str(component.get("unknown_role", "other")).strip().lower()
    if unknown_role not in {"water", "guest", "additive", "environment", "other"}:
        raise ValueError(
            "component.unknown_role must be water, guest, additive, "
            "environment, or other."
        )
    component["unknown_role"] = unknown_role
    unknown_action = str(component.get("unknown_action", "warn")).strip().lower()
    if unknown_action not in {"warn", "ignore", "error"}:
        raise ValueError("component.unknown_action must be warn, ignore, or error.")
    component["unknown_action"] = unknown_action
    if not isinstance(component.get("role_map", {}), dict):
        raise ValueError("component.role_map must be a mapping.")

    additive = config.setdefault("additive", {})
    additive["resnames"] = string_list(
        additive.get("resnames", []),
        "additive.resnames",
        allow_empty=True,
    )
    environment = config.setdefault("environment", {})
    environment["resnames"] = string_list(
        environment.get("resnames", []),
        "environment.resnames",
        allow_empty=True,
    )

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

    half = config.setdefault("half_cage", {})
    half["enabled"] = parse_auto_on_off(
        half.get("enabled", "auto"),
        "half_cage.enabled",
        auto_value=not is_cpp_mode(config.get("mode", DEFAULT_MODE)),
    )

    quasi = config.setdefault("quasi_cage", {})
    quasi["enabled"] = parse_auto_on_off(
        quasi.get("enabled", "auto"),
        "quasi_cage.enabled",
        auto_value=not is_cpp_mode(config.get("mode", DEFAULT_MODE)),
    )
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
    if not quasi["enabled"] and quasi["max_layers"] != 1:
        raise ValueError("quasi_cage.max_layers requires quasi_cage.enabled=true.")

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
    cage["occupancy_radius_nm"] = finite_float(
        cage.get("occupancy_radius_nm", 0.5),
        "cage.occupancy_radius_nm",
        positive=True,
    )
    for key in ("max_states_per_seed", "max_total_states"):
        cage[key] = nonnegative_integer(cage.get(key, 0), f"cage.{key}")
    cage["max_boundary_candidates"] = positive_integer(
        cage.get("max_boundary_candidates", 8),
        "cage.max_boundary_candidates",
    )
    cage["scientific_validation"] = parse_on_off(
        cage.get("scientific_validation", False),
        "cage.scientific_validation",
    )
    cage["max_face_planarity_rms_nm"] = finite_float(
        cage.get("max_face_planarity_rms_nm", 0.06),
        "cage.max_face_planarity_rms_nm",
        nonnegative=True,
    )
    cage["max_face_edge_cv"] = finite_float(
        cage.get("max_face_edge_cv", 0.35),
        "cage.max_face_edge_cv",
        nonnegative=True,
    )
    cage["min_cage_volume_nm3"] = finite_float(
        cage.get("min_cage_volume_nm3", 1.0e-6),
        "cage.min_cage_volume_nm3",
        positive=True,
    )
    hydrate_cluster = config.setdefault("hydrate_cluster", {})
    if "detail" in hydrate_cluster:
        raise ValueError(
            "hydrate_cluster.detail is no longer supported; "
            "add cluster-detail to output.types."
        )
    hydrate_cluster["enabled"] = parse_on_off(
        hydrate_cluster.get("enabled", False), "hydrate_cluster.enabled"
    )
    hydrate_cluster["min_cage"] = positive_integer(
        hydrate_cluster.get("min_cage", 2),
        "hydrate_cluster.min_cage / --cluster-min-cage",
    )
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
    half_output_requested = "gro" in output_types or "half-gro" in output_types
    quasi_output_requested = "gro" in output_types or "quasi-gro" in output_types
    if half_output_requested and not half["enabled"]:
        raise ValueError("half-gro output requires --find-half on.")
    if quasi_output_requested and not quasi["enabled"]:
        raise ValueError("quasi-gro output requires --find-quasi on.")
    output["types"] = output_types
    output["write_empty_files"] = parse_on_off(
        output.get("write_empty_files", False),
        "output.write_empty_files",
    )
    output.pop("summary_detail_dir", None)
    key = "summary_csv_dir"
    default = "summary"
    directory = str(output.get(key, default)).strip() or default
    directory_path = Path(directory)
    if directory_path.is_absolute() or ".." in directory_path.parts:
        raise ValueError(f"output.{key} must be a relative directory inside the output folder.")
    output[key] = directory
    cage_isomer_rows = str(output.get("cage_isomer_rows", "nonzero")).strip().lower()
    if cage_isomer_rows not in {"nonzero", "all"}:
        raise ValueError("output.cage_isomer_rows / --cage-isomer-rows must be nonzero or all.")
    output["cage_isomer_rows"] = cage_isomer_rows
    structure_layout = str(output.get("structure_layout", "grouped")).strip().lower()
    if structure_layout not in {"grouped", "flat"}:
        raise ValueError("output.structure_layout must be grouped or flat.")
    output["structure_layout"] = structure_layout
    render = config.setdefault("render", {})
    atom_scope = str(render.get("atom_scope", "full")).strip().lower()
    if atom_scope not in {"full", "compact"}:
        raise ValueError("render.atom_scope must be full or compact.")
    render["atom_scope"] = atom_scope
    guest_center_mode = str(config.get("guest", {}).get("center_mode", "center_atom")).strip().lower()
    if guest_center_mode not in {"center_atom", "centroid", "auto"}:
        raise ValueError("guest.center_mode must be center_atom, centroid, or auto.")
    config["guest"]["center_mode"] = guest_center_mode
    order = config.setdefault("order", {})
    order["parameters"] = list(
        normalize_order_parameters(order.get("parameters", ["f3", "f4"]))
    )
    from ..core.order.steinhardt import normalize_q_neighbor_mode, resolve_q_neighbor_count

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
    selected_order_parameters = set(order["parameters"])
    for parameter in ("f3", "f4"):
        if output_enabled(config, f"{parameter}-gro") and parameter not in selected_order_parameters:
            raise ValueError(
                f"{parameter}-gro output requires {parameter} in --order-parameter "
                "or order_parameter.enabled."
            )
    per_water_order = any(
        name in {"f3", "f4"} or name.startswith("q")
        for name in order["parameters"]
    )
    if output_enabled(config, "order-tsv") and not per_water_order:
        python_warnings.warn(
            "output type 'order-tsv' has no per-water F3/F4/Q_l selection; "
            "no order-parameter TSV will be written.",
            UserWarning,
            stacklevel=2,
        )
    hydrate_order = config.setdefault("hydrate_order", {})
    positive_values = (
        ("mcg_guest_cutoff_nm", 0.90),
        ("mcg_water_cutoff_nm", 0.60),
        ("dhop_neighbor_cutoff_nm", 0.35),
    )
    for key, default in positive_values:
        hydrate_order[key] = finite_float(hydrate_order.get(key, default), f"hydrate_order.{key}", positive=True)
    cone_angle = finite_float(
        hydrate_order.get("mcg_cone_half_angle_deg", 45.0),
        "hydrate_order.mcg_cone_half_angle_deg",
    )
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
    hydrate_order["mcg_guest_resnames"] = string_list(
        hydrate_order.get("mcg_guest_resnames", ["CH4", "MET"]),
        "hydrate_order.mcg_guest_resnames",
    )

    ice = config.setdefault("ice", {})
    ice["enabled"] = parse_on_off(ice.get("enabled", True), "ice.enabled")
    ice_method = str(ice.get("method", "chill")).strip().lower()
    if ice_method != "chill":
        raise ValueError("ice.method must be chill.")
    ice["method"] = ice_method
    ice["min_six_rings"] = positive_integer(ice.get("min_six_rings", 2), "ice.min_six_rings")
    ice["require_four_coord_neighbors"] = parse_on_off(
        ice.get("require_four_coord_neighbors", True),
        "ice.require_four_coord_neighbors",
    )


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


def parse_auto_on_off(value: Any, key: str, *, auto_value: bool) -> bool:
    """Resolve an auto/on/off setting to one effective boolean."""
    if str(value).strip().lower() == "auto":
        return bool(auto_value)
    return parse_on_off(value, key)


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


def nonnegative_integer(value: Any, key: str) -> int:
    """Normalize one integer that may use zero for no limit."""
    if isinstance(value, bool):
        raise ValueError(f"{key} must be a non-negative integer.")
    try:
        numeric = float(value)
        number = int(numeric)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a non-negative integer.") from exc
    if not math.isfinite(numeric) or numeric != number or number < 0:
        raise ValueError(f"{key} must be a non-negative integer.")
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
    from ..core.cage import (
        CAGE_REPORT_GROUPS,
        TARGET_FACE_COUNTS,
        canonical_cage_type,
        parse_cage_face_label,
    )

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

__all__ = [
    "finite_float",
    "nonnegative_integer",
    "normalize_analysis_scopes",
    "parse_auto_on_off",
    "parse_on_off",
    "positive_integer",
    "resolve_cage_report_types",
    "resolve_size_list",
    "string_list",
]
