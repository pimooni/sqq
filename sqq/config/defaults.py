"""Configuration constants and the immutable default schema."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

DEFAULT_MODE = "py"
CONFIG_SCHEMA_VERSION = "0.5.5"
DEFAULT_ORDER_PARAMETERS = ("f3", "f4")
CPP_MODE = "cpp"
CPP_MODES = frozenset({"99", CPP_MODE})

CPP_DEFAULT_OUTPUT_TYPES = (
    "info",
    "sqq-render",
    "summary-csv",
    "summary-detail-csv",
)
CPP_ALL_OUTPUT_TYPES = (
    "info",
    "gro",
    "cage-gro",
    "f3-gro",
    "f4-gro",
    "sqq-render",
    "summary-csv",
    "summary-xlsx",
    "summary-detail-csv",
)
CPP_OUTPUT_TYPES = frozenset(CPP_ALL_OUTPUT_TYPES)

ALL_ORDER_PARAMETERS = (
    "f3",
    "f4",
    "q6",
    "q12",
    "mcg1",
    "mcg3",
    "dhop35",
    "dhop30",
)
ORDER_PARAMETER_CHOICES = (
    "f3",
    "f4",
    "qN",
    "mcg1",
    "mcg3",
    "dhop35",
    "dhop30",
    "all",
    "none",
)
ORDER_PARAMETER_ALIASES = {
    "mcg-1": "mcg1",
    "mcg_1": "mcg1",
    "mcg-3": "mcg3",
    "mcg_3": "mcg3",
    "dhop-35": "dhop35",
    "dhop_35": "dhop35",
    "dhop-30": "dhop30",
    "dhop_30": "dhop30",
}

DEFAULT_OUTPUT_TYPES = ("info", "sqq-render", "summary-xlsx")
ALL_OUTPUT_TYPES = (
    "info",
    "membership-tsv",
    "order-tsv",
    "f3-gro",
    "f4-gro",
    "gro",
    "sqq-render",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-gro",
    "cluster-detail",
)
OUTPUT_TYPE_ORDER = (
    "info",
    "membership-tsv",
    "order-tsv",
    "f3-gro",
    "f4-gro",
    "gro",
    "sqq-render",
    "ring-gro",
    "half-gro",
    "quasi-gro",
    "cage-gro",
    "ice-gro",
    "summary-xlsx",
    "summary-csv",
    "summary-detail-csv",
    "cluster-gro",
    "cluster-detail",
)
GRO_OUTPUT_TYPES = frozenset(
    {"ring-gro", "half-gro", "quasi-gro", "cage-gro", "ice-gro"}
)

MODE_PRESETS: dict[str, dict[str, Any]] = {
    "00": {
        "engine": "sqq-py",
        "profile": "rigorous",
        "label": "rigorous",
        "worker_fraction": 1.0,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
        "find_cluster": True,
        "output_types": ["info", "sqq-render", "summary-xlsx"],
    },
    "py": {
        "engine": "sqq-py",
        "profile": "standard",
        "label": "standard",
        "worker_count": 1,
        "bond_mode": "auto",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": list(DEFAULT_OUTPUT_TYPES),
    },
    "99": {
        "engine": "sqq-cpp",
        "profile": "performance",
        "label": "sqq-cpp-performance",
        "worker_fraction": 1.0,
        "bond_mode": "hbond",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": list(CPP_DEFAULT_OUTPUT_TYPES),
    },
    "cpp": {
        "engine": "sqq-cpp",
        "profile": "standard",
        "label": "sqq-cpp",
        "worker_count": 1,
        "bond_mode": "auto",
        "ring_sizes": [4, 5, 6],
        "find_cluster": False,
        "output_types": list(CPP_DEFAULT_OUTPUT_TYPES),
    },
}

DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": CONFIG_SCHEMA_VERSION,
    "mode": DEFAULT_MODE,
    "run": {"strict": False},
    "input": {
        "pattern": "*.gro",
        "recursive": False,
        "first_file_time_ps": 0.0,
        "frame_time_step_ps": 100.0,
        "delta_time_ps": None,
        "xyz_scale": 0.1,
        "lammps": {
            "units": "real",
            "timestep": 1.0,
            "atom_style": "full",
            "coordinate_convention": "auto",
            "type_map": {},
        },
    },
    "component": {
        "auto_classify": True,
        "unknown_role": "other",
        "unknown_action": "warn",
        "role_map": {},
    },
    "additive": {"resnames": []},
    "environment": {"resnames": []},
    "water": {
        "resnames": ["SOL", "TIP", "WAT", "HOH"],
        "oxygen_names": ["OW", "O", "OH2"],
        "hydrogen_names": ["HW1", "HW2", "H1", "H2", "HW", "HT1", "HT2"],
    },
    "guest": {
        "resnames": ["CH4", "CO2", "MET", "ETH"],
        "center_atoms": {"CH4": ["C"], "CO2": ["C"], "MET": ["C"]},
        "center_mode": "center_atom",
    },
    "graph": {
        "bond_mode": "auto",
        "oo_cutoff_nm": 0.35,
        "hbond_distance_nm": 0.35,
        "hbond_angle_deg": 30.0,
        "pair_file": None,
        "pair_id": "resid",
    },
    "pbc": {"box_mode": "orthorhombic"},
    "ring": {
        "sizes": [4, 5, 6],
        "report_sizes": "auto",
        "chordless": True,
        "definition": "chordless",
    },
    "half_cage": {"enabled": "auto"},
    "quasi_cage": {
        "search_policy": "bounded",
        "enabled": "auto",
        "base_sizes": "auto",
        "side_sizes": "auto",
        "max_combinations_per_base": 50000,
        "max_layers": 1,
        "max_rings_per_layer": 6,
        "max_layer_states_per_seed": 200,
        "max_candidates_per_edge": 4,
        "max_layer_candidates": 24,
    },
    "cage": {
        "report_types": "auto",
        "max_faces": 20,
        "enabled": True,
        "search_mode": "grow",
        "seed_mode": "ring",
        "max_states_per_seed": 0,
        "max_total_states": 0,
        "max_boundary_candidates": 8,
        "scientific_validation": False,
        "max_face_planarity_rms_nm": 0.06,
        "max_face_edge_cv": 0.35,
        "min_cage_volume_nm3": 1.0e-6,
        "occupancy_mode": "polyhedron",
        "occupancy_radius_nm": 0.5,
    },
    "hydrate_cluster": {"enabled": False, "min_cage": 2},
    "hydrate_order": {
        "mcg_guest_resnames": ["CH4", "MET"],
        "mcg_guest_cutoff_nm": 0.90,
        "mcg_water_cutoff_nm": 0.60,
        "mcg_cone_half_angle_deg": 45.0,
        "mcg_min_waters": 5,
        "dhop_neighbor_cutoff_nm": 0.35,
        "dhop_planar_counts": [11, 12],
        "dhop_min_qualified_neighbors": 3,
    },
    "order": {
        "parameters": ["f3", "f4"],
        "q_neighbor_mode": "graph",
        "q_cutoff_nm": 0.35,
        "q_n_neighbor": None,
        "focus_waters": [],
    },
    "ice": {
        "enabled": True,
        "method": "chill",
        "min_six_rings": 2,
        "require_four_coord_neighbors": True,
    },
    "output": {
        "types": list(DEFAULT_OUTPUT_TYPES),
        "summary_csv_dir": "summary",
        "cage_isomer_rows": "nonzero",
        "write_empty_files": False,
        "structure_layout": "grouped",
        "center_resname": "CNT",
    },
    "render": {"atom_scope": "full"},
    "parallel": {"backend": "process", "workers": "auto", "math_threads": 1},
    "track": {
        "target": "all",
        "source": None,
        "min_jaccard": 0.50,
        "min_shared_fraction": 0.60,
        "min_shared_waters": 3,
        "max_center_distance_nm": None,
        "gap_frame": 0,
        "guest_tiebreak": True,
    },
    "debug": {"use_networkx_checks": False},
}


def default_config_copy() -> dict[str, Any]:
    """Return an isolated copy of the default runtime mapping."""
    return deepcopy(DEFAULT_CONFIG)


__all__ = [
    "ALL_ORDER_PARAMETERS",
    "ALL_OUTPUT_TYPES",
    "CONFIG_SCHEMA_VERSION",
    "CPP_ALL_OUTPUT_TYPES",
    "CPP_DEFAULT_OUTPUT_TYPES",
    "CPP_MODE",
    "CPP_MODES",
    "CPP_OUTPUT_TYPES",
    "DEFAULT_CONFIG",
    "DEFAULT_MODE",
    "DEFAULT_ORDER_PARAMETERS",
    "DEFAULT_OUTPUT_TYPES",
    "GRO_OUTPUT_TYPES",
    "MODE_PRESETS",
    "ORDER_PARAMETER_ALIASES",
    "ORDER_PARAMETER_CHOICES",
    "OUTPUT_TYPE_ORDER",
    "default_config_copy",
]
