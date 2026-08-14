"""Canonical YAML serialization and default template generation."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any, TextIO

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

from .defaults import (
    CONFIG_SCHEMA_VERSION,
    DEFAULT_CONFIG,
    DEFAULT_MODE,
    ORDER_PARAMETER_CHOICES,
)
from .migrate import YAML_KEY_ALIASES, nested_mapping
from .presets import normalize_mode


def canonical_config(config: dict[str, Any]) -> dict[str, Any]:
    """Return a YAML-facing copy using canonical singular field names."""
    data = deepcopy(config)
    engine = data.pop("engine", None)
    internal = data.pop("mode", None)
    if internal is not None:
        engine = internal
    if engine is None:
        engine = DEFAULT_MODE
    graph = data.get("graph")
    if isinstance(graph, dict) and "pairs_file" in graph:
        graph.setdefault("pair_file", graph["pairs_file"])
        graph.pop("pairs_file", None)
    for path, aliases in YAML_KEY_ALIASES.items():
        section = nested_mapping(data, path)
        if section is None:
            continue
        reverse = {internal_name: canonical for canonical, internal_name in aliases.items()}
        converted = {reverse.get(key, key): value for key, value in section.items()}
        section.clear()
        section.update(converted)
    data.pop("schema_version", None)
    data = {("order_parameter" if key == "order" else key): value for key, value in data.items()}
    return {"schema_version": CONFIG_SCHEMA_VERSION, "engine": normalize_mode(engine), **data}


_SECTION_COMMENTS = {
    "run": "Run behavior", "input": "Input discovery, time sampling, and LAMMPS reader settings",
    "component": "Automatic component classification", "additive": "Additive residue names",
    "environment": "Environment or wall residue names", "water": "Water selection",
    "guest": "Guest selection", "graph": "Water-network construction",
    "pbc": "Periodic-boundary settings", "ring": "Ring search and reporting",
    "half_cage": "Standard half-cage search (SQQ-Py)",
    "quasi_cage": "Layered quasi-cage search (SQQ-Py)",
    "cage": "Complete-cage recognition", "hydrate_cluster": "Hydrate phase and cluster recognition (SQQ-Py)",
    "hydrate_order": "MCG and DHOP definitions", "order_parameter": "Order-parameter calculation",
    "ice": "Ice-like water classification (SQQ-Py)", "output": "Output selection and layout",
    "render": "VMD render topology and trajectory", "parallel": "Worker and math-thread policy",
    "track": "Cross-frame persistent cage tracking", "debug": "Developer diagnostics",
}
_INLINE_COMMENTS: dict[tuple[str, ...], str] = {
    ("schema_version",): "managed by SQQ",
    ("engine",): "choices: py, cpp; compatibility presets: 00, 99",
    ("run", "strict"): "choices: true, false",
    ("input", "pattern"): "glob used for directory input",
    ("input", "recursive"): "choices: true, false",
    ("input", "delta_time_ps"): "ps; null analyzes every stored frame",
    ("input", "lammps", "unit"): "choices: real, metal, nano",
    ("input", "lammps", "coordinate_convention"): "choices: auto, x, xs, xu, xsu, unscaled, scaled, unwrapped, scaled_unwrapped",
    ("component", "unknown_action"): "choices: warn, ignore, error",
    ("guest", "center_mode"): "choices: center_atom, centroid, auto",
    ("graph", "mode"): "choices: auto, hbond, oo, pairs",
    ("graph", "pair_file"): "required when graph.mode is pairs",
    ("graph", "pair_id"): "choices: resid, oxygen_index, atomid",
    ("pbc", "box_mode"): "currently orthorhombic only",
    ("ring", "size"): "supported sizes, normally [4, 5, 6]",
    ("ring", "report_size"): "auto or a subset of ring.size",
    ("ring", "definition"): "choices: chordless, shortest_path",
    ("half_cage", "enabled"): "choices: auto, true, false",
    ("quasi_cage", "enabled"): "choices: auto, true, false",
    ("quasi_cage", "search_policy"): "choices: bounded, exact",
    ("cage", "scientific_validation"): "choices: true, false",
    ("hydrate_cluster", "enabled"): "choices: true, false",
    ("order_parameter", "enabled"): f"choices: {', '.join(ORDER_PARAMETER_CHOICES)}",
    ("order_parameter", "q_neighbor_mode"): "choices: graph, cutoff, nearest, lammps",
    ("output", "type"): "default may be combined; all/none are exclusive",
    ("output", "structure_layout"): "choices: grouped, flat",
    ("render", "atom_scope"): "choices: full, compact",
    ("parallel", "backend"): "choices: process, thread, serial",
    ("parallel", "worker"): "auto, integer count, fraction, or percentage",
}


def default_config_template() -> str:
    payload = canonical_config(DEFAULT_CONFIG)
    if yaml is None:
        return json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    raw = yaml.safe_dump(payload, allow_unicode=True, sort_keys=False)
    raw = raw.replace(
        "order_parameter:\n  enabled:\n  - f3\n  - f4\n",
        "order_parameter:\n  enabled: [f3, f4]\n",
    )
    output: list[str] = []
    path_by_depth: dict[int, str] = {}
    pattern = re.compile(r"^(?P<indent> *)(?P<key>[A-Za-z0-9_]+):(?P<value>.*)$")
    for line in raw.splitlines():
        match = pattern.match(line)
        if match is None:
            output.append(line)
            continue
        depth = len(match.group("indent")) // 2
        key = match.group("key")
        path_by_depth[depth] = key
        for stale in tuple(path_by_depth):
            if stale > depth:
                del path_by_depth[stale]
        path = tuple(path_by_depth[index] for index in range(depth + 1))
        if depth == 0 and key in _SECTION_COMMENTS:
            if output and output[-1] != "":
                output.append("")
            output.append(f"# {_SECTION_COMMENTS[key]}")
        comment = _INLINE_COMMENTS.get(path)
        output.append(f"{line}  # {comment}" if comment else line)
    return "\n".join(output) + "\n"


def write_default_config(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Configuration file already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(default_config_template(), encoding="utf-8", newline="\n")


def dump_config(config: dict[str, Any], handle: TextIO) -> None:
    payload = canonical_config(config)
    if yaml is not None:
        yaml.safe_dump(payload, handle, allow_unicode=True, sort_keys=False)
    else:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


__all__ = ["canonical_config", "default_config_template", "dump_config", "write_default_config"]
