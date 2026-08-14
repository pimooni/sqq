"""Mapping overrides retained by the configuration resolver."""

from __future__ import annotations

from argparse import Namespace
from typing import Any

from .validation import normalize_order_parameters


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge explicit values into a configuration mapping."""
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge_config(base[key], value)
        else:
            base[key] = value
    return base


def _on_off(value: Any, option: str) -> bool:
    text = str(value).strip().lower()
    if text == "on":
        return True
    if text == "off":
        return False
    raise ValueError(f"{option} must be on or off.")


def apply_cli_overrides(config: dict[str, Any], args: Namespace) -> None:
    """Apply only the retained public Analyze/Track command-line options."""
    delta_time = getattr(args, "delta_time", None)
    if delta_time is not None:
        config["input"]["delta_time_ps"] = delta_time

    size = getattr(args, "size", None)
    if size:
        config["ring"]["sizes"] = size
        config["quasi_cage"]["base_sizes"] = size
        config["quasi_cage"]["side_sizes"] = size

    for argument, section, option in (
        ("find_half", "half_cage", "--find-half"),
        ("find_quasi", "quasi_cage", "--find-quasi"),
        ("find_cluster", "hydrate_cluster", "--find-cluster"),
    ):
        value = getattr(args, argument, None)
        if value is not None:
            config[section]["enabled"] = _on_off(value, option)

    order_parameter = getattr(args, "order_parameter", None)
    if order_parameter is not None:
        config["order"]["parameters"] = list(
            normalize_order_parameters(order_parameter)
        )

    bond_mode = getattr(args, "bond_mode", None)
    pair_file = getattr(args, "pair", None)
    if pair_file:
        if bond_mode not in (None, "pairs"):
            raise ValueError("--pair can only be combined with --bond-mode pairs.")
        config["graph"]["pair_file"] = pair_file
        config["graph"]["bond_mode"] = "pairs"
    elif bond_mode is not None:
        config["graph"]["bond_mode"] = bond_mode
    if config["graph"].get("bond_mode") == "pairs" and not config["graph"].get("pair_file"):
        raise ValueError(
            "--bond-mode pairs requires --pair PAIRS.txt or graph.pair_file "
            "in sqq_config.yaml."
        )

    worker = getattr(args, "worker", None)
    if worker is not None:
        config["parallel"]["workers"] = worker

    output_type = getattr(args, "output_type", None)
    if output_type is not None:
        config["output"]["types"] = output_type
        if not config.get("hydrate_cluster", {}).get("enabled", False):
            requested = {
                item.strip().lower() for item in str(output_type).split(",")
            }
            if "all" not in requested and {"cluster-gro", "cluster-detail"} & requested:
                raise ValueError("Cluster outputs require --find-cluster on.")


__all__ = ["apply_cli_overrides", "merge_config"]
