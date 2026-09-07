"""Coordinate, topology, trajectory, and explicit-pair input readers."""

from .pairs import PAIR_ID_CHOICES, normalize_oxygen_edges, read_pair_edges

__all__ = ["PAIR_ID_CHOICES", "normalize_oxygen_edges", "read_pair_edges"]
