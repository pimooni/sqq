"""Public namespace for the SQQ-Py scientific kernel."""

from ..cage import find_cages
from ..graph import build_water_graph
from ..half_quasi import find_cage_patches
from ..ring import find_rings
from .backend import (
    analyze_frame,
    analyze_frame_py,
    filter_free_patches,
    is_subset_of_indexed_owner,
    select_reported_cages,
    subset_owner_index,
)

__all__ = [
    "analyze_frame",
    "analyze_frame_py",
    "build_water_graph",
    "filter_free_patches",
    "find_rings",
    "find_cage_patches",
    "find_cages",
    "is_subset_of_indexed_owner",
    "select_reported_cages",
    "subset_owner_index",
]
