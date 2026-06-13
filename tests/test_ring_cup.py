"""Smoke tests for the first graph-topology primitives."""

import numpy as np
import pytest

from sqq.core.cage import build_cage_targets, counts_fit, find_cages, grow_seed_face_sets
from sqq.core.cup import canonical_sequence
from sqq.core.ring import find_rings
from sqq.models import Cup, Frame, Ring
from sqq.pipeline import resolve_size_list


def test_find_single_five_ring():
    # A clean pentagon should yield exactly one canonical 5-ring.
    adjacency = {
        1: {2, 5},
        2: {1, 3},
        3: {2, 4},
        4: {3, 5},
        5: {4, 1},
    }

    rings = find_rings(adjacency, sizes=[5])

    assert len(rings[5]) == 1
    assert rings[5][0].nodes == (1, 2, 3, 4, 5)


def test_canonical_sequence_separates_isomers():
    # These two 5/6 side-ring sequences are intentionally distinct isomers.
    assert canonical_sequence([5, 5, 5, 6, 6]) == "55566"
    assert canonical_sequence([5, 5, 6, 5, 6]) == "55656"


def test_cup_size_auto_follows_ring_sizes():
    ring_sizes = resolve_size_list([4, 5, 6, 7], fallback=[], key="ring.sizes")

    assert resolve_size_list("auto", fallback=ring_sizes, key="cup.base_sizes") == [4, 5, 6, 7]
    assert resolve_size_list(None, fallback=ring_sizes, key="cup.side_sizes") == [4, 5, 6, 7]
    assert resolve_size_list([5, 6], fallback=ring_sizes, key="cup.side_sizes") == [5, 6]


def test_cage_grow_uses_cup_seeds_by_default():
    rings = {
        f"ring5_{idx:05d}": Ring(object_id=f"ring5_{idx:05d}", nodes=(idx, idx + 1, idx + 2, idx + 3, idx + 4))
        for idx in range(1, 4)
    }
    cup = Cup(
        object_id="cup5_55555_00001",
        cup_type="cup5_55555",
        rings=tuple(rings),
        waters=tuple(range(10)),
        center=np.zeros(3),
    )

    assert grow_seed_face_sets([cup], rings, "cup") == [frozenset(rings)]
    assert grow_seed_face_sets([cup], rings, "ring") == [frozenset([ring_id]) for ring_id in sorted(rings)]


def test_default_cage_targets_exclude_four_member_faces():
    targets = build_cage_targets(["512"], {5, 6}, output_other=False, other_max_faces=20)

    assert targets == {"512": {4: 0, 5: 12, 6: 0}}
    assert not counts_fit({4: 1}, targets["512"])


def test_other_cage_targets_can_include_four_member_faces():
    targets = build_cage_targets(["512"], {4, 5, 6}, output_other=True, other_max_faces=14)

    assert "4^1-5^10-6^2" in targets
    assert targets["4^1-5^10-6^2"] == {4: 1, 5: 10, 6: 2}


def test_cage_search_rejects_seven_member_faces():
    frame = Frame(name="empty", atoms=[], box=None)

    with pytest.raises(ValueError, match="supports only 4, 5, and 6"):
        find_cages(frame, {}, [], [], enabled=True, ring_sizes=[4, 5, 6, 7])
