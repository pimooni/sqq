"""Smoke tests for pair-file graph input and CHILL-style ice classes."""

import numpy as np

from sqq.core.graph import build_water_graph
from sqq.core.ice import classify_ice_waters
from sqq.models import Atom, Ring, Water


def make_water(resid: int, oxygen: int) -> Water:
    return Water(resid=resid, resname="SOL", oxygen=oxygen, hydrogens=(), atoms=(oxygen,))


def test_pair_file_graph_uses_resids(tmp_path):
    atoms = [Atom(index=i, resid=i + 1, resname="SOL", atomname="OW", atomid=i + 1, xyz=np.zeros(3)) for i in range(3)]
    waters = [make_water(1, 0), make_water(2, 1), make_water(3, 2)]
    pair_file = tmp_path / "pairs.txt"
    pair_file.write_text("1 2\n2 3\n", encoding="utf-8")

    graph = build_water_graph(
        atoms,
        waters,
        box=None,
        bond_mode="pairs",
        oo_cutoff_nm=0.35,
        hbond_distance_nm=0.35,
        hbond_angle_deg=30.0,
        pair_file=pair_file,
        pair_id="resid",
    )

    assert graph.mode == "pairs"
    assert graph.edges == [(0, 1), (1, 2)]
    assert graph.adjacency[1] == {0, 2}


def test_chill_style_ice_classes_split_ordered_and_interfacial():
    waters = [make_water(i + 1, i) for i in range(8)]
    adjacency = {
        0: {1, 2, 3, 4},
        1: {0, 2, 4, 5},
        2: {0, 1, 3, 5},
        3: {0, 2, 4, 6},
        4: {0, 1, 3, 6},
        5: {1, 2, 6, 7},
        6: {3, 4, 5},
        7: {5},
    }
    graph = type("Graph", (), {"mode": "pairs", "edges": [], "adjacency": adjacency})()
    rings = {6: [Ring(object_id="ring6_00001", nodes=(0, 1, 2, 3, 4, 5))]}

    classes = classify_ice_waters(graph, waters, rings, enabled=True, min_six_rings=1)

    assert 0 in classes.ice_i
    assert 5 in classes.interfacial
    assert set(classes.ice_like) == set(classes.ice_i) | set(classes.interfacial)
