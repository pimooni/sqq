from __future__ import annotations

"""CHILL-style ice classification on the shared water graph."""

from dataclasses import dataclass

from ..models import GraphResult, Ring, Water


@dataclass(frozen=True)
class IceClasses:
    """Grouped oxygen-node ids for ice-like local environments."""

    ice_like: tuple[int, ...]
    ice_i: tuple[int, ...]
    interfacial: tuple[int, ...]
    liquid_like: tuple[int, ...]


def classify_ice_waters(
    graph: GraphResult,
    waters: list[Water],
    rings: dict[int, list[Ring]],
    enabled: bool = False,
    min_six_rings: int = 2,
    require_four_coord_neighbors: bool = True,
) -> IceClasses:
    """Classify waters into ice-I-like and interfacial ice classes.

    This follows the CHILL idea used by Moore and Molinero: ice nuclei contain
    both well ordered ice-I local environments and intermediate/interfacial ice.
    SQQ applies this idea to the hydrogen-bond graph already used for rings.
    """
    water_oxygens = {water.oxygen for water in waters}
    if not enabled:
        return IceClasses(ice_like=(), ice_i=(), interfacial=(), liquid_like=tuple(sorted(water_oxygens)))

    six_ring_counts = six_ring_membership_counts(rings)
    four_coord = {
        oxygen
        for oxygen in water_oxygens
        if len(graph.adjacency.get(oxygen, set()) & water_oxygens) == 4
    }
    candidates = {
        oxygen
        for oxygen in four_coord
        if six_ring_counts.get(oxygen, 0) > 0
    }

    ice_i = set()
    for oxygen in candidates:
        neighbors = graph.adjacency.get(oxygen, set()) & water_oxygens
        enough_hex_faces = six_ring_counts.get(oxygen, 0) >= min_six_rings
        ordered_neighbors = (not require_four_coord_neighbors) or all(neighbor in four_coord for neighbor in neighbors)
        if enough_hex_faces and ordered_neighbors:
            ice_i.add(oxygen)

    interfacial = candidates - ice_i
    ice_like = ice_i | interfacial
    liquid_like = water_oxygens - ice_like
    return IceClasses(
        ice_like=tuple(sorted(ice_like)),
        ice_i=tuple(sorted(ice_i)),
        interfacial=tuple(sorted(interfacial)),
        liquid_like=tuple(sorted(liquid_like)),
    )


def find_ice_like_waters(graph: GraphResult, waters: list[Water], rings: dict[int, list[Ring]], enabled: bool = False) -> tuple[int, ...]:
    """Backward-compatible helper returning all ice-like waters."""
    return classify_ice_waters(graph, waters, rings, enabled=enabled).ice_like


def six_ring_membership_counts(rings: dict[int, list[Ring]]) -> dict[int, int]:
    """Count how many six-membered rings contain each oxygen node."""
    counts: dict[int, int] = {}
    for ring in rings.get(6, []):
        for node in ring.nodes:
            counts[node] = counts.get(node, 0) + 1
    return counts
