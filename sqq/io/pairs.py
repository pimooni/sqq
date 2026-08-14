from __future__ import annotations

"""Parse pair files into canonical selected-water oxygen atom indices."""

from collections.abc import Iterable
from pathlib import Path
import re

from ..models import Atom, Guest, Water


PAIR_ID_CHOICES = ("resid", "oxygen_index", "atomid")


def read_pair_edges(
    path: str | Path,
    atoms: list[Atom],
    waters: list[Water],
    pair_id: str,
    guests: Iterable[Guest] = (),
) -> list[tuple[int, int]]:
    """Read pair IDs and convert them once to oxygen atom indices.

    ``pair_id`` controls only the public file notation.  The returned edge
    endpoints are always ``Water.oxygen`` indices and are therefore safe to
    pass to either analysis engine.
    """
    pair_path = Path(path)
    resolver = _PairIdResolver(atoms, waters, tuple(guests), pair_id)
    edges: set[tuple[int, int]] = set()
    for lineno, raw_line in enumerate(
        pair_path.read_text(encoding="utf-8-sig").splitlines(), start=1
    ):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = [part for part in re.split(r"[\s,;]+", line) if part]
        if len(parts) < 2:
            raise ValueError(
                f"Invalid pair line {lineno} in {pair_path}: expected two integer IDs."
            )
        try:
            identifiers = (int(parts[0]), int(parts[1]))
        except ValueError as exc:
            raise ValueError(
                f"Pair IDs must be integers at {pair_path}:{lineno}."
            ) from exc
        oxygen_indices = resolver.resolve_line(identifiers, pair_path, lineno)
        left, right = sorted(oxygen_indices)
        if left != right:
            edges.add((left, right))
    return normalize_oxygen_edges(edges, waters, source=f"pair file {pair_path}")


def normalize_oxygen_edges(
    edges: Iterable[tuple[int, int]],
    waters: Iterable[Water],
    *,
    source: str = "pair_edges",
) -> list[tuple[int, int]]:
    """Validate and canonically sort already-normalized oxygen-index edges."""
    oxygen_indices = {int(water.oxygen) for water in waters}
    normalized: set[tuple[int, int]] = set()
    for edge_number, edge in enumerate(edges, start=1):
        try:
            left, right = edge
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{source} edge {edge_number} must contain exactly two oxygen atom indices."
            ) from exc
        if isinstance(left, bool) or isinstance(right, bool):
            raise ValueError(
                f"{source} edge {edge_number} must contain integer oxygen atom indices."
            )
        try:
            left = int(left)
            right = int(right)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{source} edge {edge_number} must contain integer oxygen atom indices."
            ) from exc
        unknown = [value for value in (left, right) if value not in oxygen_indices]
        if unknown:
            values = ", ".join(str(value) for value in unknown)
            raise ValueError(
                f"{source} edge {edge_number} contains unknown selected-water oxygen "
                f"atom index/indices: {values}. Water ordinals, hydrogen atoms, guest "
                "atoms, and mixed identifier schemes are not accepted."
            )
        if left != right:
            normalized.add(tuple(sorted((left, right))))
    return sorted(normalized)


class _PairIdResolver:
    def __init__(
        self,
        atoms: list[Atom],
        waters: list[Water],
        guests: tuple[Guest, ...],
        pair_id: str,
    ) -> None:
        self.atoms = atoms
        self.waters = waters
        self.pair_id = str(pair_id).strip().lower()
        if self.pair_id not in PAIR_ID_CHOICES:
            choices = ", ".join(PAIR_ID_CHOICES)
            raise ValueError(f"graph.pair_id must be one of: {choices}")

        self.oxygen_indices = {int(water.oxygen) for water in waters}
        self.hydrogen_indices = {
            int(index) for water in waters for index in water.hydrogens
        }
        self.water_atom_indices = {
            int(index) for water in waters for index in water.atoms
        }
        self.guest_atom_indices = {
            int(index) for guest in guests for index in guest.atoms
        }
        self.guest_resids = {int(guest.resid) for guest in guests}
        self.by_atomid: dict[int, list[int]] = {}
        self.by_resid: dict[int, list[int]] = {}
        for position, atom in enumerate(atoms):
            self.by_atomid.setdefault(int(atom.atomid), []).append(position)
            self.by_resid.setdefault(int(atom.resid), []).append(position)

        self.maps = {
            "resid": self._unique_water_map(
                ((int(water.resid), int(water.oxygen)) for water in waters), "resid"
            ),
            "oxygen_index": {
                int(water.oxygen): int(water.oxygen) for water in waters
            },
            "atomid": self._unique_water_map(
                (
                    (int(atoms[water.oxygen].atomid), int(water.oxygen))
                    for water in waters
                ),
                "atomid",
            ),
        }

    @staticmethod
    def _unique_water_map(
        identifiers: Iterable[tuple[int, int]], pair_id: str
    ) -> dict[int, int]:
        result: dict[int, int] = {}
        for identifier, oxygen in identifiers:
            if identifier in result and result[identifier] != oxygen:
                raise ValueError(
                    f"graph.pair_id={pair_id!r} is not unique among selected waters; "
                    "use oxygen_index or a unique atomid."
                )
            result[identifier] = oxygen
        return result

    def resolve_line(
        self, identifiers: tuple[int, int], path: Path, lineno: int
    ) -> tuple[int, int]:
        resolved: list[int] = []
        errors: list[str] = []
        for identifier in identifiers:
            oxygen, error = self._resolve_one(identifier)
            if error is not None:
                errors.append(f"{identifier} {error}")
            else:
                resolved.append(oxygen)  # type: ignore[arg-type]
        if errors:
            details = "; ".join(errors)
            raise ValueError(
                f"Invalid pair ID at {path}:{lineno} for graph.pair_id="
                f"{self.pair_id!r}: {details}. Pair files may identify selected water "
                "oxygens only and must not mix identifier schemes."
            )
        return resolved[0], resolved[1]

    def _resolve_one(self, identifier: int) -> tuple[int | None, str | None]:
        oxygen = self.maps[self.pair_id].get(identifier)
        if oxygen is not None:
            ambiguity = self._selected_identity_ambiguity(identifier, oxygen)
            if ambiguity is None:
                return oxygen, None
            return None, ambiguity

        physical_identity = self._physical_identity(identifier)
        if physical_identity is not None:
            return None, f"identifies {physical_identity}, not a selected water oxygen"

        alternate_modes = [
            mode
            for mode, mapping in self.maps.items()
            if mode != self.pair_id and identifier in mapping
        ]
        if alternate_modes:
            modes = "/".join(alternate_modes)
            return None, (
                f"matches selected-water {modes}, not {self.pair_id} "
                "(mixed identifier scheme)"
            )
        if 0 <= identifier < len(self.waters):
            return None, "is a water ordinal, not an oxygen atom index"
        return None, "is unknown"

    def _selected_identity_ambiguity(
        self, identifier: int, oxygen: int
    ) -> str | None:
        if self.pair_id == "atomid":
            positions = self.by_atomid.get(identifier, [])
            if positions != [oxygen]:
                return (
                    "has a mixed or non-unique atom identity; atomid must identify "
                    "exactly one selected water oxygen"
                )
        elif self.pair_id == "resid":
            positions = self.by_resid.get(identifier, [])
            if any(position in self.guest_atom_indices for position in positions):
                return "has mixed water/guest residue identity"
            if any(position not in self.water_atom_indices for position in positions):
                return "has mixed water/non-water residue identity"
            oxygen_matches = [
                water.oxygen for water in self.waters if int(water.resid) == identifier
            ]
            if oxygen_matches != [oxygen]:
                return "has mixed or non-unique water residue identity"
        return None

    def _physical_identity(self, identifier: int) -> str | None:
        if self.pair_id == "oxygen_index":
            if 0 <= identifier < len(self.atoms):
                return self._atom_index_identity(identifier)
            return None
        if self.pair_id == "atomid":
            positions = self.by_atomid.get(identifier, [])
            if positions:
                return self._positions_identity(positions, "atomid")
            return None
        positions = self.by_resid.get(identifier, [])
        if positions:
            return self._positions_identity(positions, "resid")
        if identifier in self.guest_resids:
            return "a guest residue"
        return None

    def _atom_index_identity(self, position: int) -> str:
        identities: list[str] = []
        if position in self.hydrogen_indices:
            identities.append("a water hydrogen atom index")
        elif position in self.water_atom_indices:
            identities.append("a non-oxygen water atom index")
        if position in self.guest_atom_indices:
            identities.append("a guest atom index")
        if 0 <= position < len(self.waters):
            identities.append("a water ordinal")
        if not identities:
            identities.append("a non-water atom index")
        return " and ".join(identities)

    def _positions_identity(self, positions: list[int], label: str) -> str:
        categories: set[str] = set()
        for position in positions:
            if position in self.hydrogen_indices:
                categories.add("water hydrogen")
            elif position in self.oxygen_indices:
                categories.add("water oxygen")
            elif position in self.water_atom_indices:
                categories.add("non-oxygen water")
            if position in self.guest_atom_indices:
                categories.add("guest")
            if (
                position not in self.water_atom_indices
                and position not in self.guest_atom_indices
            ):
                categories.add("non-water")
        description = "/".join(sorted(categories)) or "unknown"
        if len(categories) > 1 or len(positions) > 1:
            return f"mixed or non-unique {description} {label}"
        return f"a {description} {label}"


__all__ = ["PAIR_ID_CHOICES", "normalize_oxygen_edges", "read_pair_edges"]
