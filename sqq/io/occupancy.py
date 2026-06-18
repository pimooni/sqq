from __future__ import annotations

"""Shared guest-occupancy labels for cage summaries and GRO output."""

from collections import Counter

from ..models import Cage, Guest


def guest_id(guest: Guest) -> str:
    """Use the same guest identifier as the cage assignment code."""
    return f"{guest.resname}{guest.resid}"


def guest_lookup(guests: list[Guest]) -> dict[str, Guest]:
    """Index guests by the id stored on Cage.guest_ids."""
    return {guest_id(guest): guest for guest in guests}


def guest_resname_order(guests: list[Guest]) -> list[str]:
    """Return guest residue names by their first atom position in the frame."""
    first_index: dict[str, int] = {}
    for guest in guests:
        if not guest.atoms:
            continue
        first_index[guest.resname] = min(first_index.get(guest.resname, guest.atoms[0]), min(guest.atoms))
    return [resname for resname, _ in sorted(first_index.items(), key=lambda item: item[1])]


def cage_guest_names(cage: Cage, lookup: dict[str, Guest]) -> list[str]:
    """Return the residue names of guests assigned to one cage."""
    return [lookup[item].resname for item in cage.guest_ids if item in lookup]


def guest_composition_label(cage: Cage, lookup: dict[str, Guest], resname_order: list[str] | None = None) -> str:
    """Label exact cage occupancy, such as CH4, CH4x2, or CH4+ETH."""
    names = cage_guest_names(cage, lookup)
    if not names:
        return ""
    counts = Counter(names)
    order_index = {name: index for index, name in enumerate(resname_order or [])}
    ordered_names = sorted(counts, key=lambda name: (order_index.get(name, 10_000), name))
    parts = []
    for name in ordered_names:
        count = counts[name]
        parts.append(name if count == 1 else f"{name}x{count}")
    return "+".join(parts)
