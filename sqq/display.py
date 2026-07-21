from __future__ import annotations

"""Shared human-facing display helpers for terminal, Markdown, and summary output."""

from collections import Counter
from typing import Any, Iterable

_GRAPH_MODE_ORDER = ("hbond", "oo", "pairs")
_EMPTY_GRAPH_MODE_VALUES = {"", "none", "null", "nan"}


def clean_graph_mode(value: Any) -> str:
    """Return a compact graph-mode identifier suitable for UI display."""
    text = str(value).strip()
    return "" if text.lower() in _EMPTY_GRAPH_MODE_VALUES else text


def ordered_unique_graph_modes(values: Iterable[Any]) -> list[str]:
    """Return stable unique effective graph modes, with known modes first."""
    seen: set[str] = set()
    modes: list[str] = []
    for value in values:
        mode = clean_graph_mode(value)
        if not mode or mode in seen:
            continue
        seen.add(mode)
        modes.append(mode)
    return sorted(modes, key=lambda item: (_GRAPH_MODE_ORDER.index(item) if item in _GRAPH_MODE_ORDER else len(_GRAPH_MODE_ORDER), item))


def graph_mode_display(requested: Any, effective_modes: Iterable[Any] | None = None) -> str:
    """Format requested and effective graph modes for reports."""
    requested_mode = clean_graph_mode(requested) or "auto"
    if requested_mode != "auto":
        return requested_mode
    values = [
        mode
        for value in ([] if effective_modes is None else effective_modes)
        if (mode := clean_graph_mode(value))
    ]
    modes = ordered_unique_graph_modes(values)
    if len(modes) == 1:
        return f"auto -> {modes[0]}"
    if len(modes) > 1:
        counts = Counter(values)
        details = ", ".join(f"{mode}: {counts[mode]}" for mode in modes)
        return f"auto -> mixed ({details})"
    return "auto -> pending"
