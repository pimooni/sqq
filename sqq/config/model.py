"""Typed configuration-resolution records."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any


@dataclass(frozen=True)
class ResolutionAdjustment:
    """One automatic preset, migration, or capability adjustment."""

    parameter: str
    requested: Any
    effective: Any
    reason: str


@dataclass(frozen=True)
class ResolutionReport:
    """Auditable changes applied while resolving configuration."""

    requested_selector: str
    engine: str
    profile: str
    adjustments: tuple[ResolutionAdjustment, ...] = ()


@dataclass(frozen=True)
class ResolvedConfig(Mapping[str, Any]):
    """Read-only resolved configuration with its resolution report."""

    data: Mapping[str, Any]
    report: ResolutionReport
    _view: Mapping[str, Any] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        frozen = _freeze(deepcopy(dict(self.data)))
        object.__setattr__(self, "data", frozen)
        object.__setattr__(self, "_view", frozen)

    def __getitem__(self, key: str) -> Any:
        return self._view[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._view)

    def __len__(self) -> int:
        return len(self._view)

    def to_mutable_dict(self) -> dict[str, Any]:
        """Return an isolated mutable copy for an analysis invocation."""
        return _thaw(self._view)


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    if isinstance(value, frozenset):
        return {_thaw(item) for item in value}
    return deepcopy(value)


__all__ = ["ResolutionAdjustment", "ResolutionReport", "ResolvedConfig"]
