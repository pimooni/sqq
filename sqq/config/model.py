"""Typed configuration-resolution records."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
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
        object.__setattr__(self, "_view", MappingProxyType(dict(self.data)))

    def __getitem__(self, key: str) -> Any:
        return self._view[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._view)

    def __len__(self) -> int:
        return len(self._view)


__all__ = ["ResolutionAdjustment", "ResolutionReport", "ResolvedConfig"]
