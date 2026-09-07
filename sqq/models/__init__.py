"""Stable public data-model API."""

from .cage_type import (
    CAGE_REPORT_GROUPS,
    KNOWN_CAGE_TYPES,
    TARGET_FACE_COUNTS,
    cage_type_for_counts,
    canonical_cage_face_label,
    canonical_cage_type,
    parse_cage_face_label,
)
from .order import ClusterOrderValue, F3F4Result, HydrateOrderResult, WaterOrder
from .phase import HydrateCluster, HydrateDomain, HydrateMotif
from .result import FrameResult
from .structure import Atom, Frame, Guest, Water, guest_id
from .topology import Cage, CagePatch, GraphResult, Ring
from .tracking import (
    CageObservation,
    CageTrack,
    EventKind,
    FrameStamp,
    MatchStatus,
    Row,
    TargetKind,
    TargetSelection,
    TargetSpec,
    TrackCageSnapshot,
    TrackEvent,
    TrackFrameSnapshot,
    TrackingConfig,
    TrackingResult,
)

__all__ = [
    "CAGE_REPORT_GROUPS",
    "KNOWN_CAGE_TYPES",
    "TARGET_FACE_COUNTS",
    "Atom",
    "Frame",
    "Water",
    "Guest",
    "Ring",
    "CagePatch",
    "Cage",
    "HydrateMotif",
    "HydrateDomain",
    "HydrateCluster",
    "WaterOrder",
    "F3F4Result",
    "ClusterOrderValue",
    "HydrateOrderResult",
    "GraphResult",
    "FrameResult",
    "TrackingConfig",
    "FrameStamp",
    "TrackCageSnapshot",
    "TrackFrameSnapshot",
    "CageObservation",
    "CageTrack",
    "TrackEvent",
    "TrackingResult",
    "TargetSpec",
    "TargetSelection",
    "EventKind",
    "MatchStatus",
    "Row",
    "TargetKind",
    "guest_id",
    "cage_type_for_counts",
    "canonical_cage_face_label",
    "canonical_cage_type",
    "parse_cage_face_label",
]
