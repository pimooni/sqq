"""Public cage-tracking API.

The implementation is organized by responsibility while this module preserves
the stable import surface used by workflows and downstream callers.
"""

from ...models.tracking import (
    CageObservation,
    CageTrack,
    EventKind,
    FrameStamp,
    TargetSelection,
    TargetSpec,
    TrackCageSnapshot,
    TrackEvent,
    TrackFrameSnapshot,
    TrackingConfig,
    TrackingResult,
)
from .engine import TrackingAccumulator, track_cages, track_snapshots
from .snapshot import snapshot_from_frame_result
from .statistics import (
    cage_lineage_rows,
    event_rows,
    guest_event_rows,
    guest_residence_lifetime_rows,
    guest_residence_rows,
    iter_observation_rows,
    iter_tracking_quality_rows,
    lifetime_distribution_rows,
    lifetime_rows,
    lifetime_survival_rows,
    observation_rows,
    occupancy_state_lifetime_rows,
    occupancy_transition_rows,
    population_rows,
    segment_row,
    segment_track,
    tracking_quality_rows,
)
from .targets import parse_targets, select_targets

__all__ = [
    "CageObservation",
    "CageTrack",
    "EventKind",
    "FrameStamp",
    "TargetSelection",
    "TargetSpec",
    "TrackCageSnapshot",
    "TrackEvent",
    "TrackFrameSnapshot",
    "TrackingAccumulator",
    "TrackingConfig",
    "TrackingResult",
    "event_rows",
    "cage_lineage_rows",
    "guest_event_rows",
    "guest_residence_lifetime_rows",
    "guest_residence_rows",
    "iter_observation_rows",
    "iter_tracking_quality_rows",
    "lifetime_survival_rows",
    "lifetime_distribution_rows",
    "lifetime_rows",
    "observation_rows",
    "occupancy_state_lifetime_rows",
    "occupancy_transition_rows",
    "parse_targets",
    "population_rows",
    "segment_row",
    "segment_track",
    "select_targets",
    "snapshot_from_frame_result",
    "track_cages",
    "track_snapshots",
    "tracking_quality_rows",
]
