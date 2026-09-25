"""Run-private disk backing for large tracking observation/event streams."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import replace
from itertools import islice
from pathlib import Path
import pickle
import sqlite3
from typing import Generic, TypeVar

from ...models.tracking import CageObservation, TrackEvent


T = TypeVar("T")


class LazyFilteredSequence(Sequence[T], Generic[T]):
    """Repeatable, bounded-memory filtered view over another sequence."""

    def __init__(self, source: Sequence[T], predicate: Callable[[T], bool]) -> None:
        self._source = source
        self._predicate = predicate
        self._count: int | None = None

    def __iter__(self) -> Iterator[T]:
        return (item for item in self._source if self._predicate(item))

    def __len__(self) -> int:
        if self._count is None:
            self._count = sum(1 for _item in self)
        return self._count

    def __getitem__(self, index: int | slice) -> T | tuple[T, ...]:
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            return tuple(islice(self, start, stop, step))
        normalized = index if index >= 0 else len(self) + index
        if normalized < 0:
            raise IndexError(index)
        try:
            return next(islice(self, normalized, normalized + 1))
        except StopIteration as exc:
            raise IndexError(index) from exc

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Sequence):
            return NotImplemented
        return len(self) == len(other) and all(
            left == right for left, right in zip(self, other)
        )


class _SqliteSequence(Sequence[T], Generic[T]):
    """Read-only repeatable view over pickled rows in one finalized SQLite file."""

    def __init__(
        self,
        path: Path,
        query: str,
        parameters: tuple[object, ...],
        count: int,
        decoder: Callable[[bytes, int], T],
    ) -> None:
        self._path = path
        self._query = query
        self._parameters = parameters
        self._count = count
        self._decoder = decoder

    def __len__(self) -> int:
        return self._count

    def __iter__(self) -> Iterator[T]:
        connection = sqlite3.connect(f"file:{self._path}?mode=ro", uri=True)
        try:
            cursor = connection.execute(self._query, self._parameters)
            for index, row in enumerate(cursor):
                yield self._decoder(bytes(row[0]), index)
        finally:
            connection.close()

    def __getitem__(self, index: int | slice) -> T | tuple[T, ...]:
        if isinstance(index, slice):
            start, stop, step = index.indices(self._count)
            return tuple(islice(self, start, stop, step))
        normalized = index if index >= 0 else self._count + index
        if normalized < 0 or normalized >= self._count:
            raise IndexError(index)
        query = f"SELECT payload FROM ({self._query}) LIMIT 1 OFFSET ?"
        connection = sqlite3.connect(f"file:{self._path}?mode=ro", uri=True)
        try:
            row = connection.execute(
                query, (*self._parameters, normalized)
            ).fetchone()
        finally:
            connection.close()
        if row is None:
            raise IndexError(index)
        return self._decoder(bytes(row[0]), normalized)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Sequence):
            return NotImplemented
        return len(self) == len(other) and all(
            left == right for left, right in zip(self, other)
        )


class SpoolObservationBucket(Sequence[CageObservation]):
    """List-like append sink that becomes a lazy per-track sequence at finish."""

    def __init__(self, spool: "TrackingSpool", track_id: str, track_number: int) -> None:
        self._spool = spool
        self.track_id = track_id
        self.track_number = track_number
        self.first_frame_index: int | None = None
        self._count = 0

    def append(self, observation: CageObservation) -> None:
        if self.first_frame_index is None:
            self.first_frame_index = observation.frame_index
        self._spool.append_observation(observation, self.track_number)
        self._count += 1

    def __len__(self) -> int:
        return self._count

    def __iter__(self) -> Iterator[CageObservation]:
        return iter(self.sequence())

    def __getitem__(
        self, index: int | slice
    ) -> CageObservation | tuple[CageObservation, ...]:
        return self.sequence()[index]

    def sequence(self) -> Sequence[CageObservation]:
        self._spool.require_finalized()
        return _SqliteSequence(
            self._spool.path,
            "SELECT payload FROM observation WHERE track_id = ? "
            "ORDER BY frame_index, local_cage_id, ordinal",
            (self.track_id,),
            self._count,
            _decode_pickle,
        )


class _SpoolEventSink:
    def __init__(self, spool: "TrackingSpool") -> None:
        self._spool = spool

    def append(self, event: TrackEvent) -> None:
        self._spool.append_event(event)


class TrackingSpool:
    """SQLite-backed append store owned by one Analyze/raw-Track run."""

    def __init__(self, directory: str | Path) -> None:
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        self.path = root / "tracking.sqlite3"
        self._connection = sqlite3.connect(self.path)
        self._connection.execute("PRAGMA synchronous = OFF")
        self._connection.execute("PRAGMA journal_mode = OFF")
        self._connection.execute("PRAGMA temp_store = FILE")
        self._connection.execute("PRAGMA cache_size = -4096")
        self._connection.executescript(
            """
            CREATE TABLE observation (
                ordinal INTEGER PRIMARY KEY,
                track_id TEXT NOT NULL,
                track_number INTEGER NOT NULL,
                frame_index INTEGER NOT NULL,
                local_cage_id TEXT NOT NULL,
                payload BLOB NOT NULL
            );
            CREATE TABLE event (
                ordinal INTEGER PRIMARY KEY,
                frame_index INTEGER NOT NULL,
                payload BLOB NOT NULL
            );
            """
        )
        self._observation_count = 0
        self._event_count = 0
        self._writes_since_commit = 0
        self._finalized = False
        self.events = _SpoolEventSink(self)

    def observation_bucket(
        self, track_id: str, track_number: int
    ) -> SpoolObservationBucket:
        return SpoolObservationBucket(self, track_id, track_number)

    def append_observation(
        self, observation: CageObservation, track_number: int
    ) -> None:
        self._require_writable()
        self._connection.execute(
            "INSERT INTO observation VALUES (?, ?, ?, ?, ?, ?)",
            (
                self._observation_count,
                observation.track_id,
                track_number,
                observation.frame_index,
                observation.local_cage_id,
                sqlite3.Binary(pickle.dumps(observation, protocol=5)),
            ),
        )
        self._observation_count += 1
        self._checkpoint()

    def append_event(self, event: TrackEvent) -> None:
        self._require_writable()
        self._connection.execute(
            "INSERT INTO event VALUES (?, ?, ?)",
            (
                self._event_count,
                event.frame_index,
                sqlite3.Binary(pickle.dumps(event, protocol=5)),
            ),
        )
        self._event_count += 1
        self._checkpoint()

    def finalize(self) -> None:
        if self._finalized:
            return
        self._connection.executescript(
            """
            CREATE INDEX observation_track_order
                ON observation(track_id, frame_index, local_cage_id, ordinal);
            CREATE INDEX observation_frame_order
                ON observation(frame_index, track_number, local_cage_id, ordinal);
            CREATE INDEX event_frame_order ON event(frame_index, ordinal);
            """
        )
        self._connection.commit()
        self._connection.close()
        self._finalized = True

    def observations(self) -> Sequence[CageObservation]:
        self.require_finalized()
        return _SqliteSequence(
            self.path,
            "SELECT payload FROM observation "
            "ORDER BY frame_index, track_number, local_cage_id, ordinal",
            (),
            self._observation_count,
            _decode_pickle,
        )

    def event_sequence(self) -> Sequence[TrackEvent]:
        self.require_finalized()
        return _SqliteSequence(
            self.path,
            "SELECT payload FROM event ORDER BY frame_index, ordinal",
            (),
            self._event_count,
            _decode_event,
        )

    def require_finalized(self) -> None:
        if not self._finalized:
            raise RuntimeError("Tracking spool must be finalized before it is read.")

    def close(self) -> None:
        if self._finalized:
            return
        self._connection.close()
        self._finalized = True

    def _checkpoint(self) -> None:
        self._writes_since_commit += 1
        if self._writes_since_commit >= 1024:
            self._connection.commit()
            self._writes_since_commit = 0

    def _require_writable(self) -> None:
        if self._finalized:
            raise RuntimeError("Tracking spool is already finalized.")


def filtered_sequence(
    source: Sequence[T], predicate: Callable[[T], bool]
) -> Sequence[T]:
    """Return a lazy filter unless the caller already wants every item."""
    return LazyFilteredSequence(source, predicate)


def _decode_pickle(payload: bytes, _index: int) -> T:
    return pickle.loads(payload)


def _decode_event(payload: bytes, index: int) -> TrackEvent:
    return replace(pickle.loads(payload), event_id=f"e{index + 1}")


__all__ = [
    "SpoolObservationBucket",
    "TrackingSpool",
    "filtered_sequence",
]
