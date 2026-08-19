"""Cross-process lock for one SQQ output root."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import socket
from typing import Any, BinaryIO, Iterator, MutableMapping
from uuid import uuid4

from ..exceptions import OutputLockError


OUTPUT_LOCK_NAME = ".sqq.lock"


@dataclass
class OutputLock:
    """One process-held lock for an SQQ output root."""

    path: Path
    handle: BinaryIO
    token: str

    def release(self, *, remove_owned_file: bool = False) -> None:
        if self.handle.closed:
            return
        # POSIX permits unlinking the directory entry while this process still
        # holds flock on the open inode.  Doing so before unlock closes the
        # acquire-between-unlock-and-unlink race.  Windows may reject deletion
        # of an open handle, in which case the token-checked retry after close
        # is the best available path; normal 0.5.5 runs never reuse the now
        # non-empty output root and therefore cannot contend for it.
        remove_after_close = False
        try:
            self.handle.seek(0)
            if remove_owned_file:
                remove_after_close = not _remove_owned_lock_file(
                    self.path, self.token
                )
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(self.handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()
        if remove_after_close:
            _remove_owned_lock_file(self.path, self.token)


@dataclass(frozen=True, slots=True)
class OutputDirectorySelection:
    """One exclusively reserved output root and its requested-path metadata."""

    requested: Path
    resolved: Path
    auto_renamed: bool
    lock: OutputLock

    def apply_to_config(self, config: MutableMapping[str, Any]) -> None:
        """Record the output resolution in the effective configuration."""
        output = config.setdefault("output", {})
        if not isinstance(output, MutableMapping):
            raise TypeError("output configuration must be a mapping")
        output.update(
            {
                "requested_path": str(self.requested),
                "resolved_path": str(self.resolved),
                "auto_renamed": self.auto_renamed,
            }
        )


@contextmanager
def output_lock(outdir: Path) -> Iterator[OutputLock]:
    """Prevent concurrent SQQ runs from sharing one output root.

    A cleanly completed context removes only the lock file carrying this
    context's token.  Failed and interrupted contexts deliberately retain the
    lock metadata for diagnosis.
    """
    lock = acquire_output_lock(Path(outdir))
    try:
        yield lock
    except BaseException:
        lock.release()
        raise
    else:
        lock.release(remove_owned_file=True)


@contextmanager
def reserve_output_directory(
    requested: str | Path,
) -> Iterator[OutputDirectorySelection]:
    """Atomically reserve a fresh output root without overwriting old results.

    The requested directory is used when it is absent or empty.  A non-empty
    directory is preserved and the first exclusively creatable ``_NNN``
    sibling is selected.  The returned lock follows the same success/failure
    cleanup policy as :func:`output_lock`.
    """
    selection = acquire_output_directory(Path(requested))
    try:
        yield selection
    except BaseException:
        selection.lock.release()
        raise
    else:
        selection.lock.release(remove_owned_file=True)


def acquire_output_directory(requested: Path) -> OutputDirectorySelection:
    """Resolve and lock one output root, retrying suffixes on contention."""
    requested = requested.expanduser().resolve(strict=False)
    parent = requested.parent
    parent.mkdir(parents=True, exist_ok=True)
    if requested.exists() and not requested.is_dir():
        raise NotADirectoryError(f"Output path is not a directory: {requested}")

    if not requested.exists():
        try:
            requested.mkdir()
        except FileExistsError:
            # Another process claimed the requested name after our check.
            pass
        else:
            return _selection_for_created_directory(requested, requested, False)

    if requested.is_dir() and _directory_is_empty(requested):
        try:
            lock = acquire_output_lock(requested)
        except OutputLockError:
            # A concurrent SQQ process reserved the empty directory first.
            pass
        else:
            other_entries = [
                item for item in requested.iterdir() if item != lock.path
            ]
            if not other_entries:
                return OutputDirectorySelection(
                    requested=requested,
                    resolved=requested,
                    auto_renamed=False,
                    lock=lock,
                )
            # Content appeared between the empty check and lock acquisition.
            lock.release(remove_owned_file=True)

    index = 1
    while True:
        candidate = requested.with_name(f"{requested.name}_{index:03d}")
        try:
            candidate.mkdir()
        except FileExistsError:
            index += 1
            continue
        return _selection_for_created_directory(requested, candidate, True)


def _selection_for_created_directory(
    requested: Path,
    resolved: Path,
    auto_renamed: bool,
) -> OutputDirectorySelection:
    try:
        lock = acquire_output_lock(resolved)
    except BaseException:
        # The directory is known to have been created exclusively by this
        # call, so remove it only when it is still empty.
        try:
            resolved.rmdir()
        except OSError:
            pass
        raise
    return OutputDirectorySelection(
        requested=requested,
        resolved=resolved,
        auto_renamed=auto_renamed,
        lock=lock,
    )


def _directory_is_empty(path: Path) -> bool:
    try:
        next(path.iterdir())
    except StopIteration:
        return True
    return False


def _remove_owned_lock_file(path: Path, token: str) -> bool:
    """Remove a lock only when its on-disk token matches; report success."""
    try:
        metadata = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    if metadata.get("token") != token:
        return False
    try:
        path.unlink()
    except FileNotFoundError:
        return True
    except OSError:
        return False
    return True


def acquire_output_lock(root: Path) -> OutputLock:
    root.mkdir(parents=True, exist_ok=True)
    path = root / OUTPUT_LOCK_NAME
    handle = path.open("a+b")
    try:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        handle.close()
        owner = output_lock_owner(path)
        detail = f" ({owner})" if owner else ""
        raise OutputLockError(
            f"SQQ output directory is already in use: {root}{detail}. "
            "Wait for the active run or choose another --output directory."
        ) from exc

    token = uuid4().hex
    metadata = {
        "format": "SQQ output lock",
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started_at": datetime.now(timezone.utc).isoformat(),
        "token": token,
    }
    handle.seek(0)
    handle.truncate()
    handle.write(
        (json.dumps(metadata, ensure_ascii=True, sort_keys=True) + "\n").encode(
            "ascii"
        )
    )
    handle.flush()
    try:
        os.fsync(handle.fileno())
    except OSError:
        pass
    handle.seek(0)
    return OutputLock(path=path, handle=handle, token=token)


def output_lock_owner(path: Path) -> str:
    try:
        metadata = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return ""
    values = []
    if metadata.get("pid") is not None:
        values.append(f"PID {metadata['pid']}")
    if metadata.get("host"):
        values.append(f"host {metadata['host']}")
    return ", ".join(values)


__all__ = [
    "OUTPUT_LOCK_NAME",
    "OutputLock",
    "OutputDirectorySelection",
    "acquire_output_directory",
    "acquire_output_lock",
    "output_lock",
    "output_lock_owner",
    "reserve_output_directory",
]
