"""Cross-process lock for one SQQ output root."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import socket
from typing import BinaryIO, Iterator
from uuid import uuid4

from ..exceptions import OutputLockError


OUTPUT_LOCK_NAME = ".sqq.lock"


@dataclass
class OutputLock:
    """One process-held lock for an SQQ output root."""

    path: Path
    handle: BinaryIO
    token: str

    def release(self) -> None:
        if self.handle.closed:
            return
        try:
            self.handle.seek(0)
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(self.handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()


@contextmanager
def output_lock(outdir: Path) -> Iterator[OutputLock]:
    """Prevent concurrent SQQ runs from sharing one output root."""
    lock = acquire_output_lock(Path(outdir))
    try:
        yield lock
    finally:
        lock.release()


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
    "acquire_output_lock",
    "output_lock",
    "output_lock_owner",
]
