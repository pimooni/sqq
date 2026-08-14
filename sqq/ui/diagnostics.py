from __future__ import annotations

"""Run-level warning capture isolated from live progress rendering."""

from contextlib import contextmanager
import sys
from threading import Lock
from typing import Any, Iterator
import warnings as python_warnings

from .formatting import write_terminal_block


class RunDiagnostics:
    """Collect, de-duplicate, and emit run warnings outside live progress."""

    def __init__(self, *, stream: Any | None = None) -> None:
        self._stream = sys.stdout if stream is None else stream
        self._messages: list[str] = []
        self._seen: set[str] = set()
        self._emitted = False
        self._lock = Lock()

    def add(self, message: Any) -> None:
        normalized = " ".join(str(message).split())
        if not normalized:
            return
        with self._lock:
            if normalized in self._seen:
                return
            self._seen.add(normalized)
            self._messages.append(normalized)

    def showwarning(
        self,
        message: Warning | str,
        category: type[Warning],
        filename: str,
        lineno: int,
        file: Any | None = None,
        line: str | None = None,
    ) -> None:
        del category, filename, lineno, file, line
        self.add(message)

    def emit(self) -> None:
        messages = self.consume()
        if not messages:
            return
        write_terminal_block(
            [
                "Diagnostics",
                f"  {'warnings':<24}: {len(messages)}",
                *(f"  Warning: {message}" for message in messages),
                "",
            ],
            stream=self._stream,
        )

    def consume(self) -> tuple[str, ...]:
        with self._lock:
            if self._emitted:
                return ()
            self._emitted = True
            return tuple(self._messages)


@contextmanager
def capture_run_warnings(diagnostics: RunDiagnostics) -> Iterator[None]:
    """Route Python warnings into the run-level diagnostics collector."""
    with python_warnings.catch_warnings():
        python_warnings.simplefilter("always", UserWarning)
        python_warnings.showwarning = diagnostics.showwarning
        yield


__all__ = ["RunDiagnostics", "capture_run_warnings"]
