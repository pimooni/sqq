"""Ownership-aware cleanup for files written by SQQ."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
from typing import Any, Mapping
from uuid import uuid4

from ..exceptions import AnalysisError
from .render import RenderSession
from .reporting import clear_previous_summary_outputs
from .reporting.csv_writer import remove_summary_csvs, remove_summary_detail_csvs


OUTPUT_MANIFEST_NAME = "sqq_output_manifest.json"
_OUTPUT_MANIFEST_FORMAT = "SQQ output ownership manifest"
_OUTPUT_MANIFEST_VERSION = 1
_EXCLUDED_NAMES = frozenset({OUTPUT_MANIFEST_NAME, ".sqq.lock"})


@dataclass(frozen=True, slots=True)
class _FileStamp:
    size: int
    modified_ns: int
    changed_ns: int
    file_id: int


class OutputOwnershipSession:
    """Record files created or replaced by one SQQ command.

    The snapshot is taken while the output lock is held. On completion, a
    same-directory atomic replace publishes the ownership list. An unchanged
    pre-existing user file is never claimed by SQQ.
    """

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self._baseline = _snapshot_files(self.root)
        self._baseline_directories = _snapshot_directories(self.root)
        self._committed = False

    def __enter__(self) -> "OutputOwnershipSession":
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if exc_type is None:
            self.commit()
            return
        try:
            self.commit()
        except OSError:
            # Preserve the original analysis failure. A successful run still
            # treats manifest publication failure as fatal.
            pass

    def commit(self) -> Path:
        if self._committed:
            return self.root / OUTPUT_MANIFEST_NAME
        current = _snapshot_files(self.root)
        current_directories = _snapshot_directories(self.root)
        owned = sorted(
            relative
            for relative, stamp in current.items()
            if self._baseline.get(relative) != stamp
        )
        manifest = self.root / OUTPUT_MANIFEST_NAME
        payload = {
            "format": _OUTPUT_MANIFEST_FORMAT,
            "version": _OUTPUT_MANIFEST_VERSION,
            "files": owned,
            "directories": sorted(current_directories - self._baseline_directories),
        }
        _atomic_write_json(manifest, payload)
        self._committed = True
        return manifest


def cleanup_previous_analyze_outputs(
    outdir: Path,
    config: Mapping[str, Any],
    *,
    grouped_roots: bool = True,
) -> None:
    """Remove a previous Analyze run without guessing from file names."""
    root = Path(outdir)
    if cleanup_owned_outputs(root):
        return

    # Pre-manifest compatibility is deliberately conservative. Only fixed SQQ
    # paths are removed; arbitrary ``*_ice*.gro``-style matching is never used.
    child_roots = (
        [
            root / f"result_{label}"
            for label in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            if (root / f"result_{label}").is_dir()
        ]
        if grouped_roots
        else []
    )
    for candidate in [*child_roots, root]:
        cleanup_generated_output_root(candidate, config)
    for candidate in child_roots:
        try:
            candidate.rmdir()
        except OSError:
            pass


def cleanup_owned_outputs(root: str | Path) -> bool:
    """Delete only paths listed by a valid ownership manifest.

    Returns ``True`` when a valid manifest was consumed. Every entry must be
    relative and must still resolve inside the locked output root.
    """
    output_root = Path(root)
    manifest = output_root / OUTPUT_MANIFEST_NAME
    if not manifest.is_file():
        return False
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    if (
        not isinstance(payload, dict)
        or payload.get("format") != _OUTPUT_MANIFEST_FORMAT
        or payload.get("version") != _OUTPUT_MANIFEST_VERSION
        or not isinstance(payload.get("files"), list)
    ):
        return False

    resolved_root = output_root.resolve()
    failures: list[tuple[Path, OSError]] = []
    for raw in payload["files"]:
        if not isinstance(raw, str):
            continue
        relative = Path(raw)
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            continue
        candidate = output_root / relative
        try:
            candidate.resolve(strict=False).relative_to(resolved_root)
        except (OSError, ValueError):
            continue
        try:
            if candidate.is_file() or candidate.is_symlink():
                candidate.unlink(missing_ok=True)
        except OSError as exc:
            failures.append((candidate, exc))
    raw_directories = payload.get("directories", [])
    if isinstance(raw_directories, list):
        directories: list[Path] = []
        for raw in raw_directories:
            if not isinstance(raw, str):
                continue
            relative = Path(raw)
            if relative.is_absolute() or not relative.parts or ".." in relative.parts:
                continue
            candidate = output_root / relative
            try:
                candidate.resolve(strict=False).relative_to(resolved_root)
            except (OSError, ValueError):
                continue
            directories.append(candidate)
        for candidate in sorted(directories, key=lambda item: len(item.parts), reverse=True):
            try:
                candidate.rmdir()
            except OSError:
                pass
    if failures:
        preview = "; ".join(f"{path}: {exc}" for path, exc in failures[:3])
        raise AnalysisError(f"Could not remove previous SQQ output(s): {preview}")
    try:
        manifest.unlink(missing_ok=True)
    except OSError as exc:
        raise AnalysisError(
            f"Could not consume SQQ output manifest {manifest}: {exc}"
        ) from exc
    return True


def cleanup_generated_output_root(
    root: Path,
    config: Mapping[str, Any],
) -> None:
    """Conservatively remove fixed outputs from a pre-manifest SQQ run."""
    path = Path(root)
    if not path.exists():
        return
    mutable_config = deepcopy(dict(config))
    RenderSession.cleanup_output(path)
    clear_previous_summary_outputs(path, mutable_config)
    for name in (
        "summary.xlsx",
        "summary.md",
        "sqq_config_resolved.yaml",
        "run_config.yaml",
    ):
        (path / name).unlink(missing_ok=True)
    remove_summary_csvs(path, mutable_config)
    legacy_config = deepcopy(mutable_config)
    legacy_config.setdefault("output", {})["summary_csv_dir"] = "summary_csv"
    remove_summary_csvs(path, legacy_config)
    remove_summary_detail_csvs(path, mutable_config)
    shutil.rmtree(path / "track", ignore_errors=True)


def _snapshot_files(root: Path) -> dict[str, _FileStamp]:
    if not root.is_dir():
        return {}
    snapshot: dict[str, _FileStamp] = {}
    for candidate in root.rglob("*"):
        if not candidate.is_file() or candidate.name in _EXCLUDED_NAMES:
            continue
        try:
            relative = candidate.relative_to(root).as_posix()
            stat = candidate.stat()
        except OSError:
            continue
        snapshot[relative] = _FileStamp(
            size=int(stat.st_size),
            modified_ns=int(stat.st_mtime_ns),
            changed_ns=int(stat.st_ctime_ns),
            file_id=int(stat.st_ino),
        )
    return snapshot


def _snapshot_directories(root: Path) -> set[str]:
    if not root.is_dir():
        return set()
    directories: set[str] = set()
    for candidate in root.rglob("*"):
        if not candidate.is_dir():
            continue
        try:
            directories.add(candidate.relative_to(root).as_posix())
        except (OSError, ValueError):
            continue
    return directories


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}-{uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


__all__ = [
    "OUTPUT_MANIFEST_NAME",
    "OutputOwnershipSession",
    "cleanup_generated_output_root",
    "cleanup_owned_outputs",
    "cleanup_previous_analyze_outputs",
]
