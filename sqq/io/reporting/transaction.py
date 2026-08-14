from __future__ import annotations

"""Low-level atomic publication shared by reporting sinks."""

import os
from pathlib import Path
import tempfile


def temporary_output_path(target: Path) -> Path:
    """Create a same-directory temporary path for atomic replacement."""
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, raw_path = tempfile.mkstemp(
        prefix=f".{target.stem}.",
        suffix=target.suffix,
        dir=target.parent,
    )
    os.close(descriptor)
    return Path(raw_path)


def commit_output_bundle(
    pending: list[tuple[Path, Path]],
    removals: list[Path] | tuple[Path, ...] = (),
) -> None:
    """Publish related outputs together and restore old files on failure."""
    targets = list(dict.fromkeys([target for _, target in pending] + list(removals)))
    backups: dict[Path, Path] = {}
    committed: list[Path] = []
    try:
        for target in targets:
            if not target.exists():
                continue
            backup = temporary_output_path(target)
            backup.unlink(missing_ok=True)
            os.replace(target, backup)
            backups[target] = backup
        for temp_path, target in pending:
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(temp_path, target)
            committed.append(target)
    except Exception:
        for target in committed:
            target.unlink(missing_ok=True)
        for target, backup in backups.items():
            if backup.exists():
                os.replace(backup, target)
        raise
    finally:
        for backup in backups.values():
            backup.unlink(missing_ok=True)
