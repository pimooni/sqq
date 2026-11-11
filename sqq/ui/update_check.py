"""Quiet, cached update checks for interactive CLI runs."""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import sys
from threading import Thread
from time import time
from urllib.request import Request, urlopen

from .. import __version__


_LATEST_RELEASE_API = "https://api.github.com/repos/pimooni/sqq/releases/latest"
_VERSION_PATTERN = re.compile(r"^v?(\d+)\.(\d+)\.(\d+)$")
_CACHE_LIFETIME_SECONDS = 24 * 60 * 60
_CHECK_TIMEOUT_SECONDS = 1.0
_thread: Thread | None = None
_latest_tag: str | None = None


def _version_tuple(value: str) -> tuple[int, int, int] | None:
    match = _VERSION_PATTERN.fullmatch(value.strip())
    return tuple(map(int, match.groups())) if match else None


def _cache_path() -> Path:
    if sys.platform == "win32":
        root = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        root = Path.home() / "Library" / "Caches"
    else:
        root = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    return root / "sqq" / "latest_release.json"


def _read_cached_tag() -> str | None:
    try:
        data = json.loads(_cache_path().read_text(encoding="utf-8"))
        if time() - float(data["checked_at"]) > _CACHE_LIFETIME_SECONDS:
            return None
        tag = str(data["tag_name"])
        return tag if _version_tuple(tag) is not None else None
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _fetch_latest_tag() -> None:
    global _latest_tag
    cached = _read_cached_tag()
    if cached is not None:
        _latest_tag = cached
        return
    request = Request(
        _LATEST_RELEASE_API,
        headers={"Accept": "application/vnd.github+json", "User-Agent": f"SQQ/{__version__}"},
    )
    try:
        with urlopen(request, timeout=_CHECK_TIMEOUT_SECONDS) as response:
            data = json.loads(response.read(1_000_000))
        tag = str(data["tag_name"])
        if _version_tuple(tag) is None:
            return
        _latest_tag = tag
        path = _cache_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps({"checked_at": time(), "tag_name": tag}), encoding="utf-8"
        )
    except (OSError, ValueError, KeyError, TypeError):
        # An offline or restricted host should run exactly as before.
        return


def start_update_check() -> None:
    """Start one background check for an interactive Analyze/Track invocation."""
    global _thread
    if _thread is not None or not sys.stdout.isatty():
        return
    if os.environ.get("SQQ_NO_UPDATE_CHECK", "").strip().lower() in {"1", "true", "yes", "on"}:
        return
    _thread = Thread(target=_fetch_latest_tag, name="sqq-update-check", daemon=True)
    _thread.start()


def available_update() -> str | None:
    """Return a safe, copyable notice if a newer official release is known."""
    if _thread is None:
        return None
    _thread.join(timeout=_CHECK_TIMEOUT_SECONDS)
    local = _version_tuple(__version__)
    latest = _version_tuple(_latest_tag or "")
    if local is None or latest is None or latest <= local:
        return None
    return f"{_latest_tag} available: https://github.com/pimooni/sqq/releases/tag/{_latest_tag}"
