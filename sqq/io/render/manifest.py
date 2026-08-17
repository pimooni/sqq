"""Embedded manifest support for SQQ VMD packages."""

from __future__ import annotations

import json
from typing import Iterable

from .models import RenderFileReference


RENDER_MANIFEST_BEGIN = "# SQQ-RENDER-MANIFEST-BEGIN"
RENDER_MANIFEST_END = "# SQQ-RENDER-MANIFEST-END"
RENDER_MANIFEST_SCHEMA = 1


def render_manifest_block(
    *,
    kind: str,
    files: Iterable[RenderFileReference],
) -> str:
    """Return an ASCII JSON manifest embedded as Tcl comments."""
    records = [
        {"path": item.path, "required": bool(item.required), "role": item.role}
        for item in files
    ]
    payload = {
        "files": records,
        "format": "SQQ VMD render package",
        "kind": str(kind),
        "schema": RENDER_MANIFEST_SCHEMA,
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )
    return f"{RENDER_MANIFEST_BEGIN}\n# {encoded}\n{RENDER_MANIFEST_END}"


def parse_manifest(
    text: str,
) -> tuple[tuple[RenderFileReference, ...], str, tuple[str, ...]]:
    """Parse one embedded manifest without resolving its file references."""
    start = text.find(RENDER_MANIFEST_BEGIN)
    end = text.find(RENDER_MANIFEST_END, start + len(RENDER_MANIFEST_BEGIN))
    if start < 0 or end < 0 or end <= start:
        return (), "unknown", ("invalid or incomplete embedded render manifest",)
    body = text[start + len(RENDER_MANIFEST_BEGIN) : end]
    payload_text = "\n".join(
        line.lstrip()[1:].lstrip() if line.lstrip().startswith("#") else line
        for line in body.splitlines()
        if line.strip()
    )
    try:
        payload = json.loads(payload_text)
    except json.JSONDecodeError as exc:
        return (), "unknown", (f"invalid embedded render manifest JSON: {exc.msg}",)
    errors: list[str] = []
    if not isinstance(payload, dict):
        return (), "unknown", ("render manifest root must be an object",)
    if payload.get("schema") != RENDER_MANIFEST_SCHEMA:
        errors.append(f"unsupported render manifest schema: {payload.get('schema')!r}")
    kind = payload.get("kind", "unknown")
    if not isinstance(kind, str) or not kind:
        errors.append("render manifest kind must be non-empty text")
        kind = "unknown"
    raw_files = payload.get("files")
    if not isinstance(raw_files, list) or not raw_files:
        return (), kind, tuple(
            errors + ["render manifest files must be a non-empty list"]
        )
    references: list[RenderFileReference] = []
    seen_roles: set[str] = set()
    for index, item in enumerate(raw_files):
        if not isinstance(item, dict):
            errors.append(f"render manifest file {index} must be an object")
            continue
        role = item.get("role")
        path = item.get("path")
        required = item.get("required", True)
        if not isinstance(role, str) or not role:
            errors.append(f"render manifest file {index} has an invalid role")
            continue
        if role in seen_roles:
            errors.append(f"duplicate render manifest role: {role}")
            continue
        if not isinstance(path, str) or not path or "\0" in path:
            errors.append(f"render manifest file {role} has an invalid path")
            continue
        if not isinstance(required, bool):
            errors.append(
                f"render manifest file {role} has a non-boolean required flag"
            )
            continue
        seen_roles.add(role)
        references.append(RenderFileReference(role, path, required))
    return tuple(references), kind, tuple(errors)


__all__ = [
    "RENDER_MANIFEST_BEGIN",
    "RENDER_MANIFEST_END",
    "RENDER_MANIFEST_SCHEMA",
    "parse_manifest",
    "render_manifest_block",
]
