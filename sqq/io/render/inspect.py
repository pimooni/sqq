"""Discovery and validation of existing SQQ render packages."""

from __future__ import annotations

import os
from pathlib import Path
import re
import shlex
import sys
from typing import TextIO

from .manifest import RENDER_MANIFEST_BEGIN, RENDER_MANIFEST_END, parse_manifest
from .models import RenderFileReference, RenderPackageInspection


def inspect_render_script(path: Path) -> RenderPackageInspection | None:
    """Inspect one SQQ Tcl script without executing it."""
    script_path = Path(path).resolve()
    try:
        text = script_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        return RenderPackageInspection(
            script_path=script_path,
            kind="unknown",
            references=(),
            source="unreadable",
            errors=(f"cannot read Tcl script: {exc}",),
        )

    if RENDER_MANIFEST_BEGIN in text or RENDER_MANIFEST_END in text:
        references, kind, errors = parse_manifest(text)
        return _validate_references(
            script_path,
            kind,
            references,
            source="manifest",
            errors=errors,
        )

    references = _parse_legacy_references(text)
    if not references or "namespace eval ::SQQ" not in text or not re.search(
        r"(?m)^\s*proc\s+sqq\b", text
    ):
        return None
    return _validate_references(
        script_path,
        "legacy",
        references,
        source="legacy",
        errors=(),
    )


def discover_render_scripts(path: Path) -> list[Path]:
    """Return recognized SQQ VMD Tcl scripts in stable path order."""
    target = Path(path).expanduser()
    if not target.exists():
        raise FileNotFoundError(f"path does not exist: {target}")
    candidates = [target] if target.is_file() else list(target.rglob("*.vmd.tcl"))
    recognized: list[Path] = []
    for candidate in candidates:
        if not candidate.is_file() or not candidate.name.endswith(".vmd.tcl"):
            continue
        inspection = inspect_render_script(candidate)
        if inspection is not None:
            recognized.append(candidate.resolve())
    return sorted(set(recognized), key=lambda item: item.as_posix().casefold())


def run_vmd_command(
    path: str | Path = ".",
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Validate render packages and print commands without launching VMD."""
    output = sys.stdout if stdout is None else stdout
    errors = sys.stderr if stderr is None else stderr
    try:
        scripts = discover_render_scripts(Path(path))
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=errors)
        return 2
    if not scripts:
        print(f"Error: no SQQ VMD render script found under: {Path(path)}", file=errors)
        return 2

    inspections = [inspect_render_script(script) for script in scripts]
    reports = [item for item in inspections if item is not None]
    for index, report in enumerate(reports):
        if index:
            print(file=output)
        prefix = f"[{index + 1}/{len(reports)}] " if len(reports) > 1 else ""
        for line in format_render_package(report, prefix=prefix):
            print(line, file=output)
    return 0 if all(report.complete for report in reports) else 1


def format_render_package(
    report: RenderPackageInspection,
    *,
    prefix: str = "",
) -> list[str]:
    """Format one compact validation report and its two launch commands."""
    status = (
        "SQQ render package: complete"
        if report.complete
        else "Warning: incomplete SQQ render package"
    )
    labels = ["Tcl script", *(item.role for item in report.references)]
    width = max((len(label) for label in labels), default=len("Tcl script"))
    lines = [prefix + status]
    lines.append(f"  {'Tcl script':<{width}} : {_display_path(report.script_path)}")
    for reference in report.references:
        lines.append(f"  {reference.role:<{width}} : {reference.path}")
    for label, values in (
        ("Missing", report.missing),
        ("Empty", report.empty),
        ("Not a file", report.not_files),
        ("Invalid", report.errors),
    ):
        if values:
            lines.append(f"  {label:<{width}} : {', '.join(values)}")
    if not report.complete:
        lines.append("The following commands may fail:")
    script = _display_path(report.script_path)
    lines.extend(
        (
            "VMD Tcl Console:",
            f"  source {_tcl_path_argument(script)}",
            "Terminal:",
            f"  vmd -e {_terminal_path_argument(script)}",
        )
    )
    return lines


def _parse_legacy_references(text: str) -> tuple[RenderFileReference, ...]:
    references: list[RenderFileReference] = []
    for variable, role in (
        ("gro", "topology"),
        ("xtc", "trajectory"),
        ("membership", "membership"),
    ):
        pattern = re.compile(
            rf"^\s*set\s+::SQQ::{variable}_path\s+"
            rf"\[file\s+join\s+\$script_dir\s+"
            rf"(?P<value>\{{[^{{}}\r\n]+\}}|\"[^\"\r\n]+\"|[^\]\s]+)\]\s*$",
            re.MULTILINE,
        )
        match = pattern.search(text)
        if match is None:
            continue
        value = match.group("value")
        if value[:1] in {"{", '"'} and value[-1:] in {"}", '"'}:
            value = value[1:-1]
        references.append(RenderFileReference(role, value, True))
    return tuple(references)


def _validate_references(
    script_path: Path,
    kind: str,
    references: tuple[RenderFileReference, ...],
    *,
    source: str,
    errors: tuple[str, ...],
) -> RenderPackageInspection:
    missing: list[str] = []
    empty: list[str] = []
    not_files: list[str] = []
    validation_errors = list(errors)
    for reference in references:
        if "__SQQ_" in reference.path:
            validation_errors.append(
                f"unresolved placeholder in {reference.role}: {reference.path}"
            )
            continue
        target = Path(reference.path)
        if target.is_absolute() or ".." in target.parts:
            validation_errors.append(
                f"unsafe path in {reference.role}: {reference.path}"
            )
            continue
        target = script_path.parent / target
        if not target.exists():
            if reference.required:
                missing.append(reference.path)
            continue
        if not target.is_file():
            if reference.required:
                not_files.append(reference.path)
            continue
        try:
            if reference.required and target.stat().st_size == 0:
                empty.append(reference.path)
        except OSError as exc:
            validation_errors.append(f"cannot inspect {reference.path}: {exc}")
    return RenderPackageInspection(
        script_path=script_path,
        kind=kind,
        references=references,
        source=source,
        errors=tuple(validation_errors),
        missing=tuple(missing),
        empty=tuple(empty),
        not_files=tuple(not_files),
    )


def _display_path(path: Path) -> str:
    return path.resolve().as_posix()


def _terminal_path_argument(path: str, *, platform: str | None = None) -> str:
    current_platform = os.name if platform is None else platform
    if current_platform == "nt":
        return '"' + path.replace('"', '\\"') + '"'
    return shlex.quote(path)


def _tcl_path_argument(path: str) -> str:
    if "{" not in path and "}" not in path:
        return "{" + path + "}"
    escaped = path.replace("\\", "\\\\")
    for character in ('"', "$", "[", "]"):
        escaped = escaped.replace(character, "\\" + character)
    return '"' + escaped + '"'


__all__ = [
    "discover_render_scripts",
    "format_render_package",
    "inspect_render_script",
    "run_vmd_command",
]
