from __future__ import annotations

"""Locate and validate SQQ VMD render packages."""

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import shlex
import sys
from typing import Iterable, TextIO


RENDER_MANIFEST_BEGIN = "# SQQ-RENDER-MANIFEST-BEGIN"
RENDER_MANIFEST_END = "# SQQ-RENDER-MANIFEST-END"
RENDER_MANIFEST_SCHEMA = 1

VMD_COMMAND_HELP_LINES = (
    "SQQ VMD commands",
    "",
    "Usage:",
    "  sqq show <family> <target...> ?<family> <target...> ...?",
    "  sqq show label ?on|off?       (lable is accepted)",
    "  sqq color <family> <target...> <VMD-color|ColorID|default>",
    "  sqq pick center|guest|off",
    "  sqq target save",
    "  sqq clear",
    "  sqq help | sqq -h | sqq --help",
    "",
    "Families:",
    "  cage      Cage topology or exact cage ID",
    "  guest     Guests assigned to a cage topology or exact cage ID",
    "  phase     sI, sII, sH, boundary, ambiguous, unclassified, isolated",
    "  cluster   Exact cluster ID",
    "  domain    Exact domain ID",
    "  component water, guest, additive, environment, other, or a residue name",
    "",
    "Interaction:",
    "  Labels are independent of pick mode and are off by default.",
    "  sqq pick center shows yellow pickpoints without enabling labels.",
    "  Pick center enables VMD labelatom mode; pick guest enables VMD pick mode.",
    "  Do not select Mouse > Query manually; Query does not send SQQ pick callbacks.",
    "  Unselected objects are transparent; selected cages are yellow and guests orange.",
    "  sqq target save writes the selected persistent cage ID beside the render files.",
    "  sqq clear restores the source-time cage-all opaque view.",
    "",
    "Examples:",
    "  sqq show cage all",
    "  sqq show cage 512 51264",
    "  sqq show cage 512 guest 512",
    "  sqq show cage all component environment",
    "  sqq show component KLN",
    "  sqq show label",
    "  sqq pick center",
    "  sqq target save",
    "  sqq color cage 512 green",
    "  sqq color component KLN gray",
    "  sqq clear",
)
VMD_COMMAND_HELP = "\n".join(VMD_COMMAND_HELP_LINES)


@dataclass(frozen=True)
class RenderFileReference:
    role: str
    path: str
    required: bool = True


@dataclass(frozen=True)
class RenderPackageInspection:
    script_path: Path
    kind: str
    references: tuple[RenderFileReference, ...]
    source: str
    errors: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()
    empty: tuple[str, ...] = ()
    not_files: tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        return not (self.errors or self.missing or self.empty or self.not_files)


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


def tcl_help_body() -> str:
    """Return Tcl statements that print the shared command help."""
    output: list[str] = []
    for line in VMD_COMMAND_HELP_LINES:
        if any(character in line for character in "{}\\\r\n\0"):
            raise ValueError("SQQ VMD help contains unsupported Tcl characters.")
        output.append(f"    puts {{{line}}}" if line else '    puts ""')
    return "\n".join(output)


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
        references, kind, errors = _parse_embedded_manifest(text)
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


def _parse_embedded_manifest(
    text: str,
) -> tuple[tuple[RenderFileReference, ...], str, tuple[str, ...]]:
    start = text.find(RENDER_MANIFEST_BEGIN)
    end = text.find(RENDER_MANIFEST_END, start + len(RENDER_MANIFEST_BEGIN))
    if start < 0 or end < 0 or end <= start:
        return (), "unknown", ("invalid or incomplete embedded render manifest",)
    body = text[start + len(RENDER_MANIFEST_BEGIN):end]
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
        return (), kind, tuple(errors + ["render manifest files must be a non-empty list"])
    references: list[RenderFileReference] = []
    seen_roles: set[str] = set()
    for index, item in enumerate(raw_files):
        if not isinstance(item, dict):
            errors.append(f"render manifest file {index} must be an object")
            continue
        role = item.get("role")
        file_path = item.get("path")
        required = item.get("required", True)
        if not isinstance(role, str) or not role:
            errors.append(f"render manifest file {index} has an invalid role")
            continue
        if role in seen_roles:
            errors.append(f"duplicate render manifest role: {role}")
            continue
        if not isinstance(file_path, str) or not file_path or "\0" in file_path:
            errors.append(f"render manifest file {role} has an invalid path")
            continue
        if not isinstance(required, bool):
            errors.append(f"render manifest file {role} has a non-boolean required flag")
            continue
        seen_roles.add(role)
        references.append(RenderFileReference(role, file_path, required))
    return tuple(references), kind, tuple(errors)


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
        if value[:1] in {'{', '"'} and value[-1:] in {'}', '"'}:
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
        if not target.is_absolute():
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
