"""Shared Tcl renderer generation for Analyze and Track."""

from __future__ import annotations

from ...banner import SQQ_BANNER
from .manifest import render_manifest_block
from .models import (
    RenderFileReference,
    SQQ_CAGE_GRO_NAME,
    SQQ_CAGE_MEMBERSHIP_NAME,
    SQQ_CAGE_XTC_NAME,
)
from .tcl_template import SQQ_CAGE_TCL


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


def tcl_banner_body() -> str:
    """Return ASCII-only Tcl statements for the shared Unicode SQQ banner."""

    def _escaped_line(line: str) -> str:
        output: list[str] = []
        for character in line:
            codepoint = ord(character)
            if character in {'\\', '"', '$', '[', ']'}:
                output.append("\\" + character)
            elif 0x20 <= codepoint <= 0x7E:
                output.append(character)
            elif codepoint <= 0xFFFF:
                output.append(f"\\u{codepoint:04x}")
            else:
                # The current artwork is BMP-only. Keep the generator explicit
                # instead of relying on Tcl-version-specific non-BMP escaping.
                raise ValueError("SQQ VMD banner contains a non-BMP character.")
        return "".join(output)

    return "\n".join(
        f'    puts "{_escaped_line(line)}"' for line in SQQ_BANNER.splitlines()
    )


def tcl_help_body() -> str:
    """Return Tcl statements that print the shared command help."""
    output: list[str] = []
    for line in VMD_COMMAND_HELP_LINES:
        if any(character in line for character in "{}\\\r\n\0"):
            raise ValueError("SQQ VMD help contains unsupported Tcl characters.")
        output.append(f"    puts {{{line}}}" if line else '    puts ""')
    return "\n".join(output)


def tcl_braced_literal(value: str, *, label: str) -> str:
    """Return a restricted, substitution-free Tcl literal."""
    text = str(value)
    if not text:
        raise ValueError(f"{label} must not be empty.")
    if any(character in text for character in "{}\\\r\n\0"):
        raise ValueError(f"{label} contains unsupported Tcl characters.")
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"{label} must contain ASCII text only.") from exc
    return "{" + text + "}"


def vmd_script_text(
    gro_filename: str = SQQ_CAGE_GRO_NAME,
    xtc_filename: str = SQQ_CAGE_XTC_NAME,
    membership_filename: str = SQQ_CAGE_MEMBERSHIP_NAME,
    molecule_name: str = "SQQ cages",
    render_kind: str = "analyze",
) -> str:
    """Return the shared, parameterized ASCII Tcl renderer."""
    replacements = {
        "__SQQ_GRO_FILENAME__": (gro_filename, "VMD GRO filename"),
        "__SQQ_XTC_FILENAME__": (xtc_filename, "VMD XTC filename"),
        "__SQQ_MEMBERSHIP_FILENAME__": (
            membership_filename,
            "VMD membership filename",
        ),
        "__SQQ_MOLECULE_NAME__": (molecule_name, "VMD molecule name"),
    }
    script = SQQ_CAGE_TCL.replace(
        "__SQQ_RENDER_MANIFEST__",
        render_manifest_block(
            kind=render_kind,
            files=(
                RenderFileReference("topology", str(gro_filename)),
                RenderFileReference("trajectory", str(xtc_filename)),
                RenderFileReference("membership", str(membership_filename)),
            ),
        ),
    ).replace("__SQQ_HELP_BODY__", tcl_help_body()).replace(
        "__SQQ_BANNER_BODY__", tcl_banner_body()
    )
    for placeholder, (filename, label) in replacements.items():
        value = str(filename)
        if value in {".", ".."} or "/" in value or "\\" in value:
            raise ValueError(f"{label} must name a file in the script directory.")
        script = script.replace(
            placeholder,
            tcl_braced_literal(value, label=label),
        )
    if "__SQQ_" in script:
        raise AssertionError("Unresolved SQQ Tcl template placeholder.")
    return script


__all__ = [
    "SQQ_CAGE_TCL",
    "VMD_COMMAND_HELP",
    "VMD_COMMAND_HELP_LINES",
    "tcl_braced_literal",
    "tcl_banner_body",
    "tcl_help_body",
    "vmd_script_text",
]
