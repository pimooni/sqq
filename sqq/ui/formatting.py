"""Atomic terminal text formatting shared by headers and progress views."""

from __future__ import annotations

from datetime import datetime
import sys
from typing import Any


TERMINAL_LABEL_WIDTH = 24

_MONTH_NAMES = (
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
)
_WEEKDAY_NAMES = ("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")


def write_terminal_block(lines: list[str], *, stream: Any | None = None) -> None:
    """Write one complete terminal block without cross-stream interleaving."""
    target = sys.stdout if stream is None else stream
    text = "\n".join(lines) + "\n"
    encoding = getattr(target, "encoding", None) or "utf-8"
    safe_text = text.encode(encoding, errors="replace").decode(
        encoding,
        errors="replace",
    )
    target.write(safe_text)
    target.flush()


def terminal_field_line(label: str, value: Any) -> str:
    """Format one aligned terminal key-value row."""
    return f"  {label:<{TERMINAL_LABEL_WIDTH}}: {safe_terminal_text(value)}"


def print_terminal_field(label: str, value: Any) -> None:
    """Print one aligned terminal key-value row."""
    write_terminal_block([terminal_field_line(label, value)])


def format_terminal_value(value: Any) -> str:
    """Format terminal values compactly without Python container brackets."""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (list, tuple, set)):
        return ",".join(str(item) for item in value)
    return str(value)


def safe_terminal_text(value: Any) -> str:
    """Avoid UnicodeEncodeError on legacy Windows consoles."""
    text = ascii_superscript_text(format_terminal_value(value))
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    return text.encode(encoding, errors="replace").decode(encoding, errors="replace")


SUPERSCRIPT_DIGITS = {
    "⁰": "0",
    "¹": "1",
    "²": "2",
    "³": "3",
    "⁴": "4",
    "⁵": "5",
    "⁶": "6",
    "⁷": "7",
    "⁸": "8",
    "⁹": "9",
    "⁻": "-",
}


def ascii_superscript_text(text: str) -> str:
    """Convert Unicode superscripts to compact ASCII exponents."""
    result: list[str] = []
    in_superscript = False
    for char in text:
        if char in SUPERSCRIPT_DIGITS:
            if not in_superscript:
                result.append("^")
            result.append(SUPERSCRIPT_DIGITS[char])
            in_superscript = True
            continue
        if in_superscript and char.isdigit():
            result.append(" ")
        in_superscript = False
        result.append(char)
    return "".join(result)


TIME_ZONE_ALIASES = {
    ("CST", 480): "China Standard Time",
    ("中国标准时间", 480): "China Standard Time",
}


def format_time_zone(value: datetime) -> str:
    """Format a time-zone name with its signed UTC offset."""
    name = value.tzname() or "UTC"
    offset = value.utcoffset()
    if offset is None:
        return name
    total_minutes = int(offset.total_seconds() / 60)
    name = TIME_ZONE_ALIASES.get((name, total_minutes), name)
    sign = "+" if total_minutes >= 0 else "-"
    hours, minutes = divmod(abs(total_minutes), 60)
    offset_text = f"{sign}{hours}" if minutes == 0 else f"{sign}{hours}:{minutes:02d}"
    return f"{name} ({offset_text})"


def format_started(value: datetime) -> str:
    """Format a locale-independent, compact run start timestamp."""
    offset = value.utcoffset()
    if offset is None:
        offset_text = "UTC"
    else:
        total_minutes = int(offset.total_seconds() / 60)
        sign = "+" if total_minutes >= 0 else "-"
        hours, minutes = divmod(abs(total_minutes), 60)
        offset_text = f"UTC{sign}{hours:02d}:{minutes:02d}"
    month = _MONTH_NAMES[value.month - 1]
    weekday = _WEEKDAY_NAMES[value.weekday()]
    return (
        f"{month} {value.day}, {value.year} ({weekday}), "
        f"{value:%H:%M:%S} ({offset_text})"
    )


def format_seconds(seconds: float) -> str:
    """Format elapsed seconds for the live terminal display."""
    return f"{max(seconds, 0.0):.1f} s"


__all__ = [
    "SUPERSCRIPT_DIGITS",
    "TERMINAL_LABEL_WIDTH",
    "TIME_ZONE_ALIASES",
    "ascii_superscript_text",
    "format_seconds",
    "format_started",
    "format_terminal_value",
    "format_time_zone",
    "print_terminal_field",
    "safe_terminal_text",
    "terminal_field_line",
    "write_terminal_block",
]
