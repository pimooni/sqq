"""Shared command-line banner text."""

from __future__ import annotations


SQQ_TITLE = "Shell  Quant  Qualifier"
SQQ_AUTHOR = "by J. PANG & Q. SUN"

_INFO_PANEL = (
    "+-----------------------------+",
    "|   Shell  Quant  Qualifier   |",
    "|     by J. PANG & Q. SUN     |",
    "+-----------------------------+",
)


def _with_info_panel(artwork: str) -> str:
    """Place the standard information panel two columns after four-line artwork."""
    return "\n".join(
        f"{art}  {info}" for art, info in zip(artwork.splitlines(), _INFO_PANEL)
    )


DEFAULT_SOO_ART = """
┏━━━━━┳━━━━━┳━━━━━┓
┃  ━━━┫  O  ┃  O  ┃
┣━━━  ┣━━━┓ ┣━━━┓ ┃
┗━━━━━┛   ┗━┛   ┗━┛
""".strip()

PRESET_S09_ART = """
┏━━━━━┳━━━━━┳━━━━━┓
┃  ━━━┫  0  ┃  9  ┃
┣━━━  ┣━━━┓ ┣━━━┓ ┃
┗━━━━━┛   ┗━┛   ┗━┛
""".strip()

SQQ_BANNER = _with_info_panel(DEFAULT_SOO_ART)
PRESET_BANNER = _with_info_panel(PRESET_S09_ART)

ENGINE_BADGES = {
    "00": """
┏━━━━━┳━━━━━┓
┃  0  ┃  0  ┃
┃     ┃     ┃
┗━━━━━┻━━━━━┛
""".strip(),
    "99": """
┏━━━━━┳━━━━━┓
┃  9  ┃  9  ┃
┣━━   ┣━━   ┃
┗━━━━━┻━━━━━┛
""".strip(),
}


def banner_for_engine(engine: str | None = None) -> str:
    """Return the base banner with an optional compatibility-preset badge."""
    badge = ENGINE_BADGES.get(str(engine or "").strip().lower())
    if badge is None:
        return SQQ_BANNER
    return "\n".join(
        f"{base}  {extra}"
        for base, extra in zip(PRESET_BANNER.splitlines(), badge.splitlines())
    )

HELP_BANNER = f"""
{SQQ_BANNER}
SQQ (Shell Quant Qualifier): Python Joint Toolkit for Water-Shell Topology Analysis.
""".strip()


# ---------------------------------------------------------------------------
# Artwork assets
# Retained design material below is not used by runtime banner output.
# ---------------------------------------------------------------------------

LEGACY_SQQ_ART = """
┏━━━━━┳━━━━━┳━━━━━┓
┃  ═══┃  ║  ┃  ║  ┃
┃═══  ┣━━━┓ ┣━━━┓ ┃
┗━━━━━┛   ┗━┛   ┗━┛
""".strip()

IRSM_BANNER = """
┏━━━┳━━━┳━━━┳━━━━━┓
┗┓ ┏┫ ┏━┫ ══┃ ┃ ┃ ┃
┏┛ ┗┫ ┃ ┃══ ┃ ┃ ┃ ┃
┗━━━┻━┛ ┗━━━┻━┻━┻━┛
""".strip()

LEGACY_ENGINE_BADGES = {
    "00": """
┏━━━━━┳━━━━━┓
┃  ║  ┃  ║  ┃
┃     ┃     ┃
┗━━━━━┻━━━━━┛
""".strip(),
    "99": """
┏━━━━━┳━━━━━┓
┃  ║  ┃  ║  ┃
┃══   ┃══   ┃
┗━━━━━┻━━━━━┛
""".strip(),
}
