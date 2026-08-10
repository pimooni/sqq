from __future__ import annotations

"""Shared command-line banner text."""


SQQ_TITLE = "Shell  Quant  Qualifier"
SQQ_AUTHOR = "by J. PANG & Q. SUN"

SQQ_BANNER = """
┏━━━━━┳━━━━━┳━━━━━┓  +-----------------------------+
┃  ═══┃  ║  ┃  ║  ┃  |   Shell  Quant  Qualifier   |
┃═══  ┣━━━┓ ┣━━━┓ ┃  |     by J. PANG & Q. SUN     |
┗━━━━━┛   ┗━┛   ┗━┛  +-----------------------------+
""".strip()

# Reserved artwork. It is intentionally not included in the public banner.
IRSM_BANNER = """
┏━━━┳━━━┳━━━┳━━━━━┓
┗┓ ┏┫ ┏━┫ ══┃ ┃ ┃ ┃
┏┛ ┗┫ ┃ ┃══ ┃ ┃ ┃ ┃
┗━━━┻━┛ ┗━━━┻━┻━┻━┛
""".strip()

ENGINE_BADGES = {
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


def banner_for_engine(engine: str | None = None) -> str:
    """Return the base banner with an optional compatibility-preset badge."""
    badge = ENGINE_BADGES.get(str(engine or "").strip().lower())
    if badge is None:
        return SQQ_BANNER
    return "\n".join(
        f"{base}  {extra}"
        for base, extra in zip(SQQ_BANNER.splitlines(), badge.splitlines())
    )

HELP_BANNER = f"""
{SQQ_BANNER}
SQQ (Shell Quant Qualifier): Python Joint Toolkit for Water-Shell Topology Analysis.
""".strip()
