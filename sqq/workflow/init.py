"""Configuration-initialization workflow."""

from pathlib import Path

from ..config import write_default_config


def initialize_config(path: str | Path = "sqq_config.yaml") -> Path:
    """Write one default configuration and return its path."""
    output = Path(path)
    write_default_config(output)
    return output


__all__ = ["initialize_config"]
