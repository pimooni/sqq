"""Stable exception hierarchy for the public SQQ Python API."""


class SQQError(Exception):
    """Base class for recoverable SQQ API failures."""


class ConfigurationError(SQQError, ValueError):
    """A configuration source or value cannot be resolved safely."""


class InputError(SQQError, ValueError):
    """Coordinate, topology, or trajectory input is invalid."""


class AnalysisError(SQQError, RuntimeError):
    """A valid frame could not be analyzed."""


class OutputLockError(SQQError, RuntimeError):
    """Another process currently owns the selected output directory."""


__all__ = [
    "AnalysisError",
    "ConfigurationError",
    "InputError",
    "OutputLockError",
    "SQQError",
]
