"""Stable exception hierarchy for the public SQQ Python API."""


class SQQError(Exception):
    """Base class for recoverable SQQ API failures."""


class ConfigurationError(SQQError, ValueError):
    """A configuration source or value cannot be resolved safely."""


class InputError(SQQError, ValueError):
    """Coordinate, topology, or trajectory input is invalid."""


class AnalysisError(SQQError, RuntimeError):
    """A valid frame could not be analyzed."""


__all__ = ["AnalysisError", "ConfigurationError", "InputError", "SQQError"]
