"""SQQ: Shell Quant Qualifier."""

from ._version import __release_date__, __version__
from .api import analyze_frame, load_config, read_frames
from .exceptions import AnalysisError, ConfigurationError, InputError, SQQError


__all__ = [
    "AnalysisError",
    "ConfigurationError",
    "InputError",
    "SQQError",
    "analyze_frame",
    "load_config",
    "read_frames",
    "__release_date__",
    "__version__",
]
