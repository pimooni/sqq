"""SQQ: Shell Quant Qualifier."""

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

__version__ = "0.5.2"
__release_date__ = "Aug 10, 2026"
