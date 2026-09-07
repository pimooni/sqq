"""SQQ-CPP backend contract.

The native extension and its Python adapter are exposed from this package.
"""

from .backend import analyze_frame, analyze_frame_cpp, native_available, require_native

__all__ = ["analyze_frame", "analyze_frame_cpp", "native_available", "require_native"]
