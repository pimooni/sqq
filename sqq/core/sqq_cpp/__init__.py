"""SQQ-CPP backend contract.

The native extension and its Python adapter are exposed from this package.
The former :mod:`sqq.core.cpp_backend` path remains a compatibility facade.
"""

from .backend import analyze_frame, analyze_frame_cpp, native_available

__all__ = ["analyze_frame", "analyze_frame_cpp", "native_available"]
