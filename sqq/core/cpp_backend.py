"""Compatibility import for the SQQ-CPP backend."""

from .sqq_cpp.backend import analyze_frame_cpp, native_available

__all__ = ["analyze_frame_cpp", "native_available"]
