"""Public adapter for the SQQ-CPP cage-analysis engine."""

from ._backend_impl import _load_native_module, analyze_frame_cpp


def native_available() -> bool:
    """Return whether the matching native extension can be loaded."""
    try:
        _load_native_module()
    except RuntimeError:
        return False
    return True


analyze_frame = analyze_frame_cpp


__all__ = ["analyze_frame", "analyze_frame_cpp", "native_available"]
