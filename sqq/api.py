"""Stable, typed Python entry points for SQQ configuration and analysis."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

from .config import normalize_analysis_scopes, refresh_resolution_report, validate_cpp_config
from .config.model import (
    ResolvedConfig,
    ResolutionAdjustment,
    ResolutionReport,
)
from .exceptions import AnalysisError, ConfigurationError, InputError, SQQError
from .models import Frame, FrameResult


ConfigSource = str | Path | Mapping[str, Any] | None


def load_config(
    source: ConfigSource = None,
    *,
    engine: str | None = None,
) -> ResolvedConfig:
    """Load, migrate, merge, normalize, and validate an SQQ configuration.

    ``source`` may be a YAML path, a partial mapping, or ``None`` for the
    documented defaults.  The returned mapping is detached from the caller's
    data and includes an auditable resolution report.
    """
    try:
        from .config.loading import resolve_config

        data = resolve_config(source, mode=engine)
        validate_cpp_config(data)
        normalize_analysis_scopes(data)
        refresh_resolution_report(data)
        return _resolved_config(data)
    except ConfigurationError:
        raise
    except (KeyError, OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"Unable to resolve SQQ configuration: {exc}") from exc


def read_frames(
    source: str | Path | Iterable[str | Path],
    *,
    topology: str | Path | None = None,
    config: ConfigSource | ResolvedConfig = None,
    frame_indexes: Sequence[int] | None = None,
) -> Iterator[Frame]:
    """Yield supported coordinate or trajectory frames with stable errors."""
    try:
        resolved = _coerce_config(config)
        input_config = resolved["input"]
        paths = _input_paths(
            source,
            pattern=str(input_config.get("pattern", "*.gro")),
            recursive=bool(input_config.get("recursive", False)),
        )
        topology_path = (
            None if topology is None else Path(topology).expanduser()
        )
        from .io.trajectory import read_frames as read_input_frames

        yield from read_input_frames(
            paths,
            topology=topology_path,
            xyz_scale=float(input_config.get("xyz_scale", 0.1)),
            frame_indexes=frame_indexes,
            lammps_config=input_config.get("lammps", {}),
        )
    except SQQError:
        raise
    except (KeyError, OSError, TypeError, ValueError, RuntimeError) as exc:
        raise InputError(f"Unable to read SQQ input frames: {exc}") from exc


def analyze_frame(
    frame: Frame,
    config: ConfigSource | ResolvedConfig = None,
    *,
    engine: str | None = None,
    stage_callback: Callable[[str], None] | None = None,
) -> FrameResult:
    """Analyze one official :class:`Frame` and return an official FrameResult."""
    if not isinstance(frame, Frame):
        raise InputError("analyze_frame requires an sqq.models.Frame instance.")
    try:
        resolved = _coerce_config(config, engine=engine)
        mutable_config = resolved.to_mutable_dict()
        from .runtime.frame import analyze_frame as analyze_input_frame

        result = analyze_input_frame(
            frame,
            mutable_config,
            stage_callback=stage_callback,
        )
    except (ConfigurationError, InputError):
        raise
    except (KeyError, OSError, TypeError, ValueError, RuntimeError) as exc:
        raise AnalysisError(f"Unable to analyze SQQ frame {frame.name!r}: {exc}") from exc
    if not isinstance(result, FrameResult):
        raise AnalysisError("The selected SQQ engine returned an invalid frame result.")
    return result


def _coerce_config(
    source: ConfigSource | ResolvedConfig,
    *,
    engine: str | None = None,
) -> ResolvedConfig:
    if isinstance(source, ResolvedConfig):
        if engine is not None:
            selected = str(source.report.requested_selector)
            if str(engine).strip().lower() != selected.lower():
                raise ConfigurationError(
                    "engine cannot override an already resolved configuration."
                )
        return source
    return load_config(source, engine=engine)


def _resolved_config(data: Mapping[str, Any]) -> ResolvedConfig:
    detached = deepcopy(dict(data))
    raw_report = detached.get("resolution_report", {})
    if not isinstance(raw_report, Mapping):
        raise ConfigurationError("Configuration resolution report is invalid.")
    raw_adjustments = raw_report.get("adjustments", ())
    if not isinstance(raw_adjustments, Iterable) or isinstance(
        raw_adjustments, (str, bytes, Mapping)
    ):
        raise ConfigurationError("Configuration adjustments must be a sequence.")
    adjustments: list[ResolutionAdjustment] = []
    for item in raw_adjustments:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Configuration adjustment is not a mapping.")
        adjustments.append(
            ResolutionAdjustment(
                parameter=str(item.get("parameter", "unknown")),
                requested=item.get("requested"),
                effective=item.get("effective"),
                reason=str(item.get("reason", "configuration resolution")),
            )
        )
    report = ResolutionReport(
        requested_selector=str(
            raw_report.get("requested_selector", detached.get("mode", "py"))
        ),
        engine=str(raw_report.get("engine", detached.get("run", {}).get("engine", ""))),
        profile=str(
            raw_report.get("profile", detached.get("run", {}).get("profile", ""))
        ),
        adjustments=tuple(adjustments),
    )
    return ResolvedConfig(detached, report)


def _input_paths(
    source: str | Path | Iterable[str | Path],
    *,
    pattern: str,
    recursive: bool,
) -> list[Path]:
    if isinstance(source, (str, Path)):
        raw = [Path(source).expanduser()]
    else:
        try:
            raw = [Path(item).expanduser() for item in source]
        except TypeError as exc:
            raise InputError("Input source must be a path or iterable of paths.") from exc
    if not raw:
        raise InputError("Input source contains no paths.")
    from .io.trajectory import expand_inputs

    if len(raw) == 1:
        return expand_inputs(raw[0], pattern, recursive)
    paths: list[Path] = []
    for path in raw:
        paths.extend(expand_inputs(path, pattern, recursive))
    return paths


__all__ = ["analyze_frame", "load_config", "read_frames"]
