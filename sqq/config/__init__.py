"""Stable public configuration API."""

from .capabilities import normalize_engine_capabilities
from .defaults import (
    ALL_ORDER_PARAMETERS, ALL_OUTPUT_TYPES, CONFIG_SCHEMA_VERSION,
    CPP_ALL_OUTPUT_TYPES, CPP_DEFAULT_OUTPUT_TYPES, CPP_MODE, CPP_MODES,
    CPP_OUTPUT_TYPES, DEFAULT_CONFIG, DEFAULT_MODE, DEFAULT_ORDER_PARAMETERS,
    DEFAULT_OUTPUT_TYPES, GRO_OUTPUT_TYPES, MODE_PRESETS, ORDER_PARAMETER_ALIASES,
    ORDER_PARAMETER_CHOICES, OUTPUT_TYPE_ORDER, default_config_copy,
)
from .loading import load_config, refresh_resolution_report, resolve_config
from .migrate import (
    legacy_enabled, migrate_legacy_order_parameters, migrate_yaml_keys,
    strip_legacy_selection_keys,
)
from .model import ResolvedConfig, ResolutionAdjustment, ResolutionReport
from .overrides import apply_cli_overrides, merge_config
from .presets import (
    apply_mode_preset, engine_display, is_cpp_mode, mode_display, mode_label,
    mode_worker_count, mode_worker_fraction, normalize_mode, profile_name,
)
from .resolve import GraphModeResolution, resolve_graph_mode
from .scopes import (
    finite_float, nonnegative_integer, normalize_analysis_scopes,
    parse_auto_on_off, parse_on_off, positive_integer,
    resolve_cage_report_types, resolve_size_list, string_list,
)
from .serialize import canonical_config, default_config_template, dump_config, write_default_config
from .validation import (
    normalize_cpp_order_parameters, normalize_cpp_output_types,
    normalize_order_parameters, normalize_output_types, normalize_parallel_backend,
    order_parameter_display,
    order_parameter_sort_key, output_enabled, output_type_display,
    q_degrees_from_order_parameters, validate_cpp_cli, validate_cpp_config,
    validate_user_config_keys,
)

__all__ = [
    "ALL_ORDER_PARAMETERS", "ALL_OUTPUT_TYPES", "CONFIG_SCHEMA_VERSION",
    "CPP_ALL_OUTPUT_TYPES", "CPP_DEFAULT_OUTPUT_TYPES", "CPP_MODE", "CPP_MODES",
    "CPP_OUTPUT_TYPES", "DEFAULT_CONFIG", "DEFAULT_MODE", "DEFAULT_ORDER_PARAMETERS",
    "DEFAULT_OUTPUT_TYPES", "GRO_OUTPUT_TYPES", "GraphModeResolution", "MODE_PRESETS",
    "ORDER_PARAMETER_ALIASES", "ORDER_PARAMETER_CHOICES", "OUTPUT_TYPE_ORDER",
    "ResolvedConfig", "ResolutionAdjustment", "ResolutionReport", "apply_cli_overrides", "apply_mode_preset",
    "canonical_config", "default_config_copy", "default_config_template", "dump_config",
    "engine_display", "finite_float", "is_cpp_mode", "legacy_enabled", "load_config", "merge_config",
    "migrate_legacy_order_parameters", "migrate_yaml_keys", "mode_display", "mode_label",
    "mode_worker_count", "mode_worker_fraction", "normalize_cpp_order_parameters",
    "normalize_cpp_output_types", "normalize_engine_capabilities", "normalize_mode",
    "normalize_order_parameters", "normalize_output_types", "normalize_parallel_backend",
    "order_parameter_display",
    "nonnegative_integer", "normalize_analysis_scopes", "order_parameter_sort_key",
    "output_enabled", "output_type_display", "parse_auto_on_off", "parse_on_off",
    "positive_integer", "profile_name", "refresh_resolution_report",
    "q_degrees_from_order_parameters", "resolve_config", "resolve_graph_mode", "strip_legacy_selection_keys",
    "resolve_cage_report_types", "resolve_size_list", "string_list", "validate_cpp_cli",
    "validate_cpp_config",
    "validate_user_config_keys", "write_default_config",
]
