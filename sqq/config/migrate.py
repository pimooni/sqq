"""Backward-compatible configuration-key migration."""

from .._config import migrate_yaml_keys, validate_user_config_keys

__all__ = ["migrate_yaml_keys", "validate_user_config_keys"]
