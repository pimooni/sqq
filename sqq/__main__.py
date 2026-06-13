"""Run SQQ with `python -m sqq`."""

from .cli import main


if __name__ == "__main__":
    # Delegate to the same entry point used by the installed console script.
    raise SystemExit(main())
