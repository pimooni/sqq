from __future__ import annotations

"""Command-line interface for SQQ."""

import argparse
from pathlib import Path

from .config import write_default_config
from .pipeline import analyze


HELP_BANNER = """
+----------------------------+
|   Shell  Quant  Qualifier  |
+----------------------------+

SQQ for MD water-shell topology analysis.
""".strip()


def build_parser() -> argparse.ArgumentParser:
    """Create the two-command CLI: init and analyze."""
    parser = argparse.ArgumentParser(
        prog="sqq",
        description=HELP_BANNER,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Write a default config.yaml file.")
    init_parser.add_argument("-o", "--output", default="config.yaml", help="Output config path.")

    analyze_parser = subparsers.add_parser("analyze", help="Analyze MD frames.")
    analyze_parser.add_argument("-i", "--input", required=True, help="Input file or directory.")
    analyze_parser.add_argument("--pattern", help="Input pattern when --input is a directory.")
    analyze_parser.add_argument("--top", "--topology", dest="topology", help="Topology file for xtc/trr.")
    analyze_parser.add_argument("-c", "--config", help="YAML config file.")
    analyze_parser.add_argument("-o", "--output", default="result_sqq", help="Output directory.")
    analyze_parser.add_argument("--recursive", action="store_true", help="Read input directory recursively.")
    analyze_parser.add_argument("--pairs", help="Pair file for bond_mode=pairs; each line contains two water ids.")
    analyze_parser.add_argument("--pair-id", choices=["resid", "oxygen_index", "atomid"], help="How ids in --pairs are interpreted.")
    analyze_parser.add_argument("--n-jobs", default=None, help="Frame-level worker count for independent GRO/XYZ files; default auto runs serially.")
    analyze_parser.add_argument("--strict", action="store_true", help="Stop on the first failed frame.")
    analyze_parser.add_argument("--no-gro", action="store_true", help="Disable GRO structure output.")
    analyze_parser.add_argument("--no-xlsx", action="store_true", help="Disable summary.xlsx output.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch the selected SQQ subcommand."""
    args = build_parser().parse_args(argv)
    if args.command == "init":
        out = Path(args.output)
        write_default_config(out)
        print(f"Wrote default SQQ config: {out}")
        return 0
    if args.command == "analyze":
        analyze(args)
        return 0
    raise AssertionError(f"Unhandled command: {args.command}")


