from __future__ import annotations

"""Command-line interface for SQQ."""

import argparse
import sys
from pathlib import Path

from .banner import HELP_BANNER
from .config import write_default_config
from .pipeline import analyze


ROOT_EPILOG = """
Quick start:
  sqq init -o config.yaml
  sqq analyze -i test.gro -o ./result_sqq
  sqq analyze -i ./gro --pattern "*.gro" -o ./result_sqq

Use `sqq analyze -h` for analysis options and examples.
""".strip()


ANALYZE_EPILOG = """
Examples:
  sqq analyze -i test.gro -o ./result_sqq
  sqq analyze -i ./gro --pattern "*.gro" -o ./result_sqq
  sqq analyze -i "./gro/*.gro" -o ./result_sqq
  sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
  sqq analyze -i ./gro -m 00 -b hbond --workers 4 -o ./result_sqq

Analysis modes:
  -m 00  Rigorous: hbond, 4/5/6 rings and cages, other cages, 25% CPU workers
  -m 50  Standard: auto graph, 5/6 rings and cages, standard cages, 50% CPU workers
  -m 99  Performance: O-O graph, 5/6 rings and cages, standard cages, 90% CPU workers

Modes do not change quasi_cage.max_layers; its default remains 1.
-b/--bond-mode overrides the graph setting supplied by the selected mode.

Output layout:
  grouped: frame/ring/, frame/half_cage/<type>/,
           frame/quasi_cage/<type>/, frame/cage/<type>/, frame/ice/
  flat:    all per-frame structure files in the frame folder
""".strip()

def build_parser() -> argparse.ArgumentParser:
    """Create the two-command CLI: init and analyze."""
    parser = argparse.ArgumentParser(
        prog="sqq",
        description=HELP_BANNER,
        epilog=ROOT_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Write a default config.yaml file.")
    init_parser.add_argument("-o", "--output", metavar="CONFIG.yaml", default="config.yaml", help="Output config path.")

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze MD frames.",
        epilog=ANALYZE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    analyze_parser.add_argument("-i", "--input", metavar="INPUT", required=True, help="Input file or directory (.gro/.xyz/.xtc/.trr).")
    analyze_parser.add_argument("--pattern", metavar="PATTERN", help='Input pattern when --input is a directory, e.g. "*.gro".')
    analyze_parser.add_argument("--top", "--topology", metavar="TOPOLOGY.gro", dest="topology", help="Topology/structure file for .xtc/.trr input.")
    analyze_parser.add_argument("-c", "--config", metavar="CONFIG.yaml", help="YAML/JSON config file, e.g. config.yaml.")
    analyze_parser.add_argument("-o", "--output", metavar="RESULT_DIR", default="result_sqq", help="Output directory.")
    analyze_parser.add_argument(
        "-m",
        "--mode",
        choices=("00", "50", "99"),
        help="Analysis preset: 00 rigorous, 50 standard (default), or 99 performance.",
    )
    analyze_parser.add_argument(
        "-b",
        "--bond-mode",
        choices=("auto", "hbond", "oo", "pairs"),
        help="Override the water-graph connection mode selected by the analysis preset.",
    )
    analyze_parser.add_argument(
        "-s",
        "--sizes",
        metavar="4,5,6",
        help="Set ring, quasi_cage base/side, and cage ring sizes together; including 4 enables other cages unless --no-other-cages is used.",
    )
    analyze_parser.add_argument("--ring-sizes", metavar="4,5,6", help="Override ring.sizes.")
    analyze_parser.add_argument("--quasi-sizes", metavar="4,5,6", help="Override quasi_cage base_sizes and side_sizes together.")
    analyze_parser.add_argument("--quasi-base-sizes", metavar="4,5,6", help="Override quasi_cage.base_sizes.")
    analyze_parser.add_argument("--quasi-side-sizes", metavar="4,5,6", help="Override quasi_cage.side_sizes.")
    analyze_parser.add_argument("--quasi-max-layers", metavar="N", type=int, help="Override quasi_cage.max_layers; default 1 reports L1 quasi_cage and standard half_cage only.")
    analyze_parser.add_argument(
        "--cage-sizes",
        dest="cage_sizes",
        metavar="4,5,6",
        help="Override cage.ring_sizes; including 4 enables other cages unless --no-other-cages is used.",
    )
    other_group = analyze_parser.add_mutually_exclusive_group()
    other_group.add_argument("--other-cages", action="store_true", help="Enable generated unconventional 4/5/6 cage targets.")
    other_group.add_argument("--no-other-cages", action="store_true", help="Disable generated unconventional cage targets.")
    analyze_parser.add_argument("--other-max-faces", metavar="N", type=int, help="Maximum face count for generated unconventional cages.")
    analyze_parser.add_argument("--recursive", action="store_true", help="Read input directory recursively.")
    analyze_parser.add_argument("--pairs", metavar="PAIRS.txt", help="Pair file for bond_mode=pairs; each line contains two water ids.")
    analyze_parser.add_argument("--pair-id", metavar="KIND", choices=["resid", "oxygen_index", "atomid"], help="How ids in --pairs are interpreted.")
    analyze_parser.add_argument("--workers", metavar="N", default=None, help="Explicit frame-level worker count for independent GRO/XYZ files; overrides the mode CPU percentage.")
    analyze_parser.add_argument("--strict", action="store_true", help="Stop on the first failed frame.")
    analyze_parser.add_argument("--output-layout", choices=["grouped", "flat"], help="GRO layout: grouped uses ring/, half_cage/<type>/, quasi_cage/<type>/, cage/<type>/, and ice/; flat keeps same-folder files.")
    analyze_parser.add_argument("--no-info", action="store_true", help="Disable per-frame *_info.md output.")
    analyze_parser.add_argument("--no-gro", action="store_true", help="Disable GRO structure output.")
    analyze_parser.add_argument("--no-ring-gro", action="store_true", help="Disable ring GRO files.")
    analyze_parser.add_argument("--no-half-cage-gro", action="store_true", help="Disable half_cage GRO files.")
    analyze_parser.add_argument("--no-quasi-cage-gro", action="store_true", help="Disable quasi_cage GRO files.")
    analyze_parser.add_argument("--no-cage-gro", action="store_true", help="Disable cage GRO files.")
    analyze_parser.add_argument("--no-ice-gro", action="store_true", help="Disable ice GRO files.")
    analyze_parser.add_argument("--no-xlsx", action="store_true", help="Disable summary.xlsx output.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch the selected SQQ subcommand."""
    parser = build_parser()
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    if args.command == "init":
        out = Path(args.output)
        write_default_config(out)
        print(f"Wrote default SQQ config: {out}")
        return 0
    if args.command == "analyze":
        analyze(args)
        return 0
    raise AssertionError(f"Unhandled command: {args.command}")


