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
  sqq analyze -i md.gro -s 4,5,6 --cage-size H -o ./result_sqq_h

Analysis modes:
  -m 00  Rigorous: hbond, 4/5/6 search, 25% CPU workers
  -m 50  Standard: auto graph, 5/6 search, 50% CPU workers
  -m 99  Performance: O-O graph, 5/6 search, 90% CPU workers

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
    analyze_parser.add_argument("--pattern", metavar="PATTERN", help='Input pattern when --input is a directory; default "*.gro".')
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
        "--size",
        metavar="4,5,6,7",
        help="Set ring/quasi-cage search sizes; cage detection uses the selected 4/5/6 sizes.",
    )
    analyze_parser.add_argument("--ring-size", metavar="4,5,6,7", help="Report a subset of searched ring sizes; default auto follows --size.")
    analyze_parser.add_argument("--quasi-sizes", metavar="4,5,6,7", help="Override quasi-cage base and side sizes together.")
    analyze_parser.add_argument("--quasi-base-sizes", metavar="4,5,6,7", help="Override quasi-cage base-ring sizes.")
    analyze_parser.add_argument("--quasi-side-sizes", metavar="4,5,6,7", help="Override quasi-cage side-ring sizes.")
    analyze_parser.add_argument("--quasi-max-layers", metavar="N", type=int, help="Override quasi_cage.max_layers; default 1 reports L1 quasi_cage and standard half_cage only.")
    analyze_parser.add_argument("--no-q", action="store_true", help="Disable Steinhardt Q_l order-parameter calculation.")
    analyze_parser.add_argument("-q", "--q-degree", metavar="4,6,8,10,12", help="Comma-separated Q_l degree list to report; default 6,12.")
    analyze_parser.add_argument("--q-neighbor-mode", choices=["graph", "cutoff", "nearest", "lammps"], help="Neighbor source for Q_l; default graph follows the active water network.")
    analyze_parser.add_argument("--q-cutoff", metavar="NM", type=float, help="Q_l neighbor cutoff in nm for cutoff/nearest/lammps modes; default 0.35.")
    analyze_parser.add_argument("--q-n-neighbor", metavar="N|NULL", help="Fixed Q_l neighbor count; lammps mode defaults to 12, NULL uses all cutoff neighbors.")
    analyze_parser.add_argument(
        "--cage-size",
        metavar="GROUP[,GROUP...]",
        help="Report cage groups I, II, H, HS-I, TS-I, or I2II; auto/all report every detected type. Default auto follows --size.",
    )
    analyze_parser.add_argument("--max-cage-faces", metavar="N", type=int, help="Maximum face count searched for Euler-compatible cages; default 20.")
    analyze_parser.add_argument("--recursive", action="store_true", help="Read input directory recursively.")
    analyze_parser.add_argument("--pairs", metavar="PAIRS.txt", help="Pair file for bond_mode=pairs; each line contains two water ids.")
    analyze_parser.add_argument("--pair-id", metavar="KIND", choices=["resid", "oxygen_index", "atomid"], help="How ids in --pairs are interpreted; default resid.")
    analyze_parser.add_argument("--workers", metavar="N|auto", default=None, help="Frame-level worker count for independent GRO/XYZ files; overrides the mode CPU percentage.")
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
    analyze_parser.add_argument("--write-order-tsv", action="store_true", help="Write per-water *_order_parameter.tsv files.")
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


