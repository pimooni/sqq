from __future__ import annotations

"""Command-line interface for SQQ."""

import argparse
import sys
from pathlib import Path

from . import __release_date__, __version__
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


VERSION_LINE = f"SQQ version: {__version__}   Release date: {__release_date__}"
ROOT_HELP_HEADER = f"{HELP_BANNER}\n\n{VERSION_LINE}"


class DescriptionFirstArgumentParser(argparse.ArgumentParser):
    """Place the root description before argparse's usage line."""

    def format_help(self) -> str:
        formatter = self._get_formatter()
        formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
        for action_group in self._action_groups:
            formatter.start_section(action_group.title)
            formatter.add_text(action_group.description)
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()
        formatter.add_text(self.epilog)
        body = formatter.format_help().lstrip()
        if self.description:
            return f"{self.description.rstrip()}\n{body}"
        return body


ANALYZE_EPILOG = """
Examples:
  sqq analyze -i test.gro -o ./result_sqq
  sqq analyze -i ./gro --pattern "*.gro" -o ./result_sqq
  sqq analyze -i "./gro/*.gro" -o ./result_sqq
  sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
  sqq analyze -i ./gro -m 00 -b hbond -w 4 --order-parameter f3,f4,q6 -o ./result_sqq
  sqq analyze -i traj.lammpstrj --top system.data -c lammps.yaml -o ./result_sqq
  sqq analyze -i md.gro --output-type info,cage-gro,summary-xlsx -o ./result_sqq
  sqq analyze -i md.gro -s 4,5,6 --cage-size H -o ./result_sqq_h
  sqq analyze -i md.gro -s 4,5,6 --find-cluster on -o ./result_sqq_cluster
  sqq analyze -m cpp -i md.gro --output-type info,cage-gro,summary-csv -o ./result_sqq_cpp

Analysis modes:
  -m 00  SQQ-Py: hbond, 4/5/6 search, 100% workers, find cluster on
  -m 50  SQQ-Py: auto graph, 4/5/6 search, 50% workers, find cluster off
  -m 99  SQQ-CPP: hbond, internal 4/5/6 search, 100% workers
  -m cpp SQQ-CPP: auto graph, internal 4/5/6 search, 50% workers

Modes 99 and cpp use the C++17 cage/F3/F4 backend without quasi-cage,
cluster, or ice analysis. Automatic 100% workers still reserve one physical
core. Explicit -w/--worker overrides the mode fraction.
-b/--bond-mode overrides the graph setting supplied by the selected mode.
--find-cluster overrides modes 00/50; SQQ-CPP rejects cluster search.

Output layout:
  grouped: frame/ring/, frame/half_cage/<type>/,
           frame/quasi_cage/<type>/, frame/cage/<type>/, frame/ice/
  flat:    all per-frame structure files in the frame folder
""".strip()

def build_parser() -> argparse.ArgumentParser:
    """Create the two-command CLI: init and analyze."""
    parser = DescriptionFirstArgumentParser(
        prog="sqq",
        description=ROOT_HELP_HEADER,
        epilog=ROOT_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION_LINE)
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Write a default config.yaml file.")
    init_parser.add_argument("-o", "--output", metavar="CONFIG.yaml", default="config.yaml", help="Output config path.")

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze MD frames.",
        epilog=ANALYZE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    analyze_parser.add_argument("-i", "--input", metavar="INPUT", required=True, help="Input file or directory (.gro/.xyz/.xtc/.trr/.dump/.lammpstrj/.dcd).")
    analyze_parser.add_argument("--pattern", metavar="PATTERN", help='Input pattern when --input is a directory; default "*.gro".')
    analyze_parser.add_argument("--top", "--topology", metavar="TOPOLOGY", dest="topology", help="GRO topology for XTC/TRR or LAMMPS DATA topology for dump/DCD input.")
    analyze_parser.add_argument("--xyz-scale", metavar="SCALE", type=float, help="Multiply XYZ coordinates by SCALE to obtain nm; default 0.1 assumes angstrom input.")
    analyze_parser.add_argument("--trajectory-stride", metavar="N", type=int, help="Read every Nth trajectory frame; default 1.")
    analyze_parser.add_argument(
        "--lammps-units",
        choices=("real", "metal", "nano"),
        help="LAMMPS physical unit style; default real.",
    )
    analyze_parser.add_argument(
        "--lammps-timestep",
        metavar="DT",
        type=float,
        help="LAMMPS integration timestep in the selected native time unit.",
    )
    analyze_parser.add_argument(
        "--lammps-atom-style",
        choices=("full", "molecular", "bond", "angle"),
        help="LAMMPS DATA atom style; default full.",
    )
    analyze_parser.add_argument("-c", "--config", metavar="CONFIG.yaml", help="YAML/JSON config file, e.g. config.yaml.")
    analyze_parser.add_argument("-o", "--output", metavar="RESULT_DIR", default="result_sqq", help="Output directory.")
    analyze_parser.add_argument(
        "-m",
        "--mode",
        choices=("00", "50", "99", "cpp"),
        help="Analysis mode: SQQ-Py 00/50 or native SQQ-CPP 99/cpp.",
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
    analyze_parser.add_argument("--quasi-size", metavar="4,5,6,7", help="Override quasi-cage base and side size lists together.")
    analyze_parser.add_argument("--quasi-base-size", metavar="4,5,6,7", help="Override quasi-cage base-ring size list.")
    analyze_parser.add_argument("--quasi-side-size", metavar="4,5,6,7", help="Override quasi-cage side-ring size list.")
    analyze_parser.add_argument("--quasi-max-layer", metavar="N", type=int, help="Override quasi_cage.max_layers; default 1 reports L1 quasi_cage and standard half_cage only.")
    analyze_parser.add_argument("--quasi-search-policy", choices=("bounded", "exact"), help="Layer growth policy: bounded preserves the established search; exact enumerates connected layer subsets.")
    analyze_parser.add_argument("--ring-definition", choices=("chordless", "shortest_path"), help="Ring definition; default chordless preserves established output.")
    analyze_parser.add_argument(
        "--order-parameter",
        metavar="NAME[,NAME...]",
        help=(
            "Select order parameters: f3, f4, qN, mcg1, mcg3, dhop35, dhop30, "
            "all, or none. Default: f3,f4."
        ),
    )
    analyze_parser.add_argument("--no-q", action="store_true", help=argparse.SUPPRESS)
    analyze_parser.add_argument("-q", "--q-degree", metavar="4,6,8,10,12", help=argparse.SUPPRESS)
    analyze_parser.add_argument("--q-neighbor-mode", choices=["graph", "cutoff", "nearest", "lammps"], help="Neighbor source for Q_l; default graph follows the active water network.")
    analyze_parser.add_argument("--q-cutoff", metavar="NM", type=float, help="Q_l neighbor cutoff in nm for cutoff/nearest/lammps modes; default 0.35.")
    analyze_parser.add_argument("--q-n-neighbor", metavar="N|NULL", help="Fixed Q_l neighbor count; lammps mode defaults to 12, NULL uses all cutoff neighbors.")
    analyze_parser.add_argument("--mcg3", choices=("on", "off"), help=argparse.SUPPRESS)
    analyze_parser.add_argument("--dhop30", choices=("on", "off"), help=argparse.SUPPRESS)
    analyze_parser.add_argument(
        "--cage-size",
        metavar="GROUP[,GROUP...]",
        help="Report cage groups I, II, H, HS-I, TS-I, or I2II; auto/all report every detected type. Default auto follows --size.",
    )
    analyze_parser.add_argument("--max-cage-face", metavar="N", type=int, help="Maximum face count searched for Euler-compatible cages; default 20.")
    analyze_parser.add_argument("--cage-fast-closure", choices=("on", "off"), help="Enable or disable indexed 2-4 half-cage fast closure; default on.")
    analyze_parser.add_argument("--cage-scientific-validation", choices=("on", "off"), help="Enable or disable optional cage face-quality/volume validation; topology validation is always on; default off.")
    analyze_parser.add_argument("--find-cluster", choices=("on", "off"), help="Override the mode preset for hydrate-cluster search; the overall default mode 50 is off.")
    analyze_parser.add_argument("--cluster-min-cage", metavar="N", type=int, help="Minimum connected cage count required for a hydrate_cluster; default 2.")
    analyze_parser.add_argument("--recursive", action="store_true", help="Read input directory recursively.")
    analyze_parser.add_argument("--pairs", metavar="PAIRS.txt", help="Pair file for bond_mode=pairs; each line contains two water ids.")
    analyze_parser.add_argument("--pair-id", metavar="KIND", choices=["resid", "oxygen_index", "atomid"], help="How ids in --pairs are interpreted; default resid.")
    analyze_parser.add_argument("--parallel-backend", choices=("process", "thread", "serial"), help="Independent-file backend; default process uses multiple CPU cores.")
    analyze_parser.add_argument("-w", "--worker", metavar="N|auto", default=None, help="Worker count or physical-core fraction; e.g. 1, 4, 0.5, 1.0, or 50%%. Reserves one physical core.")
    analyze_parser.add_argument("--workers", dest="worker", metavar="N|auto", help=argparse.SUPPRESS)
    analyze_parser.add_argument("--strict", action="store_true", help="Stop on the first failed frame.")
    analyze_parser.add_argument("--output-layout", choices=["grouped", "flat"], help="GRO layout: grouped uses ring/, half_cage/<type>/, quasi_cage/<type>/, cage/<type>/, hydrate_cluster/, and ice/; flat keeps same-folder files.")
    analyze_parser.add_argument(
        "--output-type",
        metavar="TYPE[,TYPE...]",
        help=(
            "Select outputs: info, membership-tsv, order-tsv, vmd, "
            "sqq-cage-gro, sqq-render, gro, ring-gro, half-gro, quasi-gro, "
            "cage-gro, ice-gro, cluster-gro, "
            "summary-xlsx, summary-csv, summary-detail-csv, cluster-detail, all, or none. "
            "Defaults come from the selected mode; sqq-render implies sqq-cage-gro."
        ),
    )
    analyze_parser.add_argument(
        "--cage-isomer-rows",
        choices=("nonzero", "all"),
        help="Cage-isomer rows retained in summary/detail tables; default nonzero.",
    )
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
