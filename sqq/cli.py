from __future__ import annotations

"""Command-line interface for SQQ."""

import argparse
import sys
from pathlib import Path

from . import __release_date__, __version__
from .banner import HELP_BANNER
from .config import write_default_config
from .vmd_command import VMD_COMMAND_HELP, run_vmd_command
from .pipeline import analyze


ROOT_EPILOG = """
Quick start:
  sqq init
  sqq analyze -i test.gro -o ./result_sqq
  sqq analyze -i traj.xtc -t topol.gro -c sqq_config.yaml -o ./result_sqq
  sqq vmd ./result_sqq

Use `sqq analyze -h` for analysis options and `sqq vmd -h` for VMD commands.
""".strip()


VERSION_LINE = f"SQQ version: {__version__}   Release date: {__release_date__}"
ROOT_HELP_HEADER = f"{HELP_BANNER}\n\n{VERSION_LINE}"
ENGINE_CHOICES = ("00", "py", "99", "cpp")


class DescriptionFirstArgumentParser(argparse.ArgumentParser):
    """Place the root description before argparse's usage line."""

    def parse_args(
        self,
        args: list[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        """Reject renamed options with an actionable migration error."""
        arguments = list(sys.argv[1:] if args is None else args)
        for option in arguments:
            if (
                option == "-m"
                or option.startswith("-m=")
                or option in {f"-m{engine}" for engine in ENGINE_CHOICES}
            ):
                self.error(
                    "-m has been replaced by -e.\nUse: -e py"
                )
            if option == "--mode" or option.startswith("--mode="):
                self.error(
                    "--mode has been replaced by --engine.\nUse: --engine py"
                )
            if option == "--pairs" or option.startswith("--pairs="):
                self.error(
                    "--pairs has been replaced by --pair.\nUse: --pair PAIRS.txt"
                )
        return super().parse_args(arguments, namespace)

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


# Configuration-only fields retain neutral Namespace defaults so shared runtime
# normalization stays independent of the public command-line surface.
RUNTIME_COMPATIBILITY_DEFAULTS = {
    "pattern": None,
    "xyz_scale": None,
    "lammps_units": None,
    "lammps_timestep": None,
    "lammps_atom_style": None,
    "ring_size": None,
    "quasi_size": None,
    "quasi_base_size": None,
    "quasi_side_size": None,
    "quasi_max_layer": None,
    "quasi_search_policy": None,
    "ring_definition": None,
    "no_q": False,
    "q_degree": None,
    "q_neighbor_mode": None,
    "q_cutoff": None,
    "q_n_neighbor": None,
    "mcg3": None,
    "dhop30": None,
    "cage_size": None,
    "max_cage_face": None,
    "cage_fast_closure": None,
    "cage_scientific_validation": None,
    "cluster_min_cage": None,
    "pair_id": None,
    "parallel_backend": None,
    "recursive": False,
    "strict": False,
    "output_layout": None,
    "output_type": None,
    "cage_isomer_rows": None,
}


ANALYZE_EPILOG = """
Examples:
  sqq analyze -i test.gro -o ./result_sqq
  sqq analyze -i traj.xtc --top topol.gro -c sqq_config.yaml -o ./result_sqq
  sqq analyze -i ./gro -e py -b hbond -w 4 --order-parameter f3,f4,q6
  sqq analyze -i traj.lammpstrj -t system.data -o ./result_sqq
  sqq analyze -i md.gro -b pairs --pair water_pairs.txt
  sqq analyze -e cpp -i md.gro -s 4,5,6 -o ./result_sqq_cpp

Accepted engine values are 00, py, 99, and cpp. The default is py.
""".strip()



def _add_analysis_arguments(
    command_parser: argparse.ArgumentParser,
    *,
    input_required: bool,
) -> None:
    """Add the intentionally small public analysis option surface."""
    command_parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        required=input_required,
        help=(
            "Input file or directory "
            "(.gro/.xyz/.xtc/.trr/.dump/.lammpstrj/.dcd)."
        ),
    )
    command_parser.add_argument(
        "-t",
        "--top",
        dest="topology",
        metavar="TOPOLOGY",
        help="GRO or LAMMPS DATA topology for trajectory input.",
    )
    command_parser.add_argument(
        "-c",
        "--config",
        metavar="CONFIG.yaml",
        help="YAML/JSON configuration file, e.g. sqq_config.yaml.",
    )
    command_parser.add_argument(
        "-o",
        "--output",
        metavar="RESULT_DIR",
        default="result_sqq",
        help="Output directory; default result_sqq.",
    )
    command_parser.add_argument(
        "-e",
        "--engine",
        choices=ENGINE_CHOICES,
        help="Analysis engine: py or cpp; compatibility presets 00/99 are accepted.",
    )
    command_parser.add_argument(
        "-w",
        "--worker",
        metavar="N|auto",
        default=None,
        help="Worker count or physical-core fraction, e.g. 4, 0.5, or 50%%.",
    )
    command_parser.add_argument(
        "-dt",
        "--delta-time",
        metavar="PS",
        type=float,
        help="Physical sampling interval in ps; default all stored frames.",
    )
    command_parser.add_argument(
        "-b",
        "--bond-mode",
        choices=("auto", "hbond", "oo", "pairs"),
        help="Water-graph connection mode.",
    )
    command_parser.add_argument(
        "-s",
        "--size",
        metavar="4,5,6,7",
        help="Ring and quasi-cage search sizes.",
    )
    command_parser.add_argument(
        "--find-half",
        choices=("on", "off"),
        help="Enable or disable standard half-cage search.",
    )
    command_parser.add_argument(
        "--find-quasi",
        choices=("on", "off"),
        help="Enable or disable layered quasi-cage search.",
    )
    command_parser.add_argument(
        "--find-cluster",
        choices=("on", "off"),
        help="Enable or disable hydrate-cluster search.",
    )
    command_parser.add_argument(
        "--order-parameter",
        metavar="NAME[,NAME...]",
        help="Select order parameters, e.g. f3,f4,q6, all, or none.",
    )
    command_parser.add_argument(
        "--output-type",
        metavar="TYPE[,TYPE...]",
        help=(
            "Replace the output list; use default with extra types to extend "
            "the engine defaults."
        ),
    )
    command_parser.add_argument(
        "--pair",
        dest="pair",
        metavar="PAIRS.txt",
        help="Pair file used with --bond-mode pairs.",
    )
    command_parser.set_defaults(**RUNTIME_COMPATIBILITY_DEFAULTS)


def build_parser() -> argparse.ArgumentParser:
    """Create the SQQ command-line interface."""
    parser = DescriptionFirstArgumentParser(
        prog="sqq",
        description=ROOT_HELP_HEADER,
        epilog=ROOT_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION_LINE)
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser(
        "init",
        help="Write a default sqq_config.yaml file.",
    )
    init_parser.add_argument(
        "-o",
        "--output",
        metavar="CONFIG.yaml",
        default="sqq_config.yaml",
        help="Output config path; default sqq_config.yaml.",
    )

    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze MD frames.",
        epilog=ANALYZE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_analysis_arguments(analyze_parser, input_required=True)

    vmd_parser = subparsers.add_parser(
        "vmd",
        help="Locate and validate SQQ VMD render packages.",
        epilog=VMD_COMMAND_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    vmd_parser.add_argument(
        "path",
        nargs="?",
        default=".",
        metavar="PATH",
        help="Result directory, render directory, or SQQ .vmd.tcl file.",
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
        try:
            write_default_config(out)
        except FileExistsError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 2
        print(f"Wrote default SQQ config: {out}")
        return 0
    if args.command == "analyze":
        analyze(args)
        return 0
    if args.command == "vmd":
        return run_vmd_command(args.path)
    raise AssertionError(f"Unhandled command: {args.command}")
