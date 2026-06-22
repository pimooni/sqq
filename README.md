# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

SQQ builds a water network, reports coordination diagnostics, and finds rings, standard half-cages, quasi-cages, closed cages, cage guest occupancy, F3/F4 metrics, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

## Install

Install the released package from PyPI:

```bash
pip install sqq
```

Upgrade an existing installation:

```bash
pip install -U sqq
```

For local development from a source checkout:

```bash
pip install -e .
```

Then use:

```bash
sqq -h
sqq init -o config.yaml
sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

During source-tree development without installation:

```bash
python -m sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

## Quick Start

Single GRO file:

```bash
sqq analyze -i test1.gro -o ./result_sqq
```

Directory of GRO files:

```bash
sqq analyze -i ./gro --pattern "*.gro" -o ./result_sqq
```

Glob pattern:

```bash
sqq analyze -i "./gro/*.gro" -o ./result_sqq
```

XTC/TRR trajectory with a topology file:

```bash
sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
```

## Analysis Modes

`-m` / `--mode` selects one of three base presets:

| Mode | Purpose | Water graph | Search sizes | Automatic workers |
| --- | --- | --- | --- | --- |
| `00` | Rigorous | Hydrogen bond | 4, 5, 6 | 25% of logical CPUs |
| `50` | Standard (default) | Auto | 5, 6 | 50% of logical CPUs |
| `99` | Performance screening | O-O connectivity | 5, 6 | 90% of logical CPUs |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m 50 -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_performance
```

Modes do not change `quasi_cage.max_layers`; L1 remains the default in every mode. Use `--quasi-max-layers` explicitly for L2/L3. Automatic workers are capped by the number of independent GRO/XYZ files. A single coordinate file or XTC/TRR trajectory still runs with one worker. `--workers N` overrides the mode percentage. Parallel runs display live aggregate stages and up to six active files with per-stage and per-file timings.

## Common Commands

Write a default configuration file:

```bash
sqq init -o config.yaml
```

Search 4/5/6 ring faces and report Type H cages:

```bash
sqq analyze -i md.gro -s 4,5,6 --cage-size 512,51268,435663 -o ./result_sqq_456
```

Report every detected cage composition in the selected search scope:

```bash
sqq analyze -i md.gro -s 4,5,6 --cage-size all -o ./result_sqq_all_cages
```

Enable outer quasi-cage layers:

```bash
sqq analyze -i md.gro --quasi-max-layers 3 -o ./result_sqq_l3
```

Disable per-frame info or structure GRO output:

```bash
sqq analyze -i md.gro --no-info -o ./result_sqq_no_info
sqq analyze -i md.gro --no-gro -o ./result_sqq_report
```

Parallelize independent GRO/XYZ files:

```bash
sqq analyze -i ./gro --pattern "*.gro" --workers 4 -o ./result_sqq
```

## Important Defaults

```yaml
mode: "50"

graph:
  bond_mode: auto
  oo_cutoff_nm: 0.35
  hbond_distance_nm: 0.35
  hbond_angle_deg: 30.0

ring:
  sizes: [5, 6]
  report_sizes: auto
  chordless: true

quasi_cage:
  enabled: true
  base_sizes: auto
  side_sizes: auto
  max_layers: 1

cage:
  enabled: true
  report_types: [512, 51262, 51263, 51264]
  max_faces: 20
  search_mode: grow
  seed_mode: ring
  occupancy_mode: polyhedron

output:
  write_info: true
  write_gro: true
  write_xlsx_summary: true
  structure_layout: grouped

parallel:
  workers: auto
```

Configuration priority:

```text
built-in defaults < mode preset < config.yaml < explicit command-line options
```

## Search and Report Scope

`-s` / `--size` defines the ring-face sizes used during detection. Ring and cage reporting are filtered afterward:

```bash
# Search 4/5/6, report only ring 5/6 and the two Type H cages
sqq analyze -i md.gro -s 4,5,6 --ring-size 5,6 --cage-size 51268,435663
```

All detected cages still participate in half-cage, quasi-cage, and free-ring filtering. `--cage-size` changes user-facing counts and files, not topology ownership. Cage detection supports 4/5/6 faces; ring and quasi-cage detection also support size 7.

Named cage types are:

```text
512, 51262, 51263, 51264, 51268, 435663
```

`51268` and `435663` are the two named Type H cages. Use `--cage-size all` to report every detected Euler-compatible cage composition up to `--max-cage-faces`.

## Useful Options

| Option | Possible values | Meaning |
| --- | --- | --- |
| `-i, --input INPUT` | GRO/XYZ/XTC/TRR file, directory, or glob | Input coordinate file or trajectory source |
| `-c, --config FILE` | YAML or JSON file | User configuration file |
| `-o, --output DIR` | Directory path | Output directory |
| `-m, --mode MODE` | `00`, `50`, `99` | Select rigorous, standard, or performance preset |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --size SIZES` | Comma-separated subset of `4,5,6,7` | Define the ring-face sizes searched |
| `--ring-size SIZES` | Subset of searched sizes | Report only these ring sizes |
| `--cage-size TYPES` | Named/generic types or `all` | Report exact cage types after detection |
| `--max-cage-faces N` | Positive integer; default `20` | Limit generated cage search compositions |
| `--pattern PATTERN` | Glob such as `*.gro` | File pattern for directory input |
| `--top FILE.gro` | GRO topology file | Topology/structure file for XTC/TRR input |
| `--recursive` | Flag | Search input directories recursively |
| `--quasi-sizes SIZES` | Comma-separated subset of `4,5,6,7` | Override quasi-cage base and side sizes together |
| `--quasi-base-sizes SIZES` | Comma-separated subset of `4,5,6,7` | Override quasi-cage base-ring sizes |
| `--quasi-side-sizes SIZES` | Comma-separated subset of `4,5,6,7` | Override quasi-cage side-ring sizes |
| `--quasi-max-layers N` | Positive integer | Report quasi-cage layers up to N; default is 1 |
| `--pairs FILE` | Text pair-map file | Supply explicit water-network edges and enable pairs mode |
| `--pair-id KIND` | `resid`, `atomid`, `oxygen_index` | Select the identifier type used in the pair file |
| `--workers N` | `auto` or positive integer | Override the mode-based automatic worker count |
| `--strict` | Flag | Stop on the first failed frame |
| `--output-layout LAYOUT` | `grouped`, `flat` | Select the per-frame structure-file layout |
| `--no-info` | Flag | Disable per-frame `*_info.md` files |
| `--no-gro` | Flag | Disable all structure GRO files |
| `--no-ring-gro` | Flag | Disable ring GRO files |
| `--no-half-cage-gro` | Flag | Disable half-cage GRO files |
| `--no-quasi-cage-gro` | Flag | Disable quasi-cage GRO files |
| `--no-cage-gro` | Flag | Disable cage GRO files |
| `--no-ice-gro` | Flag | Disable ice GRO files |
| `--no-xlsx` | Flag | Disable `summary.xlsx` |

### Bond Mode

Use `-b` / `--bond-mode` to override the graph setting supplied by the selected mode or `config.yaml`:

```bash
sqq analyze -i md.gro -b auto
sqq analyze -i md.gro --bond-mode hbond
sqq analyze -i md.gro -b oo
sqq analyze -i md.gro -b pairs --pairs pairs.txt
```

Available values are `auto`, `hbond`, `oo`, and `pairs`. `--pairs PAIRS.txt` used alone remains shorthand for pairs mode. Combining `--pairs` with `-b auto`, `-b hbond`, or `-b oo` is rejected. Pairs mode requires either `--pairs` or `graph.pair_file` in `config.yaml`.

## Output Structure

SQQ writes one folder per frame plus a global workbook:

```text
result_sqq/
  summary.xlsx
  run_config.yaml
  test1/
    test1_info.md
    ring/
      test1_ring_5.gro
      test1_ring_6.gro
    half_cage/
      hc_5r_5^5/
        test1_hc_5r_5^5.gro
    quasi_cage/
      qc_5r_5^3-6^2_55566/
        test1_qc_5r_5^3-6^2_55566.gro
    cage/
      5^12/
        test1_cage_5^12.gro
        test1_cage_5^12_empty.gro
        test1_cage_5^12_occupied.gro
    ice/
      test1_ice.gro
```

Each per-frame `*_info.md` report is arranged for inspection: it shows only reported ring sizes, reports final free-ring counts, includes the active network degree distribution, groups half-cage and quasi-cage isomers below composition totals, lists reported cage types vertically, and expands exact guest compositions across the cage-occupancy table. Internal `hc_` and `qc_` prefixes are omitted from report labels.

`summary.xlsx` remains plotting-oriented. Each analysis sheet keeps one input file or trajectory frame per row.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

Cage files include cage waters, CNT center atoms, and assigned guests. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
