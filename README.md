# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

SQQ builds a water network, finds rings, standard half-cages, quasi-cages, closed cages, cage guest occupancy, F3/F4 metrics, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

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

| Mode | Purpose | Water graph | Ring/cage sizes | Other cages | Automatic workers |
| --- | --- | --- | --- | --- | --- |
| `00` | Rigorous | Hydrogen bond | 4, 5, 6 | Enabled | 25% of logical CPUs |
| `50` | Standard (default) | Auto | 5, 6 | Disabled | 50% of logical CPUs |
| `99` | Performance screening | O-O connectivity | 5, 6 | Disabled | 90% of logical CPUs |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m 50 -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_performance
```

Modes do not change `quasi_cage.max_layers`; L1 remains the default in every mode. Use `--quasi-max-layers` explicitly for L2/L3. Automatic workers are capped by the number of independent GRO/XYZ files. A single coordinate file or XTC/TRR trajectory still runs with one worker. `--workers N` overrides the mode percentage.

## Common Commands

Write a default configuration file:

```bash
sqq init -o config.yaml
```

Analyze 4/5/6 rings, quasi-cages, and cages:

```bash
sqq analyze -i md.gro -s 4,5,6 -o ./result_sqq_456
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
  chordless: true

quasi_cage:
  enabled: true
  base_sizes: auto
  side_sizes: auto
  max_layers: 1

cage:
  enabled: true
  ring_sizes: [5, 6]
  target_types: [512, 51262, 51263, 51264]
  output_other: false
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

## Useful Options

| Option | Possible values | Meaning |
| --- | --- | --- |
| `-i, --input INPUT` | GRO/XYZ/XTC/TRR file, directory, or glob | Input coordinate file or trajectory source |
| `-c, --config FILE` | YAML or JSON file | User configuration file |
| `-o, --output DIR` | Directory path | Output directory |
| `-m, --mode MODE` | `00`, `50`, `99` | Select rigorous, standard, or performance preset |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --sizes SIZES` | Comma-separated `4,5,6` subset | Set ring, quasi-cage, and cage face sizes together |
| `--pattern PATTERN` | Glob such as `*.gro` | File pattern for directory input |
| `--top FILE.gro` | GRO topology file | Topology/structure file for XTC/TRR input |
| `--recursive` | Flag: present or omitted | Search input directories recursively |
| `--ring-sizes SIZES` | Comma-separated `4,5,6,7` subset | Override only ring search sizes |
| `--quasi-sizes SIZES` | Comma-separated `4,5,6,7` subset | Override quasi-cage base and side sizes together |
| `--quasi-base-sizes SIZES` | Comma-separated `4,5,6,7` subset | Override quasi-cage base-ring sizes |
| `--quasi-side-sizes SIZES` | Comma-separated `4,5,6,7` subset | Override quasi-cage side-ring sizes |
| `--quasi-max-layers N` | Positive integer; `1`-`3` documented | Report quasi-cage layers up to N; default is 1 |
| `--cage-sizes SIZES` | Comma-separated `4,5,6` subset | Override cage face sizes |
| `--other-cages` | Flag: present or omitted | Include generated unconventional cages |
| `--no-other-cages` | Flag: present or omitted | Disable generated unconventional cages |
| `--other-max-faces N` | Positive integer | Maximum face count for unconventional cages |
| `--pairs FILE` | Text pair-map file | Supply explicit water-network edges and enable pairs mode |
| `--pair-id KIND` | `resid`, `atomid`, `oxygen_index` | Select the identifier type used in the pair file |
| `--workers N` | `auto` or positive integer | Override the mode-based automatic worker count |
| `--strict` | Flag: present or omitted | Stop on the first failed frame |
| `--output-layout LAYOUT` | `grouped`, `flat` | Select the per-frame structure-file layout |
| `--no-info` | Flag: present or omitted | Disable per-frame `*_info.md` files |
| `--no-gro` | Flag: present or omitted | Disable all structure GRO files |
| `--no-ring-gro` | Flag: present or omitted | Disable ring GRO files |
| `--no-half-cage-gro` | Flag: present or omitted | Disable half-cage GRO files |
| `--no-quasi-cage-gro` | Flag: present or omitted | Disable quasi-cage GRO files |
| `--no-cage-gro` | Flag: present or omitted | Disable cage GRO files |
| `--no-ice-gro` | Flag: present or omitted | Disable ice GRO files |
| `--no-xlsx` | Flag: present or omitted | Disable `summary.xlsx` |

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

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

Cage files include cage waters, CNT center atoms, and assigned guests. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
