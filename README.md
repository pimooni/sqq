# SQQ

**SQQ (Shell Quant Qualifier)** analyzes water-shell topology in molecular dynamics frames.

It builds a water network, finds rings, standard half-cages, quasi-cages, closed cages, cage guest occupancy, F3/F4 metrics, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

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

## Common Commands

Write a default configuration file:

```bash
sqq init -o config.yaml
```

Analyze 4/5/6 rings, quasi-cages, and cages:

```bash
sqq analyze -i md.gro --sizes 4,5,6 -o ./result_sqq_456
```

Enable outer quasi-cage layers:

```bash
sqq analyze -i md.gro --quasi-max-layers 3 -o ./result_sqq_l3
```

Disable structure GRO output and keep only reports:

```bash
sqq analyze -i md.gro --no-gro -o ./result_sqq_report
```

Parallelize independent GRO/XYZ files:

```bash
sqq analyze -i ./gro --pattern "*.gro" --n-jobs 4 -o ./result_sqq
```

## Important Defaults

```yaml
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
```

Configuration priority:

```text
built-in defaults < config.yaml < command-line options
```

## Useful Options

| Option | Meaning |
| --- | --- |
| `-i, --input` | Input file, directory, or glob pattern |
| `--pattern "*.gro"` | File pattern for directory input |
| `--top topol.gro` | Topology/structure file for XTC/TRR input |
| `-c, --config config.yaml` | User configuration file |
| `-o, --output DIR` | Output directory |
| `--sizes 4,5,6` | Set ring, quasi-cage, and cage ring sizes together |
| `--ring-sizes 4,5,6` | Override only ring search sizes |
| `--quasi-sizes 4,5,6` | Override quasi-cage base and side sizes |
| `--quasi-max-layers N` | Report quasi-cage layers up to N; default is 1 |
| `--cage-sizes 4,5,6` | Override cage face sizes |
| `--other-cages` | Include generated unconventional 4/5/6 cage targets |
| `--no-gro` | Disable all structure GRO files |
| `--no-xlsx` | Disable `summary.xlsx` |
| `--output-layout flat` | Write all per-frame GRO files in one folder |

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
