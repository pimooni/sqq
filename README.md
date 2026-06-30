# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

SQQ builds a water network, reports coordination diagnostics, and finds rings, standard half-cages, quasi-cages, closed cages, reported-cage hydrate clusters, per-frame phase domains and boundaries, cage guest occupancy, F3/F4/Q_l order parameters, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

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

Modes do not change `quasi_cage.max_layers`; L1 remains the default in every mode. Use `--quasi-max-layer` explicitly for L2/L3. Automatic workers are capped by the number of independent GRO/XYZ files. A single coordinate file or XTC/TRR trajectory still runs with one worker. `--workers N` overrides the mode percentage. Parallel runs display live aggregate stages and up to six active files with per-stage and per-file timings.

## Common Commands

Write a default configuration file:

```bash
sqq init -o config.yaml
```

Search 4/5/6 ring faces and report Type H cages:

```bash
sqq analyze -i md.gro -s 4,5,6 --cage-size H -o ./result_sqq_456
```

Explicitly report every detected cage composition in the selected search scope:

```bash
sqq analyze -i md.gro -s 4,5,6 --cage-size all -o ./result_sqq_all_cages
```

Analyze connected reported-cage hydrate clusters:

```bash
sqq analyze -i md.gro -s 4,5,6 --hydrate-cluster on -o ./result_sqq_cluster
```

Enable outer quasi-cage layers:

```bash
sqq analyze -i md.gro --quasi-max-layer 3 -o ./result_sqq_l3
```

Use LAMMPS-style Q_l degree list and neighbors:

```bash
sqq analyze -i md.gro -q 4,6,8,10,12 --q-neighbor-mode lammps --q-cutoff 0.35 --q-n-neighbor 12
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
  report_types: auto
  max_faces: 20
  search_mode: grow
  seed_mode: ring
  occupancy_mode: polyhedron

hydrate_cluster:
  enabled: false
  min_cage: 2

order:
  f3f4_enabled: true
  q_enabled: true
  q_degree: [6, 12]
  q_neighbor_mode: graph
  q_cutoff_nm: 0.35
  q_n_neighbor: null

output:
  write_info: true
  write_gro: true
  write_order_tsv: false
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

`-s` / `--size` defines the ring-face sizes used during detection and, by default, reporting. With no dedicated report filter, SQQ reports all rings, quasi-cages, and cages found in that search scope. `--ring-size` and `--cage-size` can narrow the user-facing output afterward:

```bash
# Search 4/5/6, report only ring 5/6 and the Type H cage group
sqq analyze -i md.gro -s 4,5,6 --ring-size 5,6 --cage-size H
```

For example:

```bash
# Report every detected 4/5/6 ring, quasi-cage, and cage composition
sqq analyze -i md.gro -s 4,5,6

# Keep 4/5/6 rings and quasi-cages, but report only structure-I and structure-II cages
sqq analyze -i md.gro -s 4,5,6 --cage-size I,II
```

Cage report groups expand to scientific cage compositions:

```text
I     -> 5¹², 5¹²6²
II    -> 5¹², 5¹²6⁴
H     -> 5¹², 5¹²6⁸, 4³5⁶6³
HS-I  -> 5¹², 5¹²6², 5¹²6³
TS-I  -> 5¹², 5¹²6², 5¹²6³
I2II  -> 5¹²6³
```

Repeated cage types contributed by multiple groups are reported once. All detected cages still participate in half-cage, quasi-cage, and free-ring filtering. An explicit `--cage-size` changes user-facing counts and files, not topology ownership. Cage detection supports 4/5/6 faces; ring and quasi-cage detection also support size 7.

`--cage-size` accepts the comma-separated groups `I`, `II`, `H`, `HS-I`, `TS-I`, and `I2II`. The default `auto` scope follows `--size`; `all` explicitly requests the same all-detected behavior. Use `auto` or `all` alone rather than combining either keyword with a group.

## Hydrate Cluster

`--hydrate-cluster on` analyzes the final reported cage set after `--cage-size` filtering. Cages become graph nodes and are connected through complete shared ring faces. When several detected cages reference the same face, ring-plane geometry keeps at most one cage on each physical side.

`--cluster-min-cage N` sets the minimum connected-component size; the default is `2`. Smaller components are counted as isolated cages.

Within each reported cluster, SQQ builds labelled first-shell fingerprints from neighboring cage types and shared-face sizes. Strict local sI/sII/sH seeds initialize phase evidence. sI and sII then expand through compatible face-labelled edges when a candidate has at least two accepted phase contacts; sH remains conservative and expands only through overlapping strict sH seeds. Cages claimed exclusively by one phase form deterministic per-frame domains. Remaining cages are reported as single-phase boundaries, interphase boundaries, ambiguous, or unclassified.

The analysis is off by default and does not alter ring, patch, cage, occupancy, order-parameter, or ice results. Classification follows the final `--cage-size` scope and is per-frame; temporal grain tracking and crystallographic orientation matching are not implemented.

`--cluster-detail on` adds the optional one-row-per-cluster `hydrate_cluster_detail` workbook sheet. The default workbook output includes `hydrate_cluster` and `hydrate_domain`; public motif output is not generated. Hydrate clusters do not produce separate GRO files.

## Useful Options

| Option | Possible values | Meaning |
| --- | --- | --- |
| `-i, --input INPUT` | `.gro`, `.xyz`, `.xtc`, or `.trr` file; directory; or glob | Input coordinate file or trajectory source |
| `-c, --config FILE` | YAML or JSON file | User configuration file |
| `-o, --output DIR` | Directory path; default `result_sqq` | Output directory |
| `-m, --mode MODE` | `00`, `50`, `99`; default `50` | Select rigorous, standard, or performance preset |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --size SIZES` | Comma-separated subset of `4,5,6,7` | Set ring and quasi-cage search sizes; cage search uses the selected `4,5,6` sizes |
| `--ring-size SIZES` | `auto` or a comma-separated subset of `--size` | Report only these searched ring sizes |
| `--cage-size GROUPS` | `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, `I2II`; groups may be comma-separated | Restrict cage reporting; default `auto` follows `--size` |
| `--max-cage-face N` | Positive integer; default `20` | Limit generated cage search compositions |
| `--hydrate-cluster VALUE` | `on`, `off`; default `off` | Enable reported-cage hydrate_cluster analysis |
| `--cluster-min-cage N` | Positive integer; default `2` | Minimum connected cage count required for one hydrate_cluster |
| `--cluster-detail VALUE` | `on`, `off`; default `off` | Add the one-row-per-cluster `hydrate_cluster_detail` workbook sheet |
| `--pattern PATTERN` | Glob; default `*.gro` | Select files when `--input` is a directory |
| `--top, --topology FILE.gro` | GRO topology file | Supply topology/structure data for XTC/TRR input |
| `--recursive` | Flag; default off | Search input directories recursively |
| `--quasi-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base and side size lists together |
| `--quasi-base-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base-ring size list |
| `--quasi-side-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage side-ring size list |
| `--quasi-max-layer N` | Positive integer; default `1` | Report quasi-cage layers up to N |
| `--no-q` | Flag; default off | Disable Steinhardt Q_l order-parameter calculation |
| `-q, --q-degree L1,L2` | Comma-separated non-negative integers; default `6,12` | Select Q_l degree list to report, e.g. `4,6,8,10,12` |
| `--q-neighbor-mode MODE` | `graph`, `cutoff`, `nearest`, `lammps`; default `graph` | Select the neighbor source used by Q_l |
| `--q-cutoff NM` | Positive float in nm; default `0.35` | Q_l neighbor cutoff for cutoff/nearest/lammps modes |
| `--q-n-neighbor N` | Positive integer or `NULL`; default `NULL`, or `12` in lammps mode | Fixed Q_l neighbor count |
| `--pairs FILE` | Text pair-map file | Supply explicit water-network edges and enable pairs mode |
| `--pair-id KIND` | `resid`, `oxygen_index`, `atomid`; default `resid` | Select the identifier type used in the pair file |
| `--workers N` | `auto` or a positive integer | Override the mode-based automatic worker count |
| `--strict` | Flag; default off | Stop on the first failed frame |
| `--output-layout LAYOUT` | `grouped`, `flat`; default `grouped` | Select the per-frame structure-file layout |
| `--no-info` | Flag | Disable per-frame `*_info.md` files |
| `--no-gro` | Flag | Disable all structure GRO files |
| `--no-ring-gro` | Flag | Disable ring GRO files |
| `--no-half-cage-gro` | Flag | Disable half-cage GRO files |
| `--no-quasi-cage-gro` | Flag | Disable quasi-cage GRO files |
| `--no-cage-gro` | Flag | Disable cage GRO files |
| `--no-ice-gro` | Flag | Disable ice GRO files |
| `--no-xlsx` | Flag | Disable `summary.xlsx` |
| `--write-order-tsv` | Flag; default off | Write per-water `*_order_parameter.tsv` files |

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
    test1_order_parameter.tsv   # only with --write-order-tsv
```

Each per-frame `*_info.md` report is arranged for inspection: it shows only reported ring sizes, reports final free-ring counts, includes the active network degree distribution, groups half-cage and quasi-cage isomers below composition totals, and keeps cage composition totals plus cage isomers in one vertical `Cage` table. When enabled, the same report adds `Hydrate Cluster`, hierarchy, detail, domain, and boundary sections. Internal `hc_` and `qc_` prefixes are omitted from report labels.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary.xlsx` remains plotting-oriented. Each analysis sheet keeps one input file or trajectory frame per row. When hydrate_cluster is enabled, the workbook adds `hydrate_cluster` for per-frame totals and `hydrate_domain` for one row per domain; `--cluster-detail on` additionally writes `hydrate_cluster_detail` with one row per cluster. The `order_parameter` sheet contains `F3_mean`, `F3_count`, `F4_mean`, `F4_count`, and one mean/count pair for each requested Q_l degree. By default this gives `q6_mean`, `q6_count`, `q12_mean`, and `q12_count`.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

Cage files include cage waters, CNT center atoms, and assigned guests. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
