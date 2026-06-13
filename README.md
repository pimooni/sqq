# SQQ

**SQQ (Shell Quant Qualifier)** is a Python program for identifying and quantifying water-shell topologies in molecular dynamics trajectories.

It analyzes water, ice, and hydrate-like structures from MD frames by building a water network, finding rings, cups (half-cages), closed cages, guest occupancy, F3/F4 order metrics, and ice-like waters.

## Install

Development install from this directory:

```powershell
pip install -e .
```

After installation, use the unified command:

```powershell
sqq init -o config.yaml
sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

During local development without installation:

```powershell
python -m sqq init -o config.yaml
python -m sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

## Quick Start

Single GRO file:

```powershell
sqq analyze -i test1.gro -c config.yaml -o ./result_sqq
```

Directory of GRO frames:

```powershell
sqq analyze -i ./gro --pattern "*.gro" -c config.yaml -o ./result_sqq
```

Parallel standalone GRO/XYZ frames:

```powershell
sqq analyze -i ./gro --n-jobs 4 -c config.yaml -o ./result_sqq
```

XTC/TRR trajectory with a topology/structure file:

```powershell
sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
```

## Important Defaults

```yaml
input:
  xtc_stride: 1

graph:
  bond_mode: auto
  oo_cutoff_nm: 0.35
  hbond_distance_nm: 0.35
  hbond_angle_deg: 30.0
  pair_file: null
  pair_id: resid

ring:
  sizes: [5, 6]
  primitive: true
  chordless: true

cup:
  mode: general
  enabled: true
  base_sizes: auto
  side_sizes: auto

cage:
  enabled: true
  ring_sizes: [5, 6]
  target_types: [512, 51262, 51263, 51264]
  output_other: false
  other_max_faces: 20
  search_mode: grow
  seed_mode: ring
  occupancy_mode: polyhedron

order:
  focus_waters: []

parallel:
  n_jobs: auto

output:
  write_gro: true
  write_tsv: false
  write_vmd: false
  write_xlsx_summary: true
```

Configuration priority:

```text
built-in defaults < config.yaml < command-line options
```

## Output

SQQ writes HA-style per-frame folders plus global summaries:

```text
result_sqq/
  summary.xlsx
  run_config.yaml
  test1/
    test1_info.md
    test1_ring5.gro
    test1_ring6.gro
    test1_cup5_55555.gro
    test1_512.gro
    test1_512_empty.gro
    test1_512_occupied.gro
    test1_ice.gro
```

Ring and cup outputs are free objects after cage ownership is removed. Cage GRO files include cage waters, CNT center atoms, and assigned guest molecules. Summary tables include per-cage-type empty/occupied/multi/guest columns and readable cage-isomer breakdowns.
Optional `*_membership.tsv` and `*_f3f4.tsv` files can be enabled with `output.write_tsv: true`.

## Implemented Analysis

- GRO and XYZ input in the source tree; XTC/TRR input through MDAnalysis when installed.
- Orthorhombic minimum-image PBC.
- Water and guest selection by residue and atom names.
- Shared water graph for ring, cup, cage, F3/F4, and ice metrics.
- `bond_mode=auto/hbond/oo/pairs`; `pairs` reads a user-provided water-neighbor file.
- Non-recursive DFS ring search for primitive/chordless rings.
- Ring search supports 4/5/6/7-member rings, but defaults to `[5, 6]`.
- General cup search from base rings and side-ring closure; by default cup base/side ring sizes follow `ring.sizes`, so `[4, 5, 6]` searches 4/5/6 cups and `[4, 5, 6, 7]` searches 4/5/6/7 cups.
- Default cage search by ring-face grow mode from ring seeds: grow connected face patches along open boundary edges, then validate a closed polyhedron by edge degree and Euler characteristic. Cage search supports 4/5/6-member faces through `cage.ring_sizes`, but defaults to `[5, 6]`; 7-member faces are intentionally not used for cage detection.
- Optional `cage.output_other: true` adds Euler-compatible 4/5/6 unconventional cage targets up to `cage.other_max_faces`; the default only searches `512`, `51262`, `51263`, and `51264`.
- Optional `cage.search_mode=pair` for HA/GRADE-style cup-pair comparison only.
- Default guest occupancy by `cage.occupancy_mode=polyhedron`, using an oriented solid-angle point-in-polyhedron test; `center` and `auto` modes are available for comparison/fallback.
- Cage isomer labels describe adjacent 6-ring face patterns, such as `6adj`, `6chain3`, `6star3`, or `6tri3+single`, instead of opaque `iso01` names.
- VMD helper colors: `512` blue, `51262` green, `51263` orange, `51264` red; ring centers use R4 gray, R5 purple, R6 tan, and R7 black.
- F3/F4 order metrics using the shared water graph and the reference-script formulas, including optional `order.focus_waters` averages.
- CHILL-style ice output: total ice-like waters, ice-I-like waters, and interfacial/intermediate ice waters.
- Explicit frame-level parallel execution for independent `.gro`/`.xyz` files via `--n-jobs N`.

## Main Parameters

| Parameter | Purpose |
| --- | --- |
| `water.resnames` | Water molecule residue names |
| `water.oxygen_names` | Oxygen atom names used as water graph nodes |
| `water.hydrogen_names` | Hydrogen atom names used for hydrogen-bond geometry |
| `guest.resnames` | Guest molecule residue names |
| `graph.bond_mode` | `auto`, `hbond`, `oo`, or `pairs` |
| `graph.pair_file` | Pair file used when `bond_mode=pairs`; each non-comment line has two ids |
| `graph.pair_id` | How pair ids are interpreted: `resid`, `oxygen_index`, or `atomid` |
| `ring.sizes` | Ring sizes to search; default `[5, 6]` |
| `cup.base_sizes`, `cup.side_sizes` | `auto` by default, meaning use `ring.sizes`; set explicit lists to restrict cup search |
| `cage.ring_sizes` | Ring-face sizes allowed in cage search; default `[5, 6]`, optional `[4, 5, 6]` |
| `cage.output_other` | Enable unconventional 4/5/6 cage targets; default `false` |
| `cage.other_max_faces` | Maximum face count for generated unconventional cage targets; default `20` |
| `cage.search_mode` | `grow` by default; `pair` for comparison/debugging |
| `cage.seed_mode` | `ring` by default for speed; `cup` starts grow from detected cups, `auto` uses cup seeds when cups exist and falls back to ring |
| `cage.occupancy_mode` | `polyhedron` by default; `center` and `auto` are also available |
| `input.xtc_stride` | Read every Nth XTC/TRR frame |
| `order.focus_waters` | Residue ids whose mean F3/F4 should also be reported |
| `parallel.n_jobs` / `--n-jobs` | Parallel worker count for independent `.gro`/`.xyz` files; default `auto` runs serially |
| `output.write_gro` | Write structure files for visualization |

## Current Limits

- CHILL-style ice classification is implemented as a topology/coordination classifier; separate atomistic Ih/Ic stacking assignment can be refined later if needed.
- Only orthorhombic boxes are supported in the implemented PBC path.












