# SQQ

**SQQ (Shell Quant Qualifier)** is a Python program for identifying and quantifying water-shell topologies in molecular dynamics trajectories.

It analyzes water, ice, and hydrate-like structures from MD frames by building a water network, finding rings, standard half-cages, quasi-cages, closed cages, guest occupancy, F3/F4 order metrics, and ice-like waters.

Developer-level algorithm notes are kept in `docs/design.md`. Versioned update notes are kept in `docs/update.md`.

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

Direct glob of GRO frames:

```powershell
sqq analyze -i "./gro/*.gro" -c config.yaml -o ./result_sqq
```

For files named `1.gro` through `100.gro`, put them in one directory and use natural sorting:

```powershell
sqq analyze -i ./gro --pattern "*.gro" -c config.yaml -o ./result_sqq
```

The same multi-GRO selection can be controlled from `config.yaml`:

```yaml
input:
  pattern: "*.gro"
  recursive: false
```

Then run:

```powershell
sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

Parallel standalone GRO/XYZ frames:

```powershell
sqq analyze -i ./gro --n-jobs 4 -c config.yaml -o ./result_sqq
```

XTC/TRR trajectory with a topology/structure file:

```powershell
sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
```

Temporary 4/5/6 topology-size override without editing `config.yaml`:

```powershell
sqq analyze -i md.gro -c config.yaml --sizes 4,5,6 -o ./result_sqq_456
```

Use specific overrides when ring, quasi-cage, and cage search sizes should differ:

```powershell
sqq analyze -i md.gro -c config.yaml --ring-sizes 4,5,6 --quasi-sizes 4,5,6 --cage-sizes 5,6
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

quasi_cage:
  mode: general
  enabled: true
  base_sizes: auto
  side_sizes: auto
  max_combinations_per_base: 50000
  max_layers: 3
  max_rings_per_layer: 6
  max_layer_states_per_seed: 200
  max_candidates_per_edge: 4
  max_layer_candidates: 24

cage:
  enabled: true
  ring_sizes: [5, 6]
  target_types: [512, 51262, 51263, 51264]
  output_other: false
  other_max_faces: 20
  search_mode: grow
  seed_mode: ring
  max_states_per_seed: 20000
  max_total_states: 5000000
  max_boundary_candidates: 8
  occupancy_mode: polyhedron

order:
  focus_waters: []

parallel:
  n_jobs: auto

output:
  write_info: true
  write_gro: true
  write_ring_gro: true
  write_half_cage_gro: true
  write_quasi_cage_gro: true
  write_cage_gro: true
  write_ice_gro: true
  write_tsv: false
  write_vmd: false
  write_xlsx_summary: true
  structure_layout: grouped
```

Configuration priority:

```text
built-in defaults < config.yaml < command-line options
```

## Output

SQQ writes one folder per analyzed frame plus global summaries:

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
      hc_5r_5⁵/
        test1_hc_5r_5⁵.gro
    quasi_cage/
      qc_5r_5³6²₅₅₅₆₆/
        test1_qc_5r_5³6²₅₅₅₆₆.gro
    cage/
      5¹²/
        test1_cage_5¹².gro
        test1_cage_5¹²_empty.gro
        test1_cage_5¹²_occupied.gro
      5¹²6²/
        test1_cage_5¹²6²_ETH.gro
        test1_cage_5¹²6²_MET+ETH.gro
    ice/
      test1_ice.gro
```

Output ownership is `cage > quasi_cage > half_cage > ring`. A half-cage whose ring set is a true subset of a quasi-cage is consumed by that quasi-cage and is not counted or written again; nested half-cages are also reduced to the larger half-cage patch. Cage GRO files include cage waters, CNT center atoms, and assigned guest molecules. `*_occupied.gro` is the all-occupied cage set, while labels such as `*_MET.gro`, `*_METx2.gro`, or `*_MET+ETH.gro` are exact guest-composition subsets generated from the guest residue names present in the frame.
Optional `*_membership.tsv` and `*_f3f4.tsv` files can be enabled with `output.write_tsv: true`.

Per-frame folders only create `ring`, `half_cage`, `quasi_cage`, `cage`, or `ice` subfolders when that output type has at least one file to write. In the default grouped layout, `half_cage`, `quasi_cage`, and `cage` add one more type folder before the GRO files; `ring` stays flat because it only has a few size classes. Use `output.structure_layout: flat` or `--output-layout flat` to keep the same-folder file layout.

The workbook `summary.xlsx` starts with `summary` for run/output information and molecule counts, then `frame` for per-frame core counts, then a graph tab named by the active mode (`hbond`, `oo_connection`, or `pair_connection`), followed by `ring`, `half_cage`, `quasi_cage`, `cage`, `cage_occupancy`, `cage_isomer`, `f3f4`, `ice`, and `config`.

Output type control:

```yaml
output:
  write_info: true          # per-frame *_info.md
  write_gro: true           # master switch for all structure GRO files
  write_ring_gro: true
  write_half_cage_gro: true
  write_quasi_cage_gro: true
  write_cage_gro: true
  write_ice_gro: true
  write_xlsx_summary: true  # global summary.xlsx
  structure_layout: grouped # grouped or flat
```

Equivalent command-line shortcuts:

```powershell
sqq analyze -i ./gro -o ./result_sqq --no-info
sqq analyze -i ./gro -o ./result_sqq --no-gro
sqq analyze -i ./gro -o ./result_sqq --no-ring-gro --no-half-cage-gro --no-quasi-cage-gro --no-ice-gro
sqq analyze -i ./gro -o ./result_sqq --no-xlsx
sqq analyze -i ./gro -o ./result_sqq --output-layout flat
```

## Implemented Analysis

- GRO and XYZ input in the source tree; XTC/TRR input through MDAnalysis when installed.
- Orthorhombic minimum-image PBC.
- Water and guest selection by residue and atom names.
- Shared water graph for ring, open cage patches, closed cages, F3/F4, and ice metrics.
- `bond_mode=auto/hbond/oo/pairs`; `pairs` reads a user-provided water-neighbor file. Use `bond_mode: oo` when reproducing pure O-O topology counts.
- Non-recursive DFS ring search for primitive/chordless rings.
- Ring search supports 4/5/6/7-member rings, but defaults to `[5, 6]`.
- Layered open-patch search for `half_cage` and `quasi_cage`. L1 must be a closed side-ring circle around the base ring; L2 and L3 may be dangling rings or connected dangling ring chains. Standard half-cages are `hc_5r_5⁵`, `hc_6r_5⁶`, and `hc_6r_5⁶_6¹`; other non-cage layered patches are reported as `quasi_cage`.
- Quasi-cage base/side ring sizes follow `ring.sizes` by default, so `[4, 5, 6]` searches 4/5/6 open patches and `[4, 5, 6, 7]` searches 4/5/6/7 open patches.
- Layer assignment uses the lowest possible layer: if an outer ring touches both L1 and L2, it is treated as L2 rather than L3.
- Quasi-cage candidates are pruned locally before combination: rings must first share the relevant boundary edge, then the nearest ring-center candidates are kept with `quasi_cage.max_candidates_per_edge` and `quasi_cage.max_layer_candidates`.
- Default cage search by ring-face grow mode from ring seeds: grow connected face patches along open boundary edges, then validate a closed polyhedron by edge degree and Euler characteristic. Cage search supports 4/5/6-member faces through `cage.ring_sizes`, but defaults to `[5, 6]`; 7-member faces are intentionally not used for cage detection.
- Cage growth also works from topological boundary edges first, then orders and optionally truncates addable boundary rings by ring-center distance with `cage.max_boundary_candidates`.
- Optional `cage.output_other: true` adds Euler-compatible 4/5/6 unconventional cage targets up to `cage.other_max_faces`; the default only searches `512`, `51262`, `51263`, and `51264`.
- Optional `cage.search_mode=pair` for patch-pair comparison/debugging only.
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
| `quasi_cage.base_sizes`, `quasi_cage.side_sizes` | `auto` by default, meaning use `ring.sizes`; set explicit lists to restrict open-patch search |
| `quasi_cage.max_layers` | Maximum open-patch layer index after the base ring; default `3` reports L1/L2/L3 |
| `quasi_cage.max_rings_per_layer` | Maximum number of dangling rings in one L2/L3 connected layer subset |
| `quasi_cage.max_layer_states_per_seed` | Per-seed cap for connected L2/L3 layer subsets to keep exhaustive growth finite |
| `quasi_cage.max_candidates_per_edge` | Number of nearest side-ring candidates kept per base-ring edge before L1 combinations |
| `quasi_cage.max_layer_candidates` | Number of nearest L2/L3 frontier candidates kept before dangling-layer growth |
| `cage.ring_sizes` | Ring-face sizes allowed in cage search; default `[5, 6]`, optional `[4, 5, 6]` |
| `cage.output_other` | Enable unconventional 4/5/6 cage targets; default `false` |
| `cage.other_max_faces` | Maximum face count for generated unconventional cage targets; default `20` |
| `cage.search_mode` | `grow` by default; `pair` for comparison/debugging |
| `cage.seed_mode` | `ring` by default for speed; `patch` starts grow from half_cage/quasi_cage patches, `auto` uses patch seeds when patches exist and falls back to ring |
| `cage.max_states_per_seed`, `cage.max_total_states` | Cage-grow search limits; if reached, the frame info file writes a warning because cage counts may be incomplete |
| `cage.max_boundary_candidates` | Number of nearest addable boundary rings kept at each cage-growth step; set `0` to keep all topological boundary candidates |
| `cage.occupancy_mode` | `polyhedron` by default; `center` and `auto` are also available |
| `input.xtc_stride` | Read every Nth XTC/TRR frame |
| `order.focus_waters` | Residue ids whose mean F3/F4 should also be reported |
| `parallel.n_jobs` / `--n-jobs` | Parallel worker count for independent `.gro`/`.xyz` files; default `auto` runs serially |
| `output.write_info` | Write per-frame `_info.md` files |
| `output.write_gro` | Master switch for structure files used for visualization |
| `output.write_ring_gro`, `output.write_half_cage_gro`, `output.write_quasi_cage_gro`, `output.write_cage_gro`, `output.write_ice_gro` | Enable or disable individual structure-output classes |
| `output.write_xlsx_summary` | Write the global `summary.xlsx` |
| `output.structure_layout` / `--output-layout` | `grouped` puts GRO files under `ring/`, `half_cage/<type>/`, `quasi_cage/<type>/`, `cage/<type>/`, and `ice/`; `flat` keeps the same-folder layout |

Useful CLI shortcuts:

| Option | Effect |
| --- | --- |
| `--sizes 4,5,6` | Override `ring.sizes`, `quasi_cage.base_sizes`, `quasi_cage.side_sizes`, and `cage.ring_sizes` together; including 4 enables generated other cages unless `--no-other-cages` is used |
| `--ring-sizes 4,5,6` | Override only `ring.sizes` |
| `--quasi-sizes 4,5,6` | Override both `quasi_cage.base_sizes` and `quasi_cage.side_sizes` |
| `--quasi-base-sizes 4,5,6` | Override only `quasi_cage.base_sizes` |
| `--quasi-side-sizes 4,5,6` | Override only `quasi_cage.side_sizes` |
| `--cage-sizes 4,5,6` | Override only `cage.ring_sizes`; including 4 enables generated other cages unless `--no-other-cages` is used |
| `--other-cages` | Enable generated unconventional 4/5/6 cage targets |
| `--no-other-cages` | Disable generated unconventional cage targets |
| `--other-max-faces N` | Set the maximum face count for generated unconventional cages |

## Current Limits

- CHILL-style ice classification is implemented as a topology/coordination classifier; separate atomistic Ih/Ic stacking assignment can be refined later if needed.
- Only orthorhombic boxes are supported in the implemented PBC path.












