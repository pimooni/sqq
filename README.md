# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

Current release: **0.2.6**

SQQ builds a water network, reports coordination diagnostics, and finds rings, standard half-cages, quasi-cages, closed cages, reported-cage hydrate clusters, per-frame phase domains and boundaries, cage guest occupancy, F3/F4/Q_l order parameters, MCG/DHOP hydrate-nucleation order parameters, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

## Changed in 0.2.6

- Worker parsing is now form-based: `-w 1` means one worker, while `-w 1.0` or `-w 100%` means all detected physical cores before the reserve-one-core clamp.
- Runtime metadata is aligned with the `summary.xlsx` home sheet: terminal output now reports `SQQ version`, final effective `Graph mode`, worker policy, backend, and resolved workers using the same wording as the dashboard.
- Graph mode display now preserves both requested and effective modes: `auto -> hbond`, `auto -> oo`, or `auto -> mixed (hbond, oo)`; explicit modes display as `hbond`, `oo`, or `pairs`.
- `summary.xlsx` keeps the `quasi_cage` sheet compact by aggregating quasi-cage isomers into composition-level columns; exact quasi-cage isomers are written to `summary_detail/quasi_cage_isomer.csv`.
- `run_config.yaml` keeps the raw configuration and adds a `run` block with resolved runtime metadata. Scientific analysis algorithms, coordinates, molecule membership, and topology counts are unchanged from 0.2.5.

## Changed in 0.2.5

- Worker control: use `--worker` / `-w` with either a physical-core fraction (`50%`, `0.5`, `1.0`) or an explicit worker count (`1`, `4`). SQQ reserves one physical core for the system and clamps by task count and platform limits.
- Summary dashboard clarity: the first workbook sheet reports `SQQ version` in Configuration and uses `Analysis Results (min / mean / max)` for per-frame result metrics. `Frames total / ok / failed` remains a run-level frame count.
- Single-file progress visibility: interactive serial runs now highlight the active stage with bold bright-blue text while keeping the compact three-row stage layout.

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
| `00` | Rigorous | Hydrogen bond | 4, 5, 6 | 25% of physical cores |
| `50` | Standard (default) | Auto | 5, 6 | 50% of physical cores |
| `99` | Performance screening | O-O connectivity | 5, 6 | 90% of physical cores |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m 50 -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_performance
```

Modes do not change `quasi_cage.max_layers`; L1 remains the default in every mode. Use `--quasi-max-layer` explicitly for L2/L3. Automatic workers use the mode fraction of detected physical cores, reserve one physical core for the system, and are capped by the number of independent GRO/XYZ files or selected trajectory frames. Multiple standalone files use spawned processes by default; a single indexed XTC/TRR trajectory can also distribute frames across spawned workers. `--worker N` / `-w N` overrides the mode percentage; integer text such as `1` or `4` is a worker count, while decimal text such as `0.5` or `1.0` and percentages such as `50%` or `100%` are physical-core fractions. The old `--workers` spelling is retained as a compatibility alias.

Version 0.2.6 uses process-based parallelism for independent GRO/XYZ files and selected XTC/TRR frames, so CPU-bound ring, quasi-cage, and cage searches can run on multiple cores. The main process alone owns the terminal panel and final workbook; workers report stage events through a process queue, analyze one file or a small trajectory-frame batch, write frame directories, and return summary rows. At most `3 * workers` process tasks are kept in flight; this bounds Future and serialization overhead without reducing the worker count. `parallel.math_threads: 1` prevents nested BLAS/OpenMP oversubscription.

The default `chordless`/`bounded` path preserves the established scientific definitions while accelerating neighbor generation, incremental chord pruning, L1 forward checking, cached layer growth, integer-mask subset ownership, and cage target/edge state pruning. Cage DFS also applies exact remaining-edge incidence and parity conditions before expansion. MDAnalysis supplies orthorhombic cutoff candidates when available, but SQQ still rechecks every distance and hydrogen-bond angle with its established float64 logic. F3 and graph-mode Q_l share one graph-vector cache; all Q_l degrees share candidate lists and spherical-angle work. Optional `ring.definition: shortest_path` applies the Franzblau shortest-path criterion and reuses bounded-BFS distance maps. Optional `quasi_cage.search_policy: exact` preserves distinct frontiers and enumerates connected L2/L3 subsets; these opt-in modes can change or add results. Candidate and state truncation is reported through frame warnings.

Optional scientific cage validation adds PBC-aware face planarity and edge-variation limits, manifold vertex-link checks, positive-volume validation, and volume-centroid cage centers. It remains disabled by default. Version 0.2.5 continues to use the existing orthorhombic box representation; non-orthogonal/triclinic boxes are outside this update.

The current release uses the same compact three-row stage model for serial and parallel progress: file preparation (`reading`, `settings`, `selecting`), core topology search (`graph`, `ring`, `half/quasi`, `cage`, and optional `cluster`), and post-processing (`filtering`, `order`, `ice`, `output`). In interactive single-file runs, the active stage is highlighted with bold bright-blue ANSI text. The `cluster` stage appears only when hydrate-cluster analysis is enabled. Parallel runs also show aggregate stage counts and up to six active files with per-stage and per-file timings.

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

Enable bounded outer quasi-cage layers, or opt into exact connected-subset growth:

```bash
sqq analyze -i md.gro --quasi-max-layer 3 -o ./result_sqq_l3
sqq analyze -i md.gro --quasi-max-layer 3 --quasi-search-policy exact -o ./result_sqq_l3_exact
```

Opt into Franzblau shortest-path rings:

```bash
sqq analyze -i md.gro --ring-definition shortest_path -o ./result_sqq_sp_ring
```

Opt into stricter scientific cage validation:

```bash
sqq analyze -i md.gro --cage-scientific-validation on -o ./result_sqq_scientific
```

Use LAMMPS-style Q_l degree list and neighbors:

```bash
sqq analyze -i md.gro -q 4,6,8,10,12 --q-neighbor-mode lammps --q-cutoff 0.35 --q-n-neighbor 12
```
Enable the optional MCG-3 and DHOP30 variants in addition to default MCG-1 and DHOP35:

```bash
sqq analyze -i md.gro --mcg3 on --dhop30 on -o ./result_sqq_order
```

Disable per-frame info or structure GRO output:

```bash
sqq analyze -i md.gro --no-info -o ./result_sqq_no_info
sqq analyze -i md.gro --no-gro -o ./result_sqq_report
```

Parallelize independent GRO/XYZ files with spawned processes (the default backend):

```bash
sqq analyze -i ./gro --pattern "*.gro" --parallel-backend process -w 4 -o ./result_sqq
```

The same process backend parallelizes selected frames of one indexed trajectory:

```bash
sqq analyze -i traj.xtc --top topol.gro -w 4 -o ./result_sqq
```

Use `--parallel-backend serial` for an exact one-process comparison. `thread` is retained as a compatibility backend, but CPU-bound Python topology search should normally use `process`.

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
  definition: chordless

quasi_cage:
  enabled: true
  base_sizes: auto
  side_sizes: auto
  max_layers: 1
  search_policy: bounded

cage:
  enabled: true
  report_types: auto
  max_faces: 20
  search_mode: grow
  seed_mode: ring
  fast_closure: true
  fast_closure_max_states: 20000
  scientific_validation: false
  max_face_planarity_rms_nm: 0.06
  max_face_edge_cv: 0.35
  min_cage_volume_nm3: 1.0e-6
  occupancy_mode: polyhedron

hydrate_cluster:
  enabled: false
  min_cage: 2

hydrate_order:
  mcg1_enabled: true
  mcg3_enabled: false
  mcg_guest_resnames: [CH4, MET]
  mcg_guest_cutoff_nm: 0.90
  mcg_water_cutoff_nm: 0.60
  mcg_cone_half_angle_deg: 45.0
  mcg_min_waters: 5
  dhop35_enabled: true
  dhop30_enabled: false
  dhop_neighbor_cutoff_nm: 0.35
  dhop_planar_counts: [11, 12]
  dhop_min_qualified_neighbors: 3

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
  backend: process
  workers: auto
  math_threads: 1
```

Configuration priority:

```text
built-in defaults < mode preset < config.yaml < explicit command-line options
```

## Parallel Execution

`parallel.backend: process` is the default for two or more independent GRO/XYZ inputs. SQQ uses the `spawn` start method on every supported platform. Each worker receives run configuration once, reads and writes its own frame, and sends only small stage events plus one summary row to the main process. This avoids the Python GIL limitation of the compatibility thread backend.

Automatic workers use the mode fraction of detected physical cores, then reserve one physical core for the operating system and cap the result by the number of files or selected trajectory frames. Physical-core detection prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker` / `-w` accepts either a fraction (`50%`, `0.5`, or `1.0` for 100%) or an explicit positive integer worker count (`1` means one worker). Windows `ProcessPoolExecutor` runs are capped at 61 workers; Linux workstations can use larger explicit values such as `-w 100`, subject to the reserve-one-core rule, task count, memory, and storage throughput.

One XTC/TRR file with `--top` is frame-parallel when the process backend resolves to more than one worker. Every worker opens a private MDAnalysis Universe once and seeks small contiguous batches of selected raw frame indexes; batch size is automatically bounded from 1 to 8, and complete coordinate arrays are not serialized between processes. Parent and worker trajectory readers are explicitly closed. Multiple trajectory files and the compatibility thread backend use the serial trajectory reader.

Process submission uses a bounded rolling queue of at most `3 * workers` tasks. This is a queue-depth limit, not a CPU limit: with 100 effective workers SQQ may keep up to 300 tasks submitted while still running as many as 100 workers concurrently. Results are restored to original file/frame order before workbook writing.

The parent preserves input/frame order in `summary.xlsx`. Different standalone files must have unique case-insensitive stems because each stem is the output frame-directory name. Process runs set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, and `BLIS_NUM_THREADS` to `parallel.math_threads` while workers are spawned, then restore the parent environment.

The scheduling and search-cache refinements themselves do not change existing scientific definitions or values. Before the new hydrate descriptors were enabled, they reduced the local `1200ns.gro` serial run from about 26.6 s to 18.2 s. The current default run, including MCG-1 and DHOP35, completed in about 21.6 s on the same host; every overlapping pre-existing analysis column matched the earlier workbook, while `order_parameter` gained the new columns. Performance depends on data, configuration, CPU, memory, and storage.

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

## Cage Fast Closure and Scientific Validation

One frame-local ring topology index stores `ring_by_id`, ring centers, `edge_to_ring_ids`, ring adjacency, and the symmetric distance cache. Half/quasi and cage searches reuse this object instead of rebuilding the same incidence and geometry data.

`cage.fast_closure: true` is the default. Only when generic grow reaches a configured state limit, SQQ uses an indexed half-cage overlap graph to assemble connected combinations of two to four standard half-cage patches. Every candidate must still match one generated face composition and pass the ordinary closed-polyhedron test. Existing grow detections are retained first, so exhaustive grow output and object ids remain unchanged; fast closure only adds a cage when the bounded grow path missed it. `--cage-fast-closure off` disables this supplement for exact comparison.

`cage.scientific_validation: false` is the default. When enabled with `--cage-scientific-validation on`, every accepted cage must additionally satisfy the configured PBC-aware face-planarity RMS and edge-length coefficient-of-variation limits, an edge-connected face shell, a single cyclic face link around every vertex, and a positive minimum triangulated volume. Accepted cages then use the volume centroid instead of the mean cage-water position. Enabling it can therefore remove distorted cages and can change guest occupancy or geometry-resolved hydrate-cluster edges. Raw ring and half/quasi searches, order parameters, and ice classification are unchanged; ownership-filtered free-ring and free-patch outputs can increase when a rejected cage no longer consumes them.

## Hydrate Cluster

`--hydrate-cluster on` analyzes the final reported cage set after `--cage-size` filtering. Cages become graph nodes and are connected through complete shared ring faces. When several detected cages reference the same face, ring-plane geometry keeps at most one cage on each physical side.

The hierarchy follows the HTR+ idea of classifying hydrate type, domains, and boundaries on a cage-connection graph ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)). SQQ implements this independently with labelled shared-face fingerprints, strict local seeds, mutually compatible expansion, and exclusive per-frame domains.

`--cluster-min-cage N` sets the minimum connected-component size; the default is `2`. Smaller components are counted as isolated cages.

Within each reported cluster, SQQ builds labelled first-shell fingerprints from neighboring cage types and shared-face sizes. Strict local sI/sII/sH seeds initialize phase evidence. The sH templates cover `5^12`, `4^3 5^6 6^3`, and `5^12 6^8` cages; the earlier two-anchor sH composite is retained as supplemental high-confidence evidence. All three phases expand through mutually compatible face-labelled edges when a candidate has at least two accepted phase contacts. Cages claimed exclusively by one phase form deterministic per-frame domains. Remaining cages are reported as single-phase boundaries, interphase boundaries, ambiguous, or unclassified.

The analysis is off by default and does not alter ring, patch, cage, occupancy, order-parameter, or ice results. Classification follows the final `--cage-size` scope and is per-frame; temporal grain tracking and crystallographic orientation matching are not implemented.

`--cluster-detail on` adds optional one-row-per-cluster `summary_detail/hydrate_cluster_detail.csv`. The workbook keeps the per-frame `hydrate_cluster` sheet, while one-row-per-domain details move to `summary_detail/hydrate_domain.csv`; public motif output is not generated. Hydrate clusters do not produce separate GRO files.

## Hydrate Nucleation Order Parameters

Version 0.2.5 computes MCG-1 and DHOP35 by default. These descriptors are independent of the optional reported-cage `hydrate_cluster` classifier: MCG works on selected methane-like guest centers and surrounding waters, while DHOP works on a dedicated O-O neighbor graph. They do not change graph, ring, patch, cage, occupancy, F3/F4/Q_l, hydrate-cluster, or ice results.

MCG follows the mutually coordinated guest definition. Guest pairs within `0.90 nm` are connected when at least five waters lie within `0.60 nm` of both guests and inside both 45-degree opposing cones. The threshold is **at least five**, not exactly five. MCG-1 keeps guest nodes with at least one qualifying MCG edge; optional MCG-3 applies a one-pass degree-at-least-three filter to the same qualifying graph. Connected components are measured only through qualifying MCG edges. The default guest residue names are `CH4` and `MET`; change `hydrate_order.mcg_guest_resnames` for another methane naming convention. If no configured guest type is present, MCG is reported as `N/A`, not zero.

DHOP builds its own orthorhombic-PBC oxygen graph with `hydrate_order.dhop_neighbor_cutoff_nm: 0.35`. This 0.35 nm default follows the all-atom TIP4P/Ice implementation used by Li et al.; use `0.325` in YAML when reproducing the original mW-water definition. For each central O-O bond, SQQ counts neighboring plane-normal pairs within 35 degrees (or 30 degrees for DHOP30), selects waters with counts 11 or 12, requires at least three similarly qualified neighbors, includes their first oxygen shell, and reports the largest connected water cluster. `DHOP35` and `DHOP30` name the angular thresholds, not the O-O cutoff. No transition-state value such as DHOP35=57 is hard-coded; such values are system- and condition-dependent.

Use `--mcg3 on` and `--dhop30 on` to add the optional variants. The primary metrics can be disabled only through YAML (`mcg1_enabled` or `dhop35_enabled`). All cutoff searches use deterministic cell lists and exact float64 minimum-image rechecks; there are no fixed neighbor-array limits.

References: Barnes et al., MCG ([DOI 10.1063/1.4871898](https://doi.org/10.1063/1.4871898)); Knott et al., MCG nucleation coordinate ([DOI 10.1021/jp507959q](https://doi.org/10.1021/jp507959q)); DeFever and Sarupria, DHOP ([DOI 10.1063/1.4996132](https://doi.org/10.1063/1.4996132)); Li et al., all-atom DHOP nucleation pathway ([DOI 10.1073/pnas.2011755117](https://doi.org/10.1073/pnas.2011755117)).

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
| `--cage-fast-closure VALUE` | `on`, `off`; default `on` | Enable indexed two-to-four half-cage closure after generic grow |
| `--cage-scientific-validation VALUE` | `on`, `off`; default `off` | Enable strict face/manifold/volume validation and volume centroids |
| `--hydrate-cluster VALUE` | `on`, `off`; default `off` | Enable reported-cage hydrate_cluster analysis |
| `--cluster-min-cage N` | Positive integer; default `2` | Minimum connected cage count required for one hydrate_cluster |
| `--cluster-detail VALUE` | `on`, `off`; default `off` | Add one-row-per-cluster `summary_detail/hydrate_cluster_detail.csv` |
| `--pattern PATTERN` | Glob; default `*.gro` | Select files when `--input` is a directory |
| `--top, --topology FILE.gro` | GRO topology file | Supply topology/structure data for XTC/TRR input |
| `--recursive` | Flag; default off | Search input directories recursively |
| `--quasi-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base and side size lists together |
| `--quasi-base-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base-ring size list |
| `--quasi-side-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage side-ring size list |
| `--quasi-max-layer N` | Positive integer; default `1` | Report quasi-cage layers up to N |
| `--quasi-search-policy POLICY` | `bounded`, `exact`; default `bounded` | Preserve bounded growth or enumerate connected outer-layer subsets |
| `--ring-definition DEFINITION` | `chordless`, `shortest_path`; default `chordless` | Select the detected ring definition |
| `--no-q` | Flag; default off | Disable Steinhardt Q_l order-parameter calculation |
| `-q, --q-degree L1,L2` | Comma-separated non-negative integers; default `6,12` | Select Q_l degree list to report, e.g. `4,6,8,10,12` |
| `--q-neighbor-mode MODE` | `graph`, `cutoff`, `nearest`, `lammps`; default `graph` | Select the neighbor source used by Q_l |
| `--q-cutoff NM` | Positive float in nm; default `0.35` | Q_l neighbor cutoff for cutoff/nearest/lammps modes |
| `--q-n-neighbor N` | Positive integer or `NULL`; default `NULL`, or `12` in lammps mode | Fixed Q_l neighbor count |
| `--mcg3 VALUE` | `on`, `off`; default `off` | Add optional MCG-3 calculation; MCG-1 remains enabled by default |
| `--dhop30 VALUE` | `on`, `off`; default `off` | Add optional DHOP30 calculation; DHOP35 remains enabled by default |
| `--pairs FILE` | Text pair-map file | Supply explicit water-network edges and enable pairs mode |
| `--pair-id KIND` | `resid`, `oxygen_index`, `atomid`; default `resid` | Select the identifier type used in the pair file |
| `--parallel-backend BACKEND` | `process`, `thread`, `serial`; default `process` | Select independent-file/frame execution backend |
| `--worker, -w N` | `auto`, a fraction (`50%`, `0.5`, `1.0`), or a positive integer (`1`, `4`) | Override the mode-based worker count; one physical core is reserved. Integer `1` means one worker, while `1.0` / `100%` means all physical cores before clamping. `--workers` remains a hidden compatibility alias |
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
| `--no-summary-detail` | Flag | Disable `summary_detail/*.csv` detail tables |
| `--cage-isomer-rows MODE` | `nonzero`, `all`; default `nonzero` | Choose whether `summary_detail/cage_isomer.csv` keeps only observed isomer rows or the full zero-filled matrix |
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

SQQ writes one folder per frame, a global workbook, and CSV detail tables:

```text
result_sqq/
  summary.xlsx
  summary_detail/
    cage_occupancy.csv
    cage_isomer.csv
    quasi_cage_isomer.csv
    hydrate_domain.csv
    hydrate_cluster_detail.csv        # only with --cluster-detail on
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

GRO structure folders, filenames, and title lines use portable ASCII structure labels since version 0.2.4, for example `5^126^2` and `qc_5r_5^36^2_56566`. Markdown/Excel scientific labels retain their readable superscript notation. This avoids Windows GBK/legacy-reader failures caused by Unicode superscript or subscript characters in generated GRO paths and titles.

Each per-frame `*_info.md` report is arranged for inspection. `Frame Information` begins with the SQQ version, report-generation date and local timezone, absolute source path, frame name, and trajectory time. The report shows requested/effective `graph_mode`, the effective `bond_mode`, only reported ring sizes, final free-ring counts, the active network degree distribution, groups half-cage and quasi-cage isomers below composition totals, and keeps cage composition totals plus cage isomers in one vertical `Cage` table. When enabled, the same report adds `Hydrate Cluster`, hierarchy, detail, domain, and boundary sections. Internal `hc_` and `qc_` prefixes are omitted from report labels.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary.xlsx` remains plotting-oriented. Its first sheet is a dashboard: Configuration includes `SQQ version` and requested/effective `Graph mode` such as `auto -> hbond`, and `Analysis Results (min / mean / max)` reports per-frame min/mean/max values while `Frames total / ok / failed` stays a run-level count. The analysis sheets keep one input file or trajectory frame per row, including `frame`, connection diagnostics, `ring`, `half_cage`, compact composition-level `quasi_cage`, `cage`, optional per-frame `hydrate_cluster`, `order_parameter`, `ice`, `detail_index`, and `config`. Multi-row and isomer detail tables are written as UTF-8-SIG CSV files in `summary_detail/`: `cage_occupancy.csv`, `cage_isomer.csv`, `quasi_cage_isomer.csv`, `hydrate_domain.csv`, and, with `--cluster-detail on`, `hydrate_cluster_detail.csv`. The `quasi_cage` workbook sheet aggregates exact quasi-cage isomers into composition-level columns such as `5r_5²6³`, while `quasi_cage_isomer.csv` keeps nonzero exact isomer rows with `quasi_cage_type`, `isomer`, and `count`. `cage_isomer.csv` defaults to observed nonzero isomer rows plus per-frame totals; use `--cage-isomer-rows all` to restore the full zero-filled matrix. The `order_parameter` sheet contains `F3_mean`, `F3_count`, `F4_mean`, `F4_count`, one mean/count pair for each requested Q_l degree, and the default `MCG-1` and `DHOP35` largest-cluster columns. `MCG-3` and `DHOP30` columns appear only when their switches are enabled. By default the Q_l columns are `q6_mean`, `q6_count`, `q12_mean`, and `q12_count`. Each per-frame report places a separate `Hydrate Nucleation Order Parameters` table immediately after the F3/F4/Q_l table.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

Cage files include cage waters, CNT center atoms, and assigned guests. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
