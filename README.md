# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

Current release: **0.3.1**

SQQ builds a water network, reports coordination diagnostics, and finds rings, standard half-cages, quasi-cages, closed cages, topology-wide hydrate clusters, per-frame phase domains and boundaries, cage guest occupancy, F3/F4/Q_l order parameters, MCG/DHOP hydrate-nucleation order parameters, and ice-like waters. The default numeric modes use the complete Python backend. Version 0.3.1 also adds the focused C++17 backend selected by `-m cpp` for graph, internal ring, cage, occupancy, and F3/F4 analysis. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

## Changed in 0.3.1

- Added `sqq analyze -m cpp`, a C++17 native backend for water-graph construction, internal chordless 4/5/6-ring search, generic cage topology and cage isomers, automatic guest occupancy, and F3/F4.
- Python continues to own the CLI, input/topology readers, configuration, scheduling, structure writers, Markdown, summary CSV, and optional XLSX generation. The native extension releases the GIL while one frame is analyzed.
- Mode `cpp` defaults to `auto` graph selection, internal 4/5/6 rings, `f3,f4`, approximately 90% of physical cores with one reserved, and `info,cage-gro,summary-csv` output. It therefore avoids XLSX generation unless `summary-xlsx` is selected explicitly.
- C++ reports are intentionally compact: no public ring, half/quasi, cluster, ice, Q_l/MCG/DHOP, VMD, TSV, or detail-CSV output is generated. Occupancy is marked not evaluated when no selected guests exist.
- Output names are explicit in 0.3.1: `xlsx` is replaced by `summary-xlsx`, `summary-detail` by `summary-detail-csv`, and the new `summary-csv` writes one main-result CSV per applicable workbook sheet. The removed names have no compatibility aliases. SQQ-Py now defaults to `info,gro,summary-xlsx`, so detail CSV is opt-in.
- Unsupported explicit C++ options fail before analysis, and an unavailable or failed native extension never falls back silently to SQQ-Py.
- Release automation builds platform wheels for CPython 3.10-3.14 on Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, plus a source distribution. This describes the release workflow; it does not claim that 0.3.1 has already been published to PyPI.
- Numeric modes `00`, `09`, `50`, and `99` retain their complete SQQ-Py behavior and scientific output.
- Package and native-core versions are `0.3.1`, released Jul 19, 2026.

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

Building from source compiles the native extension and requires a C++17 compiler, CMake 3.20 or newer, Python development headers, and a platform build tool. Normal releases are intended to install a prebuilt wheel and do not compile C++ on the user's machine.

Then use:

```bash
sqq -h
sqq --version
sqq init -o config.yaml
sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```

Root help prints the SQQ version and release date immediately before the usage line. Use `sqq -v` or `sqq --version` for the version line alone.

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

### Input Units and Boxes

GRO and MDAnalysis trajectory coordinates are interpreted in nm. GRO accepts exactly one frame per file and rejects truncated atom blocks, missing or malformed box lines, extra non-empty records, and non-finite coordinates. Trajectory frames also require finite coordinates. XYZ coordinates are multiplied by `input.xyz_scale` / `--xyz-scale`; the default `0.1` assumes angstrom input, while `1.0` keeps nm values. SQQ accepts exactly one declared XYZ frame per file and rejects truncated, extra, malformed, or non-finite atom records. XYZ has no periodic box unless converted through another format.

GRO atom counts and the mandatory box line are validated. A three-value positive box is orthorhombic; an all-zero box is treated as non-periodic. Nine-value GRO boxes with nonzero tilt terms and trajectory frames with non-90-degree angles are rejected because triclinic minimum-image calculations are not implemented. Molecules are formed from contiguous residue blocks in source order, preventing wrapped or repeated residue IDs from merging distinct molecules.

## Analysis Modes

`-m` / `--mode` selects one of four complete SQQ-Py presets or the focused SQQ-CPP backend. The default remains numeric mode `50`:

| Mode | Purpose | Water graph | Search sizes | Automatic workers | Find cluster |
| --- | --- | --- | --- | --- | --- |
| `00` | Rigorous | Hydrogen bond | 4, 5, 6 | 25% of physical cores | on |
| `09` | Rigorous, high parallelism | Hydrogen bond | 4, 5, 6 | 90% of physical cores | on |
| `50` | Standard (default) | Auto | 5, 6 | 50% of physical cores | off |
| `99` | Performance screening | O-O connectivity | 5, 6 | 90% of physical cores | off |
| `cpp` | Focused native cage analysis | Auto | 4, 5, 6 (internal) | 90% of physical cores | unsupported |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m 09 -o ./result_rigorous_fast
sqq analyze -i ./gro -m 50 -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_performance
sqq analyze -i ./gro -m cpp -o ./result_cpp
```

Cluster-search precedence is independent of the other preset fields:

```text
--find-cluster > hydrate_cluster.enabled in config.yaml > mode preset
```

```bash
sqq analyze -i ./gro -m 09 --find-cluster off -o ./result_no_cluster
sqq analyze -i ./gro -m 50 --find-cluster on -o ./result_with_cluster
```

Numeric modes do not change `quasi_cage.max_layers`, `order.parameters`, or the initially selected output types; L1, `f3,f4`, and `info,gro,summary-xlsx` remain their defaults. Resolved cluster search always forces `cluster-gro`. It adds `summary-xlsx` only when neither `summary-csv` nor `summary-xlsx` is already selected, and never forces `info`. Thus a selected `summary-csv` remains CSV-only during cluster analysis. When `info` remains selected, every `*_info.md` report includes the compact `Hydrate Cluster` hierarchy; an output selection without `info` creates no Markdown report. Use `--quasi-max-layer` explicitly for L2/L3 and `--order-parameter` for another descriptor set. Automatic workers use the mode fraction of detected physical cores, reserve one physical core for the system, and are capped by the number of independent GRO/XYZ files or selected trajectory frames. Multiple standalone files use spawned processes by default; a single indexed XTC/TRR trajectory can also distribute frames across spawned workers. `--worker N` / `-w N` overrides the mode percentage; integer text such as `1` or `4` is a worker count, while decimal text such as `0.5` or `1.0` and percentages such as `50%` or `100%` are physical-core fractions. The old `--workers` spelling is retained as a compatibility alias.

SQQ uses process-based parallelism for independent GRO/XYZ files and selected XTC/TRR frames, so CPU-bound ring, quasi-cage, and cage searches can run on multiple cores. The main process alone owns the terminal panel and final summary outputs; process workers report stage events through a process queue, analyze one file or a small trajectory-frame batch, write frame directories, and return summary rows. At most `3 * workers` process or compatibility-thread tasks are kept in flight; this bounds Future and serialization overhead without reducing the worker count. `parallel.math_threads: 1` prevents nested BLAS/OpenMP oversubscription.

The default `chordless`/`bounded` path preserves the established scientific definitions while accelerating neighbor generation, incremental chord pruning, L1 forward checking, cached layer growth, integer-mask subset ownership, and cage target/edge state pruning. Cage DFS also applies exact remaining-edge incidence and parity conditions before expansion. MDAnalysis supplies orthorhombic cutoff candidates when available, but SQQ still rechecks every distance and hydrogen-bond angle with its established float64 logic. F3 and graph-mode Q_l share one graph-vector cache; all Q_l degrees share candidate lists and spherical-angle work. Optional `ring.definition: shortest_path` applies the Franzblau shortest-path criterion and reuses bounded-BFS distance maps. Optional `quasi_cage.search_policy: exact` preserves distinct frontiers and enumerates connected L2/L3 subsets; these opt-in modes can change or add results. Candidate and state truncation is reported through frame warnings.

Optional scientific cage validation adds PBC-aware face planarity and edge-variation limits, manifold vertex-link checks, positive-volume validation, and volume-centroid cage centers. It remains disabled by default. SQQ uses an orthorhombic box representation and now rejects non-orthogonal/triclinic input explicitly.

The current release uses the same compact three-row stage model for serial and parallel progress: file preparation (`reading`, `settings`, `selecting`), core topology search (`graph`, `ring`, `half/quasi`, `cage`, and optional `cluster`), and post-processing (`filtering`, `order`, `ice`, `output`). In interactive single-file runs, the active stage is highlighted with bold bright-blue ANSI text. The `cluster` stage appears only when hydrate-cluster analysis is enabled. Parallel runs also show aggregate stage counts and up to six active files with per-stage and per-file timings.

### Native SQQ-CPP Backend

`-m cpp` selects a separate native engine, not a fifth SQQ-Py scientific preset. Python still parses inputs and pair maps, selects molecules, normalizes configuration, schedules independent files/frames, and writes every report. The C++17 extension receives one normalized frame, performs the supported calculations while the GIL is released, and returns graph, internal ring, cage, occupancy, and F3/F4 records through the Python adapter.

| Setting | SQQ-CPP default |
| --- | --- |
| Graph | `auto` (`hbond` when usable water hydrogens exist, otherwise `oo`) |
| Rings | Internal chordless 4/5/6 search; no public ring output |
| Order parameters | `f3,f4` |
| Workers | Approximately 90% of detected physical cores, with one core reserved |
| Output | `info,cage-gro,summary-csv` |
| Occupancy | Automatic when selected guests exist; otherwise not evaluated |
| Periodicity | Orthorhombic or non-periodic only |

The native scope is deliberately small: generic Euler-compatible 4/5/6-face cage topology, the same cage labels/report groups, hexagonal-face cage isomers, polyhedron occupancy, and F3/F4. Ring detection is an internal cage prerequisite and is not exposed as a result.

Mode `cpp` accepts the compatible CLI subset:

- Input and run control: `-i`, `--top`, `-o`, `-c`, directory/glob controls, trajectory/input controls, `--strict`, and `--output-layout`.
- Graph: `-b auto|hbond|oo|pairs`, `--pairs`, and `--pair-id`.
- Topology: `-s` with a nonempty subset of `4,5,6`, `--cage-size`, `--max-cage-face`, and `--cage-scientific-validation`.
- Order: `--order-parameter f3`, `f4`, `f3,f4`, `all` (the supported pair), or `none`.
- Execution: `-w` / `--worker` and `--parallel-backend process|serial`.
- Output: `--output-type info|gro|cage-gro|summary-csv|summary-xlsx|all|none`, comma-separated where applicable; `--output-layout grouped|flat`; and `--cage-isomer-rows nonzero|all`. In this mode `gro` expands only to `cage-gro`, and `all` expands to `info,cage-gro,summary-csv,summary-xlsx`.

Unsupported in mode `cpp`:

- public ring tables/files, `--ring-size`, `ring-gro`, ring size 7, and `shortest_path` rings;
- half-cages, quasi-cages, and their settings or files;
- hydrate-cluster search/detail/GRO output;
- ice, Q_l, MCG, and DHOP analysis;
- VMD, membership/order TSV, `summary-detail-csv`, `cluster-detail`, and other detailed exports;
- the Python fast-closure option, thread backend, and triclinic boxes.

An explicit incompatible CLI or nondefault configuration request fails before analysis. A missing or failing native extension raises an error and never falls back to SQQ-Py.

The default compact layout is:

```text
result/
  summary_csv/
    summary.csv
    cage.csv
    cage_occupancy.csv        # only when selected guests exist
    cage_isomer.csv
    order_parameter.csv       # only when F3/F4 is selected
    config.csv
    failures.csv              # only when frames fail
  run_config.yaml
  frame_name/
    frame_name_info.md
    cage/<type>/frame_name_cage_<type>.gro
```

The default `summary_csv/` directory contains one UTF-8-SIG CSV for each applicable compact summary table. These files preserve the same columns and row order as the corresponding optional `summary.xlsx` sheets but are independent files rather than workbook tabs. Explicit `summary-xlsx` writes the dashboard, `cage`, `cage_isomer`, selected F3/F4 `order_parameter`, and `config` sheets; `failures` appears only when needed and `cage_occupancy` appears only when selected guests exist. The per-frame info report keeps frame/mode settings, molecules and connection diagnostics, cage topology/isomers, occupancy status, selected F3/F4, and warnings. It omits Ring, Half Cage, Quasi Cage, Hydrate Cluster, Hydrate Nucleation, and Ice sections. With no selected guests, occupancy is explicitly reported as not evaluated rather than empty. Native cage GRO files contain cage waters and assigned guests but do not add the synthetic `CNT` cage-center pseudoatom used by the full Python output.

Release CI is configured to build and test precompiled wheels for CPython 3.10-3.14 on Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, plus a source distribution. A wheel already contains the platform-native extension; end users installing such a wheel do not compile C++. A source install instead invokes the CMake/scikit-build-core build and therefore needs CMake 3.20 or newer and a local C++17 toolchain.

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

Analyze connected hydrate clusters from all detected cages:

```bash
sqq analyze -i md.gro -s 4,5,6 --find-cluster on -o ./result_sqq_cluster
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

Select F3/F4 plus a LAMMPS-style Q_l degree list and neighbors:

```bash
sqq analyze -i md.gro --order-parameter f3,f4,q4,q6,q8,q10,q12 --q-neighbor-mode lammps --q-cutoff 0.35 --q-n-neighbor 12
```
Select all F3/F4 and hydrate-nucleation descriptors:

```bash
sqq analyze -i md.gro --order-parameter f3,f4,mcg1,mcg3,dhop35,dhop30 -o ./result_sqq_order
```

`--order-parameter` replaces the complete selection rather than adding to the default. Use `all` for `f3,f4,q6,q12,mcg1,mcg3,dhop35,dhop30`, or `none` to skip all order-parameter calculations.

Select the complete output set with one positive option:

```bash
sqq analyze -i md.gro --output-type info,cage-gro,summary-xlsx -o ./result_sqq_report
sqq analyze -i md.gro --output-type none -o ./result_sqq_config_only
```

For SQQ-Py, `--output-type` defaults to `info,gro,summary-xlsx`. SQQ-CPP defaults to `info,cage-gro,summary-csv`. The mandatory `run_config.yaml` remains with every selection, including `none`.

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

input:
  pattern: "*.gro"
  xtc_stride: 1
  xyz_scale: 0.1

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
  mcg_guest_resnames: [CH4, MET]
  mcg_guest_cutoff_nm: 0.90
  mcg_water_cutoff_nm: 0.60
  mcg_cone_half_angle_deg: 45.0
  mcg_min_waters: 5
  dhop_neighbor_cutoff_nm: 0.35
  dhop_planar_counts: [11, 12]
  dhop_min_qualified_neighbors: 3

order:
  parameters: [f3, f4]
  q_neighbor_mode: graph
  q_cutoff_nm: 0.35
  q_n_neighbor: null

output:
  types: [info, gro, summary-xlsx]
  summary_csv_dir: summary_csv
  summary_detail_dir: summary_detail
  cage_isomer_rows: nonzero
  write_empty_files: false
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

Process submission uses a bounded rolling queue of at most `3 * workers` tasks. This is a queue-depth limit, not a CPU limit: with 100 effective workers SQQ may keep up to 300 tasks submitted while still running as many as 100 workers concurrently. Results are restored to original file/frame order before main-summary writing.

The parent preserves input/frame order in every selected main summary output. Different standalone files must have unique case-insensitive stems because each stem is the output frame-directory name. Process runs set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, and `BLIS_NUM_THREADS` to `parallel.math_threads` while workers are spawned, then restore the parent environment.

The scheduling and search-cache refinements themselves do not change existing scientific definitions or values. Before the new hydrate descriptors were enabled, they reduced the local `1200ns.gro` serial run from about 26.6 s to 18.2 s. A 0.2.3 benchmark that also selected MCG-1 and DHOP35 completed in about 21.6 s on the same host; every overlapping pre-existing analysis column matched the earlier workbook. Performance depends on data, configuration, CPU, memory, and storage.

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

Guest occupancy uses the configured center atom when available. Otherwise, guest atoms are PBC-unwrapped around one molecular anchor before calculating the centroid; the same helper is used by MCG. This correction can intentionally change occupancy counts relative to 0.2.6 or early 0.2.7 results for multi-atom guests crossing a periodic boundary.

## Hydrate Cluster

`--find-cluster on` analyzes every detected cage in the selected search scope. Cages become graph nodes and are connected through complete shared ring faces. When several detected cages reference the same face, ring-plane geometry keeps at most one cage on each physical side. `--cage-size` filters user-facing cage tables and files only; it does not remove cages from cluster connectivity or phase evidence.

The hierarchy follows the HTR+ idea of classifying hydrate type, domains, and boundaries on a cage-connection graph ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)). SQQ implements this independently with labelled shared-face fingerprints, strict local seeds, mutually compatible expansion, and exclusive per-frame domains.

`--cluster-min-cage N` sets the minimum connected-component size; the default is `2`. Smaller components are counted as isolated cages.

Within each cluster, SQQ builds labelled first-shell fingerprints from neighboring cage types and shared-face sizes. Strict local sI/sII/sH seeds initialize phase evidence. The sH templates cover `5^12`, `4^3 5^6 6^3`, and `5^12 6^8` cages; the earlier two-anchor sH composite is retained as supplemental high-confidence evidence. All three phases expand through mutually compatible face-labelled edges when a candidate has at least two accepted phase contacts. Cages claimed exclusively by one phase form deterministic per-frame domains.

After the exclusive sI/sII/sH domains are finalized, SQQ partitions the remaining cluster cages. A cage enters the generic boundary only when it is outside every phase domain and directly shares a complete cage face with at least one domain cage. Boundary search stops at this first external non-phase layer. Domain cages are never relabelled as boundary, and a direct shared-face contact between different phase domains leaves both endpoint cages in their original phases.

The resulting `classified_cage_ids`, `boundary_cage_ids`, `ambiguous_cage_ids`, and `unclassified_cage_ids` are mutually exclusive and together cover every cage in a reported cluster. Competing phase claims without boundary membership remain ambiguous; all other residual cages are unclassified. There are no `sI-boundary`, `sII-boundary`, `sH-boundary`, transition, or boundary-context categories. Neighboring cages can still share face-water coordinates in structure views, so cage ownership should be verified from cage IDs or detected cage/ring edges rather than coordinate-set overlap.

The default command uses mode `50`, so cluster search is off unless it is enabled by mode `00`/`09`, `hydrate_cluster.enabled` in the configuration, or an explicit `--find-cluster on`. Explicit `--find-cluster on|off` has highest priority. Cluster search does not alter ring, patch, cage, occupancy, order-parameter, or ice results. Classification is per-frame and independent of the cage reporting filter; temporal grain tracking and crystallographic orientation matching are not implemented.

Enabling cluster search always forces `cluster-gro`. If neither `summary-csv` nor `summary-xlsx` is selected, it also adds `summary-xlsx`; an existing `summary-csv` selection stays CSV-only. The selected main summary output gains its per-frame `hydrate_cluster` table, while native category structures are written under grouped layout as `<frame>/hydrate_cluster/<frame>_cluster_sI.gro`, `<frame>_cluster_sII.gro`, `<frame>_cluster_sH.gro`, and `<frame>_cluster_boundary.gro`. Flat layout places the same filenames directly in the frame directory. All same-category domains and clusters are aggregated into one file per frame. An absent category is omitted unless `output.write_empty_files: true`.

Cluster GRO files contain only complete water molecules belonging to the selected cage IDs; guests and CNT atoms are excluded. Ambiguous, unclassified, and isolated cages are not exported. Every atom keeps the exact wrapped coordinate from the analyzed frame, and every file keeps the original box; categories are never moved or unwrapped independently. Periodic or percolating networks may therefore still show bonds crossing a box face because no single-copy GRO representation can remove every periodic seam.

Cage IDs are mutually exclusive across sI, sII, sH, and boundary, but adjacent category files can contain the same face-water molecules because neighboring cages physically share them. When resolved cluster search is on and `info` is selected, `Frame Information` records `find_cluster` as `on` and the report adds one compact `Hydrate Cluster` hierarchy. Domain rows may be sI, sII, or sH; boundary and compact unclassified rows are subdivided by cage type. The compact unclassified count is the deduplicated unresolved set: stored ambiguous and unclassified IDs plus any uncategorized residual cluster cages. Main summary and cluster-detail output preserve the distinct scientific fields. Counts use unique cage IDs, zero-count rows are omitted, multiple clusters appear sequentially, and `isolated` appears once as the final top-level row without subtype children.

```text
## Hydrate Cluster

| item               | type         | cage_qty |
| ------------------ | ------------ | -------- |
| cluster_00001      | mixed        | 334      |
| ├ domain_00001     | sI           | ├ 66     |
|   ├ 5¹²            |              |   ├ 13   |
|   └ 5¹²6²          |              |   └ 53   |
| ├ domain_00002     | sII          | ├ 194    |
|   ├ 5¹²            |              |   ├ 131  |
|   └ 5¹²6⁴          |              |   └ 63   |
| ├ boundary         | boundary     | ├ 69     |
|   ├ 5¹²            |              |   ├ 24   |
|   └ 5¹²6³          |              |   └ 45   |
| └ unclassified     | unclassified | └ 5      |
|   ├ 5¹²6³          |              |   ├ 2    |
|   └ 4¹5¹⁰6²        |              |   └ 3    |
| isolated           | isolated     | 5        |
```

The compact table does not include exact IDs, seeds, confidence values, water/guest membership, or domain adjacency. Add `cluster-detail` to `--output-type` for `summary_detail/hydrate_domain.csv` and one-row-per-cluster `summary_detail/hydrate_cluster_detail.csv`. Explicit `cluster-detail` or `cluster-gro` selection requires cluster search. Turning search off writes neither `cluster-detail` nor `cluster-gro` and removes stale generated cluster GRO files. Public motif output is not generated.

## Hydrate Nucleation Order Parameters

MCG-1 and DHOP35 were introduced as defaults in 0.2.5. Since 0.2.7, every MCG/DHOP variant is selected explicitly through `--order-parameter`; the package default is only `f3,f4`. These descriptors are independent of the optional cage-topology `hydrate_cluster` classifier: MCG works on selected methane-like guest centers and surrounding waters, while DHOP works on a dedicated O-O neighbor graph. They do not change graph, ring, patch, cage, occupancy, F3/F4/Q_l, hydrate-cluster, or ice results.

MCG follows the mutually coordinated guest definition. Guest pairs within `0.90 nm` are connected when at least five waters lie within `0.60 nm` of both guests and inside both 45-degree opposing cones. The threshold is **at least five**, not exactly five. MCG-1 keeps guest nodes with at least one qualifying MCG edge; optional MCG-3 applies a one-pass degree-at-least-three filter to the same qualifying graph. Connected components are measured only through qualifying MCG edges. The default guest residue names are `CH4` and `MET`; change `hydrate_order.mcg_guest_resnames` for another methane naming convention. If no configured guest type is present, MCG is reported as `N/A`, not zero.

DHOP builds its own orthorhombic-PBC oxygen graph with `hydrate_order.dhop_neighbor_cutoff_nm: 0.35`. This 0.35 nm default follows the all-atom TIP4P/Ice implementation used by Li et al.; use `0.325` in YAML when reproducing the original mW-water definition. For each central O-O bond, SQQ counts neighboring plane-normal pairs within 35 degrees (or 30 degrees for DHOP30), selects waters with counts 11 or 12, requires at least three similarly qualified neighbors, includes their first oxygen shell, and reports the largest connected water cluster. `DHOP35` and `DHOP30` name the angular thresholds, not the O-O cutoff. No transition-state value such as DHOP35=57 is hard-coded; such values are system- and condition-dependent.

Select any combination with names such as `--order-parameter mcg1,mcg3,dhop35,dhop30`. Selection is separate from the numerical `hydrate_order` cutoff settings. All cutoff searches use deterministic cell lists and exact float64 minimum-image rechecks; there are no fixed neighbor-array limits.

References: Barnes et al., MCG ([DOI 10.1063/1.4871898](https://doi.org/10.1063/1.4871898)); Knott et al., MCG nucleation coordinate ([DOI 10.1021/jp507959q](https://doi.org/10.1021/jp507959q)); DeFever and Sarupria, DHOP ([DOI 10.1063/1.4996132](https://doi.org/10.1063/1.4996132)); Li et al., all-atom DHOP nucleation pathway ([DOI 10.1073/pnas.2011755117](https://doi.org/10.1073/pnas.2011755117)).

## Useful Options

| Option | Possible values | Meaning |
| --- | --- | --- |
| `-i, --input INPUT` | `.gro`, `.xyz`, `.xtc`, or `.trr` file; directory; or glob | Input coordinate file or trajectory source |
| `-c, --config FILE` | YAML or JSON file | User configuration file |
| `-o, --output DIR` | Directory path; default `result_sqq` | Output directory |
| `-m, --mode MODE` | `00`, `09`, `50`, `99`, `cpp`; default `50` | Select a complete Python preset or the focused C++17 backend |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --size SIZES` | Comma-separated subset of `4,5,6,7` | Set ring and quasi-cage search sizes; cage search uses the selected `4,5,6` sizes |
| `--ring-size SIZES` | `auto` or a comma-separated subset of `--size` | Report only these searched ring sizes |
| `--cage-size GROUPS` | `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, `I2II`; groups may be comma-separated | Restrict cage reporting; default `auto` follows `--size` |
| `--max-cage-face N` | Positive integer; default `20` | Limit generated cage search compositions |
| `--cage-fast-closure VALUE` | `on`, `off`; default `on` | Enable indexed two-to-four half-cage closure after generic grow |
| `--cage-scientific-validation VALUE` | `on`, `off`; default `off` | Enable strict face/manifold/volume validation and volume centroids |
| `--find-cluster VALUE` | `on`, `off`; default follows mode/config | Override all-detected-cage cluster search; `on` forces native cluster GRO and ensures at least one main summary format |
| `--cluster-min-cage N` | Positive integer; default `2` | Minimum connected cage count required for one hydrate_cluster |
| `--pattern PATTERN` | Glob; default `*.gro` | Select files when `--input` is a directory |
| `--top, --topology FILE.gro` | GRO topology file | Supply topology/structure data for XTC/TRR input |
| `--xyz-scale SCALE` | Positive float; default `0.1` | Multiply XYZ coordinates by this value to obtain nm; use `1.0` for XYZ already in nm |
| `--recursive` | Flag; default off | Search input directories recursively |
| `--quasi-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base and side size lists together |
| `--quasi-base-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage base-ring size list |
| `--quasi-side-size SIZES` | `auto` or a comma-separated subset of searched `4,5,6,7` | Override quasi-cage side-ring size list |
| `--quasi-max-layer N` | Positive integer; default `1` | Report quasi-cage layers up to N |
| `--quasi-search-policy POLICY` | `bounded`, `exact`; default `bounded` | Preserve bounded growth or enumerate connected outer-layer subsets |
| `--ring-definition DEFINITION` | `chordless`, `shortest_path`; default `chordless` | Select the detected ring definition |
| `--order-parameter NAMES` | `f3`, `f4`, `qN`, `mcg1`, `mcg3`, `dhop35`, `dhop30`, `all`, or `none`; comma-separated | Select the complete descriptor set; default `f3,f4`. `all` expands to `f3,f4,q6,q12,mcg1,mcg3,dhop35,dhop30` |
| `--q-neighbor-mode MODE` | `graph`, `cutoff`, `nearest`, `lammps`; default `graph` | Select the neighbor source used by Q_l |
| `--q-cutoff NM` | Positive float in nm; default `0.35` | Q_l neighbor cutoff for cutoff/nearest/lammps modes |
| `--q-n-neighbor N` | Positive integer or `NULL`; default `NULL`, or `12` in lammps mode | Fixed Q_l neighbor count |
| `--pairs FILE` | Text pair-map file | Supply explicit water-network edges and enable pairs mode |
| `--pair-id KIND` | `resid`, `oxygen_index`, `atomid`; default `resid` | Select the identifier type used in the pair file |
| `--parallel-backend BACKEND` | `process`, `thread`, `serial`; default `process` | Select independent-file/frame execution backend |
| `--worker, -w N` | `auto`, a fraction (`50%`, `0.5`, `1.0`), or a positive integer (`1`, `4`) | Override the mode-based worker count; one physical core is reserved. Integer `1` means one worker, while `1.0` / `100%` means all physical cores before clamping. `--workers` remains a hidden compatibility alias |
| `--strict` | Flag; default off | Stop on the first failed frame |
| `--output-layout LAYOUT` | `grouped`, `flat`; default `grouped` | Select the per-frame structure-file layout |
| `--output-type TYPES` | Comma-separated `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `summary-xlsx`, `summary-csv`, `summary-detail-csv`, `cluster-detail`, `all`, or `none` | Select the complete output set; SQQ-Py defaults to `info,gro,summary-xlsx`, while SQQ-CPP defaults to `info,cage-gro,summary-csv` |
| `--cage-isomer-rows MODE` | `nonzero`, `all`; default `nonzero` | SQQ-Py: control `summary_detail/cage_isomer.csv` when `summary-detail-csv` is selected; SQQ-CPP: control `summary_csv/cage_isomer.csv` and the optional `summary-xlsx` sheet |

### Bond Mode

Use `-b` / `--bond-mode` to override the graph setting supplied by the selected mode or `config.yaml`:

```bash
sqq analyze -i md.gro -b auto
sqq analyze -i md.gro --bond-mode hbond
sqq analyze -i md.gro -b oo
sqq analyze -i md.gro -b pairs --pairs pairs.txt
```

Available values are `auto`, `hbond`, `oo`, and `pairs`. `--pairs PAIRS.txt` used alone remains shorthand for pairs mode. Combining `--pairs` with `-b auto`, `-b hbond`, or `-b oo` is rejected. Pairs mode requires either `--pairs` or `graph.pair_file` in `config.yaml`.

### Output Selection

`--output-type TYPE[,TYPE...]` replaces the complete configured output list. Its SQQ-Py default is `info,gro,summary-xlsx`, equivalent to:

```yaml
output:
  types: [info, gro, summary-xlsx]
  summary_csv_dir: summary_csv
  summary_detail_dir: summary_detail
```

The SQQ-Py canonical names are `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `summary-xlsx`, `summary-csv`, `summary-detail-csv`, and `cluster-detail`. `gro` expands to the five ordinary ring/half/quasi/cage/ice GRO categories; `cluster-gro` is a separate search-dependent type. SQQ-CPP accepts only `info`, `gro`, `cage-gro`, `summary-csv`, `summary-xlsx`, `all`, and `none`; its `gro` means `cage-gro`, and `all` expands to `info,cage-gro,summary-csv,summary-xlsx`. SQQ-CPP does not support `summary-detail-csv`. `all` and `none` must appear alone. `run_config.yaml` is mandatory and is always rewritten, including with `--output-type none`.

`summary-xlsx` writes `summary.xlsx`. `summary-csv` writes each applicable main-summary table as a separate UTF-8-SIG file under `output.summary_csv_dir` (default `summary_csv/`); these are the CSV equivalents of workbook sheets, not detail records. `summary-detail-csv` selects the ordinary multi-row CSV tables under `output.summary_detail_dir` (default `summary_detail/`). Both directory settings must be different relative paths that resolve inside the selected output directory. `cluster-detail` separately selects `hydrate_domain.csv` and `hydrate_cluster_detail.csv`; explicit `cluster-detail` and `cluster-gro` selections require cluster search. Enabling cluster search always forces `cluster-gro` and ensures a main summary: if neither `summary-csv` nor `summary-xlsx` is selected, it adds `summary-xlsx`; if `summary-csv` is already selected, no workbook is added. Disabling search writes no cluster GRO files and cleans stale SQQ-generated `hydrate_cluster` directories or flat cluster filenames.

Explicit CLI selection has precedence:

```text
--output-type > output.types > engine default
```

`--hydrate-cluster`, `--cluster-detail`, `--no-output`, `--write-order-tsv`, and the individual `--no-*` output switches were removed in 0.2.8. They have no compatibility aliases. Likewise, `output.disabled_outputs` is rejected rather than migrated; configurations must use `output.types`.

When an existing output directory is reused, SQQ removes known stale files for output types outside the effective selection while preserving unrelated user files. Main-summary and detail-CSV cleanup is restricted to known SQQ-generated filenames inside the currently configured `summary_csv_dir` and `summary_detail_dir`; unknown files are preserved, and changing either setting does not make SQQ scan or clean a formerly configured directory. If no per-frame output type is selected and no unrelated file remains, the empty frame directory is removed.

## Output Structure

With the SQQ-Py default `--output-type info,gro,summary-xlsx`, SQQ writes one folder per frame and a global workbook. Detail CSV files are no longer part of the default:

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
    hydrate_cluster/                     # forced when cluster search is on
      test1_cluster_sI.gro               # omitted when the category is absent
      test1_cluster_sII.gro
      test1_cluster_sH.gro
      test1_cluster_boundary.gro
    test1_order_parameter.tsv   # only with output type order-tsv
```

Selecting `summary-csv` adds `<summary_csv_dir>/<sheet>.csv` for every applicable main-result table. Selecting `summary-detail-csv` separately adds ordinary multi-row files under `<summary_detail_dir>/`, including `failures.csv`, `cage_occupancy.csv`, `cage_isomer.csv`, and `quasi_cage_isomer.csv`; `cluster-detail` adds `hydrate_domain.csv` and `hydrate_cluster_detail.csv` there when cluster search is on.

Only selected or required outputs are written. With cluster search off, `--output-type none` leaves only `run_config.yaml` and does not create empty per-frame directories. With cluster search on, the same selector forces `summary-xlsx` and grouped or flat cluster GRO files because no main summary was selected; `--output-type summary-csv` instead remains CSV-only and adds cluster GRO. Cluster search does not force `info`, so no `*_info.md` report is created unless `info` is selected. Grouped layout uses `<frame>/hydrate_cluster/`; flat layout places the four canonical filenames at the frame root. Missing categories are omitted unless `output.write_empty_files` is true.

Without `--strict`, standalone serial/process/thread read failures become failed summary rows and analysis continues where the reader remains usable. Failed inputs appear in `summary.xlsx/failures`, `<summary_csv_dir>/failures.csv`, and `<summary_detail_dir>/failures.csv` when their respective output types are enabled, and always in the mandatory `run_config.yaml` `run.failures` list. With `--strict`, SQQ re-raises the error after updating `run_config.yaml` to `status: failed`.

GRO structure folders, filenames, and title lines use portable ASCII structure labels since version 0.2.4, for example `5^126^2` and `qc_5r_5^36^2_56566`. Markdown and main-summary scientific labels retain their readable superscript notation. This avoids Windows GBK/legacy-reader failures caused by Unicode superscript or subscript characters in generated GRO paths and titles.

Each per-frame `*_info.md` report is arranged for inspection. `Frame Information` begins with `sqq version`, `mode`, `date & time`, absolute `source`, `frame`, and `time_ps`. Numeric modes display as, for example, `09 (sqq-py)`; native mode displays as `sqq-cpp`. Requested/effective `graph_mode`, effective `bond_mode`, ring sizes, status, and molecule counts follow; `find_cluster` is present only for SQQ-Py. Python reports retain the total/free Ring table, grouped half/quasi and cage-isomer sections, optional compact Hydrate Cluster hierarchy, and the other enabled analyses. SQQ-CPP omits inapplicable sections and keeps connection diagnostics, cage topology/isomers, occupancy status, selected F3/F4, and warnings.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary-xlsx` writes the plotting-oriented `summary.xlsx` workbook. `summary-csv` uses the same applicable main-table mapping and writes one UTF-8-SIG file per sheet under `summary_csv_dir`, preserving table names, columns, row order, and values without Excel formatting or tabs. The first `summary` table is a dashboard: Configuration includes `SQQ version`, requested/effective `Graph mode` such as `auto -> hbond`, normalized `Order parameters`, `Find cluster`, and normalized `Output types`; `Analysis Results (min / mean / max)` reports per-frame min/mean/max values while `Frames total / ok / failed` stays a run-level count. Analysis tables such as connection diagnostics, `ring`, `half_cage`, compact composition-level `quasi_cage`, `cage`, optional `hydrate_cluster`, `order_parameter`, and `ice` keep one input file or trajectory frame per row. The other tables have metadata-specific row units: `summary` is a dashboard, `failures` has one failed input/frame per row, `detail_index` has one generated detail file per row, and `config` contains configuration metadata rather than frame rows. Ordinary multi-row and isomer tables are written separately under `summary_detail_dir` only when `summary-detail-csv` is selected: optional `failures.csv`, `cage_occupancy.csv`, `cage_isomer.csv`, and `quasi_cage_isomer.csv`. The separate `cluster-detail` type writes `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. The compact `quasi_cage` table aggregates exact quasi-cage isomers into composition-level columns such as `5r_5²6³`, while the detail `quasi_cage_isomer.csv` keeps nonzero exact isomer rows with `quasi_cage_type`, `isomer`, and `count`. `cage_isomer.csv` defaults to observed nonzero isomer rows plus per-frame totals; use `--cage-isomer-rows all` to restore the full zero-filled matrix. The `order_parameter` table contains only the selected F3, F4, Q_l, MCG, and DHOP columns; `--order-parameter none` omits it. Focus mean/count columns are written only when `order.focus_waters` is non-empty. Output type `order-tsv` writes only selected per-water F3/F4/Q_l values because MCG/DHOP are frame-level descriptors.

Summary construction records rows, columns, cells, bytes, CSV/XLSX write time, formatting time, and final-save time in `run_config.yaml -> run.summary_write`; the terminal prints its total seconds. Main CSV, XLSX, detail CSV, and `run_config.yaml` are written to same-directory temporary files and atomically replaced on success. XLSX sheets above 200,000 cells or 128 columns keep header styling, filter, freeze pane, and fixed column widths but skip costly body-cell formatting; scientific values and table schemas are unchanged.

The `hydrate_cluster` main-summary table reports the mutually exclusive `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`. Optional cluster-detail CSV records add the corresponding cage-id groups and `boundary_composition`; hydrate-domain CSV records expose only external boundary contacts through `external_boundary_contact_count` and `external_boundary_contact_ids`.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

SQQ-Py cage files include cage waters, `CNT` center pseudoatoms, and assigned guests. SQQ-CPP cage files omit the synthetic `CNT` center pseudoatom. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
