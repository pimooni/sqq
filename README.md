# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

Current release: **0.2.9**

SQQ builds a water network, reports coordination diagnostics, and finds rings, standard half-cages, quasi-cages, closed cages, topology-wide hydrate clusters, per-frame phase domains and boundaries, cage guest occupancy, F3/F4/Q_l order parameters, MCG/DHOP hydrate-nucleation order parameters, and ice-like waters. Detailed algorithms are documented in `docs/design.md`; version notes are documented in `docs/update.md`.

## Changed in 0.2.9

- Hydrate-cluster cage categories are now mutually exclusive: every reported cluster cage belongs to exactly one of a classified sI/sII/sH domain, the generic boundary, ambiguous phase evidence, or residual unclassified topology.
- A boundary cage must be outside every final phase domain and share a complete cage face with at least one domain cage. Only this first external non-phase layer is retained.
- Phase-domain cages are never relabelled as boundary. Direct sI-sII, sI-sH, or sII-sH contacts retain their phase identities on both sides and do not create duplicate boundary membership.
- Removed the phase-boundary, transition, and boundary-context categories, including `sI-boundary`, `sII-boundary`, and `sH-boundary`.
- The `hydrate_cluster` workbook sheet now reports `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`. Optional cluster-detail CSV records expose the corresponding cage-id groups and `boundary_composition`.
- Hydrate-domain detail now names external contacts explicitly as `external_boundary_contact_count` and `external_boundary_contact_ids`.
- Resolved cluster search now forces both `xlsx` and the native `cluster-gro` output type. Turning cluster search off writes no cluster GRO files and removes stale SQQ-generated cluster GRO output from a reused result directory.
- Each frame aggregates all domains/clusters into at most four category files: sI, sII, sH, and boundary. Files contain complete water molecules only, preserve the frame's original wrapped coordinates and box exactly, and omit guests, CNT atoms, ambiguous, unclassified, and isolated cages.
- A real mixed-phase regression partitions its 334-cage main cluster into 66 sI, 194 sII, 69 boundary, 0 ambiguous, and 5 unclassified cages with disjoint cage IDs. Five additional cages remain in isolated/subthreshold components.
- Ring, half-cage, quasi-cage, closed-cage, occupancy, order-parameter, ice, phase-seed, and phase-domain expansion definitions are unchanged. Boundary membership and related cluster reports intentionally change from 0.2.8.
- Package version and root version output are updated to `0.2.9`, released Jul 16, 2026.

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

`-m` / `--mode` selects one of four base presets. The default remains `50`:

| Mode | Purpose | Water graph | Search sizes | Automatic workers | Find cluster |
| --- | --- | --- | --- | --- | --- |
| `00` | Rigorous | Hydrogen bond | 4, 5, 6 | 25% of physical cores | on |
| `09` | Rigorous, high parallelism | Hydrogen bond | 4, 5, 6 | 90% of physical cores | on |
| `50` | Standard (default) | Auto | 5, 6 | 50% of physical cores | off |
| `99` | Performance screening | O-O connectivity | 5, 6 | 90% of physical cores | off |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m 09 -o ./result_rigorous_fast
sqq analyze -i ./gro -m 50 -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_performance
```

Modes do not change `quasi_cage.max_layers`, `order.parameters`, or the initially selected output types; L1, `f3,f4`, and `info,gro,xlsx,summary-detail` remain the defaults in every mode. They do set the initial cluster-search state shown above. After that state is resolved, cluster search on forces `xlsx` and `cluster-gro`, while cluster search off writes no cluster GRO/detail output and leaves ordinary XLSX selection independent. `--find-cluster on|off` explicitly overrides both `hydrate_cluster.enabled` in `config.yaml` and the mode preset. Use `--quasi-max-layer` explicitly for L2/L3 and `--order-parameter` for another descriptor set. Automatic workers use the mode fraction of detected physical cores, reserve one physical core for the system, and are capped by the number of independent GRO/XYZ files or selected trajectory frames. Multiple standalone files use spawned processes by default; a single indexed XTC/TRR trajectory can also distribute frames across spawned workers. `--worker N` / `-w N` overrides the mode percentage; integer text such as `1` or `4` is a worker count, while decimal text such as `0.5` or `1.0` and percentages such as `50%` or `100%` are physical-core fractions. The old `--workers` spelling is retained as a compatibility alias.

SQQ uses process-based parallelism for independent GRO/XYZ files and selected XTC/TRR frames, so CPU-bound ring, quasi-cage, and cage searches can run on multiple cores. The main process alone owns the terminal panel and final workbook; process workers report stage events through a process queue, analyze one file or a small trajectory-frame batch, write frame directories, and return summary rows. At most `3 * workers` process or compatibility-thread tasks are kept in flight; this bounds Future and serialization overhead without reducing the worker count. `parallel.math_threads: 1` prevents nested BLAS/OpenMP oversubscription.

The default `chordless`/`bounded` path preserves the established scientific definitions while accelerating neighbor generation, incremental chord pruning, L1 forward checking, cached layer growth, integer-mask subset ownership, and cage target/edge state pruning. Cage DFS also applies exact remaining-edge incidence and parity conditions before expansion. MDAnalysis supplies orthorhombic cutoff candidates when available, but SQQ still rechecks every distance and hydrogen-bond angle with its established float64 logic. F3 and graph-mode Q_l share one graph-vector cache; all Q_l degrees share candidate lists and spherical-angle work. Optional `ring.definition: shortest_path` applies the Franzblau shortest-path criterion and reuses bounded-BFS distance maps. Optional `quasi_cage.search_policy: exact` preserves distinct frontiers and enumerates connected L2/L3 subsets; these opt-in modes can change or add results. Candidate and state truncation is reported through frame warnings.

Optional scientific cage validation adds PBC-aware face planarity and edge-variation limits, manifold vertex-link checks, positive-volume validation, and volume-centroid cage centers. It remains disabled by default. SQQ uses an orthorhombic box representation and now rejects non-orthogonal/triclinic input explicitly.

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
sqq analyze -i md.gro --output-type info,cage-gro,xlsx -o ./result_sqq_report
sqq analyze -i md.gro --output-type none -o ./result_sqq_config_only
```

`--output-type` defaults to `info,gro,xlsx,summary-detail`. The mandatory `run_config.yaml` remains with every selection, including `none`.

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
  types: [info, gro, xlsx, summary-detail]
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

Process submission uses a bounded rolling queue of at most `3 * workers` tasks. This is a queue-depth limit, not a CPU limit: with 100 effective workers SQQ may keep up to 300 tasks submitted while still running as many as 100 workers concurrently. Results are restored to original file/frame order before workbook writing.

The parent preserves input/frame order in `summary.xlsx`. Different standalone files must have unique case-insensitive stems because each stem is the output frame-directory name. Process runs set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, and `BLIS_NUM_THREADS` to `parallel.math_threads` while workers are spawned, then restore the parent environment.

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

Enabling cluster search always forces both `xlsx` and `cluster-gro`, even when the explicit output list omits them or is `none`. The workbook gains its per-frame `hydrate_cluster` sheet, while native category structures are written under grouped layout as `<frame>/hydrate_cluster/<frame>_cluster_sI.gro`, `<frame>_cluster_sII.gro`, `<frame>_cluster_sH.gro`, and `<frame>_cluster_boundary.gro`. Flat layout places the same filenames directly in the frame directory. All same-category domains and clusters are aggregated into one file per frame. An absent category is omitted unless `output.write_empty_files: true`.

Cluster GRO files contain only complete water molecules belonging to the selected cage IDs; guests and CNT atoms are excluded. Ambiguous, unclassified, and isolated cages are not exported. Every atom keeps the exact wrapped coordinate from the analyzed frame, and every file keeps the original box; categories are never moved or unwrapped independently. Periodic or percolating networks may therefore still show bonds crossing a box face because no single-copy GRO representation can remove every periodic seam.

Cage IDs are mutually exclusive across sI, sII, sH, and boundary, but adjacent category files can contain the same face-water molecules because neighboring cages physically share them. Cluster data is not copied into per-frame `*_info.md` reports. Add `cluster-detail` to `--output-type` for `summary_detail/hydrate_domain.csv` and one-row-per-cluster `summary_detail/hydrate_cluster_detail.csv`. Explicit `cluster-detail` or `cluster-gro` selection requires cluster search. Turning search off writes neither `cluster-detail` nor `cluster-gro` and removes stale generated cluster GRO files. Public motif output is not generated.

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
| `-m, --mode MODE` | `00`, `09`, `50`, `99`; default `50` | Select rigorous, rigorous/high-parallelism, standard, or performance preset |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --size SIZES` | Comma-separated subset of `4,5,6,7` | Set ring and quasi-cage search sizes; cage search uses the selected `4,5,6` sizes |
| `--ring-size SIZES` | `auto` or a comma-separated subset of `--size` | Report only these searched ring sizes |
| `--cage-size GROUPS` | `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, `I2II`; groups may be comma-separated | Restrict cage reporting; default `auto` follows `--size` |
| `--max-cage-face N` | Positive integer; default `20` | Limit generated cage search compositions |
| `--cage-fast-closure VALUE` | `on`, `off`; default `on` | Enable indexed two-to-four half-cage closure after generic grow |
| `--cage-scientific-validation VALUE` | `on`, `off`; default `off` | Enable strict face/manifold/volume validation and volume centroids |
| `--find-cluster VALUE` | `on`, `off`; default follows mode/config | Override all-detected-cage cluster search; `on` forces XLSX plus native cluster GRO output |
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
| `--output-type TYPES` | Comma-separated `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `xlsx`, `summary-detail`, `cluster-detail`, `all`, or `none` | Select the complete output set; default `info,gro,xlsx,summary-detail`. `cluster-gro` requires cluster search and is forced when search is on |
| `--cage-isomer-rows MODE` | `nonzero`, `all`; default `nonzero` | Choose whether `summary_detail/cage_isomer.csv` keeps only observed isomer rows or the full zero-filled matrix |

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

`--output-type TYPE[,TYPE...]` replaces the complete configured output list. Its default is `info,gro,xlsx,summary-detail`, equivalent to:

```yaml
output:
  types: [info, gro, xlsx, summary-detail]
```

The canonical names are `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `xlsx`, `summary-detail`, and `cluster-detail`. `gro` expands to the five ordinary ring/half/quasi/cage/ice GRO categories; `cluster-gro` is a separate search-dependent type. `all` selects every output applicable to the enabled analyses; `none` selects no optional output before mandatory cluster outputs are resolved. Both keywords must appear alone. `run_config.yaml` is mandatory and is always rewritten.

`summary-detail` selects the ordinary UTF-8-SIG CSV tables. `cluster-detail` separately selects `hydrate_domain.csv` and `hydrate_cluster_detail.csv`; explicit `cluster-detail` and `cluster-gro` selections require cluster search. Enabling cluster search forces both `xlsx` and `cluster-gro`, so `--find-cluster on --output-type none` still writes `summary.xlsx`, the per-frame category GRO files, and `run_config.yaml`. Disabling search writes no cluster GRO files and cleans stale SQQ-generated `hydrate_cluster` directories or flat cluster filenames.

Explicit CLI selection has precedence:

```text
--output-type > output.types > default info,gro,xlsx,summary-detail
```

`--hydrate-cluster`, `--cluster-detail`, `--no-output`, `--write-order-tsv`, and the individual `--no-*` output switches were removed in 0.2.8. They have no compatibility aliases. Likewise, `output.disabled_outputs` is rejected rather than migrated; configurations must use `output.types`.

When an existing output directory is reused, SQQ removes known stale files for output types outside the effective selection while preserving unrelated user files. If no per-frame output type is selected and no unrelated file remains, the empty frame directory is removed.

## Output Structure

With the default `--output-type info,gro,xlsx,summary-detail`, SQQ writes one folder per frame, a global workbook, and ordinary CSV detail tables:

```text
result_sqq/
  summary.xlsx
  summary_detail/
    failures.csv                       # only when frames fail
    cage_occupancy.csv
    cage_isomer.csv
    quasi_cage_isomer.csv
    hydrate_domain.csv                # only with cluster-detail and cluster search
    hydrate_cluster_detail.csv        # only with cluster-detail and cluster search
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

Only selected or required outputs are written. With cluster search off, `--output-type none` leaves only `run_config.yaml` and does not create empty per-frame directories. With cluster search on, the same selector still forces `summary.xlsx` and grouped or flat cluster GRO files. Grouped layout uses `<frame>/hydrate_cluster/`; flat layout places the four canonical filenames at the frame root. Missing categories are omitted unless `output.write_empty_files` is true.

Without `--strict`, standalone serial/process/thread read failures become failed summary rows and analysis continues where the reader remains usable. Failed inputs appear in `summary.xlsx/failures` and `summary_detail/failures.csv` when those outputs are enabled, and always in the mandatory `run_config.yaml` `run.failures` list. With `--strict`, SQQ re-raises the error after updating `run_config.yaml` to `status: failed`.

GRO structure folders, filenames, and title lines use portable ASCII structure labels since version 0.2.4, for example `5^126^2` and `qc_5r_5^36^2_56566`. Markdown/Excel scientific labels retain their readable superscript notation. This avoids Windows GBK/legacy-reader failures caused by Unicode superscript or subscript characters in generated GRO paths and titles.

Each per-frame `*_info.md` report is arranged for inspection. `Frame Information` begins with the SQQ version, report-generation date and local timezone, absolute source path, frame name, and trajectory time. The report shows requested/effective `graph_mode`, the effective `bond_mode`, only reported ring sizes, final free-ring counts, the active network degree distribution, groups half-cage and quasi-cage isomers below composition totals, and keeps cage composition totals plus cage isomers in one vertical `Cage` table. Cluster, domain, and boundary data is kept out of per-frame info reports: compact counts go to the forced XLSX cluster summary, exact records go to optional `cluster-detail` CSV files, and native category structures go to forced `cluster-gro` files. Internal `hc_` and `qc_` prefixes are omitted from report labels.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary.xlsx` remains plotting-oriented when `xlsx` is enabled or cluster search forces it. Its first sheet is a dashboard: Configuration includes `SQQ version`, requested/effective `Graph mode` such as `auto -> hbond`, normalized `Order parameters`, `Find cluster`, and normalized `Output types`; `Analysis Results (min / mean / max)` reports per-frame min/mean/max values while `Frames total / ok / failed` stays a run-level count. Failed inputs add a compact `failures` sheet without restoring the removed redundant `frame` sheet. The remaining analysis sheets keep one input file or trajectory frame per row, including connection diagnostics, `ring`, `half_cage`, compact composition-level `quasi_cage`, `cage`, the mandatory `hydrate_cluster` sheet when cluster search is enabled, optional `order_parameter`, `ice`, `detail_index`, and `config`. Ordinary multi-row and isomer tables are written as UTF-8-SIG CSV files in `summary_detail/` when `summary-detail` is selected: optional `failures.csv`, `cage_occupancy.csv`, `cage_isomer.csv`, and `quasi_cage_isomer.csv`. The separate `cluster-detail` type writes `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. The `quasi_cage` workbook sheet aggregates exact quasi-cage isomers into composition-level columns such as `5r_5²6³`, while `quasi_cage_isomer.csv` keeps nonzero exact isomer rows with `quasi_cage_type`, `isomer`, and `count`. `cage_isomer.csv` defaults to observed nonzero isomer rows plus per-frame totals; use `--cage-isomer-rows all` to restore the full zero-filled matrix. The `order_parameter` sheet contains only the selected F3, F4, Q_l, MCG, and DHOP columns; `--order-parameter none` omits the sheet. Focus mean/count columns are written only when `order.focus_waters` is non-empty. Output type `order-tsv` writes only selected per-water F3/F4/Q_l values because MCG/DHOP are frame-level descriptors.

Summary construction records rows, columns, cells, bytes, CSV/XLSX write time, formatting time, and final-save time in `run_config.yaml -> run.summary_write`; the terminal prints its total seconds. XLSX, detail CSV, and `run_config.yaml` are written to same-directory temporary files and atomically replaced on success. Data sheets above 200,000 cells or 128 columns keep header styling, filter, freeze pane, and fixed column widths but skip costly body-cell formatting; scientific values and table schemas are unchanged.

The `hydrate_cluster` sheet reports the mutually exclusive `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`. Optional cluster-detail CSV records add the corresponding cage-id groups and `boundary_composition`; hydrate-domain CSV records expose only external boundary contacts through `external_boundary_contact_count` and `external_boundary_contact_ids`.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

Cage files include cage waters, CNT center atoms, and assigned guests. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
