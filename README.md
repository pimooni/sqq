# SQQ

**SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.**

Current release: **0.3.6**

SQQ provides the complete SQQ-Py water-shell topology workflow plus the focused SQQ-CPP cage engine. Modes `00` and `py` use SQQ-Py; modes `99` and `cpp` use C++17 for graph, internal ring, cage, occupancy, and F3/F4 analysis. Algorithms are documented in `docs/design.md` and release notes in `docs/update.md`.

## Acknowledgements

Names are listed alphabetically by family name. The order does not indicate relative contribution.

- Liwei Cheng - Wuhan Institute of Technology
- Bin Fang - Hainan University
- Jihui Jia - China University of Petroleum (Beijing)
- Wuquan Li - Beijing Huairou Laboratory
- Bo Liao - China University of Petroleum (East China)
- Yingxu Lu - Wuhan Institute of Technology
- Fengyi Mi - Southwest University of Science and Technology
- Zhengcai Zhang - Laoshan Laboratory

## Changed in 0.3.6

- Package and native-core versions are synchronized at `0.3.6`, released Jul 23, 2026.
- The VMD `sqq show` command accepts one or more family/target groups, so one command can combine objects, for example `sqq show cage 512 guest 512`.
- The source default remains `sqq show cage all`. The first `show` replaces it; later `show` commands add layers, and exact repeated selections are ignored.
- `sqq clear` removes custom layers and color overrides, restores the default cage-all view, and rearms first-show replacement.
- Cross-family representations use the fixed order `phase -> cluster -> domain -> cage -> guest`, independent of command order. The existing cage-topology priority remains separate.
- Sourcing the Tcl script prints a compact command welcome. `sqq help`, `sqq -h`, and `sqq --help` print the full command guide. The `color` command remains single-family.

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

For multiple GRO files, SQQ groups compatible frames automatically. Files with one topology share the requested output root; heterogeneous inputs are separated into `result_A`, `result_B`, and so on in first-occurrence order. The grouping affects aggregation and paths only, not per-frame analysis.

Glob pattern:

```bash
sqq analyze -i "./gro/*.gro" -o ./result_sqq
```

XTC/TRR trajectory with a topology file:

```bash
sqq analyze -i traj.xtc --top topol.gro -c config.yaml -o ./result_sqq
```

LAMMPS dump or DCD with a DATA topology; standard water/methane types are inferred automatically:

```bash
sqq analyze -i traj.lammpstrj -t system.data -o ./result_sqq
```

### Input Units and Boxes

GRO and MDAnalysis trajectory coordinates are interpreted in nm. GRO accepts exactly one frame per file and rejects truncated atom blocks, missing or malformed box lines, extra non-empty records, and non-finite coordinates. Trajectory frames also require finite coordinates. XYZ coordinates are multiplied by `input.xyz_scale` / `--xyz-scale`; the default `0.1` assumes angstrom input, while `1.0` keeps nm values. SQQ accepts exactly one declared XYZ frame per file and rejects truncated, extra, malformed, or non-finite atom records. XYZ has no periodic box unless converted through another format.

GRO atom counts and the mandatory box line are validated. A three-value positive box is orthorhombic; an all-zero box is treated as non-periodic. Nine-value GRO boxes with nonzero tilt terms and trajectory frames with non-90-degree angles are rejected because triclinic minimum-image calculations are not implemented. GRO molecules are formed from contiguous residue blocks in source order, preventing wrapped or repeated residue IDs from merging distinct molecules. LAMMPS normally uses DATA molecule IDs; automatic inference can rebuild them from unambiguous Bonds components, and dump atom rows may be interleaved.

LAMMPS trajectories require `-t system.data` (equivalent to `--top` / `--topology`). A non-empty `input.lammps.type_map` explicitly maps every numeric atom type to `resname`/`atomname` or `ignore` and always takes priority. If the map is absent or empty, SQQ uses DATA masses, type comments, and Bonds to infer only unambiguous standard water (`1 O + 2 H`, two O-H bonds), all-atom methane (`1 C + 4 H`, four C-H bonds), and labeled single-site methane mappings. If molecule IDs do not define valid molecules but Bonds do, SQQ rebuilds deterministic molecule IDs and reports that decision. Ambiguous masses, inconsistent reuse of one numeric type, unsupported molecular topology, or insufficient topology evidence fail with a request for an explicit map. The resolved mapping is recorded in `config.yaml`, per-frame info, and main-summary configuration. This normalization is shared by SQQ-Py and SQQ-CPP. Supported inputs are LAMMPS DATA with `full`, `molecular`, `bond`, or `angle` atom style, fully periodic `pp pp pp` orthorhombic dump boxes, and LAMMPS DCD. Tilted boxes, nonperiodic dump boundaries, `units lj`, duplicate atom IDs, ambiguous molecule reconstruction, and topology/trajectory ID mismatches fail before analysis. `input.trajectory_stride` applies to XTC, TRR, LAMMPS dump, and LAMMPS DCD.

## Analysis Modes

`-m` / `--mode` selects one of four presets; the default is `py`:

| Mode | Engine | Water graph | Ring search | Default workers | Find cluster | Default output types |
| --- | --- | --- | --- | --- | --- | --- |
| `00` | SQQ-Py | `hbond` | 4/5/6 | 100% | on | `info,sqq-cage-gro,sqq-render,summary-xlsx` |
| `py` | SQQ-Py | `auto` | 4/5/6 | 1 worker | off | `info,sqq-cage-gro,sqq-render,summary-xlsx` |
| `99` | SQQ-CPP | `hbond` | internal 4/5/6 | 100% | unsupported | `info,sqq-cage-gro,sqq-render,summary-csv` |
| `cpp` | SQQ-CPP | `auto` | internal 4/5/6 | 1 worker | unsupported | `info,sqq-cage-gro,sqq-render,summary-csv` |

```bash
sqq analyze -i ./gro -m 00 -o ./result_rigorous
sqq analyze -i ./gro -m py -o ./result_standard
sqq analyze -i ./gro -m 99 -o ./result_cpp_hbond
sqq analyze -i ./gro -m cpp -o ./result_cpp_auto
```

For SQQ-Py, `--find-cluster` overrides YAML, which overrides the mode preset. Mode `00` enables cluster search, while mode `py` leaves it off. Search results enter selected info/main-summary outputs, but no mode selects `cluster-gro` by default; request it explicitly when split cluster structures are required. Modes `99` and `cpp` reject `--find-cluster on`.

Modes `py` and `cpp` default to one worker. Modes `00` and `99` use 100% of detected physical cores, reserve one physical core, and are capped by independent files or selected trajectory frames. `-w` / `--worker` overrides the preset: integer text is a worker count, while `0.5`, `1.0`, `50%`, and `100%` are physical-core fractions. Process parallelism supports independent GRO/XYZ files and indexed XTC/TRR/LAMMPS trajectories. At most `3 * workers` tasks are submitted at once.

The default `chordless`/`bounded` path preserves the established scientific definitions while accelerating neighbor generation, incremental chord pruning, L1 forward checking, cached layer growth, integer-mask subset ownership, and cage target/edge state pruning. Cage DFS also applies exact remaining-edge incidence and parity conditions before expansion. MDAnalysis supplies orthorhombic cutoff candidates when available, but SQQ still rechecks every distance and hydrogen-bond angle with its established float64 logic. F3 and graph-mode Q_l share one graph-vector cache; all Q_l degrees share candidate lists and spherical-angle work. Optional `ring.definition: shortest_path` applies the Franzblau shortest-path criterion and reuses bounded-BFS distance maps. Optional `quasi_cage.search_policy: exact` preserves distinct frontiers and enumerates connected L2/L3 subsets; these opt-in modes can change or add results. Candidate and state truncation is reported through frame warnings.

Every cage now passes the same mandatory topology validation in SQQ-Py and SQQ-CPP: each edge belongs to exactly two faces, `V - E + F = 2`, the face shell is connected, every vertex link is one cycle, and every shell vertex is trivalent. Optional scientific cage validation adds PBC-aware face-planarity and edge-variation limits, nonzero projected area, positive-volume validation, and volume-centroid cage centers. It remains disabled by default, but disabling it no longer bypasses topology validation. SQQ uses an orthorhombic box representation and rejects non-orthogonal/triclinic input explicitly.

The current release uses the same compact three-row stage model for serial and parallel progress: file preparation (`reading`, `settings`, `selecting`), core topology search (`graph`, `ring`, `half/quasi`, `cage`, and optional `cluster`), and post-processing (`filtering`, `order`, `ice`, `output`). In interactive single-file runs, the active stage is highlighted with bold bright-blue ANSI text. The `cluster` stage appears only when hydrate-cluster analysis is enabled. Parallel runs also show aggregate stage counts and up to six active files with per-stage and per-file timings.

### Native SQQ-CPP Backend

Modes `99` and `cpp` select the focused native engine. Python owns input normalization, molecule selection, scheduling, annotated GRO/VMD output, Markdown, summary CSV, and optional XLSX; C++17 performs graph construction, internal chordless 4/5/6 rings, cage topology/isomers, occupancy, and F3/F4 while releasing the GIL.

| Mode | Graph | Workers | Default output |
| --- | --- | --- | --- |
| `99` | `hbond` | 100%, one physical core reserved | `info,sqq-cage-gro,sqq-render,summary-csv` |
| `cpp` | `auto` | 1 worker | `info,sqq-cage-gro,sqq-render,summary-csv` |

Both accept orthorhombic GROMACS/LAMMPS inputs, compatible graph/pair settings, `-s` within 4/5/6, cage report/validation settings, `f3`/`f4`, process or serial scheduling, and `info`, `gro`, `cage-gro`, `sqq-cage-gro`, `sqq-render`, `summary-csv`, or `summary-xlsx`. `sqq-render` implies `sqq-cage-gro`. `gro` enables the supported classified cage GRO output, but `cpp` does not select it by default.

Unsupported requests fail before analysis: public ring output, size 7, shortest-path rings, half/quasi cages, cluster, ice, Q_l/MCG/DHOP, membership/order TSV, legacy per-frame `vmd`, detail CSV, Python fast closure, thread scheduling, and triclinic boxes. A failed native extension never falls back to Python.

The mode-`cpp` default layout is:

```text
result/
  sqq-cage.gro
  sqq-render.vmd.tcl
  summary/
    summary.csv
    cage.csv
    cage_occupancy.csv
    cage_isomer.csv
    order_parameter.csv
  config.yaml
  frame_name/
    frame_name_info.md
```

Neither native mode selects ordinary/classified GRO by default. Modes `99` and `cpp` write classified cage GRO only when `gro` or `cage-gro` is selected explicitly.

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

Default outputs come from the selected mode. Modes `00` and `py` use `info,sqq-cage-gro,sqq-render,summary-xlsx`; modes `99` and `cpp` use `info,sqq-cage-gro,sqq-render,summary-csv`. No preset includes ordinary, classified, or cluster GRO. `config.yaml` remains mandatory.

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
mode: "py"

input:
  pattern: "*.gro"
  trajectory_stride: 1
  xyz_scale: 0.1
  lammps:
    units: real
    timestep: 1.0
    atom_style: full
    coordinate_convention: auto
    type_map: {}  # optional override; empty enables strict DATA inference

graph:
  bond_mode: auto
  oo_cutoff_nm: 0.35
  hbond_distance_nm: 0.35
  hbond_angle_deg: 30.0

ring:
  sizes: [4, 5, 6]
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

guest:
  resnames: [CH4, CO2, MET, ETH]
  center_atoms:
    CH4: [C]
    CO2: [C]
    MET: [C]
  center_mode: center_atom

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
  types: [info, sqq-cage-gro, sqq-render, summary-xlsx]
  summary_csv_dir: summary
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

Before dispatching two or more GRO files, SQQ reads only their topology records and assigns topology groups in first-occurrence order. The fingerprint contains the atom count and ordered contiguous residue blocks, represented by each block's residue name and ordered atom-name sequence. Titles and time labels, coordinates, velocities, boxes, and numeric atom/residue IDs do not affect grouping. A supplied GRO `--top` is checked against every input fingerprint; any mismatch fails before analysis and identifies the exact source file.

All accepted groups use one shared worker pool and one global progress index. Each task also carries a group-local frame index and output root, so group summaries and annotated bundles remain correctly ordered without running groups serially. Requested `graph.bond_mode: auto` remains recorded as `auto`, but its effective `hbond` or `oo` mode is resolved once from a representative frame in each topology group and reused by both SQQ-Py and SQQ-CPP for every frame in that group.

With `parallel.workers: auto`, modes `py` and `cpp` resolve to one worker, while modes `00` and `99` use 100% of detected physical cores before reserving one physical core and applying the file/frame cap. Physical-core detection prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker` / `-w` accepts either a fraction (`50%`, `0.5`, or `1.0` for 100%) or an explicit positive integer worker count (`1` means one worker). Windows `ProcessPoolExecutor` runs are capped at 61 workers; Linux workstations can use larger explicit values such as `-w 100`, subject to the reserve-one-core rule, task count, memory, and storage throughput.

One XTC/TRR or supported LAMMPS trajectory with `--top` is frame-parallel when the process backend resolves to more than one worker. Every worker opens a private MDAnalysis Universe once and seeks small contiguous batches of selected raw frame indexes; batch size is automatically bounded from 1 to 8, and complete coordinate arrays are not serialized between processes. Parent and worker trajectory readers are explicitly closed. Multiple trajectory files and the compatibility thread backend use the serial trajectory reader.

Process submission uses a bounded rolling queue of at most `3 * workers` tasks. This is a queue-depth limit, not a CPU limit: with 100 effective workers SQQ may keep up to 300 tasks submitted while still running as many as 100 workers concurrently. Results are restored to original file/frame order before main-summary writing.

The parent preserves original input order globally and group-local order in every selected group summary and annotated bundle. Output-name collisions are resolved deterministically within each topology group. Process runs set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, and `BLIS_NUM_THREADS` to `parallel.math_threads` while workers are spawned, then restore the parent environment.

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

Topology validation is always enabled. Every candidate must use each edge exactly twice, satisfy `V - E + F = 2`, form one edge-connected face shell, have one cyclic face link around every vertex, and have only trivalent shell vertices. These inexpensive checks reject disconnected, pinched, branched, and non-manifold false cages before type/isomer assignment in both engines.

`cage.scientific_validation: false` is the default. When enabled with `--cage-scientific-validation on`, a topologically valid cage must additionally satisfy the configured PBC-aware face-planarity RMS and edge-length coefficient-of-variation limits, nonzero projected face area, and positive minimum triangulated volume. Accepted cages then use the volume centroid instead of the mean cage-water position. Enabling it can therefore remove geometrically distorted cages and can change guest occupancy or geometry-resolved hydrate-cluster edges. The mandatory topology checks can reduce cage, isomer, occupancy, and cluster results relative to earlier 0.3.2 builds that accepted non-manifold shells. Raw ring and half/quasi searches remain unchanged; ownership-filtered free-ring and free-patch outputs can increase when a rejected cage no longer consumes them.

Guest occupancy uses the configured center atom when available. The defaults select `CH4`, `CO2`, `MET`, and `ETH` as guests and map `CH4`, `CO2`, and `MET` to atom name `C`, so these residues use their carbon atom under the default `guest.center_mode: center_atom`. Otherwise, guest atoms are PBC-unwrapped around one molecular anchor before calculating the centroid; the same helper is used by MCG. This correction can intentionally change occupancy counts relative to 0.2.6 or early 0.2.7 results for multi-atom guests crossing a periodic boundary.

## Hydrate Cluster

`--find-cluster on` analyzes every detected cage in the selected search scope. Cages become graph nodes and are connected through complete shared ring faces. When several detected cages reference the same face, ring-plane geometry keeps at most one cage on each physical side. `--cage-size` filters user-facing cage tables and files only; it does not remove cages from cluster connectivity or phase evidence.

The hierarchy follows the HTR+ idea of classifying hydrate type, domains, and boundaries on a cage-connection graph ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)). SQQ implements this independently with labelled shared-face fingerprints, strict local seeds, mutually compatible expansion, and exclusive per-frame domains.

`--cluster-min-cage N` sets the minimum connected-component size; the default is `2`. Smaller components are counted as isolated cages.

Within each cluster, SQQ builds labelled first-shell fingerprints from neighboring cage types and shared-face sizes. Strict local sI/sII/sH seeds initialize phase evidence. The sH templates cover `5^12`, `4^3 5^6 6^3`, and `5^12 6^8` cages; the earlier two-anchor sH composite is retained as supplemental high-confidence evidence. All three phases expand through mutually compatible face-labelled edges when a candidate has at least two accepted phase contacts. Cages claimed exclusively by one phase form deterministic per-frame domains.

After the exclusive sI/sII/sH domains are finalized, SQQ partitions the remaining cluster cages. A cage enters the generic boundary only when it is outside every phase domain and directly shares a complete cage face with at least one domain cage. Boundary search stops at this first external non-phase layer. Domain cages are never relabelled as boundary, and a direct shared-face contact between different phase domains leaves both endpoint cages in their original phases.

The resulting `classified_cage_ids`, `boundary_cage_ids`, `ambiguous_cage_ids`, and `unclassified_cage_ids` are mutually exclusive and together cover every cage in a reported cluster. Competing phase claims without boundary membership remain ambiguous; all other residual cages are unclassified. There are no `sI-boundary`, `sII-boundary`, `sH-boundary`, transition, or boundary-context categories. Neighboring cages can still share face-water coordinates in structure views, so cage ownership should be verified from cage IDs or detected cage/ring edges rather than coordinate-set overlap.

The default command uses mode `py`, so cluster search is off unless enabled by mode `00`, `hydrate_cluster.enabled`, or explicit `--find-cluster on`. Modes `99` and `cpp` do not support cluster search. Explicit `--find-cluster on|off` has highest priority. Cluster search does not alter ring, patch, cage, occupancy, order-parameter, or ice results. Classification is per-frame and independent of the cage reporting filter; temporal grain tracking and crystallographic orientation matching are not implemented.

Cluster search populates every selected `info` and main-summary output. Split category structures are written only when `cluster-gro` is selected explicitly; no mode includes it by default. The selected main summary output gains its per-frame `hydrate_cluster` table, while native category structures are written under grouped layout as `<frame>/hydrate_cluster/<frame>_cluster_sI.gro`, `<frame>_cluster_sII.gro`, `<frame>_cluster_sH.gro`, and `<frame>_cluster_boundary.gro`. Flat layout places the same filenames directly in the frame directory. All same-category domains and clusters are aggregated into one file per frame. An absent category is omitted unless `output.write_empty_files: true`.

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
| `-i, --input INPUT` | `.gro`, `.xyz`, `.xtc`, `.trr`, `.dump`, `.lammpstrj`, or LAMMPS `.dcd`; directory; or glob | Input source |
| `-c, --config FILE` | YAML or JSON file | User configuration file |
| `-o, --output DIR` | Directory path; default `result_sqq` | Output directory |
| `-m, --mode MODE` | `00`, `py`, `99`, `cpp`; default `py` | Select SQQ-Py (`00`, `py`) or SQQ-CPP (`99`, `cpp`) |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, `pairs` | Override the water-graph connection mode |
| `-s, --size SIZES` | Comma-separated subset of `4,5,6,7` | Set ring and quasi-cage search sizes; cage search uses the selected `4,5,6` sizes |
| `--ring-size SIZES` | `auto` or a comma-separated subset of `--size` | Report only these searched ring sizes |
| `--cage-size GROUPS` | `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, `I2II`; groups may be comma-separated | Restrict cage reporting; default `auto` follows `--size` |
| `--max-cage-face N` | Positive integer; default `20` | Limit generated cage search compositions |
| `--cage-fast-closure VALUE` | `on`, `off`; default `on` | Enable indexed two-to-four half-cage closure after generic grow |
| `--cage-scientific-validation VALUE` | `on`, `off`; default `off` | Enable additional geometric face/volume validation and volume centroids; basic topology validation is always on |
| `--find-cluster VALUE` | `on`, `off` | Override SQQ-Py cluster search; split GRO still requires `cluster-gro` |
| `--cluster-min-cage N` | Positive integer; default `2` | Minimum connected cage count required for one hydrate_cluster |
| `--pattern PATTERN` | Glob; default `*.gro` | Select files when `--input` is a directory |
| `-t, --top, --topology FILE` | GRO for XTC/TRR or multi-GRO validation; LAMMPS DATA for dump/DCD | Supply trajectory topology or validate a multiple-GRO batch |
| `--xyz-scale SCALE` | Positive float; default `0.1` | Multiply XYZ coordinates by this value to obtain nm; use `1.0` for XYZ already in nm |
| `--trajectory-stride N` | Positive integer; default `1` | Read every Nth trajectory frame |
| `--lammps-units STYLE` | `real`, `metal`, `nano` | Select LAMMPS units |
| `--lammps-timestep DT` | Positive number | Convert LAMMPS steps to ps |
| `--lammps-atom-style STYLE` | `full`, `molecular`, `bond`, `angle` | Interpret LAMMPS DATA |
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
| `--output-type TYPES` | Established outputs plus `sqq-cage-gro` and `sqq-render`; comma-separated, or `all`/`none` | Replace the complete mode-specific output set |
| `--cage-isomer-rows MODE` | `nonzero`, `all`; default `nonzero` | SQQ-Py: control `summary_detail/cage_isomer.csv` when `summary-detail-csv` is selected; SQQ-CPP: control `summary/cage_isomer.csv` and the optional `summary-xlsx` sheet |

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

`--output-type TYPE[,TYPE...]` replaces the configured list. Mode `py` defaults to:

```yaml
output:
  types: [info, sqq-cage-gro, sqq-render, summary-xlsx]
```

SQQ-Py accepts the established info/TSV/per-frame GRO/summary/detail types plus `sqq-cage-gro` and `sqq-render`. `gro` expands to ordinary ring/half/quasi/cage/ice GRO categories; `sqq-cage-gro` is a separate run-level annotated trajectory. `sqq-render` implies `sqq-cage-gro`.

SQQ-CPP accepts `info`, `gro`, `cage-gro`, `sqq-cage-gro`, `sqq-render`, `summary-csv`, and `summary-xlsx`, plus `all`/`none`. Neither `99` nor `cpp` selects `gro` or `cage-gro` by default. Cluster search fills selected info/main-summary outputs but does not add unselected output types; explicit `cluster-gro` or `cluster-detail` requires search. `config.yaml` is always written.

Explicit CLI selection has precedence:

```text
--output-type > output.types > engine default
```

`--hydrate-cluster`, `--cluster-detail`, `--no-output`, `--write-order-tsv`, and the individual `--no-*` output switches were removed in 0.2.8. They have no compatibility aliases. Likewise, `output.disabled_outputs` is rejected rather than migrated; configurations must use `output.types`.

When an existing output directory is reused, SQQ removes known stale files for output types outside the effective selection while preserving unrelated user files. Main-summary and detail-CSV cleanup is restricted to known SQQ-generated filenames inside the currently configured `summary_csv_dir` and `summary_detail_dir`; unknown files are preserved, and changing either setting does not make SQQ scan or clean a formerly configured directory. If no per-frame output type is selected and no unrelated file remains, the empty frame directory is removed.

## Output Structure

A single GRO file and trajectory inputs retain their established per-frame layout. For two or more GRO files, topology grouping controls only the aggregation root. If every GRO has one compatible topology, all selected outputs are written directly under the requested result directory:

```text
result/
  config.yaml
  summary.xlsx                 # when summary-xlsx is selected
  summary/                     # when summary-csv is selected
    summary.csv
    cage.csv
    ...
  info/
    frame_001_info.md
    frame_002_info.md
  gro/                         # only when ordinary GRO output is selected
    frame_001/
    frame_002/
  sqq-cage.gro                # when sqq-cage-gro is selected
  sqq-render.vmd.tcl          # when sqq-render is selected
```

When 2-26 distinct topologies are found, groups are assigned letters by first occurrence and each group gets a complete independent result root. No summary, GRO, or VMD bundle combines incompatible systems:

```text
result/
  config.yaml              # batch manifest and source-to-group mapping
  result_A/
    config.yaml
    summary.xlsx               # and/or summary/
    info/
    gro/                       # when selected
    sqq-cage.gro               # when selected
    sqq-render.vmd.tcl         # when selected
  result_B/
    config.yaml
    summary.xlsx               # and/or summary/
    info/
    gro/                       # when selected
    sqq-cage.gro               # when selected
    sqq-render.vmd.tcl         # when selected
```

If more than 26 topologies are found, SQQ warns and switches the whole multi-GRO run to information-only output. It still analyzes every readable GRO, but writes only the root `config.yaml` and `result/info/*_info.md`; summary XLSX/CSV/detail files, ordinary GRO files, annotated cage GRO, and VMD renderer files are all suppressed. This safety override has precedence over mode defaults and explicit output requests.

For normal multiple-GRO groups, Markdown, membership/order TSV, and legacy per-frame VMD reports are placed under `info/`; ordinary structure files are placed under `gro/<frame>/`.

`sqq-cage.gro` concatenates one complete block per successful frame; atom order and the named contiguous residue-block topology must match, while numeric atom/residue IDs may differ. Coordinates and boxes remain wrapped exactly as analyzed. Annotations begin at column 69 as `; SQQ1 m=...` and retain water-to-cage plus guest-to-cage memberships; nonmembers use `m=-`. A guest assigned to several cages keeps every assignment, and every atom of a multi-atom guest remains available to the renderer.

Keep `sqq-render.vmd.tcl` and `sqq-cage.gro` together, then source the script from the VMD Tk Console:

```tcl
source {path/to/sqq-render.vmd.tcl}
```

Sourcing prints a compact welcome, reports `SQQ graph: <effective-mode>` once, and starts from the default `sqq show cage all` view. The graph line is printed again only if the effective mode changes. Use any of these equivalent commands for the full guide:

```tcl
sqq help
sqq -h
sqq --help
```

The command grammar is explicit:

```text
sqq show <family> <target...> [<family> <target...>]...
sqq color <family> <target...> <color>
sqq clear
```

Supported families are `cage`, `guest`, `phase`, `cluster`, and `domain`. Examples:

```tcl
sqq show cage all
sqq show cage 512
sqq show cage 512 51264
sqq show cage 51262_00053
sqq show cage 512 guest 512
sqq show cage 512 51264 guest 512 phase sI

sqq show guest all
sqq show guest 512
sqq show guest 51262_00053

sqq show phase all
sqq show phase sI boundary
sqq show cluster all
sqq show cluster cluster_00001
sqq show domain all
sqq show domain domain_00001

sqq color cage 512 green
sqq color cage 51262_00053 yellow
sqq color guest 512 yellow
sqq color phase boundary orange
sqq color cluster cluster_00001 cyan
sqq color cage all default
```

The startup `sqq show cage all` view is a replaceable default. The first `sqq show ...` command after sourcing the script or after `sqq clear` replaces that default; later `show` commands add independent layers without removing earlier selections. One `show` may contain several family/target groups, and an exact repeated family/target selection is ignored rather than creating another VMD representation. `sqq clear` removes all custom show layers and color overrides, restores the initial cage-all view, and makes the next `show` replace that restored default.

Each family token in `show` starts a new group and consumes the following targets until the next family token. For `cage`, a target is `all`, a registered cage type, or an exact frame-local cage ID; generic types such as `4^1-5^10-6^2` also accept `4151062`. For `guest`, the same target identifies guests assigned to all cages, to a cage type, or to one exact cage ID. Phase targets are `all`, `sI`, `sII`, `sH`, `boundary`, `ambiguous`, `unclassified`, or `isolated`; cluster/domain targets are `all` or exact IDs. Multiple targets are accepted within each family group. The former inferred forms such as `sqq show 512` and `sqq color 512 blue` are not accepted.

Unlike `show`, `sqq color` accepts exactly one family per command. Colors accept a case-insensitive VMD color name, an in-range ColorID, or `default`. Cage and guest overrides are independent and persist across frame/selection changes until `sqq clear`, re-sourcing, or an explicit `default` reset. Cross-family layers always render as `phase -> cluster -> domain -> cage -> guest`, so guests remain last and visible regardless of `show` order. This family order is separate from the fixed cage-topology priority used for coincident cage edges and multi-cage guests. Cage networks use DynamicBonds with a 3.5 angstrom cutoff; guests use CPK and include the full molecule. A single cage layer uses a 0.125 angstrom cylinder radius (0.250 angstrom diameter); multi-type layers remain bounded from 0.125 to 0.130 angstrom.

The renderer manages representations by VMD's stable representation names, so `show`, `color`, and frame changes remove only SQQ-created representations and preserve representations added by the user. Rapid frame notifications are coalesced into one pending redraw. Fully unknown cage, cage-ID, guest-selection, cluster-ID, and domain-ID targets are rejected against the complete loaded trajectory; recognized phase names remain valid even when the current frame has no matching membership. Re-sourcing a generated script resets its selection/color state.

Exact cage, cluster, and domain IDs are assigned independently in each frame. Retaining an ID selection while changing frames therefore selects the same frame-local label, not a tracked physical object. Category selections (`phase`, `cluster`, or `domain`) and recognized phase labels simply report no membership when cluster analysis was not run; an explicit cage/type/cluster/domain target that never occurs anywhere in the loaded trajectory is rejected.

When cage or guest objects are shown, the generated VMD script uses the following stable cage-type colors; guest defaults follow the cage type that selected them. The visible shades follow the active VMD ColorID palette.

| Cage type | VMD ColorID | Default color |
| --- | ---: | --- |
| `5¹²` | 7 | Green |
| `5¹²6²` | 0 | Blue |
| `5¹²6³` | 1 | Red |
| `5¹²6⁴` | 3 | Orange |
| `5¹²6⁸` | 11 | Purple |
| `4³5⁶6³` | 10 | Cyan |
| Other cage types | 2 | Gray |

Ordinary per-frame GRO files are opt-in through `gro` or individual types. `cluster-gro` is separately opt-in and requires cluster search; no mode includes either category by default. With `--output-type none`, only `config.yaml` remains.

Without `--strict`, standalone serial/process/thread read failures become failed summary rows and analysis continues where the reader remains usable. Failed inputs appear in `summary.xlsx/failures`, `<summary_csv_dir>/failures.csv`, and `<summary_detail_dir>/failures.csv` when their respective output types are enabled, and always in the mandatory `config.yaml` `run.failures` list. With `--strict`, SQQ re-raises the error after updating `config.yaml` to `status: failed`.

GRO structure folders, filenames, and title lines use portable ASCII structure labels since version 0.2.4, for example `5^126^2` and `qc_5r_5^36^2_56566`. Markdown and main-summary scientific labels retain their readable superscript notation. This avoids Windows GBK/legacy-reader failures caused by Unicode superscript or subscript characters in generated GRO paths and titles.

Each `*_info.md` report starts with SQQ version, mode/engine, date/time, source, input format, topology when applicable, trajectory stride, frame/time, requested-to-effective graph mode, effective bond mode, ring sizes, status, and molecule counts. Modes display as `00 (sqq-py)`, `py (sqq-py)`, `99 (sqq-cpp)`, or `sqq-cpp`. LAMMPS reports also record units, timestep, atom style, and type-map source.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary-xlsx` writes the plotting-oriented `summary.xlsx` workbook. `summary-csv` uses the same applicable main-table mapping and writes one UTF-8-SIG file per sheet under `summary_csv_dir` (default `summary/`), preserving table names, columns, row order, and values without Excel formatting or tabs. The first `summary` table is a dashboard: Configuration includes `SQQ version`, requested/effective `Graph mode` such as `auto -> hbond`, normalized `Order parameters`, `Find cluster`, and normalized `Output types`; `Analysis Results (min / mean / max)` reports per-frame min/mean/max values while `Frames total / ok / failed` stays a run-level count. Analysis tables such as connection diagnostics, `ring`, `half_cage`, compact composition-level `quasi_cage`, `cage`, optional `hydrate_cluster`, `order_parameter`, and `ice` keep one input file or trajectory frame per row. The other tables have metadata-specific row units: `summary` is a dashboard, `failures` has one failed input/frame per row, and `detail_index` has one generated detail file per row. Detailed configuration tables are not written; the dashboard retains only its compact Configuration block. Ordinary multi-row and isomer tables are written separately under `summary_detail_dir` only when `summary-detail-csv` is selected: optional `failures.csv`, `cage_occupancy.csv`, `cage_isomer.csv`, and `quasi_cage_isomer.csv`. The separate `cluster-detail` type writes `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. The compact `quasi_cage` table aggregates exact quasi-cage isomers into composition-level columns such as `5r_5²6³`, while the detail `quasi_cage_isomer.csv` keeps nonzero exact isomer rows with `quasi_cage_type`, `isomer`, and `count`. `cage_isomer.csv` defaults to observed nonzero isomer rows plus per-frame totals; use `--cage-isomer-rows all` to restore the full zero-filled matrix. The `order_parameter` table contains only the selected F3, F4, Q_l, MCG, and DHOP columns; `--order-parameter none` omits it. Focus mean/count columns are written only when `order.focus_waters` is non-empty. Output type `order-tsv` writes only selected per-water F3/F4/Q_l values because MCG/DHOP are frame-level descriptors.

Summary construction records rows, columns, cells, bytes, CSV/XLSX write time, formatting time, and final-save time in `config.yaml -> run.summary_write`; the terminal prints its total seconds. The mandatory output-root `config.yaml` records final SQQ version, mode/engine, requested and effective graph modes, requested and resolved workers, normalized output types, input metadata, status/failures, and summary timing. Main CSV, XLSX, detail CSV, and `config.yaml` are written to same-directory temporary files and atomically replaced on success or failure. XLSX sheets above 200,000 cells or 128 columns keep header styling, filter, freeze pane, and fixed column widths but skip costly body-cell formatting; scientific values and table schemas are unchanged.

The `hydrate_cluster` main-summary table reports the mutually exclusive `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`. Optional cluster-detail CSV records add the corresponding cage-id groups and `boundary_composition`; hydrate-domain CSV records expose only external boundary contacts through `external_boundary_contact_count` and `external_boundary_contact_ids`.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

SQQ-Py cage files include cage waters, `CNT` center pseudoatoms, and assigned guests. SQQ-CPP cage files omit the synthetic `CNT` center pseudoatom. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
