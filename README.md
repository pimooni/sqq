# SQQ

**SQQ (Shell Quant Qualifier): Python Joint Toolkit for Water-Shell Topology Analysis.**

Current development version: **0.5.4** (Aug 19, 2026; Qixi Festival)

SQQ provides the complete SQQ-Py water-shell topology workflow and the focused SQQ-CPP cage engine. Select the Python workflow with `-e py` or the C++17 graph/ring/cage/occupancy/F3/F4 workflow with `-e cpp`. Algorithms are documented in `docs/design.md` and release notes in `docs/update.md`.

## Acknowledgements

Names are listed alphabetically by family name.

- Cao, Pinqiang    @  Wuhan University of Science and Technology
- Cheng, Liwei     @  Wuhan Institute of Technology
- Fang, Bin        @  Hainan University
- Hu, Yifei        @  Fuzhou University
- Jia, Jihui       @  China University of Petroleum (Beijing)
- Li, Wuquan       @  Beijing Huairou Laboratory
- Li, Zhenchao     @  Fuzhou University
- Liao, Bo         @  China University of Petroleum (East China)
- Lu, Yingxu       @  Wuhan Institute of Technology
- Mi, Fengyi       @  Southwest University of Science and Technology
- Sun, Yingtao     @  The Hong Kong University of Science and Technology
- Xu, Hongye       @  The University of Tokyo
- Zhang, Zhengcai  @  Laoshan Laboratory
- Zhao, Jingyuan   @  Akatsuki Games Inc.

## Changed in 0.5.4

- Made generated VMD scripts responsive while membership and trajectory data load, added transactional cleanup for repeated sourcing and failed loads, and made renderer readiness the final successful initialization step.
- Kept semantic VMD representations stable and created them lazily, inherited user-edited styles when a compatible layer first appears later, and reduced large Track target initialization to one redraw.
- Stored fixed component metadata only once at render frame 0 while validating that every later frame has the same component topology.
- Condensed Analyze and Track into one terminal page, removed the execution panel for every true single-task run, and made terminal and summary citation text share the same completed-feature evidence.

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

The first invocation in a new environment can take a few extra seconds while MDAnalysis dependencies create a local Matplotlib font cache; subsequent runs are unaffected.

Then use:

```bash
sqq -h
sqq --version
sqq vmd -h
sqq init
sqq analyze -i ./gro -c sqq_config.yaml -o ./result_sqq
```

Root help prints the SQQ version and release date immediately before the usage line. Use `sqq -v` or `sqq --version` for the version line alone.

During source-tree development without installation:

```bash
python -m sqq analyze -i ./gro -c sqq_config.yaml -o ./result_sqq
```

## Quick Start

Single GRO file:

```bash
sqq analyze -i test1.gro -o ./result_sqq
```

Directory of GRO files (the default `input.pattern` in `sqq_config.yaml` is `*.gro`):

```bash
sqq analyze -i ./gro -o ./result_sqq
```

For multiple GRO files, SQQ groups compatible frames automatically. Files with one topology share the requested output root; heterogeneous inputs are separated into `result_A`, `result_B`, and so on in first-occurrence order. The grouping affects aggregation and paths only, not per-frame analysis.

Glob pattern:

```bash
sqq analyze -i "./gro/*.gro" -o ./result_sqq
```

XTC/TRR trajectory with a topology file; add `-dt 100` to analyze an exact 100 ps interval:

```bash
sqq analyze -i traj.xtc --top topol.gro -dt 100 -c sqq_config.yaml -o ./result_sqq
```

A stacked GRO trajectory uses repeated complete GRO blocks in one file and needs no separate topology:

```bash
sqq analyze -i frames.gro -dt 100 -o ./result_sqq
```

LAMMPS dump or DCD with a DATA topology; standard water/methane types are inferred automatically:

```bash
sqq analyze -i traj.lammpstrj -t system.data -o ./result_sqq
```



### Input Units and Boxes

GRO and MDAnalysis trajectory coordinates are interpreted in nm. A GRO input may contain one frame or repeated complete GRO blocks in one stacked trajectory. SQQ streams every block, validates atom counts, ordered atom identity, box records, and finite coordinates, and rejects topology changes between stacked frames. A GRO used as `--top` must still contain exactly one frame. Other trajectory frames also require finite coordinates. XYZ coordinates are multiplied by YAML `input.xyz_scale`; the default `0.1` assumes angstrom input, while `1.0` keeps nm values. SQQ accepts exactly one declared XYZ frame per file and rejects truncated, extra, malformed, or non-finite atom records. XYZ has no periodic box unless converted through another format.

GRO atom counts and the mandatory box line are validated. A three-value positive box is orthorhombic; an all-zero box is treated as non-periodic. Nine-value GRO boxes with nonzero tilt terms and trajectory frames with non-90-degree angles are rejected because triclinic minimum-image calculations are not implemented. GRO molecules are formed from contiguous residue blocks in source order, preventing wrapped or repeated residue IDs from merging distinct molecules. LAMMPS normally uses DATA molecule IDs; automatic inference can rebuild them from unambiguous Bonds components, and dump atom rows may be interleaved.

LAMMPS trajectories require `-t system.data` (equivalent to `--top`). A non-empty `input.lammps.type_map` explicitly maps numeric atom types to `resname`/`atomname` or `ignore` and always takes priority. If the map is absent or empty, SQQ uses DATA masses, type comments, molecule IDs, and Bonds to identify unambiguous water (`1 O + 2 H`), all-atom methane (`1 C + 4 H`), and labeled single-site methane. Other bonded components are retained deterministically as environment/other context instead of being mistaken for water or guests; they do not enter the water graph or cage occupancy. Use `component.role_map`, `additive.resname`, `environment.resname`, or an explicit `type_map` when the automatic role is not the intended one. If molecule IDs do not define valid water/guest molecules but Bonds do, SQQ rebuilds deterministic molecule IDs and reports that decision. Ambiguous reuse of one atom type, insufficient evidence for a requested water/guest role, or topology/trajectory ID mismatch still fails clearly. The resolved mapping and role provenance are recorded in `sqq_config_resolved.yaml`, per-frame info, and main-summary configuration. This normalization is shared by SQQ-Py and SQQ-CPP. Supported inputs are LAMMPS DATA with `full`, `molecular`, `bond`, or `angle` atom style, fully periodic `pp pp pp` orthorhombic dump boxes, and LAMMPS DCD. Tilted boxes, nonperiodic dump boundaries, `units lj`, duplicate atom IDs, and topology/trajectory ID mismatches fail before analysis. `input.delta_time_ps` / `-dt` / `--delta-time` selects a physical interval in ps for XTC, TRR, LAMMPS dump/DCD, and stacked GRO trajectories. With no delta time, every stored frame is analyzed. For LAMMPS dumps, the native reader interval is passed explicitly as `input.lammps.timestep × units-to-ps`; this preserves physical time and prevents MDAnalysis from substituting `1 ps` without hiding unrelated warnings. The requested interval must be at least, and an integer multiple of, the regular native interval; missing or irregular time metadata is rejected instead of rounded.

## Analysis Engines

`-e` / `--engine` selects the analysis engine; the default is `py`:

| Engine | Implementation | Main scope | Default workers | Default output types |
| --- | --- | --- | --- | --- |
| `py` | SQQ-Py | Complete graph, ring, open-patch, cage, cluster, order-parameter, and ice workflow | 1 worker | `info,sqq-render,summary-xlsx` |
| `cpp` | SQQ-CPP | Native graph, internal 4/5/6 rings, cage/isomer/occupancy, and F3/F4 | 1 worker | `info,sqq-render,summary-csv,summary-detail-csv` |

For every successful Analyze sequence, both SQQ-Py and SQQ-CPP write `track/track_state.json` and the six normalized Track tables: `cage_observation.csv`, `cage_track.csv`, `cage_event.csv`, `cage_population.csv`, `guest_residence.csv`, and `lifetime_distribution.csv`. These data do not depend on selecting `sqq-render`; the engine defaults happen to select the visualization package as well.

```bash
sqq analyze -i ./gro -e py -o ./result_py
sqq analyze -i ./gro -e cpp -o ./result_cpp
```

For SQQ-Py, `--find-cluster` overrides `hydrate_cluster.enabled` in YAML. Search results enter selected info/main-summary outputs; split cluster structures still require the YAML output type `cluster-gro`. SQQ-CPP rejects cluster search.

Both documented engines default to one worker. `-w` / `--worker` overrides the preset: integer text is a worker count, while `0.5`, `1.0`, `50%`, and `100%` are physical-core fractions. Process parallelism supports independent GRO/XYZ files and indexed XTC/TRR/LAMMPS trajectories. At most `3 * workers` tasks are submitted at once.

The default `chordless`/`bounded` path preserves the established scientific definitions while accelerating neighbor generation, incremental chord pruning, L1 forward checking, cached layer growth, integer-mask subset ownership, and cage target/edge state pruning. Cage DFS also applies exact remaining-edge incidence and parity conditions before expansion. MDAnalysis supplies orthorhombic cutoff candidates when available, but SQQ still rechecks every distance and hydrogen-bond angle with its established float64 logic. F3 and graph-mode Q_l share one graph-vector cache; all Q_l degrees share candidate lists and spherical-angle work. Optional `ring.definition: shortest_path` applies the Franzblau shortest-path criterion and reuses bounded-BFS distance maps. Optional `quasi_cage.search_policy: exact` preserves distinct frontiers and enumerates connected L2/L3 subsets; these opt-in modes can change or add results. Quasi-cage candidate and layer-state truncation is reported through frame warnings; the current cage GROW search itself does not truncate candidates or return partial results.

Every cage now passes the same mandatory topology validation in SQQ-Py and SQQ-CPP: each edge belongs to exactly two faces, `V - E + F = 2`, the face shell is connected, every vertex link is one cycle, and every shell vertex is trivalent. Optional scientific cage validation adds PBC-aware face-planarity and edge-variation limits, nonzero projected area, positive-volume validation, and volume-centroid cage centers. It remains disabled by default, but disabling it no longer bypasses topology validation. SQQ uses an orthorhombic box representation and rejects non-orthogonal/triclinic input explicitly.

The current development version keeps the complete run on one compact terminal page. `Started`, resolved `Input`, `SQQ`, `Graph`, `Ring sizes`, and `Additional search` combine information that was previously repeated across several rows. Serial and parallel progress retain the same three-stage model: file preparation (`reading`, `settings`, `selecting`), core topology search (`graph`, `ring`, optional `half/quasi`, `cage`, and optional `cluster`), and post-processing (`filtering`, `order`, `ice`, `output`). An execution plan containing exactly one task, including a directory that resolves to one file, has no execution row or one-item progress bar. Multi-item runs show one execution-policy row, aggregate stages, and at most five active files; shorter terminals automatically reduce live-progress detail. The TTY panel updates in place, while redirected output prints one bounded starting progress record without ANSI or carriage-return rewrites and leaves final counts to `Analysis Results`. The normal successful completed page targets no more than 35 lines; no separate 80-column by 24-row layout is defined. Nonfatal warnings are deduplicated and reduced to one final preview with full details retained in result files.

### Native SQQ-CPP Backend

Engine `cpp` selects the focused native workflow. Python owns input normalization, molecule selection, scheduling, full-frame VMD output, Markdown, summary CSV, and optional XLSX; C++17 performs graph construction, internal chordless 4/5/6 rings, cage topology/isomers, occupancy, and F3/F4 while releasing the GIL. Native periodic and non-periodic candidate generation uses an exact cell list followed by the unchanged double-precision distance/angle checks; frame-local ring and vector caches reduce repeated work without introducing an approximate neighbor rule.


It accepts orthorhombic GROMACS/LAMMPS inputs, compatible graph/pair settings, `-s` within 4/5/6, cage report/validation settings, `f3`/`f4`, process or serial scheduling, and `info`, `gro`, `cage-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, or `summary-detail-csv`. `sqq-render` owns the complete visualization package. `gro` enables the supported classified cage GRO output, but `cpp` does not select it by default.

Unsupported requests fail before analysis: public ring output, size 7, shortest-path rings, half/quasi cages, cluster, ice, Q_l/MCG/DHOP, membership/order TSV, legacy per-frame VMD, Python fast closure, thread scheduling, and triclinic boxes. A failed native extension never falls back to Python.

The `cpp` default layout is:

```text
result/
  sqq_render/
    sqq_cage.gro
    sqq_cage.xtc
    sqq_cage.membership.tsv
    sqq_cage.vmd.tcl
  track/
    track_state.json
    cage_observation.csv
    cage_track.csv
    cage_event.csv
    cage_population.csv
    guest_residence.csv
    lifetime_distribution.csv
  summary/
    summary.csv
    cage.csv
    cage_occupancy.csv
    cage_isomer.csv
    order_parameter.csv
    detail_index.csv
  sqq_config_resolved.yaml
  sqq_output_manifest.json
  frame_name/
    frame_name_info.md
```

The native engine does not select ordinary/classified GRO by default. Set YAML `output.type` to include `gro` or `cage-gro` when that output is required.

Release CI is configured to build and test precompiled wheels for CPython 3.10-3.14 on Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, plus a source distribution. A wheel already contains the platform-native extension; end users installing such a wheel do not compile C++. A source install instead invokes the CMake/scikit-build-core build and therefore needs CMake 3.20 or newer and a local C++17 toolchain.

## Package Architecture

SQQ 0.5.4 separates public configuration and data models, parallel `core/sqq_py` and `core/sqq_cpp` scientific backends, command workflows (`init`, `analyze`, `track`, and `vmd`), runtime scheduling, I/O/reporting/rendering, and terminal UI. Retired monolithic modules are not compatibility entry points; new code should use the public API or the responsibility-specific package paths. The shared VMD Tcl template is kept as a readable Python string in `sqq/io/render/tcl_template.py`; generated render packages still contain the ordinary `*.vmd.tcl` script required by VMD.

### Python API

The supported programmatic entry points are exported directly by `sqq`:

```python
from sqq import analyze_frame, load_config, read_frames

config = load_config({"graph": {"mode": "oo"}}, engine="py")
frame = next(read_frames("frame.gro", config=config))
result = analyze_frame(frame, config)
```

`load_config` returns a resolved configuration that carries its resolution record, `read_frames` yields the public `Frame` model, and `analyze_frame` returns `FrameResult`. Invalid configuration, input, or analysis requests raise the exported typed SQQ exceptions.

## Common Commands

Write the commented default configuration to `sqq_config.yaml`. The template includes `#` section and choices comments, defaults to ring sizes 4/5/6, and refuses to overwrite an existing destination:

```bash
sqq init
```

Use `-o` only when a different configuration filename is wanted; that destination must also not already exist:

```bash
sqq init -o methane.yaml
```

SQQ does not auto-load a similarly named file from the current directory. Without `-c`, built-in defaults are used. With `-c`, the named user file is read but never rewritten.

Select ring search sizes, open-patch searches, hydrate clusters, and order parameters with the retained public overrides:

```bash
sqq analyze -i md.gro -e py -s 4,5,6 --find-half on --find-quasi on
sqq analyze -i md.gro --find-cluster on --order-parameter f3,f4,q6
```

Analyze an explicit pair map:

```bash
sqq analyze -i md.gro -b pairs --pair water_pairs.txt
```

`--pair` overrides YAML `graph.pair_file`. A CLI-relative pair path is resolved from the working directory; a YAML-relative path is resolved from the directory containing the user configuration. `graph.mode: pairs` without either source fails before frame analysis.

Parallelize independent files or indexed trajectory frames with the default process backend:

```bash
sqq analyze -i ./gro -w 4 -o ./result_sqq
sqq analyze -i traj.xtc -t topol.gro -w 50% -o ./result_sqq
```

Integer worker text is an explicit count. Decimal text and percentages are physical-core fractions, so `-w 1` is one worker while `-w 1.0` and `-w 100%` request all detected physical cores before the reserve-one-core and task-count clamps.

Track every detected cage directly from a trajectory, or reuse an Analyze result without repeating frame science:

```bash
sqq track -i traj.xtc -t topol.gro -dt 100 --target all -o ./result_track
sqq track --source ./result_sqq --target 512,51264,sI,t133 -o ./result_track
```

Targets may be `all`, one or more cage types, hydrate phases, or persistent cage IDs (`t1`, `t2`, ...). Comma-separated targets are written independently, so `--target 512,sI` creates both `type_512/` and `phase_sI/`. A phase target automatically resolves `find_cluster` to `on` and requires SQQ-Py; `-e cpp` is rejected for phase targets. Source mode cannot add phase labels retroactively, so imported state must already contain them. With neither `-i` nor `--source`, Track searches the current directory for exactly one Analyze tracking state.

Source mode inherits the Analyze engine from `sqq_config_resolved.yaml`, so the terminal and new resolved configuration retain the original `sqq-py` or `sqq-cpp` identity; a C++ source also normalizes half/quasi settings to SQQ-CPP capabilities.

Raw-input Track currently analyzes selected frames serially: any `-w` / `--worker` request is normalized to one worker and the serial backend. Track output is fixed: `--output-type` does not change the set consisting of the resolved configuration, run-level state and six CSV tables, independently filtered target directories, and target render packages. For a persistent-ID target with pre-birth frames, raw mode first builds the cross-frame `tID`, then reanalyzes only the required prefix through the birth frame to add precursor tables and precursor membership to that target's VMD package. A target already present in the first selected frame has no precursor interval. Imported Analyze state contains cage snapshots only and reports precursor history as unavailable.

## Configuration

The generated file uses YAML `#` comments and canonical singular keys. The main settings are:

```yaml
schema_version: "0.5.4"
engine: py  # choices: py, cpp

run:
  strict: false  # choices: true, false

input:
  pattern: "*.gro"
  recursive: false  # choices: true, false
  delta_time_ps: null
  xyz_scale: 0.1
  lammps:
    unit: real  # choices: real, metal, nano
    timestep: 1.0
    atom_style: full  # choices: full, molecular, bond, angle
    coordinate_convention: auto  # choices: auto, x, xs, xu, xsu
    type_map: {}

component:
  auto_classify: true
  unknown_role: other
  unknown_action: warn
  role_map: {}

water:
  resname: [SOL, TIP, WAT, HOH]
  oxygen_name: [OW, O, OH2]
  hydrogen_name: [HW1, HW2, H1, H2, HW, HT1, HT2]

guest:
  resname: [CH4, CO2, MET, ETH]
  center_atom:
    CH4: [C]
    CO2: [C]
    MET: [C]
  center_mode: center_atom

additive:
  resname: []

environment:
  resname: []

graph:
  mode: auto  # choices: auto, hbond, oo, pairs
  oo_cutoff_nm: 0.35
  hbond_distance_nm: 0.35
  hbond_angle_deg: 30.0
  pair_file: null
  pair_id: resid  # choices: resid, oxygen_index, atomid

ring:
  size: [4, 5, 6]
  report_size: auto
  definition: chordless  # choices: chordless, shortest_path

half_cage:
  enabled: auto  # choices: auto, true, false

quasi_cage:
  enabled: auto  # choices: auto, true, false
  base_size: auto
  side_size: auto
  max_layer: 1
  search_policy: bounded  # choices: bounded, exact

cage:
  enabled: true
  report_type: auto
  max_face: 20
  search_mode: grow
  seed_mode: ring
  max_state_per_seed: 0  # 0 = unlimited
  max_total_state: 0  # 0 = unlimited
  max_boundary_candidate: 8  # compatibility setting; never truncates exact search
  scientific_validation: false
  max_face_planarity_rms_nm: 0.06
  max_face_edge_cv: 0.35
  min_cage_volume_nm3: 1.0e-6
  occupancy_mode: polyhedron

hydrate_cluster:
  enabled: false
  min_cage: 2

hydrate_order:
  mcg_guest_resname: [CH4, MET]
  mcg_guest_cutoff_nm: 0.90
  mcg_water_cutoff_nm: 0.60
  mcg_cone_half_angle_deg: 45.0
  mcg_min_water: 5
  dhop_neighbor_cutoff_nm: 0.35
  dhop_planar_count: [11, 12]
  dhop_min_qualified_neighbor: 3

order_parameter:
  enabled: [f3, f4]  # choices: f3, f4, qN, mcg1, mcg3, dhop35, dhop30, all, none
  q_neighbor_mode: graph  # choices: graph, cutoff, nearest, lammps
  q_cutoff_nm: 0.35
  q_n_neighbor: null


parallel:
  backend: process  # choices: process, thread, serial
  worker: auto
  math_thread: 1

output:
  type: [info, sqq-render, summary-xlsx]
  summary_csv_dir: summary
  cage_isomer_row: nonzero  # choices: nonzero, all
  write_empty_file: false
  structure_layout: grouped  # choices: grouped, flat
  center_resname: CNT

render:
  atom_scope: full  # choices: full, compact

track:
  target: all
  source: null
  min_jaccard: 0.50
  min_shared_fraction: 0.60
  min_shared_water: 3
  max_center_distance_nm: null
  gap_frame: 0
  guest_tiebreak: true
```

Unknown keys and duplicate YAML keys are errors. Canonical public collections are singular, units appear in names such as `_ps`, `_nm`, and `_deg`, engine-related three-state switches use `auto/true/false`, and ordinary booleans use `true/false`. Legacy top-level `mode`, `graph.bond_mode`, and `order.parameter` migrate with warnings to `engine`, `graph.mode`, and `order_parameter.enabled`. Former 0.3.x plural keys also remain readable for migration; generated and resolved files use only the canonical form. Inactive legacy `output.gro_atom_mode` and `output.context_role(s)` values are ignored with a recorded migration adjustment instead of being silently accepted. `output.center_resname` is the 1-5 character ASCII residue name used for synthetic ring, half/quasi-cage, and SQQ-Py cage-center atoms in selected category GRO files.

Configuration priority is:

```text
built-in defaults < engine preset < sqq_config.yaml < retained command-line overrides
```

Every run writes the final effective state to `sqq_config_resolved.yaml` in its result root, including the requested engine, effective `sqq-py`/`sqq-cpp` backend, requested and effective graph modes, requested and resolved workers, output selection, input/LAMMPS provenance, automatic adjustments, run status, failures, and summary-write timing. This file is separate from the user-owned `sqq_config.yaml`.

Track matching is deterministic. Each water is identified across frames by its stable one-based topology atom position, not the width-limited serial stored in a GRO atom record; cages receive persistent `t1`, `t2`, ... IDs. Guest membership uses a topology-stable `g########` identifier rather than potentially repeated residue labels. Member-water overlap/Jaccard drives cage assignment; cage topology and orthorhombic-PBC center displacement constrain candidates; guest continuity is only a tie-break. `gap_frame: 0` means that one missing selected-frame observation ends a track. Positive values permit only explicitly recorded gaps and produce `gap` observations/events rather than silently joining discontinuous cages.

## Parallel Execution

Public YAML uses singular collection keys such as `water.resname`, `ring.size`, `order_parameter.enabled`, `output.type`, and `parallel.worker`. Configurations from 0.3.x that use former plural spellings are migrated on read; `sqq init` and the resolved runtime file use the singular schema.

`half_cage.enabled: auto` and `quasi_cage.enabled: auto` resolve to `on` for SQQ-Py and `off` for SQQ-CPP. If an older YAML explicitly enables either Python-only search under C++, SQQ disables the unsupported work, deactivates quasi layer controls, removes incompatible half/quasi outputs, records the adjustment in `sqq_config_resolved.yaml`, and continues. Missing required inputs or settings that prevent the native cage calculation remain hard errors.

`parallel.backend: process` is the default for two or more independent GRO/XYZ inputs. SQQ uses the `spawn` start method on every supported platform. Each worker receives run configuration once, reads and writes its own frame, and sends only small stage events plus one summary row to the main process. This avoids the Python GIL limitation of the compatibility thread backend.

Before dispatching two or more GRO files, SQQ reads only their topology records and assigns topology groups in first-occurrence order. The fingerprint contains the atom count and ordered contiguous residue blocks, represented by each block's residue name and ordered atom-name sequence. Titles and time labels, coordinates, velocities, boxes, and numeric atom/residue IDs do not affect grouping. A supplied GRO `-t` / `--top` is checked against every input fingerprint; any mismatch fails before analysis and identifies the exact source file.

All accepted groups use one shared worker pool and one global progress index. Each task also carries a group-local frame index and output root, so group summaries and annotated bundles remain correctly ordered without running groups serially. Requested `graph.mode: auto` remains recorded as `auto`, but its effective mode is resolved once from a representative frame in each topology group and reused by both SQQ-Py and SQQ-CPP for every frame in that group: every selected water with at least two hydrogens selects `hbond`, every selected water with no hydrogen selects `oo`, and mixed or incomplete hydrogen topology is rejected.

With `parallel.worker: auto`, the documented `py` and `cpp` engines resolve to one worker. Physical-core detection for explicit fractional requests prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker` / `-w` accepts either a fraction (`50%`, `0.5`, or `1.0` for 100%) or an explicit positive integer worker count (`1` means one worker). Windows `ProcessPoolExecutor` runs are capped at 61 workers; Linux workstations can use larger explicit values such as `-w 100`, subject to the reserve-one-core rule, task count, memory, and storage throughput.

One XTC/TRR or supported LAMMPS trajectory with `--top` is frame-parallel when the process backend resolves to more than one worker. Every worker opens a private MDAnalysis Universe once and seeks small contiguous batches of selected raw frame indexes; batch size is automatically bounded from 1 to 8, and complete coordinate arrays are not serialized between processes. Parent and worker trajectory readers are explicitly closed. Multiple trajectory files and the compatibility thread backend use the serial trajectory reader.

Process submission uses a bounded rolling queue of at most `3 * workers` tasks. This is a queue-depth limit, not a CPU limit: with 100 effective workers SQQ may keep up to 300 tasks submitted while still running as many as 100 workers concurrently. Results are restored to original file/frame order before main-summary writing.

The parent preserves original input order globally and group-local order in every selected group summary and annotated bundle. Output-name collisions are resolved deterministically within each topology group. Process runs set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, and `BLIS_NUM_THREADS` to `parallel.math_thread` while workers are spawned, then restore the parent environment.

The scheduling and search-cache refinements themselves do not change existing scientific definitions or values. Before the new hydrate descriptors were enabled, they reduced the local `1200ns.gro` serial run from about 26.6 s to 18.2 s. A 0.2.3 benchmark that also selected MCG-1 and DHOP35 completed in about 21.6 s on the same host; every overlapping pre-existing analysis column matched the earlier workbook. Performance depends on data, configuration, CPU, memory, and storage.

## Search and Report Scope

`-s` / `--size` defines the ring-face sizes used during detection and, by default, reporting. YAML `ring.report_size` and `cage.report_type` can narrow user-facing output without changing the shared search universe:

```yaml
ring:
  size: [4, 5, 6]
  report_size: [5, 6]

cage:
  report_type: [I, II]
```

```bash
sqq analyze -i md.gro -c sqq_config.yaml -s 4,5,6
```

`cage.report_type` accepts `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, and `I2II`; group names may be listed together. `auto` follows the selected search sizes, while `all` reports every detected cage composition in scope. Do not combine `auto` or `all` with named groups.

Repeated cage types contributed by several groups are reported once. All detected cages still participate in half-cage, quasi-cage, free-ring filtering, and hydrate-cluster topology. A report filter changes user-facing cage counts and files, not topology ownership. Cage detection supports 4/5/6 faces; ring and quasi-cage detection also support size 7 in SQQ-Py.

## Exact Sparse Cage Search and Scientific Validation

SQQ-Py and SQQ-CPP use the same exact frame-local cage-search contract. The shared topology stores compact ring-to-edge and edge-to-ring incidence. A search branch stores only its selected ring IDs, local edge-use counts, and open-edge frontier; it does not allocate a frame-wide bit mask per ring, edge, or state. Ring centers, normals, and adjacency are built only when another enabled analysis actually needs them.

GROW starts from each canonical seed ring, chooses a constrained open edge, and visits every topologically eligible neighboring ring in deterministic order. Candidate count is never used to prune the search. Per-seed duplicate-state detection is exact and released when that seed finishes. `cage.max_boundary_candidate` is retained only so older configuration files remain readable; it does not truncate candidates or change the cage set.

`cage.max_state_per_seed: 0` and `cage.max_total_state: 0` are the defaults and mean unlimited exact search. Positive values are diagnostic safety guards. Reaching either guard aborts the frame with a clear error and publishes no partial result. The former half-cage fast-closure recovery path is removed; older configuration keys are ignored during migration and recorded as adjustments in `sqq_config_resolved.yaml`.

Topology validation is always enabled. Every candidate must use each edge exactly twice, satisfy `V - E + F = 2`, form one edge-connected face shell, have one cyclic face link around every vertex, and have only trivalent shell vertices. These checks reject disconnected, pinched, branched, and non-manifold false cages before type/isomer assignment in both engines. Accepted cages receive a deterministic final order from cage type, water membership, and face topology, so SQQ-Py and SQQ-CPP assign the same frame-local IDs to the same cage set.

`cage.scientific_validation: false` is the default. When set to `true` in YAML, a topologically valid cage must additionally satisfy the configured PBC-aware face-planarity RMS and edge-length coefficient-of-variation limits, nonzero projected face area, and positive minimum triangulated volume. Each face is unwrapped independently through its cyclic edges before its quality checks; the complete shell is unwrapped separately for volume and centroid. Accepted cages then use the volume centroid instead of the mean cage-water position. Enabling it can therefore remove geometrically distorted cages and can change guest occupancy or geometry-resolved hydrate-cluster edges. Raw ring and half/quasi searches remain unchanged; ownership-filtered free-ring and free-patch outputs can increase when a rejected cage no longer consumes them.

Guest occupancy uses the configured center atom when available. The defaults select `CH4`, `CO2`, `MET`, and `ETH` as guests and map `CH4`, `CO2`, and `MET` to atom name `C`, so these residues use their carbon atom under the default `guest.center_mode: center_atom`. Otherwise, guest atoms are PBC-unwrapped around one molecular anchor before calculating the centroid; the same helper is used by MCG.

## Hydrate Cluster

`--find-cluster on` analyzes every detected cage in the selected search scope. Cages become graph nodes and are connected through complete shared ring faces. When several detected cages reference the same face, ring-plane geometry keeps at most one cage on each physical side. YAML `cage.report_type` filters user-facing cage tables and files only; it does not remove cages from cluster connectivity or phase evidence.

SQQ classifies hydrate type, domains, and boundaries on a cage-connection graph using labelled shared-face fingerprints, strict local evidence, distributed spatial cores, mutually compatible expansion, and exclusive per-frame domains.

YAML `hydrate_cluster.min_cage` sets the minimum connected-component size; the default is `2`. Smaller components are counted as isolated cages.

Within each cluster, SQQ builds labelled first-shell fingerprints from neighboring cage types and shared-face sizes. Exact sI/sII/sH fingerprints remain high-confidence seeds. In addition, partial but phase-pure fingerprints can form a distributed spatial core: candidates require at least 50% template coverage, 50% phase purity, and a harmonic support score of 0.55; the compatible cage graph must retain a degree-2 core of at least three cages, mean support of 0.60, and the phase-defining hexagonal large-cage connection. The core is anchored by phase-specific cages (`5^12 6^2` for sI, `5^12 6^4` for sII, and `4^3 5^6 6^3`/`5^12 6^8` for sH), while every validated member of that core enters the initial phase seed. All three phases then expand through mutually compatible face-labelled edges when a candidate has at least two accepted phase contacts. Cages claimed exclusively by one phase form deterministic per-frame domains.

After the exclusive sI/sII/sH domains are finalized, SQQ partitions the remaining cluster cages. A cage enters the generic boundary only when it is outside every phase domain and directly shares a complete cage face with at least one domain cage. Boundary search stops at this first external non-phase layer. Domain cages are never relabelled as boundary, and a direct shared-face contact between different phase domains leaves both endpoint cages in their original phases.

The resulting `classified_cage_ids`, `boundary_cage_ids`, `ambiguous_cage_ids`, and `unclassified_cage_ids` are mutually exclusive and together cover every cage in a reported cluster. Competing phase claims without boundary membership remain ambiguous; all other residual cages are unclassified. There are no `sI-boundary`, `sII-boundary`, `sH-boundary`, transition, or boundary-context categories. Neighboring cages can still share face-water coordinates in structure views, so cage ownership should be verified from cage IDs or detected cage/ring edges rather than coordinate-set overlap.

The default `py` engine leaves cluster search off unless YAML `hydrate_cluster.enabled` or explicit `--find-cluster on` enables it. Engine `cpp` does not support cluster search. Explicit `--find-cluster on|off` has highest priority. Cluster search does not alter ring, patch, cage, occupancy, order-parameter, or ice results. Classification is per-frame and independent of the cage reporting filter. Spatial consensus uses only the current frame: it performs no temporal smoothing, so analyzing a frame alone or inside a compatible batch gives the same phase assignment. Temporal grain tracking and crystallographic orientation matching are not implemented.

Cluster search populates every selected `info` and main-summary output. Split category structures are written only when YAML `output.type` includes `cluster-gro`; no documented engine preset includes it by default. The selected main summary output gains its per-frame `hydrate_cluster` table, while native category structures are written under grouped layout as `<frame>/hydrate_cluster/<frame>_cluster_sI.gro`, `<frame>_cluster_sII.gro`, `<frame>_cluster_sH.gro`, and `<frame>_cluster_boundary.gro`. Flat layout places the same filenames directly in the frame directory. All same-category domains and clusters are aggregated into one file per frame. An absent category is omitted unless `output.write_empty_file: true`.

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

The compact table does not include exact IDs, seeds, confidence values, water/guest membership, or domain adjacency. Add `cluster-detail` to YAML `output.type` for `summary/hydrate_domain.csv` and one-row-per-cluster `summary/hydrate_cluster_detail.csv`. Explicit `cluster-detail` or `cluster-gro` selection requires cluster search. Turning search off writes neither `cluster-detail` nor `cluster-gro` and removes stale generated cluster GRO files. Public motif output is not generated.

## Hydrate Nucleation Order Parameters

MCG-1 and DHOP35 were introduced as defaults in 0.2.5. Since 0.2.7, every MCG/DHOP variant is selected explicitly through `--order-parameter`; the package default is only `f3,f4`. These descriptors are independent of the optional cage-topology `hydrate_cluster` classifier: MCG works on selected methane-like guest centers and surrounding waters, while DHOP works on a dedicated O-O neighbor graph. They do not change graph, ring, patch, cage, occupancy, F3/F4/Q_l, hydrate-cluster, or ice results.

MCG follows the mutually coordinated guest definition. Guest pairs within `0.90 nm` are connected when at least five waters lie within `0.60 nm` of both guests and inside both 45-degree opposing cones. The threshold is **at least five**, not exactly five. MCG-1 keeps guest nodes with at least one qualifying MCG edge; optional MCG-3 applies a one-pass degree-at-least-three filter to the same qualifying graph. Connected components are measured only through qualifying MCG edges. The default guest residue names are `CH4` and `MET`; change `hydrate_order.mcg_guest_resname` for another methane naming convention. If no configured guest type is present, MCG is reported as `N/A`, not zero.

DHOP builds its own orthorhombic-PBC oxygen graph with `hydrate_order.dhop_neighbor_cutoff_nm: 0.35`. This 0.35 nm default follows the all-atom TIP4P/Ice implementation used by Li et al.; use `0.325` in YAML when reproducing the original mW-water definition. For each central O-O bond, SQQ counts neighboring plane-normal pairs within 35 degrees (or 30 degrees for DHOP30), selects waters with counts 11 or 12, requires at least three similarly qualified neighbors, includes their first oxygen shell, and reports the largest connected water cluster. `DHOP35` and `DHOP30` name the angular thresholds, not the O-O cutoff. No transition-state value such as DHOP35=57 is hard-coded; such values are system- and condition-dependent.

Select any combination with names such as `--order-parameter mcg1,mcg3,dhop35,dhop30`. Selection is separate from the numerical `hydrate_order` cutoff settings. All cutoff searches use deterministic cell lists and exact float64 minimum-image rechecks; there are no fixed neighbor-array limits.

References: Barnes et al., MCG ([DOI 10.1063/1.4871898](https://doi.org/10.1063/1.4871898)); Knott et al., MCG nucleation coordinate ([DOI 10.1021/jp507959q](https://doi.org/10.1021/jp507959q)); DeFever and Sarupria, DHOP ([DOI 10.1063/1.4996132](https://doi.org/10.1063/1.4996132)); Li et al., all-atom DHOP nucleation pathway ([DOI 10.1073/pnas.2011755117](https://doi.org/10.1073/pnas.2011755117)).

## Public CLI

Analyze exposes this intentionally compact interface:

| Option | Values / role |
| --- | --- |
| `-i, --input INPUT` | Input file, directory, or glob |
| `-t, --top FILE` | GRO topology for XTC/TRR or LAMMPS DATA for dump/DCD |
| `-c, --config FILE` | User YAML configuration; normally `sqq_config.yaml` |
| `-o, --output DIR` | Result directory |
| `-e, --engine ENGINE` | `py` or `cpp`; default `py` |
| `-w, --worker N` | `auto`, a fraction such as `50%`/`0.5`/`1.0`, or a positive integer count |
| `-dt, --delta-time PS` | Exact physical sampling interval in ps |
| `-b, --bond-mode MODE` | `auto`, `hbond`, `oo`, or `pairs` |
| `-s, --size SIZES` | Comma-separated ring search sizes |
| `--find-half on|off` | Override standard half-cage search |
| `--find-quasi on|off` | Override layered quasi-cage search |
| `--find-cluster on|off` | Override SQQ-Py hydrate-cluster search |
| `--order-parameter NAMES` | `f3`, `f4`, `qN`, `mcg1`, `mcg3`, `dhop35`, `dhop30`, `all`, or `none`; comma-separated |
| `--pair FILE` | Explicit water-network edge file; enables pairs mode unless `-b pairs` is already present |
| `--output-type TYPES` | Replace the Analyze output list; Track uses its fixed Track output set |
| `-h, --help` | Show command help |


### Migration errors

The former mode spelling is not silently converted:

```text
$ sqq analyze -i md.gro --mode py
Error: --mode has been replaced by --engine.
Use: --engine py

$ sqq analyze -i md.gro -m py
Error: -m has been replaced by -e.
Use: -e py
```

Likewise, the former plural pair option stops with an actionable error:

```text
$ sqq analyze -i md.gro --pairs water_pairs.txt
Error: --pairs has been replaced by --pair.
Use: --pair water_pairs.txt
```

These errors exit with status `2`. The deprecated spellings are not hidden aliases. Expected configuration, input, I/O, output-lock, and analysis failures likewise end with one concise `Error: ...` line instead of a Python traceback; set `SQQ_DEBUG=1` when a development traceback is required. Analyze completes its read-only input/configuration preflight before it creates or cleans the requested result directory.

Former advanced CLI settings now belong only in YAML:

| Former CLI setting | Canonical YAML key |
| --- | --- |
| `--pattern`, `--recursive`, `--strict`, `--xyz-scale` | `input.pattern`, `input.recursive`, `run.strict`, `input.xyz_scale` |
| `--lammps-units`, `--lammps-timestep`, `--lammps-atom-style` | `input.lammps.unit`, `input.lammps.timestep`, `input.lammps.atom_style` |
| `--ring-size`, `--ring-definition` | `ring.report_size`, `ring.definition` |
| `--quasi-size`, `--quasi-base-size`, `--quasi-side-size`, `--quasi-max-layer`, `--quasi-search-policy` | both size lists; `quasi_cage.base_size`; `quasi_cage.side_size`; `quasi_cage.max_layer`; `quasi_cage.search_policy` |
| `--cage-size`, `--max-cage-face`, `--cage-scientific-validation` | `cage.report_type`, `cage.max_face`, `cage.scientific_validation` |
| `--cluster-min-cage` | `hydrate_cluster.min_cage` |
| `--q-neighbor-mode`, `--q-cutoff`, `--q-n-neighbor` | `order_parameter.q_neighbor_mode`, `order_parameter.q_cutoff_nm`, `order_parameter.q_n_neighbor` |
| `--pair-id`, `--parallel-backend` | `graph.pair_id`, `parallel.backend` |
| `--output-layout`, `--cage-isomer-rows` | `output.structure_layout`, `output.cage_isomer_row` |

The former `--cage-fast-closure` option is removed. Legacy YAML fast-closure keys remain readable for migration, are ignored, and are recorded as configuration adjustments; they are not active settings in the resolved configuration.

The legacy compatibility names `--workers`, `--no-q`, `-q`, `--q-degree`, `--mcg3`, `--dhop30`, and `--topology` are removed. Use `-w` / `--worker`, `--order-parameter`, and `-t` / `--top` as applicable.

### Bond Mode and Pair Files

Use `-b` / `--bond-mode` to override YAML `graph.mode`:

```bash
sqq analyze -i md.gro -b auto
sqq analyze -i md.gro --bond-mode hbond
sqq analyze -i md.gro -b oo
sqq analyze -i md.gro -b pairs --pair pairs.txt
```

`--pair PAIRS.txt` alone is shorthand for pairs mode. Combining it with explicit `-b auto`, `-b hbond`, or `-b oo` is rejected. Pairs mode requires either `--pair` or YAML `graph.pair_file`; the identifier convention is YAML `graph.pair_id`.

### Output Selection

Output selection can be set with `--output-type` or YAML `output.type`. Engine `py` defaults to:

```yaml
output:
  type: [info, sqq-render, summary-xlsx]
```

SQQ-Py accepts `info`, `membership-tsv`, `order-tsv`, `f3-gro`, `f4-gro`, `sqq-render`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `summary-xlsx`, `summary-csv`, `summary-detail-csv`, and `cluster-detail`, plus `default`, `all`, and `none`. `gro` expands to ordinary ring/half/quasi/cage/ice GRO categories. `sqq-render` writes the complete four-file visualization package; its files cannot be selected independently.

SQQ-CPP accepts `info`, `gro`, `cage-gro`, `f3-gro`, `f4-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, and `summary-detail-csv`, plus `default`, `all`, and `none`. `default` may be combined with extra types, for example `--output-type default,summary-detail-csv`; duplicate types are removed. `all` and `none` remain exclusive. Removed types `sqq-cage-gro` and `vmd` are rejected with a message to use `sqq-render`. The `cpp` preset does not select `gro` or `cage-gro` by default. Cluster-specific output requires SQQ-Py with resolved cluster search on. `sqq_config_resolved.yaml` is always written regardless of `output.type`.

### Per-water F3/F4 GRO

`--output-type f3-gro` and `--output-type f4-gro` are explicit, non-default outputs supported by both engines. The matching parameter must also be selected, for example `--order-parameter f3,f4 --output-type f3-gro,f4-gro`. Each file contains only waters with a defined per-water value, but retains every atom of each selected water in source order. Only the oxygen record is annotated (`; SQQ F3=<value>` or `; SQQ F4=<value>`); hydrogen and virtual-site records are unannotated. GRO coordinate columns and available velocity columns remain fixed-width, and the annotation begins after the velocity field. Grouped output is written below `<frame>/order/`; in separated multi-frame/multi-file output it is below `gro/<frame>/order/`. Empty files follow `output.write_empty_file`. These files expose existing per-water values and do not change their calculation.

## Output Structure

The layout decision is based on the complete job rather than a particular reader or scheduler. A single ordinary one-frame GRO or XYZ keeps the compact frame-root layout. Every trajectory-like input, including one or more XTC/TRR/LAMMPS/DCD paths or a stacked GRO, and every multi-file GRO/XYZ job uses the separated `info/` plus `gro/<frame>/` layout in serial, thread, and process execution.

For two or more independent GRO files, topology grouping controls only the aggregation root. If every GRO has one compatible topology, all selected outputs are written directly under the requested result directory:

```text
result/
  sqq_config_resolved.yaml
  sqq_output_manifest.json      # paths owned by the completed SQQ run
  summary.xlsx                 # when summary-xlsx is selected
  summary/                     # summary-csv/detail/cluster-detail CSVs
    summary.csv
    cage.csv
    ...
  info/
    frame_001_info.md
    frame_002_info.md
  gro/                         # when any per-frame GRO output is selected
    frame_001/
      order/                   # explicit f3-gro/f4-gro
        frame_001_f3.gro
        frame_001_f4.gro
    frame_002/
  sqq_render/
    sqq_cage.gro              # stable topology and first selected frame
    sqq_cage.xtc              # every selected render frame
    sqq_cage.membership.tsv   # typed frame/center/guest/membership metadata
    sqq_cage.vmd.tcl          # VMD loader and commands
  track/
    track_state.json          # persistent IDs and normalized observations/events
    cage_observation.csv
    cage_track.csv
    cage_event.csv
    cage_population.csv
    guest_residence.csv
    lifetime_distribution.csv

```

When 2-26 distinct topologies are found, groups are assigned letters by first occurrence and each group gets a complete independent result root. No summary, GRO, or VMD bundle combines incompatible systems:

```text
result/
  sqq_config_resolved.yaml              # batch manifest and source-to-group mapping
  sqq_output_manifest.json              # ownership manifest for every group output
  result_A/
    sqq_config_resolved.yaml
    summary.xlsx               # and/or summary/
    info/
    gro/                       # when selected
    sqq_render/                # when selected
      sqq_cage.gro
      sqq_cage.xtc
      sqq_cage.membership.tsv
      sqq_cage.vmd.tcl
    track/
  result_B/
    sqq_config_resolved.yaml
    summary.xlsx               # and/or summary/
    info/
    gro/                       # when selected
    sqq_render/                # when selected
      sqq_cage.gro
      sqq_cage.xtc
      sqq_cage.membership.tsv
      sqq_cage.vmd.tcl
    track/
```

If more than 26 topologies are found, SQQ warns and switches the whole multi-GRO run to information-only output. It still analyzes every readable GRO, but writes only the root `sqq_config_resolved.yaml` and `result/info/*_info.md`; summary XLSX/CSV/detail files, ordinary GRO files, and the complete `sqq_render/` bundle are suppressed. This safety override has precedence over engine defaults and configured output requests.

For every normal multi-frame or multi-file result root, Markdown and optional membership/order TSV reports are placed under `info/`; selected per-frame structure files are placed under `gro/<frame>/`. Serial, thread, and process execution use the same placement. No per-frame directory is created merely to hold one Markdown file.

The four files in `sqq_render/` form one visualization package. With the default `render.atom_scope: full`, `sqq_cage.gro` contains the complete input-frame atom topology and first selected frame, while `sqq_cage.xtc` contains every atom coordinate and box for every selected render frame, with the original physical frame times when available. This includes water hydrogens, complete guests, additives, environment/wall components, and other retained atoms. `render.atom_scope: compact` selects the legacy water-oxygen plus complete-guest topology. `sqq_cage.membership.tsv` uses five record types: `F` maps render frames to source frames, time, and effective graph mode; `C` stores every cage center after orthorhombic PBC wrapping in explicit angstrom coordinates; `G` stores every complete guest-molecule atom group, including guests outside cages; `M` stores cage/guest membership atoms with cage type and optional phase, domain, and cluster identifiers; and `P` maps rendered atoms to component roles and residue names for context selection. Component topology is fixed, so `P` records are emitted only for render frame 0 after later frames have been validated against the same component signature. Large `P` groups are split into bounded rows, so full-system metadata remains readable for large trajectories. Analyze builds persistent IDs independently of rendering; when `sqq-render` is selected, it rewrites the `C` and `M` records with those IDs before atomic publication. “Sparse” applies only to membership metadata: the XTC retains every atom selected by `atom_scope` in every selected frame. Shared waters and guests assigned to several cages retain all memberships. Analyze and raw Track use the same behavior in SQQ-Py and SQQ-CPP; source Track inherits the topology and atom scope of its imported render bundle.

A Track command writes one independent result directory per requested target:

```text
result_track/
  sqq_config_resolved.yaml
  track/
    track_info.md
    track_state.json
    cage_observation.csv
    cage_track.csv
    cage_event.csv
    cage_population.csv
    guest_residence.csv
    lifetime_distribution.csv
    all/                       # or type_512, phase_sI, cage_t133
      track_info.md
      cage_observation.csv
      cage_track.csv
      cage_event.csv
      cage_population.csv
      guest_residence.csv
      lifetime_distribution.csv
      precursor_state.csv      # persistent-ID target only
      water_history.csv        # persistent-ID target only
      sqq_render/
        sqq_track.gro
        sqq_track.xtc
        sqq_track.membership.tsv
        sqq_track.vmd.tcl
```

Target tables preserve each selected cage's complete lifecycle, not just frames in which it has the requested type or phase. Lifetime rows report left/right censoring; population rows are per selected frame; event rows include birth, death, type/phase change, split, merge, and explicit gap events; guest residence remains non-exclusive.

SQQ permits only one active Analyze run per output root. A concurrent run stops before modifying results and asks for another `--output` directory. Worker fragments use a private run workspace; final render files are replaced atomically. Temporary-directory removal retries transient Linux/shared-filesystem `ENOTEMPTY`, `EBUSY`, and permission delays. If cleanup still fails, SQQ prints the retained temporary path but does not turn an otherwise completed analysis into a failure.

Each frame stages its Markdown, TSV, and selected GRO files privately and publishes the complete frame bundle only after every selected writer succeeds. Summary XLSX/CSV/detail files use the same recoverable publication rule. A failed write removes its staging data rather than exposing empty directories or a mixture of old and new generated files. Reusing an output root clears known SQQ artifacts from previous source names, grouped/flat layouts, and obsolete lettered topology roots, including F3/F4 output; unknown user files remain untouched.

Keep all four files together in `sqq_render/`, then source only the script from the VMD Tk Console:

```tcl
source {path/to/result/sqq_render/sqq_cage.vmd.tcl}
```

The Tcl script embeds a machine-readable render manifest naming its actual topology, trajectory, and membership files; the package still contains only four files. From a terminal, locate and validate one or more SQQ render packages and print absolute-path launch commands with:

```bash
sqq vmd path/to/result
```

`sqq vmd` searches the current directory when no path is supplied. It accepts a result directory, `sqq_render/` directory, or specific `.vmd.tcl` file, lists multiple packages in stable path order, and reports missing or empty required files without launching VMD. Commands use absolute VMD/Tcl-compatible paths and Windows or POSIX terminal quoting automatically. Existing Tcl files without an embedded manifest remain readable through their declared SQQ data paths.

Sourcing returns control to the Tk Console promptly, sets the background to white, prints a compact welcome, and loads the membership table plus GRO/XTC data through short event-loop tasks. `sqq help`, `sqq -h`, and `sqq --help` are available immediately; another SQQ command entered before initialization finishes reports `SQQ renderer is still loading`. Re-sourcing or a failed load uses one cleanup path to cancel pending callbacks, close the membership file, remove SQQ traces/graphics/representations, and delete only the molecule created by that SQQ script; help and the original error remain available after failure. The renderer becomes ready only after the default view and frame traces have been initialized successfully. It then reports the loaded frame count, displays the default opaque `sqq show cage all` view, and reports `SQQ graph: <effective-mode>` once. The graph line is printed again only if the effective mode changes. Use any of these equivalent commands for the full guide:

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
sqq show label [on|off]
sqq pick center|guest|off
sqq target save
```

Supported families are `cage`, `guest`, `phase`, `cluster`, `domain`, and `component`. Component targets are `all`, the roles `water`, `guest`, `additive`, `environment`, and `other`, or an exact residue name such as `KLN`. The default sourced view remains `sqq show cage all`; full-frame context is available but hidden until requested. Examples:

```tcl
sqq show cage all
sqq show cage 512
sqq show cage 512 51264
sqq show cage 512 guest 512
sqq show cage 512 51264 guest 512 phase sI

sqq show guest all
sqq show guest 512

sqq show phase all
sqq show phase sI boundary
sqq show cluster all
sqq show cluster cluster_00001
sqq show domain all
sqq show domain domain_00001
sqq show component environment
sqq show component KLN

sqq color cage 512 green
sqq color cage 51262_00053 yellow
sqq color guest 512 yellow
sqq color phase boundary orange
sqq color cluster cluster_00001 cyan
sqq color component KLN gray
sqq color cage all default
```

```tcl
sqq show label
sqq pick center
sqq target save
sqq pick off
```

The startup `sqq show cage all` view is a replaceable default. The first `sqq show ...` command after sourcing the script or after `sqq clear` replaces that default; later `show` commands add independent layers without removing earlier selections. One `show` may contain several family/target groups, and an exact repeated family/target selection is ignored rather than creating another VMD representation. `sqq show label` toggles labels; optional `on` or `off` sets an explicit state, and the historical misspelling `lable` is accepted. Labels are off by default and remain independent of picking. `sqq pick center` makes active objects transparent and creates yellow cage-center spheres/pick points; VMD automatically enters Atom Label mode for graphics picking; click a yellow center to make exactly that cage opaque. Do not manually select Query mode: Query only prints VMD information and does not send the callbacks SQQ needs. Water-atom clicks are ignored in center mode. `sqq pick guest` automatically enters VMD Pick mode, accepts a click on any guest atom, and highlights the complete guest plus every cage containing it. Cage highlights are yellow and guest highlights are orange. Both pick paths update persistent highlight representations instead of deleting the representation being clicked. Guests with no cage membership are reported without highlighting. `sqq target save` writes the current selected persistent cage ID or IDs to `sqq_target.txt` beside the render files. Center and guest modes are mutually exclusive. A frame change clears the transient selection and rebuilds the current-frame targets without disabling the chosen mode. `sqq pick off` exits pick mode but does not restore a previous VMD mouse mode; `sqq clear` removes custom show/color/label/pick state and restores the initial opaque cage-all view.

Each family token in `show` starts a new group and consumes the following targets until the next family token. For `cage`, a target is `all`, a registered cage type, or an exact cage ID such as `t133`; a fallback frame-local ID such as `51262_00053` is also accepted when persistent state could not be built. Generic types such as `4^1-5^10-6^2` also accept `4151062`. For `guest`, the same target identifies guests assigned to all cages, to a cage type, or to one exact cage ID. Phase targets are `all`, `sI`, `sII`, `sH`, `boundary`, `ambiguous`, `unclassified`, or `isolated`; cluster/domain targets are `all` or exact frame-local IDs. Multiple targets are accepted within each family group. The former inferred forms such as `sqq show 512` and `sqq color 512 blue` are not accepted.

Unlike `show`, `sqq color` accepts exactly one family per command. Colors accept a case-insensitive VMD color name, an in-range ColorID, or `default`. Cage and guest overrides are independent and persist across frame/selection changes until `sqq clear`, re-sourcing, or an explicit `default` reset. Cross-family layers always render as `phase -> cluster -> domain -> cage -> guest`, so guests remain last and visible regardless of `show` order. This family order is separate from the fixed cage-topology priority used for coincident cage edges and multi-cage guests. Cage networks use DynamicBonds with a 3.5 angstrom cutoff; guests use CPK and include the full molecule. A single cage layer uses a 0.125 angstrom cylinder radius (0.250 angstrom diameter); multi-type layers remain bounded from 0.125 to 0.130 angstrom.

The renderer processes only the current frame and creates a semantic display layer the first time it is needed. It then keeps that SQQ layer in one stable VMD representation across frame changes. A frame redraw recalculates only current-frame membership and updates the existing selection with `mol modselect` (or `none` when the layer is empty); it does not pre-render every frame or retain a frame-by-representation atom cache. Therefore a DynamicBonds radius, color, or material changed in VMD's Graphics Representations window on any frame remains effective when moving to earlier or later frames. A compatible semantic layer that first appears after such an edit inherits the matching family, cage-type, and radius-tier style instead of returning to the generated default. Explicit SQQ view/style changes or `sqq clear` may rebuild or reset affected SQQ layers, while representations added by the user remain untouched. Track scripts merge all initially selected target IDs before one default redraw instead of rendering once per target chunk. Rapid frame notifications are coalesced into one pending redraw. Fully unknown cage, cage-ID, guest-selection, cluster-ID, and domain-ID targets are rejected against the complete loaded trajectory; recognized phase names remain valid even when the current frame has no matching membership. Re-sourcing performs the complete SQQ-owned cleanup described above before resetting the generated script's selection/color state.

Cage identifiers are persistent `tID` values when a complete Analyze/Track state is available; cluster and domain identifiers remain deterministic frame-local classifications. Category selections (`phase`, `cluster`, or `domain`) and recognized phase labels simply report no membership when cluster analysis was not run; an explicit cage/type/cluster/domain target that never occurs anywhere in the loaded trajectory is rejected.

## Final Terminal Results

In an interactive terminal, SQQ redraws the same page after all reports and render files are safely published. It retains compact `Basic Information` and `Configuration`, replaces progress with `Analysis Results` or `Tracking Results`, combines all frame counts on one row and total/analysis/output time on another, and shows mean time when meaningful. A fully successful run omits status, and the result path is not repeated because the persistent `Output` row already identifies it. Redirected output is never cleared.

The final `Citation Recommendation` is a copy-ready manuscript sentence generated from one completed-feature evidence record, not merely from requested YAML settings. The terminal and summary dashboard consume that same record. For example, an applicable run may state `Cage types, populations, and guest occupancies were analyzed using SQQ.` Exact order parameters are named, VMD is mentioned only when `sqq-render` was published, occupancy is mentioned only when guests were actually available for evaluation, and Track mentions guest residence only when a multi-frame result contains guest observations. A one-frame Track result states only `Cage tracks were generated using SQQ.` Type transitions and lifetimes are claimed only for results containing at least two frames; without guest observations, that multi-frame sentence omits guest residence. Every completed page ends with the provisional publication and GitHub lines:

```text
Publication: J. Pang & Q. Sun, SQQ: Python Joint Toolkit for Water-Shell Topology Analysis, in submission.
GitHub     : https://github.com/pimooni/sqq
```

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

Ordinary per-frame GRO files are opt-in through YAML `output.type`. `cluster-gro` is separately opt-in and requires cluster search; no documented engine preset includes either category by default. With `output.type: [none]`, only `sqq_config_resolved.yaml` remains.

With YAML `run.strict: false`, standalone serial/process/thread read failures become failed summary rows and analysis continues where the reader remains usable. Failed inputs appear in `summary.xlsx/failures` and `<summary_csv_dir>/failures.csv` when their respective main-summary output types are enabled, and always in the mandatory `sqq_config_resolved.yaml` `run.failures` list. With `run.strict: true`, SQQ re-raises the error after updating `sqq_config_resolved.yaml` to `status: failed`.

GRO structure folders, filenames, and title lines use portable ASCII structure labels since version 0.2.4, for example `5^126^2` and `qc_5r_5^36^2_56566`. Markdown and main-summary scientific labels retain their readable superscript notation. This avoids Windows GBK/legacy-reader failures caused by Unicode superscript or subscript characters in generated GRO paths and titles.

Each `*_info.md` report starts with SQQ version, `SQQ engine: sqq-py` or `SQQ engine: sqq-cpp`, date/time, source, input format, topology when applicable, resolved sampling metadata for trajectory-like input, half/quasi search state, frame/time, requested-to-effective graph mode, effective bond mode, ring sizes, status, and molecule counts. It never formats the backend as `py (sqq-py)`. LAMMPS reports also record units, timestep, atom style, and type-map source.

When quasi-cage or cage isomers are present, the same report adds description tables:

- `Quasi Cage Isomer Description` explains each observed layered quasi-cage isomer by base ring and L1/L2/L3 ring sequence.
- `Cage Isomer Description` explains each observed closed-cage isomer by face composition and 6-ring face adjacency pattern.

`Cage Occupancy` remains a separate table because it describes guest assignment rather than cage topology. It expands exact guest compositions across dynamic columns in source guest order.

`summary-xlsx` and `summary-csv` use one shared main-table builder. SQQ-Py main output contains `summary`, optional `failures`, the effective connection table, `ring`, `half_cage`, compact composition-level `quasi_cage`, `cage`, optional `hydrate_cluster`, `order_parameter`, `ice`, and `detail_index` when detail files exist. SQQ-CPP emits the applicable subset: `summary`, optional `failures`, `cage`, `order_parameter`, and optional `detail_index`. XLSX stores these as sheets; CSV stores one UTF-8-SIG file per table under `output.summary_csv_dir` (default `summary/`). `summary-detail-csv` writes `cage_occupancy.csv` and `cage_isomer.csv` for both engines, plus `quasi_cage_isomer.csv` for SQQ-Py. `cluster-detail` adds `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. Main, detail, and cluster-detail CSV files share the same `summary/` directory; they have disjoint filenames and can be selected together. The first `summary` table is the compact dashboard, `failures` has one failed input/frame per row, and `detail_index` lists generated detail files. `cage_isomer.csv` defaults to observed nonzero isomer rows plus per-frame totals; YAML `output.cage_isomer_row: all` restores the zero-filled matrix. `order_parameter` contains only the selected F3, F4, Q_l, MCG, and DHOP columns; `--order-parameter none` omits it.

Summary construction records rows, columns, cells, bytes, CSV/XLSX write time, formatting time, and final-save time in `sqq_config_resolved.yaml -> run.summary_write`; the terminal prints its total seconds. The mandatory output-root `sqq_config_resolved.yaml` records final SQQ version, requested engine, effective SQQ engine, requested and effective graph modes, requested and resolved workers, normalized output types, input metadata, status/failures, and summary timing. Main CSV, XLSX, detail CSV, and `sqq_config_resolved.yaml` are written to same-directory temporary files and atomically replaced on success or failure. A completed run also writes `sqq_output_manifest.json`; later cleanup removes only paths recorded as SQQ-owned and leaves unrelated files in a reused result directory untouched. XLSX sheets above 200,000 cells or 128 columns keep header styling, filter, freeze pane, and fixed column widths but skip costly body-cell formatting; scientific values and table schemas are unchanged.

The `hydrate_cluster` main-summary table reports the mutually exclusive `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`. Optional cluster-detail CSV records add the corresponding cage-id groups and `boundary_composition`; hydrate-domain CSV records expose only external boundary contacts through `external_boundary_contact_count` and `external_boundary_contact_ids`.

Output ownership is:

```text
cage > quasi_cage > half_cage > ring
```

SQQ-Py cage files include cage waters, `CNT` center pseudoatoms, and assigned guests. SQQ-CPP cage files omit the synthetic `CNT` center pseudoatom. Exact guest-composition files are generated from the guest names present in the frame, such as `CH4`, `CH4x2`, or `CH4+CO2`.

See `docs/design.md` for algorithm details and `docs/update.md` for release changes.
