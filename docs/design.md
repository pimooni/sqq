# SQQ Development Design

SQQ means **Shell Quant Qualifier**. This document records the current implementation logic for developers, so the code and the scientific definitions stay aligned.

## Pipeline

Modes `00` and `50` run the complete SQQ-Py pipeline:

```text
input frames
  -> molecule selection
  -> water graph: hydrogen bond / O-O / user pair map
  -> diagnostic coordination distribution
  -> chordless rings (default) or optional shortest-path rings
  -> half_cage and quasi_cage open patches
  -> closed cage search and guest occupancy
  -> optional hydrate_cluster analysis from all detected cages
  -> F3/F4/Q_l plus MCG/DHOP order parameters and ice metrics
  -> per-frame outputs, selected main summary format, and optional detail CSV
```

Modes `99` and `cpp` use the same Python shell around a reduced native frame pipeline:

```text
Python input/config/selection
  -> C++ graph
  -> C++ internal chordless 4/5/6 rings
  -> C++ cage topology, isomer, and automatic occupancy
  -> C++ F3/F4
  -> Python annotated cage GRO/VMD, selected cage GRO, compact info, and summary CSV or optional XLSX
```

In SQQ-Py, the shared water graph is used by ring, half_cage, quasi_cage, cage, selected F3/F4, graph-mode Q_l, and ice analysis. Selected MCG and DHOP descriptors are calculated during the order stage but use dedicated guest/water cutoff graphs because their published definitions are independent of the selected SQQ bond mode. Hydrate_cluster analysis starts after cage detection and uses all detected cage-ring memberships, not the raw water graph or the report-filtered cage list. In both engines, the graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode. Coordination diagnostics read this graph without adding, removing, or capping edges.

## Analysis Modes and Workers

Modes are discrete presets rather than a continuous 00-99 scale. Modes `00` and `50` select SQQ-Py; modes `99` and `cpp` select SQQ-CPP. The command default remains `50`.

| Mode | Engine | Graph | Search sizes | Auto worker fraction | Find cluster | Default output types |
| --- | --- | --- | --- | --- | --- | --- |
| `00` | SQQ-Py | `hbond` | 4/5/6 | 100% | on | `info,gro,sqq-cage-gro,sqq-render,summary-xlsx,cluster-gro` |
| `50` | SQQ-Py | `auto` | 4/5/6 | 50% | off | `info,sqq-cage-gro,sqq-render,summary-xlsx` |
| `99` | SQQ-CPP | `hbond` | internal 4/5/6 | 100% | unsupported | `info,gro,sqq-cage-gro,sqq-render,summary-csv` |
| `cpp` | SQQ-CPP | `auto` | internal 4/5/6 | 50% | unsupported | `info,sqq-cage-gro,sqq-render,summary-csv` |

Mode application order is:

```text
built-in defaults -> mode preset -> config.yaml -> explicit CLI options
```

Each preset owns graph mode, ring sizes, automatic worker fraction, initial cluster state, and default output types. Config and explicit CLI values can override compatible settings. `quasi_cage.max_layers` and `order.parameters` remain independent. Mode `00` includes split `cluster-gro`; enabling cluster search explicitly in mode `50` populates selected info and main-summary outputs but does not add an unselected output type. A disabled search writes no cluster GRO/detail output. `sqq-render` implies `sqq-cage-gro`.

Modes `99` and `cpp` share the SQQ-CPP feature boundary. Both use internal chordless 4/5/6 rings, `f3,f4`, no cluster/quasi/ice analysis, and no silent Python fallback. Mode `99` chooses the hbond/100% performance preset and ordinary classified GRO; mode `cpp` chooses auto/50% and deliberately does not select `gro` or `cage-gro` by default. Compatible config/CLI overrides remain available; incompatible requests fail before analysis.
Cluster-search precedence is:

```text
--find-cluster > hydrate_cluster.enabled in config.yaml > mode preset
```

An explicit `-b` / `--bond-mode {auto,hbond,oo,pairs}` overrides the graph mode from both the preset and `config.yaml`. `--pairs PAIRS.txt` implies pairs mode unless `-b pairs` is already given; it cannot be combined with another explicit bond mode.

`parallel.workers: auto` calculates `floor(physical_core_count * mode_fraction)`, then reserves one physical core for the operating system and caps the result by the number of independent files or selected trajectory frames. Physical-core detection prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker N` / `-w N` overrides the mode fraction using form-based parsing: integer text such as `1`, `4`, or `100` is an explicit worker count; decimal text such as `0.5` or `1.0` and percentages such as `50%` or `100%` are physical-core fractions. Thus `-w 1` means one worker, while `-w 1.0` and `-w 100%` mean all detected physical cores before the reserve-one-core clamp. Percentages above `100%` and decimal values above `1.0` are rejected. The old `--workers` spelling is retained as a hidden compatibility alias. Worker resolution remains capped by task count and the Windows `ProcessPoolExecutor` limit. `parallel.backend` defaults to `process`; `thread` is retained for compatibility and `serial` forces one process.

### SQQ-CPP Native Backend

#### Ownership boundary

SQQ-CPP is an in-process compiled backend behind the existing CLI. It is not a standalone executable and does not duplicate file parsing or reporting:

```text
Python
  parse CLI/config
  read GRO/XYZ/XTC/TRR or orthorhombic LAMMPS DATA + dump/DCD
  select waters, guests, and other molecules
  resolve pair identifiers and worker tasks
        |
        v
pybind11 adapter (one normalized frame)
        |
        v
C++17
  graph -> chordless rings -> cages/isomers/occupancy -> F3/F4
        |
        v
Python
  rebuild SQQ models -> annotated GRO/VMD / selected cage GRO / info.md / summary CSV or optional XLSX / run_config.yaml
```

The native binding releases the GIL across the compute call. The `thread` scheduler is deliberately not exposed for SQQ-CPP. `process` and `serial` retain the established independent-file/frame scheduling semantics; one individual frame remains one native task.

#### Input data contract

The Python adapter passes only normalized, index-based data:

- all atom positions in nm as finite Cartesian triples;
- each water's oxygen index and available hydrogen indices;
- each selected guest's residue id/name, atom indices, and optional center atom;
- an optional three-length orthorhombic box, or no box for non-periodic input;
- resolved water-index pair edges for `pairs` mode;
- graph thresholds, selected 4/5/6 ring sizes, cage limits/validation thresholds, occupancy settings, and F3/F4 switches.

The native return mapping contains the effective bond mode, sorted graph edges, canonical internal rings, cages with type/waters/face-ring indices/center/guest indices/isomer, optional per-water F3/F4, and warnings. The adapter reconstructs `GraphResult`, `Ring`, `Cage`, `F3F4Result`, and `FrameResult` objects and assigns deterministic public object IDs. The normal Python writers therefore remain the single owner of file schemas.

#### Native algorithm scope

Graph construction implements the same `auto`, `hbond`, `oo`, and `pairs` definitions used by SQQ-Py. `auto` selects hydrogen bonds only when all selected waters have usable hydrogen coordinates; otherwise it selects O-O connectivity. Orthorhombic minimum images and deterministic cell candidates are followed by exact distance/angle checks and stable edge sorting.

Ring search enumerates canonical chordless cycles only for selected sizes 4, 5, and 6. Those rings are returned through the internal contract so cage face membership, isomers, and structure output can reuse the Python models, but SQQ-CPP suppresses public ring tables and ring files.

Cage growth uses the same generated face-composition constraints, Euler-compatible edge/vertex incidence, deterministic state ordering and limits, and closed-shell acceptance rules as the Python reference. Both engines always require an edge-connected, trivalent manifold shell. Optional scientific validation adds face-planarity/edge-variation, projected-area, positive-volume, and volume-centroid checks. Cage isomers use the same six-ring face-adjacency signature.

Occupancy is always part of supported cage analysis. With selected guests, candidate centers are tested by the configured polyhedron definition and assigned deterministically. With no selected guests, the report state is `not evaluated`; this is distinct from an evaluated cage set in which every cage is empty. No occupancy CLI switch is added.
Selected native cage GRO files preserve standard atom records and the source box. They contain cage waters and assigned guests without adding the synthetic `CNT` cage-center pseudoatom used by SQQ-Py. Separately, `sqq-cage-gro` writes the complete source atom set with compact membership annotations; it is shared by both engines and does not add cage centers.

F3 and F4 are independently selectable. F3 uses the active graph. F4 requires usable hydrogen geometry; a frame without it reports no valid F4 values and a warning. `all` expands only to the native supported pair `f3,f4`, and `none` skips the order result.

#### CLI and validation boundary

Accepted control groups are:

- `-i`, `-t` / `--top` / `--topology`, `-o`, `-c`, directory/glob/trajectory/input controls, `--strict`, and `--output-layout`;
- `-b auto|hbond|oo|pairs`, `--pairs`, and `--pair-id`;
- `-s` restricted to a nonempty subset of 4/5/6, `--cage-size`, `--max-cage-face`, and `--cage-scientific-validation`;
- `--order-parameter f3|f4|f3,f4|all|none`;
- `-w` / `--worker` and `--parallel-backend process|serial`;
- `--output-layout grouped|flat`, `--cage-isomer-rows nonzero|all`, and `--output-type` restricted to `info`, `gro`, `cage-gro`, `sqq-cage-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, `all`, or `none`.

The C++ output normalizer keeps `gro` as the ordinary classified-GRO umbrella and makes `sqq-render` imply `sqq-cage-gro`. `all` maps to every supported native output. Mode `99` defaults to `info,gro,sqq-cage-gro,sqq-render,summary-csv`; mode `cpp` defaults to `info,sqq-cage-gro,sqq-render,summary-csv`, so it does not create classified `cage-gro` by default. `summary-detail-csv` is unsupported. Incompatible names are not silently removed. Public ring selection/output, size 7, shortest-path rings, every half/quasi option, cluster, ice, Q_l/MCG/DHOP, legacy per-frame VMD, detail/TSV output, fast closure, thread backend, and triclinic input are unsupported. Explicit CLI requests fail before frames are analyzed. A full SQQ config may retain unrelated built-in defaults, but an incompatible nondefault request fails. There is no `-m c` alias and no Python fallback if `_sqq_cpp` cannot be imported or returns an error.

#### Build and wheel architecture

`pyproject.toml` uses scikit-build-core. CMake builds the pybind11 module `_sqq_cpp` from `sqq/core/sqq-cpp/src/` with a C++17 compiler and installs it inside `sqq/core/`. The handwritten source and CMake metadata are included in the source distribution but excluded from binary wheels. Native build directories, CMake/Ninja state, extension binaries, wheels, source archives, and local test output are ignored; source `.cpp`, `.hpp`, and `CMakeLists.txt` files remain version-controlled.

Release CI uses cibuildwheel to compile and import-test CPython 3.10-3.14 wheels for Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, and separately builds a source distribution. A user installing a matching wheel receives an already compiled `.pyd` or `.so`; compilation occurs only for an explicit source build. The source path requires CMake 3.20 or newer, Python development headers, a C++17 compiler, and an appropriate platform build tool. Packaging configuration does not imply that a distribution has already been published.

#### Scientific parity contract

SQQ-Py is the regression oracle. Deterministic discrete results must match exactly for the same effective settings: graph edges, canonical ring nodes, cage type/water/face membership, cage-isomer label, and occupancy guest assignment. F3/F4 comparisons permit a small floating-point tolerance because the C++ and Python accumulation implementations are independent.

The native parity baseline established in 0.3.1 includes random graph/ring frames, random F3/F4 frames, synthetic cage/geometry/occupancy cases, package import tests, and a real 11,104-atom frame. On `tests/100.gro`, both engines produced 4,322 edges; 2,499 internal rings (45/2,147/307 for sizes 4/5/6); 339 cages across 16 types with exact water and face membership; 339/339 matching isomers and occupancy assignments; and 315 occupied cages. F3 matched exactly for 2,176 waters, while F4 had maximum absolute difference `4.44e-16`. The native core took about 0.5673 s versus 10.0967 s for the equivalent Python path on that host. This 17.8-fold core-path ratio is a benchmark, not an end-to-end guarantee.

### Process Execution Architecture

For two or more independent GRO/XYZ files, the default execution path is:

```text
main process
  -> normalize config and resolve input order
  -> validate unique frame-directory stems
  -> create spawn ProcessPoolExecutor
  -> initialize each worker once with config/output/stage queue
  -> maintain a rolling queue of at most 3 * workers tasks
  -> submit (input_index, path)

worker process
  -> report start
  -> read one frame
  -> run the ordinary frame pipeline
  -> write its own frame directory
  -> return (input_index, summary_row)

main process
  -> consume stage events
  -> reorder rows by input_index
  -> write selected main-summary CSV/XLSX, optional detail CSV, and run_config.yaml with resolved run metadata
```

`spawn` is selected explicitly on macOS, Windows, and Linux. This avoids forking a process after the interactive progress refresh thread exists and gives the same pickling/import contract on every platform. Worker callables and initializers are module-level functions. Only paths, raw trajectory indexes, small event tuples, and summary dictionaries cross process boundaries; atoms, rings, patches, and cages stay worker-local.

While the pool exists, the parent sets the common BLAS/OpenMP thread environment variables to `parallel.math_threads` (default 1). Spawned children inherit the limits before importing NumPy-backed modules. The parent environment is restored after pool shutdown.

For one XTC/TRR or supported LAMMPS trajectory with `-t` / `--top`, the parent opens the trajectory once to resolve `input.trajectory_stride`-selected raw indexes. Each worker initializer opens a private MDAnalysis Universe once, caches its immutable atom metadata once, and tasks contain a small contiguous batch of `(ordered_frame_index, raw_frame_index)` pairs. Batch size is `ceil(selected_frames / (4 * workers))`, clamped to 1 through 8. Parent and worker readers are explicitly closed. Multiple trajectory files and non-process trajectory backends remain serial. SQQ retains an orthorhombic box representation; non-90-degree trajectory angles are detected and rejected instead of being silently approximated.

Both standalone-file and indexed-trajectory process paths, plus the compatibility thread path, maintain at most `3 * workers` submitted tasks and refill the queue as futures complete. This is a bounded submission window, not a worker cap: 100 effective workers retain 100-way execution and at most 300 submitted tasks. It avoids constructing a Future and serializing arguments for every item in a very large input set.

Standalone files whose case-insensitive stems collide are rejected because the stem is the frame output-directory name. Non-strict serial, thread, and process read failures return failed summary rows; strict failures cancel queued work where possible and propagate to the main process. Indexed trajectory reader failure records the current failed frame and stops that unusable iterator. Output order is determined by the original index rather than completion order. Configuration is normalized before dispatch, so thread workers only read the shared mapping.

Terminal and main-summary dashboard metadata share the same display helpers. The requested graph mode is preserved from config/CLI, and the effective per-frame graph modes are collected from `connection_mode` in summary rows. Explicit graph modes display as `hbond`, `oo`, or `pairs`. Automatic graph mode displays as `auto -> hbond`, `auto -> oo`, or `auto -> mixed (hbond: N, oo: N)` when different frames resolve differently. Before frame analysis completes, the terminal header may show `auto -> pending`; the final run summary and the `summary` table in either main-summary format use the resolved display value.

Root `sqq` / `sqq -h` output renders the banner and product sentence, then `SQQ version: 0.3.3   Release date: Jul 22, 2026`, then the ordinary `usage:` line. Root `sqq -v` / `sqq --version` exits successfully after printing only that version line. Subcommand help retains the standard argparse layout.

`run_config.yaml` preserves normalized configuration values such as `graph.bond_mode: auto`, `hydrate_cluster.enabled: false`, `order.parameters: [f3, f4]`, `output.types: [info, sqq-cage-gro, sqq-render, summary-xlsx]`, and `parallel.workers: 1.0`, then adds a `run` block with status/error, frame totals, failure details, `sqq_version`, requested/effective graph mode, normalized order/output selections, worker policy, resolved workers, backend, math threads, and final `summary_write` timing/table-size metadata. The metadata reports per-table rows, columns, cells, bytes, write time, format mode/time, final workbook-save time, and total summary-write time. An initial `status: running` file is written before frame analysis; strict analysis or summary-write failures update it to `status: failed`. `output.types` is the only output selector; removed `output.disabled_outputs` configurations are rejected rather than migrated. Run-config, main CSV, detail CSV, and workbook writes use same-directory temporary files followed by atomic replacement, so a failed rewrite does not truncate an existing completed artifact.

Configuration is normalized once before execution. Textual booleans use explicit on/off parsing rather than Python truthiness; enum values are canonicalized and rejected when unsupported; cutoffs and scales must be finite and positive where required; counts, strides, and state limits must be positive integers rather than booleans or fractional numbers. Residue-name lists accept comma-separated text or sequences and are normalized deterministically.
### Terminal Progress Display

SQQ-Py serial and parallel runs share the same three-row stage model:

```text
file preparation       reading -> settings -> selecting
core topology search   graph -> ring -> half/quasi -> cage -> cluster
post-processing        filtering -> order -> ice -> output
```

`cluster` is included only when the resolved `hydrate_cluster.enabled` value is true, for example through mode `00` or `--find-cluster on`. When hydrate cluster analysis is disabled, the stage is omitted rather than shown as `cluster:0`.

SQQ-CPP keeps the same aligned three-row presentation but removes Python-only work:

```text
file preparation       reading -> settings -> selecting
native topology        graph -> ring -> cage
post-processing        order -> output
```

The public ring report is absent even though the native `ring` stage is required internally for cage construction.

The serial interactive panel renders the complete workflow and highlights the active stage with ANSI bold plus bright blue (`RGB(0,0,255)`). Stage columns are sized by the longest stage name in that column, so the display stays compact while `reading`, `graph`, and `filtering` remain aligned. The continuation marker `>` for a new row is placed before the aligned stage column. The timing row remains `stage / frame / total`.

```text
stage               : reading   > settings > selecting
                    > graph     > ring     > half/quasi > cage
                    > filtering > order    > ice        > output
stage / frame / total: 3.2 s / 18.7 s / 18.7 s
```

With hydrate cluster enabled, the second row ends with `cluster`:

```text
stage               : reading   > settings > selecting
                    > graph     > ring     > half/quasi > cage   > cluster
                    > filtering > order    > ice        > output
stage / frame / total: 3.2 s / 18.7 s / 18.7 s
```

Non-interactive serial output, such as redirected logs or a `tqdm` postfix, reports the short current stage label rather than ANSI highlighted text.

Parallel GRO/XYZ and indexed trajectory runs use a main-process progress aggregator. Spawned workers never write terminal control sequences; they send `start` and stage-transition tuples through a multiprocessing queue. The main process applies those events, ignores late events from already-finished tasks, and shows completed/failed/active/queued counts, compact `stage:count` rows, total elapsed time, and up to six active-file rows. Additional active files are summarized so high-worker modes do not fill the terminal.

For `stage_summary`, each column width is recalculated from the longest current `stage:count` cell in that column. Cells are left-aligned, and each column is followed by two spaces. This keeps the summary aligned without the wide `|`-separated cells used earlier.

```text
stage_summary       : reading:0    settings:0  selecting:0
                      graph:1      ring:2      half/quasi:0  cage:1
                      filtering:0  order:0     ice:0         output:0
```

With hydrate cluster enabled, `cluster` appears at the end of the core-search row and participates in the column-width calculation:

```text
stage_summary       : reading:0    settings:0  selecting:0
                      graph:1      ring:2      half/quasi:0  cage:1   cluster:2
                      filtering:0  order:0     ice:0         output:0
```

## Modules

- `sqq/pipeline.py`: top-level analysis order, config merging, frame loop, output dispatch.
- `sqq/core/cpp_backend.py`: validated frame/config adapter between Python models and the native module.
- `sqq/core/sqq-cpp/`: CMake project, C++17 core, pybind11 bindings, and native public header.
- `sqq/parallel.py`: spawned worker initialization, file/trajectory tasks, effective CPU detection, and math-thread environment control.
- `sqq/core/graph.py`: water graph construction with orthorhombic PBC and nearby-pair search.
- `sqq/core/ring.py`: non-recursive DFS ring search.
- `sqq/core/ring_topology.py`: shared per-frame ring incidence, adjacency, centers, distances, and optional face-quality metrics.
- `sqq/core/quasi_cage.py`: layered `half_cage` and `quasi_cage` search.
- `sqq/core/cage.py`: closed cage grow/fast-closure search, mandatory shell-topology validation, optional scientific geometry validation, and guest assignment.
- `sqq/core/hydrate_cluster.py`: all-detected-cage cluster graph, phase cores, domains, and boundaries.
- `sqq/core/f3f4.py`: F3/F4/Q_l order metrics.
- `sqq/core/spatial.py`: deterministic orthorhombic-PBC self/cross cutoff pairs plus a reusable point-radius cell index for occupancy candidates.
- `sqq/core/mcg.py`: MCG guest graph and MCG-1/MCG-3 largest components.
- `sqq/core/dhop.py`: DHOP plane-normal counts, water tagging, and DHOP35/DHOP30 components.
- `sqq/core/ice.py`: CHILL-style ice classification.
- `sqq/io/trajectory.py`: common dispatch plus standalone and worker-local indexed trajectory materialization.
- `sqq/io/lammps.py`: strict LAMMPS DATA/dump/DCD validation, explicit or inferred atom mapping, unit conversion, and stable frame normalization.
- `sqq/io/summary.py`: per-frame info, global workbook/CSV tables, output timing metadata, atomic replacement, and size-aware workbook formatting.
- `sqq/io/gro_writer.py`: grouped or flat GRO structure output.
- `sqq/io/vmd.py`: run-level annotated `sqq-cage.gro`, fragment merging, membership validation, and VMD Tcl generation.

## Input Validation and Coordinate Units

GRO declares coordinates in nm. SQQ requires the declared atom count, all corresponding finite atom records, one separate 3- or 9-value box line, and no extra non-empty records after that box. Each GRO file therefore contains exactly one frame. Atom and residue names must be non-empty. An all-zero box is normalized to non-periodic `None`; mixed nonpositive lengths are invalid. A nine-value box is accepted only when all six tilt terms are zero. Nonzero tilt and non-finite values fail fast.

XTC/TRR positions and lengths supplied by MDAnalysis are converted from angstrom to nm. Coordinates and cell angles must be finite, and angles must be 90 degrees within tolerance; triclinic frames are rejected. Non-finite or unparseable trajectory time is retained as unavailable rather than emitted as a numeric time.

LAMMPS trajectories use one shared Python boundary adapter before either engine runs. A LAMMPS DATA topology supplied through `-t`, `--top`, or `--topology` is mandatory. The adapter accepts `.dump` / `.lammpstrj` and LAMMPS DCD. A non-empty YAML `input.lammps.type_map` is authoritative and must map every atom type to `resname`/`atomname` or explicit `ignore`. When the map is absent or empty, the adapter derives H/O/C candidates from DATA masses and type comments, then validates molecule graphs from DATA Bonds. It recognizes exactly one O plus two H atoms joined by two O-H bonds as water, one C plus four H atoms joined by four C-H bonds as all-atom methane, and a clearly carbon-labeled unbonded singleton as a united-atom methane guest. Valid DATA molecule IDs are retained. If their grouping is invalid but the bond-connected components are all uniquely recognized, deterministic molecule IDs are rebuilt in atom-ID order and that fallback is recorded. Every occurrence of one numeric atom type must resolve consistently; ambiguous masses, mixed roles for one type, unsupported compositions, insufficient evidence, or non-unique reconstruction fail with a request for an explicit map. The resolved map and provenance are stored in `run_config.yaml`, per-frame `*_info.md`, and main-summary configuration. DATA atom IDs become SQQ atom IDs, dump rows are restored to stable topology atom-ID order, and `real`, `metal`, or `nano` lengths and times are converted to nm/ps exactly once. Every dump frame must declare `pp pp pp` and an orthorhombic cell; DCD cell angles must be 90 degrees. Tilt, mixed/nonperiodic boundaries, `units lj`, duplicate or changing atom-ID sets, and non-finite/nonpositive cells fail before topology analysis. The process scheduler gives each worker a private reader and immutable mapped atom metadata. SQQ-CPP receives the same normalized `Frame` as SQQ-Py and contains no separate LAMMPS parser.

GRO readers preserve optional fixed-width velocities. SQQ annotations after a semicolon are stripped before the ordinary record is parsed, allowing the generated visualization GRO to be read back by SQQ without changing its coordinate, identity, velocity, or box fields.

XYZ has no standard unit or box metadata. Coordinates are multiplied by `input.xyz_scale` / `--xyz-scale`; the default `0.1` assumes angstrom input, and `1.0` preserves coordinates already expressed in nm. SQQ accepts exactly one nonnegative declared atom count, exactly that many finite coordinate records, and no extra nonempty records; multi-frame XYZ must be split into files before analysis. XYZ remains non-periodic unless converted to a format carrying a box.


GRO water and guest molecules are grouped by contiguous residue blocks in source order. This preserves separate molecules when five-digit residue numbers wrap or an input reuses a residue ID later in the file. LAMMPS frames instead use resolved molecule IDs: valid DATA molecule IDs are retained, while invalid groupings may be rebuilt deterministically from unambiguous Bonds components. Interleaved dump atom rows therefore still form the correct molecules and molecule inventory. Fixed-column GRO parsing remains primary; the whitespace fallback recognizes digit-containing residue names such as `TIP3` without folding those digits into the residue number.

## Water-Graph Candidate Search

For `hbond` and `oo` modes, graph construction first requests oxygen pairs within `cutoff + 1e-7` from MDAnalysis `self_capped_distance`. MDAnalysis may select brute force, nsgrid, or periodic KD-tree internally. SQQ treats this only as a candidate generator: every pair is sorted deterministically, recomputed with the existing float64 orthorhombic minimum-image function, compared with the exact configured cutoff, and, in `hbond` mode, checked with the established donor-angle test. If the accelerated neighbor API is unavailable or rejects the input, SQQ uses the previous orthorhombic cell list.

This two-stage design prevents float32 candidate-boundary behavior from changing scientific edges while moving the broad neighbor search into compiled code. `pairs` mode bypasses geometric candidate generation. Non-orthogonal boxes are intentionally outside the current model.
## Coordination Diagnostics

The active graph is summarized by water-node degree. Per-frame outputs report degree 0, 1, 2, 3, 4, and greater than 4 as counts and fractions, together with mean coordination, the degree <=2 fraction, the four-coordinated fraction, and the over-four fraction.

The section title follows the resolved graph mode: Hydrogen-Bond Coordination, O-O Connectivity Coordination, or Pair Connectivity Coordination. These values are diagnostic only. They do not modify graph construction, ring/cage detection, F3/F4, Q_l, or ice classification.

## Order Parameters

Order-parameter calculation is selected through one normalized list:

```yaml
order:
  parameters: [f3, f4]
```

The equivalent CLI is `--order-parameter f3,f4`. Supported names are `f3`, `f4`, any non-negative `qN` degree such as `q6` or `q12`, `mcg1`, `mcg3`, `dhop35`, and `dhop30`. `all` expands to `f3,f4,q6,q12,mcg1,mcg3,dhop35,dhop30`; `none` disables the entire order stage output. Names are deduplicated and stored in canonical order. An explicit CLI list replaces the complete configured list rather than extending it:

```text
--order-parameter > order.parameters > default [f3, f4]
```

The pre-0.2.7 selectors `--no-q`, `-q` / `--q-degree`, `--mcg3`, and `--dhop30` remain hidden compatibility options. They emit a deprecation warning, and `--order-parameter` takes precedence when both forms are present. Legacy YAML enable flags are translated only when `order.parameters` is absent.

F3 and F4 follow the project reference implementation and use the active water graph as the neighbor map.

F3 and F4 are independently selectable. If only `f3` is selected, F4 is neither calculated nor written, and conversely for `f4`.

One frame-local graph-vector cache computes each undirected PBC bond vector once and stores both orientations. F3 and graph-mode Q_l share this cache. Other Q_l neighbor modes build the candidate list once per oxygen pair, and all requested degrees reuse the normalized bond vectors and spherical angles. Cached normalization constants remove repeated factorial work. These are calculation-sharing changes only; neighbor selection, accumulation order, thresholds, and reported values are retained.

Q_l is the local Steinhardt/LAMMPS-style bond-orientational order parameter:

```text
Ybar_lm(i) = (1 / Nb(i)) * sum_j Y_lm(theta_ij, phi_ij)
Q_l(i)     = sqrt(4*pi / (2*l + 1) * sum_m |Ybar_lm(i)|^2)
```

The implementation is independent Python code and does not copy LAMMPS source. It uses unweighted oxygen-neighbor bond vectors and the same rotationally invariant normalization as LAMMPS `compute orientorder/atom`. Q_l is opt-in in 0.2.7: selecting `q6,q12` computes the former default degree pair, while `q4,q6,q8,q10,q12` selects the common LAMMPS degree list.

Neighbor modes:

- `graph` (default): use the active SQQ water graph. In `hbond` mode this gives hydrogen-bond neighbors; in `oo` mode it gives O-O neighbors; in `pairs` mode it follows the user pair map.
- `cutoff`: use all water oxygens within `order.q_cutoff_nm`.
- `nearest`: use the nearest `order.q_n_neighbor` water oxygens within `order.q_cutoff_nm`.
- `lammps`: LAMMPS-compatible cutoff plus fixed-neighbor behavior; if `order.q_n_neighbor` is null, it defaults to `12`.

The non-graph modes reuse the deterministic orthorhombic/non-periodic cell-list pair search from `sqq/core/spatial.py` rather than scanning every O-O pair. Distances and vectors are still recomputed with the shared float64 minimum-image function and sorted deterministically.

When a fixed neighbor count is active and fewer than that number of neighbors are found inside the cutoff, every requested Q_l value is set to `0.0`, matching LAMMPS behavior. Without a fixed neighbor count, waters with no Q_l neighbors are omitted from the Q_l mean and count.

Q_l is a continuous structural descriptor, not a standalone ice-count classifier. The value is sensitive to the neighbor definition, so SQQ records selected Q_l names and neighbor settings in `run_config.yaml`, the terminal header, and selected main-summary output. `--q-neighbor-mode`, `--q-cutoff`, and `--q-n-neighbor` remain calculation settings and are used only when at least one `qN` parameter is selected.

## Hydrate Nucleation Order Parameters

MCG/DHOP are opt-in frame-local descriptors and are separate from the cage-topology `hydrate_cluster` hierarchy. Selecting `mcg1`, `mcg3`, `dhop35`, or `dhop30` activates only those requested outputs in the existing `order` terminal stage. Their results are stored in `HydrateOrderResult`, where `None` means not applicable and integer zero means that the calculation was applicable but no qualifying component was found. Numerical thresholds remain under `hydrate_order`; descriptor selection belongs only to `order.parameters`.

### Shared spatial search

`sqq/core/spatial.py` supplies deterministic self- and cross-cutoff pairs. For an orthorhombic box, coordinates are wrapped into cells whose widths are at least the cutoff; only the 27 neighboring cells are inspected. Candidate distances are recomputed in float64 with the same minimum-image function used elsewhere. Pairs are deduplicated and sorted, so process/serial output and tie breaking are stable. Non-periodic input uses the same cell scheme without wrapping. No fixed atom or neighbor array is used.

### MCG-1 and MCG-3

The selected guest set is controlled by `hydrate_order.mcg_guest_resnames` (default `CH4`, `MET`). A configured center atom is used when available; otherwise the residue is unwrapped around its first atom before its centroid is calculated. Original guest-list indices are retained for output membership.

For each guest pair A,B within 0.90 nm, candidate water oxygens are the intersection of waters within 0.60 nm of A and B. A water W is mutually coordinated when both conditions hold:

```text
((W-A) . (B-A)) / (|W-A| |B-A|) >= cos(45 degrees)
((W-B) . (A-B)) / (|W-B| |A-B|) >= cos(45 degrees)
```

All vectors use orthorhombic minimum images. A guest edge is accepted when at least `mcg_min_waters` waters satisfy both cones; the default is `>= 5`, not `== 5`. The qualifying edge graph is the only graph used for clustering. MCG-1 keeps nodes whose original qualifying degree is at least one. MCG-3 keeps nodes whose original qualifying degree is at least three. This is a one-pass degree filter, matching the published MCG-N convention; it is not a recursively peeled k-core. Connected components are then measured on the induced qualifying graph. Equal-size components are resolved by lexicographically smallest original guest indices.

When no configured MCG guest residue is present, the value is not applicable (`N/A`). When guests are present but no edge qualifies, the largest cluster is zero.

### DHOP35 and DHOP30

DHOP builds a dedicated oxygen neighbor graph; it does not reuse `hbond`, `oo`, or `pairs` connectivity. The default cutoff is 0.35 nm for the all-atom TIP4P/Ice workflow of Li et al. A YAML value of 0.325 nm reproduces the original mW-water distance definition of DeFever and Sarupria. `35` and `30` denote plane-normal angle limits, not distance cutoffs.

For each undirected central oxygen bond j-k, SQQ combines every neighbor i of j other than k with every neighbor l of k other than j. The reference definition permits i and l to be the same common neighbor; retaining that case reproduces the published companion calculation. The normals are:

```text
n1 = (r_i - r_j) x (r_k - r_j)
n2 = (r_j - r_k) x (r_l - r_k)
cos(theta) = (n1 . n2) / (|n1| |n2|)
```

Zero-area planes are skipped and cosine is clamped to [-1,1]. For one central O-O bond, valid left/right normals are compared as a float64 matrix; values within `1e-12` of either angular threshold are recomputed with the scalar dot-product expression. A qualifying pair increments the planar-event count of both central endpoints, which is equivalent to traversing both directed orientations without duplicate geometric work. DHOP35 uses `theta <= 35 degrees`; optional DHOP30 uses `theta <= 30 degrees` from the same loop. This batching changes allocation and loop overhead only, not the descriptor definition or threshold behavior.

A water is initially qualified when its planar-event count belongs to `dhop_planar_counts` (default 11 or 12). A qualified water becomes a seed when at least three of its O-O neighbors are also qualified. Every seed and its complete first neighbor shell is tagged. The reported DHOP value is the largest connected component of tagged waters in the dedicated O-O graph; ties are deterministic. The transition-state value DHOP35=57 reported for one Li et al. system is not hard-coded.

References: Barnes et al. (DOI 10.1063/1.4871898), Knott et al. (DOI 10.1021/jp507959q), DeFever and Sarupria (DOI 10.1063/1.4996132), and Li et al. (DOI 10.1073/pnas.2011755117).

## Ring Search

Rings are searched on the already-built water graph; geometry is not used after the graph edge set is fixed. Sorted adjacency tuples and one node-to-bit map are built once. Every DFS state stores an immutable path plus an integer visited mask. The minimum node must remain the start, and a closing path is accepted only in one direction (`second < last`), eliminating rotational and reverse rediscovery before final canonical ordering.

With the default `ring.definition: chordless`, adding a node checks its graph neighbors already present in the partial path. A connection to any non-previous node is an immediate chord and prunes the branch. A connection back to the start closes the candidate immediately; the path is not extended beyond that closure edge. Random-graph regression compares this optimized traversal with the previous final-only chord test.

With `ring.definition: shortest_path`, every chordless candidate additionally applies the Franzblau shortest-path criterion. For each ring node, a bounded BFS is run only to `floor(size/2)`; the graph distance to every other ring node must equal the shorter distance along the cycle. A frame-local `(source, depth)` cache reuses these bounded distance maps across candidates. This opt-in definition can remove rings and therefore change patch/cage results.

Current behavior:

- supported ring sizes: 4, 5, 6, 7;
- default ring sizes: 5, 6;
- default definition: `chordless`;
- optional definition: `shortest_path`;
- ring nodes are water oxygen indices;
- `ring.sizes` / `--size` controls detection, while `ring.report_sizes` / `--ring-size` filters ring tables and GRO files after detection.

The historical `ring.primitive` default was not connected to the search implementation and is no longer emitted. Explicit definition now uses `ring.definition` / `--ring-definition`.

## Half-Cage and Quasi-Cage Terms

`patch` means a connected set of ring faces during search. It is not an output class by itself.

Layer definitions:

```text
L0 = base ring
L1 = side rings sharing every base-ring edge; L1 must close into a full side wall
L2 = rings grown outward from L1
L3 = rings grown outward from L2
```

L2 and L3 may be dangling rings or connected dangling ring chains. They do not need to close. The default configuration reports L1 quasi-cages and standard half-cages only; L2/L3 remain available by setting `quasi_cage.max_layers` or `--quasi-max-layer`.

`half_cage` is the standard subset of open patches:

```text
hc_5r_5^5
hc_6r_5^6
hc_6r_5^6_6^1
```

Any other valid non-closed layered open patch is reported as `quasi_cage`.

## Half-Cage and Quasi-Cage Algorithm

The search first precomputes:

```text
edge_to_rings: graph edge -> rings using that edge
ring_centers: locally unwrapped O-centroid for each ring
```

`edge_to_rings` is the primary topology filter. `ring_centers` are only used after topology filtering to order and limit candidates.

SQQ retains frame-local caches for symmetric ring-center distances, L2/L3 topology expansions, and patch geometry. These caches are discarded after the frame and therefore never mix topology or coordinates between trajectory frames.

Search order:

1. Choose one allowed ring as `L0`.
2. For every edge of `L0`, use `edge_to_rings` to find rings sharing exactly that base edge.
3. Sort those candidates by ring-center distance to `L0`; keep at most `quasi_cage.max_candidates_per_edge` per base edge unless the limit is disabled.
4. Build one compatibility map between each pair of adjacent L1 candidate lists. DFS accepts only compatible next rings, forward-checks that the following list still has an unused compatible ring, and requires the last side ring to connect back to the first.
5. Reject shifted or overlapped L1 choices with the expected unique-water count check.
6. Classify the L0+L1 patch:
   - if it matches a standard `half_cage`, store it as `half_cage`;
   - otherwise store it as `quasi_cage`.
7. For possible `hc_6r_5^6_6^1`, inspect L2 candidates and classify that larger standard patch as `half_cage`.
8. If `quasi_cage.max_layers >= 2`, grow L2/L3 from the current frontier:
   - find exposed frontier edges;
   - use `edge_to_rings` to get rings sharing those edges;
   - reject rings already in the patch;
   - reject rings touching lower-layer edges, so layer assignment uses the lowest possible layer;
   - sort by ring-center distance and keep at most `quasi_cage.max_layer_candidates`; repeated ring-pair distances are read from the frame-local distance cache.
9. Convert layer candidates according to `quasi_cage.search_policy`:
   - build candidate adjacency through `edge_to_rings`, so candidates are adjacent only when they share a complete graph edge;
   - `bounded` (default) uses small components whole and represents large components by single dangling rings plus deterministic local connected neighborhoods;
   - `exact` enumerates every connected subset up to `max_rings_per_layer` once by fixing the subset minimum ring id; subset, adjacency, and boundary membership use integer masks internally.
10. Cache ordered growth units by sorted patch IDs, sorted frontier IDs, and active limits. Deduplication keys encode patch/frontier membership as integer masks. Bounded mode retains patch-only state deduplication for compatibility. Exact mode deduplicates `(patch, frontier)` so the same patch reached with a different last layer remains a distinct L3 growth state.
11. Classify each new layered patch again as `half_cage` or `quasi_cage`.
12. Check the final ring-set key before constructing geometry. Only a new patch performs PBC unwrapping and center calculation; identical ring sets reuse frame-local geometry.

The algorithm avoids scanning all rings after `edge_to_rings` is built. The default `bounded` policy retains existing candidate caps, layer definitions, deterministic ordering, and `max_layers = 1`. The opt-in `exact` policy may add connected L2/L3 half-layers. Candidate ranking, wall-combination, and layer-state truncation set explicit warnings. Patch maximality and free-patch ownership use ring-to-owner inverted indexes rather than all-pairs subset scans.

Important limits:

- `quasi_cage.max_layers`: default 1 for fast routine analysis; use 2 or 3 to report outer dangling quasi-cage layers.
- `quasi_cage.max_rings_per_layer`: maximum rings in one L2/L3 growth unit.
- `quasi_cage.max_layer_states_per_seed`: per-seed cap for finite growth.
- `quasi_cage.max_candidates_per_edge`: L1 candidate cap per base edge.
- `quasi_cage.max_layer_candidates`: L2/L3 candidate cap per frontier.
- `quasi_cage.search_policy`: `bounded` by default; `exact` preserves frontiers and enumerates connected layer subsets.

## Cage Algorithm

The default cage mode is:

```text
cage.search_mode = grow
cage.seed_mode = ring
cage.fast_closure = true
cage.fast_closure_max_states = 20000
cage.scientific_validation = false
```

The code can also grow from half-cage/quasi-cage patch seeds with `cage.seed_mode = patch`, but ring seeds remain the default. Indexed two-to-four half-cage closure is enabled as a supplementary path. Basic topology validation is unconditional. The stricter scientific validation phase is opt-in and changes only geometric acceptance and the center definition.

### Search Scope and Report Scope

`ring.sizes` / `--size` defines the shared face-size search universe. Ring and quasi-cage detection support 4/5/6/7. Cage detection intentionally uses only the 4/5/6 intersection of that universe.

For cage search, SQQ generates every trivalent Euler-compatible face composition up to `cage.max_faces` / `--max-cage-face`:

```text
2*n4 + n5 = 12
n4 + n5 + n6 <= max_faces
```

All generated compositions are searched in one merged grow traversal. Reports and documentation use scientific face-count notation:

```text
5¹²
5¹²6²
5¹²6³
5¹²6⁴
5¹²6⁸       # Type H large cage
4³5⁶6³      # Type H small cage
```

Other accepted compositions use generic labels such as `4^1-5^10-6^2`.

Detection and reporting are separate:

- `all_cages` contains every accepted closed cage in the search scope;
- `cage.report_types: auto` reports every detected cage allowed by `ring.sizes` / `--size`;
- `cage.report_types` / `--cage-size` filters the user-facing cage counts, occupancy, GRO, info, and main-summary tables;
- cage report groups expand to exact compositions and duplicate types are removed:

```text
I     -> 5¹², 5¹²6²
II    -> 5¹², 5¹²6⁴
H     -> 5¹², 5¹²6⁸, 4³5⁶6³
HS-I  -> 5¹², 5¹²6², 5¹²6³
TS-I  -> 5¹², 5¹²6², 5¹²6³
I2II  -> 5¹²6³
```

- `--cage-size` accepts `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, and `I2II`; group names may be comma-separated;
- `--cage-size all` reports every detected composition;
- all detected cages, including unreported types, still remove consumed half-cages, quasi-cages, and free rings;
- an explicitly requested cage group is rejected when one of its compositions requires a face size absent from `--size` or exceeds `--max-cage-face`.

The default report scope is `auto`. Therefore, `-s 4,5,6` searches and reports every accepted cage composed of 4-, 5-, and 6-membered faces unless `--cage-size` explicitly narrows the report. For example, `-s 4,5,6 --cage-size I,II` keeps 4/5/6 ring and quasi-cage reporting while cage output is restricted to 5¹², 5¹²6², and 5¹²6⁴.

### Shared Ring Topology and Fast Closure

After ring detection, the pipeline builds one frame-local `RingTopologyIndex`. It contains the stable ring-id map, locally unwrapped ring centers, `edge_to_ring_ids` incidence, ring adjacency, and a symmetric ring-distance cache. Quasi-cage and cage search receive the same object. Scientific validation additionally stores face planarity RMS, edge-length coefficient of variation, and projected area. Hydrate-cluster analysis additionally requests cached least-squares ring normals and reuses center/normal pairs during physical shared-face resolution. When both a normal and face-quality metrics are requested, one SVD supplies both results.

When generic grow reaches a configured state limit and `cage.fast_closure = true`, SQQ builds a sparse patch-overlap graph from standard half-cages. Two patches are connected when they share a face or an exposed boundary edge. A deterministic indexed traversal joins connected combinations of two, three, or four patches, rejects face-count overflow and edge overuse incrementally, and submits only exact target compositions to the normal closed-polyhedron validator.

Generic grow runs first and remains the normal path for arbitrary generated compositions and cages not expressible by available half-cage patches. If grow finishes without reaching a state limit, fast closure is skipped entirely. Recovery results are appended only when their `(cage_type, water_set)` was not already found. Therefore an exhaustive grow run keeps the same cage set and object ids; a bounded grow run can gain a standard cage that would otherwise be missed. `cage.fast_closure_max_states` bounds only this supplementary patch-combination traversal. Reaching it produces a warning without discarding generic grow results.

### Grow Logic

1. Reuse `edge_to_ring_ids` and ring centers from the frame-local topology index; assign stable integer bits to active rings, graph edges, edge-to-ring candidate sets, and generated cage targets.
2. Start from a seed face set, usually one ring, and derive its face-membership bitset plus disjoint `edge_once` and `edge_twice` masks.
3. Compute the target bitmask whose per-size counts can still contain the seed counts. If no single target remains compatible, prune immediately.
4. Pop a state `(face_ids, face_mask, edge_once, edge_twice, (n4,n5,n6), used_face_incidence, compatible_target_mask)` from iterative DFS. The fixed 4/5/6 tuple avoids copying a face-count dictionary for every branch.
5. If `edge_once` is empty, classify by exact face counts and submit the shell to ordinary polyhedron validation.
6. Otherwise choose the open boundary edge with the fewest addable rings (minimum remaining value) using `edge_to_ring_ids`.
7. Reject candidates already in `face_mask`, below single-ring seed rank, intersecting `edge_twice`, or exceeding every compatible target. Before expansion, also reject a target when its precomputed total face incidence minus the state incidence is smaller than the number of open edges or has incompatible parity; these are necessary conditions for closing every open edge exactly once.
8. Adding a ring promotes `ring_edges & edge_once` to `edge_twice` and toggles its remaining edges into `edge_once`; no edge-count dictionary is copied.
9. Update 4/5/6 counts and intersect the compatible-target mask with precomputed count-allowance masks. An empty target mask stops the branch.
10. Calculate ring-center distance only when `cage.max_boundary_candidates` truncates a larger topology-filtered candidate list; record a warning when this occurs.
11. Continue until closure, state budget, or target infeasibility. Face masks are the global duplicate-state keys.

Single-ring seed-rank pruning avoids rediscovering the same cage from every face: later growth cannot add a ring whose stable id ranks before the seed ring.

Default acceptance criteria:

- every edge is used exactly twice;
- Euler characteristic satisfies `V - E + F = 2`;
- the face-adjacency graph is edge-connected;
- incident faces around every shell vertex form one cyclic vertex link;
- every shell vertex is trivalent;
- face counts match one generated target composition;
- the same `(cage_type, water_set)` was not already accepted.

These topology checks are always applied before cage type/isomer assignment in SQQ-Py and SQQ-CPP. They reject disconnected, pinched, branched, and non-manifold shells independently of `cage.scientific_validation`.

### Optional Scientific Geometry Validation

`cage.scientific_validation = false` is the default and is independent of analysis mode. `off` disables only the additional geometry checks; it never bypasses the mandatory topology checks above. `--cage-scientific-validation on` additionally requires:

- each ordered ring face is locally unwrapped across PBC and fitted by SVD; its planarity RMS must not exceed `cage.max_face_planarity_rms_nm`;
- its cyclic O-O edge-length coefficient of variation must not exceed `cage.max_face_edge_cv`, and its projected area must be nonzero;
- the consistently outward-triangulated shell must have volume at least `cage.min_cage_volume_nm3`.

An accepted scientific-validation cage uses the tetrahedral volume centroid of the oriented triangle shell. The default path continues to use the mean of locally unwrapped cage-water coordinates. Consequently, enabling scientific validation may remove geometrically distorted cages and may change guest occupancy or geometry-resolved hydrate-cluster edges. The now-mandatory topology checks can reduce cage, isomer, occupancy, and hydrate-cluster results relative to builds that accepted non-manifold shells. Raw ring and half/quasi searches, order parameters, and ice classification are unchanged; ownership-filtered free-ring and free-patch outputs can increase when rejected cages no longer consume them.

The built-in scientific thresholds are `0.06 nm` planarity RMS, `0.35` edge-length CV, and `1.0e-6 nm^3` minimum volume. They are explicit configuration values rather than hidden constants.

## Hydrate Cluster

Hydrate_cluster is optional. Its built-in value is off, and the overall default remains off because the default mode is `50`:

```text
hydrate_cluster.enabled = false
hydrate_cluster.min_cage = 2
```

Modes `00` and `50 --find-cluster on` can run cluster analysis; modes `99` and `cpp` reject it. The command-line controls are `--find-cluster on/off` and `--cluster-min-cage N`; the explicit find-cluster value has precedence over config and mode. Cluster search populates every selected `info` and main-summary output but does not add an unselected output type. Mode `00` includes `cluster-gro` in its preset, whereas mode `50 --find-cluster on` does not. Output type `cluster-detail` controls the optional domain and one-row-per-cluster CSV files. Explicit `cluster-detail` or `cluster-gro` requires resolved cluster search on.

Hydrate_cluster uses `result.all_cages`, the complete detected cage set in the selected ring/search scope. `cage.report_types` / `--cage-size` filters cage counts, occupancy, GRO, Markdown cage tables, and main-summary cage columns only; it does not filter the cluster graph or phase evidence. Cluster hierarchy/detail/domain records resolve cage IDs against the same complete set, so an unreported cage can remain topologically necessary without appearing in the report-scoped cage table.

The high-level hierarchy is informed by HTR+ ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)): classify hydrate type and polycrystalline boundaries on a cage-connection graph. SQQ does not copy the HTR+ implementation; it uses the explicit shared-face fingerprints and deterministic domain rules below.

### Physical shared-face cage graph

1. Use `result.all_cages` (falling back to `result.cages` only for legacy result objects).
2. Build `ring_id -> cage_ids` from each cage face list.
3. Treat one complete shared ring as a potential cage-cage edge.
4. A physical ring face can separate at most two cages. When ring geometry is available, locally unwrap the face, fit its plane, and select the best cage center on each side. When geometry is unavailable, accept only a face referenced by exactly two cages.
5. Find deterministic connected components in the resulting undirected graph.
6. Report components with cage count >= `hydrate_cluster.min_cage`; count cages in smaller components as isolated cages.

### Local phase fingerprints and strict seeds

For every graph node, SQQ counts first-shell labels of the form `(neighbor cage type, shared face size)`. Strict seeds reject unexpected label types and allow each expected count to differ by at most one.

The sI templates are:

- `5^12`: twelve `5^12 6^2` neighbors through 5-ring faces;
- `5^12 6^2`: four `5^12` neighbors through 5-ring faces, eight `5^12 6^2` neighbors through 5-ring faces, and two `5^12 6^2` neighbors through 6-ring faces.

The sII templates are:

- `5^12`: six `5^12` and six `5^12 6^4` neighbors through 5-ring faces;
- `5^12 6^4`: twelve `5^12` neighbors through 5-ring faces and four `5^12 6^4` neighbors through 6-ring faces.

The sH templates are:

- `5^12`: four `5^12`, four `4^3 5^6 6^3`, and four `5^12 6^8` neighbors through 5-ring faces;
- `4^3 5^6 6^3`: three `4^3 5^6 6^3` neighbors through 4-ring faces, six `5^12` neighbors through 5-ring faces, and three `5^12 6^8` neighbors through 6-ring faces;
- `5^12 6^8`: twelve `5^12` neighbors through 5-ring faces, six `4^3 5^6 6^3` neighbors through equatorial 6-ring faces, and two `5^12 6^8` neighbors through axial 6-ring faces.

These counts are shared-face incidences. They need not be distinct cage ids in a minimal periodic cell. Their balance is consistent with the ideal sH cell ratio of three `5^12`, two `4^3 5^6 6^3`, and one `5^12 6^8`.

The existing conservative sH composite is retained as supplemental high-confidence evidence: two nonadjacent `5^12 6^8` anchors, exactly six common `5^12` cages connected to both anchors through 5-ring faces, exactly six `4^3 5^6 6^3` cages, and at least one adjacent medium-cage bridge between the anchors. It supplements rather than replaces the per-cage fingerprints.

Seeds are internal phase evidence. `seed_count` is the number of strict seed anchors contained in a domain; `seed_cage_count` is the number of unique domain cages covered by their overlapping seed neighborhoods.

### Expansion and exclusive domains

sI, sII, and sH expand independently from the union of their strict seed members. A growth candidate must:

- use a cage type supported by the phase template;
- have at least one compatible internal fingerprint label;
- not exceed any expected fingerprint count by more than one;
- connect through a face label allowed by both endpoint templates; and
- receive at least two compatible contacts from already accepted phase cages.

For sH, the same edge check includes pentagonal `5^12` contacts, square medium-medium contacts, equatorial medium-large hexagonal contacts, and axial large-large hexagonal contacts. Composite-seed internal edges remain trusted seed evidence.

After all phases collect claims independently, SQQ forms domains from cages claimed by exactly one phase. Domain edges must remain phase-compatible, and every connected domain component must contain at least one strict seed anchor. Same-phase regions separated by non-domain cages remain separate domains. Domain ids are deterministic within a frame and are not tracked across frames.

### Boundaries and cluster type

After exclusive domains exist, SQQ partitions every reported cluster into mutually exclusive cage categories:

1. `classified_cage_ids` is the union of all final sI, sII, and sH domain members.
2. For each non-domain cage, inspect its complete shared-face neighbors in the same connected component.
3. If at least one neighbor is a final domain member, place only the non-domain cage in `boundary_cage_ids`.
4. Do not mark the contacted domain cage. A direct shared-face edge between two different phase domains leaves both endpoints classified in their original phases and creates no boundary cage by itself.
5. Stop after the first external non-phase layer. Boundary membership is not propagated through another non-domain cage.
6. A remaining non-domain cage with competing phase claims enters `ambiguous_cage_ids`; every other residual cage enters `unclassified_cage_ids`.

These four sets are pairwise disjoint and cover the complete cluster:

```text
classified | boundary | ambiguous | unclassified = all cluster cages
classified & boundary = empty
boundary & ambiguous = empty
boundary & unclassified = empty
```

Boundary is a generic topological category. SQQ does not create `sI-boundary`, `sII-boundary`, `sH-boundary`, phase-boundary context, or transition-cage categories. A non-domain cage touching more than one phase domain is still one boundary cage.

`HydrateDomain.boundary_cage_ids` stores the external non-domain boundary cages that directly contact that domain. Public domain-detail fields expose this relationship as `external_boundary_contact_count` and `external_boundary_contact_ids`; these are adjacency records, not additional cage classifications. Cluster detail records expose the four category ID sets and `boundary_composition`.

Neighboring cages can share vertices and complete face waters, so separate structure views may contain common water coordinates even though cage IDs are disjoint. Scientific ownership and regression checks therefore use cage IDs or the detected cage/ring edge graph rather than coordinate unions.

The mixed sI/sII real-GRO regression contains one 334-cage main cluster partitioned into 260 classified cages (66 sI and 194 sII), 69 boundary cages, 0 ambiguous cages, and 5 residual unclassified cages. Five additional cages occur in isolated or below-threshold components. All four main-cluster cage-ID sets are disjoint.

### Compact per-frame cluster hierarchy

Cluster reporting reuses the finalized `HydrateCluster` and `HydrateDomain` objects; it never reruns search, expands a domain, or changes cage ownership. When resolved cluster search is on and `info` is selected, `Frame Information` records `find_cluster: on` and the report adds one compact `Hydrate Cluster` table. Search off omits the section. Cluster search does not force `info` or a main-summary type, so an output selection without `info` creates no Markdown report.

Each `cluster_XXXXX` row reports the number of unique cage IDs in that connected component. Its children are deterministic sI/sII/sH domain rows, boundary, and compact unclassified topology, each subdivided by detected cage type. For this compact display only, `unclassified` is the deduplicated unresolved set: stored `ambiguous_cage_ids` and `unclassified_cage_ids` plus any uncategorized residual cluster cages. The main summary and cluster-detail CSV preserve those fields separately. Zero-count rows are omitted, multiple clusters are written sequentially, and `isolated` is one final top-level count with no cage-type children. Parent-child totals are conserved. Tree symbols decorate both `item` and `cage_qty`; Markdown columns are padded by display width.

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

The table contains no cage IDs, seed lists, confidence values, water membership, or guest records. Those remain in the existing plotting/detail outputs. The formatter changes presentation only and cannot change ring, cage, phase-domain, boundary, occupancy, order-parameter, or ice results.

### Native cluster GRO views

`cluster-gro` is a separate search-dependent output. Mode `00` selects it by default; mode `50 --find-cluster on` writes it only when selected explicitly. Search off writes no cluster GRO output; reusing an output directory removes stale grouped or flat cluster GRO files generated by SQQ. The ordinary `gro` umbrella remains limited to ring, half-cage, quasi-cage, cage, and ice files.

For each frame, output aggregates cage IDs across every reported cluster and every domain of the same category:

- sI is the union of all sI domain cage IDs;
- sII is the union of all sII domain cage IDs;
- sH is the union of all sH domain cage IDs;
- boundary is the union of all generic boundary cage IDs.

Ambiguous, residual unclassified, and isolated/below-threshold cage IDs are intentionally omitted. Each category maps its cages back to water membership, deduplicates waters, and writes complete water molecules in original frame order. Guests, CNT atoms, and other non-water molecules are never added. An absent category has no file unless `output.write_empty_files` is true.

Grouped layout writes:

```text
<frame>/hydrate_cluster/<frame>_cluster_sI.gro
<frame>/hydrate_cluster/<frame>_cluster_sII.gro
<frame>/hydrate_cluster/<frame>_cluster_sH.gro
<frame>/hydrate_cluster/<frame>_cluster_boundary.gro
```

Flat layout writes the same four canonical filenames directly under `<frame>/`.

Cluster GRO export is coordinate-preserving, not a visualization reimaging step. Every selected atom retains the exact wrapped coordinate from the analyzed frame and every category file copies the original box unchanged. Categories are never independently translated, centered, unwrapped, or made whole. Periodic and percolating networks can therefore retain apparent bonds across box faces: one single-copy GRO cannot remove every periodic seam without changing the original coordinate representation.

Scientific exclusivity applies to cage IDs. Neighboring cages from different categories share face waters physically, so their category GRO files may contain some of the same complete water molecules even though no cage ID belongs to more than one category.

A cluster is `sI`, `sII`, or `sH` when its domains contain one unique phase, `mixed` when multiple domain types occur, and `unclassified` when no domain exists. No phase or boundary label is inferred from cage composition alone.

Public motif output is not generated in the current release. The compatibility motif return slot remains empty, and neither a `Hydrate Motif` Markdown section nor a `hydrate_motif` main-summary table is written.

## Guest Occupancy

Cage centers are computed from the locally unwrapped O coordinates of the cage waters. Hydrogens are not used in the cage center.

For each accepted cage, SQQ checks all selected guest molecules. `guest.center_mode` is active:

- `center_atom` (default) and `auto` use the configured center atom when present. The built-in `guest.center_atoms` maps `CH4`, `CO2`, and `MET` to atom name `C`, so these residues use their carbon atom; a guest without a matching configured center falls back to a PBC-aware residue centroid;
- `centroid` always unwraps all guest atoms around one molecular anchor before calculating the centroid.

The same centroid helper is shared by occupancy and MCG and guest centers are precomputed once per frame. This fixes cross-boundary multi-atom guests whose raw arithmetic mean lies near the box center; occupancy values can intentionally differ from earlier releases for those molecules. A reusable orthorhombic-PBC cell index first selects only guest centers within the cage candidate radius, then performs the established exact minimum-image distance check. Guest iteration order and non-exclusive assignment to overlapping cages are unchanged.

Default `cage.occupancy_mode = polyhedron` triangulates the cage ring faces and uses a point-in-polyhedron solid-angle test. Candidate points are evaluated in float64 batches while triangle contributions are accumulated in the original face order; values numerically near the `2*pi` membership boundary fall back to the scalar implementation. Degenerate point-on-vertex triangles contribute zero solid angle rather than forcing an inside result. `center` uses only a center-distance cutoff, and `auto` accepts either method.

## Ownership Rule

Output ownership follows:

```text
cage > quasi_cage > half_cage > ring
```

Rules:

- patches consumed by any detected cage are not written as free `half_cage` or `quasi_cage`, even when that cage type is not reported;
- a `half_cage` whose ring set is a true subset of a `quasi_cage` is consumed by that `quasi_cage`;
- nested `half_cage` results are reduced to the larger `half_cage` patch;
- free rings are rings not consumed by cage, quasi_cage, or half_cage outputs;
- guest occupancy is only applied to closed cages.

## Output Layout

Output selection is positive and mode-specific. The four presets are:

| Mode | Default `output.types` |
| --- | --- |
| `00` | `info,gro,sqq-cage-gro,sqq-render,summary-xlsx,cluster-gro` |
| `50` | `info,sqq-cage-gro,sqq-render,summary-xlsx` |
| `99` | `info,gro,sqq-cage-gro,sqq-render,summary-csv` |
| `cpp` | `info,sqq-cage-gro,sqq-render,summary-csv` |

Supported SQQ-Py canonical names are `info`, `membership-tsv`, `order-tsv`, legacy per-frame `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `sqq-cage-gro`, `sqq-render`, `cluster-gro`, `summary-xlsx`, `summary-csv`, `summary-detail-csv`, and `cluster-detail`. `gro` expands to the five ordinary ring/half/quasi/cage/ice subtypes. `sqq-cage-gro` is an independent run-level visualization trajectory, and `sqq-render` implies it. `cluster-gro` is separate and requires resolved cluster search. `all` selects every applicable type; `none` selects none. Both keywords must appear alone.

SQQ-CPP accepts `info`, `gro`, `cage-gro`, `sqq-cage-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, `all`, and `none`. `gro` enables the supported ordinary classified cage GRO path; it is not collapsed in the recorded output selection. Mode `cpp` does not select either `gro` or `cage-gro` by default. `summary-detail-csv` and all cluster types are unsupported. `run_config.yaml` is mandatory in every mode.

`sqq-cage.gro` contains one complete GRO block per successful frame, ordered by the original input/frame index. All source atoms, identity/order, wrapped coordinates, box, and optional GRO velocities are preserved; no PBC reimaging is performed. Atom lines reserve the optional velocity columns and place `; SQQ1 m=...` at column 69. Membership records encode cage type and ID plus phase/domain/cluster IDs when cluster analysis exists; nonmembers use `m=-`. Multiple cage memberships remain attached to one oxygen rather than being flattened into one label. A multi-frame bundle requires identical atom identity and order across frames and is finalized atomically from worker-local fragments.

`sqq-render.vmd.tcl` is self-contained and loads the neighboring `sqq-cage.gro`. Because VMD does not treat concatenated GRO blocks as a trajectory directly, the script splits complete frames into a temporary directory, loads them in order, deletes the temporary files, and tracks frame changes. Cage topology is the default `sqq show all` view. Dynamic oxygen bonds use a 3.5 angstrom display cutoff. A single cage layer uses a 0.125 angstrom cylinder radius, giving a 0.250 angstrom displayed diameter; multi-type cage views use bounded 0.125–0.130 angstrom radii for deterministic overlap visibility.

The public interface is `sqq show <object...>`, `sqq color <object> <color>`, and `sqq help`. `sqq show all` selects every cage and is the startup view; `cage` is deliberately rejected as a `show` target. The bare names `phase`, `cluster`, and `domain` select every member of those families. Registered cage labels select one or more cage types; delimiter-free aliases map generic labels such as `4151062` to `4^1-5^10-6^2`; `<cage_type>_<five-digit-id>` and its compact alias select a frame-local cage; phase labels map `sI`, `sII`, `sH`, `boundary`, `ambiguous`, `unclassified`, and `isolated` to the annotation codes; and `cluster_<id>` / `domain_<id>` select reconstructed full IDs. Multiple explicit `show` targets must belong to one family, and a category target must be used alone. Category targets and recognized phase labels may produce an empty view when cluster analysis was not run; explicit cage/type/cluster/domain targets are validated against the complete loaded trajectory.

Color values are a case-insensitive VMD color name, an in-range ColorID, or `default`. Overrides are independent of the current show selection and persist during frame and selection changes in the current VMD session. Assigning a nondefault category color first removes all existing overrides in that family and then installs one category override; later object-specific assignments take precedence. An object-level `default` forces that object's stable built-in color, while category-level `default` removes all overrides in the family. Effective cage color precedence is exact cage ID, cage type, cage category, then the stable built-in cage palette; phase/cluster/domain use exact object, category, then stable or deterministic defaults. Cage displays expand selected cage types to frame-local cage IDs and regroup atom indexes by cage topology layer and effective ColorID. Nonstandard cage types form the lowest deterministic layers; standard types follow `512 < 51262 < 51263 < 51264 < 435663 < 51268`; explicitly selected or recolored cage IDs form the final highlight layers. The sorted layer sequence is independent of `show` argument order and color number. One layer uses radius 0.125 angstrom; multiple layers are distributed from 0.125 through 0.130 angstrom so coincident shared edges resolve in favor of the higher topology layer. Full cage, cluster, and domain IDs are reconstructed from compact GRO annotations; these IDs are assigned per frame and are not physical-object tracking identifiers.

The annotation reader stores one atom list per frame-local cage ID and derives cage-type selections from the ID prefix; it does not retain a second duplicate cage-type atom list. Phase, cluster, and domain atom lists are deduplicated after each frame. A lightweight trajectory-wide registry validates explicit cage/type/cluster/domain targets while still allowing an object that is absent from the current frame but present elsewhere in the loaded trajectory.

The renderer records VMD's stable representation names after each `mol addrep` and deletes only those names on the next redraw; the initial representation created with the SQQ molecule is adopted once, while later user-created representations are preserved. Frame traces keep at most one pending idle callback, so rapid animation updates coalesce to the final frame. Cage representations are ordered by fixed topology priority, with exact-object highlights last; changing a color or reversing `show` arguments does not change overlap ownership. Other object families retain color-override specificity ordering. Re-sourcing first removes the old frame trace, cancels its pending callback, and clears color/selection/representation state before loading the new bundle.

Per-frame output folders keep the configured grouped/flat structure for ordinary category files. Generated GRO paths and title lines use ASCII-only structure labels (`5^12`, `5^126^2`, `hc_5r_5^5`, and similar); Unicode superscripts/subscripts remain display-only in Markdown and main summaries.

When an output directory is reused, SQQ removes only known generated files that are outside the new effective selection. Temporary annotated-GRO fragments and partial visible bundles are removed after success or failure. Unknown user files remain untouched.

Mode `50` writes per-frame info, run-level annotated GRO/VMD, and `summary.xlsx`; ordinary category GRO and detail CSV files are not default output. `summary-xlsx` owns the workbook, `summary-csv` owns one file per main table under `output.summary_csv_dir`, and `summary-detail-csv` owns ordinary multi-row and isomer CSV files under `output.summary_detail_dir`. The directory settings default to `summary_csv` and `summary_detail` and must be different relative paths inside the output root.

`cluster-detail` separately owns `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. Explicit `cluster-detail` and `cluster-gro` require resolved cluster search on. Cluster search populates selected info and main-summary outputs without adding an unselected output type. Mode `00` already includes `cluster-gro`; mode `50` requires it explicitly. Search off removes stale SQQ-generated grouped directories or flat category filenames. If no per-frame type remains selected, the pipeline removes an empty frame directory. Unrelated files are preserved.

Per-frame output folders use the default grouped structure:

```text
frame_name/
  frame_name_info.md
  frame_name_order_parameter.tsv  # with order-tsv and selected F3/F4/Q_l
  ring/
  half_cage/<type>/
  quasi_cage/<type>/
  cage/<type>/
  ice/
  hydrate_cluster/                    # when cluster-gro is selected
    frame_name_cluster_sI.gro
    frame_name_cluster_sII.gro
    frame_name_cluster_sH.gro
    frame_name_cluster_boundary.gro
```

`summary-xlsx` writes the global `summary.xlsx` workbook. `summary-csv` writes the same applicable main-table mapping as separate UTF-8-SIG files under the configured `summary_csv_dir`, using sheet names as filenames and preserving columns, row order, and values without workbook formatting or tabs. The first `summary` table is a dashboard: Configuration begins with `SQQ version` and `Mode`, then includes input format/trajectory stride and applicable LAMMPS provenance, the same requested/effective `Graph mode` display used by the final terminal run summary, normalized order parameters, resolved `Find cluster`, and normalized output types; `Analysis Results (min / mean / max)` reports per-frame min/mean/max values while `Frames total / ok / failed` remains a run-level frame count. Analysis tables contain one input file or trajectory frame per row. `failures` instead has one failed input/frame per row, `detail_index` has one generated detail file per row, and `config` stores configuration metadata. Exact quasi-cage isomer rows, optional `failures.csv`, and other ordinary multi-row detail tables are written under the configured `summary_detail_dir` only when `summary-detail-csv` is selected. Failure details are always retained in `run_config.yaml`.

The SQQ-CPP main summary deliberately contains only applicable tables: `summary`, `cage`, `cage_isomer`, selected F3/F4 `order_parameter`, and `config`; `failures` is conditional and `cage_occupancy` exists only when at least one frame contains selected guests. Its default `summary-csv` writes these as independent files, while explicit `summary-xlsx` writes the same mapping as workbook sheets. Ring/connection diagnostic tables, half/quasi, cluster, ice, detail index, and detail CSV are omitted from its compact schema. If no selected guests exist, the dashboard and per-frame info state that occupancy was not evaluated.

Each per-frame `*_info.md` report is optimized for inspection rather than plotting:

- Frame Information begins `sqq version`, `mode`, `date & time`, `source`, input provenance, `frame`, and `time_ps`; modes display as `00 (sqq-py)`, `50 (sqq-py)`, `99 (sqq-cpp)`, or `sqq-cpp`;
- the Ring table shows only report-selected ring sizes, with `total` primitive-ring counts, final `free` ring counts, and a sum row for both columns;
- Half Cage and Quasi Cage omit internal `hc_`/`qc_` prefixes, aggregate each composition on a parent row, and list exact isomers on synchronized child rows;
- Cage combines composition totals and structural isomers in one vertical table: each cage composition is a parent row, and observed isomers are synchronized child rows below it;
- Quasi Cage Isomer Description gives one explanation row for each observed quasi-cage isomer, including the base ring and L1/L2/L3 ring sequence;
- Cage Isomer Description gives one explanation row for each observed cage isomer, including face composition and the 6-ring face adjacency pattern;
- Cage Occupancy remains separate because it describes guest assignment, not topology; it uses one cage type per row and dynamic exact guest-composition columns in source guest order;
- when resolved cluster search is on, one compact `Hydrate Cluster` hierarchy follows the cage sections; exact IDs remain in optional `cluster-detail` CSV files, and native sI/sII/sH/boundary structures remain in selected `cluster-gro`;
- hierarchy labels use short tree symbols, and Markdown source tables are padded using Unicode display width so their pipe columns align;
- explicit overrides change later effective-setting rows without hiding the selected engine; SQQ-Py also records resolved `find_cluster`, while SQQ-CPP omits that inapplicable row;
- SQQ-CPP retains Molecules, active connection coordination, Cage, Cage Isomer Description, occupancy status, selected F3/F4, and Warnings, but omits Ring, Half Cage, Quasi Cage, Hydrate Cluster, Hydrate Nucleation, and Ice sections.

The optional `summary-xlsx` workbook and `summary-csv` directory keep plotting-oriented frame-analysis tables with one input file or trajectory frame per row; dashboard, failure, detail-index, and configuration tables use the metadata-specific row units described above. Their `hydrate_cluster` table reports `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`, which are mutually exclusive within each reported cluster. Their `quasi_cage` table aggregates exact isomer columns into composition-level columns such as `5r_5²6³`. When `summary-detail-csv` is selected, failed rows are written to `<summary_detail_dir>/failures.csv`, exact quasi-cage isomers are split into `quasi_cage_isomer.csv`, and other ordinary multi-row detail outputs are `cage_occupancy.csv` and `cage_isomer.csv`. When cluster search is enabled and `cluster-detail` is selected, `hydrate_domain.csv` and `hydrate_cluster_detail.csv` are also written. Cluster detail includes the four category ID groups and `boundary_composition`; domain detail uses `external_boundary_contact_count` and `external_boundary_contact_ids` for direct external adjacency. `detail_index` records only generated detail CSV files. `cage_isomer.csv` defaults to nonzero isomer rows plus per-frame totals; `output.cage_isomer_rows = all` or `--cage-isomer-rows all` restores the full zero-filled matrix. Public motif output is not written. The `order_parameter` table reports only selected F3/F4 mean/count pairs, selected `qN` mean/count pairs, and selected MCG/DHOP largest clusters. Matching focus-water mean/count columns are added only when `order.focus_waters` is non-empty. MCG without a configured guest remains `N/A`. `order.parameters: []` or `--order-parameter none` omits the table. Output type `order-tsv` contains only selected per-water F3/F4/Q_l columns; MCG/DHOP remain frame-level outputs, so a hydrate-only selection does not create an otherwise empty `*_order_parameter.tsv`.
Summary generation records all table dimensions and write/format/save timings in `run.summary_write`. Exact quasi isomers are carried as sparse records until `quasi_cage_isomer.csv` is built, avoiding a summary DataFrame column for every observed isomer; composition-level quasi counts remain in the compact tables. Main CSV, `summary.xlsx`, every detail CSV, and `run_config.yaml` are written through same-directory temporary files and atomically replaced after success. Detail CSV replacements/removals commit as one recoverable bundle. Stale CSV cleanup is confined to known SQQ-generated filenames inside the currently configured `summary_csv_dir` and `summary_detail_dir`; unknown files are preserved, and a previous custom directory is not scanned after the setting changes. Before any table is handed to pandas Excel output, SQQ checks the 1,048,576-row and 16,384-column workbook limits and reports an actionable error for an unexpected oversize compact sheet. For an analysis sheet above 200,000 cells or 128 columns, Excel keeps header style, filter, freeze pane, and fixed widths but skips per-body-cell formatting and broad auto-width scans. This I/O policy does not change values, row/column schemas, or CSV selection.

## Current Limits

- The implemented PBC path remains orthorhombic. Non-orthogonal/triclinic GRO boxes and trajectory cell angles are detected and rejected; conversion must occur before SQQ analysis.
- SQQ-CPP is intentionally limited to graph, internal chordless 4/5/6 rings, cage topology/isomer/occupancy, and F3/F4. SQQ-Py modes are required for public rings, half/quasi, cluster, ice, Q_l/MCG/DHOP, legacy per-frame VMD, TSV, or `summary-detail-csv` output. The run-level annotated GRO/VMD renderer is shared by both engines.
- XYZ input has configurable coordinate scaling but no periodic box metadata.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- Hydrate domains are per-frame topological regions; temporal grain tracking and crystallographic orientation matching are not implemented.
- Hydrate phase classification depends on the detected cage/search scope but not on `--cage-size`; changing ring search sizes, cage face limits, or detection thresholds can still change the available topology.
- Boundary membership is a per-frame first-layer topological classification; transition-path kinetics, temporal domain tracking, and crystallographic orientation matching are not implemented.
- Cluster GRO output preserves original wrapped coordinates and the original box rather than reimaging each category; periodic or percolating structures may therefore retain unavoidable cross-box seams in a single-copy file.
- Default `quasi_cage.search_policy = bounded` is not exhaustive for large outer-layer components. Opt-in `exact` enumerates connected subsets but remains subject to explicit candidate and state limits.
- Automatic process workers parallelize independent GRO/XYZ files or selected frames of one indexed XTC/TRR/LAMMPS trajectory. Worker counts are based on physical cores with one physical core reserved for the system; topology search inside one individual frame remains single-process.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
- MCG is meaningful only for guest residue names selected in `hydrate_order.mcg_guest_resnames`; other guest species are not silently treated as methane.
- Published DHOP transition-state thresholds are model- and condition-dependent; SQQ reports the descriptor and does not assign a universal critical-nucleus threshold.
