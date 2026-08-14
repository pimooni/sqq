# SQQ Development Design

SQQ means **Shell Quant Qualifier**. This document records the current implementation logic for developers, so the code and the scientific definitions stay aligned.

## Pipeline

Engine values `00` and `py` run the complete SQQ-Py pipeline:

```text
input frames
  -> molecule selection
  -> water graph: hydrogen bond / O-O / user pair map
  -> diagnostic coordination distribution
  -> chordless rings (default) or optional shortest-path rings
  -> optional half_cage and/or quasi_cage open patches
  -> closed cage search and guest occupancy
  -> optional hydrate_cluster analysis from all detected cages
  -> F3/F4/Q_l plus MCG/DHOP order parameters and ice metrics
  -> per-frame outputs, selected main summary format, and optional detail CSV
```

Engine values `99` and `cpp` use the same Python shell around a reduced native frame pipeline:

```text
Python input/config/selection
  -> C++ graph
  -> C++ internal chordless 4/5/6 rings
  -> C++ cage topology, isomer, and automatic occupancy
  -> C++ F3/F4
  -> Python full-frame GRO/XTC/TSV/Tcl renderer, selected cage GRO, info, and summary CSV or optional XLSX
```

In SQQ-Py, the shared water graph is used by ring, half_cage, quasi_cage, cage, selected F3/F4, graph-mode Q_l, and ice analysis. Selected MCG and DHOP descriptors are calculated during the order stage but use dedicated guest/water cutoff graphs because their published definitions are independent of the selected SQQ bond mode. Hydrate_cluster analysis starts after cage detection and uses all detected cage-ring memberships, not the raw water graph or the report-filtered cage list. In both engines, the graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode. Coordination diagnostics read this graph without adding, removing, or capping edges.

### Stacked GRO trajectory input

A single GRO path may contain one frame or repeated complete GRO blocks. Repeated blocks are one trajectory, not a multiple-GRO batch: SQQ streams them in file order, parses physical time from each title `t=` / `time_ps=` token, and applies the same time-selection contract used by XTC/TRR/LAMMPS. Every block must preserve atom count and ordered `(resid, resname, atomname, atomid)` identity. Coordinates, optional velocities, orthorhombic boxes, and text after a semicolon are parsed frame by frame. A GRO used as `--top` must contain exactly one frame.

### Multiple-GRO topology grouping

Grouping is a pre-analysis boundary used only when one invocation contains two or more GRO files; it does not own the later output-layout decision. A single one-frame GRO stays compact, while a stacked GRO and other trajectory-like inputs use the separated layout described below. The pre-scan parses each GRO atom block and builds a deterministic fingerprint from the total atom count plus ordered contiguous residue blocks. Each block contributes its `resname` and ordered `atomname` sequence. GRO title/time text, coordinates, optional velocities, box values, and numeric atom/residue identifiers are deliberately excluded, so renumbered or coordinate-shifted frames of the same molecular topology stay together. Atom order, residue-block order, names, or atom count changes create a different topology group. Groups are assigned `A` through `Z` in first-occurrence order.

One detected topology uses the requested output root directly. Two through 26 topologies use independent `result_A` through `result_Z` roots, each with its own summary, info/GRO trees, four-file `sqq_render/` package, and group `sqq_config_resolved.yaml`; the requested root keeps a batch `sqq_config_resolved.yaml` source-to-group manifest. More than 26 topologies cannot be represented by the fixed letter namespace. In that case the entire invocation switches to information-only output at the requested root: every readable GRO is still analyzed, but all summary, detail, ordinary GRO, and render-package outputs are disabled. The limit is applied before creating partial `result_A`-`result_Z` trees and overrides explicit output selections.

A single GRO topology supplied with `--top` is validated against every GRO input during the same pre-scan. A mismatch is fatal before analysis and reports the exact incompatible source; one topology cannot silently normalize heterogeneous GRO files. For requested `auto` graph mode, SQQ selects `hbond` or `oo` once from a representative frame in each topology group, stores that effective value beside the preserved requested value, and sends the same group configuration to either engine. Thus one group cannot mix graph definitions across its frames.

## Analysis Engines and Workers

The public selector is `-e` / `--engine`. Accepted values are discrete presets rather than a continuous 00-99 scale. Values `00` and `py` select SQQ-Py; values `99` and `cpp` select SQQ-CPP. The command default is `py`; removed value `50` has no alias.

| Engine value | Effective backend | Graph | Search sizes | Auto worker policy | Find cluster | Default output types |
| --- | --- | --- | --- | --- | --- | --- |
| `00` | `sqq-py` | `hbond` | 4/5/6 | 100% | on | `info,sqq-render,summary-xlsx` |
| `py` | `sqq-py` | `auto` | 4/5/6 | 1 worker | off | `info,sqq-render,summary-xlsx` |
| `99` | `sqq-cpp` | `hbond` | internal 4/5/6 | 100% | unsupported | `info,sqq-render,summary-csv,summary-detail-csv` |
| `cpp` | `sqq-cpp` | `auto` | internal 4/5/6 | 1 worker | unsupported | `info,sqq-render,summary-csv,summary-detail-csv` |

Preset/configuration application order is:

```text
built-in defaults -> engine preset -> sqq_config.yaml -> retained CLI overrides
```

Each preset owns graph mode, ring sizes, automatic worker policy, initial cluster state, and default output types. YAML and retained CLI values override compatible settings. `quasi_cage.max_layer` and `order_parameter.enabled` remain independent. Engine value `00` enables cluster search; `py` does not. Cluster search populates selected info/main-summary outputs but never adds an unselected output type. No preset includes ordinary, classified, or cluster GRO. `sqq-render` owns one indivisible GRO/XTC/membership-TSV/Tcl package.

Engine values `99` and `cpp` share the SQQ-CPP feature boundary. Both use internal chordless 4/5/6 rings, `f3,f4`, no cluster/quasi/ice analysis, and no silent Python fallback. Value `99` chooses the hbond/100% compatibility preset; `cpp` chooses auto graph with one worker. Neither native preset selects `gro` or `cage-gro` by default. Compatible YAML/CLI overrides remain available; incompatible requests fail before analysis.

Cluster-search precedence is:

```text
--find-cluster > hydrate_cluster.enabled in sqq_config.yaml > engine preset
```

An explicit `-b` / `--bond-mode {auto,hbond,oo,pairs}` overrides YAML `graph.mode` and the preset. `--pair PAIRS.txt` implies pairs mode unless `-b pairs` is already given; it cannot be combined with another explicit bond mode. CLI pair paths resolve from the working directory, while YAML `graph.pair_file` resolves from the user-config directory. `graph.pair_id` owns the identifier convention.

`parallel.worker: auto` follows the selected engine value: `py` and `cpp` resolve to one worker, while `00` and `99` calculate from 100% of detected physical cores, reserve one physical core for the operating system, and cap the result by the number of independent files or selected trajectory frames. Physical-core detection prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker N` / `-w N` overrides the preset using form-based parsing: integer text such as `1`, `4`, or `100` is an explicit worker count; decimal text such as `0.5` or `1.0` and percentages such as `50%` or `100%` are physical-core fractions. Percentages above `100%` and decimal values above `1.0` are rejected. Worker resolution remains capped by task count and the Windows `ProcessPoolExecutor` limit. YAML `parallel.backend` selects `process`, compatibility `thread`, or `serial`.

The public Analyze surface is exactly `-i/--input`, `-t/--top`, `-c/--config`, `-o/--output`, `-e/--engine`, `-w/--worker`, `-dt/--delta-time`, `-b/--bond-mode`, `-s/--size`, `--find-half`, `--find-quasi`, `--find-cluster`, `--order-parameter`, `--pair`, `--output-type`, and `-h/--help`.

The parser reserves migration diagnostics before ordinary unknown-option handling. `--mode VALUE` errors with `--mode has been replaced by --engine` and a `Use: --engine VALUE` line; `-m VALUE` gives the corresponding `-e` prompt. `--pairs FILE` errors with `--pairs has been replaced by --pair` and a `Use: --pair FILE` line. These paths exit with status 2 and never translate the invocation silently.

Former advanced CLI controls live only in canonical YAML: input discovery/strictness/scaling, LAMMPS interpretation, ring reporting/definition, quasi sizes/layers/policy, cage reporting/limits/validation, cluster minimum, Q_l neighbors, pair ID, parallel backend, and output layout/isomer rows. Legacy compatibility spellings `--workers`, `--no-q`, `-q`, `--q-degree`, `--mcg3`, `--dhop30`, and `--topology` are removed; the retained replacements are `--worker`, `--order-parameter`, and `--top` where applicable.

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
  rebuild SQQ models -> compact GRO/XTC/TSV/Tcl package / selected cage GRO / info.md / summary CSV or optional XLSX / sqq_config_resolved.yaml
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
Selected native cage GRO files preserve standard atom records and the source box. They contain cage waters and assigned guests without adding the synthetic `CNT` cage-center pseudoatom used by SQQ-Py. Separately, `sqq-render` defaults to the complete normalized input-frame topology and coordinates, including water hydrogens, complete guests, additives, environment/wall components, and other retained atoms. YAML `render.atom_scope: compact` selects the legacy water-oxygen plus complete-guest topology. This Python-owned writer is shared by SQQ-Py and SQQ-CPP and does not change native scientific calculations.

F3 and F4 are independently selectable. F3 uses the active graph. F4 requires usable hydrogen geometry; a frame without it reports no valid F4 values and a warning. `all` expands only to the native supported pair `f3,f4`, and `none` skips the order result.

#### CLI and validation boundary

The shared public CLI remains the compact Analyze surface listed above. Under SQQ-CPP, retained overrides are additionally constrained as follows:

- `-i`, `-t` / `--top`, `-o`, and `-c` retain the shared input/provenance contract;
- `-b auto|hbond|oo|pairs` and `--pair` feed the supported native graph definitions;
- `-s` is restricted to a nonempty subset of 4/5/6;
- `--order-parameter` is restricted to `f3`, `f4`, `f3,f4`, `all`, or `none`;
- `-w` / `--worker` retains the shared worker parser;
- `--find-half`, `--find-quasi`, and `--find-cluster on` are rejected because those searches are not native features.

Advanced native-compatible settings come from YAML: `graph.pair_id`; cage report/face/validation settings; `parallel.backend` restricted to `process|serial`; `output.structure_layout`; `output.cage_isomer_row`; and `output.type` restricted to `info`, `gro`, `cage-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, `summary-detail-csv`, `default`, `all`, or `none`.

The C++ output normalizer keeps `gro` as the ordinary classified-GRO umbrella. `sqq-render` owns the complete visualization package and `all` maps to every supported native output. Engine values `99` and `cpp` both default to `info,sqq-render,summary-csv,summary-detail-csv`; ordinary/classified GRO remains explicit. Native detail output contains cage occupancy and cage isomer tables. Public ring output, size 7, shortest-path rings, cluster, ice, Q_l/MCG/DHOP, legacy per-frame VMD or TSV output, fast closure, thread backend, and triclinic input are unsupported. Explicit retained-CLI requests outside the engine-neutral half/quasi compatibility rule fail before frames are analyzed. Generic `auto` half/quasi settings resolve off under C++; legacy YAML that explicitly enables those Python-only searches is downgraded with one warning, incompatible half/quasi output types are removed, and the adjustment is recorded in `sqq_config_resolved.yaml`. There is no `-m c` alias and no Python fallback if `_sqq_cpp` cannot be imported or returns an error.

#### Build and wheel architecture

`pyproject.toml` uses scikit-build-core. CMake builds the pybind11 module `_sqq_cpp` from `sqq/core/sqq_cpp/native/src/` with a C++17 compiler and installs it inside `sqq/core/sqq_cpp/`. The handwritten source and CMake metadata are included in the source distribution but excluded from binary wheels. Native build directories, CMake/Ninja state, extension binaries, wheels, source archives, and local test output are ignored; source `.cpp`, `.hpp`, and `CMakeLists.txt` files remain version-controlled.

Release CI uses cibuildwheel to compile and import-test CPython 3.10-3.14 wheels for Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, and separately builds a source distribution. A user installing a matching wheel receives an already compiled `.pyd` or `.so`; compilation occurs only for an explicit source build. The source path requires CMake 3.20 or newer, Python development headers, a C++17 compiler, and an appropriate platform build tool. Wheel and source-distribution publication are handled separately by the release workflow.

#### Scientific parity contract

SQQ-Py is the regression oracle. Deterministic discrete results must match exactly for the same effective settings: graph edges, canonical ring nodes, cage type/water/face membership, cage-isomer label, and occupancy guest assignment. F3/F4 comparisons permit a small floating-point tolerance because the C++ and Python accumulation implementations are independent.

The native parity baseline established in 0.3.1 includes random graph/ring frames, random F3/F4 frames, synthetic cage/geometry/occupancy cases, package import tests, and a real 11,104-atom frame. On `tests/100.gro`, both engines produced 4,322 edges; 2,499 internal rings (45/2,147/307 for sizes 4/5/6); 339 cages across 16 types with exact water and face membership; 339/339 matching isomers and occupancy assignments; and 315 occupied cages. F3 matched exactly for 2,176 waters, while F4 had maximum absolute difference `4.44e-16`. The native core took about 0.5673 s versus 10.0967 s for the equivalent Python path on that host. This 17.8-fold core-path ratio is a benchmark, not an end-to-end guarantee.

### Process Execution Architecture

For two or more independent GRO/XYZ files, the default execution path is:

```text
main process
  -> normalize config and resolve input order
  -> for multiple GRO files, pre-scan and assign topology groups
  -> validate shared GRO topology and group-local output names
  -> resolve requested auto graph mode once per GRO topology group
  -> create one spawn ProcessPoolExecutor for the complete invocation
  -> initialize each worker once with group configs/output roots/stage queue
  -> maintain a rolling queue of at most 3 * workers tasks
  -> submit (global_index, group_key, group_local_index, path)

worker process
  -> report progress with the global input index
  -> read one frame
  -> select the matching group config and output root
  -> run the ordinary frame pipeline
  -> write group-local info and selected structure output
  -> return (global_index, group_key, group_local_index, summary_row)

main process
  -> consume stage events across every group
  -> reorder rows globally and within each group
  -> finalize each selected group bundle and summary independently
  -> finalize group run configs and the root batch manifest
```

Topology groups are task-routing metadata, not separate serial jobs. One shared worker pool can execute tasks from different groups concurrently. The global index keeps progress and failure reporting stable, while the group-local index controls summary rows and compact render-frame order inside one compatible system.

`spawn` is selected explicitly on macOS, Windows, and Linux. This avoids forking a process after the interactive progress refresh thread exists and gives the same pickling/import contract on every platform. Worker callables and initializers are module-level functions. Only paths, raw trajectory indexes, small event tuples, and summary dictionaries cross process boundaries; atoms, rings, patches, and cages stay worker-local.

While the pool exists, the parent sets the common BLAS/OpenMP thread environment variables to `parallel.math_thread` (default 1). Spawned children inherit the limits before importing NumPy-backed modules. The parent environment is restored after pool shutdown.

For one XTC/TRR or supported LAMMPS trajectory with `-t` / `--top`, the parent reads physical frame times once and resolves `input.delta_time_ps`-selected raw indexes. The requested interval must be at least and an integer multiple of a regular native interval. Each worker initializer opens a private MDAnalysis Universe once, caches immutable atom metadata once, and tasks contain a small contiguous batch of `(ordered_frame_index, raw_frame_index)` pairs. Batch size is `ceil(selected_frames / (4 * workers))`, clamped to 1 through 8. Parent and worker readers are explicitly closed. A stacked GRO is streamed serially; multiple trajectory files and non-process trajectory backends remain serial. SQQ retains an orthorhombic box representation; non-90-degree trajectory angles are detected and rejected instead of being silently approximated.

Both standalone-file and indexed-trajectory process paths, plus the compatibility thread path, maintain at most `3 * workers` submitted tasks and refill the queue as futures complete. This is a bounded submission window, not a worker cap: 100 effective workers retain 100-way execution and at most 300 submitted tasks. It avoids constructing a Future and serializing arguments for every item in a very large input set.

Standalone files whose case-insensitive stems collide inside the same output root are rejected because the stem supplies the per-frame output name. Non-strict serial, thread, and process read failures return failed summary rows; strict failures cancel queued work where possible and propagate to the main process. Indexed trajectory reader failure records the current failed frame and stops that unusable iterator. Output order is determined by the original index rather than completion order. Configuration is normalized before dispatch, so thread workers only read the shared mapping.

Terminal and main-summary dashboard metadata share the same display helpers. The requested graph mode is preserved from YAML or retained CLI. Explicit graph modes display as `hbond`, `oo`, or `pairs`. Automatic graph mode is resolved during preflight from the first selected frame before the terminal header or initial resolved YAML is written. A compatible topology therefore displays only `auto -> hbond` or `auto -> oo` in the terminal, per-frame reports, summaries, VMD metadata, and resolved YAML. Multiple-GRO topology groups are resolved independently; when groups differ, the root terminal and manifest list one exact mode per `result_A` ... `result_Z` instead of a pending or mixed placeholder.

Root `sqq` / `sqq -h` output renders the banner and the `SQQ (Shell Quant Qualifier)` product sentence, then `SQQ version: 0.5.2   Release date: Aug 13, 2026`, then the ordinary `usage:` line. The root command list contains `init`, `analyze`, `track`, and `vmd`. Root `sqq -v` / `sqq --version` exits successfully after printing only that version line. Subcommand help retains the standard argparse layout. Every `sqq analyze` and raw-input `sqq track` invocation prints the SQQ banner before directory creation, output locking, reader initialization, or other preflight work; the version remains in the standard `Configuration` block, so it is not duplicated. Even an early validation failure therefore has a stable first screen.

`sqq init` writes the fixed commented template to `sqq_config.yaml` by default; `sqq init -o NAME.yaml` changes only that destination. It refuses to overwrite an existing file. The template uses `#` section/choices comments and defaults `ring.size` to `[4, 5, 6]`. Analyze does not auto-read a same-named current-directory file when `-c` is omitted. User YAML is never rewritten: relative paths in it resolve from the YAML directory, while CLI paths resolve from the invocation directory. Duplicate mapping keys and unknown public keys are errors.

Terminal headers, per-frame info, and summary configuration use the normalized label `SQQ engine: sqq-py` or `SQQ engine: sqq-cpp`. They do not render requested values as `py (sqq-py)` or similar compound mode labels.

The mandatory output-root `sqq_config_resolved.yaml` is the authoritative runtime record. It preserves normalized analysis settings and adds final effective metadata: SQQ version, requested `engine`, effective `sqq-py`/`sqq-cpp` backend, input/topology provenance, requested and effective graph modes, requested worker policy and resolved workers, backend/math threads, normalized output types, automatic adjustments, status/error, frame totals, failures, and `summary_write` timing/table dimensions. It is initialized with `status: running` before frame analysis and atomically replaced with `completed` or `failed`; a failed rewrite does not truncate the previous complete file. The detailed `config` worksheet in `summary.xlsx` and `summary/config.csv` are no longer built. Main summaries retain only the compact dashboard Configuration block. `--output-type` and YAML `output.type` are the two output selectors; removed `output.disabled_outputs` configurations are rejected rather than migrated.

Configuration is normalized once before execution. Public YAML uses the canonical top-level `engine` key and singular collection keys (`resname`, `size`, `enabled`, `type`, `worker`, and related forms); order selection lives at `order_parameter.enabled`, and pair settings live at `graph.pair_file` / `graph.pair_id`. The runtime may keep established internal names. Legacy top-level `mode`, `graph.bond_mode`, and `order.parameter` migrate with warnings to `engine`, `graph.mode`, and `order_parameter.enabled`; 0.3.x plural YAML is also migrated on read. Textual booleans use explicit on/off parsing rather than Python truthiness; enum values are canonicalized and rejected when unsupported; cutoffs and scales must be finite and positive where required; counts must be positive integers, while cage state guards accept nonnegative integers with `0` meaning unlimited; physical sampling intervals must be finite and positive. Residue-name lists accept comma-separated text or sequences and are normalized deterministically.

Generic `sqq_config.yaml` files set `half_cage.enabled` and `quasi_cage.enabled` to `auto`. Engine selection occurs before capability normalization: auto resolves on for SQQ-Py and off for SQQ-CPP. With C++, legacy explicit `on` values, quasi layer controls, and half/quasi output requests are converted to inactive effective values with concise warnings and `sqq_config_resolved.yaml` adjustment records. The hard-error boundary is reserved for missing required input or a setting that prevents native cage analysis.

### Terminal Progress Display

Analyze owns one deterministic stdout interface. The banner is followed by `Basic Information`, one atomic `Configuration` block containing the SQQ version, `Analysis Progress`, one final-output status line, one optional `Diagnostics` block, and `Run Summary`. After frame analysis closes its 100% panel, SQQ prints `Writing output files; please wait and do not close SQQ...` exactly once before render, Excel/CSV, and final resolved-config publication. A multiple-GRO invocation uses the same single line for all groups. Ordinary progress and diagnostics are written to stdout. Fatal exception text remains an error-path concern rather than part of the normal dashboard.

SQQ-Py serial and parallel runs share the same three-row stage model:

```text
file preparation       reading -> settings -> selecting
core topology search   graph -> ring -> half/quasi -> cage -> cluster
post-processing        filtering -> order -> ice -> output
```

`cluster` is included only when the resolved `hydrate_cluster.enabled` value is true, for example through engine value `00` or `--find-cluster on`. When hydrate cluster analysis is disabled, the stage is omitted rather than shown as `cluster:0`. Likewise, `half/quasi` is omitted when both `half_cage.enabled` and `quasi_cage.enabled` are false.

SQQ-CPP keeps the same aligned three-row presentation but removes Python-only work:

```text
file preparation       reading -> settings -> selecting
native topology        graph -> ring -> cage
post-processing        order -> output
```

The public ring report is absent even though the native `ring` stage is required internally for cage construction.

On a TTY, serial and parallel execution use one in-place progress panel. The panel renders the complete workflow and highlights the active stage with ANSI bold plus bright blue (`RGB(0,0,255)`). Stage columns are sized by the longest stage name in that column, so the display stays compact while `reading`, `graph`, and `filtering` remain aligned. The continuation marker `>` for a new row is placed before the aligned stage column. The timing row remains `stage / frame / total`. Refreshes compare the complete rendered content and skip identical updates; close is idempotent and leaves exactly one final completed panel.

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

When stdout is not a TTY, SQQ writes bounded static checkpoints at approximately 5% increments, including the initial and final state, for at most about 21 progress records. This path emits no ANSI controls, carriage-return redraws, or `tqdm` output, so redirected logs remain stable and line-oriented.

Parallel GRO/XYZ and indexed trajectory runs use a main-process progress aggregator. Spawned workers never write terminal control sequences; they send `start` and stage-transition tuples through a multiprocessing queue. The main process applies those events, ignores late events from already-finished tasks, and shows completed/failed/active/queued counts, compact `stage:count` rows, total elapsed time, and up to six active-file rows. Additional active files are summarized so high-worker modes do not fill the terminal.

Warnings captured during a successful run are normalized, deduplicated by exact rendered text, and emitted once in a single `Diagnostics` block after progress has closed and before `Run Summary`. The block is omitted when empty. A fatal path may flush the same deduplicated block immediately before re-raising, but it cannot create repeated diagnostics sections during normal execution.

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

### Final Results Page

Progress measures frame science only. The workflow records an analysis-complete timestamp immediately before report, summary, render-bundle, tracking-state, and resolved-configuration publication. After every selected output is safely written, an interactive TTY is cleared and redrawn with the retained `Basic Information` and `Configuration` sections, followed by `Analysis Results` or `Tracking Results`. Redirected/captured streams are not cleared and remain plain append-only text.

The result block reports requested, analyzed, successful, and failed frames; total, analysis, and output-writing seconds; mean seconds per successful frame; final status; and absolute result path. Multi-GRO groups use the same final renderer and report their topology-group count rather than printing raw internal mappings. Track source mode displays imported-state and matching configuration only; it does not present unused Analyze defaults.

Citation wording is derived from an explicit map of features that actually ran and output types that actually completed. A configured feature that never executed is omitted. A completed search may be described as an analysis even when it found zero objects, but the text does not claim an object was identified unless results support that claim. Supported TTYs render the recommendation sentence in bold. The block always ends with the provisional publication statement and `https://github.com/pimooni/sqq`.

## Modules

SQQ 0.5.2 uses these responsibility boundaries:

```text
sqq/
├ config/                 defaults, migration, resolution, and validation
├ models/                 structure, topology, result, and tracking contracts
├ core/
│  ├ sqq_py/              Python frame-analysis backend
│  ├ sqq_cpp/             Python adapter and native C++17 source in native/
│  ├ order/               F3/F4, Steinhardt Q, MCG, and DHOP
│  └ tracking.py          deterministic temporal matching algorithms
├ workflow/               init, analyze, track, and vmd command flows
├ runtime/                plans, frame tasks, workers, and execution policy
├ io/
│  ├ reporting/           per-frame and aggregate reports
│  ├ render/              Analyze/Track render packages and Python Tcl template
│  ├ tracking.py          Track state, CSV, and membership conversion
│  └ trajectory/LAMMPS/GRO readers and writers
└ ui/                     progress, diagnostics, final results, and citations
```

The two cage kernels are peers. `sqq.core.sqq_py.analyze_frame` and `sqq.core.sqq_cpp.analyze_frame` share the normalized `Frame`, selected `Water`/`Guest`, configuration, callback, and `FrameResult` contract. Their documented capabilities remain different; the shared Python workflow applies engine normalization before dispatch and never falls back silently from C++ to Python.

`workflow/init.py`, `workflow/analyze.py`, `workflow/track.py`, and `workflow/vmd.py` are parallel command entry points. Workflow code owns orchestration and final publication; scientific acceptance remains in `core`, serialization remains in `io`, scheduling remains in `runtime`, and terminal-only presentation remains in `ui`. The C++17 implementation and pybind11 bindings live below `core/sqq_cpp/native`, so native and Python backends are visibly parallel while wheel builds still compile one `_sqq_cpp` extension.

The shared VMD renderer is stored as a readable raw string in `io/render/tcl_template.py`; no standalone Tcl source file is tracked or packaged. Python supplies the generated manifest, actual filenames, molecule name, and synchronized help body when publishing Analyze or Track scripts. The published render package still contains the ordinary `*.vmd.tcl` file required by VMD, with the same content and command interface as before. `workflow/vmd.py` owns the public `sqq vmd` flow; render discovery, validation, and file generation remain I/O responsibilities.

The public Python surface is `sqq.load_config`, `sqq.read_frames`, and `sqq.analyze_frame`, together with the exported models and typed exceptions. Retired monolithic modules are removed rather than kept as parallel compatibility implementations. Internal callers use the responsibility-specific package paths shown above.

Important implementation modules remain:

- `core/graph.py`, `ring.py`, and `ring_topology.py`: water graph, ring search, and sparse ring incidence.
- `core/half_quasi.py`, `cage.py`, and `phase.py`: open patches, exact closed cages, and frame-local hydrate phase/domain classification.
- `core/order/{f3f4,steinhardt,mcg,dhop}.py` and `core/ice.py`: order parameters and ice classification.
- `core/tracking.py`: persistent-ID matching, events, lifetimes, targets, and streaming accumulation.
- `io/trajectory.py` and `io/lammps.py`: normalized coordinate/topology readers.
- `io/tracking.py`: versioned Track state, normalized tables, target selection, and persistent membership conversion.
- `io/render/tracking.py`: target render discovery, validation, naming, and publication.
- `ui/final_results.py` and `ui/run_statistics.py`: final-page wording and measured run statistics.

### Stable Python API

`sqq.load_config` accepts a YAML path, a partial mapping, or `None`, performs migration and engine-capability normalization, and returns a `ResolvedConfig` carrying an auditable `ResolutionReport`. `sqq.read_frames` yields normalized public `Frame` values through the same supported readers used by the CLI. `sqq.analyze_frame` dispatches one `Frame` to the resolved Python or C++ backend and returns the public `FrameResult`. Configuration, input, and analysis failures use distinct exported SQQ exception types. These functions are the supported programmatic interface; workflow, runtime, and I/O modules remain implementation details.

## Input Validation and Coordinate Units

GRO declares coordinates in nm. Each frame requires a title, declared atom count, all corresponding finite atom records, and one separate 3- or 9-value box line. A file may contain one frame or repeated complete blocks; stacked frames must keep identical ordered atom identity. Atom and residue names must be non-empty. Optional velocities are retained, while text after `;` is treated as an annotation and ignored by numeric parsing. An all-zero box is normalized to non-periodic `None`; mixed nonpositive lengths are invalid. A nine-value box is accepted only when all six tilt terms are zero. Nonzero tilt and non-finite values fail fast.

XTC/TRR positions and lengths supplied by MDAnalysis are converted from angstrom to nm. Coordinates and cell angles must be finite, and angles must be 90 degrees within tolerance; triclinic frames are rejected. Non-finite or unparseable trajectory time is retained as unavailable rather than emitted as a numeric time.

LAMMPS trajectories use one shared Python boundary adapter before either engine runs. A LAMMPS DATA topology supplied through `-t` / `--top` is mandatory. The adapter accepts `.dump` / `.lammpstrj` and LAMMPS DCD. A non-empty YAML `input.lammps.type_map` is authoritative. When the map is absent or empty, the adapter derives element candidates from DATA masses/type comments and validates molecule graphs from DATA Bonds. It recognizes one O plus two H joined by two O-H bonds as water, one C plus four H joined by four C-H bonds as all-atom methane, and a clearly carbon-labeled unbonded singleton as a united-atom guest. Other bonded components such as surfactants or walls are retained with deterministic environment/other labels rather than rejected or admitted to the water graph. `component.role_map`, `additive.resname`, and `environment.resname` can override those roles; water and guest recognition remains strict. Valid DATA molecule IDs are retained, or rebuilt from unambiguous bond-connected components when necessary. Every occurrence of one numeric atom type must resolve consistently; ambiguous mixed water/guest roles, insufficient evidence for a requested scientific role, non-unique reconstruction, or ID mismatch fails clearly. The resolved map, component role, and provenance are stored in `sqq_config_resolved.yaml`, per-frame `*_info.md`, and main-summary configuration. DATA atom IDs become SQQ atom IDs, dump rows are restored to stable topology atom-ID order, and `real`, `metal`, or `nano` lengths and times are converted to nm/ps exactly once. The MDAnalysis LAMMPS reader receives an explicit physical frame increment, `dt = settings.timestep * settings.time_to_ps` ps, so its frame times and SQQ sampling metadata use the same unit conversion without emitting the reader-missing-`dt` warning. Every dump frame must declare `pp pp pp` and an orthorhombic cell; DCD cell angles must be 90 degrees. Tilt, mixed/nonperiodic boundaries, `units lj`, duplicate or changing atom-ID sets, and non-finite/nonpositive cells fail before topology analysis. The process scheduler gives each worker a private reader and immutable mapped atom metadata. SQQ-CPP receives the same normalized `Frame` as SQQ-Py and contains no separate LAMMPS parser.

GRO readers preserve optional fixed-width velocities. SQQ annotations after a semicolon are stripped before the ordinary record is parsed, allowing the generated visualization GRO to be read back by SQQ without changing its coordinate, identity, velocity, or box fields.

XYZ has no standard unit or box metadata. Coordinates are multiplied by YAML `input.xyz_scale`; the default `0.1` assumes angstrom input, and `1.0` preserves coordinates already expressed in nm. SQQ accepts exactly one nonnegative declared atom count, exactly that many finite coordinate records, and no extra nonempty records; multi-frame XYZ must be split into files before analysis. XYZ remains non-periodic unless converted to a format carrying a box.


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
order_parameter:
  enabled: [f3, f4]
```

The equivalent CLI is `--order-parameter f3,f4`. Supported names are `f3`, `f4`, any non-negative `qN` degree such as `q6` or `q12`, `mcg1`, `mcg3`, `dhop35`, and `dhop30`. `all` expands to `f3,f4,q6,q12,mcg1,mcg3,dhop35,dhop30`; `none` disables the entire order stage output. Names are deduplicated and stored in canonical order. An explicit CLI list replaces the complete configured list rather than extending it:

```text
--order-parameter > order_parameter.enabled > default [f3, f4]
```

The pre-0.2.7 selectors `--no-q`, `-q` / `--q-degree`, `--mcg3`, and `--dhop30` are removed rather than retained as hidden compatibility options. Descriptor selection is expressed only through retained `--order-parameter` or YAML `order_parameter.enabled`.

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
- `cutoff`: use all water oxygens within `order_parameter.q_cutoff_nm`.
- `nearest`: use the nearest `order_parameter.q_n_neighbor` water oxygens within `order_parameter.q_cutoff_nm`.
- `lammps`: LAMMPS-compatible cutoff plus fixed-neighbor behavior; if `order_parameter.q_n_neighbor` is null, it defaults to `12`.

The non-graph modes reuse the deterministic orthorhombic/non-periodic cell-list pair search from `sqq/core/spatial.py` rather than scanning every O-O pair. Distances and vectors are still recomputed with the shared float64 minimum-image function and sorted deterministically.

When a fixed neighbor count is active and fewer than that number of neighbors are found inside the cutoff, every requested Q_l value is set to `0.0`, matching LAMMPS behavior. Without a fixed neighbor count, waters with no Q_l neighbors are omitted from the Q_l mean and count.

Q_l is a continuous structural descriptor, not a standalone ice-count classifier. The value is sensitive to the neighbor definition, so SQQ records selected Q_l names and YAML `order_parameter.q_neighbor_mode`, `order_parameter.q_cutoff_nm`, and `order_parameter.q_n_neighbor` settings in `sqq_config_resolved.yaml`, the terminal header, and selected main-summary output. These settings are active only when at least one `qN` parameter is selected.

## Hydrate Nucleation Order Parameters

MCG/DHOP are opt-in frame-local descriptors and are separate from the cage-topology `hydrate_cluster` hierarchy. Selecting `mcg1`, `mcg3`, `dhop35`, or `dhop30` activates only those requested outputs in the existing `order` terminal stage. Their results are stored in `HydrateOrderResult`, where `None` means not applicable and integer zero means that the calculation was applicable but no qualifying component was found. Numerical thresholds remain under `hydrate_order`; descriptor selection belongs only to `order_parameter.enabled`.

### Shared spatial search

`sqq/core/spatial.py` supplies deterministic self- and cross-cutoff pairs. For an orthorhombic box, coordinates are wrapped into cells whose widths are at least the cutoff; only the 27 neighboring cells are inspected. Candidate distances are recomputed in float64 with the same minimum-image function used elsewhere. Pairs are deduplicated and sorted, so process/serial output and tie breaking are stable. Non-periodic input uses the same cell scheme without wrapping. No fixed atom or neighbor array is used.

### MCG-1 and MCG-3

The selected guest set is controlled by `hydrate_order.mcg_guest_resname` (default `CH4`, `MET`). A configured center atom is used when available; otherwise the residue is unwrapped around its first atom before its centroid is calculated. Original guest-list indices are retained for output membership.

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

Rings are searched on the already-built water graph; geometry is not used after the graph edge set is fixed. Sorted adjacency tuples are built once. Every bounded DFS state stores one short immutable path, and membership is checked directly against that path. The minimum node must remain the start, and a closing path is accepted only in one direction (`second < last`), eliminating rotational and reverse rediscovery before final canonical ordering.

With the default `ring.definition: chordless`, adding a node checks its graph neighbors already present in the partial path. A connection to any non-previous node is an immediate chord and prunes the branch. A connection back to the start closes the candidate immediately; the path is not extended beyond that closure edge. Random-graph regression compares this optimized traversal with the previous final-only chord test.

With `ring.definition: shortest_path`, every chordless candidate additionally applies the Franzblau shortest-path criterion. For each ring node, a bounded BFS is run only to `floor(size/2)`; the graph distance to every other ring node must equal the shorter distance along the cycle. A frame-local `(source, depth)` cache reuses these bounded distance maps across candidates. This opt-in definition can remove rings and therefore change patch/cage results.

Current behavior:

- supported ring sizes: 4, 5, 6, 7;
- default ring sizes: 4, 5, 6;
- default definition: `chordless`;
- optional definition: `shortest_path`;
- ring nodes are water oxygen indices;
- `ring.size` / `-s` / `--size` controls detection, while YAML `ring.report_size` filters ring tables and GRO files after detection.

The historical `ring.primitive` default was not connected to the search implementation and is no longer emitted. Explicit definition now uses YAML `ring.definition`.

## Half-Cage and Quasi-Cage Terms

`patch` means a connected set of ring faces during search. It is not an output class by itself.

Layer definitions:

```text
L0 = base ring
L1 = side rings sharing every base-ring edge; L1 must close into a full side wall
L2 = rings grown outward from L1
L3 = rings grown outward from L2
```

L2 and L3 may be dangling rings or connected dangling ring chains. They do not need to close. `half_cage.enabled` / `--find-half` and `quasi_cage.enabled` / `--find-quasi` independently control the two reported families. Generic `auto` resolves on in SQQ-Py and off in SQQ-CPP. The default quasi search reports L1; L2/L3 remain available through YAML `quasi_cage.max_layer`, which is inactive when quasi search is off.

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
8. If `quasi_cage.max_layer >= 2`, grow L2/L3 from the current frontier:
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

- `quasi_cage.max_layer`: default 1 for fast routine analysis; use 2 or 3 to report outer dangling quasi-cage layers.
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
cage.max_state_per_seed = 0
cage.max_total_state = 0
cage.max_boundary_candidate = 8
cage.scientific_validation = false
```

Ring seeds remain the default; patch seeds are retained for compatibility and targeted comparisons. `0` state guards mean unlimited exact search. `max_boundary_candidate` does not truncate the exact result. Legacy fast-closure keys are accepted only during configuration migration, ignored, and recorded as adjustments; they are not part of the resolved schema. Mandatory topology validation is unconditional. Optional scientific validation changes only geometric acceptance and the center definition.
### Search Scope and Report Scope

`ring.size` / `--size` defines the shared face-size search universe. Ring and quasi-cage detection support 4/5/6/7. Cage detection intentionally uses only the 4/5/6 intersection of that universe.

For cage search, SQQ generates every trivalent Euler-compatible face composition up to YAML `cage.max_face`:

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
- `cage.report_type: auto` reports every detected cage allowed by `ring.size` / `--size`;
- YAML `cage.report_type` filters the user-facing cage counts, occupancy, GRO, info, and main-summary tables;
- cage report groups expand to exact compositions and duplicate types are removed:

```text
I     -> 5¹², 5¹²6²
II    -> 5¹², 5¹²6⁴
H     -> 5¹², 5¹²6⁸, 4³5⁶6³
HS-I  -> 5¹², 5¹²6², 5¹²6³
TS-I  -> 5¹², 5¹²6², 5¹²6³
I2II  -> 5¹²6³
```

- `cage.report_type` accepts `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, and `I2II`; group names may be listed together;
- `cage.report_type: all` reports every detected composition;
- all detected cages, including unreported types, still remove consumed half-cages, quasi-cages, and free rings;
- an explicitly requested cage group is rejected when one of its compositions requires a face size absent from `ring.size` / `--size` or exceeds YAML `cage.max_face`.

The default report scope is `auto`. Therefore, `-s 4,5,6` searches and reports every accepted cage composed of 4-, 5-, and 6-membered faces unless YAML `cage.report_type` narrows the report. For example, `cage.report_type: [I, II]` with `-s 4,5,6` keeps 4/5/6 ring and quasi-cage reporting while cage output is restricted to 5¹², 5¹²6², and 5¹²6⁴.

### Sparse Ring Topology

After ring detection, the pipeline builds one frame-local `RingTopologyIndex`. Its mandatory data are the stable ring-id map, each ring's sorted graph edges, and the sparse reverse index from one graph edge to the rings using it. Cage search consumes only this incidence topology. Ring centers, ring adjacency, cached distances, face-quality metrics, and least-squares normals are constructed lazily when half/quasi search, scientific validation, or hydrate-cluster analysis actually requires them. No topology or coordinate cache crosses a frame boundary.

The Python and C++ engines implement the same scientific contract with their native containers. Python uses sparse dictionaries/sets and C++ uses compact CSR-style incidence arrays plus local vectors. They do not share source code, but they use the same ordering, state, feasibility, closure, and validation rules.

### Exact Sparse GROW Logic

1. Generate every Euler-compatible 4/5/6 target composition in the configured scope and precompute target face-count and total-incidence constraints.
2. Start from each canonical seed. A single-ring seed cannot later add a lower-ranked ring, which prevents rediscovering one shell from every face. Patch seeds use an equivalent earlier-seed subset rule.
3. Initialize a branch with a sorted face list, face-membership set, local `edge_count` mapping, and the set of edges currently used once. Edges used twice are closed; any third use is rejected.
4. Encode compatible target compositions as a compact mask. Face-count excess, insufficient remaining face incidence, or incompatible parity removes an impossible target before expansion.
5. Record the exact sorted face tuple in a per-seed visited set. The set is released as soon as that seed completes rather than accumulating for the full frame.
6. If no open edge remains, require an exact target composition and submit the shell to the mandatory topology validator.
7. Otherwise inspect the open edges and choose one having the fewest currently addable rings. This minimum-remaining-value order changes traversal cost only.
8. Visit every eligible neighboring ring in stable rank order. `cage.max_boundary_candidate` is not a branch cap and never drops a candidate.
9. Add the ring by updating only its local edge counts, open frontier, face counts, and target mask; recurse; then undo those changes in place.
10. Continue until closure or target infeasibility. `cage.max_state_per_seed` and `cage.max_total_state` default to `0` (unlimited). A positive guard raises an error and rejects the complete frame; no partial cage result is accepted or published.
11. Deduplicate accepted objects by their complete face set and scientific identity. Final deterministic ordering by cage type, water membership, and face topology assigns equal Py/CPP cage sets equal frame-local cage IDs.

The old implementation attached a frame-wide mask to every active ring/edge and copied large mask states during growth. On large crystalline fixtures, memory therefore scaled with global topology multiplied by explored state. The sparse implementation stores global incidence once and makes branch memory proportional to the local shell boundary. It does not round coordinates, alter graph construction, restrict target types, or approximate the search.

Default acceptance criteria:

- every edge is used exactly twice;
- Euler characteristic satisfies `V - E + F = 2`;
- the face-adjacency graph is edge-connected;
- incident faces around every shell vertex form one cyclic vertex link;
- every shell vertex is trivalent;
- face counts match one generated target composition;
- the same scientific cage identity was not already accepted.

These topology checks are always applied before cage type/isomer assignment in SQQ-Py and SQQ-CPP. They reject disconnected, pinched, branched, and non-manifold shells independently of `cage.scientific_validation`.

### Optional Scientific Geometry Validation

`cage.scientific_validation = false` is the default and is independent of the analysis engine. `false` disables only the additional geometry checks; it never bypasses the mandatory topology checks above. Setting it to `true` additionally requires:

- each ordered ring face is locally unwrapped across PBC and fitted by SVD; its planarity RMS must not exceed `cage.max_face_planarity_rms_nm`;
- its cyclic O-O edge-length coefficient of variation must not exceed `cage.max_face_edge_cv`, and its projected area must be nonzero;
- the consistently outward-triangulated shell must have volume at least `cage.min_cage_volume_nm3`.

An accepted scientific-validation cage uses the tetrahedral volume centroid of the oriented triangle shell. The default path continues to use the mean of locally unwrapped cage-water coordinates. Consequently, enabling scientific validation may remove geometrically distorted cages and may change guest occupancy or geometry-resolved hydrate-cluster edges. The now-mandatory topology checks can reduce cage, isomer, occupancy, and hydrate-cluster results relative to builds that accepted non-manifold shells. Raw ring and half/quasi searches, order parameters, and ice classification are unchanged; ownership-filtered free-ring and free-patch outputs can increase when rejected cages no longer consume them.

The built-in scientific thresholds are `0.06 nm` planarity RMS, `0.35` edge-length CV, and `1.0e-6 nm^3` minimum volume. They are explicit configuration values rather than hidden constants.

## Hydrate Cluster

Hydrate_cluster is optional. Its built-in value is off, and the overall default remains off because the default engine value is `py`:

```text
hydrate_cluster.enabled = false
hydrate_cluster.min_cage = 2
```

Engine values `00` and `py` with resolved cluster search can run cluster analysis; values `99` and `cpp` reject it. The retained CLI control is `--find-cluster on/off`; YAML `hydrate_cluster.min_cage` owns the component threshold. The explicit find-cluster value has precedence over YAML and the engine preset. Cluster search populates every selected `info` and main-summary output but does not add an unselected output type. No preset includes `cluster-gro`; it must be selected explicitly. Output type `cluster-detail` controls the optional domain and one-row-per-cluster CSV files. Explicit `cluster-detail` or `cluster-gro` requires resolved cluster search on.

Hydrate_cluster uses `result.all_cages`, the complete detected cage set in the selected ring/search scope. YAML `cage.report_type` filters cage counts, occupancy, GRO, Markdown cage tables, and main-summary cage columns only; it does not filter the cluster graph or phase evidence. Cluster hierarchy/detail/domain records resolve cage IDs against the same complete set, so an unreported cage can remain topologically necessary without appearing in the report-scoped cage table.

The high-level hierarchy is informed by HTR+ ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)): classify hydrate type and polycrystalline boundaries on a cage-connection graph. SQQ does not copy the HTR+ implementation; it uses the explicit shared-face fingerprints and deterministic domain rules below.

### Physical shared-face cage graph

1. Use `result.all_cages` (falling back to `result.cages` only for legacy result objects).
2. Build `ring_id -> cage_ids` from each cage face list.
3. Treat one complete shared ring as a potential cage-cage edge.
4. A physical ring face can separate at most two cages. When ring geometry is available, locally unwrap the face, fit its plane, and select the best cage center on each side. When geometry is unavailable, accept only a face referenced by exactly two cages.
5. Find deterministic connected components in the resulting undirected graph.
6. Report components with cage count >= `hydrate_cluster.min_cage`; count cages in smaller components as isolated cages.

### Local phase fingerprints and spatial consensus

For every graph node, SQQ counts first-shell labels of the form `(neighbor cage type, shared face size)`. Strict evidence rejects unexpected label types and allows each expected count to differ by at most one.

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

The conservative sH composite remains supplemental high-confidence evidence: two nonadjacent `5^12 6^8` anchors, exactly six common `5^12` cages connected to both anchors through 5-ring faces, exactly six `4^3 5^6 6^3` cages, and at least one adjacent medium-cage bridge between the anchors.

A strict fingerprint is no longer the only way to initialize a domain. SQQ also scores distributed spatial evidence in the current frame:

1. `matched` is the sum of each observed expected label capped at its template count.
2. `coverage = matched / expected labels`; `purity = matched / all observed labels`.
3. `support` is the harmonic mean of coverage and purity.
4. A spatial candidate requires coverage >= 0.50, purity >= 0.50, and support >= 0.55.
5. Candidates are connected only by face labels allowed by both endpoint templates. Iterative removal of nodes with degree below two leaves the graph 2-core.
6. A retained component requires at least three cages, mean support >= 0.60, and a phase-defining hexagonal edge: `5^12 6^2`--`5^12 6^2` for sI, `5^12 6^4`--`5^12 6^4` for sII, or `4^3 5^6 6^3`--`5^12 6^8` for sH.
7. Spatial cores are anchored only by phase-specific cages: `5^12 6^2` for sI, `5^12 6^4` for sII, and `4^3 5^6 6^3`/`5^12 6^8` for sH. Shared `5^12` cages provide spatial support without becoming mandatory anchors for competing phases.

Strict and spatial evidence are internal. `seed_count` counts the contained evidence objects, and `seed_cage_count` counts their unique domain anchors or strict neighborhoods. Neither value is used as a public phase category.

### Expansion and exclusive domains

sI, sII, and sH expand independently from strict seed members and spatial-core anchors. A growth candidate must:

- use a cage type supported by the phase template;
- have at least one compatible internal fingerprint label;
- not exceed any expected fingerprint count by more than one;
- connect through a face label allowed by both endpoint templates; and
- receive at least two compatible contacts from already accepted phase cages.

For sH, the same edge check includes pentagonal `5^12` contacts, square medium-medium contacts, equatorial medium-large hexagonal contacts, and axial large-large hexagonal contacts. Strict/composite-seed internal edges remain trusted evidence; spatial-core edges must pass the ordinary template check.

After all phases collect claims independently, SQQ forms domains from cages claimed by exactly one phase. Domain edges must remain phase-compatible, and every connected domain component must contain at least one strict or spatial evidence anchor. Same-phase regions separated by non-domain cages remain separate domains. Domain ids are deterministic within a frame and are not tracked across frames.

All spatial scoring is frame-local. No previous/next frame, persistent cage id, hysteresis, or temporal smoothing is used. Therefore one frame analyzed alone has the same phase assignment as that frame analyzed in a trajectory or compatible multi-file group. True topology changes can still cross the fixed spatial thresholds.

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

`cluster-gro` is a separate search-dependent output. Engine values `00` and `py` write it only when cluster search resolves on and YAML `output.type` selects it. Search off writes no cluster GRO output; reusing an output directory removes stale grouped or flat cluster GRO files generated by SQQ. The ordinary `gro` umbrella remains limited to ring, half-cage, quasi-cage, cage, and ice files.

For each frame, output aggregates cage IDs across every reported cluster and every domain of the same category:

- sI is the union of all sI domain cage IDs;
- sII is the union of all sII domain cage IDs;
- sH is the union of all sH domain cage IDs;
- boundary is the union of all generic boundary cage IDs.

Ambiguous, residual unclassified, and isolated/below-threshold cage IDs are intentionally omitted. Each category maps its cages back to water membership, deduplicates waters, and writes complete water molecules in original frame order. Guests, CNT atoms, and other non-water molecules are never added. An absent category has no file unless `output.write_empty_file` is true.

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

Public motif output is not generated in the current development version. The compatibility motif return slot remains empty, and neither a `Hydrate Motif` Markdown section nor a `hydrate_motif` main-summary table is written.

## Cross-Frame Cage Tracking

Tracking is a temporal identity layer over the accepted closed cages of each selected frame. It does not change the water graph, ring set, cage search, cage type, occupancy, or hydrate phase/domain assignment. Each frame is reduced to a compact snapshot containing frame/time/source, orthorhombic box, frame-local cage ID, canonical cage type, sorted face-size topology, stable member-water identities, wrapped cage center, phase labels, and non-exclusive guest IDs. A water identity is its one-based topology atom position, not the width-limited serial stored in a GRO atom record.

### Persistent identity

Candidate generation is sparse by water identity: an inverted current-frame water-to-cage map counts shared waters only for pairs that actually overlap. A pair must satisfy all configured minimums for shared-water count, Jaccard similarity, and shared fraction relative to the smaller cage. `max_center_distance_nm`, when non-null, is an additional orthorhombic minimum-image guard. Face-size multiset similarity, center proximity, and optional guest similarity rank surviving candidates; guest continuity contributes only a small tie-break term and cannot create a candidate without sufficient water overlap.

One deterministic maximum-weight bipartite assignment selects the global one-to-one continuation set rather than accepting cage pairs greedily. Stable snapshot sorting and assignment tie handling give newly born cages persistent IDs `t1`, `t2`, ... in deterministic order. A compatible cage retains its ID through cage-type or phase-label changes. Frame-local cage IDs remain stored in observation metadata for auditability.

`track.gap_frame` counts consecutive selected frames in which a cage may be absent. The default `0` expires an unmatched cage immediately. A positive value keeps a dormant state for at most that many selected frames; any later match records the exact gap on the observation and emits a `gap` event. Snapshot indexes must remain consecutive, so an analyzed frame with no cages is represented explicitly as an empty snapshot and cannot disappear from the time axis.

### Events, lifetime, and guest residence

Track events are `birth`, `death`, `type_change`, `phase_change`, `split`, `merge`, and `gap`. Births in the first selected frame and tracks still active in the last selected frame are marked left/right censored. An uncensored death ends at the first absent selected frame rather than the last frame in which the cage was seen. `cage_track.csv` stores one lifetime sample per persistent cage; `lifetime_distribution.csv` aggregates exact frame/time samples and reports censoring counts. Population rows report total, cage-type, and phase counts for every selected frame.

Guest residence is calculated as contiguous observation episodes within a persistent cage. A cage gap splits an episode even if the same guest is present before and after it. One guest may contribute residence episodes to multiple cages because occupancy itself is non-exclusive.

### Analyze integration

Analyze has a dedicated tracking sink independent of optional rendering. It converts each successful `FrameResult` immediately to a compact snapshot and feeds one `TrackingAccumulator` per topology group; it does not retain all frame results. For every complete selected-frame sequence, Analyze writes `track_state.json` and six normalized CSV tables under `track/`, whether or not `sqq-render` is selected. When rendering is enabled, the completed persistent-ID result is passed to `RenderSession`, which atomically rewrites cage-center (`C`) and cage-membership (`M`) TSV records before publication. A failed frame never gets silently compressed out of a continuous trajectory: that topology group's persistent state is skipped with a diagnostic while already valid per-frame scientific outputs remain available.

The versioned JSON state contains frame stamps, configuration, tracks/observations, and events. Schema version 2 is written; version 1 can be migrated on read. A future matching configuration cannot be retroactively applied in `--source` mode because the assignments have already been made; Track reports the stored settings.

### Track workflow and targets

`sqq track` has two mutually exclusive data paths:

- raw `-i` input reuses Analyze readers, frame selection, engine dispatch, and `-dt`, then tracks the resulting stream;
- `--source` imports an existing 0.5.1 Analyze `track/track_state.json` and its complete `sqq_render/` bundle without rerunning frame science.

Source mode reads the Analyze engine from `sqq_config_resolved.yaml`, preserves the original `sqq-py` or `sqq-cpp` identity in terminal and newly resolved metadata, and applies C++ half/quasi capability normalization when the source engine is SQQ-CPP.

```bash
sqq track --source ./result_sqq --target all -o ./result_track
```

With neither option, exactly one state is discovered from the current directory. Raw Track accepts one trajectory or one stacked GRO because its target is one physical time series. It is currently serial: `-w` / `--worker` is normalized to one and the backend to `serial`. Track has a fixed configuration/state/CSV/target-render output set; `--output-type` does not alter it. A source containing several topology groups must be narrowed to one `result_A`/`result_B` directory.

`--target` accepts `all`, canonical or compact cage types, `sI`/`sII`/`sH`/boundary categories, and persistent IDs. A phase target automatically resolves `find_cluster` to `on` and requires SQQ-Py; SQQ-CPP phase targets fail before analysis. Source mode cannot add phase labels retroactively, so imported state must already contain them. Comma-separated mixed targets are de-duplicated and written to independent `all/`, `type_<type>/`, `phase_<phase>/`, or `cage_<tID>/` directories. Type and phase targets include the full lifecycle of any track that matched at least once, permitting residence and transition analysis outside the matching frames. Persistent-ID targets validate that the requested ID exists.

The run-level `track/` directory contains state plus `cage_observation.csv`, `cage_track.csv`, `cage_event.csv`, `cage_population.csv`, `guest_residence.csv`, and `lifetime_distribution.csv`. Each target directory contains filtered copies, `track_info.md`, and `sqq_render/{sqq_track.gro,sqq_track.xtc,sqq_track.membership.tsv,sqq_track.vmd.tcl}`. A raw-input persistent-ID target with pre-birth frames uses two passes: the first establishes its cross-frame `tID`; the second reanalyzes only the required prefix through the birth frame, follows the birth-cage water set backwards through dispersed, connected, ring, half/quasi, and cage states, writes `precursor_state.csv` and `water_history.csv`, and adds precursor membership to that target's VMD package. A target already present in the first selected frame has no precursor interval. Imported state lacks full pre-cage graph/ring/patch objects and writes an explicit unavailable record instead.

The shared Tcl renderer recognizes persistent IDs in both Analyze and Track packages. `sqq target save` atomically writes the current selected ID or comma-separated IDs to `sqq_target.txt` beside the render files for reuse in a later `--target` command.

### Streaming and validation

The accumulator retains frame stamps, normalized observations, events, and only active/dormant cage states; it does not retain the snapshot sequence. The 10,000-frame regression completed in about 5 seconds with roughly 10 MiB traced Python peak allocation and 57 MiB peak process working set on the validation host. Scientific tests cover membership change, orthorhombic PBC translation, type/phase change, birth/death, split/merge, guest episodes, deterministic IDs, exact `-dt` selection, gap 0/1, state/CSV round trips, source/raw CLI, and mixed targets. Both Track paths finish through the shared final `Tracking Results` and feature-aware `Citation Recommendation` page described above.

## Guest Occupancy

Cage centers are computed from the locally unwrapped O coordinates of the cage waters. Hydrogens are not used in the cage center.

For each accepted cage, SQQ checks all selected guest molecules. `guest.center_mode` is active:

- `center_atom` (default) and `auto` use the configured center atom when present. The built-in `guest.center_atom` maps `CH4`, `CO2`, and `MET` to atom name `C`, so these residues use their carbon atom; a guest without a matching configured center falls back to a PBC-aware residue centroid;
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

One shared job-shape rule owns per-frame placement for both engines and every scheduler. A single ordinary one-frame GRO or XYZ uses the compact frame-root layout. Any trajectory-like input, including one or more XTC/TRR/LAMMPS/DCD paths or a stacked GRO, and any job containing multiple independent GRO/XYZ files uses the separated layout: report/TSV files under `info/` and selected structures under `gro/<frame>/`.

Output selection is positive and engine-preset-specific. The four accepted engine values are:

| Engine value | Default `output.type` |
| --- | --- |
| `00` | `info,sqq-render,summary-xlsx` |
| `py` | `info,sqq-render,summary-xlsx` |
| `99` | `info,sqq-render,summary-csv,summary-detail-csv` |
| `cpp` | `info,sqq-render,summary-csv,summary-detail-csv` |

Supported SQQ-Py canonical names are `info`, `membership-tsv`, `order-tsv`, `f3-gro`, `f4-gro`, `sqq-render`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `cluster-gro`, `summary-xlsx`, `summary-csv`, `summary-detail-csv`, and `cluster-detail`. `gro` expands to the five ordinary ring/half/quasi/cage/ice subtypes. `sqq-render` writes one indivisible four-file visualization package. `cluster-gro` is separate and requires resolved cluster search. `default` expands to the selected engine preset and may be combined with extra types. `all` and `none` are exclusive.

SQQ-CPP accepts `info`, `gro`, `cage-gro`, `f3-gro`, `f4-gro`, `sqq-render`, `summary-csv`, `summary-xlsx`, `summary-detail-csv`, `default`, `all`, and `none`. `gro` enables the supported ordinary classified cage GRO path. Neither `99` nor `cpp` selects `gro` or `cage-gro` by default. Cluster types are unsupported. `sqq_config_resolved.yaml` is mandatory for every engine value.

For any separated job (multiple standalone GRO/XYZ files, one or more XTC/TRR/LAMMPS/DCD trajectories, or one stacked GRO), the requested root directly owns `sqq_config_resolved.yaml`, optional `summary.xlsx`, optional `summary/`, `info/`, optional `gro/`, and optional `sqq_render/`. All per-frame Markdown reports are flat files under `info/`; selected per-frame structures are under `gro/<frame>/`. A single ordinary one-frame GRO or XYZ retains its compact frame-root layout. The render directory contains:

```text
sqq_render/
  sqq_cage.gro              # stable topology and first selected frame
  sqq_cage.xtc              # every selected frame coordinate set and box
  sqq_cage.membership.tsv   # typed frame/center/guest/membership manifest
  sqq_cage.vmd.tcl          # VMD loader and interactive commands
```

Every completed Analyze sequence owns the following tracking state and six CSV tables independently of `sqq-render`:

```text
track/
  track_state.json
  cage_observation.csv
  cage_track.csv
  cage_event.csv
  cage_population.csv
  guest_residence.csv
  lifetime_distribution.csv
```

For 2-26 topology groups, the requested root owns the batch manifest and each `result_<letter>/` owns the corresponding group files, including its own `sqq_render/` when selected. Per-frame Markdown and optional membership/order TSV files are routed to the group's `info/`; selected structure output is routed below `gro/<frame>/`. The same separation rule applies to serial and process-parallel trajectory frames. If `summary-xlsx` and `summary-csv` are both selected, `summary.xlsx` and `summary/` coexist in the same group root.

For more than 26 topology groups, output normalization is replaced by the information-only safety selection for the complete run. The root contains `sqq_config_resolved.yaml` plus `info/*_info.md`; no group directories, summaries, detail CSV, GRO, or `sqq_render/` package is produced. The root manifest records the detected group count, source assignments, warning, requested output types, and effective information-only selection.

The four files in `sqq_render/` are one indivisible visualization package when `sqq-render` is selected. `sqq_cage.gro` contains one standard GRO block: the stable atom topology and coordinates/box of the first selected render frame. YAML `render.atom_scope` defaults to `full`, so this topology and every XTC frame contain all normalized input-frame atoms in original order: water hydrogens, complete guests, additives, environment/wall components, and other retained atoms. `compact` selects only water oxygens and complete guests. Wrapped coordinates are preserved without PBC reimaging. `sqq_cage.xtc` stores coordinates and orthorhombic boxes for every selected render frame with physical times when available. `sqq_cage.membership.tsv` is a sparse typed manifest: `F` records identify each render/source frame and time; `C` records store each cage center, wrapped into the primary orthorhombic box in nm and then converted to angstrom; `G` records map every rendered guest atom to its complete molecule group; `M` records preserve atom-level cage/guest membership plus cage type and optional cage/phase/domain/cluster IDs; and `P` records map rendered atom indices to component role and residue name. Large component groups are emitted as bounded repeated `P` records rather than one unbounded TSV field. Atoms without cage membership are absent only from `M`; full-scope atoms remain present in every XTC frame and are classified through `P`. Multiple memberships and complete multi-atom guest groups are retained rather than flattened. Finalization validates stable atom identity/order, applies the already completed persistent-ID result when available, and atomically replaces the four destination files from run-private fragments. Analyze and raw Track share this writer under both engines; source Track inherits the imported bundle's topology and atom scope.

Every newly generated Analyze Tcl embeds a comment-only JSON render-package manifest between `SQQ-RENDER-MANIFEST-BEGIN` and `SQQ-RENDER-MANIFEST-END`. It records `schema`, package `kind`, and a general list of actual file `role`, relative `path`, and `required` state. The manifest is stored inside `sqq_cage.vmd.tcl`, so the package remains four files. Tcl path declarations and manifest entries are generated from the same filename arguments.

`sqq vmd [PATH]` discovers recognized SQQ `.vmd.tcl` files without executing them. `PATH` may be omitted, a result/render directory, or a specific script. Manifest paths are resolved relative to the script and checked with file metadata only; required files must exist, be regular files, and be nonempty. Unsafe parent/absolute paths and unresolved template placeholders make a package incomplete. Older Tcl scripts without a manifest use restricted static parsing of their declared SQQ paths. Multiple scripts are sorted by normalized absolute path. Displayed absolute paths use VMD/Tcl-compatible forward slashes; terminal arguments use Windows quoting on Windows and POSIX shell quoting on macOS/Linux. Exit status is `0` when every package is complete, `1` when a recognized package is incomplete, and `2` for an invalid path or no recognized SQQ script.

`sqq vmd -h` adds locator usage above the same command guide printed by Tcl `sqq help`, `sqq -h`, and `sqq --help`; one Python definition generates both outputs.

`sqq_cage.vmd.tcl` resolves its neighboring files relative to the script path, opens the topology GRO, attaches the XTC frames, reads membership records for the current render frame, and builds SQQ-owned representations. The user keeps all four files together and sources only the Tcl script:

```tcl
source {path/to/result/sqq_render/sqq_cage.vmd.tcl}
```

The script sets `color Display Background white`, initializes the default opaque cage-all view, and reports `SQQ graph: <effective-mode>` when sourced. Later redraws remain silent unless the effective graph mode actually changes. Its public grammar is:

```text
sqq show <family> <target...> [<family> <target...>]...
sqq color <family> <target...> <VMD-color|ColorID|default>
sqq show label [on|off]
sqq pick center|guest|off
sqq target save
sqq clear
sqq help | sqq -h | sqq --help
```

Families are `cage`, `guest`, `phase`, `cluster`, `domain`, and `component`. In `show`, every family token begins a group and consumes at least one following target up to the next family token; one command may therefore contain several additive groups, such as `sqq show cage 512 guest 512 phase sI`. Cage targets are `all`, a canonical/generic cage label, its delimiter-free alias, or an exact persistent ID such as `t133`; fallback frame-local IDs remain readable when a partial sequence cannot publish Track state. Guest targets reuse the cage target namespace and mean guests assigned to all cages, to a cage type, or to one exact cage ID. Phase accepts `all`, `sI`, `sII`, `sH`, `boundary`, `ambiguous`, `unclassified`, and `isolated`; cluster/domain accept `all` or exact frame-local IDs. Component accepts `all`, `water`, `guest`, `additive`, `environment`, `other`, or an exact residue name. `sqq show component ...` and `sqq color component ...` expose full-frame context without changing the startup view, which remains opaque `sqq show cage all` with context hidden. `color` deliberately remains a one-family command. Bare inferred forms from 0.3.4 remain removed.

The renderer stores the active view as deduplicated family/target state. On source or reset, it contains a synthetic, replaceable `cage all` default. The first `sqq show ...` replaces the default; each later `show` merges every supplied family group without clearing existing layers. Labels are an independent state and default off; `sqq show label` toggles them unless an explicit `on`/`off` is supplied. `sqq pick center` and `sqq pick guest` are mutually exclusive modes. Center mode draws one yellow graphics sphere and pickpoint for each current-frame cage center and enters VMD Atom Label mode, because Query mode only prints information and does not dispatch the graphics callback required by SQQ. The callback tag is read from VMD pickpoint metadata instead of assuming it equals the graphics ID. Selecting a yellow center identifies the exact cage, while water-atom clicks are ignored. Guest mode enters VMD Pick mode, maps any picked guest atom to its complete molecule, then highlights the molecule and every cage membership; a guest with no cage membership clears the transient selection and is not highlighted. Both modes make the unselected context transparent. A persistent yellow DynamicBonds layer highlights selected cages, while guest mode adds a persistent orange CPK layer for the selected guest. Pick callbacks update these layers with `mol modselect` and do not delete the representation currently owned by VMD's mouse event. `sqq target save` atomically writes the selected persistent cage ID or IDs to `sqq_target.txt`. `off` leaves pick mode. A frame change clears the transient selection and rebuilds the current-frame centers, guest map, and view while retaining the chosen pick mode. `sqq clear` removes custom selection, color, label, and pick state and restores opaque `cage all`.

Multi-atom guests retain the complete molecule, and one guest may belong to several cages. Cage representations use DynamicBonds; guest representations use CPK. Cage and guest colors have independent override maps. Cross-family rendering always follows `phase -> cluster -> domain -> cage -> guest`, so guests remain visible regardless of command order. This family order is separate from the cage-topology order: a single cage layer uses a 0.125 angstrom cylinder radius (0.250 angstrom diameter); multiple cage layers remain bounded from 0.125 through 0.130 angstrom and are ordered as nonstandard below `512 < 51262 < 51263 < 51264 < 435663 < 51268`, with exact-ID highlights last. Cage identifiers use persistent Track IDs when state publication succeeds; cluster and domain identifiers remain frame-local.

The renderer records VMD's stable representation names after each `mol addrep` and deletes only those names on the next redraw; user-created representations are preserved. Frame traces keep at most one pending idle callback, so rapid animation updates coalesce to the final frame. Re-sourcing removes old frame and pick traces, cancels any pending callback, deletes prior SQQ graphics, and installs exactly one clean callback set before resetting color/selection/representation state.

Per-frame output folders keep the configured grouped/flat structure for ordinary category files. Generated GRO paths and title lines use ASCII-only structure labels (`5^12`, `5^126^2`, `hc_5r_5^5`, and similar); Unicode superscripts/subscripts remain display-only in Markdown and main summaries.

Analyze holds a nonblocking process lock on the requested output root for the complete run. The hidden lock file may remain as metadata, but the operating-system lock is released automatically when the process exits; a concurrent Analyze process targeting the same root fails before cleanup or result writes. Multi-GRO topology groups each receive a separate run workspace under their own result root.

Each invocation writes visualization fragments to a unique private staging directory and passes that exact path to every serial, thread, and process worker. Finalization sorts and validates all manifests, writes temporary GRO/XTC/TSV/Tcl files in the destination filesystem, and publishes them with atomic replacement. Cleanup retries transient `ENOTEMPTY`, `EBUSY`, `EACCES`, and `EPERM` errors with bounded backoff. Failure to remove staging after a successful publish is nonfatal: SQQ reports the retained directory and preserves finalized output. Scientific or validation failures remain fatal, and cleanup cannot replace the original exception.

Non-render per-frame output follows the same transaction boundary. A worker writes the selected info, membership/order TSV, F3/F4 GRO, and ordinary/category GRO files below a run-private frame directory. Only a fully successful frame is moved into its final `info/` and `gro/<frame>/` destinations. A failed writer discards the staged frame, including empty directories, without exposing a partial frame bundle.

When an Analyze output directory is reused, SQQ removes only known generated files that are outside the new effective selection and abandoned legacy/run-isolated fragment workspaces. Unknown user files remain untouched.

Cleanup is not limited to the current source stem. The previous resolved run manifest supplies prior input names and layout state, and known per-frame artifacts are removed across old/new source names and both grouped/flat paths. F3/F4 GRO files participate in the same cleanup. Unknown files and directories that do not match SQQ's generated-file registry are preserved.

Engine values `00` and `py` write per-frame info, the four-file run-level `sqq_render/` package, and `summary.xlsx`; ordinary/category/cluster GRO and detail CSV files are not default output. `summary-xlsx` owns the workbook. `summary-csv`, `summary-detail-csv`, and `cluster-detail` write disjoint CSV filenames under the single `output.summary_csv_dir`, which defaults to `summary`.
`cluster-detail` separately owns `hydrate_domain.csv` and `hydrate_cluster_detail.csv`. Explicit `cluster-detail` and `cluster-gro` require resolved cluster search on. Cluster search populates selected info and main-summary outputs without adding an unselected output type. Every engine preset requires `cluster-gro` explicitly through YAML `output.type`. Search off removes stale SQQ-generated grouped directories or flat category filenames. If no per-frame type remains selected, the pipeline removes an empty frame directory. Unrelated files are preserved.

Per-frame output folders use the default grouped structure:

```text
frame_name/
  frame_name_info.md
  frame_name_order_parameter.tsv  # with order-tsv and selected F3/F4/Q_l
  order/
    frame_name_f3.gro           # with f3-gro and selected F3
    frame_name_f4.gro           # with f4-gro and selected F4
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

`summary-xlsx` and `summary-csv` use the same main-table builder. SQQ-Py includes `summary`, optional `failures`, the effective connection table, `ring`, `half_cage`, compact `quasi_cage`, `cage`, optional `hydrate_cluster`, selected `order_parameter`, `ice`, and `detail_index` when detail files exist. SQQ-CPP includes the applicable subset: `summary`, optional `failures`, `cage`, selected `order_parameter`, and optional `detail_index`. XLSX writes sheets; CSV writes one UTF-8-SIG file per table under `summary_csv_dir` (default `summary/`). The dashboard contains the compact configuration and min/mean/max analysis results; failure and detail-index tables keep their own row units. Failure details are always retained in `sqq_config_resolved.yaml`.

The explicit `f3-gro` and `f4-gro` outputs use the shared Python writer after either engine has produced `F3F4Result`. A matching order parameter must be selected. For each parameter, the writer selects exactly the waters whose per-water value is defined, expands each oxygen to the full `Water.atoms` tuple, preserves source atom order and the frame box, and annotates only the oxygen. The fixed-width atom prefix occupies columns 1-44; the existing three velocity fields occupy columns 45-68 when present and are otherwise blank; `; SQQ F3=<value>` or `; SQQ F4=<value>` begins at column 69. Values use eight decimal places. With grouped layout the files are `<frame>/order/<frame>_f3.gro` and `<frame>/order/<frame>_f4.gro`; separated multi-frame/multi-file output adds the leading `gro/`. `output.write_empty_file` controls zero-valid-water files. These writers expose existing values only and do not alter graph construction or F3/F4 calculation.

The SQQ-CPP main summary deliberately contains only applicable tables: `summary`, optional `failures`, `cage`, selected F3/F4 `order_parameter`, and `detail_index` when detail CSVs exist. Its default `summary-csv` writes these as independent files, while explicit `summary-xlsx` writes the same mapping as workbook sheets. Default `summary-detail-csv` adds `cage_occupancy.csv` and `cage_isomer.csv` in the same `summary/` directory. Ring/connection diagnostic tables, half/quasi, cluster, and ice are omitted from the native schema. If no selected guests exist, the dashboard and per-frame info state that occupancy was not evaluated.

Each per-frame `*_info.md` report is optimized for inspection rather than plotting:

- Frame Information begins `sqq version`, `SQQ engine`, `date & time`, `source`, input provenance, `frame`, and `time_ps`; the engine value is normalized to `sqq-py` or `sqq-cpp`, never `py (sqq-py)`;
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

The optional `summary-xlsx` workbook and `summary-csv` directory keep plotting-oriented frame-analysis tables with one input file or trajectory frame per row; dashboard, failure, and detail-index tables use their metadata-specific row units. Their `hydrate_cluster` table reports mutually exclusive classified, boundary, ambiguous, and unclassified cage counts. Their `quasi_cage` table aggregates exact isomers into composition-level columns such as `5r_5²6³`. `summary-detail-csv` writes `cage_occupancy.csv` and `cage_isomer.csv` for both engines, plus `quasi_cage_isomer.csv` for SQQ-Py. With cluster search and `cluster-detail`, `hydrate_domain.csv` and `hydrate_cluster_detail.csv` are added. All main/detail/cluster CSVs share `summary_csv_dir` and have disjoint filenames. Cluster detail includes the four category ID groups and `boundary_composition`; domain detail uses `external_boundary_contact_count` and `external_boundary_contact_ids` for direct external adjacency. `detail_index` records only generated detail CSV files. `cage_isomer.csv` defaults to nonzero isomer rows plus per-frame totals; YAML `output.cage_isomer_row = all` restores the full zero-filled matrix. Public motif output is not written. The `order_parameter` table reports only selected F3/F4 mean/count pairs, selected `qN` mean/count pairs, and selected MCG/DHOP largest clusters. Matching focus-water mean/count columns are added only when `order_parameter.focus_water` is non-empty. MCG without a configured guest remains `N/A`. `order_parameter.enabled: []` or `--order-parameter none` omits the table. Output type `order-tsv` contains only selected per-water F3/F4/Q_l columns; MCG/DHOP remain frame-level outputs, so a hydrate-only selection does not create an otherwise empty `*_order_parameter.tsv`.
Summary generation records all table dimensions and write/format/save timings in `run.summary_write`. Exact quasi isomers are carried as sparse records until `quasi_cage_isomer.csv` is built, avoiding a summary DataFrame column for every observed isomer; composition-level quasi counts remain in the compact tables. Main CSV, `summary.xlsx`, every detail CSV, and `sqq_config_resolved.yaml` are written through same-directory temporary files and atomically replaced after success. Detail CSV replacements/removals commit as one recoverable bundle. Stale CSV cleanup is confined to known SQQ-generated filenames inside the configured `summary_csv_dir`; unknown files are preserved. Known files in the retired default `summary_detail/` directory are removed during migration. Before any table is handed to pandas Excel output, SQQ checks the 1,048,576-row and 16,384-column workbook limits and reports an actionable error for an unexpected oversize compact sheet. For an analysis sheet above 200,000 cells or 128 columns, Excel keeps header style, filter, freeze pane, and fixed widths but skips per-body-cell formatting and broad auto-width scans. This I/O policy does not change values, row/column schemas, or CSV selection.

Summary publication is one recoverable transaction over the selected XLSX, main CSV, detail CSV, cluster-detail CSV, and stale-known-file removals. All replacements are prepared before commit; if construction or publication fails, staged files are removed and the preceding complete summary set is retained rather than exposing a mixed generation. Output-root reuse also removes known generated files from obsolete `result_A`-`result_Z` roots when the next run contains only one system; unrelated user files remain untouched.

## Current Limits

- The implemented PBC path remains orthorhombic. Non-orthogonal/triclinic GRO boxes and trajectory cell angles are detected and rejected; conversion must occur before SQQ analysis.
- SQQ-CPP is intentionally limited to graph, internal chordless 4/5/6 rings, cage topology/isomer/occupancy, and F3/F4. The SQQ-Py engine is required for public rings, half/quasi, cluster, ice, Q_l/MCG/DHOP, legacy per-frame VMD, or TSV output. The full-frame-by-default GRO/XTC/membership-TSV/Tcl renderer is shared by both engines.
- XYZ input has configurable coordinate scaling but no periodic box metadata.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- Hydrate domains are per-frame topological regions; temporal grain tracking and crystallographic orientation matching are not implemented.
- Cage Track follows cage identity and records each observation's phase labels; it does not assign a persistent identity to a whole hydrate domain, cluster, or crystallographic grain.
- `sqq track --source` cannot reconstruct pre-cage water/ring/patch history because Analyze state stores compact cage snapshots. Use raw-input Track for a persistent-ID precursor history.
- Raw Track accepts one trajectory or one stacked GRO physical system. Multiple incompatible Analyze groups must be selected through one explicit `result_A`/`result_B` source directory.
- Hydrate phase classification depends on the detected cage/search scope but not on YAML `cage.report_type`; changing ring search sizes, cage face limits, or detection thresholds can still change the available topology.
- Boundary membership is a per-frame first-layer topological classification; transition-path kinetics, temporal domain tracking, and crystallographic orientation matching are not implemented.
- Cluster GRO output preserves original wrapped coordinates and the original box rather than reimaging each category; periodic or percolating structures may therefore retain unavoidable cross-box seams in a single-copy file.
- Default `quasi_cage.search_policy = bounded` is not exhaustive for large outer-layer components. Opt-in `exact` enumerates connected subsets but remains subject to explicit candidate and state limits.
- Automatic process workers parallelize independent GRO/XYZ files or selected frames of one indexed XTC/TRR/LAMMPS trajectory. Multiple GRO topology groups share one pool rather than consuming workers group by group. Worker counts are based on physical cores with one physical core reserved for the system; topology search inside one individual frame remains single-process.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
- MCG is meaningful only for guest residue names selected in `hydrate_order.mcg_guest_resname`; other guest species are not silently treated as methane.
- Published DHOP transition-state thresholds are model- and condition-dependent; SQQ reports the descriptor and does not assign a universal critical-nucleus threshold.
