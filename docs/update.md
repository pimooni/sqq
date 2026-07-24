# SQQ Update Notes

This file records versioned update notes. New releases should be appended above older entries.

## Version 0.3.7

### Short Summary

Version 0.3.7 stabilizes frame-local sI/sII/sH phase assignment without temporal smoothing. Exact local fingerprints remain high-confidence evidence, while a conservative spatial-consensus path can initialize a coherent domain when no single cage has a perfect first shell. Phase-specific anchors, compatible graph 2-cores, and the existing two-contact expansion prevent one marginal contact from switching an entire domain on or off. Public phase and boundary categories remain unchanged, and a frame produces the same assignment alone or in a compatible batch. Package and native-core metadata are synchronized at `0.3.7`, released Jul 25, 2026.

### Main Changes

1. Distributed spatial phase evidence
   - Added coverage, purity, and harmonic-support scoring for partial sI/sII/sH first-shell fingerprints.
   - Requires coverage >= 0.50, purity >= 0.50, support >= 0.55, and component mean support >= 0.60.
   - Prunes the phase-compatible candidate graph to its degree-2 core and requires at least three cages plus the phase-defining large-cage hexagonal connection.
   - Uses `5^12 6^2`, `5^12 6^4`, and `4^3 5^6 6^3`/`5^12 6^8` as phase-specific spatial anchors; shared `5^12` cages do not become mandatory competing anchors.

2. Deterministic per-frame domains
   - Retains strict fingerprints and the supplemental sH composite as high-confidence evidence.
   - Expands strict and spatial evidence through the established mutually compatible face labels and two-contact rule.
   - Uses only the current frame. No previous/next frame, persistent cage id, hysteresis, interpolation, or moving average affects classification.
   - Keeps sI, sII, sH, boundary, ambiguous, unclassified, and isolated ownership rules mutually exclusive as before.

3. Regression and trajectory validation
   - Synthetic tests cover a coherent sI spatial core with no strict seed and sparse/mixed negative fingerprints.
   - On the supplied 815--845 ns diagnostic window, zero-sI frames decreased from 15 to 0, on/off transitions from 10 to 0, maximum adjacent change from 53 to 18 cages, and total variation from 542 to 205.
   - The specified 410 ns `51262_00003` and 1000 ns `51262_00053` counterexamples remain boundary cages rather than sI cages.
   - Current 0.3.x tests and boundary regressions pass; archived tests tied to removed 0.2.x interfaces remain historical fixtures.

4. Metadata and acknowledgements
   - Synchronized Python, C++ core, CMake, and wheel-publish checks at `0.3.7` with release date Jul 25, 2026.
   - Added Zhenchao Li and Yifei Hu to the surname-alphabetized README acknowledgements and retained only the approved ordering sentence.

### Compatibility

- No CLI option, output type, table, or public phase name is added or removed.
- Phase counts can intentionally differ from 0.3.6 when a coherent spatial domain previously disappeared because its last strict local seed was missing.
- Cluster search remains SQQ-Py-only and opt-in except in mode `00`; SQQ-CPP scientific outputs are unchanged apart from version metadata.

## Version 0.3.6

### Short Summary

Version 0.3.6 makes the generated VMD renderer additive across commands and across object families. A single `sqq show` can combine cage, guest, phase, cluster, and domain groups; the first show replaces the startup cage-all view, later shows add deduplicated selections, and `sqq clear` restores the default view and colors. Cross-family rendering is deterministic with guests last. The source welcome and full Tcl help now describe the same interface. Scientific analysis, annotations, counts, and report values are unchanged. Package and native-core metadata are synchronized at `0.3.6`, released Jul 23, 2026.

### Main Changes

1. Compound additive `show`
   - Extended the grammar to `sqq show <family> <target...> [<family> <target...>]...`.
   - Added cross-family commands such as `sqq show cage 512 guest 512` while retaining every valid single-family form.
   - The first `show` after source/reset replaces the initial `cage all` view; later `show` commands add to the current view.
   - Canonical family/target selections are deduplicated, including repeats within one command and across commands.

2. Reset and deterministic rendering
   - `sqq clear` removes custom selections and color overrides, restores the initial cage-all view, and rearms first-show replacement.
   - Cross-family layers always render as `phase -> cluster -> domain -> cage -> guest`, independent of command order.
   - The cross-family order is separate from the existing fixed cage-topology priority.
   - `sqq color <family> <target...> <color>` remains a one-family command.

3. Source welcome and full help
   - Replaced the verbose source output with a compact `SQQ VMD Renderer` welcome showing the default view, additive mode, command forms, and two examples.
   - `sqq help`, `sqq -h`, and `sqq --help` remain equivalent and now print the complete usage, family descriptions, layer lifecycle, reset behavior, and examples.

4. Version metadata
   - Synchronized Python package, native C++ core, CMake, and publish checks at `0.3.6`.
   - Retained the release date Jul 23, 2026.

### Compatibility

- Existing single-family commands such as `sqq show cage 512 51264` remain valid.
- `sqq color` does not accept compound family groups; issue separate color commands for separate families.
- Regenerate `sqq-render.vmd.tcl` to use compound `show`, additive state, `sqq clear`, deterministic family order, and the revised help.
- These changes affect VMD interaction and representation order only; cage detection, scientific classification, annotations, and report values are unchanged.

## Version 0.3.5

### Short Summary

Version 0.3.5 replaces mode `50` with the explicit single-worker `py` default, makes all four presets GRO-clean by default, adds family-qualified VMD cage/guest controls, and moves authoritative runtime metadata from `run_config.yaml` to `config.yaml`. Detailed configuration sheets are removed from XLSX/CSV summaries while the compact dashboard Configuration block remains. The release also adds project acknowledgements in surname-alphabetical order without implying contribution rank. Scientific topology and order-parameter algorithms are unchanged. Package and native-core metadata are synchronized at `0.3.5`, released Jul 23, 2026.

### Main Changes

1. Mode names and exact defaults
   - Removed mode `50` and introduced `py` without a compatibility alias; `py` is the command default.
   - Retained modes `00`, `py`, `99`, and `cpp`.
   - Modes `py` and `cpp` default to one worker; modes `00` and `99` retain the 100% automatic-core policy with one physical core reserved. Explicit `-w` remains authoritative.
   - Modes `00` and `py` default to `info,sqq-cage-gro,sqq-render,summary-xlsx`.
   - Modes `99` and `cpp` default to `info,sqq-cage-gro,sqq-render,summary-csv`.
   - No mode emits ordinary, classified, or cluster GRO files by default. The corresponding `--output-type` values remain available explicitly.
   - Mode `00` still enables cluster search; that search populates selected info/main-summary output without implicitly selecting `cluster-gro`.

2. Explicit VMD cage and guest commands
   - The generated script starts with `sqq show cage all` and prints concise help when sourced.
   - Added `sqq -h` and `sqq --help` as equivalents of `sqq help`.
   - Replaced inferred object syntax with `sqq show <family> <target...>` and `sqq color <family> <target...> <color>` for `cage`, `guest`, `phase`, `cluster`, and `domain`.
   - Cage and guest targets accept `all`, cage types, exact cage IDs, and multiple targets where applicable.
   - Cage networks use DynamicBonds; guests use CPK. Cage/guest colors are independent, and guests assigned to multiple cages use the fixed cage-type priority.
   - `sqq-cage.gro` retains guest-to-cage memberships and every atom of a multi-atom guest.
   - `SQQ graph: <mode>` is printed when the script loads and thereafter only when the effective graph mode changes.

3. Runtime `config.yaml` and compact summaries
   - Replaced mandatory output-root `run_config.yaml` with `config.yaml`.
   - The file records final SQQ version, mode/engine, input metadata, requested/effective graph modes, requested/resolved workers, normalized output types, status/failures, and summary-write metadata.
   - Runtime configuration is written atomically for running, successful, and failed states.
   - Removed the detailed `config` worksheet from `summary.xlsx` and `config.csv` from `summary/`.
   - Retained a compact Configuration block on the summary dashboard.

4. Documentation and acknowledgements
   - Updated README and design documentation for the four modes, exact output presets, Tcl cage/guest commands, and runtime configuration ownership.
   - Replaced drive-specific Tcl examples with `source {path/to/sqq-render.vmd.tcl}`.
   - Added the approved README Acknowledgements section. Contributors are listed alphabetically by family name; ordering does not indicate relative contribution.

### Compatibility

- `-m 50` is rejected; use `-m py` or omit `-m`.
- Scripts generated before 0.3.5 retain their embedded commands. Regenerate the Tcl/GRO pair to use family-qualified cage/guest controls and guest memberships.
- Bare commands such as `sqq show 512` and `sqq color 512 blue` are not aliases; use `sqq show cage 512` and `sqq color cage 512 blue`.
- Consumers of `run_config.yaml`, the XLSX `config` worksheet, or `summary/config.csv` must use output-root `config.yaml` instead.
- Ordinary and cluster GRO output remains available only through explicit `--output-type` selection.
- Mode renaming, output routing, Tcl controls, and report consolidation do not change per-frame graph, ring, cage, occupancy, cluster, or order-parameter values.

## Version 0.3.4

### Short Summary

Version 0.3.4 adds deterministic topology grouping for multiple GRO inputs, separates incompatible systems into independent result roots, changes the default main-summary CSV directory to `summary/`, and adds an information-only safety fallback above 26 topology groups. The grouped scheduler still uses one worker pool, resolves `auto` graph mode once per group, and validates a shared GRO topology against every input. It also makes `sqq show all` the generated VMD renderer's default all-cage command and gives overlapping cage types a deterministic display priority. Per-frame scientific algorithms and values are unchanged; the multi-GRO aggregation and output paths are intentionally different. Package and native-core metadata are synchronized at `0.3.4`, released Jul 23, 2026.

### Main Changes

1. Version metadata
   - Updated Python package, native C++ core, CMake project, and wheel-publish checks to `0.3.4`.
   - Updated root help and `sqq -v` / `sqq --version` to report `Release date: Jul 23, 2026`.

2. Multiple-GRO topology grouping
   - Added a pre-scan for invocations containing two or more GRO files.
   - The fingerprint contains atom count and ordered contiguous residue blocks, with each block represented by `resname` and ordered `atomname` values.
   - Titles/times, coordinates, velocities, boxes, and numeric atom/residue IDs do not affect the fingerprint; atom or residue-block order and identity do.
   - Groups are assigned in first-occurrence order. One topology uses the requested root; 2-26 topologies use independent `result_A` through `result_Z` roots.

3. Group-aware output layout
   - Each normal group owns its own run config, selected summary, `info/`, optional `gro/`, and selected annotated cage GRO/VMD renderer.
   - The requested root keeps the batch `run_config.yaml` source-to-group manifest when multiple groups exist.
   - `summary-csv` now defaults to `summary/`, with one UTF-8-SIG CSV per logical summary table. `summary.xlsx` and `summary/` may coexist when both main-summary types are selected.
   - More than 26 groups activate a whole-run information-only fallback. All readable GRO files are analyzed, but only the root run config and `info/*_info.md` are retained; no partial lettered groups, summaries, detail files, GRO, or renderer are produced.

4. Group scheduling, graph resolution, and topology validation
   - All GRO topology groups share one worker pool and bounded submission queue instead of running one group at a time.
   - Tasks carry both a global progress index and group-local frame index, preserving source order while keeping each group summary and annotated bundle internally ordered.
   - Requested `auto` graph mode is preserved in configuration but resolves once to `hbond` or `oo` for each topology group; SQQ-Py and SQQ-CPP consume the same resolved group setting.
   - A shared GRO `--top` is checked against every GRO input before analysis. A mismatch fails with the exact incompatible source rather than silently applying one topology to heterogeneous inputs.

5. VMD all-cage command
   - Replaced `sqq show cage` with `sqq show all` and made it the generated script's startup view; the removed spelling has no `show` alias.
   - Kept `sqq show phase`, `sqq show cluster`, and `sqq show domain` as whole-category selectors.
   - Added `all` as the cage-category color target, for example `sqq color all default`.

6. Deterministic shared-edge rendering
   - Separated cage topology priority from color choice and command argument order.
   - Nonstandard cages form the lowest layers, followed by `512 < 51262 < 51263 < 51264 < 435663 < 51268`; explicitly selected or recolored cage IDs form the final highlight layers.
   - A single layer retains the 0.125 angstrom radius. Multi-layer views distribute topology tiers over 0.125-0.130 angstrom so higher-priority coincident edges remain visible without unbounded line growth.
   - Generic cage IDs recover their canonical cage type when applying type-level color overrides.

7. Validation
   - Validated topology-fingerprint semantics, stable first-occurrence grouping, the exact 26-group A-Z boundary, the 27-group information-only fallback, and non-strict malformed-input scanning.
   - Retained focused Tcl tests for `show all`, rejection of `show cage`, reversed argument order, fixed standard-cage priority, bounded radii, exact-ID highlights, and generic-cage aliases.
   - Verified clean reuse transitions between multi-group, single-group, and information-only layouts; analysis or summary-write failures finalize root and group manifests as `failed`.

### Compatibility

- Single GRO files and trajectory inputs retain their established analysis/output paths.
- Multiple GRO files that previously shared one aggregate output are now combined only when their topology fingerprints match; incompatible systems receive separate result roots.
- The default `summary-csv` directory changes from `summary_csv/` to `summary/`; the public output-type name remains `summary-csv`.
- More than 26 topologies intentionally suppress every optional output except per-frame info, even when a mode or explicit request selected other files.
- Regenerate `sqq-render.vmd.tcl` to use `sqq show all` and deterministic cage layering.
- Grouping and routing do not change graph, ring, cage, cluster, occupancy, or order-parameter values for an individual frame.

## Version 0.3.3

### Short Summary

Version 0.3.3 adds `-t` as the short topology option, makes `input.lammps.type_map` optional for unambiguous standard water/methane DATA topologies, and replaces the generated VMD renderer's category commands with compact object-based `show` and `color` commands. Explicit LAMMPS mappings remain authoritative, and both SQQ-Py and SQQ-CPP consume the same resolved mapping. Package and native-core version metadata are synchronized at `0.3.3`, released Jul 22, 2026.

### Main Changes

1. Automatic LAMMPS atom-type mapping
   - When `input.lammps.type_map` is absent or empty, SQQ derives H/O/C candidates from DATA masses and type comments, then validates molecule graphs from DATA Bonds.
   - Automatic water recognition requires exactly one O plus two H atoms and two O-H bonds; automatic all-atom methane recognition requires exactly one C plus four H atoms and four C-H bonds. A clearly carbon-labeled unbonded singleton is accepted as a united-atom methane guest.
   - Valid DATA molecule IDs are retained. Invalid shared IDs are rebuilt deterministically only when all Bonds components have unique supported identities; the fallback is reported explicitly.
   - One numeric atom type must have one consistent inferred role. Ambiguous masses, conflicting roles, unsupported molecule compositions, or insufficient topology evidence fail before analysis and request an explicit map.
   - A non-empty user `type_map` always takes priority and retains the existing strict complete-map validation.
   - The inferred map and provenance are written to `run_config.yaml`, per-frame info, and main-summary configuration. The inference runs in the shared Python LAMMPS adapter, so SQQ-Py and SQQ-CPP receive identical normalized atoms and frames.
   - The default `guest.center_atoms` now maps `MET` to atom name `C`. Automatically inferred single-site and all-atom LAMMPS methane therefore use the carbon atom as the guest center; multi-atom methane no longer falls back to a whole-residue centroid for cage occupancy.

2. Topology CLI alias
   - Added `-t` as the short form of `--top` / `--topology`.
   - Unambiguous standard LAMMPS input can now run directly as `sqq analyze -i traj.lammpstrj -t system.data -m 99` without `-c`.

3. Version metadata
   - Updated the Python package and native-core version to `0.3.3`.
   - Updated root help and `sqq -v` / `sqq --version` to report `Release date: Jul 22, 2026`.
   - Updated CMake and wheel-publish native-version checks to require `0.3.3`.

4. VMD object commands
   - Replaced the former `sqq cage`, `sqq phase`, `sqq cluster`, and `sqq domain` commands in newly generated render scripts with `sqq show <object...>`; `sqq show <category>` shows an entire category.
   - Added automatic object recognition for registered cage labels, delimiter-free generic-cage aliases such as `4151062`, full frame-local cage IDs, phase labels, cluster IDs, and domain IDs. Multiple explicit selections are accepted within one family.
   - Added `sqq color <object> <color>`, accepting a VMD color name, ColorID, or `default`. Overrides do not depend on the current selection and remain active across show and frame changes in the current VMD session.
   - Reconstructed `51262_00053`, `cluster_00001`, and `domain_00001` names from compact GRO annotations. Such numeric IDs remain frame-local rather than tracked physical identities.
   - Cage objects are grouped by effective ColorID during rendering, so a single-cage override remains visible in an all-cage or cage-type view without creating one representation per cage.
   - Set the rendered DynamicBonds cylinder radius to 0.125 angstrom, corresponding to a 0.250 angstrom edge diameter.
   - Removed the duplicate cage-type atom-membership lists from Tcl state and deduplicated phase/cluster/domain atom lists per frame, reducing renderer memory without changing selected atoms.
   - Added a trajectory-wide object registry so misspelled cage, cage-ID, cluster-ID, and domain-ID targets fail instead of silently selecting nothing.
   - Coalesced rapid frame callbacks, preserved user-created VMD representations by tracking only SQQ-owned stable representation names, and reset pending/color state when a script is sourced again.
   - Ordered color layers by override specificity so exact cage highlights are drawn after default, category, and cage-type layers.
   - Kept the annotation format and scientific cage/cluster results unchanged.
   - Focused Tcl tests cover object parsing, category and multi-object selection, mixed-family rejection, ID reconstruction, named/numeric/default colors, and override precedence.

### Compatibility

- Existing explicit LAMMPS type maps and the long `--top` / `--topology` options retain their behavior.
- Automatic inference changes input setup only; a successfully inferred map represents the same atoms as the equivalent explicit map and does not alter cage algorithms, mode presets, or output defaults.
- Non-water/methane, ambiguous, or incomplete DATA topologies still require an explicit YAML mapping.
- Newly generated `sqq-render.vmd.tcl` files use the new object command interface. Automation written for the 0.3.2 category commands must use `sqq show cage`, `sqq show phase`, `sqq show cluster`, or `sqq show domain` instead.
- The Tcl command change affects visualization control only; GRO annotations, cage detection, classification, counts, and other scientific outputs are unchanged.

## Version 0.3.2

### Short Summary

Version 0.3.2 reorganizes the four analysis presets, adds strict orthorhombic LAMMPS trajectory input, adds a run-level annotated cage trajectory with a reusable VMD renderer, and strengthens cage acceptance in both engines. Modes `00` and `50` use SQQ-Py; modes `99` and `cpp` use SQQ-CPP. Mode `09` is removed. Requested `auto` graph selection is now reported together with the graph actually used.

The default command remains mode `50`. Mode-specific outputs are now explicit, and mode `cpp` does not generate classified `cage-gro` unless the user selects it. Package and native-core versions are `0.3.2`, released Jul 20, 2026.

### Main Changes

1. Modes and output presets
   - Mode `00`: SQQ-Py, `hbond`, 4/5/6 rings, 100% worker policy, cluster search on, and `info,gro,sqq-cage-gro,sqq-render,summary-xlsx,cluster-gro`.
   - Mode `50`: SQQ-Py, `auto`, 4/5/6 rings, 50% worker policy, cluster search off, and `info,sqq-cage-gro,sqq-render,summary-xlsx`.
   - Mode `99`: SQQ-CPP, `hbond`, internal 4/5/6 rings, 100% worker policy, no cluster support, and `info,gro,sqq-cage-gro,sqq-render,summary-csv`.
   - Mode `cpp`: SQQ-CPP, `auto`, internal 4/5/6 rings, 50% worker policy, no cluster support, and `info,sqq-cage-gro,sqq-render,summary-csv`.
   - Automatic percentages continue to reserve one physical core. Explicit `-w` / `--worker` overrides the preset.
   - `sqq-render` implies `sqq-cage-gro`. Enabling cluster search in mode `50` populates selected info/main-summary outputs but does not implicitly add split `cluster-gro`.

2. Orthorhombic LAMMPS input
   - Added LAMMPS DATA plus `.dump` / `.lammpstrj` and LAMMPS `.dcd` input through one Python adapter shared by SQQ-Py and SQQ-CPP.
   - Added `input.trajectory_stride` / `--trajectory-stride`; legacy `input.xtc_stride` is migrated during configuration loading.
   - Added physical `real`, `metal`, and `nano` unit conversion to the internal nm/ps representation, explicit atom-style handling, and a required YAML atom-type map.
   - DATA molecule IDs become SQQ residue IDs, DATA atom IDs become atom IDs, and every frame is reordered into stable DATA atom-ID order.
   - Fully periodic orthorhombic boxes are supported, including changing NPT lengths. Tilted/triclinic boxes, mixed or nonperiodic dump boundaries, `units lj`, missing molecule IDs, incomplete type maps, and topology/trajectory ID mismatches fail before analysis.
   - Selected LAMMPS frames use the existing process scheduler; each worker owns and closes a private MDAnalysis Universe.

3. Annotated cage trajectory and VMD renderer
   - Added output types `sqq-cage-gro` and `sqq-render`; the visible files are `sqq-cage.gro` and `sqq-render.vmd.tcl` in the run root.
   - `sqq-cage.gro` concatenates complete GRO blocks in successful frame order. It preserves source atom order, coordinates, box, and optional GRO velocities without PBC movement.
   - Compact ASCII annotations begin at column 69 after the optional velocity fields and encode all cage memberships plus optional phase, domain, and cluster assignments. Nonmembers use `m=-`.
   - The VMD script temporarily splits the run-level GRO for loading, removes its temporary files, and renders cage topology by default. Its public command is `sqq cage|phase|cluster|domain|help`.
   - Split per-frame structure GRO remains controlled by `gro` / individual structure output types. `cluster-gro` remains a separate cluster-search-dependent output.

4. Reporting and robustness
   - Terminal completion, `*_info.md`, summary configuration, run metadata, and annotated GRO titles report resolved graph selection as `auto -> hbond`, `auto -> oo`, or a counted mixed result across frames.
   - Terminal, per-frame info, and summary configuration now include normalized input format, trajectory stride, topology, and applicable LAMMPS provenance.
   - GRO input now preserves optional velocity fields and ignores SQQ semicolon annotations after the fixed-width record.
   - LAMMPS selection and per-frame molecule inventories use explicit DATA molecule IDs, so interleaved atom rows no longer split one molecule into several reported molecules.
   - `auto` now chooses hydrogen bonds only when every selected water has usable hydrogen coordinates; otherwise both engines choose O-O connectivity.
   - Fixed the mixed-graph display return path and restored worker-local XTC/TRR Universe initialization after adding the LAMMPS dispatch.
   - Interrupted or failed bundle generation removes temporary fragments and partial visible files.

5. Mandatory cage topology validation
   - SQQ-Py and SQQ-CPP now always require every candidate shell to use each edge exactly twice, satisfy `V - E + F = 2`, form one edge-connected face shell, have one cyclic face link around every vertex, and contain only trivalent shell vertices.
   - These checks run before cage type/isomer assignment and reject disconnected, pinched, branched, and non-manifold false cages even when `cage.scientific_validation` is `false`.
   - `--cage-scientific-validation on|off` is retained. It now controls only the additional PBC-aware face-planarity RMS, edge-length CV, nonzero projected-area, positive-volume, and volume-centroid geometry path; its default remains `off`.
   - Corrected false positives can reduce cage/isomer counts and their downstream occupancy, hydrate-cluster, GRO, info, and summary records. Ring and half/quasi detection are unchanged; rejected cages can return their rings/patches to ownership-filtered free outputs.

6. Verification
   - Configuration tests cover every mode's engine, graph, ring sizes, worker fraction, cluster state, exact default output order, `sqq-render` implication, and `cpp` default exclusion of `cage-gro`.
   - Synthetic LAMMPS DATA/dump tests cover shuffled atom rows, stable atom mapping, nm/ps conversion, changing boxes, frame stride, and both serial and worker adapters.
   - Annotated-GRO tests cover complete multi-frame blocks, column-69 comments, velocities, nonmembers, membership validation, temporary VMD loading, and resolved graph titles.
   - Real `tests/100.gro` smoke runs cover modes `00`, `50`, `99`, and `cpp`, including default-output presence/absence and explicit CPP `cage-gro`; native version/build checks report `0.3.2`. XTC worker initialization and legacy stride migration have dedicated regression checks.

### Compatibility

- Ring, half/quasi, F3/F4, and ice definitions are unchanged. Mandatory topology validation can remove cage false positives and consequently reduce cage/isomer, occupancy, hydrate-cluster, GRO, info, and summary results. Complete-hydrogen and oxygen-only inputs keep their prior `auto` result; partially mapped water topologies intentionally change from a mixed-quality hydrogen-bond graph to O-O connectivity.
- LAMMPS normalization can change only representation units/order before analysis; an equivalent physical frame is expected to produce the same scientific result within existing floating-point tolerance.
- Mode preset changes and default output changes are intentional. Scripts using removed mode `09` must choose `00`, `50`, `99`, or `cpp`.
- `sqq-cage.gro` requires identical atom identity and order across its successful frames. It is a visualization trajectory, while ordinary classified GRO output remains independently selectable.

## Version 0.3.1

### Short Summary

Version 0.3.1 adds `sqq analyze -m cpp`, a focused C++17 backend for the performance-critical cage workflow. Python remains responsible for the CLI, input/topology readers, configuration, process scheduling, and output writers; the native extension performs graph construction, internal chordless 4/5/6-ring search, generic cage topology and cage-isomer detection, automatic guest occupancy, and F3/F4.

The native mode intentionally does not reproduce the complete SQQ-Py feature surface. Its reports contain only connection diagnostics, cage topology/isomers, occupancy status, selected F3/F4, and the corresponding cage GRO, Markdown, per-table summary CSV, and optional XLSX output. Explicit unsupported requests fail before analysis, and the program never silently falls back to Python. Numeric modes `00`, `09`, `50`, and `99` retain the complete Python pipeline and existing scientific behavior.

The package and native-core version are `0.3.1`, released Jul 19, 2026. Release automation is configured to build platform wheels and a source distribution; this entry does not claim that the release has already been uploaded to PyPI.

### Main Changes

1. Hybrid Python/C++ architecture
   - Added the C++17 source tree under `sqq/core/sqq-cpp/` and the importable native module `sqq.core._sqq_cpp`.
   - Added a Python adapter that converts normalized SQQ frames, waters, guests, pair maps, boxes, and supported configuration into the native data contract, then reconstructs the existing graph/ring/cage/F3F4 result models.
   - The native call releases the Python GIL while one frame is analyzed. Existing process scheduling remains available for independent files and selected trajectory frames.
   - Native-module import or analysis failure is a hard error in mode `cpp`; no automatic SQQ-Py fallback is permitted.

2. Native scientific scope
   - Implemented `auto`, `hbond`, `oo`, and user `pairs` water graphs with deterministic edge ordering and orthorhombic minimum images.
   - Implemented internal canonical chordless-ring search for any nonempty subset of sizes 4, 5, and 6. Rings are retained for cage construction but are not public SQQ-CPP output.
   - Implemented generic Euler-compatible cage growth, report-group filtering, the existing cage labels, optional shell/manifold/geometry scientific validation, and the established hexagonal-face cage-isomer definition.
   - Implemented polyhedron occupancy automatically when selected guest molecules exist. Reports distinguish an analysis with no selected guests (`not evaluated`) from an evaluated cage population with zero occupancy.
   - Implemented independently selectable F3 and F4. F4 without usable water-hydrogen coordinates is unavailable with a warning rather than zero.

3. Mode defaults and supported controls
   - Added `cpp` to `-m` / `--mode`; the command default remains mode `50`.
   - SQQ-CPP defaults are `auto` graph, internal 4/5/6 rings, `f3,f4`, approximately 90% of detected physical cores with one core reserved, and `info,cage-gro,summary-csv` output. XLSX is no longer part of the native default.
   - Compatible controls include input/topology/output/configuration, graph/pair settings, `-s` within 4/5/6, cage report groups and face limit, scientific cage validation, F3/F4 selection, worker count/fraction, process or serial execution, strict/input controls, grouped/flat supported output selection, and `--cage-isomer-rows nonzero|all` for native `summary_csv/cage_isomer.csv` and an optional XLSX isomer sheet.
   - In mode `cpp`, `--order-parameter all` means `f3,f4`; `--output-type gro` means `cage-gro`; `--output-type all` means `info,cage-gro,summary-csv,summary-xlsx`.
   - Added `summary-csv`, which writes each applicable main-summary table as a separate UTF-8-SIG file under `summary_csv/`. It is the CSV equivalent of workbook sheets, not detail output. SQQ-Py supports it but leaves it off by default; SQQ-CPP enables it by default.
   - Added `output.summary_csv_dir` with default `summary_csv`, alongside the retained `output.summary_detail_dir` default `summary_detail`. The two values must be different relative paths inside the output root.
   - Stale CSV cleanup is restricted to known SQQ-generated filenames inside the currently configured summary directories; unknown files and formerly configured custom directories are not removed.
   - Renamed `xlsx` to `summary-xlsx` and `summary-detail` to `summary-detail-csv`, without compatibility aliases. SQQ-Py now defaults to `info,gro,summary-xlsx`, so ordinary detail CSV is no longer written unless selected explicitly. SQQ-CPP does not support `summary-detail-csv`.
   - Resolved cluster search always forces `cluster-gro`. It adds `summary-xlsx` only when neither `summary-csv` nor `summary-xlsx` is selected; an existing `summary-csv` selection remains CSV-only.
   - `--output-type none` disables optional native reports and cage GRO files, but the mandatory `run_config.yaml` is still written.

4. Explicitly unsupported scope
   - SQQ-CPP does not expose public ring results, ring GRO, ring size 7, or shortest-path rings.
   - It does not calculate or output half-cages, quasi-cages, hydrate clusters, ice, Q_l, MCG, or DHOP.
   - VMD, membership/order TSV, `summary-detail-csv`, `cluster-detail`, and other detailed exports remain SQQ-Py features.
   - Thread backend, Python fast closure, and triclinic boxes are unsupported.
   - Explicit incompatible CLI options and nondefault incompatible configuration values are collected into clear validation errors instead of being ignored.

5. Compact native reports and main summaries
   - Terminal, main summary, and per-frame Markdown use `SQQ version` followed by `Mode`; numeric modes display as `NN (sqq-py)` and native mode displays as `sqq-cpp`.
   - SQQ-CPP progress omits Python-only half/quasi, cluster, filtering, and ice stages.
   - Per-frame `Frame Information` begins `sqq version`, `mode`, `date & time`, `source`, `frame`, and `time_ps` and omits Python-only analysis sections.
   - The compact native main-summary mapping contains `summary`, `cage`, `cage_isomer`, selected F3/F4 `order_parameter`, and `config`; `failures` is conditional, and `cage_occupancy` is present only when selected guests exist. Default `summary-csv` writes one file per table, while explicit `summary-xlsx` writes the same mapping as workbook sheets.
   - Default structure output is cage GRO only. The Python reader/writer layer preserves standard GRO records and the source box; the SQQ-CPP cage export does not add a synthetic cage-center atom.

6. Build and release packaging
   - Replaced the pure-Python build backend with scikit-build-core plus CMake/pybind11 so source builds compile the C++17 extension.
   - A local CPython 3.12 Windows x86_64 wheel was built, installed into an isolated target, and imported successfully. The wheel contains the native `.pyd` and omits the C++ source tree.
   - The source distribution contains the CMake, header, and C++ source files required for an explicit local build.
   - Release CI builds and tests precompiled wheels for CPython 3.10-3.14 on Windows x86_64, Linux x86_64, macOS x86_64, and macOS arm64, then builds a source distribution. Installing a matching wheel does not compile C++ on the user's machine; an explicit source build requires CMake 3.20 or newer and a C++17 toolchain.
   - Native build products, CMake/Ninja state, wheels, source archives, and local test outputs remain ignored; handwritten `.cpp`, `.hpp`, and `CMakeLists.txt` files remain tracked.

7. Verification
   - GCC 16, CMake, and Ninja completed a clean native build; the native module and installed wheel reported core version `0.3.1`.
   - Random parity covered 100 graph/ring cases, 40 random frames in each of `oo`, `hbond`, and `auto`, and 50 random F3/F4 frames. Graph edges and canonical rings matched exactly; F3/F4 matched within `2e-14`.
   - Cube-cage tests covered topology, scientific geometry validation, and occupancy against the Python reference.
   - The real `tests/100.gro` frame contains 11,104 atoms, 2,176 waters, and 384 guests. C++ and Python matched exactly for 4,322 graph edges; 2,499 rings (45 four-rings, 2,147 five-rings, and 307 six-rings); 339 cages across 16 types by type, water membership, and face rings; all 339 cage-isomer labels; and all 339 occupancy assignments, including 315 occupied cages.
   - All 2,176 F3 values matched exactly. All 2,176 F4 values passed tolerance with maximum absolute difference `4.44e-16`.
   - On that frame, the supported native core took about 0.5673 s versus 10.0967 s for the equivalent Python core, approximately 17.8 times faster. This is a core-path benchmark, not a guarantee for total runs, which also include input and output time.
   - The installed-wheel end-to-end run on the same frame completed in about 0.9 s and produced compact info, main-summary, and cage GRO outputs. Output-format tests cover default per-table CSV and explicit XLSX separately.
   - The current external regression suite passed 151 tests plus five unittest subtests under the 0.3.1 output contract.

### Compatibility

- Numeric SQQ-Py modes retain their graph, ring, patch, cage, cluster, occupancy, order-parameter, and ice behavior. The output selector names and defaults intentionally change as documented above. Selecting `-m cpp` is an explicit opt-in to a reduced backend and compact schema.
- There is no `-m c` alias. Scripts must use exactly `-m cpp`.
- Scientific definitions for the supported graph/ring/cage/isomer/occupancy/F3/F4 path are aligned with SQQ-Py. Native floating-point reductions are compared with tolerance where bitwise identity is not portable.
- Python-only options are not accepted as no-ops in mode `cpp`; explicit incompatible requests fail. Full configurations may retain unrelated default sections, but incompatible nondefault requests are rejected.
- SQQ-CPP supports orthorhombic or non-periodic boxes only and never falls back to the Python engine.
- Consumers must not expect ring, half/quasi, cluster, ice, hydrate-order, VMD, TSV, or `summary-detail-csv` records from SQQ-CPP. Use a numeric mode when those outputs are required.

## Version 0.2.10

### Short Summary

Per-frame reports now record the selected mode directly and expand the Ring table to show both total primitive rings and final free rings for every reported size. These are presentation-only changes.

Version 0.2.10 restores a compact hydrate-cluster hierarchy to each selected per-frame `*_info.md` report when the resolved cluster-search state is on. Mode `09` therefore reports cluster information by default, while `--find-cluster off` suppresses both the search and the section. The report summarizes unique cluster cages, sI/sII/sH domains, boundary cages, unresolved cages, and isolated cages without restoring the verbose pre-0.2.8 detail sections. Hydrate topology, phase classification, boundary assignment, and every other scientific result are unchanged.

### Main Changes

1. Resolved cluster state in per-frame info
   - `Frame Information` records `find_cluster` as `on` or `off`.
   - `Hydrate Cluster` is written only when cluster search resolves to on and `info` belongs to the effective output selection.
   - `-m 09` and `-m 00` include the section by default. `-m 09 --find-cluster off` omits it, while `-m 50 --find-cluster on` includes it.
   - Cluster search still forces `xlsx` and `cluster-gro`, but does not force `info`; `--find-cluster on --output-type none` therefore writes no `*_info.md` file.

2. Compact hierarchy
   - Reported clusters are listed sequentially as `cluster_00001`, `cluster_00002`, and so on.
   - Each cluster row reports its unique cage count. Child rows show sI/sII/sH domains, generic boundary, and compact unclassified topology; each child is further divided by cage type.
   - The compact `unclassified` row is the deduplicated unresolved set: stored `ambiguous_cage_ids` and `unclassified_cage_ids` plus any uncategorized residual cluster cages. This is display-only; the workbook and cluster-detail CSV retain their separate ambiguous and unclassified fields.
   - Zero-count categories are omitted. `isolated` appears once as the final top-level row and is not divided by cage type.
   - Child counts conserve their parent count, and cluster counts use unique cage IDs rather than water-coordinate unions.

3. Detail ownership
   - Exact cage IDs, seeds, confidence, domain adjacency, and other verbose records remain in `summary_detail/hydrate_domain.csv` and `summary_detail/hydrate_cluster_detail.csv` when `cluster-detail` is selected.
   - `summary.xlsx` retains the separate classified, boundary, ambiguous, and unclassified cluster fields and the existing plotting-oriented schema.

4. Mode and documentation alignment
   - Every `*_info.md` Frame Information table records the selected preset after `time_ps`, using the same two-digit code and label as the terminal and workbook, for example `09 (rigorous-performance)`.
   - The per-frame Ring table now reports `ring size | total | free` and sums both numeric columns. `total` is the primitive-ring population; `free` excludes rings owned by final half-cage, quasi-cage, or closed-cage output.
   - The README mode table explicitly records graph mode, ring sizes, automatic worker fraction, and initial cluster state for modes `00`, `09`, `50`, and `99`.
   - Cluster-state precedence is documented as `--find-cluster` over `hydrate_cluster.enabled` in `config.yaml` over the mode preset. Mode `50` remains the command default.

5. Package version and verification
   - Updated package and root version output to `0.2.10`, released Jul 19, 2026.
   - Regression coverage checks mode/config/CLI precedence, Info mode ordering, total/free ring conservation, conditional cluster-info creation, multiple clusters, sI/sII/sH domains, compact unresolved counts, zero omission, isolated placement, and cage-count conservation.

### Compatibility

- Ring, half-cage, quasi-cage, closed-cage, hydrate-domain, boundary, occupancy, F3/F4/Q_l, MCG/DHOP, and ice algorithms and values are unchanged.
- Existing XLSX, detail-CSV, and cluster-GRO scientific schemas are unchanged.
- Every generated `*_info.md` gains the mode row and expanded Ring table. The compact Hydrate Cluster hierarchy remains conditional on resolved cluster search being on; output selections without `info` still create no Markdown report.

## Version 0.2.9

### Short Summary

Version 0.2.9 replaces overlapping phase/boundary labels with one mutually exclusive hydrate-cluster partition. Final sI, sII, and sH domain cages remain classified phase cages; the generic boundary contains only non-domain cages in the first complete shared-face layer outside a final phase domain. Direct phase-phase contacts no longer relabel either endpoint, boundary membership does not propagate beyond the first external layer, and the former sI/sII/sH-boundary, transition, and interface-context categories are removed. Enabled cluster search now also forces native per-frame sI/sII/sH/boundary GRO views together with XLSX output. Phase seeds and domain expansion are unchanged, but boundary IDs, cluster-summary schemas, and cluster output files intentionally change from 0.2.8.

### Main Changes

1. Mutually exclusive cluster categories
   - Every cage in a reported cluster belongs to exactly one of `classified_cage_ids`, `boundary_cage_ids`, `ambiguous_cage_ids`, or `unclassified_cage_ids`.
   - `classified_cage_ids` remains the union of exclusive final sI, sII, and sH domains.
   - Boundary, ambiguous, and residual unclassified cages are selected only after the final domains are frozen.
   - The four category counts sum to the cluster's unique `cage_count`.

2. Generic external first-layer boundary
   - A boundary cage must be outside every final phase domain and share a complete cage face with at least one domain cage.
   - Only the non-domain cage is marked. Its contacted domain neighbors retain only their sI, sII, or sH identity.
   - A direct contact between different phase domains does not create a boundary cage by itself and does not relabel either endpoint.
   - Search stops at the first external non-phase layer; another non-domain graph step is not added automatically.
   - A non-domain cage touching multiple phase domains remains one generic boundary cage.

3. Removed overlap categories
   - Removed `sI-boundary`, `sII-boundary`, `sH-boundary`, phase-boundary context labels, and transition-cage membership.
   - Removed the old phase-specific, classified/unclassified/ambiguous-boundary, transition, and interface-context report fields.
   - Competing phase claims outside the boundary remain ambiguous; all other residual cluster cages remain unclassified.

4. Workbook and cluster-detail reporting
   - The `hydrate_cluster` workbook sheet now reports `classified_cage_count`, `boundary_cage_count`, `ambiguous_cage_count`, and `unclassified_cage_count`.
   - `hydrate_cluster_detail.csv` exposes the four corresponding cage-ID groups and adds `boundary_composition`.
   - `hydrate_domain.csv` names direct external adjacency as `external_boundary_contact_count` and `external_boundary_contact_ids`.
   - Domain boundary-contact records remain relationships; they do not assign another scientific category to a cage.

5. Native cluster GRO output
   - Added the search-dependent `cluster-gro` output type. Resolved `--find-cluster on` forces both `cluster-gro` and `xlsx`, even when the requested output selection is `none`.
   - Resolved cluster search off writes no category GRO files. Reusing an output directory removes stale SQQ-generated grouped `hydrate_cluster` directories and flat cluster filenames.
   - Grouped layout writes `<frame>/hydrate_cluster/<frame>_cluster_sI.gro`, `<frame>_cluster_sII.gro`, `<frame>_cluster_sH.gro`, and `<frame>_cluster_boundary.gro`; flat layout places the same filenames at the frame root.
   - Each file aggregates every domain/cluster of its category within the frame. Ambiguous, residual unclassified, isolated, and below-threshold cage IDs are omitted.
   - Export maps exclusive cage IDs to deduplicated complete water molecules only. Guests, CNT atoms, and other non-water molecules are excluded.
   - Missing categories are omitted unless `output.write_empty_files` is enabled.
   - Every exported atom retains its exact original wrapped frame coordinate and every file retains the original box. No category is moved, centered, unwrapped, or independently made whole.
   - Periodic/percolating networks may retain cross-box lines because one single-copy GRO cannot remove every periodic seam. Different category files may share face-water molecules even though their cage IDs are exclusive.

6. Package version and verification
   - Updated `pyproject.toml`, `sqq.__version__`, and root version output to `0.2.9`, released Jul 16, 2026.
   - Unit tests cover mutually exclusive classification, direct phase-phase contact, multi-phase contact through one non-domain cage, first-layer stopping, sH behavior, report schemas, category-count conservation, forced cluster output, grouped/flat paths, stale cleanup, molecule completeness, and exact PBC coordinate/box preservation.
   - The mixed sI/sII real-GRO regression contains one 334-cage main cluster: 260 classified cages (66 sI and 194 sII), 69 boundary cages, 0 ambiguous cages, and 5 residual unclassified cages. Five additional cages are isolated or in below-threshold components. All main-cluster category ID sets are disjoint.

### Compatibility

- Ring, half-cage, quasi-cage, closed-cage, occupancy, F3/F4/Q_l, MCG/DHOP, ice, phase-seed, and phase-domain expansion results are unchanged by this release.
- Boundary cage IDs and every report derived from the former overlapping boundary model can intentionally change.
- Code must no longer expect phase cages inside `boundary_cage_ids`, phase-specific boundary categories, `transition_cage_ids`, or phase-boundary context labels.
- Consumers of `summary.xlsx`, `hydrate_cluster_detail.csv`, or `hydrate_domain.csv` must use the new mutually exclusive count/ID fields and renamed external-contact fields.
- `cluster-gro` is now a canonical output type but requires cluster search. Search on forces it together with XLSX; search off suppresses it and cleans stale generated category files.
- Cluster GRO consumers receive full water molecules in original wrapped coordinates and the original box, not independently unwrapped or centered category structures.
- Neighboring cages may still share face-water coordinates in structure files; mutual exclusivity applies to detected cage IDs, not to the union of their water atoms.

## Version 0.2.8

### Short Summary

Version 0.2.8 separates hydrate phase identity from boundary membership. The first complete cage layer on both sides of a shared-face phase interface can now retain sI, sII, or sH identity while also being recorded as boundary. Boundary propagation stops after one cage-graph step, and summary output distinguishes unique boundary totals, phase-classified boundary cages, unclassified/ambiguous boundary cages, and all single-, two-, and three-phase interface contexts. This release also adds mode `09`, makes cluster-search defaults part of the four mode presets, replaces `--hydrate-cluster` with `--find-cluster`, and replaces negative output suppression with positive `--output-type` selection. The removed CLI and YAML output selectors have no compatibility aliases or migration. Ring, patch, cage, occupancy, order-parameter, and ice definitions are unchanged; hydrate boundary membership and related reports intentionally change from 0.2.7.

### Main Changes

1. Independent phase and boundary properties
   - Removed the rule that skipped every phase-domain cage during boundary resolution.
   - Removed the definition `boundary = component - classified`.
   - `classified_cage_ids` and `boundary_cage_ids` may now overlap, while each tuple remains internally unique.
   - Cluster `cage_count` continues to count unique cage ids and is not calculated by adding phase and boundary counts.

2. One-layer shared-face boundary
   - A direct shared-face contact between different phase domains marks the complete cage on both sides.
   - A non-domain cage that directly contacts a phase domain and every directly contacted domain cage enter the same boundary layer.
   - Contacted domain phases and retained phase claims produce stable labels in sI, sII, sH order.
   - Supported interface contexts are `sI`, `sII`, `sH`, `sI+sII`, `sI+sH`, `sII+sH`, and `sI+sII+sH`.
   - Expansion stops after the directly contacted cage layer; deeper unclassified cages are not treated as boundary automatically.

3. Domain and cluster semantics
   - Phase domains remain exclusive and use the established strict-seed and compatible-expansion logic.
   - `HydrateDomain.boundary_cage_ids` records external boundary contacts rather than boundary members inside that domain.
   - `transition_cage_ids` is the subset of boundary cages whose interface context contains more than one phase.
   - Non-domain cages retain unclassified or ambiguous phase status independently of boundary membership.

4. Workbook and cluster-detail reporting
   - Added classified, unclassified, and ambiguous boundary counts alongside the unique boundary total.
   - Phase-specific boundary counts now measure boundary cages that retain sI, sII, or sH domain identity.
   - Added `sI_boundary_context_cage_count`, `sII_boundary_context_cage_count`, and `sH_boundary_context_cage_count`; existing multi-phase count columns measure the corresponding multi-phase interface context.
   - Cluster detail records add `boundary_phase_counts` and the phase-specific/unclassified/ambiguous boundary breakdown.
   - Boundary totals remain unique even when a cage retains a phase identity and also belongs to the boundary.

5. Four analysis modes
   - Added mode `09`: `hbond`, 4/5/6-ring search, 90% automatic worker fraction, and cluster search on.
   - Mode `00` remains `hbond`, 4/5/6 rings, and 25% workers, and now defaults cluster search to on.
   - Modes `50` and `99` default cluster search to off. Their graph/ring/worker presets remain `auto`, 5/6, 50% and `oo`, 5/6, 90%, respectively.
   - Mode `50` remains the command default, so an unqualified run has cluster search off.

6. Unified cluster search
   - Replaced `--hydrate-cluster on/off` with `--find-cluster on/off`.
   - Explicit `--find-cluster` overrides `hydrate_cluster.enabled` from the configuration and the selected mode.
   - Enabling cluster search forces `xlsx` into the effective output set and guarantees the per-frame `hydrate_cluster` sheet in `summary.xlsx`.
   - Cluster, domain, and boundary data is no longer copied into per-frame `*_info.md` reports.
   - Output type `cluster-detail` writes `summary_detail/hydrate_domain.csv` and `summary_detail/hydrate_cluster_detail.csv` and is invalid when cluster search is off.

7. Positive output selection
   - Added `--output-type TYPE[,TYPE...]`, which replaces the complete output list. Its default is `info,gro,xlsx,summary-detail`.
   - Canonical types are `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `xlsx`, `summary-detail`, and `cluster-detail`.
   - `gro` expands to every GRO subtype; `all` selects every type applicable to enabled analyses; `none` selects no optional output. Mandatory `run_config.yaml` remains.
   - Replaced `output.disabled_outputs` with canonical `output.types`. Old YAML is rejected rather than migrated.
   - Removed `--cluster-detail`, `--no-output`, `--write-order-tsv`, and the hidden individual `--no-*` output switches. Together with `--hydrate-cluster`, these old names are no longer accepted.

8. Visualization rule
   - Scientific records retain dual phase/boundary identity.
   - An exclusive three-color rendering should give boundary display priority and subtract boundary cage ids from the sI/sII/sH display sets.
   - Hydrate clusters still do not generate separate GRO files. Coordinate unions cannot encode cage ownership because cages share waters and faces; a custom renderer should use detected cage/ring edges rather than infer cages from rebonded coordinates.

9. Package version and verification
   - Updated `pyproject.toml`, `sqq.__version__`, and root version output to `0.2.8`, released Jul 16, 2026.
   - Added focused boundary, four-mode preset, CLI-priority, positive-output, stale-output cleanup, XLSX-forcing, cluster-detail routing, and removed-argument tests.
   - Passed 132 tests plus five unittest subtests, including all retained and updated 0.2.7 regression coverage.
   - A strict real-GRO smoke run on `sqq_example.gro` completed with zero failures. Its 154 unique boundary cages partitioned into 107 classified plus 47 unclassified cages, and independently into 9 sI-context plus 91 sI+sII-context plus 54 sII-context cages.
   - Additional strict real-GRO runs confirmed that mode `09` with only `cluster-detail` selected forces XLSX and writes exactly the two cluster CSV files, while mode `50` with XLSX only omits both the cluster stage and `hydrate_cluster` sheet.

### Compatibility

- Ring, half-cage, quasi-cage, closed-cage, occupancy, F3/F4/Q_l, MCG/DHOP, and ice results are unchanged by this release.
- Hydrate domain membership uses the same phase claims, but boundary cage ids, transition counts, domain boundary contacts, cluster-detail records, and hydrate-cluster summary columns can intentionally change.
- Code that assumed `classified_cage_ids` and `boundary_cage_ids` were disjoint must switch to set-based unique totals.
- Code that assumed `cluster.cage_count == classified_cage_count + boundary_cage_count` is no longer valid.
- Mode `50` remains the default and therefore preserves the overall cluster-search default of off. Mode `00` now enables cluster search; mode `09` is new.
- Scripts and configurations must replace `--hydrate-cluster` with `--find-cluster`, replace output switches with `--output-type`, and replace `output.disabled_outputs` with `output.types`.
- The removed cluster/output CLI names and old YAML output selector are intentionally rejected without compatibility aliases or migration.
- Enabling cluster search always creates `summary.xlsx`; select `cluster-detail` only when domain and per-cluster CSV records are required.

## Version 0.2.7

### Short Summary

Version 0.2.7 unifies order-parameter selection under `--order-parameter NAME[,NAME...]` and output suppression under `--no-output TYPE[,TYPE...]`. The default descriptor set is now F3/F4 only, while the default disabled-output set is `none`. The release also fixes PBC/input/configuration validation, repeated residue-ID handling, guest occupancy centroids, serial failure handling, hydrate-cluster report-scope coupling, summary-write robustness, and several equivalent hot paths. Ring, cage, and order-parameter definitions are unchanged; corrected guest centroids and cluster scope can intentionally change their corresponding outputs.

### Main Changes

1. Unified order-parameter selector
   - Added `--order-parameter` with comma-separated names `f3`, `f4`, `qN`, `mcg1`, `mcg3`, `dhop35`, and `dhop30`.
   - Added `all`, which expands to `f3,f4,q6,q12,mcg1,mcg3,dhop35,dhop30`.
   - Added `none`, which skips all order-parameter calculations and omits the `order_parameter` workbook sheet.
   - Names are case-normalized, deduplicated, and stored in deterministic order.
   - The explicit CLI selector overrides `order.parameters`; the built-in default is now `[f3, f4]`.

2. Independent calculation and dynamic output
   - F3 and F4 now have independent calculation switches instead of one shared `f3f4_enabled` switch.
   - Each `qN` name directly selects that non-negative Q_l degree; Q_l neighbor mode, cutoff, and fixed-neighbor settings remain separate controls.
   - MCG/DHOP numerical settings remain under `hydrate_order`, while selection belongs only to `order.parameters`.
   - The terminal header/final summary, dashboard, per-frame Markdown, and `run_config.yaml` report the normalized selection.
   - `summary.xlsx/order_parameter` contains only selected columns. Selected but inapplicable MCG values remain `N/A`.
   - Focus mean/count columns are omitted unless `order.focus_waters` contains at least one requested water residue.
   - `--write-order-tsv` writes only selected per-water F3/F4/Q_l columns. MCG/DHOP remain frame-level descriptors, so a hydrate-only selection does not create an empty TSV.

3. Compatibility migration
   - Hidden compatibility options `--no-q`, `-q` / `--q-degree`, `--mcg3`, and `--dhop30` remain accepted in 0.2.7 and emit a deprecation warning.
   - When `--order-parameter` is supplied together with a legacy selector, the unified selector wins.
   - Legacy YAML enable keys are translated when `order.parameters` is absent, allowing older generated configuration files to retain their requested descriptor set.
   - New default configuration files use `order.parameters` and no longer emit the old enable booleans.

4. Workbook cleanup
   - Removed the redundant one-row-per-frame overview sheet named `frame`.
   - Dedicated connection, ring, half-cage, quasi-cage, cage, hydrate-cluster, order-parameter, and ice sheets continue to carry their per-frame data.
   - The dashboard, detail CSV files, `detail_index`, and config sheet are unchanged by this removal.

5. Unified output suppression
   - Added `--no-output` with comma-separated names `info`, `membership-tsv`, `order-tsv`, `vmd`, `gro`, `ring-gro`, `half-gro`, `quasi-gro`, `cage-gro`, `ice-gro`, `xlsx`, and `summary-detail`.
   - The default is `none`, which preserves the established default files. `all` disables every optional output but always retains `run_config.yaml`.
   - Added `output.disabled_outputs: []` as the canonical YAML setting. Explicit `--no-output` replaces the configured list.
   - `gro` supersedes all GRO subtypes. Aliases normalize to canonical names, and repeated names are removed deterministically.
   - Individual `--no-info`, `--no-*-gro`, `--no-xlsx`, and `--no-summary-detail` options remain hidden compatibility aliases for one release and emit a warning.
   - `order-tsv` suppression overrides `--write-order-tsv`; disabling `summary-detail` suppresses optional cluster detail CSVs.
   - Reused output directories remove only known stale files for disabled categories. Empty per-frame directories are removed when every per-frame output is disabled.
   - Terminal, dashboard, and `run_config.yaml` report the normalized disabled-output set.

6. Correctness, reliability, and equivalent performance fixes
   - Fixed `minimum_image()` for zero or partially nonpositive box axes; invalid axes are treated as non-periodic instead of producing NaN.
   - GRO now validates its declared atom count, finite coordinates, mandatory 3/9-value box line, and single-frame structure. Extra non-empty records after the box are rejected. All-zero boxes become non-periodic; triclinic tilt terms are rejected.
   - XTC/TRR frames with non-finite coordinates or non-90-degree angles are rejected explicitly. Invalid trajectory time is recorded as unavailable. XYZ coordinate scaling is configurable through `input.xyz_scale` / `--xyz-scale` with default `0.1`.
   - Occupancy and MCG now share one PBC-aware guest centroid. Guest centers are precomputed once per frame, eliminating repeated cage×guest centroid calculation.
   - Degenerate solid-angle triangles contribute zero, and disconnected topology objects fail explicitly instead of mixing wrapped and unwrapped coordinates.
   - Non-strict standalone serial read errors become failed rows and processing continues. Strict analysis and summary failures update mandatory `run_config.yaml` to `status: failed`.
   - Failed rows are available through `run.failures`, an optional `failures` workbook sheet, and `summary_detail/failures.csv`.
   - Configuration is normalized once before thread dispatch; queued sibling tasks are cancelled on strict thread failure, and parent-side trajectory progress no longer depends only on queue timing.
   - Legacy selector keys are removed after migration, unrelated empty output directories are preserved, and hydrate-only order TSV requests emit a warning.
   - Single-task runs skip physical-core probing, immutable ring size/edge properties are cached, and non-graph Q_l modes reuse deterministic cell-list cutoff pairs.
   - Explicit worker values are now validated even when a run has only one task or uses the serial backend; automatic single-task resolution still avoids probing physical cores.
   - XYZ now validates a nonnegative declared count, exactly that many finite coordinate records, and no extra nonempty records. Multi-frame XYZ must be split before analysis.
   - The common full-orthorhombic `minimum_image()` path is vectorized without boolean advanced-index temporaries. Large cutoff searches can use MDAnalysis compiled candidates, followed by the established exact float64 PBC recheck; fallback cell-list behavior is retained.
   - MDAnalysis trajectory readers cache immutable atom metadata per Universe. Thread scheduling now uses the same bounded `3 * workers` in-flight window as process scheduling.
   - Occupancy now uses one reusable PBC-aware guest-center cell index per frame, exact distance rechecks, batched solid-angle evaluation, and a scalar fallback at the membership boundary. Guest order and non-exclusive overlapping-cage semantics are retained.
   - Cage-grow DFS keeps face counts as a fixed `(n4, n5, n6)` tuple and precomputes target face incidences; target masks, candidate order, closure checks, and pruning semantics are unchanged.
   - DHOP compares valid plane-normal matrices in batches per central O-O bond, with scalar rechecks near 30/35-degree thresholds. F4 reuses cached O-H and graph O-O minimum-image vectors for its dihedral geometry.
   - Water and guest selection now groups contiguous residue blocks in source order, preventing wrapped/repeated five-digit GRO residue IDs from merging distinct molecules. Whitespace GRO fallback parsing also preserves digit-containing residue names such as `TIP3`.
   - Configuration normalization now parses textual booleans explicitly and rejects unsupported modes, non-finite or nonpositive cutoffs/scales, and boolean/fractional/nonpositive integer settings before dispatch.
   - `input.xtc_stride` is validated as a positive integer instead of silently changing zero to one.
   - Hydrate-cluster topology and its detail/domain lookups now use all detected cages. `cage.report_types` / `--cage-size` remains a reporting filter and no longer changes connectivity or phase evidence.

7. Summary-write observability and robust output
   - `write_summary()` now records per-table rows, columns, cells, bytes, CSV/XLSX write time, workbook-format time, final-save time, and total time under `run.summary_write` in mandatory `run_config.yaml`.
   - The final terminal Run Summary adds `Summary write (s)`.
   - `summary.xlsx`, each detail CSV, and `run_config.yaml` are written through same-directory temporary files and atomically replaced only after a successful write. Existing final files remain intact if a new write fails before replacement.
   - Exact quasi-cage isomers are retained as sparse per-frame records until `quasi_cage_isomer.csv` is written; they no longer create a DataFrame column per observed isomer. Composition-level quasi counts and all scientific count values are unchanged.
   - Related detail CSV replacements/removals now commit as one recoverable bundle. Each workbook table is preflight-checked against Excel's 1,048,576-row and 16,384-column limits, so an actionable SQQ error replaces a late pandas traceback for unexpected oversized compact sheets.
   - Data sheets above 200,000 cells or 128 columns retain header style, filter, freeze pane, and fixed widths but skip per-body-cell style/auto-width work. Small sheets retain the established full formatting. No result value or table schema is changed.


8. Verification
   - Added and passed 75 focused local 0.2.7 tests plus five unittest subtests covering the unified selectors, root version flags, PBC, guest centroid, GRO/trajectory/XYZ validation, residue grouping, configuration normalization, hydrate-cluster report-scope independence, strict/non-strict failures, failure artifacts, output cleanup, performance equivalence, and end-to-end CLI paths.
   - The existing four Q_l reference tests and eight cage report-scope tests also pass.
   - A strict real-GRO smoke run with `--order-parameter all` completed with one input and zero failed frames; the workbook contained every selected descriptor and no `frame` sheet.
   - A strict real-GRO run with `--no-output quasi-gro,cage-gro,xlsx` completed with zero failures; info, ring/half/ice GRO, detail CSVs, and `run_config.yaml` remained, while the disabled outputs and their empty directories were absent.
   - The same input was analyzed with an isolated PyPI 0.2.6 baseline. The hbond, ring, half_cage, quasi_cage, cage, and ice sheets matched exactly, and every common F3/F4/Q6/Q12/MCG/DHOP output value matched.
   - After the PBC guest-centroid fix, the real-GRO hbond, ring, half_cage, quasi_cage, cage, order_parameter, and ice sheets still match the pre-fix run exactly. Occupied cages intentionally changed from 251 to 275, and spurious multi-guest compositions disappeared for this boundary-crossing sample.
   - The same all-parameter real run fell from about 27.0 s to 20.8 s despite also writing detail output, primarily from precomputed guest centers and skipped single-task CPU probing.
   - Two-file spawned-process and compatibility-thread smoke runs both completed with `frames_total=2`, `frames_ok=2`, `frames_failed=0`, and final `status=completed`.
   - Python compilation and `git diff --check` pass. Both DOCX files reopen through `python-docx` and pass ZIP/OOXML integrity and required-text checks.
   - Added and passed seven focused performance/output tests covering PBC spatial-index equivalence, scalar versus batched occupancy membership, scalar versus batched DHOP counts, cached versus scalar F4, fixed-count cage pruning, lightweight worksheet formatting, and summary timing/temporary-file cleanup.
   - A real all-parameter `1to2_small_100ns.gro` run matched the pre-optimization hbond, ring, half_cage, quasi_cage, cage, order_parameter, and ice sheets exactly. Its `run.summary_write` reported all output stages separately.

9. Package version and documentation
   - Updated `pyproject.toml` and `sqq.__version__` from `0.2.6` to `0.2.7`.
   - Updated README, developer design documentation, and the English/Chinese design DOCX files for the unified selector and workbook layout.
   - Root `sqq` / `sqq -h` output now places `SQQ version: 0.2.7   Release date: Jul 15, 2026` immediately before `usage:`. Added `sqq -v` and `sqq --version`; each prints only that version line and exits successfully. `analyze` and `init` help retain their prior layout.

### Compatibility

- The 0.2.7 default output intentionally differs from 0.2.6: Q6/Q12, MCG-1, and DHOP35 are no longer calculated unless selected.
- To reproduce the former default descriptor set, use `--order-parameter f3,f4,q6,q12,mcg1,dhop35`.
- Descriptor formulas and cage/ring detection are unchanged. PBC-aware occupancy is intentionally corrected for multi-atom guests crossing a box boundary; wrapped/repeated residue IDs now produce separate source-order molecules.
- Hydrate-cluster classification can change when `--cage-size` previously hid cages required by the topology. Cage report tables/files remain filtered exactly as requested.
- The new occupancy, cage-grow, DHOP, and F4 implementations are result-equivalent optimizations; `run.summary_write`, terminal timing, write robustness, and large-sheet cosmetics do not change scientific values.
- Triclinic input that was previously approximated is now rejected; convert it to orthorhombic form before analysis.
- Plotting scripts that required the removed `frame` sheet should read the dedicated analysis sheets instead.
- Invalid single-frame GRO/XYZ, non-finite coordinates, unsupported configuration modes, and invalid numeric settings now fail early instead of being accepted or coerced.
- New commands should use `--order-parameter`; legacy selection options are scheduled for removal after the compatibility period.
- New output control should use `--no-output`; old individual `--no-*` switches remain hidden for the compatibility period.
- `--no-output all` intentionally retains mandatory `run_config.yaml`.

## Version 0.2.6

### Short Summary

Version 0.2.6 tightens runtime metadata, worker control, and summary-write scalability. `--worker` / `-w` now distinguishes worker counts from CPU fractions by input form, terminal output and the `summary.xlsx` dashboard share the same `SQQ version` and requested/effective `Graph mode` wording, `run_config.yaml` records resolved run metadata while preserving raw configuration values, and `summary.xlsx` keeps quasi-cage reporting compact by moving exact quasi-cage isomers to CSV detail output. Scientific analysis algorithms, coordinates, topology counts, and molecule membership are unchanged from `0.2.5`.

### Main Changes

1. Worker parsing by input form
   - `-w 1` now means exactly one requested worker.
   - `-w 1.0` and `-w 100%` mean 100% of detected physical cores before the reserve-one-core clamp.
   - Decimal values must be in `(0, 1]`, percentages must be in `(0%, 100%]`, and positive integer text remains an explicit worker count.
   - Final workers remain capped by one reserved physical core, task count, and the Windows `ProcessPoolExecutor` limit.

2. Runtime and summary dashboard metadata alignment
   - The terminal Configuration block now includes `SQQ version` and uses labels aligned with the `summary.xlsx` dashboard.
   - After analysis finishes, the terminal prints a final run summary with finish time, duration, `SQQ version`, resolved `Graph mode`, worker policy, backend, and workers.
   - The `summary.xlsx` home sheet uses the same graph-mode display as the terminal final summary.

3. Requested/effective graph-mode display
   - Explicit modes display as `hbond`, `oo`, or `pairs`.
   - Auto mode displays as `auto -> hbond`, `auto -> oo`, or `auto -> mixed (hbond, oo)` when frames resolve differently.
   - Per-frame `*_info.md` adds `graph_mode` while retaining the effective `bond_mode` row.
   - Per-frame `connection_mode` columns are unchanged for plotting compatibility.

4. Resolved run metadata in `run_config.yaml`
   - Raw config values such as `graph.bond_mode: auto` and `parallel.workers: 1.0` are preserved.
   - A `run` block records `sqq_version`, requested/effective graph mode, display graph mode, worker request, worker policy, resolved worker count, backend, and math threads.

5. Compact quasi-cage workbook output
   - The `summary.xlsx` `quasi_cage` sheet now aggregates exact quasi-cage isomers into composition-level columns such as `5r_5²6³`, matching the compact style used by the `cage` sheet.
   - Exact nonzero quasi-cage isomer rows are written to `summary_detail/quasi_cage_isomer.csv` with `frame`, `time_ps`, `quasi_cage_type`, `isomer`, and `count`.
   - The `detail_index` sheet records `quasi_cage_isomer.csv` beside the other generated CSV detail tables.
   - This reduces workbook width and final openpyxl formatting time for long trajectories with many quasi-cage isomers.

6. Package version
   - Updated `pyproject.toml` and `sqq.__version__` from `0.2.5` to `0.2.6`.
   - Updated README, design documentation, and the English/Chinese design DOCX files.

### Compatibility

- Scientific analysis results do not change from `0.2.5`.
- Existing `0.2.5` commands remain compatible except for the clarified `-w 1` behavior: `1` is now one worker, while `1.0` or `100%` is the 100% physical-core fraction.
- `summary.xlsx` dashboard and terminal text change for clarity.
- The `summary.xlsx` `quasi_cage` sheet is intentionally more compact: columns are composition-level quasi-cage types, while exact isomer rows move to `summary_detail/quasi_cage_isomer.csv`.
- `run_config.yaml` gains a `run` metadata block but preserves existing raw config keys.

## Version 0.2.5

### Short Summary

Version 0.2.5 updates runtime usability around worker selection, the `summary.xlsx` home sheet, and single-file stage visibility. The preferred CLI spelling is now `--worker` / `-w`, summary dashboard result rows now report per-frame `min / mean / max`, and interactive serial runs highlight the active stage with bold bright-blue text. Scientific analysis algorithms, per-frame result sheets, coordinates, molecule membership, and topology counts are unchanged from `0.2.4`.

### Main Changes

1. Worker option and physical-core policy
   - Adds the preferred `--worker` / `-w` CLI option.
   - `--worker` accepts `auto`, physical-core fractions such as `50%`, `0.5`, or `1`, and explicit positive integer worker counts such as `4`.
   - Worker resolution uses detected physical cores when possible, reserves one physical core for the operating system, and remains capped by task count and platform limits. A 4-core/8-thread machine therefore resolves `-w 100%` or an overlarge worker request to at most 3 workers.
   - The existing YAML key remains `parallel.workers`, and the old `--workers` CLI spelling is retained as a hidden compatibility alias for existing scripts.

2. Summary dashboard home sheet
   - The first `summary.xlsx` sheet now labels the Configuration version row as `SQQ version` and keeps the config path under `Config file`.
   - The bottom dashboard block is renamed to `Analysis Results (min / mean / max)`.
   - Except for `Frames total / ok / failed`, result rows such as water/guest molecules, ring counts, cage counts, hydrate clusters, and ice-like waters are reported as per-frame min/mean/max rather than cross-frame sums.
   - This changes dashboard presentation only; one-row-per-frame sheets and scientific analysis values are unchanged.

3. Single-file progress highlight
   - Interactive serial progress keeps the compact three-row stage layout and now renders the active stage with ANSI bold plus bright blue (`RGB(0,0,255)`) instead of bold alone.
   - Non-interactive logs and `tqdm` postfix output remain plain text.
   - This changes terminal readability only and does not affect analysis or output files.

4. Package version
   - Updated `pyproject.toml` and `sqq.__version__` from `0.2.4` to `0.2.5`.
   - Updated README, design documentation, and the English/Chinese design DOCX files.

### Compatibility

- Scientific analysis results do not change from `0.2.4`.
- Prefer `--worker` / `-w`; existing `--workers` commands continue to run through the compatibility alias.
- Effective worker counts may resolve lower than older runs because SQQ now reserves one physical core for the system.
- The summary dashboard's result block now reports per-frame min/mean/max instead of cross-frame sums, so dashboard numbers may differ from older dashboards even though the underlying per-frame sheets are unchanged.
- Interactive single-file terminal output now uses color for the active stage; non-interactive text output remains plain.

## Version 0.2.4

### Short Summary

Version 0.2.4 fixes generated GRO compatibility on Windows and other locale-dependent readers, and moves large multi-row summary details out of `summary.xlsx` into CSV files. Structure directory names, filenames, and GRO title lines now use portable ASCII labels, while Markdown and Excel retain readable scientific superscripts. Analysis algorithms, counts, coordinates, and molecule membership are unchanged.

### Bug Fixes

1. Generated GRO files could fail to open in Windows/locale-dependent readers
   - Symptom: generated half-cage, quasi-cage, or cage GRO files could be rejected by tools such as MDAnalysis on Windows even though SQQ's fixed-column audit found valid atom counts, atom records, coordinates, and box lines.
   - Cause: Unicode superscript/subscript structure labels appeared in generated directory names, filenames, and GRO title lines; some readers decoded those bytes through the active system locale, such as GBK, before parsing the file.
   - Fix: generated structure labels and GRO title lines now use portable ASCII forms such as `5^126^2`, `hc_5r_5^5`, and `qc_5r_5^36^2_56566`. Markdown and workbook display labels keep readable scientific notation.

2. Long XTC/TRR runs could fail while writing `summary.xlsx`
   - Symptom: analysis completed all frames, then crashed at final workbook writing with an Excel size error such as `This sheet is too large`, especially for `-s 4,5,6` trajectories with many frames.
   - Cause: multi-row detail tabs, especially `cage_isomer`, expanded one frame into many rows and could exceed Excel's 1,048,576-row sheet limit.
   - Fix: `cage_occupancy`, `cage_isomer`, `hydrate_domain`, and optional `hydrate_cluster_detail` are written as UTF-8-SIG CSV files under `summary_detail/`; `summary.xlsx` keeps a lightweight `detail_index` sheet. `cage_isomer.csv` defaults to nonzero isomer rows plus per-frame totals.

### Main Changes

1. Portable GRO structure labels
   - Converts Unicode display labels that contain superscript/subscript structure notation to ASCII path/title forms, for example `5^126^2`, `hc_5r_5^5`, and `qc_5r_5^36^2_56566`.
   - Applies the conversion to half-cage, quasi-cage, and cage directories and filenames, and sanitizes every generated GRO title line.
   - Keeps atom/residue names, atom membership, coordinates, box vectors, report labels, and workbook labels unchanged.

2. Reader verification
   - The reported failure was reproduced with MDAnalysis on Windows: 10 of 43 generated `1200ns.gro` structure files failed when Unicode superscript/subscript bytes were decoded through the active GBK locale.
   - Before the fix, SQQ's strict fixed-column audit found valid atom counts, 44-character atom records, numeric coordinates, and valid box lines in all 43 files, isolating the failure to Unicode path/title compatibility.
   - After the fix, the 76-test suite passes; regression tests require ASCII-only generated paths/titles and successful GRO re-reading. All 43 files in the complete `1200ns.gro` structure set pass fixed-column validation and MDAnalysis loading.

3. Long-trajectory summary detail CSV output
   - Moves `cage_occupancy`, `cage_isomer`, `hydrate_domain`, and optional `hydrate_cluster_detail` out of the workbook body and writes them as UTF-8-SIG CSV files under `summary_detail/`.
   - Adds a lightweight `detail_index` sheet to `summary.xlsx` with each generated detail table path and dimensions.
   - Makes `cage_isomer.csv` default to observed nonzero isomer rows plus per-frame totals; `--cage-isomer-rows all` or `output.cage_isomer_rows: all` restores the full zero-filled matrix.
   - Prevents long XTC/TRR runs from failing at the final Excel-writing step when a multi-row detail sheet exceeds Excel's 1,048,576-row limit.

4. Package version
   - Updated `pyproject.toml` and `sqq.__version__` from `0.2.3` to `0.2.4`.
   - Updated README, design documentation, and the English/Chinese design DOCX files.

### Compatibility

- Scientific analysis results do not change.
- Generated GRO structure paths change from Unicode display labels to ASCII labels. Scripts that hard-coded the former Unicode paths should update those path strings.
- Source frame names are retained; the compatibility conversion applies to SQQ-generated scientific structure labels and GRO title text.
- `summary.xlsx` no longer contains the multi-row `cage_occupancy`, `cage_isomer`, `hydrate_domain`, or `hydrate_cluster_detail` sheets; these tables are available as CSV files in `summary_detail/` and indexed by `detail_index`.

## Version 0.2.3

### Short Summary

Version 0.2.3 adds true process-based multi-core analysis for independent GRO/XYZ files and selected XTC/TRR frames, accelerates the water graph, ring DFS, bounded quasi-cage growth, cage target pruning, order-parameter calculation, and shared cluster geometry, and retains the established `chordless`/`bounded` default results. It also adds MCG-1 and DHOP35 hydrate-nucleation order parameters by default, with optional MCG-3 and DHOP30 switches. Optional `shortest_path` rings and `exact` quasi-cage growth expose stricter or more complete searches explicitly. The release includes indexed half-cage closure, exact local sH cage-face fingerprints, unified sI/sII/sH domain expansion, and opt-in scientific cage validation. Non-orthogonal/triclinic boxes are not part of this update.

### Parallel and Search Performance

1. Process-based multi-core execution
   - Replaced the default independent-file `ThreadPoolExecutor` path with spawned `ProcessPoolExecutor` workers, so Python ring/quasi/cage loops can execute on separate CPU cores rather than sharing one GIL.
   - The main process alone renders progress and writes the global workbook. Workers receive configuration once, read and write one frame at a time, publish compact stage events through a process queue, and return one summary row.
   - `parallel.backend` defaults to `process`; `thread` remains a compatibility backend and `serial` provides one-process comparison. `--parallel-backend` selects the backend directly.
   - `parallel.math_threads` defaults to `1` and controls inherited OpenMP, OpenBLAS, MKL, Accelerate, NumExpr, and BLIS thread limits to prevent nested oversubscription.
   - Worker resolution uses CPUs available to the current process, the mode fraction, task count, and the Windows process-pool cap. Duplicate case-insensitive standalone file stems are rejected before concurrent output can collide.
   - A single indexed XTC/TRR trajectory can distribute selected raw frame indexes across workers. Every worker opens one private MDAnalysis Universe; coordinate arrays are not serialized from the parent.
   - Process submission uses a rolling queue capped at `3 * workers` tasks. This bounds Future/pickle overhead without reducing concurrency: 100 effective workers can still execute 100 tasks while keeping at most 300 submitted.
   - Selected trajectory indexes are grouped into contiguous batches sized automatically from 1 to 8. Parent and worker MDAnalysis readers are explicitly closed after use.

2. Water graph and ring search
   - Uses MDAnalysis `self_capped_distance` to generate orthorhombic cutoff candidates when available, then rechecks every distance and hydrogen-bond angle with the established float64 SQQ logic. The prior cell list remains the fallback.
   - Pre-sorts graph adjacency once, represents path membership with integer bits, rejects chords while extending the DFS, and removes reverse-direction duplicates before final canonical ordering.
   - `ring.definition: chordless` remains the default. `ring.definition: shortest_path` / `--ring-definition shortest_path` additionally applies the Franzblau all-pairs shortest-path criterion and can therefore change downstream patch/cage results.
   - Shortest-path validation caches each bounded BFS distance map by source and depth for the current frame; acceptance rules are unchanged.

3. Bounded and exact quasi-cage policies
   - `quasi_cage.search_policy: bounded` remains the default and preserves the established candidate limits and large-frontier representation.
   - L1 candidate lists now use precomputed adjacent-list compatibility and forward checking; repeated distance, expansion, geometry, growth-edge, and adjacency work is cached per frame.
   - Patch subset filtering uses ring-to-patch inverted indexes instead of scanning every pair.
   - `quasi_cage.search_policy: exact` / `--quasi-search-policy exact` enumerates connected outer-layer subsets up to `max_rings_per_layer` and preserves distinct `(patch, frontier)` states. It can add L2/L3 half-layers that bounded local neighborhoods do not represent.
   - Candidate, wall-combination, and layer-state limits now add explicit frame warnings instead of truncating silently.
   - Connected-subset, patch, and frontier identities use integer masks internally while preserving the existing deterministic traversal and output order.

4. Cage target and edge-state pruning
   - Each grow state carries the cage compositions still compatible with its 4/5/6 face counts; branches that fit no single target stop immediately instead of surviving under merged per-size maxima.
   - Face membership, edge-used-once, and edge-used-twice state is represented with integer bitsets. Adding a face promotes shared edges without copying a complete edge-count dictionary.
   - Edge-to-ring choices and compatible target compositions are also represented as bitmasks. Remaining face-edge incidence and parity provide exact necessary-condition pruning before a branch is expanded.
   - The minimum-remaining-value boundary-edge rule and deterministic candidate ranking are retained. Exceeding `cage.max_boundary_candidates` is reported as a frame warning.

5. Shared hydrate-cluster geometry and active guest centers
   - Ring plane centers/normals can be cached in the shared `RingTopologyIndex` and reused by hydrate-cluster shared-face resolution; sI/sII/sH fingerprints and growth rules are unchanged.
   - The cluster/domain hierarchy follows the HTR+ cage-graph concept (DOI 10.1088/1361-648X/ad52df) while retaining SQQ's explicit labelled-face fingerprints, strict seeds, and deterministic boundary rules.
   - `guest.center_mode` is now active: `center_atom`/`auto` uses a configured center atom when present and falls back to the residue centroid; `centroid` always uses the residue centroid.
   - Historical no-op defaults `ring.primitive` and `quasi_cage.mode` are no longer emitted. Explicit behavior is controlled by `ring.definition` and `quasi_cage.search_policy`.

6. Equivalent order and geometry calculation
   - F3 and graph-mode Q_l reuse one PBC-aware graph-neighbor vector cache. Cutoff/nearest/LAMMPS Q_l candidates are built once per pair, and every requested degree shares normalized vectors and spherical angles.
   - Spherical-harmonic normalization constants are cached. When both ring normals and scientific face-quality values are required, one SVD supplies both rather than fitting the face twice.

7. Verification
   - Added process policy/spawn integration tests, serial/process row-equivalence coverage, incremental ring and shortest-path tests, exact connected-layer tests, target/edge-budget tests, Q_l cache reference tests, trajectory batching/cleanup tests, and existing cluster/quasi regression coverage.
   - All 70 tests pass. The search/cache-only code completed the local `1200ns.gro` serial run in 18.2 s versus the 26.6 s pre-refinement baseline. The current default run including MCG-1 and DHOP35 completed in about 21.6 s; every overlapping pre-existing analysis column matched the earlier workbook.
   - The scheduling and search-cache refinements add no scientific-value change. The MCG/DHOP addition below intentionally adds `hydrate_order`, two optional CLI switches, and frame-level order-parameter output columns without changing pre-existing graph, ring, half/quasi, cage, cluster, occupancy, F3/F4/Q_l, or ice values.
   - The current implementation intentionally retains the existing orthorhombic box representation. GRO nine-term boxes and non-orthogonal/triclinic minimum-image support remain out of scope.

### MCG and DHOP Hydrate-Nucleation Order Parameters

1. Defaults and controls
   - Added `hydrate_order` configuration. MCG-1 and DHOP35 are enabled by default; MCG-3 and DHOP30 remain off unless enabled with `--mcg3 on` and `--dhop30 on` or YAML.
   - MCG defaults to methane-like residue names `CH4` and `MET`, a 0.90 nm guest-pair cutoff, a 0.60 nm shared-water cutoff, opposing 45-degree cones, and at least five coordinated waters.
   - DHOP uses a dedicated O-O graph with a 0.35 nm default cutoff for the all-atom TIP4P/Ice workflow. The cutoff is configurable, including 0.325 nm for the original mW-water definition. The DHOP35/DHOP30 names refer to angular thresholds.

2. Corrected and bounded implementation
   - MCG accepts five **or more** coordinating waters, builds connectivity only from qualifying MCG edges, and applies MCG-1/MCG-3 as one-pass degree filters before deterministic connected components.
   - DHOP evaluates each undirected central O-O bond once, accumulates the equivalent directed-center counts for both endpoints, accepts counts 11 or 12, applies the three-qualified-neighbor criterion, includes the first neighbor shell, and reports the largest tagged-water component.
   - Orthorhombic-PBC cutoff candidates use deterministic cell lists followed by exact float64 minimum-image checks. Dynamic adjacency sets remove the fixed neighbor and atom-index array limits in the reference companion program.
   - The implementation intentionally does not copy the companion code's hard-coded residue names, O(N^2) loops, uninitialized arrays, exactly-five-water MCG test, truncated Perl neighbor readers, or non-qualifying guest-distance links.

3. Output and compatibility
   - Per-frame `*_info.md` adds `Hydrate Nucleation Order Parameters` immediately after the existing F3/F4/Q_l section, reporting largest cluster size and member type. An unavailable MCG guest selection is `N/A`, distinct from a valid zero-sized cluster.
   - `summary.xlsx` adds `MCG-1` and `DHOP35` to `order_parameter`; `MCG-3` and `DHOP30` columns appear only when enabled.
   - These descriptors are read-only analyses. They do not change the active water graph, rings, half/quasi cages, closed cages, occupancy, cage-based hydrate clusters, F3/F4/Q_l, or ice classification. Existing scientific outputs retain their values; the order-parameter output schema gains the new columns and section.

4. Reference validation
   - On the published companion `panding.gro` sample (300 methane and 3487 water molecules), corrected SQQ geometry gives MCG-1=227, MCG-3=167, DHOP35=1112, and DHOP30=1009 with the 0.35 nm all-atom cutoff.
   - References: Barnes et al., MCG (DOI 10.1063/1.4871898); Knott et al., MCG nucleation coordinate (DOI 10.1021/jp507959q); DeFever and Sarupria, DHOP (DOI 10.1063/1.4996132); Li et al., all-atom hydrate nucleation (DOI 10.1073/pnas.2011755117).

### Main Changes

1. Equivalent quasi-cage search acceleration
   - Reuses the existing graph-edge-to-ring reverse index to construct candidate-ring adjacency from shared bonds instead of repeatedly comparing every candidate-ring pair.
   - Caches ring-center distances used for deterministic candidate ranking.
   - Caches L2/L3 expansion units by the sorted patch and frontier ring identifiers while retaining the existing per-seed state budget and ordered growth rules.
   - Rejects already-seen half-cage or quasi-cage ring sets before PBC unwrapping and patch-center construction, and reuses geometry for identical ring sets within one frame.
   - The compatible default remains `quasi_cage.search_policy: bounded` with `max_layers: 1`; opt-in `exact` is explicit, and no hidden frame-global state cap is introduced.
   - Added focused regression coverage for indexed adjacency, connected growth units, open L2 half-layers, distance/geometry caches, and max-layer behavior. Exact comparison against tag `v0.2.2` matched half/quasi results for `max_layers = 1/2/3`.

2. sH phase topology and domain growth
   - Added strict local sH fingerprints for `5^12`, `4^3 5^6 6^3`, and `5^12 6^8`, keyed by neighboring cage type and shared-face size.
   - The ideal fingerprints distinguish pentagonal `5^12` contacts, square medium-medium contacts, equatorial medium-large hexagonal contacts, and axial large-large hexagonal contacts. Shared-face incidence counts follow the ideal sH cell ratio `3:2:1`.
   - Retained the existing two-anchor composite sH seed as supplemental high-confidence evidence rather than the only sH recognition path.
   - Extended the existing compatible-edge, two-contact growth rule to sH, so all three hydrate phases use the same deterministic expansion and exclusive-domain construction.
   - Added regression coverage for a complete 24-cage synthetic sH face graph while preserving the prior composite-seed test.

3. Shared ring topology and indexed half-cage closure
   - Added one frame-local ring topology index containing `ring_by_id`, locally unwrapped ring centers, `edge_to_ring_ids`, ring adjacency, and a shared symmetric distance cache.
   - Half/quasi and cage analysis now reuse the same topology object instead of rebuilding ring centers and edge incidence independently.
   - Added deterministic indexed closure of connected two-, three-, and four-half-cage combinations. Face-count overflow and edge overuse are pruned before ordinary closed-polyhedron validation.
   - Generic grow remains the normal path and its detections are retained first. Fast closure is skipped when grow finishes within its budgets; if grow reaches a state limit, the recovery path can add a missed standard shell without discarding grow results.
   - `cage.fast_closure` defaults to `true`, `cage.fast_closure_max_states` defaults to `20000`, and `--cage-fast-closure on/off` provides direct comparison control.
   - On `1200ns.gro`, fast closure off/on completed in 32.25/32.35 s and all ten analysis sheets matched exactly; the recovery traversal was skipped because grow finished within budget.

4. Optional scientific cage validation
   - Added PBC-aware SVD face-planarity RMS, cyclic edge-length coefficient of variation, and projected-area measurements.
   - Added edge-connected face-shell and cyclic vertex-link checks to reject disconnected, pinched, or non-manifold shell topology.
   - Added oriented triangle-shell volume validation and tetrahedral volume-centroid calculation.
   - `cage.scientific_validation` defaults to `false` and is enabled with `--cage-scientific-validation on`. Default runs retain the prior topological acceptance and mean-water cage center.
   - When enabled, the explicit defaults are `max_face_planarity_rms_nm: 0.06`, `max_face_edge_cv: 0.35`, and `min_cage_volume_nm3: 1.0e-6`. Enabling the option can remove distorted cages and can change occupancy, ownership-filtered free-ring/free-patch output, or geometry-resolved cluster edges; raw ring/patch search, order, and ice calculations are unchanged.
   - On `1200ns.gro`, enabling the default scientific thresholds reduced accepted cages from 301 to 290 while leaving order-parameter and ice sheets unchanged.

5. Package version
   - Updated `pyproject.toml` and `sqq.__version__` from `0.2.2` to `0.2.3`.
   - Updated current-release references in the README and design documentation.

## Version 0.2.2

### Short Summary

Version 0.2.2 adds optional reported-cage hydrate clusters, geometry-resolved shared-face connectivity, strict local sI/sII/sH phase seeds, per-frame phase-domain expansion and boundary classification, cluster/domain workbook output, compact cluster-aware terminal stage displays, expanded per-frame report metadata, and singular CLI option names.

Comparison baseline: GitHub `pimooni/sqq` tag `v0.1.6` (commit `97996b2`) -> local `0.2.2`.

### Main Changes

1. Hydrate cluster graph
   - Added optional hydrate cluster analysis controlled by `--hydrate-cluster on/off` or `hydrate_cluster.enabled`; it remains off by default.
   - Builds the graph from the final reported cage set, so `--cage-size` controls which cages can enter cluster and phase classification.
   - Connects cages through a complete shared ring face. When more than two detected cages reference one face, ring-plane geometry selects at most one cage on each physical side; without geometry, only an unambiguous two-cage face is accepted.
   - Reports connected components with at least `--cluster-min-cage` cages and counts smaller components as isolated cages.

2. Phase seeds, expansion, and domains
   - Added labelled first-shell fingerprints keyed by neighboring cage type and shared-face size.
   - Added strict sI and sII seed templates for `5^12`, `5^12 6^2`, and `5^12 6^4` environments with a count tolerance of one and no unexpected fingerprint labels.
   - Added conservative composite sH seeds made from two separated `5^12 6^8` anchors, six common `5^12` cages, six `4^3 5^6 6^3` cages, and an adjacent medium-cage bridge between the anchors.
   - Expands sI and sII from strict seed members through compatible labelled edges. A candidate must have a partial compatible fingerprint and at least two accepted phase contacts. sH currently grows only by the union of overlapping strict sH seeds.
   - Collects phase claims independently, then forms deterministic per-frame domains only from cages claimed exclusively by one phase and connected through compatible phase edges. Every domain must contain at least one strict seed anchor.

3. Boundaries and cluster labels
   - Keeps cages with competing phase claims out of exclusive domains and labels supported single-phase or multi-phase boundary cages after domains are formed.
   - Separates multi-phase/interphase boundaries, single-phase boundaries, ambiguous cages, and unclassified cages without forcing a label from cage composition alone.
   - Labels a cluster `sI`, `sII`, or `sH` when all of its domains have one phase, `mixed` when multiple domain types occur, and `unclassified` when no phase domain is found.
   - Domain and cluster identifiers are deterministic within a frame; temporal domain tracking is not implemented.

4. Hydrate cluster output
   - Adds per-frame `Hydrate Cluster`, `Hydrate Cluster Hierarchy`, `Hydrate Cluster Detail`, `Hydrate Domain`, and `Hydrate Boundary` sections when analysis is enabled.
   - The hierarchy table reports cluster, domain, and boundary cage quantities by cage composition; seed counts remain in domain detail rather than a dedicated hierarchy column.
   - Adds `hydrate_cluster` and `hydrate_domain` workbook sheets whenever hydrate analysis is enabled.
   - Adds `--cluster-detail on/off` and `hydrate_cluster.detail`; when enabled, `hydrate_cluster_detail` adds one row per cluster.
   - Public motif output is not generated in 0.2.2. The compatibility return slot and model remain internal/empty, and no `hydrate_motif` sheet or Markdown section is written.
   - Hydrate clusters do not produce additional GRO structure files.

5. Terminal progress display
   - Serial progress now shows the full three-row workflow instead of a single long current-stage phrase: file preparation (`reading`, `settings`, `selecting`), core search (`graph`, `ring`, `half/quasi`, `cage`, optional `cluster`), and post-processing/output (`filtering`, `order`, `ice`, `output`).
   - Interactive serial terminals bold the active stage while keeping `stage / frame / total` timings unchanged; non-interactive postfix output uses the same short stage labels without ANSI styling.
   - Parallel `stage_summary` now uses compact `stage:count` cells instead of wide `|`-separated columns. Each column width is recalculated from the longest current `stage:count` cell in that column, then followed by two spaces.
   - The `cluster` stage is shown only when hydrate cluster analysis is enabled, and it remains at the end of the core-search row. With `--hydrate-cluster off`, the summary omits `cluster` entirely instead of showing `cluster:0`.

6. Frame metadata
   - Per-frame `info.md` Frame Information now starts with `sqq version`, `date & time`, `source`, `frame`, and `time_ps`, followed by bond mode, ring sizes, status, and molecule counts.
   - `date & time` is generated at report-writing time using the local timezone label, and `source` is written as an absolute path without Markdown code quotes.

7. CLI and package metadata
   - Renamed `--quasi-sizes`, `--quasi-base-sizes`, `--quasi-side-sizes`, `--quasi-max-layers`, and `--max-cage-faces` to `--quasi-size`, `--quasi-base-size`, `--quasi-side-size`, `--quasi-max-layer`, and `--max-cage-face`.
   - Existing YAML keys, including `quasi_cage.max_layers` and `cage.max_faces`, remain unchanged.
   - Updated package metadata from `0.1.6` to `0.2.2`.
   - Added hydrate-cluster unit tests covering physical shared-face resolution, ideal phase seeds, incomplete or incorrect fingerprints, expansion, mixed boundaries, separated same-phase domains, output sheets, hierarchy rendering, and CLI settings.

### Compatibility

- With hydrate cluster analysis disabled, the 0.1.6 graph, ring, patch, cage, occupancy, order-parameter, and ice workflow is unchanged.
- Hydrate classification depends on the final cage report scope; excluding a required cage type can prevent the corresponding phase seed or domain from being recognized.
- Phase and boundary assignments are per-frame topological classifications, not temporal grain tracking or crystallographic orientation matching.
- New commands should use the singular long option names. Historical notes retain the option spelling used by the corresponding release.
- Terminal progress text layout changed; analysis results and output file schemas are unaffected by the stage-display update.

## Version 0.1.6

### Short Summary

Version 0.1.6 makes cage reporting follow the selected search scope, adds Q_l order parameters, and improves per-frame cage/isomer readability while keeping the scientific cage acceptance rules unchanged.

### Main Changes

1. Cage report scope now follows `--size` by default
   - Changed the built-in `cage.report_types` value from a fixed cage list to `auto`.
   - With no explicit cage report filter, every detected cage composition allowed by `-s` / `--size` is included in counts, occupancy tables, GRO files, per-frame info, and `summary.xlsx`.
   - For example, `-s 4,5,6` now reports accepted cages containing 4-, 5-, and 6-membered faces instead of silently limiting output to the former 5/6 named set.

2. Cage report groups
   - The documented `--cage-size` values are `auto`, `all`, `I`, `II`, `H`, `HS-I`, `TS-I`, and `I2II`.
   - Group names expand to exact cage compositions and overlapping types are reported once.
   - `-s 4,5,6 --cage-size I,II` reports 4/5/6 rings and quasi-cages while cage output is restricted to 5¹², 5¹²6², and 5¹²6⁴.

3. CLI and documentation correction
   - Removed obsolete concatenated cage-label examples from README and CLI help.
   - Added the complete `--cage-size` value domain and corrected option defaults, aliases, and size-scope descriptions.

4. Q_l order parameters
   - Added local Steinhardt/LAMMPS-style Q_l calculation for water oxygens.
   - Default reported degree list is Q6 and Q12 through `order.q_degree: [6, 12]`.
   - The same interface can report common LAMMPS degrees such as `4,6,8,10,12`.
   - Added Q_l neighbor modes: `graph`, `cutoff`, `nearest`, and `lammps`.
   - Added Q_l settings: `order.q_enabled`, `order.q_degree`, `order.q_neighbor_mode`, `order.q_cutoff_nm`, and `order.q_n_neighbor`.
   - Added CLI overrides: `--no-q`, `-q` / `--q-degree`, `--q-neighbor-mode`, `--q-cutoff`, and `--q-n-neighbor`.

5. Order-parameter output
   - Renamed the workbook `f3f4` sheet to `order_parameter`.
   - The sheet now reports `F3_mean`, `F3_count`, `F4_mean`, `F4_count`, and one mean/count pair for each requested Q_l degree.
   - Per-frame Markdown reports now use an `Order Parameters` section for F3, F4, and Q_l.
   - Added optional per-water `*_order_parameter.tsv` output through `--write-order-tsv` or `output.write_order_tsv`.

6. Per-frame info report readability
   - Merged the previous separate `Cage` and `Cage Isomer` Markdown sections into one `Cage` section.
   - Cage rows now show cage composition as the parent row and structural isomers as synchronized child rows.
   - Added `Quasi Cage Isomer Description` and `Cage Isomer Description` sections to each per-frame `*_info.md` report.
   - Each observed quasi-cage or cage isomer receives its own description row with the observed count.
   - `Cage Occupancy` remains separate because it describes guest assignment rather than cage topology.
   - The global `summary.xlsx` remains plotting-oriented and keeps separate workbook sheets where useful.

7. Documentation and tests
   - Updated README and design documentation for cage report scope, Q_l definitions, neighbor modes, and output layout.
   - Added Q_l reference tests for simple-cubic and tetrahedral bond-vector arrangements plus degree parsing.

8. Package version
   - Updated package metadata from `0.1.5` to `0.1.6`.

### Compatibility

- Cage detection, closure validation, ownership filtering, and search limits are unchanged.
- Use cage groups rather than implementation-specific concatenated labels in new commands and configuration files.
- Runs that relied on the old implicit four-type cage filter should now specify `--cage-size I,II,I2II`.
- The old internal F3/F4 data model is retained for compatibility, but user-facing workbook output now uses `order_parameter`.
- Existing optional `*_f3f4.tsv` output is replaced by `*_order_parameter.tsv` when order TSV output is enabled.
- Q_l is a diagnostic order parameter and does not alter graph construction, ring/cage detection, ownership filtering, guest occupancy, F3/F4, or ice classification.
- Per-frame `*_info.md` cage formatting changed for readability; `summary.xlsx` remains the stable plotting-oriented output.

## Version 0.1.5

### Short Summary

Version 0.1.5 adds diagnostic coordination distributions, named Type H cages, and explicit search/report scopes for rings and cages.

### Main Changes

1. Hydrogen-bond and connectivity coordination diagnostics
   - Added per-frame degree distributions for water nodes: 0, 1, 2, 3, 4, and greater than 4.
   - Added counts, fractions, mean coordination, degree <=2 fraction, four-coordinated fraction, and over-four fraction.
   - The active `hbond`, `oo`, or `pairs` workbook tab now contains the same plotting-friendly diagnostics.
   - Diagnostics read the existing graph only; they do not alter edges or any ring, patch, cage, F3/F4, or ice result.

2. Type H named cages
   - Added named reporting for 5¹²6⁸ and 4³5⁶6³.
   - Added superscript report labels and distinct cage-center atom names for both Type H cages.

3. Separate cage detection and reporting
   - `-s` / `--size` now defines the ring-face sizes searched by ring, quasi-cage, and cage analysis.
   - Cage search generates all trivalent Euler-compatible 4/5/6 face compositions up to `--max-cage-faces` (default 20).
   - `--cage-size` filters reported cage groups; `--cage-size all` reports every detected composition.
   - All detected cages continue to control half-cage, quasi-cage, and free-ring ownership even when their type is not reported.

4. Separate ring detection and reporting
   - Added `--ring-size` and `ring.report_sizes` to filter ring tables and GRO files after detection.
   - The report-size list must be a subset of the search sizes.

5. CLI and configuration cleanup
   - Renamed `--sizes` to `--size`.
   - Replaced the old search-oriented `--ring-sizes` and `--cage-sizes` options with report-oriented `--ring-size` and `--cage-size`.
   - Replaced `--other-cages`, `--no-other-cages`, and `--other-max-faces` with `--cage-size all` and `--max-cage-faces`.
   - Replaced `cage.ring_sizes`, `cage.target_types`, `cage.output_other`, and `cage.other_max_faces` with `cage.report_types` and `cage.max_faces`.

6. Package metadata and terminal labels
   - Added `J. PANG` and `Q. SUN` as package authors and retained `pimooni@gmail.com` as contact email.
   - Terminal graph modes now display readable names while CLI/config identifiers remain `auto`, `hbond`, `oo`, and `pairs`.
   - Retained `-m` as the short form of `--mode`.
   - Updated package version metadata from `0.1.4` to `0.1.5`.

### Compatibility

- Existing graph construction and topology acceptance criteria are unchanged.
- Commands and YAML files using the removed plural/other-cage options must be updated to the new search/report options.
- Cage-network or crystal-domain analysis is intentionally not included in this release.

## Version 0.1.4

### Short Summary

Version 0.1.4 improves CLI help, live per-worker progress, and per-frame Markdown readability. It also broadens the local DOCX ignore rule without changing topology-analysis results or workbook schemas.

### Main Changes

1. Clearer command-line help hierarchy
   - Updated the help description to `SQQ (Shell Qualification Quantifier): Python Joint Toolkit for Water-Shell Topology Analysis.`
   - Reduced `sqq` and `sqq -h` to a concise top-level overview with subcommands and quick-start examples.
   - Moved analysis examples, mode descriptions, and output-layout guidance to `sqq analyze -h`.
   - Removed the duplicated manually listed analysis options from the top-level help; argparse remains the authoritative detailed option reference.

2. Live progress for parallel file analysis
   - Added thread-safe per-worker stage reporting for multi-file GRO/XYZ runs.
   - The live panel reports completed, failed, active, and queued files plus total elapsed time.
   - Run metadata now formats time zones with their UTC offset, for example `China Standard Time (+8)`; the full China name is applied only to UTC+8.
   - `stage_summary` always shows all 11 work stages in three logical rows: preparation, topology search, and post-processing/output.
   - Up to six active files show their current stage, stage elapsed time, and file elapsed time; additional active files are summarized.
   - The serial one-worker progress display is unchanged.

3. Clearer per-frame Markdown reports
   - Ring counts now show only configured sizes and use the final free-ring population; cage counts list one topology type per row.
   - Half-cage and quasi-cage labels omit internal prefixes and group isomers below composition totals.
   - The Molecules table now reports both molecule and atom counts in source order.
   - The active network section reports the selected connection type, connection count, and mean connections per water.
   - Cage occupancy now lists cage types by row and exact guest compositions by column in source guest order.
   - Cage isomers now use nested readable labels instead of a wide type-by-isomer matrix.
   - Frame, molecule, active-connection, F3/F4, and ice information is separated into compact sections.
   - `summary.xlsx` remains unchanged with one input file per row in analysis sheets.

4. Code and documentation cleanup
   - Removed obsolete Markdown-report helpers and simplified canonical cage ordering without changing counts.
   - Tightened time-zone alias handling so non-China `CST` offsets retain their original name.
   - Updated README, design documentation, and the English/Chinese DOCX design guides.

5. Broader local-document ignore rule
   - Replaced two filename-specific DOCX ignore entries with the global `*.docx` rule.
   - Design DOCX files remain local and are not included in repository uploads.

6. Version metadata
   - Updated `pyproject.toml` and `sqq.__version__` from `0.1.3` to `0.1.4`.

### Compatibility

- Analysis algorithms, scientific defaults, `summary.xlsx` schemas, mode presets, worker allocation, and bond-mode behavior are unchanged from `0.1.3`. The per-frame Markdown layout and parallel progress reporting changed.
- Existing `sqq analyze` commands and configuration files remain compatible.

## Version 0.1.3

### Short Summary

Version 0.1.3 adds reproducible analysis presets and mode-based file-level worker allocation while keeping quasi-cage layer depth and output selection under explicit user control.

### Main Changes

1. Added analysis mode presets
   - `-m 00`: rigorous hydrogen-bond analysis with 4/5/6 ring and cage faces, unconventional cages enabled, and 25% automatic workers.
   - `-m 50`: standard default using `auto` graph selection, 5/6 faces, standard cages, and 50% automatic workers.
   - `-m 99`: performance screening using O-O connectivity, 5/6 faces, standard cages, and 90% automatic workers.
   - Modes do not change `quasi_cage.max_layers`; all modes retain the L1 default.

2. Unified parallel terminology
   - Replaced `--n-jobs` and `parallel.n_jobs` with `--workers` and `parallel.workers`.
   - Automatic worker counts use the mode CPU fraction and are capped by the number of independent GRO/XYZ files.
   - Single-file and XTC/TRR runs remain serial.

3. CLI and reporting updates
   - Added `-s` as the short alias for `--sizes`.
   - Added `-b` / `--bond-mode {auto,hbond,oo,pairs}` as an explicit graph-mode override.
   - Documented `--no-info` beside `--no-gro` and `--no-xlsx`.
   - Terminal and `summary.xlsx` configuration sections now report the active mode, worker policy, and actual worker count.

4. Project description
   - Standardized the README description as `SQQ: Python Joint Toolkit for Water-Shell Topology Analysis.`

## Version 0.1.2

### Short Summary

Version 0.1.2 is mainly a speed and usability update. The closed-cage rules are unchanged: accepted cages still require every edge to be used exactly twice, `V - E + F = 2`, and target face-count matching.

### Main Changes

1. Faster default quasi-cage search
   - `quasi_cage.max_layers` now defaults to `1`.
   - Routine runs report L1 quasi-cages and standard half-cages by default.
   - L2/L3 dangling outer layers remain available with `quasi_cage.max_layers: 2/3` or `--quasi-max-layers 2/3`.

2. New CLI option
   - Added `--quasi-max-layers N` for temporary layer-depth control without editing `config.yaml`.

3. Faster cage grow search
   - Standard target cages now share one grow pass instead of repeating the DFS separately for `512`, `51262`, `51263`, and `51264`.
   - Boundary growth first uses exact shared-edge topology.
   - The grow step chooses the boundary edge with the fewest addable rings.
   - Ring-center distances are calculated only when candidate truncation is needed.

4. Documentation cleanup
   - Simplified `README.md` for installation, quick start, common options, defaults, and output layout.
   - Kept algorithm details in `docs/design.md`.

5. Report and terminal usability
   - Reworked the first `summary.xlsx` sheet into a dashboard-style run overview.
   - Added live per-frame stage text in the default terminal progress display, such as graph building, ring search, cage search, and output writing.

6. Strict L1 quasi-cage fix
   - Default quasi-cage output now remains strict L1 when `quasi_cage.max_layers = 1`.
   - The special L2 `6^1` check is restricted to the standard `hc_6r_5^6_6^1` half-cage and is no longer reported as a quasi-cage.

## Version 0.1.1

### Short Summary

The core 0.1.1 change is replacing the old `cup` workflow with `half_cage` and `quasi_cage`, while unifying ring, open-patch, and closed-cage search around shared-edge topology lookups.

Comparison baseline: GitHub `pimooni/sqq` 0.1.0 -> local 0.1.1.

### Main Changes

1. Terminology update
   - Removed the old `cup` concept.
   - Current structure hierarchy:

```text
ring -> half_cage -> quasi_cage -> cage
```

2. Added `half_cage` and `quasi_cage`
   - `half_cage` is reserved for standard half-cage patches.
   - `quasi_cage` is used for non-closed cage-like patches.
   - Layer labels now follow L0/L1/L2/L3.
   - L1 must be closed; L2/L3 may be dangling rings.

3. Improved quasi-cage search
   - Candidate rings are found with a shared-edge reverse lookup.
   - `ring CNT` is used only for candidate ordering and pruning.
   - This avoids repeated full-ring scans and makes runtime more stable.

4. Improved cage search
   - Closed cages grow from ring faces along boundary edges.
   - Candidate faces are selected by shared boundary edges first, then ordered by ring-center distance.
   - Accepted cages must satisfy:

```text
each edge is used by exactly two faces
V - E + F = 2
face counts match the target cage type
```

5. Clearer size control
   - `ring` supports 4/5/6/7 and defaults to 5/6.
   - `half_cage` and `quasi_cage` can follow ring sizes, including 4/5/6/7.
   - `cage` supports 4/5/6 and defaults to 5/6. It does not use 7-ring faces.

6. Output layout update
   - Structure files are grouped by class:

```text
ring/
half_cage/
quasi_cage/
cage/
ice/
```

   - `summary.xlsx` has clearer sheet separation.
   - Each frame still writes an `_info.md` file by default.

7. More detailed cage occupancy
   - Empty, occupied, guest-specific, and multi-guest counts are reported.
   - Exact guest-composition labels are supported, for example:

```text
CH4
CH4x2
CH4+CO2
```

8. CLI and config updates
   - Added or improved:

```text
--sizes
--ring-sizes
--quasi-sizes
--quasi-base-sizes
--quasi-side-sizes
--cage-sizes
--other-cages
--no-other-cages
--output-layout
--no-ring-gro
--no-half-cage-gro
--no-quasi-cage-gro
--no-cage-gro
```

9. Documentation updates
   - README now describes the 0.1.1 workflow.
   - Developer design notes are available at `docs/design.md`.
