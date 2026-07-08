# SQQ Update Notes

This file records versioned update notes. New releases should be appended above older entries.

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
