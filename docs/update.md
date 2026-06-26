# SQQ Update Notes

This file records versioned update notes. New releases should be appended above older entries.

## Version 0.1.6

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

### Short Summary

Version 0.1.6 makes cage reporting follow the selected search scope, adds Q_l order parameters, and improves per-frame cage/isomer readability while keeping the scientific cage acceptance rules unchanged.

## Version 0.1.5

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

### Short Summary

Version 0.1.5 adds diagnostic coordination distributions, named Type H cages, and explicit search/report scopes for rings and cages.

## Version 0.1.4

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

### Short Summary

Version 0.1.4 improves CLI help, live per-worker progress, and per-frame Markdown readability. It also broadens the local DOCX ignore rule without changing topology-analysis results or workbook schemas.

## Version 0.1.3

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

### Short Summary

Version 0.1.3 adds reproducible analysis presets and mode-based file-level worker allocation while keeping quasi-cage layer depth and output selection under explicit user control.

## Version 0.1.2

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

### Short Summary

Version 0.1.2 is mainly a speed and usability update. The closed-cage rules are unchanged: accepted cages still require every edge to be used exactly twice, `V - E + F = 2`, and target face-count matching.

## Version 0.1.1

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

### Short Summary

The core 0.1.1 change is replacing the old `cup` workflow with `half_cage` and `quasi_cage`, while unifying ring, open-patch, and closed-cage search around shared-edge topology lookups.
