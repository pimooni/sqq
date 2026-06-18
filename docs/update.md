# SQQ Update Notes

This file records versioned update notes. New releases should be appended above older entries.

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
