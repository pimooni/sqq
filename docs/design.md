# SQQ Design Notes

SQQ means **Shell Quant Qualifier**.

SQQ analyzes water-shell topology in MD trajectories with this pipeline:

```text
input frames
  -> molecule selection
  -> hydrogen-bond / O-O graph
  -> primitive chordless rings
  -> cups (half-cages)
  -> cages and guest occupancy
  -> F3/F4 and ice metrics
  -> per-frame and global outputs
```

## Current Coding Scope

Implemented in the current source tree:

- `sqq init`
- `sqq analyze`
- YAML config loading and command-line overrides
- `.gro` and `.xyz` frame input; `.xtc` and `.trr` through MDAnalysis when installed, with `input.xtc_stride`
- orthorhombic minimum-image PBC
- water and guest selection
- external pair-file graph input through `bond_mode=pairs`
- `bond_mode=auto/hbond/oo/pairs`
- non-recursive DFS ring search; 4/5/6/7 rings are supported, with `[5, 6]` as the default
- general cup search from base-ring and side-ring closure; 4/5/6/7 cups are supported when enabled through `ring.sizes`, with `[5, 6]` as the default
- ring-face grow cage search with closed-polyhedron validation; cage faces support 4/5/6, with `[5, 6]` as the default
- optional cup-pair cage comparison mode
- cage guest assignment by point-in-polyhedron, with center-distance comparison/fallback mode
- per-cage-type empty/occupied/multi/guest and cage-isomer summary columns
- F3/F4 metrics based on the reference-script formulas, including optional `order.focus_waters` averages
- CHILL-style ice-like / ice-I-like / interfacial ice classification
- per-frame folder output
- `*_info.md`, `*_view.vmd.tcl`
- optional `*_membership.tsv` and `*_f3f4.tsv` when `output.write_tsv: true`
- free ring, free cup, cage, and ice GRO output
- global `summary.xlsx`, `run_config.yaml`
- explicit frame-level parallel execution for independent `.gro`/`.xyz` files

## Cage Principle

The default cage mode is `cage.search_mode = grow` with `cage.seed_mode = ring`.

A cage is detected by starting from ring faces, then growing along open boundary edges until the shell closes. A cup-seeded grow path is also available with `cage.seed_mode = cup`, where detected cups (half-cages) are used as the initial connected face patches. A candidate is accepted only when:

- every edge is used by exactly two ring faces;
- Euler characteristic satisfies `V - E + F = 2`;
- face counts match a target cage type such as `512`, `51262`, `51263`, or `51264`.

By default, cage search uses only 5- and 6-member faces and reports the conventional hydrate cage targets. Setting `cage.ring_sizes = [4, 5, 6]` allows 4-member faces. Setting `cage.output_other = true` additionally generates Euler-compatible 4/5/6 unconventional cage targets up to `cage.other_max_faces`. Seven-member faces are intentionally excluded from cage detection in the current design.

Cage isomer labels are based on the adjacency pattern among 6-ring faces, using readable labels such as `6adj`, `6chain3`, `6star3`, or `6tri3+single`.

`cage.search_mode = pair` is retained only as a compatibility/debug mode for HA/GRADE-style cup-pair comparison. It is not the default cage definition.

`cage.seed_mode = cup` is retained for the cup-to-cage logic discussed during design. In current HA example tests, `seed_mode = ring` is faster, so it remains the default. `cage.seed_mode = auto` uses cup seeds when cups exist and falls back to ring seeds otherwise.

## Ownership Rule

Output ownership follows:

```text
cage > cup > ring
```

Rings or cups consumed by a cage are not counted or written as free cup/free ring objects. Ring and cup outputs do not use guest occupancy; cage outputs do. The default cage occupancy mode is `polyhedron`, which tests the configured guest center against the closed ring-face polyhedron. `center` reproduces the simpler distance mode, and `auto` accepts either test.

## Remaining Refinements

- refine atomistic Ih/Ic stacking assignment if separate hexagonal/cubic ice counts are needed.

## Naming

Package name, import name, and command name are all:

```text
sqq
```

Examples:

```powershell
sqq init -o config.yaml
sqq analyze -i ./gro -c config.yaml -o ./result_sqq
```











