# SQQ Development Design

SQQ means **Shell Quant Qualifier**. This document records the current implementation logic for developers, so the code and the scientific definitions stay aligned.

## Pipeline

```text
input frames
  -> molecule selection
  -> water graph: hydrogen bond / O-O / user pair map
  -> primitive chordless rings
  -> half_cage and quasi_cage open patches
  -> closed cage search and guest occupancy
  -> F3/F4 and ice metrics
  -> per-frame outputs and summary.xlsx
```

The shared water graph is used by ring, half_cage, quasi_cage, cage, F3/F4, and ice analysis. The graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode.

## Modules

- `sqq/pipeline.py`: top-level analysis order, config merging, frame loop, output dispatch.
- `sqq/core/graph.py`: water graph construction with orthorhombic PBC and nearby-pair search.
- `sqq/core/ring.py`: non-recursive DFS ring search.
- `sqq/core/quasi_cage.py`: layered `half_cage` and `quasi_cage` search.
- `sqq/core/cage.py`: closed cage grow search, polyhedron validation, guest assignment.
- `sqq/core/f3f4.py`: F3/F4 order metrics.
- `sqq/core/ice.py`: CHILL-style ice classification.
- `sqq/io/summary.py`: per-frame info and global workbook tables.
- `sqq/io/gro_writer.py`: grouped or flat GRO structure output.

## Ring Search

Rings are searched on the already-built water graph. The algorithm does not use geometric distance at this stage; it follows graph adjacency. The implementation uses an iterative DFS instead of recursive DFS.

Current behavior:

- supported ring sizes: 4, 5, 6, 7;
- default ring sizes: 5, 6;
- default ring filter: chordless primitive rings;
- ring nodes are water oxygen indices.

## Half-Cage and Quasi-Cage Terms

`patch` means a connected set of ring faces during search. It is not an output class by itself.

Layer definitions:

```text
L0 = base ring
L1 = side rings sharing every base-ring edge; L1 must close into a full side wall
L2 = rings grown outward from L1
L3 = rings grown outward from L2
```

L2 and L3 may be dangling rings or connected dangling ring chains. They do not need to close.

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

Search order:

1. Choose one allowed ring as `L0`.
2. For every edge of `L0`, use `edge_to_rings` to find rings sharing exactly that base edge.
3. Sort those candidates by ring-center distance to `L0`; keep at most `quasi_cage.max_candidates_per_edge` per base edge unless the limit is disabled.
4. Build L1 with early pruning: each neighboring side ring must share one non-base edge, and the last side ring must connect back to the first.
5. Reject shifted or overlapped L1 choices with the expected unique-water count check.
6. Classify the L0+L1 patch:
   - if it matches a standard `half_cage`, store it as `half_cage`;
   - otherwise store it as `quasi_cage`.
7. For possible `hc_6r_5^6_6^1`, inspect L2 candidates and classify that larger standard patch as `half_cage`.
8. Grow L2/L3 from the current frontier:
   - find exposed frontier edges;
   - use `edge_to_rings` to get rings sharing those edges;
   - reject rings already in the patch;
   - reject rings touching lower-layer edges, so layer assignment uses the lowest possible layer;
   - sort by ring-center distance and keep at most `quasi_cage.max_layer_candidates`.
9. Convert layer candidates into bounded connected growth units. Small connected components are used whole; large components are represented by single dangling rings and local connected neighborhoods.
10. Each new layered patch is classified again as `half_cage` or `quasi_cage`.

The algorithm intentionally avoids scanning all rings for L1/L2/L3 once `edge_to_rings` is built.

Important limits:

- `quasi_cage.max_layers`: default 3, so L1/L2/L3 are reported.
- `quasi_cage.max_rings_per_layer`: maximum rings in one L2/L3 growth unit.
- `quasi_cage.max_layer_states_per_seed`: per-seed cap for finite growth.
- `quasi_cage.max_candidates_per_edge`: L1 candidate cap per base edge.
- `quasi_cage.max_layer_candidates`: L2/L3 candidate cap per frontier.

## Cage Algorithm

The default cage mode is:

```text
cage.search_mode = grow
cage.seed_mode = ring
```

The code can also grow from `half_cage`/`quasi_cage` patch seeds with `cage.seed_mode = patch`, but ring seeds are the default because they are faster in the tested frames.

Target cages:

```text
512   = 5^12
51262 = 5^12 6^2
51263 = 5^12 6^3
51264 = 5^12 6^4
```

Default cage search uses 5- and 6-member faces. 4-member faces can be enabled with `cage.ring_sizes = [4, 5, 6]`. 7-member faces are intentionally not used for cage detection.

Grow logic:

1. Build `edge_to_ring_ids` from all allowed cage rings.
2. Start from a seed face set, usually one ring.
3. Count how many ring faces use each edge in the current patch.
4. Edges used once are open boundary edges; edges used twice are closed; edges used more than twice are invalid.
5. Choose the most constrained boundary edge, meaning the boundary edge with the fewest addable rings.
6. For that edge, use `edge_to_ring_ids` to find rings sharing the boundary edge.
7. Reject a candidate ring if:
   - it is already in the patch;
   - it would exceed the target face count;
   - it would make any edge used more than twice;
   - it violates the single-ring seed rank pruning.
8. Sort remaining candidates by ring-center distance to the current patch and keep at most `cage.max_boundary_candidates`. Set this parameter to `0` to keep all topology-valid candidates.
9. Continue DFS growth until either the shell closes or the branch can no longer match the target face counts.

Single-ring seed rank pruning avoids rediscovering the same cage from every face: when the seed is one ring, later growth cannot add a ring whose stable id ranks before the seed ring.

Acceptance criteria:

- all edges are used exactly twice;
- Euler characteristic satisfies `V - E + F = 2`;
- face counts match the target cage type;
- the same `(cage_type, water_set)` was not already accepted.

## Guest Occupancy

Cage centers are computed from the locally unwrapped O coordinates of the cage waters. Hydrogens are not used in the cage center.

For each accepted cage, SQQ checks all selected guest molecules. The guest point is:

- the configured center atom if present, such as carbon in methane or carbon dioxide;
- otherwise the geometric center of all atoms in the guest molecule.

Default `cage.occupancy_mode = polyhedron` triangulates the cage ring faces and uses a point-in-polyhedron solid-angle test. `center` uses only a center-distance cutoff, and `auto` accepts either method.

## Ownership Rule

Output ownership follows:

```text
cage > quasi_cage > half_cage > ring
```

Rules:

- patches consumed by a cage are not written as free `half_cage` or `quasi_cage`;
- a `half_cage` whose ring set is a true subset of a `quasi_cage` is consumed by that `quasi_cage`;
- nested `half_cage` results are reduced to the larger `half_cage` patch;
- free rings are rings not consumed by cage, quasi_cage, or half_cage outputs;
- guest occupancy is only applied to closed cages.

## Output Layout

Per-frame output folders use the default grouped structure:

```text
frame_name/
  frame_name_info.md
  ring/
  half_cage/<type>/
  quasi_cage/<type>/
  cage/<type>/
  ice/
```

The global workbook is `summary.xlsx`. It contains run metadata, per-frame counts, graph statistics, ring/half_cage/quasi_cage/cage tables, occupancy tables, F3/F4, ice, and config sheets.

## Current Limits

- Orthorhombic boxes are supported in the implemented PBC path.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- L2/L3 quasi_cage growth is bounded for speed and is not an exhaustive enumeration of all possible outer-layer subsets.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
