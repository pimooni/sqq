# SQQ Development Design

SQQ means **Shell Quant Qualifier**. This document records the current implementation logic for developers, so the code and the scientific definitions stay aligned.

## Pipeline

```text
input frames
  -> molecule selection
  -> water graph: hydrogen bond / O-O / user pair map
  -> diagnostic coordination distribution
  -> primitive chordless rings
  -> half_cage and quasi_cage open patches
  -> closed cage search and guest occupancy
  -> F3/F4 and ice metrics
  -> per-frame outputs and summary.xlsx
```

The shared water graph is used by ring, half_cage, quasi_cage, cage, F3/F4, and ice analysis. The graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode. Coordination diagnostics read this graph without adding, removing, or capping edges.

## Analysis Modes and Workers

Modes are discrete base presets, not a continuous 00-99 scale. The default mode is `50`.

| Mode | Label | Graph | Search sizes | Auto worker fraction |
| --- | --- | --- | --- | --- |
| `00` | rigorous | `hbond` | 4/5/6 | 25% |
| `50` | standard | `auto` | 5/6 | 50% |
| `99` | performance | `oo` | 5/6 | 90% |

Mode application order is:

```text
built-in defaults -> mode preset -> config.yaml -> explicit CLI options
```

The mode preset controls graph mode, the shared ring-face search sizes, and the automatic worker fraction. It does not control `quasi_cage.max_layers` or output switches. L1 is therefore the default quasi-cage depth in every mode; L2/L3 require `--quasi-max-layers` or an explicit config value.

An explicit `-b` / `--bond-mode {auto,hbond,oo,pairs}` overrides the graph mode from both the preset and `config.yaml`. `--pairs PAIRS.txt` implies pairs mode unless `-b pairs` is already given; it cannot be combined with another explicit bond mode.

`parallel.workers: auto` calculates `floor(logical_cpu_count * mode_fraction)`, with a minimum of one and a maximum equal to the number of independent input files. `--workers N` overrides that calculation. Parallel execution is file-level and currently uses `ThreadPoolExecutor` for standalone GRO/XYZ files. A single file and XTC/TRR input run with one worker.
Parallel GRO/XYZ runs use a thread-safe progress aggregator. Every worker reports its active file and current pipeline stage. The interactive panel shows completed/failed/active/queued counts, a fixed 11-stage summary in three logical rows, total elapsed time, and up to six active-file rows with stage/file timings. Additional active files are summarized so high-worker modes do not fill the terminal. The serial progress panel remains unchanged.

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

## Coordination Diagnostics

The active graph is summarized by water-node degree. Per-frame outputs report degree 0, 1, 2, 3, 4, and greater than 4 as counts and fractions, together with mean coordination, the degree <=2 fraction, the four-coordinated fraction, and the over-four fraction.

The section title follows the resolved graph mode: Hydrogen-Bond Coordination, O-O Connectivity Coordination, or Pair Connectivity Coordination. These values are diagnostic only. They do not modify graph construction, ring/cage detection, F3/F4, or ice classification.

## Ring Search

Rings are searched on the already-built water graph. The algorithm does not use geometric distance at this stage; it follows graph adjacency. The implementation uses an iterative DFS instead of recursive DFS.

Current behavior:

- supported ring sizes: 4, 5, 6, 7;
- default ring sizes: 5, 6;
- default ring filter: chordless primitive rings;
- ring nodes are water oxygen indices;
- `ring.sizes` / `--size` controls detection, while `ring.report_sizes` / `--ring-size` filters ring tables and GRO files after detection.

## Half-Cage and Quasi-Cage Terms

`patch` means a connected set of ring faces during search. It is not an output class by itself.

Layer definitions:

```text
L0 = base ring
L1 = side rings sharing every base-ring edge; L1 must close into a full side wall
L2 = rings grown outward from L1
L3 = rings grown outward from L2
```

L2 and L3 may be dangling rings or connected dangling ring chains. They do not need to close. The default configuration reports L1 quasi-cages and standard half-cages only; L2/L3 remain available by setting `quasi_cage.max_layers` or `--quasi-max-layers`.

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
8. If `quasi_cage.max_layers >= 2`, grow L2/L3 from the current frontier:
   - find exposed frontier edges;
   - use `edge_to_rings` to get rings sharing those edges;
   - reject rings already in the patch;
   - reject rings touching lower-layer edges, so layer assignment uses the lowest possible layer;
   - sort by ring-center distance and keep at most `quasi_cage.max_layer_candidates`.
9. Convert layer candidates into bounded connected growth units. Small connected components are used whole; large components are represented by single dangling rings and local connected neighborhoods.
10. Each new layered patch is classified again as `half_cage` or `quasi_cage`.

The algorithm intentionally avoids scanning all rings for L1/L2/L3 once `edge_to_rings` is built.

Important limits:

- `quasi_cage.max_layers`: default 1 for fast routine analysis; use 2 or 3 to report outer dangling quasi-cage layers.
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

The code can also grow from half-cage/quasi-cage patch seeds with `cage.seed_mode = patch`, but ring seeds remain the default.

### Search Scope and Report Scope

`ring.sizes` / `--size` defines the shared face-size search universe. Ring and quasi-cage detection support 4/5/6/7. Cage detection intentionally uses only the 4/5/6 intersection of that universe.

For cage search, SQQ generates every trivalent Euler-compatible face composition up to `cage.max_faces` / `--max-cage-faces`:

```text
2*n4 + n5 = 12
n4 + n5 + n6 <= max_faces
```

All generated compositions are searched in one merged grow traversal. Named cage labels are retained when a composition matches:

```text
512    = 5^12
51262  = 5^12 6^2
51263  = 5^12 6^3
51264  = 5^12 6^4
51268  = 5^12 6^8       # Type H large cage
435663 = 4^3 5^6 6^3    # Type H small cage
```

Other accepted compositions use generic labels such as `4^1-5^10-6^2`.

Detection and reporting are separate:

- `all_cages` contains every accepted closed cage in the search scope;
- `cage.report_types` / `--cage-size` filters the user-facing cage counts, occupancy, GRO, info, and workbook tables;
- `--cage-size all` reports every detected composition;
- all detected cages, including unreported types, still remove consumed half-cages, quasi-cages, and free rings;
- an explicitly requested cage type is rejected when it requires a face size absent from `--size` or exceeds `--max-cage-faces`.

The default report set remains `512,51262,51263,51264`.

### Grow Logic

1. Build `edge_to_ring_ids` from all allowed cage rings.
2. Start from a seed face set, usually one ring.
3. Merge all generated target face-count limits into one grow pass so compatible cage types reuse DFS branches.
4. Count how many ring faces use each edge in the current patch.
5. Treat edges used once as open boundaries, twice as closed, and more than twice as invalid.
6. Choose the boundary edge with the fewest addable rings.
7. Use `edge_to_ring_ids` to obtain only rings sharing that boundary edge.
8. Reject a candidate if it is already present, exceeds every compatible target, overuses an edge, or violates seed-rank pruning.
9. Calculate ring-center distance only when `cage.max_boundary_candidates` requires candidate truncation.
10. Continue iterative DFS until the shell closes or no target composition remains feasible.
11. Classify a closed shell by exact face counts, then apply full polyhedron validation.

Single-ring seed-rank pruning avoids rediscovering the same cage from every face: later growth cannot add a ring whose stable id ranks before the seed ring.

Acceptance criteria:

- every edge is used exactly twice;
- Euler characteristic satisfies `V - E + F = 2`;
- face counts match one generated target composition;
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

- patches consumed by any detected cage are not written as free `half_cage` or `quasi_cage`, even when that cage type is not reported;
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

The global workbook is `summary.xlsx`. It contains run metadata, per-frame counts, connection and coordination diagnostics, report-scoped ring/cage tables, half_cage/quasi_cage tables, occupancy tables, F3/F4, ice, and config sheets.

Each per-frame `*_info.md` report is optimized for inspection rather than plotting:

- the Ring table shows only report-selected ring sizes and reports final free-ring counts;
- Half Cage and Quasi Cage omit internal `hc_`/`qc_` prefixes, aggregate each composition on a parent row, and list exact isomers on synchronized child rows;
- Cage and Cage Isomer use one topology type per row;
- Cage Occupancy uses one cage type per row and dynamic exact guest-composition columns in source guest order;
- Frame Information, Molecules, active connection coordination, F3/F4, and Ice are separated into compact sections.

The global `summary.xlsx` workbook is intentionally unchanged: plotting-oriented analysis sheets retain one input file or trajectory frame per row.

## Current Limits

- Orthorhombic boxes are supported in the implemented PBC path.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- Cage-network or crystal-domain classification is not implemented in this version.
- L2/L3 quasi_cage growth is bounded for speed and is not an exhaustive enumeration of all possible outer-layer subsets.
- Automatic workers parallelize independent GRO/XYZ files only; they do not parallelize topology search inside one frame.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
