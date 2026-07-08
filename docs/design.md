# SQQ Development Design

SQQ means **Shell Quant Qualifier**. This document records the current implementation logic for developers, so the code and the scientific definitions stay aligned.

## Pipeline

```text
input frames
  -> molecule selection
  -> water graph: hydrogen bond / O-O / user pair map
  -> diagnostic coordination distribution
  -> chordless rings (default) or optional shortest-path rings
  -> half_cage and quasi_cage open patches
  -> closed cage search and guest occupancy
  -> optional hydrate_cluster analysis from reported cages
  -> F3/F4/Q_l plus MCG/DHOP order parameters and ice metrics
  -> per-frame outputs, summary.xlsx, and summary_detail CSV
```

The shared water graph is used by ring, half_cage, quasi_cage, cage, F3/F4, default Q_l, and ice analysis. MCG and DHOP are calculated during the order stage but use dedicated guest/water cutoff graphs because their published definitions are independent of the selected SQQ bond mode. Hydrate_cluster analysis starts after cage reporting and uses reported cage-ring membership, not the raw water graph. The graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode. Coordination diagnostics read this graph without adding, removing, or capping edges.

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

The mode preset controls graph mode, the shared ring-face search sizes, and the automatic worker fraction. It does not control `quasi_cage.max_layers` or output switches. L1 is therefore the default quasi-cage depth in every mode; L2/L3 require `--quasi-max-layer` or an explicit config value.

An explicit `-b` / `--bond-mode {auto,hbond,oo,pairs}` overrides the graph mode from both the preset and `config.yaml`. `--pairs PAIRS.txt` implies pairs mode unless `-b pairs` is already given; it cannot be combined with another explicit bond mode.

`parallel.workers: auto` calculates `floor(physical_core_count * mode_fraction)`, then reserves one physical core for the operating system and caps the result by the number of independent files or selected trajectory frames. Physical-core detection prefers optional `psutil`, then platform probes such as Windows CIM, macOS `sysctl`, or Linux `/proc/cpuinfo`; if physical cores cannot be detected, SQQ falls back to the CPU count visible to the process. `--worker N` / `-w N` overrides the mode fraction: values such as `50%`, `0.5`, or `1` are physical-core fractions, while values greater than one are explicit integer worker counts. The old `--workers` spelling is retained as a hidden compatibility alias. Worker resolution remains capped by task count and the Windows `ProcessPoolExecutor` limit. `parallel.backend` defaults to `process`; `thread` is retained for compatibility and `serial` forces one process.

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
  -> write summary.xlsx, summary_detail CSV, and run_config.yaml
```

`spawn` is selected explicitly on macOS, Windows, and Linux. This avoids forking a process after the interactive progress refresh thread exists and gives the same pickling/import contract on every platform. Worker callables and initializers are module-level functions. Only paths, raw trajectory indexes, small event tuples, and summary dictionaries cross process boundaries; atoms, rings, patches, and cages stay worker-local.

While the pool exists, the parent sets the common BLAS/OpenMP thread environment variables to `parallel.math_threads` (default 1). Spawned children inherit the limits before importing NumPy-backed modules. The parent environment is restored after pool shutdown.

For one XTC/TRR file with `--top`, the parent opens the trajectory once to resolve stride-selected raw indexes. Each worker initializer opens a private MDAnalysis Universe once, and tasks contain a small contiguous batch of `(ordered_frame_index, raw_frame_index)` pairs. Batch size is `ceil(selected_frames / (4 * workers))`, clamped to 1 through 8. Parent and worker readers are explicitly closed. Multiple trajectory files and non-process trajectory backends remain serial. The existing first-three-length orthorhombic box representation is retained; this design does not add triclinic support.

Both standalone-file and indexed-trajectory process paths maintain at most `3 * workers` submitted tasks and refill the queue as futures complete. This is a bounded submission window, not a worker cap: 100 effective workers retain 100-way execution and at most 300 submitted tasks. It avoids constructing a Future and serializing arguments for every item in a very large input set.

Standalone files whose case-insensitive stems collide are rejected because the stem is the frame output-directory name. Non-strict worker-side read failures return a failed summary row; strict failures cancel queued futures and propagate to the main process. Output order is determined by the original index rather than completion order.
### Terminal Progress Display

Serial and parallel runs share the same three-row stage model:

```text
file preparation       reading -> settings -> selecting
core topology search   graph -> ring -> half/quasi -> cage -> cluster
post-processing        filtering -> order -> ice -> output
```

`cluster` is included only when `hydrate_cluster.enabled` is true, for example through `--hydrate-cluster on`. When hydrate cluster analysis is disabled, the stage is omitted rather than shown as `cluster:0`.

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
- `sqq/parallel.py`: spawned worker initialization, file/trajectory tasks, effective CPU detection, and math-thread environment control.
- `sqq/core/graph.py`: water graph construction with orthorhombic PBC and nearby-pair search.
- `sqq/core/ring.py`: non-recursive DFS ring search.
- `sqq/core/ring_topology.py`: shared per-frame ring incidence, adjacency, centers, distances, and optional face-quality metrics.
- `sqq/core/quasi_cage.py`: layered `half_cage` and `quasi_cage` search.
- `sqq/core/cage.py`: closed cage grow/fast-closure search, optional scientific polyhedron validation, and guest assignment.
- `sqq/core/hydrate_cluster.py`: reported-cage cluster graph, phase cores, domains, and boundaries.
- `sqq/core/f3f4.py`: F3/F4/Q_l order metrics.
- `sqq/core/spatial.py`: deterministic orthorhombic-PBC self/cross cutoff pairs for hydrate order metrics.
- `sqq/core/mcg.py`: MCG guest graph and MCG-1/MCG-3 largest components.
- `sqq/core/dhop.py`: DHOP plane-normal counts, water tagging, and DHOP35/DHOP30 components.
- `sqq/core/ice.py`: CHILL-style ice classification.
- `sqq/io/trajectory.py`: standalone readers plus worker-local indexed XTC/TRR frame materialization.
- `sqq/io/summary.py`: per-frame info and global workbook tables.
- `sqq/io/gro_writer.py`: grouped or flat GRO structure output.

## Water-Graph Candidate Search

For `hbond` and `oo` modes, graph construction first requests oxygen pairs within `cutoff + 1e-7` from MDAnalysis `self_capped_distance`. MDAnalysis may select brute force, nsgrid, or periodic KD-tree internally. SQQ treats this only as a candidate generator: every pair is sorted deterministically, recomputed with the existing float64 orthorhombic minimum-image function, compared with the exact configured cutoff, and, in `hbond` mode, checked with the established donor-angle test. If the accelerated neighbor API is unavailable or rejects the input, SQQ uses the previous orthorhombic cell list.

This two-stage design prevents float32 candidate-boundary behavior from changing scientific edges while moving the broad neighbor search into compiled code. `pairs` mode bypasses geometric candidate generation. Non-orthogonal boxes are intentionally outside the current model.
## Coordination Diagnostics

The active graph is summarized by water-node degree. Per-frame outputs report degree 0, 1, 2, 3, 4, and greater than 4 as counts and fractions, together with mean coordination, the degree <=2 fraction, the four-coordinated fraction, and the over-four fraction.

The section title follows the resolved graph mode: Hydrogen-Bond Coordination, O-O Connectivity Coordination, or Pair Connectivity Coordination. These values are diagnostic only. They do not modify graph construction, ring/cage detection, F3/F4, Q_l, or ice classification.

## Order Parameters

F3 and F4 follow the project reference implementation and use the active water graph as the neighbor map.

One frame-local graph-vector cache computes each undirected PBC bond vector once and stores both orientations. F3 and graph-mode Q_l share this cache. Other Q_l neighbor modes build the candidate list once per oxygen pair, and all requested degrees reuse the normalized bond vectors and spherical angles. Cached normalization constants remove repeated factorial work. These are calculation-sharing changes only; neighbor selection, accumulation order, thresholds, and reported values are retained.

Q_l is the local Steinhardt/LAMMPS-style bond-orientational order parameter:

```text
Ybar_lm(i) = (1 / Nb(i)) * sum_j Y_lm(theta_ij, phi_ij)
Q_l(i)     = sqrt(4*pi / (2*l + 1) * sum_m |Ybar_lm(i)|^2)
```

The implementation is independent Python code and does not copy LAMMPS source. It uses unweighted oxygen-neighbor bond vectors and the same rotationally invariant normalization as LAMMPS `compute orientorder/atom`. The default reported degree list is `order.q_degree: [6, 12]`. The same interface can report the common LAMMPS degree list `[4, 6, 8, 10, 12]`.

Neighbor modes:

- `graph` (default): use the active SQQ water graph. In `hbond` mode this gives hydrogen-bond neighbors; in `oo` mode it gives O-O neighbors; in `pairs` mode it follows the user pair map.
- `cutoff`: use all water oxygens within `order.q_cutoff_nm`.
- `nearest`: use the nearest `order.q_n_neighbor` water oxygens within `order.q_cutoff_nm`.
- `lammps`: LAMMPS-compatible cutoff plus fixed-neighbor behavior; if `order.q_n_neighbor` is null, it defaults to `12`.

When a fixed neighbor count is active and fewer than that number of neighbors are found inside the cutoff, every requested Q_l value is set to `0.0`, matching LAMMPS behavior. Without a fixed neighbor count, waters with no Q_l neighbors are omitted from the Q_l mean and count.

Q_l is a continuous structural descriptor, not a standalone ice-count classifier. The value is sensitive to the neighbor definition, so SQQ records the active Q_l degree list and neighbor settings in `run_config.yaml`, the terminal header, and `summary.xlsx`.

## Hydrate Nucleation Order Parameters

MCG/DHOP are frame-local descriptors and are separate from the reported-cage `hydrate_cluster` hierarchy. They are computed in the existing `order` terminal stage and do not introduce another stage. Their results are stored in `HydrateOrderResult`, where `None` means not applicable or disabled and integer zero means that the calculation was applicable but no qualifying component was found.

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

Zero-area planes are skipped and cosine is clamped to [-1,1]. A qualifying pair increments the planar-event count of both central endpoints, which is equivalent to traversing both directed orientations without duplicate geometric work. DHOP35 uses `theta <= 35 degrees`; optional DHOP30 uses `theta <= 30 degrees` from the same loop.

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

Version 0.2.5 retains frame-local caches for symmetric ring-center distances, L2/L3 topology expansions, and patch geometry. These caches are discarded after the frame and therefore never mix topology or coordinates between trajectory frames.

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

The code can also grow from half-cage/quasi-cage patch seeds with `cage.seed_mode = patch`, but ring seeds remain the default. Indexed two-to-four half-cage closure is enabled as a supplementary path. The stricter scientific validation phase is opt-in and therefore does not change the default cage acceptance or center definition.

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
- `cage.report_types` / `--cage-size` filters the user-facing cage counts, occupancy, GRO, info, and workbook tables;
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
4. Pop a state `(face_ids, face_mask, edge_once, edge_twice, face_counts, compatible_target_mask)` from iterative DFS.
5. If `edge_once` is empty, classify by exact face counts and submit the shell to ordinary polyhedron validation.
6. Otherwise choose the open boundary edge with the fewest addable rings (minimum remaining value) using `edge_to_ring_ids`.
7. Reject candidates already in `face_mask`, below single-ring seed rank, intersecting `edge_twice`, or exceeding every compatible target. Before expansion, also reject a target when the remaining face-edge incidence is smaller than the number of open edges or has incompatible parity; these are necessary conditions for closing every open edge exactly once.
8. Adding a ring promotes `ring_edges & edge_once` to `edge_twice` and toggles its remaining edges into `edge_once`; no edge-count dictionary is copied.
9. Update 4/5/6 counts and intersect the compatible-target mask with precomputed count-allowance masks. An empty target mask stops the branch.
10. Calculate ring-center distance only when `cage.max_boundary_candidates` truncates a larger topology-filtered candidate list; record a warning when this occurs.
11. Continue until closure, state budget, or target infeasibility. Face masks are the global duplicate-state keys.

Single-ring seed-rank pruning avoids rediscovering the same cage from every face: later growth cannot add a ring whose stable id ranks before the seed ring.

Default acceptance criteria:

- every edge is used exactly twice;
- Euler characteristic satisfies `V - E + F = 2`;
- face counts match one generated target composition;
- the same `(cage_type, water_set)` was not already accepted.

### Optional Scientific Validation

`cage.scientific_validation = false` is the default and is independent of analysis mode. `--cage-scientific-validation on` enables all of the following additional checks:

- each ordered ring face is locally unwrapped across PBC and fitted by SVD; its planarity RMS must not exceed `cage.max_face_planarity_rms_nm`;
- its cyclic O-O edge-length coefficient of variation must not exceed `cage.max_face_edge_cv`, and its projected area must be nonzero;
- the face-adjacency graph must be edge-connected;
- incident faces around every shell vertex must form one cyclic link, excluding pinched/non-manifold shells;
- the consistently outward-triangulated shell must have volume at least `cage.min_cage_volume_nm3`.

An accepted scientific-validation cage uses the tetrahedral volume centroid of the oriented triangle shell. The default path continues to use the mean of locally unwrapped cage-water coordinates. Consequently, enabling scientific validation may remove distorted cages and may change guest occupancy or geometry-resolved hydrate-cluster edges. Raw ring and half/quasi searches, order parameters, and ice classification are unchanged; ownership-filtered free-ring and free-patch outputs can increase when rejected cages no longer consume them.

The built-in scientific thresholds are `0.06 nm` planarity RMS, `0.35` edge-length CV, and `1.0e-6 nm^3` minimum volume. They are explicit configuration values rather than hidden constants.

## Hydrate Cluster

Hydrate_cluster is optional and defaults to off:

```text
hydrate_cluster.enabled = false
hydrate_cluster.min_cage = 2
hydrate_cluster.detail = false
```

The command-line controls are `--hydrate-cluster on/off`, `--cluster-min-cage N`, and `--cluster-detail on/off`. The detail switch controls the optional per-cluster workbook sheet only; the per-frame Markdown report always includes the compact cluster, hierarchy, domain, and boundary sections when hydrate analysis is enabled.

Hydrate_cluster follows the final reported cage scope. Therefore, `--cage-size I,II --hydrate-cluster on` builds clusters and phase evidence only from reported I/II cages, while `--cage-size all --hydrate-cluster on` uses every reported cage composition in the selected search scope.

The high-level hierarchy is informed by HTR+ ([DOI 10.1088/1361-648X/ad52df](https://doi.org/10.1088/1361-648X/ad52df)): classify hydrate type and polycrystalline boundaries on a cage-connection graph. SQQ does not copy the HTR+ implementation; it uses the explicit shared-face fingerprints and deterministic domain rules below.

### Physical shared-face cage graph

1. Use the reported `result.cages` after `--cage-size` filtering.
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

After all phases collect claims independently, SQQ forms domains from cages claimed by exactly one phase. Domain edges must remain phase-compatible, and every connected domain component must contain at least one strict seed anchor. Same-phase regions separated by boundary cages remain separate domains. Domain ids are deterministic within a frame and are not tracked across frames.

### Boundaries and cluster type

Cages claimed by more than one phase are excluded from exclusive domains. After domains exist, SQQ compares each remaining claimed cage with adjacent domain phases:

- one supported phase produces a single-phase boundary label;
- two or more supported phases produce an interphase label such as `sI+sII`;
- unresolved competing claims are `ambiguous`;
- cages without accepted phase evidence remain `unclassified`.

No cage is forced into a transition role from composition alone. Every non-domain cage belongs to the cluster boundary exactly once. A cluster is `sI`, `sII`, or `sH` when its domains contain one unique phase, `mixed` when multiple domain types occur, and `unclassified` when no domain exists.

Public motif output is not generated in the current release. The compatibility motif return slot remains empty, and neither a `Hydrate Motif` Markdown section nor a `hydrate_motif` workbook sheet is written.

## Guest Occupancy

Cage centers are computed from the locally unwrapped O coordinates of the cage waters. Hydrogens are not used in the cage center.

For each accepted cage, SQQ checks all selected guest molecules. `guest.center_mode` is active:

- `center_atom` (default) and `auto` use the configured center atom when present, such as carbon in methane or carbon dioxide, and otherwise fall back to the geometric residue centroid;
- `centroid` always uses the geometric center of all atoms in the guest molecule.

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

Per-frame output folders use the default grouped structure. Generated GRO paths and title lines use ASCII-only structure labels (`5^12`, `5^126^2`, `hc_5r_5^5`, and similar); Unicode superscripts/subscripts remain display-only in Markdown and Excel. This prevents Windows locale-dependent GRO parsers from failing before they read otherwise valid fixed-width atom records.

Per-frame output folders use the default grouped structure:

```text
frame_name/
  frame_name_info.md
  frame_name_order_parameter.tsv  # only when enabled
  ring/
  half_cage/<type>/
  quasi_cage/<type>/
  cage/<type>/
  ice/
```

The global workbook is `summary.xlsx`. Its first sheet is a dashboard: Configuration includes `SQQ version`, and `Analysis Results (min / mean / max)` reports per-frame min/mean/max values for result metrics while `Frames total / ok / failed` remains a run-level frame count. The workbook also contains per-frame counts, connection and coordination diagnostics, report-scoped ring/cage tables, half_cage/quasi_cage tables, optional per-frame hydrate_cluster totals, order parameters, ice, `detail_index`, and config sheets. Multi-row detail tables are written as UTF-8-SIG CSV files in `summary_detail/`.

Each per-frame `*_info.md` report is optimized for inspection rather than plotting:

- the Ring table shows only report-selected ring sizes and reports final free-ring counts;
- Half Cage and Quasi Cage omit internal `hc_`/`qc_` prefixes, aggregate each composition on a parent row, and list exact isomers on synchronized child rows;
- Cage combines composition totals and structural isomers in one vertical table: each cage composition is a parent row, and observed isomers are synchronized child rows below it;
- Quasi Cage Isomer Description gives one explanation row for each observed quasi-cage isomer, including the base ring and L1/L2/L3 ring sequence;
- Cage Isomer Description gives one explanation row for each observed cage isomer, including face composition and the 6-ring face adjacency pattern;
- Cage Occupancy remains separate because it describes guest assignment, not topology; it uses one cage type per row and dynamic exact guest-composition columns in source guest order;
- Hydrate Cluster appears only when enabled and reports frame totals, a cluster/domain/boundary hierarchy, compact cluster and domain tables, and boundary evidence;
- cluster hierarchy labels include total and nonzero cage-type counts, while cage composition is embedded in the corresponding cluster detail table;
- domain detail separates strict seed counts, unique seed-covered cages, expanded cages, and adjacent boundary contacts; the hierarchy omits a dedicated seed column;
- boundary output separates single-phase, interphase, ambiguous, and unclassified cages without assigning phase from cage composition alone;
- per-frame info reports omit long cage-id and shared-face-id expansions; exact ids remain available in `summary_detail/*.csv`;
- all hierarchy labels use the same short `├`, `└`, and `│` symbols, and Markdown source tables are padded using Unicode display width so their pipe columns align;
- Frame Information starts with `sqq version`, `date & time`, `source`, `frame`, and `time_ps`, then reports bond mode, ring sizes, status, and molecule counts; Molecules, active connection coordination, Order Parameters, Hydrate Nucleation Order Parameters, and Ice are separated into compact sections.

The global `summary.xlsx` workbook keeps plotting-oriented analysis sheets with one input file or trajectory frame per row. Multi-row detail outputs are split into `summary_detail/cage_occupancy.csv`, `summary_detail/cage_isomer.csv`, `summary_detail/hydrate_domain.csv`, and, when `hydrate_cluster.detail = true` or `--cluster-detail on`, `summary_detail/hydrate_cluster_detail.csv`. `detail_index` records the generated CSV file paths and table dimensions. `cage_isomer.csv` defaults to nonzero isomer rows plus per-frame totals; `output.cage_isomer_rows = all` or `--cage-isomer-rows all` restores the full zero-filled matrix. Public motif output is not written. The `order_parameter` sheet reports `F3_mean`, `F3_count`, `F4_mean`, `F4_count`, one mean/count pair for each requested Q_l degree, plus focus-water columns when configured. It also reports default `MCG-1` and `DHOP35` largest clusters; optional `MCG-3` and `DHOP30` columns appear only when enabled. Per-frame Markdown places these four hydrate descriptors in a separate table immediately after F3/F4/Q_l. Optional per-water F3/F4/Q_l output is written as `*_order_parameter.tsv`; MCG/DHOP currently remain frame-level outputs.

## Current Limits

- The implemented PBC path remains orthorhombic. Non-orthogonal/triclinic GRO boxes and trajectory cell angles are not handled by this update.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- Hydrate domains are per-frame topological regions; temporal grain tracking and crystallographic orientation matching are not implemented.
- Hydrate phase classification uses the final reported cage scope, so excluding required cage types or shared-face context with `--cage-size` can prevent the corresponding strict seed or domain from being recognized.
- Boundary labels are per-frame topological evidence; transition-path kinetics, temporal domain tracking, and crystallographic orientation matching are not implemented.
- Default `quasi_cage.search_policy = bounded` is not exhaustive for large outer-layer components. Opt-in `exact` enumerates connected subsets but remains subject to explicit candidate and state limits.
- Automatic process workers parallelize independent GRO/XYZ files or selected frames of one indexed XTC/TRR trajectory. Worker counts are based on physical cores with one physical core reserved for the system; topology search inside one individual frame remains single-process.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
- MCG is meaningful only for guest residue names selected in `hydrate_order.mcg_guest_resnames`; other guest species are not silently treated as methane.
- Published DHOP transition-state thresholds are model- and condition-dependent; SQQ reports the descriptor and does not assign a universal critical-nucleus threshold.
