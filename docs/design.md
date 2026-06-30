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
  -> optional hydrate_cluster analysis from reported cages
  -> F3/F4/Q_l order parameters and ice metrics
  -> per-frame outputs and summary.xlsx
```

The shared water graph is used by ring, half_cage, quasi_cage, cage, F3/F4, default Q_l, and ice analysis. Hydrate_cluster analysis starts after cage reporting and uses reported cage-ring membership, not the raw water graph. The graph node is the water oxygen. A graph edge is an O-H...O hydrogen bond in `hbond` mode, an O-O neighbor in `oo` mode, or a user-supplied pair in `pairs` mode. Coordination diagnostics read this graph without adding, removing, or capping edges.

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

`parallel.workers: auto` calculates `floor(logical_cpu_count * mode_fraction)`, with a minimum of one and a maximum equal to the number of independent input files. `--workers N` overrides that calculation. Parallel execution is file-level and currently uses `ThreadPoolExecutor` for standalone GRO/XYZ files. A single file and XTC/TRR input run with one worker.

### Terminal Progress Display

Serial and parallel runs share the same three-row stage model:

```text
file preparation       reading -> settings -> selecting
core topology search   graph -> ring -> half/quasi -> cage -> cluster
post-processing        filtering -> order -> ice -> output
```

`cluster` is included only when `hydrate_cluster.enabled` is true, for example through `--hydrate-cluster on`. When hydrate cluster analysis is disabled, the stage is omitted rather than shown as `cluster:0`.

The serial interactive panel renders the complete workflow and bolds the active stage. Stage columns are sized by the longest stage name in that column, so the display stays compact while `reading`, `graph`, and `filtering` remain aligned. The continuation marker `>` for a new row is placed before the aligned stage column. The timing row remains `stage / frame / total`.

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

Non-interactive serial output, such as redirected logs or a `tqdm` postfix, reports the short current stage label rather than ANSI bold text.

Parallel GRO/XYZ runs use a thread-safe progress aggregator. Every worker reports its active file and current pipeline stage. The interactive parallel panel shows completed/failed/active/queued counts, a compact `stage:count` summary, total elapsed time, and up to six active-file rows with stage/file timings. Additional active files are summarized so high-worker modes do not fill the terminal.

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
- `sqq/core/graph.py`: water graph construction with orthorhombic PBC and nearby-pair search.
- `sqq/core/ring.py`: non-recursive DFS ring search.
- `sqq/core/quasi_cage.py`: layered `half_cage` and `quasi_cage` search.
- `sqq/core/cage.py`: closed cage grow search, polyhedron validation, guest assignment.
- `sqq/core/hydrate_cluster.py`: reported-cage cluster graph, phase cores, domains, and boundaries.
- `sqq/core/f3f4.py`: F3/F4/Q_l order metrics.
- `sqq/core/ice.py`: CHILL-style ice classification.
- `sqq/io/summary.py`: per-frame info and global workbook tables.
- `sqq/io/gro_writer.py`: grouped or flat GRO structure output.

## Coordination Diagnostics

The active graph is summarized by water-node degree. Per-frame outputs report degree 0, 1, 2, 3, 4, and greater than 4 as counts and fractions, together with mean coordination, the degree <=2 fraction, the four-coordinated fraction, and the over-four fraction.

The section title follows the resolved graph mode: Hydrogen-Bond Coordination, O-O Connectivity Coordination, or Pair Connectivity Coordination. These values are diagnostic only. They do not modify graph construction, ring/cage detection, F3/F4, Q_l, or ice classification.

## Order Parameters

F3 and F4 follow the project reference implementation and use the active water graph as the neighbor map.

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

## Hydrate Cluster

Hydrate_cluster is optional and defaults to off:

```text
hydrate_cluster.enabled = false
hydrate_cluster.min_cage = 2
hydrate_cluster.detail = false
```

The command-line controls are `--hydrate-cluster on/off`, `--cluster-min-cage N`, and `--cluster-detail on/off`. The detail switch controls the optional per-cluster workbook sheet only; the per-frame Markdown report always includes the compact cluster, hierarchy, domain, and boundary sections when hydrate analysis is enabled.

Hydrate_cluster follows the final reported cage scope. Therefore, `--cage-size I,II --hydrate-cluster on` builds clusters and phase evidence only from reported I/II cages, while `--cage-size all --hydrate-cluster on` uses every reported cage composition in the selected search scope.

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

The conservative sH seed is a composite of two nonadjacent `5^12 6^8` anchors, exactly six common `5^12` cages connected to both anchors through 5-ring faces, and exactly six `4^3 5^6 6^3` cages. At least one adjacent pair of medium cages must bridge the two anchors.

Seeds are internal phase evidence. `seed_count` is the number of strict seed anchors contained in a domain; `seed_cage_count` is the number of unique domain cages covered by their overlapping seed neighborhoods.

### Expansion and exclusive domains

sI and sII expand independently from the union of their strict seed members. A growth candidate must:

- use a cage type supported by the phase template;
- have at least one compatible internal fingerprint label;
- not exceed any expected fingerprint count by more than one;
- connect through a face label allowed by both endpoint templates; and
- receive at least two compatible contacts from already accepted phase cages.

The characteristic sH composite is already broad, so sH currently expands only by the union of overlapping strict sH seeds.

After all phases collect claims independently, SQQ forms domains from cages claimed by exactly one phase. Domain edges must remain phase-compatible, and every connected domain component must contain at least one strict seed anchor. Same-phase regions separated by boundary cages remain separate domains. Domain ids are deterministic within a frame and are not tracked across frames.

### Boundaries and cluster type

Cages claimed by more than one phase are excluded from exclusive domains. After domains exist, SQQ compares each remaining claimed cage with adjacent domain phases:

- one supported phase produces a single-phase boundary label;
- two or more supported phases produce an interphase label such as `sI+sII`;
- unresolved competing claims are `ambiguous`;
- cages without accepted phase evidence remain `unclassified`.

No cage is forced into a transition role from composition alone. Every non-domain cage belongs to the cluster boundary exactly once. A cluster is `sI`, `sII`, or `sH` when its domains contain one unique phase, `mixed` when multiple domain types occur, and `unclassified` when no domain exists.

Public motif output is not generated in 0.2.2. The compatibility motif return slot remains empty, and neither a `Hydrate Motif` Markdown section nor a `hydrate_motif` workbook sheet is written.

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
  frame_name_order_parameter.tsv  # only when enabled
  ring/
  half_cage/<type>/
  quasi_cage/<type>/
  cage/<type>/
  ice/
```

The global workbook is `summary.xlsx`. It contains run metadata, per-frame counts, connection and coordination diagnostics, report-scoped ring/cage tables, half_cage/quasi_cage tables, occupancy tables, optional hydrate_cluster tables, order parameters, ice, and config sheets.

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
- per-frame info reports omit long cage-id and shared-face-id expansions; exact ids remain available in the optional workbook detail sheets;
- all hierarchy labels use the same short `├`, `└`, and `│` symbols, and Markdown source tables are padded using Unicode display width so their pipe columns align;
- Frame Information starts with `sqq version`, `date & time`, `source`, `frame`, and `time_ps`, then reports bond mode, ring sizes, status, and molecule counts; Molecules, active connection coordination, Order Parameters, and Ice are separated into compact sections.

The global `summary.xlsx` workbook keeps plotting-oriented analysis sheets with one input file or trajectory frame per row. When hydrate_cluster is enabled, `hydrate_cluster` keeps one row per frame and `hydrate_domain` keeps one row per domain. With `hydrate_cluster.detail = true` or `--cluster-detail on`, `hydrate_cluster_detail` adds one row per cluster. Public motif output is not written. The `order_parameter` sheet reports `F3_mean`, `F3_count`, `F4_mean`, `F4_count`, and one mean/count pair for each requested Q_l degree, plus focus-water columns when configured. Optional per-water output is written as `*_order_parameter.tsv`.

## Current Limits

- Orthorhombic boxes are supported in the implemented PBC path.
- Cage detection supports 4/5/6 faces only; 7-member rings remain available for ring and quasi_cage analysis.
- Hydrate domains are per-frame topological regions; temporal grain tracking and crystallographic orientation matching are not implemented.
- Hydrate phase classification uses the final reported cage scope, so excluding required cage types or shared-face context with `--cage-size` can prevent the corresponding strict seed or domain from being recognized.
- Boundary labels are per-frame topological evidence; transition-path kinetics, temporal domain tracking, and crystallographic orientation matching are not implemented.
- L2/L3 quasi_cage growth is bounded for speed and is not an exhaustive enumeration of all possible outer-layer subsets.
- Automatic workers parallelize independent GRO/XYZ files only; they do not parallelize topology search inside one frame.
- CHILL-style ice classification is implemented, but separate atomistic Ih/Ic stacking assignment can be refined later.
