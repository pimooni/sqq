import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

from sqq.cli import build_parser
from sqq.config import DEFAULT_CONFIG
from sqq.core.hydrate_cluster import analyze_hydrate_clusters, build_cage_graph, find_hydrate_clusters
from sqq.io.summary import display_width, hydrate_cluster_info_section, result_row, summary_sheet_tables
from sqq.models import Cage, Frame, FrameResult, GraphResult, Guest
from sqq.pipeline import apply_cli_overrides, normalize_analysis_scopes


def cage(object_id, cage_type, rings, waters, guests=()):
    return Cage(
        object_id=object_id,
        cage_type=cage_type,
        rings=tuple(rings),
        waters=tuple(waters),
        center=np.zeros(3),
        guest_ids=tuple(guests),
    )


def graph_cages(cage_types, edges):
    ring_lists = [[] for _ in cage_types]
    ring_sizes = {}
    for number, (left, right, size) in enumerate(edges, start=1):
        ring_id = f"ring{size}_{number:05d}"
        ring_lists[left].append(ring_id)
        ring_lists[right].append(ring_id)
        ring_sizes[ring_id] = size
    cages = [
        cage(f"cage_{index + 1:05d}", cage_type, ring_lists[index], [index + 1])
        for index, cage_type in enumerate(cage_types)
    ]
    return cages, ring_sizes


def s_i_graph():
    cage_types = ["512", *(["51262"] * 12)]
    edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
    return graph_cages(cage_types, edges)


def s_ii_graph():
    cage_types = ["512", *(["512"] * 6), *(["51264"] * 6)]
    edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
    return graph_cages(cage_types, edges)


def s_h_graph():
    cage_types = ["51268", "51268", *(["512"] * 6), *(["435663"] * 6)]
    edges = []
    for small in range(2, 8):
        edges.extend([(0, small, 5), (1, small, 5)])
    for offset, medium in enumerate(range(8, 14)):
        anchor = 0 if offset % 2 == 0 else 1
        edges.extend([(anchor, medium, 5), (medium, 2 + offset, 5)])
    edges.append((8, 9, 5))
    return graph_cages(cage_types, edges)


def mixed_graph_with_unclassified_bridge():
    cage_types = [
        "512", *(["51262"] * 12),
        "512", *(["512"] * 6), *(["51264"] * 6),
        "51263",
    ]
    edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
    edges.extend((13, neighbor, 5) for neighbor in range(14, 26))
    edges.extend([(1, 26, 5), (20, 26, 5)])
    return graph_cages(cage_types, edges)


def mixed_graph_with_shared_standard_cage():
    # sI large-cage seed: anchor 0, four small cages (1..4), ten large cages.
    cage_types = ["51262", *(["512"] * 4), *(["51262"] * 10)]
    edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
    edges.extend((0, neighbor, 6) for neighbor in range(13, 15))

    # sII large-cage seed shares small cage 1 with sI.
    s_ii_anchor = len(cage_types)
    cage_types.append("51264")
    additional_small = list(range(len(cage_types), len(cage_types) + 11))
    cage_types.extend(["512"] * 11)
    large_neighbors = list(range(len(cage_types), len(cage_types) + 4))
    cage_types.extend(["51264"] * 4)
    for small in [1, *additional_small]:
        edges.append((s_ii_anchor, small, 5))
    for large in large_neighbors:
        edges.append((s_ii_anchor, large, 6))
    return graph_cages(cage_types, edges)


def frame_result(cages, analysis, *, detail=False):
    clusters, motifs, domains, isolated = analysis
    return FrameResult(
        frame=Frame(name="test", atoms=[], time_ps=0.0, source=Path("test.gro")),
        waters=[],
        guests=[Guest(resid=1, resname="MET", atoms=(1,))],
        graph=GraphResult(mode="hbond", edges=[], adjacency={}),
        rings={},
        cages=cages,
        all_cages=cages,
        cage_report_types=tuple(sorted({item.cage_type for item in cages})),
        hydrate_cluster_enabled=True,
        hydrate_cluster_detail=detail,
        hydrate_clusters=clusters,
        hydrate_motifs=motifs,
        hydrate_domains=domains,
        isolated_cage_ids=isolated,
    )


class HydrateClusterTests(unittest.TestCase):
    def test_shared_ring_face_builds_cluster_and_isolated_cage(self):
        cages = [
            cage("a", "512", ["r1", "shared"], [1, 2, 3], ["MET1"]),
            cage("b", "51262", ["shared", "r3"], [3, 4, 5], ["MET2"]),
            cage("c", "51264", ["r4"], [10, 11, 12]),
        ]
        clusters, isolated = find_hydrate_clusters(cages, min_cage=2)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].cage_ids, ("a", "b"))
        self.assertEqual(clusters[0].water_count, 5)
        self.assertEqual(clusters[0].guest_ids, ("MET1", "MET2"))
        self.assertEqual(clusters[0].hydrate_type, "unclassified")
        self.assertEqual(isolated, ("c",))

    def test_min_cage_threshold_moves_small_components_to_isolated(self):
        cages = [
            cage("a", "512", ["shared"], [1]),
            cage("b", "51262", ["shared"], [2]),
        ]
        clusters, isolated = find_hydrate_clusters(cages, min_cage=3)
        self.assertEqual(clusters, [])
        self.assertEqual(isolated, ("a", "b"))

    def test_multi_owner_ring_selects_one_cage_on_each_face_side(self):
        centers = (
            np.array([0.0, 0.0, 1.0]),
            np.array([0.4, 0.0, 1.0]),
            np.array([0.0, 0.0, -1.0]),
            np.array([0.4, 0.0, -1.0]),
        )
        cages = [
            Cage(
                object_id=f"cage_{index + 1:05d}",
                cage_type="512",
                rings=("ring5_shared",),
                waters=(index,),
                center=center,
            )
            for index, center in enumerate(centers)
        ]
        adjacency, shared_faces, _ = build_cage_graph(
            cages,
            {"ring5_shared": 5},
            face_geometries={
                "ring5_shared": (np.zeros(3), np.array([0.0, 0.0, 1.0]))
            },
        )
        self.assertEqual(adjacency[0], {2})
        self.assertEqual(adjacency[2], {0})
        self.assertEqual(adjacency[1], set())
        self.assertEqual(adjacency[3], set())
        self.assertEqual(shared_faces, {(0, 2): {"ring5_shared"}})

    def test_ideal_phase_graphs_create_domains_without_public_motifs(self):
        for expected_phase, builder, expected_count in (
            ("sI", s_i_graph, 13),
            ("sII", s_ii_graph, 13),
            ("sH", s_h_graph, 14),
        ):
            with self.subTest(phase=expected_phase):
                cages, ring_sizes = builder()
                clusters, motifs, domains, isolated = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
                self.assertEqual(isolated, ())
                self.assertEqual(motifs, [])
                self.assertEqual(clusters[0].hydrate_type, expected_phase)
                self.assertEqual(len(domains), 1)
                self.assertEqual(domains[0].hydrate_type, expected_phase)
                self.assertEqual(domains[0].cage_count, expected_count)
                self.assertEqual(domains[0].seed_count, 1)
                self.assertEqual(domains[0].status, "complete")

    def test_large_cage_network_without_full_seed_is_unclassified(self):
        cages, ring_sizes = graph_cages(
            ["51262", "51262", "51262", "512"],
            [(0, 1, 5), (0, 2, 5), (0, 3, 5), (1, 3, 5), (2, 3, 5)],
        )
        clusters, motifs, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(clusters[0].hydrate_type, "unclassified")
        self.assertEqual(motifs, [])
        self.assertEqual(domains, [])

    def test_two_hexagon_connected_51264_cages_are_not_an_s_ii_domain(self):
        cages, ring_sizes = graph_cages(["51264", "51264"], [(0, 1, 6)])
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(clusters[0].hydrate_type, "unclassified")
        self.assertEqual(domains, [])

    def test_wrong_face_size_is_not_a_domain_seed(self):
        cage_types = ["512", *(["51262"] * 12)]
        edges = [(0, neighbor, 5 if neighbor < 12 else 6) for neighbor in range(1, 13)]
        cages, ring_sizes = graph_cages(cage_types, edges)
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(domains, [])
        self.assertEqual(clusters[0].hydrate_type, "unclassified")

    def test_phase_expansion_requires_seed_and_two_compatible_contacts(self):
        cages, ring_sizes = s_i_graph()
        cage_types = [item.cage_type for item in cages] + ["51262"]
        edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
        edges.extend([(13, 1, 5), (13, 2, 5)])
        cages, ring_sizes = graph_cages(cage_types, edges)
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(clusters[0].hydrate_type, "sI")
        self.assertEqual(domains[0].cage_count, 14)
        self.assertEqual(domains[0].status, "expanded")
        self.assertEqual(len(domains[0].seed_cage_ids), 13)

    def test_s_i_edge_cage_can_contact_51263_boundary(self):
        cage_types = ["512", *(["51262"] * 12), "51262", "51263"]
        edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
        edges.extend([(13, 1, 5), (13, 2, 5), (13, 14, 5)])
        cages, ring_sizes = graph_cages(cage_types, edges)
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(clusters[0].hydrate_type, "sI")
        self.assertEqual(domains[0].cage_count, 14)
        self.assertIn("cage_00014", domains[0].cage_ids)
        self.assertEqual(clusters[0].unclassified_cage_ids, ("cage_00015",))
        self.assertEqual(domains[0].boundary_cage_ids, ("cage_00015",))

    def test_mixed_cluster_keeps_nonstandard_bridge_unclassified(self):
        cages, ring_sizes = mixed_graph_with_unclassified_bridge()
        clusters, motifs, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(motifs, [])
        self.assertEqual(clusters[0].hydrate_type, "mixed")
        self.assertEqual({domain.hydrate_type for domain in domains}, {"sI", "sII"})
        self.assertEqual(clusters[0].unclassified_cage_ids, ("cage_00027",))
        self.assertEqual(clusters[0].phase_boundary_labels, ())

    def test_shared_standard_cage_is_an_s_i_s_ii_boundary(self):
        cages, ring_sizes = mixed_graph_with_shared_standard_cage()
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        cluster = clusters[0]
        self.assertEqual(cluster.hydrate_type, "mixed")
        self.assertEqual({domain.hydrate_type for domain in domains}, {"sI", "sII"})
        self.assertEqual(cluster.phase_boundary_labels, (("cage_00002", "sI+sII"),))
        self.assertEqual(cluster.transition_cage_ids, ("cage_00002",))
        self.assertNotIn("cage_00002", cluster.classified_cage_ids)
        self.assertEqual(
            cluster.cage_count,
            len(cluster.classified_cage_ids) + len(cluster.boundary_cage_ids),
        )
        self.assertEqual(sum(domain.cage_count for domain in domains) + 1, cluster.cage_count)

    def test_same_phase_regions_separated_by_bridge_form_two_domains(self):
        one_types = ["512", *(["51262"] * 12)]
        cage_types = [*one_types, *one_types, "51263"]
        edges = [(0, neighbor, 5) for neighbor in range(1, 13)]
        edges.extend((13, neighbor, 5) for neighbor in range(14, 26))
        edges.extend([(1, 26, 5), (14, 26, 5)])
        cages, ring_sizes = graph_cages(cage_types, edges)
        clusters, motifs, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(motifs, [])
        self.assertEqual(clusters[0].hydrate_type, "sI")
        self.assertEqual(len(domains), 2)
        self.assertEqual(clusters[0].unclassified_cage_ids, ("cage_00027",))

    def test_s_h_requires_adjacent_medium_bridge(self):
        cages, ring_sizes = s_h_graph()
        bridge_id = next(
            ring_id
            for ring_id in ring_sizes
            if ring_id in cages[8].rings and ring_id in cages[9].rings
        )
        cages = [
            Cage(
                object_id=item.object_id,
                cage_type=item.cage_type,
                rings=tuple(ring_id for ring_id in item.rings if ring_id != bridge_id),
                waters=item.waters,
                center=item.center,
                guest_ids=item.guest_ids,
            )
            for item in cages
        ]
        clusters, _, domains, _ = analyze_hydrate_clusters(cages, ring_sizes=ring_sizes)
        self.assertEqual(domains, [])
        self.assertEqual(clusters[0].hydrate_type, "unclassified")

    def test_summary_never_outputs_hydrate_motif_sheet(self):
        cages, ring_sizes = s_i_graph()
        result = frame_result(cages, analyze_hydrate_clusters(cages, ring_sizes=ring_sizes))
        row = result_row(result)
        for detail in (False, True):
            config = deepcopy(DEFAULT_CONFIG)
            config["hydrate_cluster"]["detail"] = detail
            tables = summary_sheet_tables(pd.DataFrame([row]), config)
            self.assertIn("hydrate_cluster", tables)
            self.assertIn("hydrate_domain", tables)
            self.assertNotIn("hydrate_motif", tables)
            if detail:
                self.assertIn("hydrate_cluster_detail", tables)
        domain_table = summary_sheet_tables(pd.DataFrame([row]), deepcopy(DEFAULT_CONFIG))["hydrate_domain"]
        self.assertEqual(domain_table.iloc[0]["seed_count"], 1)
        self.assertEqual(domain_table.iloc[0]["seed_cage_count"], 13)

    def test_info_hierarchy_is_an_aligned_table_without_motif_level(self):
        cages, ring_sizes = s_i_graph()
        result = frame_result(cages, analyze_hydrate_clusters(cages, ring_sizes=ring_sizes))
        text = "\n".join(hydrate_cluster_info_section(result, result_row(result)))
        hierarchy = text.split("## Hydrate Cluster Hierarchy", 1)[1].split("## Hydrate Cluster Detail", 1)[0]
        self.assertIn("| item", hierarchy)
        self.assertIn("| type", hierarchy)
        self.assertIn("| cage_qty", hierarchy)
        self.assertNotIn("| seeds", hierarchy)
        self.assertIn("| cluster_00001", hierarchy)
        self.assertIn("| └ domain_00001", hierarchy)
        self.assertIn("| sI", hierarchy)
        self.assertNotIn("status", hierarchy)
        self.assertNotIn("seed_cages", hierarchy)
        self.assertNotIn("expanded", hierarchy)
        self.assertNotIn("motif_", text)
        self.assertNotIn("Hydrate Motif", text)
        self.assertNotIn("```", text)

        table = []
        for line in [*text.splitlines(), ""]:
            if line.startswith("|"):
                table.append(line)
                continue
            if table:
                self.assertEqual(len({display_width(item) for item in table}), 1)
                table = []

    def test_cli_on_off_min_cage_and_detail_update_config(self):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "-i", "test.gro", "--hydrate-cluster", "on",
            "--cluster-min-cage", "4", "--cluster-detail", "on",
        ])
        config = deepcopy(DEFAULT_CONFIG)
        apply_cli_overrides(config, args)
        normalize_analysis_scopes(config)
        self.assertTrue(config["hydrate_cluster"]["enabled"])
        self.assertEqual(config["hydrate_cluster"]["min_cage"], 4)
        self.assertTrue(config["hydrate_cluster"]["detail"])


if __name__ == "__main__":
    unittest.main()
