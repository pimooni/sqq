import unittest

import numpy as np

from sqq.core.f3f4 import normalize_q_degree, q_l_from_vectors, resolve_q_neighbor_count


class Q6Tests(unittest.TestCase):
    def test_simple_cubic_q6_reference(self):
        vectors = [
            np.array([1.0, 0.0, 0.0]),
            np.array([-1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            np.array([0.0, 0.0, -1.0]),
        ]

        self.assertAlmostEqual(q_l_from_vectors(vectors, degree=6), 0.3535533905932738, places=12)

    def test_tetrahedral_q6_reference(self):
        vectors = [
            np.array([1.0, 1.0, 1.0]),
            np.array([1.0, -1.0, -1.0]),
            np.array([-1.0, 1.0, -1.0]),
            np.array([-1.0, -1.0, 1.0]),
        ]

        self.assertAlmostEqual(q_l_from_vectors(vectors, degree=6), 0.6285393610547089, places=12)

    def test_lammps_mode_defaults_to_twelve_neighbors(self):
        self.assertEqual(resolve_q_neighbor_count("lammps", None), 12)

    def test_q_degree_default_and_lammps_list(self):
        self.assertEqual(normalize_q_degree(), (6, 12))
        self.assertEqual(normalize_q_degree("4,6,8,10,12"), (4, 6, 8, 10, 12))


if __name__ == "__main__":
    unittest.main()
