#!/usr/bin/env python3
"""Deterministic unit tests for similarity-stratified invariant coverage."""

import itertools
import math
import unittest

import numpy as np

from run_invariant_coverage import (
    group_average_weights,
    invariant_weights,
    similarity_bins,
    top_k_jaccard,
)


class InvariantCoverageTests(unittest.TestCase):
    def test_group_average_preserves_bin_mass_and_total_mass(self):
        raw = np.asarray([0.05, 0.15, 0.30, 0.10, 0.40])
        bins = [np.asarray([0, 1]), np.asarray([2, 3, 4])]
        averaged = group_average_weights(raw, bins)
        self.assertTrue(math.isclose(float(averaged.sum()), 1.0))
        self.assertTrue(math.isclose(float(averaged[bins[0]].sum()), 0.20))
        self.assertTrue(math.isclose(float(averaged[bins[1]].sum()), 0.80))
        self.assertTrue(math.isclose(averaged[0], averaged[1]))
        self.assertTrue(math.isclose(averaged[2], averaged[3]))
        self.assertTrue(math.isclose(averaged[3], averaged[4]))

    def test_group_average_rejects_overlap_and_incomplete_partition(self):
        raw = np.asarray([0.2, 0.3, 0.5])
        with self.assertRaises(ValueError):
            group_average_weights(raw, [np.asarray([0, 1]), np.asarray([1, 2])])
        with self.assertRaises(ValueError):
            group_average_weights(raw, [np.asarray([0, 1])])

    def test_averaged_objective_equals_exact_permutation_expectation(self):
        raw = np.asarray([0.1, 0.3, 0.2, 0.4])
        bins = [np.asarray([0, 1]), np.asarray([2, 3])]
        hit = np.asarray([1.0, 0.0, 1.0, 0.0])
        averaged = float(group_average_weights(raw, bins) @ hit)
        objectives = []
        for first in itertools.permutations(raw[bins[0]].tolist()):
            for second in itertools.permutations(raw[bins[1]].tolist()):
                permuted = np.asarray(first + second)
                objectives.append(float(permuted @ hit))
        self.assertTrue(math.isclose(averaged, float(np.mean(objectives))))

    def test_invariant_weights_are_constant_inside_similarity_strata(self):
        similarities = np.asarray([0.1, 0.9, 0.3, 0.8, 0.4, 0.7, 0.2])
        weights = invariant_weights(similarities)
        for rows in similarity_bins(similarities, 5):
            self.assertTrue(np.allclose(weights[rows], weights[rows[0]]))
        self.assertTrue(math.isclose(float(weights.sum()), 1.0))

    def test_top_k_jaccard(self):
        self.assertTrue(math.isclose(top_k_jaccard([1, 2, 3], [2, 3, 4], 3), 0.5))
        self.assertTrue(math.isclose(top_k_jaccard([], [], 10), 1.0))


if __name__ == "__main__":
    unittest.main(verbosity=2)
