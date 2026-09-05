#!/usr/bin/env python3
"""Deterministic unit tests for cross-fitted safe coverage."""

import math
import unittest

import numpy as np

from run_safe_coverage import LCB_Z, selected_order, top_k_jaccard, weighted_lcb90


class SafeCoverageTests(unittest.TestCase):
    def test_weighted_lcb_matches_equal_weight_formula(self):
        values = np.asarray([0.0, 0.1, 0.2, 0.3])
        weights = np.ones(4)
        result = weighted_lcb90(values, weights)
        expected_variance = float(np.var(values, ddof=1))
        expected_se = math.sqrt(expected_variance / 4)
        self.assertTrue(math.isclose(result["estimated_local_uplift"], 0.15))
        self.assertTrue(math.isclose(result["effective_pseudo_holdout_count"], 4.0))
        self.assertTrue(math.isclose(result["weighted_variance"], expected_variance))
        self.assertTrue(math.isclose(result["weighted_standard_error"], expected_se))
        self.assertTrue(
            math.isclose(result["lcb90"], 0.15 - LCB_Z * expected_se)
        )

    def test_weight_scaling_does_not_change_result(self):
        values = np.asarray([-0.1, 0.2, 0.4])
        first = weighted_lcb90(values, np.asarray([1.0, 2.0, 3.0]))
        second = weighted_lcb90(values, np.asarray([10.0, 20.0, 30.0]))
        for key in first:
            self.assertTrue(math.isclose(first[key], second[key]))

    def test_invalid_weights_are_rejected(self):
        with self.assertRaises(ValueError):
            weighted_lcb90(np.asarray([0.1]), np.asarray([-1.0]))
        with self.assertRaises(ValueError):
            weighted_lcb90(np.asarray([0.1, 0.2]), np.asarray([1.0]))

    def test_selector_returns_exact_frozen_policy(self):
        orders = {
            "ORIGINAL_MARGINAL": [1, 2, 3],
            "ORIGINAL_COVERAGE": [3, 2, 1],
        }
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": False}), [1, 2, 3]
        )
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": True}), [3, 2, 1]
        )

    def test_top_k_jaccard(self):
        self.assertTrue(
            math.isclose(top_k_jaccard([1, 2, 3], [2, 3, 4], 3), 0.5)
        )
        self.assertTrue(math.isclose(top_k_jaccard([], [], 10), 1.0))


if __name__ == "__main__":
    unittest.main(verbosity=2)
