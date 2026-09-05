#!/usr/bin/env python3
"""Tests for conformal displacement-bounded profile coverage."""

import unittest

import numpy as np

from run_compiled_coverage import metrics_from_sequence, stable_rank
from run_risk_controlled_coverage import bounded_coverage_order, selected_order


class RiskControlledCoverageTests(unittest.TestCase):
    def test_displacement_and_first_hit_guarantees_random_instances(self):
        for seed in range(50):
            rng = np.random.default_rng(seed)
            labels = (rng.random((30, 40)) < 0.12).astype(float)
            weights = rng.random(30)
            weights /= weights.sum()
            candidates = np.arange(40)
            scores = rng.random(40)
            keys = [f"T{index:03d}" for index in range(40)]
            marginal = stable_rank(scores, candidates, keys)
            bounded = bounded_coverage_order(
                labels, weights, candidates, keys, 20, scores, max_displacement=9
            )
            ranks = {target: rank for rank, target in enumerate(marginal)}
            for position, target in enumerate(bounded):
                self.assertLessEqual(abs(position - ranks[target]), 9)
            truth = (rng.random(40) < 0.15).astype(float)
            marginal_first = metrics_from_sequence(marginal[:20], truth, 20)[
                "first_hit"
            ]
            bounded_first = metrics_from_sequence(bounded, truth, 20)["first_hit"]
            self.assertLessEqual(bounded_first - marginal_first, 9)

    def test_zero_displacement_reproduces_marginal_order(self):
        labels = np.asarray([[1, 0, 1], [0, 1, 0]], dtype=float)
        weights = np.asarray([0.5, 0.5])
        candidates = np.arange(3)
        scores = np.asarray([0.2, 0.9, 0.4])
        keys = ["A", "B", "C"]
        observed = bounded_coverage_order(
            labels, weights, candidates, keys, 3, scores, max_displacement=0
        )
        self.assertEqual(observed, stable_rank(scores, candidates, keys))

    def test_negative_displacement_rejected(self):
        with self.assertRaises(ValueError):
            bounded_coverage_order(
                np.zeros((2, 2)),
                np.asarray([0.5, 0.5]),
                np.arange(2),
                ["A", "B"],
                2,
                np.asarray([0.2, 0.1]),
                max_displacement=-1,
            )

    def test_selector_returns_exact_bounded_or_marginal_order(self):
        orders = {
            "ORIGINAL_MARGINAL": [0, 1, 2],
            "BOUNDED_COVERAGE": [2, 1, 0],
        }
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": False}), [0, 1, 2]
        )
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": True}), [2, 1, 0]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
