#!/usr/bin/env python3
"""Unit tests for the frozen PKIS1 external-validation runner."""

from __future__ import annotations

import unittest

import numpy as np

from run_compiled_coverage import Benchmark, metrics_from_sequence
from run_pkis1_external_validation import (
    MAX_BUDGET,
    MAX_DISPLACEMENT,
    policy_orders,
)


def synthetic_benchmark(labels: np.ndarray) -> Benchmark:
    n, targets = labels.shape
    return Benchmark(
        name="SYNTHETIC",
        compound_keys=[f"C{i:03d}" for i in range(n)],
        canonical_smiles=["C" for _ in range(n)],
        target_keys=[f"T{i:03d}" for i in range(targets)],
        target_names=[f"target-{i}" for i in range(targets)],
        labels={"primary": labels},
        candidate_masks=~np.isnan(labels),
        chemotypes=[frozenset() for _ in range(n)],
        source_ids=[f"S{i:03d}" for i in range(n)],
        metadata={},
    )


class ExternalValidationTests(unittest.TestCase):
    def test_heldout_labels_do_not_change_orders(self) -> None:
        rng = np.random.default_rng(7)
        labels = rng.binomial(1, 0.15, size=(24, 30)).astype(float)
        benchmark = synthetic_benchmark(labels)
        similarities = rng.random((24, 24))
        similarities = (similarities + similarities.T) / 2.0
        np.fill_diagonal(similarities, 1.0)
        refs = np.arange(1, 24)
        candidates = np.arange(30)
        before = policy_orders(benchmark, labels, similarities, 0, refs, candidates)
        changed = labels.copy()
        changed[0] = 1.0 - changed[0]
        after = policy_orders(benchmark, changed, similarities, 0, refs, candidates)
        self.assertEqual(before, after)

    def test_policy_orders_are_deterministic(self) -> None:
        rng = np.random.default_rng(19)
        labels = rng.binomial(1, 0.2, size=(25, 35)).astype(float)
        benchmark = synthetic_benchmark(labels)
        similarities = rng.random((25, 25))
        similarities = (similarities + similarities.T) / 2.0
        np.fill_diagonal(similarities, 1.0)
        refs = np.arange(1, 25)
        candidates = np.arange(35)
        first = policy_orders(benchmark, labels, similarities, 0, refs, candidates)
        second = policy_orders(benchmark, labels, similarities, 0, refs, candidates)
        self.assertEqual(first, second)
        for order in first.values():
            self.assertEqual(len(order), MAX_BUDGET)
            self.assertEqual(len(set(order)), MAX_BUDGET)

    def test_random_first_hit_delay_bound(self) -> None:
        rng = np.random.default_rng(23)
        for _ in range(100):
            references = int(rng.integers(5, 45))
            targets = int(rng.integers(MAX_BUDGET, 75))
            labels = rng.binomial(1, rng.uniform(0.02, 0.35), size=(references + 1, targets)).astype(float)
            benchmark = synthetic_benchmark(labels)
            similarities = rng.random((references + 1, references + 1))
            similarities = (similarities + similarities.T) / 2.0
            np.fill_diagonal(similarities, 1.0)
            refs = np.arange(1, references + 1)
            candidates = np.arange(targets)
            orders = policy_orders(benchmark, labels, similarities, 0, refs, candidates)
            marginal = metrics_from_sequence(orders["MARGINAL"], labels[0], MAX_BUDGET)
            bounded = metrics_from_sequence(
                orders["BOUNDED_COVERAGE"], labels[0], MAX_BUDGET
            )
            self.assertLessEqual(
                bounded["first_hit"] - marginal["first_hit"], MAX_DISPLACEMENT
            )

    def test_missing_reference_values_are_not_positives(self) -> None:
        labels = np.zeros((8, 24), dtype=float)
        labels[1:, 2] = np.nan
        labels[1, 3] = 1.0
        benchmark = synthetic_benchmark(labels)
        similarities = np.eye(8)
        similarities[0, 1:] = 0.5
        similarities[1:, 0] = 0.5
        orders = policy_orders(
            benchmark,
            labels,
            similarities,
            0,
            np.arange(1, 8),
            np.arange(24),
        )
        # The smoothed marginal treats a wholly unobserved target as uncertain
        # (probability 0.5), but the coverage gain must not count its missing
        # reference cells as positives.
        self.assertEqual(orders["MARGINAL"][0], 2)
        self.assertEqual(orders["UNBOUNDED_COVERAGE"][0], 3)


if __name__ == "__main__":
    unittest.main()
