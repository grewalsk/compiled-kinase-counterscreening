#!/usr/bin/env python3
"""Deterministic unit tests for the compiled-coverage method."""

from __future__ import annotations

import math

import numpy as np

from run_compiled_coverage import (
    coverage_order,
    exact_restricted_order,
    expected_cover_time,
    metrics_from_sequence,
    optimal_batch_partition,
    toy_compiler_boundary_check,
)


def main() -> None:
    # A and B are redundant; C covers a distinct plausible profile.
    labels = np.array(
        [
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ]
    )
    weights = np.array([0.5, 0.3, 0.2])
    candidates = np.arange(3, dtype=int)
    keys = ["A", "B", "C"]
    order = coverage_order(
        labels, weights, candidates, keys, 3, fallback_scores=np.array([0.5, 0.5, 0.3])
    )
    assert order == [0, 2, 1], order

    exact_order, exact_cost = exact_restricted_order(order, labels, weights)
    greedy_cost = expected_cover_time(order, labels, weights)
    assert greedy_cost >= exact_cost - 1e-12
    assert sorted(exact_order) == [0, 1, 2]

    truth = np.array([0.0, 1.0, 0.0])
    metrics = metrics_from_sequence([0, 1, 2], truth, 3)
    assert metrics["first_hit"] == 2
    assert math.isclose(metrics["audc"], 2.0 / 3.0)

    sequential = optimal_batch_partition(order, labels, weights, 0.0)
    assert sequential["batch_sizes"] == [1, 1, 1], sequential
    one_batch = optimal_batch_partition(order, labels, weights, 100.0)
    assert one_batch["batch_sizes"] == [3], one_batch

    boundary = toy_compiler_boundary_check()
    assert boundary["binary_profiles_checked"] == 16
    assert boundary["first_hit_identity_pass"] is True
    assert boundary["multi_hit_counterexample_exists"] is True
    print("PASS: compiled coverage deterministic unit tests")


if __name__ == "__main__":
    main()
