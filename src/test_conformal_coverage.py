#!/usr/bin/env python3
"""Deterministic unit tests for conformal support-gated coverage."""

import math
import unittest

import numpy as np

from run_conformal_coverage import conformal_support_gate, selected_order


class ConformalCoverageTests(unittest.TestCase):
    def setUp(self):
        # Compound 0 is the query; 1--4 are references. Reference compounds
        # support one another at 0.8, while the query support is configurable.
        self.similarities = np.eye(5, dtype=float)
        self.similarities[1:, 1:] = 0.8
        np.fill_diagonal(self.similarities, 1.0)
        self.refs = np.asarray([1, 2, 3, 4])

    def test_low_support_invokes_coverage(self):
        self.similarities[0, 1:] = 0.2
        self.similarities[1:, 0] = 0.2
        gate = conformal_support_gate(self.similarities, 0, self.refs, alpha=0.25)
        self.assertTrue(gate["gate_used_coverage"])
        self.assertTrue(math.isclose(gate["support_pvalue"], 0.2))
        self.assertTrue(math.isclose(gate["query_max_similarity"], 0.2))

    def test_supported_query_abstains_to_marginal(self):
        self.similarities[0, 1:] = 0.9
        self.similarities[1:, 0] = 0.9
        gate = conformal_support_gate(self.similarities, 0, self.refs, alpha=0.25)
        self.assertFalse(gate["gate_used_coverage"])
        self.assertTrue(math.isclose(gate["support_pvalue"], 1.0))

    def test_too_few_references_abstains(self):
        gate = conformal_support_gate(
            self.similarities, 0, np.asarray([1]), alpha=0.10
        )
        self.assertFalse(gate["gate_used_coverage"])
        self.assertIsNone(gate["support_pvalue"])

    def test_gate_does_not_accept_labels(self):
        # The public signature itself enforces label-free selection.
        self.assertEqual(
            set(conformal_support_gate.__annotations__),
            {"similarities", "heldout", "refs", "alpha", "return"},
        )

    def test_selector_returns_exact_base_order(self):
        orders = {
            "ORIGINAL_MARGINAL": [0, 1, 2],
            "ORIGINAL_COVERAGE": [2, 1, 0],
        }
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": False}), [0, 1, 2]
        )
        self.assertEqual(
            selected_order(orders, {"gate_used_coverage": True}), [2, 1, 0]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
