#!/usr/bin/env python3
"""Detailed post-run comparisons for the frozen bounded-coverage component."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from run_compiled_coverage import paired_bootstrap, sign_flip_test, write_csv


MASTER_SEED = 20260905
METRICS = ["audc_1_20", "hit_at_1", "hit_at_3", "hit_at_5", "hit_at_10", "hit_at_20"]
COMPARISONS = {
    "bounded_minus_marginal": ("BOUNDED_COVERAGE", "ORIGINAL_MARGINAL"),
    "bounded_minus_unbounded": ("BOUNDED_COVERAGE", "ORIGINAL_COVERAGE"),
}


def load_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def holm_adjust(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    adjusted = np.empty(len(pvalues), dtype=float)
    running = 0.0
    total = len(pvalues)
    for rank, index in enumerate(order):
        value = min(1.0, (total - rank) * pvalues[index])
        running = max(running, value)
        adjusted[index] = running
    return adjusted.tolist()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-metrics", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    args = parser.parse_args()

    rows = load_rows(args.case_metrics)
    indexed: dict[tuple[str, str], dict[str, dict]] = defaultdict(dict)
    datasets = {}
    for row in rows:
        indexed[(row["condition"], row["method"])][row["compound_key"]] = row
        datasets[row["condition"]] = row["dataset"]

    results = {}
    flat = []
    seed_counter = 0
    for condition in sorted(datasets):
        results[condition] = {}
        for comparison_name, (method, baseline) in COMPARISONS.items():
            method_rows = indexed[(condition, method)]
            baseline_rows = indexed[(condition, baseline)]
            keys = sorted(set(method_rows) & set(baseline_rows))
            metric_results = {}
            for metric in METRICS:
                difference = np.asarray(
                    [
                        float(method_rows[key][metric])
                        - float(baseline_rows[key][metric])
                        for key in keys
                    ]
                )
                result = {
                    "dataset": datasets[condition],
                    "comparison": f"{method} - {baseline}",
                    "metric": metric,
                    "n": len(keys),
                    "estimate": float(difference.mean()),
                    "ci95": paired_bootstrap(
                        difference,
                        args.bootstrap,
                        MASTER_SEED + 10_000 + seed_counter,
                    )["ci95"],
                    "permutation_p_two_sided": sign_flip_test(
                        difference,
                        args.permutations,
                        MASTER_SEED + 20_000 + seed_counter,
                    ),
                }
                metric_results[metric] = result
                flat.append(
                    {
                        "condition": condition,
                        "comparison_name": comparison_name,
                        **result,
                    }
                )
                seed_counter += 1
            results[condition][comparison_name] = metric_results

    for comparison_name in COMPARISONS:
        for metric in METRICS:
            selected = [
                row
                for row in flat
                if row["comparison_name"] == comparison_name and row["metric"] == metric
            ]
            adjusted = holm_adjust(
                [float(row["permutation_p_two_sided"]) for row in selected]
            )
            for row, value in zip(selected, adjusted):
                row["holm_p_across_8_conditions"] = value
                results[row["condition"]][comparison_name][metric][
                    "holm_p_across_8_conditions"
                ] = value

    output = {
        "status": "POST_RUN_DETAILED_COMPONENT_AUDIT",
        "multiplicity": (
            "Holm adjustment is reported separately for each comparison and metric "
            "across the eight conditions. Primary frozen decisions remain AUDC-based."
        ),
        "results": results,
    }
    args.output_json.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(args.output_csv, flat)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
