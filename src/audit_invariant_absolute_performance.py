#!/usr/bin/env python3
"""Post-run validity audit of absolute invariant-model performance.

The frozen modeling protocol compares invariant coverage with invariant
marginal ranking. This audit checks whether that contrast was caused by an
improved coverage model or a degraded invariant marginal comparator.
"""

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
ORIGINAL_CHEMOTYPE_GAIN = 0.053671875

COMPARISONS = {
    "invariant_marginal_minus_original_marginal": (
        "INVARIANT_MARGINAL",
        "ORIGINAL_MARGINAL",
    ),
    "invariant_coverage_minus_original_marginal": (
        "INVARIANT_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "selective_invariant_minus_original_marginal": (
        "SELECTIVE_INVARIANT_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "selective_invariant_minus_original_coverage": (
        "SELECTIVE_INVARIANT_COVERAGE",
        "ORIGINAL_COVERAGE",
    ),
}


def load_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


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
    flat_rows = []
    seed_counter = 0
    for condition in sorted(datasets):
        results[condition] = {}
        for name, (method, baseline) in COMPARISONS.items():
            method_rows = indexed[(condition, method)]
            baseline_rows = indexed[(condition, baseline)]
            keys = sorted(set(method_rows) & set(baseline_rows))
            audc = np.asarray(
                [
                    float(method_rows[key]["audc_1_20"])
                    - float(baseline_rows[key]["audc_1_20"])
                    for key in keys
                ]
            )
            assay = np.asarray(
                [
                    int(float(baseline_rows[key]["cost_to_first_censored"]))
                    - int(float(method_rows[key]["cost_to_first_censored"]))
                    for key in keys
                ],
                dtype=int,
            )
            boot = paired_bootstrap(
                audc, args.bootstrap, MASTER_SEED + 3000 + seed_counter
            )
            result = {
                "dataset": datasets[condition],
                "comparison": f"{method} - {baseline}",
                "n": len(keys),
                "estimate": float(audc.mean()),
                "ci95": boot["ci95"],
                "permutation_p_two_sided": sign_flip_test(
                    audc, args.permutations, MASTER_SEED + 4000 + seed_counter
                ),
                "mean_assays_earlier": float(assay.mean()),
                "wins": int(np.sum(assay > 0)),
                "losses": int(np.sum(assay < 0)),
                "ties": int(np.sum(assay == 0)),
                "large_advance_ge_10": int(np.sum(assay >= 10)),
                "large_delay_ge_10": int(np.sum(assay <= -10)),
            }
            results[condition][name] = result
            flat_rows.append(
                {"condition": condition, "comparison_name": name, **result}
            )
            seed_counter += 1

    chemotype = results["pkis2_gt90_chemotype"]
    standard = results["pkis2_gt90_standard"]
    sensitivity = results["pkis2_gt80_sensitivity"]
    absolute_checks = {
        "invariant_coverage_retains_75_percent_vs_original_marginal": bool(
            chemotype["invariant_coverage_minus_original_marginal"]["estimate"]
            >= 0.75 * ORIGINAL_CHEMOTYPE_GAIN
        ),
        "selective_invariant_retains_half_vs_original_marginal": bool(
            chemotype["selective_invariant_minus_original_marginal"]["estimate"]
            >= 0.5 * ORIGINAL_CHEMOTYPE_GAIN
        ),
        "selective_invariant_large_delays_vs_original_marginal_at_most_20": bool(
            chemotype["selective_invariant_minus_original_marginal"][
                "large_delay_ge_10"
            ]
            <= 20
        ),
        "pkis2_standard_reversal_eliminated_vs_original_marginal": bool(
            standard["invariant_coverage_minus_original_marginal"]["estimate"]
            >= 0.0
        ),
        "pkis2_gt80_reversal_eliminated_vs_original_marginal": bool(
            sensitivity["invariant_coverage_minus_original_marginal"]["estimate"]
            >= 0.0
        ),
        "selective_model_beats_original_coverage_on_chemotype": bool(
            chemotype["selective_invariant_minus_original_coverage"]["estimate"]
            > 0.0
        ),
    }
    absolute_checks["all_absolute_checks_pass"] = all(absolute_checks.values())
    output = {
        "status": "POST_RUN_VALIDITY_AUDIT",
        "reason": (
            "The frozen within-invariant-family contrast can increase if invariant "
            "weighting harms the invariant marginal baseline. Absolute comparisons "
            "are required before interpreting model improvement."
        ),
        "comparisons": results,
        "absolute_checks": absolute_checks,
        "interpretation_rule": (
            "Do not describe the invariant model as an absolute improvement when it "
            "only improves relative to a degraded matched comparator."
        ),
    }
    args.output_json.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(args.output_csv, flat_rows)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
