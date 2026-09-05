#!/usr/bin/env python3
"""Post-run absolute audit of the prespecified bounded-coverage component."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from run_compiled_coverage import paired_bootstrap


MASTER_SEED = 20260905
CHEMOTYPE_RETENTION_THRESHOLD = 0.04025390625
NONINFERIORITY_MARGIN = 0.015
STABILITY_MARGIN = 0.01


def read_gzip_csv(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--stability", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    args = parser.parse_args()

    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    comparisons = summary["comparisons"]
    bounded_key = "bounded_coverage_minus_original_marginal"
    chemotype = comparisons["pkis2_gt90_chemotype"][bounded_key]
    standard = comparisons["pkis2_gt90_standard"][bounded_key]
    sensitivity = comparisons["pkis2_gt80_sensitivity"][bounded_key]
    klaeger = comparisons["klaeger_primary_standard"][bounded_key]

    stability_rows = read_gzip_csv(args.stability)
    stability_difference = np.asarray(
        [
            float(row["bounded_coverage_top10_jaccard"])
            - float(row["original_coverage_top10_jaccard"])
            for row in stability_rows
        ]
    )
    stability = {
        "n": len(stability_difference),
        "subsamples_per_case": int(stability_rows[0]["replicates"]),
        "mean_bounded_top10_jaccard": float(
            np.mean(
                [float(row["bounded_coverage_top10_jaccard"]) for row in stability_rows]
            )
        ),
        "mean_original_top10_jaccard": float(
            np.mean(
                [float(row["original_coverage_top10_jaccard"]) for row in stability_rows]
            )
        ),
        "mean_difference": float(stability_difference.mean()),
        "difference_ci95": paired_bootstrap(
            stability_difference, args.bootstrap, MASTER_SEED + 9901
        )["ci95"],
    }

    bounded_delays = {
        condition: result[bounded_key]["large_delay_ge_10"]
        for condition, result in comparisons.items()
    }
    checks = {
        "chemotype_retains_75_percent": bool(
            chemotype["estimate"] >= CHEMOTYPE_RETENTION_THRESHOLD
            and chemotype["ci95"][0] > 0.0
        ),
        "pkis2_standard_nonnegative_and_noninferior": bool(
            standard["estimate"] >= 0.0
            and standard["ci95"][0] >= -NONINFERIORITY_MARGIN
        ),
        "pkis2_gt80_nonnegative_and_noninferior": bool(
            sensitivity["estimate"] >= 0.0
            and sensitivity["ci95"][0] >= -NONINFERIORITY_MARGIN
        ),
        "klaeger_primary_noninferior": bool(
            klaeger["ci95"][0] >= -NONINFERIORITY_MARGIN
        ),
        "zero_large_delays_all_conditions": bool(
            all(value == 0 for value in bounded_delays.values())
        ),
        "bounded_order_stable": bool(
            stability["mean_difference"] >= 0.0
            and stability["difference_ci95"][0] >= -STABILITY_MARGIN
        ),
    }
    checks["composite_success"] = all(checks.values())
    output = {
        "status": "POST_RUN_PRESPECIFIED_COMPONENT_AUDIT",
        "reason": (
            "BOUNDED_COVERAGE was frozen as a component ablation of the final "
            "risk-controlled model. This audit applies the same absolute criteria "
            "without altering its order or parameters."
        ),
        "comparisons": {
            condition: result[bounded_key]
            for condition, result in comparisons.items()
        },
        "stability": stability,
        "large_delays_by_condition": bounded_delays,
        "absolute_checks": checks,
        "chemotype_fraction_of_original_gain": float(
            chemotype["estimate"]
            / comparisons["pkis2_gt90_chemotype"][
                "original_coverage_minus_original_marginal"
            ]["estimate"]
        ),
        "interpretation_rule": (
            "The component may be described as risk-bounded and empirically robust, "
            "but not as more reference-resampling stable when that criterion fails."
        ),
    }
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
