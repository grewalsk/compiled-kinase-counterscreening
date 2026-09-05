#!/usr/bin/env python3
"""Run the frozen exploratory similarity-stratified invariant coverage study."""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import platform
import shutil
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from rdkit import rdBase

from run_compiled_coverage import (
    PKIS2_SHA256,
    coverage_order,
    load_klaeger,
    load_pkis2,
    metrics_from_sequence,
    molecular_arrays,
    normalized_similarity_weights,
    paired_bootstrap,
    reference_indices,
    sha256,
    sign_flip_test,
    stable_rank,
    subset_benchmark,
    weighted_probabilities,
    write_csv,
)
MASTER_SEED = 20260905
MAX_BUDGET = 20
STRATA = 5
PSEUDO_HOLDOUTS = 32
MIN_PSEUDO_HOLDOUTS = 8
STABILITY_REPLICATES = 50
STABILITY_FRACTION = 0.8
ORIGINAL_CHEMOTYPE_GAIN = 0.053671875
CHEMOTYPE_RETENTION_THRESHOLD = 0.04025390625
GATED_RETENTION_THRESHOLD = 0.0268359375
KLAEGER_NONINFERIORITY_MARGIN = 0.005

METHODS = [
    "ORIGINAL_MARGINAL",
    "ORIGINAL_COVERAGE",
    "INVARIANT_MARGINAL",
    "INVARIANT_COVERAGE",
    "SELECTIVE_INVARIANT_COVERAGE",
]

CONDITIONS = {
    "KLAEGER_CHEMBL30": [
        ("klaeger_primary_standard", "primary", "standard"),
        ("klaeger_delta_0_5", "sensitivity", "standard"),
        ("klaeger_scaffold_exclusion", "primary", "scaffold"),
        ("klaeger_leave10_nearest", "primary", "nearest10"),
    ],
    "PKIS2": [
        ("pkis2_gt90_standard", "primary", "standard"),
        ("pkis2_gt80_sensitivity", "sensitivity", "standard"),
        ("pkis2_gt90_scaffold", "primary", "scaffold"),
        ("pkis2_gt90_chemotype", "primary", "chemotype"),
    ],
}

CONDITION_IDS = {
    condition: index + 1
    for index, condition in enumerate(
        [condition for values in CONDITIONS.values() for condition, _, _ in values]
    )
}

COMPARISONS = {
    "original_coverage_minus_original_marginal": (
        "ORIGINAL_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "invariant_coverage_minus_invariant_marginal": (
        "INVARIANT_COVERAGE",
        "INVARIANT_MARGINAL",
    ),
    "invariant_coverage_minus_original_coverage": (
        "INVARIANT_COVERAGE",
        "ORIGINAL_COVERAGE",
    ),
    "selective_invariant_minus_invariant_marginal": (
        "SELECTIVE_INVARIANT_COVERAGE",
        "INVARIANT_MARGINAL",
    ),
}

_WORKER_STATE: tuple | None = None


def similarity_bins(similarities: np.ndarray, strata: int = STRATA) -> list[np.ndarray]:
    if len(similarities) == 0:
        return []
    order = np.argsort(similarities, kind="stable")
    return [
        part.astype(int)
        for part in np.array_split(order, min(strata, len(order)))
        if len(part)
    ]


def group_average_weights(raw_weights: np.ndarray, bins: list[np.ndarray]) -> np.ndarray:
    """Reynolds average of weights under within-bin symmetric groups."""
    averaged = np.zeros_like(raw_weights, dtype=float)
    seen = np.zeros(len(raw_weights), dtype=bool)
    for rows in bins:
        if np.any(seen[rows]):
            raise ValueError("Similarity strata overlap")
        seen[rows] = True
        averaged[rows] = float(raw_weights[rows].sum()) / len(rows)
    if not np.all(seen):
        raise ValueError("Similarity strata do not cover all references")
    if not math.isclose(float(averaged.sum()), 1.0, abs_tol=1e-12):
        raise AssertionError("Invariant weights do not sum to one")
    return averaged


def invariant_weights(similarities: np.ndarray) -> np.ndarray:
    raw = normalized_similarity_weights(similarities)
    return group_average_weights(raw, similarity_bins(similarities, STRATA))


def top_k_jaccard(left: list[int], right: list[int], k: int) -> float:
    left_set = set(left[:k])
    right_set = set(right[:k])
    union = left_set | right_set
    return float(len(left_set & right_set) / len(union)) if union else 1.0


def policy_orders(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
    heldout: int,
    refs: np.ndarray,
    candidates: np.ndarray,
) -> dict[str, list[int]]:
    train = labels[refs]
    sims = similarities[heldout, refs]
    raw_weights = normalized_similarity_weights(sims)
    averaged_weights = group_average_weights(
        raw_weights, similarity_bins(sims, STRATA)
    )
    raw_marginals = weighted_probabilities(raw_weights, train)
    averaged_marginals = weighted_probabilities(averaged_weights, train)
    original_marginal = stable_rank(
        raw_marginals, candidates, benchmark.target_keys
    )[:MAX_BUDGET]
    invariant_marginal = stable_rank(
        averaged_marginals, candidates, benchmark.target_keys
    )[:MAX_BUDGET]
    return {
        "ORIGINAL_MARGINAL": original_marginal,
        "ORIGINAL_COVERAGE": coverage_order(
            train,
            raw_weights,
            candidates,
            benchmark.target_keys,
            MAX_BUDGET,
            raw_marginals,
        ),
        "INVARIANT_MARGINAL": invariant_marginal,
        "INVARIANT_COVERAGE": coverage_order(
            train,
            averaged_weights,
            candidates,
            benchmark.target_keys,
            MAX_BUDGET,
            averaged_marginals,
        ),
    }


def pseudo_holdout_gate(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
    scaffolds: list[str],
    heldout: int,
    outer_refs: np.ndarray,
    exclusion: str,
) -> dict:
    ranked = sorted(
        outer_refs.tolist(), key=lambda index: (-similarities[heldout, index], index)
    )
    pseudo_indices: list[int] = []
    pseudo_deltas: list[float] = []
    for pseudo in ranked:
        allowed = reference_indices(
            pseudo,
            similarities,
            scaffolds,
            benchmark.chemotypes,
            exclusion,
        )
        inner_refs = np.intersect1d(outer_refs, allowed, assume_unique=True)
        candidates = np.flatnonzero(
            benchmark.candidate_masks[pseudo] & ~np.isnan(labels[pseudo])
        )
        if len(inner_refs) < 5 or len(candidates) < MAX_BUDGET:
            continue
        orders = policy_orders(
            benchmark, labels, similarities, pseudo, inner_refs, candidates
        )
        truth = labels[pseudo]
        marginal = metrics_from_sequence(
            orders["INVARIANT_MARGINAL"], truth, MAX_BUDGET
        )
        coverage = metrics_from_sequence(
            orders["INVARIANT_COVERAGE"], truth, MAX_BUDGET
        )
        pseudo_indices.append(pseudo)
        pseudo_deltas.append(float(coverage["audc"] - marginal["audc"]))
        if len(pseudo_indices) == PSEUDO_HOLDOUTS:
            break

    if len(pseudo_indices) < MIN_PSEUDO_HOLDOUTS:
        return {
            "gate_used_coverage": False,
            "pseudo_holdout_count": len(pseudo_indices),
            "estimated_local_uplift": None,
        }
    local_weights = normalized_similarity_weights(
        similarities[heldout, np.asarray(pseudo_indices, dtype=int)]
    )
    estimate = float(local_weights @ np.asarray(pseudo_deltas, dtype=float))
    return {
        "gate_used_coverage": bool(estimate > 0.0),
        "pseudo_holdout_count": len(pseudo_indices),
        "estimated_local_uplift": estimate,
    }


def stability_diagnostic(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
    heldout: int,
    refs: np.ndarray,
    candidates: np.ndarray,
    full_orders: dict[str, list[int]],
    condition: str,
    replicates: int,
) -> dict | None:
    if replicates <= 0:
        return None
    sample_size = max(5, int(math.floor(STABILITY_FRACTION * len(refs))))
    original: list[float] = []
    invariant: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [MASTER_SEED, 71, CONDITION_IDS[condition], heldout, replicate]
            )
        )
        selected = np.sort(rng.choice(refs, size=sample_size, replace=False)).astype(int)
        orders = policy_orders(
            benchmark, labels, similarities, heldout, selected, candidates
        )
        original.append(
            top_k_jaccard(
                full_orders["ORIGINAL_COVERAGE"], orders["ORIGINAL_COVERAGE"], 10
            )
        )
        invariant.append(
            top_k_jaccard(
                full_orders["INVARIANT_COVERAGE"], orders["INVARIANT_COVERAGE"], 10
            )
        )
    return {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        "replicates": replicates,
        "reference_fraction": STABILITY_FRACTION,
        "original_coverage_top10_jaccard": float(np.mean(original)),
        "invariant_coverage_top10_jaccard": float(np.mean(invariant)),
        "invariant_minus_original": float(np.mean(invariant) - np.mean(original)),
    }


def evaluate_outer_case(
    benchmark,
    similarities: np.ndarray,
    scaffolds: list[str],
    condition: str,
    label_name: str,
    exclusion: str,
    heldout: int,
    stability_replicates: int,
) -> tuple[list[dict], dict, dict | None] | None:
    labels = benchmark.labels[label_name]
    refs = reference_indices(
        heldout, similarities, scaffolds, benchmark.chemotypes, exclusion
    )
    candidates = np.flatnonzero(
        benchmark.candidate_masks[heldout] & ~np.isnan(labels[heldout])
    )
    if len(refs) < 5 or len(candidates) < MAX_BUDGET:
        return None
    orders = policy_orders(
        benchmark, labels, similarities, heldout, refs, candidates
    )
    gate = pseudo_holdout_gate(
        benchmark,
        labels,
        similarities,
        scaffolds,
        heldout,
        refs,
        exclusion,
    )
    orders["SELECTIVE_INVARIANT_COVERAGE"] = (
        orders["INVARIANT_COVERAGE"]
        if gate["gate_used_coverage"]
        else orders["INVARIANT_MARGINAL"]
    )
    truth = labels[heldout]
    raw_weights = normalized_similarity_weights(similarities[heldout, refs])
    averaged_weights = invariant_weights(similarities[heldout, refs])
    rows: list[dict] = []
    for method in METHODS:
        metrics = metrics_from_sequence(orders[method], truth, MAX_BUDGET)
        row = {
            "dataset": benchmark.name,
            "condition": condition,
            "compound_key": benchmark.compound_keys[heldout],
            "source_ids": benchmark.source_ids[heldout],
            "method": method,
            "reference_count": len(refs),
            "candidate_count": len(candidates),
            "hidden_positive_count": int(np.nansum(truth[candidates])),
            "audc_1_20": metrics["audc"],
            "cost_to_first_censored": metrics["first_hit"],
            "hit_at_1": metrics["any_curve"][0],
            "hit_at_3": metrics["any_curve"][2],
            "hit_at_5": metrics["any_curve"][4],
            "hit_at_10": metrics["any_curve"][9],
            "hit_at_20": metrics["any_curve"][19],
            "effective_reference_count_original": float(
                1.0 / np.sum(np.square(raw_weights))
            ),
            "effective_reference_count_invariant": float(
                1.0 / np.sum(np.square(averaged_weights))
            ),
            "gate_used_coverage": int(gate["gate_used_coverage"]),
            "gate_pseudo_holdout_count": gate["pseudo_holdout_count"],
            "gate_estimated_local_uplift": gate["estimated_local_uplift"],
            "order": ";".join(benchmark.target_keys[index] for index in orders[method]),
        }
        rows.append(row)
    gate_row = {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        **gate,
    }
    stability = stability_diagnostic(
        benchmark,
        labels,
        similarities,
        heldout,
        refs,
        candidates,
        orders,
        condition,
        stability_replicates,
    )
    return rows, gate_row, stability


def initialize_worker(benchmark, similarities, scaffolds) -> None:
    global _WORKER_STATE
    _WORKER_STATE = (benchmark, similarities, scaffolds)


def evaluate_outer_case_worker(task: tuple) -> tuple[int, object]:
    if _WORKER_STATE is None:
        raise RuntimeError("Worker was not initialized")
    benchmark, similarities, scaffolds = _WORKER_STATE
    condition, label_name, exclusion, heldout, stability_replicates = task
    return heldout, evaluate_outer_case(
        benchmark,
        similarities,
        scaffolds,
        condition,
        label_name,
        exclusion,
        heldout,
        stability_replicates,
    )


def run_condition(
    benchmark,
    similarities: np.ndarray,
    scaffolds: list[str],
    condition: str,
    label_name: str,
    exclusion: str,
    workers: int,
    stability_replicates: int,
) -> tuple[list[dict], list[dict], list[dict]]:
    print(
        f"Running {benchmark.name}/{condition}: {len(benchmark.compound_keys)} cases, "
        f"{workers} workers",
        flush=True,
    )
    start_method = "fork" if "fork" in multiprocessing.get_all_start_methods() else "spawn"
    context = multiprocessing.get_context(start_method)
    tasks = [
        (condition, label_name, exclusion, heldout, stability_replicates)
        for heldout in range(len(benchmark.compound_keys))
    ]
    by_heldout = {}
    with ProcessPoolExecutor(
        max_workers=min(workers, len(tasks)),
        mp_context=context,
        initializer=initialize_worker,
        initargs=(benchmark, similarities, scaffolds),
    ) as pool:
        futures = {pool.submit(evaluate_outer_case_worker, task): task[3] for task in tasks}
        completed = 0
        for future in as_completed(futures):
            heldout, result = future.result()
            by_heldout[heldout] = result
            completed += 1
            if completed % 50 == 0 or completed == len(tasks):
                print(f"  {condition}: {completed}/{len(tasks)}", flush=True)
    case_rows: list[dict] = []
    gate_rows: list[dict] = []
    stability_rows: list[dict] = []
    for heldout in sorted(by_heldout):
        result = by_heldout[heldout]
        if result is None:
            continue
        rows, gate, stability = result
        case_rows.extend(rows)
        gate_rows.append(gate)
        if stability is not None:
            stability_rows.append(stability)
    return case_rows, gate_rows, stability_rows


def summarize(
    case_rows: list[dict],
    gate_rows: list[dict],
    stability_rows: list[dict],
    bootstrap: int,
    permutations: int,
    enforce_reproduction: bool,
) -> dict:
    indexed: dict[tuple[str, str, str], dict[str, dict]] = defaultdict(dict)
    for row in case_rows:
        indexed[(row["dataset"], row["condition"], row["method"])][
            row["compound_key"]
        ] = row
    conditions = sorted({(row["dataset"], row["condition"]) for row in case_rows})
    comparisons: dict[str, dict] = {}
    comparison_rows: list[dict] = []
    seed_counter = 0
    for dataset, condition in conditions:
        comparisons[condition] = {}
        for name, (method, baseline) in COMPARISONS.items():
            method_rows = indexed[(dataset, condition, method)]
            baseline_rows = indexed[(dataset, condition, baseline)]
            keys = sorted(set(method_rows) & set(baseline_rows))
            audc = np.asarray(
                [method_rows[key]["audc_1_20"] - baseline_rows[key]["audc_1_20"] for key in keys]
            )
            assay = np.asarray(
                [
                    baseline_rows[key]["cost_to_first_censored"]
                    - method_rows[key]["cost_to_first_censored"]
                    for key in keys
                ],
                dtype=int,
            )
            boot = paired_bootstrap(audc, bootstrap, MASTER_SEED + 100 + seed_counter)
            pvalue = sign_flip_test(audc, permutations, MASTER_SEED + 200 + seed_counter)
            result = {
                "dataset": dataset,
                "comparison": f"{method} - {baseline}",
                "n": len(keys),
                "estimate": float(audc.mean()),
                "ci95": boot["ci95"],
                "permutation_p_two_sided": pvalue,
                "mean_assays_earlier": float(assay.mean()),
                "wins": int(np.sum(assay > 0)),
                "losses": int(np.sum(assay < 0)),
                "ties": int(np.sum(assay == 0)),
                "large_advance_ge_10": int(np.sum(assay >= 10)),
                "large_delay_ge_10": int(np.sum(assay <= -10)),
            }
            comparisons[condition][name] = result
            comparison_rows.append(
                {"condition": condition, "comparison_name": name, **result}
            )
            seed_counter += 1

    gates: dict[str, dict] = {}
    for condition in sorted({row["condition"] for row in gate_rows}):
        selected = [row for row in gate_rows if row["condition"] == condition]
        estimates = [
            row["estimated_local_uplift"]
            for row in selected
            if row["estimated_local_uplift"] is not None
        ]
        gates[condition] = {
            "n": len(selected),
            "coverage_count": int(sum(row["gate_used_coverage"] for row in selected)),
            "coverage_rate": float(np.mean([row["gate_used_coverage"] for row in selected])),
            "mean_pseudo_holdout_count": float(
                np.mean([row["pseudo_holdout_count"] for row in selected])
            ),
            "mean_estimated_local_uplift": float(np.mean(estimates)) if estimates else None,
        }

    stability = None
    if stability_rows:
        original = np.asarray(
            [row["original_coverage_top10_jaccard"] for row in stability_rows]
        )
        invariant = np.asarray(
            [row["invariant_coverage_top10_jaccard"] for row in stability_rows]
        )
        difference = invariant - original
        stability = {
            "condition": stability_rows[0]["condition"],
            "n": len(stability_rows),
            "subsamples_per_case": stability_rows[0]["replicates"],
            "mean_original_top10_jaccard": float(original.mean()),
            "mean_invariant_top10_jaccard": float(invariant.mean()),
            "mean_difference": float(difference.mean()),
            "difference_ci95": paired_bootstrap(
                difference, bootstrap, MASTER_SEED + 901
            )["ci95"],
        }

    chemotype = comparisons["pkis2_gt90_chemotype"]
    pkis_standard = comparisons["pkis2_gt90_standard"]
    pkis_sensitivity = comparisons["pkis2_gt80_sensitivity"]
    klaeger = comparisons["klaeger_primary_standard"]
    original_chemotype = chemotype["original_coverage_minus_original_marginal"]["estimate"]
    if enforce_reproduction and not math.isclose(
        original_chemotype, ORIGINAL_CHEMOTYPE_GAIN, abs_tol=1e-12
    ):
        raise RuntimeError(
            f"Original chemotype result failed reproduction: {original_chemotype}"
        )
    invariant_chemotype = chemotype["invariant_coverage_minus_invariant_marginal"]["estimate"]
    gated_chemotype = chemotype["selective_invariant_minus_invariant_marginal"]
    standard_original = pkis_standard["original_coverage_minus_original_marginal"]["estimate"]
    standard_invariant = pkis_standard["invariant_coverage_minus_invariant_marginal"]["estimate"]
    sensitivity_original = pkis_sensitivity["original_coverage_minus_original_marginal"]["estimate"]
    sensitivity_invariant = pkis_sensitivity["invariant_coverage_minus_invariant_marginal"]["estimate"]
    klaeger_revision = klaeger["invariant_coverage_minus_original_coverage"]["estimate"]

    success = {
        "chemotype_retains_at_least_75_percent": bool(
            invariant_chemotype >= CHEMOTYPE_RETENTION_THRESHOLD
        ),
        "pkis2_standard_improves_over_original_contrast": bool(
            standard_invariant > standard_original
        ),
        "pkis2_gt80_improves_over_original_contrast": bool(
            sensitivity_invariant > sensitivity_original
        ),
        "pkis2_standard_reversal_eliminated": bool(standard_invariant >= 0.0),
        "pkis2_gt80_reversal_eliminated": bool(sensitivity_invariant >= 0.0),
        "klaeger_primary_noninferior_within_0_005": bool(
            klaeger_revision >= -KLAEGER_NONINFERIORITY_MARGIN
        ),
        "gate_large_delays_at_most_20": bool(gated_chemotype["large_delay_ge_10"] <= 20),
        "gate_retains_at_least_half_original_gain": bool(
            gated_chemotype["estimate"] >= GATED_RETENTION_THRESHOLD
        ),
        "invariant_order_more_stable": bool(
            stability is not None and stability["mean_difference"] > 0.0
        ),
    }
    required = [
        "chemotype_retains_at_least_75_percent",
        "pkis2_standard_improves_over_original_contrast",
        "pkis2_gt80_improves_over_original_contrast",
        "klaeger_primary_noninferior_within_0_005",
        "gate_large_delays_at_most_20",
        "gate_retains_at_least_half_original_gain",
        "invariant_order_more_stable",
    ]
    success["composite_exploratory_success"] = all(success[key] for key in required)
    return {
        "comparisons": comparisons,
        "comparison_rows": comparison_rows,
        "gates": gates,
        "stability": stability,
        "frozen_success_criteria": success,
    }


def make_figures(summary: dict, case_rows: list[dict], stability_rows: list[dict], output: Path) -> None:
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    condition_order = [condition for values in CONDITIONS.values() for condition, _, _ in values]
    labels = [
        "K primary",
        "K delta 0.5",
        "K scaffold",
        "K remove-10",
        "P >90 standard",
        "P >80 standard",
        "P >90 scaffold",
        "P >90 chemotype",
    ]
    original = [
        summary["comparisons"][condition]["original_coverage_minus_original_marginal"]
        for condition in condition_order
    ]
    invariant = [
        summary["comparisons"][condition]["invariant_coverage_minus_invariant_marginal"]
        for condition in condition_order
    ]
    x = np.arange(len(condition_order))
    figure, axis = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    for offset, values, name, color in [
        (-0.14, original, "original", "#8c8c8c"),
        (0.14, invariant, "invariant", "#276678"),
    ]:
        estimates = np.asarray([value["estimate"] for value in values])
        lower = estimates - np.asarray([value["ci95"][0] for value in values])
        upper = np.asarray([value["ci95"][1] for value in values]) - estimates
        axis.errorbar(
            x + offset,
            estimates,
            yerr=np.vstack([lower, upper]),
            fmt="o",
            capsize=3,
            label=name,
            color=color,
        )
    axis.axhline(0.0, color="#222222", linewidth=1)
    axis.set_xticks(x, labels, rotation=30, ha="right")
    axis.set_ylabel("coverage minus matched marginal AUDC")
    axis.set_title("Original and similarity-stratified invariant coverage")
    axis.legend(frameon=False)
    figure.savefig(figures / "condition_effects.png", dpi=180)
    plt.close(figure)

    if stability_rows:
        original_s = [row["original_coverage_top10_jaccard"] for row in stability_rows]
        invariant_s = [row["invariant_coverage_top10_jaccard"] for row in stability_rows]
        figure, axis = plt.subplots(figsize=(6.5, 6), constrained_layout=True)
        axis.scatter(original_s, invariant_s, alpha=0.4, s=20, color="#276678")
        axis.plot([0, 1], [0, 1], color="#333333", linestyle="--")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("original coverage top-10 stability")
        axis.set_ylabel("invariant coverage top-10 stability")
        axis.set_title("PKIS2 chemotype: 80% reference-subsample stability")
        figure.savefig(figures / "stability_comparison.png", dpi=180)
        plt.close(figure)

    selected = [row for row in case_rows if row["condition"] == "pkis2_gt90_chemotype"]
    by_method = defaultdict(dict)
    for row in selected:
        by_method[row["method"]][row["compound_key"]] = row
    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    datasets = []
    display_names = []
    for method, baseline_method, name in [
        ("ORIGINAL_COVERAGE", "ORIGINAL_MARGINAL", "original coverage"),
        ("INVARIANT_COVERAGE", "INVARIANT_MARGINAL", "invariant coverage"),
        (
            "SELECTIVE_INVARIANT_COVERAGE",
            "INVARIANT_MARGINAL",
            "selective invariant",
        ),
    ]:
        baseline = by_method[baseline_method]
        values = [
            baseline[key]["cost_to_first_censored"]
            - by_method[method][key]["cost_to_first_censored"]
            for key in sorted(baseline)
        ]
        datasets.append(values)
        display_names.append(name)
    axis.boxplot(datasets, tick_labels=display_names, showfliers=False)
    axis.axhline(0.0, color="#333333", linewidth=1)
    axis.set_ylabel("assays earlier than invariant marginal")
    axis.set_title("PKIS2 chemotype case-level first-hit differences")
    figure.savefig(figures / "chemotype_tail_risk.png", dpi=180)
    plt.close(figure)


def output_manifest(output: Path) -> dict:
    rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "output_manifest.json":
            rows.append(
                {
                    "path": path.relative_to(output).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    return {"files": rows}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--pkis2-xlsx", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--smoke-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    if sha256(args.pkis2_xlsx) != PKIS2_SHA256:
        raise RuntimeError("PKIS2 SHA-256 mismatch")
    if args.workers < 1:
        raise ValueError("Worker count must be positive")
    bootstrap = min(args.bootstrap, 250) if args.smoke_test else args.bootstrap
    permutations = min(args.permutations, 1_000) if args.smoke_test else args.permutations
    stability_replicates = 3 if args.smoke_test else STABILITY_REPLICATES

    benchmarks = [load_klaeger(args.data_dir), load_pkis2(args.pkis2_xlsx)]
    if args.smoke_test:
        benchmarks = [
            subset_benchmark(
                benchmark, 20 if benchmark.name == "KLAEGER_CHEMBL30" else 40
            )
            for benchmark in benchmarks
        ]

    all_case_rows: list[dict] = []
    all_gate_rows: list[dict] = []
    all_stability_rows: list[dict] = []
    dataset_summary = {}
    for benchmark in benchmarks:
        similarities, scaffolds = molecular_arrays(benchmark)
        dataset_summary[benchmark.name] = {
            **benchmark.metadata,
            "compounds": len(benchmark.compound_keys),
            "targets": len(benchmark.target_keys),
        }
        for condition, label_name, exclusion in CONDITIONS[benchmark.name]:
            condition_stability = (
                stability_replicates if condition == "pkis2_gt90_chemotype" else 0
            )
            case_rows, gate_rows, stability_rows = run_condition(
                benchmark,
                similarities,
                scaffolds,
                condition,
                label_name,
                exclusion,
                args.workers,
                condition_stability,
            )
            all_case_rows.extend(case_rows)
            all_gate_rows.extend(gate_rows)
            all_stability_rows.extend(stability_rows)

    statistical = summarize(
        all_case_rows,
        all_gate_rows,
        all_stability_rows,
        bootstrap,
        permutations,
        enforce_reproduction=not args.smoke_test,
    )
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "case_metrics.csv.gz", all_case_rows)
    write_csv(output / "gate_diagnostics.csv.gz", all_gate_rows)
    write_csv(output / "stability_diagnostics.csv.gz", all_stability_rows)
    write_csv(output / "comparison_summary.csv", statistical.pop("comparison_rows"))
    make_figures(statistical, all_case_rows, all_stability_rows, output)
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "rdkit": rdBase.rdkitVersion,
        "cpu_only": True,
        "paid_api_cost_usd": 0,
        "workers": args.workers,
    }
    (output / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "status": "SMOKE_TEST_ONLY" if args.smoke_test else "EXPLORATORY_MODELING_COMPLETE",
        "interpretation_scope": (
            "Post hoc model-development study motivated by the completed mechanism audit; "
            "not a revised confirmatory test or evidence of therapeutic suitability."
        ),
        "datasets": dataset_summary,
        "models": METHODS,
        "parameters": {
            "strata": STRATA,
            "chemical_weight_beta": 8.0,
            "pseudo_holdouts": PSEUDO_HOLDOUTS,
            "minimum_pseudo_holdouts": MIN_PSEUDO_HOLDOUTS,
            "stability_replicates": stability_replicates,
            "stability_reference_fraction": STABILITY_FRACTION,
            "master_seed": MASTER_SEED,
            "bootstrap": bootstrap,
            "permutations": permutations,
        },
        **statistical,
        "runtime_seconds": time.time() - started,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    protocol = Path(__file__).resolve().parents[1] / "INVARIANT_COVERAGE_PROTOCOL.md"
    if protocol.exists():
        shutil.copy2(protocol, output / protocol.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
