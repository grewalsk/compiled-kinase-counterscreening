#!/usr/bin/env python3
"""Run the frozen cross-fitted safe-coverage study.

SAFE_COVERAGE selects between the original weighted-marginal and weighted-
coverage policies using only outcomes from eligible pseudo-held-out reference
compounds. The real held-out compound is used only for final evaluation.
"""

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
PSEUDO_HOLDOUTS = 32
MIN_PSEUDO_HOLDOUTS = 8
LCB_Z = 1.2815515655446004
STABILITY_REPLICATES = 10
STABILITY_FRACTION = 0.8
ORIGINAL_CHEMOTYPE_GAIN = 0.053671875
CHEMOTYPE_RETENTION_THRESHOLD = 0.04025390625
NONINFERIORITY_MARGIN = 0.015
STABILITY_MARGIN = 0.01

METHODS = ["ORIGINAL_MARGINAL", "ORIGINAL_COVERAGE", "SAFE_COVERAGE"]

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
    "safe_coverage_minus_original_marginal": (
        "SAFE_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "safe_coverage_minus_original_coverage": (
        "SAFE_COVERAGE",
        "ORIGINAL_COVERAGE",
    ),
}

_WORKER_STATE: tuple | None = None


def top_k_jaccard(left: list[int], right: list[int], k: int) -> float:
    left_set = set(left[:k])
    right_set = set(right[:k])
    union = left_set | right_set
    return float(len(left_set & right_set) / len(union)) if union else 1.0


def base_policy_orders(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
    heldout: int,
    refs: np.ndarray,
    candidates: np.ndarray,
) -> dict[str, list[int]]:
    train = labels[refs]
    weights = normalized_similarity_weights(similarities[heldout, refs])
    marginals = weighted_probabilities(weights, train)
    marginal_order = stable_rank(
        marginals, candidates, benchmark.target_keys
    )[:MAX_BUDGET]
    return {
        "ORIGINAL_MARGINAL": marginal_order,
        "ORIGINAL_COVERAGE": coverage_order(
            train,
            weights,
            candidates,
            benchmark.target_keys,
            MAX_BUDGET,
            marginals,
        ),
    }


def weighted_lcb90(values: np.ndarray, weights: np.ndarray) -> dict:
    """Frozen weighted mean, effective sample size, SE, and one-sided LCB."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if len(values) == 0 or len(values) != len(weights):
        raise ValueError("Values and weights must have the same nonzero length")
    if np.any(weights < 0) or not np.isfinite(weights).all():
        raise ValueError("Weights must be finite and nonnegative")
    total = float(weights.sum())
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Weights must have positive mass")
    normalized = weights / total
    sum_squares = float(np.sum(np.square(normalized)))
    effective_n = 1.0 / sum_squares
    mean = float(normalized @ values)
    if sum_squares >= 1.0 - 1e-15:
        variance = 0.0 if np.allclose(values, mean) else math.nan
    else:
        variance = float(
            np.sum(normalized * np.square(values - mean)) / (1.0 - sum_squares)
        )
    standard_error = (
        math.sqrt(max(0.0, variance) / effective_n)
        if math.isfinite(variance) and effective_n > 1.0
        else math.nan
    )
    lower = mean - LCB_Z * standard_error if math.isfinite(standard_error) else math.nan
    return {
        "estimated_local_uplift": mean,
        "effective_pseudo_holdout_count": effective_n,
        "weighted_variance": variance,
        "weighted_standard_error": standard_error,
        "lcb90": lower,
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
        orders = base_policy_orders(
            benchmark, labels, similarities, pseudo, inner_refs, candidates
        )
        truth = labels[pseudo]
        marginal = metrics_from_sequence(
            orders["ORIGINAL_MARGINAL"], truth, MAX_BUDGET
        )
        coverage = metrics_from_sequence(
            orders["ORIGINAL_COVERAGE"], truth, MAX_BUDGET
        )
        pseudo_indices.append(pseudo)
        pseudo_deltas.append(float(coverage["audc"] - marginal["audc"]))
        if len(pseudo_indices) == PSEUDO_HOLDOUTS:
            break

    empty = {
        "gate_used_coverage": False,
        "pseudo_holdout_count": len(pseudo_indices),
        "effective_pseudo_holdout_count": None,
        "estimated_local_uplift": None,
        "weighted_variance": None,
        "weighted_standard_error": None,
        "lcb90": None,
    }
    if len(pseudo_indices) < MIN_PSEUDO_HOLDOUTS:
        return empty
    local_weights = normalized_similarity_weights(
        similarities[heldout, np.asarray(pseudo_indices, dtype=int)]
    )
    diagnostic = weighted_lcb90(np.asarray(pseudo_deltas), local_weights)
    lower = diagnostic["lcb90"]
    use_coverage = bool(lower is not None and math.isfinite(lower) and lower > 0.0)
    return {
        "gate_used_coverage": use_coverage,
        "pseudo_holdout_count": len(pseudo_indices),
        **diagnostic,
    }


def selected_order(orders: dict[str, list[int]], gate: dict) -> list[int]:
    method = "ORIGINAL_COVERAGE" if gate["gate_used_coverage"] else "ORIGINAL_MARGINAL"
    return orders[method]


def stability_diagnostic(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
    scaffolds: list[str],
    heldout: int,
    refs: np.ndarray,
    candidates: np.ndarray,
    full_orders: dict[str, list[int]],
    full_gate: dict,
    condition: str,
    replicates: int,
) -> dict | None:
    if replicates <= 0:
        return None
    sample_size = max(5, int(math.floor(STABILITY_FRACTION * len(refs))))
    original: list[float] = []
    safe: list[float] = []
    gate_agreement: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [MASTER_SEED, 97, CONDITION_IDS[condition], heldout, replicate]
            )
        )
        selected = np.sort(rng.choice(refs, size=sample_size, replace=False)).astype(int)
        orders = base_policy_orders(
            benchmark, labels, similarities, heldout, selected, candidates
        )
        gate = pseudo_holdout_gate(
            benchmark,
            labels,
            similarities,
            scaffolds,
            heldout,
            selected,
            "chemotype",
        )
        orders["SAFE_COVERAGE"] = selected_order(orders, gate)
        original.append(
            top_k_jaccard(
                full_orders["ORIGINAL_COVERAGE"], orders["ORIGINAL_COVERAGE"], 10
            )
        )
        safe.append(
            top_k_jaccard(full_orders["SAFE_COVERAGE"], orders["SAFE_COVERAGE"], 10)
        )
        gate_agreement.append(
            float(gate["gate_used_coverage"] == full_gate["gate_used_coverage"])
        )
    return {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        "replicates": replicates,
        "reference_fraction": STABILITY_FRACTION,
        "original_coverage_top10_jaccard": float(np.mean(original)),
        "safe_coverage_top10_jaccard": float(np.mean(safe)),
        "safe_minus_original": float(np.mean(safe) - np.mean(original)),
        "gate_decision_agreement": float(np.mean(gate_agreement)),
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
    orders = base_policy_orders(
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
    orders["SAFE_COVERAGE"] = selected_order(orders, gate)
    truth = labels[heldout]
    marginal_metrics = metrics_from_sequence(
        orders["ORIGINAL_MARGINAL"], truth, MAX_BUDGET
    )
    coverage_metrics = metrics_from_sequence(
        orders["ORIGINAL_COVERAGE"], truth, MAX_BUDGET
    )
    realized_delta = float(coverage_metrics["audc"] - marginal_metrics["audc"])
    rows: list[dict] = []
    for method in METHODS:
        metrics = metrics_from_sequence(orders[method], truth, MAX_BUDGET)
        rows.append(
            {
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
                "gate_used_coverage": int(gate["gate_used_coverage"]),
                "order": ";".join(
                    benchmark.target_keys[index] for index in orders[method]
                ),
            }
        )
    predicted = gate["estimated_local_uplift"]
    gate_row = {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        **gate,
        "realized_coverage_minus_marginal": realized_delta,
        "prediction_sign_correct": (
            None
            if predicted is None
            else int((predicted > 0.0) == (realized_delta > 0.0))
        ),
    }
    stability = stability_diagnostic(
        benchmark,
        labels,
        similarities,
        scaffolds,
        heldout,
        refs,
        candidates,
        orders,
        gate,
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


def comparison_result(
    method_rows: dict[str, dict],
    baseline_rows: dict[str, dict],
    bootstrap: int,
    permutations: int,
    seed_offset: int,
    dataset: str,
    method: str,
    baseline: str,
) -> dict:
    keys = sorted(set(method_rows) & set(baseline_rows))
    audc = np.asarray(
        [
            method_rows[key]["audc_1_20"] - baseline_rows[key]["audc_1_20"]
            for key in keys
        ]
    )
    assay = np.asarray(
        [
            baseline_rows[key]["cost_to_first_censored"]
            - method_rows[key]["cost_to_first_censored"]
            for key in keys
        ],
        dtype=int,
    )
    return {
        "dataset": dataset,
        "comparison": f"{method} - {baseline}",
        "n": len(keys),
        "estimate": float(audc.mean()),
        "ci95": paired_bootstrap(
            audc, bootstrap, MASTER_SEED + 1000 + seed_offset
        )["ci95"],
        "permutation_p_two_sided": sign_flip_test(
            audc, permutations, MASTER_SEED + 2000 + seed_offset
        ),
        "mean_assays_earlier": float(assay.mean()),
        "wins": int(np.sum(assay > 0)),
        "losses": int(np.sum(assay < 0)),
        "ties": int(np.sum(assay == 0)),
        "large_advance_ge_10": int(np.sum(assay >= 10)),
        "large_delay_ge_10": int(np.sum(assay <= -10)),
    }


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
            result = comparison_result(
                indexed[(dataset, condition, method)],
                indexed[(dataset, condition, baseline)],
                bootstrap,
                permutations,
                seed_counter,
                dataset,
                method,
                baseline,
            )
            comparisons[condition][name] = result
            comparison_rows.append(
                {"condition": condition, "comparison_name": name, **result}
            )
            seed_counter += 1

    gates: dict[str, dict] = {}
    for condition in sorted({row["condition"] for row in gate_rows}):
        selected = [row for row in gate_rows if row["condition"] == condition]
        valid = [row for row in selected if row["estimated_local_uplift"] is not None]
        used = [row for row in selected if row["gate_used_coverage"]]
        predictions = np.asarray(
            [row["estimated_local_uplift"] for row in valid], dtype=float
        )
        realized = np.asarray(
            [row["realized_coverage_minus_marginal"] for row in valid], dtype=float
        )
        correlation = (
            float(np.corrcoef(predictions, realized)[0, 1])
            if len(valid) > 1
            and np.std(predictions) > 0
            and np.std(realized) > 0
            else None
        )
        gates[condition] = {
            "n": len(selected),
            "valid_gate_count": len(valid),
            "coverage_count": len(used),
            "coverage_rate": float(len(used) / len(selected)),
            "mean_pseudo_holdout_count": float(
                np.mean([row["pseudo_holdout_count"] for row in selected])
            ),
            "mean_lcb90": (
                float(np.mean([row["lcb90"] for row in valid])) if valid else None
            ),
            "prediction_realized_pearson": correlation,
            "prediction_sign_accuracy": (
                float(np.mean([row["prediction_sign_correct"] for row in valid]))
                if valid
                else None
            ),
            "selected_realized_mean_delta": (
                float(np.mean([row["realized_coverage_minus_marginal"] for row in used]))
                if used
                else None
            ),
        }

    stability = None
    if stability_rows:
        original = np.asarray(
            [row["original_coverage_top10_jaccard"] for row in stability_rows]
        )
        safe = np.asarray(
            [row["safe_coverage_top10_jaccard"] for row in stability_rows]
        )
        difference = safe - original
        stability = {
            "condition": stability_rows[0]["condition"],
            "n": len(stability_rows),
            "subsamples_per_case": stability_rows[0]["replicates"],
            "mean_original_top10_jaccard": float(original.mean()),
            "mean_safe_top10_jaccard": float(safe.mean()),
            "mean_difference": float(difference.mean()),
            "difference_ci95": paired_bootstrap(
                difference, bootstrap, MASTER_SEED + 9001
            )["ci95"],
            "mean_gate_decision_agreement": float(
                np.mean([row["gate_decision_agreement"] for row in stability_rows])
            ),
        }

    chemotype = comparisons["pkis2_gt90_chemotype"]
    original_chemotype = chemotype[
        "original_coverage_minus_original_marginal"
    ]["estimate"]
    if enforce_reproduction and not math.isclose(
        original_chemotype, ORIGINAL_CHEMOTYPE_GAIN, abs_tol=1e-12
    ):
        raise RuntimeError(
            f"Original chemotype result failed reproduction: {original_chemotype}"
        )
    safe_chemotype = chemotype["safe_coverage_minus_original_marginal"]
    standard = comparisons["pkis2_gt90_standard"][
        "safe_coverage_minus_original_marginal"
    ]
    sensitivity = comparisons["pkis2_gt80_sensitivity"][
        "safe_coverage_minus_original_marginal"
    ]
    klaeger = comparisons["klaeger_primary_standard"][
        "safe_coverage_minus_original_marginal"
    ]
    success = {
        "chemotype_retains_75_percent": bool(
            safe_chemotype["estimate"] >= CHEMOTYPE_RETENTION_THRESHOLD
            and safe_chemotype["ci95"][0] > 0.0
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
        "chemotype_large_delays_at_most_20": bool(
            safe_chemotype["large_delay_ge_10"] <= 20
        ),
        "safe_order_stable": bool(
            stability is not None
            and stability["mean_difference"] >= 0.0
            and stability["difference_ci95"][0] >= -STABILITY_MARGIN
        ),
    }
    success["composite_success"] = all(success.values())
    return {
        "comparisons": comparisons,
        "comparison_rows": comparison_rows,
        "gates": gates,
        "stability": stability,
        "frozen_absolute_success_criteria": success,
    }


def make_figures(
    summary: dict,
    gate_rows: list[dict],
    stability_rows: list[dict],
    output: Path,
) -> None:
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
        summary["comparisons"][condition][
            "original_coverage_minus_original_marginal"
        ]
        for condition in condition_order
    ]
    safe = [
        summary["comparisons"][condition]["safe_coverage_minus_original_marginal"]
        for condition in condition_order
    ]
    x = np.arange(len(condition_order))
    figure, axis = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    for offset, values, name, color in [
        (-0.14, original, "original coverage", "#8c8c8c"),
        (0.14, safe, "safe coverage", "#1b6f5a"),
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
    axis.set_ylabel("AUDC difference versus original marginal")
    axis.set_title("Absolute performance of cross-fitted safe coverage")
    axis.legend(frameon=False)
    figure.savefig(figures / "condition_effects.png", dpi=180)
    plt.close(figure)

    valid = [
        row for row in gate_rows if row["estimated_local_uplift"] is not None
    ]
    figure, axis = plt.subplots(figsize=(7, 6), constrained_layout=True)
    axis.scatter(
        [row["estimated_local_uplift"] for row in valid],
        [row["realized_coverage_minus_marginal"] for row in valid],
        alpha=0.25,
        s=16,
        color="#1b6f5a",
    )
    axis.axhline(0.0, color="#333333", linewidth=1)
    axis.axvline(0.0, color="#333333", linewidth=1)
    axis.set_xlabel("reference-only predicted coverage uplift")
    axis.set_ylabel("held-out realized coverage uplift")
    axis.set_title("Gate calibration diagnostic (not used for selection)")
    figure.savefig(figures / "gate_calibration.png", dpi=180)
    plt.close(figure)

    if stability_rows:
        figure, axis = plt.subplots(figsize=(6.5, 6), constrained_layout=True)
        axis.scatter(
            [row["original_coverage_top10_jaccard"] for row in stability_rows],
            [row["safe_coverage_top10_jaccard"] for row in stability_rows],
            alpha=0.4,
            s=20,
            color="#1b6f5a",
        )
        axis.plot([0, 1], [0, 1], color="#333333", linestyle="--")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("original coverage top-10 stability")
        axis.set_ylabel("safe coverage top-10 stability")
        axis.set_title("PKIS2 chemotype reference-subsample stability")
        figure.savefig(figures / "stability_comparison.png", dpi=180)
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
    stability_replicates = 2 if args.smoke_test else STABILITY_REPLICATES

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
    make_figures(statistical, all_gate_rows, all_stability_rows, output)
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
        "status": "SMOKE_TEST_ONLY" if args.smoke_test else "FROZEN_MODEL_EVALUATION_COMPLETE",
        "interpretation_scope": (
            "Prospectively frozen post hoc model-development experiment; retrospective "
            "within-panel binding-liability discovery only."
        ),
        "datasets": dataset_summary,
        "models": METHODS,
        "parameters": {
            "chemical_weight_beta": 8.0,
            "pseudo_holdouts": PSEUDO_HOLDOUTS,
            "minimum_pseudo_holdouts": MIN_PSEUDO_HOLDOUTS,
            "lcb_one_sided_confidence": 0.90,
            "lcb_z": LCB_Z,
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
    protocol = Path(__file__).resolve().parents[1] / "SAFE_COVERAGE_PROTOCOL.md"
    if protocol.exists():
        shutil.copy2(protocol, output / protocol.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
