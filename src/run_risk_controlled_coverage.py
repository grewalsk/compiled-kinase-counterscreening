#!/usr/bin/env python3
"""Run frozen conformal displacement-bounded profile coverage."""

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
from run_conformal_coverage import conformal_support_gate, top_k_jaccard


MASTER_SEED = 20260905
MAX_BUDGET = 20
MAX_DISPLACEMENT = 9
STABILITY_REPLICATES = 50
STABILITY_FRACTION = 0.8
ORIGINAL_CHEMOTYPE_GAIN = 0.053671875
CHEMOTYPE_RETENTION_THRESHOLD = 0.04025390625
NONINFERIORITY_MARGIN = 0.015
STABILITY_MARGIN = 0.01

METHODS = [
    "ORIGINAL_MARGINAL",
    "ORIGINAL_COVERAGE",
    "BOUNDED_COVERAGE",
    "RISK_CONTROLLED_COVERAGE",
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
    "bounded_coverage_minus_original_marginal": (
        "BOUNDED_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "risk_controlled_minus_original_marginal": (
        "RISK_CONTROLLED_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "risk_controlled_minus_original_coverage": (
        "RISK_CONTROLLED_COVERAGE",
        "ORIGINAL_COVERAGE",
    ),
}

_WORKER_STATE: tuple | None = None


def bounded_coverage_order(
    labels: np.ndarray,
    weights: np.ndarray,
    candidates: np.ndarray,
    target_keys: list[str],
    max_budget: int,
    marginal_scores: np.ndarray,
    max_displacement: int = MAX_DISPLACEMENT,
) -> list[int]:
    """Greedy profile coverage subject to fixed marginal-rank displacement."""
    if max_displacement < 0:
        raise ValueError("Maximum displacement must be nonnegative")
    full_marginal = stable_rank(marginal_scores, candidates, target_keys)
    marginal_rank = {target: rank for rank, target in enumerate(full_marginal)}
    positive = np.nan_to_num(labels, nan=0.0) == 1.0
    remaining = candidates.copy()
    uncovered = np.ones(len(labels), dtype=bool)
    gains_all = weights @ positive
    order: list[int] = []
    for position in range(min(max_budget, len(remaining))):
        urgent = [
            target
            for target in remaining.tolist()
            if marginal_rank[target] + max_displacement <= position
        ]
        if urgent:
            chosen = min(urgent, key=lambda target: (marginal_rank[target], target_keys[target]))
        else:
            eligible = np.asarray(
                [
                    target
                    for target in remaining.tolist()
                    if marginal_rank[target] <= position + max_displacement
                ],
                dtype=int,
            )
            if len(eligible) == 0:
                raise AssertionError("Displacement window has no eligible target")
            ranked = sorted(
                zip(
                    gains_all[eligible].tolist(),
                    marginal_scores[eligible].tolist(),
                    eligible.tolist(),
                ),
                key=lambda item: (-round(item[0], 14), -item[1], target_keys[item[2]]),
            )
            chosen = ranked[0][2]
        chosen_rank = marginal_rank[chosen]
        if chosen_rank > position + max_displacement:
            raise AssertionError("Target advanced beyond displacement bound")
        if position > chosen_rank + max_displacement:
            raise AssertionError("Target delayed beyond displacement bound")
        order.append(chosen)
        newly_covered = uncovered & positive[:, chosen]
        if np.any(newly_covered):
            gains_all -= weights[newly_covered] @ positive[newly_covered]
            gains_all[np.abs(gains_all) < 1e-14] = 0.0
        uncovered[newly_covered] = False
        remaining = remaining[remaining != chosen]
    return order


def policy_orders(
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
    marginal = stable_rank(marginals, candidates, benchmark.target_keys)[:MAX_BUDGET]
    original = coverage_order(
        train,
        weights,
        candidates,
        benchmark.target_keys,
        MAX_BUDGET,
        marginals,
    )
    bounded = bounded_coverage_order(
        train,
        weights,
        candidates,
        benchmark.target_keys,
        MAX_BUDGET,
        marginals,
        MAX_DISPLACEMENT,
    )
    return {
        "ORIGINAL_MARGINAL": marginal,
        "ORIGINAL_COVERAGE": original,
        "BOUNDED_COVERAGE": bounded,
    }


def selected_order(orders: dict[str, list[int]], gate: dict) -> list[int]:
    method = "BOUNDED_COVERAGE" if gate["gate_used_coverage"] else "ORIGINAL_MARGINAL"
    return orders[method]


def stability_diagnostic(
    benchmark,
    labels: np.ndarray,
    similarities: np.ndarray,
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
    bounded: list[float] = []
    controlled: list[float] = []
    gate_agreement: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [MASTER_SEED, 131, CONDITION_IDS[condition], heldout, replicate]
            )
        )
        selected = np.sort(rng.choice(refs, size=sample_size, replace=False)).astype(int)
        orders = policy_orders(
            benchmark, labels, similarities, heldout, selected, candidates
        )
        gate = conformal_support_gate(similarities, heldout, selected)
        orders["RISK_CONTROLLED_COVERAGE"] = selected_order(orders, gate)
        original.append(
            top_k_jaccard(
                full_orders["ORIGINAL_COVERAGE"], orders["ORIGINAL_COVERAGE"], 10
            )
        )
        bounded.append(
            top_k_jaccard(
                full_orders["BOUNDED_COVERAGE"], orders["BOUNDED_COVERAGE"], 10
            )
        )
        controlled.append(
            top_k_jaccard(
                full_orders["RISK_CONTROLLED_COVERAGE"],
                orders["RISK_CONTROLLED_COVERAGE"],
                10,
            )
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
        "bounded_coverage_top10_jaccard": float(np.mean(bounded)),
        "risk_controlled_top10_jaccard": float(np.mean(controlled)),
        "controlled_minus_original": float(np.mean(controlled) - np.mean(original)),
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
    orders = policy_orders(
        benchmark, labels, similarities, heldout, refs, candidates
    )
    gate = conformal_support_gate(similarities, heldout, refs)
    orders["RISK_CONTROLLED_COVERAGE"] = selected_order(orders, gate)
    truth = labels[heldout]
    metrics_by_method = {
        method: metrics_from_sequence(orders[method], truth, MAX_BUDGET)
        for method in METHODS
    }
    marginal_first = metrics_by_method["ORIGINAL_MARGINAL"]["first_hit"]
    for method in ["BOUNDED_COVERAGE", "RISK_CONTROLLED_COVERAGE"]:
        observed_delay = metrics_by_method[method]["first_hit"] - marginal_first
        if observed_delay > MAX_DISPLACEMENT:
            raise AssertionError(
                f"{method} violated first-hit delay bound: {observed_delay}"
            )
    rows: list[dict] = []
    for method in METHODS:
        metrics = metrics_by_method[method]
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
                "support_pvalue": gate["support_pvalue"],
                "query_max_similarity": gate["query_max_similarity"],
                "order": ";".join(
                    benchmark.target_keys[index] for index in orders[method]
                ),
            }
        )
    gate_row = {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        **gate,
        "realized_bounded_minus_marginal": float(
            metrics_by_method["BOUNDED_COVERAGE"]["audc"]
            - metrics_by_method["ORIGINAL_MARGINAL"]["audc"]
        ),
    }
    stability = stability_diagnostic(
        benchmark,
        labels,
        similarities,
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


def compare(
    method_rows: dict[str, dict],
    baseline_rows: dict[str, dict],
    dataset: str,
    method: str,
    baseline: str,
    bootstrap: int,
    permutations: int,
    seed_offset: int,
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
            audc, bootstrap, MASTER_SEED + 5000 + seed_offset
        )["ci95"],
        "permutation_p_two_sided": sign_flip_test(
            audc, permutations, MASTER_SEED + 6000 + seed_offset
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
            result = compare(
                indexed[(dataset, condition, method)],
                indexed[(dataset, condition, baseline)],
                dataset,
                method,
                baseline,
                bootstrap,
                permutations,
                seed_counter,
            )
            comparisons[condition][name] = result
            comparison_rows.append(
                {"condition": condition, "comparison_name": name, **result}
            )
            seed_counter += 1

    gates: dict[str, dict] = {}
    for condition in sorted({row["condition"] for row in gate_rows}):
        selected = [row for row in gate_rows if row["condition"] == condition]
        valid = [row for row in selected if row["support_pvalue"] is not None]
        used = [row for row in valid if row["gate_used_coverage"]]
        abstained = [row for row in valid if not row["gate_used_coverage"]]
        pvalues = np.asarray([row["support_pvalue"] for row in valid], dtype=float)
        maximums = np.asarray(
            [row["query_max_similarity"] for row in valid], dtype=float
        )
        gates[condition] = {
            "n": len(selected),
            "coverage_count": len(used),
            "coverage_rate": float(len(used) / len(selected)),
            "support_pvalue_quantiles": np.quantile(
                pvalues, [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
            ).tolist(),
            "query_max_similarity_quantiles": np.quantile(
                maximums, [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
            ).tolist(),
            "selected_realized_bounded_mean_delta": (
                float(np.mean([row["realized_bounded_minus_marginal"] for row in used]))
                if used
                else None
            ),
            "abstained_realized_bounded_mean_delta": (
                float(
                    np.mean(
                        [row["realized_bounded_minus_marginal"] for row in abstained]
                    )
                )
                if abstained
                else None
            ),
        }

    stability = None
    if stability_rows:
        original = np.asarray(
            [row["original_coverage_top10_jaccard"] for row in stability_rows]
        )
        bounded = np.asarray(
            [row["bounded_coverage_top10_jaccard"] for row in stability_rows]
        )
        controlled = np.asarray(
            [row["risk_controlled_top10_jaccard"] for row in stability_rows]
        )
        difference = controlled - original
        stability = {
            "condition": stability_rows[0]["condition"],
            "n": len(stability_rows),
            "subsamples_per_case": stability_rows[0]["replicates"],
            "mean_original_top10_jaccard": float(original.mean()),
            "mean_bounded_top10_jaccard": float(bounded.mean()),
            "mean_risk_controlled_top10_jaccard": float(controlled.mean()),
            "mean_difference_controlled_minus_original": float(difference.mean()),
            "difference_ci95": paired_bootstrap(
                difference, bootstrap, MASTER_SEED + 9201
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
    controlled_chemotype = chemotype["risk_controlled_minus_original_marginal"]
    standard = comparisons["pkis2_gt90_standard"][
        "risk_controlled_minus_original_marginal"
    ]
    sensitivity = comparisons["pkis2_gt80_sensitivity"][
        "risk_controlled_minus_original_marginal"
    ]
    klaeger = comparisons["klaeger_primary_standard"][
        "risk_controlled_minus_original_marginal"
    ]
    all_bounded_delays = [
        comparisons[condition]["bounded_coverage_minus_original_marginal"][
            "large_delay_ge_10"
        ]
        for condition in comparisons
    ]
    all_controlled_delays = [
        comparisons[condition]["risk_controlled_minus_original_marginal"][
            "large_delay_ge_10"
        ]
        for condition in comparisons
    ]
    if any(all_bounded_delays) or any(all_controlled_delays):
        raise AssertionError("Observed output violates the displacement guarantee")
    success = {
        "chemotype_retains_75_percent": bool(
            controlled_chemotype["estimate"] >= CHEMOTYPE_RETENTION_THRESHOLD
            and controlled_chemotype["ci95"][0] > 0.0
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
            not any(all_bounded_delays) and not any(all_controlled_delays)
        ),
        "risk_controlled_order_stable": bool(
            stability is not None
            and stability["mean_difference_controlled_minus_original"] >= 0.0
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
        "displacement_guarantee_verified": True,
    }


def make_figures(summary: dict, stability_rows: list[dict], output: Path) -> None:
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
    series = [
        (
            "original_coverage_minus_original_marginal",
            "original coverage",
            "#8c8c8c",
            -0.22,
        ),
        (
            "bounded_coverage_minus_original_marginal",
            "bounded coverage",
            "#bf7c2f",
            0.0,
        ),
        (
            "risk_controlled_minus_original_marginal",
            "risk-controlled coverage",
            "#136f63",
            0.22,
        ),
    ]
    x = np.arange(len(condition_order))
    figure, axis = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    for key, name, color, offset in series:
        values = [summary["comparisons"][condition][key] for condition in condition_order]
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
    axis.set_title("Absolute performance with bounded first-hit risk")
    axis.legend(frameon=False)
    figure.savefig(figures / "condition_effects.png", dpi=180)
    plt.close(figure)

    if stability_rows:
        figure, axis = plt.subplots(figsize=(6.5, 6), constrained_layout=True)
        axis.scatter(
            [row["original_coverage_top10_jaccard"] for row in stability_rows],
            [row["risk_controlled_top10_jaccard"] for row in stability_rows],
            alpha=0.4,
            s=20,
            color="#136f63",
        )
        axis.plot([0, 1], [0, 1], color="#333333", linestyle="--")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("original coverage top-10 stability")
        axis.set_ylabel("risk-controlled top-10 stability")
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
    make_figures(statistical, all_stability_rows, output)
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
            "Final prospectively specified post hoc model-development experiment; "
            "retrospective within-panel binding-liability discovery only."
        ),
        "datasets": dataset_summary,
        "models": METHODS,
        "parameters": {
            "chemical_weight_beta": 8.0,
            "conformal_alpha": 0.10,
            "maximum_marginal_rank_displacement": MAX_DISPLACEMENT,
            "support_score": "nearest eligible-reference Morgan Tanimoto",
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
    protocol = Path(__file__).resolve().parents[1] / "RISK_CONTROLLED_COVERAGE_PROTOCOL.md"
    if protocol.exists():
        shutil.copy2(protocol, output / protocol.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
