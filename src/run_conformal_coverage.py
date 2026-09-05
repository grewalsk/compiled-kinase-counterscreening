#!/usr/bin/env python3
"""Run the frozen label-free conformal support-gated coverage study."""

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
    load_klaeger,
    load_pkis2,
    metrics_from_sequence,
    molecular_arrays,
    paired_bootstrap,
    reference_indices,
    sha256,
    sign_flip_test,
    subset_benchmark,
    write_csv,
)
from run_safe_coverage import base_policy_orders, top_k_jaccard


MASTER_SEED = 20260905
MAX_BUDGET = 20
CONFORMAL_ALPHA = 0.10
STABILITY_REPLICATES = 50
STABILITY_FRACTION = 0.8
ORIGINAL_CHEMOTYPE_GAIN = 0.053671875
CHEMOTYPE_RETENTION_THRESHOLD = 0.04025390625
NONINFERIORITY_MARGIN = 0.015
STABILITY_MARGIN = 0.01

METHODS = ["ORIGINAL_MARGINAL", "ORIGINAL_COVERAGE", "CONFORMAL_COVERAGE"]

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
    "conformal_coverage_minus_original_marginal": (
        "CONFORMAL_COVERAGE",
        "ORIGINAL_MARGINAL",
    ),
    "conformal_coverage_minus_original_coverage": (
        "CONFORMAL_COVERAGE",
        "ORIGINAL_COVERAGE",
    ),
}

_WORKER_STATE: tuple | None = None


def conformal_support_gate(
    similarities: np.ndarray,
    heldout: int,
    refs: np.ndarray,
    alpha: float = CONFORMAL_ALPHA,
) -> dict:
    """Return the frozen lower-tail conformal chemical-support decision."""
    refs = np.asarray(refs, dtype=int)
    if len(refs) < 2:
        return {
            "gate_used_coverage": False,
            "reference_count": len(refs),
            "query_max_similarity": None,
            "support_pvalue": None,
            "calibration_support_min": None,
            "calibration_support_median": None,
            "calibration_support_max": None,
        }
    query_support = float(np.max(similarities[heldout, refs]))
    calibration_matrix = similarities[np.ix_(refs, refs)].copy()
    np.fill_diagonal(calibration_matrix, -np.inf)
    calibration = np.max(calibration_matrix, axis=1)
    if not math.isfinite(query_support) or not np.isfinite(calibration).all():
        return {
            "gate_used_coverage": False,
            "reference_count": len(refs),
            "query_max_similarity": query_support,
            "support_pvalue": None,
            "calibration_support_min": None,
            "calibration_support_median": None,
            "calibration_support_max": None,
        }
    pvalue = float((1 + np.sum(calibration <= query_support)) / (len(refs) + 1))
    return {
        "gate_used_coverage": bool(pvalue <= alpha),
        "reference_count": len(refs),
        "query_max_similarity": query_support,
        "support_pvalue": pvalue,
        "calibration_support_min": float(np.min(calibration)),
        "calibration_support_median": float(np.median(calibration)),
        "calibration_support_max": float(np.max(calibration)),
    }


def selected_order(orders: dict[str, list[int]], gate: dict) -> list[int]:
    method = "ORIGINAL_COVERAGE" if gate["gate_used_coverage"] else "ORIGINAL_MARGINAL"
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
    conformal: list[float] = []
    gate_agreement: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [MASTER_SEED, 113, CONDITION_IDS[condition], heldout, replicate]
            )
        )
        selected = np.sort(rng.choice(refs, size=sample_size, replace=False)).astype(int)
        orders = base_policy_orders(
            benchmark, labels, similarities, heldout, selected, candidates
        )
        gate = conformal_support_gate(similarities, heldout, selected)
        orders["CONFORMAL_COVERAGE"] = selected_order(orders, gate)
        original.append(
            top_k_jaccard(
                full_orders["ORIGINAL_COVERAGE"], orders["ORIGINAL_COVERAGE"], 10
            )
        )
        conformal.append(
            top_k_jaccard(
                full_orders["CONFORMAL_COVERAGE"],
                orders["CONFORMAL_COVERAGE"],
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
        "conformal_coverage_top10_jaccard": float(np.mean(conformal)),
        "conformal_minus_original": float(np.mean(conformal) - np.mean(original)),
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
    gate = conformal_support_gate(similarities, heldout, refs)
    orders["CONFORMAL_COVERAGE"] = selected_order(orders, gate)
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
        "realized_coverage_minus_marginal": realized_delta,
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
            audc, bootstrap, MASTER_SEED + 3000 + seed_offset
        )["ci95"],
        "permutation_p_two_sided": sign_flip_test(
            audc, permutations, MASTER_SEED + 4000 + seed_offset
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
            "valid_gate_count": len(valid),
            "coverage_count": len(used),
            "coverage_rate": float(len(used) / len(selected)),
            "support_pvalue_quantiles": (
                np.quantile(pvalues, [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]).tolist()
                if valid
                else None
            ),
            "query_max_similarity_quantiles": (
                np.quantile(maximums, [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]).tolist()
                if valid
                else None
            ),
            "selected_realized_mean_delta": (
                float(np.mean([row["realized_coverage_minus_marginal"] for row in used]))
                if used
                else None
            ),
            "abstained_realized_mean_delta": (
                float(
                    np.mean(
                        [row["realized_coverage_minus_marginal"] for row in abstained]
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
        conformal = np.asarray(
            [row["conformal_coverage_top10_jaccard"] for row in stability_rows]
        )
        difference = conformal - original
        stability = {
            "condition": stability_rows[0]["condition"],
            "n": len(stability_rows),
            "subsamples_per_case": stability_rows[0]["replicates"],
            "mean_original_top10_jaccard": float(original.mean()),
            "mean_conformal_top10_jaccard": float(conformal.mean()),
            "mean_difference": float(difference.mean()),
            "difference_ci95": paired_bootstrap(
                difference, bootstrap, MASTER_SEED + 9101
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
    conformal_chemotype = chemotype[
        "conformal_coverage_minus_original_marginal"
    ]
    standard = comparisons["pkis2_gt90_standard"][
        "conformal_coverage_minus_original_marginal"
    ]
    sensitivity = comparisons["pkis2_gt80_sensitivity"][
        "conformal_coverage_minus_original_marginal"
    ]
    klaeger = comparisons["klaeger_primary_standard"][
        "conformal_coverage_minus_original_marginal"
    ]
    success = {
        "chemotype_retains_75_percent": bool(
            conformal_chemotype["estimate"] >= CHEMOTYPE_RETENTION_THRESHOLD
            and conformal_chemotype["ci95"][0] > 0.0
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
            conformal_chemotype["large_delay_ge_10"] <= 20
        ),
        "conformal_order_stable": bool(
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
    conformal = [
        summary["comparisons"][condition][
            "conformal_coverage_minus_original_marginal"
        ]
        for condition in condition_order
    ]
    x = np.arange(len(condition_order))
    figure, axis = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    for offset, values, name, color in [
        (-0.14, original, "original coverage", "#8c8c8c"),
        (0.14, conformal, "conformal coverage", "#7a3e9d"),
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
    axis.set_title("Absolute performance of conformal support-gated coverage")
    axis.legend(frameon=False)
    figure.savefig(figures / "condition_effects.png", dpi=180)
    plt.close(figure)

    valid = [row for row in gate_rows if row["support_pvalue"] is not None]
    figure, axis = plt.subplots(figsize=(7, 6), constrained_layout=True)
    colors = ["#7a3e9d" if row["gate_used_coverage"] else "#aaaaaa" for row in valid]
    axis.scatter(
        [row["support_pvalue"] for row in valid],
        [row["realized_coverage_minus_marginal"] for row in valid],
        alpha=0.25,
        s=16,
        c=colors,
    )
    axis.axhline(0.0, color="#333333", linewidth=1)
    axis.axvline(CONFORMAL_ALPHA, color="#7a3e9d", linestyle="--")
    axis.set_xlabel("label-free conformal support p-value")
    axis.set_ylabel("held-out realized coverage uplift")
    axis.set_title("Chemical support and realized coverage benefit")
    figure.savefig(figures / "support_diagnostic.png", dpi=180)
    plt.close(figure)

    if stability_rows:
        figure, axis = plt.subplots(figsize=(6.5, 6), constrained_layout=True)
        axis.scatter(
            [row["original_coverage_top10_jaccard"] for row in stability_rows],
            [row["conformal_coverage_top10_jaccard"] for row in stability_rows],
            alpha=0.4,
            s=20,
            color="#7a3e9d",
        )
        axis.plot([0, 1], [0, 1], color="#333333", linestyle="--")
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("original coverage top-10 stability")
        axis.set_ylabel("conformal coverage top-10 stability")
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
            "Prospectively specified post hoc model-development experiment; "
            "retrospective within-panel binding-liability discovery only."
        ),
        "datasets": dataset_summary,
        "models": METHODS,
        "parameters": {
            "chemical_weight_beta": 8.0,
            "conformal_alpha": CONFORMAL_ALPHA,
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
    protocol = Path(__file__).resolve().parents[1] / "CONFORMAL_COVERAGE_PROTOCOL.md"
    if protocol.exists():
        shutil.copy2(protocol, output / protocol.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
