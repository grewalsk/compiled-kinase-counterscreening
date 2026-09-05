#!/usr/bin/env python3
"""Run the frozen external PKIS1 validation of displacement-bounded coverage."""

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
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from run_compiled_coverage import (
    Benchmark,
    PKIS2_SHA256,
    coverage_order,
    holm_adjust,
    largest_fragment_smiles,
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
    weighted_probabilities,
    write_csv,
)
from run_conformal_coverage import top_k_jaccard
from run_risk_controlled_coverage import bounded_coverage_order


PKIS1_SHA256 = "81d7f9f82f7ee8e6b0f38dafe523da7254aaabe9449758a056372511d7868ad0"
PKIS1_SOURCE_COMMIT = "5fd3934f5789c371026fc9eece1846ff1294122b"
MASTER_SEED = 20260905
MAX_BUDGET = 20
MAX_DISPLACEMENT = 9
NONINFERIORITY_MARGIN = 0.015
STABILITY_REPLICATES = 50
STABILITY_FRACTION = 0.8
METHODS = ["MARGINAL", "UNBOUNDED_COVERAGE", "BOUNDED_COVERAGE"]
CONDITIONS = [
    ("pkis1_gt90_standard", "primary", "standard"),
    ("pkis1_gt80_sensitivity", "sensitivity", "standard"),
    ("pkis1_gt90_scaffold", "primary", "scaffold"),
    ("pkis1_gt90_leave10_nearest", "primary", "nearest10"),
]
CONDITION_IDS = {condition: index + 1 for index, (condition, _, _) in enumerate(CONDITIONS)}
COMPARISONS = [
    ("bounded_minus_marginal", "BOUNDED_COVERAGE", "MARGINAL"),
    ("unbounded_minus_marginal", "UNBOUNDED_COVERAGE", "MARGINAL"),
    ("bounded_minus_unbounded", "BOUNDED_COVERAGE", "UNBOUNDED_COVERAGE"),
]
METRICS = [
    "audc_1_20",
    "hit_at_1",
    "hit_at_3",
    "hit_at_5",
    "hit_at_10",
    "hit_at_20",
]

_WORKER_STATE: tuple | None = None


def _require_single_value(frame: pd.DataFrame, column: str, expected: str) -> None:
    observed = set(frame[column].dropna().astype(str))
    if observed != {expected}:
        raise RuntimeError(f"Unexpected PKIS1 {column}: {sorted(observed)}")


def load_pkis1(
    raw_path: Path,
    pkis2_path: Path,
    klaeger_data_dir: Path,
    enforce_expected_counts: bool = True,
) -> Benchmark:
    """Construct the frozen parent-structure by parent-target PKIS1 matrix."""
    observed_hash = sha256(raw_path)
    if observed_hash != PKIS1_SHA256:
        raise RuntimeError(
            f"Unexpected PKIS1 SHA-256: {observed_hash}; expected {PKIS1_SHA256}"
        )
    if sha256(pkis2_path) != PKIS2_SHA256:
        raise RuntimeError("PKIS2 SHA-256 mismatch")

    frame = pd.read_csv(raw_path, dtype={"CHEMBL_ID": str, "ASSAY_CHEMBL_ID": str})
    required = {
        "CHEMBL_ID",
        "SMILES",
        "ASSAY_CHEMBL_ID",
        "ASSAY",
        "ENDPOINT",
        "RELATION",
        "VALUE",
        "UNITS",
        "TARGET_CHEMBL_ID",
        "TARGET",
        "SPECIES",
        "SOURCE",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"PKIS1 columns missing: {missing}")
    assay = frame["ASSAY"].fillna("").astype(str)
    selected = frame.loc[
        assay.str.contains(" 1 uM", regex=False)
        & assay.str.contains("[Nanosyn]", regex=False)
    ].copy()
    _require_single_value(selected, "ENDPOINT", "Inhibition")
    _require_single_value(selected, "RELATION", "=")
    _require_single_value(selected, "SPECIES", "Homo sapiens")
    _require_single_value(selected, "SOURCE", "GSK Published Kinase Inhibitor Set")
    units = set(selected["UNITS"].dropna().astype(str))
    if units != {"%"}:
        raise RuntimeError(f"Unexpected PKIS1 units: {sorted(units)}")
    selected["VALUE"] = pd.to_numeric(selected["VALUE"], errors="raise")
    observed_values = selected["VALUE"].dropna()
    if not np.isfinite(observed_values).all():
        raise RuntimeError("PKIS1 contains nonfinite observed inhibition values")
    if selected[["CHEMBL_ID", "SMILES", "ASSAY_CHEMBL_ID", "TARGET_CHEMBL_ID"]].isna().any().any():
        raise RuntimeError("PKIS1 contains missing identifiers or structures")

    structure_rows = selected[["CHEMBL_ID", "SMILES"]].drop_duplicates()
    if structure_rows.groupby("CHEMBL_ID")["SMILES"].nunique().max() != 1:
        raise RuntimeError("A PKIS1 ChEMBL compound ID maps to multiple SMILES")
    chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    structure_rows["parent_smiles"] = [
        largest_fragment_smiles(str(value), chooser) for value in structure_rows["SMILES"]
    ]
    smiles_by_id = dict(zip(structure_rows["CHEMBL_ID"], structure_rows["parent_smiles"]))

    assay_parent_counts = selected.groupby("ASSAY_CHEMBL_ID")["TARGET_CHEMBL_ID"].nunique()
    if int(assay_parent_counts.max()) != 1:
        raise RuntimeError("A PKIS1 assay maps to multiple parent targets")
    target_names_by_id = (
        selected.groupby("TARGET_CHEMBL_ID")["TARGET"]
        .agg(lambda values: ";".join(sorted(set(map(str, values)))))
        .to_dict()
    )
    repeat_averaged = (
        selected.groupby(
            ["CHEMBL_ID", "ASSAY_CHEMBL_ID", "TARGET_CHEMBL_ID"], sort=True
        )["VALUE"]
        .mean()
        .reset_index()
    )
    parent_collapsed = (
        repeat_averaged.groupby(["CHEMBL_ID", "TARGET_CHEMBL_ID"], sort=True)["VALUE"]
        .max()
        .unstack("TARGET_CHEMBL_ID")
    )
    target_ids = sorted(parent_collapsed.columns.astype(str).tolist())
    parent_collapsed = parent_collapsed[target_ids]

    development_structures = set(load_pkis2(pkis2_path).canonical_smiles)
    development_structures.update(load_klaeger(klaeger_data_dir).canonical_smiles)
    structure_groups: list[tuple[str, list[str]]] = []
    for parent_smiles, rows in structure_rows.groupby("parent_smiles", sort=True):
        identifiers = sorted(set(rows["CHEMBL_ID"].astype(str)))
        structure_groups.append((str(parent_smiles), identifiers))
    development_overlap = [
        (smiles, identifiers)
        for smiles, identifiers in structure_groups
        if smiles in development_structures
    ]
    retained_groups = [
        (smiles, identifiers)
        for smiles, identifiers in structure_groups
        if smiles not in development_structures
    ]

    if enforce_expected_counts:
        counts = {
            "filtered_rows": len(selected),
            "compound_ids": selected["CHEMBL_ID"].nunique(),
            "raw_smiles": selected["SMILES"].nunique(),
            "parent_structures": len(structure_groups),
            "assays": selected["ASSAY_CHEMBL_ID"].nunique(),
            "parent_targets": selected["TARGET_CHEMBL_ID"].nunique(),
            "overlap_structures": len(development_overlap),
            "retained_structures": len(retained_groups),
        }
        expected = {
            "filtered_rows": 82656,
            "compound_ids": 366,
            "raw_smiles": 366,
            "parent_structures": 364,
            "assays": 224,
            "parent_targets": 200,
            "overlap_structures": 8,
            "retained_structures": 356,
        }
        if counts != expected:
            raise RuntimeError(f"Unexpected PKIS1 construction counts: {counts}; expected {expected}")

    values = np.vstack(
        [
            parent_collapsed.reindex(identifiers).median(axis=0, skipna=True).to_numpy(float)
            for _, identifiers in retained_groups
        ]
    )
    primary = np.where(np.isnan(values), np.nan, (values > 90.0).astype(float))
    sensitivity = np.where(np.isnan(values), np.nan, (values > 80.0).astype(float))
    compound_keys = [f"V{index:04d}" for index in range(1, len(retained_groups) + 1)]
    target_keys = [f"P{index:04d}" for index in range(1, len(target_ids) + 1)]
    return Benchmark(
        name="PKIS1_EXTERNAL",
        compound_keys=compound_keys,
        canonical_smiles=[smiles for smiles, _ in retained_groups],
        target_keys=target_keys,
        target_names=[target_names_by_id[target] for target in target_ids],
        labels={"primary": primary, "sensitivity": sensitivity},
        candidate_masks=~np.isnan(primary),
        chemotypes=[frozenset() for _ in retained_groups],
        source_ids=[";".join(identifiers) for _, identifiers in retained_groups],
        metadata={
            "source": "Zhang et al. 2019 archived raw ChEMBL PKIS1 Nanosyn panel",
            "source_commit": PKIS1_SOURCE_COMMIT,
            "source_sha256": observed_hash,
            "license": "ChEMBL CC BY-SA 3.0",
            "raw_rows": len(frame),
            "filtered_1uM_nanosyn_rows": len(selected),
            "source_compound_ids": selected["CHEMBL_ID"].nunique(),
            "source_assays": selected["ASSAY_CHEMBL_ID"].nunique(),
            "parent_targets": len(target_ids),
            "standardized_parent_structures_before_overlap_exclusion": len(structure_groups),
            "excluded_development_overlap_structures": len(development_overlap),
            "retained_parent_structures": len(retained_groups),
            "blank_units_retained": int(selected["UNITS"].isna().sum()),
            "raw_missing_values": int(selected["VALUE"].isna().sum()),
            "raw_values_below_zero": int((observed_values < 0.0).sum()),
            "raw_values_above_100": int((observed_values > 100.0).sum()),
            "collapsed_missing_values": int(np.isnan(values).sum()),
            "task_language": "first strong biochemical kinase profile-activity discovery",
        },
    )


def policy_orders(
    benchmark: Benchmark,
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
    unbounded = coverage_order(
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
        "MARGINAL": marginal,
        "UNBOUNDED_COVERAGE": unbounded,
        "BOUNDED_COVERAGE": bounded,
    }


def stability_diagnostic(
    benchmark: Benchmark,
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
    unbounded: list[float] = []
    bounded: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [MASTER_SEED, 311, CONDITION_IDS[condition], heldout, replicate]
            )
        )
        selected = np.sort(rng.choice(refs, size=sample_size, replace=False)).astype(int)
        orders = policy_orders(
            benchmark, labels, similarities, heldout, selected, candidates
        )
        unbounded.append(
            top_k_jaccard(
                full_orders["UNBOUNDED_COVERAGE"], orders["UNBOUNDED_COVERAGE"], 10
            )
        )
        bounded.append(
            top_k_jaccard(
                full_orders["BOUNDED_COVERAGE"], orders["BOUNDED_COVERAGE"], 10
            )
        )
    return {
        "dataset": benchmark.name,
        "condition": condition,
        "compound_key": benchmark.compound_keys[heldout],
        "replicates": replicates,
        "reference_fraction": STABILITY_FRACTION,
        "unbounded_coverage_top10_jaccard": float(np.mean(unbounded)),
        "bounded_coverage_top10_jaccard": float(np.mean(bounded)),
        "bounded_minus_unbounded": float(np.mean(bounded) - np.mean(unbounded)),
    }


def evaluate_outer_case(
    benchmark: Benchmark,
    similarities: np.ndarray,
    scaffolds: list[str],
    condition: str,
    label_name: str,
    exclusion: str,
    heldout: int,
    stability_replicates: int,
) -> tuple[list[dict], dict | None] | None:
    labels = benchmark.labels[label_name]
    refs = reference_indices(
        heldout, similarities, scaffolds, benchmark.chemotypes, exclusion
    )
    candidates = np.flatnonzero(
        benchmark.candidate_masks[heldout] & ~np.isnan(labels[heldout])
    )
    if len(refs) < 5 or len(candidates) < MAX_BUDGET:
        return None
    orders = policy_orders(benchmark, labels, similarities, heldout, refs, candidates)
    truth = labels[heldout]
    metrics_by_method = {
        method: metrics_from_sequence(orders[method], truth, MAX_BUDGET)
        for method in METHODS
    }
    marginal_first = metrics_by_method["MARGINAL"]["first_hit"]
    observed_delay = metrics_by_method["BOUNDED_COVERAGE"]["first_hit"] - marginal_first
    if observed_delay > MAX_DISPLACEMENT:
        raise AssertionError(f"Bounded coverage violated first-hit delay bound: {observed_delay}")

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
                "order": ";".join(benchmark.target_keys[index] for index in orders[method]),
            }
        )
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
    return rows, stability


def initialize_worker(benchmark: Benchmark, similarities: np.ndarray, scaffolds: list[str]) -> None:
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
    benchmark: Benchmark,
    similarities: np.ndarray,
    scaffolds: list[str],
    condition: str,
    label_name: str,
    exclusion: str,
    workers: int,
    stability_replicates: int,
) -> tuple[list[dict], list[dict]]:
    print(f"Running {condition}: {len(benchmark.compound_keys)} cases, {workers} workers", flush=True)
    start_method = "fork" if "fork" in multiprocessing.get_all_start_methods() else "spawn"
    context = multiprocessing.get_context(start_method)
    tasks = [
        (condition, label_name, exclusion, heldout, stability_replicates)
        for heldout in range(len(benchmark.compound_keys))
    ]
    by_heldout: dict[int, object] = {}
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
    stability_rows: list[dict] = []
    for heldout in sorted(by_heldout):
        result = by_heldout[heldout]
        if result is None:
            continue
        rows, stability = result
        case_rows.extend(rows)
        if stability is not None:
            stability_rows.append(stability)
    return case_rows, stability_rows


def _indexed_rows(case_rows: list[dict]) -> dict[tuple[str, str], dict[str, dict]]:
    indexed: dict[tuple[str, str], dict[str, dict]] = defaultdict(dict)
    for row in case_rows:
        indexed[(row["condition"], row["method"])][row["compound_key"]] = row
    return indexed


def _comparison(
    first: dict[str, dict],
    second: dict[str, dict],
    condition: str,
    comparison_name: str,
    bootstrap: int,
    permutations: int,
    seed_offset: int,
) -> tuple[dict, list[dict]]:
    keys = sorted(set(first) & set(second))
    audc = np.asarray(
        [first[key]["audc_1_20"] - second[key]["audc_1_20"] for key in keys],
        dtype=float,
    )
    assay = np.asarray(
        [second[key]["cost_to_first_censored"] - first[key]["cost_to_first_censored"] for key in keys],
        dtype=int,
    )
    pvalue = sign_flip_test(audc, permutations, MASTER_SEED + 6000 + seed_offset)
    summary = {
        "condition": condition,
        "comparison": comparison_name,
        "n": len(keys),
        "estimate": float(audc.mean()),
        "ci95": paired_bootstrap(audc, bootstrap, MASTER_SEED + 5000 + seed_offset)["ci95"],
        "permutation_p_two_sided": pvalue,
        "mean_assays_earlier": float(assay.mean()),
        "wins": int(np.sum(assay > 0)),
        "losses": int(np.sum(assay < 0)),
        "ties": int(np.sum(assay == 0)),
        "large_advance_ge_10": int(np.sum(assay >= 10)),
        "large_delay_ge_10": int(np.sum(assay <= -10)),
    }
    detail: list[dict] = []
    for metric_index, metric in enumerate(METRICS):
        differences = np.asarray(
            [first[key][metric] - second[key][metric] for key in keys], dtype=float
        )
        interval = paired_bootstrap(
            differences,
            bootstrap,
            MASTER_SEED + 12000 + seed_offset * 20 + metric_index,
        )["ci95"]
        detail.append(
            {
                "condition": condition,
                "comparison": comparison_name,
                "metric": metric,
                "n": len(keys),
                "estimate": float(differences.mean()),
                "ci95_low": interval[0],
                "ci95_high": interval[1],
            }
        )
    return summary, detail


def summarize(
    case_rows: list[dict],
    stability_rows: list[dict],
    bootstrap: int,
    permutations: int,
) -> tuple[dict, list[dict], list[dict]]:
    indexed = _indexed_rows(case_rows)
    comparisons: dict[str, dict] = {}
    comparison_rows: list[dict] = []
    metric_rows: list[dict] = []
    seed_offset = 0
    for condition, _, _ in CONDITIONS:
        comparisons[condition] = {}
        for name, method, baseline in COMPARISONS:
            result, detail = _comparison(
                indexed[(condition, method)],
                indexed[(condition, baseline)],
                condition,
                name,
                bootstrap,
                permutations,
                seed_offset,
            )
            comparisons[condition][name] = result
            comparison_rows.append(
                {
                    **result,
                    "ci95_low": result["ci95"][0],
                    "ci95_high": result["ci95"][1],
                }
            )
            metric_rows.extend(detail)
            seed_offset += 1

    bounded_pvalues = {
        condition: comparisons[condition]["bounded_minus_marginal"][
            "permutation_p_two_sided"
        ]
        for condition, _, _ in CONDITIONS
    }
    adjusted = holm_adjust(bounded_pvalues)
    for condition in bounded_pvalues:
        comparisons[condition]["bounded_minus_marginal"][
            "permutation_p_holm_four_conditions"
        ] = adjusted[condition]
    for row in comparison_rows:
        if row["comparison"] == "bounded_minus_marginal":
            row["permutation_p_holm_four_conditions"] = adjusted[row["condition"]]
        else:
            row["permutation_p_holm_four_conditions"] = ""

    scaffold = indexed[("pkis1_gt90_scaffold", "BOUNDED_COVERAGE")]
    scaffold_base = indexed[("pkis1_gt90_scaffold", "MARGINAL")]
    nearest = indexed[("pkis1_gt90_leave10_nearest", "BOUNDED_COVERAGE")]
    nearest_base = indexed[("pkis1_gt90_leave10_nearest", "MARGINAL")]
    macro_keys = sorted(set(scaffold) & set(scaffold_base) & set(nearest) & set(nearest_base))
    macro = np.asarray(
        [
            0.5
            * (
                scaffold[key]["audc_1_20"]
                - scaffold_base[key]["audc_1_20"]
                + nearest[key]["audc_1_20"]
                - nearest_base[key]["audc_1_20"]
            )
            for key in macro_keys
        ],
        dtype=float,
    )
    macro_result = {
        "n": len(macro_keys),
        "estimate": float(macro.mean()),
        "ci95": paired_bootstrap(macro, bootstrap, MASTER_SEED + 19001)["ci95"],
        "permutation_p_two_sided": sign_flip_test(
            macro, permutations, MASTER_SEED + 19002
        ),
    }

    stability = None
    if stability_rows:
        unbounded = np.asarray(
            [row["unbounded_coverage_top10_jaccard"] for row in stability_rows]
        )
        bounded = np.asarray(
            [row["bounded_coverage_top10_jaccard"] for row in stability_rows]
        )
        difference = bounded - unbounded
        stability = {
            "condition": "pkis1_gt90_scaffold",
            "n": len(stability_rows),
            "subsamples_per_case": stability_rows[0]["replicates"],
            "mean_unbounded_top10_jaccard": float(unbounded.mean()),
            "mean_bounded_top10_jaccard": float(bounded.mean()),
            "mean_difference_bounded_minus_unbounded": float(difference.mean()),
            "difference_ci95": paired_bootstrap(
                difference, bootstrap, MASTER_SEED + 19501
            )["ci95"],
        }

    bm = {
        condition: comparisons[condition]["bounded_minus_marginal"]
        for condition, _, _ in CONDITIONS
    }
    no_large_delays = all(result["large_delay_ge_10"] == 0 for result in bm.values())
    if not no_large_delays:
        raise AssertionError("Observed PKIS1 output violates the displacement guarantee")
    success = {
        "gt90_standard_nonnegative_and_noninferior": bool(
            bm["pkis1_gt90_standard"]["estimate"] >= 0.0
            and bm["pkis1_gt90_standard"]["ci95"][0] >= -NONINFERIORITY_MARGIN
        ),
        "gt80_sensitivity_nonnegative_and_noninferior": bool(
            bm["pkis1_gt80_sensitivity"]["estimate"] >= 0.0
            and bm["pkis1_gt80_sensitivity"]["ci95"][0] >= -NONINFERIORITY_MARGIN
        ),
        "each_shift_nonnegative_and_noninferior": bool(
            all(
                bm[condition]["estimate"] >= 0.0
                and bm[condition]["ci95"][0] >= -NONINFERIORITY_MARGIN
                for condition in [
                    "pkis1_gt90_scaffold",
                    "pkis1_gt90_leave10_nearest",
                ]
            )
        ),
        "positive_shift_macro_ci_excludes_zero": bool(
            macro_result["estimate"] > 0.0 and macro_result["ci95"][0] > 0.0
        ),
        "zero_large_delays_all_conditions": bool(no_large_delays),
    }
    success["composite_success"] = all(success.values())
    return (
        {
            "comparisons": comparisons,
            "chemical_shift_macro_effect": macro_result,
            "stability": stability,
            "frozen_external_success_criteria": success,
            "displacement_guarantee_verified": True,
        },
        comparison_rows,
        metric_rows,
    )


def make_figures(summary: dict, output: Path) -> None:
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    conditions = [condition for condition, _, _ in CONDITIONS]
    labels = [">90 standard", ">80 standard", ">90 scaffold", ">90 remove-10"]
    series = [
        ("bounded_minus_marginal", "bounded", "#136f63", -0.12),
        ("unbounded_minus_marginal", "unbounded", "#8c8c8c", 0.12),
    ]
    x = np.arange(len(conditions))
    figure, axis = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    for key, name, color, offset in series:
        values = [summary["comparisons"][condition][key] for condition in conditions]
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
    axis.axhline(-NONINFERIORITY_MARGIN, color="#9a4f2d", linewidth=1, linestyle="--")
    axis.set_xticks(x, labels, rotation=22, ha="right")
    axis.set_ylabel("AUDC difference versus marginal")
    axis.set_title("Untuned external validation on PKIS1")
    axis.legend(frameon=False)
    figure.savefig(figures / "pkis1_external_effects.png", dpi=180)
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
    parser.add_argument("--pkis1-raw", type=Path, required=True)
    parser.add_argument("--pkis2-xlsx", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--validate-input-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    if args.workers < 1:
        raise ValueError("Worker count must be positive")
    benchmark = load_pkis1(args.pkis1_raw, args.pkis2_xlsx, args.data_dir)
    if args.validate_input_only:
        print(json.dumps({"metadata": benchmark.metadata}, indent=2, sort_keys=True))
        return

    bootstrap = min(args.bootstrap, 250) if args.smoke_test else args.bootstrap
    permutations = min(args.permutations, 1_000) if args.smoke_test else args.permutations
    stability_replicates = 3 if args.smoke_test else STABILITY_REPLICATES
    if args.smoke_test:
        keep = np.arange(min(40, len(benchmark.compound_keys)))
        benchmark = Benchmark(
            name=benchmark.name,
            compound_keys=[benchmark.compound_keys[index] for index in keep],
            canonical_smiles=[benchmark.canonical_smiles[index] for index in keep],
            target_keys=benchmark.target_keys,
            target_names=benchmark.target_names,
            labels={name: values[keep] for name, values in benchmark.labels.items()},
            candidate_masks=benchmark.candidate_masks[keep],
            chemotypes=[benchmark.chemotypes[index] for index in keep],
            source_ids=[benchmark.source_ids[index] for index in keep],
            metadata={**benchmark.metadata, "smoke_test_limit": len(keep)},
        )

    similarities, scaffolds = molecular_arrays(benchmark)
    all_case_rows: list[dict] = []
    all_stability_rows: list[dict] = []
    for condition, label_name, exclusion in CONDITIONS:
        condition_stability = (
            stability_replicates if condition == "pkis1_gt90_scaffold" else 0
        )
        case_rows, stability_rows = run_condition(
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
        all_stability_rows.extend(stability_rows)

    statistical, comparison_rows, metric_rows = summarize(
        all_case_rows, all_stability_rows, bootstrap, permutations
    )
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "case_metrics.csv.gz", all_case_rows)
    write_csv(output / "stability_diagnostics.csv.gz", all_stability_rows)
    write_csv(output / "comparison_summary.csv", comparison_rows)
    write_csv(output / "metric_effects.csv", metric_rows)
    make_figures(statistical, output)
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
        "status": "SMOKE_TEST_ONLY" if args.smoke_test else "FROZEN_EXTERNAL_VALIDATION_COMPLETE",
        "interpretation_scope": (
            "Retrospective external validation on a separate PKIS1 compound panel; "
            "biochemical within-panel kinase profile-activity discovery only."
        ),
        "dataset": {
            **benchmark.metadata,
            "compounds": len(benchmark.compound_keys),
            "targets": len(benchmark.target_keys),
        },
        "models": METHODS,
        "parameters": {
            "chemical_weight_beta": 8.0,
            "maximum_marginal_rank_displacement": MAX_DISPLACEMENT,
            "max_budget": MAX_BUDGET,
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
    protocol = Path(__file__).resolve().parents[1] / "PKIS1_EXTERNAL_VALIDATION_PROTOCOL.md"
    if protocol.exists():
        shutil.copy2(protocol, output / protocol.name)
    amendment = Path(__file__).resolve().parents[1] / "PKIS1_EXTERNAL_VALIDATION_AMENDMENT.md"
    if amendment.exists():
        shutil.copy2(amendment, output / amendment.name)
    amendment_2 = Path(__file__).resolve().parents[1] / "PKIS1_EXTERNAL_VALIDATION_AMENDMENT_2.md"
    if amendment_2.exists():
        shutil.copy2(amendment_2, output / amendment_2.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
