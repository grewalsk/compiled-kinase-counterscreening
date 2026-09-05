#!/usr/bin/env python3
"""Exploratory threshold/target-resolution audit after PKIS1 validation."""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from rdkit import rdBase

from run_compiled_coverage import Benchmark, holm_adjust, molecular_arrays, sha256, write_csv
from run_pkis1_external_validation import (
    MASTER_SEED,
    PKIS1_SHA256,
    _comparison,
    load_pkis1,
    run_condition,
)


PARENT_THRESHOLDS = [50, 70, 80, 90]
CONSTRUCT_THRESHOLDS = [80, 90]
EXCLUSIONS = ["standard", "scaffold", "nearest10"]


def continuous_matrices(
    raw_path: Path, base: Benchmark
) -> tuple[np.ndarray, list[str], list[str], np.ndarray, list[str], list[str]]:
    if sha256(raw_path) != PKIS1_SHA256:
        raise RuntimeError("PKIS1 SHA-256 mismatch")
    frame = pd.read_csv(raw_path, dtype={"CHEMBL_ID": str, "ASSAY_CHEMBL_ID": str})
    assay_text = frame["ASSAY"].fillna("").astype(str)
    selected = frame.loc[
        assay_text.str.contains(" 1 uM", regex=False)
        & assay_text.str.contains("[Nanosyn]", regex=False)
    ].copy()
    selected["VALUE"] = pd.to_numeric(selected["VALUE"], errors="raise")
    repeat = (
        selected.groupby(
            ["CHEMBL_ID", "ASSAY_CHEMBL_ID", "TARGET_CHEMBL_ID"], sort=True
        )["VALUE"]
        .mean()
        .reset_index()
    )
    parent_ids = sorted(repeat["TARGET_CHEMBL_ID"].astype(str).unique())
    assay_ids = sorted(repeat["ASSAY_CHEMBL_ID"].astype(str).unique())
    parent = (
        repeat.groupby(["CHEMBL_ID", "TARGET_CHEMBL_ID"], sort=True)["VALUE"]
        .max()
        .unstack("TARGET_CHEMBL_ID")
        .reindex(columns=parent_ids)
    )
    construct = (
        repeat.groupby(["CHEMBL_ID", "ASSAY_CHEMBL_ID"], sort=True)["VALUE"]
        .mean()
        .unstack("ASSAY_CHEMBL_ID")
        .reindex(columns=assay_ids)
    )
    source_groups = [source_ids.split(";") for source_ids in base.source_ids]
    parent_values = np.vstack(
        [parent.reindex(identifiers).median(axis=0, skipna=True).to_numpy(float) for identifiers in source_groups]
    )
    construct_values = np.vstack(
        [construct.reindex(identifiers).median(axis=0, skipna=True).to_numpy(float) for identifiers in source_groups]
    )
    if parent_values.shape != (356, 200) or construct_values.shape != (356, 224):
        raise RuntimeError(
            f"Unexpected audit matrices: {parent_values.shape}/{construct_values.shape}"
        )
    parent_name_map = (
        selected.groupby("TARGET_CHEMBL_ID")["TARGET"]
        .agg(lambda values: ";".join(sorted(set(map(str, values)))))
        .to_dict()
    )
    assay_name_map = (
        selected.groupby("ASSAY_CHEMBL_ID")["ASSAY"]
        .agg(lambda values: ";".join(sorted(set(map(str, values)))))
        .to_dict()
    )
    return (
        parent_values,
        parent_ids,
        [parent_name_map[target] for target in parent_ids],
        construct_values,
        assay_ids,
        [assay_name_map[target] for target in assay_ids],
    )


def threshold_benchmark(
    base: Benchmark,
    values: np.ndarray,
    target_ids: list[str],
    target_names: list[str],
    resolution: str,
    threshold: int,
) -> Benchmark:
    labels = np.where(np.isnan(values), np.nan, (values > threshold).astype(float))
    prefix = "P" if resolution == "parent" else "A"
    return Benchmark(
        name=f"PKIS1_{resolution.upper()}_GT{threshold}",
        compound_keys=base.compound_keys,
        canonical_smiles=base.canonical_smiles,
        target_keys=[f"{prefix}{index:04d}" for index in range(1, len(target_ids) + 1)],
        target_names=target_names,
        labels={"activity": labels},
        candidate_masks=~np.isnan(labels),
        chemotypes=base.chemotypes,
        source_ids=base.source_ids,
        metadata={**base.metadata, "target_resolution": resolution, "threshold_gt": threshold},
    )


def parse_condition(condition: str) -> tuple[str, int, str]:
    resolution, threshold_text, exclusion = condition.split("__")
    return resolution, int(threshold_text.removeprefix("gt")), exclusion


def make_figure(comparison_rows: list[dict], output: Path) -> None:
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    bounded = [row for row in comparison_rows if row["comparison"] == "bounded_minus_marginal"]
    figure, axis = plt.subplots(figsize=(8.2, 5.3), constrained_layout=True)
    colors = {"standard": "#136f63", "scaffold": "#bf7c2f", "nearest10": "#5e548e"}
    for exclusion in EXCLUSIONS:
        rows = [
            row
            for row in bounded
            if parse_condition(row["condition"])[0] == "parent"
            and parse_condition(row["condition"])[2] == exclusion
        ]
        rows.sort(key=lambda row: parse_condition(row["condition"])[1])
        axis.plot(
            [parse_condition(row["condition"])[1] for row in rows],
            [row["estimate"] for row in rows],
            marker="o",
            label=f"parent / {exclusion}",
            color=colors[exclusion],
        )
    for exclusion in EXCLUSIONS:
        rows = [
            row
            for row in bounded
            if parse_condition(row["condition"])[0] == "construct"
            and parse_condition(row["condition"])[2] == exclusion
        ]
        rows.sort(key=lambda row: parse_condition(row["condition"])[1])
        axis.scatter(
            [parse_condition(row["condition"])[1] for row in rows],
            [row["estimate"] for row in rows],
            marker="x",
            s=65,
            color=colors[exclusion],
            label=f"construct / {exclusion}",
        )
    axis.axhline(0.0, color="#222222", linewidth=1)
    axis.set_xlabel("strict percent-inhibition threshold")
    axis.set_ylabel("bounded minus marginal AUDC")
    axis.set_title("Exploratory PKIS1 threshold and target-resolution audit")
    axis.legend(frameon=False, ncol=2, fontsize=8)
    figure.savefig(figures / "threshold_resolution_audit.png", dpi=180)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkis1-raw", type=Path, required=True)
    parser.add_argument("--pkis2-xlsx", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    base = load_pkis1(args.pkis1_raw, args.pkis2_xlsx, args.data_dir)
    (
        parent_values,
        parent_ids,
        parent_names,
        construct_values,
        construct_ids,
        construct_names,
    ) = continuous_matrices(args.pkis1_raw, base)
    specifications = [
        ("parent", threshold, parent_values, parent_ids, parent_names)
        for threshold in PARENT_THRESHOLDS
    ] + [
        ("construct", threshold, construct_values, construct_ids, construct_names)
        for threshold in CONSTRUCT_THRESHOLDS
    ]
    all_case_rows: list[dict] = []
    for resolution, threshold, values, target_ids, target_names in specifications:
        benchmark = threshold_benchmark(
            base, values, target_ids, target_names, resolution, threshold
        )
        similarities, scaffolds = molecular_arrays(benchmark)
        for exclusion in EXCLUSIONS:
            condition = f"{resolution}__gt{threshold}__{exclusion}"
            rows, stability = run_condition(
                benchmark,
                similarities,
                scaffolds,
                condition,
                "activity",
                exclusion,
                args.workers,
                0,
            )
            if stability:
                raise AssertionError("Robustness audit unexpectedly ran stability")
            all_case_rows.extend(rows)

    by_key: dict[tuple[str, str], dict[str, dict]] = {}
    for row in all_case_rows:
        by_key.setdefault((row["condition"], row["method"]), {})[row["compound_key"]] = row
    comparison_rows: list[dict] = []
    metric_rows: list[dict] = []
    conditions = sorted({row["condition"] for row in all_case_rows})
    seed = 0
    for condition in conditions:
        for name, method in [
            ("bounded_minus_marginal", "BOUNDED_COVERAGE"),
            ("unbounded_minus_marginal", "UNBOUNDED_COVERAGE"),
        ]:
            summary, details = _comparison(
                by_key[(condition, method)],
                by_key[(condition, "MARGINAL")],
                condition,
                name,
                args.bootstrap,
                args.permutations,
                300 + seed,
            )
            comparison_rows.append(
                {
                    **summary,
                    "ci95_low": summary["ci95"][0],
                    "ci95_high": summary["ci95"][1],
                }
            )
            metric_rows.extend(details)
            seed += 1
    bounded_p = {
        row["condition"]: row["permutation_p_two_sided"]
        for row in comparison_rows
        if row["comparison"] == "bounded_minus_marginal"
    }
    adjusted = holm_adjust(bounded_p)
    for row in comparison_rows:
        row["p_holm_18_bounded_conditions"] = (
            adjusted[row["condition"]]
            if row["comparison"] == "bounded_minus_marginal"
            else ""
        )

    prevalence_rows = []
    for condition in conditions:
        rows = list(by_key[(condition, "MARGINAL")].values())
        positives = np.asarray([row["hidden_positive_count"] for row in rows], dtype=float)
        candidates = np.asarray([row["candidate_count"] for row in rows], dtype=float)
        resolution, threshold, exclusion = parse_condition(condition)
        prevalence_rows.append(
            {
                "condition": condition,
                "resolution": resolution,
                "threshold_gt": threshold,
                "exclusion": exclusion,
                "n": len(rows),
                "mean_positive_targets": float(positives.mean()),
                "median_positive_targets": float(np.median(positives)),
                "fraction_compounds_with_any_positive": float(np.mean(positives > 0)),
                "mean_positive_fraction_among_observed": float(np.mean(positives / candidates)),
            }
        )
    bounded_effects = [
        {
            "condition": row["condition"],
            "estimate": row["estimate"],
            "ci95": row["ci95"],
            "p_two_sided": row["permutation_p_two_sided"],
            "p_holm_18": row["p_holm_18_bounded_conditions"],
            "large_delay_ge_10": row["large_delay_ge_10"],
        }
        for row in comparison_rows
        if row["comparison"] == "bounded_minus_marginal"
    ]
    if any(row["large_delay_ge_10"] for row in bounded_effects):
        raise AssertionError("Bounded audit violated delay guarantee")
    summary = {
        "status": "EXPLORATORY_POSTHOC_ROBUSTNESS_AUDIT",
        "external_validation_verdict_unchanged": True,
        "dataset": {
            "compounds": 356,
            "parent_targets": 200,
            "assay_constructs": 224,
            "source_sha256": PKIS1_SHA256,
            "license": "ChEMBL CC BY-SA 3.0",
        },
        "grid": {
            "parent_thresholds": PARENT_THRESHOLDS,
            "construct_thresholds": CONSTRUCT_THRESHOLDS,
            "exclusions": EXCLUSIONS,
        },
        "bounded_effects": bounded_effects,
        "runtime_seconds": time.time() - started,
    }
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "case_metrics.csv.gz", all_case_rows)
    write_csv(output / "comparison_summary.csv", comparison_rows)
    write_csv(output / "metric_effects.csv", metric_rows)
    write_csv(output / "prevalence_summary.csv", prevalence_rows)
    make_figure(comparison_rows, output)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
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
        "master_seed": MASTER_SEED,
        "bootstrap": args.bootstrap,
        "permutations": args.permutations,
    }
    (output / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    protocol = Path(__file__).resolve().parents[1] / "PKIS1_ROBUSTNESS_AUDIT_PROTOCOL.md"
    shutil.copy2(protocol, output / protocol.name)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
