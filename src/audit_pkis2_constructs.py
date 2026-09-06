#!/usr/bin/env python3
"""Frozen post-result audit of PKIS2 at assay-construct resolution."""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from rdkit import rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from run_compiled_coverage import (
    Benchmark,
    PKIS2_METADATA,
    PKIS2_SHA256,
    holm_adjust,
    largest_fragment_smiles,
    molecular_arrays,
    sha256,
    write_csv,
)
from run_pkis1_external_validation import MASTER_SEED, _comparison, run_condition


CONDITIONS = [
    ("construct__gt90__standard", "gt90", "standard"),
    ("construct__gt80__standard", "gt80", "standard"),
    ("construct__gt90__scaffold", "gt90", "scaffold"),
    ("construct__gt90__chemotype", "gt90", "chemotype"),
]


def load_construct_benchmark(path: Path) -> Benchmark:
    observed_hash = sha256(path)
    if observed_hash != PKIS2_SHA256:
        raise RuntimeError(f"PKIS2 SHA-256 mismatch: {observed_hash}")
    frame = pd.read_excel(path, sheet_name="Table 4 - PKIS2 %Inh")
    missing_metadata = sorted(set(PKIS2_METADATA) - set(frame.columns))
    if missing_metadata:
        raise RuntimeError(f"PKIS2 metadata columns missing: {missing_metadata}")
    assays = sorted(str(column) for column in frame.columns if column not in PKIS2_METADATA)
    if len(assays) != 406:
        raise RuntimeError(f"Expected 406 PKIS2 assays, observed {len(assays)}")
    empty = frame[["Regno", "Compound", "Chemotype", "Smiles"]].isna().all(axis=1)
    excluded = frame.loc[empty, assays].stack()
    if (
        int(empty.sum()) != 1
        or len(excluded) != 1
        or str(excluded.index[0][1]) != "TEC"
        or not np.isclose(float(excluded.iloc[0]), 13.0)
    ):
        raise RuntimeError("Unexpected identifier-empty row or stray assay cells")
    frame = frame.loc[~empty].copy()
    numeric = frame[assays].apply(pd.to_numeric, errors="raise")
    if ((numeric < 0) | (numeric > 100)).any().any():
        raise RuntimeError("PKIS2 inhibition values outside [0,100]")

    chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    frame["parent_smiles"] = [
        largest_fragment_smiles(str(value), chooser) for value in frame["Smiles"]
    ]
    groups = [
        (str(smiles), group.index.tolist())
        for smiles, group in frame.groupby("parent_smiles", sort=True)
    ]
    values = np.vstack(
        [numeric.loc[indices].median(axis=0, skipna=True).to_numpy(float) for _, indices in groups]
    )
    if values.shape != (640, 406):
        raise RuntimeError(f"Unexpected construct matrix shape: {values.shape}")
    labels = {
        "gt90": np.where(np.isnan(values), np.nan, (values > 90.0).astype(float)),
        "gt80": np.where(np.isnan(values), np.nan, (values > 80.0).astype(float)),
    }
    source_ids = [
        ";".join(sorted(str(value) for value in frame.loc[indices, "Compound"]))
        for _, indices in groups
    ]
    chemotypes = [
        frozenset(str(value) for value in frame.loc[indices, "Chemotype"].dropna().unique())
        for _, indices in groups
    ]
    return Benchmark(
        name="PKIS2_CONSTRUCT_AUDIT",
        compound_keys=[f"X{index:04d}" for index in range(1, len(groups) + 1)],
        canonical_smiles=[smiles for smiles, _ in groups],
        target_keys=[f"A{index:04d}" for index in range(1, len(assays) + 1)],
        target_names=assays,
        labels=labels,
        candidate_masks=~np.isnan(labels["gt90"]),
        chemotypes=chemotypes,
        source_ids=source_ids,
        metadata={
            "source": "Drewry et al. 2017 PLOS ONE Supporting Table S4",
            "source_sha256": observed_hash,
            "license": "CC BY 4.0",
            "parent_structures": len(groups),
            "assay_constructs": len(assays),
            "target_resolution": "source assay construct",
        },
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkis2-xlsx", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


def output_manifest(output: Path) -> dict:
    return {
        "files": [
            {
                "path": path.relative_to(output).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in sorted(output.rglob("*"))
            if path.is_file() and path.name != "output_manifest.json"
        ]
    }


def main() -> None:
    args = parse_args()
    started = time.time()
    benchmark = load_construct_benchmark(args.pkis2_xlsx)
    similarities, scaffolds = molecular_arrays(benchmark)
    all_rows: list[dict] = []
    for condition, label, exclusion in CONDITIONS:
        rows, stability = run_condition(
            benchmark, similarities, scaffolds, condition, label, exclusion,
            args.workers, 0,
        )
        if stability:
            raise AssertionError("Construct audit unexpectedly ran stability")
        all_rows.extend(rows)

    indexed: dict[tuple[str, str], dict[str, dict]] = {}
    for row in all_rows:
        indexed.setdefault((row["condition"], row["method"]), {})[row["compound_key"]] = row
    comparisons: list[dict] = []
    metric_rows: list[dict] = []
    for offset, (condition, _, _) in enumerate(CONDITIONS):
        for comparison_index, (name, method) in enumerate([
            ("bounded_minus_marginal", "BOUNDED_COVERAGE"),
            ("unbounded_minus_marginal", "UNBOUNDED_COVERAGE"),
        ]):
            summary, details = _comparison(
                indexed[(condition, method)], indexed[(condition, "MARGINAL")],
                condition, name, args.bootstrap, args.permutations,
                700 + offset * 2 + comparison_index,
            )
            comparisons.append({
                **summary,
                "ci95_low": summary["ci95"][0],
                "ci95_high": summary["ci95"][1],
            })
            metric_rows.extend(details)
    bounded_p = {
        row["condition"]: row["permutation_p_two_sided"]
        for row in comparisons if row["comparison"] == "bounded_minus_marginal"
    }
    adjusted = holm_adjust(bounded_p)
    for row in comparisons:
        row["p_holm_4_bounded_conditions"] = (
            adjusted[row["condition"]]
            if row["comparison"] == "bounded_minus_marginal" else ""
        )
    prevalence = []
    for condition, _, exclusion in CONDITIONS:
        rows = list(indexed[(condition, "MARGINAL")].values())
        positives = np.asarray([row["hidden_positive_count"] for row in rows], float)
        candidates = np.asarray([row["candidate_count"] for row in rows], float)
        prevalence.append({
            "condition": condition,
            "exclusion": exclusion,
            "n": len(rows),
            "mean_positive_constructs": float(positives.mean()),
            "median_positive_constructs": float(np.median(positives)),
            "fraction_compounds_with_any_positive": float(np.mean(positives > 0)),
            "mean_positive_fraction_among_observed": float(np.mean(positives / candidates)),
        })
    if any(
        row["large_delay_ge_10"]
        for row in comparisons if row["comparison"] == "bounded_minus_marginal"
    ):
        raise AssertionError("Bounded method violated the rank-deadline guarantee")

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "case_metrics.csv.gz", all_rows)
    write_csv(output / "comparison_summary.csv", comparisons)
    write_csv(output / "metric_effects.csv", metric_rows)
    write_csv(output / "prevalence_summary.csv", prevalence)
    summary = {
        "status": "FROZEN_POSTRESULT_BIOLOGICAL_VALIDITY_AUDIT_COMPLETE",
        "cannot_change_external_validation_verdict": True,
        "dataset": benchmark.metadata,
        "conditions": [condition for condition, _, _ in CONDITIONS],
        "bounded_effects": [row for row in comparisons if row["comparison"] == "bounded_minus_marginal"],
        "runtime_seconds": time.time() - started,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "environment.json").write_text(
        json.dumps({
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "rdkit": rdBase.rdkitVersion,
            "cpu_only": True,
            "paid_api_cost_usd": 0,
            "workers": args.workers,
            "master_seed": MASTER_SEED,
            "bootstrap": args.bootstrap,
            "permutations": args.permutations,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    protocol = Path(__file__).resolve().parents[1] / "PKIS2_CONSTRUCT_AUDIT_PROTOCOL.md"
    shutil.copy2(protocol, output / protocol.name)
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
