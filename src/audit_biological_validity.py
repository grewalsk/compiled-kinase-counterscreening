#!/usr/bin/env python3
"""Reproducible source-semantics and case-level kinase-panel audit."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from audit_pkis1_thresholds import continuous_matrices
from run_compiled_coverage import (
    PKIS2_METADATA,
    PKIS2_SHA256,
    pkis2_parent_target,
    sha256,
    write_csv,
)
from run_pkis1_external_validation import PKIS1_SHA256, load_pkis1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkis1-raw", type=Path, required=True)
    parser.add_argument("--pkis2-xlsx", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--external-cases", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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


def selected_pkis1(path: Path) -> pd.DataFrame:
    if sha256(path) != PKIS1_SHA256:
        raise RuntimeError("PKIS1 SHA-256 mismatch")
    frame = pd.read_csv(path, dtype={
        "CHEMBL_ID": str, "ASSAY_CHEMBL_ID": str, "TARGET_CHEMBL_ID": str,
    })
    assay = frame["ASSAY"].fillna("").astype(str)
    selected = frame.loc[
        assay.str.contains(" 1 uM", regex=False)
        & assay.str.contains("[Nanosyn]", regex=False)
    ].copy()
    selected["VALUE"] = pd.to_numeric(selected["VALUE"], errors="raise")
    return selected


def construct_groups(pkis1: pd.DataFrame, pkis2_path: Path) -> list[dict]:
    rows: list[dict] = []
    for parent, group in pkis1.groupby("TARGET_CHEMBL_ID", sort=True):
        assays = group[["ASSAY_CHEMBL_ID", "ASSAY"]].drop_duplicates().sort_values("ASSAY_CHEMBL_ID")
        if len(assays) > 1:
            rows.append({
                "panel": "PKIS1",
                "parent_target": str(parent),
                "parent_display_name": ";".join(sorted(set(map(str, group["TARGET"])))),
                "construct_count": len(assays),
                "construct_ids": ";".join(assays["ASSAY_CHEMBL_ID"].astype(str)),
                "construct_names": " || ".join(assays["ASSAY"].astype(str)),
                "collapse_semantics": "maximum inhibition: any assayed construct of parent",
            })
    if sha256(pkis2_path) != PKIS2_SHA256:
        raise RuntimeError("PKIS2 SHA-256 mismatch")
    frame = pd.read_excel(pkis2_path, sheet_name="Table 4 - PKIS2 %Inh")
    assays = [str(column) for column in frame.columns if column not in PKIS2_METADATA]
    groups: dict[str, list[str]] = defaultdict(list)
    for assay in assays:
        groups[pkis2_parent_target(assay)].append(assay)
    for parent in sorted(groups):
        if len(groups[parent]) > 1:
            rows.append({
                "panel": "PKIS2",
                "parent_target": parent,
                "parent_display_name": parent,
                "construct_count": len(groups[parent]),
                "construct_ids": ";".join(sorted(groups[parent])),
                "construct_names": " || ".join(sorted(groups[parent])),
                "collapse_semantics": "maximum inhibition: any profiled state/domain/construct of parent",
            })
    return rows


def discordance_rows(
    base, parent_values: np.ndarray, parent_ids: list[str],
    construct_values: np.ndarray, construct_ids: list[str], pkis1: pd.DataFrame,
) -> list[dict]:
    assay_parent = (
        pkis1[["ASSAY_CHEMBL_ID", "TARGET_CHEMBL_ID"]]
        .drop_duplicates().set_index("ASSAY_CHEMBL_ID")["TARGET_CHEMBL_ID"].to_dict()
    )
    multi = {
        parent for parent, count in
        pd.Series(list(assay_parent.values())).value_counts().items() if count > 1
    }
    construct_indices = defaultdict(list)
    for index, assay in enumerate(construct_ids):
        construct_indices[assay_parent[assay]].append(index)
    rows = []
    for threshold in [80, 90]:
        measured = discordant = any_positive = all_positive = 0
        for compound_index in range(len(base.compound_keys)):
            for parent in sorted(multi):
                indices = construct_indices[parent]
                observed = construct_values[compound_index, indices]
                observed = observed[~np.isnan(observed)]
                if len(observed) < 2:
                    continue
                measured += 1
                labels = observed > threshold
                any_positive += int(labels.any())
                all_positive += int(labels.all())
                discordant += int(labels.any() and not labels.all())
        rows.append({
            "panel": "PKIS1 retained external set",
            "threshold": f">{threshold}",
            "multi_construct_parents": len(multi),
            "compound_parent_cells_with_at_least_two_constructs": measured,
            "construct_discordant_cells": discordant,
            "construct_discordant_fraction": discordant / measured,
            "any_construct_positive_cells": any_positive,
            "all_observed_constructs_positive_cells": all_positive,
        })
    return rows


def boundary_rows(
    parent_values: np.ndarray, construct_values: np.ndarray, pkis1: pd.DataFrame,
    pkis2_path: Path,
) -> list[dict]:
    workbook = pd.read_excel(pkis2_path, sheet_name="Table 4 - PKIS2 %Inh")
    pkis2_assays = [str(column) for column in workbook.columns if column not in PKIS2_METADATA]
    pkis2_numeric = workbook[pkis2_assays].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    datasets = [
        ("PKIS1 selected raw rows", pkis1["VALUE"].to_numpy(float)),
        ("PKIS1 retained parent matrix", parent_values.ravel()),
        ("PKIS1 retained construct matrix", construct_values.ravel()),
        ("PKIS2 raw workbook cells", pkis2_numeric.ravel()),
    ]
    rows = []
    for name, values in datasets:
        values = values[~np.isnan(values)]
        for threshold in [80, 90]:
            rows.append({
                "dataset": name,
                "boundary": threshold,
                "observed_cells": len(values),
                "exact_boundary_cells": int(np.sum(values == threshold)),
                "fraction_exact_boundary": float(np.mean(values == threshold)),
            })
    return rows


def first_hit_target(row: pd.Series) -> str:
    cost = int(row["cost_to_first_censored"])
    order = str(row["order"]).split(";")
    return order[cost - 1] if 1 <= cost <= len(order) else "CENSORED"


def case_review_rows(base, parent_ids: list[str], cases: pd.DataFrame) -> list[dict]:
    target_name = {
        key: str(name).strip() for key, name in zip(base.target_keys, base.target_names)
    }
    target_id = dict(zip(base.target_keys, parent_ids))
    output = []
    for condition in sorted(cases["condition"].unique()):
        part = cases.loc[cases["condition"] == condition]
        pivot = part.pivot(index="compound_key", columns="method")
        for compound in sorted(pivot.index):
            marginal = part[(part["compound_key"] == compound) & (part["method"] == "MARGINAL")].iloc[0]
            bounded = part[(part["compound_key"] == compound) & (part["method"] == "BOUNDED_COVERAGE")].iloc[0]
            unbounded = part[(part["compound_key"] == compound) & (part["method"] == "UNBOUNDED_COVERAGE")].iloc[0]
            bounded_delay = int(bounded["cost_to_first_censored"] - marginal["cost_to_first_censored"])
            unbounded_delay = int(unbounded["cost_to_first_censored"] - marginal["cost_to_first_censored"])
            bounded_audc = float(bounded["audc_1_20"] - marginal["audc_1_20"])
            reasons = []
            if unbounded_delay >= 10:
                reasons.append("unbounded_delay_ge_10")
            if int(marginal["hidden_positive_count"]) >= 20:
                reasons.append("broad_profile_ge_20_positives")
            if bounded_audc >= 0.30:
                reasons.append("large_bounded_gain")
            if bounded_audc <= -0.30:
                reasons.append("large_bounded_loss")
            if not reasons:
                continue
            hit_keys = {
                "marginal": first_hit_target(marginal),
                "bounded": first_hit_target(bounded),
                "unbounded": first_hit_target(unbounded),
            }
            output.append({
                "condition": condition,
                "compound_key": compound,
                "source_ids": marginal["source_ids"],
                "review_reasons": ";".join(reasons),
                "candidate_count": int(marginal["candidate_count"]),
                "hidden_positive_count": int(marginal["hidden_positive_count"]),
                "marginal_first_hit": int(marginal["cost_to_first_censored"]),
                "bounded_first_hit": int(bounded["cost_to_first_censored"]),
                "unbounded_first_hit": int(unbounded["cost_to_first_censored"]),
                "bounded_delay_vs_marginal": bounded_delay,
                "unbounded_delay_vs_marginal": unbounded_delay,
                "bounded_audc_difference": bounded_audc,
                "unbounded_audc_difference": float(unbounded["audc_1_20"] - marginal["audc_1_20"]),
                "marginal_first_target_key": hit_keys["marginal"],
                "marginal_first_target_id": target_id.get(hit_keys["marginal"], "CENSORED"),
                "marginal_first_target_name": target_name.get(hit_keys["marginal"], "CENSORED"),
                "bounded_first_target_key": hit_keys["bounded"],
                "bounded_first_target_id": target_id.get(hit_keys["bounded"], "CENSORED"),
                "bounded_first_target_name": target_name.get(hit_keys["bounded"], "CENSORED"),
                "unbounded_first_target_key": hit_keys["unbounded"],
                "unbounded_first_target_id": target_id.get(hit_keys["unbounded"], "CENSORED"),
                "unbounded_first_target_name": target_name.get(hit_keys["unbounded"], "CENSORED"),
            })
    return output


def flagged_source_records(pkis1: pd.DataFrame, review: list[dict]) -> list[dict]:
    """Expose every raw first-hit construct record for every flagged case."""
    requested: set[tuple[str, str, str, str]] = set()
    for row in review:
        for method in ["marginal", "bounded", "unbounded"]:
            target = row[f"{method}_first_target_id"]
            if target == "CENSORED":
                continue
            for compound_id in str(row["source_ids"]).split(";"):
                requested.add((row["condition"], compound_id, method, target))
    output = []
    for condition, compound_id, method, target in sorted(requested):
        records = pkis1.loc[
            (pkis1["CHEMBL_ID"] == compound_id)
            & (pkis1["TARGET_CHEMBL_ID"] == target)
        ].sort_values("ASSAY_CHEMBL_ID")
        if records.empty:
            raise RuntimeError(f"Missing requested source record: {compound_id}/{target}")
        for _, record in records.iterrows():
            output.append({
                "condition": condition,
                "compound_id": compound_id,
                "method_whose_first_hit": method,
                "target_chembl_id": target,
                "target_name": record["TARGET"],
                "assay_chembl_id": record["ASSAY_CHEMBL_ID"],
                "assay_description": record["ASSAY"],
                "endpoint": record["ENDPOINT"],
                "relation": record["RELATION"],
                "value": record["VALUE"],
                "units": record["UNITS"],
                "species": record["SPECIES"],
                "source": record["SOURCE"],
            })
    return output


def main() -> None:
    args = parse_args()
    base = load_pkis1(args.pkis1_raw, args.pkis2_xlsx, args.data_dir)
    pkis1 = selected_pkis1(args.pkis1_raw)
    (
        parent_values, parent_ids, _, construct_values, construct_ids, _,
    ) = continuous_matrices(args.pkis1_raw, base)
    cases = pd.read_csv(args.external_cases)
    if set(cases["compound_key"]) != set(base.compound_keys):
        raise RuntimeError("External case file does not match retained PKIS1 compounds")

    groups = construct_groups(pkis1, args.pkis2_xlsx)
    discordance = discordance_rows(
        base, parent_values, parent_ids, construct_values, construct_ids, pkis1
    )
    boundaries = boundary_rows(parent_values, construct_values, pkis1, args.pkis2_xlsx)
    review = case_review_rows(base, parent_ids, cases)
    source_records = flagged_source_records(pkis1, review)
    candidates = (
        cases.loc[cases["method"] == "MARGINAL"]
        .groupby("condition")["candidate_count"]
        .agg(["count", "min", "max", "mean"]).reset_index().to_dict("records")
    )
    all_unbounded_large = sum("unbounded_delay_ge_10" in row["review_reasons"] for row in review)
    summary = {
        "status": "COMPUTATIONAL_AND_SOURCE_SEMANTICS_REVIEW_COMPLETE",
        "independent_human_kinase_expert_attestation": False,
        "pkis1_multi_construct_parent_groups": sum(row["panel"] == "PKIS1" for row in groups),
        "pkis2_multi_construct_parent_groups": sum(row["panel"] == "PKIS2" for row in groups),
        "pkis1_construct_discordance": discordance,
        "threshold_boundary_counts": boundaries,
        "external_candidate_completeness": candidates,
        "flagged_case_rows": len(review),
        "flagged_first_hit_source_record_rows": len(source_records),
        "unbounded_large_delay_case_condition_rows": all_unbounded_large,
        "limitations": [
            "This is a reproducible literature/data audit, not independent human expert sign-off.",
            "Parent maximum collapse defines an any-profiled-construct endpoint and does not establish construct interchangeability.",
            "Single-concentration percent inhibition and relative-affinity events are operational labels, not potency, cellular engagement, toxicity, or therapeutic-window claims.",
        ],
    }
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "construct_groups.csv", groups)
    write_csv(output / "pkis1_construct_discordance.csv", discordance)
    write_csv(output / "threshold_boundaries.csv", boundaries)
    write_csv(output / "external_flagged_cases.csv", review)
    write_csv(output / "flagged_case_first_hit_source_records.csv", source_records)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
