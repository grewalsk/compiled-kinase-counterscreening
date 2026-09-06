#!/usr/bin/env python3
"""Complete post hoc PKIS2 strict-boundary inclusivity audit."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.Chem.MolStandardize import rdMolStandardize

from run_compiled_coverage import (
    Benchmark, PKIS2_METADATA, PKIS2_SHA256, holm_adjust,
    largest_fragment_smiles, molecular_arrays, pkis2_parent_target, sha256, write_csv,
)
from run_pkis1_external_validation import _comparison, run_condition


def matrices(path: Path):
    if sha256(path) != PKIS2_SHA256:
        raise RuntimeError("PKIS2 SHA-256 mismatch")
    frame = pd.read_excel(path, sheet_name="Table 4 - PKIS2 %Inh")
    assays = sorted(str(column) for column in frame.columns if column not in PKIS2_METADATA)
    empty = frame[["Regno", "Compound", "Chemotype", "Smiles"]].isna().all(axis=1)
    frame = frame.loc[~empty].copy()
    numeric = frame[assays].apply(pd.to_numeric, errors="raise")
    chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    frame["parent_smiles"] = [largest_fragment_smiles(str(x), chooser) for x in frame["Smiles"]]
    groups = [(str(s), g.index.tolist()) for s, g in frame.groupby("parent_smiles", sort=True)]
    construct = np.vstack([
        numeric.loc[indices].median(axis=0, skipna=True).to_numpy(float)
        for _, indices in groups
    ])
    parent_constructs: dict[str, list[str]] = defaultdict(list)
    for assay in assays:
        parent_constructs[pkis2_parent_target(assay)].append(assay)
    parents = sorted(parent_constructs)
    collapsed = pd.DataFrame({
        parent: numeric[columns].max(axis=1, skipna=True)
        for parent, columns in parent_constructs.items()
    })[parents]
    parent = np.vstack([
        collapsed.loc[indices].median(axis=0, skipna=True).to_numpy(float)
        for _, indices in groups
    ])
    smiles = [s for s, _ in groups]
    source_ids = [
        ";".join(sorted(str(x) for x in frame.loc[indices, "Compound"]))
        for _, indices in groups
    ]
    chemotypes = [
        frozenset(str(x) for x in frame.loc[indices, "Chemotype"].dropna().unique())
        for _, indices in groups
    ]
    if parent.shape != (640, 392) or construct.shape != (640, 406):
        raise RuntimeError(f"Unexpected matrices {parent.shape}/{construct.shape}")
    return smiles, source_ids, chemotypes, parent, parents, construct, assays


def benchmark(smiles, source_ids, chemotypes, values, names, resolution, threshold):
    labels = np.where(np.isnan(values), np.nan, (values >= threshold).astype(float))
    prefix = "P" if resolution == "parent" else "A"
    return Benchmark(
        name=f"PKIS2_{resolution.upper()}_GE{threshold}",
        compound_keys=[f"X{i:04d}" for i in range(1, 641)],
        canonical_smiles=smiles,
        target_keys=[f"{prefix}{i:04d}" for i in range(1, len(names) + 1)],
        target_names=names,
        labels={"activity": labels},
        candidate_masks=~np.isnan(labels),
        chemotypes=chemotypes,
        source_ids=source_ids,
        metadata={"resolution": resolution, "inclusive_threshold": threshold},
    )


def parse_args():
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


def main():
    args = parse_args()
    smiles, source_ids, chemotypes, parent, parents, construct, assays = matrices(args.pkis2_xlsx)
    rows = []
    for resolution, values, names in [("parent", parent, parents), ("construct", construct, assays)]:
        for threshold in [80, 90]:
            data = benchmark(smiles, source_ids, chemotypes, values, names, resolution, threshold)
            similarities, scaffolds = molecular_arrays(data)
            for exclusion in ["standard", "chemotype"]:
                condition = f"{resolution}__ge{threshold}__{exclusion}"
                result, stability = run_condition(
                    data, similarities, scaffolds, condition, "activity", exclusion,
                    args.workers, 0,
                )
                if stability:
                    raise AssertionError("Unexpected stability rows")
                rows.extend(result)
    indexed = {}
    for row in rows:
        indexed.setdefault((row["condition"], row["method"]), {})[row["compound_key"]] = row
    comparisons = []
    for offset, condition in enumerate(sorted({row["condition"] for row in rows})):
        summary, _ = _comparison(
            indexed[(condition, "BOUNDED_COVERAGE")], indexed[(condition, "MARGINAL")],
            condition, "bounded_minus_marginal", args.bootstrap, args.permutations,
            900 + offset,
        )
        comparisons.append({
            **summary,
            "ci95_low": summary["ci95"][0],
            "ci95_high": summary["ci95"][1],
        })
    adjusted = holm_adjust({r["condition"]: r["permutation_p_two_sided"] for r in comparisons})
    for row in comparisons:
        row["p_holm_8"] = adjusted[row["condition"]]
    if any(row["large_delay_ge_10"] for row in comparisons):
        raise AssertionError("Bounded policy violated guarantee")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "case_metrics.csv.gz", rows)
    write_csv(output / "comparison_summary.csv", comparisons)
    (output / "summary.json").write_text(json.dumps({
        "status": "EXPLORATORY_POSTHOC_BOUNDARY_AUDIT_COMPLETE",
        "label_rule": "greater than or equal to boundary",
        "complete_grid": True,
        "bounded_effects": comparisons,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "output_manifest.json").write_text(
        json.dumps(output_manifest(output), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
