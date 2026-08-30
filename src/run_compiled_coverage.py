#!/usr/bin/env python3
"""Run the compiled profile-coverage kinase-screen study.

The runner is intentionally CPU-only. It evaluates a lossless all-negative
compilation of stop-on-first-hit policies, correlation-aware profile coverage,
strong chemical-series exclusions, exact small-instance optima, support
abstention diagnostics, and fixed-order latency-aware batching.

It keeps the Klaeger/ChEMBL Kd task and PKIS2 single-concentration task
strictly separate. Nothing in this program estimates clinical safety,
toxicity, efficacy, or therapeutic suitability.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold


METHODS = [
    "RANDOM_EXPECTED",
    "PREVALENCE",
    "CHEMICAL_KNN",
    "WEIGHTED_MARGINAL",
    "GLOBAL_COVERAGE",
    "WEIGHTED_COVERAGE",
    "COMPILED_PARTICLE_MYOPIC",
    "COMPILED_PARTICLE_ROLLOUT_2",
    "ORACLE",
]
FROZEN_BUDGETS = [1, 3, 5, 10, 20]
SUPPORT_THRESHOLDS = [0.2, 0.3, 0.4, 0.5]
LAMBDA_COST_RATIOS = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
PKIS2_SHA256 = "48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c"
PKIS2_METADATA = ["Regno", "Compound", "Chemotype", "Smiles", ">90", ">80", ">70"]
PRIMARY_METHOD = "WEIGHTED_COVERAGE"
PRIMARY_BASELINE = "WEIGHTED_MARGINAL"


@dataclass
class Benchmark:
    name: str
    compound_keys: list[str]
    canonical_smiles: list[str]
    target_keys: list[str]
    target_names: list[str]
    labels: dict[str, np.ndarray]
    candidate_masks: np.ndarray
    chemotypes: list[frozenset[str]]
    source_ids: list[str]
    metadata: dict


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: Iterable[dict], fields: list[str] | None = None) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fields is None:
        fields = list(rows[0].keys())
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stable_rank(scores: np.ndarray, candidates: np.ndarray, target_keys: list[str]) -> list[int]:
    return sorted(candidates.tolist(), key=lambda target: (-float(scores[target]), target_keys[target]))


def normalized_similarity_weights(similarities: np.ndarray, beta: float = 8.0) -> np.ndarray:
    if len(similarities) == 0:
        return np.array([], dtype=float)
    shifted = beta * (similarities - float(np.max(similarities)))
    raw = np.exp(shifted)
    return raw / raw.sum()


def weighted_probabilities(weights: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Beta(1,1)-smoothed marginals, matching the completed particle study."""
    if len(weights) == 0:
        return np.full(labels.shape[1], 0.5, dtype=float)
    normalized = weights / weights.sum()
    effective = normalized * len(weights)
    observed = ~np.isnan(labels)
    positive = np.nan_to_num(labels, nan=0.0)
    numerator = 1.0 + effective @ positive
    denominator = 2.0 + effective @ observed.astype(float)
    return numerator / denominator


def prevalence_scores(labels: np.ndarray) -> np.ndarray:
    observed = ~np.isnan(labels)
    return (1.0 + np.nansum(labels, axis=0)) / (2.0 + observed.sum(axis=0))


def chemical_knn_scores(
    labels: np.ndarray, similarities: np.ndarray, k: int = 10
) -> np.ndarray:
    order = np.argsort(-similarities, kind="stable")[: min(k, len(similarities))]
    selected = labels[order]
    weights = similarities[order]
    observed = ~np.isnan(selected)
    numerator = 1.0 + np.nansum(selected * weights[:, None], axis=0)
    denominator = 2.0 + np.sum(observed * weights[:, None], axis=0)
    return numerator / denominator


def coverage_order(
    labels: np.ndarray,
    weights: np.ndarray,
    candidates: np.ndarray,
    target_keys: list[str],
    max_budget: int,
    fallback_scores: np.ndarray | None = None,
) -> list[int]:
    """Natural weighted min-sum/max-coverage greedy order.

    A reference profile leaves the uncovered pool after the first selected
    target on which it has a measured positive. Missing values never cover a
    profile. The fixed fallback only resolves zero/equal-gain optima.
    """
    positive = np.nan_to_num(labels, nan=0.0) == 1.0
    remaining = candidates.copy()
    uncovered = np.ones(len(labels), dtype=bool)
    fallback = (
        weighted_probabilities(weights, labels)
        if fallback_scores is None
        else np.asarray(fallback_scores, dtype=float)
    )
    # Maintain marginal gains incrementally. Each profile is subtracted once,
    # when it is first covered, reducing a 20-pass matrix scan to roughly two
    # profile-matrix passes without changing the greedy rule.
    gains_all = weights @ positive
    order: list[int] = []
    for _ in range(min(max_budget, len(remaining))):
        gains = gains_all[remaining]
        ranked = sorted(
            zip(gains.tolist(), fallback[remaining].tolist(), remaining.tolist()),
            key=lambda item: (-round(item[0], 14), -item[1], target_keys[item[2]]),
        )
        chosen = ranked[0][2]
        order.append(chosen)
        newly_covered = uncovered & positive[:, chosen]
        if np.any(newly_covered):
            gains_all -= weights[newly_covered] @ positive[newly_covered]
            gains_all[np.abs(gains_all) < 1e-14] = 0.0
        uncovered[newly_covered] = False
        remaining = remaining[remaining != chosen]
    return order


def update_particle(weights: np.ndarray, label_column: np.ndarray, observation: int) -> np.ndarray:
    missing = np.isnan(label_column)
    agreement = np.equal(label_column, observation, where=~missing)
    likelihood = np.where(missing, 1.0, np.where(agreement, 0.9, 0.1))
    updated = weights * likelihood
    total = updated.sum()
    if not np.isfinite(total) or total <= 0:
        return np.full(len(weights), 1.0 / len(weights))
    return updated / total


def particle_sequence(
    labels: np.ndarray,
    initial_weights: np.ndarray,
    candidates: np.ndarray,
    target_keys: list[str],
    max_budget: int,
    rollout: bool,
    true_labels: np.ndarray | None,
) -> list[int]:
    """Run online with true labels, or compile by setting true_labels=None."""
    weights = initial_weights.copy()
    remaining = candidates.copy()
    sequence: list[int] = []
    for _ in range(min(max_budget, len(remaining))):
        probabilities = weighted_probabilities(weights, labels)
        if not rollout or len(remaining) == 1:
            chosen = stable_rank(probabilities, remaining, target_keys)[0]
        else:
            sublabels = labels[:, remaining]
            missing = np.isnan(sublabels)
            negative_likelihood = np.where(
                missing, 1.0, np.where(sublabels == 0.0, 0.9, 0.1)
            )
            negative_weights = negative_likelihood.T * weights[None, :]
            totals = negative_weights.sum(axis=1, keepdims=True)
            negative_weights = np.divide(
                negative_weights,
                totals,
                out=np.full_like(negative_weights, 1.0 / len(weights)),
                where=totals > 0,
            )
            effective = negative_weights * len(weights)
            observed = (~np.isnan(sublabels)).astype(float)
            positive = np.nan_to_num(sublabels, nan=0.0)
            second_probabilities = (1.0 + effective @ positive) / (2.0 + effective @ observed)
            np.fill_diagonal(second_probabilities, -np.inf)
            best_second = np.max(second_probabilities, axis=1)
            first_probability = probabilities[remaining]
            utility = first_probability + (1.0 - first_probability) * np.maximum(best_second, 0.0)
            chosen = sorted(
                zip(utility.tolist(), remaining.tolist()),
                key=lambda item: (-item[0], target_keys[item[1]]),
            )[0][1]
        sequence.append(chosen)
        observation = 0 if true_labels is None else int(true_labels[chosen])
        weights = update_particle(weights, labels[:, chosen], observation)
        remaining = remaining[remaining != chosen]
    return sequence


def metrics_from_sequence(sequence: list[int], truth: np.ndarray, max_budget: int) -> dict:
    observed = [int(truth[target]) for target in sequence[:max_budget]]
    cumulative = np.cumsum(observed) if observed else np.array([], dtype=int)
    any_curve = [
        int(len(cumulative) > 0 and cumulative[min(b, len(cumulative)) - 1] > 0)
        if b > 0
        else 0
        for b in range(1, max_budget + 1)
    ]
    first_hit = next((step + 1 for step, value in enumerate(observed) if value), max_budget + 1)
    audc = float(np.mean(any_curve))
    identity = (max_budget + 1 - first_hit) / max_budget
    if not math.isclose(audc, identity, abs_tol=1e-12):
        raise AssertionError((audc, identity, first_hit))
    return {
        "audc": audc,
        "first_hit": int(first_hit),
        "any_curve": any_curve,
    }


def random_any_hit(n: int, positives: int, budget: int) -> float:
    b = min(budget, n)
    if positives <= 0 or b <= 0:
        return 0.0
    if b > n - positives:
        return 1.0
    return 1.0 - math.comb(n - positives, b) / math.comb(n, b)


def random_metrics(n: int, positives: int, max_budget: int) -> dict:
    curve = [random_any_hit(n, positives, budget) for budget in range(1, max_budget + 1)]
    expected_first = float(sum(1.0 - random_any_hit(n, positives, b) for b in range(max_budget + 1)))
    audc = float(np.mean(curve))
    if not math.isclose(audc, (max_budget + 1 - expected_first) / max_budget, abs_tol=1e-12):
        raise AssertionError("Random expectation violates the AUDC/first-hit identity")
    return {"audc": audc, "first_hit": expected_first, "any_curve": curve}


def load_klaeger(data_dir: Path) -> Benchmark:
    cases = [
        json.loads(line)
        for line in (data_dir / "cases.jsonl").read_text(encoding="utf-8").splitlines()
        if line
    ]
    with gzip.open(data_dir / "gold.csv.gz", "rt", encoding="utf-8", newline="") as stream:
        gold = list(csv.DictReader(stream))
    targets = pd.read_csv(data_dir / "target_mapping.csv")
    target_keys = sorted(targets.loc[targets["in_common_panel"] == 1, "target_key"].tolist())
    target_names_by_key = dict(zip(targets["target_key"], targets["target_name"]))
    target_names = [str(target_names_by_key[key]) for key in target_keys]
    compound_keys = [case["compound_key"] for case in cases]
    compound_index = {key: index for index, key in enumerate(compound_keys)}
    target_index = {key: index for index, key in enumerate(target_keys)}
    shape = (len(cases), len(target_keys))
    primary = np.full(shape, np.nan, dtype=float)
    sensitivity = np.full(shape, np.nan, dtype=float)
    candidate_masks = np.zeros(shape, dtype=bool)
    for case in cases:
        i = compound_index[case["compound_key"]]
        for key in case["candidate_target_keys"]:
            if key in target_index:
                candidate_masks[i, target_index[key]] = True
    for row in gold:
        i = compound_index[row["compound_key"]]
        t = target_index[row["target_key"]]
        primary[i, t] = int(row["liability_delta_1"])
        sensitivity[i, t] = int(row["liability_delta_0_5"])
    if np.any(candidate_masks & np.isnan(primary)):
        raise RuntimeError("Klaeger candidate mask contains missing primary labels")
    return Benchmark(
        name="KLAEGER_CHEMBL30",
        compound_keys=compound_keys,
        canonical_smiles=[case["canonical_smiles"] for case in cases],
        target_keys=target_keys,
        target_names=target_names,
        labels={"primary": primary, "sensitivity": sensitivity},
        candidate_masks=candidate_masks,
        chemotypes=[frozenset() for _ in cases],
        source_ids=compound_keys,
        metadata={
            "source": "Klaeger et al. via ChEMBL 30",
            "license": "CC BY-SA 3.0",
            "task_language": "within-panel binding liability",
        },
    )


def pkis2_parent_target(name: str) -> str:
    name = re.sub(
        r"-(?:nonphosphorylated|phosphorylated|autoinhibited|cyclinD[13])$", "", name
    )
    return re.sub(r"\(.*\)$", "", name)


def largest_fragment_smiles(smiles: str, chooser: rdMolStandardize.LargestFragmentChooser) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Unparseable PKIS2 SMILES: {smiles}")
    parent = chooser.choose(molecule)
    return Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)


def load_pkis2(path: Path) -> Benchmark:
    observed_hash = sha256(path)
    if observed_hash != PKIS2_SHA256:
        raise RuntimeError(
            f"Unexpected PKIS2 workbook SHA-256: {observed_hash}; expected {PKIS2_SHA256}"
        )
    frame = pd.read_excel(path, sheet_name="Table 4 - PKIS2 %Inh")
    missing_metadata = sorted(set(PKIS2_METADATA) - set(frame.columns))
    if missing_metadata:
        raise RuntimeError(f"PKIS2 metadata columns missing: {missing_metadata}")
    assay_columns = [str(column) for column in frame.columns if column not in PKIS2_METADATA]
    empty_identifier_mask = frame[["Regno", "Compound", "Chemotype", "Smiles"]].isna().all(axis=1)
    excluded_cells = frame.loc[empty_identifier_mask, assay_columns].stack()
    if (
        int(empty_identifier_mask.sum()) != 1
        or len(excluded_cells) != 1
        or str(excluded_cells.index[0][1]) != "TEC"
        or not math.isclose(float(excluded_cells.iloc[0]), 13.0)
    ):
        raise RuntimeError("Unexpected PKIS2 identifier-empty row or stray assay cells")
    frame = frame.loc[~empty_identifier_mask].copy()
    if frame["Compound"].isna().any() or frame["Smiles"].isna().any():
        raise RuntimeError("PKIS2 contains a partially populated compound row")
    parent_to_constructs: dict[str, list[str]] = defaultdict(list)
    for column in assay_columns:
        parent_to_constructs[pkis2_parent_target(column)].append(column)
    parent_names = sorted(parent_to_constructs)
    if len(assay_columns) != 406 or len(parent_names) != 392:
        raise RuntimeError(
            f"Unexpected PKIS2 assay/parent counts: {len(assay_columns)}/{len(parent_names)}"
        )
    numeric = frame[assay_columns].apply(pd.to_numeric, errors="raise")
    if ((numeric < 0) | (numeric > 100)).any().any():
        raise RuntimeError("PKIS2 percent-inhibition values fall outside [0,100]")
    collapsed = pd.DataFrame(
        {
            parent: numeric[constructs].max(axis=1, skipna=True)
            for parent, constructs in parent_to_constructs.items()
        }
    )[parent_names]

    chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    frame["parent_smiles"] = [largest_fragment_smiles(str(value), chooser) for value in frame["Smiles"]]
    grouped_indices = [
        (parent_smiles, group.index.tolist())
        for parent_smiles, group in frame.groupby("parent_smiles", sort=True)
    ]
    values = np.vstack(
        [collapsed.loc[indices].median(axis=0, skipna=True).to_numpy(float) for _, indices in grouped_indices]
    )
    primary = np.where(np.isnan(values), np.nan, (values > 90.0).astype(float))
    sensitivity = np.where(np.isnan(values), np.nan, (values > 80.0).astype(float))
    candidate_masks = ~np.isnan(primary)
    target_keys = [f"K{index:04d}" for index in range(1, len(parent_names) + 1)]
    compound_keys = [f"X{index:04d}" for index in range(1, len(grouped_indices) + 1)]
    source_ids: list[str] = []
    chemotypes: list[frozenset[str]] = []
    parent_smiles: list[str] = []
    for smiles, indices in grouped_indices:
        parent_smiles.append(smiles)
        source_ids.append(";".join(sorted(str(value) for value in frame.loc[indices, "Compound"])))
        chemotypes.append(
            frozenset(str(value) for value in frame.loc[indices, "Chemotype"].dropna().unique())
        )
    return Benchmark(
        name="PKIS2",
        compound_keys=compound_keys,
        canonical_smiles=parent_smiles,
        target_keys=target_keys,
        target_names=parent_names,
        labels={"primary": primary, "sensitivity": sensitivity},
        candidate_masks=candidate_masks,
        chemotypes=chemotypes,
        source_ids=source_ids,
        metadata={
            "source": "Drewry et al. 2017 PLOS ONE S4",
            "source_sha256": observed_hash,
            "license": "CC BY 4.0",
            "raw_rows_after_empty_exclusion": int(len(frame)),
            "excluded_identifier_empty_rows": int(empty_identifier_mask.sum()),
            "excluded_stray_assay_cell": "TEC=13",
            "parent_compounds": int(len(grouped_indices)),
            "raw_assay_columns": len(assay_columns),
            "parent_targets": len(parent_names),
            "task_language": "first strong profile-activity discovery; intended targets unavailable",
        },
    )


def subset_benchmark(benchmark: Benchmark, limit: int) -> Benchmark:
    indices = np.arange(min(limit, len(benchmark.compound_keys)))
    return Benchmark(
        name=benchmark.name,
        compound_keys=[benchmark.compound_keys[i] for i in indices],
        canonical_smiles=[benchmark.canonical_smiles[i] for i in indices],
        target_keys=benchmark.target_keys,
        target_names=benchmark.target_names,
        labels={key: value[indices] for key, value in benchmark.labels.items()},
        candidate_masks=benchmark.candidate_masks[indices],
        chemotypes=[benchmark.chemotypes[i] for i in indices],
        source_ids=[benchmark.source_ids[i] for i in indices],
        metadata={**benchmark.metadata, "smoke_test_limit": len(indices)},
    )


def molecular_arrays(benchmark: Benchmark):
    molecules = [Chem.MolFromSmiles(smiles) for smiles in benchmark.canonical_smiles]
    if any(molecule is None for molecule in molecules):
        raise RuntimeError(f"RDKit failed to parse a {benchmark.name} canonical SMILES")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = [generator.GetFingerprint(molecule) for molecule in molecules]
    n = len(fingerprints)
    similarities = np.eye(n, dtype=float)
    for i in range(n):
        values = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[i + 1 :])
        similarities[i, i + 1 :] = values
        similarities[i + 1 :, i] = values
    scaffolds = [
        MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=False)
        for molecule in molecules
    ]
    return similarities, scaffolds


def reference_indices(
    heldout: int,
    similarities: np.ndarray,
    scaffolds: list[str],
    chemotypes: list[frozenset[str]],
    exclusion: str,
) -> np.ndarray:
    refs = np.array(
        [
            index
            for index in range(len(similarities))
            if index != heldout and similarities[heldout, index] < 1.0 - 1e-12
        ],
        dtype=int,
    )
    if exclusion == "scaffold":
        refs = np.array([index for index in refs if scaffolds[index] != scaffolds[heldout]], dtype=int)
    elif exclusion == "nearest10":
        ranked = sorted(refs.tolist(), key=lambda index: (-similarities[heldout, index], index))
        refs = np.array(ranked[10:], dtype=int)
    elif exclusion == "chemotype":
        heldout_groups = chemotypes[heldout]
        refs = np.array(
            [index for index in refs if not (heldout_groups & chemotypes[index])], dtype=int
        )
    elif exclusion != "standard":
        raise ValueError(exclusion)
    return refs


def expected_cover_time(order: list[int], labels: np.ndarray, weights: np.ndarray) -> float:
    uncovered = np.ones(len(labels), dtype=bool)
    total = 1.0
    for target in order:
        uncovered &= ~(np.nan_to_num(labels[:, target], nan=0.0) == 1.0)
        total += float(weights[uncovered].sum())
    return total


def exact_restricted_order(
    restricted_targets: list[int], labels: np.ndarray, weights: np.ndarray
) -> tuple[list[int], float]:
    """Exact min-sum order on a small target restriction via subset DP."""
    k = len(restricted_targets)
    positive = np.nan_to_num(labels[:, restricted_targets], nan=0.0) == 1.0
    survival = np.empty(1 << k, dtype=float)
    survival[0] = 1.0
    for mask in range(1, 1 << k):
        selected = [bit for bit in range(k) if mask & (1 << bit)]
        hit = np.any(positive[:, selected], axis=1)
        survival[mask] = float(weights[~hit].sum())
    dp = np.full(1 << k, np.inf, dtype=float)
    previous = np.full(1 << k, -1, dtype=int)
    added = np.full(1 << k, -1, dtype=int)
    dp[0] = 0.0
    for mask in range(1 << k):
        for bit in range(k):
            if mask & (1 << bit):
                continue
            new_mask = mask | (1 << bit)
            candidate = dp[mask] + survival[new_mask]
            if candidate < dp[new_mask] - 1e-15:
                dp[new_mask] = candidate
                previous[new_mask] = mask
                added[new_mask] = bit
    mask = (1 << k) - 1
    reverse_bits: list[int] = []
    while mask:
        reverse_bits.append(int(added[mask]))
        mask = int(previous[mask])
    order = [restricted_targets[bit] for bit in reversed(reverse_bits)]
    return order, float(1.0 + dp[-1])


def survival_curve(order: list[int], labels: np.ndarray, weights: np.ndarray) -> np.ndarray:
    uncovered = np.ones(len(labels), dtype=bool)
    survival = [1.0]
    for target in order:
        uncovered &= ~(np.nan_to_num(labels[:, target], nan=0.0) == 1.0)
        survival.append(float(weights[uncovered].sum()))
    return np.asarray(survival, dtype=float)


def optimal_batch_partition(
    order: list[int], labels: np.ndarray, weights: np.ndarray, lambda_ratio: float
) -> dict:
    q = survival_curve(order, labels, weights)
    b = len(order)
    prefix_cost = np.arange(b + 1, dtype=float)
    dp = np.full(b + 1, np.inf, dtype=float)
    previous = np.full(b + 1, -1, dtype=int)
    dp[0] = 0.0
    for j in range(1, b + 1):
        for i in range(j):
            candidate = dp[i] + q[i] * (lambda_ratio + prefix_cost[j] - prefix_cost[i])
            if candidate < dp[j] - 1e-15 or (
                math.isclose(candidate, dp[j], abs_tol=1e-15) and i > previous[j]
            ):
                dp[j] = candidate
                previous[j] = i
    endpoints: list[int] = []
    cursor = b
    while cursor > 0:
        endpoints.append(cursor)
        cursor = int(previous[cursor])
    endpoints.reverse()
    starts = [0] + endpoints[:-1]
    sizes = [end - start for start, end in zip(starts, endpoints)]
    return {
        "modeled_expected_cost": float(dp[b]),
        "batch_endpoints": endpoints,
        "batch_sizes": sizes,
        "survival": q.tolist(),
    }


def evaluate_condition(
    benchmark: Benchmark,
    condition: str,
    label_name: str,
    exclusion: str,
    similarities: np.ndarray,
    scaffolds: list[str],
    max_budget: int,
    run_particle: bool,
) -> tuple[list[dict], list[dict], list[dict], list[dict], int]:
    labels = benchmark.labels[label_name]
    case_rows: list[dict] = []
    traces: list[dict] = []
    exact_rows: list[dict] = []
    batch_rows: list[dict] = []
    compiler_checks = 0
    for heldout, compound_key in enumerate(benchmark.compound_keys):
        refs = reference_indices(
            heldout, similarities, scaffolds, benchmark.chemotypes, exclusion
        )
        if len(refs) < 5:
            continue
        candidates = np.flatnonzero(benchmark.candidate_masks[heldout] & ~np.isnan(labels[heldout]))
        if len(candidates) < max_budget:
            continue
        train = labels[refs]
        truth = labels[heldout]
        sims = similarities[heldout, refs]
        chemical_weights = normalized_similarity_weights(sims)
        uniform_weights = np.full(len(refs), 1.0 / len(refs))
        weighted_marginals = weighted_probabilities(chemical_weights, train)
        sequences: dict[str, list[int]] = {
            "PREVALENCE": stable_rank(prevalence_scores(train), candidates, benchmark.target_keys)[
                :max_budget
            ],
            "CHEMICAL_KNN": stable_rank(
                chemical_knn_scores(train, sims), candidates, benchmark.target_keys
            )[:max_budget],
            "WEIGHTED_MARGINAL": stable_rank(
                weighted_marginals, candidates, benchmark.target_keys
            )[:max_budget],
            "GLOBAL_COVERAGE": coverage_order(
                train,
                uniform_weights,
                candidates,
                benchmark.target_keys,
                max_budget,
                prevalence_scores(train),
            ),
            "WEIGHTED_COVERAGE": coverage_order(
                train,
                chemical_weights,
                candidates,
                benchmark.target_keys,
                max_budget,
                weighted_marginals,
            ),
            "ORACLE": sorted(
                candidates.tolist(), key=lambda target: (-int(truth[target]), benchmark.target_keys[target])
            )[:max_budget],
        }
        if run_particle:
            for method, rollout in [
                ("COMPILED_PARTICLE_MYOPIC", False),
                ("COMPILED_PARTICLE_ROLLOUT_2", True),
            ]:
                compiled = particle_sequence(
                    train,
                    chemical_weights,
                    candidates,
                    benchmark.target_keys,
                    max_budget,
                    rollout,
                    None,
                )
                online = particle_sequence(
                    train,
                    chemical_weights,
                    candidates,
                    benchmark.target_keys,
                    max_budget,
                    rollout,
                    truth,
                )
                compiled_metrics = metrics_from_sequence(compiled, truth, max_budget)
                online_metrics = metrics_from_sequence(online, truth, max_budget)
                prefix_length = min(compiled_metrics["first_hit"], max_budget)
                if compiled[:prefix_length] != online[:prefix_length]:
                    raise AssertionError(
                        f"Compiler prefix mismatch: {condition}/{compound_key}/{method}"
                    )
                if compiled_metrics != online_metrics:
                    raise AssertionError(
                        f"Compiler metric mismatch: {condition}/{compound_key}/{method}"
                    )
                sequences[method] = compiled
                compiler_checks += 1

        support = {
            "reference_count": len(refs),
            "nearest_reference_similarity": float(np.max(sims)),
            "effective_reference_count": float(1.0 / np.sum(np.square(chemical_weights))),
            "positive_reference_weight": float(
                chemical_weights[
                    np.any(np.nan_to_num(train[:, candidates], nan=0.0) == 1.0, axis=1)
                ].sum()
            ),
        }
        for method, sequence in sequences.items():
            metrics = metrics_from_sequence(sequence, truth, max_budget)
            row = {
                "dataset": benchmark.name,
                "condition": condition,
                "compound_key": compound_key,
                "source_ids": benchmark.source_ids[heldout],
                "policy": method,
                "candidate_count": len(candidates),
                "total_hidden_positives": int(np.nansum(truth[candidates])),
                "audc_1_20": metrics["audc"],
                "cost_to_first_censored": metrics["first_hit"],
                **support,
            }
            for budget in FROZEN_BUDGETS:
                row[f"any_hit_b{budget}"] = metrics["any_curve"][budget - 1]
            case_rows.append(row)
            for step, target in enumerate(sequence, start=1):
                traces.append(
                    {
                        "dataset": benchmark.name,
                        "condition": condition,
                        "compound_key": compound_key,
                        "policy": method,
                        "step": step,
                        "target_key": benchmark.target_keys[target],
                        "target_name": benchmark.target_names[target],
                        "hidden_positive": int(truth[target]),
                    }
                )

        positives = int(np.nansum(truth[candidates]))
        random_result = random_metrics(len(candidates), positives, max_budget)
        random_row = {
            "dataset": benchmark.name,
            "condition": condition,
            "compound_key": compound_key,
            "source_ids": benchmark.source_ids[heldout],
            "policy": "RANDOM_EXPECTED",
            "candidate_count": len(candidates),
            "total_hidden_positives": positives,
            "audc_1_20": random_result["audc"],
            "cost_to_first_censored": random_result["first_hit"],
            **support,
        }
        for budget in FROZEN_BUDGETS:
            random_row[f"any_hit_b{budget}"] = random_result["any_curve"][budget - 1]
        case_rows.append(random_row)

        if benchmark.name == "KLAEGER_CHEMBL30" and condition == "klaeger_primary_standard":
            restricted = stable_rank(weighted_marginals, candidates, benchmark.target_keys)[:8]
            greedy = coverage_order(
                train,
                chemical_weights,
                np.asarray(restricted, dtype=int),
                benchmark.target_keys,
                8,
                weighted_marginals,
            )
            exact, exact_cost = exact_restricted_order(restricted, train, chemical_weights)
            greedy_cost = expected_cover_time(greedy, train, chemical_weights)
            if greedy_cost < exact_cost - 1e-10:
                raise AssertionError("Greedy restricted cost beat exact optimum")
            exact_rows.append(
                {
                    "compound_key": compound_key,
                    "restricted_target_count": 8,
                    "greedy_expected_cover_time": greedy_cost,
                    "exact_expected_cover_time": exact_cost,
                    "absolute_gap": greedy_cost - exact_cost,
                    "relative_gap": (greedy_cost / exact_cost - 1.0) if exact_cost else 0.0,
                    "greedy_order": ";".join(benchmark.target_keys[target] for target in greedy),
                    "exact_order": ";".join(benchmark.target_keys[target] for target in exact),
                }
            )

        if label_name == "primary" and exclusion == "standard":
            weighted_order = sequences["WEIGHTED_COVERAGE"]
            for lambda_ratio in LAMBDA_COST_RATIOS:
                schedule = optimal_batch_partition(
                    weighted_order, train, chemical_weights, lambda_ratio
                )
                batch_rows.append(
                    {
                        "dataset": benchmark.name,
                        "condition": condition,
                        "compound_key": compound_key,
                        "lambda_over_assay_cost": lambda_ratio,
                        "modeled_expected_cost": schedule["modeled_expected_cost"],
                        "batch_count": len(schedule["batch_sizes"]),
                        "batch_sizes": ";".join(map(str, schedule["batch_sizes"])),
                        "batch_endpoints": ";".join(map(str, schedule["batch_endpoints"])),
                        "q_after_20": schedule["survival"][-1],
                    }
                )
    return case_rows, traces, exact_rows, batch_rows, compiler_checks


def paired_bootstrap(values: np.ndarray, replicates: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n = len(values)
    estimates = np.empty(replicates, dtype=float)
    for start in range(0, replicates, 1000):
        size = min(1000, replicates - start)
        indices = rng.integers(0, n, size=(size, n))
        estimates[start : start + size] = values[indices].mean(axis=1)
    return {
        "estimate": float(values.mean()),
        "ci95": [float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))],
        "replicates": replicates,
        "seed": seed,
    }


def sign_flip_test(values: np.ndarray, permutations: int, seed: int) -> float:
    nonzero = values[values != 0]
    if len(nonzero) == 0:
        return 1.0
    observed = abs(float(nonzero.mean()))
    rng = np.random.default_rng(seed)
    exceed = 0
    completed = 0
    while completed < permutations:
        size = min(2000, permutations - completed)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(size, len(nonzero)))
        exceed += int(np.sum(np.abs((signs * nonzero).mean(axis=1)) >= observed - 1e-15))
        completed += size
    return float((exceed + 1) / (permutations + 1))


def holm_adjust(pvalues: dict[str, float]) -> dict[str, float]:
    ordered = sorted(pvalues.items(), key=lambda item: item[1])
    adjusted: dict[str, float] = {}
    running = 0.0
    m = len(ordered)
    for rank, (name, pvalue) in enumerate(ordered):
        running = max(running, min(1.0, (m - rank) * pvalue))
        adjusted[name] = running
    return adjusted


def summarize_results(
    case_rows: list[dict], bootstrap: int, permutations: int
) -> tuple[dict, list[dict], list[dict]]:
    by_key: dict[tuple[str, str, str], dict[str, dict]] = defaultdict(dict)
    for row in case_rows:
        by_key[(row["dataset"], row["condition"], row["policy"])][row["compound_key"]] = row
    conditions = sorted({(row["dataset"], row["condition"]) for row in case_rows})
    tests: dict[str, dict] = {}
    for test_index, (dataset, condition) in enumerate(conditions):
        first = by_key[(dataset, condition, PRIMARY_METHOD)]
        second = by_key[(dataset, condition, PRIMARY_BASELINE)]
        shared = sorted(set(first) & set(second))
        differences = np.array(
            [first[key]["audc_1_20"] - second[key]["audc_1_20"] for key in shared],
            dtype=float,
        )
        name = f"{dataset}/{condition}/{PRIMARY_METHOD}_minus_{PRIMARY_BASELINE}"
        tests[name] = {
            **paired_bootstrap(differences, bootstrap, 20260830 + test_index * 2),
            "permutation_p_two_sided": sign_flip_test(
                differences, permutations, 20260831 + test_index * 2
            ),
            "n": len(shared),
            "comparison": f"{PRIMARY_METHOD} - {PRIMARY_BASELINE}",
        }

    secondary_baselines = [
        "CHEMICAL_KNN",
        "PREVALENCE",
        "GLOBAL_COVERAGE",
        "COMPILED_PARTICLE_MYOPIC",
        "COMPILED_PARTICLE_ROLLOUT_2",
    ]
    secondary: dict[str, dict] = {}
    raw_p: dict[str, float] = {}
    primary_conditions = [
        (dataset, condition)
        for dataset, condition in conditions
        if condition in {"klaeger_primary_standard", "pkis2_gt90_standard"}
    ]
    counter = 0
    for dataset, condition in primary_conditions:
        first = by_key[(dataset, condition, PRIMARY_METHOD)]
        for baseline in secondary_baselines:
            second = by_key.get((dataset, condition, baseline), {})
            shared = sorted(set(first) & set(second))
            if not shared:
                continue
            differences = np.array(
                [first[key]["audc_1_20"] - second[key]["audc_1_20"] for key in shared],
                dtype=float,
            )
            name = f"{dataset}/{condition}/{PRIMARY_METHOD}_minus_{baseline}"
            entry = paired_bootstrap(differences, bootstrap, 20260910 + counter * 2)
            pvalue = sign_flip_test(differences, permutations, 20260911 + counter * 2)
            entry.update(
                {"permutation_p_two_sided": pvalue, "n": len(shared), "comparison": name}
            )
            secondary[name] = entry
            raw_p[name] = pvalue
            counter += 1
    adjusted = holm_adjust(raw_p)
    for name, entry in secondary.items():
        entry["holm_adjusted_p"] = adjusted[name]

    summary_rows: list[dict] = []
    curves: list[dict] = []
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in case_rows:
        grouped[(row["dataset"], row["condition"], row["policy"])].append(row)
    for (dataset, condition, policy), rows in sorted(grouped.items()):
        summary_rows.append(
            {
                "dataset": dataset,
                "condition": condition,
                "policy": policy,
                "n": len(rows),
                "mean_audc_1_20": float(np.mean([row["audc_1_20"] for row in rows])),
                "mean_cost_to_first_censored": float(
                    np.mean([row["cost_to_first_censored"] for row in rows])
                ),
                "positive_case_fraction": float(
                    np.mean([row["total_hidden_positives"] > 0 for row in rows])
                ),
            }
        )
        for budget in FROZEN_BUDGETS:
            curves.append(
                {
                    "dataset": dataset,
                    "condition": condition,
                    "policy": policy,
                    "budget": budget,
                    "n": len(rows),
                    "mean_any_hit": float(np.mean([row[f"any_hit_b{budget}"] for row in rows])),
                }
            )
    return {"primary_tests": tests, "secondary_tests": secondary}, summary_rows, curves


def support_analysis(case_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    primary_conditions = {"klaeger_primary_standard", "pkis2_gt90_standard"}
    for dataset in sorted({row["dataset"] for row in case_rows}):
        data_rows = [
            row
            for row in case_rows
            if row["dataset"] == dataset and row["condition"] in primary_conditions
        ]
        by_policy: dict[str, dict[str, dict]] = defaultdict(dict)
        for row in data_rows:
            by_policy[row["policy"]][row["compound_key"]] = row
        first = by_policy[PRIMARY_METHOD]
        second = by_policy[PRIMARY_BASELINE]
        all_shared = sorted(set(first) & set(second))
        for threshold in SUPPORT_THRESHOLDS:
            shared = [
                key
                for key in all_shared
                if first[key]["nearest_reference_similarity"] >= threshold
            ]
            differences = [
                first[key]["audc_1_20"] - second[key]["audc_1_20"] for key in shared
            ]
            output.append(
                {
                    "dataset": dataset,
                    "minimum_nearest_reference_similarity": threshold,
                    "retained_n": len(shared),
                    "total_n": len(all_shared),
                    "retained_fraction": len(shared) / len(all_shared) if all_shared else None,
                    "mean_audc_difference": float(np.mean(differences)) if differences else None,
                }
            )
    return output


def make_figures(curves: list[dict], statistics: dict, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    primary_conditions = {
        "KLAEGER_CHEMBL30": "klaeger_primary_standard",
        "PKIS2": "pkis2_gt90_standard",
    }
    display_methods = [
        "RANDOM_EXPECTED",
        "PREVALENCE",
        "CHEMICAL_KNN",
        "WEIGHTED_MARGINAL",
        "GLOBAL_COVERAGE",
        "WEIGHTED_COVERAGE",
    ]
    for dataset, condition in primary_conditions.items():
        selected = [
            row for row in curves if row["dataset"] == dataset and row["condition"] == condition
        ]
        if not selected:
            continue
        plt.figure(figsize=(7.2, 4.8))
        for policy in display_methods:
            rows = sorted(
                [row for row in selected if row["policy"] == policy],
                key=lambda row: row["budget"],
            )
            if rows:
                plt.plot(
                    [row["budget"] for row in rows],
                    [row["mean_any_hit"] for row in rows],
                    marker="o",
                    label=policy.replace("_", " ").title(),
                )
        plt.xlabel("Assay budget")
        plt.ylabel("Fraction with at least one measured hit")
        plt.ylim(0, 1)
        plt.xticks(FROZEN_BUDGETS)
        plt.legend(fontsize=7, ncol=2)
        plt.tight_layout()
        plt.savefig(figure_dir / f"{dataset.lower()}_discovery_curve.png", dpi=220)
        plt.close()

    primary_tests = statistics["primary_tests"]
    if primary_tests:
        names = list(primary_tests)
        estimates = [primary_tests[name]["estimate"] for name in names]
        lower = [primary_tests[name]["ci95"][0] for name in names]
        upper = [primary_tests[name]["ci95"][1] for name in names]
        y = np.arange(len(names))
        plt.figure(figsize=(8.2, max(4.2, 0.45 * len(names))))
        plt.errorbar(
            estimates,
            y,
            xerr=[np.array(estimates) - np.array(lower), np.array(upper) - np.array(estimates)],
            fmt="o",
            capsize=3,
        )
        plt.axvline(0.0, color="black", linewidth=1, linestyle="--")
        plt.yticks(y, [name.split("/")[1] for name in names], fontsize=8)
        plt.xlabel("AUDC difference: weighted coverage − weighted marginal")
        plt.tight_layout()
        plt.savefig(figure_dir / "coverage_effects.png", dpi=220)
        plt.close()


def toy_compiler_boundary_check() -> dict:
    """Exhaustive first-hit identity plus a multi-hit branching counterexample."""
    actions = [0, 1, 2, 3]

    def policy(history: list[tuple[int, int]]) -> int:
        remaining = [action for action in actions if action not in {item[0] for item in history}]
        if not history:
            return 0
        if history[-1][1] == 1 and 2 in remaining:
            return 2
        if history[-1][1] == 0 and 1 in remaining:
            return 1
        return remaining[0]

    compiled: list[int] = []
    history: list[tuple[int, int]] = []
    for _ in actions:
        action = policy(history)
        compiled.append(action)
        history.append((action, 0))
    first_hit_equal = 0
    multi_hit_difference_found = False
    for mask in range(1 << len(actions)):
        outcomes = [(mask >> action) & 1 for action in actions]
        online: list[int] = []
        history = []
        for _ in actions:
            action = policy(history)
            online.append(action)
            history.append((action, outcomes[action]))
        online_first = next((i + 1 for i, action in enumerate(online) if outcomes[action]), 5)
        compiled_first = next((i + 1 for i, action in enumerate(compiled) if outcomes[action]), 5)
        if online_first != compiled_first:
            raise AssertionError("Toy compiler first-hit identity failed")
        first_hit_equal += 1
        online_two_hit = sum(outcomes[action] for action in online[:2])
        compiled_two_hit = sum(outcomes[action] for action in compiled[:2])
        multi_hit_difference_found |= online_two_hit != compiled_two_hit
    if not multi_hit_difference_found:
        raise AssertionError("Toy policy failed to demonstrate the multi-hit boundary")
    return {
        "binary_profiles_checked": first_hit_equal,
        "first_hit_identity_pass": True,
        "multi_hit_counterexample_exists": True,
        "compiled_negative_path": compiled,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pkis2-xlsx", type=Path)
    parser.add_argument("--max-budget", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    if args.max_budget != 20:
        raise ValueError("The frozen study requires --max-budget 20")
    started = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    benchmarks = [load_klaeger(args.data_dir)]
    if args.pkis2_xlsx is not None:
        benchmarks.append(load_pkis2(args.pkis2_xlsx))
    if args.smoke_test:
        benchmarks = [
            subset_benchmark(benchmark, 20 if benchmark.name == "KLAEGER_CHEMBL30" else 30)
            for benchmark in benchmarks
        ]
        args.bootstrap = min(args.bootstrap, 250)
        args.permutations = min(args.permutations, 1_000)

    all_case_rows: list[dict] = []
    all_traces: list[dict] = []
    all_exact_rows: list[dict] = []
    all_batch_rows: list[dict] = []
    compiler_checks = 0
    dataset_manifests: dict[str, dict] = {}

    for benchmark in benchmarks:
        similarities, scaffolds = molecular_arrays(benchmark)
        dataset_manifests[benchmark.name] = {
            **benchmark.metadata,
            "compounds": len(benchmark.compound_keys),
            "targets": len(benchmark.target_keys),
            "primary_positive_compounds": int(np.sum(np.nansum(benchmark.labels["primary"], axis=1) > 0)),
            "fingerprint_duplicate_pairs": int(
                np.sum(np.triu(similarities >= 1.0 - 1e-12, k=1))
            ),
        }
        if benchmark.name == "KLAEGER_CHEMBL30":
            conditions = [
                ("klaeger_primary_standard", "primary", "standard", True),
                ("klaeger_delta_0_5", "sensitivity", "standard", True),
                ("klaeger_scaffold_exclusion", "primary", "scaffold", True),
                ("klaeger_leave10_nearest", "primary", "nearest10", True),
            ]
        else:
            conditions = [
                ("pkis2_gt90_standard", "primary", "standard", False),
                ("pkis2_gt80_sensitivity", "sensitivity", "standard", False),
                ("pkis2_gt90_scaffold", "primary", "scaffold", False),
                ("pkis2_gt90_chemotype", "primary", "chemotype", False),
            ]
        for condition, label_name, exclusion, run_particle in conditions:
            rows, traces, exact_rows, batch_rows, checks = evaluate_condition(
                benchmark,
                condition,
                label_name,
                exclusion,
                similarities,
                scaffolds,
                args.max_budget,
                run_particle,
            )
            all_case_rows.extend(rows)
            all_traces.extend(traces)
            all_exact_rows.extend(exact_rows)
            all_batch_rows.extend(batch_rows)
            compiler_checks += checks

    statistics, method_summary, curves = summarize_results(
        all_case_rows, args.bootstrap, args.permutations
    )
    support_rows = support_analysis(all_case_rows)
    toy_check = toy_compiler_boundary_check()
    exact_summary = {
        "n": len(all_exact_rows),
        "mean_absolute_gap": float(np.mean([row["absolute_gap"] for row in all_exact_rows]))
        if all_exact_rows
        else None,
        "maximum_relative_gap": float(np.max([row["relative_gap"] for row in all_exact_rows]))
        if all_exact_rows
        else None,
    }
    summary = {
        "status": "SMOKE_TEST" if args.smoke_test else "COMPLETE",
        "datasets": dataset_manifests,
        "compiler": {
            "case_policy_condition_checks": compiler_checks,
            "all_checks_passed": True,
            "toy_boundary": toy_check,
        },
        "statistics": statistics,
        "exact_restricted_diagnostic": exact_summary,
        "claim_boundary": {
            "klaeger": "retrospective within-panel Kd liability discovery",
            "pkis2": "retrospective first strong profile-activity discovery; not off-target annotated",
            "not_supported": [
                "cellular target engagement",
                "toxicity",
                "clinical safety",
                "efficacy",
                "therapeutic recommendation",
            ],
        },
        "runtime_seconds": time.time() - started,
    }

    write_csv(args.output_dir / "case_metrics.csv.gz", all_case_rows)
    write_csv(args.output_dir / "query_traces.csv.gz", all_traces)
    write_csv(args.output_dir / "method_summary.csv", method_summary)
    write_csv(args.output_dir / "curves.csv", curves)
    write_csv(args.output_dir / "support_abstention.csv", support_rows)
    write_csv(args.output_dir / "exact_restricted_optima.csv", all_exact_rows)
    write_csv(args.output_dir / "batch_schedules.csv.gz", all_batch_rows)
    for benchmark in benchmarks:
        write_csv(
            args.output_dir / f"{benchmark.name.lower()}_target_mapping.csv",
            [
                {"target_key": key, "target_name": name}
                for key, name in zip(benchmark.target_keys, benchmark.target_names)
            ],
        )
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "rdkit": rdBase.rdkitVersion,
        },
        "parameters": {
            "max_budget": args.max_budget,
            "bootstrap": args.bootstrap,
            "permutations": args.permutations,
            "chemical_weight_beta": 8.0,
            "fingerprint": "Morgan radius 2, 2048 bits",
            "support_thresholds": SUPPORT_THRESHOLDS,
            "lambda_cost_ratios": LAMBDA_COST_RATIOS,
        },
    }
    (args.output_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    make_figures(curves, statistics, args.output_dir)
    output_hashes = {
        str(path.relative_to(args.output_dir)): sha256(path)
        for path in sorted(args.output_dir.rglob("*"))
        if path.is_file() and path.name != "output_manifest.json"
    }
    (args.output_dir / "output_manifest.json").write_text(
        json.dumps({"files": output_hashes}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
