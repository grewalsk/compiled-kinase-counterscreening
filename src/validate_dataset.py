#!/usr/bin/env python3
"""Independently validate the compact counter-screen benchmark.

This validator deliberately does not import the extraction code.  It checks
the derived artifact hashes, relational invariants, label arithmetic,
censor-safety condition, mechanism exclusions, and provenance coverage.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


COMPOUND_RE = re.compile(r"^C\d{4}$")
TARGET_RE = re.compile(r"^T\d{4}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    data_dir = args.data_dir
    manifest = json.loads((data_dir / "manifest.json").read_text(encoding="utf-8"))
    for name, metadata in manifest["outputs"].items():
        path = data_dir / name
        require(path.is_file(), f"Missing manifest output: {name}")
        require(path.stat().st_size == metadata["bytes"], f"Size mismatch: {name}")
        require(sha256(path) == metadata["sha256"], f"SHA-256 mismatch: {name}")

    cases = [json.loads(line) for line in (data_dir / "cases.jsonl").read_text(
        encoding="utf-8"
    ).splitlines() if line]
    gold = read_csv(data_dir / "gold.csv.gz")
    compounds = read_csv(data_dir / "compound_mapping.csv")
    targets = read_csv(data_dir / "target_mapping.csv")
    mechanisms = read_csv(data_dir / "mechanism_provenance.csv")
    provenance = read_csv(data_dir / "activity_provenance.csv.gz")

    case_by_key = {case["compound_key"]: case for case in cases}
    require(len(case_by_key) == len(cases), "Duplicate case compound keys")
    require(all(COMPOUND_RE.fullmatch(key) for key in case_by_key), "Bad compound key")
    require(len(compounds) == len(cases), "Compound mapping/case count mismatch")
    require({row["compound_key"] for row in compounds} == set(case_by_key),
            "Compound mapping keys do not equal case keys")

    target_by_key = {row["target_key"]: row for row in targets}
    require(len(target_by_key) == len(targets), "Duplicate target mapping keys")
    require(all(TARGET_RE.fullmatch(key) for key in target_by_key), "Bad target key")
    common_targets = {
        row["target_key"] for row in targets if int(row["in_common_panel"]) == 1
    }
    require(len(common_targets) == manifest["counts"]["common_targets"],
            "Common target count mismatch")

    mechanism_pairs = {(row["compound_key"], row["target_key"]) for row in mechanisms}
    require(len(mechanism_pairs) > 0, "No mechanism provenance")
    require(all(c in case_by_key and t in target_by_key for c, t in mechanism_pairs),
            "Unknown key in mechanism provenance")

    gold_pairs = [(row["compound_key"], row["target_key"]) for row in gold]
    require(len(set(gold_pairs)) == len(gold_pairs), "Duplicate gold compound-target cell")
    require(len(gold) == manifest["counts"]["gold_candidate_rows"],
            "Gold row count mismatch")
    gold_by_compound: dict[str, set[str]] = defaultdict(set)
    primary_counts: Counter[str] = Counter()
    sensitivity_counts: Counter[str] = Counter()
    activity_ids: list[int] = []

    forbidden_gold_fields = {
        "molecule_chembl_id", "target_chembl_id", "target_name", "assay_chembl_id"
    }
    require(not (set(gold[0]) & forbidden_gold_fields), "Identifiers leaked into gold")

    for row in gold:
        compound = row["compound_key"]
        target = row["target_key"]
        require(compound in case_by_key and target in target_by_key, "Unknown gold key")
        require((compound, target) not in mechanism_pairs,
                f"Mechanism target leaked into candidate gold: {compound}, {target}")
        require(target in common_targets, f"Non-common candidate target: {target}")
        relation = row["standard_relation"]
        require(relation in {"=", ">"}, f"Unexpected relation: {relation}")
        kd_nm = float(row["kd_nm"])
        theta_nm = float(case_by_key[compound]["theta_nm"])
        require(0 < theta_nm <= 1000.0, f"Bad theta: {compound}")
        expected_primary = int(relation == "=" and kd_nm <= 10.0 * theta_nm)
        expected_sensitivity = int(
            relation == "=" and kd_nm <= (10.0 ** 0.5) * theta_nm
        )
        require(int(row["liability_delta_1"]) == expected_primary,
                f"Primary label mismatch: {compound}, {target}")
        require(int(row["liability_delta_0_5"]) == expected_sensitivity,
                f"Sensitivity label mismatch: {compound}, {target}")
        if relation == ">":
            require(kd_nm > 10.0 * theta_nm,
                    f"Ambiguous primary censor: {compound}, {target}")
        gold_by_compound[compound].add(target)
        primary_counts[compound] += expected_primary
        sensitivity_counts[compound] += expected_sensitivity
        activity_ids.append(int(row["activity_id"]))

    require(len(set(activity_ids)) == len(activity_ids), "Reused activity ID in gold")

    for compound, case in case_by_key.items():
        documented = set(case["documented_target_keys"])
        anchors = set(case["anchor_target_keys"])
        candidates = set(case["candidate_target_keys"])
        require(bool(anchors), f"No anchor target: {compound}")
        require(anchors <= documented, f"Anchor outside documented set: {compound}")
        require(candidates == common_targets - documented,
                f"Candidate universe mismatch: {compound}")
        require(candidates == gold_by_compound[compound],
                f"Case/gold candidate mismatch: {compound}")
        require(all((compound, target) in mechanism_pairs for target in documented),
                f"Missing mechanism provenance: {compound}")
        require(isinstance(case["canonical_smiles"], str) and case["canonical_smiles"],
                f"Missing canonical SMILES: {compound}")

    provenance_pairs = [(row["compound_key"], row["target_key"]) for row in provenance]
    require(len(provenance_pairs) == len(gold_pairs), "Activity provenance row mismatch")
    require(set(provenance_pairs) == set(gold_pairs), "Activity provenance coverage mismatch")
    require(len({int(row["activity_id"]) for row in provenance}) == len(provenance),
            "Duplicate provenance activity ID")
    require({int(row["activity_id"]) for row in provenance} == set(activity_ids),
            "Gold/provenance activity IDs differ")
    require({row["document_chembl_id"] for row in provenance} == {"CHEMBL3991601"},
            "Unexpected source document")

    counts = {
        "eligible_compounds": len(cases),
        "common_targets": len(common_targets),
        "gold_candidate_rows": len(gold),
        "primary_liability_positive_compounds": sum(
            primary_counts[key] > 0 for key in case_by_key
        ),
        "sensitivity_liability_positive_compounds": sum(
            sensitivity_counts[key] > 0 for key in case_by_key
        ),
    }
    for key, value in counts.items():
        require(value == manifest["counts"][key], f"Manifest count mismatch: {key}")

    report = {
        "status": "PASS",
        "checks": {
            "artifact_hashes": "PASS",
            "opaque_key_integrity": "PASS",
            "complete_common_panel": "PASS",
            "mechanism_exclusion": "PASS",
            "label_recalculation": "PASS",
            "primary_censor_safety": "PASS",
            "activity_provenance": "PASS",
        },
        "counts": counts,
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
