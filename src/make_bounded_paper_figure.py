#!/usr/bin/env python3
"""Make the paper's cross-panel effect and tail-risk figure from frozen JSON."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEV = json.loads((ROOT / "risk_controlled_coverage_output" / "summary.json").read_text())
EXT = json.loads((ROOT / "pkis1_external_validation_output" / "summary.json").read_text())

ROWS = [
    ("Klaeger", "primary", "klaeger_primary_standard"),
    ("Klaeger", "delta-0.5", "klaeger_delta_0_5"),
    ("Klaeger", "scaffold", "klaeger_scaffold_exclusion"),
    ("Klaeger", "remove-10", "klaeger_leave10_nearest"),
    ("PKIS2", ">90 standard", "pkis2_gt90_standard"),
    ("PKIS2", ">80 standard", "pkis2_gt80_sensitivity"),
    ("PKIS2", ">90 scaffold", "pkis2_gt90_scaffold"),
    ("PKIS2", ">90 chemotype", "pkis2_gt90_chemotype"),
    ("PKIS1 external", ">90 standard", "pkis1_gt90_standard"),
    ("PKIS1 external", ">80 standard", "pkis1_gt80_sensitivity"),
    ("PKIS1 external", ">90 scaffold", "pkis1_gt90_scaffold"),
    ("PKIS1 external", ">90 remove-10", "pkis1_gt90_leave10_nearest"),
]


def values(dataset: str, condition: str) -> tuple[dict, dict]:
    if dataset == "PKIS1 external":
        result = EXT["comparisons"][condition]
        return result["bounded_minus_marginal"], result["unbounded_minus_marginal"]
    result = DEV["comparisons"][condition]
    return (
        result["bounded_coverage_minus_original_marginal"],
        result["original_coverage_minus_original_marginal"],
    )


def main() -> None:
    figure, (effect_axis, risk_axis) = plt.subplots(
        1,
        2,
        figsize=(11.2, 6.2),
        gridspec_kw={"width_ratios": [2.2, 1]},
        constrained_layout=True,
    )
    y = np.arange(len(ROWS))[::-1]
    bounded_values = [values(dataset, condition)[0] for dataset, _, condition in ROWS]
    unbounded_values = [values(dataset, condition)[1] for dataset, _, condition in ROWS]
    for offset, entries, label, color in [
        (0.12, unbounded_values, "unbounded coverage", "#999999"),
        (-0.12, bounded_values, "bounded coverage", "#136f63"),
    ]:
        estimate = np.asarray([entry["estimate"] for entry in entries])
        lower = estimate - np.asarray([entry["ci95"][0] for entry in entries])
        upper = np.asarray([entry["ci95"][1] for entry in entries]) - estimate
        effect_axis.errorbar(
            estimate,
            y + offset,
            xerr=np.vstack([lower, upper]),
            fmt="o",
            capsize=2.5,
            markersize=4.5,
            color=color,
            label=label,
        )
    effect_axis.axvline(0.0, color="#222222", linewidth=1)
    effect_axis.set_yticks(
        y,
        [f"{dataset}: {label}" for dataset, label, _ in ROWS],
        fontsize=8,
    )
    effect_axis.set_xlabel("AUDC difference versus weighted marginal")
    effect_axis.set_title("A  Mean first-hit performance")
    effect_axis.legend(frameon=False, fontsize=8, loc="lower right")
    effect_axis.axhspan(-0.5, 3.5, color="#f5e7e3", alpha=0.55, zorder=-5)

    unbounded_delays = np.asarray([entry["large_delay_ge_10"] for entry in unbounded_values])
    bounded_delays = np.asarray([entry["large_delay_ge_10"] for entry in bounded_values])
    risk_axis.barh(y + 0.16, unbounded_delays, height=0.3, color="#999999")
    risk_axis.barh(y - 0.16, bounded_delays, height=0.3, color="#136f63")
    risk_axis.set_yticks([])
    risk_axis.set_xlabel(r"case-conditions delayed by $\geq$10 assays")
    risk_axis.set_title("B  Tail risk")
    risk_axis.set_xlim(0, max(40, int(unbounded_delays.max()) + 3))
    for yi, count in zip(y, unbounded_delays):
        if count:
            risk_axis.text(count + 0.6, yi + 0.16, str(int(count)), va="center", fontsize=8)
    risk_axis.text(
        0.5,
        -0.75,
        f"bounded: {int(bounded_delays.sum())} in all {sum(entry['n'] for entry in bounded_values):,} evaluations",
        fontsize=8,
        color="#136f63",
    )
    output = ROOT / "paper" / "figures" / "bounded_cross_panel.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220)
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
