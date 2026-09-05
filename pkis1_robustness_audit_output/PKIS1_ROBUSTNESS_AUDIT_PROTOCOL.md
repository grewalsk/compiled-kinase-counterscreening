# PKIS1 threshold and target-resolution robustness audit

Version: 1.0, frozen after the external validation and before this audit

Date frozen: 2026-09-05

## Status and purpose

The frozen external PKIS1 validation failed. In particular,
`BOUNDED_COVERAGE - MARGINAL` was negative at the prespecified `>80` and
`>90` thresholds, although the formal no-large-delay guarantee held.

This is an explicitly exploratory robustness audit. It asks whether the
negative transfer is an artifact of either (i) the chosen activity threshold
or (ii) collapsing 224 assay constructs to 200 parent kinase targets. It does
not modify the method, rescue the failed validation, select a new threshold,
or create a new confirmatory claim. All conditions will be reported.

## Frozen data and preprocessing

Use exactly the raw input, source hash, compound standardization, eight-
structure development-overlap exclusion, missing-data rule, repeat averaging,
chemical representation, weights, tie breaks, budget 20, and three policies
from `PKIS1_EXTERNAL_VALIDATION_PROTOCOL.md` plus its two pre-outcome
amendments.

Construct two outcome matrices:

1. **Parent target:** repeat-average each compound--assay pair and take the
   maximum across assays sharing a `TARGET_CHEMBL_ID`, producing 200 targets.
2. **Assay construct:** retain the repeat-averaged 224
   `ASSAY_CHEMBL_ID` columns without parent collapse.

For multiple source compound IDs that standardize to one retained parent
structure, take the per-column median. Preserve missing values.

## Frozen audit grid

For the 200-parent-target matrix, evaluate strict percent-inhibition
thresholds `>50`, `>70`, `>80`, and `>90` under:

- standard references;
- exact Bemis--Murcko scaffold exclusion;
- removal of the ten nearest eligible references.

For the 224-assay-construct matrix, evaluate `>80` and `>90` under the same
three reference conditions.

In every one of the 18 cells compare fixed `BOUNDED_COVERAGE` (`K=9`) and
fixed `UNBOUNDED_COVERAGE` with `MARGINAL`. Report AUDC over budgets 1--20,
paired 95% intervals from 10,000 compound bootstraps, two-sided p-values from
100,000 paired sign flips, hit rates at budgets 1, 3, 5, 10, and 20, activity
prevalence, and large-delay counts. Use master seed 20260905.

## Interpretation rule

- If bounded-minus-marginal stays nonpositive across most of the grid, the
  negative external result is robust to threshold and target resolution.
- If only selected cells turn positive, they are hypothesis-generating regime
  indicators, not evidence that the failed primary criterion was misspecified.
- If the construct-level sign reverses consistently at both thresholds, parent
  collapse is a plausible explanatory artifact requiring a new independent
  dataset.
- The frozen external-validation verdict remains unchanged under every audit
  outcome.
