# Cross-fitted safe-coverage protocol

Version: 1.1, frozen before evaluating `SAFE_COVERAGE` outcomes

Date frozen: 2026-09-05

## Status and question

This is a prospective, post hoc model-development experiment motivated by the
completed invariant-coverage audit. It is not a revision of the earlier frozen
confirmatory study. The question is whether a reference-only selector can use
profile coverage only when closely related pseudo-held-out reference compounds
support it, and otherwise abstain to the original weighted-marginal order.

The previous hard-invariance endpoint is not an input to this method. Its
absolute audit showed that interpolation toward that endpoint cannot repair the
observed PKIS2 reversals. No result from `SAFE_COVERAGE` may be inspected before
this protocol and its implementation freeze are written.

## Frozen policies

For held-out compound `h`, let `R_h` be the eligible reference set under the
condition-specific exclusion and let the original chemical weights be

`w_i(h) proportional to exp(8 * Tanimoto(h, i))`, for `i` in `R_h`.

All weights, fingerprints, labels, candidates, tie breaks, and the 20-assay
budget are inherited unchanged from `src/run_compiled_coverage.py`.

1. **M0 ORIGINAL_MARGINAL:** rank targets by the Beta(1,1)-smoothed weighted
   marginal probability under `w(h)`.
2. **M1 ORIGINAL_COVERAGE:** greedily maximize newly covered reference-profile
   mass under `w(h)`, with M0 as the fixed tie-break order.
3. **M2 SAFE_COVERAGE:** select M1 only when the cross-fitted lower confidence
   bound defined below is strictly positive; otherwise return M0 exactly.

M2 is a selector over two frozen policies, not a third scoring heuristic.

## Cross-fitted local evidence

Rank the compounds in `R_h` by decreasing similarity to `h`, breaking ties by
their existing integer index. Visit this order until 32 valid pseudo-held-out
compounds have been obtained or the list is exhausted. For pseudo-held-out
compound `j`:

1. remove `j` completely;
2. recompute the references allowed when `j` is held out under the same outer
   exclusion rule;
3. intersect those references with `R_h` so that the gate cannot access a
   compound excluded from the real held-out problem;
4. require at least five inner references and at least 20 measured candidate
   assays;
5. construct M0 and M1 for `j` using only the inner references;
6. record `d_j = AUDC(M1,j) - AUDC(M0,j)` from `j`'s known reference profile.

The true activity values of `h` never enter model selection.

Weight the valid pseudo-holdouts using the same original chemical weighting
between `h` and each `j`, normalized over the selected pseudo-holdouts. Let

`mu = sum_j v_j d_j`, `n_eff = 1 / sum_j v_j^2`, and

`s2 = sum_j v_j (d_j - mu)^2 / (1 - sum_j v_j^2)`.

For `n_eff > 1`, define `se = sqrt(s2 / n_eff)` and the fixed one-sided 90%
normal lower bound

`LCB90 = mu - 1.2815515655446004 * se`.

Use M1 only if there are at least eight valid pseudo-holdouts and `LCB90 > 0`.
Otherwise use M0. Zero variance with positive `mu` yields `LCB90 = mu`; invalid
or non-finite uncertainty yields abstention to M0. The threshold, confidence
level, pseudo-holdout count, and minimum count may not change after outcomes
are inspected.

## Benchmarks and conditions

Use the exact existing parsers and all eight existing conditions:

- Klaeger/ChEMBL 30: primary standard, delta-0.5 sensitivity, scaffold
  exclusion, and removal of the ten nearest references.
- PKIS2: greater-than-90 standard, greater-than-80 sensitivity,
  greater-than-90 scaffold exclusion, and greater-than-90 source-chemotype
  exclusion.

PKIS2 source SHA-256:
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
Klaeger/ChEMBL-derived records remain CC BY-SA 3.0; PKIS2 remains CC BY 4.0.
The measurement types are never pooled.

## Outcomes and comparisons

For every condition, report M1-M0, M2-M0, and M2-M1 for AUDC over budgets
1--20, censored cost to first recorded activity, hit at budgets 1, 3, 5, 10,
and 20, wins/losses/ties, mean assays earlier, and counts at least ten assays
earlier or later. Use 10,000 paired compound bootstraps and 100,000 paired sign
flips with master seed 20260905.

Report gate coverage, effective pseudo-holdout count, `mu`, `LCB90`, weighted
standard error, and the held-out realized M1-M0 difference for diagnostic
calibration. These diagnostics do not alter decisions.

For PKIS2 chemotype only, create 10 fixed 80% reference subsamples per held-out
compound. Recompute the complete selector, including its inner gate, and report
top-10 Jaccard agreement between each subsample order and the corresponding
full-reference M1 and M2 orders.

Version 1.1 changes only this diagnostic from 50 to 10 subsamples before any
implementation or outcome was run. A fully nested selector requires up to 32
inner policy comparisons per subsample; ten replicates across 640 held-out
PKIS2 compounds still provide 6,400 paired stability observations while
keeping the CPU-only experiment bounded. No model decision or success
threshold changed.

## Frozen absolute success criteria

`SAFE_COVERAGE` is a successful model only if all conditions below hold:

1. On PKIS2 chemotype exclusion, M2-M0 AUDC is at least `0.04025390625`
   (75% of the completed M1-M0 gain `0.053671875`) and its two-sided 95%
   bootstrap interval has lower endpoint greater than zero.
2. On PKIS2 standard, M2-M0 has point estimate at least zero and 95% lower
   endpoint at least `-0.015`.
3. On PKIS2 greater-than-80 sensitivity, M2-M0 has point estimate at least zero
   and 95% lower endpoint at least `-0.015`.
4. On Klaeger primary, M2-M0 has 95% lower endpoint at least `-0.015`.
5. On PKIS2 chemotype, M2 has no more than 20 cases delayed by at least ten
   assays relative to M0.
6. On PKIS2 chemotype, mean M2 top-10 reference-subsample stability is at least
   that of M1 and the lower endpoint of the paired 95% interval for M2-M1 is at
   least `-0.01`.

Passing a subset is not overall success. Failure does not authorize changing
the gate, confidence level, candidate policies, or thresholds. All eight
conditions must be reported regardless of direction.

## Leakage and interpretation controls

- Held-out outcomes may be used only for final evaluation and diagnostic
  calibration after the order is fixed.
- Every pseudo-held-out compound is removed from its own inner training pool.
- Outer and inner reference sets both enforce identical-structure and
  condition-specific exclusions.
- Missing assays are never negatives. Missingness may determine assay
  availability but not the gate outcome.
- The existing deterministic molecular standardization and Morgan fingerprints
  are unchanged.
- This retrospective task concerns discovery of recorded within-panel binding
  liabilities. It does not establish cellular engagement, clinical safety,
  toxicity, efficacy, or therapeutic suitability.
