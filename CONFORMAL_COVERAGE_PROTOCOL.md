# Conformal support-gated coverage protocol

Version: 1.0, frozen before evaluating `CONFORMAL_COVERAGE` outcomes

Date frozen: 2026-09-05

## Status and motivation

This is a prospective specification of a second post hoc model-development
experiment. The previously frozen pseudo-outcome `SAFE_COVERAGE` selector
retained a significant PKIS2 chemotype gain and reduced large delays, but
failed its composite criterion because its compound-level uplift estimates
were weak and its decisions were unstable under reference resampling.

The new method tests one mechanistically narrower claim: profile coverage is
useful when the query compound has unusually weak chemical support in the
available reference library, whereas weighted marginal ranking is the safer
default for chemically supported queries. The gate is label-free and has one
fixed conformal level. It does not use the failed invariant model or tune a
threshold against activity outcomes.

This remains post hoc development on datasets whose earlier outcomes have been
analyzed; it is not independent confirmation. All eight conditions will be
reported even if the composite fails.

## Frozen policies

For held-out compound `h`, let `R_h` be the eligible references under the
condition-specific exclusion, and let

`w_i(h) proportional to exp(8 * Tanimoto(h, i))`.

1. **M0 ORIGINAL_MARGINAL:** the unchanged Beta(1,1)-smoothed weighted
   marginal ranking.
2. **M1 ORIGINAL_COVERAGE:** the unchanged greedy weighted profile-coverage
   order with M0 as its fixed tie break.
3. **M2 CONFORMAL_COVERAGE:** use M1 only when the label-free chemical-support
   test below declares `h` under-supported; otherwise return M0 exactly.

No target-activity label, candidate activity magnitude, assay outcome, or
pseudo-held-out utility enters the M2 gate.

## Label-free conformal support gate

Let the query support be its nearest eligible-reference similarity:

`q_h = max_{i in R_h} Tanimoto(h, i)`.

For every reference compound `i` in the fixed outer pool `R_h`, compute its
leave-one-reference-out support inside that same pool:

`q_i = max_{j in R_h, j != i} Tanimoto(i, j)`.

Define the deterministic lower-tail conformal rank

`p_support = (1 + count_i[q_i <= q_h]) / (|R_h| + 1)`.

Use M1 exactly when `p_support <= 0.10`; otherwise use M0. The `0.10` level is
fixed as a conventional 90% conformal outlier rule before M2 outcomes are
computed. Equality invokes coverage. Similarity and compound-index ties are
handled by the displayed non-randomized formula. Fewer than two references or
any non-finite support value forces abstention to M0.

The calibration pool is condition-specific but label-free. This tests whether
the query is unusually poorly represented relative to the chemical support
that the available references provide one another.

## Benchmarks and conditions

Use the existing molecular standardization, Morgan fingerprints, Tanimoto
similarity, parsers, labels, candidate masks, tie breaks, exclusions, and
20-assay budget without change.

- Klaeger/ChEMBL 30: primary standard, delta-0.5 sensitivity, scaffold
  exclusion, and removal of the ten nearest references.
- PKIS2: greater-than-90 standard, greater-than-80 sensitivity,
  greater-than-90 scaffold exclusion, and greater-than-90 source-chemotype
  exclusion.

PKIS2 source SHA-256:
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
Klaeger/ChEMBL-derived records remain CC BY-SA 3.0; PKIS2 remains CC BY 4.0.
The two measurement types are never pooled.

## Outcomes and statistics

For each condition, report M1-M0, M2-M0, and M2-M1 for:

- AUDC for first recorded activity over assay budgets 1--20;
- censored cost to first recorded activity;
- hit at budgets 1, 3, 5, 10, and 20;
- wins, losses, ties, mean assays earlier, and counts at least ten assays
  earlier or later.

Use 10,000 paired compound bootstraps and 100,000 paired sign flips with master
seed 20260905. Report gate rate, support-p-value distribution, query maximum
similarity, selected and abstained realized M1-M0 differences, and all
compound-level orders. Realized differences are diagnostic only and cannot
alter the gate.

For PKIS2 chemotype exclusion, create 50 fixed 80% reference subsamples per
held-out compound. Recompute M1 and the complete label-free M2 gate. Report
mean top-10 Jaccard agreement with the corresponding full-reference order,
paired M2-M1 stability difference and interval, and gate-decision agreement.

## Frozen absolute success criteria

M2 is successful only if every criterion passes:

1. PKIS2 chemotype M2-M0 is at least `0.04025390625` (75% of the completed
   M1-M0 gain `0.053671875`) and its two-sided 95% bootstrap lower endpoint is
   greater than zero.
2. PKIS2 standard M2-M0 has point estimate at least zero and 95% lower endpoint
   at least `-0.015`.
3. PKIS2 greater-than-80 M2-M0 has point estimate at least zero and 95% lower
   endpoint at least `-0.015`.
4. Klaeger primary M2-M0 has 95% lower endpoint at least `-0.015`.
5. PKIS2 chemotype M2 has no more than 20 cases delayed by at least ten assays
   relative to M0.
6. PKIS2 chemotype mean M2 top-10 stability is at least M1 stability, and the
   lower endpoint of the paired 95% interval for M2-M1 is at least `-0.01`.

Failure does not authorize changing `0.10`, the conformal direction, support
score, policies, exclusions, success thresholds, or stability design.

## Leakage and scientific scope

- The gate uses structures and eligible-reference membership only.
- Held-out activity outcomes enter only final evaluation.
- Identical-structure and condition-specific exclusions are enforced before
  support calibration.
- Missing assays are never negatives and do not enter the support gate.
- The study evaluates retrospective discovery of recorded within-panel binding
  liabilities. It does not establish cellular engagement, selectivity,
  toxicity, clinical safety, efficacy, or therapeutic suitability.
