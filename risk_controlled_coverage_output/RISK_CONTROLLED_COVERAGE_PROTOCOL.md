# Conformal displacement-bounded coverage protocol

Version: 1.0, frozen before evaluating `RISK_CONTROLLED_COVERAGE` outcomes

Date frozen: 2026-09-05

## Status and single hypothesis

This is the final prospective specification in the current post hoc
model-development sequence. The frozen conformal-support model retained 78.31%
of the PKIS2 chemotype gain and largely removed the standard and threshold
reversals, but missed its tail criterion by one compound and slightly reduced
order stability. Its unbounded coverage component could move targets
arbitrarily far from the marginal order.

The single new hypothesis is that a bounded-displacement coverage order can
preserve the chemotype benefit of profile diversity while guaranteeing that a
recorded first activity is never delayed by ten or more assays relative to
weighted marginal ranking. The conformal support gate, alpha, chemical weights,
datasets, and all statistical criteria remain unchanged. No risk-controlled
outcome may be inspected before this protocol and implementation are frozen.

This is post hoc method development on previously analyzed datasets, not an
independent confirmatory study. A failed composite will end this modeling
sequence rather than trigger further tuning.

## Frozen policies

For held-out compound `h`, construct the unchanged original chemical weights,
weighted marginal probabilities, and full deterministic marginal target order.

1. **M0 ORIGINAL_MARGINAL:** the first 20 targets in the weighted marginal
   order.
2. **M1 ORIGINAL_COVERAGE:** the existing unconstrained weighted profile-
   coverage policy.
3. **M2 BOUNDED_COVERAGE:** the displacement-bounded greedy coverage order
   below, using fixed bound `K=9`.
4. **M3 RISK_CONTROLLED_COVERAGE:** use M2 when the already frozen label-free
   conformal support p-value is at most `0.10`; otherwise return M0 exactly.

M3 is the only proposed final method. M1 and M2 are component ablations.

## Displacement-bounded coverage

Index the complete marginal order from zero, and let `r(a)` be target `a`'s
marginal rank. At output position `t`, among unselected candidates:

1. If a target has reached its deadline, `r(a) + K <= t`, select the one with
   smallest marginal rank.
2. Otherwise restrict eligibility to targets satisfying `r(a) <= t + K` and
   select the target with greatest current uncovered weighted profile mass.
3. Break equal coverage gains by larger frozen marginal probability and then
   opaque target key, exactly as in original coverage.
4. Remove newly covered reference profiles and continue through budget 20.

Set `K=9`. The early-eligibility and deadline rules imply

`|position_M2(a) - r(a)| <= 9`

for every target reached by its deadline. If the first positive under M0 has
rank at most 10, M2 must encounter it no more than nine assays later. If its
rank exceeds 10, censoring at assay 21 itself limits the observed delay to at
most nine. Therefore the evaluated censored cost to first activity for M2, and
for M3 whenever it invokes M2, cannot be ten or more assays worse than M0. The
implementation must assert this property for every evaluated case.

`K=9` is fixed directly from the previously declared large-delay threshold of
ten assays; it is not selected from an outcome sweep.

## Label-free conformal gate

Reuse the frozen conformal gate without modification. For eligible references
`R_h`, let `q_h` be the maximum query-to-reference Morgan Tanimoto similarity
and `q_i` each reference compound's maximum similarity to another compound in
`R_h`. Define

`p_support = (1 + count_i[q_i <= q_h]) / (|R_h| + 1)`.

Invoke M2 if and only if `p_support <= 0.10`; otherwise return M0. No activity
label or observed target outcome enters this decision.

## Data, exclusions, and conditions

Use all eight unchanged conditions:

- Klaeger/ChEMBL 30: primary standard, delta-0.5 sensitivity, scaffold
  exclusion, and removal of the ten nearest references.
- PKIS2: greater-than-90 standard, greater-than-80 sensitivity,
  greater-than-90 scaffold exclusion, and greater-than-90 source-chemotype
  exclusion.

The parsers, molecular standardization, Morgan fingerprints, Tanimoto
similarity, weight temperature 8, candidate masks, activity definitions,
missing-value handling, exclusions, target-key tie breaks, and budget 20 are
unchanged. PKIS2 source SHA-256 is
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
Klaeger/ChEMBL-derived records remain CC BY-SA 3.0 and PKIS2 remains CC BY 4.0.

## Outcomes and statistics

Report M1-M0, M2-M0, M3-M0, and M3-M1 for AUDC over budgets 1--20, censored
cost to first activity, hits at budgets 1, 3, 5, 10, and 20, paired wins,
losses and ties, mean assays earlier, and at-least-ten-assay advances or delays.
Use 10,000 paired compound bootstraps and 100,000 paired sign flips with master
seed 20260905.

Report conformal gate rates and support distributions. For PKIS2 chemotype,
perform 50 fixed 80% reference subsamples per held-out compound, recomputing
M1, M2, M3, and the conformal gate. Report top-10 Jaccard stability for all
three coverage policies, the paired M3-M1 difference and interval, and gate
decision agreement.

## Frozen absolute success criteria

M3 succeeds only if every condition passes:

1. PKIS2 chemotype M3-M0 is at least `0.04025390625` and its two-sided 95%
   bootstrap lower endpoint is greater than zero.
2. PKIS2 standard M3-M0 has point estimate at least zero and 95% lower endpoint
   at least `-0.015`.
3. PKIS2 greater-than-80 M3-M0 has point estimate at least zero and 95% lower
   endpoint at least `-0.015`.
4. Klaeger primary M3-M0 has 95% lower endpoint at least `-0.015`.
5. PKIS2 chemotype has zero M3 cases delayed by at least ten assays relative
   to M0. This is an implementation check of the stated guarantee.
6. PKIS2 chemotype mean M3 top-10 stability is at least M1 stability, and the
   lower endpoint of the paired 95% interval for M3-M1 is at least `-0.01`.

Failure does not authorize changing `K`, alpha, the gate, the success criteria,
or the greedy rule.

## Leakage and scope

- The gate uses structures and eligible-reference membership only.
- Held-out activities are used only after all four orders are fixed.
- Missing assays are never negatives.
- Outer identical-structure and condition-specific exclusions are applied
  before weights or support calibration.
- Results concern retrospective discovery of recorded within-panel binding
  liabilities, not cellular engagement, toxicity, clinical safety, efficacy,
  or therapeutic suitability.
