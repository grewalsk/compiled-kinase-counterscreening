# Exploratory invariant-coverage modeling protocol

Version: 1.0, frozen before evaluating the new models

## Status

This is a prospective protocol for an exploratory model derived after the
completed PKIS2 mechanism audit. It does not revise the failed confirmatory
Klaeger decision and must not be described as preregistered confirmation. The
model is motivated by the frozen finding that targetwise shuffling removed the
PKIS2 chemotype gain while whole-profile shuffling inside similarity strata
increased it.

## Modeling question

Can a deterministic model retain coarse chemical locality while removing the
unsupported fine-grained assignment of individual activity profiles to
chemical-similarity weights, and can a reference-only gate reduce harmful
coverage decisions?

## Frozen models

For every held-out compound, let `R` be the reference set produced by the
existing exclusion, `s_i` its Morgan-fingerprint Tanimoto similarities, and
`w_i` the existing normalized weights proportional to `exp(8 s_i)`.

1. **M0 ORIGINAL_MARGINAL:** original weights and marginal activity ranking.
2. **M1 ORIGINAL_COVERAGE:** original weights and greedy profile coverage.
3. **M2 INVARIANT_MARGINAL:** group-averaged weights and marginal ranking.
4. **M3 INVARIANT_COVERAGE:** group-averaged weights and greedy profile
   coverage. This is the primary new model.
5. **M4 SELECTIVE_INVARIANT_COVERAGE:** use M3 only when a reference-only local
   pseudo-holdout estimate is positive; otherwise abstain to M2. This is the
   secondary safety model.

No outcome-dependent tuning is permitted.

## Group-averaged weights

Sort the reference similarities and split them into five deterministic,
equally sized strata, exactly as in the completed mechanism audit. Preserve
the original total weight of each stratum and distribute it uniformly inside
the stratum:

`M_b = sum_{i in B_b} w_i` and `wbar_i = M_b / |B_b|` for `i in B_b`.

This is the Reynolds average over the product of within-stratum symmetric
groups. It is invariant to any permutation of profiles inside a stratum. Both
M2 and M3 use the same averaged weights.

## Reference-only selective gate

For held-out compound `h`, rank its eligible reference compounds by chemical
similarity. Visit that ranking until 32 valid pseudo-held-out compounds are
obtained. For every pseudo-held-out compound `j`:

1. remove `j` completely;
2. intersect the outer reference pool with the references allowed by the same
   exclusion rule when `j` is held out;
3. require at least 5 remaining references and 20 measured candidate assays;
4. construct M2 and M3 using only the remaining reference data;
5. evaluate their AUDC difference on the already known profile of `j`.

Weight the pseudo-holdout differences by the existing similarity weighting
between `h` and the pseudo-held-out compounds. Use M3 only if at least eight
valid pseudo-holdouts exist and the weighted mean difference is strictly
positive. Otherwise use M2. The true activity values of `h` must never enter
the gate.

## Frozen benchmarks and conditions

Use the exact existing parsers, molecular representation, exclusions,
candidate masks, thresholds, tie-breaking, and 20-assay budget.

- Klaeger/ChEMBL 30: primary standard, delta-0.5 sensitivity, scaffold
  exclusion, and removal of ten nearest references.
- PKIS2: greater-than-90 standard, greater-than-80 sensitivity,
  greater-than-90 scaffold exclusion, and greater-than-90 source-chemotype
  exclusion.

PKIS2 source SHA-256:
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
Klaeger/ChEMBL derived records remain CC BY-SA 3.0; PKIS2 remains CC BY 4.0.
The two measurement types are never pooled.

## Outcomes and comparisons

For every condition report:

- M1 minus M0;
- M3 minus M2, the primary new within-model comparison;
- M3 minus M1, the direct model-revision comparison;
- M4 minus M2, the selective-policy comparison.

Report AUDC, censored cost to first recorded activity, budgets 1, 3, 5, 10,
and 20, 10,000 paired bootstraps, 100,000 paired sign flips, wins/losses/ties,
mean assays earlier, and counts at least 10 assays earlier or later.

## Stability diagnostic

For the PKIS2 greater-than-90 source-chemotype condition only, create 50 fixed
80-percent reference subsamples per held-out compound. Recompute M1 and M3 and
report mean top-10 set Jaccard agreement with each corresponding full-reference
order. Seeds are derived from NumPy `SeedSequence` with master seed 20260905.

## Frozen success criteria

The composite exploratory success flag requires all of the following:

1. M3 minus M2 on PKIS2 chemotype is at least 75 percent of the completed
   original gain: `0.04025390625`.
2. M3 minus M2 is greater than M1 minus M0 for both PKIS2 standard and the
   greater-than-80 sensitivity condition. Nonnegative values are separately
   reported as reversal-elimination flags.
3. M3 does not reduce Klaeger primary AUDC by more than `0.005` relative to M1.
4. On PKIS2 chemotype, M4 has no more than 20 cases delayed by at least 10
   assays relative to M2 and retains at least half the completed original gain:
   `0.0268359375`.
5. Mean top-10 stability of M3 is strictly greater than that of M1.

Failure of the composite does not authorize changing stratum count, gate
threshold, pseudo-holdout count, or success thresholds.

## Leakage, validity, and interpretation controls

- Held-out activity values are used only for final evaluation.
- Pseudo-held-out activity values are permitted because they belong to the
  reference pool, but the pseudo-held-out compound is removed from its own
  training pool.
- Assay missingness may define available candidates; activity magnitudes or
  labels may not define the gate before evaluation.
- Identical structures and the condition-specific exclusion are applied in
  every outer and inner split.
- Results from all eight conditions are reported, regardless of direction.
- This remains retrospective first-activity or binding-liability discovery.
  It does not establish cellular engagement, selectivity, toxicity, safety,
  efficacy, or therapeutic suitability.
