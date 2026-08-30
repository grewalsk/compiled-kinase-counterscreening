# Frozen protocol: compiled profile-coverage counter-screening

Protocol date: 2026-08-30 (America/New_York)

Status: written after the completed active-search study was inspected, but
before the correlation-aware coverage method was run on either the Klaeger or
PKIS2 benchmark. The earlier results motivated this method and are not treated
as confirmatory evidence for it.

Working title: **When Sequential Is Not Adaptive: Compiled Profile-Coverage
Kinase Counter-Screening**

## 1. Claim and task

For a binary screen that stops at its first measured activity, compile a
sequential policy into its all-negative path and choose a compound-specific
target order that covers chemically weighted historical activity profiles
early. The study is retrospective assay prioritization. It does not establish
biochemical inhibition outside the source assay, cellular engagement,
toxicity, safety, efficacy, or therapeutic suitability.

## 2. Exact compilation check

For a deterministic policy `pi`, construct a fixed sequence by repeatedly
feeding `pi` a negative result. For every hidden binary profile, the compiled
sequence and the online policy must select the same prefix through the first
positive result or budget exhaustion. This implies identical first-hit time,
Hit@b, and any-hit AUDC. The implementation must assert equality for every
eligible case for both the existing particle-myopic and two-step rollout
policies. A failure is an implementation failure, not a statistical result.

The claim is restricted to stop-on-first-positive utility. A toy test must
show that branching can matter when screening continues after a positive or
rewards multiple hits.

## 3. Development data

Use the already frozen Klaeger/ChEMBL 30 compact benchmark without changing
its compounds, targets, mechanism exclusions, activity interpretations, or
provenance. The primary endpoint uses `liability_delta_1`; the frozen
`liability_delta_0_5` label is a sensitivity analysis. Missing or excluded
measurements remain missing and are never converted into negatives.

Every test compound is excluded from its reference set. Reference compounds
with an identical 2,048-bit radius-2 Morgan fingerprint are also excluded in
all conditions. Prespecified robustness conditions are:

1. remove all references sharing the test compound's Bemis--Murcko scaffold;
2. remove the ten most similar remaining reference compounds; and
3. use the frozen `delta=0.5` liability label.

Conditions with fewer than five references are unavailable rather than
imputed.

## 4. Compound-conditioned profile weights

For query compound `x` and eligible reference compound `j`, use

`w_j(x) proportional to exp(8 * Tanimoto(x,j))`.

Normalize weights to sum to one. The exponent 8, Morgan radius 2, and 2,048
bits are inherited unchanged from the completed particle method. They are not
tuned on the new method's results. Report nearest-reference similarity and
effective reference count `1/sum(w^2)` for every case.

## 5. Methods and baselines

All ties are resolved by opaque target identifier.

- `PREVALENCE`: Beta(1,1)-smoothed reference prevalence.
- `CHEMICAL_KNN`: the completed study's similarity-weighted ten-neighbour
  Beta-smoothed marginal ranker.
- `WEIGHTED_MARGINAL`: rank each target independently by its weighted
  Beta-smoothed activity probability under all eligible profiles.
- `GLOBAL_COVERAGE`: greedy profile coverage with uniform reference weights.
- `WEIGHTED_COVERAGE`: repeatedly select the target covering the largest
  remaining chemical-weight mass of not-yet-covered reference profiles. Ties
  use the initial weighted marginal score and then target identifier.
- `COMPILED_PARTICLE_MYOPIC` and `COMPILED_PARTICLE_ROLLOUT_2`: exact
  all-negative compilations of the existing online policies.
- exact hypergeometric random expectation and a hidden-outcome oracle are
  reported only as lower/upper references.

For unit assay costs, the weighted-coverage greedy sequence is both the
natural weighted min-sum set-cover greedy order and the cardinality maximum-
coverage greedy prefix. No new approximation theorem is claimed. On each
development case, solve an eight-target restriction exactly by subset dynamic
programming and report the greedy-to-optimal objective gap as a diagnostic;
this is not a global optimality claim.

## 6. Outcomes and statistics

Maximum budget is 20. First-hit time is censored to 21 when no hit is found.
Primary performance outcome is any-hit `AUDC_1_20`. The implementation must
verify exactly that `AUDC_1_20 = (21 - first_hit_time) / 20` for every fixed
sequence.

Primary contrast:

`WEIGHTED_COVERAGE - WEIGHTED_MARGINAL` on Klaeger `liability_delta_1`.

This isolates the contribution of correlated profile coverage from the same
chemical weights and marginal evidence. Secondary contrasts compare weighted
coverage with chemical-kNN, prevalence, global coverage, and compiled
particle policies. Report Hit@1/3/5/10/20 and censored assays to first hit.

Use 10,000 compound-paired bootstrap replicates for 95% intervals and 100,000
paired sign flips for two-sided p-values. The primary contrast is unadjusted;
Holm-adjust all prespecified secondary performance contrasts. Report every
effect and interval regardless of direction. Exact compiler invariants receive
no p-value.

Support analyses report performance and retained-case coverage over fixed
nearest-reference thresholds `{0.2, 0.3, 0.4, 0.5}`. These are model-support
abstentions, not biological clearance rules.

## 7. External PKIS2 transfer analysis

Download the Drewry et al. PLOS ONE S4 workbook from
`https://doi.org/10.1371/journal.pone.0181585.s004`. The expected SHA-256 is
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
The article and supplement are CC BY 4.0.

PKIS2 is analyzed separately from the Klaeger `Kd` benchmark. Drop the final
row whose four core compound identifiers are all empty; the implementation
must separately assert and record its otherwise stray `TEC = 13` cell rather
than silently treating it as an ordinary compound. Require the published
metadata columns; parse every SMILES with RDKit; keep the largest molecular
fragment without neutralizing it; and group rows with the same resulting
canonical parent SMILES.

Collapse the 406 assay columns to the 392 published parent kinase labels by
removing explicit phosphorylation, autoinhibited, cyclin-partner, and kinase-
domain qualifiers. Within each raw compound row and parent label, first take
the maximum percent inhibition across constructs. Then aggregate duplicate
parent-structure rows by their median collapsed value before thresholding. The
primary screen hit is strictly
greater than 90% inhibition at 1 micromolar; strictly greater than 80% is a
threshold sensitivity. These cutoffs correspond to the source study's stated
screening/follow-up rules.

Because the workbook does not provide a complete per-compound intended-target
annotation, this analysis is called **first strong profile-activity
discovery**, not off-target-liability discovery. No PKIS2 result may be called
a safety result. In addition to leave-one-parent-structure-out evaluation,
report scaffold exclusion and exclusion of the source-provided chemotype.

The external primary contrast remains `WEIGHTED_COVERAGE -
WEIGHTED_MARGINAL`; the method and beta are not retuned.

## 8. Batching

For each weighted-coverage order, estimate survival `q_i` as the chemical
weight of historical profiles with no covered hit in the first `i` positions.
With unit assay cost, solve the exact fixed-order batch partition dynamic
program for dimensionless delay/setup penalties
`lambda/c in {0, 0.25, 0.5, 1, 2, 5, 10, 20}`. Report expected modeled cost,
batch boundaries, and batch sizes. These are sensitivity scenarios, not
claimed laboratory prices.

## 9. Decision rule

The method is a strong replacement-paper result only if weighted coverage
improves over weighted marginal ranking in the primary Klaeger analysis and
the direction transfers to PKIS2, with no reversal under the strongest
chemical-series exclusions. An interval overlapping zero is reported as
inconclusive. A reversal or disappearance under exclusions is a negative
result. No dataset, threshold, fingerprint, or primary comparison is changed
after result access to rescue the method.

## 10. Required expert gates

Before submission, a kinase/medicinal-chemistry reviewer must approve the
Klaeger liability language, PKIS2 construct collapse and thresholds, intended-
target limitation, and operational interpretation of the proposed panels. An
optimization reviewer must check the compilation proposition, AUDC identity,
coverage objective, exact small-instance dynamic program, and batching
recurrence.
