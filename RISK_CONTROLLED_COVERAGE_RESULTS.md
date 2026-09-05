# Displacement-bounded coverage results

Date completed: 2026-09-05

Headline status: **strong prespecified component result; five of six absolute
criteria pass; reference-resampling stability does not improve**

## Main result

`BOUNDED_COVERAGE` greedily diversifies kinase assays while constraining every
target to remain within nine positions of its chemical weighted-marginal rank.
The constraint gives a deterministic first-hit guarantee: relative to marginal
ranking, no recorded first activity can be delayed by ten or more assays.

The full evaluation verified zero at-least-ten-assay delays in all eight
conditions and all 3,004 held-out compound-condition cases.

| Condition | Bounded coverage minus marginal AUDC | 95% bootstrap interval |
|---|---:|---:|
| Klaeger primary | +0.01306 | [0.00045, 0.02568] |
| Klaeger delta-0.5 | +0.01802 | [-0.00180, 0.03829] |
| Klaeger scaffold exclusion | +0.01351 | [0.00090, 0.02658] |
| Klaeger remove-10-nearest | +0.04009 | [0.02117, 0.06036] |
| PKIS2 >90 standard | +0.00680 | [-0.00133, 0.01469] |
| PKIS2 >80 sensitivity | +0.00062 | [-0.00750, 0.00836] |
| PKIS2 >90 scaffold exclusion | +0.01266 | [0.00422, 0.02125] |
| PKIS2 >90 chemotype exclusion | +0.04328 | [0.03211, 0.05445] |

Every point estimate is nonnegative. On PKIS2 chemotype exclusion, bounded
coverage retains 80.64% of the original unconstrained gain (`0.04328` versus
`0.05367`) while changing the observed large-delay count from 25 to zero.

## Direct comparison with unconstrained coverage

Bounded coverage fixes the regimes where unconstrained coverage failed:

- PKIS2 standard: `+0.02141` AUDC, interval `[0.00813, 0.03477]`, nominal
  `p=0.00171`, Holm-adjusted across eight conditions `p=0.01197`.
- PKIS2 greater-than-80: `+0.02742`, interval `[0.01430, 0.04117]`, nominal
  `p=0.00007`, Holm-adjusted `p=0.00056`.
- PKIS2 scaffold exclusion: `+0.01516`, interval `[0.00055, 0.03000]`, nominal
  `p=0.04519`; the eight-condition Holm adjustment is not significant.
- PKIS2 chemotype exclusion: `-0.01039`, interval
  `[-0.02555, 0.00477]`, nominal `p=0.18487`. Thus the method trades a modest,
  statistically unresolved amount of unconstrained chemotype gain for the
  deterministic delay bound.

At the individual assay budgets, the detailed audit reports hit-rate effects,
paired bootstrap intervals, nominal sign-flip tests, and Holm adjustments for
budgets 1, 3, 5, 10, and 20.

## Frozen component audit

The bounded component passes five of six absolute criteria:

| Criterion | Result |
|---|---:|
| Chemotype retains at least 75%, CI lower > 0 | PASS |
| PKIS2 standard nonnegative and noninferior | PASS |
| PKIS2 >80 nonnegative and noninferior | PASS |
| Klaeger primary noninferior | PASS |
| Zero large delays in every condition | PASS |
| More stable than original coverage | FAIL |

Mean top-10 Jaccard under 50 fixed 80% reference subsamples was `0.76328` for
bounded coverage and `0.78835` for unconstrained coverage. The paired difference
was `-0.02507`, interval `[-0.03033, -0.01978]`. We therefore do not claim that
the bounded method improves reference-resampling stability.

## Role of the conformal wrapper

The separately frozen `RISK_CONTROLLED_COVERAGE` wrapper applied bounded
coverage only to label-free conformal chemical-support outliers. It passed the
standard, greater-than-80, Klaeger, and zero-delay criteria, but retained only
60.12% of the original chemotype gain and remained less stable than
unconstrained coverage. Its composite failed. It is an ablation, not the
headline method.

`BOUNDED_COVERAGE` was itself explicitly prespecified as M2 before the final
run, but choosing it as the headline over the failed M3 wrapper occurred after
observing this experiment. The paper must disclose that model-selection step
and should seek a fresh external validation dataset before making a
confirmatory generalization claim.

## Scientific interpretation

The evidence supports a narrow, decision-relevant claim: unrestricted profile
diversification can improve chemotype-shifted first-liability discovery but
occasionally pushes high-marginal assays dangerously far down the queue.
Marginal-rank deadlines preserve most of the diversity benefit and remove that
specific tail risk by construction. They do not make the assay order more
stable to reference resampling.

All results concern retrospective discovery of recorded within-panel binding
liabilities. They do not establish cellular engagement, toxicity, clinical
safety, efficacy, or therapeutic suitability.
