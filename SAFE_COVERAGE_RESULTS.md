# Cross-fitted safe-coverage results

Date completed: 2026-09-05

Status: **frozen composite criterion failed; useful risk--reward trade-off**

## Result

The outcome-only pseudo-holdout selector passed three of six frozen absolute
criteria. It should not be presented as the final successful method.

On PKIS2 chemotype exclusion, `SAFE_COVERAGE` improved AUDC over original
marginal ranking by `0.03078125` (95% bootstrap interval
`[0.01695117, 0.04515625]`, paired sign-flip `p=0.000020`). It retained 57.35%
of the completed original-coverage gain rather than the required 75%.

The selector nevertheless repaired most of the operational failures:

- PKIS2 greater-than-80 changed from `-0.02679688` for original coverage to
  `+0.00156250` for safe coverage versus marginal ranking.
- PKIS2 standard changed from `-0.01460938` to `-0.00187500`; the interval
  `[-0.00750000, 0.00289062]` met the noninferiority margin, but the frozen
  nonnegative point-estimate rule failed.
- Chemotype cases delayed by at least ten assays versus marginal ranking fell
  from 25 to 10.
- Klaeger primary was effectively neutral at `-0.00045045`, with interval
  `[-0.00270270, 0.00135135]`.

The selector's stability criterion failed. Mean top-10 Jaccard under 80%
reference resampling was `0.7341`, compared with `0.7856` for original
coverage; the paired difference was `-0.05147` with interval
`[-0.06546, -0.03785]`. Gate decisions agreed in 81.0% of subsamples.

## Mechanistic diagnostic

The reference-only lower-confidence gate selected coverage for 34.69% of
PKIS2 chemotype cases, but only 5.63% of standard and 2.81% of greater-than-80
cases. Thus it detected the evaluation regime and sharply reduced tail risk.
Its compound-level uplift prediction was weak, however: predicted versus
realized coverage uplift had Pearson correlation `0.068` on chemotype cases.
The unstable, low-resolution casewise outcome gate explains both the lost
chemotype gain and the reduced reference-resampling stability.

## Frozen decisions

| Criterion | Result |
|---|---:|
| Chemotype retention at least 75%, CI lower > 0 | FAIL |
| PKIS2 standard nonnegative and noninferior | FAIL |
| PKIS2 >80 nonnegative and noninferior | PASS |
| Klaeger primary noninferior | PASS |
| Chemotype large delays at most 20 | PASS |
| Safe order stable | FAIL |

## Interpretation

This experiment supports abstention as a tail-risk control, but does not
support the pseudo-outcome gate as the final paper method. The next frozen
model replaces the noisy outcome-derived gate with a label-free conformal
chemical-support test. That follow-up is motivated by the observed regime
separation and must be reported separately from this failed test.

All statements concern retrospective discovery of recorded within-panel
binding liabilities, not cellular engagement, clinical safety, toxicity,
efficacy, or therapeutic suitability.
