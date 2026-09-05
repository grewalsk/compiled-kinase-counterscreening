# Conformal support-gated coverage results

Date completed: 2026-09-05

Status: **near-success; frozen composite criterion failed**

## Primary result

The label-free conformal gate retained 78.31% of the original PKIS2 chemotype
gain. `CONFORMAL_COVERAGE` improved AUDC over original marginal ranking by
`0.04203125` (95% bootstrap interval `[0.02367188, 0.05984375]`, paired
sign-flip `p=0.000010`) compared with `0.05367188` for unbounded original
coverage.

The gate invoked coverage for 471/640 (73.59%) chemotype-exclusion cases but
only 64/640 (10.00%) standard cases, exactly matching the intended
chemical-support distinction.

## Absolute robustness

- PKIS2 standard changed from `-0.01460938` for original coverage to
  `+0.00109375` for conformal coverage versus marginal ranking. The direct
  improvement over original coverage was `+0.01570313`, interval
  `[0.00117188, 0.03015625]`, `p=0.03494`.
- PKIS2 greater-than-80 changed from `-0.02679688` to `-0.00390625`. The direct
  improvement over original coverage was `+0.02289062`, interval
  `[0.00773437, 0.03859375]`, `p=0.00366`. The frozen nonnegative point rule
  nevertheless failed.
- Klaeger primary was `+0.00540541` versus marginal ranking, interval
  `[0.00135135, 0.01081081]`.
- PKIS2 scaffold exclusion was `+0.00875000`, interval
  `[-0.00164063, 0.01867187]`, rather than the original `-0.00250000`.

## Tail risk and stability

Chemotype cases delayed by at least ten assays fell from 25 to 21, missing the
frozen maximum of 20 by one case. Gate decisions agreed in 95.03% of 50 fixed
80% reference subsamples per compound. Mean top-10 order stability was
`0.77737` versus `0.78842` for original coverage; difference `-0.01105`,
interval `[-0.01899, -0.00339]`.

## Frozen decisions

| Criterion | Result |
|---|---:|
| Chemotype retention at least 75%, CI lower > 0 | PASS |
| PKIS2 standard nonnegative and noninferior | PASS |
| PKIS2 >80 nonnegative and noninferior | FAIL |
| Klaeger primary noninferior | PASS |
| Chemotype large delays at most 20 | FAIL (21) |
| Conformal order stable | FAIL |

## Interpretation

Chemical-support gating successfully identifies when profile diversity is
useful and largely removes the regime reversals without inspecting activity
labels. The remaining failures are consistent with unconstrained coverage
being allowed to move a target arbitrarily far from the marginal order. A
separately frozen displacement-bounded model therefore addresses tail risk by
construction rather than changing the conformal threshold.

This is retrospective within-panel binding-liability discovery, not evidence
of cellular engagement, clinical safety, toxicity, efficacy, or therapeutic
suitability.
