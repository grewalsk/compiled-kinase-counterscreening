# PKIS2 construct-resolution audit protocol

Status: **frozen before execution**

This is a post-result biological-validity audit. It cannot alter the frozen
PKIS1 external-validation verdict or be used to tune the rank-deadline method.

## Question

Does treating each of the 406 PKIS2 assay constructs as a separate candidate,
rather than max-collapsing them to 392 parent labels, materially change the
reported development-panel effects?

## Data and construction

- Source: Drewry et al. (2017), PLOS ONE, Supporting Table S4, worksheet
  `Table 4 - PKIS2 %Inh`.
- Required workbook SHA-256:
  `48ead22a6616488b51a35060d65e40f100f9a59af86e11f0a7251e3b333e2a2c`.
- Apply the existing empty-identifier-row exclusion, largest-organic-fragment
  standardization, and median collapse of duplicate parent structures.
- Keep all 406 assay columns separately and in deterministic lexical order.
  Do not apply `pkis2_parent_target`.
- Use strict `>90%` and `>80%` single-concentration inhibition labels. Missing
  cells remain unavailable candidates. No inclusive-threshold substitution is
  permitted after seeing the results.

## Fixed policies and conditions

Use the already-reported Morgan/Tanimoto reference weighting, marginal order,
unbounded coverage order, and rank-deadline bounded coverage order with
maximum displacement K=9 and assay budget 20. Evaluate:

1. `>90`, standard leave-one-compound-out references;
2. `>80`, standard references;
3. `>90`, Bemis--Murcko scaffold exclusion;
4. `>90`, source-chemotype exclusion.

## Outcomes and interpretation

For bounded and unbounded coverage versus marginal ranking, report paired mean
AUDC differences, 95% compound-bootstrap confidence intervals, two-sided
paired sign-flip tests, Holm correction across the four bounded contrasts,
first-hit wins/losses/ties, and counts of censored delays of at least ten.
Also report label prevalence and compare the construct-resolved bounded effect
with the previously reported parent-collapsed effect.

The audit succeeds scientifically regardless of sign. A sign reversal or loss
of the PKIS2 chemotype gain must be reported as a limitation. Agreement does
not validate interchangeability of constructs; it only shows that the
aggregate algorithmic result is not created by the parent max-collapse.

