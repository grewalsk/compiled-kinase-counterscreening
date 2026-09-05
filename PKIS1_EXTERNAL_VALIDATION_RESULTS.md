# External PKIS1 validation results

Date completed: 2026-09-05

Headline status: **the frozen composite failed; the mean-gain hypothesis did
not transfer, while the deterministic delay guarantee held exactly**

## What was tested

After the displacement-bounded method had been selected from the development
sequence, its parameters were fixed and tested once on the separate PKIS1
compound panel. The raw ChEMBL-derived Nanosyn file was pinned to repository
commit `5fd3934f5789c371026fc9eece1846ff1294122b` and SHA-256
`81d7f9f82f7ee8e6b0f38dafe523da7254aaabe9449758a056372511d7868ad0`.
Eight parent structures shared with PKIS2 were excluded before labels or
reference sets were constructed, leaving 356 compounds and 200 parent kinase
targets.

The protocol, two input-validation amendments, implementation, and tests were
separately frozen before any method outcome was computed. The amendments
preserved genuinely missing values as missing and retained finite assay
readouts outside `[0, 100]` without clipping; neither change alters a strict
`>80` or `>90` binary label. The audit trail discloses that a schema command
inadvertently printed two processed rows before the protocol freeze. No
aggregate outcome, order, or method comparison was viewed.

## Frozen result

All effects below are paired AUDC differences over assay budgets 1--20.
Intervals are 10,000 compound-bootstrap 95% intervals.

| Condition | Bounded minus marginal | 95% interval | Unbounded minus marginal | Bounded delays >=10 | Unbounded delays >=10 |
|---|---:|---:|---:|---:|---:|
| `>90`, standard | -0.00351 | [-0.01152, 0.00435] | -0.00969 | 0 | 5 |
| `>80`, standard | -0.02065 | [-0.03062, -0.01166] | -0.03876 | 0 | 14 |
| `>90`, scaffold exclusion | -0.00197 | [-0.01039, 0.00632] | -0.00716 | 0 | 5 |
| `>90`, remove 10 nearest | -0.00323 | [-0.01447, 0.00801] | -0.01039 | 0 | 11 |

The four-condition chemical-shift macro-effect was -0.00260, with 95%
interval [-0.01131, 0.00611]. At the `>80` endpoint, the two-sided paired
sign-flip p-value was `9.9999e-6` and the four-condition Holm-adjusted value
was `3.99996e-5`. The prespecified external composite therefore failed, and
the evidence rejects a universal mean-improvement claim.

The risk statement transferred. Bounded coverage produced no delay of ten or
more assays in any of the 1,424 external compound-condition evaluations,
whereas unconstrained coverage produced 35. Bounded-versus-unbounded point
estimates were positive in all four external conditions, including +0.00618
at `>90` standard (95% interval [0.00126, 0.01236]) and +0.01812 at `>80`
([0.00955, 0.02795]). Reference-resampling stability was similar externally:
-0.0023 bounded minus unbounded, with interval approximately
[-0.0097, 0.0052].

Across development and external evaluations together, the rank-deadline rule
had zero delays of ten or more assays in 4,428 compound-condition evaluations;
unconstrained coverage had 163. These are repeated condition-level evaluations,
not 4,428 independent compounds.

## Post-result robustness audit

The separately frozen, explicitly exploratory 18-cell audit varied the
activity threshold and retained all 224 assay constructs instead of collapsing
to 200 parent targets. It did not reverse the external result:

- construct-level standard effects were -0.01320 at `>80` (95% interval
  [-0.02219, -0.00492]) and -0.00899 at `>90`
  ([-0.01868, approximately 0]);
- parent-level standard effects were negative at every tested threshold
  (`>50`, `>70`, `>80`, and `>90`);
- two remove-ten-nearest cells became positive at lower thresholds. The
  `>50` effect was +0.02809 ([0.01110, 0.04495], Holm over 18 p=0.01408), and
  the `>70` effect was +0.02051 ([0.00590, 0.03539], Holm p=0.09684).

Those lower-threshold cells are hypothesis-generating regime indicators. They
do not rescue the frozen test or license selecting a new endpoint after seeing
the result.

## Defensible paper claim

The result separates a mathematical guarantee from a biological-performance
prediction. Rank deadlines cannot make an empirical profile-coverage model
transport across platforms or chemistry; they can guarantee that, when that
model is wrong, no target and no first recorded activity is displaced by more
than the declared budget. The method is therefore a backstop for an already
validated diversity objective, not a replacement for marginal ranking in an
unvalidated setting.

All endpoints are retrospective biochemical profile labels. They do not
establish cellular engagement, toxicity, efficacy, clinical safety, or a
therapeutic recommendation.
