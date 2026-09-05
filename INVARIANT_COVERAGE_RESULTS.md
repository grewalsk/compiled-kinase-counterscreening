# Similarity-stratified invariant coverage: run report

## Bottom line

The complete frozen experiment ran successfully, but the proposed invariant
model is **not an absolute improvement** over the existing method. All frozen
within-family success flags were true because invariant coverage (M3) strongly
outperformed invariant marginal ranking (M2). A post-run validity audit showed
that this contrast was inflated by a large loss in M2, not by M3 surpassing the
original coverage model (M1).

This is a useful negative modeling result, not a successful replacement model.
It should not be presented as having fixed the original reversals.

## Execution and scope

- Status: `EXPLORATORY_MODELING_COMPLETE`
- Wall time: 285.63 seconds (4.76 minutes)
- Hardware requirement: CPU only; a GPU gives no material advantage
- Paid API use: none; cost USD 0
- Conditions: all four frozen Klaeger/ChEMBL 30 conditions and all four frozen
  PKIS2 conditions
- Cases: 111 Klaeger compounds and 640 PKIS2 compounds per condition
- Inference: 10,000 paired bootstrap replicates and 100,000 paired sign flips
- Stability: 50 fixed 80% reference subsamples for each of 640 PKIS2 chemotype
  cases
- Leakage control: every outer held-out profile was unavailable to ranking and
  to the reference-only gate; pseudo-held-out references were removed from
  their own inner training pools

This was exploratory, post hoc model development motivated by the completed
mechanism audit. It does not revise the failed confirmatory Klaeger decision.

The completed full run used Python 3.12.8, RDKit 2025.03.5, SciPy 1.16.1,
pandas 2.2.2, matplotlib 3.10.0, and NumPy 2.2.6. The Colab bundle pins NumPy
2.0.2; all unit tests, an all-condition smoke run, and the complete post-run
statistical audit passed under the pinned environment. The audit JSON and CSV
were byte-identical across the two NumPy versions. A second full pinned local
run reached PKIS2 but was terminated by the host process-pool resource limit;
it was not substituted for the completed run. This distinction is recorded
rather than claiming an unperformed byte-identical full rerun.

## Models

- M0: original similarity-weighted marginal ranking
- M1: original similarity-weighted greedy profile coverage
- M2: similarity-stratified invariant marginal ranking
- M3: similarity-stratified invariant greedy profile coverage
- M4: reference-only gated choice between M2 and M3

M2 and M3 replace fine-grained weights inside five chemical-similarity strata
with the group average while preserving each stratum's total mass.

## Frozen result versus absolute validity result

The frozen primary contrast was M3 minus M2. That is a valid matched comparison
of acquisition rules under invariant weights, but it is insufficient to claim
that the new model improved absolute performance. The required diagnostic is
M3 minus M0 and, more importantly for replacement, M3 minus M1.

| Condition | Frozen M3-M2 AUDC | Absolute M3-M0 AUDC | Absolute M4-M0 AUDC |
|---|---:|---:|---:|
| Klaeger primary | +0.0212 [-0.0036, +0.0442] | +0.0018 [-0.0225, +0.0239] | -0.0014 [-0.0252, +0.0216] |
| Klaeger delta 0.5 | +0.0207 [-0.0086, +0.0532] | +0.0194 [-0.0099, +0.0518] | 0.0000 [-0.0239, +0.0266] |
| Klaeger scaffold exclusion | +0.0207 [-0.0041, +0.0441] | +0.0014 [-0.0216, +0.0243] | -0.0018 [-0.0252, +0.0207] |
| Klaeger remove 10 nearest | +0.0149 [-0.0108, +0.0392] | +0.0167 [-0.0095, +0.0414] | +0.0207 [+0.0027, +0.0374] |
| PKIS2 >90 standard | +0.0792 [+0.0580, +0.1009] | -0.0711 [-0.0906, -0.0520] | -0.0659 [-0.0834, -0.0483] |
| PKIS2 >80 sensitivity | +0.0817 [+0.0591, +0.1046] | -0.0755 [-0.0942, -0.0568] | -0.0710 [-0.0884, -0.0543] |
| PKIS2 >90 scaffold | +0.0832 [+0.0602, +0.1059] | -0.0544 [-0.0748, -0.0343] | -0.0507 [-0.0688, -0.0327] |
| PKIS2 >90 source chemotype | +0.0805 [+0.0588, +0.1027] | +0.0309 [+0.0105, +0.0513] | +0.0270 [+0.0073, +0.0465] |

Values are mean paired differences in AUDC through assay budget 20; brackets
are paired-bootstrap 95% confidence intervals. Positive values favor the first
method named in the contrast.

## Why the frozen criteria passed but the model did not

On PKIS2 source-chemotype exclusion, M2 was worse than M0 by -0.0496 AUDC
[-0.0637, -0.0365]. M3 recovered +0.0805 relative to this weakened comparator,
but its absolute gain over M0 was only +0.0309. This retains 57.6% of the
original M1-M0 gain of +0.0537, below the frozen intended threshold of 75%.
M3 was worse than M1 by -0.0227 AUDC.

M4 retained 50.4% of the original chemotype gain but was worse than M1 by
-0.0266 AUDC [-0.0378, -0.0157]. Relative to M0 it produced 41 cases at least
10 assays earlier and 26 at least 10 assays later. The original M1 comparison
had 55 such advances and 25 such delays. The gate therefore did not reduce the
absolute tail-risk count while preserving performance.

For PKIS2 standard, >80 sensitivity, and scaffold exclusion, M3 remained
significantly worse than M0. The apparent reversal elimination existed only
against M2. It did not eliminate the real absolute reversals.

On Klaeger primary, M3-M0 was +0.0018 with a confidence interval crossing zero;
M4-M0 was -0.0014 with a confidence interval crossing zero. There is no evidence
of a meaningful improvement or degradation there.

## Stability and gate diagnostics

Mean top-10 Jaccard stability under 80% reference resampling was 0.7890 for M1
and 0.7904 for M3. The difference was +0.00144 with 95% CI
[-0.00207, +0.00493]. The frozen strict-inequality flag passed, but the effect
is negligible and statistically inconclusive.

The gate selected M3 for 79.2% of PKIS2 chemotype cases and between 55.9% and
78.4% of Klaeger cases. It reduced some losses relative to M2, but could not
recover the chemical ranking information removed by group averaging.

## Scientific interpretation

The result supports three narrower conclusions:

1. Fine-grained chemical-similarity weighting contains useful information for
   marginal assay ranking; averaging it away is destructive on PKIS2.
2. Greedy joint-profile coverage is more robust than marginal ranking to this
   destruction, which creates a large conditional M3-M2 contrast.
3. A reference-only local gate cannot rescue a representation that has already
   discarded important ranking information.

The result does not establish that group invariance improves counter-screen
selection. It also changes the interpretation of the earlier whole-profile
shuffle: a larger coverage-minus-marginal contrast can arise because the
marginal comparator degrades, so absolute performance must accompany every
relative mechanism contrast.

## Paper decision

Do not promote M3 or M4 as the main new algorithm. The original coverage method
remains the best tested absolute method for PKIS2 chemotype exclusion. If this
experiment appears in the paper, use it as a falsification/ablation showing why
coarse group averaging and local gating are insufficient. Its strongest value
is methodological: matched relative gains can be misleading when a
representation perturbation damages the baseline.

Any next algorithm should preserve M0's fine-grained marginal scores and alter
only the redundancy correction, with model selection evaluated against both M0
and M1 on absolute endpoints. No further parameter search should reuse these
eight reported conditions as if they were untouched test data.

## Reproduction

The frozen Colab notebook is `notebooks/invariant_coverage_colab.ipynb` and its
upload bundle is `invariant_coverage_colab_bundle.zip`. The notebook verifies
the prospective freeze, installs pinned dependencies, validates data, downloads
and checks the exact PKIS2 source, runs deterministic tests, executes all eight
conditions, and exports the result archive.

The original frozen outputs are in `invariant_coverage_output/summary.json` and
`comparison_summary.csv`. The explicitly post-run audit is in
`absolute_validity_audit.json` and `absolute_comparisons.csv`; it is separately
hashed rather than retroactively inserted into the prospective freeze.
