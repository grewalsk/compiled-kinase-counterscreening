# Which Counter-Screen Next?

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/grewalsk/compiled-kinase-counterscreening/blob/main/notebooks/compiled_profile_coverage_colab.ipynb)

[![Open invariant-model audit in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/grewalsk/compiled-kinase-counterscreening/blob/main/notebooks/invariant_coverage_colab.ipynb)

Reproducible CPU-only benchmark of sequential kinase counter-screen selection.
The scientific task is retrospective: given a compound structure, documented
direct kinase targets, and a limited assay budget, rank which other kinase in
a common Kinobeads panel should be measured next to discover a binding
liability relative to the compound's measured on-target affinity.

This project does **not** infer biochemical inhibition, cellular target
engagement, therapeutic benefit, toxicity, or clinical safety. It does not
recommend compounds or treatments.

## New compiled profile-coverage experiment

The post-benchmark method reframe is implemented in
`notebooks/compiled_profile_coverage_colab.ipynb`. Open it in a standard Colab
CPU runtime, download `compiled_coverage_colab_bundle.zip` from this repository,
and upload it when prompted. The notebook:

1. verifies the separate method freeze;
2. installs pinned CPU dependencies and runs deterministic unit tests;
3. validates the unchanged compact Klaeger/ChEMBL benchmark;
4. downloads and hash-checks the CC BY PKIS2 S4 workbook;
5. runs the exact first-hit compiler checks, matched marginal and coverage
   methods, chemical-series exclusions, external transfer evaluation,
   statistics, exact small-instance diagnostics, support abstention, and
   batching; and
6. downloads a complete result archive with logs and manifests.

The new analysis is frozen in `COMPILED_COVERAGE_PROTOCOL.md` and
`COMPILED_COVERAGE_FREEZE.sha256`. Its implementation is
`src/run_compiled_coverage.py`, with deterministic tests in
`src/test_compiled_coverage.py`. PKIS2 is not pooled with the Klaeger `Kd`
data and is not described as an off-target benchmark because its public
workbook does not provide complete per-compound intended-target annotations.

The full method run uses no GPU and no paid API. Use `SMOKE_TEST=True` only to
check execution; smoke-test estimates are not paper results.

## Exploratory invariant-model audit

`notebooks/invariant_coverage_colab.ipynb` evaluates similarity-stratified
group-averaged weights and a reference-only pseudo-holdout gate on all eight
existing conditions. Upload `invariant_coverage_colab_bundle.zip` when
prompted. The complete CPU run takes about five minutes after dependency
installation and costs USD 0.

This experiment is a documented negative result. Its frozen within-family
criteria passed, but an absolute-baseline audit showed that group averaging
damaged the marginal comparator and did not outperform the original coverage
model. Read `INVARIANT_COVERAGE_RESULTS.md` before interpreting the frozen
success flags. The original prospective freeze and the post-run validity audit
have separate SHA-256 manifests, and `invariant_coverage_results_full.zip`
contains the complete code, results, diagnostics, and report.

## Fast path: Google Colab

1. Open `notebooks/counterscreen_active_search_colab.ipynb` in Colab with a
   standard CPU runtime.
2. Run every cell in order.
3. When prompted, upload `counterscreen_colab_bundle.zip`.
4. The notebook verifies all three chronological freeze manifests, installs
   pinned dependencies, validates every dataset invariant, executes the frozen
   analysis and the explicitly post-result stress tests, verifies every
   reported numerical value, and downloads `counterscreen_results.zip`.

No API key, Google Drive mount, GPU, model training, or source ChEMBL database
is required. The expected API cost is USD 0. On the recorded standard Colab
CPU runtime, pinned environment setup took about 29 minutes, the frozen study
took 29 seconds, and the post-result audit took 28 seconds. Colab resource
availability is not guaranteed by Google; the results archive and console
logs make interrupted runs obvious.

## Artifact map

- `PROTOCOL.md`: frozen scientific question, case construction, policies,
  outcomes, statistics, and go/no-go gates.
- `PREANALYSIS_FREEZE.sha256`: hashes fixed before adaptive implementation.
- `IMPLEMENTATION_CLARIFICATIONS.md`: pre-run resolutions of underspecified
  numerical details.
- `IMPLEMENTATION_FREEZE.sha256`: hashes of the adaptive runner and its
  implementation clarifications before result access.
- `POSTHOC_ROBUSTNESS_PROTOCOL.md` and `POSTHOC_FREEZE.sha256`: added tests
  fixed after the primary result but before their implementation was run.
- `EXPECTED_RESULTS.json` and `src/verify_outputs.py`: post-run reporting
  values and an exact rerun/transfer checker.
- `POSTRUN_DEVIATIONS.md`: transparent account of the repaired Colab logging
  wrapper; no scientific output was changed.
- `BACKGROUND.md`: closest primary literature and narrow novelty boundary.
- `DATA_CARD.md`: exact dataset construction, exclusions, licensing, and known
  validity limitations.
- `ETHICS_AND_RISKS.md`: biosafety, misuse, licensing, and scientific-validity
  risk assessment.
- `EXPERT_REVIEW.md`: mandatory kinase-assay review checklist.
- `src/extract_panel.py`: read-only extraction from ChEMBL 30.
- `src/validate_dataset.py`: independent compact-dataset validator.
- `src/run_study.py`: frozen policies, statistics, plots, and go/no-go report.
- `src/run_robustness.py`: static-prior attribution, stronger nearest-neighbor
  exclusion, kinase-taxonomy, and absolute-affinity stress tests.
- `data/derived/`: compact benchmark and row-level provenance.
- `notebooks/counterscreen_active_search_colab.ipynb`: complete compute entry
  point.
- `paper/main.tex`: five-page AI4DD workshop manuscript source.
- `output/pdf/which_counter_screen_next_ai4dd_2026.pdf`: rendered manuscript.

## Result interpretation

`output/summary.json` is the frozen decision source;
`posthoc_output/summary.json` contains only explicitly post-result attribution
and stress tests. A computational pass still is not a biological-validation
pass: `EXPERT_REVIEW.md` must be completed, and the study must be framed as
measured within-panel binding-liability discovery. The frozen primary contrast
passes, but its attribution and safe-clearance claims do not. That mixed result
is the paper rather than a failure to be hidden.

## Licensing and attribution

The derived database records are an adaptation of ChEMBL 30 and are provided
under the [Creative Commons Attribution-ShareAlike 3.0 Unported
License](https://creativecommons.org/licenses/by-sa/3.0/). ChEMBL's official
[licensing FAQ](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions)
requires attribution and share-alike redistribution of adaptations. Cite both
ChEMBL and Klaeger et al., *Science* 2017, DOI
[`10.1126/science.aan4368`](https://doi.org/10.1126/science.aan4368).
RDKit is distributed under the [BSD 3-Clause
License](https://github.com/rdkit/rdkit/blob/master/license.txt).
