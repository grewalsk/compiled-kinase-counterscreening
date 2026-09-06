# Which Counter-Screen Next?

[![Open rank-deadline study in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/grewalsk/compiled-kinase-counterscreening/blob/main/notebooks/rank_deadline_coverage_colab.ipynb)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/grewalsk/compiled-kinase-counterscreening/blob/main/notebooks/compiled_profile_coverage_colab.ipynb)

[![Open invariant-model audit in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/grewalsk/compiled-kinase-counterscreening/blob/main/notebooks/invariant_coverage_colab_public.ipynb)

Reproducible CPU-only benchmark of sequential kinase counter-screen selection.
The scientific task is retrospective: given a compound structure, documented
direct kinase targets, and a limited assay budget, rank which other target in
a common kinase-drug Kinobeads panel should be measured next to discover an
operational relative-affinity event. This study-defined event is not a safety
margin.

This project does **not** infer biochemical inhibition, cellular target
engagement, therapeutic benefit, toxicity, or clinical safety. It does not
recommend compounds or treatments.

## Headline: rank-deadline coverage and external validation

The submission-facing method is **rank-deadline profile coverage**. It greedily
diversifies kinase assays while constraining every assay to remain within nine
positions of a chemical-similarity-weighted marginal order. This gives an
exact first-hit backstop: the method cannot delay the first recorded activity
by ten or more assays relative to marginal ranking.

The result is deliberately not presented as universal ranking superiority.
All eight development-condition point estimates were nonnegative and the
method retained 80.6% of the unconstrained PKIS2 chemotype-shift gain, but a
frozen external PKIS1 evaluation rejected mean transfer. The external `>80%`
effect was -0.02065 AUDC (95% interval [-0.03062, -0.01166]). The formal risk
result nevertheless held: zero large delays across all 4,428 development and
external compound-condition evaluations, versus 163 for unconstrained
coverage.

Start with:

- `paper/bounded_coverage.pdf`: submission-formatted manuscript;
- `PKIS1_EXTERNAL_VALIDATION_RESULTS.md`: external result and robustness;
- `RISK_CONTROLLED_COVERAGE_RESULTS.md`: development result;
- `notebooks/rank_deadline_coverage_colab.ipynb`: end-to-end CPU rerun;
- `BIOLOGICAL_VALIDITY_REVIEW.md`: completed source/data semantics and case
  review, with the remaining independent-human-attestation boundary;
- `BOUNDED_COVERAGE_EXPERT_REVIEW.md`: original review checklist and status.

No GPU, paid API, or model training is needed. A standard Colab CPU is the
recommended runtime.

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

`notebooks/invariant_coverage_colab_public.ipynb` is the one-click public
launcher for the similarity-stratified group-averaging experiment and its
reference-only pseudo-holdout gate. It downloads and verifies the separately
frozen bundle, runs all eight conditions, performs the explicitly post-run
absolute-baseline audit, and downloads the complete artifacts. The full run is
CPU-only and costs USD 0. `notebooks/invariant_coverage_colab.ipynb` is retained
unchanged as the prospectively frozen pre-result entry point.

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
- `BIOLOGICAL_VALIDITY_REVIEW.md`: completed literature/data-semantics review
  and exact boundary on remaining independent human attestation.
- `src/audit_biological_validity.py`: construct map, threshold-boundary census,
  flagged-case selection, and source-record export.
- `src/audit_pkis2_constructs.py`: frozen 406-assay PKIS2 resolution audit.
- `src/audit_pkis2_boundaries.py`: complete post hoc inclusive-boundary grid.
- `biological_validity_review_output/`, `pkis2_construct_audit_output/`, and
  `pkis2_boundary_audit_output/`: machine-readable audit results and manifests.
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
and stress tests. The submission-facing audit is
`BIOLOGICAL_VALIDITY_REVIEW.md`: it completes the reproducible source/data
checks but is not independent human attestation or biological validation. The
study is framed as measured within-panel operational-event discovery. The
frozen primary contrast passes, but its attribution and safe-clearance claims
do not. That mixed result is the paper rather than a failure to be hidden.

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
