# External PKIS1 validation protocol

Version: 1.0, frozen before computing any PKIS1 method order, aggregate
activity statistic, or performance outcome

Date frozen: 2026-09-05

## Status and confirmatory question

`BOUNDED_COVERAGE` was a prespecified component in the final PKIS2/Klaeger
experiment, but it was promoted to the headline method only after the
conformal wrapper failed. This protocol therefore evaluates that now-fixed
method, with no tuning, on PKIS1. PKIS1 was not used to design the method or
select its displacement bound.

The confirmatory question is whether displacement-bounded profile coverage
transfers to a separate 1 uM kinase inhibition matrix: it should remain
noninferior to chemical weighted-marginal ranking in ordinary conditions,
improve average first-hit discovery when close chemical references are
removed, and retain its deterministic protection against delays of ten or
more assays.

This is a retrospective external validation, not a prospective wet-lab test.
PKIS1 and PKIS2 are separate compound panels from the same broader Published
Kinase Inhibitor Set ecosystem, so the validation is stronger than another
split of PKIS2 but not fully independent across laboratory, assay family, or
target class.

## Audit trail before freezing

Before this file was frozen, only data provenance, file hashes, field names,
non-outcome metadata, row/compound/assay counts, structure counts, and exact
structure overlap with the development datasets were inspected. During a
schema command, the first two rows of the authors' processed continuous
matrix were inadvertently printed. No aggregate PKIS1 activity distribution,
method order, condition-level outcome, or performance comparison was viewed.
This deviation prevents describing the run as perfectly blinded, although the
two example rows were not used to alter any rule or criterion below.

## Source, license, and immutable inputs

Use the raw compressed ChEMBL export archived by the authors of Zhang et al.
(2019) in `SpencerEricksen/informers`:

- repository commit: `5fd3934f5789c371026fc9eece1846ff1294122b`
- relative path:
  `data/original_data/pkis1/PKIS_screening_data.csv.gz`
- raw-file SHA-256:
  `81d7f9f82f7ee8e6b0f38dafe523da7254aaabe9449758a056372511d7868ad0`
- source article/data archive: DOI `10.1371/journal.pcbi.1006813` and
  Zenodo DOI `10.5281/zenodo.3354432`
- original PKIS1 description: Drewry et al. (2014), DOI
  `10.2174/1568026613666131127160819`

The file is a ChEMBL-derived dataset. ChEMBL data are available under Creative
Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0). Any redistributed
derived data must preserve attribution, the ChEMBL identifiers, the release or
archived source information when available, and share-alike terms. The
authors' PLOS article is CC BY 4.0. No authentication or paid access is
required.

The PKIS2 development workbook remains pinned by SHA-256
`48ead22a1f860cd0d5096fa87d5acd329f722fe8d65e693bb0be682a333e2a2c`.
The locally derived Klaeger cases are used only to audit exact structure
overlap; no PKIS1 outcome is joined to either development dataset.

## Exact PKIS1 construction and exclusions

1. Read the raw CSV with strings preserved for identifiers.
2. Retain rows whose `ASSAY` contains the literal strings `" 1 uM"` and
   `"[Nanosyn]"`.
3. Require `ENDPOINT == "Inhibition"`, `RELATION == "="`,
   `SPECIES == "Homo sapiens"`, and `SOURCE == "GSK Published Kinase
   Inhibitor Set"`. `UNITS` is expected to be `%` or blank; blank units are
   retained because assay name, endpoint, relation, source, and the archived
   preprocessing identify the same percent-inhibition panel.
4. Parse `VALUE` numerically and require finite values in `[0, 100]`.
5. For repeat rows having the same `CHEMBL_ID` and `ASSAY_CHEMBL_ID`, take the
   arithmetic mean, matching the archived authors' preprocessing.
6. Collapse mutant/construct assays sharing `TARGET_CHEMBL_ID` by the maximum
   of their repeat-averaged percent inhibition, mirroring the development
   benchmark's parent-target collapse. Sort target keys by
   `TARGET_CHEMBL_ID`; target names are metadata only.
7. Standardize each compound by RDKit largest-fragment selection with
   `preferOrganic=True`, then canonical isomeric SMILES. If multiple
   `CHEMBL_ID` values map to one parent structure, take the per-target median.
8. Remove every standardized parent structure appearing in either development
   benchmark before constructing labels or reference sets. The frozen
   structure-only audit found eight PKIS1 structures in PKIS2 and none in the
   Klaeger cases.
9. Expected structure-only counts are 366 ChEMBL compound IDs, 364
   standardized PKIS1 parent structures, 224 assay identifiers, 200 parent
   target ChEMBL IDs, and 356 parent structures after development-overlap
   exclusion. Any mismatch stops the run.
10. Do not impute missing activities. A target is a candidate for a held-out
    compound only when its collapsed activity is observed. Identical
    structures are excluded from reference sets even after grouping as a
    defensive check.

Define primary activity as percent inhibition strictly greater than 90 and
sensitivity activity as strictly greater than 80. These fixed thresholds
match the development experiment. They indicate strong biochemical profile
activity in this assay, not cellular target engagement, toxicity, efficacy,
or therapeutic suitability.

## Fixed representations and policies

Use the development pipeline without modification:

- RDKit Morgan bit fingerprints, radius 2 and 2,048 bits;
- Tanimoto similarities;
- similarity weights proportional to `exp(8 * Tanimoto)`;
- deterministic opaque-key tie breaks;
- assay budget 20;
- `MARGINAL`: rank parent targets by weighted activity probability;
- `UNBOUNDED_COVERAGE`: the original weighted greedy profile-coverage order;
- `BOUNDED_COVERAGE`: the frozen displacement-bounded greedy coverage order
  with `K=9`.

At output position `t`, `BOUNDED_COVERAGE` selects the earliest target whose
marginal-rank deadline `r(a)+9` has arrived. If none is urgent, it maximizes
uncovered weighted profile mass among targets with `r(a) <= t+9`, breaking
ties by marginal probability and then opaque target key. No parameter is fit
on PKIS1.

## Fixed validation conditions

Run all 356 held-out parent structures under four conditions:

1. `pkis1_gt90_standard`: primary labels; all eligible nonidentical
   references.
2. `pkis1_gt80_sensitivity`: sensitivity labels; same references.
3. `pkis1_gt90_scaffold`: primary labels; additionally remove references with
   the same exact Bemis--Murcko scaffold string as the held-out compound.
4. `pkis1_gt90_leave10_nearest`: primary labels; additionally remove the ten
   most Tanimoto-similar eligible references, breaking ties by compound index.

The scaffold and nearest-reference conditions form the prespecified
chemical-shift pair. PKIS1 has no source chemotype field, so no chemotype
condition is invented.

## Outcomes and statistical analysis

The primary scalar outcome is area under the any-hit discovery curve (AUDC)
over budgets 1--20. Also report censored cost to first activity, hit at budgets
1, 3, 5, 10, and 20, paired wins/losses/ties, mean assays earlier, and counts
of advances or delays of at least ten assays.

For each condition report paired `BOUNDED_COVERAGE - MARGINAL` and
`UNBOUNDED_COVERAGE - MARGINAL` effects. Use 10,000 paired compound
bootstraps and 100,000 paired sign flips with master seed 20260905. Report
two-sided 95% percentile intervals, two-sided nominal sign-flip p-values, and
Holm-adjusted p-values across the four bounded-versus-marginal condition
comparisons. No test will be removed because of zero-hit compounds.

For the chemical-shift macro-effect, average each compound's bounded-minus-
marginal AUDC difference across the scaffold and leave-10-nearest conditions,
then bootstrap/sign-flip those 356 paired compound-level averages. This keeps
the resampling unit at the compound rather than treating the two conditions as
independent.

As a descriptive limitation analysis, perform 50 fixed 80% reference
subsamples in the scaffold condition and report top-10 Jaccard stability for
bounded and unbounded coverage. Stability is not a success criterion because
the development result already established that the bounded rule does not
improve it.

## Frozen external-validation success criteria

The external validation succeeds only if all of the following hold:

1. In `pkis1_gt90_standard`, bounded-minus-marginal AUDC has point estimate
   at least zero and 95% lower endpoint at least `-0.015`.
2. In `pkis1_gt80_sensitivity`, bounded-minus-marginal AUDC has point estimate
   at least zero and 95% lower endpoint at least `-0.015`.
3. Both chemical-shift conditions have nonnegative point estimates and 95%
   lower endpoints at least `-0.015`.
4. The paired chemical-shift macro-effect is positive and its 95% lower
   endpoint is greater than zero.
5. Every condition has zero cases in which bounded coverage delays the
   censored first activity by ten or more assays relative to marginal ranking.

The noninferiority margin was fixed at 1.5 percentage points of AUDC, as in
the development sequence. With 356 held-out compounds and the development
standard deviations (`0.104`--`0.108`) for the comparable PKIS2 conditions,
the approximate two-sided 5%, 80%-power minimum detectable mean difference is
about `0.015`--`0.016`. The sample size is fixed by the public panel; this
calculation is a sensitivity justification, not a guarantee about the unseen
PKIS1 variance.

Failure does not authorize changing `K`, the representation, thresholds,
filters, conditions, success criteria, or baselines. All outcomes will be
retained and reported.

## Leakage, validity, and scope controls

- PKIS1 structures overlapping either development benchmark are removed
  before labels, weights, or references are constructed.
- A held-out compound's activities never enter weights, marginal scores,
  coverage gains, or method selection.
- Condition-specific exclusions are applied before similarities are converted
  to weights.
- Target names and therapeutic annotations are never features.
- Missing measurements are never encoded as negatives.
- All orders must be deterministic and all bounded orders must pass the
  displacement and first-hit-delay assertions.
- The biological claim is limited to prioritizing biochemical kinase
  counter-screens for discovery of recorded within-panel activity. No medical
  or therapeutic recommendation is made.
