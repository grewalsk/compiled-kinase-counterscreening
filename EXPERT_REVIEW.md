# Mandatory subject-matter review before submission

Reviewer requested: a scientist with practical kinase binding/selectivity
profiling experience. This sign-off cannot be replaced by statistical or
language-model review.

## Case-construction review

- [ ] Confirm that the Klaeger/Kinobeads Kd rows used here are comparable
  enough for within-source retrospective ranking.
- [ ] Inspect a random sample of at least 15 compounds plus every case with an
  unusual documented-target set; compare ChEMBL mechanism rows against the
  source and known polypharmacology without changing the frozen analysis.
- [ ] Decide whether ChEMBL `direct_interaction=1` is defensible as a reference
  intended-target exclusion, and document counterexamples.
- [ ] Confirm that choosing the best exact documented-target Kd as `theta_i`
  does not create a biologically misleading comparison for multi-target drugs.
- [ ] Judge whether the tenfold primary and sqrt(10) sensitivity margins are
  credible operational profiling thresholds. They must not be described as
  toxicity or clinical-safety thresholds.

## Assay and interpretation review

- [ ] Verify the interpretation of exact versus `>30,000 nM` measurements and
  the proof that primary/sensitivity labels are censor-safe.
- [ ] Review whether competition binding in lysate may be distorted by protein
  abundance, probe competition, isoforms, or assay-specific missingness.
- [ ] Ensure paper language consistently says “measured within-panel binding
  liability,” not off-target effect, functional inhibition, or adverse event.
- [ ] Review the top discovered target families for assay artifacts or trivial
  panel-wide promiscuity before any biological narrative is written.
- [ ] Confirm that no result is used to recommend a molecule, indication,
  dosing decision, or therapeutic action.

## Submission decision

- Reviewer name/role: ______________________________
- Cases examined and date: _________________________
- Required corrections: ____________________________
- Approve benchmark construction as operationally meaningful? YES / NO
- Approve claim language and limitations? YES / NO
- Signature or documented email reference: ______________________________

If either final answer is NO, the computational result may be released as a
negative/diagnostic audit but must not replace the fallback submission.
