# Ethics, biosafety, and scientific-validity statement

## Scope

This is a retrospective methods study of already-public small-molecule kinase
binding measurements. It performs no wet-lab work, synthesis, pathogen work,
human-subject research, animal research, dosing analysis, or therapeutic
recommendation. It does not establish functional inhibition, cellular target
engagement, toxicity, efficacy, or clinical safety.

## Foreseeable misuse

The main foreseeable harm is epistemic rather than biosynthetic: a user could
mistake a low posterior residual-risk score for evidence that a compound is
selective or safe. The experiment directly tests that use and finds it
unsupported. The code and paper therefore use `clear-within-panel`, never
global or clinical clearance; require exact false-clear auditing; and recommend
abstention when the finite-sample criterion is not met.

The artifact contains public compound structures and binding records but no
new synthesis route, pathogen target prioritization, organism engineering, or
operational protocol that materially increases biological capability beyond
the source databases. No compound is ranked for therapeutic progression.

## Scientific-validity risks and mitigations

- **Assay proxy:** lysate competition binding is not biochemical inhibition or
  cellular pharmacology. Claims are restricted to measured within-panel
  binding and require kinase-assay expert review.
- **Mechanism annotation:** ChEMBL mechanism rows may be incomplete or reflect
  clinically broad polypharmacology. Every documented direct human-kinase
  target represented in the panel is excluded, and sampled mappings must be
  manually checked before submission.
- **Censoring and missingness:** missing cells are never made negative;
  right-censored values are used as negatives only where the censor bound is
  provably above the label threshold.
- **Chemical leakage:** held-out parents are excluded; exact-scaffold and
  leave-ten-nearest-reference audits are reported. The stronger audit is
  inconclusive, so extrapolation to new chemical series is not claimed.
- **Single-source external validity:** the benchmark is a fixed census from one
  standardized source. Heterogeneous panels are not pooled merely to create an
  appearance of replication.
- **Multiple analysis stages:** primary artifacts were frozen before result
  access; later attribution and robustness tests were separately frozen and
  marked post-result.
- **Finite sample:** the 111 cases are not presented as a power-designed random
  sample. Paired uncertainty intervals describe resolution for this census.

## Data and software terms

The compact derived ChEMBL records are redistributed as an adaptation under
CC BY-SA 3.0 with ChEMBL and Klaeger et al. attribution. The original 22 GB
database and any access credentials are not redistributed. RDKit is used under
its BSD 3-Clause license. The artifact contains no private, patient, or
personally identifying information.
