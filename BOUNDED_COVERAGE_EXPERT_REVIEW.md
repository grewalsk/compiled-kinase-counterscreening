# Required kinase-assay expert review

Status: **not yet completed**

This checklist is a mandatory scientific-validity handoff for the
rank-deadline coverage paper. It is not an invitation to tune the algorithm or
select favorable outcomes. The reviewer should record their name, relevant
expertise, review date, evidence consulted, and any requested corrections.

## 1. Parent-target collapse

Review the PKIS1 rule that averages repeat compound--assay rows and then takes
the maximum across Nanosyn assays sharing a `TARGET_CHEMBL_ID`.

- Are any merged assays mutant, fusion, complex, or construct measurements
  whose biological interpretation should remain distinct?
- Is taking the maximum defensible for a first-strong-profile-activity task?
- Does the construct-level robustness analysis adequately delimit the effect
  of this choice?

Required record: identify representative merged target IDs and source rows,
state whether the parent-level endpoint is biologically interpretable, and
request wording changes if needed.

## 2. Cross-panel comparability

Review the decision to report, but never pool, three modalities:

- Nanosyn PKIS1 percent inhibition at 1 micromolar;
- PKIS2 percent inhibition at 1 micromolar;
- Klaeger/Kinobeads binding `K_d` values.

Required record: state which cross-panel comparisons are legitimate at the
level of algorithmic transfer and which biochemical interpretations must not
be made. In particular, confirm or reject the manuscript's use of
“profile-activity discovery” rather than potency, cellular engagement, or
clinical selectivity.

## 3. Operational activity labels

Review strict `>80%` and `>90%` inhibition labels for PKIS1/PKIS2 and the
Klaeger label defined relative to a compound's best measured documented-target
`K_d`.

Required record: judge whether these are defensible operational labels for a
retrospective first-hit counter-screen task, while remaining unsuitable as
standalone potency or therapeutic-window claims. Flag assay-specific controls
or concentration-response evidence that would be required for stronger
language.

## 4. Case-level sanity audit

Inspect representative assay orders and source records for:

- well-known promiscuous kinase inhibitors;
- compounds with unusually sparse measurements;
- targets prone to assay interference or construct-dependent behavior;
- the largest bounded-versus-marginal improvements and losses;
- the 35 external large delays made by unconstrained coverage.

Required record: determine whether any headline pattern is explained by data
artifacts, target naming errors, duplicated chemistry, or an implausible
counter-screen decision. Every flagged case must be reported; favorable cases
must not be selected alone.

## Sign-off boundary

The computational paper may truthfully report that the above review is
outstanding. It must not claim biological validation until all four sections
are completed. Expert feedback may correct data semantics and interpretation;
any algorithm or endpoint change prompted after results must be labeled
post hoc and re-evaluated on a genuinely new validation set.
