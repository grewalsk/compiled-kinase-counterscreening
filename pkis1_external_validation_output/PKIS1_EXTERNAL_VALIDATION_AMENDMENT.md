# PKIS1 external-validation pre-outcome amendment

Version: 1.0

Date frozen: 2026-09-05

Status: issued after input validation and before any PKIS1 method order,
aggregate activity statistic, or performance outcome was computed

The raw PKIS1 archive contains missing `VALUE` cells. This was detected when
the input-only validator stopped; no method was run. Step 4 of the original
protocol said to require finite `VALUE` fields, whereas step 10 explicitly
said not to impute missing activities and to restrict candidates to observed
collapsed activities. Those statements conflict for the raw file.

The effective rule is therefore clarified as follows:

1. Parse `VALUE` numerically, preserving blank cells as missing.
2. Require every nonmissing value to be finite and in `[0, 100]`.
3. Average repeats and collapse constructs with missing values skipped.
4. If an entire compound--parent-target cell remains missing, retain it as
   missing. Never convert it to an inactive label.
5. Include a target in a held-out compound's candidate set only if that
   held-out collapsed value is observed.
6. Report the raw and collapsed missing-value counts in output metadata.

This amendment implements the protocol's already stated missing-data and
candidate-mask policy. It changes no activity threshold, representation,
ranking method, parameter, condition, statistical test, or success criterion.
