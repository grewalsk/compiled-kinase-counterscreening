# PKIS1 external-validation pre-outcome amendment 2

Version: 1.0

Date frozen: 2026-09-05

Status: issued after the second input-validation stop and before any PKIS1
method order, aggregate activity statistic, or performance outcome was
computed

The raw PKIS1 archive contains finite percent-inhibition measurements outside
the nominal interval `[0, 100]`, a possible feature of unbounded assay readouts
and experimental noise. The second input-only validation stopped at this
check; no method was run.

The effective rule is clarified as follows:

1. Retain every finite observed raw inhibition value without clipping.
2. Apply the already frozen strict thresholds `>90` and `>80` to the retained
   values after repeat averaging and construct collapse.
3. Report counts below zero and above 100 in output metadata.

Clipping would not change binary labels below zero or above 100, but retaining
the archived measurement is more faithful to the source. This amendment
changes no label at either decision threshold, representation, ranking method,
parameter, condition, statistical test, or success criterion.
