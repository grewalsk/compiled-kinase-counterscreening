# Kinase-assay semantics and biological-validity review

Date: **2026-09-05**
Status: **evidence-backed computational/literature review complete; independent
human kinase-assay attestation not obtained**

This record closes the four concrete data/interpretation checks in
`BOUNDED_COVERAGE_EXPERT_REVIEW.md`. It was produced as an AI-assisted
adversarial audit of the primary papers, source records, code, and all flagged
cases. It must not be represented as a named human expert's sign-off. No
algorithm or frozen external endpoint was changed in response to this review.

## Evidence reviewed

- Drewry et al. 2014, the original PKIS description and Nanosyn/Caliper panel:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4435035/>
- Elkins et al. 2016, comprehensive PKIS characterization:
  <https://doi.org/10.1038/nbt.3374>
- Drewry et al. 2017, PKIS2 KINOMEscan experiment and public S4 workbook:
  <https://doi.org/10.1371/journal.pone.0181585>
- Zhang et al. 2019, archived PKIS matrices and target-wise modeling:
  <https://doi.org/10.1371/journal.pcbi.1006813>
- Klaeger et al. 2017, Kinobeads chemical-proteomics panel:
  <https://doi.org/10.1126/science.aan4368>
- ChEMBL licensing and current target records:
  <https://chembl.github.io/chembl-licensing/>,
  <https://www.ebi.ac.uk/chembl/explore/target/CHEMBL5550>, and
  <https://www.ebi.ac.uk/chembl/explore/target/CHEMBL5146>.

Reproducible evidence is in `biological_validity_review_output/`,
`pkis2_construct_audit_output/`, and `pkis2_boundary_audit_output/`. The audit
scripts are `src/audit_biological_validity.py`,
`src/audit_pkis2_constructs.py`, and `src/audit_pkis2_boundaries.py`.

## 1. Parent-target collapse

### Finding

The operation is not ordinary duplicate removal. In PKIS1, 11 of 200 parent
target IDs have multiple assays. These include ABL1 (wild type and six
mutants), EGFR (wild type and four mutants), KIT and PDGFRA (wild type and
three mutants each), RET, FLT3, and LRRK2 variants, CDK2/cyclin complexes, and
PKC-beta and LYN splice variants. In PKIS2, 13 of 392 parent labels combine
phosphorylation or autoinhibited states, cyclin complexes, catalytic and
pseudokinase domains, or the two kinase domains of RSK-family proteins.

Taking the maximum therefore defines a logical endpoint: **activity above the
threshold in any profiled construct assigned to that parent**. It is
interpretable for broad profile-activity discovery only if stated exactly that
way. It is not evidence that the constructs are interchangeable, that the
wild-type protein is active, or that one physical follow-up assay represents
the whole parent label.

Among 3,912 retained PKIS1 compound--multi-construct-parent cells with at
least two observed constructs, 230 (5.88%) are construct-discordant at `>80`
and 134 (3.43%) at `>90`. The parent endpoint is thus biologically meaningful
but materially broader than a wild-type or construct-specific endpoint.

### Aggregate robustness

- PKIS1's negative external result persists at 224-assay resolution:
  bounded-minus-marginal AUDC is -0.0132 (95% CI -0.0222 to -0.0049) at
  `>80`, and -0.0090 (-0.0187 to approximately 0) at `>90` under standard
  references. Collapse is not an explanation for the failed transfer.
- A separately frozen post-result PKIS2 audit retains all 406 assay columns.
  The central chemotype-shift effect is +0.0442 (0.0336 to 0.0554;
  Holm-adjusted p approximately 4e-5), versus +0.0433 (0.0321 to 0.0545)
  after parent collapse. The construct-resolved scaffold effect is +0.0109
  (0.0026 to 0.0193); standard `>90` and `>80` effects remain small and
  interval-compatible with zero. The reported PKIS2 gain is not created by
  max-collapse.

### Required paper boundary

Use “parent-level any-profiled-construct event” when discussing a collapsed
PKIS endpoint. Do not turn the returned parent label into a laboratory assay
order without first selecting a specific construct and platform.

## 2. Cross-platform comparability

The manuscript does not pool modalities, which is essential. PKIS1 reports
Nanosyn single-concentration biochemical percent inhibition; PKIS2 uses a
DiscoverX competitive-binding displacement assay in singlicate at 1 micromolar
and converts percent control to percent inhibition; the Klaeger-derived rows
are exact `Kd` records from a chemical-proteomics source represented in
ChEMBL 30. Drewry et al. explicitly used different selection criteria because
PKIS1 and PKIS2 used different platforms and panel sizes.

Legitimate comparisons are limited to **algorithmic transfer**: within each
panel, did one ordering policy find that panel's operational event sooner than
another policy? It is also legitimate to ask whether the K=9 rank-deadline
property held, because that property concerns order positions rather than
biochemical scale.

The following cross-panel claims are not legitimate and are not made:

- equating a percent-inhibition event with a `Kd` event;
- comparing absolute AUDC magnitudes as biochemical efficacy or selectivity;
- inferring potency, intact-cell target engagement, toxicity, therapeutic
  window, or clinical benefit;
- attributing a sign difference solely to biology rather than platform,
  chemistry, panel composition, labeling, and reference-support differences.

“Profile-activity discovery” is defensible. “Potency discovery,” “safety
liability,” and “clinical selectivity” are not. The Klaeger label should be
called an **operational relative-affinity event**, not a universal binding
liability.

## 3. Operational activity labels

The `80` and `90` cutoffs have historical experimental context but no unique
biological status. Drewry et al. used 90% inhibition in PKIS selection, and
used PKIS2 single-concentration hits above 80% for `Kd` follow-up among a
selected subset. The same article explicitly calls inhibition-threshold
selection arbitrary and warns against direct selectivity-index comparison
across differently sized panels. Because PKIS2 was run in singlicate, a
single-concentration event requires concentration-response or `Kd`
confirmation before any potency language.

The frozen code uses strict `>` although the source paper uses both verbal
`>` criteria and a `>=90%` SI definition. This convention is disclosed. PKIS1
has only one retained matrix value exactly 80 and none exactly 90. PKIS2 has
275 raw values exactly 80 and 169 exactly 90, so an inclusive-boundary audit
was necessary. With `>=90`, the chemotype effect remains positive at both
parent (+0.0353, 0.0245 to 0.0465) and construct (+0.0367, 0.0261 to 0.0478)
resolution; standard effects remain near zero. At `>=80`, standard effects
also remain near zero. Thus the qualitative PKIS2 conclusion survives, while
the effect magnitude is correctly treated as threshold-convention dependent.

The Klaeger tenfold and half-log labels are study-defined transformations of
each compound's best measured documented-target `Kd`. They are useful for a
retrospective first-event benchmark but are not established safety margins.
Any stronger interpretation would require assay-specific replication,
concentration-response measurements, appropriate controls, a prespecified
construct, and context-relevant cellular target-engagement evidence.

## 4. Case-level sanity audit

The case audit is exhaustive under fixed flags rather than a favorable-case
sample. `external_flagged_cases.csv` contains every one of the 35
unconstrained delays of at least ten assays, every case-condition row with at
least 20 observed positive targets, and all bounded AUDC changes of magnitude
at least 0.30. Its 63 rows link to 288 underlying first-hit assay records in
`flagged_case_first_hit_source_records.csv`.

- Sparse measurement is not driving the external result. Candidate counts are
  198--200: 1,276 of 1,424 case-condition rows have 200 candidates, 140 have
  199, and eight have 198.
- Broad within-panel profiles can make first-hit discovery easy. At `>80`,
  CHEMBL1909396, CHEMBL1909349, and CHEMBL237347 have 30/200, 27/200, and
  26/199 positive parent targets. These are empirical panel descriptions, not
  general claims that the compounds are promiscuous in every assay context.
- Construct dependence is real. CHEMBL1909371 is `>90` for PDGFRA only through
  the archived T674I-mutant record (93.90; wild type 66.64). CHEMBL237347 is
  `>80` for KIT through D816V (84.47; wild type 48.21). These examples justify
  the any-construct wording and the construct-resolution analyses.
- Source-normalized percentage values can exceed 100. CHEMBL467581 has an
  archived INSR value of 107.4. It was retained without clipping as frozen;
  this does not change its binary threshold status, but the raw number must not
  be interpreted as a literal fraction.
- The largest repeated favorable example is CHEMBL66004 at `>90`, where
  bounded coverage advances the only positive parent (IKKE) by nine positions
  and gains 0.45 AUDC. Representative losses include CHEMBL1909362/Aurora-C
  and CHEMBL467581/INSR, each delayed nine positions with -0.45 AUDC. Both
  directions and every large unbounded delay are released.
- Source nomenclature includes `KIT (T6701 mutant)` and historical
  `BRAF (V599E mutant)` strings. They are preserved verbatim with assay IDs;
  no silent correction was made. Current ChEMBL calls CHEMBL5550 “Atypical
  kinase COQ8A, mitochondrial,” whereas the pinned ChEMBL 30 name is older.
  The executable filter relies on a ChEMBL protein-kinase classification path,
  not name substring matching, so these display strings do not alter target
  inclusion or ranking.
- Parent structures are canonicalized and duplicate source structures are
  median-collapsed; only one retained PKIS1 parent structure represents more
  than one source compound ID. Exact-structure references are excluded, and
  eight structures overlapping PKIS2 were removed before the external run.

No aggregate headline is explained away by missingness, duplicate chemistry,
display-name parsing, or one favorable compound. The case audit does expose a
real endpoint limitation: first-hit AUDC depends strongly on within-panel hit
prevalence and can represent mutant-only activity after parent collapse.

## Review decision

The paper is scientifically supportable as a **retrospective algorithmic study
of operational within-panel profile-activity events**, with the wording
changes above and with the negative PKIS1 transfer result retained. The review
does not support claims of prospective assay utility, biological validation,
potency, safety, or therapeutic recommendation.

An independent named kinase-assay scientist should still inspect this record
if the authors want to claim human expert review. That person should focus on
(1) whether the any-construct endpoint matches a realistic assay-ordering use
case, (2) whether a specific wild-type-only or domain-specific endpoint is more
appropriate, and (3) whether the relative-affinity benchmark is useful to an
actual screening team. Their future feedback may correct interpretation, but
any outcome or algorithm change must be labeled post hoc and validated anew.
