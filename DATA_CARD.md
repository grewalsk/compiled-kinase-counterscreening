# Data card: compound-specific kinase counter-screen benchmark

## Intended use

Retrospective evaluation of low-compute target-acquisition policies that seek
measured binding liabilities in one internally consistent kinase panel. The
dataset is suitable for methods/benchmark research, not clinical or compound
progression decisions.

## Source and access

- Database release: ChEMBL 30 SQLite, released 2022-02-22.
- Source document: `CHEMBL3991601`, Klaeger et al., *The target landscape of
  clinical kinase drugs*, *Science* 2017.
- Source access: official downloadable ChEMBL release; the extraction opens
  it read-only and records the release metadata and byte size.
- License: ChEMBL data are CC BY-SA 3.0. The compact derived records therefore
  retain attribution and share-alike conditions. No commercial calculated
  ChEMBL property is extracted.

## Exact construction

The executable SQL and all filters are in `src/extract_panel.py`. The source
document is restricted to parent ChEMBL small molecules with a canonical
parent SMILES, human single-protein kinase targets, binding assays with
confidence 9, Kd in nM, relation `=` or `>`, no validity comment, and no
potential-duplicate flag. Missing measurements are never made negative.

A case is one parent compound. `D_i` contains all ChEMBL-documented direct
human-kinase mechanism targets represented in the source panel. `M_i` contains
the subset with an exact panel Kd at most 1,000 nM. A compound is eligible only
when `M_i` is nonempty; its reference affinity is the minimum Kd in `M_i`.

The common panel is the intersection of qualifying measured kinase targets
over every eligible compound. Every target in `D_i` is removed from compound
`i`'s candidate set, including documented targets too weak or censored to set
the reference affinity.

Primary label: exact candidate Kd no greater than ten times the reference
affinity. Sensitivity label: exact candidate Kd no greater than sqrt(10) times
the reference. Right-censored values are negative only when the censor bound
is above the liability boundary; the validator independently enforces this.

## Frozen census

- 49,603 qualifying source compound–target measurements.
- 242 source small molecules and 226 human kinase targets.
- 111 eligible parent-compound cases.
- 146 targets in the common panel.
- 16,034 non-mechanism candidate measurements.
- 63 compounds with at least one primary tenfold liability.
- 53 compounds with at least one half-log sensitivity liability.

These are a fixed census of the source under the rules above, not a prospective
sample-size target.

## Files and disclosure layers

- `cases.jsonl`: opaque case/target IDs, canonical parent structure, reference
  affinity, documented and anchor target sets, and candidate set.
- `gold.csv.gz`: hidden exact/censored assay interpretation and binary labels.
- `compound_mapping.csv`, `target_mapping.csv`: auditable mapping to ChEMBL.
- `mechanism_provenance.csv`: ChEMBL mechanism record IDs and annotations.
- `activity_provenance.csv.gz`: activity, assay, molecule, target, and source
  document IDs for every gold cell.
- `manifest.json`: source metadata, parameters, counts, bytes, and SHA-256.

Policies receive only the fields specified in the frozen protocol. Mapping and
gold files are used by the simulator/evaluator, not exposed as policy inputs.

## Known limitations

- One 2017 Kinobeads competition source dominates the benchmark; within-source
  comparability is strong, external validity is narrow.
- Binding in lysate does not establish functional inhibition or intact-cell
  engagement; later cellular profiling can differ.
- ChEMBL mechanism annotations are imperfect and time-bound to release 30.
- Clinical kinase drugs are not representative of arbitrary discovery-stage
  chemistry.
- Canonical SMILES and Bemis--Murcko scaffolds are representation choices.
- The relative tenfold label is operational, not a universal safety margin.
- Censored high values support nonliability only within the frozen boundaries.

## Validation

`src/validate_dataset.py` independently checks hashes, key opacity, common-panel
completeness, mechanism exclusion, both label calculations, censor safety,
activity uniqueness, and provenance coverage. Any failed assertion aborts the
Colab study.

