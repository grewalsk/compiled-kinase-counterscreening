# Research background and novelty boundary

## Drug-discovery motivation

Large kinase profiling studies show that clinical and experimental kinase
inhibitors vary greatly in the breadth of their measured interactions. Karaman
et al. also showed that small panels do not robustly characterize selectivity.
The operational decision studied here is therefore not whether a molecule is
globally or clinically safe, but which measured counter-screen is most likely
to reveal a target-level binding liability before a limited panel is treated
as reassuring.

Primary biological sources:

- Fabian et al. (2005), *A small molecule-kinase interaction map for clinical
  kinase inhibitors*, doi:10.1038/nbt1068.
- Karaman et al. (2008), *A quantitative analysis of kinase inhibitor
  selectivity*, doi:10.1038/nbt1358.
- Klaeger et al. (2017), *The target landscape of clinical kinase drugs*,
  doi:10.1126/science.aan4368.
- Binder et al. (2026), *Cellular Context Influences Kinase Inhibitor
  Selectivity*, doi:10.1021/acs.jmedchem.5c02916.

The last paper is an important claim boundary: a lysate binding panel is not a
substitute for intact-cell target engagement or functional pharmacology.

## Methodological neighbors

- Garnett et al. (ICML 2012) formalized Bayesian active search, whose utility
  is the number of positives discovered under a query budget.
- Jiang et al. (ICML 2017) showed that nonmyopic active search can materially
  outperform myopic search and that exact optimization is hard.
- Mohamed et al. (RECOMB 2015) combined drug-target matrix factorization with
  active learning to reduce the number of experiments needed for accurate
  global interaction prediction.
- Irwin et al. (2021) evaluated deep imputation on sparse drug-discovery and
  public kinase-assay matrices. Related commercial work has discussed kinase
  assay prioritization.

## What is and is not novel

Not novel:

- kinase selectivity profiling;
- active learning or active search;
- compound-target prediction or matrix completion;
- chemical-similarity and global-prevalence baselines;
- the general idea of prioritizing experiments.

Combination-level gap sought here:

- one compound-specific counter-screen search episode per independently
  documented mechanism-target set;
- measured exact and censored negatives from a common public panel;
- explicit budgets and discovery-cost curves;
- unsafe-clear risk and abstention rather than prediction accuracy alone;
- adaptive posterior updates after each negative counter-screen;
- complete provenance and a Colab-scale reproducibility artifact.

The paper must say that no identical public benchmark was located within the
recorded search boundary. It must not claim exhaustive priority or a new
general active-search theory.

## Access and licensing

ChEMBL 30 is used under CC BY-SA 3.0 with attribution and share-alike
redistribution of derived records. RDKit is used for open-source molecular
fingerprints. The compact Colab bundle contains only the derived rows required
for this study; it excludes the full ChEMBL database and all credentials.

