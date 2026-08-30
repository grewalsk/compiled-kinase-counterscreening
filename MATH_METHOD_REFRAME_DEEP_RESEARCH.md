# Mathematical Reframe for a Constructive AI-for-Drug-Discovery Paper

**Research decision — 30 August 2026**

## Recommendation

Replace the evaluation-paper thesis with:

> **When Sequential Is Not Adaptive: Compiled Submodular Kinase Counter-Screening**

The method would take a compound, a weighted collection of historical kinase-liability profiles, assay costs, and a turnaround-time penalty. It would return:

1. an ordered, compound-specific kinase mini-panel;
2. an exact compilation of any stop-on-first-hit adaptive policy into a fixed order;
3. optimal batch boundaries trading assay use against feedback latency; and
4. provenance showing which historical liability profiles each selected kinase newly covers.

This is a concrete planning method rather than another benchmark audit. It can reuse almost all of the completed dataset, leakage controls, Colab infrastructure, and statistical code.

## The mathematical insight that changes the project

The present benchmark asks whether at least one liability has been found and permits stopping at the first positive assay. That objective cannot support a genuinely branching binary-feedback policy.

Let \(Y_t\in\{0,1\}\) be the hidden result for kinase \(t\), and let policy \(\pi\) choose the next kinase from the observed history. Define a fixed order by simulating the policy under an all-negative history:

\[
a_1=\pi(\varnothing),\qquad
a_k=\pi((a_1,0),\ldots,(a_{k-1},0)).
\]

For every possible compound profile, the adaptive policy and this fixed order query exactly the same targets until the first positive or budget exhaustion. Before a hit, every observed outcome must be negative; after a hit, the first-hit task has already terminated.

Therefore the two procedures have identical, realization by realization:

- first-hit time;
- Hit@1/3/5/10/20;
- any-hit discovery curve and AUDC;
- sequential assay cost; and
- probability of finding at least one liability.

No conditional-independence assumption is needed. Kinase outcomes may be arbitrarily correlated, and the posterior may be wrong. For a randomized policy, fix its random seed before compilation and the same statement holds conditionally and hence in distribution.

This explains an important weakness in the current comparison. The existing “static particle” baseline freezes the initial probability ranking. It is not the fixed order produced by rolling the adaptive policy forward under negative outcomes. Feedback can improve over a naive frozen marginal ranking while still adding no runtime branching value for first-hit discovery.

General stochastic-probing research studies much broader branching objectives and obtains nontrivial adaptivity gaps ([Gupta, Nagarajan, and Singla](https://arxiv.org/abs/1608.00673)). General policy compilation is also established ([Grześ, Poupart, and Hoey](https://cs.uwaterloo.ca/~ppoupart/publications/bfsc/GrzesADT2013.pdf)). We should claim neither field as new. Our proposed contribution is the exact first-hit reduction and the counter-screen method it enables.

## Method 1: negative-path policy compilation

The first algorithm is extremely small:

```text
COMPILE-FIRST-HIT(policy π, compound x, budget B)
    history ← empty
    order ← empty
    repeat B times
        target ← π(history)
        append target to order
        append (target, NEGATIVE) to history
    return order
```

We would apply it to the existing particle-myopic and rollout policies and verify exact prefix equivalence for every compound, policy, random seed, and budget. These are algebraic invariant tests, not statistical comparisons.

The theorem has clear limits. Genuine adaptivity can return when:

- screening continues after a positive;
- utility rewards multiple liabilities or kinase-family coverage;
- the policy observes graded affinity rather than a binary threshold;
- multiple distinct nonterminal outcomes can change the next action;
- assays are noisy and repeated; or
- stopping depends on evidence accumulated across several assays.

Those limits make the paper more useful: it tells researchers which assay-planning objectives actually justify MCTS, POMDPs, or other sequential planners.

## Method 2: weighted min-sum set cover

After compilation, the drug-discovery problem becomes an ordering problem.

For query compound \(x\), let historical profile \(j\) have:

- a liability set \(L_j\) containing its measured non-mechanism kinase liabilities; and
- a chemical-support weight \(w_j(x)\ge0\), with weights summing to one.

For target order \(\sigma\), let \(\tau_j(\sigma)\) be the position of the first selected target that belongs to \(L_j\), or \(B+1\) if no selected target covers that profile. Optimize

\[
C_B(\sigma\mid x)=\sum_j w_j(x)\tau_j(\sigma).
\]

This is a weighted, truncated **min-sum set-cover** problem: targets are sets, and plausible compound profiles are the elements that should be covered early.

The connection to the current outcome is exact. If \(\tau\) is first-hit time censored to \(B+1\),

\[
\mathrm{AUDC}_{1:B}
=\frac{1}{B}\sum_{b=1}^B\mathbf 1[\tau\le b]
=\frac{B+1-\tau}{B}.
\]

Minimizing expected cover time is therefore exactly equivalent to maximizing expected any-hit AUDC.

The natural greedy rule selects

\[
\arg\max_t
\sum_{j:L_j\text{ not yet covered}}
w_j(x)\mathbf 1[t\in L_j].
\]

In plain language: choose the kinase that detects the largest remaining probability mass of plausible liability profiles. The selected target is interpretable because the system can list the analogue profiles it newly covers.

Classical work proves a factor-4 guarantee for natural min-sum set-cover greedy under the standard formulation ([Feige, Lovász, and Tetali](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/tr-2003-21.pdf)). For a fixed budget, the probability-of-any-hit objective

\[
F_x(S)=\sum_jw_j(x)\mathbf 1[S\cap L_j\ne\varnothing]
\]

is monotone submodular for any joint target distribution; cardinality-greedy has the classical \(1-1/e\) guarantee ([Nemhauser, Wolsey, and Fisher](https://thibaut.horel.org/submodularity/papers/nemhauser1978.pdf)). We must check the exact theorem statement for our weighted, censored, and potentially nonuniform-cost formulation before putting a guarantee in the paper.

## Method 3: latency-aware batch optimization

Labs commonly run assays in panels or plates rather than waiting for each individual result. After finding the order, choose how much feedback is operationally worth waiting for.

Let:

- \(c_k\) be the cost of assay \(a_k\);
- \(\lambda\) be a per-batch setup or turnaround penalty; and
- \(q_i=P(Y_{a_1}=\cdots=Y_{a_i}=0\mid x)\), the probability of reaching position \(i+1\).

For batch endpoints \(0=i_0<i_1<\cdots<i_R=B\), expected cost is

\[
\sum_{r=1}^R q_{i_{r-1}}
\left(\lambda+\sum_{k=i_{r-1}+1}^{i_r}c_k\right).
\]

For a fixed order, the optimal partition is found exactly in \(O(B^2)\):

\[
D(j)=\min_{0\le i<j}
\left[D(i)+q_i\left(\lambda+\sum_{k=i+1}^{j}c_k\right)\right].
\]

This returns an actionable schedule: for example, “run four assays now; if all are negative, run the next seven together.” Sequential batch testing already has a substantial operations-research literature, including dynamic programs with fixed batch and test costs ([Al-Turki et al.](https://www.sciencedirect.com/science/article/pii/S0957417424017251)). Our contribution would be its integration with a correlated, compound-specific kinase-profile posterior, not the general batching problem.

## Why the other mathematical directions are weaker

### Information theory

A total-correlation or conditional-mutual-information certificate could bound how much one assay teaches us about other assays. It is elegant, but it still makes the paper principally about evaluating whether feedback is useful. It is also statistically fragile with 111 high-dimensional profiles. Retain it only as a secondary diagnostic once the task allows more than one continuing outcome.

### Graph signal processing

Modeling a kinase profile as a graph signal and choosing nodes by graph-sampling or D-optimal design would be mathematically rich. However, graph-regularized drug–target matrix completion is established ([DLGRMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6304580/)), as is active graph-signal sampling ([Lin et al.](https://arxiv.org/abs/1902.04265)). More importantly, the current kinase-taxonomy baseline added essentially nothing. We do not yet have evidence that this liability signal is smooth on a defensible kinase graph.

### Distributionally robust optimization

Optimizing worst-case or CVaR coverage across chemical-support strata is attractive under series shift, but robust and risk-averse submodular optimization are established and technically harder ([Staib, Wilder, and Jegelka](https://proceedings.mlr.press/v89/staib19a.html)). Add a robust sensitivity only after the core method works.

### Group theory

Group actions are natural for atom relabeling and randomized-SMILES invariance. They are not central here because the current Morgan fingerprints are already graph-invariant. Group theory would be more honest in a separate paper on whether representation-equivalent SMILES cause different predictions or assay actions.

## Novelty relative to the closest drug-discovery papers

- Bayesian active search and nonmyopic planning are established ([Garnett et al.](https://arxiv.org/abs/1206.6406); [Jiang et al.](https://proceedings.mlr.press/v70/jiang17d.html)). We exploit a special absorbing objective where runtime branching disappears.
- A recent assay-planning preprint uses similarity-weighted historical cases and ensemble MCTS ([Chen et al.](https://arxiv.org/abs/2601.14710)). Its setting is broader and includes continuous, nonterminal outcomes; we should not suggest that compilation applies to all of it.
- BATCHIE uses Bayesian and submodular sequential design for drug-combination screens ([Tansey et al.](https://www.nature.com/articles/s41467-024-55287-7)). It optimizes global model information, not first-liability cover time for a single compound.
- Focused kinase mini-panels have already been optimized to approximate whole-panel selectivity ([Sutherland et al.](https://pubmed.ncbi.nlm.nih.gov/29792797/)). Our method must be distinguished by compound-specific scenario weights, ordered first-hit utility, and batch scheduling.

A defensible novelty statement is:

> Within a bounded literature search, we did not locate a public drug-discovery study combining lossless compilation of stop-on-first-liability policies, compound-specific weighted min-sum cover over analogue profiles, and latency-aware assay batching.

That is a combination-level claim, not a claim that the underlying mathematics is new.

## Required experiments

### Data plan

- **Development:** retain Klaeger/ChEMBL 30, the frozen 111-compound benchmark.
- **Primary external validation:** PKIS2, with 645 compounds and public profiling data. Its source article is CC BY and makes the associated data available ([Drewry et al.](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0181585)). Analyze its single-concentration percent-displacement outcome separately from \(K_d\).
- **Secondary only if access is resolved:** Davis, with 72 inhibitors and 442 kinase assays ([Davis et al.](https://www.nature.com/articles/nbt.1990)). Check the original supplement's redistribution terms before bundling it.
- ChEMBL-derived data remain CC BY-SA 3.0 ([official ChEMBL licensing](https://chembl.github.io/chembl-licensing/)).

### Baselines

Compare random, prevalence, chemical-\(k\)NN, static initial-particle ranking, current adaptive particle and rollout, their compiled negative paths, weighted min-sum greedy, fixed-budget maximum-coverage greedy, a global mini-panel, and exact small-instance optimization.

### Primary checks

1. Adaptive-versus-compiled first-hit metrics must match exactly for every case.
2. Weighted min-sum greedy should reduce censored assays to first liability relative to chemical-\(k\)NN and the global panel.
3. The direction should survive scaffold and strong nearest-analogue exclusions.
4. The method should reproduce on PKIS2 under a frozen, assay-appropriate threshold.
5. Batch schedules should be reported over a dimensionless \(\lambda/c\) sensitivity grid; do not invent universal assay prices.
6. Multi-hit and graded-outcome counterexamples should demonstrate exactly where adaptivity becomes genuine.

Use compound-paired bootstrap intervals and sign-flip tests for performance contrasts. Do not attach p-values to compiler equivalence. Keep all \(K_d\), percent-displacement, and biochemical-inhibition analyses separate, and never convert missing measurements into negatives.

## Biological and validity gates

A kinase profiling expert must review:

- the liability and PKIS2 hit thresholds;
- intended-target exclusions;
- target and variant mappings;
- whether a panel interaction is legitimately called a counter-screen liability;
- whether proposed batch sizes correspond to plausible workflows; and
- every interpretation beyond the measured assay universe.

The paper must not infer cellular engagement, toxicity, clinical safety, or therapeutic value. A no-hit result is not a “safe” or “clear” compound.

An optimization expert should independently check the compiler proposition, the AUDC identity, the correct min-sum/maximum-coverage guarantee for our exact formulation, and the batching recurrence.

## Deadline and feasibility

The official workshop page gives a **5 September 2026, 11:59 PM AoE** deadline and explicitly welcomes resource-constrained, decision-aware experiment prioritization ([AI4DD call for papers](https://ai4dd-neurips2026.github.io/)). The method is CPU-only:

- compiler: \(O(B)\) policy calls;
- min-sum greedy: approximately \(O(Bmn)\);
- batch partition: \(O(B^2)\);
- API cost: USD 0.

The credible six-day sequence is:

1. Freeze the theorem, algorithm, data thresholds, and primary contrast.
2. Implement the compiler, invariant tests, greedy cover, and batching DP in the existing Colab.
3. Run and debug on Klaeger, then freeze the external protocol.
4. Run PKIS2 once as external validation.
5. Complete robustness, figures, expert checks, and the five-page rewrite.

If PKIS2 cannot be constructed cleanly with traceable licensing and assay semantics by September 2, preserve the existing frozen paper or submit this as a short theory/method paper. Do not create a heterogeneous last-minute dataset.

## Bottom line

This reframe gives the project a memorable method and a useful output:

> **Compile the all-negative path, optimize the resulting order as weighted min-sum cover, and batch it according to assay turnaround cost.**

It also turns the current negative result into motivation rather than the contribution. The paper would no longer say merely that adaptivity failed. It would show why first-hit adaptivity is structurally unnecessary, provide a lossless replacement, and identify the richer objectives under which genuine adaptive drug-discovery planning is justified.
