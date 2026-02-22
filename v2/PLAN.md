# SSA v2: Cancer Driver Gene Discovery with Structural Priors and Ablation-Validated Evaluation

**Working name**: SSA v2 (final name TBD — candidates: SSA-RANK, STRIDE, SSA-SP)

## Thesis

> **Cancer driver gene prediction faces five open challenges: PPI study bias,
> gold standard circularity, non-binary driver spectrum, mutual exclusivity
> confounding, and untested graph contributions. We address all five with a
> framework that combines structure-predicted protein networks (AlphaFold),
> positive-unlabeled ranking, biological priors as constraints, pluggable
> subnetwork validation, and bias-stratified evaluation against functional
> genomics data (DepMap).**

---

## Open Challenges in Cancer Driver Prediction

### Challenge 1: PPI network study bias
- Literature-curated PPIs (STRING, BioGRID) over-represent well-studied genes
- Well-studied genes have higher degree → GNNs preferentially select them
- [EMGNN reported 29-63% variance](https://academic.oup.com/bioinformatics/article/39/11/btad643/7325352)
  across different PPI networks — if PPI were robust signal, this wouldn't happen
- [Brouard et al. (2024)](https://academic.oup.com/bib/article/25/2/bbad027/7606638):
  MLP without network matched GNN performance in transcriptomic prediction

**Our approach**: Use [AlphaFold-predicted interactome (Science 2025)](https://www.bakerlab.org/wp-content/uploads/2025/11/Zhang-Science-Predicting-protein-protein-interactions-in-the-human-proteome-3.pdf) —
17,849 interactions at 90% precision from 200M protein pairs, structurally derived,
**no literature bias**. Supplemented by [PIONEER oncoPPIs](https://www.nature.com/articles/s41587-024-02428-4)
(586 cancer-enriched structural interactions across 33 cancer types).

### Challenge 2: Gold standard circularity
- CGC is embedded in NCG → using both is [data leakage](https://www.pnas.org/doi/10.1073/pnas.1616440113)
- Models learn to distinguish driver GENES from non-driver GENES (easy task),
  not functional drivers (hard task) — [28% performance drop](https://academic.oup.com/bib/article/25/2/bbad519/7584784)
  when negatives come from driver genes
- Ascertainment bias: kinases and well-studied genes are over-represented

**Our approach**: Validate against [DepMap CRISPR essentiality](https://depmap.org/portal/)
(1865 cell lines, functional evidence) and [AACR GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)
(110k+ tumors). Time-split evaluation: train on pre-2020 CGC, validate on 2020-2025 additions.

### Challenge 3: Non-binary driver spectrum
- [Mini-driver model](https://www.nature.com/articles/nrc3999): many mutations are weak/conditional drivers
- [Latent drivers](https://pmc.ncbi.nlm.nih.gov/articles/PMC9941481/): 155 same-gene double mutations
  where individual components only become oncogenic when paired
- [Passengers accumulate](https://www.biorxiv.org/content/10.1101/2025.11.13.688369v1.full)
  to mimic driver-like survival effects

**Our approach**: Continuous driverness score via PU-ranking, not binary classification.

### Challenge 4: ME confounding
- [Subtype confounding](https://pmc.ncbi.nlm.nih.gov/articles/PMC8664367/): different subtypes have
  different mutations → apparent ME without functional relationship
- [Low power for rare genes](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003503):
  at 1% frequency, ME is statistically undetectable in ~1000 samples
- [Co-occurrence is also informative](https://www.biorxiv.org/content/10.1101/2024.05.01.592039v1.full):
  cooperative driver pairs (PIK3CA+CDH1) are penalized by ME-only models

**Our approach**: ME as a biological constraint (not loss), corrected for subtype + TMB,
paired with co-occurrence signal for cooperative modules.

### Challenge 5: Untested graph contribution
- XGBoost baselines with identical features are rarely included in GNN cancer driver papers
- [SGCD (2025)](https://academic.oup.com/bib/article/26/1/bbae691/7942242) showed replacing GNN
  message passing with MLP improves results on heterophilic PPIs

**Our approach**: Mandatory ablation ladder — each component must demonstrably improve
over the previous level, or it is dropped.

---

## Architecture: Three-Component Design (Search + Score + Learn)

### The Key Architectural Insight from v1

SSA-ME v1 has a brilliant separation that NO GNN method preserves:

```
Search (ant colony)  ←→  Evaluation (scoreSubnetwork)
    DECOUPLED               PLUGGABLE — ANY LOGIC
```

The `scoreSubnetwork` is a **strategy pattern** — a black-box interface where you can plug
in ANY validation logic. v1 already has 6+ implementations:

| Implementation | What it evaluates |
|---|---|
| `MutualExclusivityNetworkManager` | ME (Dendrix/HyperGeometric) |
| `FunSeq2NetworkManager` | FunSeq2 functional impact scores |
| `ExpressionNetworkManager` | Differential expression |
| `PancancerNetworkManager` | Pan-cancer scoring |
| `QtlEpistasisNetworkManager` | QTL epistasis |
| `TimeExpressionNM` | Time-series expression |

**This capability is lost in end-to-end GNNs.** In a GNN, search and evaluation collapse
into a single differentiable forward pass. Non-differentiable domain tools (AlphaFold
structure analysis, metabolic flux simulation, molecular dynamics) cannot participate
in gradient-based optimization.

**Our v2 preserves this power AND adds neural acceleration.**

### Three-Component Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  COMPONENT 1: NEURAL ENCODER (differentiable, fast)             │
│                                                                 │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────┐    │
│  │ Multi-omics  │ │ Foundation   │ │ Functional           │    │
│  │ (TCGA)       │ │ Embeddings   │ │ Genomics             │    │
│  │              │ │              │ │                      │    │
│  │ • Mutations  │ │ • GenePT     │ │ • DepMap CRISPR      │    │
│  │ • CNA        │ │   (lit-based,│ │   essentiality       │    │
│  │ • Methylation│ │   best in    │ │ • Evolutionary       │    │
│  │ • Expression │ │   benchmark) │ │   constraint         │    │
│  │ • Mut impact │ │ • or scGPT   │ │   (pLI, LOEUF)      │    │
│  │   (CADD etc) │ │   (expr-     │ │ • AlphaFold          │    │
│  │              │ │   based)     │ │   interface scores   │    │
│  └──────┬───────┘ └──────┬───────┘ └──────────┬───────────┘    │
│         └────────────────┼────────────────────┘                 │
│                          │                                      │
│         ┌────────────────▼────────────────┐                     │
│         │  Encoder (Tabular OR Graph)     │                     │
│         │                                 │                     │
│         │  Level 0: XGBoost (no graph)    │                     │
│         │  Level 1: Graph Transformer     │                     │
│         │    • Multi-view: AlphaFold PPI, │                     │
│         │      PIONEER, HI-II-14, Reactome│                     │
│         │    • Bimodal: omics + topology   │                     │
│         │    • PPR diffusion              │                     │
│         └────────────────┬────────────────┘                     │
│                          │                                      │
│              Per-gene driverness prior                           │
│              s_g ∈ [0, 1] (PU-ranking)                          │
│              + per-gene embedding h_g ∈ R^d                     │
│                                                                 │
│  Trained via: PU-ranking loss + biological prior regularizers   │
│  + (optionally) RL reward from Component 3                      │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           │  gene scores guide exploration
                           │  (= neural pheromone initialization)
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  COMPONENT 2: SUBNETWORK SEARCH (graph exploration)             │
│                                                                 │
│  Search strategy (interchangeable):                             │
│                                                                 │
│  ┌─────────────────────┐   ┌──────────────────────────┐        │
│  │ Ant Colony (v1)     │   │ Neural-Guided Search     │        │
│  │ • Pheromones init   │   │ • Beam search on graph   │        │
│  │   from Component 1  │   │ • Use h_g embeddings to  │        │
│  │ • Random walks      │   │   predict promising      │        │
│  │ • Evaporation       │   │   extensions             │        │
│  │ • Proven algorithm  │   │ • Faster convergence     │        │
│  └─────────┬───────────┘   └─────────────┬────────────┘        │
│            └──────────────┬──────────────┘                      │
│                           │                                     │
│              Candidate subnetworks                               │
│              {module_1, module_2, ..., module_k}                 │
│              Each: set of genes + interactions                   │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           │  candidate modules
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  COMPONENT 3: PLUGGABLE SUBNETWORK SCORER (arbitrary logic)     │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │  SubnetworkScorer Protocol (Python ABC)                   │  │
│  │                                                           │  │
│  │  def score(module: Subnetwork) -> float                   │  │
│  │  def explain(module: Subnetwork) -> Explanation           │  │
│  │                                                           │  │
│  │  Can be ANYTHING. Does NOT need to be differentiable.     │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  Built-in scorers:                                              │
│  ┌─────────────────────┐  ┌──────────────────────────────┐     │
│  │ MEScorer            │  │ StructuralDisruptionScorer   │     │
│  │ (v1 Dendrix/HyperG) │  │ • AlphaFold interface ΔΔG   │     │
│  │ + subtype-corrected │  │ • FoldX/Rosetta binding est. │     │
│  │ + co-occurrence     │  │ • Which interactions break?  │     │
│  └─────────────────────┘  └──────────────────────────────┘     │
│  ┌─────────────────────┐  ┌──────────────────────────────┐     │
│  │ MetabolicImpactScore│  │ FunctionalGenomicsScorer     │     │
│  │ • Flux balance      │  │ • DepMap essentiality        │     │
│  │ • COBRA simulation  │  │ • Synthetic lethality        │     │
│  │ • Pathway cascade   │  │ • Drug response (GDSC)       │     │
│  └─────────────────────┘  └──────────────────────────────┘     │
│  ┌─────────────────────┐  ┌──────────────────────────────┐     │
│  │ CompositeScorer     │  │ UserDefinedScorer            │     │
│  │ • Chain any scorers │  │ • YOUR custom logic          │     │
│  │ • Weighted sum      │  │ • Import external tools      │     │
│  │ • Min/max/vote      │  │ • Call APIs, simulations     │     │
│  └─────────────────────┘  └──────────────────────────────┘     │
│                                                                 │
│  Output: scored & ranked modules with explanations              │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           │  (optional) scores as RL reward
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  FEEDBACK LOOP (optional — connects Components 1 ↔ 3)          │
│                                                                 │
│  Mode A: One-shot (simplest)                                    │
│    Component 1 trains independently → feeds Component 2 →       │
│    Component 3 scores → final ranking. No feedback.             │
│                                                                 │
│  Mode B: Score-as-label                                         │
│    Pre-compute Component 3 scores for known modules →           │
│    Use as additional training signal for Component 1             │
│                                                                 │
│  Mode C: RL loop (most powerful)                                │
│    Component 1 policy → Component 2 explores → Component 3     │
│    rewards → REINFORCE/PPO updates Component 1                  │
│    (= "neural ant colony" — GNN replaces pheromones,            │
│     scoring function stays arbitrary)                           │
└─────────────────────────────────────────────────────────────────┘
```

### What This Architecture Enables

Existing GNN cancer driver tools (EMOGI, EMGNN, CGMega, SGCD, DGCL_DWA, DGHNN,
GRAFT, CancerGATE) use end-to-end differentiable pipelines. This gains gradient
efficiency but prevents incorporating non-differentiable domain knowledge.

**Our architecture integrates domain tools that require non-differentiable computation:**

| Scenario | What happens | End-to-end GNN? | Our approach |
|---|---|---|---|
| "Does G12D at KRAS disrupt RAF binding?" | AlphaFold + FoldX ΔΔG estimation | Requires differentiable proxy | Direct via StructuralDisruptionScorer |
| "If KRAS-RAF is disrupted, what happens to MAPK flux?" | COBRA metabolic simulation | Not feasible | Direct via MetabolicImpactScorer |
| "Is this 4-gene module synthetically lethal?" | DepMap CRISPR analysis per module | Requires feature engineering | Direct via FunctionalGenomicsScorer |
| "Score this module using a new metric from 2027" | User implements new scorer | Requires retraining | Plug in new scorer, no retraining |

**v1's `scoreSubnetwork` already supported this. v2 preserves it and adds neural acceleration.**

### The Pluggable Scorer Protocol

```python
from abc import ABC, abstractmethod
from dataclasses import dataclass

@dataclass
class Subnetwork:
    """A candidate module: genes + their interactions."""
    genes: list[str]                    # gene symbols
    interactions: list[tuple[str, str]] # gene pairs (edges)
    mutations: dict[str, list]          # gene → list of mutations in module
    samples: list[str]                  # samples where module is altered

@dataclass
class ScorerExplanation:
    """Structured explanation of WHY a module scored high/low."""
    score: float
    components: dict[str, float]        # sub-score breakdown
    narrative: str                       # human-readable explanation
    evidence: dict                       # supporting data

class SubnetworkScorer(ABC):
    """
    The v1 architectural power, preserved.

    ANY scoring logic goes here. Does NOT need to be differentiable.
    v1's MutualExclusivityNetworkManager.scoreSubnetwork → this.
    """

    @abstractmethod
    def score(self, module: Subnetwork) -> float:
        """Score a candidate subnetwork. Any logic. Any complexity."""
        ...

    @abstractmethod
    def explain(self, module: Subnetwork) -> ScorerExplanation:
        """Explain why this module scored the way it did."""
        ...


class MEScorer(SubnetworkScorer):
    """v1 behavior: mutual exclusivity scoring (Dendrix + HyperGeometric)."""

    def __init__(self, alteration_matrix, subtypes=None):
        self.A = alteration_matrix
        self.subtypes = subtypes

    def score(self, module: Subnetwork) -> float:
        # Exact v1 logic: subtype-corrected Dendrix score
        if self.subtypes is not None:
            return self._subtype_corrected_me(module)
        return self._dendrix_score(module)

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        me_pairs = self._pairwise_me(module)
        return ScorerExplanation(
            score=self.score(module),
            components={f"{g1}-{g2}": s for (g1, g2), s in me_pairs.items()},
            narrative=f"Module of {len(module.genes)} genes with ME pattern",
            evidence={"coverage": self._coverage(module)}
        )


class StructuralDisruptionScorer(SubnetworkScorer):
    """
    NEW: Score by structural impact of mutations on protein interactions.

    Not feasible as a differentiable GNN loss function.
    Example: "Does G12D at KRAS disrupt the RAS-RAF interface?"
    """

    def __init__(self, alphafold_interfaces, binding_predictor):
        self.interfaces = alphafold_interfaces   # PIONEER data
        self.predictor = binding_predictor        # FoldX, Rosetta, or DL model

    def score(self, module: Subnetwork) -> float:
        disruption_score = 0.0
        for gene in module.genes:
            for mut in module.mutations.get(gene, []):
                # Which interactions does this mutation affect?
                affected = self.interfaces.get_interfaces(gene, mut.position)
                for interface in affected:
                    # Estimate binding affinity change
                    ddg = self.predictor.predict_ddg(
                        gene, mut, interface.partner
                    )
                    # Large |ΔΔG| = more disruptive = higher driver signal
                    disruption_score += abs(ddg) * interface.confidence
        return disruption_score

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        disrupted = self._get_disrupted_interfaces(module)
        return ScorerExplanation(
            score=self.score(module),
            components={f"{d.gene}:{d.mutation}→{d.partner}": d.ddg
                       for d in disrupted},
            narrative=f"Disrupts {len(disrupted)} protein interfaces",
            evidence={"interfaces": disrupted}
        )


class MetabolicImpactScorer(SubnetworkScorer):
    """
    NEW: Score by downstream metabolic impact.

    Example: "If KRAS is constitutively active and STK11 is lost,
    what happens to mTOR signaling and glycolytic flux?"

    Uses COBRA/flux balance analysis or pathway simulation.
    """

    def __init__(self, metabolic_model, pathway_map):
        self.model = metabolic_model      # e.g., Recon3D via COBRApy
        self.pathways = pathway_map

    def score(self, module: Subnetwork) -> float:
        # 1. Map mutations to enzyme perturbations
        perturbations = self._mutations_to_enzyme_effects(module)

        # 2. Simulate metabolic impact (FBA with constraints)
        baseline_flux = self.model.optimize()
        self.model.apply_perturbations(perturbations)
        perturbed_flux = self.model.optimize()

        # 3. Score = magnitude of flux change in cancer-relevant pathways
        return self._flux_disruption_score(baseline_flux, perturbed_flux)


class CompositeScorer(SubnetworkScorer):
    """Chain multiple scorers with configurable aggregation."""

    def __init__(self, scorers: list[tuple[SubnetworkScorer, float]],
                 aggregation: str = "weighted_sum"):
        self.scorers = scorers  # [(scorer, weight), ...]
        self.aggregation = aggregation

    def score(self, module: Subnetwork) -> float:
        scores = [(s.score(module), w) for s, w in self.scorers]
        if self.aggregation == "weighted_sum":
            return sum(s * w for s, w in scores)
        elif self.aggregation == "min":
            return min(s for s, _ in scores)
        elif self.aggregation == "vote":
            threshold = np.median([s for s, _ in scores])
            return sum(1 for s, _ in scores if s > threshold)


# Usage: compose any combination
scorer = CompositeScorer([
    (MEScorer(alteration_matrix, subtypes=subtypes), 0.4),
    (StructuralDisruptionScorer(alphafold, foldx), 0.3),
    (MetabolicImpactScorer(recon3d, kegg), 0.2),
    (FunctionalGenomicsScorer(depmap), 0.1),
], aggregation="weighted_sum")
```

### Neural Encoder Details

The encoder (Component 1) uses the same ablation ladder as before:

```
Level 0: XGBoost/LightGBM (tabular, no graph) — mandatory baseline
Level 1: Graph Transformer (multi-view, bimodal) — only if L0 beaten
Level 2: + biological prior regularizers — only if they help
```

Multi-view networks for Level 1:
- AlphaFold PPI (Baker 2025) — UNBIASED (structure-predicted)
- PIONEER oncoPPIs (NatBiotech 2024) — cancer-specific structural
- HI-II-14 — unbiased (proteome-scale screen)
- Reactome pathways — curated but functional

Loss for Component 1:
```
Loss = L_PU_ranking
     + α(t) · L_ME_regularizer      (subtype-corrected)
     + β(t) · L_CoOcc_constraint    (cooperative modules)
     + γ(t) · L_link_prediction     (graph regularization)
     + δ(t) · L_contrastive         (cross-view + PPR)
     + ε(t) · L_scorer_reward       (optional: from Component 3 via RL)

     α, β, γ, δ, ε: dynamic weights (uncertainty-based, Kendall)
```

### Why AlphaFold-Predicted Networks Change Everything

| Network | Source | Bias | Coverage | Cancer-specific |
|---|---|---|---|---|
| STRING | Literature mining + experiments | **HIGH** (study bias) | High (known genes) | No |
| BioGRID | Literature curation | **HIGH** | Medium | No |
| HI-II-14 | Proteome-scale Y2H screen | **LOW** | Low (~14k interactions) | No |
| **Baker 2025** | **AlphaFold structure prediction** | **NONE** (computational) | **17,849 at 90% precision** | No |
| **PIONEER oncoPPIs** | **AlphaFold + somatic mutations** | **LOW** | **586 cancer interactions** | **YES** |
| Reactome | Expert pathway curation | Medium | Medium | Partially |

**The AlphaFold-predicted interactome solves the core PPI bias problem:**
1. It covers understudied genes equally (no literature bias)
2. It's structurally grounded (physical contact, not co-citation)
3. 3,631 interactions are completely novel (never seen experimentally)
4. PIONEER adds cancer-specific structural interfaces where mutations cluster

**Novel angle**: We can create a feature per gene: "fraction of mutations that fall on
AlphaFold-predicted protein-protein interaction interfaces." This is a structural
driver signal that no GNN paper has used. If a gene's mutations cluster at interaction
interfaces, it's more likely to be a driver that disrupts specific interactions.

---

## The ME/Co-occurrence Biological Prior: Refined Design

### Why ME is a constraint, not a loss

ME as a primary loss has three problems:
1. Subtype confounding creates false ME signal → trains on noise
2. Low power for rare genes → no gradient signal where we need it most
3. Penalizes legitimate co-occurrence → suppresses cooperative driver pairs

As a **constraint/regularizer** (small weight, corrected signal), ME provides
biological inductive bias without dominating training.

### Formulation: Dual-Signal Biological Prior

```python
def biological_prior_loss(
    scores,           # predicted driverness scores [G]
    A,                # alteration matrix [G x S]
    subtypes,         # sample subtype labels [S]
    pathway_pairs,    # gene pairs in same pathway
    cooperative_pairs # known cooperative driver pairs
):
    """
    Dual-signal: ME between pathways, co-occurrence within cooperative modules.
    Both subtype-corrected.
    """
    # 1. Subtype-corrected ME (between-pathway)
    me_loss = 0
    for subtype in subtypes.unique():
        mask = subtypes == subtype
        A_sub = A[:, mask]

        for g_i, g_j in pathway_pairs:
            # Per-subtype co-occurrence
            co_obs = (A_sub[g_i] * A_sub[g_j]).mean()
            # Expected under independence (DISCOVER-style)
            co_exp = A_sub[g_i].mean() * A_sub[g_j].mean()
            # Weighted by prediction confidence
            weight = scores[g_i] * scores[g_j]
            # Penalize excess co-occurrence (ME violation)
            me_loss += weight * F.relu(co_obs - co_exp)

    # 2. Co-occurrence reward (within cooperative modules)
    coop_loss = 0
    for g_i, g_j in cooperative_pairs:
        co_obs = (A[g_i] * A[g_j]).mean()
        co_exp = A[g_i].mean() * A[g_j].mean()
        weight = scores[g_i] * scores[g_j]
        # Reward genuine co-occurrence above expected
        coop_loss -= weight * F.relu(co_obs - co_exp)

    # 3. Interface mutation enrichment (AlphaFold-derived)
    # Genes with mutations at PPI interfaces should score higher
    interface_loss = -torch.mean(
        scores * interface_mutation_fraction  # encourage interface-mutated genes
    )

    return me_loss + coop_loss + interface_loss
```

### Cooperative pair sources
- PIK3CA + CDH1 (breast cancer)
- KRAS + STK11 (lung adenocarcinoma)
- TP53 + RB1 (small cell lung cancer)
- [BeWith](https://www.researchgate.net/publication/316598763) for systematic detection
- [DIALECT (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11071465/) for latent-variable ME/co-occurrence

---

## PU-Ranking: The Correct Problem Formulation

### Why not binary classification

We have:
- **Positives**: Known drivers (CGC, OncoKB) — ~700 genes
- **Unlabeled**: Everything else — ~19,000 genes (some are undiscovered drivers)
- **NO confirmed negatives**

Every method that uses "negative examples" (genes assumed to be passengers)
is making an assumption that could be wrong. This is a textbook
[Positive-Unlabeled (PU) learning](https://pubmed.ncbi.nlm.nih.gov/22923290/) problem.

### Formulation

```python
class PURankingLoss(nn.Module):
    """
    Combines PU learning with ranking objective.

    - PU: treats unlabeled as mixture of hidden positives + true negatives
    - Ranking: known drivers should score higher than unlabeled genes
    """
    def __init__(self, prior_pi=0.05):
        super().__init__()
        self.prior_pi = prior_pi  # estimated fraction of true drivers among unlabeled

    def forward(self, scores, labels):
        # labels: 1 = known driver, 0 = unlabeled (NOT negative)
        pos_mask = labels == 1
        unl_mask = labels == 0

        # PU risk estimator (unbiased)
        # R_PU = pi * R+(f) + max(0, R-(f) - pi * R-(f_on_positives))
        pos_risk = F.binary_cross_entropy(scores[pos_mask], torch.ones_like(scores[pos_mask]))
        unl_risk = F.binary_cross_entropy(scores[unl_mask], torch.zeros_like(scores[unl_mask]))
        neg_on_pos = F.binary_cross_entropy(scores[pos_mask], torch.zeros_like(scores[pos_mask]))

        pu_loss = self.prior_pi * pos_risk + F.relu(unl_risk - self.prior_pi * neg_on_pos)

        # Pairwise ranking: each positive should score higher than random unlabeled
        pos_scores = scores[pos_mask]
        unl_scores = scores[unl_mask]
        # Sample pairs for efficiency
        n_pairs = min(len(pos_scores) * 10, len(unl_scores))
        idx = torch.randint(0, len(unl_scores), (n_pairs,))
        rank_loss = F.margin_ranking_loss(
            pos_scores.repeat_interleave(n_pairs // len(pos_scores)),
            unl_scores[idx[:len(pos_scores.repeat_interleave(n_pairs // len(pos_scores)))]],
            torch.ones(len(pos_scores.repeat_interleave(n_pairs // len(pos_scores)))),
            margin=0.1
        )

        return pu_loss + rank_loss
```

---

## Phases

### Phase 1: Data Infrastructure

- [ ] **1.1** Python project setup (pyproject.toml, DDD structure, pytest)
- [ ] **1.2** TCGA data ingestion via GDC API
  - Somatic mutations (MAF files)
  - Copy Number Alterations (GISTIC2)
  - DNA methylation (450k/EPIC arrays)
  - Gene expression (RNA-Seq, TPM normalized)
  - **Molecular subtype labels** per sample (PAM50, etc.)
- [ ] **1.3** Binary alteration matrix (v1 compatible)
- [ ] **1.4** Network construction: UNBIASED networks
  - **AlphaFold-predicted PPI** (Baker Lab Science 2025, 17,849 interactions)
  - **PIONEER oncoPPIs** (Nature Biotechnology 2024, 586 cancer interactions)
  - **HI-II-14** (proteome-scale screen, unbiased)
  - **Reactome pathways** (curated functional relationships)
  - For comparison only: STRING, CPDB (literature-biased baselines)
  - Interface mutation annotation: per-gene fraction of mutations at predicted
    interaction interfaces (from PIONEER structural models)
- [ ] **1.5** Ground truth — anti-circular design
  - **Training positives**: CGC Tier1 (pre-2020 freeze)
  - **Held-out temporal validation**: CGC genes added 2020-2025
  - **Functional validation**: DepMap CRISPR essentiality (1865 cell lines)
  - **Independent cohort**: AACR GENIE (110k+ tumors)
  - **NO combined CGC+NCG** (data leakage)
  - **PU framing**: unlabeled genes are UNLABELED, not negative
- [ ] **1.6** DepMap integration
  - Gene dependency scores (CERES) per cell line
  - Aggregate essentiality score per gene across cancer-relevant cell lines
  - Use as: (a) feature, (b) independent validation metric
- [ ] **1.7** Foundation model embeddings
  - GenePT embeddings (literature-based, [best in 2025 benchmark](https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2))
  - scGPT embeddings (expression-based, alternative)
  - Evaluate which adds more value in ablation

### Phase 2: Feature Engineering

- [ ] **2.1** Mutation features per gene
  - Frequency across samples
  - Type distribution (missense, nonsense, frameshift, splice)
  - Functional impact (CADD, PolyPhen-2, SIFT aggregated)
  - Hotspot clustering (OncodriveCLUSTL)
  - Background mutation rate correction (gene length, replication timing, expression level)
- [ ] **2.2** CNA features
  - Amplification/deletion frequency
  - Focal vs arm-level amplitude
  - CNA-expression correlation
- [ ] **2.3** Methylation features
  - Mean promoter methylation
  - Variance across tumors
  - Methylation-expression silencing score
- [ ] **2.4** Expression features
  - Differential expression tumor vs normal
  - Variance across tumors
  - Co-expression with network neighbors (WGCNA-style)
- [ ] **2.5** Functional genomics features (underexplored in GNN approaches)
  - DepMap CRISPR essentiality score (mean across cancer cell lines)
  - DepMap selective dependency (cancer-specific vs pan-essential)
  - Evolutionary constraint (pLI, LOEUF from gnomAD)
- [ ] **2.6** Structural features (NEW — AlphaFold-derived)
  - Fraction of mutations at predicted PPI interfaces (from PIONEER)
  - Interface disruption score: do mutations break specific interactions?
  - Number of predicted interaction partners (AlphaFold degree)
  - Compare: AlphaFold degree vs STRING degree (study bias test)
- [ ] **2.7** Foundation model embeddings
  - GenePT: d-dimensional embedding from NCBI gene descriptions via GPT
  - Captures functional relationships without network bias
- [ ] **2.8** Topology features (per network view)
  - Degree, betweenness, closeness (degree-corrected)
  - Node2Vec embeddings
  - Spectral positional encoding (Laplacian eigenvectors)
- [ ] **2.9** Bimodal separation (SGCD-inspired)
  - Omics features (2.1-2.5) + foundation embeddings (2.7): encoded via MLP_omics
  - Structural features (2.6) + topology (2.8): encoded via MLP_topo
  - Concatenated AFTER independent encoding

### Phase 3: Component 1 — Neural Encoder (Ablation Ladder)

The encoder is built as a stack of increasing complexity. Each level must
demonstrably improve over the previous, or it's dropped.

- [ ] **3.0** LEVEL 0: Tabular baseline (XGBoost/LightGBM)
  - All features from Phase 2 (no graph structure)
  - PU-aware training (treat unlabeled as unlabeled, not negative)
  - This is the floor. Everything else must beat this.
  - Report: AUPRC, precision@k, DepMap correlation
- [ ] **3.1** LEVEL 1: Graph model (only if L0 beaten)
  - Bimodal encoder (omics MLP + topology MLP)
  - Multi-view Graph Transformer (8 heads, 3 layers, residual + LayerNorm)
  - Views: AlphaFold PPI, PIONEER oncoPPIs, HI-II-14, Reactome
  - Shared-weight encoding + cross-view attention fusion
  - PPR diffusion + contrastive learning (DGCL_DWA-validated)
  - PU-ranking head (not binary classification)
- [ ] **3.2** LEVEL 2: Biological priors (only if they help)
  - Subtype-corrected ME regularizer (within pathways)
  - Co-occurrence constraint (cooperative modules)
  - Interface mutation enrichment (AlphaFold-derived)
  - Dynamic weight adjustment (Kendall uncertainty-based)
- [ ] **3.3** Ablation matrix (CRITICAL for paper credibility)

  | Config | Features | Graph | ME prior | CoOcc | Interface | PU-rank |
  |--------|----------|-------|----------|-------|-----------|---------|
  | L0     | All      | No    | No       | No    | No        | Yes     |
  | L1a    | Omics only | Yes | No      | No    | No        | Yes     |
  | L1b    | All      | Yes   | No       | No    | No        | Yes     |
  | L1c    | All      | Yes (STRING) | No | No  | No        | Yes     |
  | L1d    | All      | Yes (AlphaFold) | No | No | No     | Yes     |
  | L2a    | All      | Yes   | Yes (naive) | No | No      | Yes     |
  | L2b    | All      | Yes   | Yes (subtype-corr) | No | No | Yes  |
  | L2c    | All      | Yes   | Yes      | Yes   | No        | Yes     |
  | L2d    | All      | Yes   | Yes      | Yes   | Yes       | Yes     |
  | Binary | All      | Yes   | Yes      | Yes   | Yes       | No (BCE)|

  **Key comparisons:**
  - L0 vs L1b: Does the graph help at all?
  - L1c vs L1d: AlphaFold PPI vs STRING (bias test)
  - L2a vs L2b: Does subtype correction matter for ME?
  - L2c vs L2d: Does interface enrichment help?
  - L2d vs Binary: PU-ranking vs binary classification

### Phase 4: Component 2 — Subnetwork Search

- [ ] **4.1** Subnetwork search framework
  - `SearchStrategy` protocol: explores PPI graph, returns candidate modules
  - Module extraction: given gene scores from Component 1, find connected subnetworks
- [ ] **4.2** Ant colony search (v1 lineage)
  - Port v1's `NodeCostNetworkManager` exploration to Python
  - Initialize pheromones from Component 1's per-gene scores (neural pheromones)
  - Random walks + evaporation on multi-view networks
  - v1 compatibility: same exploration algorithm, better initialization
- [ ] **4.3** Neural-guided beam search (alternative)
  - Use per-gene embeddings from Component 1 to predict promising extensions
  - Beam search: start from seed gene, extend by highest-scoring neighbors
  - Faster convergence than ant colony for large networks
- [ ] **4.4** Hybrid search: ant colony warm-started by neural scores
  - Best of both: ant colony exploration diversity + neural score guidance
  - Compare all three strategies in ablation

### Phase 5: Component 3 — Pluggable Scorers

- [ ] **5.1** SubnetworkScorer protocol (Python ABC)
  - `score(module) → float` — can be any logic, non-differentiable OK
  - `explain(module) → ScorerExplanation` — structured explanation
  - `CompositeScorer` — chain multiple scorers
- [ ] **5.2** MEScorer (v1 behavior, modernized)
  - Dendrix score + HyperGeometric (exact v1 logic)
  - Subtype-corrected variant
  - Co-occurrence scoring for cooperative modules
- [ ] **5.3** StructuralDisruptionScorer (NEW)
  - Map mutations to AlphaFold-predicted PPI interfaces (PIONEER)
  - Estimate ΔΔG binding affinity change per mutation
  - Score = cumulative interface disruption in module
  - Answer: "Does this mutation break THIS specific interaction?"
- [ ] **5.4** MetabolicImpactScorer (NEW — the v1 power unleashed)
  - Map module perturbations to metabolic model (Recon3D via COBRApy)
  - Flux balance analysis: what happens to cancer-relevant pathways?
  - Score = magnitude of metabolic disruption
  - Example: KRAS activation + STK11 loss → mTOR/glycolysis impact
- [ ] **5.5** FunctionalGenomicsScorer
  - DepMap: is this module enriched for CRISPR-essential genes?
  - Synthetic lethality: do module genes have SL partners?
  - Drug response: GDSC sensitivity for module-altered cell lines
- [ ] **5.6** Feedback loop (optional)
  - Mode A (default): one-shot — encoder scores → search → scorer → done
  - Mode B: score-as-label — pre-compute scorer results, add as training signal
  - Mode C: RL — REINFORCE updates encoder using scorer rewards

### Phase 6: Interpretability

- [ ] **6.1** Multi-level explanations
  - **Gene level**: Integrated Gradients (Captum) from Component 1
  - **Module level**: SubnetworkScorer explanations from Component 3
  - **Network level**: which network view (AlphaFold vs pathway vs oncoPPI) is most
    informative per gene — enabled by cross-view attention
  - **Structural level** (NEW): for each predicted driver, which protein-protein
    interfaces are disrupted by its mutations? (from PIONEER + Component 3)
- [ ] **6.2** ME/Co-occurrence module discovery
  - For each predicted driver, extract ME partners (same pathway, exclusive alterations)
  - For each predicted driver, extract cooperative partners (co-altered)
  - Cross-reference with Reactome, KEGG pathway enrichment
  - Flag modules with therapeutic implications
- [ ] **6.3** Structural interpretation (enabled by AlphaFold networks)
  - Which mutations disrupt which specific protein-protein interactions?
  - Visualize: protein structure with mutated interface residues highlighted
  - Connect: disrupted interaction → pathway disruption → driver mechanism
  - This is not possible with literature-curated PPIs (no structural detail)
- [ ] **6.4** Visualization
  - Interactive network (Cytoscape.js, v1 compatible)
  - Oncoplot of ME + co-occurrence patterns per module
  - Sankey: gene → disrupted interface → pathway → drug target
  - Study bias analysis: prediction confidence vs gene "well-studiedness"

### Phase 7: Comprehensive Evaluation

- [ ] **7.1** Cross-validation
  - 5-fold stratified CV (PU-aware: stratify on positives only)
  - Leave-one-cancer-type-out CV
  - Leave-one-network-out CV (robustness to network choice)
- [ ] **7.2** Temporal validation (anti-circular)
  - Train: CGC pre-2020
  - Validate: CGC additions 2020-2025
  - Measure: how many new CGC genes were in our top-k predictions?
- [ ] **7.3** Functional validation (DepMap)
  - Correlation between predicted driverness score and DepMap essentiality
  - Are our top-k predictions enriched for CRISPR-essential genes?
  - Stratify by: known drivers, novel predictions, passengers
- [ ] **7.4** Independent cohort validation
  - AACR GENIE: are predicted drivers enriched in recurrently mutated genes
    across 110k+ tumors?
  - Hartwig Medical Foundation: do predictions hold in metastatic tumors?
- [ ] **7.5** Bias-stratified evaluation
  - Divide genes into quartiles by "well-studiedness" (PubMed citation count)
  - Report precision@k per quartile
  - Does the model discover drivers in the LEAST studied quartile?
  - Compare: AlphaFold network vs STRING (which is less biased in practice?)
- [ ] **7.6** Benchmark vs 15+ methods
  - **Frequency-based**: MutSigCV, MutSig2CV, OncodriveFML
  - **ME-based**: DISCOVER, DIALECT (2024)
  - **GNN methods**: EMOGI, MTGCN, EMGNN, CGMega, SGCD, DGCL_DWA,
    deepCDG, DGHNN, SEFGNN, GRAFT, CancerGATE
  - **Classical ML**: XGBoost (our own L0 baseline)
  - **Our own**: SSA-ME v1
  - **Metrics**: AUPRC, AUROC, precision@k, recall@k,
    DepMap essentiality correlation, temporal recall
  - **Stratified precision**: by mutation frequency (<1%, 1-2%, 2-5%, >5%)
  - **Study bias-adjusted metrics**: precision in understudied genes
- [ ] **7.7** Ablation study (Phase 3.3 table — the backbone of the paper)
- [ ] **7.8** Component architecture ablation (NEW)
  - Encoder-only (no search, no scorers) vs full pipeline
  - Search strategy comparison: ant colony vs beam search vs hybrid
  - Scorer comparison: ME-only vs structural vs composite
  - Feedback loop comparison: one-shot vs score-as-label vs RL

### Phase 8: Pan-Cancer and Discovery

- [ ] **8.1** Pan-cancer model
  - Cancer-type as learned embedding (not one-hot)
  - Shared encoder + cancer-specific heads
- [ ] **8.2** Transfer learning for rare cancers
  - Pre-train pan-cancer, fine-tune on small cohorts
  - Test: LAML (small) benefits from BRCA (large)?
- [ ] **8.3** Biological discovery
  - Novel candidate drivers: high score, NOT in any gold standard,
    BUT validated by DepMap essentiality
  - AlphaFold structural mechanism: which interfaces do they disrupt?
  - ME/co-occurrence modules: novel pathway relationships
  - Clinical translation: DGIdb druggability, CIViC actionability
- [ ] **8.4** Rare driver detection — capability analysis
  - Can the model actually find drivers at <2% frequency?
  - If yes: what features drive the prediction? (not ME — power is too low)
  - Hypothesis: DepMap essentiality + AlphaFold interface disruption + foundation
    embeddings may be the actual signals for rare drivers, not ME
  - ME helps for moderate-frequency (2-10%) drivers in shared pathways
  - Report per-frequency-band precision transparently

### Phase 9: Paper Writing & bioRxiv Publication (Paper 1)

- [ ] **9.1** Figures and tables
  - Figure 1: Architecture diagram (three-component design)
  - Figure 2: Ablation ladder results (L0 vs L1 vs L2, network comparison)
  - Figure 3: AlphaFold PPI vs STRING bias analysis (degree vs well-studiedness)
  - Figure 4: Phenotype decomposition overview / structural mechanism examples
  - Figure 5: Pan-cancer results (Borda aggregation, cancer-type transfer)
  - Table 1: Benchmark vs 15+ methods (AUPRC, precision@k, DepMap correlation)
  - Table 2: Ablation matrix (full Phase 3.3 table with all configs)
  - Table 3: Novel driver predictions with DepMap validation + structural mechanisms
  - Supplementary: frequency-stratified precision, bias quartile analysis, scorer comparison
- [ ] **9.2** Manuscript draft
  - Abstract (structured: background, methods, results, conclusion)
  - Introduction: five open challenges + our contributions
  - Results: ablation-first (show what works and what doesn't)
  - Methods: three-component architecture + PU-ranking + biological priors
  - Discussion: honest limitations, where graph doesn't help, ME power analysis
  - Follow Nature Methods / Nature Communications formatting guidelines
- [ ] **9.3** Supplementary materials
  - Detailed ablation results (all configs from Phase 3.3)
  - Per-cancer-type results breakdown
  - Full gene rankings (downloadable as supplementary table)
  - Code availability statement + GitHub repo link
  - Reproducibility: exact commands to regenerate all results
- [ ] **9.4** Code release preparation
  - Clean up repository for public release
  - Add installation guide (pip install, conda environment)
  - Add tutorial notebook (end-to-end example on small dataset)
  - Add CLI documentation (ssame train / ssame predict / ssame evaluate)
  - License selection (MIT or Apache 2.0)
  - Zenodo DOI for code snapshot
- [ ] **9.5** bioRxiv preprint submission
  - Format manuscript as PDF (LaTeX or Word template)
  - Submit to bioRxiv under Bioinformatics or Genomics category
  - Share preprint link on social media / relevant communities
  - Submit to target journal (Nature Methods or Nature Communications)
- [ ] **9.6** Peer review response
  - Address reviewer comments
  - Run additional experiments if requested
  - Resubmit revised manuscript

### Phase 10: Paper 2 — Phenotype-Directed Subnetwork Decomposition

Full spec in `v2/docs/PAPER2_PHENOTYPE_DECOMPOSITION.md`.

- [ ] **10.1** Heterogeneous network assembly
  - Load OmniPath directed signaling network (pypath, 100K+ directed signed interactions)
  - Load Reactome biochemical reactions (enzyme-substrate)
  - Load miRTarBase miRNA-target interactions (500K+ validated pairs)
  - Load DoRothEA / ENCODE TF-target interactions
  - Load KEGG / HMDB / Recon3D metabolic network
  - Load LncBase lncRNA interactions
  - Load ENCODE / EWAS Atlas methylation data
  - Integrate all layers into PyG HeteroData with typed nodes and edges
- [ ] **10.2** Phenotype effector curation
  - Download MSigDB Hallmarks (50 gene sets)
  - Download GO Biological Process terms (apoptosis, cell migration, angiogenesis, etc.)
  - Curate terminal effector node sets for each hallmark phenotype (multi-layer: proteins, metabolites, miRNAs)
  - Define phenotype labels from TCGA clinical data (metastasis, survival, drug response)
- [ ] **10.3** Variant impact pipeline
  - Download precomputed CADD scores for all SNVs
  - Download precomputed SpliceAI scores
  - Download GTEx eQTL data (variant-expression associations)
  - Integrate FoldX or ESM2 for protein stability ΔΔG predictions
  - Integrate AlphaFold-Multimer for binding affinity changes
  - Load UniProt functional domain annotations
  - Load PhosphoSitePlus PTM site annotations
  - Load codon usage tables
  - Build unified `variant_impact(alteration) -> float` function across all mechanism layers
- [ ] **10.4** PhenotypeScorer implementation
  - Implement Component 1: pathway connectivity scoring (directed shortest paths on heterogeneous graph)
  - Implement Component 2: patient-specific variant disruption scoring (multi-layer impact)
  - Implement Component 3: phenotype association scoring (odds ratio / Fisher exact test)
  - Implement composite score: connectivity × disruption × association
  - Implement `explain()` method: extract narrative from scoring internals (no extra computation)
  - Implement all four scorer subtypes: InformationFlow, SignalingDisruption, ExpressionCoherence, MultiOmicsPerturbation
  - Wire into existing CompositeScorer
- [ ] **10.5** Heterogeneous ACO
  - Extend ACO to walk on heterogeneous graph (multi-type node selection)
  - Ants can traverse cross-layer edges (protein → miRNA → mRNA → protein)
  - Modules can contain nodes of any type (proteins, metabolites, miRNAs, methylation sites)
  - Pheromone initialization from Component 1 adapted for heterogeneous nodes
- [ ] **10.6** Heterogeneous encoder
  - Extend Graph Transformer for cross-layer message passing (PyG HeteroData)
  - Different edge-type-specific attention weights
  - Node-type-specific feature encoders
  - Unified embedding space for cross-type scoring
- [ ] **10.7** Narrative generator
  - Template-based explain() that formats scoring internals into human-readable narratives
  - Per-patient narrative: which variants disrupt which cross-layer paths to which effectors
  - Per-module narrative: overall mechanistic summary of how the module drives the phenotype
  - Exportable as structured JSON + human-readable text
- [ ] **10.8** Validation suite
  - MSigDB overlap: predicted phenotype-gene associations recover known hallmark genes
  - Expression validation: invasion module genes correlate with EMT signatures
  - DepMap functional validation: CRISPR screens for specific phenotypes
  - Clinical validation: modules associate with patient outcomes (survival, drug response)
  - Cross-phenotype specificity: apoptosis module nodes NOT in invasion module
  - Variant-level validation: ClinVar pathogenic variants have high disruption scores
  - Cross-layer validation: recover known circuits (TP53 → miR-34a → BCL2, VHL → HIF1A → VEGFA)
  - Narrative coherence: auto-generated narratives consistent with published literature
- [ ] **10.9** Case studies & paper
  - BRCA phenotype decomposition: apoptosis, invasion, proliferation, angiogenesis, drug resistance
  - LUAD phenotype decomposition
  - COAD phenotype decomposition
  - Cross-cancer comparison: shared vs cancer-specific phenotype modules
  - Write Paper 2 manuscript following same structure as Paper 1
  - Submit to bioRxiv + target journal (Nature Methods / Genome Biology / Cell Systems)

---

## Project Structure

```
v2/
├── pyproject.toml
├── PLAN.md
│
├── ssame/
│   ├── __init__.py
│   │
│   ├── data/                           # Domain: data ingestion
│   │   ├── __init__.py
│   │   ├── tcga_client.py              # GDC API client
│   │   ├── alteration_matrix.py        # Binary alteration matrix
│   │   ├── omics_features.py           # Multi-omic features
│   │   ├── depmap_client.py            # DepMap essentiality data
│   │   ├── foundation_embeddings.py    # GenePT / scGPT embeddings
│   │   └── ground_truth.py             # CGC loader with temporal split
│   │
│   ├── network/                        # Domain: biological networks
│   │   ├── __init__.py
│   │   ├── alphafold_ppi.py            # Baker 2025 predicted PPI
│   │   ├── pioneer_oncoppi.py          # PIONEER cancer interfaces
│   │   ├── experimental_ppi.py         # HI-II-14 (unbiased screen)
│   │   ├── pathway_loader.py           # Reactome, KEGG
│   │   ├── biased_ppi.py              # STRING, CPDB (for comparison only)
│   │   ├── interface_features.py       # Per-gene interface mutation scores
│   │   └── multi_view_graph.py         # Multi-view graph construction
│   │
│   ├── encoder/                        # Domain: Component 1 (neural encoder)
│   │   ├── __init__.py
│   │   ├── tabular_baseline.py         # XGBoost/LightGBM (Level 0)
│   │   ├── bimodal_encoder.py          # Separate omics/topology MLPs
│   │   ├── graph_transformer.py        # Multi-view Graph Transformer
│   │   ├── cross_view_fusion.py        # Attention fusion across views
│   │   ├── pu_ranking_loss.py          # PU learning + ranking loss
│   │   ├── biological_priors.py        # ME + CoOcc + interface (differentiable)
│   │   ├── dynamic_weights.py          # Uncertainty-based loss balancing
│   │   ├── graph_diffusion.py          # PPR augmentation
│   │   └── encoder.py                  # Full encoder assembly
│   │
│   ├── search/                         # Domain: Component 2 (subnetwork search)
│   │   ├── __init__.py
│   │   ├── strategy.py                 # SearchStrategy protocol (ABC)
│   │   ├── ant_colony.py               # v1 ant colony (ported from Scala)
│   │   ├── beam_search.py              # Neural-guided beam search
│   │   ├── hybrid_search.py            # Ant colony + neural pheromones
│   │   └── subnetwork.py               # Subnetwork data structures
│   │
│   ├── scorer/                         # Domain: Component 3 (pluggable scorers)
│   │   ├── __init__.py
│   │   ├── protocol.py                 # SubnetworkScorer ABC + Subnetwork types
│   │   ├── me_scorer.py                # v1 ME scoring (Dendrix + HyperGeometric)
│   │   ├── structural_scorer.py        # AlphaFold interface disruption (ΔΔG)
│   │   ├── metabolic_scorer.py         # COBRA flux balance analysis
│   │   ├── functional_scorer.py        # DepMap essentiality + synthetic lethality
│   │   ├── composite_scorer.py         # Chain/combine multiple scorers
│   │   └── feedback.py                 # RL reward loop (REINFORCE/score-as-label)
│   │
│   ├── explain/                        # Domain: interpretability
│   │   ├── __init__.py
│   │   ├── integrated_gradients.py     # Per-feature attribution (Component 1)
│   │   ├── attention_explainer.py      # Multi-level attention analysis
│   │   ├── module_explainer.py         # Module-level explanations (Component 3)
│   │   ├── structural_explainer.py     # Interface disruption analysis
│   │   └── study_bias_auditor.py       # Well-studiedness analysis
│   │
│   ├── evaluation/                     # Domain: evaluation
│   │   ├── __init__.py
│   │   ├── pu_cross_validation.py      # PU-aware CV
│   │   ├── temporal_validation.py      # Time-split evaluation
│   │   ├── depmap_validation.py        # Functional validation
│   │   ├── independent_cohort.py       # GENIE / Hartwig validation
│   │   ├── study_bias_test.py          # Bias audit by gene studiedness
│   │   ├── ablation.py                 # Ablation ladder runner
│   │   ├── benchmark.py                # Multi-method comparison
│   │   └── rare_driver_analysis.py     # Frequency-stratified evaluation
│   │
│   ├── pipeline/                       # Domain: full pipeline orchestration
│   │   ├── __init__.py
│   │   ├── pipeline.py                 # Encoder → Search → Score → (Feedback)
│   │   └── config.py                   # Pipeline configuration
│   │
│   └── visualization/
│       ├── __init__.py
│       ├── network_viz.py              # Interactive network (Cytoscape.js)
│       ├── oncoplot.py                 # ME + co-occurrence patterns
│       ├── structure_viz.py            # Interface mutation visualization
│       ├── sankey.py                   # Gene → interface → pathway → drug
│       └── bias_plots.py              # Study bias analysis plots
│
├── tests/
│   ├── data/
│   ├── network/
│   ├── encoder/
│   │   ├── test_pu_ranking_loss.py
│   │   ├── test_biological_priors.py
│   │   ├── test_tabular_baseline.py
│   │   └── test_graph_transformer.py
│   ├── search/
│   │   ├── test_ant_colony.py          # v1 parity tests
│   │   ├── test_beam_search.py
│   │   └── test_hybrid_search.py
│   ├── scorer/
│   │   ├── test_me_scorer.py           # v1 scoring parity tests
│   │   ├── test_structural_scorer.py
│   │   ├── test_metabolic_scorer.py
│   │   ├── test_composite_scorer.py
│   │   └── test_scorer_protocol.py     # Any scorer satisfies protocol
│   ├── pipeline/
│   │   └── test_pipeline.py
│   ├── explain/
│   └── evaluation/
│       ├── test_temporal_validation.py
│       └── test_study_bias.py
│
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_network_bias_analysis.ipynb   # STRING vs AlphaFold bias comparison
│   ├── 03_feature_importance.ipynb
│   ├── 04_ablation_ladder.ipynb
│   ├── 05_me_confounding_analysis.ipynb # Subtype confounding in ME
│   ├── 06_search_strategy_comparison.ipynb  # Ant colony vs beam search
│   ├── 07_scorer_comparison.ipynb           # ME vs structural vs composite
│   ├── 08_depmap_validation.ipynb
│   ├── 09_structural_mechanisms.ipynb   # AlphaFold interface disruption
│   ├── 10_metabolic_impact.ipynb        # COBRA flux analysis of modules
│   ├── 11_rare_driver_analysis.ipynb     # Per-frequency rare driver analysis
│   └── 12_biological_discovery.ipynb
│
└── scripts/
    ├── train_encoder.py               # Train Component 1
    ├── run_search.py                  # Run Component 2
    ├── run_scorers.py                 # Run Component 3
    ├── run_pipeline.py                # Full pipeline: encode → search → score
    ├── ablation.py
    ├── benchmark.py
    └── visualize.py
```

---

## Technology Stack

| Component | Technology | Justification |
|---|---|---|
| **Tabular ML** | XGBoost + LightGBM | Mandatory baseline; often competitive with GNNs |
| **GNN Framework** | PyTorch Geometric (PyG) | Graph Transformer, multi-view, PPR diffusion |
| **Deep Learning** | PyTorch >=2.2 | Autograd for PU loss + biological priors |
| **Data** | pandas + polars | polars for large MAF/DepMap files |
| **Graphs** | NetworkX + PyG | Topology analysis + GPU tensors |
| **Statistics** | scipy + statsmodels | ME tests, FDR correction |
| **Interpretability** | Captum | Integrated Gradients |
| **Foundation Models** | transformers (HuggingFace) | GenePT / scGPT embeddings |
| **Structure** | BioPython + PyMOL (optional) | Interface visualization |
| **Metabolic Modeling** | COBRApy | Flux balance analysis for MetabolicImpactScorer |
| **Binding Affinity** | FoldX or DL predictor | ΔΔG estimation for StructuralDisruptionScorer |
| **Visualization** | matplotlib + plotly + Cytoscape.js | All plot types |
| **Experiment Tracking** | Weights & Biases (wandb) | Ablation management |
| **RL (optional)** | stable-baselines3 | REINFORCE for feedback loop (Mode C) |
| **Testing** | pytest | TDD |

---

## Success Metrics

### Performance (must demonstrate)
- [ ] Level 1 (graph) beats Level 0 (XGBoost) by >3% AUPRC
- [ ] AlphaFold PPI outperforms STRING on at least precision in understudied genes
- [ ] Biological priors (L2) improve rare driver precision over L1
- [ ] PU-ranking outperforms binary classification on temporal validation
- [ ] Full pipeline (encoder + search + scorer) beats encoder-only
- [ ] Composite scorer beats any single scorer

### Novelty (must show)
- [ ] **Pluggable scorer architecture**: first cancer driver method where subnetwork
  validation is decoupled from the learning algorithm — users can plug in ANY
  scoring logic (metabolic simulation, structural analysis, etc.) without retraining
- [ ] AlphaFold-predicted PPI network applied to cancer driver prediction
- [ ] PU-ranking for cancer driver prediction with biological priors
- [ ] Subtype-corrected ME as constraint (not primary loss) in GNN
- [ ] Bias-stratified evaluation with per-quartile metrics

### Biological discovery (must find)
- [ ] At least 10 novel drivers validated by DepMap but absent from CGC
- [ ] At least 3 structural mechanisms (which interface is disrupted)
- [ ] At least 2 novel ME modules with pathway coherence
- [ ] At least 1 therapeutic implication (druggable target)
- [ ] At least 1 metabolic mechanism (flux disruption) discovered via MetabolicImpactScorer

### Transparency (must report)
- [ ] Where does the graph NOT help? (ablation L0 vs L1)
- [ ] Where does ME NOT help? (ablation L2a vs L1)
- [ ] What is the actual power of ME for rare genes? (frequency-stratified analysis)
- [ ] How does prediction quality correlate with gene "well-studiedness"?
- [ ] Does the search+scorer pipeline add value over encoder-only? (7.8)

---

## Impact Strategy

### Publication target: Nature Methods or Nature Communications

**Why Nature Methods**: This paper proposes:
1. A new **evaluation framework** (temporal validation, bias-stratified metrics, DepMap functional validation) — reusable by the community
2. A new **method** (PU-ranking + biological priors on structure-predicted networks)
3. **Biological insight** (novel drivers with structural mechanisms)
Nature Methods values methodological rigor, comprehensive evaluation, and broad utility.

**Why Nature Communications**: RL-GenRisk (GCN + DQN for cancer risk genes) was
published here in 2025 with less methodological innovation. Realistic target.

### Paper structure (contributions in order of strength)

1. **The evaluation gap**: Five open challenges in the field, supported by
   quantitative evidence (bias-stratified metrics, temporal splits)
2. **The method**: PU-ranking on structure-predicted networks (AlphaFold) with
   biological priors as constraints — addresses challenges 1-4
3. **The ablation**: Each component's contribution, transparently — including
   graph vs no-graph, AlphaFold vs STRING, PU vs binary, scorer pipeline vs encoder-only
4. **The biology**: Novel drivers with structural AND metabolic mechanisms,
   validated by DepMap essentiality (not circular gene lists)
5. **The architecture**: Three-component design enabling pluggable subnetwork
   validation — presented as the vehicle, not the headline

### Three-paper strategy (if results support it)
1. **Main paper** (NMI/NatMethods): Method + architecture + evaluation + pan-cancer
2. **Resource paper** (Nucleic Acids Research): Web server / database of predictions
   with AlphaFold structural mechanisms per driver + pluggable scorer API
3. **Biological follow-up** (Cancer Discovery): Deep dive into novel rare drivers
   with structural + metabolic mechanisms and therapeutic implications

---

## v1 Compatibility

1. **Input compatible**: Accept same .m2 files and networks as v1
2. **Output compatible**: Generate same output files (edges, nodes, pattern, network.html)
3. **Direct benchmark**: Run v1 and v2 on same 12 TCGA datasets from the paper
4. **v1 wrapper**: Python script that calls v1 JAR for automated comparison
5. **Architectural lineage**: v1's `scoreSubnetwork` → v2's `SubnetworkScorer.score()`
   v1's `AntLikeOptimization` → v2's `ant_colony.py` in search/
   Same philosophy: search is decoupled from evaluation. v2 adds neural acceleration.

---

## References

### Challenges we address
1. Tokheim & Karchin (2016) Evaluating the evaluation. PNAS. doi:10.1073/pnas.1616440113
2. Brouard et al. (2024) Should we use GNNs for transcriptomic prediction? Brief. Bioinf. doi:10.1093/bib/bbad027
3. PPI study bias (2015). doi:10.1093/bioinformatics/btv036
4. Castro-Giner et al. (2015) Mini-driver model. Nat Rev Cancer. doi:10.1038/nrc3999
5. ME confounding. PMC8664367

### Our innovation sources
6. Baker Lab (2025) Human proteome PPI prediction. Science. doi:10.1126/science.adt1630
7. PIONEER (2024) Structurally informed interactome. Nat Biotechnol. doi:10.1038/s41587-024-02428-4
8. Predictomes (2025) AF-Multimer curated PPIs. Mol Cell. doi:10.1016/j.molcel.2025.01.034
9. DISCOVER (2016) ME correction. Genome Biology. doi:10.1186/s13059-016-1114-x
10. DIALECT (2024) Latent ME/co-occurrence. PMC11071465
11. MAGPIE (2024) Pathway-aware ME. Am J Hum Genet. doi:10.1016/j.ajhg.2023.12.009
12. GenePT (2023) Gene embeddings from GPT. PMC10614824
13. Gene embedding benchmark (2025). doi:10.1101/2025.01.29.635607

### GNN baselines
14. EMOGI (2021) Nat Machine Intelligence. doi:10.1038/s42256-021-00325-y
15. EMGNN (2023) Bioinformatics. doi:10.1093/bioinformatics/btad643
16. CGMega (2024) Nature Communications. doi:10.1038/s41467-024-50426-6
17. SGCD (2025) Briefings in Bioinformatics. doi:10.1093/bib/bbae691
18. DGCL_DWA (2025) JCST. doi:10.1007/s11390-025-4999-6
19. DGHNN (2025) Bioinformatics. doi:10.1093/bioinformatics/btaf379
20. GRAFT (2026) Briefings in Bioinformatics. doi:10.1093/bib/bbaf706
21. CancerGATE (2024) Comp Biol Med. doi:10.1016/j.compbiomed.2024.108367

### Validation resources
22. DepMap / Cancer Dependency Map. depmap.org
23. AACR Project GENIE. doi:10.1158/2159-8290.CD-22-0609
24. Hartwig Medical Foundation (2019) Nature. doi:10.1038/s41586-019-1689-y

### Our prior work
25. Pulido-Tamayo et al. (2016) SSA-ME. Scientific Reports. doi:10.1038/srep36257

### Foundational methods
26. Kendall et al. (2018) Multi-task uncertainty weighting. CVPR.
27. Klicpera et al. (2019) APPNP. ICLR.
28. Lin et al. (2017) Focal Loss. ICCV.
