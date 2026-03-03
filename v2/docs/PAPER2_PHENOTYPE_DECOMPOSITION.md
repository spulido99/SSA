# Paper 2: Phenotype-Directed Subnetwork Decomposition

**Working title**: "Phenotype-directed subnetwork decomposition reveals novel mechanistic drivers of cancer hallmarks"

**Prerequisite**: Paper 1 (SSA-ME v2) published or submitted.

**Last updated**: 2026-02-21 (post critical review)

---

## Critical Review Summary

### Novelty Assessment

The core concept — **"narrative-as-score"** (mechanistic explanation as optimization objective, not post-hoc) — has **no direct precedent in the literature** as of Feb 2026. Existing methods either find drivers then explain post-hoc, or find phenotype associations without mechanistic scoring during search.

### Competitive Landscape (2024-2026)

| Method | Multi-omics | Heterogeneous graph | Phenotype-specific | Patient-specific | Mechanistic explanation | Search |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| PhenoDriver (Briefings Bioinf 2023) | Partial | No | Yes | Yes | Partial (subnetworks) | No |
| CGMega (Nat Comm 2024) | Yes (3D genome, epigenome, PPI) | Partial | No | No | Modules + attention | No |
| GRAFT (Briefings Bioinf 2025) | Multi-view | Multi-view | No | No | Edge-attention | No |
| HGTDG (Methods 2024) | Yes (TCGA 16 types) | Yes | No | No | Limited | No |
| HallmarkGraph (Bioinformatics 2025) | Expression | GRN | Yes (hallmarks) | No | Classification only | No |
| SEFGNN (BIBM 2025) | Multi-view | Multi-view | No | No | DS uncertainty | No |
| Driver2Comm (PLoS CB 2026) | scRNA + spatial | No | Partial (comm) | No | Pathways | No |
| **This work** | **Yes (3+ layers)** | **Yes** | **Yes** | **Yes (variant)** | **Yes (narrative-as-score)** | **ACO** |

**Key differentiator**: No existing method combines all six columns. The unique combination of heterogeneous multi-layer graph + phenotype-specific scoring + patient-specific variant-level narratives + ACO search is the contribution.

### Primary Competitor: PhenoDriver

PhenoDriver (Briefings in Bioinformatics, 2023) identifies personalized phenotype-associated driver genes with subnetwork construction in breast cancer (988 patients). The framing overlaps significantly. **Explicit differentiation required in introduction:**

1. PhenoDriver uses homogeneous PPI — we use heterogeneous multi-layer graph
2. PhenoDriver works at gene level — we quantify at specific variant level (ΔΔG, SpliceAI, eQTL)
3. PhenoDriver scores post-hoc — our mechanistic quality IS the optimization objective
4. PhenoDriver has no search optimization — we use ACO guided by mechanistic score
5. PhenoDriver narratives are implicit — ours are explicit, quantified, causally reasoned

### Key Risks Identified

1. **Computational feasibility**: Pre-computation of all variant impacts is MANDATORY (not during ACO)
2. **Phenotype label circularity**: Molecular surrogates for phenotypes risk circular reasoning
3. **Effector set subjectivity**: Results depend on curated effector sets
4. **Cross-phenotype specificity**: Pleiotropic drivers (TP53) appear in multiple hallmarks — this is biologically correct
5. **No experimental validation**: Limits journal tier without wet-lab collaboration
6. **Scope ambition**: 6 molecular layers + 7 edge types + all hallmarks is too much for one paper

---

## Core Idea

Use SSA-ME v2's pluggable SubnetworkScorer architecture to create **phenotype-specific scorers** that guide the ACO search toward subnetworks that mechanistically explain specific cancer phenotypes (apoptosis resistance, motility, drug resistance, etc.).

The same cancer genome gets **decomposed** into phenotype-specific driver modules. Unknown genes appearing in these modules are predictions of novel phenotype-associated genes, and the subnetwork itself is the mechanistic explanation.

**Key insight**: The mechanistic narrative is NOT a post-hoc explanation — it IS the optimization objective. The PhenotypeScorer evaluates how well a candidate subnetwork explains the phenotype at the variant level for real patients. The ACO converges toward the subnetwork with the best mechanistic explanation.

## Key Difference from Paper 1

- **Paper 1** asks: "What subnetworks are drivers?"
- **Paper 2** asks: "What subnetworks mechanistically explain *this* phenotype?"

Paper 1 scores subnetworks by driver-ness (mutation recurrence, network centrality, PU learning).
Paper 2 scores subnetworks by **mechanistic explanation quality**: does this set of molecular entities, with the specific variants observed in patients, form a coherent causal chain to the phenotype across multiple molecular layers?

Different phenotype scorers produce different modules from the same graph. The ACO converges toward phenotype-relevant subnetworks instead of generic driver modules.

---

## Multi-Omics Heterogeneous Network

The search operates on a **heterogeneous multi-layer network** (not single-layer PPI). This is what enables the molecular-level explanations.

### Scope Decision (post-review)

**Phase 1 (this paper)**: Start with 3 molecular layers to prove the concept:
- **Proteins** (gene products) — core layer
- **miRNAs** (post-transcriptional regulators) — regulatory layer
- **Methylation sites** (CpG islands, promoter regions) — epigenetic layer

**Phase 2 (future work)**: Expand to full 6 layers:
- Metabolites (small molecules, lipids, amino acids)
- mRNAs / lncRNAs (transcripts)
- Chromatin regions (enhancers, super-enhancers)

**Rationale**: 3 layers are sufficient to demonstrate cross-layer mechanistic narratives (protein → miRNA → expression, methylation → silencing → protein loss) while keeping the system tractable. Each additional layer is an independent extension.

### Node Types (Phase 1)

- **Proteins** (gene products)
- **miRNAs / small RNAs** (post-transcriptional regulators)
- **Methylation sites** (CpG islands, promoter regions)

### Edge Types (Phase 1, biochemical relationships)

- Protein → Protein: phosphorylation, ubiquitination, complex formation (OmniPath, Reactome)
- miRNA → mRNA: post-transcriptional silencing (miRTarBase)
- Transcription factor → Gene: transcriptional activation/repression (DoRothEA)
- Methylation → Gene: epigenetic silencing (ENCODE, EWAS Atlas)

### Heterogeneous Graph Implementation

- PyG `HeteroData` natively supports multiple node/edge types
- Cross-layer message passing via heterogeneous Graph Transformer (see Encoder section)
- ACO walks on this heterogeneous graph — modules can span multiple molecular layers
- Phenotype effectors can be any node type (proteins, miRNAs, etc.)

### Edge Confidence Refinement with Interventional Data

**Improvement**: Use Perturb-seq (CRISPR perturbation + scRNA-seq) data to validate or refine edge directionality and confidence:
- If OmniPath says A → B (activation) but Perturb-seq shows knockout(A) has no effect on B → reduce edge confidence
- If Perturb-seq confirms knockout(A) → B drops significantly → reinforce edge with interventional evidence
- Source: Large-scale Perturb-seq datasets (K562, RPE1) are publicly available

This transforms the network from pure prior knowledge to **prior knowledge validated with interventional evidence**, strengthening Component 1 (Causal Pathway Reasoning) and all downstream narratives.

---

## Scoring Architecture: The Narrative IS the Score

### Evidence-Fused Composite Score (revised)

The original design used simple multiplication: `connectivity × disruption × association`. This is clean but **does not handle uncertainty**. A SpliceAI score of 0.6 and one of 0.99 contribute identically in a binary threshold scheme.

**Revised design**: Dempster-Shafer evidence fusion with calibrated uncertainty.

For a candidate subnetwork S and a target phenotype (e.g., "apoptosis resistance"):

#### Component 1: Causal Pathway Reasoning (CARNIVAL-based, static per phenotype)

**Improvement over original**: Replace shortest-path computation with CARNIVAL causal reasoning.

CARNIVAL (CAusal Reasoning for Network identification using Integer VALue programming) reconstructs the minimal causal subgraph connecting module nodes to phenotype effectors using Integer Linear Programming over OmniPath's signed, directed network. Benchmarked as the best algorithm for recovering signaling pathways containing compound targets (BMC Bioinformatics 2023).

```
Effectors_apoptosis = {
  Proteins: CASP3, CASP9, BAX, APAF1, BCL2, PUMA, ...
  miRNAs: miR-34a, miR-15a, ...
}

For each node Ni in S:
  causal_subgraph = CARNIVAL(
    source=Ni,
    targets=effectors,
    network=OmniPath_signed_directed,
    edge_confidence=interventional_validated  # Perturb-seq refined
  )
  # Reconstructs minimal causal paths respecting:
  # - Signed edges (activation vs. inhibition)
  # - AND/OR logic gates
  # - Cross-layer transitions (protein → miRNA → target)

score_connectivity = f(causal_subgraph_quality, effector_coverage, path_confidence)
```

**Why CARNIVAL over shortest-path**:
- Respects signaling logic: blocking an inhibitor = activation
- Produces **causally correct** narratives, not just topologically connected ones
- Minimizes the causal subgraph (Occam's razor for pathways)
- Already integrated with OmniPath (pypath → CARNIVAL pipeline exists)

**Pre-computation**: CARNIVAL solutions for all (gene, phenotype) pairs are computed ONCE and cached. During ACO, only lookup is needed.

#### Component 2: Patient-Specific Variant Disruption (Foundation Model pipeline, pre-computed)

**Improvement over original**: Unified foundation model pipeline replacing fragmented tools.

For each patient P who has alterations in nodes of S, evaluates whether those specific alterations actually disrupt the causal paths identified in Component 1.

**Tiered Variant Impact Pipeline** (ALL pre-computed before ACO):

```
Tier 1 — Fast screening (pre-computed for all possible SNVs):
  AlphaMissense pathogenicity score (71M missense variants pre-scored)
  CADD score (all possible SNVs pre-scored)
  → Binary filter: skip variants with AlphaMissense < 0.3 AND CADD < 15

Tier 2 — Mechanistic classification (pre-computed per observed variant):
  ESM2 fine-tuned (650M parameters):
    → Classifies 20 protein features at amino acid resolution
    → Infers WHY variant is pathogenic: loss of DNA_binding, loss of PTM_site,
      disruption of hydrophobic_core, loss of interface_contact, etc.
    → Output: feature_loss_vector with confidence per feature

Tier 3 — Structural quantification (pre-computed per missense in key genes):
  Rhapsody-2 (PNAS 2025):
    → Combines structure + evolution + dynamics
    → ΔΔG stability, ΔΔG binding, dynamic flexibility change
    → Mechanistic evaluation: "R248W loses DNA contact at position 248,
      ΔΔG_DNA_binding = +4.1 kcal/mol"

Tier 4 — Non-coding / regulatory (pre-computed):
  SpliceAI → splice disruption score
  GTEx eQTLs → expression effect size per tissue
  TCGA methylation → promoter methylation beta values
  TCGA CNA → segment-level copy number
  TCGA miRNA-seq → miRNA expression levels
```

```
For each patient P with alterations in S:
  For each altered node Ni in P:
    alteration = P.alterations[Ni]
    impact = tiered_variant_impact(alteration)  # Pre-computed lookup
    causal_importance = CARNIVAL_path_weight(Ni → effectors)

    disruption[P][Ni] = DS_evidence(
      belief=impact.pathogenicity,
      plausibility=impact.max_possible_effect,
      uncertainty=1 - impact.confidence
    )

  score_disruption(P) = DS_fuse(disruption[P])  # Dempster-Shafer fusion

score_variant_disruption = DS_fuse(score_disruption(P) for P with alterations in S)
```

**Variant impact across molecular layers (updated with foundation models):**

| Mechanism Layer | Effect | Tools / Data (open access) | Improvement |
|---|---|---|---|
| **Pathogenicity screening** | Binary pathogenic/benign | AlphaMissense (pre-computed, 71M variants) | **NEW**: replaces manual thresholds |
| **Mechanistic classification** | WHY variant is damaging | ESM2 fine-tuned (20 protein features) | **NEW**: feature-level explanation |
| **Protein stability** | Mutation destabilizes folding (ΔΔG) | Rhapsody-2 (structure + evolution + dynamics) | **UPGRADED**: from FoldX |
| **Binding affinity** | Mutation at interface reduces partner binding | Rhapsody-2 interface mode | **UPGRADED**: from AlphaFold-Multimer |
| **Catalytic function** | Mutation in active site abolishes enzymatic activity | UniProt active sites + ESM2 feature loss | **UPGRADED**: adds ESM2 |
| **PTM site disruption** | Mutation destroys phosphorylation/ubiquitination site | PhosphoSitePlus (open) + ESM2 PTM features | **UPGRADED**: adds ESM2 |
| **Splice disruption** | Variant alters splicing → exon skipping | SpliceAI (open, Illumina) | Unchanged |
| **mRNA stability** | UTR variant changes secondary structure | RBPmap, UTRdb | Unchanged |
| **Expression (eQTL)** | Promoter/enhancer variant reduces transcription | GTEx (open), ENCODE | Unchanged |
| **Methylation** | Promoter hypermethylation silences gene | TCGA methylation arrays | Unchanged |
| **Copy number** | Amplification/deletion alters dosage | TCGA CNA segments | Unchanged |
| **miRNA dysregulation** | miRNA over/underexpression | TCGA miRNA-seq | Unchanged |

#### Component 3: Phenotype Association (statistical, with circularity mitigation)

**Improvement over original**: Address circularity risk by prioritizing clinical endpoints.

Tests whether patients with disruptive alterations in S actually show the phenotype more than patients without.

```
Group A: patients with disruptive alterations in S
Group B: patients without alterations in S

score_association = odds_ratio(phenotype ~ group) or Fisher exact test
```

**Phenotype Label Strategy (revised to avoid circularity)**:

| Hallmark | Primary label (clinical, preferred) | Fallback label (molecular, use with caution) |
|---|---|---|
| Apoptosis resistance | Response to BH3 mimetics (if available), overall survival | Apoptosis pathway activity score (PROGENy) |
| Invasion/motility | Metastasis status (M0 vs M1, clinical data) | EMT score (signature-based) |
| Proliferation | Ki-67 index, mitotic count, tumor grade | Proliferation pathway activity (PROGENy) |
| Angiogenesis | VEGF serum levels (if available) | Hypoxia/angiogenesis signature |
| Drug resistance | Drug response data (TCGA clinical + GDSC) | ABC transporter expression |
| Immune evasion | Response to immunotherapy (if available) | Immune cell infiltration (ESTIMATE) |

**Circularity mitigation rules**:
1. ALWAYS prefer clinical endpoints (metastasis, drug response, survival) over molecular surrogates
2. When using molecular surrogates, ensure effector genes are EXCLUDED from the surrogate computation
3. Report which label type was used for each result — clinical labels get higher confidence
4. Sensitivity analysis: run with both clinical and molecular labels, compare results

This closes the causal loop: it is not enough for S to connect to effectors (Component 1) with damaging alterations (Component 2) — patients with disrupted S **must actually show the phenotype** (Component 3).

#### Composite Score: Dempster-Shafer Evidence Fusion

**Improvement over original**: Replace naive multiplication with uncertainty-aware evidence fusion.

```python
class PhenotypeScorer(SubnetworkScorer):
    def __init__(self, effector_nodes, heterogeneous_network, patients,
                 phenotype_labels, carnival_cache, variant_impact_cache, ...):
        self.effector_nodes = effector_nodes
        self.network = heterogeneous_network
        self.patients = patients
        self.phenotype_labels = phenotype_labels
        self.carnival_cache = carnival_cache          # Pre-computed CARNIVAL solutions
        self.variant_cache = variant_impact_cache      # Pre-computed variant impacts

    def score(self, module: Subnetwork) -> ScoredResult:
        # Each component returns (belief, plausibility, uncertainty)
        conn_evidence = self._causal_path_score(module, self.effector_nodes)
        disr_evidence = self._patient_disruption(module, self.patients)
        assoc_evidence = self._phenotype_association(module, self.phenotype_labels)

        # Dempster-Shafer fusion of three evidence sources
        fused = dempster_shafer_fuse([conn_evidence, disr_evidence, assoc_evidence])

        return ScoredResult(
            score=fused.belief,          # Point estimate for ACO pheromone
            confidence=fused.plausibility - fused.uncertainty,  # Confidence interval
            uncertainty=fused.uncertainty,
            evidence=fused                # Full evidence for explain()
        )

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        # Narrative extracted from scoring computation — no extra work
        # CARNIVAL causal subgraph provides causally correct paths
        # ESM2/Rhapsody-2 provide mechanistic variant impact details
        # DS fusion provides confidence intervals per claim
        # explain() formats the narrative with uncertainty bounds
        ...
```

**Why Dempster-Shafer over multiplication**:
- Each prediction comes with a **calibrated confidence interval**, not just a point score
- A module with score 0.7 ± 0.05 is qualitatively different from 0.7 ± 0.3
- Uncertainty propagates correctly: uncertain variant impacts → uncertain module scores
- Reviewers value uncertainty quantification (see SEFGNN, BIBM 2025)

**Conformal Prediction layer** (optional, for final predictions):
- After ACO converges, apply conformal prediction to generate calibrated prediction sets
- Output: "With 95% confidence, the apoptosis resistance drivers in this patient are {TP53, BAX, NOVEL_X}"
- Distribution-free, model-agnostic — works with any scorer combination

### Why This Works with ACO

The ACO explores candidate subnetworks on the heterogeneous graph. With the PhenotypeScorer:

| Candidate subnetwork | connectivity | disruption | association | Result |
|---|---|---|---|---|
| No causal paths to effectors | ≈ 0 | — | — | Pheromone evaporates |
| Connected but no patient alterations | high | ≈ 0 | — | Discarded |
| Connected + alterations but patients don't show phenotype | high | high | ≈ 0 | Discarded |
| Connected + disruptive alterations + patients show phenotype | high | high | high | **Pheromone reinforced** |

Ants converge toward the multi-layer subnetwork that **best mechanistically explains** why certain patients exhibit the phenotype.

---

## Phenotype Decomposition

The same dataset produces different modules with different scorers:

```
Same BRCA dataset, same ACO, same pipeline, same heterogeneous network:

ApoptosisScorer(effectors={CASP3, BAX, miR-34a, ...})
  → {TP53, MDM2, BAX, NOVEL_X, miR-34a↓}
  → confidence: 0.87 ± 0.04

InvasionScorer(effectors={MMP2, CDH1, miR-200, ...})
  → {CDH1, CTNNB1, NOVEL_Y, miR-200↓, DAPK1_methylated}
  → confidence: 0.79 ± 0.08

ProliferationScorer(effectors={CCND1, CDK4, RB1, ...})
  → {KRAS, BRAF, NF1, NOVEL_Z, let-7↓}
  → confidence: 0.82 ± 0.06

DrugResistanceScorer(effectors={ABCB1, ABCG2, ...})
  → {PIK3CA, PTEN, NOVEL_V}
  → confidence: 0.71 ± 0.12   ← higher uncertainty = less confident
```

Each module spans multiple molecular layers. Unknown nodes in a module are predictions of novel phenotype-associated entities, with a built-in mechanistic rationale (the causal paths and disruptions computed during scoring).

### Cross-Phenotype Specificity (revised)

**Original assumption (incorrect)**: "Apoptosis module nodes should NOT appear in invasion module."

**Revised understanding**: Pleiotropic drivers like TP53 legitimately participate in multiple hallmarks. Demanding total specificity is biologically wrong.

**Correct formulation**: Each module should have a **phenotype-specific core** of genes unique to that hallmark, plus an **expected overlap** in pleiotropic hubs. Metrics:
- Jaccard index between module pairs (expect < 0.3 for specific modules)
- Core specificity ratio: fraction of module nodes unique to this phenotype
- Pleiotropic hub identification: nodes appearing in 3+ modules are flagged as pleiotropic

---

## Example: Patient-Level Multi-Layer Narrative (updated with foundation models)

```
Module found by ACO: {TP53, MDM2, BAX, GENE_NOVEL_X, miR-34a, DAPK1_methylation}
Phenotype: apoptosis resistance
Confidence: 0.87 ± 0.04 (DS fusion of 3 evidence components)

Patient TCGA-A2-A0T3:

  TP53 R248W:
    AlphaMissense: 0.997 (pathogenic)
    ESM2 feature loss: DNA_binding=0.02, transactivation=0.05 (both lost)
    Rhapsody-2: ΔΔG_folding = +1.2 kcal/mol (mild destabilization)
                ΔΔG_DNA_binding = +4.1 kcal/mol (total loss of DNA contact)
    CARNIVAL causal path (signed, directed):
      TP53 --[activates]-→ BAX --[triggers]-→ CASP9 (BLOCKED: TP53 loss-of-function)
      TP53 --[activates]-→ miR-34a --[silences]-→ BCL2 (BLOCKED: TP53 loss-of-function)
    Disruption evidence: (belief=0.98, plausibility=0.99, uncertainty=0.01)

  MDM2 amplified (CN=8):
    Layer: Copy number → protein dosage
    → Overexpression → degrades residual wild-type p53
    → CARNIVAL path: MDM2 --[degrades]-→ TP53 (REINFORCES TP53 loss)
    Disruption evidence: (belief=0.90, plausibility=0.95, uncertainty=0.05)

  BAX promoter variant (-500 C>T, eQTL effect=-0.4):
    Layer: Regulatory → expression
    → GTEx: reduces BAX expression 40% in breast tissue
    → CARNIVAL path: BAX --[triggers]-→ CASP9 (WEAKENED)
    Disruption evidence: (belief=0.65, plausibility=0.80, uncertainty=0.15)

  miR-34a suppressed (by TP53 loss):
    Layer: miRNA → mRNA regulation
    → De-represses BCL2 (anti-apoptotic)
    → CARNIVAL path: miR-34a --[silences]-→ BCL2 mRNA (LOST → BCL2 overexpressed)
    Disruption evidence: (belief=0.85, plausibility=0.92, uncertainty=0.08)

  GENE_NOVEL_X splice variant (c.451+1G>A, SpliceAI=0.95):
    Layer: Splice → protein truncation
    → Exon 5 skipping → truncated protein, kinase domain lost
    → ESM2: loss of kinase_activity feature (score=0.03)
    → CARNIVAL path: NOVEL_X --[phosphorylates]-→ BAD --[releases]-→ BCL-XL (BLOCKED)
    Disruption evidence: (belief=0.88, plausibility=0.93, uncertainty=0.07)

  DAPK1 promoter hypermethylated:
    Layer: Methylation → gene silencing
    → TCGA methylation beta = 0.85 (heavily methylated)
    → DAPK1 silenced → parallel apoptosis trigger removed
    Disruption evidence: (belief=0.75, plausibility=0.88, uncertainty=0.12)

  Fused narrative (with confidence):
    "Apoptosis resistance in this tumor is a multi-hit, multi-layer event
    (overall confidence: 0.87 ± 0.04):

    TP53 R248W eliminates DNA binding (ΔΔG=+4.1 kcal/mol, ESM2 DNA_binding
    feature=0.02, AlphaMissense=0.997), removing transcriptional activation
    of both BAX and miR-34a [HIGH confidence: 0.98].

    Loss of miR-34a de-represses BCL2 at the post-transcriptional level
    [HIGH confidence: 0.85].

    MDM2 amplification (CN=8) degrades residual p53 [HIGH confidence: 0.90].

    The eQTL variant in BAX reduces available effector protein by ~40%
    [MODERATE confidence: 0.65 — eQTL effect is tissue-averaged].

    The splice variant in GENE_NOVEL_X (SpliceAI=0.95) truncates its kinase
    domain (ESM2 kinase_activity=0.03), blocking BAD phosphorylation → BCL-XL
    pathway [HIGH confidence: 0.88].

    DAPK1 promoter hypermethylation (beta=0.85) silences the last parallel
    apoptosis trigger [MODERATE-HIGH confidence: 0.75].

    Five independent mechanisms across three molecular layers converge on
    apoptosis resistance."
```

This narrative is not generated separately — it is extracted directly from the scoring computation. Each disruption was quantified during scoring with calibrated uncertainty; `explain()` formats it.

---

## Proposed Scorer Types

### 1. Causal Information Flow Scorer (CARNIVAL-based)

**Improvement over original "Information Flow Scorer"**: Uses CARNIVAL ILP instead of random walks.

- CARNIVAL reconstructs minimal causal subgraph from module nodes to phenotype effectors
- Respects signed, directed edges: activation vs. inhibition, AND/OR gates
- Edge confidence refined with Perturb-seq interventional data where available
- Module scores high if CARNIVAL finds short, high-confidence, causally coherent paths
- Pre-computed for all (gene, phenotype) pairs, cached

### 2. Foundation Model Disruption Scorer

**Improvement over original "Signaling Disruption Scorer"**: Uses unified foundation model pipeline.

- Tiered variant impact: AlphaMissense → ESM2 fine-tuned → Rhapsody-2 → SpliceAI/GTEx
- ESM2 provides MECHANISTIC classification: which protein feature is lost
- Rhapsody-2 provides STRUCTURAL quantification: ΔΔG values
- All pre-computed — ACO only does lookups
- Returns evidence triplet (belief, plausibility, uncertainty) per variant

### 3. Expression-Phenotype Coherence Scorer

- Requires multi-omics expression data (TCGA RNA-seq + miRNA-seq + methylation)
- Checks whether module nodes are consistently altered in tumors showing the phenotype
- Example: for "invasion/motility" scorer, module genes should be differentially expressed, relevant miRNAs dysregulated, and methylation sites altered in metastatic vs. primary tumors

### 4. Multi-Omics Perturbation Scorer

- Integrates mutation + expression + copy number + methylation + miRNA at variant level
- Scores modules by how consistently their specific multi-layer alterations correlate with the phenotype across patients
- Patient-specific: "in *this* tumor, apoptosis resistance is driven by TP53 R248W + MDM2 amp + miR-34a loss + BAX eQTL + DAPK1 methylation"
- Uses Dempster-Shafer fusion of per-variant evidence

These scorers can be composed via the existing `CompositeScorer` with learned or fixed weights.

---

## Drug Response Integration (NEW section)

**Rationale**: Connecting phenotype modules to drug response transforms the paper from "computational method" to "precision medicine tool." This is the difference between Genome Biology and Nature Communications.

### DrugResponseValidator

For each phenotype-specific module found by ACO:

```
1. Drug target mapping:
   Map module genes → known drug targets (DrugBank open, ChEMBL open)
   Example: PIK3CA in drug resistance module → PI3K inhibitors (alpelisib)

2. In silico drug response validation:
   Patients with disrupted module → response to drugs targeting that pathway?
   Source: GDSC (Genomics of Drug Sensitivity in Cancer, open)
   Source: TCGA clinical drug response data

3. Combination therapy rationale:
   If 3 mechanisms converge on apoptosis resistance:
     TP53 loss → BH3 mimetics unlikely to work alone
     MDM2 amp → MDM2 inhibitor (idasanutlin)
     BCL2 overexpression → venetoclax
     → Predicted combination: MDM2 inhibitor + venetoclax

4. Drug resistance prediction:
   DrugResistanceScorer modules identify which molecular mechanisms
   enable resistance BEFORE clinical failure
   → Proactive therapy switching
```

### Clinical Utility Metrics

- **Target actionability**: Fraction of module genes with approved or trial-stage drugs
- **Combination coverage**: Minimum drug combination to hit all convergent mechanisms
- **Resistance anticipation**: Early warning of resistance mechanisms from multi-omics profile

---

## Phenotype-Specific DepMap Validation (NEW section)

**Rationale**: DepMap CRISPR screens in hundreds of cell lines provide phenotype-specific functional validation WITHOUT requiring a wet lab. This should be CENTRAL to validation, not accessory.

### Strategy

```
For each phenotype module (e.g., apoptosis resistance):

1. Identify phenotype-relevant cell lines in DepMap:
   Apoptosis resistance:
     → Lines with high BCL2 expression, low BH3 mimetic response
     → Lines classified as "apoptosis-resistant" by drug sensitivity data

   Invasion/motility:
     → Lines with high EMT score, known metastatic origin
     → Lines with high migration/invasion assay data (CCLE annotations)

2. Test phenotype-specific essentiality:
   For each gene G in module:
     dependency_resistant = DepMap_score(G, resistant_lines)
     dependency_sensitive = DepMap_score(G, sensitive_lines)

     PASS if: G is MORE essential in phenotype-relevant lines
              (dependency_resistant < dependency_sensitive, p < 0.05)
     FAIL if: G is equally essential everywhere (generic essential)
     FAIL if: G is not essential anywhere

3. Novel gene validation:
   For predicted novel phenotype-associated genes:
     Check DepMap dependency in phenotype-relevant context
     High dependency + phenotype-specificity = strong functional evidence

4. Context Explorer (DepMap 25Q3):
   Use DepMap's Context Explorer to find enriched genetic dependencies
   across lineage and subtype contexts, cross-referencing with modules
```

### Expected Metrics

- Phenotype-specific essentiality enrichment (odds ratio, Fisher test)
- Novel gene validation rate (fraction of novel predictions confirmed by DepMap)
- Specificity: module genes should NOT be essential in phenotype-irrelevant lines

---

## Encoder Architecture (revised)

### Heterogeneous Graph Transformer (informed by GRAFT/HGTDG)

**Improvement**: Adopt edge-attention bias from GRAFT (Briefings Bioinf 2025) and heterogeneous message passing from HGTDG (Methods 2024).

```python
class PhenotypeHGT(torch.nn.Module):
    """
    Heterogeneous Graph Transformer for phenotype-directed encoding.
    Informed by GRAFT edge-attention bias + HGTDG multi-omics design.
    """
    def __init__(self, node_types, edge_types, hidden_dim, num_heads, num_layers):
        # Per-type linear projections (different feature spaces per molecular layer)
        self.node_projections = {nt: Linear(feat_dim[nt], hidden_dim) for nt in node_types}

        # Edge-attention bias: different attention patterns per edge type
        # (PPI edges vs miRNA-target edges vs TF-target edges have different semantics)
        self.edge_attention_bias = {et: Parameter(num_heads) for et in edge_types}

        # Structural encodings per molecular layer
        self.structural_encodings = {nt: StructuralEncoding(nt) for nt in node_types}

        # Cross-layer transformer layers
        self.layers = [HGTLayer(hidden_dim, num_heads, edge_types) for _ in range(num_layers)]

    def forward(self, hetero_data):
        # 1. Project each node type to shared dimension
        # 2. Add structural encodings
        # 3. Cross-layer message passing with edge-type-specific attention
        # 4. Output: node embeddings aware of cross-layer context
        ...
```

**Key design decisions**:
- Edge-attention bias differentiates PPI from miRNA-target from TF-target interactions
- Structural encodings capture layer-specific topology (degree, betweenness, etc.)
- Cross-layer message passing learns which cross-layer interactions are most informative
- Builds on proven architectures (GRAFT, HGTDG) rather than reinventing

---

## Data Requirements (all open access)

### Phenotype Definitions
- **MSigDB Hallmarks**: 50 gene sets (apoptosis, EMT, angiogenesis, etc.), freely downloadable
- **GO Biological Process**: Apoptosis, cell migration, angiogenesis, etc.
- **Phenotype effector nodes**: Curated from Hallmarks + GO (terminal effectors of each phenotype, any molecular layer)
- **Sensitivity analysis**: Vary effector sets (±20% genes) and report robustness

### Heterogeneous Network (Phase 1)

| Layer | Source | Access | Coverage |
|-------|--------|--------|----------|
| Directed signaling (PPI) | OmniPath | Open (pypath) | 100K+ interactions, directed, signed |
| Biochemical reactions | Reactome | Open | 14K+ reactions, enzyme-substrate |
| miRNA targets | miRTarBase | Open | 500K+ validated miRNA-target pairs |
| TF-target | DoRothEA | Open | 1M+ TF binding events |
| Methylation | ENCODE / EWAS Atlas | Open | Genome-wide CpG methylation |

### Variant Impact Tools (updated with foundation models)

**Tier 1 — Screening (pre-computed)**:
- **AlphaMissense**: Pathogenicity for 71M missense variants (open, DeepMind)
- **CADD scores**: Precomputed for all possible SNVs (open, download)

**Tier 2 — Mechanistic (pre-computed per observed variant)**:
- **ESM2 fine-tuned**: 20 protein features at amino acid resolution (open, Meta)
- **Rhapsody-2**: Structure + evolution + dynamics variant impact (open, PNAS 2025)

**Tier 3 — Regulatory/splicing (pre-computed)**:
- **SpliceAI**: Precomputed splice impact scores (open, Illumina)
- **GTEx eQTLs**: Variant-expression associations across tissues (open)

**Tier 4 — Annotations**:
- **UniProt domains**: Functional domain annotations (open)
- **PhosphoSitePlus**: PTM site annotations (open for non-commercial)

### Interventional Data (network refinement)
- **Perturb-seq K562**: Large-scale CRISPR + scRNA-seq (open, Replogle et al.)
- **Perturb-seq RPE1**: Additional cell type for cross-validation

### Patient Multi-Omics Data — Pan-Cancer Cohort (10 cancer types)

**Strategy**: Full pan-cancer analysis across the 10 TCGA cancer types with the most data and highest therapeutic relevance. This enables discovery of conserved vs. cancer-specific hallmark driver modules — a contribution no existing method has made.

| # | TCGA Code | Cancer Type | n patients | Selection rationale |
|:-:|:---------:|-------------|:----------:|---------------------|
| 1 | **BRCA** | Breast invasive carcinoma | 1,097 | Largest cohort; benchmark vs PhenoDriver; well-defined molecular subtypes (Luminal A/B, HER2+, Basal) |
| 2 | **KIRC** | Kidney renal clear cell carcinoma | 537 | VHL→HIF1A→VEGFA is the ideal case for multi-layer causal narrative; resistant to conventional chemo |
| 3 | **HNSC** | Head & neck squamous cell carcinoma | 528 | HPV+ vs HPV- are biologically distinct diseases within the same cohort — natural internal control |
| 4 | **LUAD** | Lung adenocarcinoma | 522 | #1 cancer mortality globally; many targeted therapies (EGFR, ALK, KRAS G12C) for drug response integration |
| 5 | **LUSC** | Lung squamous cell carcinoma | 504 | No effective targeted therapies — greatest clinical need for novel target discovery |
| 6 | **SKCM** | Skin cutaneous melanoma | 470 | Highest tumor mutational burden; paradigm for immunotherapy (checkpoint inhibitors) |
| 7 | **COAD** | Colon adenocarcinoma | 459 | MSI-high vs MSS dichotomy directly maps to immunotherapy response; well-studied hallmarks |
| 8 | **STAD** | Stomach adenocarcinoma | 443 | 4 molecular subtypes (EBV, MSI, GS, CIN) — maximum molecular diversity within one cohort |
| 9 | **BLCA** | Bladder urothelial carcinoma | 412 | Responds to checkpoint inhibitors; luminal vs basal subtypes affect treatment |
| 10 | **UCEC** | Uterine corpus endometrial carcinoma | 548 | 4 molecular subtypes (POLE ultramutated, MSI, CN-low, CN-high); POLE has exceptionally high TMB |

**Total**: ~5,520 patients across 10 cancer types.

**Excluded cancer types and rationale**:
- GBM/LGG (brain): Blood-brain barrier limits drug response integration; fundamentally different biology
- OV (ovarian): Heavily studied but few novel therapeutic options
- THCA (thyroid): Generally excellent prognosis; little need for phenotype decomposition
- PRAD (prostate): Many indolent cases make phenotype labeling ambiguous
- LIHC (liver): Viral etiology (HBV/HCV) confounds driver analysis

**Multi-omics data per patient** (all from GDC, open access):
- **Mutations**: MAF files (specific variants per patient)
- **RNA-seq expression**: Gene-level FPKM/TPM
- **miRNA-seq**: miRNA expression quantification
- **DNA methylation**: 450K/EPIC arrays (promoter CpG beta values)
- **Copy number**: Segment-level CNA from SNP arrays
- **Clinical data**: Stage, metastasis status, survival, drug response (where available)

**Pan-cancer analysis dimensions**:

```
For each of the 10 cancer types × each of 5 hallmark phenotypes:
  → Run ACO with PhenotypeScorer
  → Produce phenotype-specific module with confidence interval
  → Generate patient-level mechanistic narratives

Cross-cancer aggregation:
  1. Conserved hallmark cores: genes appearing in hallmark H
     across ≥7/10 cancer types → universal hallmark drivers
  2. Cancer-specific hallmark drivers: genes appearing only
     in 1-2 cancer types for hallmark H → tissue-specific mechanisms
  3. Pleiotropic hubs: genes appearing in ≥3 hallmarks across
     ≥5 cancer types → master regulators (e.g., TP53, PIK3CA)
  4. Cross-cancer drug target map: module genes × approved drugs
     × cancer types → pan-cancer actionability matrix
```

### Drug Response Data
- **GDSC**: Genomics of Drug Sensitivity in Cancer (open, Sanger)
- **DrugBank**: Drug-target mapping (open access portion)
- **ChEMBL**: Bioactivity data (open, EMBL-EBI)

### Functional Validation Data
- **DepMap 25Q3**: CRISPR screens in 1000+ cell lines (open, Broad)
- **ClinVar**: Known pathogenic variants for variant-level validation (open, NCBI)

---

## Architecture (builds on SSA-ME v2)

All phenotype scorers implement the existing `SubnetworkScorer` protocol:

```python
class PhenotypeScorer(SubnetworkScorer):
    def __init__(self, effector_nodes, heterogeneous_network, patients,
                 phenotype_labels, carnival_cache, variant_impact_cache, ...):
        ...
    def score(self, module: Subnetwork) -> ScoredResult:
        conn_evidence = self._causal_path_score(module, self.effector_nodes)   # CARNIVAL
        disr_evidence = self._patient_disruption(module, self.patients)         # Foundation models
        assoc_evidence = self._phenotype_association(module, self.phenotype_labels)
        fused = dempster_shafer_fuse([conn_evidence, disr_evidence, assoc_evidence])
        return ScoredResult(score=fused.belief, confidence=..., uncertainty=fused.uncertainty)
    def explain(self, module: Subnetwork) -> ScorerExplanation:
        # CARNIVAL causal subgraph → causally correct paths
        # ESM2/Rhapsody-2 → mechanistic variant details
        # DS fusion → confidence intervals per claim
        # explain() formats narrative with uncertainty bounds
        ...
```

**Changes from Paper 1:**
- **Network**: Heterogeneous multi-layer graph (not single-layer PPI)
- **ACO**: Walks on heterogeneous graph (modules can span molecular layers)
- **Encoder**: Heterogeneous Graph Transformer with edge-attention bias (GRAFT-informed)
- **Scorers**: New PhenotypeScorer family with DS fusion (causal × disruption × association)
- **Variant impact**: Foundation model pipeline (AlphaMissense → ESM2 → Rhapsody-2)
- **Pathway reasoning**: CARNIVAL ILP (not shortest-path)
- **Uncertainty**: Dempster-Shafer fusion + optional conformal prediction
- **Pipeline, feedback loop, CLI**: No changes

The `explain()` method is essentially free: all causal paths, variant impacts, and patient associations were already computed during `score()`. The narrative is just a formatted view of the scoring internals with calibrated confidence.

---

## Validation Strategy (revised)

### Tier 1: Recovery of Known Biology (sanity checks)

1. **Known associations**: Predicted phenotype-gene associations should recover known hallmark genes (MSigDB overlap). Expected: >70% recall of top hallmark genes.
2. **Known regulatory circuits**: Multi-layer modules should recover known circuits (e.g., TP53 → miR-34a → BCL2, VHL → HIF1A → VEGFA).
3. **Variant-level validation**: Known pathogenic variants (ClinVar) in module genes should have high disruption scores.

### Tier 2: Functional Validation (DepMap, no wet lab needed)

4. **Phenotype-specific essentiality**: Module genes should be preferentially essential in phenotype-relevant DepMap cell lines (see DepMap Validation section).
5. **Novel gene functional evidence**: Predicted novel phenotype genes should show dependency in relevant DepMap context.

### Tier 3: Clinical Validation

6. **Patient outcome association**: Modules should associate with clinical outcomes (survival, drug response) using TCGA clinical data.
7. **Drug response prediction**: DrugResistanceScorer modules should predict poor response to specific drugs (cross-reference with GDSC/TCGA drug response).

### Tier 4: Methodological Validation

8. **Cross-phenotype specificity**: Modules should have phenotype-specific cores (Jaccard < 0.3) with expected overlap at pleiotropic hubs.
9. **Expression validation**: Predicted invasion module genes should correlate with EMT signatures independently.
10. **Narrative coherence**: Auto-generated narratives should be consistent with published literature for known modules.
11. **Effector set sensitivity**: Results should be robust to ±20% variation in effector gene sets.
12. **Uncertainty calibration**: DS-fused confidence intervals should be well-calibrated (coverage probability matches nominal level).

### Tier 5: Pan-Cancer Validation (NEW)

13. **Conserved hallmark core recovery**: Genes in hallmark H modules across ≥7/10 cancer types should be enriched in known hallmark genes (MSigDB).
14. **Cancer-specific driver validation**: Tissue-specific module genes should match known tissue-specific drivers (e.g., VHL in KIRC, CDH1 in BRCA/STAD).
15. **Cross-cancer actionability**: Pan-cancer drug target matrix should recover known pan-cancer targets (e.g., PIK3CA, KRAS) and identify cancer-specific targets with DepMap support.
16. **Subtype concordance**: Modules should differ between molecular subtypes within the same cancer type (e.g., MSI-high vs MSS in COAD, HPV+ vs HPV- in HNSC).
17. **Pan-cancer pleiotropic hubs**: Genes in ≥3 hallmarks across ≥5 types should be enriched in known master regulators (e.g., TP53, MYC, PIK3CA).

### Tier 5b: Low-Mutation Cancer Ablation (motivation from Paper 1)

**Key insight from Paper 1 pan-cancer results**: Low-mutation cancers (THCA: 1,411 genes, PRAD: 4,772 genes, LGG: 5,944 genes) show lower P@50 (36%, 54%, 68% respectively) compared to high-mutation types (UCEC: 98%, LUSC: 98%). Mutation frequency alone is insufficient for these cancers.

Per-type feature ablation should demonstrate that mechanistic features (gnomAD constraint, DepMap essentiality, and especially the multi-layer molecular features from Paper 2) provide the largest relative AUPRC lift in precisely these low-mutation cancers:

- **Hypothesis**: lift_full_vs_mut_only ∝ 1/mutation_load (inverse relationship)
- **Script**: `scripts/run_per_type_ablation.py` (already implemented)
- **Expected figure**: scatter plot of mutation load (x) vs AUPRC lift (y), showing mechanistic features matter most where mutations are few
- **Narrative**: "In mutationally quiet cancers, Paper 2's multi-layer molecular features (variant-level ΔΔG, splice disruption, methylation silencing) fill the information gap that mutation frequency cannot, enabling driver discovery even in THCA/PRAD"

This directly motivates Paper 2's raison d'être: when mutations alone don't discriminate, you need mechanistic explanation.

### Tier 6: Comparative Benchmarks

18. **Benchmark vs. PhenoDriver**: Direct comparison on BRCA dataset using PhenoDriver's evaluation framework (988 patients).
19. **Benchmark vs. CGMega**: Compare module discovery capability.
20. **Benchmark vs. HallmarkGraph**: Compare hallmark-informed classification vs. our hallmark-informed module discovery.
21. **Benchmark vs. generic driver methods**: Compare against non-phenotype-specific methods (MutSigCV, OncodriveFML, etc.) on phenotype recovery across all 10 cancer types.

---

## Target Journals (revised assessment, post pan-cancer decision)

| Journal | Fit | Requirements | Realistic? |
|---------|-----|-------------|:-:|
| Nature Methods | Methodological: narrative-as-score + DS fusion + CARNIVAL integration | Needs experimental validation or very strong DepMap evidence | Stretch |
| **Genome Biology** | Computational biology audience, strong methods focus | Strong computational validation sufficient | **Best fit** |
| **Nature Communications** | Broader audience, biological discoveries | Pan-cancer scope + novel gene discoveries + clinical validation | **Now realistic** |
| Cell Systems | Systems biology focus | Needs strong systems-level biological insight | Good fit |

**Recommendation**: With 10-cancer pan-cancer scope, **Nature Communications** becomes realistic as primary target if biological discoveries are strong (conserved vs. cancer-specific hallmark drivers, novel drug targets). **Genome Biology** as solid backup — guaranteed fit for methods audience. The pan-cancer analysis elevates the biological contribution significantly: no existing method has done phenotype decomposition across 10 cancer types with patient-level mechanistic narratives.

---

## Implementation Order (revised)

The following is the internal implementation order. The paper presents the full system as a single contribution.

### Phase 0: Pre-computation (can start immediately, parallelize)
1. **TCGA multi-omics download** — Mutations, RNA-seq, miRNA-seq, methylation, CNA, clinical for all 10 cancer types (~5,520 patients) from GDC
2. **Variant impact pre-computation** — Download AlphaMissense, run ESM2 fine-tuned + Rhapsody-2 on TCGA observed variants across all 10 types, cache all SpliceAI + GTEx eQTL scores
3. **CARNIVAL cache** — Pre-compute CARNIVAL solutions for all (gene, phenotype-effector) pairs on OmniPath

### Phase 1: Core System
4. **Heterogeneous network assembly** — Load and integrate OmniPath + Reactome + miRTarBase + DoRothEA into PyG HeteroData (3 node types, 4 edge types)
5. **Phenotype effector curation** — Define effector node sets for 5 hallmarks from MSigDB + GO (apoptosis, invasion, proliferation, immune evasion, drug resistance)
6. **Dempster-Shafer fusion module** — Implement evidence fusion with (belief, plausibility, uncertainty) triplets
7. **PhenotypeScorer** — Three-component scoring with DS fusion (CARNIVAL connectivity × foundation model disruption × phenotype association)
8. **Heterogeneous ACO** — Extend ACO to walk heterogeneous graph (multi-type node selection)
9. **Heterogeneous encoder** — Graph Transformer with edge-attention bias (GRAFT-informed)

### Phase 2: Single-cancer validation (BRCA pilot)
10. **BRCA pilot run** — End-to-end pipeline on BRCA (largest cohort) with all 5 hallmarks to validate system before scaling
11. **Narrative generator** — Template-based explain() from CARNIVAL paths + ESM2 features + DS confidence
12. **DepMap validation suite** — Phenotype-specific essentiality analysis on BRCA modules
13. **Benchmark vs PhenoDriver** — Direct comparison on BRCA dataset (their evaluation framework, 988 patients)

### Phase 3: Pan-cancer execution
14. **Run pipeline on 9 remaining cancer types** — Parallelizable across cancer types (each type is independent)
15. **Cross-cancer aggregation** — Conserved hallmark cores, cancer-specific drivers, pleiotropic hubs, pan-cancer actionability matrix
16. **Drug response integration** — DrugBank/ChEMBL mapping + GDSC validation across all 10 types
17. **Clinical validation** — Survival analysis, drug response prediction per cancer type

### Phase 4: Paper
18. **Pan-cancer figures** — Hallmark conservation heatmaps, cancer-specific driver maps, actionability matrix
19. **Case studies** — Deep narratives for 3-4 showcase cancers (BRCA, KIRC, SKCM, COAD)
20. **Sensitivity analyses** — Effector set variation, uncertainty calibration, cross-phenotype specificity
21. **Benchmark suite** — vs. PhenoDriver, CGMega, HallmarkGraph, generic methods (MutSigCV, OncodriveFML)
22. **Paper writing** — Methods, results, figures, supplementary (per-cancer-type results)

---

## Future Work (explicitly deferred)

These items are valuable but out of scope for this paper. Mention in Discussion:

1. **Full 6-layer heterogeneous network** — Add metabolites, lncRNAs, chromatin regions (Phase 2 of network)
2. **Perturb-seq edge refinement** — Use interventional data to validate/refine network edges at scale
3. **Driver2Comm microenvironment integration** — Connect intracellular driver modules to cell-cell communication (immune evasion, angiogenesis signaling to microenvironment)
4. **Conformal prediction sets** — Formal statistical guarantees on module membership
5. **Prospective clinical validation** — Partner with clinical trial for drug response prediction validation
6. **CGC/OncoKB integration** — Additional ground truth when academic access is available

---

## Timeline

Start after Paper 1 is submitted. Estimated **10-12 months** (revised from 8-10 due to pan-cancer scope) for the complete system and paper.

| Month | Phase | Key deliverables |
|:-----:|-------|------------------|
| 1-2 | Phase 0 | TCGA data for 10 cancer types downloaded; variant impacts pre-computed; CARNIVAL cache built |
| 2-4 | Phase 1 | Heterogeneous network, DS fusion, PhenotypeScorer, heterogeneous ACO |
| 4-5 | Phase 1 | Heterogeneous encoder, end-to-end pipeline working |
| 5-6 | Phase 2 | BRCA pilot: all 5 hallmarks, narratives, DepMap validation, PhenoDriver benchmark |
| 6-8 | Phase 3 | Pan-cancer execution (9 remaining types), cross-cancer aggregation, drug response |
| 8-9 | Phase 3 | Clinical validation, full benchmark suite |
| 9-10 | Phase 4 | Pan-cancer figures, case studies, sensitivity analyses |
| 10-12 | Phase 4 | Paper writing, revision, submission |

**Parallelization opportunities**:
- Phase 0 can overlap with Paper 1 finalization
- Phase 3 cancer types run in parallel (10 independent jobs)
- DepMap validation runs in parallel with pan-cancer execution

---

## Key References

- **PhenoDriver**: Zhang et al. Briefings in Bioinformatics 2023. bbad291. (Primary competitor)
- **CGMega**: Nat Comm 2024. s41467-024-50426-6. (Module dissection with attention)
- **GRAFT**: Briefings Bioinf 2025. bbaf706. (Graph-aware fusion transformer)
- **HGTDG**: Methods 2024. (Heterogeneous graph transformer for driver genes)
- **HallmarkGraph**: Bioinformatics 2025. btaf444. (Cancer hallmark GNN)
- **SEFGNN**: BIBM 2025. arxiv 2510.06290. (Dempster-Shafer evidence fusion)
- **Driver2Comm**: PLoS Comp Bio 2026. pcbi.1013973. (Driver-communication bridge)
- **CARNIVAL**: npj Sys Bio App 2019. s41540-019-0118-z. (Causal reasoning with ILP)
- **ESM2 fine-tuned**: arXiv 2410.10919. (Protein feature classification)
- **Rhapsody-2**: PNAS 2025. pnas.2418100122. (Mechanistic variant evaluation)
- **AlphaMissense**: Science 2023. adj8672. (71M missense variant predictions)
- **Conformal prediction genomics**: Frontiers Bioinf 2025. fbinf.2025.1507448.
