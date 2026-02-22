# Paper 2: Critical Review & Competitive Analysis

**Date**: 2026-02-21
**Purpose**: Reference document for implementation decisions. Summarizes critical review, competitive landscape, and rationale for improvements.

---

## Competitive Landscape Detail

### PhenoDriver (Briefings in Bioinformatics, Sep 2023)
- **What it does**: Personalized phenotype-associated driver gene identification + subnetwork construction
- **Dataset**: 988 breast cancer patients
- **Strengths**: Identifies recurrent AND rare driver mutations; discovers cancer subtypes from personalized driver profiles
- **Weaknesses**: Homogeneous PPI graph only; gene-level (not variant-level); no search optimization; no uncertainty quantification
- **Our differentiation**: Heterogeneous multi-layer graph; variant-level mechanistic quantification (ESM2, Rhapsody-2); narrative-as-score paradigm; ACO search; Dempster-Shafer uncertainty
- **Action**: Explicit comparison table in Introduction; benchmark on BRCA dataset using their evaluation framework

### CGMega (Nature Communications, Jul 2024)
- **What it does**: Explainable graph attention framework for cancer gene module dissection
- **Dataset**: Breast cancer cell lines, AML patients
- **Strengths**: Multi-omics (3D genome, epigenome, PPI); GNNExplainer for module interpretation; 396 candidate AML genes; AUPRC=0.914
- **Weaknesses**: Not phenotype-specific; not patient-specific; no causal reasoning; no variant-level analysis
- **Our differentiation**: Phenotype-directed (not generic); patient-specific variant narratives; CARNIVAL causal paths
- **Action**: Benchmark module discovery on AML dataset if feasible

### GRAFT (Briefings in Bioinformatics, Jan 2025)
- **What it does**: Graph-Aware Fusion Transformer for cancer driver gene prediction
- **Dataset**: Pan-cancer + individual cancer types
- **Strengths**: Multi-view graph encoder (PPI, pathway co-occurrence, gene semantic similarity); edge-attention bias; structural + functional embeddings
- **Weaknesses**: Not phenotype-specific; no patient-level analysis; classification only (not module discovery)
- **Our differentiation**: Phenotype decomposition; module discovery with narratives; ACO search
- **Action**: Adopt edge-attention bias architecture for our heterogeneous encoder

### HGTDG (Methods, Oct 2024)
- **What it does**: Heterogeneous Graph Transformer for cancer driver gene prediction
- **Dataset**: TCGA 16 cancer types with multi-omics
- **Strengths**: True heterogeneous graph; multi-omics integration; transformer architecture
- **Weaknesses**: Not phenotype-specific; no variant-level mechanistic analysis; no module discovery
- **Our differentiation**: Phenotype-directed scoring; variant-level narratives; ACO module discovery
- **Action**: Use as reference for heterogeneous graph construction patterns

### HallmarkGraph (Bioinformatics, Sep 2025)
- **What it does**: Cancer hallmark-informed GNN for hierarchical tumor subtype classification
- **Dataset**: 10,000+ samples, 26 cancer types, 405 subtypes
- **Strengths**: First hallmark-informed GNN; hierarchical classification (up to 8 levels); 85-99% accuracy
- **Weaknesses**: Classification only (not discovery); expression data only; no variant analysis; no module output
- **Our differentiation**: Module DISCOVERY with novel gene predictions; multi-omics; patient-specific narratives
- **Action**: Use their hallmark gene sets (6,448 genes from Reactome) as reference for effector curation

### SEFGNN (IEEE BIBM, Oct 2025)
- **What it does**: Dempster-Shafer evidence fusion across multiple PPI networks for driver gene identification
- **Dataset**: 3 cancer types
- **Strengths**: Uncertainty-aware fusion; Soft Evidence Smoothing; principled evidence combination
- **Weaknesses**: PPI networks only; not phenotype-specific; not patient-specific; no mechanistic explanation
- **Our differentiation**: Multi-layer heterogeneous graph; phenotype-specific; variant-level; mechanistic narratives
- **Action**: Adopt DS fusion approach for our composite score; cite as methodological precedent

### Driver2Comm (PLoS Computational Biology, Feb 2026)
- **What it does**: Bridges cancer driver genes with cell-cell communication events
- **Dataset**: Pancreatic, breast, neuroblastoma, oral squamous cell carcinoma (scRNA-seq + spatial)
- **Strengths**: Links intracellular drivers to microenvironment; predicts immunotherapy response and survival
- **Weaknesses**: No heterogeneous graph; no variant-level analysis; no mechanistic scoring
- **Our differentiation**: Mechanistic variant-level narratives; heterogeneous graph; ACO search
- **Action**: Mention as future integration (microenvironment layer); cite in Discussion

---

## Key Technical References for Implementation

### CARNIVAL (Causal Reasoning)
- Paper: Liu et al. npj Systems Biology and Applications 2019
- Code: https://saezlab.github.io/CARNIVAL/
- R package; need Python wrapper or port
- Uses OmniPath as prior knowledge network (9,306 signed directed edges, 3,610 nodes)
- ILP solver: CPLEX or CBC (open source)
- Benchmarked best for pathway recovery with OmniPath (BMC Bioinformatics 2023)
- **Implementation note**: May need to implement core ILP formulation in Python using PuLP or OR-Tools rather than depending on R package

### ESM2 Fine-tuned (Variant Mechanistic Classification)
- Paper: arXiv 2410.10919
- Model: ESM2-650M fine-tuned on 20 protein features
- Classifies at amino acid resolution: DNA_binding, kinase_activity, PTM_site, etc.
- Can run on RTX 3060 (650M params ~2.5GB VRAM in fp16)
- **Implementation note**: Pre-compute for all TCGA observed missense variants; cache feature_loss_vector per variant

### Rhapsody-2 (Structural Variant Evaluation)
- Paper: PNAS 2025, pnas.2418100122
- Combines AlphaFold2 structure + evolutionary conservation + protein dynamics
- Outputs: ΔΔG stability, ΔΔG binding, flexibility change, mechanistic category
- **Implementation note**: May be computationally expensive; prioritize key genes in modules

### AlphaMissense (Pathogenicity Screening)
- Paper: Science 2023
- Pre-computed scores for 71M possible human missense variants
- Download: https://zenodo.org/records/8208688
- Simple lookup table — no computation needed
- **Implementation note**: ~3GB download, index by (chromosome, position, ref, alt)

### Dempster-Shafer Evidence Fusion
- Theory: Well-established in decision theory
- Python libraries: pyds, dst (lightweight)
- SEFGNN implementation: https://github.com/hanzo2020/SEFGNN (reference for Soft Evidence Smoothing)
- **Implementation note**: Need custom implementation for 3-component fusion with heterogeneous evidence types

### Conformal Prediction (Optional)
- Paper: Frontiers Bioinformatics 2025
- Python library: MAPIE, crepes
- Distribution-free, model-agnostic
- **Implementation note**: Apply post-ACO for calibrated prediction sets; low implementation effort

---

## Risk Mitigation Strategies

| Risk | Mitigation | Status |
|------|-----------|--------|
| PhenoDriver overlap | Explicit differentiation in intro; benchmark on BRCA | Planned |
| Computational cost | ALL variant impacts pre-computed; CARNIVAL cached | Design decision |
| Phenotype label circularity | Prefer clinical endpoints; exclude effectors from surrogates | Design decision |
| Effector set subjectivity | Sensitivity analysis ±20% | Planned |
| Cross-phenotype specificity | Revised: expect pleiotropic overlap, measure specificity ratio | Design decision |
| No wet-lab validation | DepMap phenotype-specific essentiality as central validation | Design decision |
| 6GB VRAM limit | ESM2-650M fits in fp16; batch processing for large datasets | Verified |
| CARNIVAL Python port | Implement ILP core with PuLP/OR-Tools; OmniPath has Python API (pypath) | To verify |
| Timeline 10-12 months | Phase 0 parallelizable with Paper 1; pan-cancer types run in parallel | Planned |
| Pan-cancer compute cost | 10 types × 5 hallmarks = 50 ACO runs; each independent, parallelizable | Design decision |
| Data download volume | ~5,520 patients × 6 omics types from GDC; estimate 200-500 GB total | To verify |

---

## Improvement Priority Matrix

| # | Improvement | Scientific impact | Implementation effort | Priority |
|---|-----------|:-:|:-:|:-:|
| 1 | Dempster-Shafer fusion + uncertainty | HIGH | Medium (2-3 weeks) | **P1** |
| 2 | CARNIVAL causal pathway reasoning | HIGH | Medium (3-4 weeks) | **P1** |
| 3 | Foundation models (AlphaMissense/ESM2/Rhapsody-2) | HIGH | Low (pre-computed) | **P1** |
| 4 | DepMap phenotype-specific validation | HIGH | Low (data analysis) | **P1** |
| 5 | Drug response integration | HIGH | Medium (2-3 weeks) | **P1** |
| 6 | GRAFT/HGTDG-informed encoder | Medium | Medium (3-4 weeks) | P2 |
| 7 | Perturb-seq edge refinement | Medium | High (4-6 weeks) | P3 |
| 8 | Driver2Comm microenvironment | Medium | High (6-8 weeks) | P3 (future) |

---

## Pan-Cancer Cohort Selection (decided 2026-02-21)

10 TCGA cancer types selected for maximum data × therapeutic relevance:

| # | Code | Cancer | n | Key biology for phenotype decomposition |
|:-:|:----:|--------|:-:|----------------------------------------|
| 1 | BRCA | Breast | 1,097 | PhenoDriver benchmark; Luminal/Basal subtypes; PIK3CA/TP53/CDH1 |
| 2 | KIRC | Kidney cc | 537 | VHL→HIF1A→VEGFA perfect for multi-layer narrative; chemo-resistant |
| 3 | HNSC | Head/neck | 528 | HPV+/HPV- = two diseases = natural internal control |
| 4 | LUAD | Lung adeno | 522 | EGFR/ALK/KRAS targets; #1 mortality; drug response data rich |
| 5 | LUSC | Lung squamous | 504 | No targeted therapy → need for novel target discovery |
| 6 | SKCM | Melanoma | 470 | Highest TMB; immunotherapy paradigm; PD-L1/CTLA4 |
| 7 | COAD | Colon | 459 | MSI-H vs MSS → immunotherapy response; WNT/RAS pathways |
| 8 | STAD | Stomach | 443 | 4 subtypes (EBV/MSI/GS/CIN); max molecular diversity |
| 9 | BLCA | Bladder | 412 | Checkpoint inhibitor responsive; luminal vs basal |
| 10 | UCEC | Endometrium | 548 | POLE ultramutated/MSI/CN-low/CN-high subtypes |

**Total**: ~5,520 patients. **Excluded**: GBM/LGG (BBB), OV (few novel options), THCA (good prognosis), PRAD (indolent), LIHC (viral confound).

**Expected pan-cancer contributions**:
1. Conserved hallmark cores (genes in hallmark H across ≥7/10 types)
2. Cancer-specific hallmark drivers (unique to 1-2 types)
3. Pleiotropic master regulators (≥3 hallmarks, ≥5 types)
4. Pan-cancer actionability matrix (module genes × drugs × cancer types)
5. Subtype-aware decomposition (MSI vs MSS, HPV+ vs HPV-, etc.)
