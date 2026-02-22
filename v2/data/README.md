# Data Directory

All data is excluded from git. Download using scripts or manually as described below.

## Directory Structure

```
data/
├── raw/
│   ├── tcga/                  # TCGA somatic mutations, clinical, CNA
│   │   ├── mc3.maf.gz        # MC3 unified MAF (GDC API, public)
│   │   └── clinical/         # Clinical + subtype data
│   ├── networks/
│   │   ├── reactome/          # Reactome pathways (public)
│   │   ├── string/            # STRING PPI (public, biased baseline)
│   │   ├── hi_union/          # HI-II-14 proteome screen (public)
│   │   ├── alphafold_baker/   # Baker 2025 predicted PPI (manual)
│   │   └── pioneer/           # PIONEER oncoPPI (manual)
│   ├── ground_truth/
│   │   ├── cancer_gene_census.csv  # COSMIC CGC (manual, requires license)
│   │   ├── oncokb_cancer_genes.tsv # OncoKB (manual)
│   │   └── genie/             # AACR GENIE (manual, requires Synapse)
│   ├── depmap/                # DepMap CRISPR essentiality (public)
│   ├── embeddings/
│   │   └── genept/            # GenePT gene embeddings (manual)
│   └── annotations/
│       └── gnomad_constraint.tsv  # gnomAD pLI/LOEUF (public)
├── processed/                 # Processed/transformed data
└── cache/                     # Intermediate computation cache
```

## Data Sources

| Source | Access | Script |
|--------|--------|--------|
| TCGA MC3 MAF | Public (GDC API) | `scripts/download_tcga.py` |
| DepMap | Public (FTP) | `scripts/download_depmap.py` |
| Reactome | Public | `scripts/download_reactome.py` |
| STRING | Public | `scripts/download_string.py` |
| HI-II-14 | Public | `scripts/download_hi_union.py` |
| gnomAD | Public | `scripts/download_gnomad.py` |
| COSMIC CGC | Academic license required | Manual download |
| OncoKB | Free registration | Manual download |
| Baker AlphaFold PPI | Paper supplementary | Manual download |
| PIONEER oncoPPI | Paper supplementary | Manual download |
| AACR GENIE | Synapse DUA required | Manual download |
| GenePT | GitHub repo | Manual download |
