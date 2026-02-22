"""DepMap CRISPR essentiality data loader."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_GENE_ENTREZ_RE = re.compile(r"^(.+)\s+\((\d+)\)$")


@dataclass
class DepMapResult:
    """Parsed DepMap CRISPR essentiality data."""

    mean_essentiality: dict[str, float]
    selective_dependency: dict[str, float]
    entrez_ids: dict[str, int]
    matrix: np.ndarray
    gene_names: list[str]
    cell_line_ids: list[str]

    def get_essential_genes(self, threshold: float = -0.5) -> set[str]:
        """Return genes with mean essentiality below threshold."""
        return {
            gene
            for gene, score in self.mean_essentiality.items()
            if score < threshold
        }


def load_depmap_essentiality(
    path: Path,
    threshold: float | None = None,
) -> DepMapResult:
    """Load DepMap CRISPRGeneEffect CSV.

    Column format: 'GENE (ENTREZ_ID)' or just 'GENE'.
    Rows: cell lines (ACH-XXXXXX).

    Args:
        path: Path to CRISPRGeneEffect.csv.
        threshold: Optional essentiality threshold for filtering.

    Returns:
        DepMapResult with per-gene essentiality scores.
    """
    df = pd.read_csv(path, index_col=0)

    # Parse column names: "TP53 (7157)" -> ("TP53", 7157)
    gene_names: list[str] = []
    entrez_ids: dict[str, int] = {}

    for col in df.columns:
        match = _GENE_ENTREZ_RE.match(col)
        if match:
            gene = match.group(1)
            entrez = int(match.group(2))
            gene_names.append(gene)
            entrez_ids[gene] = entrez
        else:
            gene_names.append(col)

    df.columns = gene_names

    matrix = df.values.astype(np.float64)
    cell_line_ids = df.index.tolist()

    # Per-gene statistics
    mean_essentiality = {
        gene: float(np.nanmean(matrix[:, i])) for i, gene in enumerate(gene_names)
    }

    # Selective dependency: variance across cell lines (high variance = selective)
    selective_dependency = {
        gene: float(np.nanvar(matrix[:, i])) for i, gene in enumerate(gene_names)
    }

    return DepMapResult(
        mean_essentiality=mean_essentiality,
        selective_dependency=selective_dependency,
        entrez_ids=entrez_ids,
        matrix=matrix,
        gene_names=gene_names,
        cell_line_ids=cell_line_ids,
    )
