"""gnomAD constraint score loader (pLI, LOEUF)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class GnomADResult:
    """Parsed gnomAD gene constraint data."""

    pli: dict[str, float]
    loeuf: dict[str, float]

    def get_constrained_genes(
        self,
        pli_threshold: float = 0.9,
        loeuf_threshold: float | None = None,
    ) -> set[str]:
        """Return constrained genes by pLI and/or LOEUF threshold.

        Args:
            pli_threshold: Genes with pLI >= this are constrained.
            loeuf_threshold: If set, genes with LOEUF <= this are also included.

        Returns:
            Set of constrained gene symbols.
        """
        constrained = {g for g, v in self.pli.items() if v >= pli_threshold}
        if loeuf_threshold is not None:
            constrained |= {g for g, v in self.loeuf.items() if v <= loeuf_threshold}
        return constrained

    def as_features(self, gene_list: list[str]) -> np.ndarray:
        """Return (n_genes, 2) feature matrix [pLI, LOEUF] for given genes.

        Missing genes get (0.0, 1.0) defaults (no constraint assumed).
        """
        features = np.zeros((len(gene_list), 2))
        for i, gene in enumerate(gene_list):
            features[i, 0] = self.pli.get(gene, 0.0)
            features[i, 1] = self.loeuf.get(gene, 1.0)
        return features


def load_gnomad_constraints(path: Path) -> GnomADResult:
    """Load gnomAD gene constraint metrics.

    Supports both v4.1 (constraint_metrics.tsv) and v2.1.1 (lof_metrics.by_gene.txt).

    Args:
        path: Path to gnomAD constraint TSV.

    Returns:
        GnomADResult with per-gene pLI and LOEUF scores.
    """
    df = pd.read_csv(path, sep="\t")

    # Normalize column names for different gnomAD versions
    col_map = {c.lower(): c for c in df.columns}

    gene_col = col_map.get("gene", col_map.get("gene_symbol", "gene"))
    pli_col = col_map.get("pli", "pLI")
    loeuf_col = col_map.get("loeuf", col_map.get("oe_lof_upper", "loeuf"))

    pli: dict[str, float] = {}
    loeuf: dict[str, float] = {}

    for _, row in df.iterrows():
        gene = str(row[gene_col])
        try:
            pli_val = float(row[pli_col])
        except (ValueError, TypeError):
            pli_val = 0.0
        try:
            loeuf_val = float(row[loeuf_col])
        except (ValueError, TypeError):
            loeuf_val = 1.0

        # Keep first transcript per gene (usually canonical)
        if gene not in pli:
            pli[gene] = pli_val
            loeuf[gene] = loeuf_val

    return GnomADResult(pli=pli, loeuf=loeuf)
