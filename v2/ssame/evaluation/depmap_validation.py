"""DepMap functional validation: correlate predictions with CRISPR essentiality."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class DepMapValidationResult:
    """Results from DepMap functional validation."""

    correlation: float
    p_value: float
    enrichment_top_k: float
    mean_essentiality_top: float
    mean_essentiality_bottom: float
    n_genes_evaluated: int


class DepMapValidator:
    """Validate predicted driverness scores against DepMap CRISPR essentiality.

    Key insight: if predictions are biologically meaningful, top-ranked
    genes should be enriched for CRISPR-essential genes (more negative
    essentiality scores in DepMap).
    """

    def __init__(
        self,
        gene_scores: dict[str, float],
        essentiality_scores: dict[str, float],
        top_k_fraction: float = 0.25,
    ) -> None:
        self.gene_scores = gene_scores
        self.essentiality_scores = essentiality_scores
        self.top_k_fraction = top_k_fraction

    def validate(self) -> DepMapValidationResult:
        """Run DepMap validation.

        Returns:
            DepMapValidationResult with correlation and enrichment metrics.
        """
        # Common genes
        common = sorted(
            set(self.gene_scores) & set(self.essentiality_scores)
        )

        if len(common) < 3:
            return DepMapValidationResult(
                correlation=0.0,
                p_value=1.0,
                enrichment_top_k=0.0,
                mean_essentiality_top=0.0,
                mean_essentiality_bottom=0.0,
                n_genes_evaluated=len(common),
            )

        scores = np.array([self.gene_scores[g] for g in common])
        essentiality = np.array([self.essentiality_scores[g] for g in common])

        # Spearman correlation (predicted score vs essentiality)
        corr, p_val = stats.spearmanr(scores, essentiality)

        # Top-k enrichment
        k = max(1, int(len(common) * self.top_k_fraction))
        ranking = np.argsort(-scores)

        top_k_idx = ranking[:k]
        bottom_k_idx = ranking[-k:]

        mean_top = float(np.mean(essentiality[top_k_idx]))
        mean_bottom = float(np.mean(essentiality[bottom_k_idx]))

        # Enrichment: how much more essential are top predictions vs bottom
        enrichment = mean_bottom - mean_top  # Positive if top is more essential

        return DepMapValidationResult(
            correlation=float(corr),
            p_value=float(p_val),
            enrichment_top_k=float(enrichment),
            mean_essentiality_top=mean_top,
            mean_essentiality_bottom=mean_bottom,
            n_genes_evaluated=len(common),
        )
