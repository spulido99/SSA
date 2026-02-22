"""Functional genomics scorer: scores modules by CRISPR essentiality enrichment."""

from __future__ import annotations

from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import ScorerExplanation, Subnetwork


class FunctionalGenomicsScorer(SubnetworkScorer):
    """Score subnetworks by enrichment for functionally essential genes.

    Uses DepMap CRISPR essentiality scores (more negative = more essential).
    A module enriched for essential genes is more likely to be a real driver module.
    """

    def __init__(
        self,
        essentiality_scores: dict[str, float],
        threshold: float = -0.5,
    ) -> None:
        self.essentiality_scores = essentiality_scores
        self.threshold = threshold

    def score(self, module: Subnetwork) -> float:
        """Score module by mean essentiality of its genes.

        Transforms negative essentiality scores to positive module scores:
        more essential genes (more negative) produce higher module scores.
        """
        gene_scores = []
        for gene in module.genes:
            if gene.symbol in self.essentiality_scores:
                gene_scores.append(self.essentiality_scores[gene.symbol])

        if not gene_scores:
            return 0.0

        # Transform: more negative essentiality -> higher score
        mean_essentiality = sum(gene_scores) / len(gene_scores)
        return max(0.0, -mean_essentiality)

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        """Explain the functional genomics score."""
        components: dict[str, float] = {}
        for gene in module.genes:
            if gene.symbol in self.essentiality_scores:
                components[gene.symbol] = self.essentiality_scores[gene.symbol]

        n_essential = sum(
            1 for v in components.values() if v < self.threshold
        )

        return ScorerExplanation(
            score=self.score(module),
            components=components,
            narrative=(
                f"Module contains {n_essential}/{len(module.genes)} genes with "
                f"CRISPR essentiality below {self.threshold} threshold."
            ),
        )
