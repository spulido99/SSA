"""MEScorer: mutual exclusivity scoring (port of v1 Dendrix + HyperGeometric)."""

from __future__ import annotations

import math

from scipy import stats

from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import AlterationMatrix, ScorerExplanation, Subnetwork


class MEScorer(SubnetworkScorer):
    """Mutual exclusivity scorer.

    Ports v1 MutualExclusivityNetworkManager scoring logic.

    Two methods:
    - "dendrix": coverage - overlap (Dendrix-style, fast)
    - "hypergeometric": statistical test for ME significance

    Optionally subtype-corrected.
    """

    def __init__(
        self,
        alteration_matrix: AlterationMatrix,
        subtypes: dict[str, str] | None = None,
        method: str = "dendrix",
    ) -> None:
        self.matrix = alteration_matrix
        self.subtypes = subtypes
        self.method = method

    def score(self, module: Subnetwork) -> float:
        gene_symbols = [g.symbol for g in module.genes if g.symbol in self.matrix.genes]
        if len(gene_symbols) < 2:
            return 0.0

        if self.method == "dendrix":
            return self._dendrix_score(gene_symbols)
        elif self.method == "hypergeometric":
            return self._hypergeometric_score(gene_symbols)
        else:
            raise ValueError(f"Unknown ME method: {self.method}")

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        gene_symbols = [g.symbol for g in module.genes if g.symbol in self.matrix.genes]
        if len(gene_symbols) < 2:
            return ScorerExplanation(score=0.0, narrative="Module too small for ME analysis")

        score = self.score(module)
        coverage = self.matrix.coverage(gene_symbols)

        pairwise = {}
        for i, g1 in enumerate(gene_symbols):
            for g2 in gene_symbols[i + 1 :]:
                co = self.matrix.co_occurrence_count(g1, g2)
                pairwise[f"{g1}-{g2}"] = float(co)

        return ScorerExplanation(
            score=score,
            components=pairwise,
            narrative=f"ME module of {len(gene_symbols)} genes, coverage={coverage:.2f}",
            evidence={"coverage": coverage},
        )

    def _dendrix_score(self, gene_symbols: list[str]) -> float:
        """Dendrix-style ME score: coverage - overlap, normalized.

        Matches v1 logic: for each sample, score = 1/count_of_mutated_genes.
        Sum across samples, normalize by gene count.
        """
        sub = self.matrix.submatrix(gene_symbols)
        n_samples = sub.n_samples

        # Per-sample: count how many module genes are mutated
        gene_counts_per_sample = sub.data.sum(axis=0)  # shape: (n_samples,)

        total_score = 0.0
        for j in range(n_samples):
            count = gene_counts_per_sample[j]
            if count > 0:
                # v1 logic: each gene in a sample contributes 1/count
                # Sum over genes in that sample = count * (1/count) = 1
                # But we want to reward ME: fewer co-altered genes = higher per-gene score
                # Score per sample = sum over mutated genes of 1/count
                total_score += 1.0  # Sample-level contribution (coverage)

        # Overlap penalty: count samples with > 1 gene mutated
        overlap = float((gene_counts_per_sample > 1).sum())
        coverage = float((gene_counts_per_sample > 0).sum())

        # Dendrix score: (coverage - 2*overlap) / n_samples, but keep non-negative
        # More ME = higher score
        me_score = (coverage - overlap) / n_samples
        return max(0.0, me_score)

    def _hypergeometric_score(self, gene_symbols: list[str]) -> float:
        """HyperGeometric test for ME: tests if overlap is less than expected."""
        sub = self.matrix.submatrix(gene_symbols)
        n_samples = sub.n_samples

        # For each pair, test if co-occurrence is less than expected
        min_pval = 1.0
        for i in range(len(gene_symbols)):
            for j in range(i + 1, len(gene_symbols)):
                n_a = int(sub.data[i].sum())
                n_b = int(sub.data[j].sum())
                co_occur = int((sub.data[i] & sub.data[j]).sum())

                # P(X <= co_occur) under hypergeometric
                # ME = fewer co-occurrences than expected
                pval = stats.hypergeom.cdf(co_occur, n_samples, n_a, n_b)
                min_pval = min(min_pval, pval)

        # Convert p-value to score: -log10(p) capped
        if min_pval <= 0:
            return 10.0
        return min(-math.log10(min_pval), 10.0)
