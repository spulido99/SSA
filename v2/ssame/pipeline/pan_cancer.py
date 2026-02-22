"""Pan-Cancer framework (WP10).

Runs the SSA-ME pipeline independently per cancer type, then
aggregates results using Borda rank aggregation to identify
recurrent driver genes across cancers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ssame.network.multi_view_graph import MultiViewGraph
from ssame.pipeline.pipeline import Pipeline, PipelineResult
from ssame.types import AlterationMatrix


@dataclass
class PanCancerResult:
    """Output of a pan-cancer analysis."""

    per_cancer: dict[str, PipelineResult]
    aggregated_ranking: list[tuple[str, float]]
    recurrence_count: dict[str, int]


class PanCancerRunner:
    """Run SSA-ME independently per cancer type and aggregate.

    Args:
        alteration_matrices: cancer_type → AlterationMatrix.
        labels: Shared known driver labels (gene → 1/0).
        encoder_level: 0 (XGBoost) or 1 (Graph Transformer).
        encoder_n_estimators: Number of XGBoost estimators (L0).
        search_n_iterations: ACO iterations.
        seed: Random seed.
        top_k_fraction: Fraction of genes flagged as "hit" per cancer for recurrence.
    """

    def __init__(
        self,
        alteration_matrices: dict[str, AlterationMatrix],
        labels: dict[str, int],
        *,
        encoder_level: int = 0,
        encoder_n_estimators: int = 100,
        search_n_iterations: int = 1000,
        seed: int = 42,
        top_k_fraction: float = 0.1,
    ) -> None:
        self._matrices = alteration_matrices
        self._labels = labels
        self._encoder_level = encoder_level
        self._encoder_n_estimators = encoder_n_estimators
        self._search_n_iterations = search_n_iterations
        self._seed = seed
        self._top_k_fraction = top_k_fraction

    def run(self) -> PanCancerResult:
        """Run pipeline for each cancer type and aggregate."""
        per_cancer: dict[str, PipelineResult] = {}

        for cancer_type, alt_matrix in self._matrices.items():
            # Build a minimal single-view graph from the alteration matrix genes
            multi_view = self._build_coalteration_graph(alt_matrix)
            pipeline = Pipeline(
                alteration_matrix=alt_matrix,
                multi_view=multi_view,
                labels=self._labels,
                encoder_level=self._encoder_level,
                encoder_n_estimators=self._encoder_n_estimators,
                search_n_iterations=self._search_n_iterations,
                seed=self._seed,
            )
            per_cancer[cancer_type] = pipeline.run()

        aggregated = self._borda_aggregate(per_cancer)
        recurrence = self._compute_recurrence(per_cancer)

        return PanCancerResult(
            per_cancer=per_cancer,
            aggregated_ranking=aggregated,
            recurrence_count=recurrence,
        )

    def _build_coalteration_graph(self, alt_matrix: AlterationMatrix) -> MultiViewGraph:
        """Build a co-alteration network from the alteration matrix.

        Two genes are connected if they are co-altered in >= 2 samples.
        """
        data = alt_matrix.data.astype(np.float32)
        co_occur = data @ data.T
        np.fill_diagonal(co_occur, 0)
        adj = (co_occur >= 2).astype(np.float32)

        node_to_index = {g: i for i, g in enumerate(alt_matrix.genes)}
        return MultiViewGraph(
            node_list=alt_matrix.genes,
            node_to_index=node_to_index,
            view_names=["coalteration"],
            _adjacencies={"coalteration": adj},
        )

    def _borda_aggregate(
        self, per_cancer: dict[str, PipelineResult]
    ) -> list[tuple[str, float]]:
        """Borda count rank aggregation across cancer types."""
        # Collect all genes
        all_genes: set[str] = set()
        for result in per_cancer.values():
            all_genes.update(result.gene_scores.keys())

        n_cancers = len(per_cancer)
        n_genes = len(all_genes)

        # For each cancer, compute rank (1 = best)
        gene_borda: dict[str, float] = {g: 0.0 for g in all_genes}

        for result in per_cancer.values():
            # Sort by score descending → rank 1 = highest
            ranked = sorted(result.gene_scores.items(), key=lambda x: -x[1])
            gene_set = {g for g, _ in ranked}
            for rank, (gene, _score) in enumerate(ranked, 1):
                # Borda score: (n_genes - rank) / n_genes
                gene_borda[gene] += (n_genes - rank) / n_genes

            # Genes not present in this cancer get rank = n_genes (worst)
            for gene in all_genes - gene_set:
                gene_borda[gene] += 0.0

        # Normalize by number of cancers
        for gene in gene_borda:
            gene_borda[gene] /= n_cancers

        return sorted(gene_borda.items(), key=lambda x: -x[1])

    def _compute_recurrence(
        self, per_cancer: dict[str, PipelineResult]
    ) -> dict[str, int]:
        """Count in how many cancer types each gene is in the top K."""
        recurrence: dict[str, int] = {}

        all_genes: set[str] = set()
        for result in per_cancer.values():
            all_genes.update(result.gene_scores.keys())

        for gene in all_genes:
            recurrence[gene] = 0

        for result in per_cancer.values():
            ranked = sorted(result.gene_scores.items(), key=lambda x: -x[1])
            top_k = max(1, int(len(ranked) * self._top_k_fraction))
            top_genes = {g for g, _ in ranked[:top_k]}
            for gene in top_genes:
                recurrence[gene] = recurrence.get(gene, 0) + 1

        return recurrence
