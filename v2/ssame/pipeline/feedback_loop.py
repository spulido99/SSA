"""Score-as-label feedback loop: module scores feed back as soft labels."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ssame.network.multi_view_graph import MultiViewGraph
from ssame.pipeline.pipeline import Pipeline, PipelineResult
from ssame.types import AlterationMatrix


@dataclass
class FeedbackResult:
    """Results from the full feedback pipeline."""

    round_results: list[PipelineResult]
    final_ranked_genes: list[tuple[str, float]]


class FeedbackPipeline:
    """Iterative pipeline with score-as-label feedback.

    After each round:
    1. Run full pipeline (encode -> search -> score)
    2. Genes in high-scoring modules get soft label boost
    3. Updated soft labels feed the next round's encoder training

    This implements the feedback loop from v1's Fungus where
    pheromone levels (posterior probabilities) are reinforced by
    subnetwork quality scores.
    """

    def __init__(
        self,
        alteration_matrix: AlterationMatrix,
        multi_view: MultiViewGraph,
        labels: dict[str, int],
        *,
        n_rounds: int = 3,
        feedback_weight: float = 0.3,
        encoder_n_estimators: int = 100,
        search_n_ants: int = 50,
        search_n_iterations: int = 1000,
        search_subnetwork_size: int = 3,
        seed: int = 42,
    ) -> None:
        self.alteration_matrix = alteration_matrix
        self.multi_view = multi_view
        self.original_labels = labels
        self.n_rounds = n_rounds
        self.feedback_weight = feedback_weight
        self.encoder_n_estimators = encoder_n_estimators
        self.search_n_ants = search_n_ants
        self.search_n_iterations = search_n_iterations
        self.search_subnetwork_size = search_subnetwork_size
        self.seed = seed

    def run(self) -> FeedbackResult:
        """Run iterative feedback pipeline."""
        round_results: list[PipelineResult] = []
        current_labels = dict(self.original_labels)

        for round_idx in range(self.n_rounds):
            pipeline = Pipeline(
                alteration_matrix=self.alteration_matrix,
                multi_view=self.multi_view,
                labels=current_labels,
                encoder_n_estimators=self.encoder_n_estimators,
                search_n_ants=self.search_n_ants,
                search_n_iterations=self.search_n_iterations,
                search_subnetwork_size=self.search_subnetwork_size,
                seed=self.seed + round_idx,
            )
            result = pipeline.run()
            round_results.append(result)

            # Update labels for next round (except last)
            if round_idx < self.n_rounds - 1:
                current_labels = self._update_labels(
                    self.original_labels, result
                )

        # Final ranking from last round
        final_ranked = round_results[-1].ranked_genes

        return FeedbackResult(
            round_results=round_results,
            final_ranked_genes=final_ranked,
        )

    def _update_labels(
        self,
        original_labels: dict[str, int],
        result: PipelineResult,
    ) -> dict[str, int]:
        """Create updated soft labels by boosting genes in high-scoring modules.

        Known positives (label=1) remain 1. Unlabeled genes that appear
        in high-scoring modules get promoted to positive for the next round.
        """
        updated = dict(original_labels)

        if not result.module_scores:
            return updated

        # Find score threshold (top 25% of module scores)
        scores_arr = np.array(result.module_scores)
        if len(scores_arr) == 0:
            return updated
        threshold = np.percentile(scores_arr, 75)

        # Genes in top-scoring modules get label boost
        for module, score in zip(result.modules, result.module_scores, strict=True):
            if score >= threshold:
                for gene in module.genes:
                    if gene.symbol not in original_labels or original_labels[gene.symbol] == 0:
                        # Soft promotion: only if feedback_weight roll passes
                        # Deterministic: use gene score as proxy
                        gene_score = result.gene_scores.get(gene.symbol, 0.0)
                        if gene_score > self.feedback_weight:
                            updated[gene.symbol] = 1

        return updated
