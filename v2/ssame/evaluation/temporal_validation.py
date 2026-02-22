"""Temporal validation: train on pre-2020 labels, evaluate on 2020-2025 holdout."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ssame.encoder.tabular_baseline import TabularBaseline


@dataclass
class TemporalResult:
    """Results from temporal validation."""

    metrics: dict[str, float]
    n_holdout_in_top_k: int
    ranked_scores: np.ndarray
    holdout_ranks: list[int] = field(default_factory=list)


class TemporalValidator:
    """Temporal validation: train on historical labels, evaluate on new discoveries.

    Simulates the real prediction scenario: train on CGC genes known before
    a cutoff year, then check how many post-cutoff CGC additions appear in
    the top predictions.
    """

    def __init__(
        self,
        n_estimators: int = 100,
        seed: int = 42,
        top_k: int = 50,
    ) -> None:
        self.n_estimators = n_estimators
        self.seed = seed
        self.top_k = top_k

    def run(
        self,
        X: np.ndarray,
        train_labels: np.ndarray,
        holdout_labels: np.ndarray,
    ) -> TemporalResult:
        """Run temporal validation.

        Args:
            X: Feature matrix (n_genes, n_features).
            train_labels: PU labels for training (1=pre-cutoff known positive).
            holdout_labels: Labels for holdout (1=post-cutoff new addition).

        Returns:
            TemporalResult with recall metrics on holdout set.
        """
        # Train on historical labels
        model = TabularBaseline(
            n_estimators=self.n_estimators,
            seed=self.seed,
        )
        model.fit(X, train_labels)
        scores = model.predict(X)

        # Rank genes by predicted score
        ranking = np.argsort(-scores)

        # Find holdout genes in top-k predictions
        holdout_idx = np.where(holdout_labels == 1)[0]
        n_holdout = len(holdout_idx)

        top_k_50 = set(ranking[: min(50, len(ranking))])
        top_k_100 = set(ranking[: min(100, len(ranking))])

        holdout_in_top_50 = sum(1 for idx in holdout_idx if idx in top_k_50)
        holdout_in_top_100 = sum(1 for idx in holdout_idx if idx in top_k_100)

        # Compute ranks of holdout genes
        rank_map = {gene_idx: rank for rank, gene_idx in enumerate(ranking)}
        holdout_ranks = [rank_map[idx] for idx in holdout_idx]

        metrics = {
            "recall_at_50": holdout_in_top_50 / n_holdout if n_holdout > 0 else 0.0,
            "recall_at_100": holdout_in_top_100 / n_holdout if n_holdout > 0 else 0.0,
            "n_holdout_total": float(n_holdout),
            "n_holdout_in_top_50": float(holdout_in_top_50),
            "n_holdout_in_top_100": float(holdout_in_top_100),
            "median_holdout_rank": float(np.median(holdout_ranks)) if holdout_ranks else 0.0,
        }

        return TemporalResult(
            metrics=metrics,
            n_holdout_in_top_k=holdout_in_top_50,
            ranked_scores=scores,
            holdout_ranks=holdout_ranks,
        )
