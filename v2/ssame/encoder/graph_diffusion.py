"""Personalized PageRank graph diffusion."""

from __future__ import annotations

import numpy as np


class PPRDiffusion:
    """Personalized PageRank diffusion on a graph.

    Computes PPR scores using the power iteration method.
    Used to smooth node features or scores across the graph topology.
    """

    def __init__(
        self,
        adj: np.ndarray,
        alpha: float = 0.15,
        n_iterations: int = 50,
    ) -> None:
        self.alpha = alpha
        self.n_iterations = n_iterations
        n = adj.shape[0]

        # Build transition matrix (row-stochastic)
        degree = adj.sum(axis=1)
        # Handle isolated nodes
        degree[degree == 0] = 1.0
        self._transition = adj / degree[:, np.newaxis]
        self._n = n

    def diffuse_from_seed(self, seed_idx: int) -> np.ndarray:
        """Compute PPR scores from a single seed node.

        Args:
            seed_idx: Index of the seed node.

        Returns:
            PPR score vector (n_nodes,) summing to ~1.
        """
        personalization = np.zeros(self._n)
        personalization[seed_idx] = 1.0
        return self._power_iteration(personalization)

    def diffuse_signal(self, signal: np.ndarray) -> np.ndarray:
        """Diffuse a signal vector across the graph.

        Args:
            signal: Signal vector (n_nodes,).

        Returns:
            Diffused signal (n_nodes,).
        """
        # Normalize signal as personalization
        s = signal.sum()
        personalization = signal / s if s > 0 else np.ones(self._n) / self._n
        ppr = self._power_iteration(personalization)
        # Scale back to original signal magnitude
        return ppr * s

    def diffuse_features(self, features: np.ndarray) -> np.ndarray:
        """Diffuse a feature matrix (n_nodes, d) column by column.

        Args:
            features: Feature matrix (n_nodes, d).

        Returns:
            Diffused features (n_nodes, d).
        """
        result = np.zeros_like(features)
        for col in range(features.shape[1]):
            result[:, col] = self.diffuse_signal(features[:, col])
        return result

    def _power_iteration(self, personalization: np.ndarray) -> np.ndarray:
        """Run power iteration for PPR computation."""
        scores = personalization.copy()
        for _ in range(self.n_iterations):
            scores = (
                (1 - self.alpha) * (self._transition.T @ scores)
                + self.alpha * personalization
            )
        return scores
