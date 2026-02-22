"""Biological prior regularizers (differentiable)."""

from __future__ import annotations

import torch
from torch import nn


class BiologicalPriorLoss(nn.Module):
    """Mutual exclusivity regularizer for predicted gene scores.

    Penalizes pairs of genes that co-occur in the same samples (overlap)
    when both have high predicted scores. This encourages the encoder to
    favor mutually exclusive gene sets, consistent with pathway-level
    driver mutation patterns.
    """

    def __init__(
        self,
        alteration_matrix: torch.Tensor,
        weight: float = 1.0,
    ) -> None:
        super().__init__()
        # alteration_matrix: (n_genes, n_samples), binary
        self.register_buffer("alt", alteration_matrix.float())
        self.weight = weight

        # Pre-compute co-occurrence matrix: (n_genes, n_genes)
        # co_occur[i,j] = number of samples where both gene i and j are altered
        self.register_buffer("co_occur", self.alt @ self.alt.T)

    def forward(self, scores: torch.Tensor) -> torch.Tensor:
        """Compute ME regularization loss.

        Args:
            scores: Predicted gene driverness scores (n_genes,).

        Returns:
            Scalar loss penalizing co-occurrence among high-scoring genes.
        """
        # Outer product of scores: weight co-occurrence by both genes' scores
        score_product = scores.unsqueeze(0) * scores.unsqueeze(1)

        # Penalty: sum of (score_i * score_j * co_occur_ij) for i != j
        # This penalizes high scores on gene pairs that co-occur
        n = scores.shape[0]
        mask = 1.0 - torch.eye(n, device=scores.device)
        penalty = (score_product * self.co_occur * mask).sum()

        # Normalize by number of pairs and number of samples
        n_samples = self.alt.shape[1]
        n_pairs = n * (n - 1)
        if n_pairs > 0 and n_samples > 0:
            penalty = penalty / (n_pairs * n_samples)

        return self.weight * penalty
