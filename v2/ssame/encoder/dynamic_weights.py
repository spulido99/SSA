"""Kendall uncertainty-based dynamic multi-task weighting."""

from __future__ import annotations

import torch
from torch import nn


class DynamicWeights(nn.Module):
    """Learnable task weights based on homoscedastic uncertainty.

    Based on Kendall et al. "Multi-Task Learning Using Uncertainty to
    Weigh Losses for Scene Geometry and Semantics" (CVPR 2018).

    Each task i has a learnable log-variance parameter log_sigma_i.
    The effective weight is 1/(2*sigma_i^2), plus a regularization term
    log(sigma_i) that prevents all weights from going to zero.
    """

    def __init__(self, n_tasks: int) -> None:
        super().__init__()
        # Initialize log-variance parameters to 0 (equal weights)
        self.log_sigmas = nn.Parameter(torch.zeros(n_tasks))

    def forward(self, losses: list[torch.Tensor]) -> torch.Tensor:
        """Compute weighted sum of losses using learned uncertainties.

        Args:
            losses: List of scalar loss tensors, one per task.

        Returns:
            Weighted total loss.
        """
        total = torch.tensor(0.0, device=self.log_sigmas.device)
        for i, loss in enumerate(losses):
            # weight_i = 1 / (2 * sigma_i^2) = exp(-2 * log_sigma_i) / 2
            precision = torch.exp(-2 * self.log_sigmas[i])
            total = total + precision * loss + self.log_sigmas[i]
        return total

    def get_weights(self) -> list[float]:
        """Return current effective weights (for logging)."""
        with torch.no_grad():
            precisions = torch.exp(-2 * self.log_sigmas)
            return precisions.tolist()
