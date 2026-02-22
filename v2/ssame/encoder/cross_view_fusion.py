"""Cross-view attention fusion for multi-network representations."""

from __future__ import annotations

import torch
from torch import nn


class CrossViewFusion(nn.Module):
    """Fuses node representations from multiple network views using attention.

    For each gene, computes attention weights over views, then produces
    a weighted combination of view-specific embeddings.
    """

    def __init__(
        self,
        node_dim: int,
        n_views: int,
    ) -> None:
        super().__init__()
        self.node_dim = node_dim
        self.n_views = n_views

        # Attention scoring: project each view embedding to a scalar
        self.attn_proj = nn.Sequential(
            nn.Linear(node_dim, node_dim // 2),
            nn.Tanh(),
            nn.Linear(node_dim // 2, 1),
        )

        # Output projection
        self.out_proj = nn.Linear(node_dim, node_dim)

    def forward(
        self,
        view_embeddings: list[torch.Tensor],
        return_attention: bool = False,
    ) -> torch.Tensor | tuple[torch.Tensor, torch.Tensor]:
        """Fuse multi-view embeddings.

        Args:
            view_embeddings: List of (n_nodes, node_dim) tensors, one per view.
            return_attention: If True, also return attention weights.

        Returns:
            Fused embeddings (n_nodes, node_dim).
            If return_attention: tuple of (fused, attention_weights).
        """
        # Stack: (n_nodes, n_views, node_dim)
        stacked = torch.stack(view_embeddings, dim=1)

        # Compute attention scores per node per view
        attn_scores = self.attn_proj(stacked).squeeze(-1)  # (n_nodes, n_views)
        attn_weights = torch.softmax(attn_scores, dim=1)  # (n_nodes, n_views)

        # Weighted combination
        # (n_nodes, n_views, 1) * (n_nodes, n_views, node_dim) -> sum over views
        fused = (attn_weights.unsqueeze(-1) * stacked).sum(dim=1)  # (n_nodes, node_dim)

        fused = self.out_proj(fused)

        if return_attention:
            return fused, attn_weights
        return fused
