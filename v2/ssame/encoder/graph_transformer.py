"""Graph Transformer layer with multi-head attention."""

from __future__ import annotations

import math

import torch
from torch import nn


class GraphTransformerLayer(nn.Module):
    """Single Graph Transformer layer.

    Implements multi-head attention over graph neighborhoods with
    residual connections and LayerNorm, similar to the GPS framework.
    Uses sparse message passing (only attends to neighbors in edge_index).
    """

    def __init__(
        self,
        in_dim: int,
        out_dim: int,
        n_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.in_dim = in_dim
        self.out_dim = out_dim
        self.n_heads = n_heads
        self.head_dim = out_dim // n_heads

        # Query, Key, Value projections
        self.q_proj = nn.Linear(in_dim, out_dim)
        self.k_proj = nn.Linear(in_dim, out_dim)
        self.v_proj = nn.Linear(in_dim, out_dim)
        self.out_proj = nn.Linear(out_dim, out_dim)

        # Feed-forward network
        self.ffn = nn.Sequential(
            nn.Linear(out_dim, out_dim * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(out_dim * 2, out_dim),
            nn.Dropout(dropout),
        )

        self.norm1 = nn.LayerNorm(out_dim)
        self.norm2 = nn.LayerNorm(out_dim)
        self.dropout = nn.Dropout(dropout)

        # Residual projection if dimensions differ
        self.residual_proj = (
            nn.Linear(in_dim, out_dim) if in_dim != out_dim else nn.Identity()
        )

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        """Forward pass with sparse multi-head attention.

        Args:
            x: Node features (n_nodes, in_dim).
            edge_index: Edge indices in COO format (2, n_edges).

        Returns:
            Updated node features (n_nodes, out_dim).
        """
        n_nodes = x.size(0)
        residual = self.residual_proj(x)

        # Compute Q, K, V
        q = self.q_proj(x)  # (n_nodes, out_dim)
        k = self.k_proj(x)
        v = self.v_proj(x)

        # Reshape for multi-head: (n_nodes, n_heads, head_dim)
        q = q.view(n_nodes, self.n_heads, self.head_dim)
        k = k.view(n_nodes, self.n_heads, self.head_dim)
        v = v.view(n_nodes, self.n_heads, self.head_dim)

        # Sparse attention: compute attention only between connected nodes
        src, dst = edge_index[0], edge_index[1]

        # Attention scores for each edge
        q_dst = q[dst]  # (n_edges, n_heads, head_dim)
        k_src = k[src]  # (n_edges, n_heads, head_dim)
        attn_scores = (q_dst * k_src).sum(dim=-1) / math.sqrt(self.head_dim)
        # (n_edges, n_heads)

        # Softmax over incoming edges per destination node
        attn_weights = self._sparse_softmax(attn_scores, dst, n_nodes)
        attn_weights = self.dropout(attn_weights)

        # Weighted aggregation of values
        v_src = v[src]  # (n_edges, n_heads, head_dim)
        weighted_v = v_src * attn_weights.unsqueeze(-1)  # (n_edges, n_heads, head_dim)

        # Scatter-add to destination nodes
        out = torch.zeros(n_nodes, self.n_heads, self.head_dim, device=x.device)
        out.scatter_add_(0, dst.unsqueeze(-1).unsqueeze(-1).expand_as(weighted_v), weighted_v)

        # Reshape back: (n_nodes, out_dim)
        out = out.view(n_nodes, self.out_dim)
        out = self.out_proj(out)

        # Add & Norm (attention block)
        out = self.norm1(residual + out)

        # FFN block with residual
        out = self.norm2(out + self.ffn(out))

        return out

    @staticmethod
    def _sparse_softmax(
        scores: torch.Tensor,
        indices: torch.Tensor,
        n_nodes: int,
    ) -> torch.Tensor:
        """Compute softmax over groups defined by indices.

        Args:
            scores: Edge attention scores (n_edges, n_heads).
            indices: Destination node indices (n_edges,).
            n_nodes: Total number of nodes.

        Returns:
            Softmax-normalized scores (n_edges, n_heads).
        """
        # Subtract max for numerical stability (per destination node)
        max_scores = torch.zeros(n_nodes, scores.size(1), device=scores.device)
        max_scores.scatter_reduce_(
            0,
            indices.unsqueeze(-1).expand_as(scores),
            scores,
            reduce="amax",
            include_self=True,
        )
        scores = scores - max_scores[indices]

        exp_scores = scores.exp()

        # Sum of exp scores per destination node
        sum_exp = torch.zeros(n_nodes, scores.size(1), device=scores.device)
        sum_exp.scatter_add_(
            0,
            indices.unsqueeze(-1).expand_as(exp_scores),
            exp_scores,
        )

        # Normalize
        return exp_scores / (sum_exp[indices] + 1e-10)
