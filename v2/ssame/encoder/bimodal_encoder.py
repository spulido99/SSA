"""Bimodal encoder: separate MLPs for omics and topology features."""

from __future__ import annotations

import torch
from torch import nn


class BimodalEncoder(nn.Module):
    """Encodes omics and topology features through independent MLPs, then concatenates.

    SGCD-inspired: omics features (mutation freq, essentiality, etc.) and topology
    features (degree, betweenness, node2vec) are encoded separately to prevent
    one modality from dominating the other.
    """

    def __init__(
        self,
        omics_dim: int,
        topo_dim: int,
        hidden_dim: int = 64,
        output_dim: int = 32,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.omics_dim = omics_dim
        self.topo_dim = topo_dim

        if omics_dim > 0:
            self.omics_mlp = nn.Sequential(
                nn.Linear(omics_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
            )
            omics_out = hidden_dim
        else:
            self.omics_mlp = None
            omics_out = 0

        if topo_dim > 0:
            self.topo_mlp = nn.Sequential(
                nn.Linear(topo_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
            )
            topo_out = hidden_dim
        else:
            self.topo_mlp = None
            topo_out = 0

        self.projection = nn.Linear(omics_out + topo_out, output_dim)

    def forward(
        self,
        omics: torch.Tensor | None,
        topo: torch.Tensor | None,
    ) -> torch.Tensor:
        """Encode features through bimodal MLPs.

        Args:
            omics: Omics features (n_genes, omics_dim) or None.
            topo: Topology features (n_genes, topo_dim) or None.

        Returns:
            Encoded representation (n_genes, output_dim).
        """
        parts = []
        if self.omics_mlp is not None and omics is not None:
            parts.append(self.omics_mlp(omics))
        if self.topo_mlp is not None and topo is not None:
            parts.append(self.topo_mlp(topo))

        combined = torch.cat(parts, dim=-1)
        return self.projection(combined)
