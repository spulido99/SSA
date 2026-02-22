"""Full SSA-ME v2 encoder (Level 1): bimodal encoding -> graph transformer -> PU ranking head."""

from __future__ import annotations

import torch
from torch import nn

from ssame.encoder.bimodal_encoder import BimodalEncoder
from ssame.encoder.biological_priors import BiologicalPriorLoss
from ssame.encoder.cross_view_fusion import CrossViewFusion
from ssame.encoder.dynamic_weights import DynamicWeights
from ssame.encoder.graph_transformer import GraphTransformerLayer
from ssame.encoder.pu_ranking_loss import PURankingLoss


class SSAEncoder(nn.Module):
    """Full encoder: features -> bimodal MLP -> graph transformer -> fusion -> scores.

    Level 1 encoder that combines:
    1. Bimodal encoder (omics MLP + topology MLP)
    2. Per-view Graph Transformer layers
    3. Cross-view attention fusion
    4. PU-ranking scoring head
    """

    def __init__(
        self,
        omics_dim: int,
        topo_dim: int,
        n_views: int = 1,
        hidden_dim: int = 64,
        n_transformer_layers: int = 2,
        n_heads: int = 4,
        dropout: float = 0.1,
        pu_prior: float = 0.05,
        alteration_matrix: torch.Tensor | None = None,
        bio_prior_weight: float = 0.1,
    ) -> None:
        super().__init__()
        self.n_views = n_views
        self.hidden_dim = hidden_dim

        # Step 1: Bimodal encoder
        self.bimodal = BimodalEncoder(
            omics_dim=omics_dim,
            topo_dim=topo_dim,
            hidden_dim=hidden_dim,
            output_dim=hidden_dim,
            dropout=dropout,
        )

        # Step 2: Per-view Graph Transformer stacks
        self.view_transformers = nn.ModuleList()
        for _ in range(n_views):
            layers = nn.ModuleList()
            for _ in range(n_transformer_layers):
                layers.append(
                    GraphTransformerLayer(
                        in_dim=hidden_dim,
                        out_dim=hidden_dim,
                        n_heads=n_heads,
                        dropout=dropout,
                    )
                )
            self.view_transformers.append(layers)

        # Step 3: Cross-view fusion
        self.fusion = CrossViewFusion(
            node_dim=hidden_dim,
            n_views=n_views,
        )

        # Step 4: Scoring head
        self.score_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

        # Losses
        self.pu_loss = PURankingLoss(prior_pi=pu_prior)
        self.dynamic_weights = DynamicWeights(
            n_tasks=2 if alteration_matrix is not None else 1
        )

        if alteration_matrix is not None:
            self.bio_prior_loss = BiologicalPriorLoss(
                alteration_matrix=alteration_matrix,
                weight=bio_prior_weight,
            )
        else:
            self.bio_prior_loss = None

    def forward(
        self,
        omics_features: torch.Tensor,
        topo_features: torch.Tensor,
        edge_indices: list[torch.Tensor],
    ) -> torch.Tensor:
        """Forward pass producing per-gene driverness scores.

        Args:
            omics_features: (n_genes, omics_dim).
            topo_features: (n_genes, topo_dim).
            edge_indices: List of edge_index tensors (2, n_edges), one per view.

        Returns:
            Scores in [0, 1] with shape (n_genes,).
        """
        # Bimodal encoding
        h = self.bimodal(omics_features, topo_features)

        # Per-view graph transformer
        view_outputs = []
        for view_idx, transformer_layers in enumerate(self.view_transformers):
            h_view = h
            edge_index = edge_indices[view_idx]
            for layer in transformer_layers:
                h_view = layer(h_view, edge_index)
            view_outputs.append(h_view)

        # Cross-view fusion
        fused = self.fusion(view_outputs)

        # Scoring head
        logits = self.score_head(fused).squeeze(-1)
        scores = torch.sigmoid(logits)

        return scores

    def compute_loss(
        self,
        scores: torch.Tensor,
        labels: torch.Tensor,
    ) -> torch.Tensor:
        """Compute combined loss (PU ranking + optional biological priors).

        Args:
            scores: Predicted scores (n_genes,).
            labels: PU labels (1=known positive, 0=unlabeled).

        Returns:
            Scalar loss.
        """
        losses = [self.pu_loss(scores, labels)]

        if self.bio_prior_loss is not None:
            losses.append(self.bio_prior_loss(scores))

        return self.dynamic_weights(losses)
