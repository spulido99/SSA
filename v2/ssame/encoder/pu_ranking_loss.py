"""PU-Ranking loss: combines PU learning with ranking objective."""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class PURankingLoss(nn.Module):
    """Positive-Unlabeled ranking loss for cancer driver prediction.

    Combines:
    - PU risk estimator (unbiased, with known prior pi)
    - Pairwise ranking: known drivers should score higher than unlabeled genes

    Labels convention: 1 = known driver, 0 = unlabeled (NOT negative).
    """

    def __init__(self, prior_pi: float = 0.05, ranking_weight: float = 1.0) -> None:
        super().__init__()
        self.prior_pi = prior_pi
        self.ranking_weight = ranking_weight

    def forward(self, scores: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
        """Compute combined PU + ranking loss.

        Args:
            scores: Predicted scores in [0, 1], shape (n_genes,).
            labels: 1 for known positives, 0 for unlabeled.

        Returns:
            Scalar loss value.
        """
        pos_mask = labels == 1
        unl_mask = labels == 0

        n_pos = pos_mask.sum().float()
        n_unl = unl_mask.sum().float()

        loss = torch.tensor(0.0, device=scores.device)

        # --- PU risk estimator ---
        if n_pos > 0:
            pos_scores = scores[pos_mask]
            # R+(f): risk on positives
            pos_risk = F.binary_cross_entropy(
                pos_scores, torch.ones_like(pos_scores), reduction="mean"
            )
            # R-(f) on positives: treating positives as negative
            neg_on_pos = F.binary_cross_entropy(
                pos_scores, torch.zeros_like(pos_scores), reduction="mean"
            )

            pu_loss = self.prior_pi * pos_risk

            if n_unl > 0:
                unl_scores = scores[unl_mask]
                # R-(f) on unlabeled
                unl_risk = F.binary_cross_entropy(
                    unl_scores, torch.zeros_like(unl_scores), reduction="mean"
                )
                # Non-negative PU correction
                pu_loss = pu_loss + F.relu(unl_risk - self.prior_pi * neg_on_pos)

            loss = loss + pu_loss

        # --- Pairwise ranking loss ---
        if n_pos > 0 and n_unl > 0:
            pos_scores = scores[pos_mask]
            unl_scores = scores[unl_mask]

            # Sample pairs for efficiency
            n_pairs = min(int(n_pos.item()) * 10, int(n_unl.item()))
            if n_pairs > 0:
                idx_pos = torch.randint(0, int(n_pos.item()), (n_pairs,), device=scores.device)
                idx_unl = torch.randint(0, int(n_unl.item()), (n_pairs,), device=scores.device)

                rank_loss = F.margin_ranking_loss(
                    pos_scores[idx_pos],
                    unl_scores[idx_unl],
                    torch.ones(n_pairs, device=scores.device),
                    margin=0.1,
                )
                loss = loss + self.ranking_weight * rank_loss

        return loss
