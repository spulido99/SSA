"""Tests for biological prior regularizers (WP4.2). Written FIRST per TDD."""

import numpy as np
import pytest
import torch

from ssame.encoder.biological_priors import BiologicalPriorLoss


@pytest.fixture
def toy_alteration_data() -> tuple[torch.Tensor, torch.Tensor]:
    """Alteration matrix (10 genes, 30 samples) and predicted scores."""
    # ME pattern: genes 0-1 are mutually exclusive
    alt = np.zeros((10, 30), dtype=np.float32)
    alt[0, :15] = 1  # Gene 0 in first 15 samples
    alt[1, 15:25] = 1  # Gene 1 in next 10 samples (ME with gene 0)
    alt[2, :5] = 1  # Gene 2 overlaps with gene 0
    alt[3, 12:18] = 1  # Gene 3 overlaps with both

    scores = torch.tensor(
        [0.9, 0.8, 0.5, 0.4, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        requires_grad=True,
    )
    return torch.from_numpy(alt), scores


class TestBiologicalPriorLoss:
    def test_creates_loss(self) -> None:
        alt = torch.zeros(10, 30)
        loss_fn = BiologicalPriorLoss(alteration_matrix=alt)
        assert loss_fn is not None

    def test_me_regularizer_gradient_flows(
        self, toy_alteration_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        alt, scores = toy_alteration_data
        loss_fn = BiologicalPriorLoss(alteration_matrix=alt)
        loss = loss_fn(scores)
        loss.backward()
        assert scores.grad is not None
        assert scores.grad.abs().sum() > 0

    def test_me_penalty_for_overlapping_high_scores(
        self, toy_alteration_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        """Co-occurring high-scoring genes should incur penalty."""
        alt, _ = toy_alteration_data
        # Scores where co-occurring genes both score high
        bad_scores = torch.tensor(
            [0.9, 0.9, 0.9, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            requires_grad=True,
        )
        # Scores where ME genes score high, co-occurring genes score low
        good_scores = torch.tensor(
            [0.9, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            requires_grad=True,
        )
        loss_fn = BiologicalPriorLoss(alteration_matrix=alt)
        loss_bad = loss_fn(good_scores)
        loss_good = loss_fn(bad_scores)
        # ME-violating pattern should have higher loss
        assert loss_good > loss_bad

    def test_loss_is_non_negative(
        self, toy_alteration_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        alt, scores = toy_alteration_data
        loss_fn = BiologicalPriorLoss(alteration_matrix=alt)
        loss = loss_fn(scores)
        assert loss.item() >= 0

    def test_weight_parameter(
        self, toy_alteration_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        """Weight should scale the loss."""
        alt, scores = toy_alteration_data
        loss_1 = BiologicalPriorLoss(alteration_matrix=alt, weight=1.0)(scores)
        scores2 = scores.detach().clone().requires_grad_(True)
        loss_2 = BiologicalPriorLoss(alteration_matrix=alt, weight=2.0)(scores2)
        assert abs(loss_2.item() - 2 * loss_1.item()) < 1e-5

    def test_empty_alteration_returns_zero(self) -> None:
        """No alterations should produce zero regularizer."""
        alt = torch.zeros(5, 20)
        scores = torch.tensor([0.5, 0.5, 0.5, 0.5, 0.5], requires_grad=True)
        loss_fn = BiologicalPriorLoss(alteration_matrix=alt)
        loss = loss_fn(scores)
        assert loss.item() == pytest.approx(0.0, abs=1e-6)
