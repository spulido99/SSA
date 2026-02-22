"""Tests for PU-Ranking loss (WP4.1). Written FIRST per TDD."""

import torch

from ssame.encoder.pu_ranking_loss import PURankingLoss


class TestPURankingLoss:
    def test_creates_loss(self) -> None:
        loss_fn = PURankingLoss(prior_pi=0.05)
        assert loss_fn.prior_pi == 0.05

    def test_forward_returns_scalar(self) -> None:
        loss_fn = PURankingLoss(prior_pi=0.05)
        scores = torch.sigmoid(torch.randn(100))
        labels = torch.zeros(100)
        labels[:10] = 1  # 10 known positives
        loss = loss_fn(scores, labels)
        assert loss.dim() == 0  # scalar
        assert loss.item() > 0

    def test_gradient_flows(self) -> None:
        loss_fn = PURankingLoss(prior_pi=0.05)
        raw = torch.randn(100, requires_grad=True)
        scores = torch.sigmoid(raw)
        labels = torch.zeros(100)
        labels[:10] = 1
        loss = loss_fn(scores, labels)
        loss.backward()
        # Check gradient on the leaf tensor (raw), not the non-leaf (scores)
        assert raw.grad is not None
        assert not torch.isnan(raw.grad).any()

    def test_known_drivers_score_higher_after_optimization(self) -> None:
        """After gradient steps, known positives should have higher scores."""
        loss_fn = PURankingLoss(prior_pi=0.05)
        # 50 genes: 5 known positives, 45 unlabeled
        raw_scores = torch.randn(50, requires_grad=True)
        labels = torch.zeros(50)
        labels[:5] = 1

        optimizer = torch.optim.Adam([raw_scores], lr=0.1)
        for _ in range(100):
            optimizer.zero_grad()
            scores = torch.sigmoid(raw_scores)
            loss = loss_fn(scores, labels)
            loss.backward()
            optimizer.step()

        final_scores = torch.sigmoid(raw_scores).detach()
        pos_mean = final_scores[:5].mean()
        unl_mean = final_scores[5:].mean()
        assert pos_mean > unl_mean, (
            f"Positive mean {pos_mean:.3f} should be > unlabeled mean {unl_mean:.3f}"
        )

    def test_all_positive_does_not_crash(self) -> None:
        loss_fn = PURankingLoss(prior_pi=0.05)
        scores = torch.sigmoid(torch.randn(10))
        labels = torch.ones(10)
        loss = loss_fn(scores, labels)
        assert not torch.isnan(loss)

    def test_all_unlabeled_does_not_crash(self) -> None:
        loss_fn = PURankingLoss(prior_pi=0.05)
        scores = torch.sigmoid(torch.randn(10))
        labels = torch.zeros(10)
        loss = loss_fn(scores, labels)
        assert not torch.isnan(loss)

    def test_pu_risk_is_nonnegative(self) -> None:
        """PU risk estimator should be non-negative (due to max(0, ...) correction)."""
        loss_fn = PURankingLoss(prior_pi=0.05)
        torch.manual_seed(42)
        for _ in range(10):
            scores = torch.sigmoid(torch.randn(100))
            labels = torch.zeros(100)
            labels[:10] = 1
            loss = loss_fn(scores, labels)
            assert loss.item() >= 0, f"Loss should be non-negative, got {loss.item()}"

    def test_custom_prior(self) -> None:
        loss1 = PURankingLoss(prior_pi=0.01)
        loss2 = PURankingLoss(prior_pi=0.1)
        scores = torch.sigmoid(torch.randn(100))
        labels = torch.zeros(100)
        labels[:10] = 1
        l1 = loss1(scores, labels)
        l2 = loss2(scores, labels)
        # Different priors should give different loss values
        assert l1.item() != l2.item()
