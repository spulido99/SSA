"""Tests for dynamic weight adjustment (WP4.3). Written FIRST per TDD."""

import pytest
import torch

from ssame.encoder.dynamic_weights import DynamicWeights


class TestDynamicWeights:
    def test_creates_dynamic_weights(self) -> None:
        dw = DynamicWeights(n_tasks=3)
        assert dw is not None

    def test_initial_weights_equal(self) -> None:
        dw = DynamicWeights(n_tasks=4)
        weights = dw.get_weights()
        assert len(weights) == 4
        # All initial weights should be equal
        for w in weights:
            assert w == pytest.approx(weights[0], abs=1e-6)

    def test_weighted_loss_computation(self) -> None:
        """Weighted sum of losses should work."""
        dw = DynamicWeights(n_tasks=3)
        losses = [torch.tensor(1.0), torch.tensor(2.0), torch.tensor(3.0)]
        total = dw(losses)
        assert total.item() > 0

    def test_gradient_flows_through_weights(self) -> None:
        dw = DynamicWeights(n_tasks=2)
        x = torch.tensor(1.0, requires_grad=True)
        losses = [x * 2, x * 3]
        total = dw(losses)
        total.backward()
        assert x.grad is not None

    def test_weights_are_learnable(self) -> None:
        """Log-variance parameters should be in module parameters."""
        dw = DynamicWeights(n_tasks=3)
        params = list(dw.parameters())
        assert len(params) > 0

    def test_all_losses_receive_gradient(self) -> None:
        """Each loss component should contribute to the total."""
        dw = DynamicWeights(n_tasks=3)
        xs = [torch.tensor(float(i), requires_grad=True) for i in range(1, 4)]
        losses = [x**2 for x in xs]
        total = dw(losses)
        total.backward()
        for x in xs:
            assert x.grad is not None
            assert x.grad.abs() > 0

    def test_weights_adjust_during_training(self) -> None:
        """After optimizer steps, weights should change."""
        dw = DynamicWeights(n_tasks=2)
        optimizer = torch.optim.Adam(dw.parameters(), lr=0.1)
        initial_weights = dw.get_weights().copy()

        for _ in range(10):
            # Loss 0 is always large, loss 1 is small
            losses = [torch.tensor(10.0), torch.tensor(0.1)]
            total = dw(losses)
            optimizer.zero_grad()
            total.backward()
            optimizer.step()

        final_weights = dw.get_weights()
        # Weights should have changed
        assert initial_weights != pytest.approx(final_weights, abs=1e-3)
