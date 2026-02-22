"""Tests for Integrated Gradients explanations (WP9.1). Written FIRST per TDD."""

import numpy as np
import pytest
import torch

from ssame.explain.integrated_gradients import IntegratedGradients


@pytest.fixture
def simple_model():
    """A simple linear model where feature 0 matters most."""
    model = torch.nn.Linear(5, 1, bias=False)
    with torch.no_grad():
        model.weight.copy_(torch.tensor([[3.0, 0.5, 0.1, 0.0, -0.2]]))
    return model


@pytest.fixture
def ig(simple_model):
    return IntegratedGradients(model=simple_model)


class TestIntegratedGradients:
    def test_returns_attributions(self, ig) -> None:
        x = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0]])
        attrs = ig.attribute(x)
        assert isinstance(attrs, np.ndarray)
        assert attrs.shape == (5,)

    def test_feature_0_highest_attribution(self, ig) -> None:
        """Feature 0 has weight 3.0 — should have highest attribution."""
        x = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0]])
        attrs = ig.attribute(x)
        assert np.argmax(np.abs(attrs)) == 0

    def test_zero_input_zero_attributions(self, ig) -> None:
        """Zero input should give zero attributions."""
        x = torch.zeros(1, 5)
        attrs = ig.attribute(x)
        assert np.allclose(attrs, 0.0, atol=1e-6)

    def test_custom_baseline(self, ig) -> None:
        """Should accept a custom baseline."""
        x = torch.tensor([[2.0, 1.0, 0.0, 0.0, 0.0]])
        baseline = torch.tensor([[1.0, 0.0, 0.0, 0.0, 0.0]])
        attrs = ig.attribute(x, baseline=baseline)
        assert isinstance(attrs, np.ndarray)

    def test_n_steps_parameter(self, ig) -> None:
        """More steps should give similar but more precise attributions."""
        x = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0]])
        attrs_50 = ig.attribute(x, n_steps=50)
        attrs_200 = ig.attribute(x, n_steps=200)
        # Both should agree on which feature is most important
        assert np.argmax(np.abs(attrs_50)) == np.argmax(np.abs(attrs_200))

    def test_completeness_axiom(self, ig) -> None:
        """Sum of attributions should approximate (f(x) - f(baseline))."""
        x = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0]])
        baseline = torch.zeros(1, 5)
        attrs = ig.attribute(x, baseline=baseline, n_steps=300)

        with torch.no_grad():
            fx = ig.model(x).item()
            fb = ig.model(baseline).item()

        assert np.sum(attrs) == pytest.approx(fx - fb, abs=0.1)

    def test_batch_attributions(self, ig) -> None:
        """Should handle batch of inputs — return per-input attributions."""
        x = torch.tensor([
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
        ])
        attrs = ig.attribute_batch(x)
        assert attrs.shape == (2, 5)
        # First input: feature 0 should dominate
        assert np.argmax(np.abs(attrs[0])) == 0
        # Second input: feature 1 should dominate
        assert np.argmax(np.abs(attrs[1])) == 1
