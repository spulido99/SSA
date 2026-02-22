"""Integrated Gradients for feature attribution (WP9.1).

Implementation of the Integrated Gradients method (Sundararajan et al., 2017)
for explaining neural network predictions at the gene level.
"""

from __future__ import annotations

import numpy as np
import torch


class IntegratedGradients:
    """Compute feature attributions using Integrated Gradients.

    For each input feature, computes the integral of gradients along
    the path from a baseline (default: zero) to the input.

    Args:
        model: A PyTorch model (callable) that takes input tensor and
            returns scalar or (batch, 1) output.
    """

    def __init__(self, model: torch.nn.Module) -> None:
        self.model = model

    def attribute(
        self,
        x: torch.Tensor,
        *,
        baseline: torch.Tensor | None = None,
        n_steps: int = 50,
        target_idx: int = 0,
    ) -> np.ndarray:
        """Compute attributions for a single input.

        Args:
            x: Input tensor of shape (1, n_features).
            baseline: Baseline tensor (default: zeros).
            n_steps: Number of interpolation steps.
            target_idx: Output index if model returns multi-dim output.

        Returns:
            Attributions array of shape (n_features,).
        """
        if baseline is None:
            baseline = torch.zeros_like(x)

        # Generate interpolated inputs along the path
        alphas = torch.linspace(0, 1, n_steps + 1, device=x.device)
        # Shape: (n_steps+1, n_features)
        delta = x - baseline
        interpolated = baseline + alphas.unsqueeze(1) * delta

        interpolated.requires_grad_(True)

        # Forward pass
        self.model.eval()
        outputs = self.model(interpolated)
        if outputs.dim() > 1:
            outputs = outputs[:, target_idx]

        # Backward pass to get gradients
        grads = torch.autograd.grad(
            outputs=outputs.sum(),
            inputs=interpolated,
            create_graph=False,
        )[0]

        # Trapezoidal integration
        avg_grads = (grads[:-1] + grads[1:]).mean(dim=0) / 2.0
        # But simpler: average all gradients (Riemann sum)
        avg_grads = grads.mean(dim=0)

        # Scale by (input - baseline)
        attributions = (avg_grads * delta.squeeze(0)).detach().cpu().numpy()

        return attributions

    def attribute_batch(
        self,
        x: torch.Tensor,
        *,
        baseline: torch.Tensor | None = None,
        n_steps: int = 50,
    ) -> np.ndarray:
        """Compute attributions for a batch of inputs.

        Args:
            x: Input tensor of shape (batch, n_features).
            baseline: Baseline tensor (default: zeros).
            n_steps: Number of interpolation steps.

        Returns:
            Attributions array of shape (batch, n_features).
        """
        results = []
        for i in range(x.shape[0]):
            xi = x[i : i + 1]
            bi = baseline[i : i + 1] if baseline is not None else None
            attrs = self.attribute(xi, baseline=bi, n_steps=n_steps)
            results.append(attrs)
        return np.stack(results)
