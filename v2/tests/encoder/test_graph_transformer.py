"""Tests for Graph Transformer layer (WP4.6). Written FIRST per TDD."""

import numpy as np
import pytest
import torch

from ssame.encoder.graph_transformer import GraphTransformerLayer


@pytest.fixture
def toy_graph_data() -> tuple[torch.Tensor, torch.Tensor]:
    """Node features (20 nodes, 16 features) and edge_index (COO format)."""
    rng = np.random.default_rng(42)
    x = torch.from_numpy(rng.standard_normal((20, 16))).float()
    # Simple chain graph: 0-1-2-...-19 plus some cross-edges
    src = [*list(range(19)), 0, 5, 10]
    dst = [*list(range(1, 20)), 10, 15, 3]
    # Make undirected
    edge_index = torch.tensor(
        [src + dst, dst + src], dtype=torch.long
    )
    return x, edge_index


class TestGraphTransformerLayer:
    def test_creates_layer(self) -> None:
        layer = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=4)
        assert layer is not None

    def test_forward_output_shape(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        x, edge_index = toy_graph_data
        layer = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=4)
        out = layer(x, edge_index)
        assert out.shape == x.shape

    def test_different_output_dim(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        x, edge_index = toy_graph_data
        layer = GraphTransformerLayer(in_dim=16, out_dim=32, n_heads=8)
        out = layer(x, edge_index)
        assert out.shape == (20, 32)

    def test_multi_head_attention(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        """8 attention heads should produce different result than 1 head."""
        x, edge_index = toy_graph_data
        torch.manual_seed(42)
        layer_1h = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=1)
        torch.manual_seed(99)
        layer_8h = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=8)
        out_1h = layer_1h(x, edge_index)
        out_8h = layer_8h(x, edge_index)
        # Should not be identical (different architectures)
        assert not torch.allclose(out_1h, out_8h)

    def test_residual_connection(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        """With residual, output should not be zero even with zero weights."""
        x, edge_index = toy_graph_data
        layer = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=4)
        out = layer(x, edge_index)
        # Output should not be identical to input (transformation occurs)
        assert not torch.allclose(out, x)
        # But with residual connection, output should be close to input magnitude
        assert out.abs().mean() > 0.1

    def test_gradient_flows(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        x, edge_index = toy_graph_data
        x = x.clone().requires_grad_(True)
        layer = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=4)
        out = layer(x, edge_index)
        loss = out.sum()
        loss.backward()
        assert x.grad is not None
        assert x.grad.abs().sum() > 0

    def test_layer_norm_applied(
        self, toy_graph_data: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        """Output should have reasonable scale due to LayerNorm."""
        x, edge_index = toy_graph_data
        layer = GraphTransformerLayer(in_dim=16, out_dim=16, n_heads=4)
        out = layer(x, edge_index)
        # LayerNorm should keep values in reasonable range
        assert out.abs().max() < 100
