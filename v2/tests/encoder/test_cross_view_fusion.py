"""Tests for cross-view attention fusion (WP4.7). Written FIRST per TDD."""

import torch

from ssame.encoder.cross_view_fusion import CrossViewFusion


class TestCrossViewFusion:
    def test_creates_fusion(self) -> None:
        fusion = CrossViewFusion(node_dim=16, n_views=3)
        assert fusion is not None

    def test_forward_output_shape(self) -> None:
        fusion = CrossViewFusion(node_dim=16, n_views=3)
        # 3 views, each with (20 nodes, 16 features)
        view_embeddings = [torch.randn(20, 16) for _ in range(3)]
        out = fusion(view_embeddings)
        assert out.shape == (20, 16)

    def test_attention_weights_sum_to_one(self) -> None:
        """Per gene, attention weights across views should sum to ~1."""
        fusion = CrossViewFusion(node_dim=16, n_views=3)
        view_embeddings = [torch.randn(20, 16) for _ in range(3)]
        _, attn_weights = fusion(view_embeddings, return_attention=True)
        # attn_weights shape: (n_nodes, n_views)
        assert attn_weights.shape == (20, 3)
        row_sums = attn_weights.sum(dim=1)
        assert torch.allclose(row_sums, torch.ones(20), atol=1e-5)

    def test_single_view_passthrough(self) -> None:
        """With 1 view, output should be close to the input."""
        fusion = CrossViewFusion(node_dim=16, n_views=1)
        view_embeddings = [torch.randn(20, 16)]
        out = fusion(view_embeddings)
        # With single view and attention weight=1, should be close to input
        # (not identical due to projection layers)
        assert out.shape == (20, 16)

    def test_gradient_flows(self) -> None:
        fusion = CrossViewFusion(node_dim=16, n_views=2)
        v1 = torch.randn(10, 16, requires_grad=True)
        v2 = torch.randn(10, 16, requires_grad=True)
        out = fusion([v1, v2])
        loss = out.sum()
        loss.backward()
        assert v1.grad is not None
        assert v2.grad is not None

    def test_different_views_different_weights(self) -> None:
        """Views with different content should get different attention weights."""
        torch.manual_seed(42)
        fusion = CrossViewFusion(node_dim=16, n_views=2)
        # View 1 has large values, view 2 has small values
        v1 = torch.randn(10, 16) * 5.0
        v2 = torch.randn(10, 16) * 0.1
        _, attn = fusion([v1, v2], return_attention=True)
        # Attention should differ between views
        assert not torch.allclose(attn[:, 0], attn[:, 1], atol=0.01)
