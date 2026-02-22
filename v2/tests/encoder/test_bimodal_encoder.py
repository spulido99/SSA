"""Tests for bimodal encoder (WP4.5). Written FIRST per TDD."""

import numpy as np
import pytest
import torch

from ssame.encoder.bimodal_encoder import BimodalEncoder


@pytest.fixture
def toy_features() -> tuple[torch.Tensor, torch.Tensor]:
    """Omics features (200 genes, 4 features) and topo features (200 genes, 8 features)."""
    rng = np.random.default_rng(42)
    omics = torch.from_numpy(rng.standard_normal((200, 4))).float()
    topo = torch.from_numpy(rng.standard_normal((200, 8))).float()
    return omics, topo


class TestBimodalEncoder:
    def test_creates_encoder(self) -> None:
        enc = BimodalEncoder(omics_dim=4, topo_dim=8, hidden_dim=32, output_dim=16)
        assert enc is not None

    def test_forward_output_shape(
        self, toy_features: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        omics, topo = toy_features
        enc = BimodalEncoder(omics_dim=4, topo_dim=8, hidden_dim=32, output_dim=16)
        out = enc(omics, topo)
        assert out.shape == (200, 16)

    def test_omics_only(self) -> None:
        """Should work with only omics features (topo=None)."""
        omics = torch.randn(50, 4)
        enc = BimodalEncoder(omics_dim=4, topo_dim=0, hidden_dim=32, output_dim=16)
        out = enc(omics, None)
        assert out.shape == (50, 16)

    def test_topo_only(self) -> None:
        """Should work with only topology features (omics=None)."""
        topo = torch.randn(50, 8)
        enc = BimodalEncoder(omics_dim=0, topo_dim=8, hidden_dim=32, output_dim=16)
        out = enc(None, topo)
        assert out.shape == (50, 16)

    def test_gradient_flows(
        self, toy_features: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        omics, topo = toy_features
        omics.requires_grad_(True)
        topo.requires_grad_(True)
        enc = BimodalEncoder(omics_dim=4, topo_dim=8, hidden_dim=32, output_dim=16)
        out = enc(omics, topo)
        loss = out.sum()
        loss.backward()
        assert omics.grad is not None
        assert topo.grad is not None

    def test_output_dimensionality_independent_of_input(self) -> None:
        """Output dim should be fixed regardless of input feature dimensions."""
        enc_small = BimodalEncoder(omics_dim=2, topo_dim=3, hidden_dim=16, output_dim=8)
        enc_large = BimodalEncoder(omics_dim=20, topo_dim=30, hidden_dim=16, output_dim=8)
        out_small = enc_small(torch.randn(10, 2), torch.randn(10, 3))
        out_large = enc_large(torch.randn(10, 20), torch.randn(10, 30))
        assert out_small.shape == out_large.shape == (10, 8)

    def test_dropout_affects_training(
        self, toy_features: tuple[torch.Tensor, torch.Tensor]
    ) -> None:
        omics, topo = toy_features
        enc = BimodalEncoder(
            omics_dim=4, topo_dim=8, hidden_dim=32, output_dim=16, dropout=0.5
        )
        enc.train()
        out1 = enc(omics, topo)
        out2 = enc(omics, topo)
        # With dropout=0.5, outputs should differ in training mode
        assert not torch.allclose(out1, out2)
