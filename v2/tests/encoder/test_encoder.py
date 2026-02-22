"""Tests for full encoder assembly (WP4.9). Written FIRST per TDD."""

import numpy as np
import pytest
import torch

from ssame.encoder.encoder import SSAEncoder


@pytest.fixture
def toy_encoder_data() -> dict:
    """Feature matrix, edge indices, and labels for a toy 20-gene problem."""
    rng = np.random.default_rng(42)
    n_genes = 20
    omics_dim = 4
    topo_dim = 6

    omics_features = torch.from_numpy(
        rng.standard_normal((n_genes, omics_dim))
    ).float()
    topo_features = torch.from_numpy(
        rng.standard_normal((n_genes, topo_dim))
    ).float()

    # Add signal: first 5 genes are positives, make them distinguishable
    omics_features[:5] += 2.0
    topo_features[:5] += 1.5

    # Simple chain graph for single view
    src = list(range(n_genes - 1))
    dst = list(range(1, n_genes))
    edge_index = torch.tensor([src + dst, dst + src], dtype=torch.long)

    labels = torch.zeros(n_genes)
    labels[:5] = 1.0

    # Alteration matrix for biological priors
    alt = np.zeros((n_genes, 30), dtype=np.float32)
    alt[0, :15] = 1
    alt[1, 15:25] = 1
    alt[2, :5] = 1

    return {
        "omics_features": omics_features,
        "topo_features": topo_features,
        "edge_indices": [edge_index],  # Single view
        "labels": labels,
        "alteration_matrix": torch.from_numpy(alt),
        "n_genes": n_genes,
        "omics_dim": omics_dim,
        "topo_dim": topo_dim,
    }


class TestSSAEncoder:
    def test_creates_encoder(self, toy_encoder_data: dict) -> None:
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
        )
        assert enc is not None

    def test_forward_produces_scores(self, toy_encoder_data: dict) -> None:
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
        )
        scores = enc(d["omics_features"], d["topo_features"], d["edge_indices"])
        assert scores.shape == (d["n_genes"],)
        # Scores should be in [0, 1] (sigmoid output)
        assert (scores >= 0).all()
        assert (scores <= 1).all()

    def test_loss_computation(self, toy_encoder_data: dict) -> None:
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
            alteration_matrix=d["alteration_matrix"],
        )
        scores = enc(d["omics_features"], d["topo_features"], d["edge_indices"])
        loss = enc.compute_loss(scores, d["labels"])
        assert loss.item() > 0

    def test_backward_pass(self, toy_encoder_data: dict) -> None:
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
        )
        scores = enc(d["omics_features"], d["topo_features"], d["edge_indices"])
        loss = enc.compute_loss(scores, d["labels"])
        loss.backward()
        # Check that gradients exist on encoder parameters
        for p in enc.parameters():
            if p.requires_grad:
                assert p.grad is not None

    def test_training_loop_converges(self, toy_encoder_data: dict) -> None:
        """Loss should decrease over training iterations."""
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
        )
        optimizer = torch.optim.Adam(enc.parameters(), lr=0.01)

        losses = []
        for _ in range(50):
            scores = enc(d["omics_features"], d["topo_features"], d["edge_indices"])
            loss = enc.compute_loss(scores, d["labels"])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            losses.append(loss.item())

        # Loss should decrease from start to end
        assert losses[-1] < losses[0]

    def test_known_drivers_score_higher_after_training(
        self, toy_encoder_data: dict
    ) -> None:
        """After training, known positives should score higher than unlabeled."""
        d = toy_encoder_data
        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=1,
            hidden_dim=16,
        )
        optimizer = torch.optim.Adam(enc.parameters(), lr=0.01)

        for _ in range(100):
            scores = enc(d["omics_features"], d["topo_features"], d["edge_indices"])
            loss = enc.compute_loss(scores, d["labels"])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        with torch.no_grad():
            final_scores = enc(
                d["omics_features"], d["topo_features"], d["edge_indices"]
            )
        pos_mean = final_scores[:5].mean()
        neg_mean = final_scores[5:].mean()
        assert pos_mean > neg_mean

    def test_multi_view(self, toy_encoder_data: dict) -> None:
        """Should work with multiple graph views."""
        d = toy_encoder_data
        # Add a second view (different random edges)
        src2 = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
        dst2 = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        edge_index2 = torch.tensor([src2 + dst2, dst2 + src2], dtype=torch.long)
        edge_indices = [d["edge_indices"][0], edge_index2]

        enc = SSAEncoder(
            omics_dim=d["omics_dim"],
            topo_dim=d["topo_dim"],
            n_views=2,
            hidden_dim=16,
        )
        scores = enc(d["omics_features"], d["topo_features"], edge_indices)
        assert scores.shape == (d["n_genes"],)
