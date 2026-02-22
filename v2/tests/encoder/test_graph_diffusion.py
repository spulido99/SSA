"""Tests for PPR graph diffusion (WP4.8). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.encoder.graph_diffusion import PPRDiffusion


@pytest.fixture
def toy_adj() -> np.ndarray:
    """5-node graph adjacency matrix (chain: 0-1-2-3-4)."""
    adj = np.zeros((5, 5))
    for i in range(4):
        adj[i, i + 1] = 1
        adj[i + 1, i] = 1
    return adj


@pytest.fixture
def star_adj() -> np.ndarray:
    """Star graph: node 0 connected to all others."""
    n = 6
    adj = np.zeros((n, n))
    for i in range(1, n):
        adj[0, i] = 1
        adj[i, 0] = 1
    return adj


class TestPPRDiffusion:
    def test_creates_diffusion(self, toy_adj: np.ndarray) -> None:
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        assert ppr is not None

    def test_seed_node_highest_score(self, toy_adj: np.ndarray) -> None:
        """Seed node should have highest PPR score to itself (center of chain)."""
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        # Use center node (index 2) — endpoint nodes can lose to higher-degree
        # neighbors at low alpha, which is correct PPR behavior
        scores = ppr.diffuse_from_seed(seed_idx=2)
        assert scores[2] == scores.max()

    def test_nearby_nodes_higher(self, toy_adj: np.ndarray) -> None:
        """In a chain, closer nodes should have higher PPR scores."""
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        scores = ppr.diffuse_from_seed(seed_idx=0)
        # Node 1 should score higher than node 4 (further away)
        assert scores[1] > scores[4]

    def test_all_scores_positive(self, toy_adj: np.ndarray) -> None:
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        scores = ppr.diffuse_from_seed(seed_idx=2)
        assert (scores >= 0).all()

    def test_scores_sum_to_one(self, toy_adj: np.ndarray) -> None:
        """PPR scores should sum to approximately 1."""
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        scores = ppr.diffuse_from_seed(seed_idx=0)
        assert abs(scores.sum() - 1.0) < 1e-6

    def test_diffuse_signal(self, toy_adj: np.ndarray) -> None:
        """Diffuse a signal vector across the graph."""
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        signal = np.array([0.0, 0.0, 1.0, 0.0, 0.0])
        diffused = ppr.diffuse_signal(signal)
        assert diffused.shape == signal.shape
        # Center node should still have highest value
        assert diffused[2] == diffused.max()
        # Signal should spread to neighbors symmetrically
        assert diffused[1] > diffused[4]  # node 1 closer to center than node 4

    def test_star_graph_symmetry(self, star_adj: np.ndarray) -> None:
        """In star graph, all leaves should have equal PPR from center."""
        ppr = PPRDiffusion(adj=star_adj, alpha=0.15)
        scores = ppr.diffuse_from_seed(seed_idx=0)
        # All leaf nodes should have same score
        leaf_scores = scores[1:]
        assert np.allclose(leaf_scores, leaf_scores[0])

    def test_high_alpha_concentrates(self, toy_adj: np.ndarray) -> None:
        """High alpha (restart probability) should concentrate score on seed."""
        ppr_low = PPRDiffusion(adj=toy_adj, alpha=0.05)
        ppr_high = PPRDiffusion(adj=toy_adj, alpha=0.85)
        scores_low = ppr_low.diffuse_from_seed(seed_idx=0)
        scores_high = ppr_high.diffuse_from_seed(seed_idx=0)
        # High alpha should give more weight to seed node
        assert scores_high[0] > scores_low[0]

    def test_diffuse_features(self, toy_adj: np.ndarray) -> None:
        """Diffuse a feature matrix (n_nodes, d)."""
        ppr = PPRDiffusion(adj=toy_adj, alpha=0.15)
        features = np.eye(5)  # One-hot features
        diffused = ppr.diffuse_features(features)
        assert diffused.shape == features.shape
        # After diffusion, features should be smoothed (no longer one-hot)
        assert not np.allclose(diffused, features)
