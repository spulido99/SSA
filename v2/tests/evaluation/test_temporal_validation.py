"""Tests for temporal validation (WP8.2). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.evaluation.temporal_validation import TemporalResult, TemporalValidator


@pytest.fixture
def toy_temporal_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Features, training labels (pre-2020), and holdout labels (2020-2025)."""
    rng = np.random.default_rng(42)
    n_genes = 100
    n_features = 6
    X = rng.standard_normal((n_genes, n_features))

    # 10 pre-2020 known drivers
    train_labels = np.zeros(n_genes)
    train_labels[:10] = 1
    X[:10, 0] += 2.0

    # 5 post-2020 additions (holdout)
    holdout_labels = np.zeros(n_genes)
    holdout_labels[10:15] = 1
    X[10:15, 0] += 1.5  # weaker signal

    return X, train_labels, holdout_labels


class TestTemporalValidator:
    def test_creates_validator(self) -> None:
        tv = TemporalValidator(n_estimators=10)
        assert tv is not None

    def test_run_returns_result(
        self, toy_temporal_data: tuple[np.ndarray, np.ndarray, np.ndarray]
    ) -> None:
        X, train_labels, holdout_labels = toy_temporal_data
        tv = TemporalValidator(n_estimators=10)
        result = tv.run(X, train_labels, holdout_labels)
        assert isinstance(result, TemporalResult)

    def test_result_has_recall_at_k(
        self, toy_temporal_data: tuple[np.ndarray, np.ndarray, np.ndarray]
    ) -> None:
        X, train_labels, holdout_labels = toy_temporal_data
        tv = TemporalValidator(n_estimators=20)
        result = tv.run(X, train_labels, holdout_labels)
        assert "recall_at_50" in result.metrics
        assert "recall_at_100" in result.metrics
        assert result.metrics["recall_at_50"] >= 0

    def test_holdout_genes_found(
        self, toy_temporal_data: tuple[np.ndarray, np.ndarray, np.ndarray]
    ) -> None:
        """Some holdout genes should appear in top predictions."""
        X, train_labels, holdout_labels = toy_temporal_data
        tv = TemporalValidator(n_estimators=50, seed=42)
        result = tv.run(X, train_labels, holdout_labels)
        # At least some holdout genes should be recovered in top 50
        assert result.n_holdout_in_top_k > 0

    def test_result_has_ranked_predictions(
        self, toy_temporal_data: tuple[np.ndarray, np.ndarray, np.ndarray]
    ) -> None:
        X, train_labels, holdout_labels = toy_temporal_data
        tv = TemporalValidator(n_estimators=10)
        result = tv.run(X, train_labels, holdout_labels)
        assert len(result.ranked_scores) > 0
