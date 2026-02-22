"""Tests for PU-aware cross-validation (WP8.1). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.evaluation.pu_cross_validation import (
    CVResult,
    PUCrossValidator,
    precision_at_k,
    recall_at_k,
)


@pytest.fixture
def toy_pu_data() -> tuple[np.ndarray, np.ndarray]:
    """Features and PU labels (200 genes, 20 known positives)."""
    rng = np.random.default_rng(42)
    n_genes = 200
    n_features = 6
    X = rng.standard_normal((n_genes, n_features))
    labels = np.zeros(n_genes)
    labels[:20] = 1
    # Add signal so that positives are distinguishable
    X[:20, 0] += 2.0
    X[:20, 1] += 1.5
    return X, labels


class TestPUCrossValidator:
    def test_creates_validator(self) -> None:
        cv = PUCrossValidator(n_folds=5, seed=42)
        assert cv.n_folds == 5

    def test_run_returns_result(
        self, toy_pu_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_pu_data
        cv = PUCrossValidator(n_folds=3, seed=42)
        result = cv.run(X, labels, n_estimators=20)
        assert isinstance(result, CVResult)

    def test_result_has_metrics(
        self, toy_pu_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_pu_data
        cv = PUCrossValidator(n_folds=3, seed=42)
        result = cv.run(X, labels, n_estimators=20)
        assert "auprc" in result.metrics
        assert "auroc" in result.metrics
        assert result.metrics["auprc"] > 0
        assert result.metrics["auroc"] > 0.5  # better than random

    def test_stratifies_on_positives(
        self, toy_pu_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """PU-aware: folds should stratify on positive labels only."""
        _X, labels = toy_pu_data
        cv = PUCrossValidator(n_folds=5, seed=42)
        folds = cv._make_folds(labels)
        for _train_idx, test_idx in folds:
            # Each test fold should have some positives
            test_pos = labels[test_idx].sum()
            assert test_pos > 0

    def test_per_fold_scores(
        self, toy_pu_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_pu_data
        cv = PUCrossValidator(n_folds=3, seed=42)
        result = cv.run(X, labels, n_estimators=20)
        assert len(result.per_fold_auprc) == 3


class TestMetricFunctions:
    def test_precision_at_k(self) -> None:
        # 5 genes ranked: first 3 are true positives
        scores = np.array([0.9, 0.8, 0.7, 0.3, 0.1])
        labels = np.array([1, 1, 1, 0, 0])
        assert precision_at_k(scores, labels, k=3) == pytest.approx(1.0)
        assert precision_at_k(scores, labels, k=5) == pytest.approx(3 / 5)

    def test_recall_at_k(self) -> None:
        scores = np.array([0.9, 0.8, 0.7, 0.3, 0.1])
        labels = np.array([1, 1, 1, 0, 0])
        assert recall_at_k(scores, labels, k=2) == pytest.approx(2 / 3)
        assert recall_at_k(scores, labels, k=3) == pytest.approx(1.0)
