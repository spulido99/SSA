"""Tests for tabular baseline encoder (WP4.4). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.encoder.tabular_baseline import TabularBaseline


@pytest.fixture
def toy_data() -> tuple[np.ndarray, np.ndarray]:
    """Create toy feature matrix and PU labels.

    Features are designed so that known drivers (labels=1) have higher values
    in feature 0 (mutation frequency proxy).
    """
    rng = np.random.default_rng(42)
    n_genes = 200
    n_features = 6

    X = rng.standard_normal((n_genes, n_features))
    labels = np.zeros(n_genes)
    # 20 known positives with higher feature values
    labels[:20] = 1
    X[:20, 0] += 2.0  # mutation frequency signal
    X[:20, 1] += 1.5  # another signal

    return X, labels


class TestTabularBaseline:
    def test_creates_model(self) -> None:
        model = TabularBaseline()
        assert model is not None

    def test_fit_predict(self, toy_data: tuple[np.ndarray, np.ndarray]) -> None:
        X, labels = toy_data
        model = TabularBaseline(n_estimators=50)
        model.fit(X, labels)
        scores = model.predict(X)
        assert isinstance(scores, np.ndarray)
        assert scores.shape == (X.shape[0],)
        assert np.all(scores >= 0) and np.all(scores <= 1)

    def test_known_drivers_score_higher(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        model = TabularBaseline(n_estimators=100)
        model.fit(X, labels)
        scores = model.predict(X)
        pos_mean = scores[labels == 1].mean()
        unl_mean = scores[labels == 0].mean()
        assert pos_mean > unl_mean

    def test_feature_importance(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        model = TabularBaseline(n_estimators=50)
        model.fit(X, labels)
        importance = model.feature_importance()
        assert isinstance(importance, np.ndarray)
        assert importance.shape == (X.shape[1],)
        # Feature 0 should be among the most important (strongest signal)
        assert importance[0] > importance.mean()

    def test_pu_bagging(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """PU bagging: train on subsamples where unlabeled are treated as negative."""
        X, labels = toy_data
        model = TabularBaseline(n_estimators=50, n_pu_bags=5)
        model.fit(X, labels)
        scores = model.predict(X)
        pos_mean = scores[labels == 1].mean()
        unl_mean = scores[labels == 0].mean()
        assert pos_mean > unl_mean

    def test_reproducible_with_seed(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        model1 = TabularBaseline(n_estimators=50, seed=42)
        model1.fit(X, labels)
        scores1 = model1.predict(X)

        model2 = TabularBaseline(n_estimators=50, seed=42)
        model2.fit(X, labels)
        scores2 = model2.predict(X)

        np.testing.assert_array_almost_equal(scores1, scores2)
