"""Tests for ablation runner (WP8.5). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.evaluation.ablation import AblationConfig, AblationResult, AblationRunner


@pytest.fixture
def toy_data() -> tuple[np.ndarray, np.ndarray]:
    """Features and labels for ablation testing."""
    rng = np.random.default_rng(42)
    n_genes = 100
    n_features = 6
    X = rng.standard_normal((n_genes, n_features))
    labels = np.zeros(n_genes)
    labels[:15] = 1
    # Signal
    X[:15, 0] += 2.0
    X[:15, 1] += 1.5
    return X, labels


class TestAblationConfig:
    def test_creates_config(self) -> None:
        config = AblationConfig(name="baseline", encoder_level=0, n_estimators=10)
        assert config.name == "baseline"

    def test_default_config(self) -> None:
        config = AblationConfig(name="test", encoder_level=0)
        assert config.n_folds == 5


class TestAblationRunner:
    def test_creates_runner(self, toy_data: tuple[np.ndarray, np.ndarray]) -> None:
        X, labels = toy_data
        runner = AblationRunner(X=X, labels=labels)
        assert runner is not None

    def test_run_single_config(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        runner = AblationRunner(X=X, labels=labels)
        config = AblationConfig(
            name="baseline", encoder_level=0, n_estimators=10, n_folds=3,
        )
        result = runner.run_config(config)
        assert isinstance(result, AblationResult)
        assert result.name == "baseline"
        assert result.auprc > 0
        assert result.auroc > 0.5

    def test_run_multiple_configs(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        runner = AblationRunner(X=X, labels=labels)
        configs = [
            AblationConfig(name="small", encoder_level=0, n_estimators=5, n_folds=3),
            AblationConfig(name="large", encoder_level=0, n_estimators=20, n_folds=3),
        ]
        results = runner.run_all(configs)
        assert len(results) == 2
        assert results[0].name == "small"
        assert results[1].name == "large"

    def test_comparison_table(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        runner = AblationRunner(X=X, labels=labels)
        configs = [
            AblationConfig(name="A", encoder_level=0, n_estimators=10, n_folds=3),
            AblationConfig(name="B", encoder_level=0, n_estimators=20, n_folds=3),
        ]
        results = runner.run_all(configs)
        table = AblationRunner.comparison_table(results)
        assert "A" in table
        assert "B" in table
        assert "auprc" in table.lower() or "AUPRC" in table

    def test_result_has_per_fold_metrics(
        self, toy_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        X, labels = toy_data
        runner = AblationRunner(X=X, labels=labels)
        config = AblationConfig(
            name="test", encoder_level=0, n_estimators=10, n_folds=3,
        )
        result = runner.run_config(config)
        assert len(result.per_fold_auprc) == 3
        assert len(result.per_fold_auroc) == 3
