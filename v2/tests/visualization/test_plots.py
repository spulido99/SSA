"""Tests for visualization module (WP9). Written FIRST per TDD."""

import numpy as np

from ssame.visualization.plots import (
    prepare_bias_analysis_data,
    prepare_precision_recall_data,
    prepare_rank_comparison_data,
)


class TestPrepareRankComparisonData:
    def test_returns_dict(self) -> None:
        scores = {"TP53": 0.9, "PIK3CA": 0.7, "BRAF": 0.3}
        positives = {"TP53"}
        result = prepare_rank_comparison_data(scores, positives)
        assert isinstance(result, dict)
        assert "ranks" in result
        assert "is_positive" in result
        assert "genes" in result

    def test_positives_flagged(self) -> None:
        scores = {"TP53": 0.9, "PIK3CA": 0.7, "BRAF": 0.3}
        positives = {"TP53", "PIK3CA"}
        result = prepare_rank_comparison_data(scores, positives)
        # TP53 is rank 0 (highest score), PIK3CA rank 1
        assert result["is_positive"][0] is True
        assert result["is_positive"][1] is True
        assert result["is_positive"][2] is False

    def test_ranks_are_sorted(self) -> None:
        scores = {"A": 0.1, "B": 0.9, "C": 0.5}
        result = prepare_rank_comparison_data(scores, set())
        assert result["genes"] == ["B", "C", "A"]


class TestPreparePrecisionRecallData:
    def test_returns_curves(self) -> None:
        scores = np.array([0.9, 0.8, 0.7, 0.3, 0.1])
        labels = np.array([1, 1, 0, 0, 0])
        result = prepare_precision_recall_data(scores, labels)
        assert "precision" in result
        assert "recall" in result
        assert len(result["precision"]) == len(result["recall"])

    def test_perfect_classifier(self) -> None:
        scores = np.array([1.0, 0.9, 0.0, 0.0])
        labels = np.array([1, 1, 0, 0])
        result = prepare_precision_recall_data(scores, labels)
        # Perfect separation should have high AUPRC
        assert result["auprc"] > 0.9


class TestPrepareBiasAnalysisData:
    def test_returns_study_bias_metrics(self) -> None:
        gene_scores = {"TP53": 0.9, "UNKNOWN1": 0.5, "UNKNOWN2": 0.3}
        string_degree = {"TP53": 500, "UNKNOWN1": 5, "UNKNOWN2": 2}
        alphafold_degree = {"TP53": 10, "UNKNOWN1": 8, "UNKNOWN2": 3}
        result = prepare_bias_analysis_data(
            gene_scores, string_degree, alphafold_degree
        )
        assert "string_correlation" in result
        assert "alphafold_correlation" in result

    def test_string_more_biased(self) -> None:
        """STRING degree should correlate more with scores than AlphaFold."""
        # Setup: scores correlate perfectly with STRING degree
        gene_scores = {f"G{i}": float(i) / 10 for i in range(20)}
        string_degree = {f"G{i}": i * 50 for i in range(20)}
        # AlphaFold degree is random
        rng = np.random.default_rng(42)
        alphafold_degree = {f"G{i}": int(rng.integers(1, 20)) for i in range(20)}
        result = prepare_bias_analysis_data(
            gene_scores, string_degree, alphafold_degree
        )
        assert abs(result["string_correlation"]) > abs(result["alphafold_correlation"])
