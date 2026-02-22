"""Tests for bias-stratified evaluation (WP8.4). Written FIRST per TDD."""

import pytest

from ssame.evaluation.bias_stratified import BiasStratifiedEvaluator, StratifiedResult


@pytest.fixture
def toy_stratified_data():
    """Gene scores with study-bias strata (publication count)."""
    gene_scores = {
        # Well-studied genes
        "TP53": 0.95, "KRAS": 0.90, "PIK3CA": 0.85, "BRAF": 0.80,
        # Moderately studied
        "CDH1": 0.70, "MDM2": 0.65, "RB1": 0.60, "PTEN": 0.55,
        # Under-studied
        "NOVEL1": 0.50, "NOVEL2": 0.45, "NOVEL3": 0.40, "NOVEL4": 0.35,
    }
    publication_counts = {
        "TP53": 100000, "KRAS": 50000, "PIK3CA": 30000, "BRAF": 25000,
        "CDH1": 5000, "MDM2": 8000, "RB1": 6000, "PTEN": 7000,
        "NOVEL1": 50, "NOVEL2": 30, "NOVEL3": 20, "NOVEL4": 10,
    }
    labels = {
        "TP53": 1, "KRAS": 1, "PIK3CA": 1, "BRAF": 1,
        "PTEN": 1, "RB1": 1,
    }
    return gene_scores, publication_counts, labels


class TestBiasStratifiedEvaluator:
    def test_creates_evaluator(self, toy_stratified_data) -> None:
        scores, pub_counts, labels = toy_stratified_data
        evaluator = BiasStratifiedEvaluator(
            gene_scores=scores,
            study_bias_proxy=pub_counts,
            labels=labels,
        )
        assert evaluator is not None

    def test_returns_result(self, toy_stratified_data) -> None:
        scores, pub_counts, labels = toy_stratified_data
        evaluator = BiasStratifiedEvaluator(
            gene_scores=scores,
            study_bias_proxy=pub_counts,
            labels=labels,
        )
        result = evaluator.evaluate()
        assert isinstance(result, StratifiedResult)

    def test_has_per_stratum_metrics(self, toy_stratified_data) -> None:
        scores, pub_counts, labels = toy_stratified_data
        evaluator = BiasStratifiedEvaluator(
            gene_scores=scores,
            study_bias_proxy=pub_counts,
            labels=labels,
            n_strata=3,
        )
        result = evaluator.evaluate()
        assert len(result.strata_metrics) == 3

    def test_strata_are_ordered(self, toy_stratified_data) -> None:
        """Strata should be ordered by study bias (low to high)."""
        scores, pub_counts, labels = toy_stratified_data
        evaluator = BiasStratifiedEvaluator(
            gene_scores=scores,
            study_bias_proxy=pub_counts,
            labels=labels,
            n_strata=3,
        )
        result = evaluator.evaluate()
        mean_biases = [s["mean_study_bias"] for s in result.strata_metrics]
        assert mean_biases == sorted(mean_biases)

    def test_bias_correlation(self, toy_stratified_data) -> None:
        """Should report correlation between scores and study bias."""
        scores, pub_counts, labels = toy_stratified_data
        evaluator = BiasStratifiedEvaluator(
            gene_scores=scores,
            study_bias_proxy=pub_counts,
            labels=labels,
        )
        result = evaluator.evaluate()
        assert "score_bias_correlation" in result.summary
