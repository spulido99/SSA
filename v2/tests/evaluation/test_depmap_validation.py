"""Tests for DepMap functional validation (WP8.3). Written FIRST per TDD."""

import pytest

from ssame.evaluation.depmap_validation import DepMapValidationResult, DepMapValidator


@pytest.fixture
def toy_validation_data():
    gene_scores = {
        "TP53": 0.9, "PIK3CA": 0.8, "KRAS": 0.7,
        "CDH1": 0.5, "BRAF": 0.4, "MDM2": 0.3,
        "RANDOM1": 0.1, "RANDOM2": 0.05,
    }
    essentiality = {
        "TP53": -1.5, "PIK3CA": -1.2, "KRAS": -1.8,
        "CDH1": -0.8, "BRAF": -0.3, "MDM2": -0.9,
        "RANDOM1": 0.1, "RANDOM2": 0.2,
    }
    return gene_scores, essentiality


class TestDepMapValidator:
    def test_creates_validator(self, toy_validation_data) -> None:
        gene_scores, essentiality = toy_validation_data
        validator = DepMapValidator(
            gene_scores=gene_scores,
            essentiality_scores=essentiality,
        )
        assert validator is not None

    def test_returns_result(self, toy_validation_data) -> None:
        gene_scores, essentiality = toy_validation_data
        validator = DepMapValidator(
            gene_scores=gene_scores,
            essentiality_scores=essentiality,
        )
        result = validator.validate()
        assert isinstance(result, DepMapValidationResult)

    def test_correlation_exists(self, toy_validation_data) -> None:
        """Predicted scores should correlate with essentiality."""
        gene_scores, essentiality = toy_validation_data
        validator = DepMapValidator(
            gene_scores=gene_scores,
            essentiality_scores=essentiality,
        )
        result = validator.validate()
        # Higher predicted scores should correspond to more negative essentiality
        # (i.e., more essential = more likely driver)
        assert result.correlation < 0  # Negative correlation expected

    def test_enrichment_at_top_k(self, toy_validation_data) -> None:
        """Top predicted genes should be enriched for essential genes."""
        gene_scores, essentiality = toy_validation_data
        validator = DepMapValidator(
            gene_scores=gene_scores,
            essentiality_scores=essentiality,
        )
        result = validator.validate()
        assert result.enrichment_top_k > 0
        # Mean essentiality of top-k should be more negative than bottom-k
        assert result.mean_essentiality_top < result.mean_essentiality_bottom

    def test_handles_missing_genes(self) -> None:
        """Should handle genes missing from essentiality data."""
        gene_scores = {"TP53": 0.9, "UNKNOWN": 0.5}
        essentiality = {"TP53": -1.5}
        validator = DepMapValidator(
            gene_scores=gene_scores,
            essentiality_scores=essentiality,
        )
        result = validator.validate()
        assert result is not None
