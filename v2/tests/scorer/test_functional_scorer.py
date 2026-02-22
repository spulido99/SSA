"""Tests for FunctionalGenomicsScorer (WP6.3). Written FIRST per TDD."""

import pytest

from ssame.scorer.functional_scorer import FunctionalGenomicsScorer
from ssame.types import Gene, Interaction, ScorerExplanation, Subnetwork


@pytest.fixture
def essentiality_scores() -> dict[str, float]:
    """Per-gene CRISPR essentiality scores (more negative = more essential)."""
    return {
        "TP53": -1.5,   # Highly essential
        "PIK3CA": -1.2,  # Essential
        "CDH1": -0.8,   # Moderately essential
        "KRAS": -1.8,   # Very essential
        "BRAF": -0.3,   # Not very essential
        "MDM2": -0.9,
        "CDKN2A": -0.5,
        "RANDOM1": 0.1,  # Non-essential
        "RANDOM2": 0.2,
    }


@pytest.fixture
def driver_module() -> Subnetwork:
    genes = [Gene(symbol="TP53"), Gene(symbol="PIK3CA"), Gene(symbol="KRAS")]
    interactions = [
        Interaction(gene_a=Gene(symbol="TP53"), gene_b=Gene(symbol="PIK3CA"), source="ppi"),
    ]
    return Subnetwork(genes=genes, interactions=interactions)


@pytest.fixture
def random_module() -> Subnetwork:
    genes = [Gene(symbol="RANDOM1"), Gene(symbol="RANDOM2"), Gene(symbol="BRAF")]
    return Subnetwork(genes=genes, interactions=[])


class TestFunctionalGenomicsScorer:
    def test_creates_scorer(self, essentiality_scores: dict[str, float]) -> None:
        scorer = FunctionalGenomicsScorer(essentiality_scores=essentiality_scores)
        assert scorer is not None

    def test_driver_module_scores_higher(
        self,
        essentiality_scores: dict[str, float],
        driver_module: Subnetwork,
        random_module: Subnetwork,
    ) -> None:
        """Module with essential genes should score higher."""
        scorer = FunctionalGenomicsScorer(essentiality_scores=essentiality_scores)
        driver_score = scorer.score(driver_module)
        random_score = scorer.score(random_module)
        assert driver_score > random_score

    def test_score_is_non_negative(
        self,
        essentiality_scores: dict[str, float],
        driver_module: Subnetwork,
    ) -> None:
        scorer = FunctionalGenomicsScorer(essentiality_scores=essentiality_scores)
        assert scorer.score(driver_module) >= 0

    def test_explain_returns_explanation(
        self,
        essentiality_scores: dict[str, float],
        driver_module: Subnetwork,
    ) -> None:
        scorer = FunctionalGenomicsScorer(essentiality_scores=essentiality_scores)
        explanation = scorer.explain(driver_module)
        assert isinstance(explanation, ScorerExplanation)
        assert explanation.score > 0
        assert "essentiality" in explanation.narrative.lower() or "essential" in explanation.narrative.lower()

    def test_missing_genes_handled(self) -> None:
        """Genes not in essentiality data should be handled gracefully."""
        scorer = FunctionalGenomicsScorer(essentiality_scores={"TP53": -1.5})
        module = Subnetwork(
            genes=[Gene(symbol="TP53"), Gene(symbol="UNKNOWN")],
            interactions=[],
        )
        score = scorer.score(module)
        assert score >= 0
