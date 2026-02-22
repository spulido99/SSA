"""Tests for MetabolicScorer (WP6.7). Written FIRST per TDD.

Scores subnetworks based on metabolic pathway disruption using
Reactome pathway membership data.
"""

import pytest

from ssame.scorer.metabolic_scorer import MetabolicScorer
from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import Gene, Interaction, ScorerExplanation, Subnetwork


@pytest.fixture
def pathway_data() -> dict[str, set[str]]:
    """Pathway name → set of member gene symbols."""
    return {
        "PI3K-AKT signaling": {"PIK3CA", "PTEN", "AKT1", "MTOR", "TSC1"},
        "RAS-MAPK signaling": {"KRAS", "BRAF", "MEK1", "ERK1", "NF1"},
        "p53 pathway": {"TP53", "MDM2", "CDKN2A", "CDKN1A"},
        "Cell cycle": {"RB1", "CDKN2A", "CDK4", "CCND1"},
    }


@pytest.fixture
def scorer(pathway_data) -> MetabolicScorer:
    return MetabolicScorer(pathways=pathway_data)


@pytest.fixture
def pi3k_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene("PIK3CA"), Gene("PTEN"), Gene("AKT1")],
        interactions=[
            Interaction(Gene("PIK3CA"), Gene("PTEN"), source="test"),
            Interaction(Gene("PIK3CA"), Gene("AKT1"), source="test"),
        ],
    )


@pytest.fixture
def cross_pathway_module() -> Subnetwork:
    """Module spanning two pathways."""
    return Subnetwork(
        genes=[Gene("TP53"), Gene("CDKN2A"), Gene("RB1")],
        interactions=[
            Interaction(Gene("TP53"), Gene("CDKN2A"), source="test"),
            Interaction(Gene("CDKN2A"), Gene("RB1"), source="test"),
        ],
    )


class TestMetabolicScorer:
    def test_implements_protocol(self, scorer) -> None:
        assert isinstance(scorer, SubnetworkScorer)

    def test_scores_module(self, scorer, pi3k_module) -> None:
        score = scorer.score(pi3k_module)
        assert isinstance(score, float)
        assert score > 0.0

    def test_high_pathway_coverage_high_score(self, scorer, pi3k_module) -> None:
        """3 out of 5 genes in PI3K-AKT → high coverage score."""
        score = scorer.score(pi3k_module)
        assert score > 0.4  # 3/5 = 0.6 coverage

    def test_cross_pathway_bonus(self, scorer, cross_pathway_module) -> None:
        """Module spanning multiple pathways should score well."""
        score = scorer.score(cross_pathway_module)
        assert score > 0.0

    def test_no_pathway_genes_zero_score(self, scorer) -> None:
        module = Subnetwork(
            genes=[Gene("NOVEL1"), Gene("NOVEL2")],
            interactions=[Interaction(Gene("NOVEL1"), Gene("NOVEL2"), source="test")],
        )
        assert scorer.score(module) == pytest.approx(0.0)

    def test_explain_returns_explanation(self, scorer, pi3k_module) -> None:
        explanation = scorer.explain(pi3k_module)
        assert isinstance(explanation, ScorerExplanation)
        assert explanation.narrative != ""

    def test_explain_lists_pathways(self, scorer, pi3k_module) -> None:
        explanation = scorer.explain(pi3k_module)
        assert "PI3K-AKT signaling" in explanation.narrative

    def test_empty_module(self, scorer) -> None:
        module = Subnetwork(genes=[], interactions=[])
        assert scorer.score(module) == 0.0
