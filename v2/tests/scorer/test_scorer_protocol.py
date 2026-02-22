"""Tests for SubnetworkScorer protocol (WP6.1). Written FIRST per TDD."""

import pytest

from ssame.scorer.protocol import CompositeScorer, SubnetworkScorer
from ssame.types import Gene, ScorerExplanation, Subnetwork


class DummyScorer(SubnetworkScorer):
    """Simple scorer returning module size as score."""

    def score(self, module: Subnetwork) -> float:
        return float(module.size)

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        return ScorerExplanation(
            score=self.score(module),
            narrative=f"Module has {module.size} genes",
        )


class ConstantScorer(SubnetworkScorer):
    def __init__(self, value: float) -> None:
        self.value = value

    def score(self, module: Subnetwork) -> float:
        return self.value

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        return ScorerExplanation(score=self.value)


@pytest.fixture
def sample_module() -> Subnetwork:
    genes = [Gene(symbol="TP53"), Gene(symbol="MDM2"), Gene(symbol="CDKN2A")]
    return Subnetwork(genes=genes, interactions=[])


class TestScorerProtocol:
    def test_protocol_compliance(self, sample_module: Subnetwork) -> None:
        scorer = DummyScorer()
        assert isinstance(scorer, SubnetworkScorer)
        score = scorer.score(sample_module)
        assert isinstance(score, float)
        explanation = scorer.explain(sample_module)
        assert isinstance(explanation, ScorerExplanation)

    def test_scorer_abc_not_instantiable(self) -> None:
        with pytest.raises(TypeError):
            SubnetworkScorer()  # type: ignore[abstract]


class TestCompositeScorer:
    def test_weighted_sum(self, sample_module: Subnetwork) -> None:
        s1 = ConstantScorer(0.8)
        s2 = ConstantScorer(0.6)
        composite = CompositeScorer(
            scorers=[(s1, 0.5), (s2, 0.5)],
            aggregation="weighted_sum",
        )
        score = composite.score(sample_module)
        assert score == pytest.approx(0.7)

    def test_min_aggregation(self, sample_module: Subnetwork) -> None:
        s1 = ConstantScorer(0.8)
        s2 = ConstantScorer(0.3)
        composite = CompositeScorer(
            scorers=[(s1, 1.0), (s2, 1.0)],
            aggregation="min",
        )
        score = composite.score(sample_module)
        assert score == pytest.approx(0.3)

    def test_vote_aggregation(self, sample_module: Subnetwork) -> None:
        s1 = ConstantScorer(0.8)
        s2 = ConstantScorer(0.3)
        s3 = ConstantScorer(0.9)
        composite = CompositeScorer(
            scorers=[(s1, 1.0), (s2, 1.0), (s3, 1.0)],
            aggregation="vote",
        )
        score = composite.score(sample_module)
        # Median of [0.3, 0.8, 0.9] = 0.8; 1 score strictly > 0.8
        assert score == 1

    def test_explain_combines_components(self, sample_module: Subnetwork) -> None:
        s1 = ConstantScorer(0.8)
        s2 = ConstantScorer(0.6)
        composite = CompositeScorer(
            scorers=[(s1, 0.5), (s2, 0.5)],
            aggregation="weighted_sum",
        )
        explanation = composite.explain(sample_module)
        assert len(explanation.components) == 2
