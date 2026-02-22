"""Tests for CompositeScorer (WP6.5). Written FIRST per TDD."""

import pytest

from ssame.scorer.protocol import CompositeScorer, SubnetworkScorer
from ssame.types import Gene, Interaction, ScorerExplanation, Subnetwork


class ConstantScorer(SubnetworkScorer):
    """A scorer that always returns the same value, for testing."""

    def __init__(self, value: float) -> None:
        self.value = value

    def score(self, module: Subnetwork) -> float:
        return self.value

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        return ScorerExplanation(
            score=self.value,
            narrative=f"Constant: {self.value}",
        )


@pytest.fixture
def toy_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene("TP53"), Gene("MDM2")],
        interactions=[Interaction(Gene("TP53"), Gene("MDM2"), source="test")],
    )


class TestCompositeScorer:
    def test_weighted_sum(self, toy_module) -> None:
        scorer = CompositeScorer(
            scorers=[
                (ConstantScorer(0.8), 0.5),
                (ConstantScorer(0.6), 0.5),
            ],
            aggregation="weighted_sum",
        )
        assert scorer.score(toy_module) == pytest.approx(0.7)

    def test_min_aggregation(self, toy_module) -> None:
        scorer = CompositeScorer(
            scorers=[
                (ConstantScorer(0.8), 1.0),
                (ConstantScorer(0.2), 1.0),
            ],
            aggregation="min",
        )
        assert scorer.score(toy_module) == pytest.approx(0.2)

    def test_vote_aggregation(self, toy_module) -> None:
        scorer = CompositeScorer(
            scorers=[
                (ConstantScorer(0.9), 1.0),
                (ConstantScorer(0.1), 1.0),
                (ConstantScorer(0.8), 1.0),
            ],
            aggregation="vote",
        )
        # median is 0.8, votes above: 0.9 → 1 vote
        score = scorer.score(toy_module)
        assert score == 1

    def test_explain_lists_components(self, toy_module) -> None:
        scorer = CompositeScorer(
            scorers=[
                (ConstantScorer(0.8), 0.5),
                (ConstantScorer(0.6), 0.5),
            ],
        )
        explanation = scorer.explain(toy_module)
        assert isinstance(explanation, ScorerExplanation)
        assert len(explanation.components) == 2

    def test_implements_protocol(self, toy_module) -> None:
        scorer = CompositeScorer(scorers=[(ConstantScorer(0.5), 1.0)])
        assert isinstance(scorer, SubnetworkScorer)

    def test_unknown_aggregation_raises(self, toy_module) -> None:
        scorer = CompositeScorer(
            scorers=[(ConstantScorer(0.5), 1.0)],
            aggregation="invalid",
        )
        with pytest.raises(ValueError, match="Unknown aggregation"):
            scorer.score(toy_module)
