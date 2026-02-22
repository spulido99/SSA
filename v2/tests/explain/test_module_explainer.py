"""Tests for module explainer (WP9.2). Written FIRST per TDD."""

import pytest

from ssame.explain.module_explainer import ModuleExplainer
from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import Gene, Interaction, ScorerExplanation, Subnetwork


class FakeScorer(SubnetworkScorer):
    def __init__(self, name: str, value: float, narrative: str) -> None:
        self._name = name
        self._value = value
        self._narrative = narrative

    def score(self, module: Subnetwork) -> float:
        return self._value

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        return ScorerExplanation(
            score=self._value,
            components={self._name: self._value},
            narrative=self._narrative,
        )


@pytest.fixture
def toy_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene("TP53"), Gene("MDM2"), Gene("KRAS")],
        interactions=[
            Interaction(Gene("TP53"), Gene("MDM2"), source="test"),
            Interaction(Gene("KRAS"), Gene("MDM2"), source="test"),
        ],
    )


@pytest.fixture
def scorers() -> list[tuple[SubnetworkScorer, float]]:
    return [
        (FakeScorer("ME", 0.85, "Strong mutual exclusivity in module"), 0.5),
        (FakeScorer("Functional", 0.70, "2/3 genes essential in DepMap"), 0.3),
        (FakeScorer("Structural", 0.60, "Interface mutations in TP53-MDM2"), 0.2),
    ]


class TestModuleExplainer:
    def test_creates_explainer(self, scorers) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        assert explainer is not None

    def test_explain_module(self, scorers, toy_module) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(toy_module)
        assert isinstance(result, dict)

    def test_has_per_scorer_explanations(self, scorers, toy_module) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(toy_module)
        assert "scorer_explanations" in result
        assert len(result["scorer_explanations"]) == 3

    def test_has_composite_score(self, scorers, toy_module) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(toy_module)
        assert "composite_score" in result
        assert isinstance(result["composite_score"], float)

    def test_has_narrative(self, scorers, toy_module) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(toy_module)
        assert "narrative" in result
        assert len(result["narrative"]) > 0

    def test_narrative_mentions_genes(self, scorers, toy_module) -> None:
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(toy_module)
        # Narrative should mention the module genes
        assert "TP53" in result["narrative"] or "3 genes" in result["narrative"]

    def test_empty_module(self, scorers) -> None:
        module = Subnetwork(genes=[], interactions=[])
        explainer = ModuleExplainer(scorers=scorers)
        result = explainer.explain(module)
        assert "composite_score" in result
