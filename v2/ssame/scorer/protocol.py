"""SubnetworkScorer protocol (ABC) and CompositeScorer."""

from __future__ import annotations

from abc import ABC, abstractmethod
from statistics import median

from ssame.types import ScorerExplanation, Subnetwork


class SubnetworkScorer(ABC):
    """Abstract base class for subnetwork scoring.

    The v1 architectural power, preserved. ANY scoring logic goes here.
    Does NOT need to be differentiable.
    """

    @abstractmethod
    def score(self, module: Subnetwork) -> float:
        """Score a candidate subnetwork. Any logic. Any complexity."""
        ...

    @abstractmethod
    def explain(self, module: Subnetwork) -> ScorerExplanation:
        """Explain why this module scored the way it did."""
        ...


class CompositeScorer(SubnetworkScorer):
    """Chain multiple scorers with configurable aggregation."""

    def __init__(
        self,
        scorers: list[tuple[SubnetworkScorer, float]],
        aggregation: str = "weighted_sum",
    ) -> None:
        self.scorers = scorers
        self.aggregation = aggregation

    def score(self, module: Subnetwork) -> float:
        scores = [(s.score(module), w) for s, w in self.scorers]

        if self.aggregation == "weighted_sum":
            return sum(s * w for s, w in scores)
        elif self.aggregation == "min":
            return min(s for s, _ in scores)
        elif self.aggregation == "vote":
            raw = [s for s, _ in scores]
            threshold = median(raw)
            return sum(1 for s in raw if s > threshold)
        else:
            raise ValueError(f"Unknown aggregation: {self.aggregation}")

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        components: dict[str, float] = {}
        for i, (scorer, _weight) in enumerate(self.scorers):
            name = type(scorer).__name__
            components[f"{name}_{i}"] = scorer.score(module)

        return ScorerExplanation(
            score=self.score(module),
            components=components,
            narrative=f"Composite score from {len(self.scorers)} scorers ({self.aggregation})",
        )
