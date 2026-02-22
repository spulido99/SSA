"""Module Explainer — human-readable explanations for scored subnetworks (WP9.2).

Combines explanations from multiple scorers into a unified narrative
explaining why a candidate driver module scored high or low.
"""

from __future__ import annotations

from ssame.scorer.protocol import CompositeScorer, SubnetworkScorer
from ssame.types import Subnetwork


class ModuleExplainer:
    """Generate unified explanations for scored modules.

    Takes the same (scorer, weight) pairs used by CompositeScorer
    and produces a structured explanation dict with per-scorer
    breakdowns and a combined narrative.

    Args:
        scorers: List of (scorer, weight) tuples.
        aggregation: Aggregation method (passed to CompositeScorer).
    """

    def __init__(
        self,
        scorers: list[tuple[SubnetworkScorer, float]],
        aggregation: str = "weighted_sum",
    ) -> None:
        self._scorers = scorers
        self._composite = CompositeScorer(scorers, aggregation=aggregation)

    def explain(self, module: Subnetwork) -> dict:
        """Generate a full explanation for a module.

        Returns:
            Dict with keys:
            - composite_score: float
            - scorer_explanations: list of per-scorer explanation dicts
            - narrative: str combining all explanations
        """
        composite_score = self._composite.score(module)

        scorer_explanations = []
        narratives = []

        for scorer, weight in self._scorers:
            explanation = scorer.explain(module)
            scorer_explanations.append({
                "scorer": type(scorer).__name__,
                "weight": weight,
                "score": explanation.score,
                "components": explanation.components,
                "narrative": explanation.narrative,
            })
            if explanation.narrative:
                narratives.append(explanation.narrative)

        gene_symbols = sorted(module.gene_symbols) if module.genes else []
        gene_str = ", ".join(gene_symbols) if gene_symbols else "empty"
        n_genes = len(gene_symbols)

        narrative = (
            f"Module [{gene_str}] ({n_genes} genes) — "
            f"composite score: {composite_score:.3f}. "
        )
        if narratives:
            narrative += " ".join(narratives)

        return {
            "composite_score": composite_score,
            "scorer_explanations": scorer_explanations,
            "narrative": narrative,
            "genes": gene_symbols,
        }
