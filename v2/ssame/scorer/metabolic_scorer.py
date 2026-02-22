"""Metabolic/Pathway Disruption Scorer (WP6.7).

Scores subnetworks based on how well they cover known biological pathways.
A module that disrupts a significant portion of a single pathway (or bridges
multiple pathways) is more likely to contain driver genes.

Uses Reactome GMT or similar pathway membership data.
"""

from __future__ import annotations

from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import ScorerExplanation, Subnetwork


class MetabolicScorer(SubnetworkScorer):
    """Score modules by pathway coverage and cross-pathway bridging.

    For each pathway, compute the fraction of pathway genes that appear
    in the module. The score is the maximum single-pathway coverage,
    with a bonus for modules that bridge multiple pathways.

    Args:
        pathways: Dict of pathway_name → set of gene symbols.
    """

    def __init__(self, pathways: dict[str, set[str]]) -> None:
        self._pathways = pathways
        # Reverse index: gene → list of pathways
        self._gene_to_pathways: dict[str, list[str]] = {}
        for pathway, genes in pathways.items():
            for gene in genes:
                self._gene_to_pathways.setdefault(gene, []).append(pathway)

    def score(self, module: Subnetwork) -> float:
        if not module.genes:
            return 0.0

        symbols = module.gene_symbols

        # Per-pathway coverage
        coverages: dict[str, float] = {}
        for pathway, members in self._pathways.items():
            overlap = symbols & members
            if overlap:
                coverages[pathway] = len(overlap) / len(members)

        if not coverages:
            return 0.0

        max_coverage = max(coverages.values())

        # Cross-pathway bonus: count distinct pathways hit
        n_pathways_hit = len(coverages)
        cross_bonus = min(0.1 * (n_pathways_hit - 1), 0.3) if n_pathways_hit > 1 else 0.0

        return min(max_coverage + cross_bonus, 1.0)

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        symbols = module.gene_symbols
        per_pathway: dict[str, float] = {}
        evidence: dict[str, object] = {}

        for pathway, members in self._pathways.items():
            overlap = symbols & members
            if overlap:
                per_pathway[pathway] = len(overlap) / len(members)
                evidence[f"{pathway}_genes"] = sorted(overlap)

        total_score = self.score(module)
        gene_list = ", ".join(sorted(symbols))

        if per_pathway:
            top_pathway = max(per_pathway, key=per_pathway.get)  # type: ignore[arg-type]
            narrative = (
                f"Metabolic score {total_score:.2f} for module [{gene_list}]. "
                f"Best coverage in {top_pathway} ({per_pathway[top_pathway]:.0%}). "
                f"Hits {len(per_pathway)} pathway(s)."
            )
        else:
            narrative = (
                f"Metabolic score 0.00 for module [{gene_list}]: "
                "no genes in known pathways."
            )

        return ScorerExplanation(
            score=total_score,
            components=per_pathway,
            narrative=narrative,
            evidence=evidence,
        )
