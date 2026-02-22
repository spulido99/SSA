"""Structural Disruption Scorer (WP6.4).

Scores subnetworks based on how many mutations fall at protein-protein
interface residues. Mutations at binding interfaces are more likely to
disrupt function and act as driver events.

Uses PIONEER or AlphaFold PPI interface residue annotations.
"""

from __future__ import annotations

from ssame.scorer.protocol import SubnetworkScorer
from ssame.types import ScorerExplanation, Subnetwork


class StructuralDisruptionScorer(SubnetworkScorer):
    """Score modules by interface mutation enrichment.

    For each edge in the module, compute the fraction of mutations
    that fall at known interface residue positions, then aggregate
    across all edges.

    Args:
        interface_residues: Per-edge interface positions.
            Keys are (gene_a, gene_b) tuples, values are lists of
            residue positions on gene_a that form the interface.
        mutation_positions: Per-gene mutation positions.
            Keys are gene symbols, values are lists of mutated residue
            positions (e.g., from TCGA MAF data).
    """

    def __init__(
        self,
        interface_residues: dict[tuple[str, str], list[int]],
        mutation_positions: dict[str, list[int]],
    ) -> None:
        # Build lookup: gene → {(partner, interface_positions)}
        self._interfaces: dict[str, list[tuple[str, set[int]]]] = {}
        for (a, b), positions in interface_residues.items():
            self._interfaces.setdefault(a, []).append((b, set(positions)))
        self._mutations = {g: set(pos) for g, pos in mutation_positions.items()}

    def score(self, module: Subnetwork) -> float:
        if not module.genes:
            return 0.0

        symbols = module.gene_symbols
        edge_scores: list[float] = []

        for gene in symbols:
            mut_positions = self._mutations.get(gene, set())
            if not mut_positions:
                continue

            partners = self._interfaces.get(gene, [])
            for partner, iface_positions in partners:
                if partner not in symbols:
                    continue
                overlap = mut_positions & iface_positions
                if iface_positions:
                    frac = len(overlap) / len(mut_positions) if mut_positions else 0.0
                    edge_scores.append(frac)

        if not edge_scores:
            return 0.0

        return sum(edge_scores) / len(edge_scores)

    def explain(self, module: Subnetwork) -> ScorerExplanation:
        symbols = module.gene_symbols
        per_gene: dict[str, float] = {}
        details: dict[str, object] = {}

        for gene in sorted(symbols):
            mut_positions = self._mutations.get(gene, set())
            partners = self._interfaces.get(gene, [])
            gene_hits = 0
            gene_total = len(mut_positions)

            for partner, iface_positions in partners:
                if partner not in symbols:
                    continue
                overlap = mut_positions & iface_positions
                gene_hits += len(overlap)
                if overlap:
                    details[f"{gene}-{partner}_interface_hits"] = sorted(overlap)

            frac = gene_hits / gene_total if gene_total > 0 else 0.0
            per_gene[gene] = frac

        total_score = self.score(module)
        gene_list = ", ".join(sorted(symbols))
        narrative = (
            f"Structural disruption score {total_score:.2f} for module [{gene_list}]. "
        )
        if per_gene:
            top = max(per_gene, key=per_gene.get)  # type: ignore[arg-type]
            narrative += f"{top} has highest interface mutation fraction ({per_gene[top]:.2f})."

        return ScorerExplanation(
            score=total_score,
            components=per_gene,
            narrative=narrative,
            evidence=details,
        )
