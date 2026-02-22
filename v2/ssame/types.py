"""Core domain types used across all SSA-ME v2 modules."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Gene:
    """A gene identified by its Hugo symbol."""

    symbol: str
    entrez_id: int | None = None

    def __post_init__(self) -> None:
        if not self.symbol:
            raise ValueError("Gene symbol must not be empty")

    def __str__(self) -> str:
        return self.symbol

    def __hash__(self) -> int:
        return hash(self.symbol)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Gene):
            return NotImplemented
        return self.symbol == other.symbol


@dataclass(frozen=True)
class Interaction:
    """An interaction between two genes from a specific network source."""

    gene_a: Gene
    gene_b: Gene
    source: str
    confidence: float = 1.0

    def __post_init__(self) -> None:
        if self.gene_a == self.gene_b:
            raise ValueError(f"Interaction cannot be a self-loop: {self.gene_a}")

    def as_tuple(self) -> tuple[str, str]:
        return (self.gene_a.symbol, self.gene_b.symbol)


@dataclass
class Subnetwork:
    """A candidate module: genes + their interactions."""

    genes: list[Gene]
    interactions: list[Interaction]
    mutations: dict[str, list[str]] = field(default_factory=dict)
    samples: list[str] = field(default_factory=list)

    @property
    def gene_symbols(self) -> set[str]:
        return {g.symbol for g in self.genes}

    @property
    def size(self) -> int:
        return len(self.genes)

    def contains(self, symbol: str) -> bool:
        return symbol in self.gene_symbols

    def merge(self, other: Subnetwork) -> Subnetwork:
        seen = set()
        merged_genes = []
        for g in self.genes + other.genes:
            if g.symbol not in seen:
                seen.add(g.symbol)
                merged_genes.append(g)
        merged_symbols = {g.symbol for g in merged_genes}
        merged_interactions = []
        seen_pairs: set[tuple[str, str]] = set()
        for ix in self.interactions + other.interactions:
            pair = ix.as_tuple()
            rev = (pair[1], pair[0])
            if (
                pair not in seen_pairs
                and rev not in seen_pairs
                and ix.gene_a.symbol in merged_symbols
                and ix.gene_b.symbol in merged_symbols
            ):
                seen_pairs.add(pair)
                merged_interactions.append(ix)
        return Subnetwork(genes=merged_genes, interactions=merged_interactions)

    def intersection(self, other: Subnetwork) -> Subnetwork:
        common_symbols = self.gene_symbols & other.gene_symbols
        common_genes = [g for g in self.genes if g.symbol in common_symbols]
        common_interactions = [
            ix
            for ix in self.interactions
            if ix.gene_a.symbol in common_symbols and ix.gene_b.symbol in common_symbols
        ]
        return Subnetwork(genes=common_genes, interactions=common_interactions)


@dataclass
class AlterationMatrix:
    """Binary gene x sample alteration matrix."""

    data: np.ndarray  # shape: (n_genes, n_samples), dtype: int8
    genes: list[str]
    samples: list[str]

    def __post_init__(self) -> None:
        if self.data.shape != (len(self.genes), len(self.samples)):
            raise ValueError(
                f"Data shape {self.data.shape} does not match "
                f"genes ({len(self.genes)}) x samples ({len(self.samples)})"
            )
        self._gene_index: dict[str, int] = {g: i for i, g in enumerate(self.genes)}

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def gene_frequency(self, gene: str) -> float:
        idx = self._gene_index[gene]
        return float(self.data[idx].sum()) / self.n_samples

    def coverage(self, gene_list: list[str]) -> float:
        indices = [self._gene_index[g] for g in gene_list]
        altered = self.data[indices].any(axis=0)
        return float(altered.sum()) / self.n_samples

    def co_occurrence_count(self, gene_a: str, gene_b: str) -> int:
        idx_a = self._gene_index[gene_a]
        idx_b = self._gene_index[gene_b]
        return int((self.data[idx_a] & self.data[idx_b]).sum())

    def filter_by_min_frequency(self, min_freq: float) -> AlterationMatrix:
        freqs = self.data.sum(axis=1) / self.n_samples
        mask = freqs >= min_freq
        return AlterationMatrix(
            data=self.data[mask],
            genes=[g for g, m in zip(self.genes, mask, strict=True) if m],
            samples=self.samples,
        )

    def submatrix(self, gene_subset: list[str]) -> AlterationMatrix:
        indices = [self._gene_index[g] for g in gene_subset]
        return AlterationMatrix(
            data=self.data[indices],
            genes=gene_subset,
            samples=self.samples,
        )


@dataclass
class ScorerExplanation:
    """Structured explanation of WHY a module scored high/low."""

    score: float
    components: dict[str, float] = field(default_factory=dict)
    narrative: str = ""
    evidence: dict[str, object] = field(default_factory=dict)
