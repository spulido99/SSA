"""Open-access ground truth loaders: IntOGen, NCG, CancerMine.

Replaces CGC/OncoKB dependency with three freely available databases:
- IntOGen: ~600 driver genes from 28K tumors (intogen.org)
- NCG 7.1: 3,347 cancer drivers (network-cancer-genes.org)
- CancerMine: text-mined oncogenes/TSGs (cancermine.ai4biomed.org)

A gene is "known driver" if it appears in >= min_sources databases.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ssame.data.ground_truth import GroundTruthResult

# ── Role normalization ────────────────────────────────────────────────

_ROLE_MAP_INTOGEN = {"lof": "TSG", "act": "oncogene"}
_ROLE_MAP_CANCERMINE = {
    "tumor_suppressor": "TSG",
    "oncogene": "oncogene",
    "driver": "driver",
}


def _normalize_role(raw: str, mapping: dict[str, str]) -> str:
    return mapping.get(raw.strip().lower(), raw.strip())


# ── IntOGen ───────────────────────────────────────────────────────────

def load_intogen(
    path: Path,
    *,
    max_qvalue: float = 0.05,
) -> GroundTruthResult:
    """Load IntOGen driver gene catalog.

    Args:
        path: Path to drivers.tsv (from IntOGen-Drivers-*.zip).
        max_qvalue: Maximum q-value threshold (default 0.05).

    Returns:
        GroundTruthResult with driver genes as positives.
    """
    if not path.exists():
        raise FileNotFoundError(f"IntOGen file not found: {path}")

    df = pd.read_csv(path, sep="\t")

    # Filter to actual drivers (Unfiltered_drivers.tsv has DRIVER column,
    # Compendium_Cancer_Genes.tsv has IS_DRIVER column)
    if "DRIVER" in df.columns:
        df = df[df["DRIVER"] == True]  # noqa: E712
    elif "IS_DRIVER" in df.columns:
        df = df[df["IS_DRIVER"] == True]  # noqa: E712

    # Filter by q-value
    qval_col = "QVALUE_COMBINATION"
    if qval_col in df.columns:
        df = df[df[qval_col] <= max_qvalue]

    # Deduplicate: same gene across multiple cohorts
    metadata: dict[str, dict] = {}
    positives: set[str] = set()

    for _, row in df.iterrows():
        symbol = str(row["SYMBOL"]).strip()
        positives.add(symbol)

        if symbol not in metadata:
            role_raw = str(row.get("ROLE", "")).strip()
            metadata[symbol] = {
                "role": _normalize_role(role_raw, _ROLE_MAP_INTOGEN),
                "source": "IntOGen",
                "qvalue": float(row[qval_col]) if qval_col in df.columns else None,
            }

    return GroundTruthResult(positives=positives, _metadata=metadata)


# ── NCG ───────────────────────────────────────────────────────────────

def load_ncg(
    path: Path,
) -> GroundTruthResult:
    """Load Network of Cancer Genes (NCG) driver list.

    Args:
        path: Path to NCG cancer drivers annotation TSV.

    Returns:
        GroundTruthResult with driver genes as positives.
    """
    if not path.exists():
        raise FileNotFoundError(f"NCG file not found: {path}")

    df = pd.read_csv(path, sep="\t")

    # Auto-detect gene symbol column
    sym_col = None
    for candidate in ("ncg_symbol", "symbol", "gene_symbol", "SYMBOL"):
        if candidate in df.columns:
            sym_col = candidate
            break
    if sym_col is None:
        sym_col = df.columns[0]

    # Role normalization map (handles both original NCG format and PMC supplementary)
    _ncg_role_map = {
        "oncogene": "oncogene",
        "tumour suppressor": "TSG",
        "tumor suppressor": "TSG",
        "dual role": "both",
        "conflicting annotations": "both",
        "unclassified": "unknown",
    }

    metadata: dict[str, dict] = {}
    positives: set[str] = set()

    for _, row in df.iterrows():
        symbol = str(row[sym_col]).strip()
        positives.add(symbol)

        if symbol not in metadata:
            # Try new format: 'role' column from PMC supplementary
            if "role" in df.columns:
                raw_role = str(row.get("role", "")).strip().lower()
                role = _ncg_role_map.get(raw_role, "unknown")
            else:
                # Legacy format: NCG_oncogene / NCG_tsg boolean columns
                is_oncogene = bool(row.get("NCG_oncogene", 0))
                is_tsg = bool(row.get("NCG_tsg", 0))
                if is_oncogene and is_tsg:
                    role = "both"
                elif is_oncogene:
                    role = "oncogene"
                elif is_tsg:
                    role = "TSG"
                else:
                    role = "unknown"

            entrez = row.get("entrez")
            metadata[symbol] = {
                "role": role,
                "source": "NCG",
                "entrez_id": int(entrez) if pd.notna(entrez) else None,
            }

    return GroundTruthResult(positives=positives, _metadata=metadata)


# ── CancerMine ────────────────────────────────────────────────────────

def load_cancermine(
    path: Path,
    *,
    min_citations: int = 5,
) -> GroundTruthResult:
    """Load CancerMine text-mined cancer gene annotations.

    Args:
        path: Path to cancermine_collated.tsv from Zenodo.
        min_citations: Minimum total citation count to include a gene.

    Returns:
        GroundTruthResult with driver genes as positives.
    """
    if not path.exists():
        raise FileNotFoundError(f"CancerMine file not found: {path}")

    df = pd.read_csv(path, sep="\t")

    # Auto-detect gene column (prefer symbol over HGNC ID)
    gene_col = None
    for candidate in ("gene_normalized", "hugo_symbol", "gene_hugo_id"):
        if candidate in df.columns:
            gene_col = candidate
            break
    if gene_col is None:
        gene_col = df.columns[1]

    # Aggregate citations per gene
    cite_col = "citation_count" if "citation_count" in df.columns else None

    gene_citations: dict[str, int] = {}
    gene_roles: dict[str, list[str]] = {}

    for _, row in df.iterrows():
        symbol = str(row[gene_col]).strip()
        citations = int(row[cite_col]) if cite_col else 1
        gene_citations[symbol] = gene_citations.get(symbol, 0) + citations

        role_raw = str(row.get("role", "")).strip()
        gene_roles.setdefault(symbol, []).append(role_raw)

    # Filter by citation threshold
    metadata: dict[str, dict] = {}
    positives: set[str] = set()

    for symbol, total_cites in gene_citations.items():
        if total_cites < min_citations:
            continue

        positives.add(symbol)

        # Determine dominant role
        role_counts = Counter(
            _normalize_role(r, _ROLE_MAP_CANCERMINE) for r in gene_roles[symbol]
        )
        dominant_role = role_counts.most_common(1)[0][0]

        metadata[symbol] = {
            "role": dominant_role,
            "source": "CancerMine",
            "total_citations": total_cites,
            "roles": dict(role_counts),
        }

    return GroundTruthResult(positives=positives, _metadata=metadata)


# ── Consensus ─────────────────────────────────────────────────────────

@dataclass
class ConsensusGroundTruth:
    """Combined ground truth from multiple open databases.

    Attributes:
        positives: Genes in >= min_sources databases.
        high_confidence: Genes in all provided databases.
        source_count: Gene → number of sources it appears in.
        consensus_role: Gene → majority-vote role (oncogene/TSG).
        per_source: Dict of source name → GroundTruthResult.
    """

    positives: set[str]
    high_confidence: set[str]
    source_count: dict[str, int]
    consensus_role: dict[str, str]
    per_source: dict[str, GroundTruthResult] = field(repr=False)

    def as_labels(self) -> dict[str, int]:
        """Return dict[str, int] suitable for pipeline consumption."""
        return {gene: 1 for gene in self.positives}


def build_consensus_ground_truth(
    *,
    intogen_path: Path | None = None,
    ncg_path: Path | None = None,
    cancermine_path: Path | None = None,
    min_sources: int = 2,
    intogen_max_qvalue: float = 0.05,
    cancermine_min_citations: int = 5,
) -> ConsensusGroundTruth:
    """Build consensus ground truth from open-access databases.

    Args:
        intogen_path: Path to IntOGen drivers.tsv.
        ncg_path: Path to NCG annotation TSV.
        cancermine_path: Path to cancermine_collated.tsv.
        min_sources: Minimum number of databases a gene must appear in.
        intogen_max_qvalue: Q-value threshold for IntOGen.
        cancermine_min_citations: Citation threshold for CancerMine.

    Returns:
        ConsensusGroundTruth with gene sets at different confidence levels.
    """
    sources: dict[str, GroundTruthResult] = {}

    if intogen_path is not None:
        sources["IntOGen"] = load_intogen(intogen_path, max_qvalue=intogen_max_qvalue)
    if ncg_path is not None:
        sources["NCG"] = load_ncg(ncg_path)
    if cancermine_path is not None:
        sources["CancerMine"] = load_cancermine(
            cancermine_path, min_citations=cancermine_min_citations
        )

    # Count sources per gene
    source_count: dict[str, int] = {}
    for result in sources.values():
        for gene in result.positives:
            source_count[gene] = source_count.get(gene, 0) + 1

    n_sources = len(sources)
    positives = {g for g, c in source_count.items() if c >= min_sources}
    high_confidence = {g for g, c in source_count.items() if c >= n_sources}

    # Consensus role via majority vote
    consensus_role: dict[str, str] = {}
    for gene in source_count:
        roles: list[str] = []
        for result in sources.values():
            meta = result.get_metadata(gene)
            if meta and meta.get("role") and meta["role"] not in ("unknown", "driver"):
                roles.append(meta["role"])
        if roles:
            role_counts = Counter(roles)
            consensus_role[gene] = role_counts.most_common(1)[0][0]
        else:
            consensus_role[gene] = "unknown"

    return ConsensusGroundTruth(
        positives=positives,
        high_confidence=high_confidence,
        source_count=source_count,
        consensus_role=consensus_role,
        per_source=sources,
    )
