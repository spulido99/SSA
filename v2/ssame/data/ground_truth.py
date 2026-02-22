"""CGC ground truth loader with temporal split and PU framing."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class GroundTruthResult:
    """Result of loading ground truth data.

    PU framing: positives are known drivers, everything else is UNLABELED (not negative).
    """

    positives: set[str]
    temporal_holdout: set[str] = field(default_factory=set)
    _metadata: dict[str, dict] = field(default_factory=dict, repr=False)

    def get_metadata(self, gene: str) -> dict | None:
        return self._metadata.get(gene)


def load_cgc(
    path: Path,
    *,
    tier: int | None = None,
    before_year: int | None = None,
    role: str | None = None,
) -> GroundTruthResult:
    """Load COSMIC Cancer Gene Census and return gene sets.

    Args:
        path: Path to CGC CSV file.
        tier: If set, filter to this tier only (1 or 2).
        before_year: If set, only include genes curated before this year.
            Genes curated after this year go into temporal_holdout.
        role: If set, filter to genes with this role ("TSG", "oncogene").

    Returns:
        GroundTruthResult with positives and optional temporal holdout.
    """
    df = pd.read_csv(path)

    # Standardize column names
    df.columns = df.columns.str.strip()

    # Parse curation date
    date_col = "First Curated Date"
    if date_col in df.columns:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["curated_year"] = df[date_col].dt.year
    else:
        df["curated_year"] = pd.NA

    # Filter by tier
    if tier is not None:
        df = df[df["Tier"] == tier]

    # Filter by role
    if role is not None:
        df = df[df["Role in Cancer"].str.contains(role, case=False, na=False)]

    # Build metadata
    metadata: dict[str, dict] = {}
    for _, row in df.iterrows():
        symbol = row["Gene Symbol"]
        metadata[symbol] = {
            "entrez_id": int(row["Entrez GeneId"]) if pd.notna(row["Entrez GeneId"]) else None,
            "tier": int(row["Tier"]) if pd.notna(row["Tier"]) else None,
            "role": row.get("Role in Cancer", ""),
            "curated_year": int(row["curated_year"]) if pd.notna(row["curated_year"]) else None,
        }

    # Temporal split
    if before_year is not None:
        pre_mask = df["curated_year"] < before_year
        post_mask = df["curated_year"] >= before_year
        positives = set(df.loc[pre_mask, "Gene Symbol"])
        temporal_holdout = set(df.loc[post_mask, "Gene Symbol"])
    else:
        positives = set(df["Gene Symbol"])
        temporal_holdout = set()

    return GroundTruthResult(
        positives=positives,
        temporal_holdout=temporal_holdout,
        _metadata=metadata,
    )
