#!/usr/bin/env python
"""DepMap functional validation: assess predicted driver scores against CRISPR essentiality.

Uses existing DepMapValidator and load_depmap_essentiality to:
1. Load ranked gene scores from a prior pipeline run (or re-run pipeline)
2. Load DepMap CRISPRGeneEffect data
3. Compute Spearman correlation + enrichment at multiple k values
4. Stratify by known_drivers vs novel_predictions
5. Save results to v2/results/

Usage:
    python scripts/run_depmap_validation.py
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.depmap_loader import load_depmap_essentiality
from ssame.evaluation.depmap_validation import DepMapValidator

RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "raw"


def load_ranked_genes(path: Path) -> dict[str, float]:
    """Load ranked_genes.tsv -> dict[gene, score]."""
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["gene"], df["score"].astype(float)))


def load_known_drivers(path: Path) -> set[str]:
    """Load ground truth genes from ground_truth_genes.tsv."""
    df = pd.read_csv(path, sep="\t")
    return set(df[df["is_consensus"] == True]["gene"])  # noqa: E712


def enrichment_at_k(
    gene_scores: dict[str, float],
    essentiality: dict[str, float],
    k: int,
    threshold: float = -0.5,
) -> dict[str, float]:
    """Compute fraction of CRISPR-essential genes in top-k predicted."""
    common = sorted(set(gene_scores) & set(essentiality))
    if not common:
        return {"k": k, "n_essential_in_top_k": 0, "fraction_essential": 0.0}

    ranked = sorted(common, key=lambda g: -gene_scores[g])
    top_k = ranked[:k]
    n_essential = sum(1 for g in top_k if essentiality[g] < threshold)

    return {
        "k": k,
        "n_essential_in_top_k": n_essential,
        "fraction_essential": n_essential / k if k > 0 else 0.0,
        "n_genes_evaluated": len(common),
    }


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 60)
    print("DEPMAP FUNCTIONAL VALIDATION")
    print("=" * 60)

    # ── Load gene scores ──────────────────────────────────────────
    ranked_path = RESULTS_DIR / "ranked_genes.tsv"
    if not ranked_path.exists():
        print(f"SKIP: ranked_genes.tsv not found at {ranked_path}")
        print("Run scripts/save_results.py first.")
        return

    gene_scores = load_ranked_genes(ranked_path)
    print(f"Loaded {len(gene_scores)} gene scores")

    # ── Load DepMap ───────────────────────────────────────────────
    depmap_path = DATA_DIR / "depmap" / "CRISPRGeneEffect.csv"
    if not depmap_path.exists():
        print(f"SKIP: DepMap data not found at {depmap_path}")
        return

    print("Loading DepMap CRISPRGeneEffect...")
    depmap = load_depmap_essentiality(depmap_path)
    print(f"  {len(depmap.mean_essentiality)} genes, {len(depmap.cell_line_ids)} cell lines")

    # ── Run validation ────────────────────────────────────────────
    print("\nRunning DepMap validation...")
    validator = DepMapValidator(
        gene_scores=gene_scores,
        essentiality_scores=depmap.mean_essentiality,
    )
    result = validator.validate()
    print(f"  Spearman correlation: {result.correlation:.4f} (p={result.p_value:.2e})")
    print(f"  Enrichment (top vs bottom 25%): {result.enrichment_top_k:.4f}")
    print(f"  Mean essentiality (top 25%): {result.mean_essentiality_top:.4f}")
    print(f"  Mean essentiality (bottom 25%): {result.mean_essentiality_bottom:.4f}")
    print(f"  Genes evaluated: {result.n_genes_evaluated}")

    # ── Enrichment at multiple k ──────────────────────────────────
    print("\nEnrichment at multiple k values:")
    enrichment_rows = []
    for k in [10, 25, 50, 100, 200]:
        enr = enrichment_at_k(gene_scores, depmap.mean_essentiality, k)
        enrichment_rows.append(enr)
        print(f"  Top-{k}: {enr['fraction_essential']:.2%} CRISPR-essential")

    df_enr = pd.DataFrame(enrichment_rows)
    df_enr.to_csv(RESULTS_DIR / "depmap_enrichment_curve.tsv", sep="\t", index=False)
    print(f"  Saved: depmap_enrichment_curve.tsv")

    # ── Stratify: known drivers vs novel predictions ──────────────
    gt_path = RESULTS_DIR / "ground_truth_genes.tsv"
    if gt_path.exists():
        known_drivers = load_known_drivers(gt_path)
        print(f"\nStratification: {len(known_drivers)} known drivers")

        common_genes = sorted(set(gene_scores) & set(depmap.mean_essentiality))
        rows = []
        for gene in common_genes:
            is_known = gene in known_drivers
            rows.append({
                "gene": gene,
                "predicted_score": gene_scores[gene],
                "essentiality": depmap.mean_essentiality[gene],
                "is_known_driver": is_known,
                "category": "known_driver" if is_known else "novel",
            })
        df_strat = pd.DataFrame(rows)
        df_strat.to_csv(RESULTS_DIR / "depmap_validation.tsv", sep="\t", index=False)

        # Summary by category
        for cat in ["known_driver", "novel"]:
            subset = df_strat[df_strat["category"] == cat]
            n_genes = len(subset)
            n_essential = (subset["essentiality"] < -0.5).sum()
            mean_ess = subset["essentiality"].mean()
            mean_score = subset["predicted_score"].mean()
            print(f"  {cat}: {n_genes} genes, {n_essential} essential, "
                  f"mean_ess={mean_ess:.3f}, mean_score={mean_score:.3f}")

        # Novel driver candidates: high score + essential + NOT known
        novel_essential = df_strat[
            (df_strat["category"] == "novel")
            & (df_strat["essentiality"] < -0.5)
        ].sort_values("predicted_score", ascending=False)
        if len(novel_essential) > 0:
            print(f"\n  Top novel driver candidates (essential in DepMap, not in ground truth):")
            for _, row in novel_essential.head(20).iterrows():
                print(f"    {row['gene']}: score={row['predicted_score']:.4f}, "
                      f"essentiality={row['essentiality']:.3f}")
    else:
        # Save basic results without stratification
        common_genes = sorted(set(gene_scores) & set(depmap.mean_essentiality))
        rows = []
        for gene in common_genes:
            rows.append({
                "gene": gene,
                "predicted_score": gene_scores[gene],
                "essentiality": depmap.mean_essentiality[gene],
            })
        df_basic = pd.DataFrame(rows)
        df_basic.to_csv(RESULTS_DIR / "depmap_validation.tsv", sep="\t", index=False)

    # ── Save summary ──────────────────────────────────────────────
    summary = {
        "spearman_correlation": f"{result.correlation:.6f}",
        "p_value": f"{result.p_value:.2e}",
        "enrichment_top_25pct": f"{result.enrichment_top_k:.6f}",
        "mean_essentiality_top_25pct": f"{result.mean_essentiality_top:.6f}",
        "mean_essentiality_bottom_25pct": f"{result.mean_essentiality_bottom:.6f}",
        "n_genes_evaluated": str(result.n_genes_evaluated),
    }
    with open(RESULTS_DIR / "depmap_summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        for k, v in summary.items():
            writer.writerow([k, v])
    print(f"\n  Saved: depmap_summary.tsv")

    t_total = time.time() - t_start
    print(f"\nDepMap validation completed in {t_total:.1f}s")


if __name__ == "__main__":
    main()
