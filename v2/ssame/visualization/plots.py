"""Data preparation for visualization plots.

Prepares structured data for plotting (matplotlib/seaborn agnostic).
Actual plotting is deferred to notebooks or scripts.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from sklearn.metrics import average_precision_score, precision_recall_curve


def prepare_rank_comparison_data(
    gene_scores: dict[str, float],
    positives: set[str],
) -> dict:
    """Prepare ranked gene list with positive/unlabeled flags.

    Returns dict with keys: genes, ranks, scores, is_positive.
    """
    sorted_genes = sorted(gene_scores.items(), key=lambda x: -x[1])

    return {
        "genes": [g for g, _ in sorted_genes],
        "scores": [s for _, s in sorted_genes],
        "ranks": list(range(len(sorted_genes))),
        "is_positive": [g in positives for g, _ in sorted_genes],
    }


def prepare_precision_recall_data(
    scores: np.ndarray,
    labels: np.ndarray,
) -> dict:
    """Compute precision-recall curve data.

    Returns dict with keys: precision, recall, thresholds, auprc.
    """
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    auprc = average_precision_score(labels, scores)

    return {
        "precision": precision.tolist(),
        "recall": recall.tolist(),
        "thresholds": thresholds.tolist(),
        "auprc": float(auprc),
    }


def prepare_bias_analysis_data(
    gene_scores: dict[str, float],
    string_degree: dict[str, int],
    alphafold_degree: dict[str, int],
) -> dict:
    """Analyze study bias: correlation between scores and network degree.

    High correlation with STRING degree = study bias.
    Low correlation with AlphaFold degree = less biased predictions.
    """
    common_genes = (
        set(gene_scores) & set(string_degree) & set(alphafold_degree)
    )
    genes = sorted(common_genes)

    if len(genes) < 3:
        return {
            "string_correlation": 0.0,
            "alphafold_correlation": 0.0,
            "n_genes": len(genes),
        }

    scores = np.array([gene_scores[g] for g in genes])
    s_deg = np.array([string_degree[g] for g in genes], dtype=float)
    a_deg = np.array([alphafold_degree[g] for g in genes], dtype=float)

    string_corr, string_pval = stats.spearmanr(scores, s_deg)
    alpha_corr, alpha_pval = stats.spearmanr(scores, a_deg)

    return {
        "string_correlation": float(string_corr),
        "string_pvalue": float(string_pval),
        "alphafold_correlation": float(alpha_corr),
        "alphafold_pvalue": float(alpha_pval),
        "n_genes": len(genes),
        "genes": genes,
    }
