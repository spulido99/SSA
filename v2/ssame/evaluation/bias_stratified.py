"""Bias-stratified evaluation: measure performance across study-bias strata."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class StratifiedResult:
    """Results from bias-stratified evaluation."""

    strata_metrics: list[dict]
    summary: dict[str, float]


class BiasStratifiedEvaluator:
    """Evaluate predictions stratified by study bias proxy.

    Key idea: if a method is biased, it will perform well on well-studied
    genes but poorly on under-studied ones. An unbiased method should
    perform uniformly across strata.

    Study bias proxy: publication count, STRING degree, or PubMed hits.
    """

    def __init__(
        self,
        gene_scores: dict[str, float],
        study_bias_proxy: dict[str, float | int],
        labels: dict[str, int],
        n_strata: int = 3,
    ) -> None:
        self.gene_scores = gene_scores
        self.study_bias_proxy = study_bias_proxy
        self.labels = labels
        self.n_strata = n_strata

    def evaluate(self) -> StratifiedResult:
        """Run bias-stratified evaluation."""
        common = sorted(
            set(self.gene_scores) & set(self.study_bias_proxy)
        )

        scores = np.array([self.gene_scores[g] for g in common])
        bias = np.array([float(self.study_bias_proxy[g]) for g in common])
        is_positive = np.array([self.labels.get(g, 0) for g in common])

        # Score-bias correlation
        corr, p_val = stats.spearmanr(scores, bias)

        # Stratify by bias proxy quantiles
        quantiles = np.linspace(0, 100, self.n_strata + 1)
        thresholds = np.percentile(bias, quantiles)

        strata_metrics = []
        for i in range(self.n_strata):
            low = thresholds[i]
            high = thresholds[i + 1]
            if i == self.n_strata - 1:
                mask = (bias >= low) & (bias <= high)
            else:
                mask = (bias >= low) & (bias < high)

            stratum_scores = scores[mask]
            stratum_labels = is_positive[mask]
            stratum_bias = bias[mask]
            stratum_genes = [g for g, m in zip(common, mask, strict=True) if m]

            n_pos = int(stratum_labels.sum())
            n_total = len(stratum_scores)
            mean_score = float(np.mean(stratum_scores)) if n_total > 0 else 0.0

            # Precision at k (k = number of positives in stratum)
            p_at_k = 0.0
            if n_pos > 0 and n_total > 0:
                top_k_idx = np.argsort(-stratum_scores)[:n_pos]
                p_at_k = float(stratum_labels[top_k_idx].sum() / n_pos)

            strata_metrics.append({
                "stratum": i,
                "n_genes": n_total,
                "n_positives": n_pos,
                "mean_study_bias": float(np.mean(stratum_bias)) if n_total > 0 else 0.0,
                "mean_score": mean_score,
                "precision_at_k": p_at_k,
                "genes": stratum_genes,
            })

        summary = {
            "score_bias_correlation": float(corr),
            "score_bias_pvalue": float(p_val),
            "n_genes": len(common),
            "n_strata": self.n_strata,
        }

        return StratifiedResult(
            strata_metrics=strata_metrics,
            summary=summary,
        )
