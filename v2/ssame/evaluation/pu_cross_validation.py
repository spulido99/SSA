"""PU-aware cross-validation for cancer driver prediction."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from ssame.encoder.tabular_baseline import TabularBaseline


@dataclass
class CVResult:
    """Cross-validation results."""

    metrics: dict[str, float]
    per_fold_auprc: list[float] = field(default_factory=list)
    per_fold_auroc: list[float] = field(default_factory=list)
    gene_scores: dict[int, float] = field(default_factory=dict)


def precision_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Precision at k: fraction of top-k predictions that are true positives."""
    top_k_idx = np.argsort(-scores)[:k]
    return float(labels[top_k_idx].sum()) / k


def recall_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Recall at k: fraction of true positives found in top-k predictions."""
    n_pos = labels.sum()
    if n_pos == 0:
        return 0.0
    top_k_idx = np.argsort(-scores)[:k]
    return float(labels[top_k_idx].sum()) / n_pos


class PUCrossValidator:
    """PU-aware cross-validation.

    Stratifies on positive labels only (PU framing: unlabeled != negative).
    Positives are split across folds; all unlabeled go into training for
    each fold (since we don't know which unlabeled are truly negative).
    """

    def __init__(self, n_folds: int = 5, seed: int = 42) -> None:
        self.n_folds = n_folds
        self.seed = seed

    def _make_folds(
        self, labels: np.ndarray
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        """Create PU-aware fold splits.

        Positives are split into n_folds groups.
        For each fold: test = fold's positives + all unlabeled.
        Training = remaining positives + all unlabeled.
        """
        rng = np.random.default_rng(self.seed)
        pos_idx = np.where(labels == 1)[0]
        unl_idx = np.where(labels == 0)[0]

        # Shuffle positives
        shuffled_pos = rng.permutation(pos_idx)
        fold_sizes = np.array_split(shuffled_pos, self.n_folds)

        folds = []
        for i in range(self.n_folds):
            test_pos = fold_sizes[i]
            train_pos = np.concatenate([fold_sizes[j] for j in range(self.n_folds) if j != i])

            # Training: remaining positives + all unlabeled
            train_idx = np.concatenate([train_pos, unl_idx])
            # Test: held-out positives + all unlabeled (for ranking evaluation)
            test_idx = np.concatenate([test_pos, unl_idx])

            folds.append((train_idx, test_idx))

        return folds

    def run(
        self,
        X: np.ndarray,
        labels: np.ndarray,
        *,
        n_estimators: int = 100,
    ) -> CVResult:
        """Run PU-aware cross-validation.

        Args:
            X: Feature matrix (n_genes, n_features).
            labels: PU labels (1=known positive, 0=unlabeled).
            n_estimators: Number of trees for XGBoost.

        Returns:
            CVResult with aggregated metrics.
        """
        folds = self._make_folds(labels)
        all_scores = np.zeros(len(labels))
        score_counts = np.zeros(len(labels))

        fold_auprcs = []
        fold_aurocs = []

        for fold_idx, (train_idx, test_idx) in enumerate(folds):
            model = TabularBaseline(
                n_estimators=n_estimators,
                seed=self.seed + fold_idx,
            )
            model.fit(X[train_idx], labels[train_idx])
            test_scores = model.predict(X[test_idx])

            # Accumulate scores
            for idx, score in zip(test_idx, test_scores, strict=True):
                all_scores[idx] += score
                score_counts[idx] += 1

            # Compute per-fold metrics on test set
            test_labels = labels[test_idx]
            if test_labels.sum() > 0 and test_labels.sum() < len(test_labels):
                auprc = average_precision_score(test_labels, test_scores)
                auroc = roc_auc_score(test_labels, test_scores)
                fold_auprcs.append(auprc)
                fold_aurocs.append(auroc)

        # Average scores for genes seen in multiple folds
        mask = score_counts > 0
        all_scores[mask] /= score_counts[mask]

        # Overall metrics
        pos_mask = labels == 1
        if pos_mask.sum() > 0 and (~pos_mask).sum() > 0:
            overall_auprc = average_precision_score(labels, all_scores)
            overall_auroc = roc_auc_score(labels, all_scores)
        else:
            overall_auprc = 0.0
            overall_auroc = 0.0

        return CVResult(
            metrics={
                "auprc": float(overall_auprc),
                "auroc": float(overall_auroc),
                "precision_at_10": precision_at_k(all_scores, labels, k=10),
                "precision_at_50": precision_at_k(all_scores, labels, k=min(50, len(labels))),
                "recall_at_50": recall_at_k(all_scores, labels, k=min(50, len(labels))),
            },
            per_fold_auprc=fold_auprcs,
            per_fold_auroc=fold_aurocs,
            gene_scores={i: float(s) for i, s in enumerate(all_scores)},
        )
