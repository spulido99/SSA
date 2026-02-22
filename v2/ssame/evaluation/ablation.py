"""Ablation runner for systematic evaluation of configuration variants."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ssame.evaluation.pu_cross_validation import PUCrossValidator


@dataclass
class AblationConfig:
    """Configuration for a single ablation experiment."""

    name: str
    encoder_level: int
    n_estimators: int = 100
    n_folds: int = 5
    seed: int = 42


@dataclass
class AblationResult:
    """Results from a single ablation experiment."""

    name: str
    auprc: float
    auroc: float
    precision_at_10: float
    precision_at_50: float
    recall_at_50: float
    per_fold_auprc: list[float] = field(default_factory=list)
    per_fold_auroc: list[float] = field(default_factory=list)


class AblationRunner:
    """Runs ablation experiments across configuration variants.

    Currently supports Level 0 (XGBoost) ablations. Level 1 (graph encoder)
    ablations require graph data and will be added when graph features are
    integrated into the evaluation pipeline.
    """

    def __init__(self, X: np.ndarray, labels: np.ndarray) -> None:
        self.X = X
        self.labels = labels

    def run_config(self, config: AblationConfig) -> AblationResult:
        """Run a single ablation configuration.

        Args:
            config: Ablation configuration to evaluate.

        Returns:
            AblationResult with cross-validated metrics.
        """
        cv = PUCrossValidator(n_folds=config.n_folds, seed=config.seed)
        cv_result = cv.run(self.X, self.labels, n_estimators=config.n_estimators)

        return AblationResult(
            name=config.name,
            auprc=cv_result.metrics["auprc"],
            auroc=cv_result.metrics["auroc"],
            precision_at_10=cv_result.metrics["precision_at_10"],
            precision_at_50=cv_result.metrics["precision_at_50"],
            recall_at_50=cv_result.metrics["recall_at_50"],
            per_fold_auprc=cv_result.per_fold_auprc,
            per_fold_auroc=cv_result.per_fold_auroc,
        )

    def run_all(self, configs: list[AblationConfig]) -> list[AblationResult]:
        """Run all ablation configurations.

        Args:
            configs: List of configurations to evaluate.

        Returns:
            List of results in same order as configs.
        """
        return [self.run_config(c) for c in configs]

    @staticmethod
    def comparison_table(results: list[AblationResult]) -> str:
        """Format results as a comparison table.

        Args:
            results: List of ablation results.

        Returns:
            Formatted string table.
        """
        header = f"{'Name':<20} {'AUPRC':>8} {'AUROC':>8} {'P@10':>8} {'P@50':>8} {'R@50':>8}"
        separator = "-" * len(header)
        rows = [header, separator]

        for r in results:
            row = (
                f"{r.name:<20} "
                f"{r.auprc:>8.4f} "
                f"{r.auroc:>8.4f} "
                f"{r.precision_at_10:>8.4f} "
                f"{r.precision_at_50:>8.4f} "
                f"{r.recall_at_50:>8.4f}"
            )
            rows.append(row)

        return "\n".join(rows)
