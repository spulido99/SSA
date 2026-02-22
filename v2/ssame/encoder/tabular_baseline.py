"""Level 0 tabular baseline: PU-aware XGBoost/LightGBM encoder."""

from __future__ import annotations

import numpy as np
import xgboost as xgb


class TabularBaseline:
    """PU-aware XGBoost baseline for gene driverness scoring.

    Uses PU bagging: trains multiple models on subsamples where
    unlabeled genes are treated as negative, then averages predictions.
    This reduces false negative bias from the PU assumption.
    """

    def __init__(
        self,
        n_estimators: int = 100,
        max_depth: int = 6,
        learning_rate: float = 0.1,
        n_pu_bags: int = 1,
        seed: int = 42,
    ) -> None:
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.n_pu_bags = n_pu_bags
        self.seed = seed
        self._models: list[xgb.XGBClassifier] = []

    def fit(self, X: np.ndarray, labels: np.ndarray) -> None:
        """Fit the model using PU bagging.

        Args:
            X: Feature matrix (n_genes, n_features).
            labels: 1 for known positives, 0 for unlabeled.
        """
        rng = np.random.default_rng(self.seed)
        pos_mask = labels == 1
        unl_mask = labels == 0
        pos_idx = np.where(pos_mask)[0]
        unl_idx = np.where(unl_mask)[0]

        self._models = []
        for i in range(self.n_pu_bags):
            model = xgb.XGBClassifier(
                n_estimators=self.n_estimators,
                max_depth=self.max_depth,
                learning_rate=self.learning_rate,
                random_state=self.seed + i,
                eval_metric="logloss",
                verbosity=0,
            )

            if self.n_pu_bags > 1:
                # PU bagging: subsample unlabeled set
                n_sample = min(len(pos_idx) * 5, len(unl_idx))
                sampled_unl = rng.choice(unl_idx, size=n_sample, replace=False)
                train_idx = np.concatenate([pos_idx, sampled_unl])
            else:
                train_idx = np.arange(len(labels))

            # Weight positives more heavily to compensate for label noise
            sample_weight = np.ones(len(train_idx))
            pos_in_train = np.isin(train_idx, pos_idx)
            n_pos = pos_in_train.sum()
            n_neg = (~pos_in_train).sum()
            if n_pos > 0 and n_neg > 0:
                sample_weight[pos_in_train] = n_neg / n_pos

            model.fit(
                X[train_idx],
                labels[train_idx],
                sample_weight=sample_weight,
            )
            self._models.append(model)

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict driverness scores in [0, 1].

        Args:
            X: Feature matrix (n_genes, n_features).

        Returns:
            Per-gene driverness scores averaged across PU bags.
        """
        if not self._models:
            raise RuntimeError("Model not fitted. Call fit() first.")

        all_probs = np.stack([
            m.predict_proba(X)[:, 1] for m in self._models
        ])
        return all_probs.mean(axis=0)

    def feature_importance(self) -> np.ndarray:
        """Return averaged feature importance across PU bags."""
        if not self._models:
            raise RuntimeError("Model not fitted. Call fit() first.")

        importances = np.stack([
            m.feature_importances_ for m in self._models
        ])
        return importances.mean(axis=0)
