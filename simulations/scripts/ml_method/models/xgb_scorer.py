"""
XGBoost-based bipartition scorer.

Learns to predict whether a bipartition is in the true MUL-tree
based on features computed from gene trees.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from pathlib import Path
import pickle

try:
    import xgboost as xgb
except ImportError:
    raise ImportError("Please install xgboost: pip install xgboost")

from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import (
    precision_recall_curve,
    average_precision_score,
    roc_auc_score,
    classification_report,
    f1_score
)


class BipartitionScorer:
    """
    XGBoost classifier for scoring bipartitions.

    Predicts P(bipartition ∈ true MUL-tree | features).
    """

    def __init__(
        self,
        n_estimators: int = 100,
        max_depth: int = 6,
        learning_rate: float = 0.1,
        min_child_weight: int = 1,
        subsample: float = 0.8,
        colsample_bytree: float = 0.8,
        random_state: int = 42,
        **kwargs
    ):
        """
        Initialize the scorer.

        Args:
            n_estimators: Number of boosting rounds
            max_depth: Maximum tree depth
            learning_rate: Boosting learning rate
            min_child_weight: Minimum sum of instance weight in a child
            subsample: Subsample ratio of training instances
            colsample_bytree: Subsample ratio of columns
            random_state: Random seed
            **kwargs: Additional arguments passed to XGBClassifier
        """
        self.model = xgb.XGBClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            min_child_weight=min_child_weight,
            subsample=subsample,
            colsample_bytree=colsample_bytree,
            random_state=random_state,
            use_label_encoder=False,
            eval_metric='logloss',
            **kwargs
        )
        self.feature_names: Optional[List[str]] = None
        self.is_fitted = False

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        feature_names: Optional[List[str]] = None,
        sample_weight: Optional[np.ndarray] = None,
        verbose: bool = True
    ) -> 'BipartitionScorer':
        """
        Train the scorer.

        Args:
            X: Feature matrix [n_samples, n_features]
            y: Labels [n_samples] (1 = in true tree, 0 = not)
            feature_names: Names of features (for interpretability)
            sample_weight: Sample weights (optional)
            verbose: Print training info

        Returns:
            self
        """
        self.feature_names = feature_names

        # Handle class imbalance via scale_pos_weight
        n_pos = y.sum()
        n_neg = len(y) - n_pos
        scale_pos_weight = n_neg / n_pos if n_pos > 0 else 1.0

        if verbose:
            print(f"Training on {len(y)} samples ({n_pos:.0f} positive, {n_neg:.0f} negative)")
            print(f"Using scale_pos_weight={scale_pos_weight:.2f}")

        self.model.set_params(scale_pos_weight=scale_pos_weight)
        self.model.fit(X, y, sample_weight=sample_weight)
        self.is_fitted = True

        if verbose:
            print("Training complete.")

        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """
        Predict probability of bipartition being in true tree.

        Args:
            X: Feature matrix [n_samples, n_features]

        Returns:
            Probability of positive class [n_samples]
        """
        if not self.is_fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        return self.model.predict_proba(X)[:, 1]

    def predict(self, X: np.ndarray, threshold: float = 0.5) -> np.ndarray:
        """
        Predict binary labels.

        Args:
            X: Feature matrix
            threshold: Classification threshold

        Returns:
            Binary predictions [n_samples]
        """
        probs = self.predict_proba(X)
        return (probs >= threshold).astype(int)

    def get_feature_importance(
        self,
        importance_type: str = 'weight'
    ) -> Dict[str, float]:
        """
        Get feature importance scores.

        Args:
            importance_type: 'weight', 'gain', or 'cover'

        Returns:
            Dict mapping feature name to importance score
        """
        if not self.is_fitted:
            raise RuntimeError("Model not fitted.")

        importance = self.model.get_booster().get_score(importance_type=importance_type)

        # Map to feature names
        if self.feature_names:
            # XGBoost uses f0, f1, ... as default names
            result = {}
            for key, value in importance.items():
                if key.startswith('f'):
                    idx = int(key[1:])
                    if idx < len(self.feature_names):
                        result[self.feature_names[idx]] = value
                    else:
                        result[key] = value
                else:
                    result[key] = value
            return result
        else:
            return importance

    def cross_validate(
        self,
        X: np.ndarray,
        y: np.ndarray,
        n_splits: int = 5,
        verbose: bool = True
    ) -> Dict[str, float]:
        """
        Perform cross-validation and return metrics.

        Args:
            X: Feature matrix
            y: Labels
            n_splits: Number of CV folds
            verbose: Print results

        Returns:
            Dict of metric name to value
        """
        # Clone model for CV
        cv_model = xgb.XGBClassifier(**self.model.get_params())

        # Handle class imbalance
        n_pos = y.sum()
        n_neg = len(y) - n_pos
        cv_model.set_params(scale_pos_weight=n_neg / n_pos if n_pos > 0 else 1.0)

        # Stratified k-fold
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

        # Get cross-val predictions
        y_pred_proba = cross_val_predict(cv_model, X, y, cv=cv, method='predict_proba')[:, 1]
        y_pred = (y_pred_proba >= 0.5).astype(int)

        # Compute metrics
        metrics = {
            'roc_auc': roc_auc_score(y, y_pred_proba),
            'average_precision': average_precision_score(y, y_pred_proba),
            'f1': f1_score(y, y_pred),
            'precision': (y_pred & y.astype(bool)).sum() / y_pred.sum() if y_pred.sum() > 0 else 0,
            'recall': (y_pred & y.astype(bool)).sum() / y.sum() if y.sum() > 0 else 0,
        }

        if verbose:
            print(f"\n{n_splits}-Fold Cross-Validation Results:")
            for name, value in metrics.items():
                print(f"  {name}: {value:.4f}")

        return metrics

    def find_optimal_threshold(
        self,
        X: np.ndarray,
        y: np.ndarray,
        metric: str = 'f1'
    ) -> Tuple[float, float]:
        """
        Find optimal classification threshold.

        Args:
            X: Feature matrix
            y: Labels
            metric: Metric to optimize ('f1' or 'precision' or 'recall')

        Returns:
            (optimal_threshold, metric_value)
        """
        probs = self.predict_proba(X)
        precisions, recalls, thresholds = precision_recall_curve(y, probs)

        if metric == 'f1':
            f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)
            best_idx = np.argmax(f1_scores)
            best_threshold = thresholds[best_idx] if best_idx < len(thresholds) else 0.5
            return best_threshold, f1_scores[best_idx]

        elif metric == 'precision':
            # Find threshold that gives precision > 0.8 with highest recall
            valid_idx = np.where(precisions >= 0.8)[0]
            if len(valid_idx) == 0:
                return 0.9, precisions.max()
            best_idx = valid_idx[np.argmax(recalls[valid_idx])]
            return thresholds[best_idx] if best_idx < len(thresholds) else 0.5, precisions[best_idx]

        else:
            raise ValueError(f"Unknown metric: {metric}")

    def save(self, filepath: str):
        """Save model to file."""
        with open(filepath, 'wb') as f:
            pickle.dump({
                'model': self.model,
                'feature_names': self.feature_names,
                'is_fitted': self.is_fitted,
            }, f)

    @classmethod
    def load(cls, filepath: str) -> 'BipartitionScorer':
        """Load model from file."""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)

        scorer = cls()
        scorer.model = data['model']
        scorer.feature_names = data['feature_names']
        scorer.is_fitted = data['is_fitted']
        return scorer


def train_scorer(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: Optional[np.ndarray] = None,
    y_val: Optional[np.ndarray] = None,
    feature_names: Optional[List[str]] = None,
    **model_kwargs
) -> Tuple[BipartitionScorer, Dict[str, Any]]:
    """
    Train a bipartition scorer with optional validation.

    Args:
        X_train: Training features
        y_train: Training labels
        X_val: Validation features (optional)
        y_val: Validation labels (optional)
        feature_names: Feature names
        **model_kwargs: Arguments for BipartitionScorer

    Returns:
        (trained_scorer, results_dict)
    """
    scorer = BipartitionScorer(**model_kwargs)
    scorer.fit(X_train, y_train, feature_names=feature_names)

    results = {
        'train_metrics': {},
        'val_metrics': {},
        'feature_importance': scorer.get_feature_importance()
    }

    # Evaluate on training data
    train_probs = scorer.predict_proba(X_train)
    results['train_metrics'] = {
        'roc_auc': roc_auc_score(y_train, train_probs),
        'average_precision': average_precision_score(y_train, train_probs),
    }

    # Evaluate on validation data if provided
    if X_val is not None and y_val is not None:
        val_probs = scorer.predict_proba(X_val)
        results['val_metrics'] = {
            'roc_auc': roc_auc_score(y_val, val_probs),
            'average_precision': average_precision_score(y_val, val_probs),
        }

    return scorer, results


# Command-line interface for testing
if __name__ == "__main__":
    # Simple test with synthetic data
    print("Testing BipartitionScorer with synthetic data...")

    np.random.seed(42)

    # Generate synthetic data
    n_samples = 1000
    n_features = 18

    X = np.random.randn(n_samples, n_features)
    # Make first feature (frequency) most predictive
    y = ((X[:, 0] > 0) & (X[:, 1] > -0.5)).astype(int)

    print(f"Data shape: {X.shape}")
    print(f"Positive rate: {y.mean():.1%}")

    # Test cross-validation
    scorer = BipartitionScorer(n_estimators=50, max_depth=4)
    metrics = scorer.cross_validate(X, y)

    # Test full training
    scorer.fit(X, y)

    # Test feature importance
    print("\nFeature importance (top 5):")
    importance = scorer.get_feature_importance()
    sorted_imp = sorted(importance.items(), key=lambda x: -x[1])[:5]
    for name, value in sorted_imp:
        print(f"  {name}: {value:.2f}")

    # Test optimal threshold
    best_thresh, best_f1 = scorer.find_optimal_threshold(X, y)
    print(f"\nOptimal threshold: {best_thresh:.3f} (F1={best_f1:.3f})")

    print("\nTest passed!")
