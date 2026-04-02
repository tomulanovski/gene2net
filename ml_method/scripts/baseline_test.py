"""Simple baseline: logistic regression on raw features to predict WGD edges.

If this works well → the GNN architecture has a problem.
If this also fails → features have too much overlap despite different means.

Run on cluster:
    python scripts/baseline_test.py --dir /path/to/mul_trees_2k/training/ils_low
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", required=True)
    parser.add_argument("--max-samples", type=int, default=500)
    args = parser.parse_args()

    dataset = Gene2NetDataset(args.dir)
    n_total = min(len(dataset), args.max_samples)
    print(f"Loading {n_total} samples...")

    X_all = []
    y_all = []

    for i in range(n_total):
        try:
            sample = dataset[i]
        except Exception:
            continue

        if sample.labels is None:
            continue

        node_feats = sample.species_tree_node_features.numpy()   # [N, 13]
        edge_feats = sample.species_tree_edge_features.numpy()   # [E, 4]
        edge_index = sample.species_tree_edge_index              # [2, 2*E]
        wgd_counts = np.array(sample.labels.wgd_counts)
        n_edges = sample.labels.n_edges

        if edge_feats.shape[0] != n_edges:
            continue

        for e_idx in range(n_edges):
            # Parent and child node indices (from directed edge pairs)
            parent_idx = int(edge_index[0, 2 * e_idx])
            child_idx = int(edge_index[1, 2 * e_idx])

            parent_feat = node_feats[parent_idx]  # 13-dim
            child_feat = node_feats[child_idx]     # 13-dim
            e_feat = edge_feats[e_idx]             # 4-dim

            # Concatenate: parent_node(13) + child_node(13) + edge(4) = 30 features
            feat = np.concatenate([parent_feat, child_feat, e_feat])
            X_all.append(feat)
            y_all.append(1 if wgd_counts[e_idx] > 0 else 0)

        if (i + 1) % 100 == 0:
            print(f"  Loaded {i+1}/{n_total}")

    X = np.array(X_all)
    y = np.array(y_all)

    print(f"\nTotal edges: {len(y)}")
    print(f"  Positive (WGD): {y.sum()} ({100*y.mean():.1f}%)")
    print(f"  Negative: {(1-y).sum()} ({100*(1-y.mean()):.1f}%)")

    # Split train/test (80/20)
    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.metrics import f1_score, precision_score, recall_score, classification_report
    from sklearn.preprocessing import StandardScaler

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    print(f"\nTrain: {len(y_train)} edges ({y_train.sum()} positive)")
    print(f"Test:  {len(y_test)} edges ({y_test.sum()} positive)")

    # --- Baseline 1: Logistic Regression ---
    print(f"\n{'='*60}")
    print("BASELINE 1: Logistic Regression (class_weight='balanced')")
    print(f"{'='*60}")
    lr = LogisticRegression(class_weight="balanced", max_iter=1000, random_state=42)
    lr.fit(X_train_s, y_train)
    y_pred_lr = lr.predict(X_test_s)
    print(classification_report(y_test, y_pred_lr, target_names=["no_WGD", "WGD"]))

    # Feature importances from logistic regression
    feat_names = (
        [f"parent_{n}" for n in ["mean_cp", "var_cp", "mode_cp", "p_absent",
         "p_1cp", "p_2cp", "p_3+cp", "max_cp",
         "cl_mean", "cl_std", "cl_max", "cl_min", "cl_med"]] +
        [f"child_{n}" for n in ["mean_cp", "var_cp", "mode_cp", "p_absent",
         "p_1cp", "p_2cp", "p_3+cp", "max_cp",
         "cl_mean", "cl_std", "cl_max", "cl_min", "cl_med"]] +
        ["concordance", "branch_len", "clade_size", "depth"]
    )
    coefs = lr.coef_[0]
    sorted_idx = np.argsort(np.abs(coefs))[::-1]
    print("Top 10 features (by |coefficient|):")
    for rank, idx in enumerate(sorted_idx[:10]):
        print(f"  {rank+1}. {feat_names[idx]:<20} coef={coefs[idx]:+.4f}")

    # --- Baseline 2: Random Forest ---
    print(f"\n{'='*60}")
    print("BASELINE 2: Random Forest (class_weight='balanced')")
    print(f"{'='*60}")
    rf = RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    y_pred_rf = rf.predict(X_test)
    print(classification_report(y_test, y_pred_rf, target_names=["no_WGD", "WGD"]))

    # Feature importances
    importances = rf.feature_importances_
    sorted_idx = np.argsort(importances)[::-1]
    print("Top 10 features (by importance):")
    for rank, idx in enumerate(sorted_idx[:10]):
        print(f"  {rank+1}. {feat_names[idx]:<20} importance={importances[idx]:.4f}")

    # --- Baseline 3: Gradient Boosting ---
    print(f"\n{'='*60}")
    print("BASELINE 3: Gradient Boosting")
    print(f"{'='*60}")
    # Compute scale_pos_weight for imbalanced classes
    n_neg = (y_train == 0).sum()
    n_pos = (y_train == 1).sum()
    gb = GradientBoostingClassifier(
        n_estimators=200, max_depth=5, random_state=42,
        # No built-in class_weight, but sample_weight equivalent
    )
    sample_weights = np.where(y_train == 1, n_neg / n_pos, 1.0)
    gb.fit(X_train, y_train, sample_weight=sample_weights)
    y_pred_gb = gb.predict(X_test)
    print(classification_report(y_test, y_pred_gb, target_names=["no_WGD", "WGD"]))

    # --- Summary ---
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for name, y_pred in [("Logistic Regression", y_pred_lr),
                          ("Random Forest", y_pred_rf),
                          ("Gradient Boosting", y_pred_gb)]:
        f1 = f1_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred)
        print(f"  {name:<25} F1={f1:.3f}  Precision={prec:.3f}  Recall={rec:.3f}")

    print(f"\nIf baselines achieve F1 > 0.3, the GNN architecture needs fixing.")
    print(f"If baselines also fail (F1 < 0.1), features have too much class overlap.")

    # Save
    out_path = os.path.join(args.dir, "baseline_results.txt")
    # Re-run and capture (simple approach: just note the scores)
    with open(out_path, "w") as f:
        f.write(f"Baseline results on {len(y_test)} test edges\n")
        for name, y_pred in [("Logistic Regression", y_pred_lr),
                              ("Random Forest", y_pred_rf),
                              ("Gradient Boosting", y_pred_gb)]:
            f1 = f1_score(y_test, y_pred)
            prec = precision_score(y_test, y_pred, zero_division=0)
            rec = recall_score(y_test, y_pred)
            f.write(f"{name}: F1={f1:.3f} Precision={prec:.3f} Recall={rec:.3f}\n")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
