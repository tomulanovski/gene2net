#!/usr/bin/env python3
"""
compute_pairwise_methods.py - Pairwise method-vs-method comparison for simulations

Computes pairwise disagreement metrics between Polyphest and grandma_split
per (network, replicate), per config. Mirrors the empirical pipeline in
scripts/compute_comparisons.py so simulation disagreement levels are directly
comparable to the empirical pairwise numbers.

Output (per config):
    simulations/analysis/summary/{config}/pairwise_polyphest_vs_grandma.csv
        Columns: config, network, replicate, metric, value, status

Usage:
    python compute_pairwise_methods.py CONFIG [CONFIG ...]
    python compute_pairwise_methods.py conf_ils_low_10M conf_ils_medium_10M
    python compute_pairwise_methods.py conf_ils_low_10M --force-recompute
"""

import argparse
import hashlib
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from reticulate_tree import ReticulateTree
from compare_reticulations import pairwise_compare
from create_analysis_figures import merge_polyphest_inventory

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"

METHOD_A = "polyphest"
METHOD_B = "grandma_split"

# Output filename is derived at runtime from the method pair


def file_hash(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except Exception:
        return "missing"


def load_tree(path):
    with open(path, "r") as f:
        newick = f.read().strip()
    if newick and not newick.endswith(";") and not newick.endswith("nexml>"):
        newick += ";"
    return ReticulateTree(newick)


def compare_pair(path_a, path_b):
    """Run pairwise_compare on two network files. Returns dict with status/error/metrics."""
    try:
        tree_a = load_tree(path_a)
    except Exception as e:
        return {"status": "ERROR", "error": f"load {METHOD_A}: {e}", "metrics": None}
    try:
        tree_b = load_tree(path_b)
    except Exception as e:
        return {"status": "ERROR", "error": f"load {METHOD_B}: {e}", "metrics": None}
    try:
        # Neither method is in SINGLE_RETICULATION_METHODS, so partial_match=False
        metrics = pairwise_compare(tree_a, tree_b, partial_match=False)
        return {"status": "SUCCESS", "error": None, "metrics": metrics}
    except Exception as e:
        return {"status": "ERROR", "error": str(e), "metrics": None}


def flatten_metrics(metrics):
    flat = {}
    for k, v in metrics.items():
        if isinstance(v, dict):
            for sk, sv in v.items():
                flat[f"{k}.{sk}"] = sv
        else:
            flat[k] = v
    return flat


def process_config(config, summary_base, method_a=None, method_b=None, force_recompute=False):
    method_a = method_a or METHOD_A
    method_b = method_b or METHOD_B
    config_dir = summary_base / config
    inv_path = config_dir / "inventory.csv"
    if not inv_path.exists():
        print(f"  [SKIP] {config}: inventory.csv not found")
        return None

    inventory = pd.read_csv(inv_path)
    inventory = merge_polyphest_inventory(inventory)
    inventory = inventory[inventory["method"].isin([method_a, method_b])].copy()

    cache_dir = config_dir / "pairwise_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    n_total = n_success = n_failed = n_cached = 0

    # Group by (network, replicate); keep pairs where both methods completed
    for (network, replicate), group in inventory.groupby(["network", "replicate"]):
        a_row = group[(group["method"] == method_a) & (group["inferred_exists"])]
        b_row = group[(group["method"] == method_b) & (group["inferred_exists"])]
        if a_row.empty or b_row.empty:
            continue

        path_a = a_row.iloc[0]["inferred_path"]
        path_b = b_row.iloc[0]["inferred_path"]
        n_total += 1

        cache_path = cache_dir / f"{network}_rep{replicate}_{method_a}_vs_{method_b}.pkl"
        comparison = None

        if not force_recompute and cache_path.exists():
            try:
                with open(cache_path, "rb") as f:
                    cached = pickle.load(f)
                if (
                    cached.get("hash_a") == file_hash(path_a)
                    and cached.get("hash_b") == file_hash(path_b)
                ):
                    comparison = cached["comparison"]
                    n_cached += 1
            except Exception:
                comparison = None

        if comparison is None:
            comparison = compare_pair(path_a, path_b)
            if comparison["status"] == "SUCCESS":
                try:
                    with open(cache_path, "wb") as f:
                        pickle.dump(
                            {
                                "network": network,
                                "replicate": replicate,
                                "hash_a": file_hash(path_a),
                                "hash_b": file_hash(path_b),
                                "timestamp": datetime.now().isoformat(),
                                "comparison": comparison,
                            },
                            f,
                        )
                except Exception:
                    pass

        if comparison["status"] == "SUCCESS":
            n_success += 1
            for metric, value in flatten_metrics(comparison["metrics"]).items():
                rows.append(
                    {
                        "config": config,
                        "network": network,
                        "replicate": replicate,
                        "metric": metric,
                        "value": value,
                        "status": "SUCCESS",
                    }
                )
        else:
            n_failed += 1
            rows.append(
                {
                    "config": config,
                    "network": network,
                    "replicate": replicate,
                    "metric": "FAILED",
                    "value": np.nan,
                    "status": "FAILED",
                }
            )
            print(f"    [FAIL] {network} rep{replicate}: {comparison['error']}")

    df = pd.DataFrame(rows)
    out_name = f"pairwise_{method_a}_vs_{method_b}.csv"
    out_path = config_dir / out_name
    df.to_csv(out_path, index=False)
    print(
        f"  [{config}] pairs={n_total} success={n_success} failed={n_failed} "
        f"cached={n_cached} -> {out_path}"
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Pairwise method comparison for simulations"
    )
    parser.add_argument("configs", nargs="+", help="One or more config names")
    parser.add_argument("--method-a", default=METHOD_A, help=f"First method (default: {METHOD_A})")
    parser.add_argument("--method-b", default=METHOD_B, help=f"Second method (default: {METHOD_B})")
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Ignore cache and recompute all comparisons",
    )
    args = parser.parse_args()

    print(f"Pairwise comparison: {args.method_a} vs {args.method_b}")
    print(f"Configs: {', '.join(args.configs)}\n")

    for config in args.configs:
        process_config(config, SUMMARY_BASE, method_a=args.method_a, method_b=args.method_b,
                       force_recompute=args.force_recompute)


if __name__ == "__main__":
    main()
