#!/usr/bin/env python
"""Evaluate Gene2Net-GNN predictions against ground truth MUL-trees.

Usage:
    # Evaluate on a single biological network
    python evaluate.py --model model.pt --network Shahrestani_2015 \
        --gene-trees /path/to/gene_trees.nwk --species-tree /path/to/astral.nwk

    # Evaluate on all biological networks (requires gene trees for each)
    python evaluate.py --model model.pt --data-dir /path/to/bio_data/ \
        --networks-dir ../simulations/networks/

    # Quick test with mock data (randomly initialized model)
    python evaluate.py --mock --networks-dir ../simulations/networks/
"""
import argparse
import csv
import json
import os
import sys
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree
from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.inference.predict import predict_mul_tree
from gene2net_gnn.data.tree_io import load_gene_trees_from_file, get_species_set
from gene2net_gnn.data.mul_tree_generator import get_polyploid_species
from gene2net_gnn.data.label_extractor import extract_backbone


def load_ground_truth_stats(stats_csv):
    """Load network statistics from mul_tree_final_stats.csv."""
    stats = {}
    with open(stats_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("Network", row.get("network", ""))
            stats[name] = {
                "H_Strict": int(row.get("H_Strict", 0)),
                "Num_Polyploids": int(row.get("Num_Polyploids", 0)),
                "Total_WGD": int(row.get("Total_WGD", 0)),
            }
    return stats


def evaluate_prediction(pred_mul_tree, true_mul_tree, species_tree):
    """Compare predicted MUL-tree against ground truth.

    Returns dict of metrics.
    """
    metrics = {}

    # Ploidy comparison
    pred_polyploids = get_polyploid_species(pred_mul_tree)
    true_polyploids = get_polyploid_species(true_mul_tree)

    all_species = set(species_tree.get_leaf_names())

    # Per-species ploidy accuracy
    correct = 0
    total = len(all_species)
    for sp in all_species:
        pred_ploidy = pred_polyploids.get(sp, 1)
        true_ploidy = true_polyploids.get(sp, 1)
        if pred_ploidy == true_ploidy:
            correct += 1
    metrics["ploidy_accuracy"] = correct / total if total > 0 else 0.0

    # Event count comparison
    pred_n_events = len([sp for sp, c in pred_polyploids.items() if c > 1])
    true_n_events = len([sp for sp, c in true_polyploids.items() if c > 1])
    metrics["pred_n_polyploids"] = len(pred_polyploids)
    metrics["true_n_polyploids"] = len(true_polyploids)
    metrics["polyploid_count_error"] = abs(len(pred_polyploids) - len(true_polyploids))

    # Polyploid detection F1
    pred_set = set(pred_polyploids.keys())
    true_set = set(true_polyploids.keys())
    tp = len(pred_set & true_set)
    fp = len(pred_set - true_set)
    fn = len(true_set - pred_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    metrics["polyploid_precision"] = precision
    metrics["polyploid_recall"] = recall
    metrics["polyploid_f1"] = f1

    # Leaf count comparison
    metrics["pred_n_leaves"] = len(pred_mul_tree.get_leaves())
    metrics["true_n_leaves"] = len(true_mul_tree.get_leaves())

    return metrics


def evaluate_single_network(
    model, network_name, gene_trees_path, species_tree_path,
    ground_truth_path, device, threshold=0.5, max_events=None
):
    """Evaluate on a single biological network."""
    # Load data
    gene_trees = load_gene_trees_from_file(gene_trees_path)
    species_tree = Tree(open(species_tree_path).read(), format=1)
    true_mul_tree = Tree(open(ground_truth_path).read(), format=1)

    # Run inference
    pred_mul_tree, events, raw = predict_mul_tree(
        model, gene_trees, species_tree, device=device,
        threshold=threshold, max_events=max_events,
    )

    # Evaluate
    metrics = evaluate_prediction(pred_mul_tree, true_mul_tree, species_tree)
    metrics["network"] = network_name
    metrics["n_gene_trees"] = len(gene_trees)
    metrics["n_events_predicted"] = len(events)

    return metrics, pred_mul_tree, events


def run_mock_evaluation(networks_dir):
    """Run evaluation with a randomly initialized model on mock data.

    Useful for testing the pipeline end-to-end.
    """
    print("=== Mock Evaluation (untrained model) ===\n")

    # Find network files
    network_files = [f for f in os.listdir(networks_dir) if f.endswith(".tre")]
    if not network_files:
        print(f"No .tre files found in {networks_dir}")
        return

    # Use smallest network
    network_file = sorted(network_files)[0]
    network_name = network_file.replace(".tre", "")

    print(f"Testing on: {network_name}")
    true_mul_tree = Tree(open(os.path.join(networks_dir, network_file)).read(), format=1)

    # Extract backbone as species tree
    species_tree = extract_backbone(true_mul_tree)
    species = sorted(set(species_tree.get_leaf_names()))
    n_species = len(species)

    print(f"  Species: {n_species}")
    print(f"  True polyploids: {get_polyploid_species(true_mul_tree)}")

    # Create mock gene trees (just copies of the MUL-tree)
    gene_trees = [true_mul_tree.copy("deepcopy") for _ in range(10)]

    # Create model (random weights)
    model = Gene2NetModel(
        n_species=n_species,
        embed_dim=16,
        n_attention_heads=2,
        n_gat_layers=1,
        handcrafted_dim=13,  # 8 copy count + 5 clustering summary (fixed)
    )

    # Run inference
    pred_mul_tree, events, raw = predict_mul_tree(
        model, gene_trees, species_tree,
        threshold=0.3,  # lower threshold for untrained model
    )

    # Evaluate
    metrics = evaluate_prediction(pred_mul_tree, true_mul_tree, species_tree)

    print(f"\n  Results (untrained model - expect poor):")
    print(f"    Ploidy accuracy: {metrics['ploidy_accuracy']:.3f}")
    print(f"    Polyploid F1: {metrics['polyploid_f1']:.3f}")
    print(f"    Predicted polyploids: {metrics['pred_n_polyploids']}")
    print(f"    True polyploids: {metrics['true_n_polyploids']}")
    print(f"    Events predicted: {len(events)}")
    print(f"\n  Pipeline test: PASSED (inference completed without errors)")


def main():
    parser = argparse.ArgumentParser(description="Evaluate Gene2Net-GNN")
    parser.add_argument("--model", help="Path to trained model (.pt)")
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--network", help="Single network name to evaluate")
    parser.add_argument("--gene-trees", help="Gene trees file for single network")
    parser.add_argument("--species-tree", help="Species tree file for single network")
    parser.add_argument("--networks-dir", default="../simulations/networks/", help="Ground truth networks directory")
    parser.add_argument("--data-dir", help="Directory with gene trees/species trees per network")
    parser.add_argument("--output", help="Output CSV file for results")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--mock", action="store_true", help="Run mock evaluation with random model")
    args = parser.parse_args()

    if args.mock:
        run_mock_evaluation(args.networks_dir)
        return

    if not args.model:
        parser.error("--model is required (or use --mock for testing)")

    device = torch.device(args.device)

    # Single network evaluation
    if args.network and args.gene_trees and args.species_tree:
        ground_truth = os.path.join(args.networks_dir, f"{args.network}.tre")

        import yaml
        with open(args.config) as f:
            config = yaml.safe_load(f)

        # Load model
        data_config = config.get("data", {})
        max_species = data_config.get("max_species", 80)

        model_config = config.get("model", {})
        model = Gene2NetModel(
            n_species=max_species,
            embed_dim=model_config.get("embed_dim", 64),
            n_attention_heads=model_config.get("num_attention_heads", 8),
            n_gat_layers=model_config.get("num_gat_layers", 3),
            handcrafted_dim=13,  # 8 copy count + 5 clustering summary (fixed)
            conv_type=model_config.get("conv_type", "gat"),
        )
        model.load_state_dict(torch.load(args.model, map_location=device, weights_only=True))
        model.to(device)

        metrics, pred_tree, events = evaluate_single_network(
            model, args.network, args.gene_trees, args.species_tree,
            ground_truth, device, args.threshold, args.max_events,
        )

        print(f"\n=== {args.network} ===")
        for k, v in metrics.items():
            if isinstance(v, float):
                print(f"  {k}: {v:.3f}")
            else:
                print(f"  {k}: {v}")
    else:
        print("Please provide either:")
        print("  --mock (for pipeline testing)")
        print("  --network NAME --gene-trees FILE --species-tree FILE --model FILE")
        print("  --data-dir DIR --model FILE (for batch evaluation)")


if __name__ == "__main__":
    main()
