#!/usr/bin/env python
"""Generate one training example using SimPhy for both species and gene trees.

Pipeline:
    1. SimPhy generates a species tree under birth-death
    2. Python adds random WGD events → MUL-tree
    3. SimPhy generates gene trees on the MUL-tree (with ILS)
    4. ASTRAL 4 infers species tree from gene trees
    5. Package as Gene2NetSample (features + labels)

Usage:
    python generate_one_example.py --output-dir /path/to/example_0001 \
        --n-species 20 --n-events 3 --ils-level medium --seed 42

    # For local testing without SimPhy/ASTRAL:
    python generate_one_example.py --output-dir /path/to/example_0001 \
        --n-species 10 --n-events 1 --seed 42 --mock
"""
import argparse
import glob
import json
import os
import random
import subprocess
import sys
import tempfile

from ete3 import Tree

# Add parent dir to path so gene2net_gnn is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.mul_tree_generator import (
    add_wgd_event,
    get_polyploid_species,
)
from gene2net_gnn.data.label_extractor import extract_backbone
from gene2net_gnn.data.tree_io import load_gene_trees_from_file
from gene2net_gnn.data.dataset import Gene2NetSample


# ---------------------------------------------------------------------------
# SimPhy helpers
# ---------------------------------------------------------------------------

def run_simphy_species_tree(n_species, output_dir, seed=42):
    """Run SimPhy to generate a species tree under birth-death.

    Returns the path to the species tree Newick file.
    """
    cmd = [
        "simphy",
        "-rs", "1",                    # 1 replicate
        "-rl", "f:1",                  # 1 locus (minimum, we only want the species tree)
        "-rg", "1",                    # 1 gene tree (minimum)
        "-sb", "f:0.000001",           # speciation rate (events/time unit)
        "-sd", "f:0.0000005",          # extinction rate
        "-sl", f"f:{n_species}",       # number of taxa
        "-sp", "f:10000",              # effective population size
        "-su", "f:0.00001",            # substitution rate
        "-sg", "f:1",                  # generation time
        "-cs", str(seed),              # random seed
        "-v", "0",                     # minimal verbosity
        "-o", output_dir,
    ]
    subprocess.run(cmd, check=True)
    return os.path.join(output_dir, "1", "s_tree.trees")


def suffixed_newick(mul_tree):
    """Convert MUL-tree to Newick with unique leaf names (sp_0_1, sp_0_2 for duplicates).

    SimPhy requires unique leaf labels, but MUL-trees have repeated species names.
    Returns the Newick string.
    """
    tree = mul_tree.copy("deepcopy")
    name_counts = {}
    for leaf in tree.get_leaves():
        name_counts[leaf.name] = name_counts.get(leaf.name, 0) + 1

    seen = {}
    for leaf in tree.get_leaves():
        name = leaf.name
        if name_counts[name] > 1:
            seen[name] = seen.get(name, 0) + 1
            leaf.name = f"{name}_{seen[name]}"

    return tree.write(format=5)


def run_simphy_gene_trees(mul_tree_newick, output_dir, n_gene_trees=1000,
                          ils_level="medium", seed=42):
    """Run SimPhy to generate gene trees on a MUL-tree via coalescent simulation.

    Uses -s to pass a fixed species tree (the MUL-tree with suffixed leaves).
    Uses -rl for number of locus trees (= independent gene families).
    """
    ne_map = {"low": 100000, "medium": 500000, "high": 1000000}
    ne = ne_map.get(ils_level, 500000)

    cmd = [
        "simphy",
        "-rs", "1",                         # 1 replicate
        "-rl", f"f:{n_gene_trees}",          # N locus trees = N independent gene families
        "-rg", "1",                          # 1 gene tree per locus
        "-s", mul_tree_newick,               # fixed species tree (MUL-tree)
        "-sp", f"f:{ne}",                    # effective population size (controls ILS)
        "-su", "f:0.00001",                  # substitution rate
        "-sg", "f:1",                        # generation time
        "-cs", str(seed),                    # random seed
        "-v", "0",                           # minimal verbosity
        "-o", output_dir,
    ]
    subprocess.run(cmd, check=True)


def collect_gene_trees(simphy_output_dir, output_path):
    """Collect all gene tree files from SimPhy output into one file, stripping suffixes.

    SimPhy with -rl N creates N gene tree files: g_trees1.trees, g_trees2.trees, etc.
    Each contains one gene tree. We combine them and strip _N suffixes from leaf names.
    """
    replicate_dir = os.path.join(simphy_output_dir, "1")
    gene_tree_files = sorted(glob.glob(os.path.join(replicate_dir, "g_trees*.trees")))

    all_trees = []
    for gt_file in gene_tree_files:
        trees = load_gene_trees_from_file(gt_file)
        all_trees.extend(trees)

    # Strip _N suffixes (sp_0_1 → sp_0)
    with open(output_path, "w") as f:
        for tree in all_trees:
            for leaf in tree.get_leaves():
                parts = leaf.name.rsplit("_", 1)
                if len(parts) == 2 and parts[1].isdigit():
                    leaf.name = parts[0]
            f.write(tree.write(format=5) + "\n")

    return all_trees


def run_astral(gene_trees_path, output_path):
    """Run ASTRAL 4 on multi-copy gene trees to infer species tree."""
    cmd = [
        "astral",             # adjust to your ASTRAL 4 binary name/path
        "-i", gene_trees_path,
        "-o", output_path,
    ]
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# Random WGD event generation
# ---------------------------------------------------------------------------

def add_random_wgd_events(species_tree, n_events, seed=42):
    """Add n_events random WGD events to a species tree.

    Returns (mul_tree, events_list).
    """
    rng = random.Random(seed)
    tree = species_tree.copy("deepcopy")
    events = []

    for _ in range(n_events):
        non_root_nodes = [n for n in tree.traverse("postorder") if n.up is not None]
        if not non_root_nodes:
            break

        target_node = rng.choice(non_root_nodes)
        target_clade = set(target_node.get_leaf_names())

        event_type = rng.choice(["auto", "allo"])
        partner_clade = None

        if event_type == "allo":
            candidate_partners = [
                n for n in non_root_nodes
                if n is not target_node and not (set(n.get_leaf_names()) & target_clade)
            ]
            if not candidate_partners:
                event_type = "auto"
            else:
                partner_node = rng.choice(candidate_partners)
                partner_clade = set(partner_node.get_leaf_names())

        tree = add_wgd_event(
            tree,
            target_clade=target_clade,
            partner_clade=partner_clade,
            event_type=event_type,
        )

        events.append({
            "target_clade": target_clade,
            "partner_clade": partner_clade,
            "event_type": event_type,
        })

    return tree, events


# ---------------------------------------------------------------------------
# Mock mode (local testing without SimPhy/ASTRAL)
# ---------------------------------------------------------------------------

def generate_mock_gene_trees(mul_tree, n_trees=100, seed=42):
    """Generate mock gene trees for local testing (no SimPhy needed).

    Copies the MUL-tree topology with random gene loss.
    """
    rng = random.Random(seed)
    trees = []
    for _ in range(n_trees):
        t = mul_tree.copy("deepcopy")
        polyploids = get_polyploid_species(t)
        for sp, count in polyploids.items():
            if count > 1 and rng.random() < 0.3:
                leaves = [l for l in t.get_leaves() if l.name == sp]
                if len(leaves) > 1:
                    leaf_to_remove = rng.choice(leaves)
                    parent = leaf_to_remove.up
                    leaf_to_remove.detach()
                    if parent and len(parent.children) == 1:
                        child = parent.children[0]
                        if parent.up:
                            grandparent = parent.up
                            parent.detach()
                            grandparent.add_child(child)
        trees.append(t)
    return trees


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate one training example")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-species", type=int, default=20)
    parser.add_argument("--n-events", type=int, default=2)
    parser.add_argument("--ils-level", choices=["low", "medium", "high"], default="medium")
    parser.add_argument("--n-gene-trees", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--mock", action="store_true", help="Use mock gene trees (no SimPhy/ASTRAL)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.mock:
        # ---- Mock mode (local testing) ----
        from gene2net_gnn.data.mul_tree_generator import generate_random_mul_tree
        print(f"Mock mode: generating MUL-tree with ETE3 ({args.n_species} species, {args.n_events} events)...")
        mul_tree, events = generate_random_mul_tree(args.n_species, args.n_events, seed=args.seed)
        backbone = extract_backbone(mul_tree)
        gene_trees = generate_mock_gene_trees(mul_tree, n_trees=args.n_gene_trees, seed=args.seed)
        species_tree = backbone.copy("deepcopy")
    else:
        # ---- Real mode: SimPhy + ASTRAL 4 ----
        with tempfile.TemporaryDirectory() as tmpdir:
            # Step 1: SimPhy generates species tree under birth-death
            print(f"Step 1: SimPhy generating species tree ({args.n_species} taxa, birth-death)...")
            simphy_sp_dir = os.path.join(tmpdir, "simphy_sp")
            sp_tree_path = run_simphy_species_tree(args.n_species, simphy_sp_dir, seed=args.seed)
            species_tree_bd = Tree(open(sp_tree_path).read().strip(), format=1)
            print(f"  Generated species tree with {len(species_tree_bd.get_leaves())} leaves")

            # Step 2: Add random WGD events → MUL-tree
            print(f"Step 2: Adding {args.n_events} WGD events...")
            mul_tree, events = add_random_wgd_events(species_tree_bd, args.n_events, seed=args.seed)
            backbone = extract_backbone(mul_tree)
            print(f"  Polyploids: {get_polyploid_species(mul_tree)}")

            # Step 3: SimPhy generates gene trees on MUL-tree
            print(f"Step 3: SimPhy generating {args.n_gene_trees} gene trees (ILS={args.ils_level})...")
            mul_newick = suffixed_newick(mul_tree)
            simphy_gt_dir = os.path.join(tmpdir, "simphy_gt")
            run_simphy_gene_trees(mul_newick, simphy_gt_dir, args.n_gene_trees, args.ils_level, args.seed + 1)

            # Collect and strip suffixes
            stripped_path = os.path.join(tmpdir, "gene_trees.nwk")
            collect_gene_trees(simphy_gt_dir, stripped_path)
            gene_trees = load_gene_trees_from_file(stripped_path)
            print(f"  Collected {len(gene_trees)} gene trees")

            # Step 4: ASTRAL 4 infers species tree from gene trees
            print("Step 4: Running ASTRAL 4...")
            astral_out = os.path.join(tmpdir, "astral.nwk")
            run_astral(stripped_path, astral_out)
            species_tree = Tree(open(astral_out).read().strip(), format=1)

    # Save MUL-tree
    with open(os.path.join(args.output_dir, "mul_tree.nwk"), "w") as f:
        f.write(mul_tree.write(format=5))

    # Step 5: Build training sample (features + labels)
    print("Step 5: Building training sample (features + labels)...")
    species_list = sorted(set(backbone.get_leaf_names()))
    sample = Gene2NetSample.from_trees(
        species_tree, gene_trees, species_list, mul_tree=mul_tree
    )
    sample.save(args.output_dir)

    # Save metadata
    metadata = {
        "n_species": args.n_species,
        "n_events": args.n_events,
        "ils_level": args.ils_level,
        "n_gene_trees": len(gene_trees),
        "seed": args.seed,
        "events": [
            {"target_clade": list(e.get("target_clade", set())),
             "partner_clade": list(e["partner_clade"]) if e.get("partner_clade") else [],
             "event_type": e.get("event_type", "unknown")}
            for e in events
        ],
        "polyploid_species": dict(get_polyploid_species(mul_tree)),
    }
    with open(os.path.join(args.output_dir, "metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"\nDone! Sample saved to {args.output_dir}")
    print(f"  Species: {len(species_list)}, Events: {args.n_events}")
    print(f"  Polyploids: {get_polyploid_species(mul_tree)}")
    print(f"  Gene trees: {len(gene_trees)}")


if __name__ == "__main__":
    main()
