#!/usr/bin/env python
"""Generate random MUL-trees for training data.

Uses SimPhy to generate ultrametric species trees via birth-death, adds random
WGD events, and saves in NEXUS format so the existing SimPhy gene tree pipeline
can be reused.

Usage:
    # Generate one MUL-tree (called by SLURM array job)
    python generate_mul_trees.py --index 42 --output-dir /path/to/mul_trees

    # Generate a batch locally for testing
    python generate_mul_trees.py --index 0 --n-batch 10 --output-dir /path/to/mul_trees

    # Local testing without SimPhy (uses Python birth-death fallback)
    python generate_mul_trees.py --index 0 --n-batch 5 --output-dir /tmp/test --mock
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

from ete3 import Tree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.mul_tree_generator import (
    generate_birth_death_tree,
    add_wgd_event,
    get_polyploid_species,
    _subtree_height,
    _force_ultrametric,
)
from gene2net_gnn.data.label_extractor import extract_backbone

# SimPhy binary on the cluster
SIMPHY_BIN = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulators/simphy/SimPhy_1.0.2/bin/simphy_lnx64"


def run_simphy_species_tree(n_species, output_dir, tree_height=10_000_000.0, seed=42):
    """Run SimPhy to generate a species tree under birth-death.

    Returns an ETE3 Tree object (ultrametric, from SimPhy's birth-death model).
    """
    # SimPhy birth/death rates tuned for reasonable tree sizes
    # Same as generate_one_example.py reference implementation
    cmd = [
        SIMPHY_BIN,
        "-rs", "1",                    # 1 replicate
        "-rl", "f:1",                  # 1 locus (minimum, we only want the species tree)
        "-rg", "1",                    # 1 gene tree (minimum)
        "-sb", "f:0.000001",           # speciation rate (events/generation)
        "-sd", "f:0.0000005",          # extinction rate
        "-sl", f"f:{n_species}",       # number of taxa
        "-sp", "f:10000",              # effective population size
        "-su", "f:0.00001",            # substitution rate
        "-sg", "f:1",                  # generation time
        "-cs", str(seed),              # random seed
        "-v", "0",                     # minimal verbosity
        "-o", output_dir,
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    # SimPhy writes species tree to <output_dir>/1/s_tree.trees
    sp_tree_path = os.path.join(output_dir, "1", "s_tree.trees")
    with open(sp_tree_path) as f:
        newick = f.read().strip()
    return Tree(newick, format=1)


def generate_one_mul_tree(index, tree_height=10_000_000.0, mock=False):
    """Generate one random MUL-tree with deterministic parameters based on index.

    Uses SimPhy to generate the species tree (birth-death model), then adds
    random WGD events in Python.

    Parameter ranges (designed to cover the diversity seen in real networks):
        n_species: 5-80 (varying)
        n_events:  0-15 (varying, including 0 for negative examples)
        event_type: mix of auto and allo

    Args:
        mock: If True, use Python birth-death fallback (for local testing without SimPhy).
    """
    # Deterministic parameter assignment from index
    rng = random.Random(index)

    # Species count: weighted toward 10-40 (most common in real networks)
    n_species = rng.choice(
        [5, 5, 10, 10, 10, 15, 15, 15, 20, 20, 20, 20,
         25, 25, 30, 30, 35, 40, 45, 50, 60, 80]
    )

    # Number of WGD events: weighted toward low counts (realistic)
    # ~10% chance of 0 events (negative examples)
    n_events = rng.choices(
        population=[0, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 8, 10, 12, 15],
        weights=[3, 5, 5, 5, 5, 5, 4, 4, 3, 3, 2, 1, 1, 1, 1],
        k=1,
    )[0]

    # Cap events at n_species - 1 (can't have more events than edges)
    n_events = min(n_events, n_species - 1)

    seed = index * 7919 + 42  # deterministic seed per index

    # Generate ultrametric species tree
    if mock:
        # Python fallback for local testing
        species_tree = generate_birth_death_tree(
            n_species=n_species,
            birth_rate=1.0,
            death_rate=0.3,
            tree_height=tree_height,
            seed=seed,
        )
    else:
        # Use SimPhy birth-death (published simulator — Mallo et al. 2016)
        with tempfile.TemporaryDirectory() as tmpdir:
            simphy_dir = os.path.join(tmpdir, f"simphy_sp_{index}")
            print(f"  Running SimPhy for species tree (n={n_species}, seed={seed})...")
            species_tree = run_simphy_species_tree(
                n_species=n_species,
                output_dir=simphy_dir,
                tree_height=tree_height,
                seed=seed,
            )
            print(f"  SimPhy generated tree with {len(species_tree.get_leaves())} leaves")

    # Rename leaves to sp_0, sp_1, ... for consistency
    for i, leaf in enumerate(sorted(species_tree.get_leaves(), key=lambda l: l.name)):
        leaf.name = f"sp_{i}"

    # Add random WGD events
    mul_tree = species_tree.copy("deepcopy")
    events = []
    event_rng = random.Random(seed + 1)

    for _ in range(n_events):
        non_root_nodes = [n for n in mul_tree.traverse("postorder") if n.up is not None]
        if not non_root_nodes:
            break

        target_node = event_rng.choice(non_root_nodes)
        target_clade = set(target_node.get_leaf_names())

        event_type = event_rng.choice(["auto", "allo"])
        partner_clade = None

        if event_type == "allo":
            candidate_partners = [
                n for n in non_root_nodes
                if n is not target_node and not (set(n.get_leaf_names()) & target_clade)
            ]
            if not candidate_partners:
                event_type = "auto"
            else:
                partner_node = event_rng.choice(candidate_partners)
                partner_clade = set(partner_node.get_leaf_names())

        try:
            mul_tree = add_wgd_event(
                mul_tree,
                target_clade=target_clade,
                partner_clade=partner_clade,
                event_type=event_type,
            )
            events.append({
                "target_clade": sorted(target_clade),
                "partner_clade": sorted(partner_clade) if partner_clade else None,
                "event_type": event_type,
            })
        except ValueError:
            # Skip invalid events (e.g., clade not found after previous modifications)
            continue

    # Fix any floating-point drift from WGD additions
    actual_height = _subtree_height(mul_tree)
    _force_ultrametric(mul_tree, actual_height)

    return species_tree, mul_tree, events, {
        "n_species": n_species,
        "n_events_requested": n_events,
        "n_events_applied": len(events),
        "seed": seed,
    }


def _tree_to_nexus_newick(tree):
    """Convert ETE3 tree to Newick string with 12-decimal branch lengths.

    Same format as simulations/scripts/rescale_and_keep_ultrametric.py —
    this is what SimPhy expects via -SR flag.
    """
    def format_node(node):
        name = node.name if node.is_leaf() and node.name else ""
        if node.is_leaf():
            result = name
        else:
            children = ",".join(format_node(child) for child in node.children)
            result = f"({children})"
        if not node.is_root():
            result += f":{node.dist:.12f}"
        return result
    return format_node(tree) + ";"


def _write_nexus(tree, filepath, tree_name="1"):
    """Write tree in NEXUS format (what SimPhy reads via -SR)."""
    newick_str = _tree_to_nexus_newick(tree)
    with open(filepath, "w") as f:
        f.write("#NEXUS\n")
        f.write("begin trees;\n")
        f.write(f"tree {tree_name}={newick_str}\n")
        f.write("end;\n")


def save_mul_tree(output_dir, index, species_tree, mul_tree, events, params):
    """Save a MUL-tree in NEXUS format expected by SimPhy (-SR flag).

    Directory structure:
        output_dir/
            mul_tree_NNNN.nex          # MUL-tree in NEXUS (SimPhy input)
            species_tree_NNNN.nex      # Original species tree (backbone)
            metadata_NNNN.json         # Parameters and event details
    """
    name = f"{index:04d}"

    # Save MUL-tree in NEXUS format (SimPhy reads this via -SR)
    mul_path = os.path.join(output_dir, f"mul_tree_{name}.nex")
    _write_nexus(mul_tree, mul_path)

    # Save original species tree (backbone) in NEXUS too
    sp_path = os.path.join(output_dir, f"species_tree_{name}.nex")
    _write_nexus(species_tree, sp_path)

    # Save metadata
    polyploids = get_polyploid_species(mul_tree)
    meta = {
        **params,
        "index": index,
        "n_leaves_mul_tree": len(mul_tree.get_leaves()),
        "n_polyploid_species": len(polyploids),
        "polyploid_species": {k: v for k, v in sorted(polyploids.items())},
        "events": events,
        "tree_height": _subtree_height(mul_tree),
    }
    meta_path = os.path.join(output_dir, f"metadata_{name}.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return mul_path


def main():
    parser = argparse.ArgumentParser(description="Generate random MUL-trees for training")
    parser.add_argument("--index", type=int, required=True,
                        help="Starting index (or SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--n-batch", type=int, default=1,
                        help="Number of MUL-trees to generate (for local batch mode)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for MUL-tree files")
    parser.add_argument("--tree-height", type=float, default=10_000_000.0,
                        help="Root-to-tip distance in generations (default: 10M)")
    parser.add_argument("--mock", action="store_true",
                        help="Use Python birth-death fallback (local testing without SimPhy)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for i in range(args.n_batch):
        idx = args.index + i
        species_tree, mul_tree, events, params = generate_one_mul_tree(
            idx, tree_height=args.tree_height, mock=args.mock
        )
        path = save_mul_tree(args.output_dir, idx, species_tree, mul_tree, events, params)

        polyploids = get_polyploid_species(mul_tree)
        print(f"[{idx:04d}] species={params['n_species']}, "
              f"events={params['n_events_applied']}/{params['n_events_requested']}, "
              f"polyploids={len(polyploids)}, "
              f"leaves={len(mul_tree.get_leaves())} -> {path}")


if __name__ == "__main__":
    main()
