"""Validate rooting the ASTRAL species tree from the gene-tree root consensus.

ASTRAL outputs an unrooted tree; its arbitrary root is what inflates the MUL-tree
edit distance. But the gene trees ARE rooted, so we can recover a root for the
species tree from where the gene trees consistently put theirs.

Method per sample:
  - For each rooted gene tree, take the species set on the smaller side of its
    root (its basal split).
  - The most common basal split across gene trees = the consensus outgroup.
  - Re-root the ASTRAL tree on that outgroup.

We then check, against the known true backbone, how often each rooting recovers
the true root:
  - original ASTRAL (arbitrary)
  - consensus-from-gene-trees
  - midpoint (baseline; no gene-tree info)

If consensus rooting recovers the true root well, it is worth wiring into the
pipeline (root species tree -> recompute features -> re-run reconstruction).

Run in the gene2net env (ete3). No model needed.

Usage:
    python scripts/root_from_gene_trees.py \
        --mul-trees-dir data/mul_trees_2k --config ils_low --n 50
"""
import argparse
import os
from collections import Counter

from ete3 import Tree


def load_nexus_tree(path):
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree in {path}")


def root_smaller_side(tree):
    """Species set on the smaller side of the tree's root (its basal split)."""
    kids = tree.children
    if len(kids) < 2:
        return None
    sides = [frozenset(c.get_leaf_names()) for c in kids]
    # if >2 children (unrooted trifurcation), take the smallest child
    return min(sides, key=len)


def reroot_on(tree, outgroup_species):
    """Re-root a copy of `tree` on the given species set; None if not possible."""
    t = tree.copy()
    present = [s for s in outgroup_species if s in set(t.get_leaf_names())]
    if not present:
        return None
    try:
        if len(present) == 1:
            t.set_outgroup(present[0])
        else:
            t.set_outgroup(t.get_common_ancestor(present))
        return t
    except Exception:
        return None


def root_bip(tree):
    """Canonical root bipartition (smaller side species set) for comparison."""
    s = root_smaller_side(tree)
    return s


def freq_reroot(astral, basal_freq, n_gt):
    """Root ASTRAL using per-species basal frequency from the gene trees.

    basal_freq[s] = how often species s sat on the basal (smaller-root) side
    across gene trees. We root at the ASTRAL edge whose clade is most enriched
    for high-basal-frequency species, i.e. the most 'outgroup-like' clade. This
    uses the full distribution of gene-tree roots and always picks a real edge."""
    species = set(astral.get_leaf_names())
    freq = {s: basal_freq.get(s, 0) / max(n_gt, 1) for s in species}
    overall = sum(freq.values()) / max(len(freq), 1)
    best_node, best_score = None, -1e9
    for node in astral.traverse():
        if node.is_root() or node.is_leaf():
            continue
        clade = node.get_leaf_names()
        if len(clade) >= len(species):
            continue
        inside = sum(freq[s] for s in clade) / len(clade)
        # prefer smaller, strongly-basal clades
        score = inside - overall
        if score > best_score:
            best_score, best_node = score, node
    if best_node is None:
        return None
    t = astral.copy()
    try:
        t.set_outgroup(t.get_common_ancestor(best_node.get_leaf_names())
                       if len(best_node.get_leaf_names()) > 1
                       else best_node.get_leaf_names()[0])
        return t
    except Exception:
        return None


def matches_true(tree, true_bip, all_species):
    """Does this tree's root bipartition equal the true root bipartition?"""
    b = root_bip(tree)
    if b is None or true_bip is None:
        return False
    # a bipartition is the same whether named by side A or its complement
    return b == true_bip or b == frozenset(all_species - true_bip)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--max-gene-trees", type=int, default=500)
    args = parser.parse_args()

    orig_ok = cons_ok = freq_ok = mid_ok = hybrid_ok = total = 0
    cons_failed = freq_failed = 0

    for idx in range(args.start, args.start + args.n):
        idx_str = f"{idx:04d}"
        rep = os.path.join(args.mul_trees_dir, "simphy", args.config, idx_str,
                           f"replicate_{args.replicate}")
        gt_path = os.path.join(rep, "gene_trees.tre")
        astral_path = os.path.join(rep, "astral_species.tre")
        true_path = os.path.join(args.mul_trees_dir, f"species_tree_{idx_str}.nex")
        if not (os.path.exists(gt_path) and os.path.exists(astral_path)
                and os.path.exists(true_path)):
            continue

        try:
            astral = Tree(open(astral_path).read().strip(), format=1)
            true = load_nexus_tree(true_path)
        except Exception:
            continue
        all_species = set(astral.get_leaf_names())
        true_bip = root_bip(true)

        # consensus basal split across gene trees + per-species basal frequency
        counter = Counter()
        basal_freq = Counter()
        n_read = 0
        for line in open(gt_path):
            line = line.strip()
            if not line:
                continue
            try:
                gt = Tree(line, format=1)
            except Exception:
                continue
            b = root_smaller_side(gt)
            if b:
                counter[b] += 1
                for s in b:
                    basal_freq[s] += 1
            n_read += 1
            if n_read >= args.max_gene_trees:
                break
        if not counter:
            continue

        consensus = counter.most_common(1)[0][0]
        rerooted = reroot_on(astral, consensus)
        freq_rooted = freq_reroot(astral, basal_freq, n_read)
        midpoint = astral.copy()
        try:
            midpoint.set_outgroup(midpoint.get_midpoint_outgroup())
        except Exception:
            midpoint = None

        total += 1
        if matches_true(astral, true_bip, all_species):
            orig_ok += 1
        if rerooted is not None:
            if matches_true(rerooted, true_bip, all_species):
                cons_ok += 1
        else:
            cons_failed += 1
        if freq_rooted is not None:
            if matches_true(freq_rooted, true_bip, all_species):
                freq_ok += 1
        else:
            freq_failed += 1
        if midpoint is not None and matches_true(midpoint, true_bip, all_species):
            mid_ok += 1
        # hybrid: gene-tree consensus when it produced a clean reroot, else midpoint
        hybrid = rerooted if rerooted is not None else midpoint
        if hybrid is not None and matches_true(hybrid, true_bip, all_species):
            hybrid_ok += 1

    print(f"\nRoot-recovery on {total} samples ({args.config}):")
    print(f"  original ASTRAL (arbitrary) recovers true root: {orig_ok}/{total} "
          f"({100*orig_ok/max(total,1):.0f}%)")
    print(f"  gene-tree CONSENSUS (mode split) recovers root: {cons_ok}/{total} "
          f"({100*cons_ok/max(total,1):.0f}%)   [reroot failed on {cons_failed}]")
    print(f"  gene-tree FREQUENCY rooting recovers true root: {freq_ok}/{total} "
          f"({100*freq_ok/max(total,1):.0f}%)   [reroot failed on {freq_failed}]")
    print(f"  MIDPOINT rooting recovers true root:            {mid_ok}/{total} "
          f"({100*mid_ok/max(total,1):.0f}%)")
    print(f"  HYBRID (consensus-if-clean else midpoint):      {hybrid_ok}/{total} "
          f"({100*hybrid_ok/max(total,1):.0f}%)")
    print("\nIf consensus (or midpoint) >> original, that rooting method is the fix:")
    print("wire it in, recompute features on the rooted tree, re-run reconstruction.")


if __name__ == "__main__":
    main()
