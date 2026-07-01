"""Root an (unrooted) ASTRAL species tree for MUL-tree reconstruction.

ASTRAL outputs an unrooted species tree; the arbitrary root it gets when read
inflates the MUL-tree edit distance badly (only ~4% of arbitrary roots match the
true root). We root it with a hybrid of two signals, validated against the true
root on simulated data (scripts/root_from_gene_trees.py):

    rooting method          ils_low   dup_loss_high   ils_high
    arbitrary ASTRAL          4%          6%            2%
    midpoint                 72%         64%           34%
    gene-tree consensus      48%*        62%*          44%*   (*~89% when it commits)
    HYBRID (this module)     82%         82%           64%

Hybrid = use the gene-tree consensus root when it gives a clean clade, else fall
back to midpoint. The gene trees carry the rooting signal where the clock breaks
down (high ILS), midpoint covers the rest.
"""
from collections import Counter
from typing import List

from ete3 import Tree


def robust_set_outgroup(tree: Tree, outgroup_leaves) -> bool:
    """Re-root `tree` in-place so `outgroup_leaves` is the outgroup clade.

    Robust to the common ete3 case where the target outgroup spans the tree's
    current (arbitrary) root — then get_common_ancestor returns the whole tree and
    set_outgroup throws. We detect that, re-root at a leaf OUTSIDE the outgroup so
    the outgroup becomes monophyletic, and retry.

    Returns True if the tree was rooted on the outgroup, False if the outgroup is
    genuinely not a clade in this (unrooted) topology.
    """
    leaves = set(tree.get_leaf_names())
    og = set(l for l in outgroup_leaves if l in leaves)
    if not og or len(og) >= len(leaves):
        return False
    try:
        if len(og) == 1:
            tree.set_outgroup(next(iter(og)))
            return True
        mrca = tree.get_common_ancestor(list(og))
        if set(mrca.get_leaf_names()) != og:
            # Outgroup spans the current root: re-root at a non-outgroup leaf first.
            anchor = next((l for l in leaves if l not in og), None)
            if anchor is None:
                return False
            tree.set_outgroup(anchor)
            mrca = tree.get_common_ancestor(list(og))
            if set(mrca.get_leaf_names()) != og:
                return False  # genuinely not monophyletic (unrooted topologies differ)
        tree.set_outgroup(mrca)
        return True
    except Exception:
        return False


def _root_smaller_side(tree: Tree):
    """Species set on the smaller side of a rooted tree's basal split."""
    kids = tree.children
    if len(kids) < 2:
        return None
    sides = [frozenset(c.get_leaf_names()) for c in kids]
    return min(sides, key=len)


def consensus_outgroup(gene_trees: List[Tree], max_gene_trees: int = 500):
    """Most common basal split (smaller side) across the rooted gene trees."""
    counter = Counter()
    for gt in gene_trees[:max_gene_trees]:
        b = _root_smaller_side(gt)
        if b:
            counter[b] += 1
    if not counter:
        return None
    return counter.most_common(1)[0][0]


def hybrid_root(species_tree: Tree, gene_trees: List[Tree],
                max_gene_trees: int = 500) -> Tree:
    """Return a rooted copy of `species_tree`.

    Tries the gene-tree consensus outgroup; if that set re-roots cleanly (it is
    monophyletic in the species tree) we use it, otherwise we fall back to
    midpoint rooting. Matches the validated 82%/82%/64% hybrid.
    """
    cand = consensus_outgroup(gene_trees, max_gene_trees)
    if cand:
        t = species_tree.copy()
        # robust_set_outgroup handles the case where the consensus outgroup spans
        # the arbitrary root (previously this silently fell back to midpoint ~42%
        # of the time, discarding the more-accurate consensus signal).
        if robust_set_outgroup(t, cand):
            return t

    # midpoint fallback (only when the consensus set is genuinely not a clade)
    t = species_tree.copy()
    try:
        t.set_outgroup(t.get_midpoint_outgroup())
    except Exception:
        pass
    return t
