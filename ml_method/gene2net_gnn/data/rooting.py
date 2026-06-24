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
    species = set(species_tree.get_leaf_names())

    cand = consensus_outgroup(gene_trees, max_gene_trees)
    if cand:
        present = [s for s in cand if s in species]
        if present:
            t = species_tree.copy()
            try:
                if len(present) == 1:
                    t.set_outgroup(present[0])
                else:
                    t.set_outgroup(t.get_common_ancestor(present))
                return t
            except Exception:
                pass  # consensus set not monophyletic -> fall back

    # midpoint fallback
    t = species_tree.copy()
    try:
        t.set_outgroup(t.get_midpoint_outgroup())
    except Exception:
        pass
    return t
