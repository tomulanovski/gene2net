"""Characterize the KIND of ASTRAL backbone error at polyploids: how FAR is each
misplaced, not just whether it is misplaced.

Companion to backbone_polyploid_placement_accuracy.py (which only asks hit/miss:
is ASTRAL next to ONE true parent or neither). Here we measure the *topological
displacement* between where ASTRAL attaches each polyploid and where it belongs on
the true tree, then classify it. The thesis needs this to know whether a bounded,
local backbone repair (NNI-scale re-attachment) can fix the misplacements or
whether they are deep, long-range errors that need a full rebuild.

---------------------------------------------------------------------------
DISTANCE METRIC (chosen, well-defined — read this before trusting the numbers)
---------------------------------------------------------------------------
Coordinate system: the DIPLOID SKELETON. We build it by pruning the true MUL-tree
(mul_tree_{idx}.nex) to the single-copy (diploid) species only. Diploids appear
exactly once, so this backbone is unambiguous, and prior diagnostics show ASTRAL
reconstructs this diploid skeleton reliably (~86% exact). It is the stable ruler
against which we measure where a polyploid sits.

Attachment node of a polyploid copy: walk UP from that copy's leaf until the
enclosing clade contains at least one diploid species, then take the MRCA of those
diploid species *in the skeleton*. This node is the polyploid's diploid-context
anchor — effectively its parent node projected onto the diploid backbone. Walking
up handles the case where a polyploid's immediate sibling is another polyploid
(no diploid neighbor yet).

  - ASTRAL attachment node  : anchor of X's single copy in astral_species.tre.
  - True attachment node(s) : anchor of EACH copy of X in the MUL-tree. An
    allopolyploid appears twice, giving two true positions (its two parental
    origins); either is a legitimate "home" (matching the hit/miss logic of
    backbone_polyploid_placement_accuracy.py, which takes the BEST of the two).

Displacement = min over the true attachment nodes of the topological path length
(number of edges, ete3 get_distance(..., topology_only=True)) between the ASTRAL
attachment node and that true attachment node, measured on the skeleton. We take
the MIN so a polyploid placed next to EITHER true parent scores as correct.

Classification:
    displacement == 0  -> CORRECT      (ASTRAL anchors to a true parent's context)
    displacement 1-2   -> LOCAL        (NNI-scale error, a bounded re-attach fixes it)
    displacement  > 2  -> LONG_RANGE   (deep error, local repair insufficient)

NOTE — this is a defensible ANCHOR-BASED PROXY, not an exact node-to-node identity
between two different trees (ASTRAL and MUL are distinct trees; there is no shared
labelling of internal nodes). We project both attachments onto a common diploid
skeleton and measure distance there. The MRCA anchor can under-report distance when
a polyploid's diploid context is polyphyletic in the skeleton (MRCA floats up
toward the root). Treat the LOCAL/LONG_RANGE split as a well-grounded estimate, and
review before quoting exact percentages. FLAGGED as a proxy.

Read-only on all inputs. Run in any env with ete3 + pandas.

Usage:
    python scripts/backbone_displacement.py --config ils_low --max-samples 300 \
        --out output/backbone_displacement_ils_low.csv
    python scripts/backbone_displacement.py --config ils_low ils_med ils_high \
        --out output/backbone_displacement.csv
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd
from ete3 import Tree

from gene2net_gnn.data.mul_tree_generator import get_polyploid_species


def load_tree_any(path):
    """Load a Newick or NEXUS tree (mirrors the reference backbone scripts)."""
    txt = open(path).read()
    if "#nexus" in txt.lower() or "begin trees" in txt.lower():
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(txt.strip(), format=1)


def prune_to(tree, keep_names):
    """Copy `tree` and prune to `keep_names`; None if fewer than 3 survive.

    Reused from backbone_polyploid_localization.py.
    """
    t = tree.copy()
    keep = [l for l in t.get_leaf_names() if l in keep_names]
    if len(keep) < 3:
        return None
    t.prune(keep)
    return t


def skeleton_anchor(skeleton, skel_leaves, copy_leaf):
    """Project a polyploid copy's attachment onto the diploid skeleton.

    Walk up from `copy_leaf` until the enclosing clade holds >=1 diploid species,
    then return the MRCA of those diploids *in the skeleton* (a single leaf node if
    there is only one). None if no diploid context is ever found.
    """
    node = copy_leaf
    while node.up is not None:
        node = node.up
        dips = set(node.get_leaf_names()) & skel_leaves
        if dips:
            dips = sorted(dips)
            if len(dips) == 1:
                hits = skeleton.search_nodes(name=dips[0])
                return hits[0] if hits else None
            try:
                return skeleton.get_common_ancestor(dips)
            except Exception:
                return None
    return None


def displacement_for(skeleton, skel_leaves, astral, mul, x):
    """Topological displacement (edges on the skeleton) of polyploid `x`.

    min over the true copies of the path length between ASTRAL's attachment node
    and each true attachment node. None if it cannot be measured.
    """
    astral_copies = [l for l in astral.get_leaves() if l.name == x]
    mul_copies = [l for l in mul.get_leaves() if l.name == x]
    if not astral_copies or not mul_copies:
        return None

    astral_anchor = skeleton_anchor(skeleton, skel_leaves, astral_copies[0])
    if astral_anchor is None:
        return None

    true_anchors = [skeleton_anchor(skeleton, skel_leaves, c) for c in mul_copies]
    true_anchors = [a for a in true_anchors if a is not None]
    if not true_anchors:
        return None

    dists = []
    for ta in true_anchors:
        try:
            dists.append(int(round(skeleton.get_distance(astral_anchor, ta,
                                                          topology_only=True))))
        except Exception:
            continue
    if not dists:
        return None
    return min(dists)


def classify(d):
    if d == 0:
        return "correct"
    if d <= 2:
        return "local"
    return "long_range"


def rows_for_config(config, data_root, replicate, max_samples):
    sim_dir = os.path.join(data_root, "simphy", config)
    if not os.path.isdir(sim_dir):
        print(f"[{config}] simphy dir not found: {sim_dir} — skipping")
        return []
    idxs = sorted(d for d in os.listdir(sim_dir)
                  if os.path.isdir(os.path.join(sim_dir, d)))[:max_samples]

    rows = []
    for idx in idxs:
        astral_path = os.path.join(sim_dir, idx, f"replicate_{replicate}",
                                   "astral_species.tre")
        mul_path = os.path.join(data_root, f"mul_tree_{idx}.nex")
        if not (os.path.exists(astral_path) and os.path.exists(mul_path)):
            continue

        astral = load_tree_any(astral_path)
        mul = load_tree_any(mul_path)
        polys = get_polyploid_species(mul)      # {species: count} for count > 1
        if not polys:
            continue

        diploids = set(mul.get_leaf_names()) - set(polys)
        skeleton = prune_to(mul, diploids)
        if skeleton is None:
            continue
        skel_leaves = set(skeleton.get_leaf_names())

        for x in polys:
            d = displacement_for(skeleton, skel_leaves, astral, mul, x)
            if d is None:
                continue
            rows.append({
                "config": config,
                "sample": idx,
                "species": x,
                "displacement": d,
                "class": classify(d),
            })
    return rows


def print_aggregate(df, label):
    n = len(df)
    if n == 0:
        print(f"  {label}: no polyploids measured")
        return
    counts = df["class"].value_counts()
    n_correct = int(counts.get("correct", 0))
    n_local = int(counts.get("local", 0))
    n_long = int(counts.get("long_range", 0))
    # of the MISplacements (displacement > 0), what fraction is local vs long-range
    n_mis = n_local + n_long
    mis_local = (100 * n_local / n_mis) if n_mis else 0.0
    mis_long = (100 * n_long / n_mis) if n_mis else 0.0
    print(f"  {label}: {n} polyploids")
    print(f"    correct     (d=0):  {n_correct:5d} ({100 * n_correct / n:5.1f}%)")
    print(f"    local     (d=1-2):  {n_local:5d} ({100 * n_local / n:5.1f}%)")
    print(f"    long_range (d>2):   {n_long:5d} ({100 * n_long / n:5.1f}%)")
    print(f"    median displacement: {df['displacement'].median():.1f}   "
          f"mean {df['displacement'].mean():.2f}")
    if n_mis:
        print(f"    of misplacements: {mis_local:.1f}% local (repairable), "
              f"{mis_long:.1f}% long-range")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", nargs="+", default=["ils_low"],
                    help="one or more ILS configs")
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-samples", type=int, default=300)
    ap.add_argument("--out", default="output/backbone_displacement.csv",
                    help="per-polyploid CSV output path")
    args = ap.parse_args()

    all_rows = []
    for config in args.config:
        all_rows.extend(rows_for_config(config, args.data_root,
                                        args.replicate, args.max_samples))

    df = pd.DataFrame(all_rows,
                      columns=["config", "sample", "species", "displacement", "class"])

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    df.to_csv(args.out, index=False)

    print(f"\nWrote {len(df)} per-polyploid rows -> {args.out}\n")
    if df.empty:
        print("No polyploids measured — check --data-root / --config / paths.")
        return

    print("Displacement of ASTRAL polyploid placements "
          "(edges on the diploid skeleton):\n")
    print_aggregate(df, "OVERALL")
    if len(args.config) > 1 or df["config"].nunique() > 1:
        print("\nPer ILS config:")
        for config in sorted(df["config"].unique()):
            print_aggregate(df[df["config"] == config], config)

    print("\nReading:")
    print("  mostly correct/local -> polyploid misplacements are NNI-scale ->")
    print("    a bounded/local backbone repair can fix them.")
    print("  lots of long_range   -> misplacements are deep -> local repair is")
    print("    insufficient; the backbone needs a rebuild (phasing / joint).")


if __name__ == "__main__":
    main()
