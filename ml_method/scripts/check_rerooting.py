"""Does reroot_to_match actually re-root ASTRAL to the true root on RF=0 samples?

The astral_rerooted build (oracle rooting to the true root) scored 0.556 on RF=0
samples, between astral (0.699) and true (0.192). That blend suggests
reroot_to_match succeeds on some samples and falls back to the arbitrary root on
others. On RF=0 samples (ASTRAL unrooted topology == true), a correct re-root
should reproduce the true tree exactly (rooted RF = 0). This measures how often it
actually does.

ete3 only (no torch). Run in gene2net.

Usage:
    python scripts/check_rerooting.py --config ils_low --max-samples 300
"""
import argparse
import os

from ete3 import Tree


def load_astral(path):
    return Tree(open(path).read().strip(), format=1)


def load_nexus(path):
    for line in open(path).read().split("\n"):
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(open(path).read().strip(), format=1)


def reroot_to_match(astral_tree, true_tree):
    """Copy of reconstruct_mul_tree.reroot_to_match (inlined to avoid torch import)."""
    kids = true_tree.children
    if len(kids) < 2:
        return astral_tree, "true_no_root"
    sides = [set(c.get_leaf_names()) for c in kids]
    outgroup = min(sides, key=len)
    try:
        og = list(outgroup)
        if len(og) == 1:
            astral_tree.set_outgroup(og[0])
        else:
            astral_tree.set_outgroup(astral_tree.get_common_ancestor(og))
        return astral_tree, "ok"
    except Exception:
        return astral_tree, "fallback"


def norm_rf(t1, t2, unrooted):
    try:
        v = float(t1.compare(t2, unrooted=unrooted)["norm_rf"])
    except Exception:
        return None
    return None if v != v else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-samples", type=int, default=300)
    args = ap.parse_args()

    sim_dir = os.path.join(args.data_root, "simphy", args.config)
    idxs = sorted(d for d in os.listdir(sim_dir)
                  if os.path.isdir(os.path.join(sim_dir, d)))[:args.max_samples]

    n_rf0 = 0
    n_reroot_ok = n_fallback = 0
    n_rooted_match = 0   # rerooted rooted-topology == true (the thing that matters)
    for idx in idxs:
        ap_path = os.path.join(sim_dir, idx, f"replicate_{args.replicate}", "astral_species.tre")
        tp_path = os.path.join(args.data_root, f"species_tree_{idx}.nex")
        if not (os.path.exists(ap_path) and os.path.exists(tp_path)):
            continue
        astral = load_astral(ap_path)
        true = load_nexus(tp_path)
        if norm_rf(astral, true, unrooted=True) != 0:
            continue  # only RF=0 (exact unrooted topology)
        n_rf0 += 1
        rr, status = reroot_to_match(astral.copy(), true)
        if status == "ok":
            n_reroot_ok += 1
        elif status == "fallback":
            n_fallback += 1
        rooted = norm_rf(rr, true, unrooted=False)
        if rooted == 0:
            n_rooted_match += 1

    if n_rf0 == 0:
        print("No RF=0 samples found.")
        return

    print(f"\nReroot check on {n_rf0} RF=0 samples (config {args.config})\n")
    print(f"  reroot_to_match status: ok={n_reroot_ok}, fallback(exception)={n_fallback}")
    print(f"  rerooted rooted-topology == true: {n_rooted_match}/{n_rf0} "
          f"({100 * n_rooted_match / n_rf0:.0f}%)")
    print("\nReading:")
    print("  high match% -> rerooting works; the 0.556 build score is then a BUILD bug")
    print("    (events mis-placed on the rerooted tree), not a rooting failure.")
    print("  low match%  -> reroot_to_match is fragile; fixing it pushes RF=0 toward true (0.192).")


if __name__ == "__main__":
    main()
