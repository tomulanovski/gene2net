"""DEPRECATED / DO NOT TRUST. This standalone re-scorer builds trees OUTSIDE
ReticulateTree and gave results inconsistent with the pipeline's own scorer, which
already computes canonical (order-invariant) edit distance via
ReticulateTree.get_edit_distance_multree(canonical=True). Use score_reconstructions.py
for edit distance. Kept only as a record of the mistake.

Re-score a benchmark output directory's MUL-tree edit distance BOTH ways:
current (Newick-child-order-sensitive, as the pipeline does it) and canonical
(child order removed). Used to test whether an apparent edit-distance improvement
is a real reconstruction gain or an artifact of the greedy-GED child-order confound.

Reads output/<bench-dir>/<config>/<strategy>/<net>/output.tre and ground_truth.nex.
Self-contained copy of the canonical GED from simulations/check_edit_distance_ordering.py.

Run on two models and compare the CANONICAL columns:
  python scripts/rescore_canonical_edit.py --bench-dir output/bench_away --configs conf_ils_low_10M ...
  python scripts/rescore_canonical_edit.py --bench-dir output/bench_away_clade --configs conf_ils_low_10M ...
"""
import argparse
import os

import networkx as nx
from ete3 import Tree

MAX_NODES_FOR_GED = 500
DEFAULT_TIMEOUT = 90
NODE_MATCH = lambda u, v: u.get("label") == v.get("label")


def load_nexus_or_newick(path):
    txt = open(path).read()
    for line in txt.splitlines():
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(txt.strip(), format=1)


def build_graph(tree, canonical):
    G = nx.DiGraph()
    if not canonical:
        for node in tree.traverse():
            G.add_node(id(node), label=node.name if node.is_leaf() else None)
            if not node.is_root():
                G.add_edge(id(node.up), id(node))
        return G
    canon = {}

    def form(node):
        f = node.name if node.is_leaf() else "(" + ",".join(sorted(form(c) for c in node.children)) + ")"
        canon[node] = f
        return f

    form(tree)

    def add(node):
        G.add_node(id(node), label=node.name if node.is_leaf() else None)
        if not node.is_root():
            G.add_edge(id(node.up), id(node))
        for c in sorted(node.children, key=lambda c: canon[c]):
            add(c)

    add(tree)
    return G


def edit_distance(t1, t2, canonical, timeout=DEFAULT_TIMEOUT):
    g1, g2 = build_graph(t1, canonical), build_graph(t2, canonical)
    if max(len(g1.nodes), len(g2.nodes)) > MAX_NODES_FOR_GED:
        return float("nan")
    import signal
    use_alarm = hasattr(signal, "SIGALRM")
    if use_alarm:
        def _h(s, f):
            raise TimeoutError()
        old = signal.signal(signal.SIGALRM, _h)
        signal.alarm(timeout)
    try:
        d = next(nx.optimize_graph_edit_distance(g1, g2, node_match=NODE_MATCH))
    except Exception:
        return float("nan")
    finally:
        if use_alarm:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
    norm = max(len(g1.nodes) + len(g1.edges), len(g2.nodes) + len(g2.edges))
    return d / norm if norm else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-dir", required=True, help="e.g. output/bench_away")
    ap.add_argument("--strategy", default="cap")
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--max-per-config", type=int, default=21)
    args = ap.parse_args()

    print(f"{'config':<32}{'ed_current':>12}{'ed_canonical':>14}{'n':>5}")
    grand_cur, grand_can = [], []
    for cfg in args.configs:
        base = os.path.join(args.bench_dir, cfg, args.strategy)
        if not os.path.isdir(base):
            print(f"{cfg:<32}  (missing {base})"); continue
        cur, can = [], []
        nets = sorted(os.listdir(base))[:args.max_per_config]
        for net in nets:
            outp = os.path.join(base, net, "output.tre")
            gtp = os.path.join(base, net, "ground_truth.nex")
            if not (os.path.exists(outp) and os.path.exists(gtp)):
                continue
            try:
                rec = Tree(open(outp).read().strip(), format=1)
                tru = load_nexus_or_newick(gtp)
            except Exception:
                continue
            c = edit_distance(rec, tru, canonical=False)
            k = edit_distance(rec, tru, canonical=True)
            if c == c:
                cur.append(c)
            if k == k:
                can.append(k)
        mc = sum(cur) / len(cur) if cur else float("nan")
        mk = sum(can) / len(can) if can else float("nan")
        grand_cur += cur; grand_can += can
        print(f"{cfg:<32}{mc:>12.3f}{mk:>14.3f}{len(cur):>5}")
    if grand_cur:
        print(f"{'MEAN':<32}{sum(grand_cur)/len(grand_cur):>12.3f}"
              f"{sum(grand_can)/len(grand_can):>14.3f}{len(grand_cur):>5}")
    print("\nCompare the ed_canonical column between the two models. If the gap that")
    print("ed_current showed shrinks or vanishes under ed_canonical, it was the child-")
    print("order confound, not a real reconstruction improvement.")


if __name__ == "__main__":
    main()
