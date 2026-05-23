"""
metric_justification_panel.py - Build a thesis figure that demonstrates why
edit distance alone is insufficient for evaluating phylogenetic network
inference, and how the complementary metrics in `pairwise_compare`
decompose disagreement into three error categories: identity, count, parents.

Outputs (under scripts/figures/output/metric_justification/):
  - metric_justification_panel.{pdf,png}  (3x2 figure of folded networks)
  - metric_justification_table.{csv,tex}  (combined metric table)

Run from the repo root:
    python scripts/figures/metric_justification_panel.py
"""

import sys
from pathlib import Path

# Make the comparison engine importable without installing anything.
REPO_ROOT = Path(__file__).resolve().parents[2]
SIMS_SCRIPTS = REPO_ROOT / 'simulations' / 'scripts'
if str(SIMS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SIMS_SCRIPTS))

import matplotlib
matplotlib.use('Agg')  # non-interactive backend (cluster-safe, no X11)
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

from reticulate_tree import ReticulateTree
from compare_reticulations import (
    compare_num_rets,
    compare_ploidy_diff,
    compare_polyploid_species,
    match_and_compare,
    run_hungarian_on_groups,
)


def _tree_to_simple_graph(tree_obj: ReticulateTree) -> nx.DiGraph:
    """Convert the ete3 tree inside a ReticulateTree to a labelled DiGraph.

    Mirrors `ReticulateTree.get_edit_distance_multree`'s inner helper but is
    extracted here so we can compute GED without the signal-based timeout
    (which is POSIX-only and breaks the script on Windows). Safe for tiny
    toy trees.
    """
    G = nx.DiGraph()
    for node in tree_obj.tree.traverse():
        nid = id(node)
        G.add_node(nid, label=node.name if node.is_leaf() else None)
        if not node.is_root():
            G.add_edge(id(node.up), nid)
    return G


def safe_edit_distance_multree(t1: ReticulateTree, t2: ReticulateTree,
                                normalize: bool = True) -> float:
    """Cross-platform graph edit distance on the underlying MUL-trees.

    Calls `networkx.optimize_graph_edit_distance` directly; no SIGALRM.
    Intended for use only on the toy trees defined in this module.
    """
    g1 = _tree_to_simple_graph(t1)
    g2 = _tree_to_simple_graph(t2)
    distance = next(nx.optimize_graph_edit_distance(
        g1, g2,
        node_match=lambda u, v: u.get('label') == v.get('label'),
    ))
    if normalize:
        norm = max(len(g1.nodes) + len(g1.edges),
                   len(g2.nodes) + len(g2.edges))
        if norm > 0:
            distance /= norm
    return distance


# ───────────────────────────────────────────────────────────────────────────
# Toy networks. Each is a MUL-tree Newick string; duplicated leaf labels mark
# polyploid species. ReticulateTree auto-folds MUL-trees to DAGs on load.
# ───────────────────────────────────────────────────────────────────────────

# Single 6-species truth: A = (((q,p),(r,p)),(s,(t,u)));
# p is the polyploid; its parental clades are q and r.
A_NEWICK = "(((q,p),(r,p)),(s,(t,u)));"

# Three inferred networks, each introducing exactly one error category.
B_NEWICK = "(((q,p),(s,p)),(r,(t,u)));"          # Parents error: p hybrid of q,s instead of q,r
C_NEWICK = "(((q,u),(r,u)),(s,(t,p)));"          # Identity error: u is polyploid instead of p
D_NEWICK = "(((q,p),(r,p)),((s,u),(t,u)));"      # Count error: extra spurious reticulation on u


EXAMPLES = [
    {
        'label': 'Ex.1 parents error',
        'category': 'parents',
        'truth_name': 'A',
        'inferred_name': 'B',
        'truth_newick': A_NEWICK,
        'inferred_newick': B_NEWICK,
        'caption': 'wrong reticulation sister clade',
    },
    {
        'label': 'Ex.2 identity error',
        'category': 'identity',
        'truth_name': 'A',
        'inferred_name': 'C',
        'truth_newick': A_NEWICK,
        'inferred_newick': C_NEWICK,
        'caption': 'wrong species is polyploid',
    },
    {
        'label': 'Ex.3 count error',
        'category': 'count',
        'truth_name': 'A',
        'inferred_name': 'D',
        'truth_newick': A_NEWICK,
        'inferred_newick': D_NEWICK,
        'caption': 'spurious extra reticulation',
    },
]


OUTPUT_DIR = REPO_ROOT / 'scripts' / 'figures' / 'output' / 'metric_justification'


def draw_reticulate_tree(tree: ReticulateTree, ax, title: str = '') -> None:
    """Draw the folded DAG of a ReticulateTree onto a given matplotlib Axes.

    Mirrors ReticulateTree.visualize() but targets `ax` and never calls
    plt.show(), so it composes inside subplots.
    """
    G = tree.dag
    depths = ReticulateTree.compute_depths(G)
    for node, d in depths.items():
        G.nodes[node]['layer'] = d
    pos = nx.multipartite_layout(G, subset_key='layer')

    node_colors = []
    node_shape_of = {}
    for n in G.nodes():
        if ReticulateTree.is_reticulation_node(G, n):
            node_shape_of[n] = 'D'
            node_colors.append('orange')
        elif G.out_degree(n) == 0:
            node_shape_of[n] = 'o'
            node_colors.append('lightgreen')
        else:
            node_shape_of[n] = 'o'
            node_colors.append('lightblue')

    nodes_list = list(G.nodes())
    color_of = {n: c for n, c in zip(nodes_list, node_colors)}
    for shape in set(node_shape_of.values()):
        shaped = [n for n in nodes_list if node_shape_of[n] == shape]
        nx.draw_networkx_nodes(
            G, pos, nodelist=shaped,
            node_shape=shape,
            node_color=[color_of[n] for n in shaped],
            node_size=600,
            ax=ax,
        )
    nx.draw_networkx_edges(G, pos, arrows=True, ax=ax)

    labels = {n: G.nodes[n].get('label', '') or '' for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, ax=ax)

    if title:
        ax.set_title(title, fontsize=11)
    ax.set_axis_off()


METRIC_ORDER = [
    'edit_distance_multree',
    'num_rets_diff',
    'num_rets_bias',
    'polyploid_species_jaccard',
    'ploidy_diff.dist',
    'ploidy_diff.FP',
    'ploidy_diff.FN',
    'ploidy_diff.TP',
    'ret_leaf_jaccard.dist',
    'ret_sisters_jaccard.dist',
]


def compute_pair_metrics(truth: ReticulateTree,
                         inferred: ReticulateTree) -> dict:
    """Build the full metric vector for one (truth, inferred) pair.

    Mirrors `pairwise_compare` from compare_reticulations.py minus
    rf_distance (excluded per spec) and using the Windows-safe edit
    distance helper.
    """
    num_rets = compare_num_rets(
        truth.get_reticulation_count(),
        inferred.get_reticulation_count(),
    )
    ploidy = compare_ploidy_diff(
        truth.get_leaf_counts(),
        inferred.get_leaf_counts(),
    )
    poly_jac = compare_polyploid_species(
        truth.get_leaf_counts(),
        inferred.get_leaf_counts(),
    )

    truth_ret_leaves = truth.get_reticulation_leaves()
    inferred_ret_leaves = inferred.get_reticulation_leaves()
    precomputed = run_hungarian_on_groups(
        truth_ret_leaves.values(),
        inferred_ret_leaves.values(),
    )
    ret_leaf = match_and_compare(
        truth_ret_leaves, inferred_ret_leaves,
        precomputed_match=precomputed,
    )
    ret_sisters = match_and_compare(
        truth_ret_leaves, inferred_ret_leaves,
        truth.get_reticulation_sisters(),
        inferred.get_reticulation_sisters(),
        precomputed_match=precomputed,
    )

    return {
        'edit_distance_multree': safe_edit_distance_multree(truth, inferred),
        'num_rets_diff': num_rets['abs'],
        'num_rets_bias': num_rets['signed'],
        'polyploid_species_jaccard': poly_jac,
        'ploidy_diff.dist': ploidy['dist'],
        'ploidy_diff.FP': ploidy['FP'],
        'ploidy_diff.FN': ploidy['FN'],
        'ploidy_diff.TP': ploidy['TP'],
        'ret_leaf_jaccard.dist': ret_leaf['dist'],
        'ret_sisters_jaccard.dist': ret_sisters['dist'],
    }


def build_table_dataframe(examples: list) -> pd.DataFrame:
    """Return a long-form DataFrame with one row per (example, metric)."""
    rows = []
    for ex in examples:
        truth = ReticulateTree(ex['truth_newick'])
        inferred = ReticulateTree(ex['inferred_newick'])
        metrics = compute_pair_metrics(truth, inferred)
        for metric_name in METRIC_ORDER:
            rows.append({
                'example': ex['label'],
                'pair': f"{ex['truth_name']} vs {ex['inferred_name']}",
                'metric': metric_name,
                'value': metrics[metric_name],
            })
    return pd.DataFrame(rows)


def export_table(df: pd.DataFrame, output_dir: Path) -> tuple:
    """Pivot to wide form, write CSV + LaTeX tabular block.

    Returns (csv_path, tex_path).
    """
    wide = df.pivot(index='metric', columns='example', values='value')
    # preserve the metric order from METRIC_ORDER for both files
    wide = wide.loc[[m for m in METRIC_ORDER if m in wide.index]]

    csv_path = output_dir / 'metric_justification_table.csv'
    wide.to_csv(csv_path, float_format='%.4f')

    tex_path = output_dir / 'metric_justification_table.tex'
    with open(tex_path, 'w', encoding='utf-8') as fh:
        # to_latex on the wide frame; caller can drop into a \begin{table} env.
        fh.write(wide.to_latex(float_format='%.4f', na_rep='--'))
    return csv_path, tex_path


def render_panel(examples: list, output_dir: Path) -> tuple:
    """Render the 3-row x 2-column figure and save PDF + PNG.

    Returns (pdf_path, png_path).
    """
    fig, axes = plt.subplots(
        nrows=len(examples), ncols=2,
        figsize=(9, 3.2 * len(examples)),
        constrained_layout=True,
    )
    # axes is shape (n_rows, 2); guarantee 2D even for n_rows == 1
    if len(examples) == 1:
        axes = [axes]

    for row, ex in enumerate(examples):
        truth = ReticulateTree(ex['truth_newick'])
        inferred = ReticulateTree(ex['inferred_newick'])

        draw_reticulate_tree(
            truth, axes[row][0],
            title=f"{ex['truth_name']}: truth",
        )
        draw_reticulate_tree(
            inferred, axes[row][1],
            title=f"{ex['inferred_name']}: inferred — {ex['caption']}",
        )

    # row labels on the left for the three error categories
    for row, ex in enumerate(examples):
        axes[row][0].annotate(
            ex['label'],
            xy=(-0.08, 0.5), xycoords='axes fraction',
            ha='right', va='center',
            fontsize=12, fontweight='bold', rotation=90,
        )

    pdf_path = output_dir / 'metric_justification_panel.pdf'
    png_path = output_dir / 'metric_justification_panel.png'
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return pdf_path, png_path


def main():
    print(f"[metric_justification_panel] repo root: {REPO_ROOT}")
    print(f"[metric_justification_panel] output dir: {OUTPUT_DIR}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for ex in EXAMPLES:
        truth = ReticulateTree(ex['truth_newick'])
        inferred = ReticulateTree(ex['inferred_newick'])
        metrics = compute_pair_metrics(truth, inferred)
        print(f"\n{ex['label']} ({ex['truth_name']} vs {ex['inferred_name']})")
        for k in METRIC_ORDER:
            print(f"  {k:<32} {metrics[k]}")

    # ─── enforce the spec's qualitative signatures so future edits to
    # ─── upstream comparison functions can't silently break the figure
    cached = {ex['label']: compute_pair_metrics(
                  ReticulateTree(ex['truth_newick']),
                  ReticulateTree(ex['inferred_newick']))
              for ex in EXAMPLES}

    m1 = cached['Ex.1 parents error']
    assert m1['num_rets_diff'] == 0, m1
    assert m1['polyploid_species_jaccard'] == 0.0, m1
    assert m1['ret_leaf_jaccard.dist'] == 0.0, m1
    assert m1['ret_sisters_jaccard.dist'] > 0.0, m1

    m2 = cached['Ex.2 identity error']
    assert m2['num_rets_diff'] == 0, m2
    assert m2['num_rets_bias'] == 0, m2
    assert m2['polyploid_species_jaccard'] == 1.0, m2
    assert m2['ret_leaf_jaccard.dist'] > 0.5, m2

    m3 = cached['Ex.3 count error']
    assert m3['num_rets_diff'] == 1, m3
    assert m3['num_rets_bias'] == 1, m3
    assert 0.4 <= m3['polyploid_species_jaccard'] <= 0.6, m3

    print("\n[assertions] all metric signatures match the spec.")

    pdf_path, png_path = render_panel(EXAMPLES, OUTPUT_DIR)
    print(f"\n[panel] wrote {pdf_path}")
    print(f"[panel] wrote {png_path}")

    df_long = build_table_dataframe(EXAMPLES)
    csv_path, tex_path = export_table(df_long, OUTPUT_DIR)
    print(f"[table] wrote {csv_path}")
    print(f"[table] wrote {tex_path}")


if __name__ == '__main__':
    main()
