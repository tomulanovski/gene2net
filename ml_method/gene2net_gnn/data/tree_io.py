"""Tree I/O: parse Newick, extract species, convert to PyG format."""
from typing import Dict, List, Set, Tuple
from ete3 import Tree
import torch

def load_gene_trees_from_file(filepath: str, newick_format: int = 1) -> List[Tree]:
    """Load gene trees from a file (one Newick per line or semicolon-separated).

    Args:
        filepath: path to Newick file
        newick_format: ETE3 newick format code.
            0 = flexible (name, dist, support), 1 = name:dist (default),
            5 = internal+leaf names, 9 = leaf names only.
            See ETE3 docs for all options.
    """
    trees = []
    with open(filepath) as f:
        content = f.read()
    for part in content.split(";"):
        part = part.strip()
        if part:
            trees.append(Tree(part + ";", format=newick_format))
    return trees

def get_species_set(trees: List[Tree]) -> Set[str]:
    """Get the set of unique species labels across all trees."""
    species = set()
    for tree in trees:
        for leaf in tree.get_leaves():
            species.add(leaf.name)
    return species

def get_leaf_species_map(tree: Tree) -> Dict[str, str]:
    """Map leaf node names to species. For MUL-trees, multiple leaves share a species.
    Returns dict of {unique_leaf_id: species_name}.
    Since ETE3 leaves may share names, we use enumerated id.
    """
    mapping = {}
    for i, leaf in enumerate(tree.get_leaves()):
        mapping[f"{leaf.name}_{i}"] = leaf.name
    return mapping

def tree_to_edge_index(tree: Tree) -> Tuple[torch.Tensor, List[str]]:
    """Convert an ETE3 tree to a PyG-style edge_index tensor (undirected).
    Returns (edge_index [2, num_edges], node_names [num_nodes]).
    Nodes are indexed in traversal order.
    """
    node_to_idx = {}
    node_names = []

    for i, node in enumerate(tree.traverse("preorder")):
        node_to_idx[id(node)] = i
        node_names.append(node.name if node.name else f"internal_{i}")

    edges = []
    for node in tree.traverse("preorder"):
        for child in node.children:
            u = node_to_idx[id(node)]
            v = node_to_idx[id(child)]
            edges.append([u, v])
            edges.append([v, u])

    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    return edge_index, node_names


def reorder_edge_index_preorder(edge_index: torch.Tensor) -> torch.Tensor:
    """Reorder undirected edge pairs so the k-th pair is the k-th non-root node
    in preorder — matching the order used by the edge features and labels.

    ``tree_to_edge_index`` emits edges grouped by parent, which is a *different*
    order than the preorder-by-child order used in compute_species_tree_edge_*
    and the label extractor. That mismatch pairs each prediction's node
    embeddings with the wrong edge's features/label. This realigns them.

    Node 0 is the tree root (preorder index 0). Children order is recovered from
    the even-indexed (parent->child) edges, which preserves the original order.
    """
    children = {}
    for k in range(0, edge_index.shape[1], 2):
        p = int(edge_index[0, k])
        c = int(edge_index[1, k])
        children.setdefault(p, []).append(c)

    # Iterative preorder: emit (parent, child) the moment the child is visited.
    order = []
    stack = [(None, 0)]
    while stack:
        parent, node = stack.pop()
        if parent is not None:
            order.append((parent, node))
        for c in reversed(children.get(node, [])):
            stack.append((node, c))

    new_edges = []
    for p, c in order:
        new_edges.append([p, c])
        new_edges.append([c, p])
    return torch.tensor(new_edges, dtype=edge_index.dtype).t().contiguous()
