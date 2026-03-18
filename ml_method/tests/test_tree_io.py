from gene2net_gnn.data.tree_io import (
    load_gene_trees_from_file,
    get_species_set,
    get_leaf_species_map,
    tree_to_edge_index,
)

def test_load_gene_trees(tmp_path):
    f = tmp_path / "trees.nwk"
    f.write_text("((A:1,B:1):1,C:1);\n((A:1,C:1):1,B:1);\n")
    trees = load_gene_trees_from_file(str(f))
    assert len(trees) == 2

def test_get_species_set(sample_gene_trees):
    species = get_species_set(sample_gene_trees)
    assert species == {"A", "B", "C", "D", "E"}

def test_get_leaf_species_map(sample_gene_trees):
    mapping = get_leaf_species_map(sample_gene_trees[0])
    d_leaves = [name for name, sp in mapping.items() if sp == "D"]
    assert len(d_leaves) == 2

def test_tree_to_edge_index(simple_species_tree):
    edge_index, node_names = tree_to_edge_index(simple_species_tree)
    assert edge_index.shape[0] == 2
    assert edge_index.shape[1] == 16  # 5 leaves + 4 internal = 9 nodes, 8 edges undirected = 16
