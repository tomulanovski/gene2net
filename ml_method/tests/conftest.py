import pytest
from ete3 import Tree

@pytest.fixture
def simple_species_tree():
    """5-species rooted binary tree: ((A,B),(C,(D,E)));"""
    return Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)

@pytest.fixture
def simple_mul_tree():
    """MUL-tree: species D is allotetraploid (copy near A/B clade).
    ((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);
    """
    return Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)

@pytest.fixture
def auto_mul_tree():
    """MUL-tree: autopolyploidy of clade (D,E).
    (((D:1,E:1):1,(D:1,E:1):1):1,(A:1,(B:1,C:1):1):1);
    """
    return Tree("(((D:1,E:1):1,(D:1,E:1):1):1,(A:1,(B:1,C:1):1):1);", format=1)

@pytest.fixture
def sample_gene_trees():
    """10 simple gene trees with some copies of D."""
    trees = []
    for _ in range(7):
        trees.append(Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1))
    for _ in range(3):
        trees.append(Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1))
    return trees
