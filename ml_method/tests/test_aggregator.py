import torch
from gene2net_gnn.model.aggregator import CrossTreeAggregator

def test_aggregator_output_shape():
    agg = CrossTreeAggregator(embed_dim=32, n_heads=8)
    embeddings = torch.randn(15, 32)
    output = agg(embeddings)
    assert output.shape == (8 * 32,)  # K heads × embed_dim

def test_aggregator_single_copy():
    """Diploid species: single copy per tree."""
    agg = CrossTreeAggregator(embed_dim=32, n_heads=8)
    embeddings = torch.randn(10, 32)
    output = agg(embeddings)
    assert output.shape == (8 * 32,)

def test_aggregator_empty():
    """Species absent from all gene trees -> zero vector."""
    agg = CrossTreeAggregator(embed_dim=32, n_heads=8)
    embeddings = torch.zeros(0, 32)
    output = agg(embeddings)
    assert output.shape == (8 * 32,)
    assert torch.all(output == 0)

def test_aggregator_permutation_invariant():
    """Output should not depend on order of input embeddings."""
    agg = CrossTreeAggregator(embed_dim=32, n_heads=4)
    agg.eval()
    embeddings = torch.randn(20, 32)
    perm = torch.randperm(20)
    out1 = agg(embeddings)
    out2 = agg(embeddings[perm])
    assert torch.allclose(out1, out2, atol=1e-5)

def test_aggregator_gradient():
    agg = CrossTreeAggregator(embed_dim=16, n_heads=4)
    embeddings = torch.randn(10, 16, requires_grad=True)
    output = agg(embeddings)
    output.sum().backward()
    assert embeddings.grad is not None
