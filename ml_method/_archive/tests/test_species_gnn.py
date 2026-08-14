import torch
from gene2net_gnn.model.species_gnn import SpeciesTreeGNN
from gene2net_gnn.model.heads import WGDHead, PartnerHead

def test_species_gnn_output_shape():
    gnn = SpeciesTreeGNN(input_dim=64, hidden_dim=64, num_layers=3)
    # 9 nodes, undirected edges
    x = torch.randn(9, 64)
    edge_index = torch.tensor([
        [0,1,0,2,1,3,1,4,2,5,5,6,5,7,7,8],
        [1,0,2,0,3,1,4,1,5,2,6,5,7,5,8,7]
    ], dtype=torch.long)
    node_out, edge_out = gnn(x, edge_index)
    assert node_out.shape == (9, 64)
    assert edge_out.shape == (8, 64)  # one per undirected edge

def test_wgd_head():
    head = WGDHead(edge_dim=64, max_events=3)
    edge_embeddings = torch.randn(8, 64)
    logits = head(edge_embeddings)
    assert logits.shape == (8, 4)  # 0, 1, 2, 3 events

def test_partner_head():
    head = PartnerHead(edge_dim=64)
    edge_embeddings = torch.randn(8, 64)
    query_edge = edge_embeddings[2:3]
    scores = head(query_edge, edge_embeddings)
    assert scores.shape == (1, 8)

def test_partner_head_batch():
    head = PartnerHead(edge_dim=64)
    edge_embeddings = torch.randn(8, 64)
    query_edges = edge_embeddings[:3]  # 3 WGD edges
    scores = head(query_edges, edge_embeddings)
    assert scores.shape == (3, 8)

def test_gnn_gradient():
    gnn = SpeciesTreeGNN(input_dim=32, hidden_dim=32, num_layers=2)
    x = torch.randn(5, 32, requires_grad=True)
    edge_index = torch.tensor([[0,1,1,2,2,3,3,4,1,0,2,1,3,2,4,3]], dtype=torch.long).reshape(2,-1)
    node_out, edge_out = gnn(x, edge_index)
    (node_out.sum() + edge_out.sum()).backward()
    assert x.grad is not None
