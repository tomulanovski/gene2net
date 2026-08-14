"""Full Gene2Net model combining Modules 1-3 and prediction heads."""
import torch
import torch.nn as nn
from collections import defaultdict

from gene2net_gnn.model.tree_encoder import GeneTreeEncoder
from gene2net_gnn.model.aggregator import CrossTreeAggregator
from gene2net_gnn.model.species_gnn import SpeciesTreeGNN
from gene2net_gnn.model.heads import WGDHead, PartnerHead


class Gene2NetModel(nn.Module):
    def __init__(
        self,
        n_species: int,
        embed_dim: int = 64,
        n_attention_heads: int = 8,
        n_gat_layers: int = 3,
        handcrafted_dim: int = 13,
        max_wgd_events: int = 3,
        dropout: float = 0.2,
        conv_type: str = "gat",
        edge_feat_dim: int = 4,
    ):
        super().__init__()
        self.n_species = n_species
        self.embed_dim = embed_dim
        self.n_attention_heads = n_attention_heads

        # Module 1: Gene tree encoder
        self.tree_encoder = GeneTreeEncoder(n_species, embed_dim)

        # Module 2: Cross-tree aggregator (one shared instance, applied per species)
        self.aggregator = CrossTreeAggregator(embed_dim, n_attention_heads)

        # Species tree input dim = aggregator output + hand-crafted features
        aggregator_out_dim = n_attention_heads * embed_dim
        species_input_dim = aggregator_out_dim + handcrafted_dim

        # Module 3: Species tree GNN
        self.species_gnn = SpeciesTreeGNN(
            input_dim=species_input_dim,
            hidden_dim=embed_dim,
            num_layers=n_gat_layers,
            conv_type=conv_type,
            dropout=dropout,
            edge_feat_dim=edge_feat_dim,
        )

        # Prediction heads
        self.wgd_head = WGDHead(embed_dim, max_events=max_wgd_events)
        self.partner_head = PartnerHead(embed_dim)

    def forward(self, gene_tree_data_list, species_tree_data):
        """
        Args:
            gene_tree_data_list: list of dicts, each with:
                - edge_index: [2, 2E]
                - species_ids: [N] (-1 for internal)
                - branch_lengths: [2E]
                - leaf_mask: [N] bool
            species_tree_data: dict with:
                - edge_index: [2, 2E_sp]
                - node_features: [N_sp, handcrafted_dim]
                - edge_features: [E_sp, edge_feat_dim] (optional, concordance/branch_len/etc.)
                - is_leaf: [N_sp] bool
                - species_ids: [N_sp] (-1 for internal)

        Returns:
            wgd_logits: [E_sp, max_events+1]
            partner_logits: [E_sp, E_sp] (all edges as potential partners)
            edge_embeddings: [E_sp, embed_dim]
        """
        device = species_tree_data["node_features"].device

        # Step 1: Encode all gene trees
        # Collect per-species embeddings across all gene trees
        species_embeddings = defaultdict(list)  # species_idx -> list of [embed_dim] tensors

        for gt_data in gene_tree_data_list:
            leaf_emb = self.tree_encoder(
                gt_data["edge_index"],
                gt_data["species_ids"],
                gt_data["branch_lengths"],
            )
            # leaf_emb: [n_leaves_in_this_tree, embed_dim]
            # Map back to species
            leaf_mask = gt_data["leaf_mask"]
            leaf_species = gt_data["species_ids"][leaf_mask]
            for i, sp_idx in enumerate(leaf_species):
                sp_idx_val = sp_idx.item()
                if sp_idx_val >= 0:
                    species_embeddings[sp_idx_val].append(leaf_emb[i])

        # Step 2: Aggregate per species
        n_sp_nodes = species_tree_data["node_features"].shape[0]
        aggregator_out_dim = self.n_attention_heads * self.embed_dim
        aggregated = torch.zeros(n_sp_nodes, aggregator_out_dim, device=device)

        sp_is_leaf = species_tree_data["is_leaf"]
        sp_species_ids = species_tree_data["species_ids"]

        for node_idx in range(n_sp_nodes):
            if sp_is_leaf[node_idx]:
                sp_idx = sp_species_ids[node_idx].item()
                if sp_idx in species_embeddings and len(species_embeddings[sp_idx]) > 0:
                    stacked = torch.stack(species_embeddings[sp_idx])
                    aggregated[node_idx] = self.aggregator(stacked)
                else:
                    aggregated[node_idx] = self.aggregator(torch.zeros(0, self.embed_dim, device=device))

        # Step 3: Build species tree node features
        # Concatenate aggregated embeddings + hand-crafted features
        node_features = torch.cat([aggregated, species_tree_data["node_features"]], dim=-1)

        # Step 4: Run species tree GNN (with optional edge features)
        edge_features = species_tree_data.get("edge_features", None)
        node_emb, edge_emb = self.species_gnn(node_features, species_tree_data["edge_index"], edge_features)

        # Step 5: Prediction heads
        wgd_logits = self.wgd_head(edge_emb)          # [E, max_events+1]
        partner_logits = self.partner_head(edge_emb, edge_emb)  # [E, E]

        return wgd_logits, partner_logits, edge_emb
