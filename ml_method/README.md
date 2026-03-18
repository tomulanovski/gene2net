# Gene2Net-GNN

A GNN-based method for phylogenetic network inference. Given a collection of gene trees
and a species tree (e.g., from ASTRAL 4 or any species tree method), predicts Whole Genome
Duplication (WGD) events and outputs a MUL-tree (multiply-labeled tree).

## The Idea

Existing methods (GRAMPA, Polyphest, PADRE) use combinatorial approaches to infer MUL-trees.
This method uses Graph Neural Networks to **learn** which edges in a species tree had WGD events
and where each duplicated clade was placed.

**Three-module architecture:**

```
Gene trees (hundreds)          Species tree (e.g., ASTRAL 4)
       |                                    |
  [Module 1]                                |
  Gene Tree Encoder                         |
  (bottom-up + top-down GNN)                |
       |                                    |
  per-leaf embeddings                       |
       |                                    |
  [Module 2]                                |
  Cross-Tree Aggregator                     |
  (multi-head attention pooling             |
   per species across all trees)            |
       |                                    |
  per-species embeddings -------> [Module 3]
                                  Species Tree GNN (GAT)
                                       |
                              per-edge embeddings
                                    /       \
                          [WGD Head]         [Partner Head]
                          How many events    Where does the
                          on this edge?      duplicate attach?
                                    \       /
                              [MUL-tree Builder]
                                       |
                                  Output MUL-tree
```

**Key predictions per species tree edge:**
- WGD count: 0, 1, 2, or 3 events per edge (focal loss for class imbalance). `max_wgd_events=3` is configurable.
- Partner edge: which other edge the duplicated clade attaches to (pointer network)
- Auto vs allo: implicit from partner direction (self = auto, distant = allo)
- No ploidy prior required — the model learns WGD patterns from data. Prior knowledge can be added later via extended features.

## How It Works

### Forward pass (step by step)

**Module 1 — Gene Tree Encoder** (`tree_encoder.py`):
Each gene tree is processed independently. Every leaf gets an initial embedding from a learned
species embedding table (all copies of species A start identical). Internal nodes get a shared
learned "internal" embedding. Then two rounds of message passing:
- **Bottom-up** (leaves → root): aggregates "what species are below me?"
- **Top-down** (root → leaves): propagates "what's the evolutionary context above me?"

Each leaf's final embedding is the concatenation of its bottom-up and top-down representations,
projected to `embed_dim`. Branch lengths are used as edge features in message passing.

Gene loss is handled naturally — if a species is absent from a gene tree, it simply contributes
no embeddings from that tree. The aggregator handles variable counts.

**Module 2 — Cross-Tree Aggregator** (`aggregator.py`):
For each species, all its leaf embeddings across all gene trees are collected (e.g., species A
appeared in 280 of 300 gene trees → 280 embeddings). Multi-head attention with K learned query
vectors pools these into a fixed-size output (`K × embed_dim`). This captures how consistently
a species appears and how its gene tree context varies — key signals for WGD.

Species absent from all gene trees get a zero vector (handled gracefully).

**Module 3 — Species Tree GNN** (`species_gnn.py`):
Species tree nodes get their aggregated embeddings concatenated with hand-crafted features
(13-dim, see below). A GNN (configurable: GAT, GCN, GIN, or GraphSAGE) propagates information
along the species tree topology with residual connections. Then edge embeddings are computed by
concatenating parent + child node embeddings and passing through an MLP, optionally with
edge-level features (concordance, branch length, clade size, depth).

**Prediction Heads** (`heads.py`):
- **WGD Head**: MLP mapping edge embedding → 4-class logits (0, 1, 2, or 3 events)
- **Partner Head**: Bilinear attention between all edge pairs → logits matrix [E × E].
  For each WGD edge, the partner is the argmax of its row.

**MUL-tree Builder** (`mul_tree_builder.py`):
Takes the species tree + predicted WGD events and constructs the output MUL-tree by duplicating
clades and attaching them at partner locations.

### Features

**Node features (13-dim, per species tree node):**

Copy count statistics (8 features) — computed per species across all gene trees:
| Feature | What it captures |
|---------|-----------------|
| `mean_copies` | Average copy count. Diploids ≈ 1, polyploids > 1 |
| `var_copies` | Variance in copy count across trees. High = noisy signal |
| `mode_copies` | Most common copy count |
| `p_absent` | Fraction of trees where species is missing (gene loss signal) |
| `p_1_copy` | Fraction with exactly 1 copy (diploid signal) |
| `p_2_copies` | Fraction with exactly 2 copies (WGD signal) |
| `p_3plus_copies` | Fraction with 3+ copies (higher ploidy signal) |
| `max_copies` | Maximum observed copies in any tree |

Clustering summary (5 features) — for each species, compute co-clustering frequency
with every other species (how often they're siblings in gene trees), then summarize:
| Feature | What it captures |
|---------|-----------------|
| `mean` | Average co-clustering across all other species |
| `std` | How variable the co-clustering is |
| `max` | Strongest co-clustering partner |
| `min` | Weakest co-clustering partner |
| `median` | Typical co-clustering level |

Internal (non-leaf) nodes get zeros — the GNN fills them via message passing from leaves.

**Edge features (4-dim, per species tree edge):**
| Feature | What it captures |
|---------|-----------------|
| `concordance_factor` | Fraction of gene trees supporting this bipartition. Low concordance near WGD events |
| `branch_length` | Branch length in the species tree |
| `clade_size` | Number of leaves below this edge |
| `depth` | Distance from root (number of edges) |

### Loss function

Three components combined: `L = λ_wgd × L_wgd + λ_partner × L_partner + λ_count × L_count`

- **L_wgd**: Focal loss on WGD count prediction (handles ~95% of edges having 0 events)
- **L_partner**: Cross-entropy on partner edge prediction (only computed for edges with WGD)
- **L_count**: MSE between predicted and true total event count (soft regularizer)

### Variable species count

The model handles datasets with varying numbers of species across samples:
- Species embedding table is sized by `max_species` (config, default 100) — an upper bound
- Node features are always 13-dim (fixed), not dependent on species count
- Species absent from gene trees (gene loss) get zero aggregated embeddings
- The `max_species` parameter only affects embedding table size — increase if needed

## Quick Start

### Run tests
```bash
cd ml_method
PYTHONPATH=$(pwd) python -m pytest tests/ -v
```

### Inference (once you have a trained model)
```bash
python scripts/infer.py \
    --gene-trees my_gene_trees.nwk \
    --species-tree my_astral.nwk \
    --model output/best_model.pt \
    --output predicted_mul_tree.nwk
```

### Training
```bash
# 1. Generate training data on the cluster
sbatch scripts/generate_training_data.sh /path/to/output

# 2. Train
python scripts/train.py \
    --data-dir /path/to/training_data \
    --config configs/default.yaml \
    --output-dir output/ \
    --device cuda
```

### Evaluate against ground truth
```bash
# Quick pipeline test (random model, no training needed)
python scripts/evaluate.py --mock --networks-dir ../simulations/networks/

# Evaluate trained model on a specific network
python scripts/evaluate.py \
    --model output/best_model.pt \
    --network Shahrestani_2015 \
    --gene-trees /path/to/gene_trees.nwk \
    --species-tree /path/to/astral.nwk
```

## File Guide

### Core Model (`gene2net_gnn/model/`)
| File | What it does |
|------|-------------|
| `gene2net_model.py` | **Start here.** Full model forward pass connecting all modules |
| `tree_encoder.py` | Module 1: GNN on individual gene trees (bottom-up + top-down message passing) |
| `aggregator.py` | Module 2: Multi-head attention pooling across gene trees per species |
| `species_gnn.py` | Module 3: GNN on species tree (GAT/GCN/GIN/SAGE), produces edge embeddings |
| `heads.py` | WGD count head (ordinal classification) + partner pointer head |

### Data Pipeline (`gene2net_gnn/data/`)
| File | What it does |
|------|-------------|
| `tree_io.py` | Parse Newick trees, convert to PyG edge_index tensors |
| `mul_tree_generator.py` | Generate random MUL-trees for training (birth-death + WGD) |
| `label_extractor.py` | Decompose MUL-tree into WGD events, map to ASTRAL edges (training labels) |
| `features.py` | Hand-crafted features: copy counts, clustering profiles, concordance |
| `dataset.py` | PyG Dataset wrapping one training example |
| `collate.py` | Custom batching for variable-size gene tree sets |

### Training (`gene2net_gnn/training/`)
| File | What it does |
|------|-------------|
| `loss.py` | Focal loss (class imbalance), partner CE, event count MSE |
| `trainer.py` | Training loop with early stopping, LR scheduling, checkpointing |
| `metrics.py` | WGD F1, partner accuracy, event count error |

### Inference (`gene2net_gnn/inference/`)
| File | What it does |
|------|-------------|
| `predict.py` | End-to-end: gene trees + species tree + model -> MUL-tree |
| `mul_tree_builder.py` | Convert predicted WGD events into an actual MUL-tree |

### Scripts (`scripts/`)
| File | What it does |
|------|-------------|
| `train.py` | Training CLI entry point |
| `infer.py` | Inference CLI: gene trees in, MUL-tree out |
| `evaluate.py` | Compare predictions against ground truth networks |
| `generate_one_example.py` | Generate one training example (SimPhy + species tree inference) |
| `generate_training_data.sh` | SLURM array job for bulk training data generation |

## Configuration

All hyperparameters in `configs/default.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Model** | | |
| `embed_dim` | 64 | Embedding dimension throughout the model |
| `conv_type` | gat | GNN architecture: gat, gcn, gin, sage |
| `num_gat_layers` | 3 | GNN layers in Module 3 |
| `num_attention_heads` | 8 | Attention heads in Module 2 (K) |
| `dropout` | 0.2 | Dropout rate |
| **Training** | | |
| `lr` | 1e-3 | Learning rate |
| `weight_decay` | 1e-4 | L2 regularization |
| `batch_size` | 4 | Gradient accumulation steps (effective batch size) |
| `max_epochs` | 200 | Maximum training epochs |
| `patience` | 20 | Early stopping patience (epochs without improvement) |
| `num_gene_trees` | 300 | Gene trees subsampled per training step |
| `focal_alpha` | 0.25 | Focal loss alpha (class balance) |
| `focal_gamma` | 2.0 | Focal loss gamma (hard example focus) |
| `lambda_wgd` | 1.0 | WGD loss weight |
| `lambda_partner` | 1.0 | Partner loss weight |
| `lambda_count` | 0.1 | Count loss weight |
| **Data** | | |
| `max_species` | 100 | Upper bound on species count (sizes embedding table) |
| `max_ploidy` | 10 | Maximum ploidy level for synthetic data generation |

## Training Data

Training data is generated synthetically:
1. Generate random MUL-trees with `mul_tree_generator.py`
2. Simulate gene trees using SimPhy (adds ILS + gene duplication/loss noise)
3. Infer species tree (e.g., ASTRAL 4 or any species tree method)
4. Extract ground truth labels by decomposing the known MUL-tree

The plan targets 500 base MUL-trees x 3 ILS levels = 1500 training examples.

## Experiment Plan

See `docs/superpowers/plans/2026-03-17-gnn-network-inference.md` (Chunk 4) for:
- Synthetic data validation targets (WGD F1 > 0.5)
- Ablation studies (features, attention heads, gene tree count, GNN depth)
- Leave-one-out evaluation on 21 biological networks
- Comparison against GRAMPA, Polyphest, PADRE, MPSUGAR

## Design Document

Full technical design with biological motivation, loss function derivations,
and risk analysis: `docs/superpowers/specs/2026-03-17-gnn-network-inference-design.md`

## Dependencies

- Python 3.10+
- PyTorch
- PyTorch Geometric
- ETE3
- PyYAML
- SimPhy (cluster, for training data generation)
- ASTRAL 4 or equivalent (cluster, for species tree inference)
