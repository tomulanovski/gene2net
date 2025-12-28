# ML Algorithm Design for Polyploidy Inference using GNNs

## Task Definition

**Problem Type**: Supervised Graph-to-Graph Prediction / Tree Reconstruction

**Input**: 
- ~1000 gene trees (Newick format) per replicate
- Each gene tree is a phylogenetic tree with species as leaves

**Output**: 
- Single MUL tree (Multi-labeled tree) representing the polyploidy network
- MUL tree has duplicate leaf labels indicating polyploid species

**Ground Truth**: 
- Known MUL trees for 21 networks
- 5 replicates per network = 105 training examples (can be split train/val/test)

**Task Category**: 
- **Graph-to-Graph Generation**: Multiple input graphs → Single output graph
- **Structured Prediction**: Output has complex structure (tree topology + copy numbers)
- **Multi-instance Learning**: Multiple gene trees (instances) → One MUL tree (bag label)

---

## Architecture Overview

### High-Level Pipeline

```
1000 Gene Trees (Newick)
    ↓
[Tree Encoder] → 1000 Graph Embeddings [each: d_hidden]
    ↓
[Multi-Tree Aggregator] → Single Aggregated Embedding [d_hidden]
    ↓
[MUL Tree Decoder] → MUL Tree Structure + Copy Numbers
    ↓
Output: MUL Tree (Newick format)
```

### Key Design Decisions

1. **Encode each gene tree independently** → More efficient, parallelizable
2. **Aggregate 1000 embeddings** → Learn consensus across gene trees
3. **Decode to MUL tree** → Generate tree structure + copy numbers

---

## Part 1: Gene Tree Encoding (1000 trees → 1000 embeddings)

### Step 1: Tree to Graph Conversion

For each of the 1000 gene trees:

```
Gene Tree (Newick) → Parse with ETE3 → Directed Graph
```

**Graph Structure**:
- **Nodes**: All tree nodes (leaves + internal nodes)
- **Edges**: Parent → Child relationships (directed)
- **Node Features**: See below
- **Edge Features**: Branch lengths, edge types

### Step 2: Node Feature Engineering

Each node gets a feature vector:

**Leaf Nodes**:
- Species embedding (learned, d_species=32): One-hot species ID → Embedding layer
- Node type: [1, 0] (is_leaf, is_internal)
- Depth: Normalized distance from root
- Branch length: Normalized distance to parent
- Number of descendants: 0 (always)

**Internal Nodes**:
- Node type: [0, 1]
- Depth: Normalized distance from root
- Branch length: Normalized distance to parent
- Number of descendants: Normalized count
- Number of children: Typically 2 (binary trees)
- Species embedding: Zero vector (no species)

**Total Node Feature Dimension**: d_node = 32 (species) + 2 (type) + 1 (depth) + 1 (branch) + 1 (descendants) + 1 (children) = **38 dimensions**

### Step 3: Graph Neural Network Encoding

**Architecture**: Graph Attention Network (GAT) with multiple layers

```python
# For each gene tree:
node_features [N, 38] 
    ↓
Linear Projection [N, 38] → [N, 128]
    ↓
GAT Layer 1: Message passing, attention weights
    ↓
GAT Layer 2: Deeper relationships
    ↓
GAT Layer 3: Global tree structure
    ↓
Node Embeddings [N, 128]
```

**GAT Layer Details**:
- **Input**: Node features [N, d_hidden], Edge index [2, E], Edge features [E, d_edge]
- **Attention Mechanism**: 
  - For each node, compute attention weights over neighbors
  - Weighted aggregation of neighbor features
  - Multi-head attention (4-8 heads) for richer representations
- **Output**: Updated node embeddings [N, d_hidden]

**Why GAT?**
- Attention learns which neighbors are most important
- Handles variable tree sizes naturally
- Captures both local (immediate neighbors) and global (through layers) structure

### Step 4: Graph-Level Pooling

Convert node embeddings → single graph embedding:

**Options**:
1. **Mean Pooling**: Average all node embeddings (simple, permutation-invariant)
2. **Max Pooling**: Take maximum across nodes (captures extreme features)
3. **Root Node**: Use root node embedding (captures global structure)
4. **Attention Pooling** (Recommended): Learn which nodes are most important

```python
# Attention pooling:
gate = tanh(Linear(node_embeddings))  # [N, 128]
attention = softmax(Linear(gate))     # [N, 1]
graph_embedding = sum(attention * node_embeddings)  # [128]
```

**Result**: Each of 1000 gene trees → One embedding vector [128 dimensions]

---

## Part 2: Multi-Tree Aggregation (1000 embeddings → 1 embedding)

### Challenge

We have 1000 gene tree embeddings, each [128]. Need to combine into single [128] embedding that captures:
- Consensus topology across gene trees
- Which gene trees are most informative
- Conflicts/agreements between trees

### Aggregation Methods

#### Option 1: Simple Mean Pooling
```python
aggregated = mean([emb_1, emb_2, ..., emb_1000])  # [128]
```
**Pros**: Simple, permutation-invariant
**Cons**: Treats all trees equally, loses information about tree importance

#### Option 2: Attention-Based Aggregation (Recommended)
```python
# Learn attention weights for each tree
attention_weights = softmax(MLP(tree_embeddings))  # [1000, 1]
aggregated = sum(attention_weights * tree_embeddings)  # [128]
```
**Pros**: Learns which trees are most informative, flexible
**Cons**: Requires learning additional parameters

#### Option 3: Set Transformer
```python
# Permutation-invariant transformer
aggregated = SetTransformer(tree_embeddings)  # [128]
```
**Pros**: State-of-the-art for set aggregation, handles interactions
**Cons**: More complex, slower training

#### Option 4: Hierarchical Aggregation
```python
# First aggregate in groups, then aggregate groups
group_1 = mean(embeddings[0:100])
group_2 = mean(embeddings[100:200])
...
group_10 = mean(embeddings[900:1000])
final = attention_pool([group_1, ..., group_10])
```
**Pros**: Can capture different tree "clusters"
**Cons**: More complex, requires grouping strategy

**Recommended**: **Attention-based aggregation** or **Set Transformer**

---

## Part 3: MUL Tree Decoding (Embedding → MUL Tree)

### Challenge

Convert single embedding [128] → MUL tree structure:
- Tree topology (branching pattern)
- Copy numbers (how many times each species appears)

### Decoding Approaches

#### Approach A: Autoregressive Tree Generation

Generate tree node-by-node in a sequence:

```python
# Start with root
current_node = root
embedding = aggregated_embedding

while not complete:
    # Predict: split or stop?
    action = MLP(embedding, current_node_state)  # [split, stop]
    
    if action == "split":
        # Predict: which species/clade to add?
        left_child = predict_child(embedding, current_node)
        right_child = predict_child(embedding, current_node)
        # Recursively continue
    else:
        # Leaf node: predict species and copy number
        species = predict_species(embedding, current_node)
        copy_number = predict_copy_number(embedding, species)
```

**Pros**: Flexible, can generate any tree structure
**Cons**: Sequential generation is slow, error propagation

#### Approach B: Graph Generation Network

Generate tree as a graph directly:

```python
# Step 1: Predict all node types and species
node_predictions = MLP(embedding)  # [max_nodes, num_species + 1]
# For each potential node: which species? (or internal)

# Step 2: Predict edge connectivity
edge_predictions = MLP(embedding, node_embeddings)  # [max_nodes, max_nodes]
# Probability of edge between each pair

# Step 3: Construct tree from predictions
tree = construct_tree(node_predictions, edge_predictions)
```

**Pros**: Parallel generation, faster
**Cons**: Fixed maximum tree size, harder to ensure valid tree structure

#### Approach C: Template-Based Generation (Recommended)

Use a template tree structure and predict modifications:

```python
# Step 1: Predict species set and copy numbers
species_copy_numbers = MLP(embedding)  # [num_species, max_copies]
# For each species: how many copies?

# Step 2: Use existing method (e.g., ASTRAL) to get base topology
base_topology = ASTRAL(gene_trees)  # Species tree

# Step 3: Predict where to add duplicate leaves
duplication_points = MLP(embedding, base_topology)  # Which nodes get duplicates?

# Step 4: Construct MUL tree
mul_tree = add_duplicates(base_topology, species_copy_numbers, duplication_points)
```

**Pros**: Leverages existing methods, ensures valid tree structure
**Cons**: Depends on base topology method

#### Approach D: Direct Newick String Prediction

Treat as sequence-to-sequence problem:

```python
# Predict Newick string directly
newick_tokens = TransformerDecoder(aggregated_embedding)  
# Output: "((A,A,B),C);"
```

**Pros**: Simple, direct
**Cons**: Hard to ensure valid tree syntax, doesn't leverage tree structure

**Recommended**: **Approach C (Template-Based)** or **Approach A (Autoregressive)** with careful design

---

## Training Strategy

### Loss Functions

#### 1. Topology Loss

Compare predicted MUL tree topology with ground truth:

**Option A: Edit Distance on MUL-trees**
```python
loss_topology = edit_distance_multree(predicted_tree, ground_truth_tree)
```
- Uses existing `edit_distance_multree` function from codebase
- Measures structural differences
- Normalized to [0, 1]

**Option B: Robinson-Foulds Distance**
```python
loss_topology = rf_distance(predicted_tree, ground_truth_tree)
```
- Uses existing `rf_distance` function
- Measures bipartition differences
- Normalized to [0, 1]

**Option C: Tree Kernel Distance**
```python
# Use tree kernels (e.g., subtree kernel) to measure similarity
loss_topology = 1 - tree_kernel_similarity(predicted, ground_truth)
```

#### 2. Copy Number Loss

Compare predicted copy numbers with ground truth:

```python
# For each species:
predicted_copies = predicted_tree.get_copy_numbers()  # Dict: {species: count}
ground_truth_copies = ground_truth_tree.get_copy_numbers()

# Mean Squared Error
loss_copy = MSE(predicted_copies, ground_truth_copies)

# Or Cross-Entropy (if treating as classification)
loss_copy = CrossEntropy(predicted_copies, ground_truth_copies)
```

#### 3. Combined Loss

```python
total_loss = α * loss_topology + β * loss_copy
# α = 0.7, β = 0.3 (tune hyperparameters)
```

### Training Procedure

```python
# For each epoch:
for network in training_networks:
    for replicate in [1, 2, 3, 4, 5]:
        # Load 1000 gene trees
        gene_trees = load_gene_trees(network, replicate)  # List of 1000 Newick strings
        
        # Load ground truth MUL tree
        ground_truth = load_ground_truth(network)  # Newick MUL tree
        
        # Forward pass
        tree_embeddings = []
        for tree in gene_trees:
            graph = tree_to_graph(tree)
            embedding = tree_encoder(graph)  # [128]
            tree_embeddings.append(embedding)
        
        tree_embeddings = stack(tree_embeddings)  # [1000, 128]
        aggregated = multi_tree_aggregator(tree_embeddings)  # [128]
        predicted_tree = mul_tree_decoder(aggregated)  # MUL tree
        
        # Compute loss
        loss = compute_loss(predicted_tree, ground_truth)
        
        # Backward pass
        loss.backward()
        optimizer.step()
```

### Data Split

- **Training**: 16 networks × 5 replicates = 80 examples
- **Validation**: 3 networks × 5 replicates = 15 examples
- **Test**: 2 networks × 5 replicates = 10 examples

Or use **cross-validation** across networks.

### Optimization

- **Optimizer**: Adam or AdamW
- **Learning Rate**: 1e-4 to 1e-3 (with learning rate scheduling)
- **Batch Size**: 1 (one network replicate at a time) or small batches
- **Regularization**: Dropout (0.1-0.3), Weight decay (1e-5)
- **Early Stopping**: Monitor validation loss

---

## Architecture Details

### Complete Model Architecture

```python
class PolyploidyInferenceModel(nn.Module):
    def __init__(self):
        # Tree encoder (shared across all 1000 trees)
        self.tree_encoder = GeneTreeEncoder(
            d_node=38,
            d_hidden=128,
            num_layers=3,
            num_heads=4
        )
        
        # Multi-tree aggregator
        self.aggregator = MultiTreeAggregator(
            d_input=128,
            d_output=128,
            method='attention'  # or 'set_transformer'
        )
        
        # MUL tree decoder
        self.decoder = MULTreeDecoder(
            d_input=128,
            max_species=100,
            max_copies=10
        )
    
    def forward(self, gene_trees):
        # gene_trees: List of 1000 Newick strings
        
        # Encode each tree
        embeddings = []
        for tree in gene_trees:
            graph = tree_to_graph(tree)
            emb = self.tree_encoder(graph)  # [128]
            embeddings.append(emb)
        
        embeddings = torch.stack(embeddings)  # [1000, 128]
        
        # Aggregate
        aggregated = self.aggregator(embeddings)  # [128]
        
        # Decode
        mul_tree = self.decoder(aggregated)  # MUL tree structure
        
        return mul_tree
```

### GeneTreeEncoder

```python
class GeneTreeEncoder(nn.Module):
    def __init__(self, d_node, d_hidden, num_layers, num_heads):
        # Species embedding
        self.species_embedding = nn.Embedding(vocab_size, 32)
        
        # Node feature projection
        self.node_proj = nn.Linear(d_node, d_hidden)
        
        # GAT layers
        self.gat_layers = nn.ModuleList([
            GATConv(d_hidden, d_hidden // num_heads, heads=num_heads, edge_dim=5)
            for _ in range(num_layers)
        ])
        
        # Graph pooling
        self.pool_gate = nn.Sequential(nn.Linear(d_hidden, d_hidden), nn.Tanh())
        self.pool_lin = nn.Linear(d_hidden, 1)
    
    def forward(self, graph_data):
        x, edge_index, edge_attr = graph_data.x, graph_data.edge_index, graph_data.edge_attr
        
        # Project
        x = self.node_proj(x)
        
        # GAT layers
        for gat in self.gat_layers:
            x = gat(x, edge_index, edge_attr)
            x = F.relu(x)
        
        # Pool
        gate = self.pool_gate(x)
        attention = F.softmax(self.pool_lin(gate), dim=0)
        graph_emb = (attention * x).sum(dim=0)
        
        return graph_emb
```

### MultiTreeAggregator

```python
class MultiTreeAggregator(nn.Module):
    def __init__(self, d_input, d_output, method='attention'):
        self.method = method
        if method == 'attention':
            self.attention = nn.Sequential(
                nn.Linear(d_input, d_input),
                nn.Tanh(),
                nn.Linear(d_input, 1)
            )
        elif method == 'set_transformer':
            from set_transformer import SetTransformer
            self.set_transformer = SetTransformer(d_input, d_output)
    
    def forward(self, tree_embeddings):
        # tree_embeddings: [num_trees, d_input]
        if self.method == 'attention':
            weights = F.softmax(self.attention(tree_embeddings).squeeze(-1), dim=0)
            return (weights.unsqueeze(-1) * tree_embeddings).sum(dim=0)
        elif self.method == 'set_transformer':
            return self.set_transformer(tree_embeddings)
```

### MULTreeDecoder

```python
class MULTreeDecoder(nn.Module):
    def __init__(self, d_input, max_species, max_copies):
        # Predict copy numbers for each species
        self.copy_number_predictor = nn.Sequential(
            nn.Linear(d_input, 256),
            nn.ReLU(),
            nn.Linear(256, max_species * max_copies)  # [num_species, max_copies]
        )
        
        # Predict tree topology modifications
        self.topology_predictor = nn.Sequential(
            nn.Linear(d_input, 256),
            nn.ReLU(),
            nn.Linear(256, max_species * max_species)  # Adjacency matrix
        )
    
    def forward(self, embedding):
        # Predict copy numbers
        copy_logits = self.copy_number_predictor(embedding)
        copy_numbers = F.softmax(copy_logits.view(-1, max_copies), dim=-1)
        
        # Predict topology (or use template-based approach)
        topology = self.topology_predictor(embedding)
        
        # Construct MUL tree
        mul_tree = construct_mul_tree(copy_numbers, topology)
        
        return mul_tree
```

---

## What Happens with 1000 Gene Trees?

### Step-by-Step Process

1. **Input**: 1000 Newick gene tree strings
   - Each tree has different topology (due to ILS, gene duplication/loss)
   - All trees share the same set of species (leaves)

2. **Encoding Phase** (Parallel):
   - Convert each tree → graph (nodes, edges, features)
   - Pass through GNN encoder → embedding [128]
   - Result: 1000 embeddings, each [128]

3. **Aggregation Phase**:
   - Stack embeddings: [1000, 128]
   - Apply attention/aggregation → [128]
   - This single embedding captures:
     - Consensus topology across 1000 trees
     - Which trees agree/disagree
     - Overall species relationships

4. **Decoding Phase**:
   - Take aggregated embedding [128]
   - Predict:
     - Which species are polyploids (copy numbers > 1)
     - Where duplicates should be placed in tree
     - Final MUL tree structure

5. **Output**: Single MUL tree (Newick format)
   - Has duplicate leaf labels for polyploid species
   - Represents the polyploidy network

### Why This Works

- **Gene trees are noisy**: Each gene tree has its own evolutionary history
- **Consensus is key**: Aggregating 1000 trees finds the common signal
- **MUL tree captures polyploidy**: Duplicate labels indicate polyploid species
- **GNNs learn patterns**: Model learns which tree features predict correct MUL tree

---

## Task Classification

**Primary Task**: **Graph-to-Graph Generation / Structured Prediction**

**Sub-tasks**:
1. **Multi-instance Learning**: 1000 instances (gene trees) → 1 label (MUL tree)
2. **Set-to-Graph**: Set of graphs → Single graph
3. **Tree Reconstruction**: Reconstruct species tree from gene trees (with polyploidy)

**Related Tasks**:
- Species tree inference (ASTRAL, MP-EST)
- Phylogenetic network inference (existing methods)
- Graph generation (GraphRNN, GraphVAE)

**Novel Aspects**:
- Multiple input graphs (not just one)
- Output is MUL tree (not standard tree)
- Supervised learning (has ground truth)

---

## Key Challenges & Solutions

### Challenge 1: Variable Tree Sizes
**Problem**: Different gene trees have different numbers of nodes
**Solution**: PyTorch Geometric handles variable-sized graphs with batching

### Challenge 2: 1000 Trees is Large
**Problem**: Processing 1000 trees is computationally expensive
**Solution**: 
- Encode trees in parallel (batch processing)
- Use efficient GNN implementations
- Consider sampling subset of trees during training

### Challenge 3: Tree Structure Generation
**Problem**: Generating valid tree structures is hard
**Solution**: 
- Template-based approach (use ASTRAL for base topology)
- Autoregressive generation with constraints
- Post-processing to ensure valid tree

### Challenge 4: Limited Training Data
**Problem**: Only 21 networks × 5 replicates = 105 examples
**Solution**:
- Data augmentation (perturb gene trees, add noise)
- Transfer learning (pre-train on larger tree datasets)
- Regularization (dropout, weight decay)
- Cross-validation

### Challenge 5: Evaluation
**Problem**: How to measure success?
**Solution**: Use existing metrics:
- `edit_distance_multree`: Structural similarity
- `rf_distance`: Topological similarity
- Copy number accuracy: How many species have correct copy numbers?

---

## Implementation Considerations

### Computational Requirements

- **GPU**: Recommended (GNNs are GPU-friendly)
- **Memory**: ~1000 trees × ~50 nodes/tree × 128 dims ≈ 25MB per batch
- **Training Time**: Hours to days depending on architecture

### Libraries

- **PyTorch**: Deep learning framework
- **PyTorch Geometric**: GNN implementations
- **ETE3**: Tree parsing and manipulation (already in codebase)
- **NetworkX**: Graph utilities (optional)

### Integration with Existing Pipeline

- Input: Gene trees from `replicate_N/1/g_*` files
- Output: `ml_method_result.tre` (MUL tree in Newick format)
- Evaluation: Use existing `ComparisonEngine` from `compute_comparisons.py`
- Add to `summary_config.yaml` methods list

---

## Summary

**Task**: Supervised graph-to-graph generation using GNNs

**Input**: 1000 gene trees → Encode each → Aggregate → Decode → Output: 1 MUL tree

**Architecture**: 
- Tree Encoder (GAT) → Graph embeddings
- Multi-Tree Aggregator (Attention/Set Transformer) → Single embedding
- MUL Tree Decoder → MUL tree structure

**Training**: Supervised learning with topology + copy number loss

**Key Innovation**: Learning to aggregate information from 1000 noisy gene trees into a single accurate MUL tree representation

