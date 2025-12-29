# PhyloMUL: Machine Learning Method for Polyploid Species Tree Inference

## Executive Summary

**Goal**: Given N gene trees (with multi-labels from polyploidy), output a consensus MUL tree.

**Core Innovation**: Learn to weight gene trees by informativeness using attention, rather than treating all trees equally.

**Comparison Methods**: GRAMPA, Polyphest, PADRE, MPAllopp/MPSUGAR, AlloppNET

---

## Part 1: Theoretical Foundation

### 1.1 Problem Formulation

**Input**: Set of N gene trees G = {g₁, g₂, ..., gₙ}
- Each tree may have species appearing multiple times (polyploidy)
- Typical N = 1000 trees

**Output**: MUL tree M* that best explains the observed gene trees
- Topology showing evolutionary relationships
- Copy numbers for each species (k=1 for diploids, k≥2 for polyploids)

### 1.2 Why Learning Should Outperform Existing Methods

#### Justification 1: Optimal Weighting (Information Theory)

**Current methods**: Treat all gene trees equally
```
M* = argmax_M Σᵢ f(M, gᵢ)
```

**Problem**: Gene trees have varying informativeness:
- Trees from genes under strong selection → more reliable signal
- Trees with high ILS → noisy, conflicting
- Trees with gene loss → missing taxa, incomplete information
- Trees with short internal branches → uncertain relationships

**Theorem (informal)**: The optimal estimator weights observations by their informativeness. Uniform weighting is optimal ONLY when all observations are equally informative.

**Our approach**: Learn optimal weights
```
M* = argmax_M Σᵢ wᵢ · f_θ(M, gᵢ)
where wᵢ = α(gᵢ, G) is learned
```

#### Justification 2: Feature Learning (Representation Learning)

**Current methods use hand-crafted features**:
- GRAMPA: Parsimony reconciliation cost
- Polyphest: Bipartition frequencies with fixed thresholds
- PADRE: Reconciliation events
- MPSUGAR: Bayesian scoring with fixed likelihood model

These capture known patterns but may miss complex relationships.

**Our approach**: Learn what features matter from data
- GNN learns tree representations
- Attention learns which features indicate reliability
- End-to-end optimization discovers patterns humans didn't code

#### Justification 3: Joint Optimization (End-to-End Learning)

**Current methods separate**:
1. Step 1: Infer topology
2. Step 2: Infer copy numbers

**Problem**: Errors in step 1 propagate to step 2. Decisions are made without full information.

**Our approach**: End-to-end learning
- Topology and copy number predictions inform each other
- Shared representation captures both aspects
- Joint optimization finds globally better solutions

### 1.3 Mathematical Framework

Let:
- G = {g₁, ..., gₙ} be the input gene trees
- M be the target MUL tree
- φ(gᵢ) be a learned representation of tree gᵢ
- α(gᵢ, G) be a learned attention weight

**Our model**:
```
z = Σᵢ α(gᵢ, G) · φ(gᵢ)          # Weighted aggregation
M = decode(z)                      # Decode to MUL tree
```

**Attention weights** satisfy:
- α(gᵢ, G) ≥ 0 (non-negative)
- Σᵢ α(gᵢ, G) = 1 (sum to 1)
- Computed via softmax over learned scores

**Deep Sets guarantee**: Any permutation-invariant function can be represented in this form (Zaheer et al., 2017). Our architecture is a universal approximator for set functions.

### 1.4 Why This Architecture? (Design Rationale)

#### Why Not Merge All Trees Into One Graph?

An alternative approach might be:
> "Merge all 1000 gene trees into one big graph → do edge classification"

**This doesn't work well because:**

**1. Gene trees are INDEPENDENT SAMPLES, not one structure**

Each gene tree is a separate "observation" of the true species tree with noise:
```
True MUL Tree → Gene Tree 1 (with ILS noise)
             → Gene Tree 2 (with gene loss)
             → Gene Tree 3 (with duplication noise)
             → ... 1000 independent observations
```

If you merge them into one graph, you **lose the sample structure**. Which edges came from which tree? You can't tell.

**2. Consensus requires COUNTING agreement**

The key insight of consensus methods:
> "If 800 out of 1000 trees have the split {A,B}|{C,D}, it's probably real"

Edge classification on a merged graph doesn't naturally count "how many trees agree."

**3. Copy numbers aren't edges**

"Species D has 2 copies" is NOT an edge question. It's a higher-level inference about polyploidy patterns across all trees.

**4. Permutation invariance is required**

The order of input trees shouldn't matter. Our Set Transformer guarantees this mathematically. A merged graph doesn't have this property naturally.

#### Why Each Component Is Needed

| Challenge | Component | Solution |
|-----------|-----------|----------|
| Trees are complex variable-size structures | **GNN Encoder** | Compress each tree to fixed-size vector |
| Some trees are noisy/unreliable | **Set Transformer** | Learn attention weights to trust good trees more |
| Need to output structured tree | **Bipartition Predictor** | Predict which splits exist |
| Need to identify polyploids | **Copy Number Predictor** | Predict copies per species |

#### What Each Component Learns (Training Signal)

The entire model is trained **end-to-end**. The loss (bipartition prediction + copy number prediction) backpropagates through ALL components. Each component learns what helps minimize the final loss.

**GNN Encoder learns:**
- Which topology patterns are informative
- How branch lengths indicate confidence
- What it means when a species appears multiple times in a tree
- How to compress a full tree into a 256-dimensional vector

**Set Transformer (Attention) learns:**
- Which trees AGREE with each other → weight them higher
- Which trees are OUTLIERS → weight them lower
- The "consensus query" - what pattern indicates a reliable tree
- How to combine information from 1000 trees efficiently

**Bipartition Predictor learns:**
- Given the consensus representation, which splits are supported?
- How species groupings translate to tree structure

**Copy Number Predictor learns:**
- Which species tend to be polyploid based on gene tree patterns
- How the consensus embedding indicates duplication events

#### Is the Attention a Transformer?

**Yes!** The Set Transformer is a specialized Transformer architecture:

| Regular Transformer (GPT, BERT) | Set Transformer (Our Choice) |
|--------------------------------|------------------------------|
| For sequences (order matters) | For sets (order doesn't matter) |
| Uses positional encoding | No positional encoding |
| Self-attention O(N²) | Inducing points O(N·m), m << N |
| Outputs sequence | Outputs single "consensus" vector |

The Set Transformer uses the **same attention mechanism** as GPT/BERT, but modified for sets:
- No positional encodings (gene tree order is arbitrary)
- Inducing points for efficiency (can't do O(N²) with 1000 trees)
- Pooling by Multihead Attention (PMA) to extract consensus

When we say "attention weights the trees," it's literally **transformer attention** - the same mathematical operation, just in the Set Transformer variant designed for unordered sets.

---

## Part 2: Architecture Overview

```
┌──────────────────────────────────────────────────────────────────────────┐
│                              PHYLOMUL                                    │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│   INPUT: N Gene Trees (Newick format)                                    │
│       │                                                                  │
│       ▼                                                                  │
│   ┌────────────────────────────────────────────────────────────┐         │
│   │                    TREE ENCODER                            │         │
│   │  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    │         │
│   │  │ Parse Tree  │───▶│  Extract    │───▶│    GNN      │    │         │
│   │  │ to Graph    │    │  Features   │    │  (3 layers) │    │         │
│   │  └─────────────┘    └─────────────┘    └─────────────┘    │         │
│   │                                              │             │         │
│   │                                              ▼             │         │
│   │                                    ┌─────────────────┐     │         │
│   │                                    │  Graph Pooling  │     │         │
│   │                                    │  → Tree Vector  │     │         │
│   │                                    │    (256-dim)    │     │         │
│   │                                    └─────────────────┘     │         │
│   └────────────────────────────────────────────────────────────┘         │
│       │                                                                  │
│       │  N tree vectors: h₁, h₂, ..., hₙ                                 │
│       ▼                                                                  │
│   ┌────────────────────────────────────────────────────────────┐         │
│   │                  SET TRANSFORMER                           │         │
│   │                                                            │         │
│   │  ┌──────────────────────────────────────────────────┐     │         │
│   │  │  Induced Set Attention Block (ISAB)              │     │         │
│   │  │  - Reduces O(N²) to O(N·m) via inducing points   │     │         │
│   │  │  - Self-attention over tree embeddings           │     │         │
│   │  └──────────────────────────────────────────────────┘     │         │
│   │                          │                                 │         │
│   │                          ▼                                 │         │
│   │  ┌──────────────────────────────────────────────────┐     │         │
│   │  │  Pooling by Multihead Attention (PMA)            │     │         │
│   │  │  - Learns "consensus query"                      │     │         │
│   │  │  - Outputs attention weights w₁, ..., wₙ         │     │         │
│   │  │  - Produces combined vector z                    │     │         │
│   │  └──────────────────────────────────────────────────┘     │         │
│   │                                                            │         │
│   └────────────────────────────────────────────────────────────┘         │
│       │                                                                  │
│       │  Combined representation z (256-dim)                             │
│       │  Attention weights w₁, ..., wₙ (for interpretability)            │
│       ▼                                                                  │
│   ┌────────────────────────────────────────────────────────────┐         │
│   │                      DECODER                               │         │
│   │                                                            │         │
│   │    ┌─────────────────────┐    ┌─────────────────────┐     │         │
│   │    │   Bipartition Head  │    │  Copy Number Head   │     │         │
│   │    │   P(split s ∈ M)    │    │  P(k copies of sp)  │     │         │
│   │    └──────────┬──────────┘    └──────────┬──────────┘     │         │
│   │               │                          │                 │         │
│   │               └──────────┬───────────────┘                 │         │
│   │                          ▼                                 │         │
│   │              ┌─────────────────────┐                       │         │
│   │              │  Tree Reconstruction │                       │         │
│   │              │  (greedy algorithm)  │                       │         │
│   │              └─────────────────────┘                       │         │
│   │                                                            │         │
│   └────────────────────────────────────────────────────────────┘         │
│       │                                                                  │
│       ▼                                                                  │
│   OUTPUT: MUL Tree + Attention Weights (interpretability)                │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## Part 3: Tree Encoder (Detailed)

### 3.1 Tree Parsing

Convert Newick string to graph representation:
- **Nodes**: All internal nodes + leaves
- **Edges**: Parent-child relationships (bidirectional for GNN message passing)

```python
def parse_newick_to_graph(newick_string):
    """
    Convert Newick to PyTorch Geometric Data object.
    
    Example:
    "((A:0.1,B:0.2):0.3,(C:0.1,(D:0.05,D:0.05):0.1):0.2);"
    
    Returns:
    - node_features: [num_nodes, feature_dim]
    - edge_index: [2, num_edges]
    """
    tree = parse_newick(newick_string)
    nodes = list(tree.traverse())
    
    # Build adjacency
    edges = []
    for node in nodes:
        if node.up:  # has parent
            edges.append((node.id, node.up.id))
            edges.append((node.up.id, node.id))  # bidirectional
    
    edge_index = torch.tensor(edges).T
    
    # Extract features
    node_features = extract_node_features(tree, nodes)
    
    return Data(x=node_features, edge_index=edge_index)
```

### 3.2 Node Features (Complete Set)

**Category 1: Basic Node Properties**

| Feature | Dim | Description | Computation |
|---------|-----|-------------|-------------|
| `species_embedding` | 32 | Learned embedding per species | `nn.Embedding(num_species, 32)` for leaves; zeros for internal |
| `is_leaf` | 1 | Binary indicator | 1 if no children, else 0 |
| `branch_length` | 1 | Edge length to parent | Normalized by tree height |
| `depth` | 1 | Depth from root | depth / max_depth in tree |
| `distance_from_root` | 1 | Sum of branch lengths to root | Normalized by tree height |

**Category 2: Subtree Properties**

| Feature | Dim | Description | Computation |
|---------|-----|-------------|-------------|
| `num_descendants` | 1 | Number of leaves below | log(count + 1) for stability |
| `subtree_height` | 1 | Max distance to any descendant leaf | Normalized |
| `num_children` | 1 | Direct children count | Detects polytomies (>2) |
| `subtree_balance` | 1 | Balance of subtree | 1 - abs(left_size - right_size) / total_size |

**Category 3: Polyploidy-Specific Features** (Crucial for this task)

| Feature | Dim | Description | Computation |
|---------|-----|-------------|-------------|
| `species_copy_count` | 1 | Times this species appears in tree | Count all leaves with same species |
| `is_duplicated_species` | 1 | Does this species appear >1 time? | Binary |
| `copies_are_monophyletic` | 1 | Are all copies in one clade? | Check if MRCA of copies contains only this species |
| `sister_is_same_species` | 1 | Sister leaf is same species? | Suggests recent duplication |
| `mrca_depth_of_copies` | 1 | Depth of MRCA of all copies | Normalized; indicates when duplication happened |

**Category 4: Tree-Level Context** (Same for all nodes in a tree)

| Feature | Dim | Description | Computation |
|---------|-----|-------------|-------------|
| `total_taxa` | 1 | Number of leaves in tree | log(count) |
| `num_unique_species` | 1 | Unique species count | Indicates completeness |
| `num_duplicated_species` | 1 | Species appearing >1 time | Overall polyploidy signal |
| `tree_height` | 1 | Root to deepest leaf distance | Normalization reference |
| `mean_branch_length` | 1 | Average branch length | Tree-level signal quality |

**Total: 32 + 17 = 49 dimensions per node**

### 3.3 GNN Architecture

```python
class TreeEncoder(nn.Module):
    """
    Encode a gene tree into a fixed-size vector using Graph Attention Network.
    """
    def __init__(
        self,
        node_input_dim=49,
        hidden_dim=128,
        output_dim=256,
        num_layers=3,
        num_heads=4,
        dropout=0.1
    ):
        super().__init__()
        
        # Project node features to hidden dimension
        self.input_projection = nn.Sequential(
            nn.Linear(node_input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Graph Attention layers
        self.gnn_layers = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        
        for _ in range(num_layers):
            self.gnn_layers.append(
                GATConv(
                    in_channels=hidden_dim,
                    out_channels=hidden_dim // num_heads,
                    heads=num_heads,
                    dropout=dropout,
                    concat=True  # Output: hidden_dim
                )
            )
            self.layer_norms.append(nn.LayerNorm(hidden_dim))
        
        # Output projection
        self.output_projection = nn.Sequential(
            nn.Linear(hidden_dim, output_dim),
            nn.LayerNorm(output_dim)
        )
        
        # Pooling: attention-weighted mean
        self.pool_attention = nn.Sequential(
            nn.Linear(hidden_dim, 1),
            nn.Softmax(dim=0)
        )
    
    def forward(self, data):
        """
        Args:
            data: PyTorch Geometric Data object with:
                - x: node features [num_nodes, node_input_dim]
                - edge_index: adjacency [2, num_edges]
                - batch: batch assignment [num_nodes]
        
        Returns:
            tree_embedding: [batch_size, output_dim]
        """
        x = self.input_projection(data.x)
        
        # Message passing with residual connections
        for gnn, norm in zip(self.gnn_layers, self.layer_norms):
            x_new = gnn(x, data.edge_index)
            x_new = F.elu(x_new)
            x = norm(x + x_new)  # Residual connection
        
        # Attention-weighted pooling per graph
        attn_weights = self.pool_attention(x)
        x_pooled = scatter_sum(attn_weights * x, data.batch, dim=0)
        
        return self.output_projection(x_pooled)
```

**Design choices explained**:
- **GATConv**: Graph Attention lets model learn which neighbors matter more (parent vs children)
- **3 layers**: Sufficient for information to flow across typical tree depth
- **Residual connections**: Stabilize training, prevent vanishing gradients
- **Attention pooling**: Better than mean pooling; learns which nodes are most informative

---

## Part 4: Set Transformer Aggregator (Detailed)

### 4.1 Why Set Transformer?

**Problem**: Given N=1000 tree embeddings, compute weighted combination.

**Simple attention** (what I wrote before): O(N) but limited expressiveness

**Full self-attention**: O(N²) - too expensive for N=1000

**Set Transformer** (Lee et al., 2019): O(N·m) where m << N
- Uses "inducing points" to summarize the set
- Full expressiveness with computational efficiency
- Proven universal approximator for set functions

#### What Are Inducing Points? (Plain English)

**The Problem**: Standard self-attention compares every element to every other element:
```
1000 trees × 1000 trees = 1,000,000 attention computations (O(N²))
```
This is too slow for N=1000 gene trees.

**The Solution**: Instead of comparing trees directly, use a small set of **learned "summary" vectors** (m=32) as intermediaries.

**Analogy - The Committee**:

Think of inducing points as "committee members" that summarize the full set:
1. Each tree "reports" to the committee (trees → inducing points)
2. The committee aggregates information from all trees
3. Each tree gets the committee's summary (inducing points → trees)

```
WITHOUT inducing points (O(N²)):
  Tree₁ ←→ Tree₂ ←→ Tree₃ ←→ ... ←→ Tree₁₀₀₀
  (Every tree talks to every other tree - expensive!)

WITH inducing points (O(N·m)):
  Step 1: All trees → [I₁, I₂, ..., I₃₂]    (N × m operations)
  Step 2: [I₁, I₂, ..., I₃₂] → All trees    (m × N operations)
  (Trees talk through committee - much cheaper!)
```

**Key insight**: The inducing points are **LEARNED** end-to-end. The model figures out what "summaries" are most useful for the task. For our problem, they might learn to represent:
- "Typical well-resolved tree"
- "Tree with lots of polyploidy signals"
- "Outlier/noisy tree pattern"

**Complexity reduction**:
- Without: O(N²) = O(1,000,000) 
- With: O(2·N·m) = O(64,000) where m=32
- **~15x faster** with same expressiveness

### 4.2 Architecture

```python
class SetTransformerAggregator(nn.Module):
    """
    Aggregate N tree embeddings into single consensus representation.
    Based on: Lee et al., "Set Transformer" (ICML 2019)
    """
    def __init__(
        self,
        input_dim=256,
        hidden_dim=256,
        num_heads=4,
        num_inducing_points=32,  # m in O(N·m)
        num_isab_blocks=2,
        num_pma_seeds=1
    ):
        super().__init__()
        
        # Induced Set Attention Blocks (ISAB)
        # Reduces O(N²) to O(N·m) using inducing points
        self.isab_blocks = nn.ModuleList([
            ISAB(
                dim_input=input_dim if i == 0 else hidden_dim,
                dim_output=hidden_dim,
                num_heads=num_heads,
                num_inducing_points=num_inducing_points
            )
            for i in range(num_isab_blocks)
        ])
        
        # Pooling by Multihead Attention (PMA)
        # Learns to extract k summary vectors (k=1 for our case)
        self.pma = PMA(
            dim_input=hidden_dim,
            dim_output=hidden_dim,
            num_heads=num_heads,
            num_seeds=num_pma_seeds
        )
        
        # Output projection
        self.output_projection = nn.Linear(hidden_dim, input_dim)
    
    def forward(self, tree_embeddings, return_attention=True):
        """
        Args:
            tree_embeddings: [batch_size, N_trees, embed_dim]
            return_attention: whether to return attention weights
        
        Returns:
            combined: [batch_size, embed_dim]
            attention_weights: [batch_size, N_trees] (if return_attention)
        """
        # Apply ISAB blocks: trees interact through inducing points
        x = tree_embeddings
        for isab in self.isab_blocks:
            x = isab(x)
        
        # PMA: pool to single vector
        combined, attention_weights = self.pma(x, return_attention=True)
        combined = combined.squeeze(1)  # [batch, 1, dim] -> [batch, dim]
        
        combined = self.output_projection(combined)
        
        if return_attention:
            return combined, attention_weights.squeeze(1)
        return combined


class ISAB(nn.Module):
    """Induced Set Attention Block"""
    def __init__(self, dim_input, dim_output, num_heads, num_inducing_points):
        super().__init__()
        self.inducing_points = nn.Parameter(
            torch.randn(1, num_inducing_points, dim_output)
        )
        self.mab1 = MAB(dim_output, dim_input, dim_output, num_heads)
        self.mab2 = MAB(dim_input, dim_output, dim_output, num_heads)
    
    def forward(self, x):
        # x: [batch, N, dim]
        batch_size = x.size(0)
        I = self.inducing_points.expand(batch_size, -1, -1)
        
        # Inducing points attend to input
        H = self.mab1(I, x)  # [batch, m, dim]
        
        # Input attends to inducing points
        return self.mab2(x, H)  # [batch, N, dim]


class PMA(nn.Module):
    """Pooling by Multihead Attention"""
    def __init__(self, dim_input, dim_output, num_heads, num_seeds):
        super().__init__()
        self.seeds = nn.Parameter(torch.randn(1, num_seeds, dim_output))
        self.mab = MAB(dim_output, dim_input, dim_output, num_heads)
    
    def forward(self, x, return_attention=False):
        batch_size = x.size(0)
        S = self.seeds.expand(batch_size, -1, -1)
        output, attention = self.mab(S, x, return_attention=True)
        
        if return_attention:
            return output, attention
        return output


class MAB(nn.Module):
    """Multihead Attention Block"""
    def __init__(self, dim_Q, dim_KV, dim_output, num_heads):
        super().__init__()
        self.attention = nn.MultiheadAttention(
            embed_dim=dim_output,
            num_heads=num_heads,
            kdim=dim_KV,
            vdim=dim_KV,
            batch_first=True
        )
        self.layer_norm = nn.LayerNorm(dim_output)
        self.ffn = nn.Sequential(
            nn.Linear(dim_output, dim_output * 4),
            nn.GELU(),
            nn.Linear(dim_output * 4, dim_output)
        )
        self.layer_norm2 = nn.LayerNorm(dim_output)
        self.query_proj = nn.Linear(dim_Q, dim_output)
    
    def forward(self, Q, KV, return_attention=False):
        Q = self.query_proj(Q)
        attn_output, attn_weights = self.attention(Q, KV, KV)
        x = self.layer_norm(Q + attn_output)
        x = self.layer_norm2(x + self.ffn(x))
        
        if return_attention:
            return x, attn_weights
        return x
```

### 4.3 Interpretation

The Set Transformer learns:
1. **ISAB blocks**: How trees relate to each other (conflicting? agreeing? complementary?)
2. **PMA**: What "consensus" looks like (learned seed vector acts as query for "ideal tree")
3. **Attention weights**: Which trees contributed most (interpretable!)

---

## Part 5: Decoder (Detailed)

### 5.1 Overview

Two prediction heads that share the aggregated representation:
1. **Bipartition Head**: Which splits appear in the output tree
2. **Copy Number Head**: How many copies of each species

### 5.2 Bipartition Representation

A bipartition splits all taxa into two groups:
- Example: {A, B, C} | {D, E} means there's an internal edge separating these groups

**Encoding**: For S unique species, a bipartition is a binary vector of length S.

**Candidate generation**: Instead of predicting all 2^S possible bipartitions:
1. Extract all bipartitions from input gene trees
2. Union these to get candidate set (typically << 2^S)
3. Only predict on candidates

```python
def extract_candidate_bipartitions(gene_trees, species_list):
    """
    Extract all bipartitions that appear in any input gene tree.
    These become the candidates for prediction.
    """
    candidates = set()
    
    for tree in gene_trees:
        for node in tree.traverse():
            if not node.is_leaf() and not node.is_root():
                # Get species on each side of this edge
                left_species = frozenset(leaf.species for leaf in node.get_leaves())
                right_species = frozenset(species_list) - left_species
                
                # Canonical form: smaller set first
                if len(left_species) <= len(right_species):
                    bipart = (left_species, right_species)
                else:
                    bipart = (right_species, left_species)
                
                candidates.add(bipart)
    
    return list(candidates)
```

### 5.3 Bipartition Predictor

```python
class BipartitionPredictor(nn.Module):
    """
    Predict probability of each candidate bipartition appearing in output tree.
    """
    def __init__(
        self,
        consensus_dim=256,
        species_embed_dim=32,
        hidden_dim=128,
        num_species=50
    ):
        super().__init__()
        
        # Learnable species embeddings
        self.species_embeddings = nn.Embedding(num_species, species_embed_dim)
        
        # Bipartition encoder: encode the two sides of the split
        self.bipart_encoder = nn.Sequential(
            nn.Linear(species_embed_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        
        # Scorer: does this bipartition match the consensus?
        self.scorer = nn.Sequential(
            nn.Linear(consensus_dim + hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)
        )
    
    def encode_bipartition(self, left_species_ids, right_species_ids):
        """Encode a bipartition as a vector."""
        # Mean pooling of species embeddings on each side
        left_emb = self.species_embeddings(left_species_ids).mean(dim=0)
        right_emb = self.species_embeddings(right_species_ids).mean(dim=0)
        
        # Concatenate (order-invariant via sorting during candidate generation)
        bipart_emb = torch.cat([left_emb, right_emb])
        return self.bipart_encoder(bipart_emb)
    
    def forward(self, consensus_embedding, candidate_bipartitions):
        """
        Args:
            consensus_embedding: [batch_size, consensus_dim]
            candidate_bipartitions: list of (left_species_ids, right_species_ids) tuples
        
        Returns:
            probabilities: [batch_size, num_candidates]
        """
        batch_size = consensus_embedding.size(0)
        num_candidates = len(candidate_bipartitions)
        
        # Encode all bipartitions
        bipart_embeddings = torch.stack([
            self.encode_bipartition(left, right)
            for left, right in candidate_bipartitions
        ])  # [num_candidates, hidden_dim]
        
        # Expand for batch
        bipart_embeddings = bipart_embeddings.unsqueeze(0).expand(batch_size, -1, -1)
        consensus_expanded = consensus_embedding.unsqueeze(1).expand(-1, num_candidates, -1)
        
        # Concatenate and score
        combined = torch.cat([consensus_expanded, bipart_embeddings], dim=-1)
        logits = self.scorer(combined).squeeze(-1)  # [batch_size, num_candidates]
        
        return torch.sigmoid(logits)
```

### 5.4 Copy Number Predictor

```python
class CopyNumberPredictor(nn.Module):
    """
    Predict number of copies for each species (1, 2, 3, ..., max_copies).
    """
    def __init__(
        self,
        consensus_dim=256,
        species_embed_dim=32,
        hidden_dim=128,
        num_species=50,
        max_copies=5
    ):
        super().__init__()
        
        self.species_embeddings = nn.Embedding(num_species, species_embed_dim)
        
        self.predictor = nn.Sequential(
            nn.Linear(consensus_dim + species_embed_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, max_copies)
        )
        
        self.max_copies = max_copies
    
    def forward(self, consensus_embedding, species_ids):
        """
        Args:
            consensus_embedding: [batch_size, consensus_dim]
            species_ids: [num_species] tensor of species indices
        
        Returns:
            copy_probs: [batch_size, num_species, max_copies]
        """
        batch_size = consensus_embedding.size(0)
        num_species = len(species_ids)
        
        # Get species embeddings
        species_embs = self.species_embeddings(species_ids)  # [num_species, embed_dim]
        
        # Expand for batch
        species_embs = species_embs.unsqueeze(0).expand(batch_size, -1, -1)
        consensus_expanded = consensus_embedding.unsqueeze(1).expand(-1, num_species, -1)
        
        # Concatenate and predict
        combined = torch.cat([consensus_expanded, species_embs], dim=-1)
        logits = self.predictor(combined)  # [batch_size, num_species, max_copies]
        
        return F.softmax(logits, dim=-1)
```

### 5.5 Tree Reconstruction

Given predicted bipartitions and copy numbers, reconstruct the MUL tree:

```python
def reconstruct_mul_tree(bipartition_probs, copy_number_probs, 
                         candidate_bipartitions, species_list,
                         threshold=0.5):
    """
    Reconstruct MUL tree from predictions.
    
    Args:
        bipartition_probs: [num_candidates] probabilities
        copy_number_probs: [num_species, max_copies] probabilities
        candidate_bipartitions: list of bipartition tuples
        species_list: list of species names
        threshold: probability threshold for including bipartition
    
    Returns:
        mul_tree: Newick string of reconstructed tree
    """
    # Step 1: Select bipartitions above threshold
    selected_biparts = [
        bipart for bipart, prob in zip(candidate_bipartitions, bipartition_probs)
        if prob > threshold
    ]
    
    # Step 2: Determine copy numbers
    copy_numbers = {
        species: probs.argmax().item() + 1  # +1 because index 0 = 1 copy
        for species, probs in zip(species_list, copy_number_probs)
    }
    
    # Step 3: Expand species with copies
    expanded_taxa = []
    for species, num_copies in copy_numbers.items():
        if num_copies == 1:
            expanded_taxa.append(species)
        else:
            for i in range(num_copies):
                expanded_taxa.append(f"{species}_copy{i+1}")
    
    # Step 4: Build tree from compatible bipartitions (greedy)
    tree = build_tree_from_bipartitions(selected_biparts, expanded_taxa)
    
    return tree.write(format=1)  # Newick string
```

---

## Part 6: Full Model

```python
class PhyloMUL(nn.Module):
    """
    Complete PhyloMUL model for MUL tree inference from gene trees.
    """
    def __init__(
        self,
        num_species=50,
        max_copies=5,
        node_feature_dim=49,
        tree_embed_dim=256,
        hidden_dim=128,
        gnn_layers=3,
        num_inducing_points=32,
        num_attention_heads=4
    ):
        super().__init__()
        
        # Tree encoder
        self.encoder = TreeEncoder(
            node_input_dim=node_feature_dim,
            hidden_dim=hidden_dim,
            output_dim=tree_embed_dim,
            num_layers=gnn_layers,
            num_heads=num_attention_heads
        )
        
        # Set transformer aggregator
        self.aggregator = SetTransformerAggregator(
            input_dim=tree_embed_dim,
            hidden_dim=tree_embed_dim,
            num_heads=num_attention_heads,
            num_inducing_points=num_inducing_points
        )
        
        # Decoder heads
        self.bipartition_predictor = BipartitionPredictor(
            consensus_dim=tree_embed_dim,
            num_species=num_species,
            hidden_dim=hidden_dim
        )
        
        self.copy_number_predictor = CopyNumberPredictor(
            consensus_dim=tree_embed_dim,
            num_species=num_species,
            hidden_dim=hidden_dim,
            max_copies=max_copies
        )
    
    def forward(self, gene_trees, candidate_bipartitions, species_ids):
        """
        Args:
            gene_trees: list of N PyG Data objects (parsed gene trees)
            candidate_bipartitions: list of bipartition tuples
            species_ids: tensor of species indices
        
        Returns:
            dict with:
                - bipartition_probs: [num_candidates]
                - copy_number_probs: [num_species, max_copies]
                - attention_weights: [N] (for interpretability)
        """
        # Encode each tree
        tree_embeddings = torch.stack([
            self.encoder(tree) for tree in gene_trees
        ])  # [N, tree_embed_dim]
        
        # Add batch dimension
        tree_embeddings = tree_embeddings.unsqueeze(0)  # [1, N, tree_embed_dim]
        
        # Aggregate with attention
        consensus, attention_weights = self.aggregator(tree_embeddings)
        
        # Decode
        bipartition_probs = self.bipartition_predictor(
            consensus, candidate_bipartitions
        )
        copy_number_probs = self.copy_number_predictor(
            consensus, species_ids
        )
        
        return {
            'bipartition_probs': bipartition_probs.squeeze(0),
            'copy_number_probs': copy_number_probs.squeeze(0),
            'attention_weights': attention_weights.squeeze(0)
        }
```

---

## Part 7: Training

### 7.1 Loss Functions

```python
class PhyloMULLoss(nn.Module):
    """Combined loss for bipartition and copy number prediction."""
    
    def __init__(self, lambda_bipart=1.0, lambda_copy=1.0, pos_weight=2.0):
        super().__init__()
        self.lambda_bipart = lambda_bipart
        self.lambda_copy = lambda_copy
        
        # Positive weight for bipartition BCE (handle class imbalance)
        self.bipart_loss = nn.BCELoss(reduction='none')
        self.pos_weight = pos_weight
    
    def forward(self, predictions, targets):
        """
        Args:
            predictions: dict from model forward pass
            targets: dict with 'bipartition_labels', 'copy_numbers'
        
        Returns:
            total_loss, dict of individual losses
        """
        # Bipartition loss (weighted BCE)
        bipart_pred = predictions['bipartition_probs']
        bipart_target = targets['bipartition_labels'].float()
        
        bce = self.bipart_loss(bipart_pred, bipart_target)
        # Weight positive examples more (fewer positive bipartitions)
        weights = torch.where(bipart_target == 1, self.pos_weight, 1.0)
        bipart_loss = (bce * weights).mean()
        
        # Copy number loss (cross-entropy)
        copy_pred = predictions['copy_number_probs']
        copy_target = targets['copy_numbers']  # [num_species] integers
        
        copy_loss = F.cross_entropy(
            copy_pred.view(-1, copy_pred.size(-1)),
            copy_target.view(-1)
        )
        
        # Combined
        total_loss = self.lambda_bipart * bipart_loss + self.lambda_copy * copy_loss
        
        return total_loss, {
            'bipartition_loss': bipart_loss.item(),
            'copy_number_loss': copy_loss.item()
        }
```

### 7.2 Training Configuration

```yaml
# config.yaml
training:
  batch_size: 16  # Number of (gene_tree_set, MUL_tree) examples per batch
  learning_rate: 1e-4
  weight_decay: 1e-5
  num_epochs: 100
  early_stopping_patience: 10
  
optimizer:
  type: AdamW
  betas: [0.9, 0.999]
  
scheduler:
  type: CosineAnnealingLR
  T_max: 100
  eta_min: 1e-6

loss:
  lambda_bipart: 1.0
  lambda_copy: 1.0
  pos_weight: 2.0

model:
  tree_embed_dim: 256
  hidden_dim: 128
  gnn_layers: 3
  num_inducing_points: 32
  num_attention_heads: 4
  max_copies: 5

data:
  num_gene_trees: 1000  # N trees per example
  train_split: 0.8
  val_split: 0.1
  test_split: 0.1
```

### 7.3 Training Loop

```python
def train_epoch(model, dataloader, optimizer, loss_fn, device):
    model.train()
    total_loss = 0
    
    for batch in dataloader:
        gene_trees = [tree.to(device) for tree in batch['gene_trees']]
        candidate_biparts = batch['candidate_bipartitions']
        species_ids = batch['species_ids'].to(device)
        targets = {k: v.to(device) for k, v in batch['targets'].items()}
        
        # Forward
        predictions = model(gene_trees, candidate_biparts, species_ids)
        
        # Loss
        loss, loss_dict = loss_fn(predictions, targets)
        
        # Backward
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        total_loss += loss.item()
    
    return total_loss / len(dataloader)
```

---

## Part 8: Evaluation

### 8.1 Metrics

```python
def evaluate_model(model, dataloader, device):
    """Comprehensive evaluation metrics."""
    model.eval()
    
    all_metrics = {
        'bipartition_precision': [],
        'bipartition_recall': [],
        'bipartition_f1': [],
        'copy_number_accuracy': [],
        'edit_distance': [],
        'robinson_foulds': []
    }
    
    with torch.no_grad():
        for batch in dataloader:
            predictions = model(...)
            
            # Bipartition metrics
            pred_biparts = (predictions['bipartition_probs'] > 0.5).float()
            true_biparts = batch['targets']['bipartition_labels']
            
            precision = (pred_biparts * true_biparts).sum() / pred_biparts.sum()
            recall = (pred_biparts * true_biparts).sum() / true_biparts.sum()
            f1 = 2 * precision * recall / (precision + recall + 1e-8)
            
            # Copy number accuracy
            pred_copies = predictions['copy_number_probs'].argmax(dim=-1)
            true_copies = batch['targets']['copy_numbers']
            copy_acc = (pred_copies == true_copies).float().mean()
            
            # Tree-level metrics (reconstruct tree first)
            pred_tree = reconstruct_mul_tree(predictions, ...)
            true_tree = batch['ground_truth_tree']
            
            edit_dist = compute_edit_distance(pred_tree, true_tree)
            rf_dist = compute_robinson_foulds(pred_tree, true_tree)
            
            # Store
            all_metrics['bipartition_f1'].append(f1.item())
            all_metrics['copy_number_accuracy'].append(copy_acc.item())
            all_metrics['edit_distance'].append(edit_dist)
            all_metrics['robinson_foulds'].append(rf_dist)
    
    return {k: np.mean(v) for k, v in all_metrics.items()}
```

### 8.2 Baseline Comparison

Compare against:
1. **GRAMPA** - your existing benchmark
2. **Polyphest** - your existing benchmark  
3. **PADRE** - your existing benchmark
4. **MPAllopp/MPSUGAR** - your existing benchmark
5. **AlloppNET** - where applicable

Use your existing evaluation pipeline in `simulations/scripts/`.

### 8.3 Ablation Studies

| Ablation | What it tests |
|----------|---------------|
| No attention (uniform weights) | Does learned weighting help? |
| No polyploidy features | Are copy-specific features important? |
| Fewer GNN layers (1, 2) | How much message passing needed? |
| Simple aggregation (mean) vs Set Transformer | Is Set Transformer worth it? |
| Fewer gene trees (100, 500) | How many trees needed? |

---

## Part 9: Timeline (8 Months)

### Phase 1: Foundation (Month 1-2)

**Month 1: Data Pipeline**
- Week 1-2: Newick parser → PyG graph with features
- Week 3-4: Dataset class, data loading, train/val/test splits
- Deliverable: Can load your existing simulation data

**Month 2: Tree Encoder**
- Week 1-2: Implement TreeEncoder with all features
- Week 3-4: Unit tests, embedding visualization
- Deliverable: Working encoder, validated on toy examples

### Phase 2: Core Model (Month 3-4)

**Month 3: Aggregator**
- Week 1-2: Implement Set Transformer
- Week 3-4: Test attention weights are meaningful
- Deliverable: Encoder + Aggregator pipeline

**Month 4: Decoder**
- Week 1-2: Bipartition predictor
- Week 3-4: Copy number predictor, tree reconstruction
- Deliverable: Full model, end-to-end trainable

### Phase 3: Training & Evaluation (Month 5-6)

**Month 5: Training**
- Week 1-2: Training loop, hyperparameter tuning
- Week 3-4: Full-scale training, model selection
- Deliverable: Trained model

**Month 6: Evaluation**
- Week 1-2: Benchmark vs GRAMPA, Polyphest, etc.
- Week 3-4: Ablation studies
- Deliverable: Complete results

### Phase 4: Analysis & Paper (Month 7-8)

**Month 7: Analysis**
- Week 1-2: Attention analysis (which trees matter?)
- Week 3-4: Error analysis, biological insights
- Deliverable: Interpretability section

**Month 8: Paper**
- Week 1-2: Write paper
- Week 3-4: Code cleanup, documentation, release
- Deliverable: Submitted paper + GitHub repo

---

## Part 10: Fallback Strategies

### If Model Doesn't Beat Baselines

Focus on **interpretability angle**:
> "Competitive accuracy with interpretable gene tree weighting"

The attention weights are novel and useful even if accuracy is similar.

### If Set Transformer is Too Complex

Fall back to **simple attention**:
```python
# Instead of Set Transformer
scores = self.scorer(tree_embeddings)  # [N, 1]
weights = F.softmax(scores, dim=0)
combined = (weights * tree_embeddings).sum(dim=0)
```

Still learns weighting, just simpler.

### If GNN Training Fails

Use **bipartition-based encoder** instead:
```python
# Extract bipartitions from tree
# Embed each bipartition
# Aggregate bipartition embeddings
```

More interpretable, easier to train.

---

## Summary

**PhyloMUL** is a principled ML approach to MUL tree inference:

1. **Theoretical foundation**: Optimal weighting, feature learning, joint optimization
2. **Architecture**: GNN encoder → Set Transformer → Two-head decoder
3. **Interpretability**: Attention weights show which trees contributed
4. **Comprehensive evaluation**: vs 5 existing methods + ablations

**Key innovation**: Learn which gene trees to trust, rather than treating all equally.

---

## Appendix: Glossary of ML Terms

Quick reference for machine learning terminology used in this document.

### Neural Network Basics

| Term | Definition |
|------|------------|
| **Embedding** | Converting discrete items (species names, nodes) into continuous vectors that capture semantic meaning |
| **Hidden dimension** | Size of internal vector representations (e.g., 256-dim means vectors with 256 numbers) |
| **Forward pass** | Running input through the model to get output |
| **Backpropagation** | Computing gradients of loss with respect to all parameters |
| **End-to-end learning** | Training the entire pipeline jointly, rather than separate steps |

### Attention and Transformers

| Term | Definition |
|------|------------|
| **Attention** | Mechanism that computes weighted combinations, where weights depend on content (not fixed) |
| **Self-attention** | Each element attends to all other elements in the same set/sequence |
| **Query, Key, Value (Q, K, V)** | Three projections used in attention: Query asks "what am I looking for?", Keys are "what do I contain?", Values are "what information do I provide?" |
| **Multi-head attention** | Running multiple attention operations in parallel, each learning different patterns |
| **Transformer** | Architecture based on self-attention, used in GPT, BERT, etc. |
| **Set Transformer** | Transformer variant for unordered sets (no positional encoding) |
| **Inducing points** | Learned summary vectors that reduce O(N²) attention to O(N·m) |

### Graph Neural Networks

| Term | Definition |
|------|------------|
| **GNN** | Neural network that operates on graph-structured data |
| **Message passing** | Nodes update their representations by aggregating information from neighbors |
| **GATConv** | Graph Attention Convolution - message passing where neighbor contributions are weighted by attention |
| **Graph pooling** | Aggregating node representations into a single graph-level vector |
| **Node features** | Input attributes for each node (e.g., species ID, branch length) |

### Training

| Term | Definition |
|------|------------|
| **Loss function** | Measures how wrong the model's predictions are (lower is better) |
| **Cross-entropy loss** | Standard loss for classification (predicting categories) |
| **Binary cross-entropy (BCE)** | Loss for yes/no predictions (e.g., "does this bipartition exist?") |
| **Optimizer** | Algorithm that updates model parameters to reduce loss (e.g., Adam, SGD) |
| **Learning rate** | How big of a step to take when updating parameters |
| **Epoch** | One complete pass through all training data |
| **Batch size** | Number of examples processed together before updating parameters |
| **Overfitting** | Model memorizes training data but fails on new data |
| **Early stopping** | Stop training when validation performance stops improving |

### Model Components

| Term | Definition |
|------|------------|
| **Encoder** | Converts input (trees) into internal representation (vectors) |
| **Decoder** | Converts internal representation into output (predictions) |
| **Linear layer** | Simple transformation: output = W·input + b |
| **ReLU / ELU / GELU** | Activation functions that add non-linearity |
| **LayerNorm** | Normalizes activations to stabilize training |
| **Dropout** | Randomly zeros some values during training to prevent overfitting |
| **Softmax** | Converts numbers into probabilities that sum to 1 |
| **Sigmoid** | Converts number into probability between 0 and 1 |

### Evaluation

| Term | Definition |
|------|------------|
| **Precision** | Of predicted positives, how many are correct? |
| **Recall** | Of actual positives, how many did we find? |
| **F1 score** | Harmonic mean of precision and recall |
| **Ablation study** | Removing components to measure their contribution |
| **Hyperparameter** | Settings chosen before training (learning rate, hidden dim, etc.) |

### Set Learning Specific

| Term | Definition |
|------|------------|
| **Permutation invariance** | Output doesn't change if input order changes |
| **Deep Sets** | Framework proving any permutation-invariant function can be written as ρ(Σφ(xᵢ)) |
| **ISAB** | Induced Set Attention Block - efficient attention using inducing points |
| **PMA** | Pooling by Multihead Attention - aggregates set into fixed-size output |
