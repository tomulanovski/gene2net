# Consensus MUL Tree from Gene Trees: Approach & Research

## Your Goal

**Input**: 1000 gene trees (Newick format)  
**Output**: Single consensus MUL tree (with duplicate labels for polyploids)

This is essentially **MUL tree consensus inference** - finding the "best" MUL tree that represents the common signal across all gene trees.

---

## How Consensus Works in Phylogenetics

### Traditional Species Tree Consensus

**Standard approach** (used by ASTRAL, MP-EST, etc.):

1. **Extract bipartitions** (splits) from each gene tree
   - Each internal node creates a bipartition: splits species into two groups
   - Example: Node splits `{A, B} | {C, D}`

2. **Count frequency** of each bipartition across all gene trees
   - Bipartition `{A, B} | {C, D}` appears in 800/1000 trees → support = 0.8

3. **Build consensus tree** from most frequent bipartitions
   - Use bipartitions with high support (e.g., >50%)
   - Combine them into a single tree structure

**Why this works**: 
- Gene trees disagree due to ILS, gene duplication/loss
- Common bipartitions represent true species relationships
- Consensus tree captures the shared signal

### Example

```
Gene Tree 1: ((A,B),C);
Gene Tree 2: ((A,B),C);
Gene Tree 3: ((A,C),B);
Gene Tree 4: ((A,B),C);
Gene Tree 5: ((A,B),C);

Bipartition {A,B} | {C} appears in 4/5 trees (80% support)
→ Consensus: ((A,B),C);
```

---

## Extending to MUL Trees

### The Challenge

**MUL trees have duplicate labels**, so standard consensus methods don't work directly:

- Standard consensus: Each species appears once
- MUL tree consensus: Some species appear multiple times (polyploids)

**Example**:
```
Gene Tree 1: ((A,B),C);
Gene Tree 2: ((A,A),B);  # A appears twice (polyploid)
Gene Tree 3: ((A,B),C);
```

How do we decide:
- Is A polyploid? (appears twice in tree 2)
- How many copies? (2? 3?)
- Where in the tree should duplicates go?

### Existing Approaches

#### 1. **Two-Step Approach** (What GRAMPA does)

```
Step 1: Run ASTRAL on gene trees → Species tree (no duplicates)
Step 2: Use GRAMPA to add polyploidy information → MUL tree
```

**Pros**: Leverages proven ASTRAL method  
**Cons**: Separates topology from polyploidy inference

#### 2. **Direct MUL Tree Consensus** (Research Area)

**Limited research exists**, but some approaches:

**a) Bipartition-based with multisets**:
- Extract bipartitions from gene trees
- Allow multisets in bipartitions: `{A, A, B} | {C}` (A appears twice)
- Count frequency of multiset bipartitions
- Build MUL tree from frequent multiset bipartitions

**b) Copy number estimation**:
- Count how many times each species appears across gene trees
- Estimate copy number from frequency
- Build tree topology, then add duplicates

**c) Network-based**:
- Some methods (Polyphest, MPSUGAR) infer networks directly
- Networks can be converted to MUL trees

---

## Research Landscape

### Well-Studied: Species Tree Consensus

**Established methods** (extensive research):
- **ASTRAL** (Mirarab et al. 2014): Quartet-based consensus
- **MP-EST** (Liu et al. 2010): Maximum pseudo-likelihood
- **STAR** (Liu et al. 2009): Species tree from average ranks
- **NJst** (Liu & Yu 2011): Neighbor-joining on average distances

**Research status**: ✅ **Mature field** with many methods and benchmarks

### Emerging: MUL Tree Consensus

**Limited but growing research**:

1. **GRAMPA** (Jones et al. 2013):
   - Uses ASTRAL for species tree
   - Then infers polyploidy via optimization
   - **Not pure consensus** - uses optimization, not just frequency

2. **Polyphest** (Jones 2017):
   - Infers MUL trees directly from gene trees
   - Uses percentile-based consensus
   - **Closer to consensus approach**

3. **MPSUGAR**:
   - Maximum parsimony on gene trees
   - Can output MUL trees
   - **Optimization-based**, not pure consensus

4. **AlloppNET**:
   - Bayesian approach with MCMC
   - Uses sequence data, not just gene trees
   - **Probabilistic**, not consensus

**Research status**: ⚠️ **Emerging field** - fewer methods, less benchmarking

### Machine Learning Approaches

**Very limited research** on ML for MUL tree consensus:

- Most ML work focuses on **species tree inference** (not MUL trees)
- Some work on **phylogenetic network inference** (related but different)
- **No published ML methods** specifically for MUL tree consensus from gene trees

**Your approach would be novel** if it:
- Uses ML/GNNs for MUL tree consensus
- Learns directly from gene trees to MUL trees
- Uses supervised learning (you have ground truth)

---

## Practical Approaches for Your Method

### Option 1: Learn Consensus Directly (Recommended)

**Idea**: Model learns to find consensus MUL tree from gene trees

```
1000 Gene Trees
    ↓
[GNN Encoder] → Extract features from each tree
    ↓
[Consensus Aggregator] → Find common patterns
    ↓
[MUL Tree Decoder] → Output consensus MUL tree
```

**How consensus is learned**:
- GNNs learn which tree features are most informative
- Aggregator learns which trees agree/disagree
- Decoder learns how to combine information into MUL tree
- **Supervised learning** teaches the model what "good consensus" looks like

**Advantages**:
- Learns optimal consensus strategy (not just frequency-based)
- Can handle complex patterns ML methods might miss
- End-to-end learning

### Option 2: Hybrid Approach

**Idea**: Combine traditional consensus with ML

```
Step 1: Use ASTRAL to get base species tree topology
Step 2: Use ML to predict:
    - Which species are polyploids (copy numbers)
    - Where to place duplicates in the tree
```

**Advantages**:
- Leverages proven ASTRAL method
- ML only handles the "hard part" (polyploidy)
- More likely to produce valid trees

### Option 3: Frequency-Based with ML Refinement

**Idea**: Use traditional consensus, refine with ML

```
Step 1: Extract bipartitions from gene trees (with multisets)
Step 2: Count frequencies, build initial MUL tree
Step 3: Use ML to refine/improve the consensus
```

**Advantages**:
- Starts with reasonable consensus
- ML adds learned improvements
- Interpretable (can see what ML changed)

---

## How Your GNN Approach Fits

### Your Method = Learned Consensus

Your GNN approach is essentially **learning how to do consensus**:

1. **Encoding phase**: Extract information from each gene tree
   - Similar to extracting bipartitions, but learned features

2. **Aggregation phase**: Find consensus across trees
   - Similar to counting frequencies, but learned weighting

3. **Decoding phase**: Build MUL tree from consensus
   - Similar to combining bipartitions, but learned construction

### Why This Could Work Well

**Advantages over traditional consensus**:

1. **Learns optimal strategy**: Not limited to frequency-based consensus
2. **Handles complex patterns**: Can learn non-obvious relationships
3. **Supervised learning**: Uses ground truth to learn what "good" means
4. **Flexible**: Can adapt to different network types

**Challenges**:

1. **Limited training data**: 21 networks × 5 replicates = 105 examples
2. **Complex output space**: MUL trees are complex structures
3. **Computational cost**: Processing 1000 trees is expensive

---

## Research Quality Assessment

### Is This Good Research?

**Yes, for several reasons**:

1. **Novel application**: ML/GNNs for MUL tree consensus is underexplored
2. **Practical problem**: Real-world need (polyploidy inference)
3. **Clear evaluation**: You have ground truth and existing methods to compare
4. **Reasonable approach**: GNNs are well-suited for tree-structured data

### Potential Contributions

1. **First ML method** for MUL tree consensus from gene trees
2. **Learned consensus** vs. frequency-based (could be better)
3. **End-to-end learning** (could capture complex patterns)
4. **Benchmarking** (compare with existing methods)

### Publication Potential

**Good potential** if:
- Method works reasonably well (doesn't need to beat all methods)
- Provides insights (which trees are most informative, etc.)
- Well-evaluated (comprehensive comparison with existing methods)

**Journals**:
- Bioinformatics
- Systematic Biology
- Molecular Biology and Evolution
- Journal of Computational Biology

---

## Recommended Approach

### Start with Hybrid Method

**Phase 1**: Use ASTRAL + ML for copy numbers
```
ASTRAL(gene_trees) → Species tree topology
ML(gene_trees) → Copy numbers per species
Combine → MUL tree
```

**Why start here**:
- Easier to implement
- More likely to work
- Validates the ML approach

**Phase 2**: Full end-to-end learning
```
ML(gene_trees) → MUL tree (topology + copy numbers)
```

**Why move here**:
- More novel
- Potentially better performance
- True learned consensus

---

## Summary

**Your goal**: Consensus MUL tree from gene trees

**Research status**:
- ✅ Species tree consensus: Well-studied
- ⚠️ MUL tree consensus: Emerging field
- ❌ ML for MUL tree consensus: Very limited (your opportunity!)

**Your approach**: 
- Novel application of GNNs to MUL tree consensus
- Supervised learning (you have ground truth)
- Could learn better consensus than frequency-based methods

**Recommendation**: 
- Start with hybrid approach (ASTRAL + ML)
- Move to full end-to-end if successful
- Good research potential with practical value

