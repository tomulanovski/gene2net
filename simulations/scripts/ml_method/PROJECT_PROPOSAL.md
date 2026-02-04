# Machine Learning for MUL-Tree Inference from Gene Trees

## A Novel Approach to Polyploid Species Tree Reconstruction

---

## Executive Summary

We propose developing a machine learning method for inferring MUL-trees (multi-labeled trees representing polyploid species relationships) from gene trees. Unlike existing methods that rely on fixed heuristics and counting rules, our approach **learns optimal decision functions from data**, enabling context-aware inference that adapts to varying levels of incomplete lineage sorting (ILS), gene tree estimation error, and polyploidy complexity.

**Key Innovation**: Learn which bipartitions to include in the species tree based on multiple features, rather than using fixed frequency thresholds.

**Timeline**: 6 months to publication-ready method

---

## 1. The Problem

### 1.1 Biological Context

Polyploidy (whole genome duplication) is a major driver of plant evolution, occurring in ~15% of angiosperm speciation events and ~35% of fern speciation events. Reconstructing evolutionary relationships among polyploid species requires **MUL-trees** (multi-labeled trees) where polyploid species appear multiple times, representing their duplicated genomes.

### 1.2 Current Methods: Two Paradigms

Existing methods fall into two algorithmic paradigms:

#### Paradigm 1: Consensus-Based (Polyphest)

**Polyphest** builds a MUL-tree by analyzing **cluster frequencies** across gene trees:
- Extract clusters (clades/monophyletic groups) from each gene tree
- Count frequency of each cluster across all trees
- Include clusters above a **fixed percentile threshold** (e.g., 50%, 70%, 90%)
- Handle polyploidy via **multisets** (species can appear multiple times)

```
Gene Trees → Extract Clusters → Filter by Frequency → Build MUL-tree
```

**Limitation**: The threshold is fixed and doesn't adapt to data quality or context.

#### Paradigm 2: Reconciliation-Based (GRAMPA, PADRE, MPSUGAR)

These methods score candidate trees by how well gene trees "reconcile" to them:

**GRAMPA** (Gene Reconciliation and Polyploidy Analyzer):
- First infers species tree using ASTRAL (quartet-based)
- Maps each gene tree to species tree via **reconciliation**
- Counts events: duplications, losses, extra lineages
- Minimizes total **parsimony cost**

**PADRE** (Parsimony-based Deduction of Reticulate Evolution):
- Uses **maximum parsimony** to infer networks
- Identifies **clusters** (clades) that should exist
- Minimizes duplications, losses, deep coalescences
- Doesn't require explicit copy number input

**MPSUGAR** (Multispecies Phylogenetic SUGAR):
- **Bayesian MCMC** approach
- Samples networks weighted by posterior probability
- Uses **multispecies coalescent** model for likelihood
- Requires taxon map (copy → species mapping)

```
Candidate MUL-tree → Reconcile Gene Trees → Score by Events → Optimize
```

**Limitation**: Fixed cost functions; all gene trees contribute equally to the score.

### 1.3 Summary of Current Methods

| Method | Paradigm | Core Mechanism | Key Limitation |
|--------|----------|----------------|----------------|
| **Polyphest** | Consensus | Cluster frequency thresholds | Fixed thresholds |
| **GRAMPA** | Reconciliation | Parsimony (dup/loss/extra lineages) | Fixed costs; treats trees equally |
| **PADRE** | Reconciliation | Parsimony (clusters + events) | Fixed costs; no weighting |
| **MPSUGAR** | Reconciliation | Bayesian coalescent likelihood | Computationally expensive |

### 1.4 The Gap We Address

**Current methods**: Use **fixed rules** (thresholds, costs) that don't adapt to data.

**Our question**: Can we **learn** optimal decision rules from data with known ground truth?

| Paradigm | Current | Our Extension |
|----------|---------|---------------|
| Consensus | Fixed frequency threshold | **Learned** threshold based on features |
| Reconciliation | Fixed event costs | **Learned** cost function (future work) |

---

## 2. Our Approach: Method Alternatives

### 2.1 Key Design Principles

**Topology-only features**: Like Polyphest, GRAMPA, and other existing methods, we use **topology** (tree structure) rather than branch lengths. Branch support values (annotations) can be used, but not branch lengths.

**Bipartitions vs Clusters**: We use bipartitions (splits) rather than clusters (clades). Bipartitions are standard in phylogenetics (ASTRAL, consensus methods) with well-defined compatibility checking. This is a representation choice, not a contribution.

### 2.2 Method Alternatives Overview

We propose three method variants, from simplest to most novel:

| Method | Description | Novelty | Risk |
|--------|-------------|---------|------|
| **Option A** | Learned bipartition scoring | Low-Medium | Low |
| **Option B** | Learned tree weighting + bipartition scoring | **High** | Low |
| **Option C** | GNN-based approach | High | High |

---

### 2.3 Option A: Learned Bipartition Scoring (Baseline)

**What it does**: Replace Polyphest's frequency threshold with a learned scoring function.

```
Polyphest:  Include cluster if frequency > percentile_threshold
Us:         Include bipartition if f(features; θ) > 0.5
```

**Features** (topology-based):
- Frequency, weighted frequency, presence variance
- Branch support statistics (mean, min, max, std)
- Conflict rate, consistency score
- Polyploidy indicators (copy ratios, monophyly)
- Structural features (split size, depth)

**Novelty**: Low-Medium (learned scoring, but similar structure to Polyphest)

**Pros**: Simple, interpretable, fast to implement
**Cons**: May not be sufficiently novel for top venues

---

### 2.4 Option B: Learned Tree Weighting (RECOMMENDED - Main Novelty)

**Key insight**: Current methods treat all gene trees equally. But gene trees vary in reliability:
- Some have low branch support → uncertain
- Some are outliers → possibly erroneous
- Some have missing taxa → incomplete

#### Bootstrap: The Foundation for Tree Reliability

**Bootstrap support** (0-100% per branch) already measures branch reliability:
- High bootstrap (>70%): Branch is well-supported
- Low bootstrap (<50%): Branch is uncertain

**Simple approaches** (not novel):
```python
tree_weight = mean(bootstrap_values)  # Simple average
tree_weight = 1 if mean(bootstrap) > 0.7 else 0  # Threshold
```

**Our contribution** (novel): Learn the **optimal aggregation**:
```python
tree_weight = learned_function(
    mean_bootstrap,
    min_bootstrap,
    bootstrap_variance,
    n_low_support_branches,
    outlier_score,
    missing_taxa_rate,
    ...
)
```

**Why learning is better than simple aggregation:**
- Maybe min_bootstrap matters more than mean?
- High mean but one very low branch → should downweight?
- Outlier trees should be downweighted regardless of bootstrap?
- Optimal combination is data-dependent

**No existing method learns the optimal way to weight trees.**

**Architecture**:

```
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Compute Tree-Level Features (per gene tree)        │
├─────────────────────────────────────────────────────────────┤
│  Tree 1 → [mean_support, n_taxa, outlier_score, ...]       │
│  Tree 2 → [mean_support, n_taxa, outlier_score, ...]       │
│  ...                                                        │
│  Tree N → features                                          │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: Learn Tree Weights (attention mechanism)           │
├─────────────────────────────────────────────────────────────┤
│  tree_features → learned scorer → weight w_i                │
│  w = softmax([w_1, w_2, ..., w_N])                         │
│  "Trust tree 1 at 0.02, tree 47 at 0.05, ..."              │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: Weighted Bipartition Scoring                       │
├─────────────────────────────────────────────────────────────┤
│  weighted_freq(B) = Σ w_i × 𝟙(B in tree_i)                 │
│  score(B) = f(weighted_freq, support, consistency, ...)     │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 4: ILP Selection (like Polyphest) + Tree Building     │
└─────────────────────────────────────────────────────────────┘
```

**Tree-Level Features** (topology-based):

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `mean_branch_support` | Average support across internal nodes | Low = uncertain tree |
| `min_branch_support` | Minimum support | Weakest link |
| `support_variance` | How variable is support | Consistency |
| `n_taxa` | Number of leaves | Completeness |
| `missing_taxa_rate` | Fraction of expected taxa missing | Missing data |
| `resolution` | Fraction of nodes that are bifurcating | Tree quality |
| `outlier_score` | Distance to other trees (RF or similar) | Anomalous tree? |
| `diameter` | Maximum leaf-to-leaf distance | Tree depth/shape |

**Novelty**: **HIGH** - No existing method learns tree weights

**Justification**: Optimal estimation theory says we should weight observations by reliability. Gene trees are not equally reliable.

**Interpretability**: "Tree 47 was downweighted because it had low support and was an outlier"

---

### 2.5 Option C: GNN-Based Approach (Future/High-Risk)

**Idea**: Instead of hand-crafted features, learn tree representations with a Graph Neural Network.

```
Gene Tree → GNN Encoder → Embedding
All embeddings → Attention Aggregation → Combined representation
Combined → Decoder → Bipartition scores or MUL-tree
```

**Challenges**:
- **Node features**: Unclear what to use (species one-hot? support?)
- **Compression**: Tree → vector loses structural information
- **Decoding**: Vector → tree is non-trivial
- **Data**: May need more training data than we have

**Novelty**: High
**Risk**: High (may not work well)

**Recommendation**: Keep as future work. Start with Option B.

---

### 2.6 Selection Algorithm: ILP (Not Greedy)

All options should use **ILP** for selecting compatible bipartitions (like Polyphest), not greedy:

```
Maximize:   Σ score(bp_i) × x_i
Subject to: x_i ∈ {0, 1}
            x_i + x_j ≤ 1  for all incompatible pairs
```

This finds the **maximum weight compatible set** - optimal, not heuristic.

---

### 2.7 Positioning Summary

```
                    CONSENSUS-BASED                 RECONCILIATION-BASED
                          │                                  │
              ┌───────────┴───────────┐                      │
              │                       │                      │
         Polyphest              Standard Consensus    GRAMPA / PADRE / MPSUGAR
         (clusters,            (bipartitions,        (reconciliation scoring,
          frequency)            frequency)            fixed costs)
              │                       │                      │
              │                       │                      │
              └───────────┬───────────┘                      │
                          │                                  │
                  ┌───────┴───────┐                          │
                  │               │                          │
              OPTION A        OPTION B                       │
           (learned          (learned tree weights           │
            bipartition    + learned bipartition             │
            scoring)           scoring)                      │
                                  │                          │
                          ┌───────┴───────┐                  │
                          │               │                  │
                     OPTION C        LEARNED COSTS           │
                     (GNN-based)    (future work) ───────────┘
```

---

## 3. Ploidy Handling Options

### 3.1 The Ploidy Problem

MUL-trees require knowing how many copies of each species to include. For example:
- Species A: 2 copies (tetraploid)
- Species B: 1 copy (diploid)
- Species C: 4 copies (octoploid)

### 3.2 How Existing Methods Handle Ploidy

| Method | Ploidy Input | Notes |
|--------|--------------|-------|
| **Polyphest** | User provides for ALL species | Required "consensus multiset" |
| **GRAMPA** | Searches over possibilities | Tries different mappings |
| **MPSUGAR** | User provides taxon map | Maps gene copies to species |

### 3.3 Our Options

#### Option P1: User Provides Ploidy for All Species (Like Polyphest)

**Input**: Gene trees + ploidy levels for each species

```python
ploidy = {
    "SpeciesA": 2,
    "SpeciesB": 1,
    "SpeciesC": 4
}
```

**Pros**: Simplest, uses expert knowledge
**Cons**: Requires user knowledge, may not always be available

#### Option P2: User Provides Ploidy for Some, We Infer Rest

**Input**: Gene trees + ploidy for known species

```python
ploidy = {
    "SpeciesA": 2,   # known tetraploid
    "SpeciesB": 1,   # known diploid
    "SpeciesC": None # unknown - infer
}
```

**How to infer**:
- Learn from patterns in gene trees (copy counts, clustering)
- Train predictor on species with known ploidy

**Pros**: Flexible, reduces user burden
**Cons**: Inference may be uncertain

#### Option P3: Infer All Ploidy Levels (Most Ambitious)

**Input**: Gene trees only (no ploidy info)

**How**: Predict ploidy per species from gene tree patterns:

| Feature | Description |
|---------|-------------|
| `mode_copies` | Most common copy count in gene trees |
| `mean_copies` | Average copies per tree |
| `max_copies` | Maximum observed |
| `copy_variance` | Variability |
| `monophyly_rate` | How often copies cluster together |
| `sister_rate` | How often copies are sisters |

**Pros**: No user input needed, fully automated
**Cons**: May be inaccurate, adds complexity

### 3.4 Recommendation

**Start with Option P1** (user provides all ploidy), then explore P2/P3 as extensions.

Why:
- Matches Polyphest (fair comparison)
- Focuses novelty on tree weighting / bipartition scoring
- Ploidy inference can be added later if time permits

### 3.5 Ploidy in the Pipeline

```
                         Option P1           Option P2/P3
                              │                    │
                              ▼                    ▼
                         User provides      Learn/Infer
                         all ploidy           ploidy
                              │                    │
                              └────────┬───────────┘
                                       │
                                       ▼
                          ┌─────────────────────────┐
                          │  Filter bipartitions    │
                          │  by ploidy constraints  │
                          │  (like Polyphest)       │
                          └─────────────────────────┘
                                       │
                                       ▼
                               ILP Selection
```

---

## 4. Theoretical Foundation

### 4.1 Why Machine Learning Should Work

#### Argument 1: Optimal Weighting Theory

**Theorem (informal)**: The optimal estimator weights observations by their informativeness. Uniform weighting is optimal **only** when all observations are equally reliable.

**Reality**: Gene trees vary dramatically in quality:
- Trees from genes under strong selection → more reliable
- Trees with high ILS → noisy, conflicting
- Trees with gene loss → incomplete
- Trees with short internal branches → uncertain

**Current methods**: Treat all gene trees equally
```
score(bipartition) = count(bipartition in gene trees) / N
```

**Our approach**: Learn optimal weights
```
score(bipartition) = f(frequency, support, consistency, polyploidy_indicators, ...)
```

#### Argument 2: Non-Linear Feature Interactions

A bipartition at 60% frequency might be:
- **True** if supported by high-quality trees with strong branch support
- **False** if driven by systematic gene tree estimation error

The "right" decision depends on **multiple features interacting non-linearly**. Fixed thresholds cannot capture this.

#### Argument 3: Context Awareness

The optimal threshold varies by:
- ILS level (high ILS → lower frequency threshold acceptable)
- Tree quality distribution
- Polyploidy patterns
- Network complexity

Learning enables **context-specific decisions** that fixed rules cannot make.

### 4.2 Mathematical Framework

Let:
- G = {g₁, g₂, ..., gₙ} be the input gene trees
- B = {b₁, b₂, ..., bₘ} be candidate bipartitions extracted from G
- φ(bᵢ, G) be a feature vector for bipartition bᵢ

**Classic consensus**:
```
P(bᵢ ∈ true tree) = 1[frequency(bᵢ) > θ]
```

**Our learned approach**:
```
P(bᵢ ∈ true tree) = σ(f_θ(φ(bᵢ, G)))
```

Where f_θ is learned from training data with ground truth.

---

## 5. Proposed Method: Option B (Recommended)

### 5.1 Overview: Learned Tree Weighting + Bipartition Scoring

```
Input: N Gene Trees (Newick format) + Ploidy Info (user-provided)
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   1. Extract Candidate Bipartitions                         │
│   - All bipartitions from trees (multisets for polyploidy)  │
│   - Filter by ploidy constraints                            │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   2. Compute Tree-Level Features (NOVEL)                    │
│   - Per gene tree: support stats, completeness, outlier     │
│   - Learn tree weights via attention                        │
│   - "Which trees to trust"                                  │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   3. Compute Bipartition Features (weighted by tree trust)  │
│   - Weighted frequency (using learned tree weights)         │
│   - Support, consistency, polyploidy, structural features   │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   4. Learned Bipartition Scorer                             │
│   - P(bipartition ∈ true MUL-tree | features)               │
│   - XGBoost or neural network                               │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   5. ILP Selection (like Polyphest)                         │
│   - Maximize: Σ score(bp) × selected(bp)                    │
│   - Subject to: compatibility constraints                   │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│   6. MUL-Tree Construction                                  │
│   - Build tree from selected bipartitions                   │
│   - Fold to network (optional, like Polyphest)              │
└─────────────────────────────────────────────────────────────┘
         │
         ▼
Output: MUL-Tree + Interpretability Report
        (which trees were trusted, which features mattered)
```

### 5.2 Feature Engineering (Topology-Based Only)

**Important**: Like Polyphest, GRAMPA, and other methods, we use **topology only** (no branch lengths). Branch support values are annotations, not lengths.

#### Tree-Level Features (NEW - for learning tree weights)

| Feature | Description | Intuition |
|---------|-------------|-----------|
| `mean_support` | Average internal node support | Low = uncertain tree |
| `min_support` | Minimum support | Weakest link |
| `support_variance` | Standard deviation of support | Consistency |
| `n_taxa` | Number of leaves | Completeness |
| `missing_taxa_rate` | Fraction of expected taxa missing | Missing data |
| `resolution` | Fraction bifurcating nodes | Tree quality |
| `outlier_score` | RF distance to median tree | Anomaly detection |
| `diameter` | Max path length (in nodes) | Tree shape |

#### Bipartition Features (weighted by tree trust)

**Frequency Features** (now weighted):
| Feature | Description |
|---------|-------------|
| `frequency` | Raw count / N |
| `weighted_frequency` | Σ tree_weight × 𝟙(bp in tree) |
| `presence_variance` | Variance across trees |

**Support Features**:
| Feature | Description |
|---------|-------------|
| `mean_support` | Average support when bipartition present |
| `min_support` | Minimum support |
| `max_support` | Maximum support |
| `support_std` | Standard deviation |

**Consistency Features**:
| Feature | Description |
|---------|-------------|
| `conflict_rate` | Fraction of bipartitions that conflict |
| `consistency_score` | Co-occurrence with other high-freq bipartitions |

**Polyploidy Features**:
| Feature | Description |
|---------|-------------|
| `involves_polyploid` | Does bipartition contain polyploid species? |
| `copy_ratio` | Balance of copies between sides |
| `monophyly_score` | Do copies cluster together? |
| `left_polyploid_count` | Number of polyploid species on left |
| `right_polyploid_count` | Number of polyploid species on right |

**Structural Features**:
| Feature | Description |
|---------|-------------|
| `split_size` | Size of smaller side |
| `split_size_ratio` | Ratio of split sizes |
| `depth_indicator` | Estimated depth in tree |
| `total_species` | Total unique species in bipartition |

### 5.3 Model Architecture

**Component 1: Tree Weight Learner**
```python
# Input: tree-level features for each gene tree
# Output: weight for each tree (softmax normalized)

tree_features → MLP/XGBoost → raw_score
weights = softmax([score_1, ..., score_N])
```

**Component 2: Bipartition Scorer**
```python
# Input: bipartition features (using weighted frequency from Component 1)
# Output: P(bipartition in true MUL-tree)

bipartition_features → XGBoost → probability
```

**Training**: Joint or two-stage
- Two-stage: First train tree weighter, then bipartition scorer
- Joint: End-to-end backpropagation (requires neural network)

### 5.4 Why This Design?

| Design Choice | Justification |
|---------------|---------------|
| Topology-only | Matches existing methods (fair comparison) |
| Tree weighting | Novel contribution (no one does this) |
| ILP selection | Optimal (like Polyphest), not greedy |
| XGBoost | Interpretable, works with limited data |
| User-provided ploidy | Matches Polyphest (fair comparison) |

---

## 6. Current Implementation Status

### 6.1 Completed Components (Option A - Baseline)

| Component | Status | Description |
|-----------|--------|-------------|
| **Bipartition Extractor** | ✅ Complete | Extracts bipartitions with polyploidy support |
| **Bipartition Feature Builder** | ✅ Complete | Computes 18 features per bipartition |
| **XGBoost Scorer** | ✅ Complete | Classifier with CV and feature importance |
| **Tree Builder** | ✅ Complete | Greedy reconstruction from bipartitions |
| **Training Pipeline** | ✅ Complete | Multi-network training with cross-validation |

### 6.2 TODO: Components for Option B (Main Method)

| Component | Status | Description |
|-----------|--------|-------------|
| **Tree-Level Feature Extractor** | ❌ TODO | Compute features per gene tree |
| **Tree Weight Learner** | ❌ TODO | Learn attention weights over trees |
| **ILP Selector** | ❌ TODO | Replace greedy with ILP (like Polyphest) |
| **Ploidy Constraint Handler** | ❌ TODO | Filter bipartitions by ploidy |
| **Joint Training Pipeline** | ❌ TODO | Train tree weighter + bipartition scorer |

### 6.3 Code Structure

```
simulations/scripts/ml_method/
├── data/
│   ├── bipartition_extractor.py   # ✅ Parse trees, extract bipartitions
│   ├── feature_builder.py         # ✅ Bipartition features
│   └── tree_features.py           # ❌ TODO: Tree-level features
├── models/
│   ├── xgb_scorer.py              # ✅ Bipartition scorer
│   ├── tree_weighter.py           # ❌ TODO: Learn tree weights
│   ├── tree_builder.py            # ✅ MUL-tree reconstruction
│   └── ilp_selector.py            # ❌ TODO: ILP optimization
├── training/
│   └── train.py                   # ✅ Training (needs update for Option B)
└── analysis/
    └── (to be implemented)        # Interpretability tools
```

### 6.4 What's Ready to Run (Option A)

- End-to-end pipeline from gene trees to MUL-tree (greedy, not ILP)
- Leave-one-out cross-validation
- Feature importance analysis
- Model persistence (save/load)

---

## 7. Development Roadmap

### 7.1 Timeline

| Month | Focus | Deliverables |
|-------|-------|--------------|
| **1** | Complete Option B Implementation | Tree weighter, ILP selector, ploidy handling |
| **2** | Synthetic Data Generation | Generate 100-500 synthetic networks via SimPhy |
| **3** | Training & Iteration | Train model on synthetic data; iterate on features |
| **4** | Evaluation | Test on 21 biological networks; compare to methods |
| **5** | Analysis | Benchmarking, interpretability, ablation studies |
| **6** | Publication | Write paper; code release |

### 7.2 Immediate Next Steps

1. **Implement tree-level features**: `tree_features.py`
2. **Implement tree weight learner**: `tree_weighter.py`
3. **Implement ILP selector**: `ilp_selector.py` (using PuLP)
4. **Add ploidy constraint handling**
5. **Generate synthetic training data**: SimPhy
6. **Train and evaluate**

### 7.3 Implementation Priority

```
Week 1-2: Tree-level features + tree weighter
Week 3:   ILP selector (replace greedy)
Week 4:   Ploidy handling + integration
Week 5-6: Synthetic data generation
Week 7-8: Training on synthetic data
Week 9+:  Evaluation on biological data
```

### 7.4 Future Extensions (If Time Permits)

| Extension | Description | Priority |
|-----------|-------------|----------|
| **Ploidy Inference (P2/P3)** | Infer ploidy from gene trees | Medium |
| **GNN Encoder (Option C)** | Learn tree representations | Low |
| **Network Output** | Fold MUL-tree to network | Medium |
| **Learned Reconciliation** | Learn GRAMPA-style costs | Low |

---

## 8. Expected Contributions

### 8.1 Scientific Contributions

1. **First ML method** for MUL-tree inference from gene trees
2. **Learned tree weighting** - novel, no existing method does this
3. **Learned bipartition scoring** - replaces fixed thresholds
4. **Interpretable predictions** - which trees were trusted, which features mattered
5. **Comprehensive evaluation** on biological networks

### 8.2 Practical Impact

- **For researchers**: Better MUL-tree inference for polyploid studies
- **For the field**: Understanding of what makes gene trees/bipartitions informative
- **For methods**: Framework for applying ML to phylogenetics

### 8.3 Interpretability (Key Differentiator)

Unlike black-box methods, we can explain:

1. **Feature importance**: "Frequency and consistency are most predictive"
2. **Per-bipartition**: "This split was selected because: frequency=0.85, support=0.92"
3. **Per-network**: "Network X was challenging due to low-frequency true bipartitions"

---

## 9. Data Strategy & Evaluation Plan

### 9.1 Training vs Testing Split

**Critical**: We do NOT train on the biological networks. They are held out for final evaluation.

```
┌─────────────────────────────┐     ┌─────────────────────────────┐
│     SYNTHETIC DATA          │     │     BIOLOGICAL DATA         │
│     (Training)              │     │     (Testing)               │
├─────────────────────────────┤     ├─────────────────────────────┤
│ - Generate via SimPhy       │     │ - 21 biological networks    │
│ - Many networks (100-500+)  │     │ - 5 replicates × 1000 trees │
│ - Varied ILS, polyploidy    │     │ - HELD OUT for evaluation   │
│ - Known ground truth        │     │ - Compare to other methods  │
└─────────────────────────────┘     └─────────────────────────────┘
              │                                   │
              ▼                                   ▼
         Train Model                      Final Evaluation
```

### 9.2 Data Summary

| Dataset | Size | Use |
|---------|------|-----|
| **Synthetic networks** | 100-500 networks | **Training** |
| **21 biological networks** | 21 × 5 replicates × 1000 trees | **Testing only** |

### 9.3 Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| **Bipartition F1** | Precision/recall on bipartition prediction | > 0.8 |
| **Edit Distance** | MUL-tree similarity to ground truth | < existing methods |
| **Copy Number Accuracy** | Polyploid detection | > 0.9 |
| **Completion Rate** | Successful inference | 100% |

### 9.4 Baselines

**Existing methods:**
- GRAMPA
- Polyphest (p50, p70, p90)
- PADRE
- MPSUGAR

**Our baselines:**
- Frequency-only (no ML, no weighting)
- Mean-bootstrap weighting (simple, not learned)
- Option A (learned bipartition scoring, no tree weighting)
- Option B (learned tree weighting + bipartition scoring) ← **Our method**

### 9.5 Ablation Studies

| Ablation | Question |
|----------|----------|
| Learned vs mean-bootstrap weighting | Does learning improve over simple aggregation? |
| Learned vs no tree weighting | Does tree weighting help at all? |
| Remove frequency features | How much does raw frequency contribute? |
| Remove support features | Does branch support matter? |
| Remove polyploidy features | Are polyploidy-specific features important? |
| Greedy vs ILP | Does optimal selection matter? |

---

## 10. Publication Strategy

### 10.1 Target Venues

| Venue | Type | Fit |
|-------|------|-----|
| **Bioinformatics** | Journal | Methods-focused, good fit |
| **Systematic Biology** | Journal | Strong methodology focus |
| **MBE** | Journal | If biological insights are strong |
| **RECOMB/ISMB** | Conference | Faster review, ML audience |

### 10.2 Paper Narrative

**Title Options**:
- "Learning to Weight Gene Trees for MUL-Tree Inference"
- "Beyond Counting: Machine Learning for Polyploid Species Tree Reconstruction"
- "Attention-Weighted Consensus for MUL-Tree Inference"

**Story Arc**:
1. **Problem**: MUL-tree inference treats all gene trees equally, uses fixed rules
2. **Insight**: Gene trees vary in reliability; optimal selection depends on multiple features
3. **Method**: Learn which trees to trust + learn bipartition scoring
4. **Results**: Improved accuracy + interpretable decisions
5. **Impact**: Understanding what makes gene trees and bipartitions reliable

### 10.3 Key Figures

1. Method overview diagram (tree weighting + bipartition scoring)
2. Tree weight distribution (which trees were trusted)
3. Feature importance analysis
4. Performance comparison with existing methods
5. Performance by network complexity (H, polyploid count)
6. Example case study with explanation

---

## 11. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Not enough training data | Medium | High | Many bipartitions per network; generate synthetic |
| Tree weighting doesn't help | Medium | Medium | Still have bipartition scoring (Option A fallback) |
| Doesn't beat existing methods | Medium | Medium | Focus on interpretability angle |
| Overfitting | Low | High | Cross-validation; regularization; simple models |
| Feature engineering inadequate | Medium | Medium | Iterate; try different features |
| Time overrun | Low | Medium | Option A already complete as baseline |

---

## 12. Why This Will Work

### 12.1 Theoretical Soundness

- **Optimal weighting theory** supports learned over uniform weighting
- **Supervised learning** with clear ground truth signal
- **Feature design** captures known phylogenetic principles
- **Tree weighting** is theoretically justified (not all observations equally reliable)

### 12.2 Practical Advantages

- **Interpretable**: Which trees were trusted, which features mattered
- **Scalable**: Once trained, inference is fast
- **Robust**: Handles varying tree quality by downweighting bad trees
- **Extensible**: Framework allows adding new features/components

### 12.3 Feasibility

- **Baseline implementation complete**: Option A runs end-to-end
- **Data available**: 21 biological networks + simulation capability
- **Evaluation framework**: Existing benchmarking infrastructure
- **Timeline realistic**: 6 months with buffer
- **Fallback**: If Option B doesn't work, Option A is still publishable

### 12.4 Novelty

- **Learned optimal tree weighting**: Bootstrap measures reliability, but we learn the optimal way to aggregate it (not just mean)
- **Learned bipartition scoring**: Replaces fixed frequency thresholds with multi-feature scoring
- **First ML method** for MUL-tree inference from gene trees
- **Interpretable**: Which trees were trusted and why, which features mattered

**Important distinction**: The novelty is not "using bootstrap" (that's obvious) but **learning the optimal aggregation function** from data.

---

## 13. Summary

We propose a machine learning approach to MUL-tree inference that:

1. **Learns which gene trees to trust** (main novelty - no existing method does this)
2. **Learns optimal bipartition scoring** instead of fixed frequency thresholds
3. **Uses ILP for selection** (like Polyphest, optimal not greedy)
4. **Provides interpretability** through tree weights and feature importance
5. **Has fallback options** (Option A if tree weighting doesn't help)
6. **Is implementable** in the proposed timeline

**The key insight**: Not all gene trees are equally reliable, and not all bipartitions at the same frequency are equally likely to be true. Learning enables context-aware decisions that fixed heuristics cannot make.

---

## 14. Method Comparison Summary

| Aspect | Polyphest | Simple Bootstrap | Our Method (Option B) |
|--------|-----------|------------------|----------------------|
| Representation | Clusters | Bipartitions | Bipartitions |
| Tree weighting | Equal (1/N) | mean(bootstrap) | **Learned (optimal aggregation)** |
| Scoring | Frequency | Weighted frequency | **18 features, learned** |
| Selection | MWACC (ILP) | Threshold | ILP |
| Ploidy | User-provided | - | User-provided |
| Interpretability | Limited | Limited | **High** (tree weights + feature importance) |

### Key Innovation

Bootstrap values measure branch reliability, but:
- **Polyphest**: Ignores bootstrap for tree weighting (all trees equal)
- **Simple baseline**: Uses mean(bootstrap) as tree weight
- **Our method**: **Learns optimal combination** of bootstrap + other features

The novelty is not "using bootstrap" but **learning how to optimally use it**.

---

## Appendix: Technical Details

### A.1 Bipartition Representation

A bipartition splits taxa into two groups. For MUL-trees with polyploidy, we use **multisets**:
- Standard: `{A, B} | {C, D}`
- With polyploidy: `{A, A, B} | {C, D}` (species A appears twice)

### A.2 Compatibility

Two bipartitions are **compatible** if they can coexist in the same tree. We use greedy selection: sort by score, add if compatible with all selected.

### A.3 Dependencies

- Python 3.8+
- ete3 (tree parsing)
- xgboost (classification)
- numpy, scikit-learn (utilities)

### A.4 Computational Requirements

- Training: ~1 hour on CPU for 21 networks
- Inference: ~1 minute per network
- GPU: Not required for current approach
