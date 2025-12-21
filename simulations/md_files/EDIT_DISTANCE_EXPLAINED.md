# Edit Distance Metric - Detailed Explanation

## What is Edit Distance?

**Edit distance** (also called graph edit distance) measures the minimum number of operations needed to transform one network into another.

### Algorithm

Located in `simulations/scripts/reticulate_tree.py:798-822`:

```python
def get_edit_distance(self, other: 'ReticulateTree', normalize=True):
    # Uses NetworkX optimize_graph_edit_distance
    distance = next(nx.optimize_graph_edit_distance(
        self.dag, other.dag,
        node_match=lambda u, v: u.get('label') == v.get('label')
    ))

    if normalize:
        normalization = max(
            self.dag.number_of_nodes() + self.dag.number_of_edges(),
            other.dag.number_of_nodes() + other.dag.number_of_edges()
        )
        distance = distance / normalization

    return distance
```

### Key Properties

1. **Operations counted**:
   - Node insertions
   - Node deletions
   - Edge insertions
   - Edge deletions

2. **Node matching**: Nodes are considered equal if they have the same label

3. **Normalization**:
   ```
   normalized_distance = raw_distance / max(nodes₁ + edges₁, nodes₂ + edges₂)
   ```

4. **Range**:
   - **0** = Identical networks (perfect inference)
   - **~1.0** = Completely different networks
   - Typical values: 0.0 - 0.8

### Why This Metric?

**Appropriate for phylogenetic networks because:**
- ✅ Handles reticulation nodes properly
- ✅ Considers both topology (nodes/edges) and structure
- ✅ Works for networks of different sizes
- ✅ Normalized for fair comparison

**Better than Robinson-Foulds because:**
- ❌ RF only works for trees (not networks with reticulations)
- ❌ RF requires unique leaf labels (not MUL-trees)
- ❌ RF doesn't handle polyploidy well

### Interpretation

| Edit Distance | Interpretation |
|--------------|----------------|
| 0.00 - 0.10  | Excellent - nearly perfect reconstruction |
| 0.10 - 0.30  | Good - minor differences |
| 0.30 - 0.50  | Moderate - noticeable differences |
| 0.50 - 0.70  | Poor - substantial differences |
| > 0.70       | Very poor - major structural errors |

### Example

**True network**: 5 nodes, 4 edges (total = 9)
**Inferred network**: 6 nodes, 5 edges (total = 11)
**Operations needed**: 3 (1 extra node, 2 edge differences)

```
normalized_distance = 3 / max(9, 11) = 3 / 11 = 0.27
```

This indicates a **good** reconstruction with minor differences.

## How It's Used in Your Analysis

Edit distance is calculated for each:
- **Configuration** (ILS low/medium/high)
- **Method** (GRAMPA, Polyphest, MPSUGAR, PADRE)
- **Network** (21 empirical networks)
- **Replicate** (5 replicates per network)

**Aggregation**:
- Per network: Mean across 5 replicates
- Per method: Mean/median across all networks
- Used to rank methods and identify difficult networks

## Related Metrics

Your pipeline also computes:

1. **num_rets_diff**: Difference in reticulation count (inferred - true)
2. **ploidy_diff**: Polyploid species identification accuracy (TP/FP/FN)
3. **ret_leaf_jaccard**: Jaccard similarity of reticulation leaf sets
4. **ret_sisters_jaccard**: Jaccard similarity of sister relationships

All metrics together provide comprehensive evaluation!
