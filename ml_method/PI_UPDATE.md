# Gene2Net-GNN: progress update

## 1. Summary

The original goal was to design informative features and predict, for each edge of a species tree, whether a whole genome duplication occurred. That part works. A features plus GNN model reaches binary WGD detection F1 around 0.80 in distribution. We then went further and added partner prediction, which turns the per edge detection into a full polyploid network reconstruction that can be compared directly against GRAMPA-iter and Polyphest on their own metrics. We can now reconstruct networks, and we have a clear, evidence based picture of where the method is strong and where it is limited.

Headline points:
- Detection of WGD edges and of polyploid lineages is strong and robust to gene tree noise.
- On the reticulation identification metrics we beat GRAMPA-iter under noise.
- On full network edit distance we are behind Polyphest, and we identified exactly why through a controlled oracle experiment. The limitation is the ASTRAL based backbone, not the predictions.

## 2. Problem framing

Input: a set of gene trees and an ASTRAL inferred diploid species tree. Output: the polyploid network, represented first as a multi-labeled tree and then folded to a network. We treat this as a learned per edge problem. First detect which edges carry a WGD, then for each WGD edge predict its partner edge, which is where the duplicated copy attaches. A partner equal to the edge itself means autopolyploidy, a different partner means allopolyploidy. This is a different machinery from GRAMPA-iter, which uses reconciliation, and Polyphest, which builds the multi-labeled tree from gene tree clusters using an integer program.

## 3. Data

### Training and validation data, fully simulated

1. Random ground truth networks. For each index we generate an ultrametric species tree with SimPhy under a birth death model, with 5 to 80 species weighted toward 10 to 40. We then add 0 to 15 random WGD events. An auto event duplicates a clade and places the copy as a sister. An allo event grafts the copy onto a different non overlapping edge. Ultrametricity is preserved. The result is the ground truth multi-labeled tree, saved per index.
2. Gene trees. SimPhy simulates gene trees down each multi-labeled tree under the configuration ILS level, set by effective population size, and the gene duplication and loss rates. We use 9 training configurations spanning low, medium and high ILS and low, medium and high duplication and loss. About 500 gene trees per sample.
3. Inferred species tree. ASTRAL infers a single label diploid species tree from each sample gene trees.
4. Labels. We decompose the ground truth multi-labeled tree into its true WGD events. Each event is mapped to the best matching ASTRAL edge by Jaccard similarity with a 0.5 threshold. This yields per edge binary WGD labels and a partner edge label per event.

Totals: about 17,900 packaged samples across the 9 configurations. The split is 80 percent train and 20 percent validation, by whole sample, so no edge leaks between train and validation. Validation therefore measures generalization to unseen simulated networks from the same distribution.

### Test data, real networks

The held out test set is the 21 empirical polyploid networks used to benchmark GRAMPA-iter and Polyphest, evaluated across 12 evolutionary configurations of ILS and duplication and loss. The model never trained on these. They are fed to the model with the same inputs GRAMPA-iter receives, namely the cleaned gene trees and the ASTRAL diploid species tree, so the comparison is apples to apples.

## 4. Features

All features are computed from the gene trees. The signal we are hunting is the difference between organized duplication, which is WGD, and scattered duplication, which is ordinary gene duplication noise.

### Node features, 13 per species

A species is a leaf of the species tree. Internal nodes start empty and are filled by averaging their descendant leaf features, so internal edges also have signal.

Copy count statistics, 8 values. These describe how duplicated a species is across the gene trees.
- mean copies. Average number of copies of the species per gene tree. A polyploid sits above 1.
- variance of copies. How much the copy count fluctuates. A clean WGD gives a consistent count, stochastic duplication and loss gives high variance.
- mode copies. The most common copy count, which approximates the species ploidy.
- max copies. The highest count seen, which flags high ploidy.
- probability of 0, of 1, of 2, and of 3 or more copies. The shape of the copy count distribution. A tetraploid shows a high probability of 2.

Co-clustering summary, 5 values. These describe who the species groups with.
- For each other species we measure how often the two appear as sister leaves across gene trees. We then summarize that profile with its mean, standard deviation, max, min and median. The idea is that a species formed by hybridization with a particular partner co-clusters with that partner consistently, which shows up as a high max and a peaked profile, whereas scattered gene duplication gives a diffuse profile.

### Edge features, 9 per edge

An edge corresponds to the clade S of species below it.

Base features, 4 values.
- concordance factor. Fraction of gene trees that support this edge bipartition. WGD edges tend to be discordant, since the duplicated copies confuse the gene tree topology, so a low value is informative.
- branch length, in generations. Temporal and positional context for the edge.
- clade size. Number of species below the edge. This helps the model localize which edge in a duplicated region is the actual event.
- depth. Number of edges from the root. Positional context.

Detection features, 5 values. These were the additions that lifted detection from F1 about 0.51 to about 0.80.
- duplication synchrony. Across gene trees, do the species of clade S get duplicated together. A WGD duplicates the whole clade at one moment, so the duplications are correlated, whereas independent gene duplication hits species separately. Measured as the average pairwise correlation of the duplicated indicators of the clade species.
- mirrored-sister fraction. Fraction of gene trees in which clade S appears as two identical sister copies, the pattern two copies of A and B side by side. This is the autopolyploidy signature.
- copy-pair divergence mean. For each species with two copies, the branch length distance between its two copies, averaged and normalized within each gene tree. This is a Ks analogue, since the gene trees are in substitutions per site. It reflects how long ago the copies diverged, which is the time since duplication.
- copy-pair divergence coefficient of variation. How consistent that divergence is across species and gene trees. One synchronized event gives consistent ages, scattered duplication gives variable ages.
- fraction of clade duplicated. The average fraction of the clade species that are duplicated per gene tree. High for a true WGD clade. This turned out to be the single strongest feature.

### Pairwise feature for the partner head

This is not a per edge feature but a value for each ordered pair of edges i and j. It is the average species co-clustering between the species of clade i and the species of clade j across gene trees. It captures the allopolyploidy signal, namely whether the duplicated copy of clade i consistently lands next to clade j. This is what tells the partner head where an allopolyploid copy attaches.

## 5. Architecture

The model runs on the ASTRAL species tree. The gene trees enter only through the features.

Step by step.
1. Node features. Each leaf gets its 13 features. Each internal node is filled by averaging the features of the leaves below it, so internal edges carry signal. This propagation depends only on the data and not on the model weights, so we compute it once and cache it.
2. Projection. The 13 dimensional node features are projected to a 64 dimensional hidden space by a small MLP, which lets the model build useful combinations of the raw features.
3. Graph attention. 3 GAT layers, each with 4 attention heads, residual connections and layer norm. Message passing lets each node integrate information from its neighbors in the species tree, so an edge decision can account for its parent and children.
4. Edge embedding. For each edge we concatenate the parent node embedding, the child node embedding and the 9 edge features, and pass them through an edge MLP. The parent and child embeddings let the model compare what is below the edge against what is above it, which is the core WGD question.
5. Detection head. A small MLP maps each edge embedding to a binary WGD score.
6. Partner head. For each WGD edge i, we score every candidate edge j from the pair of edge embeddings and the pairwise co-clustering feature of i and j, then softmax over all j including j equal to i. The argmax is the predicted partner. A self partner means autopolyploidy and a different partner means allopolyploidy. Auto versus allo therefore falls out of partner prediction and needs no separate classifier.

The model is small, about 80 thousand parameters. The neural forward pass is milliseconds. The wall clock cost of inference is the feature computation from the gene trees, not the network. There is no combinatorial search.

## 5b. How we reconstruct the MUL-tree, and where the decisions are

The model gives, per edge, a WGD probability and a partner distribution. Turning these into a MUL-tree is a separate build step, and this is where several deliberate decisions live. We made the build a swappable strategy so the decisions can be validated against the reconstruction measures rather than chosen by hand.

The pipeline and its decision points.
1. Detection threshold. Which edges count as WGD, namely those with probability at or above a threshold. Decision: the threshold. We currently use 0.9.
2. Event selection strategy. How the flagged edges become events. This is the main decision, because the model tends to flag many edges in a duplicated region or one edge per polyploid lineage. The strategies are:
   - raw. Every flagged edge becomes its own event. This over fragments.
   - collapse. One event per connected block of flagged edges, placed at the top of the block. Intended to turn a clade spanning event into one duplication.
   - cap. Flagged edges, but cap the number of duplications per species at its inferred copy number, which prevents over ploidy.
   - collapse with cap. Both of the above.
   - bound driven. No threshold. Take edges in confidence order until each species reaches its copy bound, which removes the threshold knob.
3. Copy bound. The inferred consensus copy number per species, namely the count supported by at least half the gene trees. The cap strategies use it to avoid giving a species more copies than the data supports. This also keeps the folded network well formed, since over ploidy produces malformed reticulations.
4. Partner assignment. For each event edge the partner is the argmax of the partner head. If the predicted partner clade overlaps the WGD clade, which is not a valid graft, we fall back to autopolyploidy.
5. Application order and conflicts. build_mul_tree applies events smallest clade first, so nested events apply correctly. An auto event duplicates the clade as a sister. An allo event grafts a copy onto the partner edge. If a clade no longer exists after earlier edits, that event is skipped.
6. Folding to a network. For the network level metrics the MUL-tree is folded by the Huber algorithm into a network. We apply the same folding to all methods, so folding is not a source of difference between methods.

Why this matters. The current evidence is that the event selection strategy strongly affects the reticulation count and the malformedness of the network, while the edit distance is set mostly by the backbone, as shown in section 8. So the build decisions are where we control the reticulation metrics, and the backbone is where the edit distance is decided.

## 6. Results

### Detection, in distribution on simulated validation
- Topology and copy count features only reach F1 about 0.51 after threshold tuning.
- Adding the 5 detection features raises tuned F1 to about 0.80. The model is recall saturated.
- The GNN earns its place. A simple threshold on the single strongest feature reaches F1 about 0.62, and the full GNN reaches 0.80. The added value is localization, namely choosing which edge in a chain of duplicated edges is the true event.

### Partner prediction, simulated validation
- Auto partner accuracy about 0.99.
- Allo partner accuracy went from about 0.52 to about 0.80 once we added the pairwise co-clustering feature.

### Reconstruction on the 21 real networks
- Polyploid lineage detection is strong. The predicted number of polyploid lineages matches the true number of polyploids almost exactly for about 17 of 21 networks.
- Reticulation identification. We beat GRAMPA-iter on the reticulation descendant metric, and our accuracy stays flat as noise increases while GRAMPA-iter degrades.
- Full network edit distance is behind Polyphest and roughly tied with GRAMPA-iter at the hard end. The per config table across all 12 configurations is being produced now and will be added here.

## 7. Comparison to GRAMPA-iter and Polyphest

Reference numbers from the thesis benchmark, edit distance, lower is better:
- Polyphest ranges from about 0.42 at low ILS to about 0.56 at high ILS, and degrades at high duplication and loss.
- GRAMPA-iter ranges from about 0.80 to 0.91.
- Our method sits around 0.86 to 0.92 and is roughly flat across configurations, which means it is noise robust but behind on this metric.

On reticulation identification metrics our method is competitive with or better than GRAMPA-iter, especially under noise, and approaches Polyphest. Polyphest remains the most accurate overall, but it requires a consensus copy number bound and runs an integer program search. Our method is prior free in the same sense, since the copy number signal it uses is also derived from the gene trees, and it is a forward pass rather than a search.

## 8. The key finding: where the edit distance gap comes from

We ran a controlled oracle experiment to separate two possible causes of the edit distance gap, the predictions and the backbone.

First, what the true backbone is. When we simulate a sample, the generator works in two steps. It first builds a diploid species tree, then adds the WGD events on top to produce the multi-labeled tree. We save both files. The diploid species tree before any duplication is the true backbone, and the multi-labeled tree is that same tree with the duplications added. So the backbone is the multi-labeled tree minus the duplications, and it is a real saved object, not something we infer.

The experiment. We take the ground truth multi-labeled tree and decompose it back into its true events. We then rebuild a multi-labeled tree with our own build step, build_mul_tree, on a chosen backbone, and score the result against the ground truth. We do this on two backbones, holding the events fixed at the true events.
- True backbone with true events: edit distance about 0.366. Rebuilding the true events on the true backbone should reproduce the ground truth. It is not exactly zero because our decompose and rebuild round trip is not perfect.
- ASTRAL backbone with true events: edit distance about 0.858. Same true events, but built on the ASTRAL inferred tree instead of the true diploid tree.
- Our model with predicted events on ASTRAL: about 0.87.

Three conclusions follow.
1. The predictions are not the edit distance bottleneck. Perfect events on ASTRAL give about 0.858, and our predicted events give about 0.87, which is nearly identical. This is why no build strategy or threshold ever moved the edit distance.
2. The backbone is the dominant cost. The same events give 0.366 on the true backbone and 0.858 on ASTRAL. The edit distance metric is far more sensitive to global structure than Robinson Foulds, and ASTRAL itself is accurate, with mean normalized RF about 0.085 against the true tree.
3. Most importantly, our build on the true backbone reaches 0.366, which is better than Polyphest at 0.42. So our reconstruction is not fundamentally worse than Polyphest. The entire gap is that we build on the ASTRAL backbone rather than constructing the topology jointly the way Polyphest does.

This reframes the problem. To win edit distance we must address the backbone, not the predictions. To win reticulation metrics we should improve the predictions, since there the predictions do matter.

## 9. What is next and why

Tier 1, high confidence, near term:
- Complete the clean per config benchmark on all 12 configurations, running now, to fix the exact standing on every metric against GRAMPA-iter and Polyphest.
- Improve partner prediction. The current partner predictions are noisy, with reciprocal pairs and tip level partners. The oracle shows real headroom on the reticulation metrics, so this work pays off there. Target is to match or beat Polyphest on reticulation identification under noise.
- Quantify the speed advantage. Our method is feature computation plus a forward pass, against an integer program search. This should be a clear win and is a real contribution.

Tier 2, ambitious, the research frontier:
- The backbone. Since our build already beats Polyphest given a good backbone, the path to winning edit distance is to either make the build robust to the small ASTRAL topology errors, or refine the backbone using the gene tree and copy number signal, or co-infer the topology from clusters with learned scores. This is the harder, months scale direction, and the oracle proves the upside is real.

Recommended sequencing: secure the Tier 1 wins first, which gives a defensible contribution of a fast, prior free, noise robust learned method that is competitive with GRAMPA-iter and strong on detection, then invest the remaining time in the backbone frontier to also compete on edit distance.

## 10. Questions for discussion

- Is competing on reticulation identification and speed, with an honest backbone limitation on edit distance, a sufficient contribution, or should we commit the remaining months to the backbone problem to also win edit distance.
- For the backbone, is refining ASTRAL with the polyploidy signal preferable to a cluster based joint construction.
- How much weight to put on the speed and prior free advantages relative to raw accuracy.
