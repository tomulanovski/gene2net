# Method: a graph neural network for detect-then-place polyploid reconstruction

DRAFT for the thesis method section. Prose style follows the thesis convention of no
semicolons, no non-mathematical parentheses, and no em-dashes. Numbers and mechanics here
are code-verified. Empirical results belong in the results and diagnostic sections.

This section describes PlaceNet in detail, its input representation, its architecture, and its
training. The output of PlaceNet is a MUL tree. Comparing that output to a phylogenetic network,
whether for scoring against the true network or for reporting, uses the same MUL-tree folding based
on the Holm algorithm that the evaluation applies to every method, so the folding is a property of
the comparison and not part of PlaceNet.

## Input representation

The species tree is represented as a graph whose nodes are the species tree nodes and whose
edges are the tree branches. Every leaf node carries a feature vector of 13 dimensions. Internal
nodes carry no observed features and are filled deterministically before the network runs, as
described below. Every branch carries a feature vector of 9 dimensions. A small number of
pairwise features, described with the placement head, are attached to ordered pairs of branches.

### Node features

The 13 node features are computed per species from the gene trees. Eight of them summarise copy
number. For each species we count, in each gene tree, how many leaves carry that species label,
and we summarise the resulting distribution by its mean, variance, mode, maximum, and the
fraction of gene trees in which the species is absent, present once, present twice, or present
three or more times. A species that underwent a whole genome duplication tends to appear as two
copies, so these features carry the primary duplication signal.

The remaining five node features summarise co-clustering. For a target species and every other
species we measure the fraction of gene trees in which a copy of the target is a direct
leaf-sister of that other species. This yields one value per other species. Because the number
of species varies, we reduce this vector to a fixed summary of five statistics, namely its mean,
standard deviation, maximum, minimum, and median. This summary describes how concentrated a
species co-clusters, not with whom, so it behaves as a status descriptor that separates a clean
single-ancestry lineage from one with divided ancestry.

### Edge features

The nine edge features are computed per backbone edge, where each edge corresponds to the clade of
species below it. Four of them are structural. The concordance factor is the fraction of the
informative gene trees that recover the species-level bipartition the edge induces, where a gene
tree is informative for an edge when it has species on both sides of that bipartition. The branch
length is the length of the edge in the ASTRAL species tree. The clade size is the number of species
below the edge, and the depth is the number of edges from the root down to it. The clade size and
depth give the model context for where in the tree an edge sits.

The other five edge features capture the signal of a whole genome duplication on the edge. Like the
other features they feed both prediction heads through the shared trunk, but these are the ones that
measure the duplication itself. Because a whole genome duplication copies an entire clade at once,
its trace in the gene trees is that the species below the edge are duplicated together and to a
common age. Duplication synchrony measures whether the clade species tend to be duplicated in the
same gene trees, which a single shared event produces and scattered duplications do not. The
clade-duplicated fraction is the mean, per gene tree, of the fraction of the clade species present in
that tree that are duplicated, which measures how pervasive the duplication is. The mirrored-sister
fraction is the fraction of gene trees in which part of the clade appears as two identical sister
subtrees, the signature of an autopolyploidy. The copy-pair divergence is the mean distance between a
species' two closest copies, taken over the duplicated clade species and gene trees and normalised by
the tree scale. It is a tree-based analogue of the synonymous-site divergence, or Ks, used to date
duplications, so it estimates the typical age of the copies. The copy-pair divergence coefficient of
variation measures how consistent those ages are, low for a single shared event and high for
independent gene duplications.

## Model architecture

The network has about 1.5 million parameters and proceeds in five stages. Its width, depth, and the
other hyperparameters below are the values selected by the hyperparameter search described in the
training section.

First, the leaf feature vectors are propagated to the internal nodes. Each internal node receives
the mean of the feature vectors of the leaves below it. This step has no learnable parameters and
is computed once per sample, so it is a fixed preprocessing operation rather than a learned
aggregation.

Second, the node features are projected to a hidden dimension of 256 by a two-layer MLP.

Third, four graph attention layers refine the node representations. Each layer uses four attention
heads, a dropout of 0.1, a residual connection, and layer normalisation. The attention passes
messages between nodes, so after these layers each node representation carries information from
across the whole tree. Because the edge embeddings built in the next stage are formed from their two
endpoint nodes, the per-edge predictions inherit that global context indirectly, through the nodes,
rather than resting on the local edge features alone.

Fourth, each branch is given an embedding. For a branch with parent node representation and child
node representation, the embedding is produced by a two-layer MLP applied to the concatenation of the
parent representation, the child representation, and the 9 edge features of that branch.

Fifth, two prediction heads read the branch embeddings. The detection head is a two-layer MLP that
maps each branch embedding to a binary decision, whether a whole genome duplication occurred on
that branch. The placement head scores, for each detected branch, every other branch as a
candidate second parent. For a source branch i and a candidate branch j, the score is produced by
a two-layer MLP applied to the concatenation of the embedding of i, the embedding of j, and the
pairwise features of the ordered pair i and j. Applying a softmax over the candidates gives a
distribution over parents. When the highest scoring candidate is the source branch itself the
event is an autopolyploidy, and otherwise it is an allopolyploidy whose second parent is the
chosen branch.

### Pairwise placement features

For a clade the detection head has flagged as duplicated, the placement head scores every other edge
as the candidate second parent, the edge onto which a copy of that clade is grafted. This choice is
relational. It depends on how the duplicated clade sits relative to a candidate parent rather than on
either edge alone, so it needs a signal defined on the ordered pair of the source clade and the
candidate, which no single-edge feature can provide. That signal is co-clustering, the frequency
across the gene trees with which a species from the source clade appears as a sister to a species from
the candidate. It is informative because an
allopolyploid grafts a copy of its clade next to its second parent, so the two co-cluster there. There
are four pairwise features, in two pairs. The first pair is the mean and the maximum, across the
species of the two clades, of this co-clustering. The second pair is a copy-aware cluster-support
signal. When a species is duplicated, one copy tends to stay home among its own clade while the other
lands away near its allopolyploid partner. For each such away copy we measure the fraction of its
local neighborhood in the gene tree that belongs to the candidate partner clade, and we summarise that
fraction over the gene trees as an intensity and as a per-copy peak. This second pair is the sharper
allopolyploid signal, because it reads the neighborhood of the copy that actually moved rather than
symmetric sisterhood alone.

## Training objective

Detection and placement are trained jointly on synthetic species networks with known events. The
total loss is a weighted sum of a detection loss and a placement loss, and a single backward pass
sends gradients from both terms into the shared trunk of the network, so the refined
representations serve both tasks at once.

Detection uses a focal loss with class weighting, because whole genome duplication branches are
rare relative to non-event branches. The focal term downweights easy branches by a factor that
grows with the model confidence, so the effective difficulty of each branch is derived from the
prediction rather than being labeled in advance.

Placement uses a cross-entropy loss over the candidate branches, and it is applied only on
branches that carry a true event. The target for an autopolyploidy is the branch itself, and the
target for an allopolyploidy is the second parent branch.

Optimisation uses Adam with a learning rate of 0.0005 and weight decay of 0.00001. Each species
tree is one graph, and gradients are accumulated over eight graphs before each optimiser step. The
hyperparameters given here and in the architecture, namely the hidden dimension, the number of
layers, the learning rate, and the weight decay, were selected by a two-stage hyperparameter search
on the validation split.

## From predictions to a network

The predicted events are turned into a network in three steps. First, detection produces a
probability per branch. Branches are selected in order of confidence, and a branch is accepted
only while it keeps every species within its inferred copy bound. The copy bound is the largest
copy count that at least half of the gene trees support for that species, which is a consensus
estimate rather than a supplied ploidy. This selection removes the surplus branches that the
recall-oriented detector produces.

Second, each selected event is stamped onto a copy of the ASTRAL backbone. For an autopolyploidy
the clade is duplicated as an identical sibling. For an allopolyploidy a copy of the clade is
grafted at the predicted second parent, while the original copy stays at its backbone position.
After all events are stamped, some species labels appear more than once, which is what makes the
result a multi-labeled tree.

Third, the multi-labeled tree is folded into a network. Identical sibling duplicates collapse and
contribute ploidy without a reticulation, which matches autopolyploidy. Duplicates that sit under
two different parents become reticulation edges, which matches allopolyploidy.

A structural consequence of the stamping step is worth stating here, because it frames the
limitations discussed later. The build always starts from the ASTRAL backbone and only adds
duplicated clades. It never re-arranges the existing branches. Every topological error in the
backbone is therefore inherited by the multi-labeled tree and by the final network. For an
allopolyploidy in particular, the original copy remains at the backbone position that ASTRAL
assigned, so only the second parent is ever chosen by the model.
