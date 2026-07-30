# Method: a graph neural network for detect-then-place polyploid reconstruction

DRAFT for the thesis method section. Prose style follows the thesis convention of no
semicolons, no non-mathematical parentheses, and no em-dashes. Numbers and mechanics here
are code-verified. Empirical results belong in the results and diagnostic sections.

## Overview and rationale

The method reconstructs a polyploid phylogenetic network from a set of gene trees. It follows
a detect-then-place strategy. First it builds a single-copy species tree backbone from the gene
trees with ASTRAL. Then a graph neural network reads that backbone together with features
summarising the gene trees, and it makes two predictions on every backbone edge. The first
prediction is detection, namely whether a whole genome duplication occurred on that edge. The
second is placement, namely which other lineage is the second parent of the duplicated lineage.
The predicted events are then stamped onto the backbone to produce a multi-labeled tree, and
that tree is folded into a network.

Two properties motivate this design in the context of a benchmarking study. The method is fast,
because it reasons over one fixed backbone rather than searching tree space. It is also
prior-free, because it infers its own per-species copy number from the gene trees rather than
being supplied a ploidy vector. This makes it directly comparable to the prior-free iterative
GRAMPA baseline, and it removes an input that some competitors require.

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

The 9 edge features are computed per backbone branch, where each branch corresponds to the clade
of species below it. Four are structural. Concordance is the gene-tree support for the
species-level bipartition induced by the branch. The remaining three are the branch length from
the species tree, the clade size, and the depth.

The other five edge features are designed to detect whole genome duplication, which is one event
that duplicates an entire lineage at a single time. Duplication synchrony is the mean pairwise
correlation, across the clade species, of a per gene tree indicator of whether the species is
duplicated. A shared event makes the clade species duplicate together, so synchrony is high. The
mirrored-sister fraction is the fraction of gene trees in which a subset of the clade appears as
two identical sister subtrees, which is the signature of autopolyploidy. The copy-pair divergence
mean is the mean, over duplicated clade species and gene trees, of the branch-length distance
between the two closest copies of a species, normalised by the tree scale. This is a tree-based
analogue of the synonymous-site divergence used to date duplication events, so it estimates the
typical age of the duplication. The copy-pair divergence coefficient of variation measures
whether those ages are consistent, which is low for one shared event and high for scattered gene
duplications. The clade-duplicated fraction is the mean fraction of clade species that are
duplicated per gene tree.

## Model architecture

The network has around 82 thousand parameters and proceeds in five stages.

First, the leaf feature vectors are propagated to the internal nodes. Each internal node receives
the mean of the feature vectors of the leaves below it. This step has no learnable parameters and
is computed once per sample, so it is a fixed preprocessing operation rather than a learned
aggregation.

Second, the node features are projected to a hidden dimension of 64 by a two-layer perceptron.

Third, three graph attention layers refine the node representations. Each layer uses four
attention heads, a residual connection, and layer normalisation. Message passing over the tree
lets the representation of each edge endpoint be informed by the whole tree, so the per-edge
predictions that follow are made in global context rather than from local features alone.

Fourth, each branch is given an embedding. For a branch with parent node representation and child
node representation, the embedding is produced by a perceptron applied to the concatenation of the
parent representation, the child representation, and the 9 edge features of that branch.

Fifth, two prediction heads read the branch embeddings. The detection head is a perceptron that
maps each branch embedding to a binary decision, whether a whole genome duplication occurred on
that branch. The placement head scores, for each detected branch, every other branch as a
candidate second parent. For a source branch i and a candidate branch j, the score is produced by
a perceptron applied to the concatenation of the embedding of i, the embedding of j, and the
pairwise feature of the ordered pair i and j. Applying a softmax over the candidates gives a
distribution over parents. When the highest scoring candidate is the source branch itself the
event is an autopolyploidy, and otherwise it is an allopolyploidy whose second parent is the
chosen branch.

### Pairwise placement features

The pairwise feature for an ordered pair of branches carries the relational co-clustering signal
that the placement head needs and that no per-branch feature can express. For clade i and clade j
it is the mean and maximum, over the species in i crossed with the species in j, of the leaf-sister
co-clustering matrix. For a candidate that is a single species this reduces to the exact
co-clustering value between the two species. A copy-aware variant based on local neighbourhoods
was also implemented, but it did not improve accuracy, so the reported model uses the mean and
maximum co-clustering only.

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

Optimisation uses Adam with a learning rate of 0.001 and weight decay of 0.0001. Each species
tree is one graph, and gradients are accumulated over eight graphs before each optimiser step.

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
