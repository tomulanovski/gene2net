# PlaceNet

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Numbers and mechanics here are code-verified.
Empirical results belong in the results and diagnostic sections.

PlaceNet reconstructs a polyploid network in two stages. A single-copy species tree is first
inferred from the gene trees with ASTRAL-4 and taken as a fixed backbone. A graph neural network
then reads that backbone together with features summarizing the gene trees and predicts, for every
branch, whether a polyploidization event occurred on it and, for the branches it accepts, which
branch contributed the second subgenome. The accepted events are grafted onto the backbone,
yielding a MUL-tree, which is the output of the method. Bringing that MUL-tree to a network,
whether for scoring or for reporting, uses the same HOLM folding applied to every other method, so
folding is part of the comparison rather than part of PlaceNet.

## Input Representation

The backbone is represented as a graph whose nodes are the nodes of the species tree and whose
edges are its branches. Each leaf carries a vector of 13 features, each branch a vector of 9
features, and each ordered pair of branches a vector of 4 features used only by the placement head.
Internal nodes carry no observed features and are filled deterministically before inference, as
described below. All features were computed from the gene trees and the inferred backbone alone,
with no reference to ploidy levels, chromosome counts or any other external information.

Two terms recur throughout. A gene copy is a leaf of a gene tree, labeled by the taxon it belongs
to, so a taxon retaining two subgenomes contributes up to two copies per gene tree. The clade of a
branch is the set of taxa descended from that branch in the backbone, and because each branch is
identified with its clade, predicting on branches and predicting on clades are the same operation.

### Node Features

Eight of the 13 node features summarize copy number. For each taxon we counted its copies in every
gene tree, including the gene trees from which it is absent, which contribute a count of zero, and
summarized the resulting distribution by its mean, variance, mode and maximum and by the fraction
of gene trees in which the taxon was represented by zero, one, two, or three or more copies. A
taxon retaining both subgenomes of a polyploidization event is represented by two copies in most
gene trees, so these features carry the primary duplication signal.

The remaining five features summarize co-clustering. For a pair of taxa, the co-clustering fraction
is the fraction of all gene trees in which the two appear as sibling leaves, meaning that a copy of
each descends directly from a common parent node. The sibling must itself be a leaf, so a copy
whose sister is an internal node contributes nothing even when every taxon in that subtree is the
second taxon of the pair. Each gene tree contributes at most one to a pair however many copies
satisfy the condition, and the denominator is the complete set of gene trees rather than only those
containing both taxa. For a given taxon this yields one fraction per other taxon. Because the
number of taxa differs among datasets, we reduced that vector to five statistics independent of its
length, its mean, standard deviation, maximum, minimum and median. These describe how concentrated
the co-clustering of a taxon is and not which taxa it clusters with. The identity of a partner
lineage is instead supplied by the pairwise features described below.

### Branch Features

Four of the 9 branch features are structural. The concordance of a branch is the fraction of
informative gene trees that recover the split the branch induces on the taxa. Gene copies are
collapsed to taxon labels, a gene tree recovers the split when some branch of the gene tree
separates exactly the taxa of the clade restricted to the taxa present in that gene tree, with
either side of the gene tree branch counting, and a gene tree holding taxa from only one side of
the split is uninformative and was excluded from the fraction. Restricting the split to the taxa
present means that gene loss relaxes the requirement rather than failing it, whereas a clade
recovered only in part never counts. The branch length is the length of the branch in the ASTRAL
tree. The clade size is the number of taxa in the clade and the depth is the number of branches
between the root and the branch, the two together locating the branch within the tree.

The remaining five branch features target the signature of a polyploidization event on the branch.
Because such an event duplicates an entire clade at once, its trace in the gene trees is that the
taxa of the clade are duplicated together and that their duplicate copies diverged at a common
time. Copy number alone cannot localize the event, because a taxon duplicated by an event is
duplicated on every branch ancestral to it as well, and these five features are what separate a
branch carrying an event from its ancestors.

The clade-duplicated fraction is the mean, over the gene trees holding at least one taxon of the
clade, of the fraction of the clade taxa present in the gene tree that are represented by two or
more copies. Only the number of copies of a taxon enters, so a taxon whose two copies lie in
unrelated parts of the gene tree counts as duplicated in exactly the same way as one whose copies
are siblings. Because the denominator is the clade's own taxa, the fraction falls on the ancestors
of a duplicated clade, which also hold taxa the event did not duplicate, and this contrast is what
localizes an event upward along the tree.

Duplication synchrony is the mean Pearson correlation, over pairs of clade taxa, of the indicator
that a taxon is represented by two or more copies in a gene tree. The indicator was recorded for
every taxon of the clade in every gene tree holding at least one of them, with absent taxa recorded
as not duplicated, and pairs involving a taxon duplicated in all or in none of those gene trees
were omitted as uninformative. A single event duplicates the whole clade in the same gene trees and
drives the correlation up, whereas gene duplications striking taxa independently leave it near
zero.

The mirrored-sister fraction is the fraction of all gene trees holding a node whose two child
subtrees carry the same set of taxa, with that set contained in the clade. Only the taxon sets are
compared and not the topologies of the two subtrees, and a single mirrored taxon, two sibling
leaves carrying the same label, satisfies the condition in the same way as a mirrored subtree of
several taxa. This is the signature of an autopolyploidization, in which the duplicated clade is
grafted beside itself.

The copy-pair divergence is the mean distance between the two closest copies of a taxon, taken over
every combination of a clade taxon and a gene tree in which that taxon has at least two copies and
normalized by the mean root-to-leaf distance of the gene tree. The closest pair is used because,
where copies derived from an event coexist with copies from later gene duplications, the most
recent split is the best single proxy for the event. The measure is a tree-based analogue of the
synonymous-site divergence, or Ks, used to date duplications, so it estimates the age of the
duplicate copies. Its coefficient of variation over the same set is the fifth feature, low when the
copies of a clade diverged at a common time and high when they arose from independent gene
duplications.

### Pairwise Placement Features

For a clade the detection head has accepted, the placement head scores every branch as the
candidate source of the second subgenome, the branch itself included. No feature of either branch
alone can answer this, because whether a candidate is correct depends on the relation between the
two clades. Both pairwise signals rest on the same fact. An allopolyploidization leaves a copy of
the duplicated clade beside each of the two lineages that contributed to it, so in the gene trees
the copies of the clade recur beside the taxa of the true partner.

The first signal reuses the co-clustering fractions of the node features, summarized over a pair of
clades instead of per taxon. Over all pairs of taxa with one taxon drawn from each clade we took
the mean and the maximum of their co-clustering fractions. The mean asks whether the two clades sit
together in general and the maximum whether some one pair of taxa sits together very often. Both
are needed because the signal is uneven. The duplicated copy is grafted at a single position inside
the candidate clade rather than beside all of it, and discordance and fractionation leave only some
cross pairs recovering the sibling relation at all, so averaging over every cross pair dilutes what
the maximum retains.

The second signal is cluster support, which relaxes the requirement that a copy be the sibling of a
single taxon. A copy frequently lands beside a group of taxa from the candidate clade, in which
case it is the sibling of the group and of no single taxon, so no pair forms and the gene tree
contributes nothing to the first signal. Cluster support instead measures how much of the immediate
neighborhood of a copy comes from the candidate clade. The neighborhood of a copy is the largest
clade holding it whose leaf count does not exceed ten, with the copy itself excluded, and the cap is
needed because the neighborhood must extend beyond the duplicated clade before anything outside it
becomes visible. For each clade under consideration we retained only the copies of its taxa whose
neighborhood is less than half composed of taxa of that same clade, since a copy buried among its
own clade carries no information about where its duplicate went, and for each retained copy we
recorded the fraction of its neighborhood belonging to the candidate clade. The two features are the
sum of those fractions over all retained copies in all gene trees divided by the number of gene
trees, which is a rate and is therefore not bounded by one, and the largest single fraction among
them.

## Model Architecture

The network holds approximately 1.5 million parameters and proceeds in five stages. Its width,
depth and remaining hyperparameters are the values selected by the search described below.

The leaf feature vectors are first propagated to the internal nodes, each internal node receiving
the mean of the feature vectors of the leaves beneath it. This step has no learnable parameters and
is computed once per sample, so it is a fixed preprocessing operation rather than a learned
aggregation. The node features are then projected to a hidden dimension of 256 by a two-layer
perceptron.

Four graph attention layers next refine the node representations, each with four attention heads, a
dropout of 0.1, a residual connection and layer normalization. Each layer passes messages across
one branch, so four layers reach four branches out. Context beyond that range comes from the
propagation step, in which every internal node already holds the mean of the leaves beneath it, so
each node carries a summary of its whole subtree before any attention is applied.

Each branch is then given an embedding, produced by a two-layer perceptron applied to the
concatenation of the representation of its parent node, the representation of its child node, and
its 9 branch features. Because a branch embedding is formed from its two endpoint nodes, the
per-branch predictions inherit the context accumulated in those nodes rather than resting on the
branch features alone.

Two prediction heads finally read the branch embeddings. The detection head is a two-layer
perceptron mapping each branch embedding to the probability that a polyploidization event occurred
on that branch. The placement head scores, for an accepted branch i, every branch j of the tree as
the candidate source of the second subgenome, by applying a two-layer perceptron to the
concatenation of the embedding of i, the embedding of j, and the 4 pairwise features of the ordered
pair. A softmax over the candidates gives a distribution over parents. When the highest scoring
candidate is branch i itself the event is an autopolyploidization, and otherwise it is an
allopolyploidization whose second subgenome derives from the chosen branch.

## Training

Detection and placement were trained jointly. The total loss is a weighted sum of a detection and a
placement loss, and a single backward pass carries gradients from both terms into the shared layers,
so one set of node and branch representations serves both tasks.

Detection used a focal loss with class weighting, because the branches carrying an event are rare
relative to the rest. The focal term downweights each branch by a factor growing with the
confidence of its own prediction, so the difficulty of a branch is derived from the prediction
rather than labeled in advance, and the class weight raises the cost of missing an event relative to
the cost of a false positive.

Placement was trained one branch at a time. For a branch carrying a true event we scored every
branch of the tree as a candidate, the branch itself included, and a softmax over those scores gave
a distribution over candidates. The loss for that branch is the negative log of the probability
assigned to the correct candidate, which is the branch itself for an autopolyploidization and the
branch of the second parental lineage for an allopolyploidization. A network with n events
contributes n such losses, which were averaged. Branches carrying no event were excluded, having no
correct candidate.

Optimization used Adam with a learning rate of 0.0005 and a weight decay of 0.00001. Each species
tree is one graph, and gradients were accumulated over eight graphs before each optimizer step. Six
settings were tuned with Optuna on the validation split, the convolution type, the hidden dimension,
the number of layers, the dropout, the learning rate and the weight decay, in two rounds. The first
round converged on the edge of its range, at a hidden dimension of 256 and a depth of four, so the
second widened the range to a hidden dimension of 512 and a depth of five and confirmed that nothing
larger performed better. The remaining settings were held fixed, among them the batch size, the loss
weights and the constants of the features.

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
