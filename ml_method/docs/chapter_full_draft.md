# A graph neural network for detect-then-place polyploid network reconstruction

CONSOLIDATED FULL DRAFT of the method chapter, assembled from the section drafts. Prose style
follows the thesis convention of no semicolons, no non-mathematical parentheses, and no
em-dashes. Numbers and mechanics are code-verified. Two kinds of editorial marks appear:

- `[PENDING HPO]` marks a number or subsection that the running hyperparameter study will fix
  (the final architecture, the tuning section, the learning curve, and the benchmark table).
- `[TODO]` marks a data-provenance value to confirm from the simulation configuration.

Assembly decisions worth flagging to the author:
- The copy-aware cluster-support pairwise feature is described as INCLUDED in the model, because
  it is (verified in `build_pairwise_feat`, and the shipped and HPO-base config use it). If the
  decision is to remove it, the HPO must restart on a two-channel base and this section changes.
- The central framing follows the corrected canonical diagnostic: the backbone is a moderate,
  condition-dependent contributor and event prediction is the larger, improvable lever. The
  earlier "backbone is a 0.54 wall" language from the older limitation drafts has been removed as
  an artifact of the uncorrected edit-distance measure.

---

## 1. Overview and rationale

The method reconstructs a polyploid phylogenetic network from a set of gene trees. It follows a
detect-then-place strategy. First it builds a single-copy species tree backbone from the gene
trees with ASTRAL. Then a graph neural network reads that backbone together with features
summarising the gene trees, and it makes two predictions on every backbone edge. The first
prediction is detection, namely whether a whole genome duplication occurred on that edge. The
second is placement, namely which other lineage is the second parent of the duplicated lineage.
The predicted events are then stamped onto the backbone to produce a multi-labeled tree, and that
tree is folded into a network.

Two properties motivate this design in the context of a benchmarking study. The method is fast,
because it reasons over one fixed backbone rather than searching tree space. It is also
prior-free, because it infers its own per-species copy number from the gene trees rather than
being supplied a ploidy vector. This makes it directly comparable to the prior-free iterative
GRAMPA baseline, and it removes an input that some competitors require.

## 2. Experimental setup

### Data

The method is trained and evaluated on simulated data, so that the true network is known and the
reconstruction can be scored exactly. Species networks are simulated across two sources of
gene-tree discordance. The first is incomplete lineage sorting, varied at low, medium, and high
levels. The second is gene duplication and loss, varied at low, medium, and high rates at a fixed
effective population size. This gives six configurations that span the regimes under which a
detect-then-place method is expected to succeed and to struggle. [TODO confirm the number of taxa
per network, the number of gene trees per network, and the replicate structure.]

Training uses a large set of synthetic networks simulated under the same six configurations. [TODO
confirm count, on the order of two thousand networks per configuration.] For each network the gene
trees are simulated, the ASTRAL backbone is inferred, and the true events are mapped to backbone
edges to produce the training labels. The species tree is re-rooted before feature extraction and
labeling, because ASTRAL output is unrooted and an arbitrary root inflates the edit distance and
misaligns the labels. A fixed random split holds out one fifth of the networks for validation, and
the same split is used across all experiments so that validation numbers are comparable.

Evaluation uses a separate benchmark of 21 networks, disjoint from the training networks, with
gene trees simulated under each of the six configurations. This benchmark is the primary test of
generalization, because it measures reconstruction of networks the method never saw during
training.

### Metrics

Reconstructions are scored against the true network by six distance measures, all oriented so that
lower is better. The multi-labeled-tree edit distance is the primary measure of overall
reconstruction quality. The Robinson-Foulds distance measures the topological error of the
backbone. The reticulation count difference measures whether the right number of reticulations was
inferred. The reticulation-leaf and reticulation-sister Jaccard distances measure whether the
reticulations involve the right lineages, and they are the measures most directly controlled by
the placement head. The ploidy difference measures whether each species is assigned the correct
copy number.

The multi-labeled-tree edit distance is computed with a canonical child ordering, described in the
diagnostic section, so that it does not depend on the arbitrary Newick child order. All edit
distances in this chapter use that canonical measure.

### Baselines

The method is compared against two published approaches with different information demands. The
first is Polyphest, which is supplied the true ploidy of each species. The second is an iterative
application of GRAMPA, which like the present method infers ploidy from the data rather than
receiving it. The comparison to GRAMPA is therefore the fair prior-free comparison, because both
methods operate from the same inputs, while the comparison to Polyphest is informative but favors
Polyphest by an input it is given and the present method is not.

## 3. Method

### 3.1 Input representation

The species tree is represented as a graph whose nodes are the species tree nodes and whose edges
are the tree branches. Every leaf node carries a feature vector of 13 dimensions. Internal nodes
carry no observed features and are filled deterministically before the network runs, as described
below. Every branch carries a feature vector of 9 dimensions. A small number of pairwise features,
described with the placement head, are attached to ordered pairs of branches.

#### Node features

The 13 node features are computed per species from the gene trees. Eight of them summarise copy
number. For each species we count, in each gene tree, how many leaves carry that species label,
and we summarise the resulting distribution by its mean, variance, mode, maximum, and the fraction
of gene trees in which the species is absent, present once, present twice, or present three or more
times. A species that underwent a whole genome duplication tends to appear as two copies, so these
features carry the primary duplication signal.

The remaining five node features summarise co-clustering. For a target species and every other
species we measure the fraction of gene trees in which a copy of the target is a direct
leaf-sister of that other species. This yields one value per other species. Because the number of
species varies, we reduce this vector to a fixed summary of five statistics, namely its mean,
standard deviation, maximum, minimum, and median. This summary describes how concentrated a
species co-clusters, not with whom, so it behaves as a status descriptor that separates a clean
single-ancestry lineage from one with divided ancestry.

#### Edge features

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
typical age of the duplication. The copy-pair divergence coefficient of variation measures whether
those ages are consistent, which is low for one shared event and high for scattered gene
duplications. The clade-duplicated fraction is the mean fraction of clade species that are
duplicated per gene tree.

### 3.2 Model architecture

The network proceeds in five stages. `[PENDING HPO: the hidden dimension, the number and type of
message-passing layers, and the total parameter count are set by the hyperparameter search in
Section 4. The description below states the structure; the selected sizes are filled from the
study winner.]`

First, the leaf feature vectors are propagated to the internal nodes. Each internal node receives
the mean of the feature vectors of the leaves below it. This step has no learnable parameters and
is computed once per sample, so it is a fixed preprocessing operation rather than a learned
aggregation.

Second, the node features are projected to the hidden dimension by a two-layer perceptron.

Third, a stack of graph message-passing layers refines the node representations. Each layer uses a
residual connection and layer normalisation. Message passing over the tree lets the representation
of each edge endpoint be informed by the whole tree, so the per-edge predictions that follow are
made in global context rather than from local features alone. The layer type, the width, and the
depth are selected by the hyperparameter search.

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

#### Pairwise placement features

The pairwise feature for an ordered pair of branches carries the relational signal that the
placement head needs and that no per-branch feature can express. It has four channels. The first
two are the mean and maximum, over the species in clade i crossed with the species in clade j, of
the leaf-sister co-clustering matrix. For a candidate that is a single species this reduces to the
exact co-clustering value between the two species. The other two are a copy-aware cluster-support
feature, namely a support intensity and a peak support computed from the local neighbourhoods of
the duplicated copies in the gene trees, which sharpens the allopolyploid signal. The feature
importance analysis in Section 8 shows the trained model draws on both the co-clustering channels
and the cluster-support channels for placement.

### 3.3 Training objective

Detection and placement are trained jointly on synthetic species networks with known events. The
total loss is a weighted sum of a detection loss and a placement loss, and a single backward pass
sends gradients from both terms into the shared trunk of the network, so the refined
representations serve both tasks at once.

Detection uses a focal loss with class weighting, because whole genome duplication branches are
rare relative to non-event branches. The focal term downweights easy branches by a factor that
grows with the model confidence, so the effective difficulty of each branch is derived from the
prediction rather than being labeled in advance.

Placement uses a cross-entropy loss over the candidate branches, and it is applied only on
branches that carry a true event. The target for an autopolyploidy is the branch itself. The
target for an allopolyploidy is the away parent, namely the true parent that is not the lineage
ASTRAL placed the polyploid beside. The reason this specific target is used, rather than the
parent drawn as the reticulation in the true network, is a labeling defect analysed in Section 6.

Optimisation uses Adam. `[PENDING HPO: the learning rate and weight decay are selected by the
hyperparameter search in Section 4.]` Each species tree is one graph, and gradients are accumulated
over eight graphs before each optimiser step.

### 3.4 From predictions to a network

The predicted events are turned into a network in three steps. First, detection produces a
probability per branch. Branches are selected in order of confidence, and a branch is accepted
only while it keeps every species within its inferred copy bound. The copy bound is the largest
copy count that at least half of the gene trees support for that species, which is a consensus
estimate rather than a supplied ploidy. This fill-to-inferred-ploidy selection reaches each
species' inferred copy number from the ranked candidate events, rather than thresholding on a
fixed detection probability, so a moderately confident true event is still used and the surplus
that a recall-oriented detector produces is removed.

Second, each selected event is stamped onto a copy of the ASTRAL backbone. For an autopolyploidy
the clade is duplicated as an identical sibling. For an allopolyploidy a copy of the clade is
grafted at the predicted second parent, while the original copy stays at its backbone position.
After all events are stamped, some species labels appear more than once, which is what makes the
result a multi-labeled tree.

Third, the multi-labeled tree is folded into a network. Identical sibling duplicates collapse and
contribute ploidy without a reticulation, which matches autopolyploidy. Duplicates that sit under
two different parents become reticulation edges, which matches allopolyploidy.

A structural consequence of the stamping step is worth stating here, because it frames the
limitations. The build always starts from the ASTRAL backbone and only adds duplicated clades. It
never re-arranges the existing branches. Every topological error in the backbone is therefore
inherited by the multi-labeled tree and by the final network. For an allopolyploidy in
particular, the original copy remains at the backbone position that ASTRAL assigned, so only the
second parent is ever chosen by the model.

## 4. Hyperparameter selection

`[PENDING HPO: this section is written from the running Optuna study. It will report the search
space (graph-convolution type in gat, gin, gcn; hidden dimension in 64, 128, 192, 256; depth in 1
to 4; learning rate; dropout; weight decay), the bounded budget of about 100 trials with a
Bayesian sampler and no pruner, the objective, which is validation allopolyploid partner accuracy
subject to a detection-F1 guard, and the selected configuration with its seed-confirmed
performance. Early trials show the selected model at the high-capacity, low-regularisation corner
of the space, with a large gain in partner accuracy over the untuned model, to be confirmed on the
benchmark.]`

## 5. Results and comparison to baselines

`[PENDING HPO: the benchmark table below is from the untuned model and will be regenerated with the
final tuned model. The qualitative findings are expected to hold or strengthen.]`

### Detection

The whole genome duplication detector reaches an F1 of about 0.80 on the held-out validation
split. [TODO confirm precision and recall from the training log.] Detection is deliberately
recall-oriented, because the fill-to-inferred-ploidy build reaches each species' inferred copy
number from the ranked candidate events rather than from a hard probability threshold, so a
moderately confident true event is still used.

### Reconstruction benchmark

We compare against the baselines on the 21 benchmark networks across all six configurations. All
methods except Polyphest infer ploidy from the gene trees. The table reports mean canonical edit
distance. GNN scores are over all 21 networks. The competitor columns are over the 15 to 18
networks each competitor completes.

| Configuration | GNN | Polyphest | GRAMPA | GRAMPA-iter |
| --- | --- | --- | --- | --- |
| ils_low | 0.573 | 0.065 | 0.682 | 0.681 |
| ils_medium | 0.571 | 0.112 | 0.695 | 0.720 |
| ils_high | 0.554 | 0.209 | 0.709 | 0.825 |
| dup_loss_low | 0.570 | 0.196 | 0.700 | 0.729 |
| dup_loss_medium | 0.564 | 0.177 | 0.707 | 0.756 |
| dup_loss_high | 0.603 | 0.636 | 0.716 | 0.812 |

Four findings.

First, the method beats both GRAMPA and iterative GRAMPA in every configuration, by roughly 0.1 to
0.25 in edit distance. Among the prior-free methods that reconstruct a full network from gene
trees, the method is clearly better than the GRAMPA family across all conditions.

Second, Polyphest is the strongest method at low and moderate difficulty, and by a large margin.
This traces to copy number. The edit distance on multi-labeled trees is dominated by getting each
species' copy number right, and Polyphest builds its start tree from the inferred consensus
multiset, which on clean simulations is essentially the true ploidy. So on the easy conditions
Polyphest is effectively handed the copy-number structure and only has to place events, which is
why it reaches an edit distance as low as 0.065.

Third, that advantage erodes as duplication and loss rise and reverses at the hardest condition.
Under high duplication and loss the copy-number signal is corrupted, the inferred multiset stops
matching the truth, and Polyphest's start tree is wrong. There the method reaches 0.603 against
Polyphest's 0.636, so it is the better method, and the reticulation-leaf distance agrees at that
condition, with the method at 0.264 against Polyphest's 0.412. The method reaches each species'
inferred ploidy from its ranked events rather than from a fixed start tree, so it degrades more
gracefully when ploidy inference is hard.

Fourth, the comparison understates the method's robustness. Polyphest is scored on 15 to 18 of the
networks per configuration because it does not complete the rest, while the method completes all
21. So Polyphest's means are taken over the subset it can finish, and the method always produces an
answer.

On the reticulation distances the ordering matches the edit distance. Polyphest leads on the
networks it completes at low and moderate difficulty, the method beats iterative GRAMPA everywhere
and GRAMPA on the reticulation-leaf distance in most configurations, and the method overtakes
Polyphest at high duplication and loss.

### Runtime

The method reconstructs a network in a few seconds. Its per-network compute, excluding the ASTRAL
step, is about 8 seconds on the benchmark, of which the graph neural network forward pass is a
fraction of a second and the rest is feature extraction from the gene trees. ASTRAL adds a few
seconds. Polyphest, by contrast, ranges from seconds on easy networks to days on hard ones, and
does not complete the networks it is missing from the table above. [TODO insert a concrete
Polyphest wall time or range from the run logs.] So the method offers a bounded, predictable
runtime of seconds per network, against a search-based competitor whose runtime is variable and
can be prohibitive.

## 6. Diagnostic: where the edit distance comes from

All edit distances here use the canonical multi-labeled-tree edit distance. Numbers are the oracle
on the validation split, n = 200 per configuration.

### The question

A reconstruction can lose edit distance for three reasons. The build and scoring themselves leave a
residual even on a perfect input. The backbone can be wrong, because the events are placed on the
ASTRAL species tree and ASTRAL misplaces lineages. And the events can be wrong, namely the wrong
branches are flagged, the wrong copy number is produced, or the wrong parent is chosen. This
section separates the three so we know which one to work on.

### The oracle experiment

We remove event error by feeding the true events through the build, and we vary only the backbone.
In the first setting the true events are stamped onto the true species tree, so everything is
correct and only the build residual remains. We call this the floor. In the second setting the same
true events are stamped onto the ASTRAL backbone, so the only added error is the backbone. We call
this the ceiling, because it is the best any predictor could reach on the ASTRAL backbone. The
difference between the two is the backbone's contribution, and whatever a real method loses beyond
the ceiling is its event-prediction error.

### Results

| Configuration | Floor (true backbone) | Ceiling (ASTRAL backbone) | Backbone gap |
| --- | --- | --- | --- |
| ils_low | 0.109 | 0.251 | 0.142 |
| ils_medium | 0.109 | 0.310 | 0.201 |
| ils_high | 0.109 | 0.391 | 0.282 |
| dup_loss_low | 0.109 | 0.311 | 0.202 |
| dup_loss_medium | 0.109 | 0.341 | 0.232 |

Two structural facts stand out. The floor is exactly 0.109 in every configuration. This is
expected, because the floor uses the true species tree and the true events, which are the same
underlying networks across configurations, and only the gene trees change. So 0.109 is the
irreducible residual of the build and the scoring convention, not a modelling error. The ceiling
rises from 0.251 at low incomplete lineage sorting to 0.391 at high, and sits near 0.31 to 0.34
under gene duplication and loss. The ceiling depends on the gene trees only through ASTRAL, so its
rise is the degradation of ASTRAL as discordance grows.

### Interpretation

The backbone is a moderate contributor, not a wall. Placing true events on the ASTRAL backbone
still reconstructs well, with a ceiling between 0.25 and 0.39 depending on conditions, against a
floor of 0.109. The backbone's contribution is the gap between them, from 0.14 at low incomplete
lineage sorting to 0.28 at high. It grows with discordance, as expected, but it does not dominate
the edit distance, and on the easier conditions half of the networks reach the floor exactly even
on the ASTRAL backbone.

The larger source of error is event prediction. A real method must also detect the right branches,
produce the right copy number, and choose the right parents, and whatever it loses beyond the
ceiling is that event error. Our method sits well above the ceiling, so event prediction, not the
backbone, is where most of its edit distance is lost. This is encouraging, because event
prediction is improvable within the present design. Filling each species to its inferred copy
number rather than capping at a detection threshold recovers events the threshold suppressed and
improves copy number and reticulation together. Better placement of the second parent is the
remaining lever, and it is the partner problem analysed next.

### A note on the metric

The numbers above required correcting the edit-distance measure. The multi-labeled-tree edit
distance is computed by a greedy graph-edit-distance search, and that search visits nodes in the
order the tree was written, so its result depended on the phylogenetically meaningless Newick child
order. Two identical trees could score far from zero when their child orders happened to misalign.
Ordering both trees canonically before the search removes this dependence. Under the uncorrected
measure the ASTRAL ceiling appeared to be about 0.54 and the backbone looked like the dominant
wall. Under the corrected canonical measure the ceiling is 0.25 to 0.39, which reversed that
reading and is why this chapter treats event prediction, not the backbone, as the dominant lever.

## 7. Partner prediction: a labeling defect and its repair

Partner accuracy is the fraction of allopolyploid events for which the predicted parent edge
matches the target. Autopolyploid events, where the parent is the event branch itself, are excluded
from this figure because they are easy and inflate the combined number.

### The question

The diagnostic showed that event prediction, not the backbone, is the dominant source of edit
distance. Within event prediction the choice of the allopolyploid parent is the weakest part on its
own terms, and it is worth locating. The untuned one-parter method reaches an allopolyploid partner
accuracy of only about 0.45. This section asks whether that number can be raised, first by better
features and then by a correction to the training target, and what the answer implies.

### Better hand-crafted features do not help

We added hand-crafted partner features to the model and retrained the one-parter method from
scratch on all six configurations. One restricts the co-clustering signal to gene trees in which
the target is actually duplicated, on the reasoning that single-copy trees carry no partner signal
and only dilute it. Another summarizes how concentrated a species co-clusters, on the reasoning
that an allopolyploid with one clear partner looks different from an autopolyploid spread over a
clade. A branch-length distance feature was also tried.

| Variant | Allopolyploid partner accuracy | Detection F1 |
| --- | --- | --- |
| Baseline | 0.453 | 0.775 |
| Duplication-conditioned co-clustering | 0.453 | 0.784 |
| Effective-partner-count | 0.466 | 0.788 |
| Both | 0.460 | 0.777 |

The duplication-conditioned feature is a flat null. The effective-partner-count feature adds 0.013,
within the sampling noise of a single run. Across these hand-crafted additions the number does not
leave the neighborhood of 0.45. Adding further hand-crafted features on top of the co-clustering
and cluster-support channels the model already uses is not the lever.

### Where the error comes from

The reason features saturate is visible when the signal is examined directly. Ranking the candidate
parents of an allopolyploid by co-clustering with the target in the gene trees places the true
partner first only 0.42 of the time, but within the top two candidates 0.86 of the time at low
incomplete lineage sorting. The information that identifies the partner is therefore largely
present. The failure is that a second candidate outranks the true partner.

We identified that second candidate. In the 0.58 of events where the true partner is not first, it
is almost always second, and the candidate that beats it is a species from a different part of the
tree in 0.89 of cases, with a mean co-clustering to the true partner of only 0.069. It is a diploid
in 0.64 of cases. This is the home parent, namely the other lineage that contributed a subgenome.
Both subgenomes leave a co-clustering trace, so the two strongest candidates are the two parents,
and the method must return one of them. The home parent, whose copy sits cleanly in its own clade,
frequently co-clusters more strongly than the away partner the label asks for.

### A labeling defect

Examining which of the two parents the label names revealed a defect in the training target itself.
The label is taken from the true network, where one parent is drawn as the base-tree attachment and
the other as a reticulation edge, and the reticulation parent becomes the label. This choice is
made without reference to where ASTRAL placed the polyploid. ASTRAL places the polyploid next to
whichever subgenome is stronger, which is not correlated with the network drawing, so the two
disagree about half the time. A direct audit against ground truth confirms it. In 0.566 of
single-species allopolyploid events at low incomplete lineage sorting and 0.539 at high duplication
and loss, the labelled partner is the same lineage ASTRAL placed the target next to. In those
events the target is asked to reticulate to the parent it is already attached to, which is
degenerate, and the correct partner is the other parent.

The audit also bounds how often the defect is fixable. ASTRAL places a single-species target next
to a true parent in about 0.94 of events, so the correct partner, the parent that is not the
target's home, is well-defined almost always. Only about 0.06 of targets are placed away from both
parents, which is a backbone error that no partner choice can repair.

### The repair and its effect

The repair is to recompute the target as the parent that is not the polyploid's ASTRAL home. When
the ASTRAL home is the labelled partner, the target is retargeted to the other true parent.
Otherwise it is left unchanged. The true tree is used only to construct the corrected label. At
inference the model sees only the ASTRAL tree, the gene trees, and the features, so it learns the
rule that the partner is the parent it was not placed beside, which is inferable from its inputs.

Training the one-parter twice under identical settings, differing only in the labels, isolates the
effect. Partner accuracy rises from 0.456 to 0.526, an increase of 0.070 that is many standard
errors beyond noise, while detection F1 is unchanged, which confirms the change is confined to the
partner target. `[PENDING HPO: with the tuned model the partner accuracy is substantially higher
again, to be reported with the final configuration.]`

The gain also carries to the reconstruction, measured against the true networks on the 21 benchmark
networks across all six configurations. Ploidy distance improves in all six configurations, from a
mean of about 0.455 to about 0.395. Edit distance improves in five of six. The reticulation-sister
and reticulation-leaf distances are mixed across configurations and flat on average. So the
corrected target recovers polyploidy structure more reliably and edit distance slightly, while the
finer reticulation metrics move less, which is consistent with the residual sitting where the
backbone is incoherent.

### Coverage and its bound

The repair as described covers single-species allopolyploid events. A clade-level event, where a
group of species shares one hybridization, is handled by the same rule applied to the clade's common
ancestor, but only when the clade is monophyletic in the ASTRAL tree, so that its home is a single
well-defined sibling. This holds for about 0.60 of clade events at low incomplete lineage sorting
and about 0.42 at high duplication and loss. The remaining clade events are scattered by ASTRAL and
left unchanged, because for a scattered clade there is no single home and therefore no well-defined
away parent to assign. This is not an unrepaired bug. It is the absence of a correct target, and it
arises because the backbone placed the clade incorrectly in the first place.

## 8. Feature importance

We rank features by permutation importance. For each feature we shuffle its values across the
validation examples, breaking its association with the labels while keeping its marginal
distribution, rerun the trained model, and measure the drop in a target metric. A large drop means
the feature matters to the trained model. We do this twice with two targets, detection F1 and
allopolyploid partner accuracy, because the method has two heads with different needs. `[PENDING
HPO: these rankings are measured on the untuned model. The qualitative story is expected to hold;
the numbers will be refreshed on the final model.]`

### Detection

Detection is carried by a single feature and a short tail. The fraction of the clade duplicated per
gene tree dominates, with an F1 drop of 0.38, roughly seven times the next feature. The rest of the
signal is copy number and clade size, namely mean copy number, clade size, mirrored-sister
fraction, and copy-pair divergence, each with a drop below 0.05. The features designed for whole
genome duplication detection are exactly the important ones, which is direct evidence for that
design choice. The co-clustering node summaries contribute essentially nothing to detection, which
is expected, because co-clustering is a partner signal, not a detection signal.

### Partner

Partner prediction draws on a broader set of features, and the ranking is different. The pairwise
co-clustering feature dominates, with partner-accuracy drops of 0.33 for its mean channel and 0.23
for its max channel. This is the direct signal the placement head consumes. Below it the copy-aware
cluster-support pairwise feature contributes a drop of about 0.15, so the model does draw on it. The
structural and copy features act through the shared trunk and the edge embeddings, so clade size,
mean copy number, concordance, depth, and copy-pair divergence each contribute drops between 0.11
and 0.19.

### What earns its place, and what is prunable

Reading the two rankings together gives a clean account. The detection features earn their place
through detection, led by the clade-duplicated fraction. The pairwise co-clustering features earn
their place through partner, as the dominant partner signal, and the cluster-support pairwise
feature earns a smaller place there. The copy-number and structural features earn their place
through both heads. Branch length and duplication synchrony contribute nothing to either head and
are prunable. The co-clustering node summaries contribute little to either head. One feature is
worth singling out. Depth is inert for detection but contributes a partner drop of 0.13, so it
earns its place through the partner head. We report these rankings but retain the full feature set,
because the low-importance features contribute little by construction and pruning them would require
retraining for a negligible expected change.

## 9. Learning curve and data sufficiency

`[PENDING HPO: this section reports the learning curve of the final model, trained on 10, 25, 50,
75, and 100 percent of the training pool with a fixed validation split. It will state whether
detection and partner accuracy have saturated in the available data or would benefit from more.
Preliminary curves on the untuned model suggest detection saturates early while partner accuracy is
still improving with data.]`

## 10. Limitations and future work

### The backbone is the binding structural limitation, at high discordance

The reconstruction is built by stamping events onto a fixed ASTRAL backbone and never rearranges
that backbone, so every backbone error is inherited by the final network. The diagnostic quantifies
this. Under an oracle that supplies the true events, placing them on the ASTRAL backbone leaves an
edit distance between 0.25 and 0.39 depending on conditions, against 0.11 on the true backbone. The
backbone contribution is therefore moderate at low discordance and grows to about 0.28 at high
incomplete lineage sorting. It is not a wall that dominates every condition, but it is the one part
of the pipeline the present design cannot repair, because the build does not rebuild the backbone.

### Placement shares the same root cause

Placement of the second parent looked like an independent weakness, but it reduces to the same
cause. Hand-crafted features do not improve it. A defect in the training target does, and repairing
the target raises partner accuracy and improves ploidy and edit distance on the benchmark, though
the improvement in the finer reticulation measures is modest. The reason is that the correct parent
is well-defined only when ASTRAL places the polyploid coherently, and undefined when it does not.
The residual placement error and the coverage bound of the clade repair both trace to backbone
coherence. When ASTRAL scatters a polyploid or a polyploid clade, there is no coherent home against
which to name the second parent, and the event cannot be repaired by any change to the target or the
features.

### The root cause and the path through it

Both the backbone gap at high discordance and the placement coverage bound share one mechanical
cause. The two copies of a polyploid are not assigned to their subgenomes. Unassigned copies pull
ASTRAL in two directions at once, which scatters the backbone, and the same unassigned copies leave
two co-clustering peaks without telling the method which peak is the home and which is the second
parent. Assigning the copies to subgenomes, which is phasing, addresses both at once. A phased copy
defines the backbone position of the polyploid, and the other phased copy names the second parent.
Phasing is therefore the principled next step, and it is a different method rather than a change to
the present one, because it operates before the backbone is built rather than decorating a backbone
after the fact. A learned two-parent placement head that predicts both parents jointly was attempted
as a lighter alternative and did not improve the benchmark, which further supports locating the fix
in the copy assignment rather than in the decoding of a single fixed backbone.

### Smaller limitations and extensions

Several narrower items remain. The clade-level target repair covers only clades that ASTRAL keeps
monophyletic, so at high discordance a portion of clade events is left uncorrected, and extending
coverage there depends on a better backbone rather than on the repair itself. The placement head
returns a single second parent, which matches the reconstruction but is a one-sided view of a
symmetric pair of parents. The method also does not model post-duplication diploidization, where one
subgenome is gradually fractionated by loss, which would make a polyploid look diploid in a growing
fraction of gene trees and is a natural refinement of both the copy-number features and the
simulation. Finally, the reported comparison uses a copy-bound event selection, and a held-out
calibration of the operating point would make the reported numbers an honest operating choice rather
than a fixed default.
