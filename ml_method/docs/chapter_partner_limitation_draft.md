# Partner prediction: a labeling defect and its repair

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Numbers are measured on the validation split
for accuracy and on the 21 benchmark networks for the reconstruction metrics, unless stated
otherwise. Partner accuracy is the fraction of allopolyploid events for which the predicted
parent edge matches the target. Autopolyploid events, where the parent is the event branch
itself, are excluded from this figure because they are easy and inflate the combined number.

## The question

The previous section showed that the backbone, not the events, is the dominant source of
mu-distance. One part of the event prediction is nonetheless weak on its own terms, namely
the choice of the allopolyploid parent, and it is worth locating. The one-parter method
reaches an allopolyploid partner accuracy of only 0.45. This section asks whether that number
can be raised, first by better features and then by a correction to the training target, and
what the answer implies.

## Better features do not help

We added two hand-crafted partner features to the shipped model and retrained the one-parter
method from scratch on all six configurations. The first restricts the co-clustering signal
to gene trees in which the target is actually duplicated, on the reasoning that single-copy
trees carry no partner signal and only dilute it. The second summarizes how concentrated a
species co-clusters, on the reasoning that an allopolyploid with one clear partner looks
different from an autopolyploid spread over a clade.

| Variant | Allopolyploid partner accuracy | Detection F1 |
| --- | --- | --- |
| Baseline | 0.453 | 0.775 |
| Duplication-conditioned co-clustering | 0.453 | 0.784 |
| Effective-partner-count | 0.466 | 0.788 |
| Both | 0.460 | 0.777 |

The duplication-conditioned feature is a flat null. The effective-partner-count feature adds
0.013, within the sampling noise of a single run. These join earlier negative results for a
cluster-support feature and a branch-length distance feature. Across five hand-crafted
features the number does not leave the neighborhood of 0.45. Better features are not the lever.

## Where the error comes from

The reason features saturate is visible when the signal is examined directly. Ranking the
candidate parents of an allopolyploid by co-clustering with the target in the gene trees
places the true partner first only 0.42 of the time, but within the top two candidates 0.86
of the time at low incomplete lineage sorting. The information that identifies the partner is
therefore largely present. The failure is that a second candidate outranks the true partner.

We identified that second candidate. In the 0.58 of events where the true partner is not
first, it is almost always second, and the candidate that beats it is a species from a
different part of the tree in 0.89 of cases, with a mean co-clustering to the true partner of
only 0.069. It is a diploid in 0.64 of cases. This is the home parent, namely the other
lineage that contributed a subgenome. Both subgenomes leave a co-clustering trace, so the two
strongest candidates are the two parents, and the method must return one of them. The home
parent, whose copy sits cleanly in its own clade, frequently co-clusters more strongly than
the away partner the label asks for.

## A labeling defect

Examining which of the two parents the label names revealed a defect in the training target
itself. The label is taken from the true network, where one parent is drawn as the base-tree
attachment and the other as a reticulation edge, and the reticulation parent becomes the
label. This choice is made without reference to where ASTRAL placed the polyploid. ASTRAL
places the polyploid next to whichever subgenome is stronger, which is not correlated with the
network drawing, so the two disagree about half the time. A direct audit against ground truth
confirms it. In 0.566 of single-species allopolyploid events at low incomplete lineage sorting
and 0.539 at high duplication and loss, the labelled partner is the same lineage ASTRAL placed
the target next to. In those events the target is asked to reticulate to the parent it is
already attached to, which is degenerate, and the correct partner is the other parent.

The audit also bounds how often the defect is fixable. ASTRAL places a single-species target
next to a true parent in about 0.94 of events, so the correct partner, the parent that is not
the target's home, is well-defined almost always. Only about 0.06 of targets are placed away
from both parents, which is a backbone error that no partner choice can repair.

## The repair and its effect

The repair is to recompute the target as the parent that is not the polyploid's ASTRAL home.
When the ASTRAL home is the labelled partner, the target is retargeted to the other true
parent. Otherwise it is left unchanged. The true tree is used only to construct the corrected
label. At inference the model sees only the ASTRAL tree, the gene trees, and the features, so
it learns the rule that the partner is the parent it was not placed beside, which is inferable
from its inputs.

Training the one-parter twice under identical settings, differing only in the labels, isolates
the effect. Partner accuracy rises from 0.456 to 0.526, an increase of 0.070 that is many
standard errors beyond noise, while detection F1 is unchanged at 0.811 against 0.819, which
confirms the change is confined to the partner target.

The gain also carries to the reconstruction, measured against the true networks on the 21
benchmark networks across all six configurations. The effect is real but modest. Ploidy
distance improves in all six configurations, from a mean of about 0.455 to about 0.395. Edit
distance improves in five of six, from a mean of about 0.715 to about 0.701, with a small
regression at high incomplete lineage sorting that is within the per-configuration variance.
The reticulation-sister and reticulation-leaf distances are mixed across configurations and
flat on average. So the corrected target recovers polyploidy structure more reliably and edit
distance slightly, but the large reticulation gain seen at the easiest configuration does not
generalize. The validation gain in partner accuracy does not fully transfer to the finer
reticulation metrics on the benchmark, which is consistent with the backbone remaining the
limiting factor.

## Coverage and its bound

The repair as described covers single-species allopolyploid events. A clade-level event, where
a group of species shares one hybridization, is handled by the same rule applied to the
clade's common ancestor, but only when the clade is monophyletic in the ASTRAL tree, so that
its home is a single well-defined sibling. This holds for about 0.60 of clade events at low
incomplete lineage sorting and about 0.42 at high duplication and loss. The remaining clade
events are scattered by ASTRAL and left unchanged, because for a scattered clade there is no
single home and therefore no well-defined away parent to assign. This is not an unrepaired
bug. It is the absence of a correct target, and it arises because the backbone placed the
clade incorrectly in the first place.

## Consequence

Partner prediction is limited by two things, and only one of them is a partner problem.
Hand-crafted features are exhausted, but the target itself was wrong in the majority of
allopolyploid events, and repairing it raises partner accuracy by 0.070 and improves ploidy
and mu-distance on the benchmark. What remains is not a partner problem. The parents are
recoverable whenever ASTRAL places the polyploid coherently, and undefined whenever it does
not, so the residual error and the coverage bound both trace to the same source as the
reconstruction wall, namely the backbone. This is why the corrected target helps but does not close the gap,
and why the remaining ceiling is a better backbone through phasing rather than any further
change to the target or the features.
