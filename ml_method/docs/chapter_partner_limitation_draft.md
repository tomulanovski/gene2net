# The partner limitation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Numbers are measured on the validation split.
Partner accuracy is the fraction of allopolyploid events for which the predicted parent edge
matches the true parent. Autopolyploid events, where the parent is the event branch itself, are
excluded from this figure because they are easy and inflate the combined number.

## The question

The previous section showed that the backbone, not the events, is the dominant source of edit
distance. One part of the event prediction is nonetheless weak on its own terms and worth
locating, namely the choice of the allopolyploid parent. The one-parter method reaches an
allopolyploid partner accuracy of only 0.45. This section asks whether that number can be raised
by better features or a better decoding rule, and if not, why.

## Better features do not help

We added two new hand-crafted partner features to the shipped model and retrained the one-parter
method from scratch on all six configurations. The first, a duplication-conditioned co-clustering
feature, restricts the co-clustering signal to the gene trees in which the target is actually
duplicated, on the reasoning that single-copy trees carry no partner signal and only dilute it.
The second, an effective-partner-count feature, summarizes how concentrated a species co-clusters,
on the reasoning that an allopolyploid with one clear partner looks different from an autopolyploid
spread over a clade.

| Variant | Allopolyploid partner accuracy | Detection F1 |
| --- | --- | --- |
| Baseline | 0.453 | 0.775 |
| Duplication-conditioned co-clustering | 0.453 | 0.784 |
| Effective-partner-count | 0.466 | 0.788 |
| Both | 0.460 | 0.777 |

The duplication-conditioned feature is a flat null. The effective-partner-count feature adds 0.013,
which is within the sampling noise of a single run given the number of allopolyploid events. These
join the earlier negative results for a cluster-support feature and a branch-length distance
feature, neither of which moved partner accuracy. Across five distinct hand-crafted features the
number does not leave the neighborhood of 0.45. Better features are not the lever.

## Where the error comes from

The reason features saturate is visible when the partner signal is examined directly. Ranking the
candidate parents of an allopolyploid by co-clustering with the target in the gene trees places the
true partner first only 0.42 of the time, but within the top two candidates 0.86 of the time at low
incomplete lineage sorting. The information that identifies the partner is therefore largely
present. The failure is that a second candidate outranks the true partner.

We identified that second candidate. In the 0.58 of events where the true partner is not first, it
is almost always second, and the candidate that beats it is a different species in 0.89 of cases,
with a mean co-clustering to the true partner of only 0.069. That is, the winning candidate is not
a relative of the true partner that a finer resolution would fix. It is a species from a different
part of the tree. It is a diploid in 0.64 of cases. This is the home parent, namely the other
lineage that contributed a subgenome to the allopolyploid. Both subgenomes leave a co-clustering
trace, one copy of the target sisters the home parent and the other sisters the away partner, so
the two strongest candidates are the two parents. The method must return one of them, and the home
parent, whose copy sits cleanly in its own clade, frequently co-clusters more strongly than the
away partner that the label asks for.

## Why the tie cannot be broken cheaply

Breaking the tie between the two parents requires knowing which is the home and which is the away
partner. We tested the three cheap ways to supply that knowledge and all three fail.

The first is to down-weight the home. This requires identifying the home, and the home is the
lineage the target sits beside on the backbone. Identifying it therefore depends on the backbone
placement of the polyploid, which is precisely what the previous section showed ASTRAL gets wrong.
A home mask built on an unreliable placement masks the wrong candidate, and in practice it did not
help.

The second is to separate the two parents by when they co-cluster, on the reasoning that the home
copy is present even in single-copy gene trees while the away partner appears only when the target
is duplicated. This is a non-starter in the data. At low incomplete lineage sorting the target is
duplicated in every gene tree, so there is no single-copy regime to contrast against.

The third is to return both parents rather than one. This is the correct construction, since the
two candidates are the two parents, and it improves edit distance by about 0.19 under an oracle.
But the learned two-parent head reached an allopolyploid set accuracy of only about 0.29 and did
not improve the benchmark, and the co-clustering heuristic that returns the top two candidates also
did not survive on the benchmark networks. The signal that names both parents is clear on the
validation data but does not transfer, and the co-clustering top-two itself degrades from 0.86 at
low incomplete lineage sorting to 0.73 under high duplication and loss.

## Consequence: a shared root cause

The partner limitation and the backbone limitation are the same problem seen from two sides. Both
arise because the two copies of a polyploid are not assigned to their subgenomes. The unassigned
copies confuse ASTRAL, which places the polyploid wrongly and inflates the backbone edit distance.
The same unassigned copies leave two co-clustering peaks, one per parent, and the method cannot say
which peak is the away partner without knowing which copy is which. Assigning the copies to
subgenomes, namely phasing, resolves both at once. One assigned copy fixes the backbone placement
of the polyploid, and the other names the away partner. This is why the improvement path for the
method is not a further feature or a further decoding rule, both of which we have shown to be
exhausted, but the phasing direction developed in the next section.
