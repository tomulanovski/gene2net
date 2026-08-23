# From predictions to a network: the decode

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. This section belongs in the Method chapter,
after the architecture and before the experimental setup. It is the central methodological
novelty of the decode, so it is written to stand on its own.

## The decode problem

The model produces two quantities on every edge of the ASTRAL species tree. The first is a
detection probability, the model's confidence that a whole genome duplication occurred on that
edge. The second is a partner distribution, which for a confirmed event names the edge the
duplicated lineage merged with. These are per-edge predictions. A reconstructed network is a
global object. The decode is the step that turns the per-edge predictions into one multi-labeled
tree, and its only real decision is which edges become events. Once the events are chosen the
partners follow from the partner head and the network is built by grafting each event onto the
backbone.

Choosing the events is not a matter of thresholding the detection probability in isolation,
because the number of events is constrained. A species that carries two gene copies is the
product of one duplication somewhere on its ancestry, and a species that carries one copy is the
product of none. This copy number is exactly the information that a ploidy method such as
Polyphest consumes. Whether to use it, and how much to trust it, is the axis along which the
decode is defined.

## Two modes on a ploidy-dependence axis

The decode is published as two modes that sit at the two ends of this axis. The user selects
between them by a single practical question. Does a trusted list of species copy numbers exist
for the data at hand.

The first mode is ploidy-informed. It fills each species up to its copy number and stops there.
The copy number is taken as an upper bound on the events that touch a species. The method can
supply this list itself, by a self-contained inference described below, or the user can supply a
list obtained from chromosome counts or from an external tool. This is the natural mode when
copy number is reliable, and it is the fair comparison to Polyphest, because feeding both methods
the same copy-number list isolates the difference to placement.

The second mode is ploidy-free. It ignores copy number entirely. An edge becomes an event if and
only if its detection probability is at least a threshold, with no copy-number bound anywhere in the
decode. This makes the mode genuinely free of ploidy, since it computes no copy-number estimate
at all. It is the natural mode when no trusted copy number exists, and it is the fair comparison
to iterative GRAMPA, which likewise searches for reticulations without a ploidy input. It is also
the mode that survives fractionation, the regime in which the copy number is not merely unknown
but actively wrong, and it is analysed in the fractionation section.

The two modes need no separate training. They are two ways of reading the same detection and
partner heads, so a single trained model serves both, and the choice is made at inference time.

## The copy number is a ceiling, not a target

In ploidy-informed mode the copy number bounds the events but does not dictate them. The decode
adds events in order of detection confidence and refuses any event that would push a species past
its copy number. The result is that each species ends with at most its copy number of copies, and
it can end with fewer when the model is not confident enough to fill it. This ceiling behavior is
a deliberate choice, and the alternative is worth stating because it is the obvious one.

The alternative is to force each species to exactly its copy number, adding events until the
count is met whether or not the model supports them. We tested this forcing. It recovers the copy
number by construction, so it scores well on the copy-number measures, but it places the forced
copies on edges where the detection head has no support, and those arbitrary placements raise the
mu-distance and the reticulation-sister error. Forcing the exact count therefore buys the
copy-number measures at the cost of the topology, and the cost is the larger one. The ceiling
declines to invent placements it cannot justify, which is why it is the published behavior.

There is a principled version of exact filling that we did not implement and note as future work.
Rather than force copies onto arbitrary edges, one could search for the set of predicted events
that meets every species copy number exactly and has the highest total detection confidence among
all such sets. This is a constrained optimization over the model's own predictions, so it invents
nothing, and the confidence-ordered ceiling is a greedy approximation of it. The two agree
whenever the confident events do not compete for a shared copy-number budget, which is the common
case when events are sparse, so the expected gain is small and confined to networks with
overlapping events, and it comes at the cost of a search. We therefore report the greedy ceiling
and leave the exact objective to future work.

## The method's own copy-number inference

When the user does not supply copy numbers, the ploidy-informed mode infers them from the gene
trees, so the method remains self-contained. For each species we count its copies in every gene
tree, including the trees in which it is absent, which contributes a count of zero. This gives a
distribution of per-tree copy numbers for the species. Occasional duplication and loss make this
distribution noisy, so we do not take its maximum or its mean. We smooth it with a triangular
kernel and take the peak, which is robust to a few gene trees that duplicate or lose a copy. A
species seen with a single copy in almost every tree is inferred diploid, and a species seen with
two copies in the bulk of trees is inferred to carry two, regardless of a handful of outlier
trees. This is the same representative-copy-number estimate that underlies an inferred-ploidy
multiset, reimplemented inside the method so that ploidy-informed mode needs no external input.

## The threshold in ploidy-free mode

Ploidy-free mode has one parameter, the detection threshold at or above which an edge becomes an
event. We fix it at 0.5, the classifier's natural decision boundary, and expose it as a configurable
option. A fixed default is the honest choice here, because the regime in which the ploidy-free mode
is used, namely corrupted copy number, is not the regime the model was trained and selected on, so a
threshold tuned on the clean training distribution would not transfer to it in a principled way.
Calibrating the threshold per regime, ideally on data that matches the corruption the mode targets,
is left as future work.

## Building the multi-labeled tree

Once the events are selected and their partners named, they are grafted onto a copy of the ASTRAL
backbone to produce the multi-labeled tree. The grafting mirrors the way the training networks are
generated, so the output has the same form as the ground truth. Events are applied smallest target
clade first, so that an event nested inside another is placed before the outer event that may
duplicate it.

An autopolyploidy, where the partner is the clade itself, is grafted by attaching an identical copy
of the target clade as a sibling. The target is detached, a new internal node is inserted in its
place under the same parent, and the target and a deep copy of it become the two children of that
node. This produces two identical sibling subtrees, which is the signature of an autopolyploidy and
folds later to ploidy without a reticulation.

An allopolyploidy, where the partner is a different clade, is grafted by leaving the target at its
backbone position and attaching a copy of it on the partner's edge. The partner edge is subdivided
by a new internal node that holds the partner and the copy of the target. The original target
therefore stays where ASTRAL placed it and only the copy moves to the second parent, which is why
the model predicts a single second parent rather than both. This is the structural root of the
placement limitation discussed later, since the home position of the polyploid is fixed by the
backbone and never chosen by the model.

An event whose target or partner clade can no longer be located, which happens when an earlier graft
has added foreign leaves inside that clade, is dropped and counted. The dropped count is small and
concentrated on networks with nested or overlapping events, and it is the source of the small build
residual reported in the diagnostic.

## Summary

The decode turns per-edge detection and partner predictions into a multi-labeled tree by choosing
which edges become events, and it is defined along how much it trusts the species copy number. The
published method offers the two ends of that axis as two modes. Ploidy-informed mode fills each
species up to a copy-number ceiling, inferred by the method or supplied by the user, and never
forces an unsupported copy. Ploidy-free mode discards copy number and keeps the edges the
detection head is confident about, at a fixed default threshold. The user selects the
mode by whether a trusted copy-number list exists, which also fixes the fair baseline for each
mode, Polyphest for the ploidy-informed mode and iterative GRAMPA for the ploidy-free mode.
