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

The ploidy level of a species is not always known. It can come from chromosome counts or from an
earlier study, and for many datasets there is nothing reliable to go on. The decode therefore works
either way. One mode uses prior ploidy information when the user has it, and the other uses none.

## Two decode modes

The first mode is ploidy-informed. It takes a ploidy list, one entry per species giving how many
subgenomes that species carries, so an entry of two marks a tetraploid with two subgenomes. Where
the list comes from is up to the user. The decode reads each entry as an upper bound on the events
that may touch that species, and fills the species up to it. This is the natural mode when the list
is trustworthy, and it is the fair comparison to Polyphest, because giving both methods the same
list isolates the difference to placement.

The second mode is ploidy-free. It needs no list and estimates nothing in its place. An edge becomes
an event if and only if its detection probability is at least a threshold, with no bound anywhere in
the decode. This makes the mode genuinely free of ploidy. It is the natural mode when no trustworthy
list exists, and it is the fair comparison
to iterative GRAMPA, which likewise searches for reticulations without a ploidy input. It is also
the mode that survives fractionation, the regime in which the copy number is not merely unknown
but actively wrong, and it is analyzed in the fractionation section.

The two modes need no separate training. They are two ways of reading the same detection and
partner heads, so a single trained model serves both, and the choice is made at inference time.

## Filling to the ploidy list

In ploidy-informed mode the list bounds the events but does not dictate them. The decode adds
events in order of detection confidence and refuses any event that would push a species past its
entry. The result is that each species ends with at most that many copies. It ends with fewer only
when no remaining edge fits, because every edge that would raise this species would push some other
species past its own entry. An exact version, which selects the set of events that meets every entry
exactly and has the highest total detection confidence among all such sets, is a constrained
optimization over the model's own predictions rather than a greedy fill. We report the greedy fill
and discuss the exact objective in the future-work section.

## Estimating the list when none is supplied

A user without a ploidy list can run the ploidy-free mode, or let the method estimate a list from
the gene trees and run ploidy-informed on that. The benchmark takes the second route, since the
simulated datasets carry no external list. For each species we count its copies in every gene
tree, including the trees in which it is absent, which contributes a count of zero. This gives a
distribution of per-tree copy numbers for the species. Occasional duplication and loss make this
distribution noisy, so we do not take its maximum or its mean. We smooth it with a triangular
kernel and take the peak, which is robust to a few gene trees that duplicate or lose a copy. A
species seen with a single copy in almost every tree is inferred diploid, and a species seen with
two copies in the bulk of trees is inferred to carry two, regardless of a handful of outlier
trees.

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
backbone to produce the multi-labeled tree. Events are applied smallest source
clade first, so that an event nested inside another is placed before the outer event that may
duplicate it.

Every event names a source clade and a partner edge, and the graft is one operation in both cases.
The partner edge is subdivided by a new node, and that node holds whatever the edge already carried
together with a copy of the source clade.

When the partner edge is the source edge itself, the copy lands beside the clade and gives two
identical sibling subtrees. This is an autopolyploidy, and it folds later to ploidy
without a reticulation. When the partner edge is elsewhere, the copy lands beside a different
lineage and folds to a reticulation, which is an allopolyploidy.

Either way the source clade stays where ASTRAL placed it and only the copy moves, which is why the
model predicts one partner rather than two. This is the structural root of the placement limitation
discussed later, since the home position of the polyploid is fixed by the backbone and never chosen
by the model.

Each graft has to find the node whose leaves are exactly the source clade, and the node for the
partner clade. An earlier graft can add copies inside one of those clades, and then no node has
exactly those leaves any more. Such an event is dropped and counted. The dropped count is small and
concentrated on networks with nested or overlapping events, and it is the source of the small build
residual reported in the diagnostic.

## Summary

The decode turns per-edge detection and partner predictions into a multi-labeled tree by choosing
which edges become events, and it is defined along how much it trusts the species copy number. The
published method offers both, as two modes. Ploidy-informed mode fills each
species up to a copy-number ceiling, inferred by the method or supplied by the user, and never
forces an unsupported copy. Ploidy-free mode discards copy number and keeps the edges the
detection head is confident about, at a fixed default threshold. The user selects the
mode by whether a trustworthy ploidy list exists, which also fixes the fair baseline for each
mode, Polyphest for the ploidy-informed mode and iterative GRAMPA for the ploidy-free mode.
