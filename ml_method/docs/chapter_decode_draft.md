# The Decode

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes.

The two prediction heads return per-branch quantities, a detection probability for every branch of
the backbone and, for every branch, a distribution over the branches that could have contributed the
second subgenome. A reconstruction is a single global object. The decode is the step turning those
per-branch quantities into one MUL-tree, and its only substantive decision is which branches become
events, since the partner of an accepted branch is read off as the highest scoring candidate.

PlaceNet decodes in either of two modes, differing only in how far they trust the gene copy number.
Both read the same trained model, so no separate training is required and the mode is selected at
inference.

## Ploidy-Informed Decode

The ploidy-informed mode takes a ploidy list, one entry per taxon giving the number of subgenomes
that taxon carries, so an entry of two marks a tetraploid. Each entry is read as an upper bound on
the number of events that may touch the taxon. The decode sorts the branches by detection
probability and accepts each in turn, refusing any event that would push a taxon past its entry, and
no probability threshold is applied. A taxon therefore ends with at most as many copies as its
entry, and with fewer only when every remaining branch that would raise it would push some other
taxon past its own entry. This is the natural mode where the list is trustworthy and the fair
comparison against Polyphest, since supplying both methods with the same list isolates the
difference to the placement of the events.

Where no external list is available it can be inferred from the gene trees, which is the route taken
for the simulated benchmark, since the simulated datasets carry no external list. For each taxon we
counted its copies in every gene tree, including the gene trees from which it is absent, which
contribute a count of zero, yielding a distribution of per-tree copy numbers. Gene duplication and
loss make that distribution noisy, so neither its maximum nor its mean reflects the number of
subgenomes. We therefore smoothed it with a triangular kernel and took its peak, by the same
procedure applied to the empirical datasets in the preceding chapter. A taxon observed with a single
copy in almost every gene tree is inferred diploid and one observed with two copies in the bulk of
gene trees is inferred to carry two, irrespective of a few outlying gene trees.

The greedy fill is not the only way to respect the list. An exact decode would select, among all
sets of events meeting every entry exactly, the set of highest total detection probability, which is
a constrained optimization over the predictions rather than a sorted sweep. We report the greedy
fill and return to the exact objective in the future work.

## Ploidy-Free Decode

The ploidy-free mode requires no list and estimates nothing in its place. A branch becomes an event
if and only if its detection probability reaches a threshold, with no bound imposed anywhere in the
decode, which makes the mode independent of copy number. It is the natural mode where no trustworthy
list exists and the fair comparison against GRAMPAIter, which likewise searches for events without a
ploidy input. It is also the mode that survives fractionation, the regime in which copy number is
not merely unknown but systematically wrong.

The threshold is the single parameter of this mode. We fixed it at 0.5, the natural decision
boundary of the classifier, and exposed it as a configurable option. A fixed default was used
because the regime the mode targets, corrupted copy number, is not the regime on which the model was
trained and selected, so a threshold tuned on the clean training distribution would not transfer to
it. Calibrating the threshold on data matching the intended corruption is left as future work.

## Building the MUL-tree

Once the events are selected and their partners named, they are grafted onto a copy of the backbone.
Events are applied smallest source clade first, so that an event nested inside another is placed
before the outer event that may duplicate it.

Every event names a source clade and a partner branch, and the graft is one operation in both cases.
The partner branch is subdivided by a new node, and that node carries whatever the branch already
held together with a copy of the source clade. Where the partner branch is the source branch itself,
the copy lands beside the clade and gives two identical sibling subtrees, which is an
autopolyploidization and folds later to a higher ploidy level without a reticulation. Where the
partner branch is elsewhere, the copy lands beside a different lineage and folds to a reticulation,
which is an allopolyploidization. Either way the source clade remains where ASTRAL placed it and
only the copy moves, which is why the model predicts one partner rather than two.

An event is recorded as a clade, a set of taxa, because that is what a branch of the backbone is.
Each graft changes the tree, so every later graft must locate its source clade and its partner clade
in the tree as it then stands. Usually the taxa still match exactly. Where an earlier graft has added
copies inside the clade they no longer do, and the builder falls back to the smallest clade still
containing all of them, which is that clade as grown by the earlier graft. An event is skipped only
when the search reaches the root. Both cases are confined to networks with nested or overlapping
events and together constitute the build residual measured in the diagnostic analysis.
