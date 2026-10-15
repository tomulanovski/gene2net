# Discussion

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. This file opens the Discussion and runs to the
limitations. The future directions follow in chapter_future_work_draft.md as a subsection of it.

PlaceNet identifies the polyploid lineages more accurately than Polyphest in all but the mildest
conditions, and more accurately than GRAMPAIter on each of the three measures in almost every
configuration, at a cost of approximately eight seconds per network beyond the ASTRAL step. Its two
decode modes allow it to use a trusted ploidy list or to work without one, and at retention 0.50,
where copy counting has begun to fail, the ploidy-free mode is the most accurate of the five methods
on both event measures. Polyphest nonetheless remains more accurate on the reticulation sister
measure outside the two more severe fractionation levels, and on the μ-distance throughout. Those two
gaps have a single cause.

That cause is the backbone, and it is the principal limitation of the method. PlaceNet grafts events
onto a fixed ASTRAL backbone and never rearranges it, so every error in that backbone reaches the
final network. The diagnostic analysis puts the backbone at 0.1291 of the 0.2017 μ-distance reached
by the model, against 0.0252 for event prediction, so better detection, a better threshold or a
better event selection can move the reconstruction by at most 0.025. The reticulation sister measure
scores the parental context of each event and the μ-distance scores every node of the network, so
both charge PlaceNet for a backbone it never rebuilds, whereas the reticulation descendants measure,
which asks only which lineages are involved, is the one measure on which it leads.

The placement of the second parental lineage traces to the same cause. The model reaches an
allopolyploid partner accuracy of 0.828 on the validation split, but raising it has only a modest
effect on the reconstruction, because the correct partner is well defined only where ASTRAL places
the polyploid coherently. Where ASTRAL scatters the copies of a polyploid across the tree, no change
to the training target or to the features can name its second parental lineage.

Both limitations share one origin, that the copies of a polyploid are not assigned to subgenomes.
Unassigned copies pull ASTRAL in two directions and leave two co-clustering peaks without
establishing which is the home position of the lineage. Phasing the copies would address both, but
phasing acts before the backbone is built, so it constitutes a different method rather than a
modification of this one. A placement head predicting both parental lineages jointly did not improve
the benchmark, which locates the remedy in the assignment of copies rather than in the head.

Two smaller limitations remain. Cluster support carries no signal where the duplicated clade is a
single taxon, because nothing can then be excluded, and it fades once the clade holds about five
taxa, because the ten-leaf neighborhood can no longer reach past the clade. The threshold of the
ploidy-free decode is a fixed default rather than a calibrated one, because the corrupted regime it
targets does not match the clean distribution on which the model was trained.
