# Results: comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The measures are those defined earlier in the
thesis. The reticulation descendants and reticulation sister measures come first, because they ask
directly whether a method recovers the polyploidization events, and the mu-distance follows as a
summary of overall structural agreement. Every measure is computed on the phylogenetic network, so
for scoring, the multi-labeled-tree output of every method, including the present one, is folded to
a network by the shared Holm-algorithm folding of the evaluation. This folding is part of the
scoring, applied identically to every method, and not part of any method. Lower is better for every
distance. PlaceNet appears in its two decode modes, ploidy-informed and ploidy-free. Every method
infers its own ploidy from the gene trees, so Polyphest here is the inferred-ploidy variant and the
comparison is prior-free on both sides. Scores are the mean over five replicates.

The competitors do not complete every network, and the networks they skip are the harder ones.
Every comparison below therefore restricts all methods to the networks that all of them completed.
That subset is limited by Polyphest and by iterative GRAMPA with the ploidy prior, since PlaceNet
completes every network, and the networks it leaves out are discussed at the end of the section.

## Where the method wins

@fig:discordance plots the three measures across the twelve discordance configurations, and
@tab:descendants, @tab:sisters and @tab:mu in the appendix give the exact values.

The comparison spans the fifteen configurations (soon 17). On the reticulation descendants measure
the ploidy-free decode is more accurate than Polyphest in eleven of the fifteen configurations. The
four losses occur under the simplest conditions across two ILS levels: low ILS with no, low, or
medium duplication and loss rates, and medium ILS with no duplication or loss. The ploidy-informed
decode wins nine of the fifteen, so both modes are ahead of Polyphest on most configurations.

Against iterative GRAMPA, the other method that needs no supplied ploidy, the ploidy-free decode is
more accurate on every measure in at least thirteen of the fifteen configurations, and it is ahead
of iterative GRAMPA with the ploidy prior in most of them. On the reticulation sister measure it is
also more accurate than Polyphest at the two more severe fractionation levels. PlaceNet gains on
Polyphest as the conditions get harder, which is where published polyploid datasets sit.

Figure: {#fig:discordance} figures/discordance_degradation.png | Reconstruction accuracy across the twelve discordance configurations. Rows are the reticulation descendants distance, the reticulation sister distance and the mu-distance, and columns are the four configuration families. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants that use ploidy information. Lower is better.

## Accuracy on each measure

The top row of @fig:discordance shows the reticulation descendants measure. The better of PlaceNet's
two modes is ahead of both GRAMPA variants in every configuration, usually by more than a factor of
two against iterative GRAMPA. Against Polyphest the lines cross as the conditions get harder.
Polyphest leads only in the mildest configurations, and from there the ploidy-free decode is ahead,
at high duplication and loss 0.284 against 0.319, 0.247 against 0.405 and 0.230 against 0.368 as ILS
rises. At those same high rates the corrupted copy number drives the ploidy-informed mode up to
between 0.776 and 0.817, while the ploidy-free mode stays between 0.230 and 0.284, which is why the
method offers both modes.

The middle row shows the reticulation sister measure. The ploidy-free decode is more accurate than
iterative GRAMPA in eleven of the twelve configurations, and the same crossover between the decode
modes appears at every high duplication and loss rate. Polyphest leads on this measure throughout,
because the partner head is the weaker of PlaceNet's two heads and inherits the error of the ASTRAL
backbone, as the diagnostic section shows.

The bottom row shows the mu-distance. The ploidy-free decode is more accurate than iterative GRAMPA
in every configuration, and than iterative GRAMPA with the ploidy prior in all but the three
duplication and loss rates at low ILS. Polyphest is ahead on this measure throughout, because the
mu-distance scores every node of the network, including the backbone that PlaceNet takes from ASTRAL
and never rebuilds.

## Completion and runtime

Restricting to the common subset removes four to five networks per configuration, all of them
reconstructed by PlaceNet and missed by at least one competitor. PlaceNet completes all 21 in every
discordance configuration, against 17 to 20 for Polyphest. These networks are harder, with a
PlaceNet mu-distance of 0.56 to 0.73 against 0.39 to 0.61 on the common subset, so what PlaceNet
offers on them is an answer where the competitors return none.

PlaceNet reconstructs a network in about eight seconds beyond the ASTRAL step, almost all of it
feature extraction from the gene trees, and ASTRAL adds a few seconds. Polyphest's runtime spans
four orders of magnitude. On one configuration's 21 networks it ran from about four minutes to more
than three days, and fourteen of the 21 never finished, five at the five-day time limit and nine out
of memory.
