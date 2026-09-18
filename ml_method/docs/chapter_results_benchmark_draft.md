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

@fig:discordance shows all three measures across the twelve discordance configurations, and
@tab:descendants, @tab:sisters and @tab:mu in the appendix give the exact values.

Figure: {#fig:discordance} figures/discordance_degradation.png | Reconstruction accuracy across the twelve discordance configurations. Rows are the reticulation descendants distance, the reticulation sister distance and the mu-distance. Columns are the duplication and loss rate, none, low, medium and high from left to right, and the x-axis within each panel is the ILS level. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants that use ploidy information. Lower is better.

At the highest duplication and loss rate the ploidy-free decode is ahead of Polyphest at every ILS
level, 0.284 against 0.319 at low ILS, 0.247 against 0.405 at medium ILS and 0.230 against 0.368 at
high ILS. The last two are its widest margins anywhere in the benchmark. Those are also the
configurations where the two decode modes part company. The inferred copy number is corrupted there, which drives the
ploidy-informed mode up to between 0.776 and 0.817 while the ploidy-free mode stays between 0.230
and 0.284, and that is the reason the method offers both modes.

Polyphest leads in every configuration on the other two measures, the parental context in the middle
row and the mu-distance in the bottom row. Both gaps have one cause. The partner head is the weaker
of PlaceNet's two heads, and the mu-distance scores every node of the network, so both charge
PlaceNet for the ASTRAL backbone it never rebuilds, which the diagnostic section measures directly.

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
