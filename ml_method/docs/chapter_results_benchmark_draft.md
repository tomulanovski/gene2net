# Results

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Lower is better for every distance.

PlaceNet is reported in both decode modes, ploidy-informed and ploidy-free, against Polyphest and
GRAMPAIter. Every output was folded to a network by the shared HOLM folding of the evaluation before
scoring, so the folding is applied identically to every method and is part of neither. Every value
is the mean over five replicates.

The baselines do not complete every network, and the networks they leave out are the harder ones.
Every comparison below therefore restricts all methods to the networks that all of them completed,
a subset limited by Polyphest and by GRAMPAIter with the ploidy prior, since PlaceNet completed
every network in every discordance configuration. Between 16 and 18 of the 21 networks remain, and
the networks left out are treated at the end of this section.

## Comparison Against the Baselines

@fig:discordance shows the three measures across the twelve discordance configurations, and
@tab:descendants, @tab:sisters and @tab:mu in the appendix give the exact values.

Figure: {#fig:discordance} figures/discordance_degradation.png | Reconstruction accuracy across the twelve discordance configurations. Rows are the reticulation descendants distance, the reticulation sister distance and the μ-distance. Columns are the duplication and loss rate, none, low, medium and high from left to right, and the horizontal axis within each panel is the ILS level. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants using ploidy information. Lower is better.

Alt text: A grid of twelve line plots in three rows and four columns, each plotting a distance against three ILS levels for five methods. In the top row PlaceNet lies below Polyphest in most panels and the gap widens from left to right. In the middle and bottom rows Polyphest lies below every other method in every panel.

On the reticulation descendants measure the ploidy-free decode is more accurate than Polyphest in
eleven of the fifteen configurations. The four exceptions lie under the mildest discordance, namely
low ILS at each of the three lowest duplication and loss rates and medium ILS with no duplication
and loss. The ploidy-informed decode is more accurate than Polyphest in nine of the fifteen, so both
modes are ahead on most configurations. The advantage widens as discordance increases, which is the
regime the published polyploid datasets occupy. At the highest duplication and loss rate the
ploidy-free decode is ahead of Polyphest at every ILS level, 0.284 against 0.319 at low ILS, 0.247
against 0.405 at medium ILS and 0.230 against 0.368 at high ILS, the last two being its widest
margins anywhere in the benchmark.

Those same configurations are where the two decode modes part company. The copy number inferred from
the gene trees is corrupted at the highest duplication and loss rate, which drives the
ploidy-informed mode to between 0.776 and 0.817 while the ploidy-free mode stays between 0.230 and
0.284. This is the reason the method offers both modes.

Against GRAMPAIter, the other method requiring no supplied ploidy, the ploidy-free decode is more
accurate in at least thirteen of the fifteen configurations on every one of the three measures,
fourteen on the reticulation descendants distance, thirteen on the reticulation sister distance and
all fifteen on the μ-distance. Supplying GRAMPAIter with the same inferred ploidy list as a prior
narrows the gap on the reticulation sister distance, where the ploidy-free decode leads in nine of
the fifteen, but not on the other two, where it leads in fourteen and twelve.

Polyphest is more accurate than both decode modes on the reticulation sister distance in thirteen of
the fifteen configurations, the two exceptions being the more severe fractionation levels, and on the
μ-distance in all fifteen. Both gaps have one cause. The reticulation sister distance scores the
parental context of each event and the μ-distance scores every node of the network, so both charge
PlaceNet for the backbone it never rebuilds, which the diagnostic analysis measures directly.

Restricting the comparison to the common subset understates what PlaceNet offers, because it is
never the method that limits the subset. On the networks left out its μ-distance is 0.56 to 0.73
against 0.39 to 0.61 on the common subset, so these are harder networks, but on them it
returns a reconstruction where the baselines return none.

## Runtime

PlaceNet reconstructs a network in approximately eight seconds beyond the ASTRAL step, almost all of
it feature extraction from the gene trees. The model itself is one forward pass and the decode is a
single sweep over the branches, so the cost follows the size of the input rather than the difficulty
of the instance. Polyphest searches instead, and its runtime on the same networks spans several
orders of magnitude, from minutes on the smallest networks to days on the largest. It returned a
result for all five replicates of between 15 and 17 of the 21 networks depending on the
configuration, and for at least one replicate of between 17 and 20, the remainder reaching the
five-day limit or exhausting the available memory. PlaceNet returned a reconstruction for all 21
networks in every discordance configuration.
