## Where the Reconstruction Error Arises

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. All distances are the normalized μ-distance.

The two event measures report how often PlaceNet recovers the right lineages and the right parental
lineages, but they cannot say where the error arises, because they score only the reticulations. The
μ-distance can, because it scores every node of the network, so this analysis works in μ-distance
throughout.

A reconstruction can lose μ-distance in three places. The build and fold can leave a residual even
on a perfect input, the ASTRAL backbone can misplace lineages, and the predicted events can be
wrong. We separated the three with an oracle. Grafting the true events, read from the simulation
metadata in clade-level form rather than decomposed from the MUL-tree, onto the true species tree
leaves only the build residual, which we call the floor. Grafting the same events onto the ASTRAL
backbone adds the backbone error and gives the ceiling, the best any event predictor could reach on
that backbone. Whatever the model loses beyond the ceiling is its own event error. The three points
were computed on a seeded random subsample of 300 samples from the held-out validation split, pooled
across the six training configurations, and are given in @tab:oracle.

Table: {#tab:oracle} μ-distance of three reconstructions on the validation split. The floor grafts the true events onto the true backbone, the ceiling grafts them onto the ASTRAL backbone, and the model grafts its own predicted events onto the ASTRAL backbone. Lower is better.
| Setting | μ-distance |
| --- | --- |
| Floor, true events on the true backbone | 0.0474 |
| Ceiling, true events on the ASTRAL backbone | 0.1765 |
| Model, predicted events on the ASTRAL backbone | 0.2017 |

The μ-distance of 0.2017 reached by the model therefore decomposes into a build residual of 0.0474,
a backbone contribution of 0.1291 and an event-prediction error of 0.0252. The backbone contributes
approximately five times the error of event prediction, so the model already extracts close to
everything the ASTRAL backbone permits, and better detection or placement could recover at most
0.025. The hyperparameter search is consistent with this, having raised the detection and placement
measures it was tuned on without moving the reconstruction measures appreciably. The lever that
remains is the species tree itself.

The build residual is small and concentrated. The floor has a median of zero, so more than half of
the networks reconstruct exactly, and its mean of 0.047 comes from a minority of networks carrying
nested or overlapping reticulations, where one event lands inside a region another has already
reshaped. On those networks the reticulate lineages themselves are still recovered, the reticulation
descendants distance being 0.076 against 0.116 for the reticulation sister distance, so the residual
sits in the parental context rather than in the identity of the polyploid clades. The floor does not
depend on discordance, because the oracle on the true species tree never reads the gene trees.

These values describe error on networks drawn from the training simulator. The benchmark networks
are published topologies rather than simulated ones, and the μ-distance of PlaceNet on them is
higher, between 0.39 and 0.61, which this decomposition does not account for.

A second question is what the detection head sees once fractionation has deleted the copies. We read
the detection probability the model assigned to every lineage in the ground truth and separated the
polyploids whose copies fractionation had deleted, which copy counting reads as diploid, from the
true diploids. On those lost events the head fires at 0.55 on average and clears the 0.5 threshold
of the decode for 59 percent of them at retention 0.50, against 0.09 and 6 percent on the true
diploids, and at retention 0.25 it still averages 0.30. The two figures are informative only as a
pair, since a detector keyed on copy number would score both groups alike. The ploidy-free decode
therefore recovers events under fractionation that the ploidy-informed decode cannot use at all,
that mode filling each taxon only up to a copy number fractionation has already reduced to one.
