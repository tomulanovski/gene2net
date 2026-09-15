# Diagnostic: where the reconstruction error comes from

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. All distances are the normalized mu-distance.
Numbers are computed on the held-out validation split, n = 300 networks pooled across the six
training configurations, and the oracle events are read from the ground-truth simulation metadata
in clade-level form rather than decomposed from the MUL-tree.

A reconstruction can lose mu-distance in three places. The build and fold can leave a residual even
on a perfect input, the ASTRAL backbone can misplace lineages, and the predicted events can be
wrong. We separate the three with an oracle. Feeding the true events, read from the simulation
metadata, onto the true species tree leaves only the build residual, which we call the floor.
Feeding the same events onto the ASTRAL backbone adds only the backbone error and gives the ceiling,
the best any event predictor could reach on that backbone. Whatever the model loses beyond the
ceiling is its own event error.

@tab:oracle gives the three points on the validation split, 300 networks pooled across the six
training configurations.

Table: {#tab:oracle} Mu-distance of three reconstructions on the pooled validation split. The floor places the true events on the true backbone, the ceiling places them on the ASTRAL backbone, and the model places its predicted events on the ASTRAL backbone. Lower is better.
| Setting | mu-distance |
| --- | --- |
| Floor, true events on the true backbone | 0.0474 |
| Ceiling, true events on the ASTRAL backbone | 0.1765 |
| Model, predicted events on the ASTRAL backbone | 0.2017 |

The model's mu-distance of 0.2017 splits into a build residual of 0.0474, a backbone contribution of
0.1291 and an event-prediction error of 0.0252. The backbone contributes about five times the error
of event prediction, so the model already extracts nearly everything the ASTRAL backbone allows, and
better detection or partner prediction could recover at most 0.025. The hyperparameter search
agrees, since it raised the detection and partner proxies it was tuned on without moving the
reconstruction measures much. The lever that remains is the species tree, for example a backbone
rebuilt after phasing.

The build residual is small and concentrated. The floor has a median of zero, so more than half of
the networks reconstruct exactly, and its mean of 0.047 comes from a minority with nested or
overlapping reticulations, where one event lands inside a region another has already reshaped. On
those networks the reticulate lineages still match, with a leaf distance of 0.076 against a sister
distance of 0.116, so the residual sits in the sister structure. The floor does not depend on
discordance, because the oracle on the true species tree never uses the gene trees.

These numbers describe error on networks drawn from the training simulator. The benchmark networks
are published topologies rather than simulated ones, and PlaceNet's mu-distance there is higher,
around 0.39 to 0.61, so this decomposition does not account for that gap.
