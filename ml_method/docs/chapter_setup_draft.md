# Experimental setup

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Items marked TODO are data-provenance
specifics to confirm from the simulation configuration before submission.

## Data

The method is trained and evaluated on simulated data, so that the true network is known and
the reconstruction can be scored exactly. Species networks are simulated across two sources of
gene-tree discordance. The first is incomplete lineage sorting, varied at low, medium, and high
levels. The second is gene duplication and loss, varied at low, medium, and high rates at a
fixed effective population size. This gives six configurations that span the regimes under
which a detect-then-place method is expected to succeed and to struggle. TODO confirm the
number of taxa per network, the number of gene trees per network, and the replicate structure.

Training uses a large set of synthetic networks simulated under the same six configurations,
[TODO confirm count, on the order of two thousand networks per configuration]. For each network
the gene trees are simulated, the ASTRAL backbone is inferred, and the true events are mapped
to backbone edges to produce the training labels. The species tree is re-rooted before feature
extraction and labeling, because ASTRAL output is unrooted and an arbitrary root inflates the
edit distance and misaligns the labels. A fixed random split holds out one fifth of the
networks for validation, and the same split is used across all experiments so that validation
numbers are comparable.

Evaluation uses a separate benchmark of 21 networks, disjoint from the training networks, with
gene trees simulated under each of the six configurations. This benchmark is the primary test
of generalization, because it measures reconstruction of networks the method never saw during
training.

## Metrics

Reconstructions are scored against the true network by six distance measures, all oriented so
that lower is better. The multi-labeled-tree edit distance is the primary measure of overall
reconstruction quality. The Robinson-Foulds distance measures the topological error of the
backbone. The reticulation count difference measures whether the right number of reticulations
was inferred. The reticulation-leaf and reticulation-sister Jaccard distances measure whether
the reticulations involve the right lineages, and they are the measures most directly
controlled by the placement head. The ploidy difference measures whether each species is
assigned the correct copy number.

## Baselines

The method is compared against two published approaches with different information demands. The
first is Polyphest, which is supplied the true ploidy of each species. The second is an
iterative application of GRAMPA, which like the present method infers ploidy from the data
rather than receiving it. The comparison to GRAMPA is therefore the fair prior-free comparison,
because both methods operate from the same inputs, while the comparison to Polyphest is
informative but favors Polyphest by an input it is given and the present method is not.
