# Limitations and future work

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes.

## The backbone is the binding limitation

The central limitation of the method is structural and is established directly by the diagnostic
section. The reconstruction is built by stamping events onto a fixed ASTRAL backbone and never
rearranges that backbone, so every backbone error is inherited by the final network. Under an
oracle that supplies the true events with their true parents, placing them on the ASTRAL
backbone still leaves an edit distance above 0.54, while placing the same events on the true
backbone reaches 0.11. No detector and no placement head can cross that gap, because the method
does not rebuild the backbone. This is why better detection, better thresholds, and better event
selection do not move edit distance, and it is the frame for everything below.

## Placement is bounded by the same backbone

Placement of the second parent looked like an independent weakness, but it reduces to the same
cause. Hand-crafted features do not improve it. A defect in the training target does, and
repairing the target raises partner accuracy and improves ploidy and edit distance on the
benchmark, but the improvement is modest and does not generalize to the finer reticulation
measures. The reason is that the correct parent is well-defined only when ASTRAL places the
polyploid coherently, and undefined when it does not. The residual placement error and the
coverage bound of the repair both trace to the backbone. When ASTRAL scatters a polyploid or a
polyploid clade, there is no coherent home against which to name the second parent, and the
event cannot be repaired by any change to the target or the features. Placement is therefore not
a separate ceiling but a second face of the backbone limitation.

## The root cause and the path through it

Both limitations share one mechanical cause. The two copies of a polyploid are not assigned to
their subgenomes. Unassigned copies pull ASTRAL in two directions at once, which is what
scatters the backbone, and the same unassigned copies leave two co-clustering peaks without
telling the method which peak is the home and which is the second parent. Assigning the copies
to subgenomes, which is phasing, addresses both faces at once. A phased copy defines the
backbone position of the polyploid, and the other phased copy names the second parent. Phasing
is therefore the principled next step, and it is a different method rather than a change to the
present one, because it operates before the backbone is built rather than decorating a backbone
after the fact. A learned two-parent placement head that predicts both parents jointly was
attempted as a lighter alternative and did not improve the benchmark, which further supports
locating the fix in the copy assignment rather than in the decoding of a single fixed backbone.

## Smaller limitations and extensions

Several narrower items remain. The clade-level target repair covers only clades that ASTRAL
keeps monophyletic, so at high discordance a portion of clade events is left uncorrected, and
extending coverage there depends on a better backbone rather than on the repair itself. The
placement head returns a single second parent, which matches the reconstruction but is a
one-sided view of a symmetric pair of parents. The method also does not model post-duplication
diploidization, where one subgenome is gradually fractionated by loss, which would make a
polyploid look diploid in a growing fraction of gene trees and is a natural refinement of both
the copy-number features and the simulation. Finally, the reported comparison uses a single
detection threshold and a copy-bound event selection, and a held-out calibration of these
operating points would make the reported numbers an honest operating choice rather than a fixed
default.
