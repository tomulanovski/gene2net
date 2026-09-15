# Limitations

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes.

The main limitation is the backbone. PlaceNet places events on a fixed ASTRAL backbone and never
rearranges it, so every backbone error reaches the final network. The diagnostic puts this at 0.1291
of the model's 0.2017 mu-distance, against 0.0252 for event prediction, so better detection,
thresholds or event selection can move the reconstruction by at most 0.025.

Placement of the second parent traces to the same cause. The model reaches an allopolyploid partner
accuracy of 0.828 on the validation split, but raising it has only a modest effect on the
reconstruction, because the correct partner is well defined only when ASTRAL places the polyploid
coherently. When ASTRAL scatters a polyploid, no change to the training target or the features can
name its second parent.

Both problems share one cause, which is that the copies of a polyploid are not assigned to their
subgenomes. Unassigned copies pull ASTRAL in two directions and leave two co-clustering peaks
without saying which is the home position. Phasing the copies would address both, but it acts before
the backbone is built, so it is a different method rather than a change to this one. A placement
head that predicted both parents jointly did not improve the benchmark, which supports placing the
fix in the copy assignment.

Smaller limitations remain. The clade-level target repair covers only clades that ASTRAL keeps
monophyletic. Cluster support carries no signal when the duplicated clade is a single species,
because nothing can be excluded there, and it fades once the clade holds about five species, because
the ten-leaf group can no longer reach past the clade. The threshold of the ploidy-free decode is a
fixed default, because the corrupted regime it targets does not match the clean distribution the
model was trained on.
