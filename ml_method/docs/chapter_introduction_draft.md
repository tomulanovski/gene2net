# Introduction

NOTE: move candidate. The first paragraph below motivates PlaceNet from the two limitations the
benchmarking chapter exposed. Once the shared introduction to both chapters is written, that
paragraph belongs there and can be lifted verbatim, leaving this chapter to open at "To address both
limitations".

The preceding chapter benchmarked the available methods for polyploid network reconstruction and
exposed two practical limitations. The first is runtime. The existing methods search a space of
candidate networks, so their cost grows with the number of candidate polyploid clades rather than
with the size of the input, and on the larger simulated networks Polyphest, the most accurate of
them, reached the five-day limit or exhausted memory instead of returning a result. The second is
the reliance on gene copy number. Polyphest and AlloppNET require the ploidy level of every taxon
as input, and where it is unknown it must be inferred from the number of gene copies observed per
taxon. Copy number is a reliable proxy on clean data, but fractionation, the gradual loss of one
duplicate after a polyploidization event, erodes it, and published polyploid complexes are
extensively fractionated. A method whose reconstruction does not rest on per-taxon copy number
should therefore be more robust on empirical data.

To address both limitations we developed PlaceNet, a graph neural network that reconstructs a
polyploid network by detecting polyploidization events on a single-copy species tree and then
placing them. PlaceNet infers a diploid backbone with ASTRAL-4, represents that backbone as a graph
whose nodes are the species tree nodes and whose edges are its branches, and annotates the graph
with features summarizing the gene trees. Two prediction heads then read every branch of the
backbone, one deciding whether a polyploidization event occurred on it and the other naming the
branch that contributed the second subgenome. The accepted events are grafted onto the backbone,
yielding a MUL-tree, which is the output of the method. Inference is a single forward pass over one
graph, so its cost follows the size of the input rather than the number of events, and the features
driving detection measure duplication shared across a clade rather than counting copies per taxon.

Because the ploidy levels are sometimes known and sometimes not, PlaceNet decodes its predictions in
either of two modes, one taking a ploidy list as an upper bound on the events touching each taxon
and one using no ploidy information at all. The two modes read the same trained model and are
selected at inference, which also fixes the fair baseline for each, Polyphest for the
ploidy-informed mode and GRAMPAIter for the ploidy-free mode. PlaceNet is useful as it stands but
remains under active development, and the directions identified below are set out in the discussion.
