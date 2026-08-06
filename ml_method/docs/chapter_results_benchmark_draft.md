# Results: detection and comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Edit distance is the multi-labeled-tree edit
distance computed with canonical child ordering, so it does not depend on the arbitrary Newick
child order. Lower is better for every distance. Items marked TODO are to confirm before
submission.

## Detection

The whole genome duplication detector reaches an F1 of about 0.80 on the held-out validation
split. [TODO confirm precision and recall from the training log.] Detection is deliberately
recall-oriented, because the downstream event selection prunes surplus events under the
inferred copy bound, so a missed event cannot be recovered later while an extra event can be
dropped.

## Reconstruction benchmark

We compare the method to three baselines on the 21 benchmark networks across all six
configurations. Polyphest and iterative GRAMPA both infer ploidy from the data, as our method
does, so all four methods operate from the same inputs. The table reports mean edit distance
per configuration.

| Configuration | GNN | Polyphest | GRAMPA | GRAMPA-iter |
| --- | --- | --- | --- | --- |
| ils_low | 0.594 | 0.387 | 0.802 | 0.789 |
| ils_medium | 0.591 | 0.400 | 0.804 | 0.830 |
| ils_high | 0.589 | 0.538 | 0.825 | 0.906 |
| dup_loss_low | 0.583 | 0.417 | 0.811 | 0.821 |
| dup_loss_medium | 0.593 | 0.442 | 0.823 | 0.838 |
| dup_loss_high | 0.621 | 0.778 | 0.826 | 0.904 |

Three findings stand out.

First, the method beats both GRAMPA and iterative GRAMPA in every configuration, by a wide and
consistent margin of roughly 0.2 to 0.3 in edit distance. Among the prior-free methods that do
not use a supplied ploidy, GRAMPA and its iterative variant are the natural comparison, and the
method is clearly better than both across all conditions.

Second, Polyphest is better than the method at low and moderate difficulty, but the gap closes
as conditions harden and reverses at the hardest one. At high duplication and loss the method
reaches 0.621 against Polyphest's 0.778, so it is the better method there. The reticulation-leaf
distance tells the same story at that configuration, where the method reaches 0.353 against
Polyphest's 0.412. Polyphest infers ploidy by a heuristic that degrades as duplication and loss
erode the copy-number signal, and the method, which detects and places events directly, does not
share that failure mode.

Third, the comparison to Polyphest understates the method's robustness. Polyphest is scored on
15 to 18 of the 19 networks per configuration, because it does not complete the rest, whereas
the method completes all 19. So Polyphest's means are taken over the subset it can finish, while
the method always produces an answer.

On the reticulation distances, Polyphest leads on the networks it completes, the method beats
iterative GRAMPA in every configuration, and the method beats GRAMPA on the reticulation-leaf
distance in most configurations. So the ordering by reticulation recovery is Polyphest first on
its completed subset, then the method, then the GRAMPA variants.

## Runtime

The method reconstructs a network in a few seconds. On the benchmark its per-network compute,
excluding the ASTRAL step, is about 8 seconds, of which the graph neural network forward pass is
a fraction of a second and the rest is feature extraction from the gene trees. ASTRAL adds a few
seconds. Polyphest, by contrast, ranges from seconds on easy networks to days on hard ones, and
does not complete the networks it is missing from the table above. [TODO insert a concrete
Polyphest wall time or range from the run logs.] So the method offers a bounded and predictable
runtime of seconds per network, against a search-based competitor whose runtime is variable and
can be prohibitive.

## Summary

The method is a fast, prior-free reconstruction that beats the GRAMPA family across all
conditions, is competitive with Polyphest and better than it at high duplication and loss, and
runs in seconds on every network where Polyphest can take days and sometimes does not finish.
Polyphest remains the more accurate method at low and moderate difficulty on the networks it
completes, which the next sections trace to the backbone.
