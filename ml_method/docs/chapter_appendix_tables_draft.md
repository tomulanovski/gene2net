# Appendix: numerical results

DRAFT for the thesis. Exact values behind the two degradation figures, one table per measure for
the discordance benchmark and one per level for fractionation. scripts/chapter_figures.py reads
these tables for its preview and checks the cluster figures against them.

@tab:descendants, @tab:sisters and @tab:mu give the values plotted in @fig:discordance, and
@tab:frac075, @tab:frac050 and @tab:frac025 give those plotted in @fig:fractionation.

Table: {#tab:descendants} Reticulation descendants distance on the twelve discordance configurations, plotted in the top row of @fig:discordance. Each value is the mean over five replicates on the networks that every method completed, and n is the number of those networks. Lower is better.
| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.112 | 0.174 | 0.022 | 0.418 | 0.223 | 17 |
| ILS medium | 0.134 | 0.205 | 0.148 | 0.452 | 0.228 | 16 |
| ILS high | 0.157 | 0.207 | 0.245 | 0.573 | 0.210 | 16 |
| dup/loss low, Ne 200k | 0.111 | 0.166 | 0.127 | 0.374 | 0.210 | 16 |
| dup/loss medium, Ne 200k | 0.122 | 0.190 | 0.135 | 0.405 | 0.188 | 16 |
| dup/loss high, Ne 200k | 0.782 | 0.284 | 0.319 | 0.527 | 0.306 | 17 |
| dup/loss low, Ne 1M | 0.154 | 0.188 | 0.194 | 0.467 | 0.222 | 16 |
| dup/loss medium, Ne 1M | 0.136 | 0.195 | 0.243 | 0.473 | 0.237 | 16 |
| dup/loss high, Ne 1M | 0.776 | 0.247 | 0.405 | 0.549 | 0.356 | 18 |
| dup/loss low, Ne 2M | 0.137 | 0.208 | 0.251 | 0.557 | 0.228 | 16 |
| dup/loss medium, Ne 2M | 0.148 | 0.195 | 0.249 | 0.550 | 0.232 | 16 |
| dup/loss high, Ne 2M | 0.817 | 0.230 | 0.368 | 0.585 | 0.316 | 16 |

Table: {#tab:sisters} Reticulation sister distance on the twelve discordance configurations, plotted in the middle row of @fig:discordance. Each value is the mean over five replicates on the networks that every method completed, and n is the number of those networks. Lower is better.
| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.382 | 0.428 | 0.052 | 0.594 | 0.471 | 17 |
| ILS medium | 0.511 | 0.555 | 0.196 | 0.656 | 0.540 | 16 |
| ILS high | 0.522 | 0.563 | 0.310 | 0.729 | 0.565 | 16 |
| dup/loss low, Ne 200k | 0.482 | 0.524 | 0.156 | 0.574 | 0.464 | 16 |
| dup/loss medium, Ne 200k | 0.505 | 0.558 | 0.171 | 0.597 | 0.485 | 16 |
| dup/loss high, Ne 200k | 0.942 | 0.696 | 0.444 | 0.680 | 0.581 | 17 |
| dup/loss low, Ne 1M | 0.516 | 0.540 | 0.244 | 0.646 | 0.543 | 16 |
| dup/loss medium, Ne 1M | 0.515 | 0.554 | 0.274 | 0.682 | 0.559 | 16 |
| dup/loss high, Ne 1M | 0.925 | 0.671 | 0.529 | 0.725 | 0.618 | 18 |
| dup/loss low, Ne 2M | 0.514 | 0.553 | 0.316 | 0.741 | 0.590 | 16 |
| dup/loss medium, Ne 2M | 0.520 | 0.556 | 0.299 | 0.708 | 0.580 | 16 |
| dup/loss high, Ne 2M | 0.925 | 0.652 | 0.495 | 0.756 | 0.604 | 16 |

Table: {#tab:mu} Mu-distance on the twelve discordance configurations, plotted in the bottom row of @fig:discordance. Each value is the mean over five replicates on the networks that every method completed, and n is the number of those networks. Lower is better.
| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.393 | 0.404 | 0.027 | 0.480 | 0.429 | 17 |
| ILS medium | 0.417 | 0.436 | 0.062 | 0.487 | 0.437 | 16 |
| ILS high | 0.413 | 0.442 | 0.150 | 0.567 | 0.465 | 16 |
| dup/loss low, Ne 200k | 0.402 | 0.423 | 0.052 | 0.447 | 0.398 | 16 |
| dup/loss medium, Ne 200k | 0.410 | 0.445 | 0.057 | 0.466 | 0.419 | 16 |
| dup/loss high, Ne 200k | 0.607 | 0.523 | 0.317 | 0.563 | 0.507 | 17 |
| dup/loss low, Ne 1M | 0.404 | 0.416 | 0.075 | 0.511 | 0.442 | 16 |
| dup/loss medium, Ne 1M | 0.414 | 0.430 | 0.098 | 0.505 | 0.447 | 16 |
| dup/loss high, Ne 1M | 0.567 | 0.477 | 0.351 | 0.567 | 0.497 | 18 |
| dup/loss low, Ne 2M | 0.408 | 0.430 | 0.144 | 0.574 | 0.475 | 16 |
| dup/loss medium, Ne 2M | 0.421 | 0.437 | 0.132 | 0.575 | 0.472 | 16 |
| dup/loss high, Ne 2M | 0.579 | 0.456 | 0.348 | 0.574 | 0.505 | 16 |

Table: {#tab:frac075} Mild fractionation, retention 0.75, plotted in @fig:fractionation. Each value is the mean over five replicates on the 18 networks that every method completed. Lower is better.
| Method | Ret. descendants | Ret. sisters | Mu-distance |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.160 | 0.474 | 0.430 |
| PlaceNet ploidy-free | 0.233 | 0.519 | 0.460 |
| Polyphest | 0.244 | 0.335 | 0.251 |
| GRAMPA-iter | 0.492 | 0.705 | 0.546 |
| GRAMPA-iter with prior | 0.249 | 0.585 | 0.491 |

Table: {#tab:frac050} Medium fractionation, retention 0.50, plotted in @fig:fractionation. Each value is the mean over five replicates on the 18 networks that every method completed. Lower is better.
| Method | Ret. descendants | Ret. sisters | Mu-distance |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.682 | 0.829 | 0.577 |
| PlaceNet ploidy-free | 0.506 | 0.720 | 0.532 |
| Polyphest | 0.609 | 0.766 | 0.420 |
| GRAMPA-iter | 0.542 | 0.733 | 0.535 |
| GRAMPA-iter with prior | 0.672 | 0.836 | 0.562 |

Table: {#tab:frac025} Severe fractionation, retention 0.25, plotted in @fig:fractionation. Each value is the mean over five replicates on the 16 networks that every method completed. Lower is better.
| Method | Ret. descendants | Ret. sisters | Mu-distance |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.922 | 0.964 | 0.548 |
| PlaceNet ploidy-free | 0.876 | 0.942 | 0.538 |
| Polyphest | 0.917 | 0.948 | 0.462 |
| GRAMPA-iter | 0.602 | 0.774 | 0.547 |
| GRAMPA-iter with prior | 0.923 | 0.962 | 0.592 |
