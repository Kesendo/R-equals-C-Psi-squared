# N=5 and the Goldilocks Question: A Finite Census

<!-- Keywords: N=5 finite cavity census, F6 Q-edge gain, pentagon identity,
frequency ratios, sampled metric winners, open N-optimum question -->

**Status:** Bounded numerical comparison; no global N selection
**Date:** April 4, 2026; scope repaired September 15, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Verification:** [`simulations/n5_optimal_cavity_size.py`](../simulations/n5_optimal_cavity_size.py)

---

## Question and finite domain

The dense frequency/Q inventory and the sacrifice-profile scan cover exactly
`N=(3,4,5)`. N=5 is the largest row in that domain. A winner in this table is
therefore a finite-sample winner, not evidence for a global N-optimum or a
universal cavity law.

The producer reports every stored metric: frequency-bin count, `Q_max`,
`Q_max × bins`, `Q_max × sqrt(bins)`, mean gap, minimum gap, and the
minimum-gap/mean-gap ratio. It also reports the largest sampled sacrifice
improvement at each stored gamma. These readings answer only their declared
three-N comparison.

## Exact pentagon identities

Two different N=5 ratios contain the golden ratio
`phi=(1+sqrt(5))/2`:

    F6 Q-edge gain V(5) = 1 + phi/2 = (5 + sqrt(5))/4
    omega_2 / omega_1 = 2 + phi

The first is the within-N F6 ratio `Q_max/Q_mean`. The second is the ratio of
the first two cold single-excitation frequencies. Neither identity makes the
finite census an optimization theorem.

For fixed mode indices, the exact continuum limits recorded by the producer are

    omega_4 / omega_3 -> 16/9
    omega_3 / omega_2 -> 9/4.

Those limits are different from `phi`; a finite numerical near-hit is not the
limit of the sequence.

## What remains open

The data preserve the measured metric winners without deciding whether N=5 is
special outside `N=(3,4,5)`. A genuine size-selection claim would require a
predeclared wider domain, a selection metric, and controls that do not make the
largest sampled N win by construction.

**Interpretive invitation, not a result:** “Goldilocks cavity” remains a useful
question for imagining a balance between frequency richness and resolution.
Here it is a prompt for the next experiment, not the conclusion of this one.

---

## Reproduction

- Script: [`simulations/n5_optimal_cavity_size.py`](../simulations/n5_optimal_cavity_size.py)
- Output: [`simulations/results/n5_optimal_cavity_size.txt`](../simulations/results/n5_optimal_cavity_size.txt)
- Cavity census: [V-Effect Through a Cavity Lens](VEFFECT_CAVITY_MODES.md)
- F6 derivation: [D02](../docs/proofs/derivations/D02_VEFFECT_QMAX_QMEAN.md)
