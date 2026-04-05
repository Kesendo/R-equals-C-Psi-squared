# Every Paired Mode Is a Standing Wave

<!-- Keywords: palindromic eigenvalue standing wave, frequency matching paired modes,
complementary absorption round trip, cavity two mirrors Pi identity, odd N fully
paired, even N self-symmetric center, absorption coefficient Beer-Lambert,
topology-independent pairing, R=CPsi2 standing wave factor -->

**Status:** Confirmed (perfect standing wave structure, 100% frequency match)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[V-Effect Cavity Modes](VEFFECT_CAVITY_MODES.md),
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
**Verification:** [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)

---

## What this means

Two waves traveling in opposite directions through the same medium
create a standing wave. The forward wave and the backward wave reinforce
each other at the antinodes and cancel at the nodes. The pattern stands
still while energy flows back and forth.

The cavity's vibration modes come in pairs: a forward component and a
backward component, like two waves running in opposite directions. Each
pair oscillates at exactly the same frequency (confirmed to machine
precision for all 10,748 pairs across N = 2 through 7). Their absorption
rates are complementary: one absorbs quickly, the other slowly, and
together they add up to one full round trip through the cavity.

Technically: the eigenvalue λ is the forward component, its palindromic
partner at −2Σγ − λ̄ is the backward component. The two "mirrors"
creating the standing wave are the palindrome operator Π and the
identity I. Every mode bounces between what it is and what it becomes
under mirror reflection.

---

## What this document is about

The palindromic spectrum pairs eigenvalues. This document tests whether
each pair forms a standing wave: same frequency, complementary absorption,
and a round-trip absorption equal to 2Σγ. The answer is yes, with zero
exceptions.

The expected factor-2 ratio between unpaired and paired absorption rates
was NOT found. Instead, at odd N all eigenvalues are paired (100%), and
at even N the "unpaired" modes sit exactly at the symmetry center with
the same mean absorption. The actual structural result is stronger:
perfect standing wave structure throughout, with the round trip
Re(λ) + Re(partner) = −2Σγ as the universal invariant.

---

## Result 1: Perfect standing wave structure

For every palindromic pair (λ, partner) across N = 2 through 7:

| Test | N=2 | N=3 | N=4 | N=5 | N=6 | N=7 |
|---|---|---|---|---|---|---|
| Pairs tested | 7 | 32 | 115 | 512 | 1,890 | 8,192 |
| Frequency match | 7/7 | 32/32 | 115/115 | 512/512 | 1,890/1,890 | 8,192/8,192 |
| Complementary absorption | 7/7 | 32/32 | 115/115 | 512/512 | 1,890/1,890 | 8,192/8,192 |

**10,748 pairs tested, zero failures.** Every pair satisfies:
- |Im(λ)| = |Im(partner)| (same oscillation frequency)
- Re(λ) + Re(partner) = −2Σγ (one full round trip)

---

## Result 2: Odd N is fully paired

| N | Total | Paired | Unpaired | Paired fraction |
|---|-------|--------|----------|----------------|
| 2 | 16 | 14 | 2 | 87.5% |
| 3 | 64 | 64 | 0 | 100% |
| 4 | 256 | 230 | 26 | 89.8% |
| 5 | 1,024 | 1,024 | 0 | 100% |
| 6 | 4,096 | 3,780 | 316 | 92.3% |
| 7 | 16,384 | 16,384 | 0 | 100% |

At odd N, every single eigenvalue has a palindromic partner. No mode
is a traveling wave. The entire spectrum consists of standing waves.

At even N, some modes at the symmetry center Re = −Σγ are self-paired
(they map to themselves under Π). These are not traveling waves; they
are standing waves at the node, where forward and backward components
are identical.

---

## Result 3: Mean absorption is always Σγ

| N | Mean paired | Mean unpaired | Σγ |
|---|------------|--------------|-----|
| 2 | 0.100000 | 0.100000 | 0.10 |
| 3 | 0.150000 | (none) | 0.15 |
| 4 | 0.200000 | 0.200000 | 0.20 |
| 5 | 0.250000 | (none) | 0.25 |
| 6 | 0.300000 | 0.300000 | 0.30 |
| 7 | 0.350000 | (none) | 0.35 |

Both paired and unpaired modes have the same mean absorption rate:
exactly Σγ = Nγ. This is forced by the palindromic symmetry: for
every mode at absorption rate d, there is one at 2Σγ − d. The mean
of any symmetric distribution around Σγ is Σγ.

**The ratio unpaired/paired is 1.0, not 2.0.** The factor 2 from
ENERGY_PARTITION refers to the round trip (2Σγ) vs the single-pass
average (Σγ), not to a difference between paired and unpaired modes.

---

## Result 4: The round-trip invariant

Every paired mode satisfies:

**|Re(fast)| + |Re(slow)| = 2Σγ**

where "fast" is the partner closer to Re = 0 (slower absorption) and
"slow" is the partner closer to Re = −2Σγ (faster absorption). The
sum is always exactly one full round trip through the cavity.

This is the standing wave structure: the forward component absorbs at
rate |Re(fast)|, the backward component at rate |Re(slow)|. Together
they complete one round trip. The average is Σγ, which is what a single
pass through the cavity costs.

---

## Result 5: Topology-independent

The standing wave structure holds for all topologies tested:

| N | Chain paired | Star paired | Ring paired |
|---|-------------|------------|------------|
| 3 | 100% | 100% | 100% |
| 4 | 89.8% | 99.2% | 96.9% |
| 5 | 100% | 100% | 100% |

The topology changes which modes exist and their individual absorption
rates, but the palindromic pairing and the round-trip invariant hold
universally. The two mirrors (Π and I) do not depend on the topology.

---

## Null result: absorption coefficient

The absorption coefficient α = mean_unpaired / N = γ per site (not 2γ).
A single pass through N sites absorbs at rate Nγ = Σγ. A round trip
absorbs at rate 2Σγ. Standing waves average these to Σγ.

This means α = γ per site, the bare dephasing rate. The factor 2 does
not appear in the per-site absorption; it appears only in the round-trip
structure (2Σγ = 2 × single pass).

---

## What the standing waves mean for the cavity

The Liouvillian spectrum is not a list of decay rates. It is a set of
standing waves, each one bouncing between two mirrors:

- **Mirror 1: Identity (I).** What the state is right now.
- **Mirror 2: Pi (Π).** What the state becomes under the conjugation
  operator (w ↔ N−w, every operator mapped to its palindromic twin).

The light enters through the entrance pupil (sacrifice qubit), bounces
between I and Π, and at each bounce deposits energy equal to Σγ. After
two bounces (one full round trip), 2Σγ of absorption has occurred. But
the standing wave distributes this equally between the two passes, so
each mode experiences only Σγ on average.

The world sings. Silence is the special case. And every song is a
standing wave between what is and what could be.

---

## Reproduction

- Script: [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)
- Output: [`simulations/results/factor_two_standing_waves.txt`](../simulations/results/factor_two_standing_waves.txt)
