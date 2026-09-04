# Every Paired Mode Is a Standing Wave

<!-- Keywords: palindromic eigenvalue standing wave, palindromic pairing fraction,
complementary absorption round trip, cavity two mirrors Pi identity, odd N fully
paired, even N self-symmetric center, absorption coefficient Beer-Lambert,
topology-independent pairing, R=CPsi2 standing wave factor -->

**Status:** Confirmed (the spectrum is completely paired, at every N)
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
pair oscillates at exactly the same frequency, and their absorption
rates are complementary: one absorbs quickly, the other slowly, and
together they add up to one full round trip through the cavity. Both of
those follow from where the partner sits, so the question this document
actually answers is the other one: how much of the spectrum is paired at
all. Across N = 2 through 7 the answer is all of it, 21,840 eigenvalues
with no exception.

Technically: the eigenvalue λ is the forward component, its palindromic
partner at −2Σγ − λ̄ is the backward component. The two "mirrors"
creating the standing wave are the palindrome operator Π and the
identity I. Every mode bounces between what it is and what it becomes
under mirror reflection.

---

## What this document is about

The palindromic spectrum pairs eigenvalues, and a pair reads as a standing
wave: same frequency, complementary absorption, a round trip of 2Σγ. Those
three properties come with the partner's position and are not in question
here (Result 1). What this document measures is how much of the spectrum is
paired, and what the rest does.

The answer is all of it, at every N. 21,840 eigenvalues across N = 2
through 7, every one carrying a palindromic partner: 9,921 pairs of
distinct partners and 1,998 modes that are their own partner, sitting on
the fixed locus Re = −Σγ. There is no traveling wave anywhere in the
range. This is the palindrome theorem read off the spectrum rather than
proved, and it agrees with it.

The factor-2 ratio this document set out to find between "unpaired" and
paired absorption is not there, because it was never a statement about
pairing. See Result 3.

---

## Result 1: The standing wave relations are the pairing map

The palindromic partner of λ sits at −2Σγ − λ̄. Two relations follow from
that expression alone, for any spectrum at all:

- |Im(λ)| = |Im(partner)|. The map conjugates, so it preserves frequency.
- Re(λ) + Re(partner) = −2Σγ. The map reflects the real axis about −Σγ.

These are the standing wave relations, and they are algebra rather than
measurement. The search that finds a partner enforces them: it accepts the
eigenvalue nearest −2Σγ − λ̄ within 1e-8, and that one distance bounds both
the frequency difference and the deviation of the real-part sum. A pair
violating either relation is not a pair the search can return, so counting
how many pairs satisfy them counts pairs found, not tests passed.

The sum rule Re(λ) + Re(partner) = −2Σγ is nonetheless a real result about
the spectrum; it is proven independently as a one-line corollary of the
absorption theorem and the palindromic weight swap, in
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §4.2.

What this document measures is whether a partner exists at all, which is
Result 2, and how many pairs there are:

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| Pairs of distinct partners | 3 | 32 | 52 | 512 | 1,130 | 8,192 |
| Self-paired modes | 10 | 0 | 152 | 0 | 1,836 | 0 |

9,921 pairs and 1,998 self-paired modes, which accounts for all 21,840
eigenvalues exactly.

---

## Result 2: Every N is fully paired

| N | Total | With a distinct partner | Self-paired at Re = −Σγ | Paired |
|---|-------|------|------|------|
| 2 | 16 | 6 | 10 | 100% |
| 3 | 64 | 64 | 0 | 100% |
| 4 | 256 | 104 | 152 | 100% |
| 5 | 1,024 | 1,024 | 0 | 100% |
| 6 | 4,096 | 2,260 | 1,836 | 100% |
| 7 | 16,384 | 16,384 | 0 | 100% |

Every eigenvalue has a palindromic partner, at even N as well as odd. No
mode is a traveling wave anywhere in this range.

Even and odd differ in how the partner is reached, not in whether one
exists. At odd N every mode has a partner other than itself. At even N a
large share of the spectrum sits on Re = −Σγ, the fixed locus of the
palindrome map, where a mode IS its own partner: 10 of 16 at N=2, 152 of
256 at N=4, 1,836 of 4,096 at N=6, each on the locus to better than 1e-13.
Those are standing waves at the node, where the forward and backward
components are identical.

A count that treats a self-paired mode as unpaired sees a deficit at even
N that is not in the spectrum. What such a count reports is the number of
odd-sized degeneracy clusters on the locus, 2 at N=2, 26 at N=4 and 316 at
N=6, a property of how ties are broken rather than of the palindrome.

---

## Result 3: Mean absorption is always Σγ

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| Mean absorption | 0.100000 | 0.150000 | 0.200000 | 0.250000 | 0.300000 | 0.350000 |
| Σγ | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 |

The mean absorption rate over the spectrum is exactly Σγ = Nγ,
forced by the palindrome: for every mode at rate d there is one at
2Σγ − d, so the distribution is symmetric about Σγ. The
self-paired modes sit on Σγ individually, being the fixed locus, so
their mean is Σγ for the same reason read at a point.

**The factor 2 is not a ratio between two sets of modes.** It is the ratio
of the full range (0 to 2Σγ) to the centre (Σγ), which is how
[Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) states it. The
all-light modes, the ones at XY-weight N, do sit at 2Σγ, twice the
mean, and there are exactly N+1 of them at each N here. They are not
unpaired: they are the fast end of pairs whose slow end is the stationary
set at Re = 0.

---

## Result 4: The round-trip invariant

Every paired mode satisfies:

**|Re(fast)| + |Re(slow)| = 2Σγ**

where "fast" is the partner closer to Re = 0 (slower absorption) and
"slow" is the partner closer to Re = −2Σγ (faster absorption). This
is Result 1's second relation with the signs written out, so it holds on
any pair by where the partner sits, and it is established as a theorem in
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) section 4.2.

This is the standing wave structure: the forward component absorbs at
rate |Re(fast)|, the backward component at rate |Re(slow)|. Together
they complete one round trip. The average is Σγ, which is what a single
pass through the cavity costs.

---

## Result 5: Topology-independent

The paired fraction was measured on three topologies and is 100% on all of
them at N = 3, 4 and 5. What the topology moves is how the pairing is
reached, not whether it holds:

| N=4 | Chain | Star | Ring |
|---|---|---|---|
| Self-paired at Re = −Σγ | 152 | 126 | 156 |
| Pairs of distinct partners | 52 | 65 | 50 |
| Paired | 100% | 100% | 100% |

The topology changes which modes exist, their individual absorption rates,
and how much of the spectrum sits on the fixed locus. The palindromic
pairing itself does not depend on it, and neither do the two mirrors
(Π and I).

---

## Null result: absorption coefficient

The absorption coefficient α = Σγ / N = γ per site (not
2γ), the mean absorption over the spectrum divided by the number of
sites. A single pass through N sites absorbs at rate Nγ = Σγ. A
round trip absorbs at rate 2Σγ. Standing waves average these to Σγ.

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
