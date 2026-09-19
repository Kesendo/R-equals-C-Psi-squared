# The V-Effect Through a Cavity Lens

<!-- Keywords: finite V-Effect census, F6 Q-edge gain, cold warm frequency
neighbourhood, cavity lens, strict nearest-frequency distance, R=CPsi2 -->

**Status:** Finite numerical inventory; interpretive cavity lens
**Date:** April 4, 2026; scope repaired September 15, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Verification:** [`simulations/veffect_cavity_modes.py`](../simulations/veffect_cavity_modes.py)

---

## Two objects that must stay separate

The **V-Effect census** is a finite classification and frequency-bin record.
Its primary N=3 two-term sample has 36 distinct pairs and the exact split
14 hard / 19 soft / 3 truly. Historical frequency tables use their stated
rounding and generator configurations.

The **F6 Q-edge gain** is the within-one-N ratio

    V(N) = Q_max / Q_mean = 1 + cos(pi/N).

“V-Effect gain” is retained only as a historical alias for this F6 ratio. The
ratio neither counts frequencies nor supplies an ancestry relation between
eigenmodes. Conversely, a finite bin census does not derive F6.

## What the cavity computation measures

For each uniform Heisenberg chain, the producer deduplicates nonzero absolute
frequencies at tolerance `1e-6`. It then treats every cold (`gamma=0`) frequency
as a target and asks for its nearest reusable warm (`gamma=0.05`) frequency.
Acceptance is the strict, one-sided condition

    nearest distance < 0.1.

A warm frequency may cover more than one cold target. This is therefore not a
one-to-one assignment and does not compare eigenvectors, spectral projectors,
or invariant subspaces.

| N | Cold bins | Warm bins | Covered at `< 0.1` | Maximum nearest distance | Covered at `< 1e-6` |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 2 | 1/1 | 0.0000000 | 1/1 |
| 3 | 3 | 5 | 3/3 | 0.0016193 | 2/3 |
| 4 | 14 | 47 | 14/14 | 0.0226557 | 6/14 |
| 5 | 43 | 112 | 43/43 | 0.0289709 | 9/43 |

The broad radius covers all cold targets for N=2..5; the `1e-6` control loses
full coverage at N=3, N=4, and N=5. The equality boundary is rejected: distance
exactly `0.1` does not count. Those controls make the result a neighbourhood
statement rather than an identity statement.

At `gamma=0`, the Liouvillian is `-i[H, ·]`, so its eigenvalues are purely
imaginary. For a uniform chain the `(0,1)` block is

    L_(0,1) = -2iJ L_path - 2 gamma Id,

and its frequencies are `4J(1-cos(pi k/N))`. This exact block formula explains
one named family inside the inventory. It does not exhaust the full spectrum or
identify a cold vector with a warm vector.

## Finite topology inventory

The regenerated result also lists chain, star, and ring frequency-bin counts
for N=3,4,5. Same N with a different graph can give a different inventory, and
graphs with the same number of bonds need not agree. Bond count alone is not a
mechanism for the bin count.

For the uniform-chain `(0,1)` block, the largest frequency gives

    Q_max = J mu_max / gamma
          = (2J/gamma) (1 + cos(pi/N)).

Dividing by `Q_mean = 2J/gamma` gives F6. This statement concerns the Q edge of
that block; it does not turn the full finite frequency inventory into a theorem
about a universal cavity.

## Interpretive boundary

**Interpretive invitation — not a result:** a flute or soundbox is a useful way
to picture how changing a finite generator changes its available frequency
bins. The computation licenses the table and the block formula above. It does
not establish that light reveals pre-existing individual modes, that a larger
system replaces a smaller physical system, or that connection causes life,
complexity, or irreversible creation.

---

## Reproduction

- Script: [`simulations/veffect_cavity_modes.py`](../simulations/veffect_cavity_modes.py)
- Raw output: [`simulations/results/veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt)
- Related finite census: [V-Effect Palindrome](V_EFFECT_PALINDROME.md)
- Exact Q ratio: [D02](../docs/proofs/derivations/D02_VEFFECT_QMAX_QMEAN.md)
