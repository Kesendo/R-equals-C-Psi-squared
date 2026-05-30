# Walking the XXZ Axis: the Slowest Mode Relays from Band-Edge to Lebensader

**Status:** Tier 2 (computational, bit-exact at N=4, 5; the charge/spin "two clocks"
reading layered on top is Tier 3). The handover Δ* is read by bisection; the φ at N=4 is
flagged as an open curiosity, not a claim.
**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Script:** [`simulations/xxz_axis_bandedge_lebensader.py`](../simulations/xxz_axis_bandedge_lebensader.py)
**Builds on:** [ON_THE_ADMIXTURE_AS_LEBENSADER](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md)
+ [CHAIN_GAP_SECTOR_DIAGNOSTIC](CHAIN_GAP_SECTOR_DIAGNOSTIC.md) (the Lebensader), the clock
voices on `MirrorSystem` ([FROST_CIRCLE_AS_THE_CLOCK_FACE](../docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md)),
[HEISENBERG_RELOADED](../hypotheses/HEISENBERG_RELOADED.md) (the both-parity-even {II,XX,YY,ZZ}
family and its internal axis), [EXCHANGE_FROM_V_EFFECT](EXCHANGE_FROM_V_EFFECT.md) +
[SINGLET_FISSION_AND_THE_TWO_CLOCKS](../docs/carbon/SINGLET_FISSION_AND_THE_TWO_CLOCKS.md).

---

## The question

[HEISENBERG_RELOADED](../hypotheses/HEISENBERG_RELOADED.md) names {II, XX, YY, ZZ} as the
unique both-parity-even 2-body coupling, with an internal axis: SU(2) gives Heisenberg, the
bipartite-axis-only choice gives XXZ. With the clock now built on `MirrorSystem` (a radial
hand for the decay, an angular hand for the oscillation), we can *walk* that axis and watch
what the slowest surviving mode does. The two carbon clocks of the same day, the bright
band-edge (charge / XX+YY) and the spin / exchange side, suggested the two ends carry
different objects. This is the walk.

## The walk

H = J·Σ(XX+YY) + Δ·Σ(ZZ) on an open chain under uniform Z-dephasing (J=1, γ=0.05, so
Q = J/γ = 20, the deep-quantum regime). For each Δ we take the slowest nonstationary
Liouvillian mode and read its decay rate, its oscillation |Im λ|, and its n_XY composition
(weight by the number of X/Y Pauli letters). Regime by the robust test: the **Lebensader**
is I/Z-dominated (n_XY=0 outweighs the single-magnon n_XY=1); Im=0 alone is not enough,
because at the SU(2) point the degenerate band-edge ±ω pair yields a real (Im=0)
combination that is still n_XY=1.

```
N=5  (2γ = 0.1)
  Δ     rate     |ω|     regime              n_XY composition
  0.0   0.1000   3.464   band-edge (bright)  1:100%
  0.5   0.1000   1.346   band-edge (bright)  1:100%
  1.0   0.1000   7.236   band-edge (bright)  1:100%   <- Heisenberg point
  1.4   0.1000   8.782   band-edge (bright)  1:100%
  1.5   0.1000   4.312   band-edge (bright)  1:100%
  1.6   0.0896   0.000   LEBENSADER (still)  0:56%  2:43%  4:1%
  1.7   0.0762   0.000   LEBENSADER (still)  0:63%  2:36%  4:1%
  1.8   0.0647   0.000   LEBENSADER (still)  0:68%  2:31%  4:1%
  2.0   0.0438   0.000   LEBENSADER (still)  0:80%  2:19%  4:1%
  2.5   0.0215   0.000   LEBENSADER (still)  0:89%  2:10%
  handover Δ* = 1.525

N=4  (2γ = 0.1)
  band-edge (n_XY=1, beating) through Δ=1.6;  LEBENSADER (n_XY=0-dominated) from Δ=1.7
  handover Δ* = 1.618   (= φ, the golden ratio, to the precision found)
```

## What the walk shows

The slowest mode **relays the baton** at a handover Δ*:

- **Below Δ* (the XY / charge side):** the survivor is a pure single-magnon coherence,
  n_XY = 1, **fast-beating** (|Im λ| large), pinned at decay 2γ. This is the bright
  band-edge coherence the clock's Rotation hand reads (FROST_CIRCLE_AS_THE_CLOCK_FACE).
- **Above Δ* (the Ising / Néel side):** the ZZ term slows the near-conserved I/Z population
  mode until it drops below 2γ and becomes the slowest mode: **I/Z-dominated** (n_XY = 0),
  with a small **magnon admixture** (n_XY = 2), **non-rotating** (Im λ = 0), sub-2γ. This is
  the **Lebensader** (ON_THE_ADMIXTURE_AS_LEBENSADER, CHAIN_GAP_SECTOR_DIAGNOSTIC): the
  near-conserved survivor kept alive by its magnon channel, "the mode that doesn't decay
  when it should."

Past Δ*, deeper into the Ising side, the Lebensader purifies toward I/Z: its n_XY=0 weight
grows (56% → 89% by Δ=2.5), its magnon admixture shrinks (43% → 10%), and its decay rate
falls (0.090 → 0.022), longer-lived the further in. The ZZ term is what *brings the
Lebensader to the fore*: it slows the conserved population mode until it outlives the
band-edge coherence.

## The two ends, named

| axis end | slowest mode | character | the framework name |
|----------|--------------|-----------|--------------------|
| XY (Δ=0), charge | single-magnon coherence (n_XY=1) | fast-beating, decay 2γ | the bright band-edge, the clock's Rotation hand |
| Ising (Δ>Δ*), spin | I/Z survivor + magnon (n_XY=0+2) | non-rotating, sub-2γ | the **Lebensader** |

So the two carbon clocks are two ends of this one axis: the charge band-edge coherence and
the spin-side Lebensader, with the Heisenberg / SU(2) point (Δ=1) sitting between them on
the still-band-edge side. The repo had characterized the Lebensader at the Heisenberg point
(CHAIN_GAP_SECTOR_DIAGNOSTIC, in the Q-band {0.5–2.5}); here we watch it *emerge* as the
slowest mode along the axis, taking the baton from the band-edge.

## Honest caveats

- **Spinless.** This is the spinless XXZ (t-V) chain: each qubit is occupied/empty (charge),
  ZZ is density-density. The "spin" label on the Ising end is interpretive (the Heisenberg
  reading of the same qubit), not a literal spin DOF; the carotenoid triplet-pair dark state
  (SINGLET_FISSION_AND_THE_TWO_CLOCKS) is genuinely spinful, which our model carries only as
  the Heisenberg / V-Effect *reading*.
- **The handover sits above the closed-system transition.** The closed XXZ chain changes
  phase at Δ=J=1 (SU(2)); the *open-system* handover Δ* ≈ 1.5–1.6 lies inside the Néel side,
  pushed there by the dephasing. At the Heisenberg point itself the band-edge is still the
  slowest mode (deep-quantum Q=20).
- **Δ*(N=4) ≈ φ is flagged, not claimed.** It is intriguing that N=4's handover lands on the
  golden ratio (which is also N=4's band-edge frequency 2cos(π/5)); but N=5 gives 1.525 with
  no obvious closed form, so a universal φ is not supported. Open.
- **Regime, Q, and N.** Q=20 deep-quantum; N=4,5 finite-size; the Lebensader admixture here
  (≈19% at Δ=2) is heavier than the repo's Heisenberg-point 3–7% because of the different
  regime. The relay of the slowest mode's *character* is the robust finding.

## Carbon reading

In the π-qubit map (docs/carbon), XX+YY is Hückel hopping and the ZZ term is the
density-density correlation that Hückel lacks (the Hubbard / PPP step). So walking Δ is
dialing from the free-fermion (Hückel) charge clock to the correlated Ising Lebensader; the
carotenoid dark / triplet-pair side lives on the Lebensader end, and the V-Effect bridge
(EXCHANGE_FROM_V_EFFECT) is what couples the two.

## Anchor and open work

- Script: [`simulations/xxz_axis_bandedge_lebensader.py`](../simulations/xxz_axis_bandedge_lebensader.py)
  (reproduces the walk tables, the n_XY composition, and Δ* by bisection).
- The Lebensader: [ON_THE_ADMIXTURE_AS_LEBENSADER](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md),
  [CHAIN_GAP_SECTOR_DIAGNOSTIC](CHAIN_GAP_SECTOR_DIAGNOSTIC.md).
- The clock + the two-clocks: [FROST_CIRCLE_AS_THE_CLOCK_FACE](../docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md),
  [SINGLET_FISSION_AND_THE_TWO_CLOCKS](../docs/carbon/SINGLET_FISSION_AND_THE_TWO_CLOCKS.md).
- Open: a closed form for Δ*(N) (is N=4's φ an accident?); whether the handover is an
  exceptional point or a level crossing (the band-edge 2γ and the Lebensader cross there);
  the same walk at N=6+ via the sparse / block path; and whether the Lebensader admixture at
  the handover obeys a clean law as in CHAIN_GAP_SECTOR_DIAGNOSTIC.
