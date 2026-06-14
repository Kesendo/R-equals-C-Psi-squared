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
- **Resolved (2026-06-14, the handover-Q unification; arc `xxz_axis_handover`):**
  - **EP or level crossing?** A **level crossing**: the frozen Lebensader (|Im|≈1e-15) meets the
    oscillating band edge (|Im|≈9) and they cross in Re, they do not coalesce. The Δ-handover is the
    *same event* as the dephasing-axis Q-handover
    ([`handover_q.py`](../simulations/carbon/handover_q.py), the typed `HandoverFloorClaim`): the
    interior/Lebensader darkness ⟨n_XY⟩ crossing 1 = the Absorption-Theorem band-edge floor 2γ.
    Verified [`simulations/xxz_handover_unification.py`](../simulations/xxz_handover_unification.py) —
    the band edge sits at *exactly* 2γ for all Δ (|vac⟩⟨magnon| is an eigenoperator of [H,·]; the ZZ
    shifts only Im, so F50's 2N *count* breaks for Δ≠1 but the *floor* persists).
  - **The N=6+ walk:** done via dead-centre sector reduction
    ([`simulations/xxz_delta_star.py`](../simulations/xxz_delta_star.py), validated bit-exact vs the
    full L): Δ*(6)=1.381, Δ*(7)=1.325.
  - **A closed form for Δ*(N)? No clean elementary form, and the limit is the critical point**
    (resolved to N=14, 2026-06-14). Pushed via the γ→0 reduction (Δ* ⟺ gap(R) = 2, where R is the
    Z-coupled classical rate matrix among the half-filling XXZ eigenstates — a Pauli / Fermi-golden-rule
    relaxation), in the new self-validating verifier
    [`simulations/xxz_delta_star_descent.py`](../simulations/xxz_delta_star_descent.py) (γ·gap(R)
    reproduces the full-Liouvillian Lebensader rate as γ→0, ratio→1). The **γ→0** sequence (the
    physical fit target; the earlier Q=20 numbers sit below it, by a drift that grows with N) is
    **monotone decreasing**: Δ*(4..14) = 1.61961, 1.52798, 1.38463, 1.33007, 1.27243, 1.24738,
    1.21578, 1.19958, 1.17933, 1.16827, 1.15389.
    - **The verdict: Δ*(N) → Δ = 1 (the SU(2)/Heisenberg point, the closed-system critical point),
      from the Néel side, consistent with EXACTLY 1.** A free-exponent fit Δ* = L + a·N^(−α) per
      parity gives L just *above* 1 (even ≈ 1.02, odd ≈ 1.05); a fixed-1/N ansatz gives L just *below*
      1 (even ≈ 1.00, odd ≈ 0.98); the two forms **bracket Δ = 1** (L ∈ [0.98, 1.05]). Every computed
      Δ*(N ≤ 14) stays > 1 (no finite-N crossing; the fixed-1/N form would dip below 1 only at
      N ≈ 100–450, far beyond the data). The earlier 4-point ambiguity (a 1/N fit → ≈ 0.85, a 1/N²
      fit → ≈ 1.15) is **collapsed onto the critical point**. (The curve_fit covariance bars, ±0.001
      to ±0.003, are statistical-only; the true uncertainty is the fit-window/fit-form spread above.)
    - **Even/odd:** the γ→0 sequence is monotone (no raw zigzag); the even and odd subsequences are
      two smooth approaches that *both* land on Δ = 1. The fitted exponent is **non-universal**
      (α ≈ 1.16–1.73, even vs odd, window-dependent) — consistent with the **marginal/logarithmic
      corrections** characteristic of the Δ = 1 SU(2) point, which is *why* no clean power law fits and
      the limit is best read as "exactly 1 with a log-corrected approach" rather than a clean exponent.
      (A rigorous Bethe-ansatz derivation of the limit and the log structure is the remaining open item.)
    - **φ refuted (in the physical γ→0 regime):** |Δ*(4) − φ| = **1.6e-3** (φ = 2cos(π/5) = 1.61803,
      Δ*(4)_{γ→0} = 1.61961). The often-quoted "1e-4" was a *finite-γ* (Q=20) artifact of the
      evaluation point; 2cos(π/(N+1)) equals φ at N=4 (that *is* the accident) but fails badly for
      N ≥ 5, and 1+1/N fails everywhere.
    - **Mechanism (against the naive intuition).** The naive reading "Néel order sharpens with N, so
      the handover needs higher Δ" predicts Δ* should *rise*; it *falls*. The handover is set by the
      slow mode's darkness (the Z-coupled rate-matrix gap), and the finite-N offset above Δ = 1 is the
      finite Néel correlation length: once N exceeds it, the dissipative handover tracks the
      closed-system quantum critical point.
- **Ring topology (2026-06-14, the periodic twin; arc `xxz_axis_handover`):** the **ring** Δ\*(N) is
  *qualitatively unlike* the chain. Computed with the same γ→0 reduction plus one wrap bond
  ([`simulations/ring_xxz_delta_star_descent.py`](../simulations/ring_xxz_delta_star_descent.py),
  self-validating). It is **non-monotone**: both parities rise to a peak ≈ 1.31–1.33 near N=9–10
  (odd 1.331 at N=9, even 1.308 at N=10) then **descend** through N=14 (odd → 1.295 at N=13, even →
  1.278 at N=14). The ring **crosses above the chain** near N=7–8 (the chain keeps descending toward
  Δ=1, the ring humps up). **N=4 has no handover**: the full half-filling block is *tangent* to the
  floor at the XY point Δ=0 (peak 0.99998·2γ) and a survivor elsewhere — the K₂,₂ special case (the
  reduction is ~1.5% off there, so N=4 is read off the full block). The **N→∞ limit is open at N ≤ 14**
  (a power-law fit to the hump degenerates, α ≈ 33). This **refutes** the chain's premise that the
  dissipative handover tracks the closed-system Δ=1 critical point: on the ring it is a dynamical,
  topology-sensitive scale. The frame ([`ON_THE_ONE_DIAGONAL`](../reflections/ON_THE_ONE_DIAGONAL.md)):
  the floor 2γ is the first rung of the one diagonal popcount(i⊕j) the light touches — universal,
  topology-free; Δ\*(N) is the Hamiltonian's argument about that fixed floor — topology-dependent.
  *The diagonal is one; the climb is many.*
- **Still open:** a rigorous (Bethe-ansatz) derivation that the *chain* N→∞ limit is *exactly* Δ = 1
  and the log-correction structure of its approach; and, for the *ring*, whether the hump returns
  toward Δ=1 or settles above (needs N ≫ 14 or a closed-form account of the hump).
