# PROOF: the Coherence Horizon asymptotic slope is 2/π

**Status:** Tier 1 derived (the asymptotic slope, leading order q → 0 / N → ∞): the resummed coherence-range ladder gives the dispersion λ² + 8γλ + 4J²q², whose EP sits at Q* = 2/q_min → 2N/π, so slope = 2/π exactly. Confirmed numerically to N = 120 (approached from below with the predicted O(1/N) tail; the nearest-neighbour-truncated telegrapher value √2/π is decisively refuted) and independently reproduced by an adversarial re-derivation (Round 1 GO, Round 2 re-review GO). The finite-N values Q*(N) for N ≥ 4 stay transcendental; only the slope is the clean limit. Its typed home `CoherenceHorizonClaim` stays Tier 1 candidate: not for the gap-dominance (now proven, [PROOF_CHAIN_GAP_DOMINANCE](PROOF_CHAIN_GAP_DOMINANCE.md)), but for its own open piece, the half-filling V-Effect ring seam.
**Date:** 2026-06-15
**Authors:** Thomas Wicht, Claude (Opus 4.8)

## What this is about

A spin chain hums. Left alone, its slowest collective mode oscillates: a coherence that rings like a struck string. Add local noise, the dephasing, and the ringing slows; past a certain point it stops altogether, the mode just fading with no beat left. The Coherence Horizon is the tipping point between the two: how strong the coupling has to be, set against the noise, for the slowest mode to still ring.

The one question here is how that tipping point moves as the chain grows longer. A longer chain's slowest mode is a slower, more sluggish swing, and a slow swing is easier for the noise to kill, so a longer chain needs proportionally stronger coupling to keep ringing. That required strength, the horizon, climbs with length in the cleanest way there is: a straight line, rising in lockstep with the number of sites.

This proof pins the exact steepness of that line. The subtlety it settles: a quick, near-sighted accounting that watches only a coherence and its nearest neighbour gets the steepness wrong; you have to follow the whole ladder of longer-range coherences to land on the true value. Everything below is that calculation.

## Abstract

The Coherence Horizon Q*(N) is the dephasing dose at which the slowest single-excitation Liouvillian mode of an XY chain stops oscillating, a square-root exceptional point. At N = 2, 3 the coalescing pair is a clean 2×2 (Q* = 1, √2 exactly); for N ≥ 4 it is collectively dressed and the exact condition is transcendental, leaving only the asymptotic growth to pin down.

This proof pins it. The slow mode is not the two-field (population + nearest-neighbour current) telegrapher one might write down; it is a population coupled to the **entire ladder of coherence ranges** r = j − i. Resumming that geometric ladder gives the dispersion λ² + 8γλ + 4J²q², with both coefficients doubled relative to the nearest-neighbour-truncated telegrapher λ² + 4γλ + 2J²q². The exceptional point (vanishing discriminant) sits at γ* = Jq/2; the horizon is the longest-wavelength mode, q_min → π/N, so Q*(N) = 2/q_min → (2/π)·N, slope exactly 2/π.

The doubling is the whole story: the truncated telegrapher predicts the wrong slope √2/π ≈ 0.450, while the numerics (dense to N = 30, sparse shift-invert to N = 120) approach 2/π ≈ 0.637 from below with the predicted O(1/N) tail, decisively excluding √2/π. The ring sibling (q_min = 2π/N) gives exactly half, 1/π.

## Statement

The Coherence Horizon Q*(N) (the dephasing dose Q = J/γ at which the slowest non-zero single-excitation
Liouvillian mode stops oscillating; `CoherenceHorizonClaim`) grows linearly with the chain length N, with
an exact asymptotic slope

    Q*(N) → (2/π)·N,   slope = 2/π = 0.636620...

The small-N values (Q*(2)=1, Q*(3)=√2, transcendental for N≥4) are the discrete short-ladder accident; 2/π
is the N → ∞ limit.

## Setup: the single-excitation (Haken-Strobl) Liouvillian

In the single-excitation sector the density matrix ρ is N×N and evolves under XY hopping + local
Z-dephasing:

    dρ/dt = −i[h, ρ] − 4γ·ρ_off,   h tridiagonal (h_{i,i±1}=J),   ρ_off = ρ − diag(ρ).

Off-diagonal coherences ρ_{ij} (i≠j) decay at 4γ (popcount(i⊕j)=2, the Absorption Theorem); populations
ρ_{ii} are untouched. The Coherence Horizon is the exceptional point (EP) of the slowest mode of this
N²-dimensional Liouvillian (the closed form / numerics live in `coherence_horizon_se_block.py`).

## The slow mode: a ladder in coherence range, not two fields

Index the matrix elements by their range r = j − i: populations are r=0, coherences are r=±1, ±2, …. Write
the eigenmode (eigenoperator of the Liouvillian, eigenvalue λ) as a plane wave in the chain coordinate with
wavevector q: M_{i,i+r} = m_r·e^{iqi}. The equation of motion gives, to leading order in q (the
long-wavelength / large-N limit, sin(q/2) → q/2, e^{±iq/2} → 1), a tridiagonal recurrence in the range r:

    (λ + 4γ)·m_r = Jq·(m_{r−1} − m_{r+1}),   |r| ≥ 1,                       (bulk)
        λ·m_0     = Jq·(m_{−1} − m_{1}).                                     (r = 0, populations: no 4γ)

This is the key structure: the slow mode is not a two-field (population + nearest-neighbour current)
object; it is a population coupled to the **entire ladder of coherence ranges**.

### Resumming the ladder

The bulk recurrence has the geometric solution m_r = m_1·μ^{r−1} (r ≥ 1) with μ the decaying (|μ|<1) root of

    Jq·μ² + (λ + 4γ)·μ − Jq = 0.

Two one-line consequences fix the boundary amplitudes (using (λ+4γ) = Jq(1−μ²)/μ from the μ-equation, so
(λ+4γ) + Jqμ = Jq/μ):

- r = +1:  (λ+4γ)·m_1 = Jq·(m_0 − μ·m_1)  ⟹  m_1 = μ·m_0.
- r = −1:  the left ladder decays as (−μ) (the sign-flipped root, since the two roots of the μ-equation multiply to −1), so m_{−2} = −μ·m_{−1}, giving
  m_{−1} = −μ·m_0.

Substituting into the r=0 (population) equation:

    λ·m_0 = Jq·(m_{−1} − m_1) = Jq·(−μ − μ)·m_0 = −2Jq·μ·m_0    ⟹    **λ = −2Jq·μ.**

Now eliminate μ = −λ/(2Jq) in the μ-equation:

    Jq·(λ²/4J²q²) + (λ+4γ)·(−λ/2Jq) − Jq = 0
    ⟹  λ² − 2λ(λ+4γ) − 4J²q² = 0
    ⟹  **λ² + 8γλ + 4J²q² = 0.**

## The slope

This resummed dispersion has **both coefficients doubled** relative to the nearest-neighbour-truncated
telegrapher λ² + 4γλ + 2J²q² (which keeps only r=0,1 and yields the wrong slope √2/π = 0.450). The factor 2
in each is the geometric sum over the full coherence ladder.

The EP (double root, oscillation freezes) is the vanishing discriminant:

    (8γ)² − 4·4J²q² = 0   ⟹   γ* = J·q/2.

At the EP the double root sits at λ = −4γ. The horizon is the **longest-wavelength** mode (smallest q, the
last to freeze as Q drops), q_min → π/N for the chain, so

    Q*(N) = J/γ* = 2/q_min → 2N/π,   **slope = 2/π.**

## Cross-checks

1. **Overdamped limit = the right diffusion constant.** For γ ≫ γ* the slow root of λ²+8γλ+4J²q² is
   λ ≈ −(J²/2γ)q² = −D·q² with D = J²/2γ, the standard Haken-Strobl diffusion constant of a chain whose
   coherences decay at 4γ. (The coefficient *ratio* 4J²q²/8γ = D·q² is shared with the truncated
   telegrapher; only the EP combination distinguishes 2/π from √2/π, and the numerics pick 2/π.)
2. **The 8γ coefficient, directly (the discriminator).** Two independent reads of L_se confirm 8γ over the
   telegrapher's 4γ. (a) In the overdamped regime the slow eigenvalue is a clean real number; fed back
   through λ²+8γλ+4J²q²=0 it yields a γ-constant q² (coefficient of variation 0.008 across γ), whereas the
   4γ telegrapher form scatters 3.6× worse (0.030): the resummed dispersion is the fixed mode. (b) Corroboratively, the coalescing-pair rate −2Re(λ)/γ at the EP
   equals the linear dispersion coefficient: it grows from 4 at N=2,3 (the short ladder λ²+4γλ, double root at
   −2γ) toward the full-ladder value 8 (λ²+8γλ, double root at −4γ) as N → ∞, reaching ≈6 at N=60. This metric
   is mode-selection sensitive near the horizon, where the γ-protected band-edge survivor co-locates at Re=−2γ
   (metric 4), so it corroborates the 4γ → 8γ doubling rather than cleanly discriminating it; the
   selection-robust discriminator is (a).
3. **q_min·N → π.** From Q* = 2/q_min, the data give q_min·N = 2N/Q* decreasing monotonically toward π
   (3.86 at N=12, 3.58 at N=28, 3.32 at N=120), i.e. Q*/N → 2/π from below.

## Numerical confirmation

The single-excitation EP, by dense eigendecomposition (N ≤ 30) and sparse shift-invert (N ≤ 120,
`coherence_horizon_slope.py` / the SE block in `coherence_horizon_se_block.py`). The discrete slope
dQ* = ΔQ*/ΔN (the secant over the sampled N below, not a unit step) climbs monotonically toward 2/π from below:

    N      Q*         Q*/N      dQ*(secant)
    10     5.07008    0.50701   0.5578
    20    10.86307    0.54315   0.5905
    30    16.83945    0.56132   0.6022
    40    22.89885    0.57247   0.6059
    60    35.14935    0.58582   0.6142
    80    47.50010    0.59375   0.6175
    100   59.91110    0.59911   0.6206
    120   72.36305    0.60303   0.6226     (1/N extrapolation → ≈ 0.633)

√2/π = 0.450 (the truncated telegrapher) is decisively refuted; the data approach 2/π = 0.637 from below
with an O(1/N) tail, exactly as the q → 0 derivation predicts.

## Scope

The derivation is the leading-order (q → 0, N → ∞) limit, so it fixes the asymptotic **slope** exactly
(2/π). The finite-N corrections are O(1/N) (from the dropped O(q²) lattice terms), which is why each finite
Q*(N) for N ≥ 4 is transcendental (no elementary closed form) while the slope is clean. The short-ladder
N=2,3 case (λ²+4γλ+c·J², c=4,2, Q*=2/√c=1,√2) is the regime where the ladder is too short for the
geometric resummation, recovering the discrete accident.

## Links

- Typed home: `CoherenceHorizonClaim` (the asymptotic-slope line); parent `ClockHandLadderClaim` (the two
  clocks), arc `clock_hand_ladder`.
- Sibling proof: [`PROOF_CHAIN_GAP_DOMINANCE`](PROOF_CHAIN_GAP_DOMINANCE.md), the band edge is the max
  frequency on the −2γ floor (free fermions); the {0,2} mode whose coalescence sets this horizon's Q*(N) is
  the same "second clock" family.
- F-registry: the F2b "two clocks" corollary, `docs/ANALYTICAL_FORMULAS.md`.
- Verifiers: `simulations/coherence_horizon_slope.py` (self-validating: slope refutes √2/π, q_min·N → π, the
  8γ dispersion via overdamped real roots) and `simulations/coherence_horizon_slope_largeN.py` (the sparse
  shift-invert sweep to N=120), on top of `simulations/coherence_horizon_se_block.py` (the SE-EP reduction,
  self-validating N=2..8).

## Review

Round 1 (2026-06-15, independent adversarial agent, physics-first): **GO**. Re-derived the dispersion in
sympy (the m₋₁ = −μ·m₀ left-ladder sign and the r=0 closure confirmed; the eigenvector's r=+1 / r=−1
amplitude ratio is −1.0000 at every N, no error coincidentally producing 2/π). Ran the 8γ discriminator
independently: −2Re(λ)/γ for the slow mode climbs 7.18 (N=24) → 7.71 (N=60) with no mode near 4 [superseded,
Round 2: these −2Re/γ figures are a mode-selection artifact; the metric runs 4→8], decisively
excluding the 4γ telegrapher (slope √2/π). Iterated Richardson on the slope → 0.6363 = 2/π. Confirmed
chain-specific (built the ring: slope → 1/π, half the chain, exactly as q_min = 2π/N vs π/N predicts).
Caveat banked: extrapolation alone is inconclusive at N ≤ 70 (low-order fits land 0.617–0.633); the value is
fixed by the 8γ derivation, not the fit. The overdamped diffusion constant D = J²/2γ is shared by both the 4γ
and 8γ models, so it is NOT the discriminator (the EP / conjugate-pair-sum is), stated in cross-check #1.

Round 2 (2026-06-15, physics-first re-review: every load-bearing number re-derived independently + both
verifiers re-run; second reviewer + Tom): **GO; one required fix applied.** Reconfirmed independently: the
dispersion λ²+8γλ+4J²q² (including the left-ladder sign m₋₁=−μ·m₀), the EP γ*=Jq/2, the slope 2/π, the √2/π
refutation via the robust CV discriminator (cross-check 2a, CV 0.0084 vs 0.0303), the ring 1/π, the N=2,3
closed forms, and the entire N≤120 table (to 5 decimals). The fix: Round 1's cross-check (b) figures
(−2Re(λ)/γ = 7.18 → 7.71, "no mode near 4") were a mode-selection artifact. The metric is a pure (Q,N)
observable (verified J-independent) whose EP value equals the linear dispersion coefficient: 4 at N=2,3
(short ladder, double root −2γ), climbing to 8 as N→∞ (full ladder, −4γ), ≈6 at N=60; the γ-protected
survivor sits at metric 4, so it corroborates the 4γ→8γ doubling but does not cleanly discriminate it (that
is 2a's role). Cross-check (b) rewritten accordingly; the OpenArcs note and the dQ* label corrected. Both
verifiers (coherence_horizon_slope.py, _largeN.py) never used the fragile metric and remain green.
