# PROOF: the ring double-excitation handover slope is √3/(2π)

**Status:** Tier 1 derived (the asymptotic slope, leading order q → 0 / N → ∞): the (2,2) double-excitation seam's slow mode obeys the resummed coherence-ladder dispersion `λ² + 8γλ + 4J²q²` (CV-confirmed, the SE proof's own discriminator), and the handover (its darkness reaching the band-edge floor `⟨n_XY⟩ = 1`) sits at `Q·q = √3`, so `Q_h = √3/q_min → N·√3/(2π)`, slope `√3/(2π)`. The sibling of [the coherence horizon slope](PROOF_COHERENCE_HORIZON_SLOPE.md): same dispersion, but the `darkness = 1` condition (`Q·q = √3`) instead of the EP (`Q·q = 2`). Confirmed numerically N = 6..14 (Q_h·2π/N → √3 monotonically from above) and by the γ-constant-q CV discriminator (CV(8γ) = 0.012/0.019 at N = 12/10, beating the 4γ telegrapher ~10×). Resolves the open piece that kept `CoherenceHorizonClaim` + `SecondClockRegimeClaim` at Tier 1 candidate (the ring V-Effect seam), and corrects its sector label (the 2-excitation doublet, not half-filling).
**Date:** 2026-06-20
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Sibling of:** [the coherence horizon slope](PROOF_COHERENCE_HORIZON_SLOPE.md) (the single-excitation coherence horizon, slope 2/π chain / 1/π ring) and the parallel ceiling result [the ring g2=1 / commutant 2(N−2)/N, commit b191df3].

## What this is about

A ring of spins is watched, and the watching fades its coherences. At weak watching (large Q = J/γ) the longest-lived note is the single-excitation band edge, the smoothest standing wave. Turn the watching up (lower Q) and a different note outlasts it: a two-excitation pattern, two ripples instead of one, that the noise happens to spare longer. The handover Q_h is the tipping point between the two, the dephasing dose at which the two-ripple survivor stops being darker than the band edge and the band edge takes over.

This proof pins how that tipping point grows with the ring's length: a straight line of slope √3/(2π). The number is the band edge's own coherence-horizon slope (1/π for the ring) scaled by √3/2, and the √3 comes from one clean place: the two-ripple mode is governed by the same coherence-ladder dispersion as the single-excitation horizon, and asking when its darkness equals exactly one (the band-edge floor) gives `(Q·q)² = 3`.

## Abstract

On the cyclic XY ring under uniform Z-dephasing, below a crossover Q the longest-lived interior mode is the 2-excitation `(2,2)/(N−2,N−2)` doublet coherence (particle-hole partners, isospectral); above it, the single-excitation band edge wins. The handover Q_h is defined by the doublet survivor's darkness reaching the F50 / Absorption floor `⟨n_XY⟩ = 1` exactly (`Re(λ) = −2γ`).

The (2,2) slow mode is a two-particle density coupled to the **entire ladder of coherence ranges**, the same structure PROOF_COHERENCE_HORIZON_SLOPE resums for one particle: dilute on the ring, the two particles act as two independent single-particle density modes, so the long-wavelength dispersion is the same `λ² + 8γλ + 4J²q²`. This is confirmed directly by the SE proof's own discriminator: the overdamped slow eigenvalue, fed back through `λ² + 8γλ + 4J²q² = 0`, yields a γ-constant `q` (CV = 0.012 at N = 12, 0.019 at N = 10; the 4γ-truncated telegrapher scatters ~10× worse), with `q → q_min = 2π/N` (q_eff/q_min = 0.986 at N = 12).

The overdamped slow root gives the darkness `⟨n_XY⟩(Q) = 2 − √(4 − (Q·q)²)`. The handover `⟨n_XY⟩ = 1` is at `Q·q = √3`; the EP (oscillation onset, `⟨n_XY⟩ = 2`) of the **same** dispersion is at `Q·q = 2` and is the SE coherence horizon `Q* = 2/q_min = N/π`. Hence `Q_h = √3/q_min → N·√3/(2π)`, and `Q_h/Q* = √3/2`. The handover's high-Q limit (the EP darkness 2) is the (2,2) commutant `2(N−2)/N → 2`, the parallel ceiling result.

## Statement

The ring double-excitation handover Q_h(N) (the dephasing dose Q = J/γ at which the slowest 2-excitation
`(2,2)`-sector mode reaches the Absorption floor `Re(λ) = −2γ`, i.e. darkness `⟨n_XY⟩ = 1`, the band-edge
takeover) grows linearly with the ring length N, with exact asymptotic slope

    Q_h(N) → (√3/(2π))·N,   slope = √3/(2π) = 0.275664...

equivalently `Q_h = (√3/2)·Q*(N)` with `Q*(N) = N/π` the ring single-excitation coherence horizon. The
small-N values are transcendental (the short-ladder accident); √3/(2π) is the N → ∞ limit.

## Setup: the (2,2) sector and the handover

The XY ring conserves excitation number and Z-dephasing is diagonal in the coherence basis, so the
Liouvillian block-decomposes by `(ket #, bra #)` sectors. The longest-lived interior mode below Q_h lives
in the 2-excitation `(2,2)` sector (full-Liouvillian-verified at N = 6: the global survivor is the
particle-hole partner `(4,4)`, isospectral to `(2,2)`; the half-filling `(N/2,N/2)` sector gives a
different, non-matching Q_h, so the survivor is the **2-excitation doublet, not half-filling**, correcting
the earlier label). On the `(2,2)` block (dimension `C(N,2)²`, sector-projected, reaching N = 14):

    Re(λ) = −2γ·⟨n_XY⟩   (the Absorption Theorem),   handover:   ⟨n_XY⟩ = 1.

## The slow mode: the SE coherence-ladder dispersion (CV-confirmed)

The `(2,2)` slow mode is a two-particle density wave. In the long-wavelength (small-q) limit the two
particles are dilute and act as two independent single-particle density modes; each is the population
coupled to the entire ladder of coherence ranges of PROOF_COHERENCE_HORIZON_SLOPE (per-coherence dephasing
`4γ` since `n_XY = 2`, hopping `J`), so the collective density mode obeys the same resummed dispersion

    λ² + 8γλ + 4J²q² = 0.

This is confirmed by the SE proof's decisive discriminator (its cross-check 2a), which is finite-N-robust
where a full-curve fit is not: in the overdamped regime the slow `(2,2)` eigenvalue is real; fed back through
`q² = −(λ² + 8γλ)/(4J²)` it is γ-constant across a γ-sweep, CV = 0.0121 (N = 12) / 0.0194 (N = 10), while
the 4γ-truncated telegrapher `q² = −(λ² + 4γλ)/(2J²)` scatters ~10× worse (CV 0.16). The fitted q tracks the
ring `q_min = 2π/N` (q_eff/q_min = 0.978 / 0.986 at N = 10 / 12, approaching 1 from below with the O(1/N)
finite-size tail).

## The handover slope

The overdamped (real) slow root of `λ² + 8γλ + 4J²q²` is

    λ = −4γ + 2√(4γ² − J²q²)    (4γ² > J²q², i.e. Q·q < 2),

so the darkness, via the Absorption Theorem, is

    ⟨n_XY⟩(Q) = −Re(λ)/(2γ) = 2 − √(4 − (Q·q)²),    Q = J/γ.

This single function carries both clocks of the same dispersion:

- **the EP** (vanishing discriminant, oscillation onset, the double root at `λ = −4γ`): `Q·q = 2`,
  `⟨n_XY⟩ = 2`. This is the single-excitation coherence horizon `Q* = 2/q_min = N/π` (ring).
- **the handover** (darkness reaches the band-edge floor): `⟨n_XY⟩ = 1`, i.e.

      2 − √(4 − (Q·q)²) = 1   ⟹   √(4 − (Q·q)²) = 1   ⟹   (Q·q)² = 3   ⟹   **Q·q = √3.**

The handover is the **longest-wavelength** mode (smallest q, the last to reach the floor as Q rises),
`q_min = 2π/N` for the ring, so

    Q_h(N) = √3/q_min → N·√3/(2π),   **slope = √3/(2π) = (√3/2)·(1/π).**

The factor √3/2 against the SE coherence horizon is exactly the darkness-1 vs EP-2 condition on the one
dispersion: `√3 = √(Qq)²|_{D=1}` against `2 = (Qq)|_{EP}`.

## Cross-checks

1. **The endpoint, numerically (the analogue of the SE proof's `q_min·N → π`).** Q_h·2π/N → √3 from above,
   monotonically: 1.846 / 1.785 / 1.759 / 1.748 at N = 8 / 10 / 12 / 14, against √3 = 1.732, an O(1/N) tail
   exactly as the q → 0 derivation predicts (`ring_handover_qh.py` + `ring_handover_extend.py`).
2. **The dispersion, directly (the discriminator).** The γ-constant-q CV test above: CV(8γ) = 0.012 / 0.019
   (N = 12 / 10) vs CV(4γ) = 0.16, the resummed 8γ dispersion is the fixed mode, the truncated telegrapher
   refuted, exactly as in the SE proof (`ring_handover_dispersion_cv.py`).
3. **The high-Q limit = the parallel ceiling.** The EP darkness of the same dispersion is 2; the finite-N
   high-Q `(2,2)` darkness is `2(N−2)/N` (1.333 / 1.500 / 1.600 / 1.667 at N = 6 / 8 / 10 / 12, exact to
   4-5 digits), which equals the parallel-session commutant closed form `2(N−2)/N → 2` (b191df3): the
   handover's overdamped curve and the high-Q ceiling are the two ends of one dispersion. One object.
4. **The sector.** The full-Liouvillian survivor at N = 6 is the `(2,2)/(4,4)` doublet (verified, 100%
   weight on the PH partner), not half-filling; the `(2,2)`-block Q_h matches the full-L Q_h, the
   half-filling block does not.

## Scope

The derivation is leading-order (q → 0, N → ∞), so it fixes the asymptotic **slope** exactly (√3/(2π)); the
finite-N Q_h(N) for N ≥ 6 are transcendental (O(1/N) corrections, the same as the SE horizon). The
governing dispersion `λ² + 8γλ + 4J²q²` is established for the `(2,2)` mode by the CV discriminator (the SE
proof's load-bearing test); the explicit two-particle EOM resummation (the dilute-limit argument that the
`(2,2)` density mode reduces to two SE coherence ladders) is the physically-transparent remaining write-out,
not yet set down term-by-term as the SE proof does for one particle. Topology: the ring (the `q_min = 2π/N`
COM density wave); the chain handover is the co-located SE-EP (filling-degenerate), a different mechanism.

Axis: the dephasing Q-axis (XY, free fermions) only. The XXZ anisotropy Δ-axis handover (`Δ*`, the same
`darkness = 1` floor driven by the anisotropy instead of Q; ANALYTICAL_FORMULAS "the handover Δ",
`Δ*(N) → 1`) is a **distinct, interacting mechanism**: turning on the ZZ term moves the handover off
`Q·q = √3` (gate-checked, `simulations/ring_handover_zz_extension.py`: at N = 6, 8 the (2,2) handover
`Q_h·q_min` shifts at Δ = 0.5 and the crossing leaves the range by Δ = 1), so √3 is XY/free-fermion-specific
and the dispersion `λ² + 8γλ + 4J²q²` does not survive the ZZ coupling. The two handovers mirror only
through the shared `darkness = 1` floor (Q_h grows ~N, Δ* descends to 1); the Δ* closed form stays the open
Bethe-ansatz problem.

## Links

- Sibling proof: [the coherence horizon slope](PROOF_COHERENCE_HORIZON_SLOPE.md) (same dispersion, the EP
  condition `Q·q = 2 → N/π`; this proof is its `darkness = 1`, `Q·q = √3` sibling).
- Parallel ceiling: the ring `g2 = 1` (no high-Q ceiling) and the `(1,1)/(2,2)/(N/2,N/2)` commutant
  `2(N−2)/N → 2` (commit b191df3, `StructuralCeilingWitness` RingNode, `ring_ceiling_commutant_sweep.py`):
  the EP-darkness limit of this proof's curve.
- Typed homes (graduation-eligible once reviewed): `CoherenceHorizonClaim`, `SecondClockRegimeClaim`
  (the ring 2-excitation doublet V-Effect seam is their Tier 1 candidate reason; this derivation addresses it),
  arc `clock_hand_ladder`.
- Verifiers: `simulations/ring_handover_qh.py` (the sector discriminator + Q_h data, gate-first),
  `simulations/ring_handover_extend.py` (N = 14 + the `2(N−2)/N` cross-session convergence),
  `simulations/ring_handover_dispersion_cv.py` (the decisive γ-constant-q dispersion discriminator),
  `simulations/ring_handover_derivation_gate.py` (the leading-order curve, with the finite-N caveat that the
  full-curve collapse is the leading-order limit, not a finite-N identity; the CV test is the robust read).

## Review

Pending. Recommended lenses: physics-first-review (the dilute-limit dispersion assumption: is the two-particle
`(2,2)` density mode really two independent SE ladders, or does the hardcore constraint / relative coordinate
change the coefficients? the CV test says the 8γ form holds, but the EOM resummation is not written out) +
mathematical-review (the darkness-1 algebra and the `Q·q = √3` condition; the q_eff → q_min finite-size
approach). The √3/(2π) endpoint is numerically solid (cross-check 1); the load-bearing claim to attack is
that the governing dispersion is the SE one (cross-check 2).
