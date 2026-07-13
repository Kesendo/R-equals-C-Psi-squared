# The Walk-Time Step: a Coupling Defect Reads as a Local Metric Perturbation in the Ballistic Channel

**Date**: 2026-07-12
**Status**: Computationally verified (pre-registered predictions, two design reviews before running; N=7, 20, 60, 120)
**Script**: `simulations/cone_defect_arrival.py`
**Question**: does a single-bond coupling defect δJ read quantitatively as a local change of walk-time, the way "Distanz ist t" would demand? And is that the same object as the PTF's per-site α_i response to the same defect?

## The two answers up front

1. **Yes, in the ballistic channel, and exactly at first order.** The arrival-time delay profile of a single-excitation front crossing a defect bond J′ = J(1+δ) is a step: zero upstream, constant −δ/(2J) downstream (first order; the broadband seed undershoots, see P2). That constant is the naive metric answer, "this one bond now costs 1/(2J′) instead of 1/(2J) of walk-time", and the scattering derivation shows the coincidence is exact at first order in δ. The front pays the bond's local walk-time and nothing else.

2. **No, it is not the α_i object.** The PTF's α_i response to the same class of defect is smooth, sign-mixed, and nonlocal (a 15% response four sites from the bond). The defect writes into two separate channels: a local, additive time-metric shift carried by the front (this experiment), and a global eigenvector rotation read by dissipative site observables (the PTF). The channels do not reduce to one another; the locality scan below makes the contrast quantitative.

## Setup

Single-excitation sector of the XY chain (Pauli-J convention, H = J·Σ ½(X_bX_{b+1} + Y_bY_{b+1})): hopping matrix h with amplitude J = 1, band E(q) = 2J cos q, front speed |v_g|max = 2J (a corollary of F2b, Tier 1). One defect bond J′ = J(1+δ). γ = 0 runs use exact eigenbasis evolution of the pure state; γ > 0 runs use RK4 on the N×N single-excitation ρ with the exact sector dissipator ρ̇ = −i[h,ρ] − 4γ(ρ − diag ρ) (every |1_a⟩⟨1_b| coherence has joint popcount distance 2, so the rate is −2γk = −4γ; the sector is closed under Z-dephasing).

Arrival time: t_arr,i = first t with P_i(t) ≥ θ·max_t P_i(t)|clean, θ = 0.2, crossing linearly interpolated on a dt = 0.001 grid (dt = 0.002 for the γ > 0 runs). The threshold is relative because the coherent front peak decays along the chain (≈ n^(−2/3), the Airy caustic, derived in the 2026-07-13 follow-up); an absolute threshold collides with that decay (see "Reconciliation" below, where exactly this produced a wrong scout verdict).

The design was pre-registered and passed two adversarial design reviews (mathematical + physics-first) before the production runs; the physics review falsified one original prediction (P4, dephasing robustness as first stated) and forced the relative threshold, the quasi-monochromatic run, and the locality scan into the design. What follows are the revised, pre-registered tests and their outcomes.

## The derivation (first, on paper)

Plane-wave scattering across the defect bond (u = J′/J, incident e^{iqn}):

    τ(q) = −2iu sin q / (e^{−iq} − u² e^{iq}),   |τ|² = 1 − O(δ²),   φ(q) := arg τ ≈ δ·cot q.

The transmitted packet is shifted by Δn = −φ′(q) = +δ/sin²q; the right-moving front is carried by q = −π/2 (v_g = +2J), and the arrival-time change is Δt = −Δn/|v_g|. At the front:

    Δt = −δ/(2J) = 1/(2J′) − 1/(2J) + O(δ²).

Two structural facts force the clean result: |τ|² = 1 − O(δ²) (no first-order reflection, so the front loses no amplitude to the defect at this order), and φ″(±π/2) = 0 (the transmission phase is locally linear at band center, so the Airy caustic at the cone edge translates rigidly instead of deforming). At second order the scattering value −(u²−1)/(2J(u²+1)) ≈ −(δ − δ²/2)/(2J) and the metric walk-time (−δ + δ²)/(2J) part ways; the coincidence is a first-order statement and is claimed as such.

## Results against the pre-registered predictions

**P1, the step law (N=60, γ=0, seed site 0, defect bond (29,30)).** Confirmed. Upstream side ≤ 10⁻⁵ in |Δt|; downstream plateau flat to σ ≤ 0.002 across 16 far sites:

| δ | plateau (sites 40–55) | first order −δ/(2J) |
|---|---|---|
| −0.10 | +0.0633 ± 0.0017 | +0.0500 |
| −0.05 | +0.0287 ± 0.0006 | +0.0250 |
| +0.05 | −0.0235 ± 0.0001 | −0.0250 |
| +0.10 | −0.0424 ± 0.0004 | −0.0500 |

A stronger bond advances the front, a weaker bond delays it; the ± asymmetry is second-order plus finite-packet structure; the single-q second-order value and the 1/sin³q dispersion weighting do not jointly account for both signs, so the residual mechanism is left unpinned here (see P2; it is pinned in the 2026-07-13 follow-up as the relative-threshold amplitude artifact, with the intrinsic single-q second order the modest +δ²/(4J)).

**P2, the magnitude (quasi-monochromatic packet, q = −π/2, σ_q = 0.1π, N=120).** Confirmed within the pre-registered ±5%: fitted slope d(Δt)/dδ = −0.518 vs −1/(2J) = −0.500 (3.7%). The broadband single-site seed gives −0.528 on the symmetric fit; the per-δ second-order residuals above are documented as measured and not chased further here.

**P3′, the locality scan (defect bond b ∈ {10, 20, 30, 40}, δ = +0.10).** Confirmed. The step edge sits at the bond (detected at b or b−1, within one site, for all four positions); the plateau height is b-independent to ~6% (−0.0450, −0.0439, −0.0430, −0.0422; the residual is a monotone drift with b, a dispersion/window effect of the measurement segment sliding downstream, not scatter); the upstream side stays ~0.0012 or below. This is the quantitative half of the channel contrast: the ballistic response is a pure downstream step whose edge tracks the bond, while the α_i response to a defect at bond (0,1) at N=7 is spread over all seven sites with mixed signs (α = 1.095, 1.182, 1.051, 0.991, 0.845, 0.923, 0.997; `hypotheses/PERSPECTIVAL_TIME_FIELD.md` §1.3). Same perturbation class, two structurally different response shapes. (Part of the difference is inherited from the probes themselves, a localized packet vs a delocalized mode's purity; the contrast is between the channels as measured, and is fenced as such, not promoted to an intrinsic-split theorem.)

**P4′, dephasing (γ = 0.05).** The near-field window test passed far better than the pre-registered factor 2: at N=20, defect (4,5), sites 6–9 read −0.0488 to −0.0491 vs −0.0500 (2%). The far field produced the one result that outran the design: at N=60, all 16 far sites (40–55) cross and carry the step at −0.0492 ± 0.0001. The pre-registration, following the physics review's scout, expected the far field to be diffusive and made no step claim there.

**Why the far field still carries the step (the surprise, checked before believed).** The step is threshold-stable for θ ∈ {0.05, 0.1, 0.2, 0.4} (−0.047 to −0.050) and only degrades at θ ≥ 0.6, where the crossing sits on the shallow late maximum. The mechanism is that at Q = J/γ = 20 the front's *timing* is still ballistic at these distances: site 50 crosses clean at t = 25.6 ≈ 50/(2J), not at a diffusive t ~ n²·4γ/(4J²) ≈ 125. Dephasing eats the front's amplitude (max P at site 50 is 0.0069, far below the γ=0 value), but the arrival structure, and with it the metric step, survives in the timing. The dose the front pays is amplitude, not schedule.

**Reconciliation of the design-review scout.** The physics review's scout run saw the far-field step "destroyed" (+0.068, 5/16 sites crossing) because it used an absolute threshold θ = 0.01, which exceeds the damped front peak at far sites (0.0069 at site 50); re-running that observable reproduces the artifact class (+0.095, 2/16 sites here; the scout's exact figures differ, almost certainly a different time grid or window, but both show the same signature: most sites never cross, and the few crossings sit near their shallow peaks and scatter positive). The same review's own Finding 4 named the collision of absolute thresholds with front decay; the relative threshold it demanded is what resolves it. The artifact observable is reproduced in the script's follow-up probe, so the disagreement is closed, not filed away.

**Controls.** dt-halving: max change 10⁻⁷. Lattice reflection (seed AND defect mirrored together, i ↦ N−1−i): profiles mirror to 1.4·10⁻¹⁴. θ-independence: above. The upstream O(δ) reflected wave (r ≈ δ(−1 + i cot q), |r|² = O(δ²)) is visible as a faint standing pattern and is not read as signal.

**The N=7 anchor (descriptive, no test).** Near-middle defect bond (2,3), seed 0: Δt = [seed, 0.000, +0.007, −0.045, −0.047, −0.055, −0.048] at γ=0 and nearly identically at γ=0.05. Even at N=7 the shape is the step: upstream flat, downstream plateau near −0.05. Printed beside the α row above for the eye; the two observables are incommensurate (localized-packet arrival vs delocalized-state purity rescaling), which is exactly why the locality scan, not this table, carries the P3′ contrast.

## What this means, and what it does not

**"Distanz ist t" gains its computation.** Distance in this world is walk-time converted by J (L = 2Jt, the F2b group velocity); this experiment shows the conversion is *locally additive*: deform one bond and precisely the walk-time of that one bond changes, for every site beyond it, by the same amount, at first order. A bond defect is a local metric perturbation in the literal sense that the total arrival time is the sum of per-bond walk-times and the defect edits one summand. The step's site-independence downstream is the additivity; its absence upstream is the locality; the θ- and γ-robustness say the front's schedule, not its amplitude, is the carrier.

**The channel split stands.** The same δJ writes a local step into the ballistic value channel and a global, smooth rotation into the dissipative vector channel (the PTF's α_i, which arise from eigenvector mixing while the slow eigenvalues are first-order protected; `hypotheses/PERSPECTIVAL_TIME_FIELD.md` §3). This is the value/vector decoupling of the felt-time arc, now with the value side measured as a metric step rather than inferred.

**No tier promotion is implied.** The PTF stays Tier 2 (EQ-014's closure verdict is untouched). The step law itself is ordinary tight-binding scattering, claimed with its forcing facts named (first-order unitarity, φ″(±π/2) = 0); the repo-new content is the measured coincidence with the metric reading, the locality contrast against α_i at the same perturbation, and the γ-robustness of the timing. The second-order residuals of the finite packet and the exact bookkeeping of how much front amplitude survives to a given distance at given γ were left open here; both are pinned in the 2026-07-13 follow-up below.

## Follow-up (2026-07-13): the two open residuals, pinned

The step law above claimed first order only and left two residuals open (the P1 note, and the closing "What this means"). Both are now pinned with the same single-excitation engine; the verification is `simulations/cone_walk_time_residuals.py` (output in `simulations/results/cone_defect_arrival/walk_time_residuals.txt`). The object both items needed is the real-space propagator: one excitation released from a single site spreads as aₙ(t) = (−i)ⁿ Jₙ(2Jt), so Pₙ(t) = Jₙ(2Jt)². The light-cone edge is the Bessel caustic, and near it Jₙ(n + ξ·n^(1/3)) ≈ (2/n)^(1/3)·Ai(−2^(1/3)·ξ); the Airy maximum 0.5357 fixes the front peak at max_t Pₙ ≈ 0.46·n^(−2/3) for interior injection, and ≈ 1.8·n^(−2/3) at the chain end, where the image sum Jₙ + Jₙ₊₂ adds in phase. This is the source of the n^(−2/3) decay used to justify the relative threshold in the Setup, with one caveat: the boundary coefficient converges slowly, so the measured Pₙ·n^(2/3) still drifts (1.32 → 1.58 across sites 15 to 55) and the local log-log slope reads shallower than 2/3 in this range.

**Item 1: the second-order residual is the relative threshold, not new dynamics.** The exact single-bond scattering gives the front delay in closed form,

    Δt_sq(δ) = −(u² − 1) / (2J·(u² + 1)),   u = J′/J = 1 + δ,

with Taylor series −δ/(2J) + δ²/(4J) − 0·δ³ + O(δ⁴). The intrinsic second order is the modest +δ²/(4J) (coefficient +0.25 at J = 1), and there is exactly no δ³ term (the even +δ²/(4J) still makes the magnitudes |Δt(+δ)| and |Δt(−δ)| differ at O(δ²), but the single-q front carries no odd asymmetry). A direct numeric Wigner delay reproduces the closed form to six digits.

The packet nonetheless measured a large even residual (fitted coefficient +1.18) and a strong odd asymmetry (−1.52). Both are the relative threshold itself. The defect reflects the fraction |r|² ≈ δ², so the transmitted front is dimmed by the same δ² (measured: the plateau peak ratio is 1 − 1.08·δ²). A front that is δ²-dimmer reaches a threshold fixed at θ·max_t Pᵢ|clean δ²-later, an even O(δ²) delay with nothing to do with walk-time. Referencing the threshold to each run's own peak removes it and collapses the residual to single-q order: the even coefficient falls +1.18 → +0.135 (comparable to the intrinsic +0.25, not equal to it), the asymmetry −1.52 → +0.45. A smaller finite-packet remainder stays, but the amplitude artifact is ≈ 89% of the even residual; the "residual mechanism left unpinned" in P1 is the measurement method. The first-order metric coincidence is untouched, and the honest second-order statement is that the intrinsic single-q term is the small +δ²/(4J) and the rest was the threshold.

**Item 2: the coherent front pays the full dose; the survival boost is incoherent.** The front's schedule stays ballistic while its amplitude decays (P4′); the bookkeeping is now measured. Let g(n,γ) = max_t Pₙ(γ) / max_t Pₙ(0) be the survival of the first-arrival front peak. It falls off exponentially in the dose K = γ·t_arr = γ·n/(2J):

    g ≈ exp(−A·K),   A ≈ 2.8   (power exponent 0.98; the window-edge-contaminated points, γ = 0.08 at n ≥ 30 where the diffusive peak has drifted past the arrival window, are excluded from the fit).

The collapse onto K alone is only approximate: the effective slope (−ln g)/K spreads 2.55 to 3.05 (about ±9%), drifting lower as γ grows.

The natural benchmark is the coherent front. The noise-averaged amplitude ⟨aₙ⟩ damps uniformly at γ_φ/2 = 2γ (each site's dephasing damps the amplitude there at half the coherence rate γ_φ = 4γ, Setup), so the coherent front population survives as g_coh = |⟨aₙ⟩|² / Pₙ^coh = e^(−4γ·t_arr) = e^(−4K), the A = 4 line. The measured front beats it: at n = 50, γ = 0.05 (K = 1.25) the coherent front gives g_coh ≈ 0.005, about 0.7× the naive point value 0.007 (the caustic peaks a little later than t_arr, at t ≈ n/(2J) + O(n^(1/3)), so it pays the dose slightly longer), while the full front survives at 0.030, six times better, decaying at an effective rate 2.8·γ = 0.7·γ_φ.

The gap is not the caustic decaying more slowly: the coherent caustic pays the full 4γ. Because ρ is the noise average, Pₙ = ⟨|aₙ|²⟩ = |⟨aₙ⟩|² + Var(aₙ); the first term is the coherent front (paying 4γ), the second is incoherent population that dephasing converts locally from the coherent front, a halo co-moving with the leading edge (diffusion is sub-ballistic and cannot overtake the front, so this is local conversion, not transport to it). The front survives on that incoherent halo, and the reduced effective A ≈ 2.8 is a property of the total arrival-window population, not of a slower caustic. This is the quantitative form of "the dose the front pays is amplitude, not schedule": the schedule (the ballistic timing) is the coherent front's, but the surviving amplitude is increasingly the incoherent halo's.

**What stays open, smaller.** Whether the effective A settles on a clean fraction of γ_φ, and the exact coherent-versus-incoherent split of the front peak as n and γ vary, are not pinned (the K-collapse is only approximate and A drifts with γ). The robust content is that the coherent caustic pays the full 4γ, the front's extra survival is incoherent refill, and the naive full-rate estimate for the total front is falsified; and, for Item 1, that the intrinsic single-q second order is the small +δ²/(4J) with no odd asymmetry while the large measured residual is the relative-threshold method, not physics.

## Links

- The frame this computes: `hypotheses/PERSPECTIVAL_TIME_FIELD.md` (α_i, §1.3 numbers, §3 mechanism); the distance pin in MirrorWorld (`compute/MirrorWorld.Tests`, SpookyActionTests: the price −2γk never sees distance; this experiment is where distance DOES live, in the way walked).
- Front speed and dispersion: F2b in `docs/ANALYTICAL_FORMULAS.md` (the single-excitation band propagated here); `docs/proofs/derivations/D10_W1_DISPERSION.md` (the weight-1 Liouvillian dispersion, the same physics one sector up, background).
- The dissipative defect law the far field would obey if the timing ever went diffusive: `SurvivorDiffusionGradientClaim` (Tier-1 Candidate, `compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs`).
- The earlier dissipative arrival table (perturbation-response arrival, a different observable): `simulations/results/n7_coupling_defect_overlay/run_log.txt`.
