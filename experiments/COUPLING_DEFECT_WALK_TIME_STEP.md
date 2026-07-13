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

**What stayed open here, closed below.** Whether the effective A settles on a clean fraction of γ_φ, and the exact coherent-versus-incoherent split of the front peak as n and γ vary, were left unpinned at this point (the K-collapse is only approximate and A drifts with γ); the second follow-up below answers both from the exact model. The robust content of this section stands unchanged: the coherent caustic pays the full 4γ, the front's extra survival is incoherent refill, and the naive full-rate estimate for the total front is falsified; and, for Item 1, the intrinsic single-q second order is the small +δ²/(4J) with no odd asymmetry while the large measured residual is the relative-threshold method, not physics.

## Follow-up (2026-07-13, later the same day): the survival exponent, a first closure (corrected in part below)

The unpinned A question is answered by deriving the front survival from the exact model; the verification is `simulations/cone_front_survival_asymptote.py` (output in `simulations/results/cone_defect_arrival/front_survival_asymptote.txt`). The derivation passed independent mathematical and physics-first referee rounds; each referee recomputed from scratch and each caught one real correction, both folded in below.

**The exact renewal representation.** Split the generator as L = (−i[h, ·] − Γ) + Γ·diag with Γ = γ_φ = 4γ. The no-refill part propagates as e^{−Γτ}·e^{−ihτ}ρe^{+ihτ} (a uniform scalar decay riding the clean evolution), and the Dyson ladder read at the diagonal gives the exact closed equation

    P_n(t) = e^{−Γt}·S_n(t),   S_n(t) = |G_{n0}(t)|² + Γ ∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S_m(s),

with G_{nm}(τ) = (−i)^{|n−m|}·J_{|n−m|}(2Jτ) the clean propagator. The reading: the population is a ballistic front that dephasing re-diagonalizes at rate Γ and that re-propagates cleanly between events. Every refill order carries the same universal prefactor e^{−Γt} (segment times sum to t), the j = 0 term is exactly the coherent front |⟨aₙ⟩|² of Item 2, and everything above it is the incoherent halo. In momentum-Laplace space the equation closes to Ŝ(p, z) = 1/(√(z² + a²) − Γ), a = 4J·sin(p/2); probability is conserved at p = 0 and γ = 0 returns the clean Bessel front. The representation is verified against direct RK4 to 1.6·10⁻⁶.

**The answer: no sub-γ_φ constant exists; the one closed form is the endpoint A_∞ = γ_φ.** Write −ln g = Γt* − ln B with B ≥ 1 the incoherent boost. The coherent piece pays A_coh = 8Jt*/n with 2Jt* = n + 0.809·n^{1/3} (the first Airy maximum sits just inside the cone), so A_coh = 4·(1 + 0.809·n^{−2/3}) → 4 from above (4.34, 4.24, 4.15 at n = 30, 50, 100). The boost piece vanishes like n^{−1/3}: the single-refill integral I₁ saturates to the n-independent constant 0.2767 while the caustic peak decays as 0.4553·n^{−2/3}, giving (ln B)/K ≈ (8·I₁/0.4553)·n^{−1/3} = 4.86·n^{−1/3} → 0 (a single-refill, small-γ asymptotic; the verifier pins it at γ = 0.002, the first-order regime). Hence

    A_eff(n) = 4 − O(n^{−1/3}),   A_∞ = γ_φ = 4.

(Corrected the same evening, third follow-up below: this law is the pre-asymptotic window n ≪ n* ≈ 6·(J/γ)^{3/2}. At fixed γ the true ceiling is A_∞(γ) = 4 − φ(2J)/γ < 4; the endpoint 4 survives only as the γ → 0 limit of that ceiling.)

The measured 2.55–3.05 is the pre-asymptotic surface of this slow climb (4 − A fits 3.13·n^{−0.26} in the window; the apparent exponent is shallowed from 1/3 by the n^{−2/3} sub-leading), reaching only ≈ 3.2 at n = 200 (extrapolated from the fit; the verifier computes through n = 55). There is no 8/3 and no fixed 2.8; A ≈ 2.8 is simply where the climb sits at n ≈ 20–55. The climb is already visible in the previous section's own numbers (at γ = 0.01, A rises 2.77 → 3.05 as n goes 15 → 55).

**The two referee corrections.** First, the caustic power-counting fixes the exponent but not the constant: the joint-caustic contribution to I₁ is only ≈ 30% (0.083 of 0.277); the bulk interior of the cone supplies the rest. Second, the γ-direction of the drift is mostly the measurement again: at fixed n = 40 the physical same-instant survival (numerator and denominator both read at the clean caustic time) drifts only −0.14 across γ = 0.002–0.08, while the ratio-of-maxima method drifts several times more (−0.54 by γ = 0.05, −0.91 by γ = 0.08; at this n the windowed peak collides with the diffusive bump from γ = 0.05 on, the EDGE class of the previous section), because the damped front peak slides in time against the fixed t_arr = n/(2J) inside K. Once that peak-time shift is removed, A is nearly a function of n alone. So the honest two-variable statement collapses to one variable: A_eff ≈ 4 − O(n^{−1/3}), physically almost γ-independent; the measured γ-spread of the K-collapse is dominated by the max-ratio method. The same lesson as Item 1, for the third time in this experiment: the residual lived in the method.

**What stayed open here, closed and corrected below.** The two smallest items this section named, the closed form of the n^{−1/3} coefficient and the all-order fate of A_eff at fixed γ, are both answered in the third follow-up: the coefficient closes exactly, and the fixed-γ conjecture is refuted (the pinned "endpoint A_∞ = γ_φ" of this section holds only as the γ → 0 limit). What survives of this section unconditionally: the exact renewal representation, the coherent front paying the full 4γ, and the re-attribution of the γ-drift to the measurement.

## Follow-up (2026-07-13, evening): the coefficient closes, and the ceiling is not 4

The two items above are answered from the same exact renewal object, and the second answer overturns the sentence the previous section pinned. Both derivations passed independent mathematical and physics-first referee rounds; both referees recomputed from scratch with their own methods (an independent Lindblad ODE reproducing the tables to the digit; an independent local-slope discriminator at γ = 0.1), and both corrections they carried are folded in. The verification sections [7] and [8] of `simulations/cone_front_survival_asymptote.py` pin the numbers.

**Item A: the n^{−1/3} coefficient in closed form.** The single-refill constant is exactly

    I₁ = 1/12 + (1/4)·∫₀^{2c} Ai(−w) dw = 0.27694424…,   2c = 2^{2/3}·α,

with α = 1.0187929716… the magnitude of the first zero of Ai′ (the front-peak location). The route is the exact momentum-Laplace form of the refill integral (the Graf identity plus the convolution ∫₀ᵗ J₀(a(t−s))·J₀(as) ds = sin(at)/a collapse the (s, m) double sum to one θ-integral), which splits at the front into a Dirichlet endpoint term, exactly 1/8, the only rational piece of that split, and an Airy caustic term, I₁ − 1/8 = 0.15194. The previous section's 0.2767 was the pre-asymptotic reading (the exact integral, extrapolated from n ≤ 6400, gives 0.276949, agreeing with the closed form to 4.4·10⁻⁶). Sharpened alongside: the caustic peak constant is 2^{2/3}·Ai(−α)² = 0.45547, so the pre-asymptotic coefficient is 8·I₁/0.45547 = 4.864.

**Item B: the fixed-γ conjecture is refuted; the ceiling sits below 4.** At every fixed γ > 0,

    A_∞(γ) = 4 − φ(2J)/γ,   φ(2J) = √(Γ(Γ + 4J)) − 4J·arcsinh√(Γ/4J),   Γ = 4γ,

with the small-γ form A_∞ = 4 − (8/3)·√(γ/J) + O(γ^{3/2}): 3.88 at γ = 0.002, 3.41 at 0.05, 3.18 at 0.10. The mechanism: the single-refill linearization behind the n^{−1/3} law breaks at large n (the boost B − 1 grows like n^{2/3}), and the full refill ladder resums to a large-deviation pole, μ(θ) = √(Γ² + 16J²·sinh²(θ/2)), whose Legendre transform at the front gives ln S_n(t*) = (φ/(2J))·n + o(n): the incoherent halo boost is exponential in n, not a vanishing correction, so the front's asymptotic decay rate stays strictly below the bare γ_φ forever. The n^{−1/3} climb toward 4 is real but pre-asymptotic, confined to n ≪ n* ≈ 6·(J/γ)^{3/2} (n* ≈ 540 at γ = 0.05, ≈ 67000 at γ = 0.002); every window of this experiment sits deep inside it, which is why the climb read as "toward 4". The discriminator, measured from below at γ = 0.10: the local slope d(ln S_n)/dn climbs monotonically 0.032 → 0.037 across n = 50–120, toward φ/2 = 0.0410 and past the single-refill ceiling (~0.031, which moreover predicts the wrong direction). Causality holds: v = 2J is an interior point of the tilted ensemble (the bare Bessel propagator carries exponentially suppressed super-cone tails), so the refund is legitimate physics, order-n collapse-and-re-propagate events riding the cone.

**What this does to the story.** The headline sharpens rather than flips: the coherent front still pays the full γ_φ exactly, and the halo still rescues the total front; what changes is that the rescue never fades. The watched front keeps, asymptotically and forever, the fraction φ(2J)/γ_φ of its rate bill as incoherent refund; "the dose the front pays is amplitude, not schedule" gains the clause "and never quite the full amplitude bill either". There is still no clean fraction of γ_φ: the ceiling is a transcendental function of γ/J, and the one clean number, 4, is only its γ → 0 endpoint.

**What stayed open here, closed below.** The o(n) prefactor of the large-deviation law and the windowed-max counterpart are both pinned in the fourth follow-up; what remains after it is only the validity edge at γ ~ J, where the ballistic caustic observable itself degenerates (all quoted values live at γ ≪ J), and that is a fence, not a question this experiment's observable can ask.

## Follow-up (2026-07-13, night): the o(n) closes, and the third reading

The last two open items, the prefactor of the large-deviation law and the windowed-max counterpart, are answered from the same saddle geometry; section [9] of `simulations/cone_front_survival_asymptote.py` pins the numbers. Both derivations passed independent mathematical and physics-first referee rounds again; the math referee re-derived every constant symbolically and added a third method (direct steepest-descent quadrature of the pole contribution alone), the physics referee cross-checked with a structurally independent direct Lindblad integration and the exact Haken-Strobl walk.

**The prefactor.** The sub-exponential structure of the front ladder is a plain Gaussian saddle, power exactly −1/2 (no Airy degeneracy: v = 2J is a regular interior velocity of the tilted ensemble), with the constant in closed form:

    S_n(t*₀) = C(γ)·n^{−1/2}·e^{(φ/2J)·n}·(1 + O(1/n)),   C(γ) = (2π)^{−1/2}·(γ/(γ+J))^{1/4},

carried by the clean saddle-curvature identity μ″(θ*) = 2J·√(Γ/(Γ+4J)). The renewal-measured constant converges to ln C within 0.012 by n = 160 at γ = 0.1 on the verifier's grid (the referee's finer-grid run reaches 0.002). The γ → 0 handover is continuous: C ~ γ^{1/4} and φ ~ γ^{3/2} vanish together (both trace to the refill vertex; no dephasing, no ladder, no halo), and for fixed n the reading reverts to the coherent caustic below the crossover n*.

**The approach law.** Against the coherent peak reference 0.45547·n^{−2/3}, the −1/2 numerator power forces the logarithm's coefficient:

    A_eff(n, γ) = A_∞(γ) + (2J/(γn))·[−(1/6)·ln n + ln(0.45547/C(γ))] + o(1/n),

an O(ln n / n) approach that reproduces the exact renewal A_eff to 0.0014–0.0064 across n = 60–160 (the verifier's window; the derivation's runs extend the agreement to n = 240). Read against this reference, A_eff crosses its own ceiling (predicted n_x = 80.8 at γ = 0.1; the verifier measures the sign change between n = 81 and 82, the referee's independent Lindblad solver between n = 80 and 100), dips shallowly below, and returns along the analytic tail (the return is the law's asymptote, too shallow to resolve as an independent measurement). The crossing itself is a reading, not physics: it exists because the numerator is read at the caustic time while the reference is the coherent peak; a same-instant reference removes it and the approach is monotone from below. Convention-robust are the ceiling A_∞, the −1/6, and the −1/2 power; the halo-refilled front decays as n^{−1/2}, genuinely slower than the coherent n^{−2/3}, which is the physical content under the crossing.

**The third reading: the peak-tracking exponent is exactly zero.** The windowed-max counterpart degenerates, and it does so exactly: the stationarity ∂_t ln S_n = Γ forces θ* = 0, v → 0, so the global maximum of P_n(t) at fixed n is not a front feature at all but the late diffusive plateau, the exact Haken-Strobl walk P_n(t) = e^{−2Dt}·I_n(2Dt) with D = 2J²/Γ: peak at t = n²/(2D), height Pₙ·n → e^{−1/2}/√(2π) = 0.242, algebraic in n with no exponential at all. So the survival exponent comes in three readings: the fixed-time (caustic) reading pays 4 − φ/γ; the global-max reading pays 0 (the halo refills the peak height to algebraic order); a finite arrival window interpolates and runs to its late edge. This is the retroactive explanation of the EDGE class: the excluded points (γ = 0.08 at n ≥ 30 in the survival fit, and the γ ≥ 0.05 rows of the artifact split) are exactly where the diffusive plateau drifts into the arrival window and drags the reading off the ceiling branch toward zero. The committed A ≈ 2.8 sits where it should: a narrow front window at the ballistic time, on the ceiling branch, pre-asymptotic.

**What stayed open, the fence, docked below.** The γ ~ J edge resolves by survey plus one honest instantiation in the last follow-up; with it the arc's ladder is closed end to end: the step −δ/(2J), the threshold artifact, the halo mechanism, the ceiling A_∞(γ) = 4 − φ(2J)/γ, the coefficient I₁ in Airy form, the prefactor C(γ), the approach law, the trichotomy of readings, and the fence.

## Follow-up (2026-07-13, the fence): docked, not derived anew

The fence sentence "the validity edge at γ ~ J, where the ballistic caustic observable degenerates" was imprecise in a useful way, and the repo already owned everything needed to sharpen it; the verification is section [10] of the same verifier, the derivation passed an adversarial referee round. Three findings:

1. **There is no n-free front EP, and the tempting one is a variable trap.** The owned single-excitation dispersion (λ² + 8γλ + 4J²q², √-EP at γ*(q) = Jq/2; `docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md`, the CoherenceHorizonClaim) is indexed by the coherence's centre-of-mass wavevector q, not by band momentum. The ballistic front rides band momenta k ≈ ±π/2, but its coherence content is q → 0 (E(π/2 + q/2) − E(π/2 − q/2) = 2Jq at leading order, the group velocity 2J), where γ*(q) → 0: no single freezing point. Substituting the band momentum instead (γ* = πJ/4, Q = 4/π ≈ 1.27, temptingly between the hardware brackets) is wrong twice over, wrong variable and outside the dispersion's q → 0 validity. Consistently, the fourth follow-up's tilted pole μ(θ) is analytic in γ: the front has no EP in γ at all.

2. **At this experiment's distances the fence is a dose, an order of magnitude below J.** The ballistic arrival reading dies where the survival law meets its own threshold, g(K) < θ: with θ = 0.2 that is K_deg ≈ 0.68, 0.63, 0.60 at n = 20, 30, 40 (section [10]), i.e. γ_deg ≈ (1.2–1.4)·J/n and Q_deg ≈ 0.8·n. The direction of the θ-dependence is K_deg ~ ln(1/θ)/A_eff, but it is not a single-A identity (ln g is concave in K; the implied A drifts ≈ 2.1 to 2.5 across θ = 0.1 to 0.4). So "γ ~ J" was never the working edge at n = 20–60; the working edge is the dose line.

3. **The n-free content the fence does have is already owned and measured.** The coherent-transport window is L* = Q/2 (the ENAQT crossover MirrorWorld's `spread` mode validates), and the n → 1 corner of the dose line is the exact small-N EP ladder, Q*(2) = 1 and Q*(3) = √2, which are the two hardware-confirmed handovers (critical damping at Q = 1 for N = 2; the single-excitation-walk revival handover at Q ≈ 1.5 for N = 3, √2 vs 1.5 bracketed, not matched; ConfirmationsRegistry, ibm_kingston). Extrapolating Q_deg ≈ 0.8·n down to n = 1–2 straddles that bracket. One family read at three lengths, all in the repo before this follow-up, though not one continuous law (the observables differ): the small-N EP corner, the dose line at this experiment's n, and the band-edge slope Q*(N) = 2N/π at the far end.

The caustic fine structure (the Airy oscillations, coherence width q ~ n^{−1/3}) overdamps at the larger scale γ ~ (J/2)·n^{−1/3}; asymptotic, not measured, and above the arrival threshold at every n ≥ 4. With this the experiment carries no open item; what lies beyond its observable belongs to the coherence-horizon arc.

## Links

- The frame this computes: `hypotheses/PERSPECTIVAL_TIME_FIELD.md` (α_i, §1.3 numbers, §3 mechanism); the distance pin in MirrorWorld (`compute/MirrorWorld.Tests`, SpookyActionTests: the price −2γk never sees distance; this experiment is where distance DOES live, in the way walked).
- The adopted C# witness (2026-07-13): `compute/MirrorWorld/WalkTime.cs` reads the step live over the Cone engine (`Cone.SetBond` is the defect knob; run mode `walk N [delta]`); `compute/MirrorWorld.Tests/WalkTimeTests.cs` pins the committed plateau numbers at both signs of δ and the γ = 0.05 near-field survival from below.
- The typed home of the exact object (2026-07-13): `docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md` (the renewal representation and its Green's function, the sector's exact solution) = F126 in `docs/ANALYTICAL_FORMULAS.md`, typed as `DephasingFrontRenewalClaim` with the live witness `inspect --root renewal`; the survival readings of the follow-ups are its Tier 1 candidate corollaries.
- The play sequel (2026-07-14): `experiments/FRONT_PEDIGREE.md` resolves the F126 ladder by catch count, the observable no stepping engine can show: at K = 1.25 only 29% of the surviving front was never caught, survivors carry a quarter of the free walker's catches, and the last catch is spread over the whole trip (the bulk-dominance of I₁ made visible).
- Front speed and dispersion: F2b in `docs/ANALYTICAL_FORMULAS.md` (the single-excitation band propagated here); `docs/proofs/derivations/D10_W1_DISPERSION.md` (the weight-1 Liouvillian dispersion, the same physics one sector up, background).
- The dissipative defect law the far field would obey if the timing ever went diffusive: `SurvivorDiffusionGradientClaim` (Tier-1 Candidate, `compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientClaim.cs`).
- The earlier dissipative arrival table (perturbation-response arrival, a different observable): `simulations/results/n7_coupling_defect_overlay/run_log.txt`.
