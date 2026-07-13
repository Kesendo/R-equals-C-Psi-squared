# PROOF: the Dephasing Front Renewal Representation (the exact solution of the watched walk)

**Status:** Tier 1 (the representation and its momentum-Laplace closed form are exact; verified against direct RK4 to 1.6·10⁻⁶ and independently re-derived by two adversarial referee rounds plus an independent Lindblad ODE reproducing the tables to the digit). The asymptotic readings built on it (the survival ceiling, the prefactor, the closed refill constant) are listed as corollaries with their own, weaker labels; they cap at Tier 1 candidate.
**Date:** 2026-07-13
**Authors:** Thomas Wicht, Claude (Fable 5)

## What this is about

Release one excitation at the end of a chain and watch it run. Unwatched, it runs as a wave: a ballistic front at the maximum band speed, with the sharp Airy caustic at its leading edge. Watched, every step is noisy, and the naive guess is that the noise simply grinds the wave down into a random walk.

What actually happens is cleaner and stranger, and this proof holds the exact bookkeeping of it. The watched walk is the unwatched wave, repeatedly caught and released: the watching collapses the excitation onto the sites at a fixed rate, and between collapses it runs free. That sentence is not a picture, it is an equation, the renewal representation, and it is exact. Everything the walk-time experiment measured about the surviving front, the halo that rescues it, the ceiling its decay rate approaches, falls out of this one representation by reading it in different limits.

## Abstract

The single-excitation sector of a chain under local Z-dephasing evolves by ρ̇ = −i[h, ρ] − Γ(ρ − diag ρ) with Γ = γ_φ = 4γ (the Absorption Theorem rate for the sector's Hamming-2 coherences). Splitting the generator into a uniformly damped free evolution plus the diagonal refill Γ·diag and resumming the Dyson ladder at the diagonal gives the exact renewal representation

    P_n(t) = e^{−Γt}·S_n(t),   S_n(t) = |G_{n0}(t)|² + Γ ∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S_m(s),

with G the clean propagator (on the infinite chain G_{nm}(τ) = (−i)^{|n−m|}·J_{|n−m|}(2Jτ)). Every refill order carries the same universal prefactor e^{−Γt}; the zeroth term is exactly the coherent front |⟨a_n⟩|². In momentum-Laplace space the equation closes to an explicit Green's function,

    Ŝ(p, z) = 1 / (√(z² + a(p)²) − Γ),   a(p) = 4J·sin(p/2),

which conserves probability at p = 0 and returns the clean Bessel wave at Γ = 0. The representation is the exact solution of the dephasing tight-binding (Haken-Strobl) sector; the experiment's survival ceiling A_∞(γ) = 4 − φ(2J)/γ, the n^{−1/2} prefactor with constant C(γ) = (2π)^{−1/2}·(γ/(γ+J))^{1/4}, the closed refill constant I₁ = 1/12 + ¼∫₀^{2c}Ai(−w)dw, and the diffusive plateau reading are corollaries.

## Statement

Let h be the nearest-neighbour hopping matrix (amplitude J) of a chain, and let the single-excitation density matrix evolve by

    ρ̇ = −i[h, ρ] − Γ·(ρ − diag ρ),   Γ = 4γ.

Then the site populations P_n(t) = ρ_nn(t), for any seed site and any chain (finite or infinite), satisfy exactly

    P_n(t) = e^{−Γt}·S_n(t),
    S_n(t) = |G_{n0}(t)|² + Γ·∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S_m(s),               (★)

where G(τ) = e^{−ihτ} is the clean single-particle propagator and the seed is site 0. On the infinite chain, Fourier transform in space and Laplace transform in time close (★) to

    Ŝ(p, z) = 1 / (√(z² + a(p)²) − Γ),   a(p) = 4J·sin(p/2).                  (☆)

## Setup and the two definitions

The sector: one excitation on N sites, ρ an N×N matrix. Local Z-dephasing at rate γ per site damps every coherence |i⟩⟨j| at −2γ·popcount(i⊕j); two distinct single-excitation basis states differ in exactly two bits, so every off-diagonal damps at the uniform Γ = 4γ = γ_φ and the diagonal is untouched. That uniformity is the Absorption Theorem read inside one sector, and it is the whole reason the sector closes.

Write the generator as L = L₀ + 𝒥 with

    L₀ ρ = −i[h, ρ] − Γρ   (uniformly damped free evolution),
    𝒥 ρ = Γ·diag ρ          (the refill: the watching returns the diagonal).

L₀ is a scalar decay riding the clean unitary evolution, so its propagator is closed: e^{L₀τ}ρ = e^{−Γτ}·G(τ) ρ G(τ)†. The refill 𝒥 is where the physics of the watching sits: at rate Γ, the coherent state is projected onto the site basis and released again.

## The ladder at the diagonal

Dyson-expand e^{Lt} in powers of 𝒥 and read the (n,n) matrix element. The j-th term is a time-ordered chain of j refills at times s₁ < s₂ < … < s_j: between refills the state propagates cleanly and decays uniformly; at each refill only the diagonal survives and re-seeds. Three structural facts assemble (★):

1. **The uniform decay factors out.** Each segment of length τ contributes e^{−Γτ}, and the segment lengths of every j-refill history sum to t. So every ladder order carries the same e^{−Γt}, and P_n(t) = e^{−Γt}·S_n(t) with S_n the undamped ladder sum.

2. **Each refill vertex is a diagonal injection with weight Γ.** A refill at site m turns the incoming state into the population at m times Γ; the next segment propagates |m⟩ cleanly, contributing |G_{nm}(τ)|² to the next population reading. The kernel of one segment is therefore |G|², the clean single-particle population propagator.

3. **The ladder is a renewal.** Summing over j with these kernels is exactly the Volterra fixed-point equation (★): the population at (n, t) is the never-refilled coherent term plus, for every intermediate (m, s), the rate Γ of a last refill at (m, s) times clean re-propagation |G_{nm}(t−s)|² times the full (already-resummed) S_m(s).

The zeroth term identifies physically: the noise-averaged amplitude ⟨a_n⟩ damps at Γ/2 = 2γ (half the coherence rate), so e^{−Γt}|G_{n0}(t)|² = |⟨a_n(t)⟩|² is exactly the coherent front. Everything above it, the j ≥ 1 terms, is the incoherent halo Var(a_n) of the experiment: population the watching converts out of the coherent front, locally, and releases to run again.

## Closing the ladder: the Green's function

On the infinite chain |G_{nm}(τ)|² = J_{|n−m|}(2Jτ)², and the spatial convolution structure of (★) diagonalizes in the momentum p conjugate to the site index. The one identity needed is Graf's addition theorem summed over the displacement,

    Σ_d e^{ipd}·J_d(x)² = J₀(2x·sin(p/2)),

so the per-momentum kernel is K̂(p, τ) = J₀(a(p)·τ) with a(p) = 4J·sin(p/2). The time convolution then Laplace-transforms with L[J₀(aτ)] = 1/√(z² + a²), and the renewal geometric series sums:

    Ŝ(p, z) = (1/√(z² + a²)) · Σ_{j≥0} (Γ/√(z² + a²))^j = 1/(√(z² + a²) − Γ).      (☆)

Two limits pin the object. At p = 0, a = 0 and Ŝ = 1/(z − Γ), so S_total(t) = e^{Γt} and P_total ≡ 1: probability is conserved, the watching moves weight around but loses none. At Γ = 0, Ŝ = 1/√(z² + a²) inverts to the clean wave S_n = J_n(2Jt)². The pole of (☆) at z = √(Γ² − a²) (present for a < Γ, the long-wavelength window) is the diffusive branch: at small p, √(Γ² − a²) ≈ Γ − a²/(2Γ) with a² ≈ 4J²p², so S carries e^{(Γ − 2J²p²/Γ)t}; against the universal e^{−Γt} this is diffusion with D = 2J²/Γ, the Haken-Strobl diffusion constant (the F123 sibling's rate object).

## Verification

- (★) versus direct RK4 Lindblad integration of the sector: maximum deviation 1.6·10⁻⁶ over displacements 5–30 and t ≤ 18 (grid-limited), probability conserved to 10⁻¹⁵. Committed: `simulations/cone_front_survival_asymptote.py`, sections [1]–[2], result file `simulations/results/cone_defect_arrival/front_survival_asymptote.txt`.
- Independently re-derived and reproduced to the digit by two adversarial referee rounds (an independent Lindblad ODE; a per-momentum Volterra solver with independent quadrature; a symbolic re-derivation of the Graf and Laplace steps), 2026-07-13.
- The j = 0 = coherent-front identity: verified as e^{−Γt}J_n(2Jt)² = |⟨a_n⟩|² to 1.4·10⁻⁹ (section [3]).

## Corollaries (each with its own, weaker label)

The experiment `experiments/COUPLING_DEFECT_WALK_TIME_STEP.md` (follow-ups two to five) reads (★)/(☆) in four limits; none of these is part of the Tier-1 statement above:

1. **The survival ceiling** (asymptotic, Gärtner-Ellis on the tilted pole μ(θ) = √(Γ² + 16J²·sinh²(θ/2)), numerically confirmed): the front-peak survival exponent approaches A_∞(γ) = 4 − φ(2J)/γ with φ(2J) = √(Γ(Γ+4J)) − 4J·arcsinh√(Γ/4J); small γ: 4 − (8/3)·√(γ/J). The naive "the front asymptotically pays the full γ_φ" is false at every fixed γ.
2. **The prefactor** (leading saddle, symbolically exact constants): S_n(t*₀) = C(γ)·n^{−1/2}·e^{(φ/2J)n} with C(γ) = (2π)^{−1/2}·(γ/(γ+J))^{1/4}, via μ″(θ*) = 2J·√(Γ/(Γ+4J)).
3. **The refill constant** (rigorous leading-order stationary phase): the single-refill front integral saturates to I₁ = 1/12 + ¼∫₀^{2c}Ai(−w)dw = 0.27694424, 2c = 2^{2/3}·α, α the first zero of Ai′; the pre-asymptotic climb of the exponent is 4 − (8·I₁/0.45547)·n^{−1/3}.
4. **The third reading** (exact): the global-in-time maximum of P_n is the diffusive plateau of the p → 0 pole (peak at t = n²/2D, height e^{−1/2}/√(2π)/n), so the peak-tracking survival exponent is exactly zero; fixed-time, windowed, and peak-tracking readings of "front survival" form the experiment's trichotomy.

## Links

- The measuring experiment and its five follow-ups: `experiments/COUPLING_DEFECT_WALK_TIME_STEP.md`.
- The rate Γ = 4γ: `PROOF_ABSORPTION_THEOREM.md` (the −2γ·k law read inside the sector).
- The clean propagator and band: F2b in `docs/ANALYTICAL_FORMULAS.md` (band 2J·cos q, front speed 2J).
- The sector and its slow-mode EP: `PROOF_COHERENCE_HORIZON_SLOPE.md` (the same Liouvillian at the band edge; this proof solves the sector the horizon claim locates the EP of).
- The diffusive sibling: `SurvivorDiffusionGradientClaim` / F123 (D = 2J²/Γ is the (☆) pole's diffusion constant).
- The typed home: `DephasingFrontRenewalClaim` (live witness `inspect --root renewal`).
