# The Coherence Horizon asymptotic slope is 2/π

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
- r = −1:  the left ladder decays as (−μ) (the sign-flipped root), so m_{−2} = −μ·m_{−1}, giving
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
   4γ telegrapher form scatters 3.6× worse (0.030) — the resummed dispersion is the fixed mode. (b) The
   coalescing-pair trace −2Re(λ)/γ climbs monotonically toward 8 (7.18 at N=24, 7.71 at N=60) with no mode
   anywhere near 4. Incidentally, placing the shift-invert shift at σ = −4γ makes the factorization exactly
   singular — the slow mode sits precisely there.
3. **q_min·N → π.** From Q* = 2/q_min, the data give q_min·N = 2N/Q* decreasing monotonically toward π
   (3.86 at N=12, 3.58 at N=28, 3.32 at N=120), i.e. Q*/N → 2/π from below.

## Numerical confirmation

The single-excitation EP, by dense eigendecomposition (N ≤ 30) and sparse shift-invert (N ≤ 120,
`coherence_horizon_slope.py` / the SE block in `coherence_horizon_se_block.py`). The discrete slope
dQ* = Q*(N) − Q*(N−1) climbs monotonically toward 2/π from below:

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
- F-registry: the F2b "two clocks" corollary, `docs/ANALYTICAL_FORMULAS.md`.
- Verifiers: `simulations/coherence_horizon_slope.py` (self-validating: slope refutes √2/π, q_min·N → π, the
  8γ dispersion via overdamped real roots) and `simulations/coherence_horizon_slope_largeN.py` (the sparse
  shift-invert sweep to N=120), on top of `simulations/coherence_horizon_se_block.py` (the SE-EP reduction,
  self-validating N=2..8).

## Review

Round 1 (2026-06-15, independent adversarial agent, physics-first): **GO**. Re-derived the dispersion in
sympy (the m₋₁ = −μ·m₀ left-ladder sign and the r=0 closure confirmed; the eigenvector's r=+1 / r=−1
amplitude ratio is −1.0000 at every N — no error coincidentally producing 2/π). Ran the 8γ discriminator
independently: −2Re(λ)/γ for the slow mode climbs 7.18 (N=24) → 7.71 (N=60) with no mode near 4, decisively
excluding the 4γ telegrapher (slope √2/π). Iterated Richardson on the slope → 0.6363 = 2/π. Confirmed
chain-specific (built the ring: slope → 1/π, half the chain, exactly as q_min = 2π/N vs π/N predicts).
Caveat banked: extrapolation alone is inconclusive at N ≤ 70 (low-order fits land 0.617–0.633); the value is
fixed by the 8γ derivation, not the fit. The overdamped diffusion constant D = J²/2γ is shared by both the 4γ
and 8γ models, so it is NOT the discriminator (the EP / conjugate-pair-sum is) — stated in cross-check #1.
