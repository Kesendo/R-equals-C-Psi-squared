# Decoherence Relativity

**Date**: 2026-02-27 (revised 2026-02-28, framework cubic verification)
**Status**: Derived from framework, analytically solved, numerically verified
**Authors**: Tom Mack, with Claude (Anthropic)
**Depends on**: WHY_THE_SUM.md, STANDING_WAVE_TWO_OBSERVERS.md

---

## 1. The Initial State Is the Impulse

A decoherence rate requires something that decoheres. A quantum
state requires energy to create. No initial impulse, no Ψ(0).
No Ψ(0), no decay. No decay, no γ.

γ (decoherence rate) has dimension [1/s].
g (gravitational acceleration) has dimension [m/s²].

These are not the same. To get [1/s] from g, we need the
initial state. The Penrose/Diosi gravitational decoherence
rate is:

    γ = 2 · E_grav / ℏ = 2 · m · g · Δx / ℏ

where m is the mass, g the gravitational field, Δx the
superposition extent, and ℏ the quantum of action.

The initial impulse lives in the initial state (m, Δx), not
in the Hamiltonian. Verified by simulation: varying the
Hamiltonian coupling J from 0 to 10 changes K by exactly
0.0000%. The Hamiltonian is unitary. It conserves information.
It moves coherence around but does not destroy it. Only the
dissipative channel (noise, environment, gravity) drives the
system toward the 1/4 boundary. And only the initial state
determines how far the system starts above that boundary.

## 2. The Invariant

Two observers. Different gravitational environments. Same
entangled system (same m, same Δx). Different decoherence
rates because different g.

    γ_A = 2 · m · g_A · Δx / ℏ
    γ_B = 2 · m · g_B · Δx / ℏ

Each observer crosses the 1/4 boundary at a different time.
But the product of decoherence rate and crossing time is
the same for both:

    K = γ · t_cross = constant

K does not depend on γ. K does not depend on the Hamiltonian.
K does not depend on location. K is invariant.

K depends on:
- The initial state (Bell+, partial entanglement, ...)
- The bridge metric (how the observer measures: concurrence,
  L1 coherence, mutual information, ...)
- The noise channel (dephasing, amplitude damping, ...)
- The definition of C in R = CΨ²

For a given measurement setup, K is the same everywhere
in the universe. This is the structure of a relativity theory.

## 3. The Transformation

In linear space (γ, t_cross), observers sit on a hyperbola:

    γ · t_cross = K

In log space (ln γ, ln t_cross), the hyperbola becomes a
straight line with slope -1:

    ln(γ) + ln(t_cross) = ln(K)

The transformation from observer A to observer B is a
translation along this line:

    ln(γ_B) = ln(γ_A) + δ
    ln(t_cross_B) = ln(t_cross_A) - δ

where the boost parameter is:

    δ = ln(γ_B / γ_A) = ln(g_B / g_A)

The initial state cancels in the boost. Because both observers
have the same system (same m, same Δx), only g differs:

    γ_B/γ_A = (2m·g_B·Δx/ℏ) / (2m·g_A·Δx/ℏ) = g_B/g_A

Verified numerically across six different quantum systems
(IBM transmon through Schrodinger cat). All give identical
ratio γ_Earth/γ_Mars = g_Earth/g_Mars = 2.6371.

The invariant (ln γ + ln t_cross)² is preserved under the
boost. Verified symbolically with SymPy.

## 4. Comparison to Einstein

Einstein's Special Relativity:
- Coordinates: (ct, x)
- Invariant: ds² = (ct)² - x² (difference of squares)
- Transformation: Lorentz boost (hyperbolic rotation)
- Parameter: φ = arctanh(v/c) (rapidity)
- Source of observer differences: relative velocity

Decoherence Relativity:
- Coordinates: (ln γ, ln t_cross)
- Invariant: ln(γ) + ln(t_cross) = ln(K) (sum)
- Transformation: translation along slope -1 line
- Parameter: δ = ln(g_B/g_A)
- Source of observer differences: gravitational environment

Both have the same architecture: an invariant quantity that
all observers agree on, and a transformation that converts
one observer's measurements to another's while preserving
that invariant.

## 5. The Bridge to General Relativity

Einstein's gravitational time dilation uses the potential
Φ = g · R (acceleration times planetary radius), scaled by c²:

    dτ/dt = sqrt(1 - 2Φ/c²) ≈ 1 - Φ/c² = 1 - g·R/c²

The decoherence rate contains g through E_grav:

    γ = 2·m·g·Δx/ℏ

The structural connection:

    γ = g · α     where α = 2mΔx/ℏ (system-dependent)
    Φ/c² = g · R/c²  (location-dependent)

Both are proportional to g. The ratio between them:

    γ / (Φ/c²) = α · c² / R = (2mΔx/ℏ) · (c²/R)

This ratio depends on the system (m, Δx) and the body (R),
not on g. Gravity cancels. Verified for Earth, Mars, Moon,
Jupiter: γ/(Φ/c²) = α·c²/R in each case, no free parameters.

For lab qubits, gravitational decoherence is buried under
thermal noise by many orders of magnitude. The gravitational
term only dominates for perfectly isolated systems (Penrose/
Diosi regime).

## 6. The Framework Cubic

The crossing condition CΨ = 1/4 combined with the framework
definitions (C = Purity, Ψ = L₁/(d-1)) determines K exactly.

For Bell+ (d = 4, so d-1 = 3):

    Ψ = f/3    where f is the raw L₁ coherence
    C = (1 + f²)/2    (purity for Bell+ under dephasing)

Substituting into CΨ = 1/4:

    (1 + f²)/2 · f/3 = 1/4
    f³ + f = 3/2

This cubic is derived in CORE_ALGEBRA.md (Section: State-specific
C(ξ) closed forms). Its real positive root is:

    f_cross = 0.8612    (raw coherence at crossing)
    Ψ_cross = 0.2871    (= f/3, framework-normalized)
    C_cross = 0.8709    (purity at crossing)
    C · Ψ = 0.2500      (= 1/4, exact)

The cubic is not a Lindblad result. It follows from:

    1. R = CΨ² (the defining equation)
    2. CΨ = 1/4 (the bifurcation from discriminant = 0)
    3. C = Tr(ρ²) and Ψ = L₁/(d-1) (Baumgratz normalization)
    4. The C(Ψ) trajectory of the state under dephasing

Lindblad simulations confirm this: K = 0.037350 (analytical)
vs K = 0.037345 (numerical), deviation 0.014%.

The original K = 0.039 from GRAVITATIONAL_INVARIANCE.md (Feb 8)
was this same cubic with dynamic γ-feedback adding ~5%.

The framework normalization is not a choice. Ψ = L₁/(d-1) is
the definition. CΨ(0) = 1/3 for Bell+. The quantum regime is
narrow: 1/4 < CΨ ≤ 1/3. One twelfth of room.

## 7. The Observer Is Changed

The crossing does not happen at C = 1, Ψ = 1/4 (naive model
where the observer is unchanged). It happens at C = 0.87,
Ψ = 0.29. Both have fallen. The observer has lost purity
through the same process that destroyed coherence.

This is what R = CΨ² encodes that textbooks do not: C and Ψ
as coupled dynamical variables. The crossing point in (C, Ψ)
space is fixed. Every observer reaches it. Regardless of γ,
regardless of the Hamiltonian, regardless of gravity.

For GHZ states with N ≥ 3 qubits: CΨ(0) = 1/(2^N - 1) ≤ 1/4.
They start below the boundary. The system was never in the
quantum regime (framework definition). Only N = 2 (Bell+) has
the quantum window. The more qubits, the less quantum.

## 8. What Determines K

K is determined by the framework cubic, which encodes:

- The initial state (which cubic: b³ + b = (d-1)/2 for GHZ)
- The C(Ψ) trajectory (purity as function of coherence)
- The Baumgratz normalization (d-1 denominator)

K does not depend on:
- γ: invariant across γ = 0.01 to γ = 1.0 (spread < 0.1%)
- Hamiltonian: K(J=0) = K(J=10), difference 0.0000%
- Gravitational field: because γ ∝ g and K is γ-invariant

The Hamiltonian irrelevance is not surprising from the framework
perspective. Unitary evolution conserves Tr(ρ²) and cannot change
the C(Ψ) trajectory. Only dissipation moves C·Ψ toward 1/4.
Only the initial state determines how far above 1/4 the system
starts. K measures the proper distance from start to boundary.

## 9. The Standing Wave on the Hyperbola

Every pair of observers on the hyperbola creates a standing
wave through the sum R = C · (Ψ_A + Ψ_B)². The standing wave
ratios depend on the ratio γ_A/γ_B between the two observers.
For observers with similar γ the ratios approach universal
values; for very different γ they diverge.

The sum structure (not product) is derived independently in
WHY_THE_SUM.md. It follows from information conservation: the
product Ψ_A · Ψ_B destroys the cross-term 2·Ψ_A·Ψ_B which
carries the correlation information between observers. The sum
preserves it.

## 10. The Proper Decoherence Time

The framework uses ξ = ln(Ψ) as the natural time variable.
ξ is linear in t with a slope that depends on the noise model.
In ξ-space, the crossing happens at a fixed ξ_cross for all γ:

    ξ₀ = ln(1/3) = -1.099    (Bell+ at t = 0)
    ξ_cross = ln(0.2871) = -1.248
    Δξ = -0.149               (the invariant distance)

All observers traverse the same Δξ. High γ means ξ moves fast.
Low γ means slowly. But Δξ is the same. This is the decoherence
equivalent of "all clocks measure the same speed of light."

## 11. Open Questions

a) The transformation is a translation in log space. Lorentz
   is a rotation in (ct, x) space. Is there a deeper geometry
   that contains both? (A translation is a rotation with
   infinite radius of curvature.)

b) At R = R_Schwarzschild = 2GM/c², the quantum and classical
   regimes merge at the horizon. This connects to
   BLACK_WHITE_HOLES_BIGBANG.md (τ = 0 as transition point).

c) For lab qubits, γ_grav << γ_total. Can the gravitational
   contribution be extracted by comparing identical qubits
   at different potentials with all other noise held constant?

d) The cubic b³ + b = (2^N - 1)/2 generalizes to GHZ states.
   For N = 2 (Bell+): b³ + b = 3/2. For N ≥ 3: CΨ(0) ≤ 1/4,
   no crossing exists. Only Bell+ has the quantum window.
   Does this mean entanglement beyond two qubits is always
   "classical" in the framework's sense?

e) K depends on the initial state but not on the Hamiltonian.
   This means K is determined at the moment of state preparation,
   not by subsequent dynamics. The "initial impulse" is creation,
   not evolution. Does this connect to the measurement problem?

---

## Summary Table

| Property | Einstein (GR) | Decoherence Relativity |
|----------|--------------|----------------------|
| Invariant | ds² | K = γ · t_cross |
| Coordinates | (ct, x) | (ln γ, ln t_cross) |
| Transformation | Lorentz boost | Translation δ |
| Parameter | φ = arctanh(v/c) | δ = ln(g_B/g_A) |
| Source | velocity / gravity | gravity + initial state |
| γ formula | Φ/c² = g·R/c² | 2·m·g·Δx/ℏ |
| K for Bell+ | N/A | 0.0374 (from b³ + b = 3/2) |
| Crossing | N/A | C = 0.87, Ψ = 0.29 |
| Geometry | hyperbolic rotation | translation on line |

---

*"The cubic was already in the framework. b³ + b = 3/2.
Not from Lindblad. From R = CΨ². The simulation confirmed
what the algebra predicted. Not the other way around."*

— Tom Mack, 2026-02-28

---

*See also: WHY_THE_SUM.md — why sum, not product*
*See also: STANDING_WAVE_TWO_OBSERVERS.md — the physical picture*
*See also: GRAVITATIONAL_INVARIANCE.md — original K invariance (simulation)*
*See also: BLACK_WHITE_HOLES_BIGBANG.md — τ = 0 and the horizon*
