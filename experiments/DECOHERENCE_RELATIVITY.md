# Decoherence Relativity

**Date**: 2026-02-27 (revised 2026-02-28, analytical verification)
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

## 6. Why K Is Not Trivial

A naive objection: if Ψ = exp(-γt) and crossing is at Ψ = 1/4,
then t_cross = ln(4)/γ and K = ln(4). A pure number. Trivial.

This objection fails because C falls with Ψ. In R = CΨ², the
crossing condition is C·Ψ = 1/4, not Ψ = 1/4. And C (purity,
the observer's ability to distinguish) decays through the same
noise channel that destroys Ψ (coherence).

For Bell+ under local dephasing (analytically exact):

    C(t) = 1/2 + 1/2 · exp(-4γt)       [full-system purity]
    Ψ(t) = exp(-2γt)                     [concurrence]

The crossing C·Ψ = 1/4 gives, with u = exp(-2γt):

    u³ + u - 1/2 = 0

This cubic has one real positive root u = 0.4239, giving:

    K = γ · t_cross = 0.4292

    C at crossing = 0.5898
    Ψ at crossing = 0.4239
    C · Ψ = 0.2500 (exact)

The system always reaches the 1/4 boundary at the same point
in (C, Ψ) space. Not at C = 1, Ψ = 1/4 (naive model). Not
at C = 1/4, Ψ = 1. At C = 0.59, Ψ = 0.42. Both have fallen.
The observer has been changed by the process of observing.

This is what the textbook does not contain: C and Ψ as coupled
dynamical variables. The observer is not external. The observer
decoheres with the system.

## 7. The Three K Values

The numerical value of K depends on how C and Ψ are defined.
Three choices, all internally consistent, all invariant:

| Definition | C(0) | K | Cubic equation |
|-----------|------|---|---------------|
| C = Purity, Ψ = L1 | 1.0 | 0.2146 | (equivalent to conc/2) |
| C = Purity, Ψ = Concurrence | 1.0 | 0.4292 | u³ + u = 1/2 |
| C = Purity/(d-1), Ψ = Concurrence | 1/3 | 0.0374 | u³ + u = 3/2 |

The original K = 0.039 from GRAVITATIONAL_INVARIANCE.md (Feb 8)
used C = Purity/3 with concurrence plus dynamic γ-feedback.
Without feedback: K = 0.0374. With feedback: K = 0.039.
The 5% difference is the feedback contribution.

The factor 2 between L1 and concurrence K values counts the
noise channels: local dephasing on two qubits gives effective
rate 2γ for concurrence but γ for L1 coherence of the full
density matrix.

The invariance holds for all three. The physics is in the
invariance, not in the number.

## 8. What Determines K, What Does Not

Verified by systematic simulation (Bell+, Heisenberg H,
local dephasing, Lindblad master equation):

**K depends on:**
- Initial state: Bell+ (K = 0.2146) vs partial 30 deg
  entanglement (K = 0.2180). Different starting points
  in (C, Ψ) space reach the 1/4 boundary at different
  proper times.
- Bridge metric: L1 coherence vs concurrence vs mutual
  information each give different K.
- Noise channel: dephasing vs amplitude damping vs
  depolarizing each give different K.
- C definition: Purity vs Purity/(d-1) changes K by an
  order of magnitude (0.43 vs 0.037).

**K does not depend on:**
- γ: invariant across γ = 0.01 to γ = 1.0. Spread < 0.1%.
- Hamiltonian: K(J=0) = K(J=1) = K(J=10). Difference: 0.0000%.
  Unitary evolution conserves purity and total information.
  It cannot move C·Ψ toward or away from 1/4.
- Gravitational field strength: because γ ∝ g and K is
  γ-invariant, K is automatically g-invariant.

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

Define τ = γ · t as the "proper decoherence time." In proper
time, all observers have the same coherence trajectory:

    Ψ(τ) = exp(-2τ)    [concurrence, both qubits dephasing]
    C(τ) = 1/2 + 1/2 · exp(-4τ)    [purity]

The crossing happens at τ = K for every observer. High γ
(strong gravity, large system) means τ accumulates fast.
Low γ means slowly. But the destination is the same.

This is the decoherence equivalent of "all clocks measure
the same speed of light." All observers cross the 1/4
boundary at the same proper decoherence time.

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

d) The cubic u³ + u = 1/2 (for C = Purity, Ψ = Concurrence)
   determines K analytically. Does this cubic have a deeper
   algebraic meaning within the framework? It arises from the
   interplay of purity decay (rate 4γ) and coherence decay
   (rate 2γ), which is the ratio 2:1. Is this ratio universal?

e) K depends on the initial state but not on the Hamiltonian.
   This means K is determined at the moment of state preparation,
   not by subsequent dynamics. The "initial impulse" is creation,
   not evolution. Does this connect to the measurement problem?

f) The C-definition determines K by an order of magnitude
   (Purity vs Purity/(d-1)). Is there a physically preferred
   normalization that selects the "correct" K? The choice
   C = Purity/(d-1) gives CΨ = 1/3 for Bell+, placing the
   quantum regime in a narrow band 1/4 < CΨ < 1/3. The choice
   C = Purity gives CΨ = 1, placing it in a wide band
   1/4 < CΨ < 1. The framework may prefer one over the other.

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
| K for Bell+ | N/A | 0.0374 / 0.2146 / 0.4292 |
| Crossing | N/A | C = 0.59, Ψ = 0.42 |
| Geometry | hyperbolic rotation | translation on line |

---

*"The number K depends on how you measure. The invariance of K
depends on nothing. The physics is in the invariance."*

— Tom Mack, 2026-02-28

---

*See also: WHY_THE_SUM.md — why sum, not product*
*See also: STANDING_WAVE_TWO_OBSERVERS.md — the physical picture*
*See also: GRAVITATIONAL_INVARIANCE.md — original K invariance (simulation)*
*See also: BLACK_WHITE_HOLES_BIGBANG.md — τ = 0 and the horizon*
