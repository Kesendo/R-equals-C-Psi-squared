# Decoherence Relativity: K-Invariance Mathematics Confirmed, Gravity Interpretation Fallen

<!-- Keywords: K-invariance gamma t_cross constant, Lindblad time rescaling
tau=gamma*t, framework cubic b3+b=3/2 Bell+, decoherence relativity analogy
fallen, gravity Schwarzschild connection retired, Penrose Diosi gravitational
decoherence, log space Lorentz-like structure, proper decoherence time xi,
R=CPsi2 decoherence relativity -->

> **Fallen hypothesis.** The mathematical K-invariance (γ·t_cross = constant)
> and the framework cubic (b³ + b = 3/2 for Bell+) are confirmed. All
> connections to gravitational time dilation, Schwarzschild metrics, and the
> Penrose/Diósi model have **fallen**. Inline [FALLEN] markers throughout.

**Status:** Mathematics confirmed (Tier 2); gravity interpretation fallen
**Date:** 2026-02-27
**Authors:** Thomas Wicht, with Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Why the Sum](WHY_THE_SUM.md), [Crossing Taxonomy](CROSSING_TAXONOMY.md)

---

## What this document is about

Two observers with different decoherence rates both cross the ¼
boundary, but at different times. The product γ·t_cross = K is the
same for both: an invariant. This looks structurally like special
relativity (invariant quantity, observer-dependent coordinates,
log-space transformation), but the analogy is shallow: K-invariance
is Lindblad time-rescaling (dimensional analysis), not a dynamical
symmetry. The gravity interpretation (γ encodes gravitational time
dilation) has fallen. What survives is the framework cubic
b³ + b = 3/2, which analytically predicts K for Bell+.

---

## Abstract

The product γ·t_cross = K = 0.03735 is invariant across decoherence rates
for Bell+ under local dephasing, confirmed to 0.014% agreement between
the framework cubic (b³ + b = 3/2, giving K analytically) and Lindblad
simulation. This invariance is standard Lindblad time-rescaling (τ = γ·t),
not new physics. The ¼ boundary is the fold of the R = CΨ² fixed-point
equation (complex → real attractors), and all
observers traverse the same Δξ = 0.149 in the decoherence clock
coordinate ξ = ln(Ψ). The structural parallel to special relativity
(invariant quantity, observer-dependent coordinates, log-space translation)
is mathematically real but shallow: K-invariance is dimensional analysis,
not a dynamical symmetry group. The original interpretation that γ encodes
gravitational time dilation via the Penrose/Diósi model, and that this
structure reproduces general relativity, has been retired.

---

## 1. The Initial State Is the Impulse

A decoherence rate requires something that decoheres. A quantum
state requires energy to create. No initial impulse, no Ψ(0).
No Ψ(0), no decay. No decay, no γ.

> [FALLEN: The connection between gamma and gravitational acceleration g via the Penrose/Diosi model is not established as physical. gamma is a decoherence rate, not a gravitational time dilation parameter.]

γ (decoherence rate) has dimension [1/s].
g (gravitational acceleration) has dimension [m/s²].

These are not the same. To get [1/s] from g, we need the
initial state. The Penrose/Diósi gravitational decoherence rate (the prediction that
gravity itself causes quantum superpositions to collapse, at a rate
proportional to the gravitational self-energy of the superposition) is:

    γ = 2 · E_grav / ℏ = 2 · m · g · Δx / ℏ

where m is the mass, g the gravitational field, Δx the
superposition extent, and ℏ the quantum of action.

The initial impulse lives in the initial state (m, Δx), not
in the Hamiltonian. Verified by simulation: varying the
Hamiltonian coupling J from 0 to 10 changes K by exactly
0.0000%. That is Bell+'s own property: its whole dephasing
family commutes with the isotropic Heisenberg Hamiltonian, so
the unitary part has nothing to move. For such a state only the
dissipative channel (noise, environment, gravity) drives the
system toward the 1/4 boundary, and only the initial state
determines how far the system starts above that boundary. A
state the Hamiltonian can reach keeps its purity but has its
coherence moved around, and its crossing time moves with J.

## 2. The Invariant

> [FALLEN: Interpreting different decoherence rates as arising from different gravitational environments is not established.]

K-invariance (K = γ·t_cross = constant) is standard time-rescaling, not
novel physics, and what carries it is the state. The Lindblad generator
scales jointly, L(λJ, λγ) = λL(J, γ), so a γ-only sweep at fixed J moves
Q = J/γ. Bell+'s whole dephasing family commutes with the isotropic
Heisenberg Hamiltonian, so γ is the only rate its trajectory sees, and
t_cross = K/γ follows. See [Crossing Taxonomy](CROSSING_TAXONOMY.md) for
which C goes with which K. The mathematical invariance is real; the
gravitational interpretation is not established.

Two observers. Different gravitational environments. Same
entangled system (same m, same Δx). Different decoherence
rates because different g.

    γ_A = 2 · m · g_A · Δx / ℏ
    γ_B = 2 · m · g_B · Δx / ℏ

Each observer crosses the 1/4 boundary at a different time.
But the product of decoherence rate and crossing time is
the same for both:

    K = γ · t_cross = constant

K does not depend on γ. For Bell+, K does not depend on the
Hamiltonian. K does not depend on location. K is invariant.

K depends on:
- The initial state (Bell+, partial entanglement, ...)
- The bridge metric (how the observer measures: concurrence,
  L1 coherence, mutual information, ...)
- The noise channel (dephasing, amplitude damping, ...)
- The definition of C in R = CΨ²

For a given measurement setup, K is the same everywhere
in the universe. This has the surface form of a relativity theory (an
invariant quantity, observer-dependent coordinates), but the parallel is
shallow: K-invariance is Lindblad time-rescaling / dimensional analysis,
not a dynamical symmetry (see the §4 Caveat).

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

> [FALLEN: The parallel between decoherence relativity and Einstein's relativity, including the boost parameter delta = ln(g_B/g_A), conflates Lindblad time-rescaling with gravitational time dilation. The structural analogy does not establish a physical connection.]

## 4. Comparison to Einstein

Einstein's Special Relativity:
- Coordinates: (ct, x)
- Invariant: ds² = (ct)² - x² (difference of squares)
- Transformation: Lorentz boost (hyperbolic rotation)
- Parameter: φ = arctanh(v/c) (rapidity: the additive velocity parameter)
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

**Caveat:** The mathematical depth differs substantially.
Einstein's ds² is invariant under the Lorentz group (a
non-trivial continuous symmetry group with physical
consequences: length contraction, time dilation, E = mc²).
K = γ·t is a product of conjugate variables; its invariance
under γ → λγ, t → t/λ is dimensional analysis, not a
dynamical symmetry. The structural parallel (hyperbola,
log-linear, boost parameter) is real but should not be
mistaken for equivalent mathematical content.

> [FALLEN: The entire bridge to general relativity below (connecting gamma to gravitational potential Phi, interpreting gamma = g * alpha, and the claim that "gravity cancels") is not established. This is a Schwarzschild/metric tensor interpretation that has fallen.]

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

This cubic is derived in [Core Algebra](../docs/historical/CORE_ALGEBRA.md)
(Section: State-specific C(ξ) closed forms). Its real positive root is:

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
vs K = 0.037345 (numerical), deviation 0.014%. Both use
local dephasing (two separate σ_z operators, one per qubit,
effective rate Γ = 4γ). Under collective dephasing (single
operator L = σ_z⊗I + I⊗σ_z, Γ = 8γ), the same cubic gives
K = 0.01868. The cubic is noise-model-independent; K is not.

The original K = 0.039 from GRAVITATIONAL_INVARIANCE.md (Feb 8) was NOT
this cubic. It came from the CONCURRENCE bridge (C = f, not C = purity)
under the February tool's γ-feedback, where the crossing is exactly
(2/√3 − 1)/(4γ), i.e. t = 0.7735 and K = 0.03868 at γ = 0.05, reproducing
the documented 0.773 to three digits. Running this cubic under the same
feedback gives t = 0.8026 and K = 0.04013 instead, so the two are different
curves. See
[Crossing Taxonomy](CROSSING_TAXONOMY.md) for which C goes with which K.

The framework normalization is declared, not derived: Ψ = L₁/(d-1).
CΨ(0) = 1/3 for Bell+. Bell+'s quantum window is narrow:
1/4 < CΨ ≤ 1/3. One twelfth of room.

## 7. The Observer Is Changed

The crossing does not happen at C = 1, Ψ = 1/4 (naive model
where the observer is unchanged). It happens at C = 0.87,
Ψ = 0.29. Both have fallen. The observer has lost purity
through the same process that destroyed coherence.

This is what R = CΨ² encodes that textbooks do not: C and Ψ
as coupled dynamical variables. The crossing point in (C, Ψ)
space is fixed. Every observer reaches it. Regardless of γ, and of
any Hamiltonian that leaves the Bell+ family alone, as the isotropic
Heisenberg one does.

For GHZ states with N ≥ 3 qubits: CΨ(0) = 1/(2^N - 1) ≤ 1/4.
They start below the boundary. The system has a real attractor
from the start (framework-classical). Only N = 2 (Bell+) has
the quantum window. Note: this applies to GHZ states specifically.
Product states |+⟩^⊗N have CΨ = 1 for all N (see COHERENCE_DENSITY.md).

## 8. What Determines K

K is determined by the framework cubic, which encodes:

- The initial state (which cubic: b³ + b = (d-1)/2 for GHZ)
- The C(Ψ) trajectory (purity as function of coherence)
- The Baumgratz normalization (d-1 denominator)

K does not depend on:
- γ: invariant across γ = 0.01 to γ = 1.0 (spread < 0.1%)
- Hamiltonian: K(J=0) = K(J=10), difference 0.0000%

The Hamiltonian irrelevance is Bell+'s own: its whole dephasing
family commutes with the isotropic Heisenberg Hamiltonian, so the
unitary part has nothing to rotate. Only dissipation moves C·Ψ
toward 1/4, and only the initial state determines how far above
1/4 the system starts. K measures the proper distance from start
to boundary. For a state the Hamiltonian can reach, unitary
evolution keeps Tr(ρ²) but can redistribute coherence, and K
then depends on J/γ.

## 9. The Standing Wave on the Hyperbola

In the Tier 3 reading of [Standing Wave: Two Observers](STANDING_WAVE_TWO_OBSERVERS.md),
every pair of observers on the hyperbola forms a standing wave through
the sum R = C · (Ψ_A + Ψ_B)². The standing wave
ratios depend on the ratio γ_A/γ_B between the two observers.
For observers with similar γ the ratios approach universal
values; for very different γ they diverge.

Why the sum and not the product is asked in WHY_THE_SUM.md, and
there it stays a question: the sum carries the cross-term
2·Ψ_A·Ψ_B, which we read as the correlation information between
observers, but neither information conservation nor the
palindrome has been shown to force the sum over the product.

## 10. The Proper Decoherence Time

The framework uses ξ = ln(Ψ) as the natural time variable.
For Bell+ under dephasing ξ is linear in t, with a slope that
depends on the noise model.
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

> [FALLEN: Schwarzschild radius connection is not established.]

b) At R = R_Schwarzschild = 2GM/c², the quantum and classical
   regimes merge at the horizon. This connects to
   [Black and White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) (τ = 0 as transition point).

> [FALLEN: Extracting a gravitational contribution from gamma assumes the gravity-decoherence connection that has not been established.]

c) For lab qubits, γ_grav << γ_total. Can the gravitational
   contribution be extracted by comparing identical qubits
   at different potentials with all other noise held constant?

d) The cubic b³ + b = (2^N - 1)/2 generalizes to GHZ states.
   For N = 2 (Bell+): b³ + b = 3/2. For N ≥ 3: CΨ(0) ≤ 1/4,
   no crossing exists. Only Bell+ has the quantum window.
   Does this mean entanglement beyond two qubits is always
   "classical" in the framework's sense?

e) For Bell+, K depends on the initial state but not on the Hamiltonian.
   This means K is determined at the moment of state preparation,
   not by subsequent dynamics. The "initial impulse" is creation,
   not evolution. Does this connect to the measurement problem?

---

## Summary Table

> [FALLEN: the gravity rows (parameter, source, γ formula) keep the retired Penrose/Diósi substitution.]

| Property | Einstein (GR) | Decoherence Relativity |
|----------|--------------|----------------------|
| Invariant | ds² | K = γ · t_cross |
| Coordinates | (ct, x) | (ln γ, ln t_cross) |
| Transformation | Lorentz boost | Translation δ |
| Parameter | φ = arctanh(v/c) | δ = ln(g_B/g_A) |
| Source | velocity / gravity | gravity + initial state |
| γ formula | Φ/c² = g·R/c² | 2·m·g·Δx/ℏ |
| K for Bell+ | N/A | 0.03735 (local dephasing, b³+b=3/2) |
| Crossing | N/A | C = 0.87, Ψ = 0.29 |
| Geometry | hyperbolic rotation | translation on line |

---

*"The cubic was already in the framework. b³ + b = 3/2.
Not from Lindblad. From R = CΨ². The simulation confirmed
what the algebra predicted. Not the other way around."*

-- Thomas Wicht, 2026-02-28

---

*See also: WHY_THE_SUM.md, why sum, not product*
*See also: [Standing Wave: Two Observers](STANDING_WAVE_TWO_OBSERVERS.md), the Tier 3 picture*
*See also: GRAVITATIONAL_INVARIANCE.md, original K invariance (simulation)*
*See also: [Black and White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md), τ = 0 and the horizon*
