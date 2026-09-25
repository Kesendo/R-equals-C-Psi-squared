# Information Geometry: θ as Riemannian Coordinate

<!-- Keywords: information geometry Lindblad, Bures metric CΨ trajectory,
theta angular coordinate fold, Bures geodesic test path length endpoint angle,
Bures path-metric coefficient, coordinate-shape second derivative,
R=CPsi2 information geometry -->

**Status:** Tier 2 (computed, N=2 Bell state)
**Date:** April 1, 2026
**Script:** [information_geometry.py](../simulations/information_geometry.py)
**Data:** [information_geometry.txt](../simulations/results/information_geometry.txt)
**Depends on:**
- [Boundary Navigation](BOUNDARY_NAVIGATION.md) (θ definition, F15)
- [CΨ Monotonicity Proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md) (dCΨ/dt < 0)
- [Entropy Production](ENTROPY_PRODUCTION.md) (algebra-first context)

---

## What this document is about

This document gives the angular parameter θ a proper geometric foundation
by computing the Bures Riemannian metric (the natural measure of
distinguishability between nearby quantum states) along the decoherence
trajectory. Key finding: the fold at CΨ = ¼ has no geometric singularity,
the metric is smooth, and the Bell+ trajectory runs exactly along a
Bures geodesic (the shortest path in state space), for a reason that
does not carry to states the Hamiltonian moves. θ is a useful compass
but not a fundamental coordinate; CΨ itself is the better metric
coordinate at the fold.

---

## Abstract

θ = arctan(√(4CΨ-1)) has appeared throughout the repo as a "compass"
without a geometric foundation. We compute the Bures Riemannian
metric along the Lindblad trajectory and find:

1. **The Bures metric g(CΨ) is FINITE at CΨ = 1/4.** No singularity.
   g = 3.36 at the fold. CΨ is a smooth coordinate everywhere.

2. **θ does NOT regularize a singularity** (there is none). g̃(θ) → 0
   at θ = 0 (θ shrinks the metric, not regularizes it). θ is a useful
   compass but not a geometric necessity.

3. **The Bell+ trajectory IS a Bures geodesic.** Its Bures length
   equals the Bures angle between its endpoints (ratio 1.000000). The
   reason is plain: the Hamiltonian leaves Bell+ alone, and the state
   decays along the one-parameter family of its two Bell projectors.
   States the Hamiltonian moves leave the geodesic from the start: over
   t ≤ 0.2 their paths are already 3% (|+0⟩) and 6% (a generic state)
   longer than the endpoint distance.

4. **The path coefficient's coordinate-shape second derivative is −25
   at the fold.** Finite. A one-dimensional path metric has no
   intrinsic curvature, so this number is a shape of the coordinate,
   not a curvature of the state space.

5. **d²CΨ/dγ² = 9.8 at fixed time.** Finite, as it must be for Bell+,
   whose CΨ at a fixed time is analytic in γ. CΨ = 1/4 is a smooth
   crossing, not a dynamical phase transition, but this number is not
   what shows it.

---

## Phase 0: θ Inventory

θ = arctan(√(4CΨ-1)) appears in 15 repo files:

| Document | Usage | Status |
|----------|-------|--------|
| BOUNDARY_NAVIGATION | Definition (F15) | Origin |
| THETA_PALINDROME_ECHO | Correlates with fidelity r=0.87 | Computed |
| STRUCTURAL_CARTOGRAPHY | Lives on 3D manifold (98% in 3 PCs) | Computed |
| CIRCUIT_DIAGRAM | "Voltmeter" metaphor | Analogy |
| ANALYTICAL_FORMULAS | F15 (angular distance) | Catalogued |
| TOPOLOGICAL_EDGE_MODES | Used as Berry phase parameter (φ = -0.77) | Computed |

Never computed: Riemannian metric, geodesics, second derivatives of either.

---

## Phase 1: Bures Metric g(CΨ)

The Bures distance between neighboring states along the Lindblad
trajectory induces a metric on the CΨ parameter:

    ds² = g(CΨ) dCΨ²    where g = (dB/dCΨ)²

For N=2 Bell+ under Z-dephasing (γ = 0.05, J = 1.0):

| CΨ | g(CΨ) | dCΨ/dt | dB/dt |
|----|-------|--------|-------|
| 0.319 | 14.3 | -0.126 | 0.478 |
| 0.273 | 4.1 | -0.104 | 0.211 |
| **0.250** | **3.36** | **-0.093** | **0.170** |
| 0.228 | 3.0 | -0.082 | 0.142 |
| 0.198 | 2.85 | -0.068 | 0.115 |

**g(CΨ) = 3.36 at CΨ = 1/4.** Finite. No singularity. The Bures
metric is smooth across the fold boundary. Both dCΨ/dt and dB/dt are
nonzero at the crossing: the trajectory passes through CΨ = 1/4 with
finite velocity in both CΨ-space and Bures-space.

The metric increases toward large CΨ (near the initial Bell state)
because the state is purer and Bures distance is more sensitive. It
plateaus toward small CΨ, as the state approaches its limit
½(|00⟩⟨00| + |11⟩⟨11|).

---

## Phase 2: θ as Coordinate

In θ coordinates: g̃(θ) = g(CΨ) × (dCΨ/dθ)² where
dCΨ/dθ = sin(θ)/(2cos³θ).

At θ = 0 (CΨ = 1/4): dCΨ/dθ = 0, so g̃(θ) → 0.

**θ SHRINKS the metric at the fold** (g̃ → 0), it does not regularize
a divergence (there is none). This makes θ a poor coordinate near the
fold: it maps a finite metric region to zero, losing resolution.

θ remains useful as a compass (angular distance from the boundary:
θ = 0 means "at 1/4", θ = π/4 means "at 1/2"). But it is not a
Riemannian normal coordinate. CΨ itself is the better metric
coordinate at the fold.

---

## Phase 3: Geodesic Analysis

Along the trajectory's own coordinate the geodesic equation
d²CΨ/ds² + Γ(dCΨ/ds)² = 0, with Γ = (1/2g) dg/dCΨ, holds for any
monotone curve: on a line, every monotone path is the shortest one, so
the equation cannot test anything. The question that can fail lives in
the state space: a path is a Bures geodesic exactly when its Bures
length equals the Bures angle arccos √F between its endpoints.

| Initial state | Ratio, t ≤ 0.2 | Ratio, t ≤ 1.5 |
|---|---|---|
| Bell+ | 1.000000 | 1.000000 |
| \|+0⟩ | 1.030 | 7.68 |
| a generic state | 1.062 | 3.43 |

(γ = 0.05, J = 1; Bures length over the endpoint Bures angle. The
angle is capped at π/2 while a path can keep growing, so the long
window's ratios grow with the window; the short window shows where the
departure starts.)

**Bell+ runs exactly along a Bures geodesic.** The Hamiltonian leaves
it alone, and the state decays along the line of mixtures of its two
Bell projectors, monotonically, so the path and the geodesic between
its endpoints are the same curve. That is a POSITIVE result for Bell+:
here decoherence takes the shortest path from the initial state
toward equilibrium. It is not the second law in general. A state the
Hamiltonian moves leaves the geodesic from the first moment, and the
longer it is followed the further its path runs beyond the distance it
ends up covering.

---

## Phase 4: The Path Coefficient's Shape

The coordinate-shape second derivative of the Bures path-metric
coefficient, S = -(1/2g) d²(ln g)/d(CΨ)²:

| CΨ | S |
|----|---|
| 0.33 | -438 |
| 0.30 | -62 |
| **0.25** | **-25** |
| 0.21 | -17 |
| 0.19 | -15 |

**S = -25 at the fold, finite.** It is not a curvature: a
one-dimensional path metric has none, and S changes under a change of
coordinate (the flat line written as g = 1/(4x) gives S = -2/x). What
it shows is the shape of g(CΨ): steep near the initial pure state,
flattening toward the state's limit, with nothing special at the fold,
which sits in the smooth middle of the curve.

---

## Phase 5: Second γ-Derivative

d²CΨ/dγ² at fixed observation time (t = 0.75):

**d²CΨ/dγ² = 9.8 at the γ where CΨ = 1/4.** Finite, and it could not
be otherwise: for Bell+ every entry of the state is analytic in γ (a
matrix exponential) and the coherence e^(−4γt)/2 never passes through
zero, so CΨ, which sums absolute values, is analytic too. This is not a fidelity
susceptibility (that would be the second derivative of the state
fidelity), and a finite value here does not by itself rule out
critical behaviour. That the fold is algebraic (discriminant = 0) and
not critical rests on the algebra, not on this number.

---

## Connection to Previous Results

**Algebra-first principle:** The Bures geometry confirms what six
prior analyses showed: the fold at CΨ = 1/4 is algebraic, not
geometric/thermodynamic/topological. The path coefficient is smooth
and finite at the fold. θ is a convenient label, not a fundamental
coordinate.

**Geodesic decoherence:** The one positive surprise, and a narrow one.
For Bell+, among all paths from the initial state to where it is at
time t, decoherence takes the shortest. That gives CΨ monotonicity
(the proven dCΨ/dt < 0) a geometric picture for this state. It does
not carry over: once the Hamiltonian moves the state, the path leaves
the geodesic.

**θ and fidelity (r = 0.87):** The correlation between θ and fidelity
([Theta Palindrome Echo](THETA_PALINDROME_ECHO.md)) is NOT because θ
is a natural Bures coordinate (it shrinks the metric at the fold).
The correlation likely arises because both θ and fidelity are monotone
functions of CΨ along this trajectory.

---

## What This Does NOT Establish

- Whether the geodesic property of H-dead states holds for N > 2 (only
  tested N=2)
- Whether different initial states give the same metric behavior
- Any curvature: the metric computed is one-dimensional
- Whether there exists a better coordinate than CΨ or θ
- Whether the metric structure extends to the full density matrix space
  (we only computed the 1D induced metric along the trajectory)

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [`simulations/information_geometry.py`](../simulations/information_geometry.py) |
| Output | [`simulations/results/information_geometry.txt`](../simulations/results/information_geometry.txt) |

---

## References

- Braunstein, S.L., Caves, C.M. (1994). "Statistical distance and the
  geometry of quantum states." PRL 72, 3439.
- Zanardi, P., Giorda, P., Cozzini, M. (2007). "Information-theoretic
  differential geometry of quantum phase transitions." PRL 99, 100603.
- Petz, D. (1996). "Monotone metrics on matrix spaces."
  Lin. Alg. Appl. 244, 81.

---

*θ is a compass, not a coordinate. The fold has no geometric
singularity. And where the Hamiltonian leaves the state alone, the path
through the fold is the shortest possible: decoherence is not just
inevitable there, it is efficient.*
