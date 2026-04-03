# Information Geometry: θ as Riemannian Coordinate

<!-- Keywords: information geometry Lindblad, Bures metric CΨ trajectory,
theta angular coordinate fold, geodesic decoherence shortest path,
Fisher susceptibility quantum, Gaussian curvature fold point,
R=CPsi2 information geometry -->

**Status:** Tier 2 (computed, N=2 Bell state)
**Date:** April 1, 2026
**Script:** [information_geometry.py](../simulations/information_geometry.py)
**Data:** [information_geometry.txt](../simulations/results/information_geometry.txt)
**Depends on:**
- [Boundary Navigation](BOUNDARY_NAVIGATION.md) (θ definition, formula 15)
- [CΨ Monotonicity Proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md) (dCΨ/dt < 0)
- [Entropy Production](ENTROPY_PRODUCTION.md) (algebra-first context)

---

## What this document is about

This document gives the angular parameter θ a proper geometric foundation
by computing the Bures Riemannian metric (the natural measure of
distinguishability between nearby quantum states) along the decoherence
trajectory. Key finding: the fold at CΨ = ¼ has no geometric singularity,
the metric is smooth, and the Lindblad trajectory is approximately
geodesic (the shortest path in state space). θ is a useful compass but
not a fundamental coordinate; CΨ itself is the better metric coordinate
at the fold.

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

3. **The Lindblad trajectory IS approximately geodesic.** Mean geodesic
   deviation 9.1 × 10⁻⁴. Decoherence follows the geometrically
   shortest path in the Bures metric. The second law is the shortest
   path to equilibrium.

4. **Curvature K = -25 at the fold.** Negative (hyperbolic), finite.
   No geometric singularity. The fold region has strong but bounded
   negative curvature.

5. **Fisher susceptibility χ_F = 9.8.** Finite. CΨ = 1/4 is a smooth
   crossing, not a dynamical phase transition.

---

## Phase 0: θ Inventory

θ = arctan(√(4CΨ-1)) appears in 15 repo files:

| Document | Usage | Status |
|----------|-------|--------|
| BOUNDARY_NAVIGATION | Definition (formula 15) | Origin |
| THETA_PALINDROME_ECHO | Correlates with fidelity r=0.87 | Computed |
| STRUCTURAL_CARTOGRAPHY | Lives on 3D manifold (98% in 3 PCs) | Computed |
| CIRCUIT_DIAGRAM | "Voltmeter" metaphor | Analogy |
| ANALYTICAL_FORMULAS | Formula 15 (angular distance) | Catalogued |
| TOPOLOGICAL_EDGE_MODES | Used as Berry phase parameter (φ = -0.77) | Computed |

Never computed: Riemannian metric, geodesics, curvature, Fisher info.

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
plateaus near the maximally mixed state (small CΨ).

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

The geodesic equation in the CΨ metric: d²CΨ/ds² + Γ(dCΨ/ds)² = 0,
where Γ = (1/2g) dg/dCΨ is the Christoffel symbol (the correction term that accounts for curvature when computing acceleration on a curved manifold).

**Mean geodesic deviation: 9.1 × 10⁻⁴.** The Lindblad trajectory is
approximately geodesic in the Bures metric. Decoherence follows the
geometrically shortest path from the initial state to equilibrium.

This is a POSITIVE result: the second law (CΨ monotonically
decreasing) corresponds to the geometric shortest path (geodesic) in
the Bures metric. Decoherence is not just monotone; it is efficient.
No "wasted" entropy production relative to the geometric optimum.

---

## Phase 4: Curvature

Gaussian curvature K = -(1/2g) d²(ln g)/d(CΨ)²:

| CΨ | K |
|----|---|
| 0.33 | -438 |
| 0.30 | -62 |
| **0.25** | **-25** |
| 0.21 | -17 |
| 0.19 | -15 |

**K = -25 at the fold.** Negative (hyperbolic), finite. The CΨ-space
has strong negative curvature near the fold: states diverge quickly
in this region (high sensitivity to perturbations). But the curvature
does not diverge: no geometric singularity.

The curvature is strongest at large CΨ (near the initial pure state)
and decays toward the maximally mixed state. The fold at 1/4 sits in
a region of moderate negative curvature, not at an extremum.

---

## Phase 5: Fisher Susceptibility

χ_F(γ) = d²CΨ/dγ² at fixed observation time (t = 0.75):

**χ_F = 9.8 at the γ where CΨ = 1/4.** Finite. CΨ = 1/4 is NOT a
dynamical phase transition in γ-space. The crossing is smooth: the
state passes through CΨ = 1/4 without any critical behavior in the
control parameter γ.

This distinguishes CΨ = 1/4 from quantum phase transitions (where
χ_F diverges at the critical point, signaling a qualitative change in the ground state). The fold is algebraic (discriminant
= 0), not critical (susceptibility finite).

---

## Connection to Previous Results

**Algebra-first principle:** The Bures geometry confirms what six
prior analyses showed: the fold at CΨ = 1/4 is algebraic, not
geometric/thermodynamic/topological. The metric is smooth, the
curvature is finite, the susceptibility doesn't diverge. θ is a
convenient label, not a fundamental coordinate.

**Geodesic decoherence:** The one positive surprise. The Lindblad
trajectory being approximately geodesic means: among all paths from
the initial state to equilibrium, decoherence chooses (approximately)
the shortest one. This gives CΨ monotonicity (the proven dCΨ/dt < 0)
a geometric interpretation: it is the gradient flow along the shortest
Bures path.

**θ and fidelity (r = 0.87):** The correlation between θ and fidelity
([Theta Palindrome Echo](THETA_PALINDROME_ECHO.md)) is NOT because θ
is a natural Bures coordinate (it shrinks the metric at the fold).
The correlation likely arises because both θ and fidelity are monotone
functions of CΨ, which itself is approximately geodesic.

---

## What This Does NOT Establish

- Whether the geodesic property holds for N > 2 (only tested N=2)
- Whether different initial states give the same metric behavior
- Whether the negative curvature has physical consequences
- Whether there exists a better coordinate than CΨ or θ
- Whether the metric structure extends to the full density matrix space
  (we only computed the 1D induced metric along the trajectory)

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | simulations/information_geometry.py |
| Output | simulations/results/information_geometry.txt |

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
singularity. But the path through the fold is the shortest possible.
Decoherence is not just inevitable; it is efficient.*
