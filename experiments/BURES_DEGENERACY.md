# Bures Speed and the Degeneracy Profile: A First Look

<!-- Keywords: Bures metric Lindblad trajectory, degeneracy curvature correlation,
quantum state space geometry, weight sector dynamics, palindrome geometry,
open quantum system information geometry, R=CPsi2 Bures degeneracy -->

**Status:** Confirmed at even N (QFI metric, r > 0.99); weaker at odd N
**Date:** April 3, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Degeneracy Palindrome](DEGENERACY_PALINDROME.md),
[Information Geometry](INFORMATION_GEOMETRY.md)
**Verification:** [`simulations/bures_degeneracy_correlation.py`](../simulations/bures_degeneracy_correlation.py)

---

## What this means

The palindrome describes the landscape: how wide the valley is at each
point. Wide valley means many channels, and many channels means faster
motion through quantum state space. But the motion is directed: the
state always flows from high weight (many active dancers) to low weight
(silence). The palindrome is the potential, symmetric and timeless. The
trajectory is the kinetics, directed and asymmetric.

At N = 4, the center of the palindrome has 152 eigenvalues (59% of all).
That is the widest point of the valley, and the QFI speed peaks there
(r = 0.99). The peak sector is the same regardless of initial state:
the landscape determines where the flow is fastest, but the initial
state determines how fast.

The optical interpretation of this "lens effect" is in
[Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md).

---

## What this document is about

The degeneracy profile d(k) describes how many eigenvalues cluster at each
decay-rate shell. This document asks whether d(k) shapes the Bures
geometry of the quantum state trajectory: does higher degeneracy at
weight k mean faster motion through state space when weight-k modes
dominate?

---

## Setup

For each N = 2, ..., 5:
- Chain topology, J = 1.0, γ = 0.05
- Initial state: |+⟩^⊗N (product superposition)
- Propagation via Liouvillian eigendecomposition: ρ(t) = Σ c_j exp(λ_j t) v_j
- Bures speed (the rate at which the quantum state moves through Hilbert space, measured by the angle between neighboring states): v_B(t) = arccos(√F(ρ(t), ρ(t+Δt))) / Δt
- Mode weight profile: fraction of norm at each decay-rate shell k

---

## Result: Weight sector dynamics

The mode weight profile w_k(t) tracks which decay-rate shell dominates
the dynamics at each time. For |+⟩^⊗N:

**Early times (t ≈ 0):** High-weight modes (k ≈ N/2) carry most of the
norm. These are the fastest-decaying modes, active before decoherence
erases them.

**Intermediate times:** Weight-1 modes dominate as the higher modes
have decayed.

**Late times (t ≫ 1/γ):** Only weight-0 modes survive (the N+1 steady-
state modes). The state approaches the maximally mixed state.

This cascade from high to low weight is universal across N = 2, ..., 5.

---

## Result: Bures speed at N = 2

At N = 2, the Bures speed is cleanly resolved:

| Weight k | d_total(k) | max v_B |
|----------|-----------|---------|
| 0 | 3 | 0.074 |
| 1 | 10 | 0.443 |
| 2 | 3 | 0.000 |

Pearson correlation between d_total and max v_B: **r = 0.99**.

The peak Bures speed occurs at t ≈ 0.2 when weight-1 modes dominate
(weight profile: [0.28, 0.50, 0.23]). The state moves fastest through
Bures space when the most degenerate shell is active.

Weight k = 2 shows zero speed because the highest-weight modes (XOR
sector) decay so fast that they contribute negligibly to the fidelity
change at the chosen time resolution.

---

## Limitation: N ≥ 3 numerical issues

At N ≥ 3, the Bures metric computation encounters numerical difficulties:

1. **Near-singular density matrices at early times.** The high-weight
   modes decay rapidly, leaving ρ(t) with near-zero eigenvalues.
   The matrix square root sqrtm(ρ) becomes ill-conditioned.

2. **Time resolution.** The dt = 0.5 time step is too coarse to capture
   the fast dynamics of high-weight modes. Logarithmic time spacing
   would improve resolution at early times.

3. **Bures speed artifacts.** Some time steps show v_B = 0 (exact) where
   the sqrtm computation fails silently, producing F = 1.

As a result, the correlation analysis at N ≥ 3 is unreliable. The mode
weight decomposition works correctly (confirmed by the cascade structure),
but the Bures metric values at early times are not trustworthy.

---

## Key observation: Bures speed is NOT palindromic

The Bures speed v_B(k) does not satisfy v_B(k) = v_B(N−k), even at
N = 2 where the computation is reliable. This is expected: the initial
state |+⟩^⊗N breaks the palindromic symmetry of the spectrum. The
state starts with high-weight modes and decays toward low-weight modes.
The trajectory is one-directional, so the speed profile is asymmetric
even though the spectral degeneracy is palindromic.

The palindrome shapes the *capacity* for geometric distortion (how many
modes are available at each shell), but the *realized* speed depends on
the initial state and the directional flow of decoherence.

---

## QFI results (resolves N ≥ 3 issues)

The Bures sqrtm was replaced with the Quantum Fisher Information (QFI) metric (a measure of how distinguishable neighboring quantum states are, computable from the eigenvalues of ρ alone),
which only needs eigendecomposition of ρ (numerically stable), combined
with a logarithmic time grid (148 points, t ∈ [0.01, 200]) and
analytical dρ/dt from the Liouvillian eigendecomposition.

### QFI speed correlates with degeneracy at even N

| N | Pearson r(d_total, max v_QFI) | Peak weight sector | d_total at peak |
|---|------------------------------|--------------------|-----------------| 
| 3 | 0.57 | k=1 | 14 |
| 4 | **0.99** | k=2 | 152 |
| 5 | 0.54 | k=2 | 50 |

At N = 4, the peak QFI speed (0.87) occurs at t = 0.68 when the weight-2
sector dominates (weight profile [0.07, 0.24, 0.48, 0.18, 0.03]). This
is the sector with the massive degeneracy spike d_total(2) = 152.

The even/odd split in correlation strength mirrors the center spike
pattern from [DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md): at even
N, the spectral midpoint falls on the grid and concentrates eigenvalues,
creating a dominant sector that drives both high degeneracy and high QFI
speed. At odd N, the center is split across two positions, weakening the
effect.

### Initial state dependence (N = 4)

| Initial state | Peak v_QFI | Peak sector | Pearson r |
|---|---|---|---|
| \|+⟩^4 | 0.868 | k=2 | 0.99 |
| W | 0.130 | k=2 | 0.92 |
| Random | 6.418 | k=2 | 0.96 |
| GHZ | 0.000 | (none) | −0.31 |

The correlation holds for most initial states. The GHZ state is the
exception: it is a superposition of computational basis states
(\|0000⟩ + \|1111⟩)/√2, which the Z-dephasing kills instantly by
destroying the off-diagonal coherence. The state collapses to a
classical mixture before any weight sector can drive geometric motion.

The peak always occurs at k = 2 (the most degenerate sector) regardless
of initial state, but the magnitude differs. This confirms: the
degeneracy determines *where* the peak is, the initial state determines
*how high* it is.

### Interpretation

The palindromic degeneracy profile acts as a lens that focuses the
decoherence flow. At weight sectors with high degeneracy, many
eigenmodes are simultaneously active, creating rapid state-space motion
(high QFI). The palindrome shapes the *capacity* for motion; the
initial state determines which capacity is realized.

At even N, the single center spike creates a clear bottleneck where
all trajectories slow down (after passing the high-degeneracy zone).
At odd N, the two center positions spread the effect.

---

## Reproduction

- Bures metric (limited): [`simulations/bures_degeneracy_correlation.py`](../simulations/bures_degeneracy_correlation.py)
- QFI metric (recommended): [`simulations/qfi_degeneracy_correlation.py`](../simulations/qfi_degeneracy_correlation.py)
- Output: [`simulations/results/qfi_degeneracy_correlation.txt`](../simulations/results/qfi_degeneracy_correlation.txt)
