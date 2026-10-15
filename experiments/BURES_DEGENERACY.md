# Bures Speed and the Degeneracy Profile: A First Look

<!-- Keywords: Bures metric Lindblad trajectory, degeneracy curvature correlation,
quantum state space geometry, weight sector dynamics, palindrome geometry,
open quantum system information geometry, R=CPsi2 Bures degeneracy -->

**Status:** Null result for the trajectories tested: the |+⟩^⊗N speed is a closed form that the degeneracy profile does not enter; for random pure and mixed states the Hamiltonian moves, the speed is carried by the Hamiltonian part, and the chain with more centre degeneracy is not faster than the star with less. Whether the profile shapes the speed in any other way is open.
**Date:** April 3, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Degeneracy Palindrome](DEGENERACY_PALINDROME.md),
[Information Geometry](INFORMATION_GEOMETRY.md)
**Verification:** [`simulations/qfi_degeneracy_correlation.py`](../simulations/qfi_degeneracy_correlation.py)

---

## What this means

The palindrome describes the landscape: how wide the valley is at each
point. The question here was whether a wide valley, many eigenvalues on
one decay-rate shell, means faster motion through quantum state space
when that shell is active. The motion itself is directed: the state
flows from the middle and high weights toward low weight (silence).
The palindrome is the potential, symmetric and timeless. The trajectory
is the kinetics, directed and asymmetric.

For the state this document follows, |+⟩^⊗N on the Heisenberg chain under
equal Z-dephasing, the answer is that the landscape cannot matter. The
Hamiltonian never touches that state ([K-Dosimetry](K_DOSIMETRY.md)): it
stays a product of identical one-site states, and its speed through state
space is a one-line formula in γt alone. It is fastest at the start and
slows down monotonically. The 152 eigenvalues at the centre of the N = 4
palindrome are there, and the centre shell even carries the most weight at
the start, but the speed does not depend on them.

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
- Initial state: |+⟩^⊗N (product superposition); at N = 4 also GHZ, W and a random pure state
- Propagation by the exact propagator expm(L·t)
- Speed: the Bures-angle speed v = √F_Q / 2, with the quantum Fisher
  information F_Q = Σ_{p_i+p_j>0} 2|⟨i|ρ̇|j⟩|²/(p_i + p_j) in the
  eigenbasis of ρ(t), population terms included, ρ̇ = L(ρ)
- Mode weight profile: the fraction of ρ(t) in each decay-rate shell k,
  through the shell projectors
- Time grid: 148 points, logarithmic from t = 0.01 to 1, then linear up to t = 200

---

## Result: Weight sector dynamics

The mode weight profile w_k(t) tracks which decay-rate shell carries
ρ(t) at each time. For |+⟩^⊗N:

**Early times (t ≈ 0):** The profile is binomial, [1, N, …, 1]/2^N
([0.063, 0.250, 0.375, 0.250, 0.062] at N = 4): each of the N sites carries
its X letter or its I letter with equal weight, so the middle shells hold
the most.

**Intermediate times:** Weight-1 modes dominate as the higher modes
have decayed.

**Late times (t ≫ 1/γ):** Only weight-0 modes survive (the N+1 steady-
state modes). The state approaches the maximally mixed state.

This cascade from the middle shells down to weight 0 is universal across
N = 2, ..., 5.

---

## Result: the speed is a closed form

ρ(t) = ⊗ (I + g·X)/2 with g = e^(−2γt), and the eigenbasis of each one-site
factor never moves; only its two populations (1 ± g)/2 change. The quantum
Fisher information is then the classical one of those populations, added
over the N sites:

  F_Q = N · 4γ²g² / (1 − g²),  v(t) = ½·√F_Q.

The propagated trajectory meets this formula at every grid point (largest
relative distance of order 10⁻¹¹ over N = 2–5). The speed diverges like 1/√t at
t → 0, where the pure state starts to mix, and falls monotonically after
that. J does not appear. The trajectory runs on the N + 1 symmetric sums of the
X strings, eigenvectors of L at the shell rates 2kγ, one mode each, so of the
spectrum only those rates enter and no degeneracy does.

So the largest sampled speed is always the first grid point, t = 0.01:
1.580, 1.936, 2.235, 2.499 at N = 2, 3, 4, 5. A correlation between d_total(k)
and the largest speed seen while shell k dominates reads only which shell
carries the most weight at that first point, the binomial middle, and says
nothing about the spectrum. It comes out large at even N (0.999 at N = 2,
0.994 at N = 4), where the binomial middle and the degeneracy centre are the
same shell, and small at odd N (0.56 at N = 3 and 5), where the two middle shells hold equal weight and the first-point maximum lands on only one of them.

---

## Other initial states (N = 4)

| Initial state | Peak v (at t = 0.01) | Shell dominating at t = 0.01 |
|---|---|---|
| \|+⟩^4 | 2.235 | k = 2 |
| W | 1.935 | k = 2 |
| Random | 3.920 | k = 2 |
| GHZ | 2.232 | k = 0 |

Every state starts pure, so every speed has its largest value at the first
grid point. The GHZ state, (|0000⟩ + |1111⟩)/√2, never leaves
span{|0000⟩, |1111⟩}: its eigenvectors (|0000⟩ ± |1111⟩)/√2 stay fixed and
only their weights (1 ± e^(−2Nγt))/2 change, so its whole speed is the
population term; its weight sits on the shells k = 0 (the populations) and
k = 4 (the one coherence), never on k = 2.

---

## Trajectories the Hamiltonian moves

A single trajectory's speed is fixed by ρ(t) and ρ̇ = L(ρ): by the
eigenvalues the trajectory populates and their amplitudes. How many other
modes sit on the same grid line cannot enter it, so the degeneracy profile can
shape the speed only through which states it is asked of. Of the states
above only the random one is moved by the Hamiltonian (W, like |+⟩^⊗N, is
a state the Hamiltonian never touches under equal dephasing), and it starts
pure, so its speed too peaks at t → 0.

Random states at later times, pure and mixed, give a first look. The N = 4
chain and star both have three bonds, so equal Frobenius norm ‖H‖ = 12
(the chain's spectral radius is the larger, 6.46 against 5.0), but different
centre degeneracies, d_total(2) = 152 on the chain and 126 on the star. Over 60
random pure and 60 random rank-4 states:

| t | chain | star | star faster |
|---|---|---|---|
| 0.5 | 2.646 | 2.678 | 58 of 120 |
| 2 | 2.199 | 2.235 | 60 of 120 |
| 5 | 1.502 | 1.556 | 73 of 120 |
| 10 | 0.780 | 0.834 | 85 of 120 |
| 20 | 0.212 | 0.242 | 86 of 120 |
| 40 | 0.024 | 0.027 | 76 of 120 |

The more degenerate chain is not faster, although its spectral radius is
the larger: the star's mean is ahead at every time, by about 1% early on,
where it wins only half the states, and by up to 14% later. The speed is carried by the Hamiltonian part: computed
with L_D(ρ) alone in place of ρ̇ it is 3 to 14% of the total, and the same
for both graphs. This is in line with the simulated Bures-velocity check of EQ-030 in
[`review/EMERGING_QUESTIONS.md`](../review/EMERGING_QUESTIONS.md), in a
different setting (N = 3, six soft Hamiltonians under T1), where the Bures
speed at the cusp was uncorrelated with the count of Π-protected observables
and its fast cases were read as fast oscillation. One comparison of two graphs is a sighting, not
a law. A sharper form of the question follows the
[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): d(k) counts
the modes exactly on the grid lines Re λ = −2kγ, while a trajectory the
Hamiltonian moves runs largely on modes between them, so the split of ρ(t)
between on-grid and off-grid modes, not d(k), is what could matter. That is
open.

---

## Reproduction

- QFI speed, shell weights and the closed-form check: [`simulations/qfi_degeneracy_correlation.py`](../simulations/qfi_degeneracy_correlation.py)
- Output: [`simulations/results/qfi_degeneracy_correlation.txt`](../simulations/results/qfi_degeneracy_correlation.txt)
- [`simulations/bures_degeneracy_correlation.py`](../simulations/bures_degeneracy_correlation.py), a finite-difference Bures estimate, pairs separately computed left and right eigenvectors by column, which is wrong inside the degenerate eigenspaces present from N = 3 on; no number here comes from it.
