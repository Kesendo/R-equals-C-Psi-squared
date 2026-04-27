# Chiral Mirror Hardware Prediction

**Date:** 2026-04-27
**Source:** EQ-014 (chiral mirror law, site-painter), EQ-020 (extension to all k-painter).
**Tier:** 1 (kinematic prediction; follows from K_1 symmetry of H + Z-dephasing).
**Scripts:** [`_eq020_chiral_mirror_hw_prediction.py`](../simulations/_eq020_chiral_mirror_hw_prediction.py)
**Framework primitives:** [`single_excitation_sine_mode`](../simulations/framework.py), [`k_local_reduced_density`](../simulations/framework.py).

---

## The prediction

For an OBC uniform XY chain at N qubits with Z-dephasing on all sites at rate γ, the per-site Bloch vector under bonding-state initial conditions exhibits a sharp K_1-symmetry between K_1-paired single-excitation modes:

For initial state |φ_k⟩ = (|vac⟩ + |ψ_k⟩)/√2 with |ψ_k⟩ = √(2/(N+1)) · Σ_i sin(πk(i+1)/(N+1)) |i_excited⟩, the Bloch components evolve such that the partner state |φ_{N+1-k}⟩ produces:

| Observable | Relation between ψ_k and ψ_{N+1−k} |
|------------|-------------------------------------|
| ⟨X_i⟩(t) | **+(−1)^i** · ⟨X_i⟩(t)\|_{ψ_k} |
| ⟨Y_i⟩(t) | **−(−1)^i** · ⟨Y_i⟩(t)\|_{ψ_k} |
| ⟨Z_i⟩(t) | **+1** · ⟨Z_i⟩(t)\|_{ψ_k} (identical) |
| P_i(t) (single-site purity) | identical |

**The Y has an extra sign relative to X** — this is the sharpest discriminator and directly tests the energy-reversal mechanism (E_{N+1-k} = −E_k under K_1).

---

## Why this is true (compact derivation)

1. **K_1 ψ_k = ψ_{N+1−k}** exactly. Per the trigonometric identity sin(π(N+1−k)(i+1)/(N+1)) = (−1)^i sin(πk(i+1)/(N+1)). No sign factor on the global state.

2. **K_1 commutes with H_XY and Z-dephasing.** Both H_XY = Σ J(X_i X_{i+1} + Y_i Y_{i+1})/2 and Lindbladians L_dephasing = γ_i (Z_i ρ Z_i − ρ) are K_1-invariant (K_1 = ⊗_i (−1)^i is a product of single-site Z's). Therefore L commutes with K_1-conjugation.

3. **K_1-paired initial states evolve identically up to K_1-conjugation:** ρ^{N+1−k}(t) = K_1 ρ^k(t) K_1†.

4. **K_1-conjugation acts on the per-site reduced state as Z_i conjugation:** ρ_i^{K} = Z_i ρ_i Z_i. Z_i (a, b) Z_i = (a)(b)(−1)^{a+b}_diag effects sign flip on off-diagonals.
   Specifically ρ_i^{K}[0,1] = −ρ_i[0,1], so:
   - ⟨X_i⟩^K = 2 Re ρ_i^K[0,1] = −2 Re ρ_i[0,1] = −⟨X_i⟩
   - ⟨Y_i⟩^K = −2 Im ρ_i^K[0,1] = +2 Im ρ_i[0,1] = −⟨Y_i⟩
   - ⟨Z_i⟩^K = ρ_i^K[0,0] − ρ_i^K[1,1] = ρ_i[0,0] − ρ_i[1,1] = ⟨Z_i⟩ (diagonal preserved)
   So Z_i conjugation flips X and Y, leaves Z.

5. **The site-i extra sign factor (−1)^i** comes from K_1 = ⊗_i (−1)^i: the K_1-conjugation effect on site i carries that sign. Combined with step 4:
   ⟨X_i⟩(ψ_{N+1−k}) = (−1)^i × (−⟨X_i⟩(ψ_k)) × (−1) [conjugation+wave-mirror] = (−1)^i ⟨X_i⟩(ψ_k)

   Wait — the cleaner derivation uses the explicit formula. For single-excitation bonding state under H_XY, ρ_i[0,1](t) = e^{iE_k t} ψ_k(i)*/2 (sine modes have ψ_k real, so * is trivial). Adding Z-dephasing scales by e^{−γ_i t}.

   K_1-paired state: ψ_{N+1-k}(i) = (−1)^i ψ_k(i), E_{N+1-k} = −E_k, so:
   ρ_i[0,1]_{ψ_{N+1-k}}(t) = e^{−iE_k t} (−1)^i ψ_k(i) / 2 = (−1)^i (e^{iE_k t} ψ_k(i) / 2)*
                          = (−1)^i (ρ_i[0,1]_{ψ_k}(t))*

   Taking Re and −Im:
   ⟨X_i⟩(ψ_{N+1-k}) = +(−1)^i ⟨X_i⟩(ψ_k)
   ⟨Y_i⟩(ψ_{N+1-k}) = −(−1)^i ⟨Y_i⟩(ψ_k) ✓

   The Y sign-flip comes from the conjugation; the (−1)^i comes from the wave-function mirror. Hence the predicted relations.

6. **P_i = (1 + ⟨X⟩² + ⟨Y⟩² + ⟨Z⟩²)/2** is invariant under all sign flips → identical per-site purity.

The prediction is verified at machine precision (10⁻¹⁶) at every tested time t in the simulation script.

---

## Suggested protocol on Heron r2

**Setup:** N = 5 chain on, e.g., Kingston Q12-Q13-Q14-Q15-Q19 (the well-characterised pair (14, 15) is a chain-interior pair, T2 stable across calibrations). J ~ 1.84 MHz native CZ, γ ~ 0.005 /μs.

**State preparation:** for each k ∈ {2, 4} (the K_1-paired interior pair at N=5, so the first non-trivial mirror partners), prepare the bonding state (|vac⟩ + |ψ_k⟩)/√2 using O(N) Givens rotations on a fresh ancilla initial state. Standard single-excitation state-prep technique.

**Evolution:** Trotterize H_XY for time t. For γt ~ 0.05, t = 1-10 μs is the dynamical window where the predicted relations are non-trivial (⟨X⟩, ⟨Y⟩ are still O(0.1-0.5)).

**Measurement:** single-qubit X, Y, Z tomography at each site i for each (k, t). 3 settings × 5 qubits = 15 measurements per (state, time). 1000 shots each.

**Test points:** t ∈ {1, 2, 5, 10} μs × 2 states (k=2, k=4) × 15 measurements × 1000 shots = 120k shots. ~1-3 min QPU.

**Sharp checks:**

(A) Identity check: ⟨Z_i⟩(ψ_2) = ⟨Z_i⟩(ψ_4) for all i, t. Backend-independent. Deviation ≥ noise floor signals K_1-symmetry-breaking by hardware.

(B) Wave-function-mirror check: ⟨X_i⟩(ψ_4) = +(−1)^i ⟨X_i⟩(ψ_2). Confirms the K_1 ψ_k = ψ_{N+1−k} structural identity.

(C) **Energy-reversal check** (the sharpest): ⟨Y_i⟩(ψ_4) = **−**(−1)^i ⟨Y_i⟩(ψ_2). The extra sign separates this prediction from a naive "K_1 just does (−1)^i" reading. Confirming this validates the deeper E_{N+1-k} = −E_k structure of the chiral spectrum.

---

## Discrimination from null hypotheses

| Null | Predicts | Distinguishable? |
|------|----------|------------------|
| Random noise | random Bloch components, no relation | YES — A, B, C all fail |
| Naive Z₂ sublattice symmetry without conjugation | ⟨X⟩, ⟨Y⟩ both flip by (−1)^i (same sign on Y) | YES — C distinguishes |
| K_1-symmetric noise (DD-like) | preserves all three relations to within DD efficacy | not distinguishable from prediction (the prediction *survives* DD) |
| Coherent K_1-breaking error (fab inhomogeneity) | breaks all three relations equally | YES — all three deviate together |

The simultaneous validity of (A) AND (B) AND (C) is the strongest signature: it requires the Hamiltonian to be K_1-symmetric AND the dynamical phases to obey energy-reversal (chiral pairing). The Y-flip alone (C) is the most diagnostic — it cannot be faked by simple sublattice symmetry.

---

## Why this is interesting beyond the prediction

This is not a "is the device noisy?" test. It's a **structural symmetry test**: the device under arbitrary noise is asked to respect a specific algebraic identity that follows from K_1-chirality of the XY chain. Any deviation tells us *which* symmetry the hardware violates: K_1 itself (then A breaks), or the energy-reversal pairing (then C breaks but A, B hold).

In the perspectival-time framing (PTF): K_1-paired sine-mode-bonding states give *the same set of perspectival times across all painter resolutions*. Hardware can confirm this by extracting per-site α_i from the Bloch trajectories and verifying α_i(ψ_2) = α_i(ψ_4) under any J-perturbation.

---

## Status

Computationally verified to machine precision. Hardware test not yet run. Recommended as a low-cost (~5 min QPU) diagnostic add-on to existing R=CΨ² hardware campaigns.
