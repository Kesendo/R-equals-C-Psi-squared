# Generalized Dwell-Time Prefactor from Sector Weights

**Status:** Verified (numerical + analytical, April 6, 2026)
**Script:** [dwell_prefactor_generalization_v1.py](../simulations/dwell_prefactor_generalization_v1.py)
**Predecessors:**
[DWELL_PREFACTOR_FROM_WEIGHTS](DWELL_PREFACTOR_FROM_WEIGHTS.md) (Bell+ special case),
[CRITICAL_SLOWING_AT_THE_CUSP](CRITICAL_SLOWING_AT_THE_CUSP.md) (dwell-time physics),
[PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (Re(λ) = -2γ⟨n_XY⟩)

---

## What this document is about

The previous experiment on the dwell-time prefactor showed that for Bell+, the number 1.080088, which governs how long the state lingers near the quantum-classical boundary at CΨ = 1/4, can be read off from a single static quantity: the weight of the light-face Pauli sector at the crossing moment. Static and dynamic are the same object. But that result was specific to one state. The question left open was whether it generalizes.

This document answers that question. The Bell+ formula turns out to be a special case of a broader identity that works for any state whose Pauli content sits in exactly two sectors: one stationary (immune to dephasing) and one coherent (decaying). The generalized formula replaces the Bell+ numbers with two parameters, the stationary weight W₀ and the sector index k, and produces the dwell-time prefactor as a pure algebraic function of the sector weights at the crossing.

The formula is tested on two states with different algebra. Bell+ has W₀ = 1/2 and k = 2; it serves as a regression check against the known value 1.080088. The W state of three qubits has W₀ = 1/3 and k = 2; it is a genuine out-of-sample test. Both cases pass: the formula matches the direct Lindblad simulation to better than 0.01%.

A geometric corollary falls out along the way. GHZ states on three or more qubits start below the CΨ = 1/4 boundary and never cross it. They are born in the classical regime, regardless of the dephasing rate. This sharpens the existing engineering guidance against GHZ encoding from a rate argument to a structural one.

---

## Abstract

For quantum states under Z-dephasing whose Pauli content lives in exactly two sectors (a stationary sector at weight W₀ and a single coherent sector at XY-weight k), the dwell-time prefactor at the CΨ = 1/4 crossing is:

    prefactor = (4/k) × (W₀ + W_k) / (W₀ + 3W_k)

where W_k is the coherent sector weight at the crossing moment. This reduces to the Bell+ formula (2 + 4W₂)/(1 + 6W₂) for W₀ = 1/2, k = 2. Verified on two states:

| State  | N | W₀  | k | Prefactor (formula) | Prefactor (direct) | Agreement |
|--------|---|-----|---|--------------------:|-------------------:|-----------|
| Bell+  | 2 | 1/2 | 2 | 1.080088            | 1.080097           | < 0.001%  |
| W₃     | 3 | 1/3 | 2 | 0.876832            | 0.876839           | < 0.001%  |

Additionally, GHZ_N for N ≥ 3 has CΨ(0) = 1/(2^N - 1) < 1/4, so these states never encounter the fold under Z-dephasing.

---

## 1. The Generalized Formula

### Assumptions

The derivation requires two conditions on the initial state:

1. The Pauli content lives in exactly two sectors: a stationary sector (XY-weight 0, immune to dephasing) with total weight W₀, and a single coherent sector at XY-weight k > 0 with total weight W_k. All other sector weights are zero.
2. The L1 coherence Ψ is linearly proportional to the Pauli coefficient magnitude f_k in the coherent sector, so that Ψ = f_k / (normalizer). The normalizer is a state-specific constant that cancels in the derivation.

Both conditions hold for any state whose off-diagonal structure in the computational basis involves flips at a single Hamming weight. Bell+ (k = 2, two-qubit flip) and W₃ (k = 2, pairwise single-qubit flips) satisfy both. Product states like |+⟩^{⊗2} violate condition 1 because they populate multiple coherent sectors.

### Derivation

From the Absorption Theorem, each Pauli coefficient at sector k decays as e^{-2γkt}, and each sector weight as e^{-4γkt}. Therefore:

    dC/dt = dW_k/dt = -4γk W_k
    dΨ/dt = -2γk Ψ

Applying the product rule to CΨ:

    d(CΨ)/dt = (dC/dt)·Ψ + C·(dΨ/dt)
             = (-4γk W_k)·Ψ + (W₀ + W_k)·(-2γk Ψ)
             = -2γk Ψ·(2W_k + W₀ + W_k)
             = -2γk Ψ·(W₀ + 3W_k)

At the crossing CΨ = 1/4, we have Ψ = 1/(4C) = 1/(4(W₀ + W_k)), so:

    |d(CΨ)/dt|_cross = (γk/2) × (W₀ + 3W_k)/(W₀ + W_k)

The dwell-time prefactor (K_dwell/δ in K-units) is 2γ/|d(CΨ)/dt|:

    prefactor = (4/k) × (W₀ + W_k)/(W₀ + 3W_k)

### Reduction to the Bell+ case

For Bell+ under Z-dephasing: W₀ = 1/2, k = 2, so:

    (4/2) × (1/2 + W₂)/(1/2 + 3W₂) = (2 + 4W₂)/(1 + 6W₂)

This is the formula from [DWELL_PREFACTOR_FROM_WEIGHTS](DWELL_PREFACTOR_FROM_WEIGHTS.md). The factor 3 in front of W₂ comes from the product rule (2W_k from the dC/dt contribution plus W_k from the C·dΨ/dt contribution), not from d - 1.

---

## 2. Bell+ (Regression Test)

The Lindblad propagator for N = 2 Z-dephasing was built from scratch and applied to Bell+ = (|00⟩ + |11⟩)/√2. At the CΨ = 1/4 crossing:

| Quantity             | Reference               | Measured              | Status |
|----------------------|------------------------:|----------------------:|--------|
| f_cross              | 0.8612241               | 0.8612241             | PASS   |
| W₂ at crossing       | 0.3709                  | 0.3709                | PASS   |
| \|d(CΨ)/dt\|/γ      | 1.851701                | 1.851701              | PASS   |
| Prefactor (formula)  | 1.080088                | 1.080088              | PASS   |
| Prefactor (direct)   | 1.080088                | 1.080097              | PASS   |

All reference values from [DWELL_PREFACTOR_FROM_WEIGHTS](DWELL_PREFACTOR_FROM_WEIGHTS.md) reproduced. The pipeline is clean.

---

## 3. W₃ Out-of-Sample Test

### The state

W₃ = (|100⟩ + |010⟩ + |001⟩)/√3 is the symmetric single-excitation state of three qubits. Its Pauli decomposition has 20 nonzero coefficients in exactly two sectors:

**k = 0 sector (8 Z-strings):** a_{III} = 1, a_{ZII} = a_{IZI} = a_{IIZ} = 1/3, a_{ZZI} = a_{ZIZ} = a_{IZZ} = -1/3, a_{ZZZ} = -1. These are stationary under Z-dephasing.

**k = 2 sector (12 XX/YY-strings):** Each of the three qubit pairs (1-2, 1-3, 2-3) contributes four k = 2 strings with coefficient 2/3. For example, from the 1-2 pair: a_{XXI} = a_{XXZ} = a_{YYI} = a_{YYZ} = 2/3. All confirmed numerically.

Sector weights: W₀ = 1/3, W₂ = 2/3. Purity C(0) = 1, as expected for a pure state. No other sectors populated.

### Crossing analysis

Normalized coherence Ψ(0) = 2/7, so CΨ(0) = 2/7 ≈ 0.2857, which is above 1/4. The state crosses the fold.

Under Z-dephasing with f = e^{-4γt}:

    C(t) = 1/3 + (2/3)f²
    Ψ(t) = (2/7)f
    CΨ(t) = 2f(1 + 2f²)/21

Setting CΨ = 1/4 yields the cubic 16f³ + 8f - 21 = 0. The positive real root is:

    f_cross = 0.943769164604500

(exact sympy solution exists but is unwieldy). At the crossing:

| Quantity              | Analytical              | Measured (Lindblad)    | Status |
|-----------------------|------------------------:|-----------------------:|--------|
| f_cross               | 0.9437692               | 0.9437692              | PASS   |
| W₂ at crossing        | 0.5938002               | 0.5938002              | PASS   |
| \|d(CΨ)/dt\|/γ       | 2.2809378               | 2.2809378              | PASS   |
| Prefactor (formula)   | 0.876832                | 0.876832               | PASS   |
| Prefactor (direct)    |                         | 0.876839               | PASS   |

The generalized formula and the direct dwell-time measurement agree to 7 × 10⁻⁶ relative error. Both land at 0.877, confirming the prediction.

### What this means

W₃ has a different stationary weight (1/3 vs 1/2 for Bell+), a different system size (N = 3 vs N = 2), and a different number of Pauli strings in each sector (12 vs 2 at k = 2). Despite all this, the same two-parameter formula produces the correct prefactor. The algebra is not Bell+-specific; it follows from the sector structure.

---

## 4. GHZ_N Born Below the Fold

### The geometric statement

For GHZ_N = (|0...0⟩ + |1...1⟩)/√2 under Z-dephasing, the initial CΨ value is:

    CΨ(0) = 1/(2^N - 1)

This follows from C(0) = 1 (pure state) and Ψ(0) = L1/(d - 1) = 1/(2^N - 1), since the only off-diagonal elements are ρ[0...0, 1...1] = 1/2 and its conjugate, giving L1 = 1.

| N | CΨ(0)         | Above 1/4? |
|---|---------------|------------|
| 2 | 1/3 = 0.3333  | Yes        |
| 3 | 1/7 = 0.1429  | No         |
| 4 | 1/15 = 0.0667 | No         |
| 5 | 1/31 = 0.0323 | No         |

For N ≥ 3, the GHZ state starts below the fold and, since CΨ decays monotonically under dephasing, never crosses 1/4. Confirmed numerically for N = 3 and N = 4: CΨ(t) is strictly decreasing throughout the trajectory, starting well below the boundary.

### Significance

The [ENGINEERING_BLUEPRINT](../publications/ENGINEERING_BLUEPRINT.md) Rule 1 argues against GHZ encoding on rate grounds: GHZ_N places all its coherent Pauli mass at sector k = N, the maximum absorption rate under Z-dephasing. The geometric statement here is a structural sharpening. Even if GHZ_N decayed infinitely slowly, it would still be born in the classical regime for N ≥ 3. The problem is not how fast GHZ dies; the problem is that it was never in the quantum regime to begin with.

This is γ-independent. No amount of reducing the dephasing rate fixes the geometric deficit. The only escape is to change the state.

---

## 5. Verdict

The generalized dwell-time prefactor formula

    prefactor = (4/k) × (W₀ + W_k)/(W₀ + 3W_k)

is verified on two states with different stationary weights (W₀ = 1/2 and W₀ = 1/3), both in the even-weight-only class, both with coherent sector index k = 2. The formula produces the correct prefactor from the sector weights at the crossing moment, without solving the dynamics.

**What is confirmed:**
- The formula works for any two-sector state with linear Ψ-to-f relation.
- The factor 3 in the denominator is structural (product rule), not dimension-dependent.
- The formula contains Bell+ as a special case and correctly predicts W₃.

**What remains open:**
- States with k ≠ 2 (e.g., GHZ_2 has k = 2, but hypothetical states at k = 3 or k = 4 would test the 4/k prefactor).
- States that populate multiple coherent sectors simultaneously (the |+⟩^{⊗2} case from [DWELL_PREFACTOR_FROM_WEIGHTS](DWELL_PREFACTOR_FROM_WEIGHTS.md), where dΨ/dt involves √W₁ and the two-sector assumption breaks down).
- Odd-weight Pauli content, where the linear Ψ-to-f relation may not hold.

---

## References

- [DWELL_PREFACTOR_FROM_WEIGHTS](DWELL_PREFACTOR_FROM_WEIGHTS.md): Bell+ special case, partial success on |+⟩^{⊗2}
- [CRITICAL_SLOWING_AT_THE_CUSP](CRITICAL_SLOWING_AT_THE_CUSP.md): Dwell-time physics, K = γt scaling
- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): Re(λ) = -2γ⟨n_XY⟩, sector decay rates
- [ENGINEERING_BLUEPRINT](../publications/ENGINEERING_BLUEPRINT.md): Rule 1, W-type encoding recommendation
- [Simulation output](../simulations/results/dwell_prefactor_generalization_v1.txt): Full numerical results
- [Simulation script](../simulations/dwell_prefactor_generalization_v1.py): Standalone verification code
