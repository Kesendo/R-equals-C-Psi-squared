# Glossary

**Date:** 2026-03-16
**Purpose:** Notation and term reference for the repository

**Tier:** Reference
**Status:** Living document, aligned with [THE_CPSI_LENS](THE_CPSI_LENS.md) as canonical description

---

## Core notation

The symbol C is used in three contexts in this repository. See the Notation section in [THE_CPSI_LENS](THE_CPSI_LENS.md) for details.

| Symbol | Algebraic usage | Experimental usage | Legacy usage (deprecated) |
|--------|----------------|-------------------|--------------------------|
| **C** | Coupling scalar in the recurrence | Wootters concurrence | "Consciousness" |
| **Ψ** | Scalar in the recurrence | Normalized l1-coherence | "Possibility" |
| **CΨ** | Product c in the Mandelbrot iteration | Concurrence x coherence (the lens) | "Consciousness x Possibility" |
| **R** | Fixed-point value of the iteration | CΨ^2 (product of concurrence, coherence squared) | "Reality" |

---

## Derived quantities

| Symbol | Definition |
|--------|------------|
| **θ** | arctan(sqrt(4CΨ - 1)). Angular distance from the 1/4 boundary in the complex regime. |
| **τ** | γ * t. Normalized decoherence time. |
| **t_cross** | Time at which CΨ crosses 1/4. |
| **δ** | Tr(ρ^2) - Tr(ρ_predicted^2). Purity residual from a simple noise model. |

---

## Palindromic symmetry (proven March 14, 2026)

| Symbol / Term | Definition |
|---------------|------------|
| **Π** | Conjugation operator. Per-site action: I→X, X→I, Y→iZ, Z→iY. Satisfies Π·L·Π⁻¹ = -L - 2Σγ·I. Maps every decay rate d to its palindromic partner 2Σγ-d. This is one member of the P1 family. A second family P4 (I↔Y, X↔Z) supports additional terms (XZ, ZX). XY/YX terms require non-uniform (alternating) operators. See [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md), [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md). |
| **Σγ** | Sum of all individual dephasing rates: Σγ = γ₁ + γ₂ + ... + γ_N. The palindromic axis sits at -Σγ. |
| **Palindromic pairing** | Every non-zero Liouvillian eigenvalue d has a partner d' such that Re(d) + Re(d') = -2Σγ. Proven for Heisenberg+dephasing, verified for all standard models (XY, Ising, XXZ, DM) under single-axis dephasing. Depolarizing noise breaks palindrome at err ≈ γ·2(N-2)/3. Verified N=2-8. |
| **Pauli weight** | Number of non-identity Pauli operators in a Pauli string. E.g., XYI has weight 2. Π maps weight k to N-k (complementarity). Equivalent to "XY-weight" in Haga et al. (2023). |
| **Incoherenton** | Term from Haga et al. (2023). Quasiparticle that counts Pauli weight in open quantum systems. Their XY-weight = our Pauli weight. Their particle-hole transformation = our Π. |

---

## XOR space (discovered March 16, 2026)

| Term | Definition |
|------|------------|
| **XOR modes** | The N+1 Liouvillian eigenmodes at λ = -2Σγ (maximum decay rate). They are not paired with other modes; their palindromic partner is the steady state (λ=0). They are purely off-diagonal (coherences). |
| **XOR fraction** | Fraction of an input state's weight that projects onto the XOR modes. GHZ: 100%. W (N≥3): 0%. Predicted by mixed XY Pauli weight at r=0.976 (N≥3). |
| **Mixed XY Pauli weight** | Fraction of Pauli decomposition terms containing both X and Y operators simultaneously (e.g., XYI, YXZ). Predicts how much of a state falls into the fastest-decaying XOR drain. |
| **Palindromic modes** | All Liouvillian modes except the N+1 XOR modes and steady state. They come in pairs at various decay rates. States in palindromic modes are more robust because some pairs decay slowly. |
| **Spectral filter** | The palindrome acts as a filter: it separates every input into a fragile XOR component (fast decay, mixed XY) and a distributable palindromic component (various rates, some survive). |

---

## Quantum state transfer (verified March 14, 2026)

| Term | Definition |
|------|------------|
| **QST** | Quantum State Transfer. Moving a quantum state from Alice (A) to Bob (B) through a mediator (S). |
| **F_avg** | Average fidelity of the transferred state, averaged over all pure input states. Our best: F_avg = 0.888 (star, 2:1 coupling). |
| **Holevo capacity** | Maximum classical information transmittable per channel use. Our star channel: χ = 0.534 bits. |
| **2:1 coupling** | Optimal coupling ratio J_SB/J_SA = 2 for star-topology QST. Asymmetric. Not intuitive. Outperforms symmetric 1:1. |
| **Echo** | Entanglement oscillation in the SA pair. Period ~π/(4J). Envelope decays at 8γ/3. Peak C_SB = 0.598 (N=3). Scales as ~1/(N-1). |

---

## Parameters

| Symbol | Definition |
|--------|------------|
| **γ** | Decoherence rate (Lindblad dephasing strength). Higher γ = faster loss of quantum coherence. In star topology: each qubit has individual γ. Receiver noise (γ_A) is more destructive than sender noise (γ_B). |
| **J** | Hamiltonian coupling strength between qubits. In star topology: J_SA (A-to-S coupling), J_SB (B-to-S coupling). |
| **h** | Transverse field strength in the Hamiltonian. |

---

## Star topology (S + A + B)

| Symbol | Definition |
|--------|------------|
| **S** | Central qubit (shared object), coupled to both A and B |
| **A, B** | Observer qubits, not directly coupled to each other |
| **CΨ_AB** | CΨ for the AB reduced pair (traced over S). Crossing 1/4 means the observer-observer connection is active. |

### Three conditions for AB crossing (at γ = 0.05)

1. **Strong sender:** J_SB/J_SA >= 1.466
2. **Quiet receiver:** γ_A < ~0.20
3. **Deep pre-existing connection:** Initial C_SA > 0.8 (Bell-like)

---

## Observer distinction (Lindblad model)

| Term | Technical definition | Status |
|------|---------------------|--------|
| **C_int** | Symmetric Hamiltonian coupling (mutual interaction) | Tier 2 formalization |
| **C_ext** | Projective measurement (one-directional intervention) | Tier 2 formalization |

Note: The original claim that C_int preserves coherence 33x longer than C_ext has been disproven (see MATHEMATICAL_FINDINGS Section 9). The formal Lindblad distinction is verified; the physical claim about different "kinds of observation" is not.

---

## Boundary values

| Value | Source |
|-------|--------|
| **1/4** | Discriminant of the fixed-point equation. Below 1/4: real fixed points exist. Above: complex. Algebraically exact within the iteration. |
| **~1.466** | J_SB/J_SA threshold for AB crossing at γ = 0.05 (star topology). |
| **1/3** | CΨ of a maximally entangled Bell pair (C = 1, Ψ = 1/3). |
| **2γ** | Decay rate of the c+ supermode (even, slow). Topology-independent for N=3. |
| **8γ/3** | Decay rate of the concurrence envelope. Topology-independent for N=3. |
| **10γ/3** | Decay rate of the c- supermode (odd, fast). Topology-independent for N=3. |
| **-2Σγ** | Location of XOR modes (maximum decay rate). Always N+1 modes here. |
| **0.888** | Best average fidelity for QST (star topology, J_SB/J_SA = 2:1, γ = 0.05). |
| **0.039/γ** | Approximate crossing time t_cross for Bell+ under Heisenberg+dephasing. |
| **0.976** | Correlation between mixed XY Pauli weight and XOR fraction (N≥3). |

---

## Epistemic tier labels

| Tier | Meaning |
|------|---------|
| **1** | Algebraically proven |
| **2** | Computationally verified |
| **3** | Proposed interpretation or speculative extension |
| **4** | Agent-generated, not independently verified |
| **5** | Speculative philosophy |

---

*See [THE_CPSI_LENS](THE_CPSI_LENS.md) for the canonical project description.*
*See [Core Algebra](historical/CORE_ALGEBRA.md) for the proven mathematics.*
*See [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md) for the palindrome theorem.*
*See [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md) for the extended palindrome analysis.*
*See [XOR Space](../experiments/XOR_SPACE.md) for the spectral filter discovery.*
