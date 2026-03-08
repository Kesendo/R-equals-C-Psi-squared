# Glossary

**Date:** 2026-03-08
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
*See [Core Algebra](CORE_ALGEBRA.md) for the proven mathematics.*
