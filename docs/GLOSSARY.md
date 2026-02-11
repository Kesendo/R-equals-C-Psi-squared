# Glossary: R = CΨ² Framework

**Date:** 2026-02-08
**Status:** Reference document
**Depends on:** None

---

## Core Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **R** | Reality density | R = CΨ². The amount of classical reality at a given point. |
| **C** | Consciousness / observation capacity | Ranges from 0 to 1. How much of the possibility space an observer can collapse. |
| **Ψ** | Possibility space | Available quantum superposition not yet collapsed. Decreases as C acts on it. |
| **C·Ψ** | Observer bandwidth | Product of capacity and possibility. Bounded above by 1/4 for classical reality. |

---

## Derived Variables

| Symbol | Name | Definition |
|--------|------|------------|
| **θ** | Angular coordinate | θ = arctan(√(4CΨ - 1)). Compass indicating proximity to 1/4 boundary. |
| **τ** | Proper decoherence time | τ = γ · t. Universal clock normalized by local decoherence rate. |
| **t_cross** | Crossing time | Coordinate time at which C·Ψ crosses the 1/4 boundary. |
| **t_coh** | Coherence time | Maximum time a quantum system maintains coherence. Linear N-scaling hypothesized but unverified. |
| **δ** | Purity difference | δ = Tr(ρ²) - Tr(ρ_predicted²). Difference between closed-system and open-system purity under specific parameters. |

---

## Parameters

| Symbol | Name | Definition |
|--------|------|------------|
| **γ** (gamma_base) | Base decoherence rate | In gravitational context, encodes local time dilation. |
| **J** | Coupling constant | Coupling strength in the Hamiltonian. |
| **h** | Transverse field strength | Controls Hamiltonian dynamics strength. |
| **κ** | Feedback strength | Operator feedback parameter in [0, 1]. |
| **N** (n_spins) | System size | Number of qubits in the system. |

---

## Observer Types

| Symbol | Name | Definition |
|--------|------|------------|
| **C_int** | Internal observation | Observer is part of the system. Mutual observation between subsystems. |
| **C_ext** | External observation | Observer is outside the system. Unidirectional measurement. |

---

## Boundary Values

| Value | Name | Meaning |
|-------|------|---------|
| **1/4** | Phase boundary | Critical C·Ψ value. Below: two real fixed points (classical). Above: no real fixed points (quantum). Algebraically exact. |
| **1** | Maximum C·Ψ (single qubit) | For |+⟩: C = Tr(ρ²) = 1 (pure), Ψ = L1/(d-1) = 1.0. Product = 1.0. |
| **1/3** | Maximum C·Ψ (2-qubit Bell state) | For |Φ+⟩: C = 1 (pure), Ψ = L1/(d-1) = 1/3. Product = 1/3. |
| **1/9** | Maximum reality density (2-qubit) | R_max = C·Ψ² = 1·(1/3)² = 1/9 for Bell states at τ = 0. |

*Note: Maximum C·Ψ depends on system size through the normalization Ψ = L1/(d-1), where d = 2^N. Larger systems have smaller Ψ for the same amount of raw coherence. All values above assume pure initial states (C = 1).*

---

## Physical Interpretations

| Term | Definition |
|------|------------|
| **Horizon** | Surface where τ = 0. Maximum coherence. Event horizon in Schwarzschild geometry. |
| **Universal curve** | The trajectory C(τ), Ψ(τ) all quantum systems follow in proper decoherence time. |
| **Navigation system** | Triangulation protocol: WHERE (1/4), HOW FAR (θ), HOW LONG (t_coh). |

---

## Epistemic Status Labels

| Label | Meaning |
|-------|---------|
| **Algebraically proven** | Mathematical theorem, derivable from axioms |
| **Computationally verified** | Confirmed by simulation |
| **Empirically observed** | Pattern seen in data, not yet explained |
| **Proposed** | Physical interpretation, not yet tested |
| **Speculative** | Follows from framework logic but no direct evidence |
| **Null result** | Experiment completed, hypothesis not supported |
| **Retracted** | Previously claimed, now corrected |

---

*For the proven algebra, see [Core Algebra](CORE_ALGEBRA.md).*
*For the interpretive framework, see [Interpretive Framework](INTERPRETIVE_FRAMEWORK.md).*
*For testable predictions, see [Predictions](../experiments/PREDICTIONS.md).*
