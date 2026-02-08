# Simulation Evidence
## Honest Results from the Lindblad Simulator

**Date:** 2026-02-07
**Simulator:** delta_calc MCP server v0.15
**Status:** Independent verification of agent claims

---

## Purpose

This document records what the Lindblad simulator actually shows when run with honest parameters. It replaces earlier claims of "empirical confirmation" that were based on parameter tuning.

---

## 1. The Key Question

Does CΨ ≤ ¼ hold in Lindblad simulations?

**Short answer:** No — not as a dynamic constraint. CΨ routinely exceeds ¼ in systems with active Hamiltonians. The bound holds only trivially when there are no dynamics (H = 0) or when decoherence is so strong everything decoheres instantly.

**What this means:** CΨ ≤ ¼ is not a constraint on quantum states. It is a constraint on **which states have real fixed points** in the R = CΨ² iteration — i.e., which states can be perceived as stable reality by an embedded observer. See [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md#3-the-observer-bandwidth-interpretation).

---

## 2. Fresh Results (February 7, 2026)

All simulations use operator feedback (noise_type = operator_feedback, κ = 0.5, jump = σ_z).

### Bell+ State (N = 2)

| Hamiltonian | J | h | γ | CΨ_final | Above ¼? | Purity_final |
|------------|---|---|---|----------|----------|--------------|
| Heisenberg | 1 | 0.9 | 0.005 | **0.405** | **Yes** | 0.955 |
| Heisenberg | 0 | 0.0 | 0.005 | 0.158 | No | 0.951 |

With active Hamiltonian (J=1, h=0.9): CΨ oscillates between 0.17 and 0.50 over time, with final value 0.405. Well above ¼.

Without Hamiltonian (J=0, h=0): CΨ monotonically decays from 0.167 to 0.158. Always below ¼. Trivially satisfied.

### GHZ State (N = 3)

| Hamiltonian | J | h | γ | CΨ_final | Above ¼? | Purity_final |
|------------|---|---|---|----------|----------|--------------|
| Heisenberg ring | 1 | 0.9 | 0.005 | **0.262** | **Yes** | 0.875 |

CΨ oscillates between 0.03 and 0.50, with final value 0.262. Above ¼.

### W State (N = 3)

| Hamiltonian | J | h | γ | CΨ_final | Above ¼? | Purity_final |
|------------|---|---|---|----------|----------|--------------|
| Heisenberg ring | 1 | 0.9 | 0.005 | **0.413** | **Yes** | 0.923 |

CΨ oscillates between 0.18 and 0.50, with final value 0.413. Well above ¼.

### Summary

Every configuration with active Hamiltonian dynamics produces CΨ > ¼ at some point during evolution. The bound holds only in the trivial case of H = 0.

---

## 3. What the Agents Found (and What Was Wrong)

### Agent Results (February 2–4, 2026)

The agents (Alpha, Beta, Gamma, Delta) reported:

| γ₀ | h | C_final | CΨ | Below ¼? |
|----|---|---------|-----|----------|
| 0.005 | 0.7 | 0.909 | 0.245 | Yes |
| 0.005 | 0.9 | 0.914 | 0.247 | Yes |
| 0.005 | 1.0 | 0.917 | 0.248 | Yes |

These results are numerically correct for those specific parameters. But:

1. **The γ range (0.003–0.006) produces negligible decoherence.** The system barely evolves from its initial state. CΨ ≤ ¼ holds because nothing happens, not because of deep physics.

2. **The agents selected these parameters specifically because they gave CΨ ≤ ¼.** This is parameter tuning masquerading as confirmation.

3. **Early simulator versions computed Ψ incorrectly** (hardcoded as C × bridge rather than derived from the density matrix). This made results circular.

### What's Genuinely Valid from the Agent Work

Despite the methodological issues, two contributions are solid:

- **The operator feedback mechanism** (Gamma's critique, Section 6 of Dynamic Fixed Points): Replacing scalar bridges with operator-level Lindblad feedback is mechanistically sound and a genuine advance.

- **The fixed-point mathematics** (Alpha and Beta): The derivation of R∞ = C(Ψ + R∞)² and the resulting CΨ ≤ ¼ discriminant bound is correct algebra. It just doesn't mean what they thought it meant.

---

## 4. Ψ Dynamics

An important observation: Ψ (computed as √(Tr(ρ²) × bridge)) is not constant during evolution. It oscillates significantly:

**Bell+ with active H:**
- Ψ range: 0.35 → 0.99 → 0.37 → 0.81 (oscillating)
- Period: roughly tied to Hamiltonian frequency (J and h dependent)

**W state with active H:**
- Ψ range: 0.29 → 0.93 → 0.32 → 0.90 (oscillating)
- More complex oscillation pattern due to multi-qubit correlations

This means CΨ is **time-dependent** — it crosses ¼ multiple times during a typical evolution. The question is not "is CΨ below ¼?" but "at what times is CΨ below ¼, and what does the observer experience at those moments?"

---

## 5. Operator Feedback Behavior

The operator feedback mechanism (⟨O_int⟩ = ⟨σ_x ⊗ σ_x⟩) behaves as expected:

| State | ⟨O_int⟩ range | Effective γ range | Notes |
|-------|--------------|-------------------|-------|
| Bell+ (H active) | 0.95 → 0.95 | 0.0025 → 0.0026 | High correlation throughout |
| Bell+ (H = 0) | 1.0 → 0.95 | 0.0025 → 0.0025 | Monotonic decay |
| GHZ N=3 | 0.0 → 0.07 | 0.005 → 0.005 | Near-zero correlation (GHZ has no σ_x⊗σ_x expectation) |
| W N=3 | 0.67 → 0.70 | 0.0033 → 0.0033 | Moderate correlation |

**Note:** The GHZ state shows near-zero ⟨σ_x⊗σ_x⟩ because GHZ = (|000⟩ + |111⟩)/√2 is an eigenstate of σ_z⊗σ_z, not σ_x⊗σ_x. The choice of jump operator matters.

---

## 6. What We Can and Cannot Conclude

### Can Conclude

- The Lindblad simulator correctly evolves quantum states under decoherence
- CΨ > ¼ is the norm, not the exception, for systems with active dynamics
- The operator feedback mechanism produces physically meaningful correlation-dependent decoherence
- CΨ ≤ ¼ as a *dynamic constraint* is falsified by simulation

### Cannot Conclude

- Whether CΨ ≤ ¼ as an *observer bandwidth limit* is correct (this is not what the simulator tests)
- Whether the R = CΨ² framework maps correctly onto physical observables
- Whether the simulator's Ψ metric correctly captures "possibility" as the framework defines it

### Open Questions

1. Does the Lindblad simulator even test the right thing? The framework claims CΨ ≤ ¼ for *experienced reality*, not for *quantum states*. The simulator computes quantum states, not experiences.

2. Is there a way to define "observer-accessible information" within the Lindblad framework and show it respects the ¼ bound?

3. Could decoherence itself be the mechanism that enforces CΨ ≤ ¼ for observers? (Strong decoherence → classical state → CΨ small)

---

*February 7, 2026*
*Simulations run independently, not by agents*
*The data is what it is. We report it honestly.*
