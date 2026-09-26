# Sector Projection Formula: p_w(t) = Tr(P_w ρ_0)

**Status:** Theorem for sector populations at every time (proved and verified for 9 initial states at N=5). A formula for the full asymptotic state requires the separate full-support hypotheses of the linked proof.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track B)
**Proof:** [the Asymptotic Sector Projection proof](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)

---

## Statement

For any initial state ρ_0 evolving under the Heisenberg + Z-dephasing Lindblad equation, the population of excitation sector w at every time is:

    p_w(t) = Tr(P_w ρ_0)

where P_w = Σ_{i: popcount(i)=w} |i⟩⟨i| is the projector onto the w-excitation sector. In words: the sector populations are conserved. This identity alone makes no assertion about coherences or the existence of a limit for ρ(t).

## Proof

The proof follows from two facts:

1. **Z-dephasing preserves diagonal elements.** The jump operators L_k = √γ_k Z_k are diagonal in the computational basis. The dissipator D[ρ] = Σ_k γ_k(Z_k ρ Z_k − ρ) gives zero contribution to diagonal elements: D[ρ][i,i] = Σ_k γ_k(Z_k[i,i]² − 1) ρ[i,i] = 0 (since Z_k[i,i]² = 1).

2. **The Hamiltonian commutator conserves total sector population.** The Heisenberg Hamiltonian conserves excitation number. Its commutator [H, ρ] can redistribute population within a sector (via off-diagonal H elements connecting states of the same weight) but the sum p_w = Σ_{i: popcount(i)=w} ρ[i,i] is invariant:

       d(p_w)/dt = −i Σ_{i∈w} [H, ρ][i,i]
                 = −i Tr([H_w, ρ_w])
                 = 0    (Tr([A, B]) = Tr(AB) − Tr(BA) = 0 by cyclicity of trace)

Therefore d(p_w)/dt = 0 for all t, so p_w(t) = p_w(0) = Tr(P_w ρ_0). QED.

## Numerical verification (N=5, uniform γ = 0.1)

| State | p_0 | p_1 | p_2 | p_3 | p_4 | p_5 | Match? |
|-------|-----|-----|-----|-----|-----|-----|--------|
| \|0⟩⊗N | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| \|1⟩⊗N | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | OK |
| W_N | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| GHZ | 0.500 | 0.000 | 0.000 | 0.000 | 0.000 | 0.500 | OK |
| Bell+(2,3)\|0⟩ | 0.500 | 0.000 | 0.500 | 0.000 | 0.000 | 0.000 | OK |
| \|+⟩⊗N | 0.031 | 0.156 | 0.312 | 0.312 | 0.156 | 0.031 | OK |
| Neel | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | OK |
| ψ_opt | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| Bell(2,3)+exc(4) | 0.000 | 0.500 | 0.000 | 0.500 | 0.000 | 0.000 | OK |

All 9 states match to machine precision (max error < 10⁻⁶). The formula value was compared against the sector populations of ρ(t=100) computed by full eigendecomposition-based time evolution.

## Physical interpretation

The weights of the excitation sectors remain fixed while the state moves. For a connected graph with a positive Z-dephasing rate at **every** site, the [separate asymptotic proof](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) establishes the maximally mixed state within each sector. With sparse dephasing, conserved sector weights do not determine a full-state limit: the [N=5 centre-only odd block](OPERATOR_PAIR_VIEW_COMPARISON.md) carries nonstationary, undamped motion. The population formula above remains exact in that case.

For |+⟩⊗N, this distribution is binomial: p_w = C(N,w)/2^N. For GHZ, it is bimodal: p_0 = p_N = ½. For all SE states (W_N, ψ_opt), it is a delta function: p_1 = 1.

---

## Files

- `simulations/three_values.py` (Track B: verification computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) (sector conservation theorem)
- [Symmetry Census](SYMMETRY_CENSUS.md) (N+1 attractor enumeration)

---

*April 12, 2026. This is not a conjecture; it is a theorem with a one-paragraph proof.*
