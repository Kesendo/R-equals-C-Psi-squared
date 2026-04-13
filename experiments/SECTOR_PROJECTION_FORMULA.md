# Sector Projection Formula: p_w(∞) = Tr(P_w ρ_0)

**Status:** Theorem (proved and verified for 9 initial states at N=5).
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track B)
**Proof:** [PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)

---

## Statement

For any initial state ρ_0 evolving under the Heisenberg + Z-dephasing Lindblad equation, the asymptotic population of excitation sector w is:

    p_w(∞) = Tr(P_w ρ_0)

where P_w = Σ_{i: popcount(i)=w} |i⟩⟨i| is the projector onto the w-excitation sector. In words: the long-time sector populations equal the initial sector populations. No information is lost about which sectors are populated; all information is lost about coherences between them.

## Proof

The proof follows from two facts:

1. **Z-dephasing preserves diagonal elements.** The jump operators L_k = √γ_k Z_k are diagonal in the computational basis. The dissipator D[ρ] = Σ_k γ_k(Z_k ρ Z_k − ρ) gives zero contribution to diagonal elements: D[ρ][i,i] = Σ_k γ_k(Z_k[i,i]² − 1) ρ[i,i] = 0 (since Z_k[i,i]² = 1).

2. **The Hamiltonian commutator conserves total sector population.** The Heisenberg Hamiltonian conserves excitation number. Its commutator [H, ρ] can redistribute population within a sector (via off-diagonal H elements connecting states of the same weight) but the sum p_w = Σ_{i: popcount(i)=w} ρ[i,i] is invariant:

       d(p_w)/dt = −i Σ_{i∈w} [H, ρ][i,i]
                 = −i Tr([H_w, ρ_w])
                 = 0    (Tr([A, B]) = Tr(AB) − Tr(BA) = 0 by cyclicity of trace)

Therefore d(p_w)/dt = 0 for all t, so p_w(∞) = p_w(0) = Tr(P_w ρ_0). QED.

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

The asymptotic state of the system is fully determined by which excitation sectors the initial state populates and with what weights. The initial coherences between sectors (the off-diagonal blocks) are destroyed by dephasing. The coherences within each sector are destroyed by the interplay of Hamiltonian dynamics and dephasing. What survives is the sector membership: a classical probability distribution over excitation numbers.

For |+⟩⊗N, this distribution is binomial: p_w = C(N,w)/2^N. For GHZ, it is bimodal: p_0 = p_N = ½. For all SE states (W_N, ψ_opt), it is a delta function: p_1 = 1.

---

## Files

- `simulations/three_values.py` (Track B: verification computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) (sector conservation theorem)
- [Symmetry Census](SYMMETRY_CENSUS.md) (N+1 attractor enumeration)

---

*April 12, 2026. This is not a conjecture; it is a theorem with a one-paragraph proof.*
