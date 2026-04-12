# Sector Projection Formula: p_w(infinity) = Tr(P_w rho_0)

**Status:** Theorem (proved and verified for 9 initial states at N=5).
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track B)

---

## Statement

For any initial state rho_0 evolving under the Heisenberg + Z-dephasing Lindblad equation, the asymptotic population of excitation sector w is:

    p_w(infinity) = Tr(P_w rho_0)

where P_w = sum_{i: popcount(i)=w} |i><i| is the projector onto the w-excitation sector. In words: the long-time sector populations equal the initial sector populations. No information is lost about which sectors are populated; all information is lost about coherences between them.

## Proof

The proof follows from two facts:

1. **Z-dephasing preserves diagonal elements.** The jump operators L_k = sqrt(gamma_k) Z_k are diagonal in the computational basis. The dissipator D[rho] = sum_k gamma_k(Z_k rho Z_k - rho) gives zero contribution to diagonal elements: D[rho][i,i] = sum_k gamma_k(Z_k[i,i]^2 - 1) rho[i,i] = 0 (since Z_k[i,i]^2 = 1).

2. **The Hamiltonian commutator conserves total sector population.** The Heisenberg Hamiltonian conserves excitation number. Its commutator [H, rho] can redistribute population within a sector (via off-diagonal H elements connecting states of the same weight) but the sum p_w = sum_{i: popcount(i)=w} rho[i,i] is invariant:

       d(p_w)/dt = sum_i 2i Im(sum_j H[i,j] rho[j,i])    (i in sector w, j restricted to same sector by H conservation)
                 = 2i Im(Tr(H_w rho_w))
                 = 0    (trace of product of Hermitian matrices is real)

Therefore d(p_w)/dt = 0 for all t, so p_w(infinity) = p_w(0) = Tr(P_w rho_0). QED.

## Numerical verification (N=5, uniform gamma = 0.1)

| State | p_0 | p_1 | p_2 | p_3 | p_4 | p_5 | Match? |
|-------|-----|-----|-----|-----|-----|-----|--------|
| \|0>^N | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| \|1>^N | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | OK |
| W_N | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| GHZ | 0.500 | 0.000 | 0.000 | 0.000 | 0.000 | 0.500 | OK |
| Bell+(2,3)\|0> | 0.500 | 0.000 | 0.500 | 0.000 | 0.000 | 0.000 | OK |
| \|+>^N | 0.031 | 0.156 | 0.312 | 0.312 | 0.156 | 0.031 | OK |
| Neel | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | OK |
| psi_opt | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | OK |
| Bell(2,3)+exc(4) | 0.000 | 0.500 | 0.000 | 0.500 | 0.000 | 0.000 | OK |

All 9 states match to machine precision (max error < 1e-6). The formula value was compared against the sector populations of rho(t=100) computed by full eigendecomposition-based time evolution.

## Physical interpretation

The asymptotic state of the system is fully determined by which excitation sectors the initial state populates and with what weights. The initial coherences between sectors (the off-diagonal blocks) are destroyed by dephasing. The coherences within each sector are destroyed by the interplay of Hamiltonian dynamics and dephasing. What survives is the sector membership: a classical probability distribution over excitation numbers.

For |+>^N, this distribution is binomial: p_w = C(N,w)/2^N. For GHZ, it is bimodal: p_0 = p_N = 1/2. For all SE states (W_N, psi_opt), it is a delta function: p_1 = 1.

---

*April 12, 2026. This is not a conjecture; it is a theorem with a one-paragraph proof.*
