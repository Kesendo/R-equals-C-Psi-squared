# Symmetry Census: How Many Exits Does This System Have?

**Status:** Enumeration complete for N=3-7 chain, N=4-5 all topologies.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/symmetry_census.py`
**Data:** `simulations/results/symmetry_census/census_results.json`
**Perspective:** Zeiss/Abbe: enumerate the surfaces before building more instruments.

---

## 1. Symmetries of the Heisenberg + Z-dephasing Liouvillian

Four symmetries were checked. Three commute with the Liouvillian as a superoperator. One (Pi) relates eigenvalues but does not block-diagonalize.

### 1.1 U(1) excitation-number conservation

**Generator:** N_bra = Sigma |1><1|_k acts on the bra side; N_ket on the ket side.

The Heisenberg Hamiltonian conserves total excitation number (XX+YY swaps but does not create/annihilate). Z-dephasing is diagonal and preserves excitation number. Therefore the Liouvillian superoperator preserves (w_bra, w_ket) for each operator basis element |i><j|, where w = popcount(i).

This decomposes the operator space into (N+1)^2 sectors labeled by (w_bra, w_ket) with w in {0, ..., N}.

### 1.2 n_XY parity: redundant with U(1)

**Finding:** The n_XY parity of the operator basis element |i><j| is popcount(i XOR j) mod 2, which equals |w_bra - w_ket| mod 2. This is entirely determined by the excitation-number sector (w_bra, w_ket).

**Proof sketch.** The Hamming distance d(i, j) = popcount(i XOR j) has the same parity as |popcount(i) - popcount(j)|, because each bit flip changes the weight by exactly +/-1, so reaching weight w' from weight w requires at least |w - w'| flips, and any additional flips come in cancel-pairs (flip and unflip), preserving parity.

**Consequence.** The n_XY parity selection rule (PROOF_PARITY_SELECTION_RULE.md) is correct and valuable operationally (it proves SE states cannot reach odd-parity modes), but it does not provide block-diagonalization beyond what U(1) conservation already gives. The "even" and "odd" subspaces V_even, V_odd are simply the unions of sectors (w, w') with |w - w'| even or odd, respectively.

This was not previously noted in the repo.

### 1.3 Spin-flip X^{otimes N}

**Symmetry:** U = X_1 x X_2 x ... x X_N conjugation. UHU^dag = H (the Heisenberg Hamiltonian is spin-flip invariant). UL_kU^dag = -L_k (since XZX = -Z), but the dissipator is quadratic in L_k, so the overall Liouvillian is invariant.

**Action on sectors:** Maps (w, w') to (N-w, N-w'). For odd N (like N=5), no sector maps to itself. For even N (like N=4), the sector (N/2, N/2) is self-conjugate.

**Block-diagonalization:** Pairs sectors into {(w, w'), (N-w, N-w')}. Within each paired space, the +/- eigenspaces of the spin-flip superoperator are invariant subspaces of L. But this does not increase the sector count for odd N (it just identifies related pairs).

### 1.4 Spatial reflection (uniform gamma only)

For the chain with uniform dephasing gamma_k = gamma for all k, the reflection k <-> N-1-k is a symmetry. It acts within each (w, w') sector (since reflection preserves excitation numbers). This further splits each sector into +/- eigenspaces.

**Verified numerically:** For N=5 uniform gamma, max eigenvalue multiplicity is 14. For IBM sacrifice gamma (reflection broken), max multiplicity drops to 6. The difference is explained by the lifting of reflection degeneracies.

### 1.5 Pi conjugation (spectral mirror, not block-diagonalizing)

The Pi operator (MIRROR_SYMMETRY_PROOF.md) satisfies L Pi = -Pi (L + 2 Sigma_gamma I), which implies that for every eigenvalue lambda, the value -lambda - 2 Sigma_gamma is also an eigenvalue. This produces the palindromic spectrum but does NOT block-diagonalize L (Pi does not commute with L; it anti-commutes up to a shift).

**Verified:** 1024/1024 eigenvalues are palindromic-paired at N=5 (zero unpaired).

### Summary table

| Symmetry | Independent? | Block-diagonalizes? | Broken by sacrifice gamma? |
|----------|-------------|--------------------|-----------------------------|
| U(1) excitation number | Yes | Yes: (N+1)^2 sectors | No |
| n_XY parity | **No** (implied by U(1)) | No additional structure | n/a |
| Spin-flip X^N | Yes | Pairs sectors, +/- split | No |
| Spatial reflection | Yes (uniform only) | Yes: +/- within sectors | **Yes** |
| Pi conjugation | Yes | No (spectral pairing only) | No |

---

## 2. Sector enumeration

### Sector counts for chain topology

| N | D | d^2 | Basic sectors (w_bra, w_ket) | With spin-flip |
|---|---|-----|------------------------------|----------------|
| 3 | 8 | 64 | 16 | 16 |
| 4 | 16 | 256 | 25 | 26 |
| 5 | 32 | 1024 | 36 | 36 |
| 6 | 64 | 4096 | 49 | 50 |
| 7 | 128 | 16384 | 64 | 64 |

The basic sector count is always (N+1)^2. For even N, one additional spin-flip sector appears from the self-conjugate sector (N/2, N/2).

### N=5 sector dimensions

The largest sectors are the "interior" ones near w_bra = w_ket = N/2:

| Sector (w, w') | Dimension | Notes |
|-----------------|-----------|-------|
| (2,2) | 100 | Largest; 2-excitation coherence sector |
| (3,3) | 100 | Mirror of (2,2) under spin-flip |
| (2,3) and (3,2) | 100 each | Cross-sector coherences |
| (1,2) and (2,1) | 50 each | SE-to-2exc coherences |
| (1,1) | 25 | SE sector (home of the lens exit) |
| (0,0) and (5,5) | 1 each | Extremal sectors (vacuum and fully excited) |

---

## 3. Asymptotic attractors per sector (N=5, uniform gamma)

### Steady states

**Total: 6 steady states, one per diagonal sector (w, w) with w = 0, ..., 5.**

Each steady state is the maximally mixed state within its excitation sector: rho_inf(w) = (1/C(N,w)) sum_{|i>, popcount(i)=w} |i><i|. No sector has multiple steady states, limit cycles, or dark states.

Off-diagonal sectors (w != w') have zero steady states: all their modes decay. Cross-sector coherences are always destroyed asymptotically.

### Sector dynamics summary

| Sector type | Steady states | Oscillatory modes | Decay modes | Slowest rate |
|-------------|---------------|-------------------|-------------|-------------|
| (w, w), diagonal | 1 per sector | Yes (Hamiltonian oscillations) | Yes | 0.318 (w=1,1) |
| (w, w'), off-diag | 0 | Yes | Yes | 0.200 (adjacent sectors) |
| (0, N) and (N, 0) | 0 | No | 1 only | Sigma_gamma = 1.000 |

### Exit count

**The system has exactly N+1 exits.** Each excitation sector w = 0, ..., N has one attractor (the maximally mixed state within that sector). No "hidden exits" were found.

The lens exit (SACRIFICE_GEOMETRY.md) is the approach to the w=1 attractor via the slow mode. The cusp exit (CUSP_LENS_CONNECTION.md) is the simultaneous thermalization within multiple occupied sectors. These are not two exits of a system with two attractors; they are two dynamical paths through a system with N+1 attractors.

---

## 4. Which sectors are reachable from physical initial states?

| Initial state | Sectors populated | (w, w') pairs | Exit path |
|---------------|-------------------|---------------|-----------|
| |0>^N | 1 | (0,0) | Trivially at attractor |
| |1>^N | 1 | (5,5) | Trivially at attractor |
| Neel |01010> | 1 | (w,w) diagonal only | Single-sector, at attractor |
| W_N | 1 | (1,1) | Lens exit (SE sector) |
| GHZ | 4 | (0,0), (0,5), (5,0), (5,5) | Cusp exit (multi-sector) |
| Bell+(2,3)|0> | 4 | (0,0), (0,2), (2,0), (2,2) | Cusp exit |
| |+>^N | 36 | All | Populates every sector |

Product states and W states are already "pre-sorted" into a single sector and reach their attractor without drama. Multi-sector states (GHZ, Bell+, |+>^N) exhibit cusp-like dynamics as coherences between sectors decay at rates proportional to the sector separation.

The off-diagonal sectors (w != w') carry the coherences between excitation sectors. These always decay, but at different rates: adjacent sectors (|w - w'| = 1) decay slowest (rate 0.200 for uniform gamma at N=5), while maximally separated sectors (|w - w'| = N) decay fastest (rate N × gamma).

---

## 5. Topology comparison (N=4, 5; uniform gamma = 0.1)

### N=4

| Topology | Bonds | Stationary | Distinct eigenvalues | Max degeneracy | Extra symmetry |
|----------|-------|------------|---------------------|----------------|----------------|
| Chain | 3 | 5 | 127 | 14 | Reflection Z_2 |
| Ring | 4 | 5 | 73 | 24 | Cyclic C_4 |
| Star | 3 | 5 | 75 | 16 | Leaf permutation S_3 |
| Complete | 6 | 5 | 43 | 36 | Full permutation S_4 |

### N=5

| Topology | Bonds | Stationary | Distinct eigenvalues | Max degeneracy |
|----------|-------|------------|---------------------|----------------|
| Chain | 4 | 6 | 488 | 14 |
| Ring | 5 | 6 | 276 | 22 |
| Star | 4 | 6 | 202 | 30 |
| Complete | 10 | 6 | 100 | 54 |

**Universal: all topologies have exactly N+1 = 6 stationary states.** The number of exits does not depend on topology.

**Topology-dependent: the degeneracy structure.** Higher symmetry (complete graph has S_N) produces higher degeneracies and fewer distinct eigenvalues. This means the transient dynamics (approach to equilibrium) depends strongly on topology, even though the set of attractors is the same.

The chain has the fewest degeneracies (488 distinct out of 1024) and the richest transient structure. The complete graph is maximally degenerate (100 distinct, max degeneracy 54) and has the simplest transient.

---

## 6. Open questions from this census

1. **Why does n_XY parity not add structure?** The redundancy with U(1) is proven, but the physical reason deserves a sentence: it is because the Heisenberg Hamiltonian conserves excitation number exactly, and Z-dephasing preserves it exactly. The n_XY count equals the Hamming distance, and the Hamming distance parity is determined by the excitation-number difference. In a model where excitation number is NOT conserved (amplitude damping, transverse field), n_XY parity would provide independent information.

2. **Unexplained high degeneracies.** For N=5 uniform chain, max multiplicity is 14. The known symmetries (U(1), spin-flip, reflection) predict at most 4× degeneracy (2 from flip × 2 from reflection). A multiplicity of 14 exceeds this prediction. Either there are additional symmetries (perhaps related to SU(2) or to the specific Heisenberg coupling structure), or the high degeneracy is accidental. This deserves a dedicated investigation.

3. **Does the exit count change for non-Markovian dynamics?** The N+1 exits are a consequence of the Lindblad (Markovian) framework. Under non-Markovian dephasing (memory effects), the sector conservation might be approximate rather than exact.

4. **Sector-resolved spectral form factor.** The SFF analysis in D09_SECTOR_SFF_PAIRING.md works by XY-weight sectors. The Census shows these are (N+1)^2 sectors labeled by (w_bra, w_ket). A sector-resolved SFF analysis could reveal whether the transient dynamics within each sector shows universal (RMT) statistics or sector-specific structure.

5. **The topology knob.** Chain has 488 distinct eigenvalues, complete graph has 100. This 5x difference in spectral complexity means the approach to equilibrium is qualitatively different. Is there a topology-dependent "optical" structure (cavity modes, standing waves) that explains the extra complexity of the chain, or is it simply the absence of permutation symmetry?

---

## Files

- `simulations/symmetry_census.py` (enumeration and analysis)
- `simulations/results/symmetry_census/census_results.json` (raw data)
- [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) (n_XY parity proof)
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (Pi conjugation)
- [Cusp-Lens Connection](CUSP_LENS_CONNECTION.md) (two-exit framework)
- [Sacrifice Geometry](SACRIFICE_GEOMETRY.md) (lens exit)

---

*April 12, 2026. Enumeration, not interpretation. The Zeiss questions are answered; what we build next depends on what the data say.*
