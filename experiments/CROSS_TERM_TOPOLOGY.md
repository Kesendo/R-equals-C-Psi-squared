# Cross-Term Topology Dependence

<!-- Keywords: cross-term topology dependence, anti-commutator L_H L_Dc,
relative orthogonality 1/sqrt(N*2^(N+1)), topology universality,
Pythagorean breaking, Heisenberg chain ring star complete graph -->

**Status:** Tier 2 (computed)
**Date:** April 13, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [cross_term_topology.py](../simulations/cross_term_topology.py)
**Output:** [cross_term_topology/](../simulations/results/cross_term_topology/)
**Depends on:**
- [Primordial Qubit Algebra](PRIMORDIAL_QUBIT_ALGEBRA.md) (parent; this resolves "What to compute next" item 1)
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (definition of the palindromic structure)
**Task:** [TASK_CROSS_TERM_TOPOLOGY.md](../ClaudeTasks/TASK_CROSS_TERM_TOPOLOGY.md)

---

## What this document is about

At N=2, the Hamiltonian superoperator L_H and the centered dissipator
L_Dc = L_D + Sg*I are exactly orthogonal: {L_H, L_Dc} = 0. This is the
Pythagorean decomposition from [PRIMORDIAL_QUBIT_ALGEBRA](PRIMORDIAL_QUBIT_ALGEBRA.md).

At N=3 on a Heisenberg chain, the relative orthogonality
||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) equals 1/sqrt(48), a
gamma-independent geometric constant. The open question was: is this
constant chain-specific, or does it hold for other topologies?

**Answer: the constant is topology-independent.** It depends only on N.

---

## The result

### N=3

| Topology | Edges | rel_ortho | Closed form |
|----------|-------|-----------|-------------|
| chain | 2 | 0.1443375673 | 1/sqrt(48) |
| star | 2 | 0.1443375673 | 1/sqrt(48) |
| ring | 3 | 0.1443375673 | 1/sqrt(48) |
| complete | 3 | 0.1443375673 | 1/sqrt(48) |

### N=4

| Topology | Edges | rel_ortho | Closed form |
|----------|-------|-----------|-------------|
| chain | 3 | 0.0883883476 | 1/sqrt(128) |
| star | 3 | 0.0883883476 | 1/sqrt(128) |
| ring | 4 | 0.0883883476 | 1/sqrt(128) |
| complete | 6 | 0.0883883476 | 1/sqrt(128) |

All eight numbers match to machine precision (deviation < 10^-15). The
cross-term magnitude is identical across chain, star, ring, and complete
graph, despite different edge counts, different degree sequences, and
different graph symmetry groups.

### Gamma-independence

Verified for all (topology, N) pairs with gamma in {0.001, 0.01, 0.05, 0.1, 0.5}.
Maximum spread across gamma values: < 2 * 10^-16 (machine zero).

The relative orthogonality is a pure geometric property of the
Heisenberg coupling type and the Z-dephasing structure, independent of
both the dephasing rate and the graph topology.

---

## Topology equivalence at N=3

Two degeneracies exist at N=3:

1. **chain = star** (graph isomorphism): both are the path graph P_3
   with degree sequence [1, 1, 2]. The Liouvillian spectra are identical
   under the site relabeling that maps one to the other, so they must
   give the same cross-term.

2. **ring = complete** (identical edge set): at N=3, the ring
   0-1-2-0 and the complete graph K_3 both have edges
   {(0,1), (0,2), (1,2)}. They are literally the same graph.

At N=4, all four topologies are genuinely distinct (different edge sets,
different degree sequences). The cross-term equality at N=4 is therefore
the non-trivial result.

---

## Conjecture: N-scaling formula

The denominators factor cleanly:

| N | rel_ortho^2 | Denominator | Factorisation |
|---|-------------|-------------|---------------|
| 2 | 0 (exact) | n/a | Pythagorean theorem holds |
| 3 | 1/48 | 48 | 3 * 2^4 = N * 2^(N+1) |
| 4 | 1/128 | 128 | 4 * 2^5 = N * 2^(N+1) |

**Conjecture (topology-universal N-scaling):**

For N >= 3 qubits with Heisenberg XXX coupling on any graph and uniform
Z-dephasing at rate gamma per site:

    ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) = 1 / sqrt(N * 2^(N+1))

This is a geometric constant: independent of gamma, independent of J,
and independent of the graph topology.

At N=2, the formula gives 1/sqrt(16) = 1/4, but the actual value is 0.
The N=2 exception is explained in PRIMORDIAL_QUBIT_ALGEBRA: at N=2 every
L_H transition satisfies w_XY(a) + w_XY(b) = N, making the
anti-commutator vanish identically. This constraint is impossible to
satisfy at odd N (the sum is always even while N is odd) and only
partially satisfied at even N >= 4.

**Status:** conjecture with evidence at N=3 and N=4, both verified across
four topologies and five gamma values. Confirmation at N=5 would require
diagonalisation of a 1024x1024 superoperator (feasible). Analytical proof
would follow from the w_XY transition statistics of the Heisenberg bond.

---

## Hypothesis tests

### Primary hypothesis: topology-specific constants

**REFUTED.** The original hypothesis (PRIMORDIAL_QUBIT_ALGEBRA "What to
compute next" item 1) asked whether each topology has its own constant.
The answer is no: the constant depends only on N. This is a stronger
result than anticipated.

### Secondary hypothesis: complete graph vanishing at N=3

**REFUTED.** The prediction was that the complete graph might restore
orthogonality at N=3 because every site participates in every bond,
eliminating "spectator sites." The complete graph at N=3 has the same
cross-term as the chain (1/sqrt(48)).

### New finding: topology universality

The cross-term relative orthogonality is a **graph-invariant** for
Heisenberg XXX coupling with uniform Z-dephasing. It depends on:
- N (number of qubits): via the formula 1/sqrt(N * 2^(N+1))
- The coupling type: Heisenberg XXX (other couplings not tested here)
- The noise type: Z-dephasing (other noise types explored in
  [pythagorean_breaking.py](../simulations/pythagorean_breaking.py))

It does NOT depend on:
- The graph topology (chain, star, ring, complete: all identical)
- The coupling strength J (cancels in the ratio)
- The dephasing rate gamma (cancels in the ratio)

---

## Why topology cancels

Each Heisenberg bond (i,j) contributes to L_H a set of transitions
between Pauli strings that differ at sites i and j. The w_XY change at
the bond is determined by the Heisenberg coupling structure; the
remaining N-2 "spectator" sites contribute their w_XY unchanged.

The anti-commutator {L_H, L_Dc} is a sum over bonds. For each bond,
the contribution to the anti-commutator depends on the w_XY distribution
of transitions, which is the same for every bond (since all bonds are
Heisenberg XXX). The norm ||{L_H, L_Dc}|| and the norm ||L_H|| both
scale linearly with the number of bonds (for orthogonal bond
contributions) or with a topology-dependent factor, but the ratio
cancels the scaling.

Specifically:
- ||L_H||^2 = sum over bonds of ||L_H^bond||^2 + cross terms
- ||{L_H, L_Dc}||^2 = sum over bonds of ||{L_H^bond, L_Dc}||^2 + cross terms

The cross terms between different bonds contribute to both numerator and
denominator in the same proportion, preserving the ratio. The underlying
reason is that L_Dc is diagonal in the Pauli string basis and each bond
generates transitions with the same w_XY statistics relative to L_Dc.

---

## Summary table

| Statement | Result |
|-----------|--------|
| Cross-term is topology-specific | **Refuted** (topology-universal) |
| Cross-term depends on edge count | **Refuted** (independent of edges) |
| Complete graph restores orthogonality at N=3 | **Refuted** (same constant) |
| Cross-term is gamma-independent | **Confirmed** (all topologies, all N) |
| N=3 value = 1/sqrt(48) | **Confirmed** (all topologies) |
| N=4 value = 1/sqrt(128) | **Confirmed** (all topologies) |
| N-scaling: 1/sqrt(N * 2^(N+1)) | **Conjecture** (verified N=3, N=4) |

---

## What this establishes

1. **The Pythagorean breaking is a property of N, not of topology.**
   The angle between L_H and L_Dc in the operator norm depends only on
   the number of qubits and the coupling/noise type. The graph structure
   is irrelevant.

2. **1/sqrt(48) is not chain-specific.** The value reported in
   PRIMORDIAL_QUBIT_ALGEBRA for the N=3 chain is universal across all
   Heisenberg topologies. The original question ("is 1/sqrt(48)
   chain-specific?") has a definitive negative answer.

3. **A candidate formula exists.** The conjecture 1/sqrt(N * 2^(N+1))
   matches all measured data. It predicts N=5: 1/sqrt(320) = 0.05590...
   and N=6: 1/sqrt(768) = 0.03608..., both testable.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [simulations/cross_term_topology.py](../simulations/cross_term_topology.py) |
| Log | [simulations/results/cross_term_topology/cross_term_topology.txt](../simulations/results/cross_term_topology/cross_term_topology.txt) |
| Data | [simulations/results/cross_term_topology/cross_term_topology.json](../simulations/results/cross_term_topology/cross_term_topology.json) |
| Plot | [simulations/results/cross_term_topology/cross_term_topology.png](../simulations/results/cross_term_topology/cross_term_topology.png) |

---

*The question was: does the graph matter? The answer: no. The angle
between oscillation and cooling is set by N alone. Whether the qubits
form a line, a star, a ring, or a complete web, the cross-term is the
same. The geometry that breaks the Pythagorean theorem is not the
geometry of the graph. It is the geometry of the Pauli algebra at N
sites.*

*Thomas Wicht, Claude (Anthropic), April 13, 2026*
