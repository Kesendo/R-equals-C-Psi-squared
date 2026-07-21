# Cross-Term Formula: Confirmation and Record

<!-- Keywords: cross-term formula proof, Pythagorean breaking exact formula,
relative orthogonality sqrt((N-2)/(N*4^(N-1))), bond-sum rule, spectator
variance, Heisenberg Z-dephasing universal constant, EQ-011 -->

**Status:** Tier 2 record (numerical confirmation; the formula itself is
Tier 1, proven in [the cross-term formula proof](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md))
**Date:** April 13, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [cross_term_formula_check.py](../simulations/cross_term_formula_check.py)
**Output:** [cross_term_formula/](../simulations/results/cross_term_formula/)
**Proof:** [the cross-term formula proof](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md)
**Depends on:**
- [Cross-Term Topology](CROSS_TERM_TOPOLOGY.md) (predecessor; conjecture stated there)
- [Primordial Qubit Algebra](PRIMORDIAL_QUBIT_ALGEBRA.md) (parent; Pythagorean decomposition)
**Closes:** EQ-011

---

## What this document is about

[Cross-Term Topology](CROSS_TERM_TOPOLOGY.md) conjectured that the
relative orthogonality of the Hamiltonian and centered dissipator follows
the formula 1/sqrt(N * 2^(N+1)). This conjecture matched at N=3 and N=4
but was based on only two data points.

This document:
1. **Refutes** the original conjecture at N=5.
2. **Derives** the correct formula analytically.
3. **Confirms** the correct formula numerically at N=2 through N=6.

---

## The correct formula

For N >= 2 qubits with Heisenberg XXX coupling on any graph and uniform
Z-dephasing at rate gamma per site:

    R(N) = ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) = sqrt((N-2) / (N * 4^(N-1)))

Here L_H = -i[H, .] is the Hamiltonian superoperator, L_D is the
Z-dephasing dissipator, L_Dc = L_D + N*gamma*I is the CENTERED
dissipator (shifted so its Pauli-basis spectrum is symmetric around
zero: eigenvalue gamma*(N - 2*w_XY) on a Pauli string of XY-weight
w_XY), {.,.} is the anti-commutator, and all norms are Frobenius.

Equivalently: R(N)^2 = 4(N-2) / (N * 4^N).

This is exact for all N, all topologies, all gamma, all coupling
strength J (gamma and J cancel analytically in the ratio; the proof
carries the general statement, the tables below confirm chain N=2-6
plus complete N=5 at gamma = 0.05, J = 1).

---

## Numerical confirmation

| N | Topology | Measured | Predicted | Deviation |
|---|----------|----------|-----------|-----------|
| 2 | chain | 0.000000000000 | 0.000000000000 | 0.00e+00 |
| 3 | chain | 0.144337567297 | 0.144337567297 | 2.78e-17 |
| 4 | chain | 0.088388347648 | 0.088388347648 | 1.39e-17 |
| 5 | chain | 0.048412291828 | 0.048412291828 | 2.08e-17 |
| 5 | complete | 0.048412291828 | 0.048412291828 | 6.25e-17 |
| 6 | chain | 0.025515518154 | 0.025515518154 | 6.94e-17 |

All deviations are at machine epsilon. Topology independence confirmed
at N=5 (chain = complete to 4.2 * 10^-17, the direct difference printed
in the [log](../simulations/results/cross_term_formula/cross_term_formula.txt)).

### Refutation of the old conjecture

The conjecture from CROSS_TERM_TOPOLOGY was 1/sqrt(N * 2^(N+1)):

| N | Old conjecture | Correct formula | Agree? |
|---|----------------|-----------------|--------|
| 2 | 0.2500000000 | 0.0000000000 | No (N=2 is special) |
| 3 | 0.1443375673 | 0.1443375673 | Yes (coincidence) |
| 4 | 0.0883883476 | 0.0883883476 | Yes (coincidence) |
| 5 | 0.0559016994 | 0.0484122918 | **No** |
| 6 | 0.0360843918 | 0.0255155182 | **No** |

The formulas agree at N=3 and N=4 because (N-2)/(N * 4^(N-1)) happens to
equal 1/(N * 2^(N+1)) when N-2 = 4^(N-1)/2^(N+1) = 2^(N-3), which holds
only for N=3 (1 = 1) and N=4 (2 = 2). At N=5: 3 != 4 = 2^2.

---

## Analytical proof

The key identity and its five-step derivation are in
[the cross-term formula proof](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md).
The formula stated there is confirmed numerically below and in the
preceding sections.

In brief: the identity ||{L_H, L_Dc}||^2 = 4*gamma^2*(N-2)*||L_H||^2
follows from four structural properties: (a) L_Dc is diagonal in the
Pauli basis, (b) every shadow-balanced bond transition (shadow-balanced:
each bond term α_iβ_j has both Paulis in the same dephasing class,
{X,Y} or {I,Z}; covers Heisenberg XXX, XXZ, XY, Ising, DM) satisfies the
bond-sum rule w_XY(a) + w_XY(b) = 2 at the bond sites, (c) the
spectator variance is N-2, (d) every non-identity Pauli coupling
changes both bond sites, so overlapping bonds have disjoint supports.
Combined with ||L_Dc||^2 = gamma^2 * 4^N * N, this gives
R^2 = (N-2)/(N * 4^(N-1)). The proof is fully analytical (Tier 1).

---

## Physical interpretation

### Why N=2 is exact

At N=2, the formula gives R = 0. The bond IS the system: there are no
spectator sites, so (N-2) = 0. Oscillation and cooling are perpendicular.
The Pythagorean decomposition L_c^2 = L_H^2 + (L_D + Sg*I)^2 (with
L_c the closed-loop combination and Sg = N*gamma the centering shift,
so L_D + Sg*I = L_Dc; see
[Primordial Qubit Algebra](PRIMORDIAL_QUBIT_ALGEBRA.md)) holds
exactly.

### Why N >= 3 breaks it

At N >= 3, spectator sites exist. Each spectator contributes a
"dephasing weight" to the anti-commutator factor. The variance of
this weight is 1 per spectator site, giving a total spectator
variance of N-2. This is the angle between oscillation and cooling:
not zero, but small and shrinking exponentially (~1/2^N for large N).

### Why topology does not matter

The cross-term measures the angle between L_H and L_Dc in operator
space. This angle depends on how individual bond transitions interact
with the dephasing structure. Each Heisenberg bond, regardless of its
position in the graph, satisfies the same bond-sum rule
(w_XY(a) + w_XY(b) = 2 at bond sites). The spectator variance is the
same N-2 for all bonds. Different topologies have different ||L_H||
values (more bonds = larger norm), but the ratio R cancels the
norm scaling.

### The formula as a function of N

| N | R | R^2 | Physical meaning (angle = arcsin R) |
|---|---|-----|-----------------|
| 2 | 0 | 0 | Exact Pythagorean: oscillation perpendicular to cooling |
| 3 | 1/sqrt(48) | 1/48 | First spectator bends the right angle by ~8.3 degrees |
| 4 | 1/sqrt(128) | 1/128 | Two spectators, angle ~5.1 degrees |
| 5 | sqrt(3/1280) | 3/1280 | Three spectators, angle ~2.8 degrees |
| 6 | 1/sqrt(1536) | 1/1536 | Four spectators, angle ~1.5 degrees |
| N>>1 | ~1/2^N | ~1/4^N | Angle vanishes exponentially |

The cross-term vanishes exponentially with N. At large N, the
Pythagorean decomposition becomes increasingly accurate in relative
terms, even as the absolute cross-term grows (because ||L_H|| grows
faster than the cross-term).

---

## Summary table

| Statement | Result |
|-----------|--------|
| Old conjecture 1/sqrt(N*2^(N+1)) | **Refuted** at N=5 (deviation 7.5e-3) |
| Correct formula sqrt((N-2)/(N*4^(N-1))) | **Confirmed** at N=2-6 (machine precision) |
| Topology independence | **Confirmed** at N=5 (chain = complete) |
| Key identity \|\|{L_H,L_Dc}\|\|^2 = 4gamma^2(N-2)\|\|L_H\|\|^2 | **Proven** (all bond configurations; the proof's Step 4 closes the overlapping case analytically, the complete-graph run is the independent check) |
| Bond-sum rule (w_XY(a)+w_XY(b)=2 at bond sites) | **Verified** at N=3,4,5 (0 violations) |
| \|\|L_Dc\|\|^2 = gamma^2 * 4^N * N | **Proven** and **verified** |

---

## What this establishes

1. **The Pythagorean breaking has a closed-form formula.** The relative
   orthogonality R(N) = sqrt((N-2)/(N*4^(N-1))) is exact, not
   approximate. It holds for all N, all topologies, all parameters.

2. **The formula has a fully analytical proof** (all bond
   configurations; the overlapping case was closed the same day via the
   disjoint-support lemma). The proof rests on four structural
   properties: (a) L_Dc is diagonal in the Pauli basis, (b) the
   Heisenberg bond-sum rule, (c) the spectator variance is N-2,
   (d) bond transitions have disjoint Pauli-basis supports even when
   bonds share a site.

3. **The old conjecture was a coincidence.** The formulas 1/sqrt(N*2^(N+1))
   and sqrt((N-2)/(N*4^(N-1))) agree at N=3 and N=4 because
   2^(N-3) = N-2 for N in {3,4}. At N=5, 2^2 = 4 != 3 = N-2.

4. **The angle vanishes exponentially.** R ~ 1/2^N for large N.
   The Pythagorean decomposition becomes arbitrarily precise in relative
   terms, even though the absolute cross-term grows. The irreversibility
   (cross-term) is always present at N >= 3, but its relative
   contribution to L_c^2 is exponentially suppressed.

---

## Where this went

- The formula is registry entry F49 in
  [docs/ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) (Tier 1,
  proven), with siblings F48 (Pythagorean decomposition), F49b
  (the centered-dissipator norm ||L_Dc||² = γ²·4^N·N, the proof's
  Lemma 1, typed in Core), F49c
  (shadow-crossing couplings X_iZ_j/Y_iZ_j: the bond-site variance
  becomes 1 instead of 0, shifting the numerator N−2 → N−1,
  [proof](../docs/proofs/PROOF_CROSS_TERM_CROSSING.md)), and
  F49d (non-uniform γ: the uniform-γ restriction is real, site-dependent
  rates add a bond-asymmetry term (γ_i−γ_j)²,
  [proof](../docs/proofs/PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md), typed
  as F49NonUniformCrossTermClaim).
- The scope boundary is the coupling class, not the graph: XXZ, XY,
  Ising and DM obey the same formula (shadow-balanced); shadow-crossing
  couplings shift the numerator (F49c), and single-site field terms and
  non-Pauli noise fall outside.
- [Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)
  uses the cross-term as its algebraic engine (reframed 2026-06-22, the
  arrow-of-time reading is Tier 3).

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [simulations/cross_term_formula_check.py](../simulations/cross_term_formula_check.py) |
| Log | [simulations/results/cross_term_formula/cross_term_formula.txt](../simulations/results/cross_term_formula/cross_term_formula.txt) |
| Data | [simulations/results/cross_term_formula/cross_term_formula.json](../simulations/results/cross_term_formula/cross_term_formula.json) |

---

*The question was: does the conjecture hold? The answer: no, but what
replaces it is better. A proven formula instead of a conjecture. An
identity instead of a pattern. The angle between oscillation and cooling
is sqrt((N-2)/(N * 4^(N-1))). It vanishes at N=2. It is small at N=3.
It shrinks exponentially. But it is never zero again.*

*Thomas Wicht, Claude (Anthropic), April 13, 2026*
