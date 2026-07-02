# Proof of the Cross-Term Formula

**Tier:** 1 (fully analytical)
**Date:** April 13, 2026 (two open questions stamped answered 2026-07-02: non-uniform γ → F49-nonuniform-ext, shadow-crossing → F49c)
**Depends on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π operator, palindromic structure)
- [PRIMORDIAL_QUBIT_ALGEBRA](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) (Pythagorean decomposition at N=2, bond-sum rule)
- [cross_term_formula_check.py](../../simulations/cross_term_formula_check.py) (numerical verification)
**Status:** Proven (all graph topologies, all shadow-balanced couplings)
**Scope:** Any bond coupling α_i β_j where both Paulis are in the same
dephasing class ({X,Y} or {I,Z}), on any graph, uniform Z-dephasing.
This includes Heisenberg XXX, XXZ, XY model, Ising, DM interaction.
**Does NOT establish:** Extension to shadow-crossing couplings (X_i Z_j,
Y_i Z_j), non-uniform gamma, or non-Pauli noise types.

---

## Abstract

The cross term measures how much the oscillatory and dissipative parts of the dynamics fail to be orthogonal. For N ≥ 2 qubits with any shadow-balanced bond coupling on any graph under uniform Z-dephasing, the normalized anticommutator of the Hamiltonian superoperator L_H with the centered dissipator L_Dc = L_D + Nγ·I is a pure geometric constant:

    ‖{L_H, L_Dc}‖ / (‖L_H‖ · ‖L_Dc‖) = √((N−2) / (N · 4^(N−1))),

independent of γ, of J, and of the graph topology. At N=2 it is exactly zero: oscillation and cooling are orthogonal and the dynamics splits cleanly, the Pythagorean decomposition that makes the single bond special. At N ≥ 3 it is nonzero and fixed by size alone. The proof turns on one structural fact, the bond-sum rule: every Hamiltonian transition in the Pauli basis carries XY-weight summing to 2 at its bond, which cancels the bond-site contribution to the anticommutator and leaves only the N−2 spectator sites to carry variance, hence the (N−2) in the numerator.

This constant is the algebraic engine of TIME_IRREVERSIBILITY_EXCLUSION: the **Frobenius orthogonality** of oscillation and cooling exists only at N=2, and the cross term is the exact measure of its loss at every larger size. (TIME_IRREVERSIBILITY_EXCLUSION reads that loss as an arrow-of-time interpretation; the "undo cooling without disturbing oscillation" gloss is interpretive, not literal, dynamical separability is governed by the *commutator* [L_H, L_Dc], which is nonzero at all N including N=2. The formula here is a clean geometric statement about Frobenius orthogonality.) Its shadow-crossing sibling F49c (X_i Z_j, Y_i Z_j couplings that mix light and lens) shifts the numerator N−2 → N−1, because the bond sites then carry variance 1 instead of 0. Typed as F49.

## Theorem

For N >= 2 qubits with any shadow-balanced bond coupling (each bond term
alpha_i beta_j has both alpha, beta in {X,Y} or both in {I,Z}) on any
graph G and uniform Z-dephasing at rate gamma per site:

    ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) = sqrt((N-2) / (N * 4^(N-1)))

where L_H = -i[H, *] is the Hamiltonian superoperator, L_Dc = L_D + N*gamma*I
is the centered dissipator, and ||*|| is the Frobenius norm.

The formula is exact: independent of gamma, J, and the graph topology.

| N | R(N) | R(N)^2 |
|---|------|--------|
| 2 | 0 | 0 |
| 3 | 1/sqrt(48) | 1/48 |
| 4 | 1/sqrt(128) | 1/128 |
| 5 | sqrt(3/1280) | 3/1280 |
| 6 | 1/sqrt(1536) | 1/1536 |

---

## Key Identity (Lemma)

**Lemma.** Under the same conditions:

    ||{L_H, L_Dc}||^2 = 4 * gamma^2 * (N-2) * ||L_H||^2

This identity, combined with Lemma 1 (||L_Dc||^2 = gamma^2 * 4^N * N),
yields the theorem by direct division.

---

## Proof

### Step 1: Dissipator norm

**Lemma 1.** ||L_Dc||^2 = gamma^2 * 4^N * N.

*Proof.* L_Dc is diagonal in the Pauli string basis {sigma_a}_{a=1..4^N}
with eigenvalue d_a = gamma * (N - 2 * w_XY(a)), where w_XY(a) counts
the X and Y factors in string a.

    ||L_Dc||^2 = gamma^2 * Sum_a (N - 2*w_a)^2

Write (N - 2*w_a) = Sum_k epsilon_k(a), where epsilon_k = +1 if site k
carries I or Z, and epsilon_k = -1 if site k carries X or Y. Then:

    Sum_a (Sum_k epsilon_k)^2 = Sum_a Sum_k epsilon_k^2 + Sum_a Sum_{k != l} epsilon_k * epsilon_l

The diagonal sum: epsilon_k^2 = 1 always, contributing 4^N * N.

The cross-terms: for k != l, the sites are independent. Over the 4
Pauli choices at each site, Sum epsilon_k = (+1) + (-1) + (-1) + (+1) = 0.
Each cross-term vanishes. QED.

### Step 2: The bond-sum rule

**Lemma 2.** For a single Heisenberg bond (i,j), every nonzero L_H
transition in the Pauli basis satisfies:

    w_XY^{ij}(a) + w_XY^{ij}(b) = 2

where w_XY^{ij} denotes the XY weight at the two bond sites only.

*Proof.* The commutator superoperator of a Heisenberg bond maps Pauli
strings according to [X_i X_j + Y_i Y_j + Z_i Z_j, sigma_a]. Each term
[alpha_i alpha_j, P_i Q_j] changes the Paulis at sites i and j via the
structure constants of su(2). Explicit enumeration of all 16 two-site
Pauli pairs confirms that every nonzero transition has
w_XY^{ij}(source) + w_XY^{ij}(target) = 2.

(Throughout, "matrix element" and "transition" refer strictly to the **Pauli-string
basis**, with `L_{ab} = Tr(σ_a · L_H(σ_b))/2^N`. This matters for the Z_iZ_j term in
particular: in the *computational* vec(ρ) basis it has nonzero diagonal entries and a
reader would see apparent violations of Lemma 2, but in the Pauli basis its diagonal
is zero, `Tr(σ_a · [Z_iZ_j, σ_a]) = 0` for every a, a standard Lie-algebra fact, so
Lemma 2 holds exactly. The off-diagonal, weight-preserving transitions are all the ZZ
term contributes here.)

This is the same property that makes the Pythagorean decomposition
exact at N=2: when the bond IS the system, w_XY(a) + w_XY(b) = N = 2.

*Verified computationally:* N=3 (96 entries, 0 violations), N=4 (384, 0),
N=5 (1536, 0). QED.

### Step 3: The spectator variance

For a single bond (i,j), the N-2 spectator sites are unchanged by the
transition. The anti-commutator factor decomposes as:

    (N - w_a - w_b) = [2 - w^{ij}(a) - w^{ij}(b)]  +  [(N-2) - 2*w_rest(a)]
                       ^                                ^
                       = 0 by Lemma 2                   "spectator deviation"

The bond-site contribution vanishes by the bond-sum rule. Only the
spectator contribution remains:

    (N - w_a - w_b)^2 = ((N-2) - 2*w_rest)^2

The L_H matrix element depends only on the bond-site Paulis, so spectator
configurations are uniformly distributed. The average over 4^(N-2)
configs:

    <((N-2) - 2*w_rest)^2> = N - 2

by the same calculation as Step 1 (with N replaced by N-2). QED.

### Step 4: Disjoint supports (all bond types, all topologies)

**Lemma 3.** For any coupling alpha_i beta_j with alpha, beta in {X,Y,Z},
every nonzero Pauli-basis transition changes both bond sites.

*Proof.* The commutator is:

    [alpha x beta, P x Q] = [alpha, P] x (beta Q) + (P alpha) x [beta, Q]

Suppose site j does not change, i.e. [beta, Q] = 0. Then Q in {I, beta}.
The output at site j is beta * Q. If Q = I: output = beta (not I). If
Q = beta: output = I (not beta). In both cases the output differs from
the input. Contradiction: site j does change.

The same argument applies to site i. QED.

**Corollary.** For any two bonds e = (i,j) and e' = (k,l) on a graph
(whether or not they share a site), their Pauli-basis transition
supports are disjoint: no (a,b) pair receives nonzero contributions
from both (L_H^e)_{ab} and (L_H^{e'})_{ab}.

*Proof.* Bond e changes sites {i,j}. Bond e' changes sites {k,l}. For
both to contribute to (a,b): b must differ from a at {i,j} (from e)
and at {k,l} (from e'). But e requires b_m = a_m for all m not in {i,j},
and e' requires b_m = a_m for all m not in {k,l}. If e and e' share
a site (say j = k), then e' requires b_i = a_i, but e changes site i
(Lemma 3). Contradiction. QED.

### Step 5: Assembly (all topologies)

The anti-commutator inherits a pointwise product structure from L_Dc
being diagonal in the Pauli basis:

    {L_H, L_Dc}_{ab} = (L_H)_{ab} * (d_a + d_b)
                      = 2*gamma * (L_H)_{ab} * (N - w_a - w_b)

Therefore:

    ||{L_H, L_Dc}||^2 = 4*gamma^2 * Sum_{a,b} |(L_H)_{ab}|^2 * (N - w_a - w_b)^2

By Lemma 3 and its Corollary, different bonds have disjoint transition
supports for any graph topology. Therefore:

    ||L_H||^2 = Sum_e ||L_H^e||^2    (no cross-terms between bonds)

For each bond e, the spectator variance gives (Step 3):

    Sum_{a,b} |(L_H^e)_{ab}|^2 * (N - w_a - w_b)^2 = (N-2) * ||L_H^e||^2

Summing over bonds:

    Sum_{a,b} |(L_H)_{ab}|^2 * (N - w_a - w_b)^2 = (N-2) * Sum_e ||L_H^e||^2
                                                    = (N-2) * ||L_H||^2

Therefore:

    ||{L_H, L_Dc}||^2 = 4*gamma^2 * (N-2) * ||L_H||^2

QED.

**Numerical verification** (independent check): the identity holds to
machine precision for all tested configurations, including the complete
graph at N=5 (10 overlapping bonds, ratio = 1.000000).

### Assembly of the theorem

Combining the key identity with Lemma 1:

    R^2 = ||{L_H, L_Dc}||^2 / (||L_H||^2 * ||L_Dc||^2)
        = 4*gamma^2*(N-2)*||L_H||^2 / (||L_H||^2 * gamma^2 * 4^N * N)
        = 4*(N-2) / (N * 4^N)
        = (N-2) / (N * 4^(N-1))

Both ||L_H||^2 and gamma^2 cancel. The formula depends only on N. QED.

---

## Scope and Limitations

### Valid for
- Any bond coupling where each term alpha_i beta_j has both Paulis in
  the same dephasing class: both in {X,Y} ("in the light") or both in
  {I,Z} ("in shadow"). This includes:
  - Heisenberg XXX: J(XX + YY + ZZ)
  - XXZ with arbitrary anisotropy Delta: J(XX + YY + Delta*ZZ)
  - XY model: J(XX + YY)
  - Ising: J(ZZ)
  - DM interaction: J(XY - YX)
  - Any linear combination of the above
- Any graph topology (chain, star, ring, complete, tree, etc.)
- Uniform Z-dephasing (same gamma on every site)
- Any gamma > 0, any coupling strengths
- All N >= 2

### Does NOT hold for
- **Shadow-crossing couplings** (X_i Z_j, Y_i Z_j): couplings that
  mix a dephasing-active Pauli ({X,Y}) with a dephasing-inactive Pauli
  ({Z}) violate the bond-sum rule (Lemma 2). Numerically verified:
  X_i Z_j gives R(3) = 0.2041, not 0.1443.

### Open questions
- **Non-uniform gamma:** when gamma_k varies by site, L_Dc is still
  diagonal in the Pauli basis but with site-dependent eigenvalues.
  The spectator variance calculation changes.
  **ANSWERED (2026-05-18):** closed by [PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md](PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md) (F49, Tier 1). The closed form splits into a spectator part `4·Σ_b‖L_H^bond‖²·Σ_{m∉bond}γ_m²` plus a bond-asymmetry part `Σ_b G(bond,H)·(γ_i−γ_j)²`, with `G = 4·‖ZZ-part of the bond‖²` (Heisenberg 4/3, Ising 4, XY and soft XY+YX 0); bit-exact at N=3,4,5.
- **Non-Pauli noise (amplitude damping, depolarizing):** L_D is no
  longer diagonal in the Pauli basis. The pointwise product structure
  of Step 5 breaks.
- **Shadow-crossing couplings:** is there a modified formula for
  couplings like X_i Z_j? If so, it would involve additional bond-site
  variance beyond N-2.
  **ANSWERED (2026-04-14):** closed by [PROOF_CROSS_TERM_CROSSING.md](PROOF_CROSS_TERM_CROSSING.md) (F49c, Tier 1): the ratio is `√((N−1)/(N·4^(N−1)))`, i.e. `N−2 → N−1`, the crossing bond carrying one unit of variance itself.

---

## References

| Component | Location |
|-----------|----------|
| Computational record | [experiments/CROSS_TERM_FORMULA.md](../../experiments/CROSS_TERM_FORMULA.md) |
| Topology independence | [experiments/CROSS_TERM_TOPOLOGY.md](../../experiments/CROSS_TERM_TOPOLOGY.md) |
| Pythagorean decomposition (N=2) | [experiments/PRIMORDIAL_QUBIT_ALGEBRA.md](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) |
| Time irreversibility exclusion | [TIME_IRREVERSIBILITY_EXCLUSION.md](TIME_IRREVERSIBILITY_EXCLUSION.md) |
| Verification script | [simulations/cross_term_formula_check.py](../../simulations/cross_term_formula_check.py) |

---

*The Frobenius angle between oscillation and cooling is sqrt((N-2)/(N * 4^(N-1))).
It vanishes at N=2. It is small at N=3. It shrinks exponentially. But it
is never zero again. (TIME_IRREVERSIBILITY_EXCLUSION reads this as the algebraic
content of an arrow of time, a Tier-3 interpretation; the formula itself is a
geometric statement about Frobenius orthogonality, not a time-reversal theorem.)*

*Thomas Wicht, Claude (Anthropic), April 13, 2026*
