# Proof of the Cross-Term Formula

**Tier:** 1 (fully analytical)
**Date:** April 13, 2026 (two open questions stamped answered 2026-07-02: non-uniform γ → F49-nonuniform-ext, shadow-crossing → F49c)
**Depends on:**
- [Mirror Symmetry](MIRROR_SYMMETRY_PROOF.md) (Π operator, palindromic structure)
- [Primordial Qubit Algebra](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) (Pythagorean decomposition at N=2, bond-sum rule)
- [cross_term_formula_check.py](../../simulations/cross_term_formula_check.py) (numerical verification)
**Status:** Proven (all graph topologies, all shadow-balanced couplings)
**Scope:** Any bond coupling α_i β_j with both letters in {X,Y} or both
equal to Z (a leg on I makes the term a field, which the section
*Per-letter rates* treats), on any graph, uniform Z-dephasing.
This includes Heisenberg XXX, XXZ, XY model, Ising, DM interaction.
**Does NOT establish:** amplitude damping, the one noise considered here
that is not diagonal in the Pauli basis. Shadow-crossing couplings (X_i Z_j,
Y_i Z_j) are F49c, non-uniform γ is F49d, and light along X, Y and Z per
site, depolarizing included, is the section *Per-letter rates* below.

---

## Abstract

The cross term measures how much the oscillatory and dissipative parts of the dynamics fail to be orthogonal. For N ≥ 2 qubits with any shadow-balanced bond coupling on any graph under uniform Z-dephasing, the normalized anticommutator of the Hamiltonian superoperator L_H with the centered dissipator L_Dc = L_D + Nγ·I is a pure geometric constant:

    ‖{L_H, L_Dc}‖ / (‖L_H‖ · ‖L_Dc‖) = √((N−2) / (N · 4^(N−1))),

independent of γ, of J, and of the graph topology. At N=2 it is exactly zero: oscillation and cooling are orthogonal and the dynamics splits cleanly, the Pythagorean decomposition that makes the single bond special. At N ≥ 3 it is nonzero and fixed by size alone. The proof turns on one structural fact, the bond-sum rule: every Hamiltonian transition in the Pauli basis carries XY-weight summing to 2 at its bond, which cancels the bond-site contribution to the anticommutator and leaves only the N−2 spectator sites to carry variance, hence the (N−2) in the numerator.

This constant is the algebraic engine of TIME_IRREVERSIBILITY_EXCLUSION: the **Frobenius orthogonality** of oscillation and cooling exists only at N=2, and the cross term is the exact measure of its loss at every larger size. (TIME_IRREVERSIBILITY_EXCLUSION reads that loss as an arrow-of-time interpretation; the "undo cooling without disturbing oscillation" gloss is interpretive, not literal, dynamical separability is governed by the *commutator* [L_H, L_Dc], which is nonzero at all N including N=2. The formula here is a clean geometric statement about Frobenius orthogonality.) Its shadow-crossing sibling F49c (X_i Z_j, Y_i Z_j couplings that mix light and lens) shifts the numerator N−2 → N−1, because the bond sites then carry variance 1 instead of 0. Typed as F49.

## Theorem

For N >= 2 qubits with any shadow-balanced bond coupling (each bond term
alpha_i beta_j has both alpha, beta in {X,Y} or both equal to Z) on any
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

## Per-letter rates: light along X, Y and Z

The proof uses one property of Z-dephasing: L_D is diagonal in the Pauli basis. Dephasing along any letter has it, so the proof carries over with one change, the rate a site charges. Let site l be lit along each letter C ∈ {X, Y, Z} at rate γ_{l,C} ≥ 0, so that L_D = Σ_{l,C} γ_{l,C}·D[C_l], and centre it as L_Dc = L_D + Γ·I with Γ = Σ_{l,C} γ_{l,C}. Depolarizing is the case of equal letter rates on each site, in either of the repo's conventions (γ/3 per letter as in F5, or γ per letter as in the absorption theorem), with the site rates free to differ. Write χ_C(P) = +1 when the letter P commutes with C and −1 when it does not.

**Step 1′ (the centred rate is a character sum).** L_Dc is diagonal on strings, d_a = Σ_l e_l(a_l), with

    e_l(P) = Σ_C γ_{l,C} · χ_C(P).

For Z light alone this is Step 1's γ·ε_k. The three characters are orthogonal over the four letters and each sums to zero there, so

    ||L_Dc||^2 = 4^N · Σ_l g_l^2,    g_l^2 = Σ_C γ_{l,C}^2.

**Step 2′ (the moving letter picks its own light).** A term τ = c·⊗_l a_l acts on a string only where it anticommutes with it, and then moves the letter on each site of its support by a_l: P ↦ a_l·P up to phase. Since χ_C(a·P) = χ_C(a)·χ_C(P), and χ_C(a) = +1 only for C = a,

    e_l(P) + e_l(a·P) = Σ_C γ_{l,C} · χ_C(P) · (1 + χ_C(a)) = 2·γ_{l,a}·χ_a(P).

On each site of the support only the light along the moving letter is felt. For Z light this is Lemma 2's bond-sum rule, and it is the A ∈ {−2, 0, +2} of [the non-uniform γ extension](PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md).

**Step 3′ (assembly).** The two strings a transition joins fix its term, since their product is the term's string, so different terms have disjoint transition supports; the Corollary to Lemma 3 is the case of two bonds. Each term's ||L_H^τ||^2 = 2·4^N·|c_τ|^2 is spread evenly over its transitions. A spectator site contributes 2·e_m(P_m), of mean zero and mean square 4·g_m^2, as in Step 3. On the support write s_l = χ_{a_l}(P_l): a transition needs an odd number of anticommuting sites, ∏_l s_l = −1, and the support contributes the mean of (2·Σ_l γ_{l,a_l}·s_l)^2 over the transitions. Therefore

    ||{L_H, L_Dc}||^2 = Σ_τ ||L_H^τ||^2 · [ 4·Σ_{m ∉ supp τ} g_m^2 + B_τ ],

    B_τ = 4·γ_{i,a}^2                  for a field a_i              (s_i = −1),
          4·(γ_{i,a} − γ_{j,b})^2      for a two-site term a_i b_j  (s_i·s_j = −1),
          4·Σ_{l ∈ supp τ} γ_{l,a_l}^2  for three or more sites      (the s_l pairwise uncorrelated).

**What it contains.**

- *Uniform Z light* (γ_{l,Z} = γ, the rest zero). A two-site term with both letters in {X, Y}, or both Z, has B = 0, and the theorem follows. X_iZ_j has B = 4γ^2, one unit of variance on the bond beside the spectators' N − 2: F49c's N − 2 → N − 1. A bond that carries terms of both kinds is the sum of its terms, which settles the mixed Hamiltonians [the shadow-crossing proof](PROOF_CROSS_TERM_CROSSING.md) left open.
- *Per-site Z rates.* ZZ terms give 4(γ_i − γ_j)^2 and the terms with both letters in {X, Y} nothing, which is F49d; a crossing term X_iZ_j gives 4γ_j^2, the case F49d left out of scope.
- *Uniform depolarizing* (the same rate for every letter on every site). Then γ_{i,a} = γ_{j,b} for any two letters and sites, so every two-site term has B = 0, ||{L_H, L_Dc}||^2 = 4(N − 2)·g^2·||L_H||^2 and ||L_Dc||^2 = 4^N·N·g^2, and

      R(N) = sqrt((N − 2) / (N · 4^(N−1)))

  for every two-site coupling on any graph, X_iZ_j included. Under equal light no letter is in the shade, so no coupling crosses anything. Per-letter site rates p_l that differ give every two-site term B = 4(p_i − p_j)^2 whatever its letters (a unit-coupling chain at N = 3 with p = (1, 2, 3): R^2 = 1/42 for XX, XZ and ZZ alike, against 1/48), and a field carries B = 4γ_{i,a}^2 > 0 in any case.
- *N = 2.* Uniform depolarizing, centred at its own total rate Γ, keeps the Pythagorean decomposition ([F48](../ANALYTICAL_FORMULAS.md)) exact for every two-site bond.

On [the letter cube](../THE_ONE_SQUARE.md) Step 1′ is the centred light formula: the cube's Q_P = Σ_l P_l(·)P_l has eigenvalue N − 2k_P on a string, so F49's L_Dc is γ·Q_Z and with the same letter rates on every site L_Dc = Σ_C γ_C·Q_C. Step 2′ reads the corner before and after a move: a move by a keeps the coordinate k_a and flips the other two, and in the sum of the two rates the flipped ones cancel, which leaves the light along a.

**What the repo held.** A sweep of the registry, the proofs, the experiments, the OpenArcs registry, the caught-errors ledger, the glossary, the Core claims (all F49 ones), the Diagnostics witnesses and MirrorWorld's formulas found no cross-term statement for light along Y, for mixed or per-letter light, for depolarizing, or for the character sum. Light along X alone is there: the Hadamard notes of the Core's F49b and F49c claims ([the bit_a twin via Hadamard](PROOF_BIT_A_TWIN_VIA_HADAMARD.md)) take Z-dephasing to X-dephasing with the norms, the uniform X case of this section; the cross term's Z-light home is F49NonUniformCrossTermClaim. (The Core's F49Pi2Inheritance types a different object under the same number, the scaling of the palindrome residual.) fw.Confirmations holds nothing on the cross term. The pieces nearest to the rest: Step 1's ε_k and the A-classification of the non-uniform extension (Z light only), the per-letter cost law of [the absorption theorem](PROOF_ABSORPTION_THEOREM.md) §2, the letter cube's centred Q_P with [the three diagonals](../THE_THREE_DIAGONALS.md)' L_D = γ·(Q_Z − N·I), and the remark in [the depolarizing palindrome](../../experiments/DEPOLARIZING_PALINDROME.md) that each site's four letter rates average to its rate sum, which is the centring used here. Keep this N − 2 apart from F5's depolarizing error, whose (N − 2) the ledger records as a mis-import from this formula. Adjacent: the commutator [L_H, L_Dc], the separability measure of [Time Irreversibility Exclusion](TIME_IRREVERSIBILITY_EXCLUSION.md), is Step 2′'s complement, since e_l(aP) − e_l(P) = −2·Σ_{C ≠ a} γ_{l,C}·χ_C(P) keeps the lights a move flips and loses the spectators; that is [the letter cube](../THE_ONE_SQUARE.md)'s account of which moves the watching feels.

**Verified** exactly, with integer data and no tolerance, by [`f49_per_letter_cross_term_gate.py`](../../simulations/f49_per_letter_cross_term_gate.py): random terms (fields, two-site, three- and four-site strings) with random integer rates per site and letter at N = 3 and 4, the anticommutator built from the Lindblad form itself; F49, F49c, F49d, mixed bonds and the per-site crossing case as special cases; all nine two-site letter pairs under uniform depolarizing on chains and rings at N = 3 to 5; per-site depolarizing rates and a field as the controls that do not reach F49's constant; and the identity of Step 2′.

---

## Scope and Limitations

### Valid for
- Any bond coupling where each term alpha_i beta_j has both Paulis in
  the same dephasing class: both in {X,Y} ("in the light") or both Z
  ("in shadow"); a leg on I makes the term a field, which the per-letter
  section treats. This includes:
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
  **ANSWERED (2026-05-18):** closed by [F49 non-uniform γ extension](PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md) (F49, Tier 1). The closed form splits into a spectator part `4·Σ_b‖L_H^bond‖²·Σ_{m∉bond}γ_m²` plus a bond-asymmetry part `Σ_b G(bond,H)·(γ_i−γ_j)²`, with `G = 4·‖ZZ-part of the bond‖²` (Heisenberg 4/3, Ising 4, XY and soft XY+YX 0); bit-exact at N=3,4,5.
- **Amplitude damping:** its L_D is not diagonal in the Pauli basis, so
  Step 5's pointwise product does not apply as it stands. On the letter
  cube it is two equal transverse lights, γ/4 along X and along Y, which
  the per-letter section covers, plus one move per site, I → Z; that move
  is what remains open. (Depolarizing is diagonal, a Pauli channel, and is
  answered in the per-letter section.)
- **Shadow-crossing couplings:** is there a modified formula for
  couplings like X_i Z_j? If so, it would involve additional bond-site
  variance beyond N-2.
  **ANSWERED (2026-04-14):** closed by [Cross-Term Crossing](PROOF_CROSS_TERM_CROSSING.md) (F49c, Tier 1): the ratio is `√((N−1)/(N·4^(N−1)))`, i.e. `N−2 → N−1`, the crossing bond carrying one unit of variance itself.

---

## References

| Component | Location |
|-----------|----------|
| Computational record | [Cross-Term Formula experiment](../../experiments/CROSS_TERM_FORMULA.md) |
| Topology independence | [Cross-Term Topology](../../experiments/CROSS_TERM_TOPOLOGY.md) |
| Pythagorean decomposition (N=2) | [Primordial Qubit Algebra](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) |
| Time irreversibility exclusion | [Time Irreversibility Exclusion](TIME_IRREVERSIBILITY_EXCLUSION.md) |
| Verification script | [simulations/cross_term_formula_check.py](../../simulations/cross_term_formula_check.py) |

---

*The Frobenius angle between oscillation and cooling is sqrt((N-2)/(N * 4^(N-1))).
It vanishes at N=2. It is small at N=3. It shrinks exponentially. But it
is never zero again. (TIME_IRREVERSIBILITY_EXCLUSION reads this as the algebraic
content of an arrow of time, a Tier-3 interpretation; the formula itself is a
geometric statement about Frobenius orthogonality, not a time-reversal theorem.)*

*Thomas Wicht, Claude (Anthropic), April 13, 2026*
