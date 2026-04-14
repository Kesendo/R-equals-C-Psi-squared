# Proof of the Cross-Term Formula for Shadow-Crossing Couplings

**Tier:** 1 (fully analytical)
**Date:** April 14, 2026
**Depends on:**
- [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md) (parent proof; Lemmas 1, 3 reused)
- [cross_term_crossing.py](../../simulations/cross_term_crossing.py) (numerical verification)
**Status:** Proven (all graph topologies, all shadow-crossing couplings)
**Scope:** Any bond coupling alpha_i beta_j with one Pauli in {X,Y} and
one in {I,Z}, on any graph, uniform Z-dephasing.
**Does NOT establish:** Mixed Hamiltonians combining shadow-balanced and
shadow-crossing terms in the same bond.

---

## Theorem

For N >= 2 qubits with any shadow-crossing bond coupling on any graph G
and uniform Z-dephasing at rate gamma per site:

    R(N) = sqrt((N-1) / (N * 4^(N-1)))

---

## Key Identity

    ||{L_H, L_Dc}||^2 = 4 * gamma^2 * (N-1) * ||L_H||^2

The only change from the balanced case (N-2) is the replacement N-2 -> N-1.

---

## Proof

Steps 1 (dissipator norm), 4 (disjoint supports), and the final assembly
are identical to [PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md).
Only Steps 2 and 3 change.

### Modified Step 2: Bond-site deviation

**Lemma 2'.** For a shadow-crossing coupling alpha_i beta_j (alpha in
{X,Y}, beta in {I,Z} or vice versa), the bond-site deviation
s = w_XY(a) + w_XY(b) - 2 satisfies:

    <s> = 0,    <s^2> = 1

weighted by |M|^2 over the nonzero Pauli-basis transitions.

*Proof (by enumeration for X_i Z_j).* The commutator [X tensor Z, P tensor Q]
has 8 nonzero transitions (out of 16 input pairs), all with |M|^2 = 4:

| Source | Target | w_src + w_tgt | s |
|--------|--------|---------------|---|
| IX | XY | 3 | +1 |
| IY | XX | 3 | +1 |
| XX | IY | 3 | +1 |
| XY | IX | 3 | +1 |
| YI | ZZ | 1 | -1 |
| YZ | ZI | 1 | -1 |
| ZI | YZ | 1 | -1 |
| ZZ | YI | 1 | -1 |

Four transitions with s = +1, four with s = -1, all equal weight.
Mean s = 0. Mean s^2 = 1. QED.

The same calculation holds for YZ, ZX, ZY by the cyclic symmetry of
the Pauli algebra. The bond-site deviation depends only on the
shadow-crossing structure, not on which specific Paulis are involved.

### Modified Step 3: Total variance

The anti-commutator factor decomposes as:

    (N - w_a - w_b) = s + t

where s is the bond-site deviation (<s> = 0, <s^2> = 1 by Lemma 2')
and t = (N-2) - 2*w_rest is the spectator deviation (<t> = 0,
<t^2> = N-2 by the parent proof's Step 3).

Since s and t are independent (bond sites and spectator sites are
disjoint, and s depends on the bond-site Paulis while t depends on
spectator-site Paulis):

    <(s + t)^2> = <s^2> + 2<s><t> + <t^2> = 1 + 0 + (N-2) = N-1

### Assembly

By the parent proof's Lemma 3 (disjoint supports, unchanged for
shadow-crossing couplings since both bond Paulis are non-identity):

    ||L_H||^2 = Sum_e ||L_H^e||^2

For each bond e, the total variance is N-1:

    ||{L_H, L_Dc}||^2 = 4*gamma^2 * (N-1) * ||L_H||^2

Combined with ||L_Dc||^2 = gamma^2 * 4^N * N:

    R^2 = (N-1) / (N * 4^(N-1))

QED.

---

## Scope

### Valid for
- X_iZ_j, Y_iZ_j, Z_iX_j, Z_iY_j, and any linear combination thereof
- Any graph topology
- Uniform Z-dephasing, any gamma > 0
- All N >= 2

### Not covered
- Mixed Hamiltonians (shadow-balanced + shadow-crossing terms on the
  same bond). These may require a weighted average of the two variances.
- Non-uniform gamma, non-Pauli noise (same limitations as parent proof).

---

## References

| Component | Location |
|-----------|----------|
| Parent proof (balanced) | [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md) |
| Experiment | [experiments/CROSS_TERM_CROSSING.md](../../experiments/CROSS_TERM_CROSSING.md) |
| Verification script | [simulations/cross_term_crossing.py](../../simulations/cross_term_crossing.py) |

---

*The balanced bond contributes zero variance: both Paulis at the same
shadow depth. The crossing bond contributes one: one in the light, one
in the shadow. The difference between the two theorems is this single
unit of asymmetry.*

*Thomas Wicht, Claude (Anthropic), April 14, 2026*
