# Proof of M's SVD Cluster Structure (F78 single-body + F79 two-body)

**Tier:** 1 (fully analytical)
**Date:** April 29, 2026
**Depends on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π operator construction, palindrome equation F1)
- [PROOF_BIT_B_PARITY_SYMMETRY.md](PROOF_BIT_B_PARITY_SYMMETRY.md) (bit_b parity, Π² eigenstructure)
- [`framework/symmetry.py`](../../simulations/framework/symmetry.py) (`pi_action`, `pi_squared_eigenvalue`, `build_pi_full`)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`palindrome_residual`)
- Numerical verification: [`_svd_active_spectator.py`](../../simulations/_svd_active_spectator.py) (single-body, chain/star/complete, N=3-5), [`_svd_two_body_pi_squared_block.py`](../../simulations/_svd_two_body_pi_squared_block.py) (Π²-block verification, N=3-5)
- Pytest locks: `test_F78_single_body_M_additive_decomposition`, `test_F79_two_body_pi_squared_block_decomposition`

**Status:** Proven (Master Lemma + F78 + F79).
**Scope:** Any Hamiltonian H with uniform Z-dephasing γ_l, on any topology (chain, ring, star, complete K_N, arbitrary graph). The Master Lemma is universal in H. F78 specializes to single-body H = Σ_l c_l·P_l. F79 covers all two-body bond-bilinears via three cases (all-even, all-odd, mixed); the structural form of M (block-diagonal vs off-diagonal) is established for all three.

**Does NOT establish:**
- Mixed-Π²-parity 2-body Hamiltonians at the cluster level: the F79 Case 3 derivation gives the structural form (M = M_diag + M_offdiag) but the specific cluster pattern resulting from the combination is not derived in closed form.
- Inhomogeneous γ_l (the Master Lemma carries through, but per-site cancellations need site-by-site γ_l indexing).
- Non-Z dissipators (T1 amplitude damping, X- or Y-dephasing): Master Lemma fails as stated; separate analysis needed.

---

## Theorems

### Theorem F78 (Single-body Additive Decomposition)

For any single-body Hamiltonian H = Σ_l c_l·P_l with P ∈ {X, Y, Z} and uniform Z-dephasing γ_l per site, the palindrome residual M = Π·L·Π⁻¹ + L + 2σ·I (with σ = Σ_l γ_l) decomposes additively:

    M = Σ_l M_l ⊗ I_(others)

where each M_l is a 4×4 matrix on per-site Pauli space. Furthermore:

1. Each M_l is **normal**.
2. For P = X: M_l = 0.
3. For P ∈ {Y, Z}: M_l has eigenvalues +2c_l·γ_l·i and −2c_l·γ_l·i, each with multiplicity 2 in the per-site basis (singular values all equal to 2|c_l|·γ_l). Y and Z give algebraically distinct M_l matrices but spectrally identical.

The full M's eigenvalues are Σ_l ε_l·2c_l·γ_l·i for ε_l ∈ {±1}, each sign-combination with multiplicity 2^N. Singular values are |Σ_l ε_l·2c_l·γ_l|, with the same multiplicities. Cluster sizes are pure sign-combination statistics on the weight vector (c_l).

### Theorem F79 (Two-body Π²-Block Decomposition)

For a two-body bond bilinear Hamiltonian H = Σ_b Σ_t c_t·(P_t ⊗ Q_t)_b (sum over bonds b and bilinear types t), define the Π²-parity of each bilinear term:

    p(P_t, Q_t) = (bit_b(P_t) + bit_b(Q_t)) mod 2 ∈ {0, 1}

where bit_b(I) = bit_b(X) = 0 and bit_b(Y) = bit_b(Z) = 1.

Let V_± denote the ±1 eigenspaces of Π² (i.e., Pauli strings with even/odd total bit_b count).

1. **All p_t = 0 (Π²-even)**: M is block-diagonal in V_+ ⊕ V_-. The off-diagonal blocks M[V_+, V_-] and M[V_-, V_+] vanish exactly.

2. **All p_t = 1 (Π²-odd)**: M is purely off-diagonal between V_+ and V_-. The diagonal blocks M[V_+, V_+] and M[V_-, V_-] vanish exactly.

3. **Mixed parities**: M = M_diag + M_offdiag, where M_diag ≠ 0 receives contributions from p_t = 0 terms and M_offdiag ≠ 0 from p_t = 1 terms.

---

## Master Lemma: Dissipator-conjugation cancellation

### Lemma 1 (γ-independence of M for pure Z-dephasing)

For any Hamiltonian H and uniform Z-dephasing dissipator L_diss = Σ_l γ_l·(Z_l·Z_l − I):

    Π·L_diss·Π⁻¹ + L_diss + 2σ·I = 0

where σ = Σ_l γ_l. Consequently:

    M = Π·L·Π⁻¹ + L + 2σ·I = Π·L_H·Π⁻¹ + L_H

with L_H = −i[H, ·]. **M depends only on H** (independent of γ_l).

**Proof.** L_diss is additive over sites:

    L_diss = Σ_l L_diss,l ⊗ I_(others), where L_diss,l acts on site l only.

In the per-site Pauli basis (I, X, Y, Z):

    L_diss,l = diag(0, −2γ_l, −2γ_l, 0).

Π factors per site as Π = ⊗_l Π_l, and consequently the conjugation T_Π = ⊗_l T_Π,l (this follows from `pi_action` in `framework/symmetry.py`: each site idx maps independently). Therefore:

    Π·L_diss·Π⁻¹ = Σ_l (Π_l · L_diss,l · Π_l⁻¹) ⊗ I_(others).

It suffices to verify the per-site identity:

    Π_l · L_diss,l · Π_l⁻¹ + L_diss,l + 2γ_l·I_l = 0.

In the per-site basis (I, X, Y, Z), the framework's Π_l acts as (`pi_action`):
- I → X (phase 1), X → I (phase 1)
- Y → i·Z, Z → i·Y

Hence Π_l⁻¹ acts as: X → I, I → X, Z → −i·Y, Y → −i·Z.

Computing Π_l · L_diss,l · Π_l⁻¹·|σ⟩ = Π_l · L_diss,l · (Π_l⁻¹|σ⟩):

| σ | Π_l⁻¹·σ | L_diss,l(Π_l⁻¹·σ) | Π_l of result | Diagonal entry |
|---|---------|--------------------|----------------|-----------------|
| I | X | −2γ_l·X | −2γ_l·I | −2γ_l |
| X | I | 0 | 0 | 0 |
| Y | −i·Z | −i·0 = 0 | 0 | 0 |
| Z | −i·Y | −i·(−2γ_l)·Y = 2γ_l·i·Y | 2γ_l·i·(i·Z) = −2γ_l·Z | −2γ_l |

So Π_l · L_diss,l · Π_l⁻¹ has diagonal (−2γ_l, 0, 0, −2γ_l) in the basis (I, X, Y, Z).

Adding L_diss,l = diag(0, −2γ_l, −2γ_l, 0) gives diag(−2γ_l, −2γ_l, −2γ_l, −2γ_l) = −2γ_l·I_l.

Adding 2γ_l·I_l gives 0. **Per-site identity proved.**

Summing over sites:

    Π · L_diss · Π⁻¹ + L_diss = Σ_l (Π_l · L_diss,l · Π_l⁻¹ + L_diss,l) ⊗ I_(others)
                              = Σ_l (−2γ_l·I_l) ⊗ I_(others)
                              = −2σ·I.

Adding 2σ·I gives 0. □

**Corollary.** For pure Z-dephasing, M = Π·L_H·Π⁻¹ + L_H. The structure of M is determined entirely by the Hamiltonian commutator superoperator and its Π-conjugate. ∎

---

## Auxiliary Lemmas

### Lemma 2 (Π preserves bit_b parity)

For any Pauli string σ = σ_0 ⊗ σ_1 ⊗ ... ⊗ σ_{N−1}, the action of Π preserves total bit_b parity: Π(σ) has the same Σ_l bit_b(σ_l) mod 2 as σ.

**Proof.** Per-site action (from `pi_action`):
- I (b=0) → X (b=0): bit_b unchanged.
- X (b=0) → I (b=0): unchanged.
- Y (b=1) → i·Z (b=1): unchanged.
- Z (b=1) → i·Y (b=1): unchanged.

Each site preserves its bit_b. The total Σ_l bit_b(σ_l) is conserved. □

**Corollary.** T_Π preserves the V_± eigenspaces of Π², i.e., T_Π is block-diagonal in the V_+ ⊕ V_- decomposition.

### Lemma 3 (Pauli commutator parity arithmetic)

For two Pauli strings A, B: the commutator [A, B] is either 0 (if A, B commute as operators) or ±2i·AB (if they anticommute). The product AB is a Pauli string (up to phase) with bit_b parity equal to bit_b(A) + bit_b(B) mod 2.

**Proof.** Per-site Pauli multiplication: σ_a · σ_b = δ_{ab}·I + i·ε_{abc}·σ_c (for a, b, c ∈ {x, y, z}, identity case a=b giving I, a≠b giving the third Pauli with phase). Bit_b summing follows since each Pauli has a fixed bit_b: I, X have b=0; Y, Z have b=1. The product of Paulis with bits b_1, b_2 is a Pauli with bit b_1 ⊕ b_2 (mod 2). For full strings, summing per-site gives total parity. □

**Corollary.** Let H = Σ_t c_t·H_t be a sum of Pauli-string Hamiltonians, each with bit_b parity p_t = Σ_l bit_b(H_t,l) mod 2. Then for any Pauli string σ:

    L_H · σ = −i[H, σ] = Σ_t c_t · (some Pauli string σ'_t)

where each σ'_t has bit_b parity equal to (p_t + bit_b(σ)) mod 2.

If all p_t are equal to some common p ∈ {0, 1}, then L_H maps V_+ to V_+ ⊕ V_- with parity shifted by p:
- If p = 0: L_H preserves V_+ and V_- separately.
- If p = 1: L_H maps V_+ → V_- and V_- → V_+.

---

## Proof of F78 (Single-body Additive Decomposition)

For H = Σ_l c_l·P_l with P ∈ {X, Y, Z}:

L_H = −i[Σ_l c_l·P_l, ·] = Σ_l c_l·(−i[P_l, ·]) = Σ_l c_l·K_P^(l).

Each K_P^(l) acts only on site l (since [P_l, σ_l'] = 0 for l ≠ l'). Hence L_H is a sum of single-site superoperators:

    L_H = Σ_l c_l·K_P^(l) ⊗ I_(others).

By the Master Lemma:

    M = T_Π·L_H·T_Π⁻¹ + L_H.

Since T_Π factors per site and L_H is per-site additive:

    T_Π · L_H · T_Π⁻¹ = Σ_l c_l · (T_Π,l · K_P^(l) · T_Π,l⁻¹) ⊗ I_(others).

Therefore:

    M = Σ_l c_l · (T_Π,l · K_P^(l) · T_Π,l⁻¹ + K_P^(l)) ⊗ I_(others) = Σ_l M_l ⊗ I_(others)

where M_l = c_l·(T_Π,l · K_P^(l) · T_Π,l⁻¹ + K_P^(l)).

### M_l for P = X: vanishing

In the per-site (I, X, Y, Z) basis, K_X (= −i·[X, ·]) is:

    K_X = ⎡0  0  0  0 ⎤
          ⎢0  0  0  0 ⎥
          ⎢0  0  0 −2 ⎥
          ⎣0  0  2  0 ⎦

(action: Y → 2Z, Z → −2Y; I, X annihilated).

Computing T_Π,l · K_X · T_Π,l⁻¹ via the |σ⟩ table (I → −i·Z basis…) yields:

    T_Π,l · K_X · T_Π,l⁻¹ = ⎡0  0  0  0 ⎤
                            ⎢0  0  0  0 ⎥
                            ⎢0  0  0  2 ⎥
                            ⎣0  0 −2  0 ⎦

i.e., −K_X. Therefore T_Π,l · K_X · T_Π,l⁻¹ + K_X = 0. **M_l = 0 for P = X.**

### M_l for P = Y

K_Y has action Y commutators: −i[Y, X] = −2Z, −i[Y, Z] = 2X.

    K_Y = ⎡0  0  0  0 ⎤
          ⎢0  0  0  2 ⎥
          ⎢0  0  0  0 ⎥
          ⎣0 −2  0  0 ⎦

Computing T_Π,l · K_Y · T_Π,l⁻¹ via the |σ⟩ table yields:

    T_Π,l · K_Y · T_Π,l⁻¹ = ⎡0  0 −2i  0 ⎤
                            ⎢0  0   0  0 ⎥
                            ⎢−2i  0  0  0 ⎥
                            ⎣0  0   0  0 ⎦

Sum:

    M_l (P=Y) / c_l = ⎡0   0  −2i  0 ⎤
                     ⎢0   0   0   2 ⎥
                     ⎢−2i  0   0   0 ⎥
                     ⎣0  −2   0   0 ⎦

This block-decomposes into two 2×2 blocks: one on (I, Y) with [[0, −2i], [−2i, 0]], one on (X, Z) with [[0, 2], [−2, 0]]. Each 2×2 block has trace 0 and:
- (I, Y) block: det = −(−2i)² = 4 (but trace²=0 so eigenvalues are ±2i).
- (X, Z) block: det = 0 − (2)(−2) = 4, eigenvalues ±2i.

Each block contributes eigenvalues ±2i (one of each). Total: M_l/c_l has eigenvalues +2i, −2i, +2i, −2i (multiplicity 2 each). Multiplied by c_l: ±2c_l·i (mult 2). The zero-trace, equal-eigenvalue-magnitude structure makes each block normal (verifiable: M_l·M_l^† = M_l^†·M_l = 4c_l²·I, all SVs equal 2|c_l|). **M_l for P=Y is normal with eigenvalues ±2c_l·i, mult 2 each.**

### M_l for P = Z

K_Z has action: −i[Z, X] = 2Y, −i[Z, Y] = −2X.

    K_Z = ⎡0  0  0  0 ⎤
          ⎢0  0 −2  0 ⎥
          ⎢0  2  0  0 ⎥
          ⎣0  0  0  0 ⎦

Computation analogous to P=Y yields (by direct application of T_Π,l):

    M_l (P=Z) / c_l = ⎡0   0   0   −2i ⎤
                     ⎢0   0  −2   0   ⎥
                     ⎢0   2   0   0   ⎥
                     ⎣−2i  0   0   0  ⎦

Block-decomposes into (I, Z) and (X, Y) blocks; same trace/det analysis: eigenvalues ±2c_l·i (mult 2 each), normal. **M_l for P=Z is spectrally identical to M_l for P=Y.**

### From M_l to full M's spectrum

Since each M_l is normal and the {M_l ⊗ I_(others)}_l mutually commute (acting on different tensor factors), M = Σ_l M_l ⊗ I_(others) is normal. Eigenvalues of M are sums Σ_l λ_l where λ_l ∈ Spec(M_l). Multiplicities multiply across sites:

    Spec(M) = { Σ_l ε_l · 2c_l·i : ε_l ∈ {+1, −1} }, multiplicity 2^N per sign-combination.

Singular values (= |eigenvalues| for normal M) are |Σ_l ε_l · 2c_l| with the same multiplicities. **F78 proved.** ∎

---

## Proof of F79 (Two-body Π²-Block Decomposition)

By the Master Lemma, M = T_Π·L_H·T_Π⁻¹ + L_H. By Lemma 2, T_Π preserves V_± (block-diagonal). The structure of M follows from how L_H interacts with V_± (Lemma 3).

### Case 1: All terms Π²-even (p_t = 0 for all t)

For each Pauli string σ ∈ V_+ (parity 0), L_H,t · σ has parity (0 + 0) mod 2 = 0, hence stays in V_+. Similarly L_H sends V_- to V_-.

Equivalently: L_H is block-diagonal in V_+ ⊕ V_-:

    L_H = ⎡L_++  0   ⎤
          ⎣ 0   L_-- ⎦   (in V_+ ⊕ V_- decomposition)

Since T_Π is also block-diagonal:

    T_Π·L_H·T_Π⁻¹ = ⎡ T_Π,+ · L_++ · T_Π,+⁻¹     0     ⎤
                    ⎣        0          T_Π,− · L_−− · T_Π,−⁻¹ ⎦

Sum M = T_Π·L_H·T_Π⁻¹ + L_H is also block-diagonal. Off-diagonal blocks vanish exactly:

    M[V_+, V_-] = 0 = M[V_-, V_+].

### Case 2: All terms Π²-odd (p_t = 1 for all t)

For σ ∈ V_+, L_H · σ has parity (1 + 0) mod 2 = 1, hence lands in V_-. Similarly L_H sends V_- to V_+.

L_H is purely off-diagonal:

    L_H = ⎡ 0     L_+- ⎤
          ⎣L_-+   0    ⎦

Conjugating by T_Π (block-diagonal):

    T_Π·L_H·T_Π⁻¹ = ⎡    0           T_Π,+ · L_+- · T_Π,−⁻¹ ⎤
                    ⎣T_Π,− · L_-+ · T_Π,+⁻¹      0          ⎦

Also off-diagonal. Sum M = T_Π·L_H·T_Π⁻¹ + L_H is purely off-diagonal:

    M[V_+, V_+] = 0 = M[V_-, V_-].

Singular values of M arise from the off-diagonal blocks. The block M[V_+, V_-] is some operator V_- → V_+; its singular values appear once from the V_+ side of the decomposition and once from the V_- side of the SVD of M, giving each singular value of M[V_+, V_-] a multiplicity of 2 in M's spectrum.

### Case 3: Mixed parities

If H has terms of both parities:

    L_H = L_H^(even) + L_H^(odd)

where L_H^(even) is block-diagonal and L_H^(odd) is block-anti-diagonal. After T_Π conjugation (which preserves block structure):

    M = (T_Π · L_H^(even) · T_Π⁻¹ + L_H^(even))^(diag) + (T_Π · L_H^(odd) · T_Π⁻¹ + L_H^(odd))^(off-diag)

So M has both diagonal and off-diagonal contributions. The specific cluster pattern requires computing each piece. **F79 proved (cases 1 and 2 explicitly; case 3 by linearity).** ∎

---

## Connection to existing framework formulas

- **F1 (Palindrome equation)**: For "truly" Hamiltonians, M = 0 exactly. Both single-body X (F78 case M_l = 0) and bond-bilinear Heisenberg XX+YY+ZZ (F79: all terms Π²-even, but additionally satisfying the truly condition) recover F1. F78 and F79 generalize to non-truly cases by giving the explicit cluster structure of M.

- **F49 (Cross-term formula, Frobenius)**: The Frobenius norm ‖M‖²_F is per-bond additive (F49). F79 refines this by giving the full SVD spectrum, not just the norm.

- **F61 (n_XY parity selection rule)**: The bit_a parity (n_XY) selection rule selects which Pauli strings can mix under L. F79's bit_b parity is a complementary Π² selection rule.

- **F63 ([L, Π²] = 0 for w_YZ-even)**: F63 proves that for Π²-even Hamiltonians, [L, Π²] = 0. F79 makes this concrete at the M-block level: even Hamiltonians give block-diagonal M, odd Hamiltonians give purely off-diagonal M.

---

## Lemma 4 (M is anti-Hermitian for any Hermitian H, Z-dephasing)

For any Hermitian Hamiltonian H and uniform Z-dephasing, the palindrome residual M = Π·L·Π⁻¹ + L + 2σ·I is **anti-Hermitian**: M^† = −M.

**Proof.** By the Master Lemma, M = T_Π·L_H·T_Π⁻¹ + L_H. We show each summand is anti-Hermitian.

(a) L_H = −i[H, ·] is anti-Hermitian on operator space. For the Frobenius inner product ⟨A, B⟩ = tr(A^† B):

    ⟨L_H(A), B⟩ = tr((−i[H, A])^† B) = i·tr((HA − AH)^† B)
                = i·tr((A^†H − HA^†) B)  [since H = H^†]
                = i·(tr(A^†HB) − tr(HA^†B))

    −⟨A, L_H(B)⟩ = −tr(A^†·(−i[H, B])) = i·tr(A^†(HB − BH))
                 = i·(tr(A^†HB) − tr(A^†BH))

By cyclic invariance, tr(HA^†B) = tr(A^†BH). Hence ⟨L_H(A), B⟩ = −⟨A, L_H(B)⟩, so L_H^† = −L_H. ✓

(b) T_Π is unitary on operator space (it permutes Pauli-string basis vectors with phases of unit modulus). Hence T_Π·L_H·T_Π⁻¹ is anti-Hermitian: (T_Π L_H T_Π⁻¹)^† = T_Π·L_H^†·T_Π⁻¹ = T_Π·(−L_H)·T_Π⁻¹ = −T_Π L_H T_Π⁻¹. ✓

(c) Sum of anti-Hermitian matrices is anti-Hermitian: M = (T_Π·L_H·T_Π⁻¹) + L_H satisfies M^† = −M. □

**Corollary 4.1.** Eigenvalues of M are purely imaginary (Re(λ) = 0 for any eigenvalue λ of an anti-Hermitian operator). Furthermore, M is normal (since anti-Hermitian matrices satisfy M·M^† = M^†·M = −M²). Singular values of M are equal to |λ| = |Im(λ)|.

**Corollary 4.2 (spectrum-equality ⇒ unitary equivalence).** If two Hamiltonians H_1, H_2 yield M_1, M_2 with identical spectra (as multi-sets), then M_1 and M_2 are unitarily equivalent: there exists a unitary U on operator space such that U·M_1·U^† = M_2.

This follows from the spectral theorem for normal operators: each anti-Hermitian M_k diagonalizes as M_k = V_k·D·V_k^†, where D is the common diagonal matrix and V_k are unitaries. Then U = V_2·V_1^† satisfies U·M_1·U^† = M_2. ∎

---

## Π²-odd Universality (Empirical observation, structurally bounded)

**Observation.** Within the Π²-odd pure 2-body bilinear class (i.e., terms (P, Q) with P, Q ∈ {X, Y, Z}, P · Q ∈ {Y, Z} parity-odd), the specific Pauli-letter choice is **M-spectrum-irrelevant**. Verified at chain N=3, 4, 5 and star N=4: all four Π²-odd pure 2-body single bilinears (X,Y), (X,Z), (Y,X), (Z,X) give bit-identical M-eigenvalues (and consequently bit-identical SVs).

By Corollary 4.2, this means M_XY, M_XZ, M_YX, M_ZX are pairwise **unitarily equivalent**.

**Frobenius accounting (necessary condition).** F49 gives ‖M‖²_F as a function of bond pattern and per-bond coefficients squared. For any Π²-odd 2-body single-bond bilinear (P, Q) with |c|² = 1, ‖M‖² is determined entirely by N and topology — not by the specific (P, Q) letters within the parity class. So Frobenius norm equality is automatic for these pairs. Spectrum equality is strictly stronger.

**Bound: where the universality breaks.** The Π²-even non-truly class (only YZ, ZY at the pure 2-body level) is universal in the same way: M_YZ ≡ M_ZY at the spectrum level. But mixing across parities (e.g., adding YZ + XY) produces a non-pure Hamiltonian whose M has both diagonal and off-diagonal contributions and a richer spectrum. The universality is strictly within-parity-class.

**Open structural question.** The full unitary U realizing M_XY ≡ M_XZ does not factor as a per-site Pauli-letter swap or per-site state-space rotation, because:
- Per-site label-swap (Y ↔ Z basis permutation at one site) commutes with the framework's Π, but does not commute with L_diss (Z-dephasing distinguishes Y and Z via bit_b at each site). Although the Master Lemma bypasses L_diss for M, the per-site swap also does not directly map L_H_XY to L_H_XZ at the chain level (it changes one bond but not the others' Y → Z).
- Per-site state-space rotation V = e^{iπ X / 4} maps Y → Z and conjugates X⊗Y → X⊗Z at the Hamiltonian level, but does NOT commute with the framework's Π (Π has phase i for Y ↔ Z; rotation has sign −1 for Z → −Y).

So the unitary equivalence exists abstractly (Corollary 4.2) but does not reduce to a single per-site or natural global symmetry. A proof of spectrum equality would need a different route — most likely identifying that M's characteristic polynomial depends only on the bilinear-class invariants (Π²-parity, Frobenius norm, bond graph), not on Pauli-letter details. We leave this as an open structural question.

---

---

## Numerical Verification

| Theorem | Verification Path | Test Case |
|---------|-------------------|-----------|
| Master Lemma | `palindrome_residual` numerical check | Any H with Z-dephasing, ‖M(γ=1) − M(γ=2)‖ = 0 (after subtracting proper σ shifts) |
| F78 (single-body, additive) | `_svd_active_spectator.py` | chain/star/complete N=3-5, P ∈ {X,Y,Z}, all match additive prediction exactly |
| F78 normality of M_l | `_svd_single_body_extension.py` | M_l direct construction, normality check returns True |
| F79 Π²-even block-diag | `_svd_two_body_pi_squared_block.py` | YZ, YZ+ZY, Heisenberg: ‖M[V_+,V_-]‖ = 0 exactly |
| F79 Π²-odd off-diag | `_svd_two_body_pi_squared_block.py` | XY, XZ, XX+XY: ‖M[V_+,V_+]‖ = ‖M[V_-,V_-]‖ = 0 exactly |
| F79 Π²-odd universality | `_svd_two_body_structure.py` | XY ≡ XZ ≡ XX+XY ≡ XX+XZ at N=5 chain, identical SV-cluster spectra |
| Pytest locks | `simulations/framework/tests/test_smoke.py` | `test_F78_single_body_M_additive_decomposition` (chain N=4 IY+YI clusters), `test_F79_two_body_pi_squared_block_decomposition` (block structure + universality) |
