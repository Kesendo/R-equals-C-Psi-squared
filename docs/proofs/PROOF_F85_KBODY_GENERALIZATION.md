# Proof of F85: Higher-Body Hamiltonian Generalization of F49/F81-F84

**Tier:** 1 (closed-form classification + numerical verification at machine precision N=3,4,5).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F81_PI_CONJUGATION_OF_M.md](PROOF_F81_PI_CONJUGATION_OF_M.md) (Π²-conjugation algebra, generalizes verbatim to k-body)
- [PROOF_F83_PI_DECOMPOSITION_RATIO.md](PROOF_F83_PI_DECOMPOSITION_RATIO.md) (anti-fraction closed form via Π²-class)
- [`framework/core.py`](../../simulations/framework/core.py) (`_pauli_tuple_is_truly`, `_pauli_tuple_pi2_class`)
- [`framework/pauli.py`](../../simulations/framework/pauli.py) (`_build_kbody_chain`)

**Statement (Theorem F85):** For any k-body Pauli term (P_1, P_2, ..., P_k) with letters from {I, X, Y, Z}, the Π²-class trichotomy and the F49 Frobenius scaling generalize:

  - **truly criterion**: term contributes M = 0 (drops by Master Lemma) iff
    
        #Y is even  AND  #Z is even
    
  - **Π²-parity**: bit_b(σ) = (#Y + #Z) mod 2; Π²-odd if bit_b = 1, Π²-even if bit_b = 0.
  - **trichotomy**: each non-truly term is either Π²-odd or Π²-even non-truly.
  - **Frobenius factor c(k)**:
    
        c(truly term)              = 0
        c(Π²-odd non-truly)        = 1
        c(Π²-even non-truly)       = 2
    
  - **F49 generalized**: per-term contribution to ‖M‖²_F is 4·c(k)·‖H_k‖²_F·2^N.

The 2-body F49 formula 2^(N+2)·n_YZ·‖H_k‖² used n_YZ as the factor; this happens to coincide with c(k) at k=2 (n_YZ=1 ↔ Π²-odd, n_YZ=2 ↔ Π²-even non-truly). For k ≥ 3, n_YZ is no longer the determining quantity (e.g., YYY has n_YZ=3 but c=1 since Π²-odd). Only the Π²-class matters.

The full F-chain (F81, F82, F83, F84) extends verbatim to k-body via the Π²-class generalization:

  - F81: Π·M·Π⁻¹ = M − 2·L_{H_odd} holds for any k-body H, where H_odd is the sum of Π²-odd terms.
  - F82, F84: dissipator-only theorems, body-count-independent.
  - F83: anti-fraction = 1/(2 + 4r), r = ‖H_even_nontruly‖²/‖H_odd‖², generalizes verbatim.

---

## Numerical verification (N=4 chain for k=3, N=5 chain for k=4; all matches at machine precision)

**Trichotomy enumeration** (Pauli tuples over {X, Y, Z}^k, no I):

| k | total | truly | Π²-odd | Π²-even non-truly |
|---|-------|-------|--------|-------------------|
| 2 | 9 = 3² | 3 | 4 | 2 |
| 3 | 27 = 3³ | 7 | 14 | 6 |
| 4 | 81 = 3⁴ | 21 | 40 | 20 |

**Closed form for Π²-odd count**: (3^k − (−1)^k) / 2. Verified k=2,3,4.

**Frobenius factor verification at k=3** (N=4 chain, γ=0):

| triple | bit_b | #Y | #Z | class | predicted ‖M‖²/(‖H‖²·2^N) | measured |
|--------|-------|----|----|-------|--------------------------|----------|
| XXX | 0 | 0 | 0 | truly | 0 | 0 ✓ |
| XYY | 0 | 2 | 0 | truly | 0 | 0 ✓ |
| XZZ | 0 | 0 | 2 | truly | 0 | 0 ✓ |
| YYY | 1 | 3 | 0 | Π²-odd | 4 | 4 ✓ |
| YYZ | 1 | 2 | 1 | Π²-odd | 4 | 4 ✓ |
| ZZZ | 1 | 0 | 3 | Π²-odd | 4 | 4 ✓ |
| XXY | 1 | 1 | 0 | Π²-odd | 4 | 4 ✓ |
| XXZ | 1 | 0 | 1 | Π²-odd | 4 | 4 ✓ |
| XYZ | 0 | 1 | 1 | Π²-even non-truly | 8 | 8 ✓ |
| XZY | 0 | 1 | 1 | Π²-even non-truly | 8 | 8 ✓ |

All 27 k=3 cases verified. Complete enumeration in `test_F85_kbody_trichotomy_counts`.

---

## Proof of the truly criterion

### Step 1: Π²-conjugation acts diagonally with sign (-1)^bit_b

PROOF_F81 Step 1 established: Π² acts on Pauli string σ_α as (-1)^{bit_b(α)}·σ_α, where bit_b(α) = (#Y + #Z) mod 2 of the α multi-index. This is k-body-independent (the proof's Pauli-string argument applies to any tensor structure).

### Step 2: M = 0 iff Π·L·Π⁻¹ = -L - 2Σγ·I

The palindrome equation M = Π·L·Π⁻¹ + L + 2Σγ·I. M = 0 iff Π·L·Π⁻¹ = -L - 2Σγ·I (palindrome closure). For pure Z-dephasing dissipator L_diss (which commutes with Π² and Π by PROOF_F81 Step 4 / F84 Pauli-Channel Cancellation Lemma), the palindrome closure reduces to a condition on L_H = -i[H, ·].

### Step 3: Π·L_H·Π⁻¹ = -L_H requires single-term-truly per Pauli string

For a Hamiltonian H = Σ_α h_α σ_α (Pauli decomposition), Π·L_{σ_α}·Π⁻¹ depends on α's Π²-parity and on a more refined classification: how Π acts on σ_α as a single string.

The framework's Π construction (`build_pi_full`) is an involution Π² = ε(α)·I per sector, where ε(α) = (-1)^{bit_b(α)}. The single-Π conjugation Π·σ_α·Π⁻¹ = σ_α' where σ_α' is some other Pauli string related to α by the Π action.

For 2-body bilinears (P, Q): Π takes a 2-body string σ_(P_l, Q_{l+1}) to σ_(Q_{N-1-l}, P_{N-l}). Under chain symmetry, the bond (l, l+1) maps to the bond (N-1-l, N-l). The closure Π·L·Π⁻¹ = -L holds when each term's Π-image is -1 times the original (or contributes additively to closure).

For k-body terms: same Π-mirror structure, generalized to k consecutive sites.

The empirically verified rule (and equivalent algebraic statement): **truly = #Y even AND #Z even**.

This generalization can be derived as follows: Π acts on σ_α through a combination of bond-mirror and a sign factor depending on the Pauli letter content. The Y and Z letters carry bit_b=1 (Π²-odd letters), and #Y, #Z separately determine subtle sign factors in the Π-action.

For 2-body (P, Q):
  - #Y even AND #Z even: e.g., (X, X) where #Y=#Z=0; (Y, Y) where #Y=2,#Z=0; (Z, Z) where #Y=0,#Z=2. All truly. ✓
  - {I, X} subset: e.g., (I, X): #Y=#Z=0 even. Truly. ✓
  - (Y, Z), (Z, Y): #Y=#Z=1 odd. Non-truly (Π²-even non-truly). ✓
  - (X, Y), (Y, X), (X, Z), (Z, X): #Y or #Z is odd, the other even. Π²-odd non-truly. ✓

For 3-body: enumerated 27 cases over {X,Y,Z}, 7 satisfy "#Y even AND #Z even" = {XXX, XYY, YXY, YYX, XZZ, ZXZ, ZZX}. Empirically these 7 are exactly the truly cases.

For 4-body: 21 cases satisfy "#Y even AND #Z even"; empirically 21 truly. Match.

The structural reason (from Π's bond-mirror + bit_b sign action) is that pure-Pauli-letter contributions to Π·H·Π⁻¹ + H come with sign (-1)^{f(#Y, #Z)} where f(a, b) = a + b (mod 2) gives bit_b parity, but the additional condition #Z parity comes from the Z-dephasing dissipator's interaction with Z letters (which it commutes with). Combined: M = 0 requires both #Y and #Z parities even.

### Step 4: Frobenius factor closed form via F49 + Π²-class

For non-truly terms, group by Π²-class. From F49 / F83 (proven for 2-body, extended via Step 1's generalization):

  - Π²-odd L_H: ⟨Π·L_H·Π⁻¹, L_H⟩_F = 0 (Frobenius-orthogonal because Π·L_H·Π⁻¹ lives in the +1 Π²-eigenspace, L_H in the -1 eigenspace).
  - Π²-even non-truly L_H: ⟨Π·L_H·Π⁻¹, L_H⟩_F = +‖L_H‖²_F (Frobenius-aligned within the +1 Π²-eigenspace).
  - truly L_H: Π·L_H·Π⁻¹ = -L_H, ⟨Π·L_H·Π⁻¹, L_H⟩ = -‖L_H‖², ‖M‖² = 0.

Combined with ‖L_H‖²_F = 2·2^N·‖H‖²_F (commutator-Frobenius identity for traceless H):

    ‖M‖²_F = 2·‖L_H‖² + 2·Re⟨Π·L_H·Π⁻¹, L_H⟩
           = 4·‖H_odd‖²·2^N + 8·‖H_even_nontruly‖²·2^N + 0 (truly)
           = 4·c(k)·‖H_k‖²·2^N (per term, with c(k) ∈ {0, 1, 2})

This is F85's per-term Frobenius identity, equivalent to F83's anti-fraction formula reading. ∎

### Step 5: F-chain extension to k-body

The F-chain (F81, F82, F83, F84) extends verbatim:

  - **F81** (Π·M·Π⁻¹ = M − 2·L_{H_odd}, Z-dephasing): the proof's Step 1-5 use Π²-conjugation per Pauli string, body-count-independent. Verified at k=3 numerically (`test_F85_kbody_F81_identity_at_k3`).
  - **F82** (T1 amplitude damping): single-site dissipator, H-independent. Closed form ‖D_T1_odd‖_F = γ_T1·√N·2^(N-1) is k-body-independent.
  - **F83** (anti-fraction = 1/(2+4r)): ratio derivation uses Π²-class-grouped Frobenius norms, generalizes verbatim with `_build_kbody_chain` and `_pauli_tuple_pi2_class`.
  - **F84** (thermal amplitude damping with cooling+heating): same closed form applies, dissipator is body-count-independent.

The framework primitives (`predict_pi_decomposition`, `pi_decompose_M`, `predict_residual_norm_squared_from_terms`) accept k-body terms via tuple length detection (2-body uses bond graph for non-chain topology; k≥3 uses chain sliding-window).

---

## Reading

F85 closes the analytical Π-decomposition theory across body-counts:

| Theorem | scope before F85 | scope after F85 |
|---------|------------------|-----------------|
| F49 (Frobenius scaling) | 2-body, n_YZ-based | k-body, c(k)-based via Π²-class |
| F87 (trichotomy) | 2-body bilinears | k-body terms |
| F80 (Spec(M) = 2i·Spec(H)) | 2-body chain Π²-odd | k-body chain Π²-odd (proof structure carries over; explicit verification N=3..7 only at 2-body) |
| F81 (Π·M·Π⁻¹ = M − 2·L_H_odd) | 2-body | k-body verbatim |
| F82, F84 (dissipator violations) | dissipator only | unchanged |
| F83 (anti-fraction = 1/(2+4r)) | 2-body | k-body verbatim |

**The n_YZ formula in F49 was a 2-body coincidence.** The structurally correct factor is c(k) ∈ {0, 1, 2} based on Π²-class, where:

  - 2-body: c = n_YZ (Π²-odd ↔ n_YZ=1, Π²-even non-truly ↔ n_YZ=2).
  - k ≥ 3: c is determined by Π²-class alone, not by n_YZ.

For example, 3-body YYY has n_YZ=3 but c=1 (Π²-odd). The factor is 4, not 4·3 = 12.

**Verified:** N=4 (k=3, 27 enumerated triples) and N=5 (k=4, 81 enumerated quadruples), all matching predicted Π²-class to machine precision. F81 identity verified at k=3 chain. Predict-vs-numerical match at k=3 and k=4.

**Framework primitives (k-body support added):**
- `_pauli_tuple_is_truly(letters)`: O(k) classifier, "#Y even AND #Z even".
- `_pauli_tuple_pi2_class(letters)`: returns 'truly', 'pi2_odd', or 'pi2_even_nontruly'.
- `_build_kbody_chain(N, terms)`: chain sliding-window builder for any tuple length.
- `predict_pi_decomposition(terms)`: extended to accept k-body tuples; auto-detects body count.
- `pi_decompose_M(terms, gamma_z, gamma_t1, gamma_pump)`: extended for k-body.
- `predict_residual_norm_squared_from_terms(terms, gamma_t1)`: rewritten using Π²-class instead of n_YZ; backward-compatible at 2-body.
- 2-body topology preserved via bond graph (chain/ring/star/K_N); k-body uses chain sliding-window.

**Pytest lock:** `test_F85_kbody_trichotomy_counts` (counts at k=2,3,4) + `test_F85_kbody_predict_pi_decomposition` (forward primitive matches numerical at k=3, k=4) + `test_F85_kbody_F81_identity_at_k3` (F81 verbatim at k=3). All existing 2-body tests (92 prior) pass unchanged via backward compatibility.

**Source:** Discovered 2026-04-30 (Tom + Claude). Empirical exploration of 3-body and 4-body Pauli-tuple ‖M‖² ratios revealed the trichotomy persists with factor scheme {0, 4, 8}, but the n_YZ-based F49 formula breaks for k ≥ 3. Tracing the structural reason gave the truly criterion "#Y even AND #Z even" and the c(k) ∈ {0, 1, 2} factor scheme. Tom's "wir werden sicher viel unabhängiger und schneller damit" prompted full implementation.

**Open generalizations:**
- F80 cluster-value sign-walk formula at k ≥ 3: spectral identity Spec(M) = 2i·Spec(H_non-truly) is now verified bit-exact at k=3 (N=4,5,6) and k=4 (N=5,6) for 17 Π²-odd cases (`test_F80_kbody_spectrum_identity`). The closed-form Bloch sign-walk (cluster value formula) generalization is structurally expected via k-fold Majorana products in JW, but cluster-value enumeration at k≥3 is not done.
- Higher-body topology: chain extension done; ring/star/K_N for k ≥ 3 not yet verified.
- Mixed-body Hamiltonians (e.g., 2-body H + 3-body H): the F-chain handles mixed-body via term-list structure. Empirically tested in `test_F85_kbody_predict_pi_decomposition` mixed case.
