# F112 Non-Hermitian Extension: Universal-N Closure (Welle 11), Basis-Enumeration Anchor at N=2, 3, 4, 5

**Status:** Tier1Derived for all N (Welle 11, 2026-05-27) via the two-lemma structural proof in [`PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`](../docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md). The basis-enumeration result at N ≤ 5 below is preserved as the **empirical anchor** that motivated the search for the structural proof; it remains the historical numerical validation (559,912 pair F-values bit-exact 0).
**Date:** 2026-05-26 (initial enumeration N≤4); 2026-05-26 Welle 10b (extended to N=5); 2026-05-27 Welle 11 (universal-N structural proof).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Scripts:** Structural proof verifier [`simulations/_f112_universal_n_proof_verify.py`](../simulations/_f112_universal_n_proof_verify.py) (Welle 11, N=1, 2, 3) + Python enumeration [`simulations/_f112_open_identity_basis_enum.py`](../simulations/_f112_open_identity_basis_enum.py) (Welle 10a, N=2..5) + C# enumeration [`compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs`](../compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs) (Welle 10b, N=2..5 via SLOW_F112-tagged test, parallelized via Parallel.For/ForEach)
**Connects:** [PROOF_F112_NONHERMITIAN_UNIVERSAL_N](../docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md) (the structural proof, this writeup's parent), [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) (Hermitian-H parent theorem), [F112 ANALYTICAL_FORMULAS entry](../docs/ANALYTICAL_FORMULAS.md), [LindbladBitBPiBalance](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs)

## The open identity

F112 typed Tier1Derived covers Hermitian H + bit_b-homogeneous c. The non-Hermitian extension allows H = H_re + i·H_im (with H_re, H_im both Hermitian); the F112 proof's Step 5 was previously rigorous only for Hermitian H. The non-Hermitian case had been observed bit-exact empirically across 20+ random configurations (probes 1-14), but the structural proof was deferred to:

> **Open identity:** F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im.

If this identity holds, then F112 holds for any non-Hermitian H.

## Reduction (2026-05-26)

Two algebraic observations:

**Observation 1 (bilinearity).** L_H = -i[H, ·] is linear in H. The Π-conjugation eigenspace projection is also linear. Hence L_{H,-i} is linear in H. The Frobenius inner product ⟨A, B⟩ is sesquilinear (anti-linear in A, linear in B). For real-linear inputs (Hermitian operators expanded in real coefficients of the Pauli basis), F(H_re, H_im) is real-bilinear in (H_re, H_im).

**Observation 2 (antisymmetry).** The Hilbert-Schmidt inner product satisfies ⟨X, Y⟩^* = ⟨Y, X⟩, so Im⟨X, Y⟩ = -Im⟨Y, X⟩. Under H_re ↔ H_im exchange, L_{H_re,-i} ↔ L_{H_im,-i}, so F(H_re, H_im) = -F(H_im, H_re). F is antisymmetric in the exchange.

**Consequence.** F is determined by its values on a basis of pairs of Hermitian operators. The Hermitian operator space at chain length N is spanned by the 4^N Pauli strings with real coefficients. If F(σ_α, σ_β) = 0 for every pair of Pauli strings (σ_α, σ_β), then F ≡ 0 on the entire Hermitian operator space.

## Numerical enumeration result

The script `_f112_open_identity_basis_enum.py` enumerates F on all upper-triangular Pauli-string pairs at N=2, 3, 4. Welle 10b adds N=5 via two independent pipelines:

| N | distinct upper-triangular pairs | max \|Im F\| | pipeline | wall time |
|---|---|---|---|---|
| 2 | 136 | 0.0000e+00 | Python + C# | < 1 sec each |
| 3 | 2,080 | 0.0000e+00 | Python + C# | < 1 sec each |
| 4 | 32,896 | 0.0000e+00 | Python + C# | ~25 sec each |
| 5 | 524,800 | 0.0000e+00 (Python); < 1e-10 (C#) | Python + C# | Python 90.7 min; C# 2 h 45 m (outerDop=6) → ~85 min (outerDop=24 forecast after Welle 10c) |

All 559,912 pair F-values across N=2, 3, 4, 5 give Im below the 1e-10 threshold (N=2..4 bit-exact 0.0e+00; N=5 verified by Python with bit-exact 0.0e+00 across 524,800 pairs, and re-verified independently in C# via SLOW_F112-tagged xUnit test with MaxImaginary < 1e-10). The bit_b-parity breakdown also gives bit-exact 0 in each of the four (α_bit_b, β_bit_b) cells separately at every N tested, confirming this is not a sum-cancellation artifact restricted to specific bit_b sectors.

The C# pipeline runs on the canonical typed-knowledge layer; both pipelines use the same algorithm (build L_α,-i per Pauli string, then Frobenius inner product per upper-triangular pair) and share the Π conjugation construction up to language-level differences.

## Conclusion

**F112 non-Hermitian extension is Tier1Derived for all N** (Welle 11, 2026-05-27) via the two-lemma structural proof. The basis enumeration at N=2, 3, 4, 5 (559,912 pairs total, all Im bit-exact 0 or < 1e-10) was the empirical anchor that motivated the proof; the structural proof now lifts the closure from N ≤ 5 to arbitrary N.

The closure of the universal-N path makes the basis-enumeration argument obsolete as a Tier-promotion mechanism but preserves its standing as historical empirical validation: 559,912 pairs at N ≤ 5 already verified the per-pair identity bit-exact, providing strong numerical evidence that motivated the search for and gave confidence in the structural proof.

## Welle-11 closure: two-lemma structural proof

The universal-N closure rests on two lemmas (see [`PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`](../docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md) for the full proof):

**Lemma N-A (Diagonal-Norm).** For any BitB-odd Pauli string σ at chain length N, ‖L_{σ,−i}‖² = 4^N exactly. Proof composes (i) the support-count identity ‖L_σ‖² = 2 · 4^N (counting Pauli strings that anticommute with σ); (ii) the F38 / F63 Π²-conjugation eigenvalue (−1)^{BitBParity(σ)} = −1 for BitB-odd σ, putting L_σ in the Π-conjugation {+i, −i} eigenspaces; (iii) the cross-term vanishing ⟨L_σ, Π L_σ Π⁻¹⟩ = 0 via Pauli-basis matrix-support disjointness (M(L_σ) and Π M(L_σ) Π⁻¹ have complementary support on the same shifted-diagonal positions); (iv) the standard "halving" identity ‖L_{σ,±i}‖² = (1/2) ‖L_σ‖² = 4^N.

**Lemma N-B (Off-Diagonal-Orthogonality).** For σ_α ≠ σ_β both BitB-odd at chain length N, ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0 exactly. Proof reduces to ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩ = 0 for all m ∈ {0, 1, 2, 3}, each established by matrix-support disjointness: M(L_α) is supported on {(σ_α ⊕ α', α')}, Π^m M(L_β) Π^{−m} is supported on {(σ_β ⊕ α', α')} (possibly with shifted "condition on α'"); overlap requires σ_α = σ_β.

Both Welle 11 lemmas (N-A and N-B) reduce to per-position checks on the 4^N × 4^N matrix of L_σ that are uniform in N; the proof is N-independent.

**Verifier:** [`simulations/_f112_universal_n_proof_verify.py`](../simulations/_f112_universal_n_proof_verify.py) confirms each step within 1e-12 numpy double-precision tolerance at N = 1, 2, 3 (42 BitB-odd strings, 1050 off-diagonal pairs, 4368 all-pair F-values, all max deviations < 1e-12, i.e. machine zero to numpy double precision). The Welle 10a Python enumeration above is genuinely bit-exact at N ≤ 4 (rational matrix entries); the Welle 11 verifier is numerical.

## Implications

- **F112 non-Hermitian extension is Tier1Derived for all N.** The algebraic argument bilinearity + Pauli-basis spanning reduces F = 0 to the per-pair identity F(σ_α, σ_β) = 0; the per-pair identity holds structurally via Lemmas N-A and N-B.
- **The polarity_coordinates_from_L diagnostic** is a structural witness for L not in the Lindblad form `−i[H, ·] + Σ γ_k np.kron(c_k, c_k^*)` with bit_b-homogeneous c, universally in N and for any H (Hermitian or non-Hermitian).
- **The bit_b Z₂-axis carries three Tier1Derived universal-N theorems**: F108 Parts 1/2/3 (palindrome closure of bit_b = 0 bilinears), F112 (Hermitian and non-Hermitian H, now universal N), F113 (T1 break-magnitude closed form). The bit_b axis description is structurally complete on the BitB-axis side.

## Reproduction

**Python (N=2..5):**
```
python -X utf8 simulations/_f112_open_identity_basis_enum.py        # N=2 + N=3 by default
PYTHONIOENCODING=utf-8 python simulations/_f112_n5_run.py            # N=5 wrapper, ~16 GB, ~90 min
```

**C# (canonical typed-knowledge layer):**
```
# Fast tests (N=2, 3, 4):
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release \
    --filter "FullyQualifiedName~F112NonHermitianBasisEnumerationTests&Category!=SLOW_F112"

# SLOW_F112 N=5 test (~85 min on 24-core, ~16 GB):
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release \
    --filter "Category=SLOW_F112"
```

Runs N=2 + N=3 in under 5 seconds (both pipelines); N=4 in ~25 sec (Python) or ~15 sec (C# with outerDop=24). N=5: Python 90.7 min, C# 2 h 45 m at outerDop=6 (Welle 10b initial) or ~85 min at outerDop=24 (Welle 10c). N=6 (8.4M pairs) requires sparse Pauli representation.

## Related

- [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md): the parent theorem; this writeup fills the "Step 5 extension to non-Hermitian H" open item listed there.
- [F112_HARDWARE_LENS_KINGSTON.md](F112_HARDWARE_LENS_KINGSTON.md): hardware-side application of F112 that surfaced the first structural counterexample to the broader empirical envelope (Z-drive + σ⁻ T1). The non-Hermitian extension proven here is orthogonal to that counterexample (the counterexample concerns bit_b-mixed c with Hermitian H, not non-Hermitian H with bit_b-homogeneous c).
- [LindbladBitBPiBalance.cs](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs): the typed Tier1Derived C# Claim; this writeup informs the NonHermitianExtension inspectable.
