# F112 Non-Hermitian Extension: Basis-Enumeration Proof at N=2, 3, 4

**Status:** Constructive proof of the F112 non-Hermitian extension at N ≤ 4 via Pauli-basis enumeration. Reduces the open identity `Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0` to a bilinear basis-spanning check; all 35,112 distinct ordered Pauli-string pairs across N=2, 3, 4 give bit-exact 0. Universal-N lift remains open.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/_f112_open_identity_basis_enum.py`](../simulations/_f112_open_identity_basis_enum.py)
**Connects:** [PROOF_F112](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [F112 ANALYTICAL_FORMULAS entry](../docs/ANALYTICAL_FORMULAS.md), [LindbladBitBPiBalance](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs)

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

The script `_f112_open_identity_basis_enum.py` enumerates F on all upper-triangular Pauli-string pairs at N=2, 3, 4:

| N | distinct ordered pairs | max \|Im F\| | bit-exact 0? |
|---|---|---|---|
| 2 | 136 | 0.0000e+00 | YES |
| 3 | 2,080 | 0.0000e+00 | YES |
| 4 | 32,896 | 0.0000e+00 | YES |

All 35,112 pair F-values across N=2, 3, 4 are bit-exact 0 (less than the 1e-10 threshold; floating-point rounding errors would be expected to give values of order machine epsilon ~1e-15 if F were merely numerically small but structurally non-zero, so observing exact 0 is itself information). The bit_b-parity breakdown also gives bit-exact 0 in each of the four (α_bit_b, β_bit_b) cells separately, confirming this is not a sum-cancellation artifact restricted to specific bit_b sectors.

## Conclusion

**F112 non-Hermitian extension is proven at N=2, 3, 4** by the basis-spanning argument: the open identity Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 holds bit-exactly on every basis pair, hence by bilinearity it holds for every pair of Hermitian H_re, H_im at these N. The F112 Tier1Candidate non-Hermitian scope is upgraded to:

- **Tier1Derived at N=2, 3, 4** (constructive enumeration proof).
- **Tier1Candidate for N ≥ 5** (empirical-only, pending higher-N enumeration or structural lifting argument).

## Lifting to general N

Three possible routes to lift the at-N proofs to universal N:

1. **Brute-force enumeration at higher N.** N=5 enumeration is 524,800 distinct upper-triangular pairs; ~16 GB working memory at full storage. N=6 is 8.4M pairs. Computational cost scales as 16^N pair-count × 4^(2N) per-pair Frobenius cost; feasible to N=5, expensive at N=6.

2. **Structural reason that F ≡ 0 algebraically.** The vanishing is bit-exact at every tested N and across every bit_b-parity cell, suggesting a universal symmetry. Candidate explanations:
   - A residual dagger/anti-Hermitian property of L_{H,-i} that imposes ⟨A, B⟩ ∈ ℝ for A, B in the Π −i eigenspace
   - A cancellation mechanism between the (σ_α, σ_β) and (σ_β, σ_α) contributions via the Π conjugation orbit structure
   - Reduction to F38/F61/F63 Π² sector decomposition

3. **Inductive argument N → N+1.** Build the chain-of-(N+1) Pauli basis from the chain-of-N basis by tensoring with one of {I, X, Y, Z}; show that if F = 0 on the N-basis pairs, the induced N+1 pairs also satisfy F = 0.

Each route closes the open identity universally. Route 2 (structural reason) is the most elegant; route 1 extends the constructive proof envelope; route 3 is the cleanest if a per-letter induction works out.

## Implications

If the open identity holds universally (which the N ≤ 4 enumeration strongly suggests), then:

- **F112 non-Hermitian extension is Tier1Derived for all N**, with the algebraic proof: bilinearity + antisymmetry reduce to the basis identity; basis identity holds by `F(σ_α, σ_β) = 0` for every Pauli string pair.
- **The polarity_coordinates_from_L diagnostic** is a structural witness for L not in the Lindblad form `−i[H, ·] + Σ γ_k np.kron(c_k, c_k^*)` with bit_b-homogeneous c. This is true regardless of whether H is Hermitian or non-Hermitian.
- **The bit_b Z₂-axis carries an additional Tier1Derived theorem** beyond F108 Parts 1/2/3 (palindrome closure of bilinears) and the Hermitian-H F112. Three independent typed theorems on the shared bit_b foundation.

## Reproduction

```
python -X utf8 simulations/_f112_open_identity_basis_enum.py
```

Runs N=2 + N=3 in under 5 seconds; N=4 takes ~25 seconds. N=5 (524k pairs) is feasible but takes several minutes; N=6 (8.4M pairs) requires algorithmic optimization.

## Related

- [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md): the parent theorem; this writeup fills the "Step 5 extension to non-Hermitian H" open item listed there.
- [F112_HARDWARE_LENS_KINGSTON.md](F112_HARDWARE_LENS_KINGSTON.md): hardware-side application of F112 that surfaced the first structural counterexample to the broader empirical envelope (Z-drive + σ⁻ T1). The non-Hermitian extension proven here is orthogonal to that counterexample (the counterexample concerns bit_b-mixed c with Hermitian H, not non-Hermitian H with bit_b-homogeneous c).
- [LindbladBitBPiBalance.cs](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs): the typed Tier1Derived C# Claim; this writeup informs the NonHermitianExtension inspectable.
