# F112 Non-Hermitian Extension: Basis-Enumeration Proof at N=2, 3, 4, 5

**Status:** Constructive proof of the F112 non-Hermitian extension at N ≤ 5 via Pauli-basis enumeration. Reduces the open identity `Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0` to a bilinear basis-spanning check; all 559,912 distinct upper-triangular Pauli-string pairs across N=2, 3, 4, 5 give Im < 1e-10 (with N=2..4 bit-exact 0.0e+00; N=5 verified by both Python and C# pipelines, see below). Universal-N lift remains open for N ≥ 6.
**Date:** 2026-05-26 (initial N≤4); Welle 10b 2026-05-26 (extended to N=5)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Scripts:** Python [`simulations/_f112_open_identity_basis_enum.py`](../simulations/_f112_open_identity_basis_enum.py) (N=2..5) + C# [`compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs`](../compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs) (N=2..5 via SLOW_F112-tagged test, parallelized via Parallel.For/ForEach)
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

**F112 non-Hermitian extension is proven at N=2, 3, 4, 5** by the basis-spanning argument: the open identity Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 holds on every basis pair (bit-exact at N=2..4, < 1e-10 at N=5), hence by bilinearity it holds for every pair of Hermitian H_re, H_im at these N. The F112 Tier1Candidate non-Hermitian scope is upgraded to:

- **Tier1Derived at N=2, 3, 4, 5** (constructive enumeration proof, 559,912 pairs total).
- **Tier1Candidate for N ≥ 6** (empirical-only, pending higher-N enumeration via sparse Pauli representation or structural lifting argument; ~8.4M pairs at N=6).

## Lifting to general N

Three possible routes to lift the at-N proofs to universal N:

1. **Brute-force enumeration at higher N.** N=5 was closed in Welle 10b (Python 90.7 min + C# 2 h 45 m, ~16 GB working memory). N=6 is 8.4M pairs × 4^12 = 16M per-pair Frobenius elements, requiring sparse Pauli representation (L_P has exactly 4^N nonzero entries per row) to fit a feasible cache. Without sparse rep, N=6 cache alone would be ~1 TB.

2. **Structural reason that F ≡ 0 algebraically.** The vanishing is bit-exact at every tested N and across every bit_b-parity cell, suggesting a universal symmetry. Candidate explanations:
   - A residual dagger/anti-Hermitian property of L_{H,-i} that imposes ⟨A, B⟩ ∈ ℝ for A, B in the Π −i eigenspace
   - A cancellation mechanism between the (σ_α, σ_β) and (σ_β, σ_α) contributions via the Π conjugation orbit structure
   - Reduction to F38/F61/F63 Π² sector decomposition

3. **Inductive argument N → N+1.** Build the chain-of-(N+1) Pauli basis from the chain-of-N basis by tensoring with one of {I, X, Y, Z}; show that if F = 0 on the N-basis pairs, the induced N+1 pairs also satisfy F = 0.

Each route closes the open identity universally. Route 2 (structural reason) is the most elegant; route 1 extends the constructive proof envelope; route 3 is the cleanest if a per-letter induction works out.

## Implications

If the open identity holds universally (which the N ≤ 5 enumeration strongly suggests, 559,912 pairs all giving Im below threshold), then:

- **F112 non-Hermitian extension is Tier1Derived for all N**, with the algebraic proof: bilinearity + antisymmetry reduce to the basis identity; basis identity holds by `F(σ_α, σ_β) = 0` for every Pauli string pair.
- **The polarity_coordinates_from_L diagnostic** is a structural witness for L not in the Lindblad form `−i[H, ·] + Σ γ_k np.kron(c_k, c_k^*)` with bit_b-homogeneous c. This is true regardless of whether H is Hermitian or non-Hermitian.
- **The bit_b Z₂-axis carries an additional Tier1Derived theorem** beyond F108 Parts 1/2/3 (palindrome closure of bilinears) and the Hermitian-H F112. Three independent typed theorems on the shared bit_b foundation.

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
