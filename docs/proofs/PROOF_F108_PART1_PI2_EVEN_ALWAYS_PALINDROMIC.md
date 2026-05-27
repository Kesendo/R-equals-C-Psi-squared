# PROOF F108 Part 1: Π²-Even Hamiltonians Always Admit an Exact Palindrome Operator

**Status:** Tier 1 derived (closed-form via Π_5bilinear phase-variant Π operator + F1-style algebra). Acts as the base claim for the F108 family: Parts 2 and 3 are Klein-V₄ corollaries via Hadamard transport and D-conjugation respectively (`PROOF_F108_KLEIN_V4_EQUIVALENCE.md`, Welle 14).
**Date:** 2026-05-25 (Part 1 direct proof); 2026-05-27 (Klein-V₄ corollaries Parts 2, 3 added).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (the original F1 palindrome theorem for truly Heisenberg H)
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (truly criterion for k-body Pauli terms)
- [NON_HEISENBERG_PALINDROME](../../experiments/NON_HEISENBERG_PALINDROME.md) (Π-family classification: P1, P4, alternating, non-local)
- [`compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs) (Π_5bilinear builder)
- [PROOF_F108_KLEIN_V4_EQUIVALENCE.md](PROOF_F108_KLEIN_V4_EQUIVALENCE.md) (Welle 14: F108 Parts 2, 3 as Klein-V₄ corollaries of Part 1; precise mechanism documentation)

## Abstract

F1 closes the palindromic spectrum for the canonical "truly Heisenberg" Hamiltonians: those where the canonical Π operator's conjugation aligns cleanly with the dynamics. The natural follow-up question was whether F1's closure extends to a wider Hamiltonian family. Specifically, to Hamiltonians built from any of the five two-site bilinears that the canonical Π operator's Π² eigenvalue rule classifies as even-parity (the bilinears XX, YY, YZ, ZY, ZZ) on any sites with any real bond couplings. The canonical Π itself does not close the palindrome for these Hamiltonians; we needed a different conjugating operator.

This proof writes down that operator. It is a phase variant of the canonical Heisenberg Π, with the same per-site letter permutation (I↔X, Y↔Z) but two sign flips on specific phase arrows. We call it Π_5bilinear because it works precisely on the five even-parity bilinears mentioned above. With this operator in hand, the operator-level palindrome identity holds bit-exactly for any Hamiltonian in the family plus any per-site Z-dephasing rates.

The structural consequence is that the F87 trichotomy collapses on this family: no Π²-even bilinear pair can be F87-hard, because F87-hardness requires the spectrum to break palindromy, and the spectrum here is palindromic by construction. Hardness comes only from Π²-odd or mixed-parity content, never from the Π²-even sector.

The proof's strategy is the F1 algebra applied to Π_5bilinear instead of canonical Π. The anti-commutation argument on the Hamiltonian commutator superoperator goes through verbatim once the phase variant is fixed. The dissipator-side identity works because the Z-dephasing operator commutes with itself, and the phase variant preserves the relevant dissipator structure. Parts 2 and 3 of the F108 trinity extend this to X- and Y-dephasing as Klein-V₄ corollaries via the companion proof.

**Statement (Theorem F108 Part 1):** For any Hamiltonian H built as a linear combination of Π²-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} on N sites with arbitrary real bond coefficients, and Z-dephasing on every site with arbitrary per-site rates γ_l, there exists a per-site Liouville-space operator Π_5bilinear such that

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.

In particular, spec(L) is palindromic around −σ, hence no pure-Π²-even Pauli pair (truly or non-truly) can be F87-hard.

## The Π_5bilinear operator

Π_5bilinear is the per-site Liouville-space automorphism whose action on the four Pauli labels is

  I → +1 · X,    X → −1 · I,    Y → +i · Z,    Z → −i · Y.

In the 4×4 label-basis matrix form on {I, X, Y, Z}:

```
        I   X   Y   Z
   I  [ 0  -1   0   0 ]
   X  [ 1   0   0   0 ]
   Y  [ 0   0   0  -i ]
   Z  [ 0   0   i   0 ]
```

It is the same I↔X, Y↔Z permutation as the canonical Heisenberg Π (the P1 family of [NON_HEISENBERG_PALINDROME](../../experiments/NON_HEISENBERG_PALINDROME.md)), with two phase flips relative to the canonical choice: the X→I arrow has phase −1 (canonical: +1) and the Z→Y arrow has phase −i (canonical: +i). Both are unit-modulus sign flips against the canonical P1.

Key per-site facts:

1. **M is a Liouville-space automorphism, not a Hilbert-space conjugation.** No 2×2 unitary U satisfies U·I·U† = X (since U·I·U† = I always); M acts on the operator (Pauli) algebra directly, not as ρ → U ρ U†.
2. **M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}.** So M⁴ = I and M is order-4. The {I, X} and {Y, Z} 2-cycles each square to a sign, with the I-X cycle squaring to −1 and the Y-Z cycle squaring to +1.
3. **Π_5bilinear is unitary on the d²-dim Liouville space.** Each column has one non-zero entry of unit modulus; columns and rows are pairwise orthogonal.

## Proof

### Step 1: anti-commutation with every Π²-even bilinear

Let Q = M^⊗N be the full N-site Π_5bilinear operator on the 4^N-dim Pauli basis. For every Π²-even 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ}, the commutator superoperator [B, ·] anti-commutes with Q:

  {Q, [B, ·]} = Q · [B, ·] + [B, ·] · Q = 0.

This is verified bit-exactly (residual = 0 at machine precision) at the 2-qubit level for the 16×16 superoperator. The 4 Π²-odd 2-body bilinears {XY, XZ, YX, ZX} produce residual = 8.00 (clean separation; Π_5bilinear does NOT anti-commute with them, as expected).

The anti-commutation result extends from 2 sites to N sites because B is a 2-body operator: [B_(l,l+1), ·] acts non-trivially only on the (l, l+1) factor and as identity on the other N−2 factors. M acts on each identity factor as a permutation with overall sign ±1 that cancels in the {·, ·} bracket.

**Consequence for the Hamiltonian part of L:** L_H = −i [H, ·]. For H = Σ_b α_b B_b a sum of Π²-even bilinears (each B_b in the set above with coefficient α_b ∈ ℝ),

  Q · L_H · Q⁻¹ = Σ_b α_b · Q · (−i [B_b, ·]) · Q⁻¹ = Σ_b α_b · (−i) · (−[B_b, ·]) = −L_H.

The Hamiltonian piece flips sign under conjugation by Q.

### Step 2: per-site identity for the Z-dephasing dissipator

The Lindblad Z-dephasing dissipator on site l is

  D[Z_l] · ρ = γ_l · (Z_l · ρ · Z_l − ρ).

In vec basis: D[Z_l] = γ_l · (Z_l ⊗ Z_l* − I_{d²}). Per site, conjugation by the single-site M satisfies

  M · D[Z] · M⁻¹ = −D[Z] − 2γ · I_4.

This is verified bit-exactly at the 1-qubit level (residual = 0). The mechanism is a diagonal permutation in the Pauli basis. The single-qubit Z-dephasing dissipator, expressed in the {I, X, Y, Z} Pauli basis, is the diagonal super-operator

  D[Z]_pauli = γ · diag(0, −2, −2, 0)

(zeros on the {I, Z} commuting sector, −2γ on the {X, Y} anti-commuting sector). M is the per-site signed permutation with permutation (I↔X, Y↔Z) and phases (+1 on I→X, −1 on X→I, +i on Y→Z, −i on Z→Y); the conjugation M · D · M⁻¹ for a diagonal D in this basis simply permutes the diagonal entries by the underlying letter permutation (the phase factors cancel pairwise: e.g. on the (Y, Z) swap, +i · (entry at Z) · (−i) gives back the original entry). Applying the swap (I↔X, Y↔Z) to diag(0, −2γ, −2γ, 0) yields

  M · D[Z]_pauli · M⁻¹ = γ · diag(−2, 0, 0, −2) = −D[Z]_pauli − 2γ · I_4.

The identity transfers from the Pauli basis to the standard vec basis by the unitary change-of-basis T, since both sides of the identity are unchanged by similarity. Bit-exact at machine precision.

**Consequence for the dissipator part of L:** L_D = Σ_l D[Z_l] is a sum of single-site terms. M acts as a per-site product Q = M^⊗N, so

  Q · L_D · Q⁻¹ = Σ_l (M · D[Z_l] · M⁻¹)
                 = Σ_l (−D[Z_l] − 2γ_l · I_{d²})
                 = −L_D − 2σ · I_{d²},

where σ = Σ_l γ_l.

### Step 3: combining Hamiltonian and dissipator

Adding Steps 1 and 2:

  Q · L · Q⁻¹ = Q · L_H · Q⁻¹ + Q · L_D · Q⁻¹
              = −L_H − L_D − 2σ · I
              = −L − 2σ · I.

This is the operator-level palindrome identity, bit-exact for every H in the Π²-even bilinear family + Z-dephasing on every site.

### Step 4: spectral palindrome and F108 Part 1 corollary

From Q · L · Q⁻¹ = −L − 2σ · I and unitarity of Q:

  spec(L) = spec(Q · L · Q⁻¹) = spec(−L − 2σ · I) = {−λ − 2σ : λ ∈ spec(L)}.

So spec(L) is invariant under λ ↦ −λ − 2σ, i.e. palindromic around −σ.

**F87 corollary:** A Π²-even Pauli pair is F87-hard iff spec(L) breaks palindromy. Since spec(L) is palindromic for every Π²-even H (truly or non-truly), no Π²-even pair can be F87-hard. ∎

## Empirical verification

Bit-exact residual ‖Π_5bilinear · L · Π_5bilinear⁻¹ + L + 2σ · I‖_F = 0 at machine precision, across:

| Setup | N range | residual |
|-------|---------|----------|
| All 9 pure-Π²-even non-truly pairs (single-bilinear + two-term) | N = 3, 4, 5 | 0 |
| 15 random non-uniform-J instances on Π²-even bilinear family | N = 3, 4, 5 | 0 |
| 9 asymmetric J_YZ ≠ J_ZY instances | N = 3, 4, 5 | 0 |
| Pure D[Z]^⊗N dissipator (no Hamiltonian) | N = 3, 4, 5 | 0 |

Reproduction: [`simulations/_f108_part1_pi_family_scan.py`](../../simulations/_f108_part1_pi_family_scan.py); algebraic anti-commutation check in [`simulations/_f108_part1_proof_algebra.py`](../../simulations/_f108_part1_proof_algebra.py); operator output log in `simulations/results/f108_part1_pi_family_scan.txt`.

## Resolution of the previously open ker(M) attempt

Earlier exploration (`F108KerMLemmaExploration.cs`, since removed from the repository) attempted to close F108 Part 1 by proving that every L-eigenvector lies in the kernel of M (the F1 palindrome residual for the canonical Heisenberg Π). Numerical verification falsified that lemma: only 2 of 64 L-eigenvectors at the N=3 YZ+ZY test case satisfied M·v_λ = 0. The follow-up Critical Isospectral Lemma reduction was circular (it restated the spectral palindrome under another name).

The actual mechanism is structurally different: a DIFFERENT per-site Π operator (Π_5bilinear, phase variant of the canonical P1 Π) gives EXACT operator-level palindrome with residual zero, without going through M at all. The canonical Heisenberg Π has M ≠ 0 for Π²-even non-truly H; the Π_5bilinear Π has its own residual (call it M_5bilinear), which is zero for the full Π²-even bilinear family by Steps 1-2 above.

## Significance

F108 Part 1 closes the long-open F87 hardness criterion for Π²-even pairs:

- **F107** (Tier 1 derived): truly classification ⟹ y_par = 0 across all dephase letters.
- **F108 Part 1** (Tier 1 derived, THIS PROOF): no Π²-even pair is F87-hard; equivalently, every Π²-even H admits an exact operator-level palindrome Π.
- **F109** (Tier 1 derived, was modulo F108 Part 1): mother sector Klein (0, 0) soft ⟹ y_par = 1. With F108 Part 1 closed, F109 is fully unconditional Tier 1 derived.

The PROOF STRATEGY is generalizable: the per-site label-permutation algebra (which Π preserves which Pauli-pair sectors with which phase signs) is a complete classification of when palindromes hold. NON_HEISENBERG_PALINDROME's 4-family taxonomy (P1, P4, alternating, non-local Π) is the catalog; F108 Part 1 identifies P1's "5-bilinear phase variant" as the Π operator that covers every Π²-even Hamiltonian.

## Sibling y_par-axis claims

Closed 2026-05-25: F108 Part 1+2+3 (Π²-even palindrome family, Tier1Derived); F109 (MotherSoftYParityOnePurity, Tier1Derived unconditional); F110 (HardCellYInversionPattern, Tier1Candidate). Together F107+F109+F110 close the y_par-axis F87 trichotomy classification.

## Open

- Closed-form derivation of F110 Aspect C exact ratios (42:8 at k=3, 228:0 at k=4) per Pauli-letter combinatorics. F103 Section 5 explicitly lists as open.
- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4.
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).

∎
