# PROOF F112 Cross-Dephase Extension: X- and Y-Dephase Versions via Klein-V₄

**Status:** Tier 1 derived, universal N. Two complementary structural routes (axis-direct re-run of the Welle-11 lemmas, and Hadamard-conjugation transport from F112-Z to F112-X via the Klein-V₄ swap group).
**Date:** 2026-05-27 (Welle 13)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- F38 (Π_d² eigenvalue formula; `docs/ANALYTICAL_FORMULAS.md` F38 entry)
- F63 ([L, Π²] = 0 for Z-dephasing; F63 entry; analogous statements for X- and Y-dephasing follow from F38)
- [F112 Lindblad Π-balance](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) (parent: Hermitian-H 5-step proof for d = Z)
- [F112 non-Hermitian extension](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md) (Welle 11: non-Hermitian H closure for d = Z via two structural lemmas)
- [Klein-V₄ dephase-letter swaps on operator space](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) (Welle 12: Klein-V₄ subgroup {I, D, Q_zx, Q_yx} of U(4^N), where Q_yx = ⊗h is the element this proof also writes **H**; with D · Π_Z · D = Π_Y, Q_zx · Π_Z · Q_zx = Π_X, Q_yx · Π_Y · Q_yx = Π_X)
- Verifier `simulations/f112_klein_v4_cross_dephase_verify.py` (Welle 13)

## Abstract

F112 was proved for Z-dephasing. The natural next question was whether the same polarity-balance identity holds for X-dephasing and Y-dephasing as well. After the Klein-V₄ dephase-swap group landed in Welle 12, the answer looked structurally suggested: the V₄ relates the palindrome operators Π_Z, Π_X, Π_Y (precisely, by the regular V₄-action on the four oriented operators {Π_X, Π_X⁻¹, Π_Y = Π_Z⁻¹, Π_Z}; see (d)). But the group-level relation only carries an *identity* across dephase letters when the swap element is a genuine Hilbert-space lift: true for the Hadamard Q_zx (Z → X), false for D (Z → Y). The job was to write the two routes that actually close it, and to spell out how the hypothesis on the collapse operators rotates with the dephase axis.

Both routes work, and they are independent. The first re-runs the F112 proof with the dephase letter replaced by X or Y throughout. The two structural inputs that change per dephase axis are the Π² eigenvalue formula on Pauli strings (which switches from bit_b to bit_a when the dephase letter switches from Y/Z to X) and the per-site palindrome-operator matrix (a 4×4 signed permutation for every dephase letter, just with different non-zero positions). Otherwise the lemmas re-run verbatim, and the conclusion lands.

The second route is the Klein-V₄ conjugation transport, and it is one-sided. The relevant element is the Hadamard lift Q_zx (a genuine Hilbert-space unitary, U_H^⊗N · ρ · U_H^⊗N in operator space), which rotates the dephasing axis Z → X. Frobenius norms are preserved by unitary conjugation, so conjugating the F112-Z identity by Q_zx gives F112-X for free, mapping the bit_b hypothesis on the collapse operators to the bit_a hypothesis. The other swap element D conjugates Π_Z into Π_Y but does NOT rotate the dephasing axis (D is the transpose superoperator; D · L_Z · D⁻¹ is still a Z-dephasing Lindbladian, namely L_Z with H → −H^T), so D-conjugation does not yield F112-Y. F112-Y comes instead from Route 1, since it shares the bit_b axis with Z.

The hypothesis on the collapse operators carries the axis dependence. For dephase letter Z or Y, the operators need to be bit_b-homogeneous (every Pauli string in the expansion shares the same bit_b parity). For dephase letter X, they need to be bit_a-homogeneous. The axis switches with the dephase letter because Π² classifies by the axis NOT flipped by Π's per-letter swap: Π_X's swap is on the bit_b axis, while Π_Y's and Π_Z's swaps are on the bit_a axis.

The diagnostic reading: F112 is one identity, not three. It holds for any of the three F1 palindrome conventions with the appropriate per-axis homogeneity hypothesis, for any Hamiltonian (Hermitian or non-Hermitian, read as the commutator generator −i[H, ·]; see Scope in (a)), for any number of qubits. The polarity-coordinates primitive reads the same diagnostic information regardless of which dephasing convention the user invokes, as long as the conventions stay matched (the M_d-anti decomposition must use the Π_d that corresponds to the L_d under analysis).

## Introduction

**The motivating question.** [F112 Lindblad Π-balance](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) closed the polarity-balance identity for d = Z, the canonical Z-dephasing convention used throughout the codebase. The natural follow-up: does the same identity hold for d = X and d = Y, or is there something Z-specific about the construction? Once the Klein-V₄ swap group on operator space landed in [Welle 12 Task 2](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), the three identities looked closely related (the V₄ conjugates the palindrome operators into one another). The remaining job was to write down the two routes that establish this rigorously; and, as it turned out, to find that only one of the two swap directions (the Hadamard Q_zx, Z → X) actually transports an identity, while Y must be re-derived directly.

**The empirical anchor.** [simulations/f112_klein_v4_cross_dephase_verify.py](../../simulations/f112_klein_v4_cross_dephase_verify.py) verified the F112-X and F112-Y identities at N = 2, 3 with the appropriate axis_d-homogeneous c construction; the asymmetry vanished to ≤ 4.3e-14 for every configuration tested (literally 0.0 for the F112-X-direct and structural-lemma checks; floating-point noise ~1e-14 for the F112-Y and Hadamard routes, which carry the −i Π-phase). The commutator identity additionally holds at N = 4 (Hermitian and non-Hermitian H, checked 2026-06-20). The structural argument lifts the per-N verification to universal N via either route.

**What this proof closes.** Two complementary closures:

1. **Direct re-run (Route 1).** The Welle 11 two-lemma proof (Lemma N-A: `‖L_{σ,−i_d}‖² = 4^N` for σ axis_d-odd; Lemma N-B: cross-pair Pauli-basis disjointness) uses two structural inputs: F38's Π_d² eigenvalue formula on Pauli strings, and the per-site π_d_local matrix structure. F38 changes axis per d (bit_b for Y/Z, bit_a for X) but otherwise has the same algebraic form; π_d_local is a 4×4 signed permutation for every d. Re-running the two lemmas with d = X or d = Y produces the per-d identity directly.

2. **Klein-V₄ conjugation transport (Route 2), one-sided.** The Welle 12 Task 2 result establishes the Klein-V₄ subgroup `{I, D, Q_zx, Q_yx}` of `U(4^N)` with `D · Π_Z · D = Π_Y` and `Q_zx · Π_Z · Q_zx = Π_X` (each operator unitary, so Frobenius norms and asymmetry values are preserved exactly). Only **Q_zx** transports an identity between dephase letters: it is the Hadamard Hilbert-space lift, so it rotates the dephasing axis Z → X, and conjugating the F112-Z identity by Q_zx yields F112-X (mapping the bit_b hypothesis on c to bit_a). **D does NOT yield F112-Y**: D is the transpose superoperator, which conjugates Π_Z → Π_Y but leaves the dephasing axis on Z (D · L_Z · D⁻¹ is still a Z-dephasing Lindbladian), so it produces only F112-Z in a crossed (L_Z, Π_Y) disguise, not F112-Y. F112-Y is supplied by Route 1 instead (it shares the bit_b axis with Z); see (d).

**Diagnostic consequence.** F112 is now a universal Lindblad-channel polarity-balance identity: for any of the three F1 palindrome conventions (Z, X, Y), any Hamiltonian H, any axis-adapted homogeneity hypothesis on c, the asymmetry vanishes bit-exactly for any N. The `polarity_coordinates_from_L` primitive can be run under any dephasing convention with the same diagnostic interpretation. Asymmetry ≠ 0 localises the failure: (i) c violates the axis_d-homogeneous hypothesis, (ii) L not in Lindblad form, or (iii) the conventions are crossed (computing M_d with Π_d' for d' ≠ d).

## (a) Statement

For each dephase letter d ∈ {X, Y, Z}, write
- Π_d for the per-d palindrome conjugation operator on operator space (linear, unitary, order 4; see `PiOperator.cs` / `framework.symmetry.build_pi_full`), and
- L_d = -i[H, ·] + Σ_k γ_k · (np.kron(c_k, c_k^*) − ½·{c_k^† c_k, ·}) for a standard Lindblad-form (GKSL) Liouvillian. The anticommutator term is Π²-even and contributes only to M_zero (the Π²-even part of M_d, orthogonal to the M_d_anti sector that carries the asymmetry), so it does not affect the asymmetry (it is present in the verifier; omitting it changes nothing).

Define the per-d polarity asymmetry of M_d := Π_d · L_d · Π_d⁻¹ + L_d + 2σ_c · I (where σ_c := Σ_k γ_k is the total dephasing rate, which recenters the F1 palindrome about the origin; this constant is distinct from the Pauli-string σ used below, and matches the `sigma` argument of `polarity_coordinates_from_L`) as

```
asymmetry_d(L_d) := ‖M_d_plus_half‖² − ‖M_d_minus_half‖²,
```

where M_d_plus_half / M_d_minus_half are the standard Π_d-eigenvalue ±i projections of M_d_anti (see `polarity_coordinates_from_L`).

**Theorem (F112 universal cross-dephase).** Let H be any Hamiltonian (Hermitian or non-Hermitian), N ≥ 1, γ_k ∈ ℂ. Define the per-d hypothesis on each c_k:

- d = Z: each c_k is **bit_b-homogeneous** (every Pauli string σ in c_k's expansion satisfies bit_b(σ) = #Y(σ) + #Z(σ) mod 2 = const).
- d = Y: each c_k is **bit_b-homogeneous** (same axis as Z).
- d = X: each c_k is **bit_a-homogeneous** (every Pauli string σ satisfies bit_a(σ) = #X(σ) + #Y(σ) mod 2 = const).

Then asymmetry_d(L_d) = 0 bit-exactly for the corresponding d.

**Equivalent compact statement.** Define axis_d := bit_b for d ∈ {Y, Z} and axis_d := bit_a for d = X. Then F112-d holds for any (H, {c_k}) with each c_k axis_d-homogeneous, both Hermitian and non-Hermitian H, universal N.

**Scope (commutator generator only).** As with the non-Hermitian parent ([F112 non-Hermitian extension](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), corrected 2026-06-20), "non-Hermitian H" here means the **commutator** superoperator −i[H, ·] for an arbitrary matrix H, a structural ‖·‖² identity. It is NOT the physical generator of PT-symmetric / gain-loss / post-selection dynamics, −i(Hρ − ρH†) = −i[A, ρ] + {B, ρ} for H = A + iB; that anticommutator part breaks the balance (asymmetry ≈ 240–1500 at N = 2, 3 for X- and Y-dephasing, mean ≈ 54/200 over random fixed-norm configs; the same boundary the parent documents for Z). The two generators coincide exactly when H is Hermitian.

## (b) Asymmetry of hypotheses: why F112-Y uses bit_b but F112-X uses bit_a

Per F38 (PiOperator.SquaredEigenvalue):

```
Π_Z² acts on Pauli string σ by scalar (−1)^bit_b(σ),
Π_Y² acts on Pauli string σ by scalar (−1)^bit_b(σ),       (same as Π_Z²)
Π_X² acts on Pauli string σ by scalar (−1)^bit_a(σ).       (different axis)
```

The asymmetry between Y/Z and X is structural: the Π_d operator preserves the axis NOT flipped by Π_d's per-letter swap. For Π_Z, the per-letter swap is I↔X (flipping bit_a), so Π_Z preserves bit_b; squaring then collapses Π_Z² to a scalar (−1)^bit_b on each Pauli string. Same for Π_Y (also flips bit_a, preserves bit_b). For Π_X the per-letter swap is I↔Z (flipping bit_b), so Π_X preserves bit_a and Π_X² scales by (−1)^bit_a.

In the Welle-11 proof's Step (b) (BitB-parity restriction), the "bit_?-odd" sector is the only sector contributing non-trivial Π-eigenvalue ±i content. For d ∈ {Y, Z} that sector is bit_b-odd; for d = X it is bit_a-odd. The natural homogeneity hypothesis on c follows: each c_k must be supported on a single bit_d-parity to keep `np.kron(c_k, c_k^*)` in the Π_d²-conjugation +1 eigenspace (Step 2 of the parent Hermitian-H proof) and equivalently in the {+1, −1} sector of Π_d (Step 3).

## (c) Proof Route 1: Direct re-run of Welle-11 lemmas for d ∈ {X, Y}

Here `L_{σ,±i_d}` denotes the Π_d-eigenvalue ±i spectral components of the single-Pauli Lindbladian generator L_σ (the σ subscript is a Pauli string; notation inherited from [`PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md)). The Welle-11 proof of F112 non-Hermitian extension closes the per-pair identity F(σ_α, σ_β) := Im⟨L_{σ_α,−i_d}, L_{σ_β,−i_d}⟩ = 0 via two lemmas (Lemma N-A: ‖L_{σ,−i_d}‖² = 4^N for σ axis_d-odd; Lemma N-B: ⟨L_{σ_α,−i_d}, L_{σ_β,−i_d}⟩ = 0 for σ_α ≠ σ_β both axis_d-odd). The proof of each lemma reduces to per-position checks on the 4^N × 4^N matrix of L_σ in the Pauli basis, using two structural inputs:

1. **F38 Π_d² eigenvalue formula on Pauli strings** (axis = bit_b for d ∈ {Y, Z}, axis = bit_a for d = X).
2. **Pauli-basis matrix-support disjointness** of L_σ and Π^m L_σ Π^{−m} for m ∈ {0, 1, 2, 3}.

Input 1 changes axis per d but otherwise has the same algebraic form. Input 2 depends on the per-site π_d_local matrix structure, which is a 4×4 signed permutation for every d ∈ {X, Y, Z} (see `PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md` Step 1 for the explicit 4×4 matrices). Specifically, π_d_local has exactly four non-zero entries per row/column with magnitudes in {1, i, −i}, the same sparse structure for all three dephase letters.

**Lemma N-A^d (Diagonal-Norm, per-d).** For any axis_d-odd Pauli string σ at chain length N, ‖L_{σ,−i_d}‖² = 4^N exactly.

**Proof (direct re-run of Welle-11 Lemma N-A).** Steps A.1, A.2, A.3, A.4 of the Welle-11 proof use only:
- Step A.1: ‖L_σ‖² = 2 · 4^N. This is independent of the dephase letter: it is the standard Pauli-anticommutation count (exactly half of the 4^N Pauli strings anticommute with σ for σ non-identity; each non-zero matrix entry has magnitude 2). The result is `2 · 4^N` regardless of d.
- Step A.2: Π_d² L_σ Π_d^{−2} = −L_σ for σ axis_d-odd. This is F38 applied to L_σ as a superoperator: per the calculation in `PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md` section (b), Π_d² L_σ Π_d^{−2} = (−1)^{axis_d(σ)} · L_σ; substituting axis_d(σ) = 1 (odd sector) gives −L_σ for every d.
- Step A.3: ⟨L_σ, Π_d L_σ Π_d^{−1}⟩_F = 0 for σ axis_d-odd. The Welle-11 calculation uses four facts (numbered 1–4 in section (c) of the parent proof) about the per-site Klein-index map π_d_local; each of these facts holds verbatim for any d ∈ {X, Y, Z}, modulo the substitution of axis_d for bit_b:
  - Fact 1 (π_d involution on Klein indices: π_d² = identity on per-site Klein), verified for all three d in `PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md` Step 5 via d² = h² = q_zx² = I.
  - Fact 2 (π_d homomorphism on XOR up to a fixed shift): for d = Z the shift is the all-X string (1, 0)_l per site; for d = Y the shift is the same (Y-dephase also flips bit_a per letter); for d = X the shift is the all-Z string (0, 1)_l per site (X-dephase flips bit_b per letter).
  - Fact 3 (π_d(σ ⊕ α') = π_d(σ) ⊕ α' ⊕ shift; π_d(σ) ⊕ shift = σ): follows from Fact 2 plus involutivity.
  - Fact 4 (anticommutation symplectic shift): for σ axis_d-odd, anticomm(σ, σ_{π_d(α')}) = COMMUTE(σ, σ_{α'}). The argument is the same Klein-symplectic-form bilinearity used in the Z case, now with the role of bit_b played by axis_d.

  By the same support-disjointness argument as in Welle-11, the non-zero entries of M(L_σ) and Π_d M(L_σ) Π_d^{−1} sit on the same shifted-diagonal {(σ ⊕ α', α')} but with mutually exclusive conditions on α' (anticomm-with-σ versus commute-with-σ for axis_d-odd σ). The entry-wise product vanishes, giving ⟨L_σ, Π_d L_σ Π_d^{−1}⟩_F = 0.
- Step A.4: From A.1–A.3 and the L_σ = L_σ^{+i_d} + L_σ^{−i_d} decomposition (since Π_d² L_σ Π_d^{−2} = −L_σ kills the {+1, −1} components), ‖L_σ^{±i_d}‖² = (1/2) · ‖L_σ‖² = 4^N each. ∎

**Lemma N-B^d (Off-Diagonal-Orthogonality, per-d).** For σ_α ≠ σ_β both axis_d-odd Pauli strings at chain length N, ⟨L_{σ_α,−i_d}, L_{σ_β,−i_d}⟩_F = 0 exactly.

**Proof (direct re-run of Welle-11 Lemma N-B).** The Welle-11 spectral-projector identity gives

```
⟨L_{σ_α,−i_d}, L_{σ_β,−i_d}⟩_F = (1/4) Σ_{m=0}^{3} (−i)^{−m} · ⟨L_{σ_α}, Π_d^m L_{σ_β} Π_d^{−m}⟩_F.
```

For each m, the inner product on the right is shown to vanish for σ_α ≠ σ_β via the same shifted-diagonal-support disjointness used in Step B.1, B.2 of the parent proof:

- m = 0: M(L_{σ_α}) supported on {(σ_α ⊕ α', α')}; M(L_{σ_β}) on {(σ_β ⊕ α', α')}. Equating β' requires σ_α = σ_β, contradicting α ≠ β.
- m = 2: Π_d² acts as a diagonal phase per Pauli string (F38), so Π_d² M Π_d^{−2} = D · M · D^{−1} (diagonal similarity). Support is identical to M, hence m = 0 conclusion applies.
- m = 1, m = 3: the support of Π_d M(L_σ) Π_d^{−1} is again {(σ ⊕ α', α')} (the shifted-diagonal structure is preserved; only the "condition on α'" changes from "anticomm" to "commute" per Fact 4 with axis_d in place of bit_b). Disjointness with M(L_{σ_α}) still requires σ_α = σ_β.

Sum over m gives 0. ∎

**Main theorem F112-d (consequence of A^d + B^d + trivial axis_d-even case).** For any Pauli pair (σ_α, σ_β) at length N,

```
F_d(σ_α, σ_β) := Im⟨L_{σ_α,−i_d}, L_{σ_β,−i_d}⟩_F = 0.
```

- axis_d-even σ_α or σ_β: Π_d² σ Π_d^{−2} = +σ, hence L_σ has no ±i_d components, so L_{σ,−i_d} = 0. F_d = 0 trivially.
- axis_d-odd σ_α = σ_β: ⟨L_{σ,−i_d}, L_{σ,−i_d}⟩_F = ‖L_{σ,−i_d}‖² = 4^N ∈ ℝ by Lemma N-A^d. Im = 0.
- axis_d-odd σ_α ≠ σ_β: ⟨ , ⟩ = 0 exactly by Lemma N-B^d.

By bilinearity + Pauli-basis spanning, F_d(H_re, H_im) = 0 for any Hermitian H_re, H_im. Combined with the Hermitian-H parent proof (Steps 1–4 of `PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md` work for any d ∈ {X, Y, Z} because they use only F38 and the spectral structure of Π_d; Step 5's dagger argument is independent of d), F112-d holds universal N for any H, Hermitian or non-Hermitian. ∎

## (d) Proof Route 2: Hadamard transport from F112-Z to F112-X via Q_zx

The Klein-V₄ swap operator Q_zx ∈ U(4^N) (Welle 12) realizes Π_Z ↔ Π_X by conjugation: Q_zx · Π_Z · Q_zx⁻¹ = Π_X. Crucially, Q_zx is the operator-space lift of the per-site Hadamard rotation U_H = (1/√2) · [[1, 1], [1, −1]] acting by U_H^⊗N · ρ · U_H^⊗N on the Hilbert space (verified in the controller's Welle-12 reasoning + this proof's verifier).

Hadamard acts on per-site Pauli operators as

```
U_H · X · U_H = Z,    U_H · Z · U_H = X,    U_H · Y · U_H = −Y,    U_H · I · U_H = I.
```

In Klein-index coordinates:

```
X = (a, b) = (1, 0)  ↔  Z = (0, 1)        (Hadamard swaps bit_a ↔ bit_b)
Y = (1, 1)  → ±Y (preserves both)
I = (0, 0)  → I (preserves both)
```

**Consequence: bit_b ↔ bit_a swap on Pauli strings under Hadamard^⊗N conjugation.** For any Pauli string σ on N sites, write σ' := U_H^⊗N · σ · (U_H^⊗N)^†. Then bit_a(σ') = bit_b(σ) and bit_b(σ') = bit_a(σ). In particular, a Pauli sum c = Σ_α a_α σ_α is bit_b-homogeneous iff c' = U_H^⊗N · c · (U_H^⊗N)^† is bit_a-homogeneous.

**Claim B (Hadamard transport).** If F112-Z holds for (H, {c_k}) with Hermitian H + bit_b-homogeneous c_k, then F112-X holds for (H', {c'_k}) := (U_H^⊗N · H · (U_H^⊗N)^†, {U_H^⊗N · c_k · (U_H^⊗N)^†}) with H' Hermitian + bit_a-homogeneous c'_k.

**Proof of Claim B.** Let U := U_H^⊗N (Hilbert-space unitary). Conjugation by U preserves Hermiticity: H' = U · H · U^† is Hermitian iff H is. The bit_b-homogeneity of c maps to bit_a-homogeneity of c' = U · c · U^† by the Klein-index swap above.

Compute the Lindbladian under (H', c'_k). The standard Lindblad form gives

```
L_X-style(H', {c'_k}) = −i[H', ·] + Σ_k γ_k · (np.kron(c'_k, c'_k^*) − (1/2)·{c'_k^† c'_k, ·}_anti).
```

In the standard vec basis, conjugation by U on operators corresponds to a left-Kronecker action by U on vec(ρ): vec(U ρ U^†) = (U ⊗ U^*) · vec(ρ). The transformation Liouvillian L_X-style(H', c'_k) is then the conjugate of L_Z-style(H, c_k) by the operator-space unitary

```
U_op := U ⊗ U^*.
```

The Π_X palindrome operator on operator space is Π_X = Q_zx · Π_Z · Q_zx^† (Welle 12). Therefore the M_d_anti decomposition under (H', c'_k, Π_X) is unitarily equivalent to the M_d_anti decomposition under (H, c_k, Π_Z) by simultaneous conjugation by Q_zx (which factors as U_op times a basis-permutation matching the (I, X, Z, Y) basis to the standard Pauli-string enumeration; the basis convention detail is spelled out in `PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md` Convention note).

The asymmetry is a Π-eigenspace Frobenius norm difference. Under simultaneous unitary conjugation of L and Π by Q_zx, the Π-eigenspaces map to Π-eigenspaces with the same Π-eigenvalues, and Frobenius norms are preserved. Hence asymmetry_X(L_X-style(H', c'_k)) = asymmetry_Z(L_Z-style(H, c_k)) = 0 by F112-Z. ∎

**Remark on the "operator-space-only" vs "Hilbert-space lift" distinction (controller's pre-dispatch concern).** The Klein-V₄ swap group {I, D, Q_zx, H} on operator space has different lift properties:

- **Q_zx** = U_op = U ⊗ U^* IS the lift of a Hilbert-space unitary U (the Hadamard ⊗N). So Q_zx · L · Q_zx^† maps Lindblad-form L to Lindblad-form L' with rotated (H', c'_k) = (U H U^†, {U c_k U^†}).
- **D** = ⊗ diag(1, 1, 1, −1) is the **transpose superoperator** ρ ↦ ρ^T (it negates Y and fixes I, X, Z). It has no Hilbert-space unitary lift (a V with V · Y · V^{-1} = −Y, V · X · V^{-1} = X, V · Z · V^{-1} = Z is contradictory by Pauli algebra, since Y = iXZ forces V · Y · V^{-1} = +Y). But "no unitary lift" does NOT mean "leaves Lindblad form": D · L_Z · D⁻¹ is the **Z-dephasing Lindbladian with H → −H^T** (a valid, completely-positive generator, verified by a Choi-positivity check). What D fails to do is rotate the dephasing axis: it conjugates Π_Z → Π_Y but keeps the dissipator on Z, so it never reaches a Y-dephasing channel.
- **H** = ⊗ h (the X↔Z basis-index permutation) is also operator-space-only. (H = Q_yx in the Welle-12 naming.)

**What the "Klein-V₄" really acts on.** The matrix group {I, D, Q_zx, Q_yx} is a genuine Klein-V₄ (abelian, every non-identity element order 2, Q_yx = D · Q_zx), but its conjugation action is *not* "three swaps of the three letters {X, Y, Z}" (three transpositions would generate S₃, and the induced bare-letter map is in fact ill-defined for Q_zx and Q_yx). Because **Π_Y = Π_Z⁻¹** exactly (the Y- and Z-palindromes are mutual inverses at all N, since their per-site phase conventions are complex-conjugate on an order-4 operator), the three letters carry only **two** independent palindrome operators. The V₄ acts as the **regular (free, transitive) action on the four oriented operators {Π_X, Π_X⁻¹, Π_Y = Π_Z⁻¹, Π_Z}**: a Klein-four torsor, with the identity-letter slot I ↦ Π_X⁻¹, X ↦ Π_X, Y ↦ Π_Y, Z ↦ Π_Z. Each generator is a fixed-point-free double transposition (e.g. D = (Π_X Π_X⁻¹)(Π_Y Π_Z)); the ambient group is the signed-permutation group B₂ ≅ D₄, not S₃. This is precisely why "Klein-V₄" is the right name: the four oriented palindrome operators *are* a copy of the Pauli Klein group {I, X, Y, Z}. (Note: the three *dephasing-diagonal* operators Q_v carry the complementary S₃ structure; the V₄-vs-S₃ difference is the Π-vs-Q distinction, cf. the `one_diagonal_mirror_group` arc.)

This means Route 2 (Hadamard transport) gives F112-Z → F112-X cleanly, but does NOT give F112-Z → F112-Y by an analogous unitary lift. F112-Y requires Route 1 (direct re-run of Welle-11 lemmas in the bit_b axis with d = Y phase substituted for d = Z phase). The two routes are independent and complementary: Route 1 covers all three d ∈ {X, Y, Z} structurally; Route 2 provides the bonus insight that the Hadamard-conjugated Z-dephase config IS a natural F112-X config.

In particular, **F112-Y does NOT follow from D-conjugation of L_Z**, but not for the reason of leaving Lindblad form. D · L_Z · D⁻¹ is a perfectly valid Z-dephasing Lindbladian (= L_Z with H → −H^T); conjugating the F112-Z identity by D therefore reproduces F112-Z in a crossed (L_Z-style, Π_Y) disguise, a Z-dephasing channel paired with the Y-palindrome, not an F112-Y instance, which requires a genuine **Y-dephasing** L_Y. The correct path for F112-Y is the direct axis re-run of Welle-11 (Route 1 above), and it is immediate there: Y shares the bit_b axis with Z, so any bit_b-homogeneous c that satisfies F112-Z satisfies F112-Y directly (verified asymmetry_Y(L_Z) ≤ 3e-14).

## (e) Verification

The verifier `simulations/f112_klein_v4_cross_dephase_verify.py` (Welle 13) confirms all four claims numerically at N = 2 and N = 3:

| Claim | N=2 max\|asymmetry\| or \|Im F\| | N=3 max\|asymmetry\| or \|Im F\| |
|---|---|---|
| F112-Y direct (Hermitian H + bit_b-homo c, Π_Y) | 3.6e-15 | 4.3e-14 |
| F112-Y direct (non-Hermitian H + bit_b-homo c, Π_Y) | 1.4e-14 | 1.3e-32 |
| F112-X direct (Hermitian H + bit_a-homo c, Π_X) | 0.0 | 0.0 |
| F112-X direct (non-Hermitian H + bit_a-homo c, Π_X) | 0.0 | 0.0 |
| F112-X via Hadamard (Z-config rotated → Π_X measurement) | 3.6e-15 | 1.9e-30 |
| Lemma N-A^Z (axis=bit_b, n_odd=8/32, n_off=56/992): max deviation | 0.0 | 0.0 |
| Lemma N-A^Y (axis=bit_b, n_odd=8/32, n_off=56/992): max deviation | 0.0 | 0.0 |
| Lemma N-A^X (axis=bit_a, n_odd=8/32, n_off=56/992): max deviation | 0.0 | 0.0 |
| Lemma N-B^Z (axis=bit_b): max \|<L_{α,-i}, L_{β,-i}>\| | 0.0 | 0.0 |
| Lemma N-B^Y (axis=bit_b): max \|<L_{α,-i}, L_{β,-i}>\| | 0.0 | 0.0 |
| Lemma N-B^X (axis=bit_a): max \|<L_{α,-i}, L_{β,-i}>\| | 0.0 | 0.0 |
| Main theorem F_d for d ∈ {Z, Y, X}: max \|Im F\| over 256/4096 pairs | 0.0 | 0.0 |

Additionally, the edge-case table at N = 2, 3 exhibits the structural axis breakage:

| c structure | F112-X (axis=bit_a) | F112-Y (axis=bit_b) | F112-Z (axis=bit_b) |
|---|---|---|---|
| Single-Pauli c (any letter, trivially homo on both axes) | 0 | 0 | 0 |
| Mixed bit_b + mixed bit_a, complex coefs (e.g. X+Z site) | non-zero | non-zero | non-zero |
| Mixed bit_b, homogeneous bit_a (e.g. X+Y site) | **0** | non-zero | non-zero |
| Homogeneous bit_b, mixed bit_a (e.g. Y+Z site) | non-zero | **0** | **0** |

Each row's pattern matches the per-d hypothesis: F112-d gives 0 iff axis_d-homogeneous; non-zero otherwise. (At N = 2, 3 the breakage magnitude is ~225 and ~901 respectively.)

## (f) Implications

1. **F112 is universal in dephase letter:** the bit_b-axis closure for d ∈ {Y, Z} and the bit_a-axis closure for d = X are all Tier 1 derived, universal N, both Hermitian and non-Hermitian H.

2. **Klein-V₄ equivariance is partial but real:** the Hadamard subgroup {I, Q_zx} of {I, D, Q_zx, Q_yx} cleanly transports F112-Z to F112-X (Route 2). The D and H = Q_yx involutions are operator-space-only (no Hilbert-space unitary lift), so they do not rotate the dephasing axis of L: D is the transpose, and D · L_Z · D⁻¹ stays a Z-dephasing Lindbladian (H → −H^T), never a Y-dephasing one. They intertwine the Π's but cannot carry L between dephase letters. Hence F112-Y must be proven via Route 1 (direct axis re-run); it is immediate there because Y shares the bit_b axis with Z.

3. **The `Pi2KleinV4DephaseSwapGroup` typed Claim's significance docstring was softened in Welle 15** (it now states D and H = Q_yx are operator-space-only and do NOT transport L between dephase letters). The reading it now carries:
   - F1-family identities that depend only on (a) the Π-eigenvalue / Π²-eigenvalue spectrum and (b) Pauli-basis support structure transfer between dephase letters via Route 1 (direct axis re-run); these include F112-d, F108-d (conjecturally), and any other F-formula whose proof reduces to F38 + matrix-support disjointness.
   - Identities depending on the Lindblad-form L itself transfer only via the Hadamard subgroup {I, Q_zx} (Route 2), and only between (Z, X), not (Z, Y) or (Y, X).
   - The D and H operators retain mathematical value as Π-intertwiners and witness the Klein-V₄ structure on dephase letters, but should not be misread as L-transporting unitaries.

4. **For typed-knowledge consumers:** F112-X and F112-Y were added as separate typed Claims in **Welle 15**: `LindbladBitAPiBalance` (F112-X, the BitA-axis twin) and `LindbladBitBPiYBalance` (F112-Y, the BitB-axis Π_Y sibling of `LindbladBitBPiBalance`). This flipped `LindbladBitBPiBalance.BitATwinStatus` from BitBSpecific to **Filled** (the code names this proof as the trigger).

## (g) Open questions and follow-ups

- **F108 Parts 2 and 3 closed by Welle 14 (see `PROOF_F108_KLEIN_V4_EQUIVALENCE.md`):** the Route 1 / Route 2 analysis applied to F108 (2026-05-27) yields a mixed result. Part 1 ↔ Part 3 follows by operator-space D-conjugation (D · Π_5b(Z) · D = Π_5b(Y) bit-exact; bilinear set fixed on bit_b axis); this is a CLEAN Klein-V₄ corollary on the operator side. Part 1 ↔ Part 2 follows by Hilbert-space Hadamard transport (Route 2) on L itself; the operator-space Q_zx does NOT swap Π_5b(Z) ↔ Π_5b(X), so Route 1 in its strict form (operator-space conjugation by Klein-V₄) is partial. Welle 14 thereby confirmed that for F108, as for F112-Y, the Route 1 / Route 2 distinction matters; Klein-V₄ equivariance on Π_5bilinear is REAL but PARTIAL ({I, D} subgroup only on the operator side).

- **Typed Claim split for F112-X and F112-Y:** ✓ completed in Welle 15 (`LindbladBitAPiBalance` for F112-X, `LindbladBitBPiYBalance` for F112-Y; `LindbladBitBPiBalance.BitATwinStatus` now Filled). No further deferred work on the typing.

## Status

**F112 cross-dephase (X- and Y-dephase versions) is Tier 1 derived for both Hermitian and non-Hermitian H, universal in N.** Two complementary structural routes:
- Route 1 (direct re-run of Welle-11 lemmas with axis_d substituted for bit_b): proves F112-X, F112-Y, and re-derives F112-Z under a unified axis_d framework.
- Route 2 (Hadamard transport via the Q_zx operator-space lift of U_H^⊗N): proves F112-X from F112-Z via the Klein-V₄ Hadamard subgroup, exhibiting the explicit Π_Z ↔ Π_X intertwiner.

The D involution (Welle 12's Z↔Y swap) is the transpose superoperator: it conjugates Π_Z → Π_Y but does not rotate the dephasing axis (D · L_Z · D⁻¹ is still Z-dephasing, = L_Z with H → −H^T, a valid Lindbladian, not a "category mistake", just not a Y-dephasing channel). So F112-Y cannot be "obtained for free" from F112-Z via D-conjugation; it follows from Route 1, sharing the bit_b axis with Z.
