# PROOF F112: Lindblad Π-Eigenvalue Balance under bit_b Homogeneity

**Status:** Tier 1 derived for Hermitian H (this proof) and for non-Hermitian H (separate proof, [PROOF_F112_NONHERMITIAN_UNIVERSAL_N](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md)). Universal in N.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- F38 (Π² = (-1)^{w_YZ} on Pauli strings; `docs/ANALYTICAL_FORMULAS.md` F38 entry)
- F63 ([L, Π²] = 0 for Z-dephasing; `docs/ANALYTICAL_FORMULAS.md` F63 entry)
- F108 Part 1 ([PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md)); shared bit_b Z₂ grading
- `polarity_coordinates_from_L` primitive (`simulations/framework/diagnostics/polarity_coordinates.py`, added 2026-05-25)
- F87 dissipator-resonance law (orthogonal axis; first observed via `simulations/_polarity_probe_f87_connection.py`, derived 2026-06-10 in the dated section below)

## Abstract

When we first ran the polarity-coordinates decomposition on the F1 residual, we expected the +1/2 and −1/2 polarity sides to carry asymmetric information; that is what a non-trivial decomposition usually delivers. The probes told us something different: across every standard Lindblad system we could construct, the two sides balanced bit-exactly. F112 asks the obvious question after such a measurement: what is the exact condition under which the balance holds?

The answer turns out to be both restrictive and natural. The Hamiltonian needs to be Hermitian (the unsurprising part), and every collapse operator needs to respect a single Z₂ parity axis on the Pauli alphabet, the bit_b parity, which counts Y and Z letters modulo 2. Single-Pauli collapse operators (Z dephasing, X dephasing, depolarizing components on a single site) satisfy this trivially. The standard T1 σ⁻ operator does not, because σ⁻ has one foot in each parity sector at once.

The proof works in five steps. The first four are bookkeeping: the polarity asymmetry decomposes through the Π-eigenspace structure, the dissipator vanishes on the relevant ±i sectors by the bit_b homogeneity, and what remains is a Frobenius identity on the Hamiltonian commutator superoperator alone. The fifth step is where the work happens: for Hermitian H, the Hamiltonian commutator is anti-Hermitian as a superoperator, and the dagger map shuffles the Π +i and Π −i eigenspaces isometrically. That isometry forces equal Frobenius weights on both sides.

The diagnostic upshot is that the polarity-coordinates primitive becomes a precise witness for "this Lindblad system lives in the well-behaved channel". When the balance breaks, the input has stepped outside that channel in a specific way (mixed-bit_b collapse, or L not in standard Lindblad form), and the F113 follow-up gives the exact magnitude of the resulting asymmetry. The non-Hermitian-H case is closed separately by a different argument.

## Introduction

**The motivating question.** [POLARITY_COORDINATES.md](../../reflections/POLARITY_COORDINATES.md) introduced the `polarity_coordinates_from_L` primitive, which decomposes the F1 residual M into three orthogonal Frobenius components: `M_zero` (the F81 M_sym part) on the 0-polarity axis, and `M_plus_half` / `M_minus_half` on the +1/2 and −1/2 polarity axes. An empirical pattern emerged across a 14-probe sweep: for every standard physical Lindblad system constructed (T1, T2, depolarizing, mixed dephase letters, single-site transverse fields), `‖M_plus_half‖² = ‖M_minus_half‖²` held bit-exactly. Outside the construction channel (Probe 6: hand-engineered Π +i eigenmode L), the balance broke. The structural question: what is the exact condition that separates preserved from broken?

**The empirical anchor.** Probes 9-14 narrowed the condition to **bit_b-homogeneity of c with Hermitian H**: every collapse operator's Pauli-string support must share a single bit_b parity. Single-Pauli c (T2 = Z, single-site depolarizing components) is trivially bit_b-homogeneous; T1 σ⁻ (with support on both bit_b = 0 and bit_b = 1 components) is bit_b-mixed and sits outside the typed scope.

**What this proof closes.** Five steps, all rigorous for Hermitian H:

1. **Reduction.** Π-eigenspace decomposition of `M_plus_half` and `M_minus_half` reduces the asymmetry to `(1/2)(‖M_{+i}‖² − ‖M_{-i}‖²)`.
2. **Dissipator vanishing.** For bit_b-homogeneous c, the dissipator `np.kron(c, c.conj())` lies entirely in the Π²-conjugation +1 eigenspace (via F38 / F63).
3. **No ±i content from dissipator.** The Π²-conj +1 eigenspace equals Π-conj {+1, −1}, so the dissipator contributes zero +i / −i content.
4. **M_{±i} from L_H only.** M_{+i} and M_{-i} therefore reduce to `‖M_{±i}‖² = 2 · ‖L_{H,±i}‖²` with L_H = −i[H, ·].
5. **Anti-Hermitian closure.** For Hermitian H, L_H is anti-Hermitian as a superoperator (Lemma B), and the dagger map sends the Π +i eigenspace bijectively onto Π −i while preserving Frobenius norm (Lemma A). Hence `‖L_{H,+i}‖² = ‖L_{H,-i}‖²`.

The non-Hermitian-H case requires a separate argument because Step 5's `L_H^† = −L_H` fails when H ≠ H^†; that case is closed in [PROOF_F112_NONHERMITIAN_UNIVERSAL_N](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md).

**Diagnostic consequence.** The `polarity_coordinates_from_L` primitive is the precise witness for "L came from a standard Lindblad channel with bit_b-homogeneous collapse operators and (Hermitian or non-Hermitian) H". Asymmetry = 0 confirms the regime; asymmetry ≠ 0 localises the failure (mixed-bit_b c, or L not in Lindblad form at all). The closed-form magnitude in the canonical break case (Z-drive Hamiltonian + σ⁻ T1) is derived in [PROOF_F113_COEFFICIENT_DERIVATION](PROOF_F113_COEFFICIENT_DERIVATION.md).

## Theorem

Let L = -i[H, ·] + Σ_k γ_k · `np.kron(c_k, c_k^*)` be a Lindblad-form Liouvillian on N qubits with Hermitian H, arbitrary γ_k ∈ ℂ, and operators c_k. If each c_k has Pauli-string support entirely within a single Π²-Z parity sector (every Pauli string P in c_k's expansion satisfies bit_b(P) = (#Y(P) + #Z(P)) mod 2 = const), then the `polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I satisfies

    ‖M_plus_half‖² = ‖M_minus_half‖²

bit-exactly (machine precision), for any choice of complex coefficients in each c_k's Pauli expansion.

**Non-Hermitian extension:** the equality also holds for non-Hermitian H, Tier 1 derived universal N (Welle 11, 2026-05-27, structural proof in [PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md)). The rigorous proof of Steps 1-4 in this writeup covers both cases; Step 5 here is the Hermitian-H specialization (dagger + anti-Hermitian L_H), while the non-Hermitian closure is a separate two-lemma argument via Pauli-basis matrix-support disjointness.

## Empirical anchors

- Probes 1-5 (2026-05-25): five candidate-breakers across various standard Lindblad L (single-Pauli or low-Pauli-rank c), all preserve balance bit-exact.
- Probe 6 (2026-05-26, `_polarity_probe_pi_eigenmode.py`): hand-engineered L outside Lindblad form, balance broken with asymmetry = ‖M‖²/2 exact.
- Probes 7-8 (`_polarity_probe_random_lindblad.py`, `_polarity_probe_real_or_imag.py`): random c (full Pauli rank, NOT bit_b-homogeneous) break balance across 240 random configurations.
- Probe 9 (`_polarity_probe_pauli_rank.py`): k_max boundary search; k_pauli = 1 always preserves, k_pauli ≥ 2 selection-dependent.
- Probe 10 (`_polarity_probe_pair_enumeration.py`): at N=2, exhaustive 136 Pauli-pair enumeration with fixed coefficients (1, i) preserves balance for every pair; with random coefficients (probe 11), the structural axis emerges.
- Probe 11 (`_polarity_probe_coefficients.py`): coefficient sweep at N=2 reveals two pair classes: same Z₂³ cell preserves for all coefficients; cross-cell conditional.
- Probe 12 (`_polarity_probe_z2cubed_scaling.py`): pattern scales to N=3, N=4. Sweep A (within Z₂³ cell, random complex coefs) gives 27/27, 72/72, 72/72 BALANCED. Sweep B (cross-cell) splits exactly by bit_b parity match: same-bit_b cells preserve, cross-bit_b cells break. The Z₂³-cell structure reduces to the single Z₂ axis of bit_b.
- Probe 13 (`_polarity_proof_verify.py`, the Step 2 numerical check): bit_b-homogeneous c gives `np.kron(c, c.conj())` entirely in Π²-conj +1 eigenspace (100.00% at N=2, 3; mixed-bit_b c splits ~50/50 between Π² eigenspaces).
- Probe 14 (`_polarity_step5_stress.py`): direct Π-eigenspace projection of L_H for 30 random H configurations (10 Hermitian Pauli + 10 non-Hermitian Pauli + 10 random complex matrix) at N=2, 3. All 30 give ‖L_H,+i‖² = ‖L_H,-i‖² bit-exact.
- Probe (F87 connection, `_polarity_probe_f87_connection.py`): all three F87 trichotomy classes (truly, soft, hard) at N=3 under standard Z-deph give asymmetry = 0 bit-exact. F112 polarity balance is insensitive to F87 trichotomy; the two are orthogonal axes on the bit_b Z₂-grading.

## Proof (Hermitian H case, rigorous)

### Step 1: Reduction to Π-conjugation ±i Frobenius content equality

Decompose any operator A on Liouville space in the Π-conjugation eigenspaces: A = A_{+1} + A_{-1} + A_{+i} + A_{-i} where Π A_{λ} Π⁻¹ = λ A_{λ}. The four eigenspaces are Frobenius-orthogonal (since Π is unitary, the eigenspaces of conjugation by Π form an orthogonal decomposition).

For M, the `polarity_coordinates_from_L` definitions give:

    M_sym  := (M + Π M Π⁻¹) / 2 = M_{+1} + (1+i)/2 M_{+i} + (1-i)/2 M_{-i}
    M_anti := (M - Π M Π⁻¹) / 2 = M_{-1} + (1-i)/2 M_{+i} + (1+i)/2 M_{-i}

(verify: action of Π conjugation on each eigenmode gives the stated coefficients).

Then:

    M_plus_half  := (M_anti - i · Π M_anti Π⁻¹) / 2 = (1+i)/2 M_{-1} + (1-i)/2 M_{+i}
    M_minus_half := (M_anti + i · Π M_anti Π⁻¹) / 2 = (1-i)/2 M_{-1} + (1+i)/2 M_{-i}

Frobenius norms (using |(1±i)/2|² = 1/2 and orthogonality):

    ‖M_plus_half‖²  = (1/2) ‖M_{-1}‖² + (1/2) ‖M_{+i}‖²
    ‖M_minus_half‖² = (1/2) ‖M_{-1}‖² + (1/2) ‖M_{-i}‖²

Hence:

    asymmetry = ‖M_plus_half‖² - ‖M_minus_half‖² = (1/2) (‖M_{+i}‖² - ‖M_{-i}‖²)

**The balance condition reduces to ‖M_{+i}‖² = ‖M_{-i}‖².**

### Step 2: bit_b-homogeneous c gives dissipator entirely in Π²-conj +1 eigenspace

Per F38 / F63: Π² acts on a Pauli string σ by conjugation as Π² σ Π⁻² = (-1)^{bit_b(σ)} σ, where bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2.

(F63 anchors this for the Z-dephasing Π; F38 is the underlying Π² eigenvalue formula. F61 covers the bit_a parity and is not invoked here.)

For c = Σ_α a_α σ_α with all σ_α sharing bit_b(σ_α) = b (some fixed value in {0, 1}):

    Π² c Π⁻² = Σ_α a_α (-1)^b σ_α = (-1)^b · c

So c is a Π²-conjugation eigenvector with eigenvalue ε = (-1)^b ∈ {+1, -1}.

For the dissipator term `np.kron(c, c.conj())` (numpy kron / Liouville-superoperator convention):

    Π² · np.kron(c, c^*) · Π⁻²
        = np.kron(Π² c Π⁻², Π² c^* Π⁻²)        [Π² real on the kron factor structure]
        = np.kron(ε c, ε c^*)                   [c is Π²-eigenmode with eigenvalue ε]
        = ε² · np.kron(c, c^*)
        = +1 · np.kron(c, c^*)                  [ε² = 1 for ε ∈ {+1, -1}]

So `np.kron(c, c^*)` lies entirely in the Π²-conjugation +1 eigenspace.

**Numerical verification** (`simulations/_polarity_proof_verify.py`):

| c structure (N=2, 3) | Π²=+1 content | Π²=−1 content |
|---|---|---|
| bit_b-homogeneous (b=0 or b=1) | **100.00%** | **0.00%** |
| mixed bit_b | ~50% | ~50% (depending on coefs) |

Bit-exact to machine precision (1e-15) across all tested configurations.

### Step 3: Π²-conj +1 eigenspace = Π-conj {+1, -1} eigenspaces

For Π unitary order-4 with eigenvalues {+1, -1, +i, -i}, the squares are {+1, +1, -1, -1}. Π-conjugation eigenvalues +1 and -1 square to +1; Π-conjugation eigenvalues +i and -i square to -1.

So Π-conjugation +1 and -1 eigenspaces are contained in Π²-conjugation +1 eigenspace. Π-conjugation +i and -i eigenspaces are contained in Π²-conjugation -1 eigenspace.

By Step 2, `np.kron(c, c^*)` has zero Π²-conj -1 content, hence zero Π-conj +i and zero Π-conj -i content.

**The dissipator part of M does NOT contribute to M_{+i} or M_{-i}.**

### Step 4: M_{+i} and M_{-i} come entirely from the Hamiltonian part L_H

M = Π L Π⁻¹ + L + 2σ·I has four contributions:
- 2σ·I: Π·I·Π⁻¹ = I (identity is Π-conjugation-fixed), so 2σ·I is in Π-conj +1 eigenspace. No +i or -i content.
- L_H = -i[H, ·]: contributes to all four eigenspaces in general.
- Σ_k γ_k `np.kron(c_k, c_k^*)`: per Step 3, no +i or -i content (when each c_k is bit_b-homogeneous).
- Π·L·Π⁻¹: same structural decomposition under Π conjugation; the Π-fixed identity and Π²-conj +1 dissipator pieces stay in {+1, -1}.

Hence M_{+i} = (Π L_H Π⁻¹ + L_H)_{+i} and M_{-i} = (Π L_H Π⁻¹ + L_H)_{-i}.

Decompose L_H = L_{H,+1} + L_{H,-1} + L_{H,+i} + L_{H,-i} in Π-conjugation eigenspaces. Then:

    (L_H + Π L_H Π⁻¹)_{+i} = L_{H,+i} + i · L_{H,+i} = (1+i) L_{H,+i}
    (L_H + Π L_H Π⁻¹)_{-i} = L_{H,-i} + (-i) · L_{H,-i} = (1-i) L_{H,-i}

So M_{+i} = (1+i) L_{H,+i} and M_{-i} = (1-i) L_{H,-i}, giving:

    ‖M_{+i}‖² = |1+i|² · ‖L_{H,+i}‖² = 2 · ‖L_{H,+i}‖²
    ‖M_{-i}‖² = |1-i|² · ‖L_{H,-i}‖² = 2 · ‖L_{H,-i}‖²

So **‖M_{+i}‖² = ‖M_{-i}‖² ⟺ ‖L_{H,+i}‖² = ‖L_{H,-i}‖²**.

### Step 5: For Hermitian H, ‖L_{H,+i}‖² = ‖L_{H,-i}‖² via dagger / anti-Hermitian argument

**Lemma A (dagger maps Π +i ↔ Π -i):** For Π unitary and any superoperator A:

    (Π A Π⁻¹)^† = (Π⁻¹)^† · A^† · Π^† = Π · A^† · Π⁻¹    [using Π^† = Π⁻¹]

So if A has Π-conjugation eigenvalue λ, then A^† has Π-conjugation eigenvalue λ^* (complex conjugate). Specifically, the Hilbert-Schmidt dagger A → A^† is an antilinear isometry mapping the Π +i eigenspace bijectively to the Π -i eigenspace.

Since dagger preserves Frobenius norm (‖A^†‖² = ‖A‖²) and the eigenspace decomposition is orthogonal, the components satisfy:

    ‖A_{+i}‖² = ‖(A^†)_{-i}‖²

**Lemma B (L_H for Hermitian H is anti-Hermitian as superoperator):** For Hermitian H, the commutator superoperator L_H = -i[H, ·] satisfies L_H^† = -L_H.

Proof: For any operators X, Y, the Hilbert-Schmidt inner product gives
    ⟨L_H(X), Y⟩ = Tr((L_H(X))^† Y) = Tr((-i[H,X])^† Y) = i Tr([X^†, H^†] Y) = i Tr([X^†, H] Y)
                                                                                     [H Hermitian ⟹ H^† = H]
                = i Tr(X^† H Y - H X^† Y) = i Tr(X^† H Y) - i Tr(X^† Y H)              [cyclic]
                = i Tr(X^† (H Y - Y H)) = i Tr(X^† [H, Y])
                = -Tr(X^† (-i[H, Y])) = -⟨X, L_H(Y)⟩

So L_H^† acts as -L_H, i.e., L_H^† = -L_H. ∎

**Combining Lemma A and Lemma B for Hermitian H:**

    ‖L_{H,+i}‖² = ‖(L_H^†)_{-i}‖²        [Lemma A, applied to L_H]
                = ‖(-L_H)_{-i}‖²           [Lemma B]
                = ‖L_{H,-i}‖²              [norm of -X equals norm of X]

So **‖L_{H,+i}‖² = ‖L_{H,-i}‖² for any Hermitian H**, completing Step 5.

Combining Steps 1-5: the balance condition ‖M_plus_half‖² = ‖M_minus_half‖² holds for any Lindblad-form L with Hermitian H and bit_b-homogeneous c_k. ∎

**Empirical verification of Step 5** (probe 14, `_polarity_step5_stress.py`): direct Π-eigenspace projection at N=2, 3 across 10 random Hermitian H gives ‖L_{H,+i}‖² = ‖L_{H,-i}‖² bit-exact (relative difference < 1e-15) in all 10 cases.

## Non-Hermitian H extension: Tier1Derived universal N (2026-05-27, Welle 11)

The non-Hermitian extension is now closed structurally for all N in [PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md).

Writing H = H_re + i H_im (Hermitian decomposition with H_re = (H + H^†)/2, H_im = (H − H^†)/(2i) both Hermitian), the equality reduces algebraically to the identity:

    F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0   for any Hermitian H_re, H_im.

Bilinearity (L_H linear in H, Π-projection linear, Frobenius sesquilinear) plus antisymmetry (under H_re ↔ H_im exchange) reduce F ≡ 0 to F(σ_α, σ_β) = 0 on every Pauli-string pair. The per-pair identity is closed structurally via two lemmas:

- **Lemma N-A (Diagonal-Norm)** (Welle 11 non-Hermitian-extension lemma, distinct from the parent Hermitian-H Lemma A above). For any BitB-odd Pauli string σ at length N, ‖L_{σ,−i}‖² = 4^N. Combines support-counting ‖L_σ‖² = 2 · 4^N + the F38 / F63 Π²-conjugation eigenvalue + the cross-term vanishing ⟨L_σ, Π L_σ Π⁻¹⟩ = 0 via Pauli-basis matrix-support disjointness (M(L_σ) and Π M(L_σ) Π⁻¹ have complementary non-zero entries on the same shifted-diagonal positions).
- **Lemma N-B (Off-Diagonal-Orthogonality)** (Welle 11 non-Hermitian-extension lemma, distinct from the parent Hermitian-H Lemma B above). For σ_α ≠ σ_β both BitB-odd, ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0. Reduces to ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩ = 0 for all m ∈ {0, 1, 2, 3}; each established by matrix-support disjointness (overlap requires σ_α = σ_β).

Both Welle 11 lemmas (N-A and N-B) reduce to per-position checks on the 4^N × 4^N matrix of L_σ that are uniform in N. The proof is N-independent.

**Welle 10 numerical anchor** (preserved as historical empirical validation): the per-pair identity F(σ_α, σ_β) = 0 was verified bit-exact across 559,912 distinct upper-triangular pairs at N = 2, 3, 4, 5 by `simulations/_f112_open_identity_basis_enum.py` (Python, Welle 10a) and `compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs` (C#, Welle 10b). All N = 2..4 bit-exact 0; N = 5 < 1e-10. See [F112_NONHERMITIAN_BASIS_ENUMERATION.md](../../experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md).

**Welle 11 structural verifier** `simulations/_f112_universal_n_proof_verify.py` confirms each step of Lemmas A and B bit-exact at N = 1, 2, 3 (42 BitB-odd strings, 1050 off-diagonal pairs, 4368 all-pair F-values, all 0.000e+00 in numpy double precision).

## The F87 orthogonality, derived (2026-06-10)

The F87 connection probe (last entry in the empirical anchor list) left one line open since 2026-05-26: all three F87 trichotomy classes give F112 asymmetry = 0, observed at N=3 but not derived. The derivation turns out to need no new machinery. It has three parts, all verified in `simulations/f112_f87_orthogonality.py` (committed verifier; supersedes the hunt drafts `_f112_f87_orthogonality_check.py` and `_f112_f87_orth_supplement.py`).

**(a) Scope inclusion.** Every F87 input is an instance of this proof's hypotheses. F87's domain is the commutator of a Hermitian Pauli Hamiltonian plus pure Z-dephasing; each collapse operator Z_l is a single Pauli string, hence trivially bit_b-homogeneous. So the theorem applies verbatim, and the F112 asymmetry is identically zero on F87's entire domain, across all three trichotomy classes (truly, soft, hard). The orthogonality is not a coincidence the probe happened to sample; it is the statement that the two functionals never co-vary because one is identically zero where the other lives. Verified: asymmetry = 0.0 exactly (float zero, not merely small) for truly Heisenberg, soft XXZ+ZXX, the hard pairs XXZ+ZZZ (m\*=5), XXZ+XZX and IXY+XIY (m\*=9), and the Z-drive XY chain, at N = 3 and N = 4.

**(b) Mechanism separation: same involution, different degree.** For H supported on bit_b-odd Pauli strings, X^⊗N-conjugation flips the sign of every term (F H F = −H with F = X^⊗N; this is the driving lemma of [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2, and it covers the diagonal Klein cell where all F87 pair-hardness lives, F110). There the dagger involution that powers Step 5 coincides with the windowed converse's first reflection, on the windowed converse's recentred object M_rec = L + σ·I:

    M_rec† = 𝓕 M_rec 𝓕,   𝓕 = F ⊗ F = X^⊗N ⊗ X^⊗N.

Derivation in one line: 𝓕 L_H 𝓕 = L_{FHF} = L_{−H} = −L_H = L_H† (Lemma B), while the Z-dephasing dissipator and the σ·I shift are Hermitian and 𝓕-invariant (F Z_l F = −Z_l, two sign flips cancelling across the two tensor legs). Verified diff 0.00e+00 for the diagonal-cell pairs XXZ+ZXX and IXY+XIY at N=3. The dagger map that 𝓕 implements here is the same L_H† = −L_H that Step 5 applies to the F1 residual. The separation: F112's functional reads exactly the degree-2 content of this shared involution, Frobenius norms of Π-eigenprojections (Steps 1-5 use nothing else). The F87 hardness decision lives at odd degree: the recentred power sums p₃, p₅, … under the second reflection R = I ⊗ F and the unsigned girth. Degree-2 norms are blind to odd-degree sign structure, so the same bit_b grading carries both functionals without coupling them. The verifier shows hard pairs firing at their true m\* (p₉ for XXZ+XZX and IXY+XIY at N=4) while the F112 asymmetry stays exactly 0.0.

**(c) Scoped one-way implication on the σ⁻/σ⁺ probe family.** Off the shared domain, with amplitude damping added (the F113 regime), the two functionals finally meet, and they meet on a single shared moment:

    t₁^(l) = Tr(Z_l H) = 2^N c_l = 2^(N−1) ω_l,

where c_l is the Z_l coefficient of H and ω_l the documented F113 drive convention (H ⊇ Σ_l (ω_l/2) Z_l). The F113 closed form reads the moment linearly, paired with the net rate:

    asymmetry = 2^N · Σ_l t₁^(l) · (γ_pump,l − γ_T1,l)   [= (4^N/2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)],

while the girth ladder's ℓ = 1 deg-1 face squares it:

    p₃ = 6γ · Σ_l (t₁^(l))²   [= 6 · 4^N · Σ_l c_l² · γ].

Hence: F112 balance broken under the σ⁻/σ⁺ family (generic rates) ⟹ some t₁^(l) ≠ 0 ⟹ p₃ > 0 ⟹ F87-hard. Strictly one-way: hard does not imply broken (the flux and K3 hard pairs have all t₁^(l) = 0 and keep the balance under T1, asymmetry = 0.0 in the verifier). Reconciled numbers at N=3, c = (1, 2, 1), γ_T1 = 0.03, γ_pump = 0: predicted −7.680000000000000 vs measured −7.680000000000007, |diff| = 7.1e-15, machine precision (not exactly equal as floats, so not "bit-exact"). The hunt draft printed pred = −asym/2; the missing chain was ω_l = 2c_l plus the standard-physics σ⁻-lowering sign convention (reconciled per the F113 entry and `LindbladBitBPiBreakMagnitude.cs`).

**Scope note.** The verifier runs at N = 3 and N = 4. The trichotomy labels of the K3 and flux pairs are N-dependent (N=3-soft; hard from N=4 where m\* = 9 first fires), which is why the original probe at N=3 could only sample soft representatives of those pairs. Part (a) is N-independent (the theorem's hypotheses do not see N); part (b)'s identity is derived for any N and verified at N=3; part (c) combines the F113 universal-N closed form with the girth ladder's cell-free m = 3 face.

## Status

**F112 is Tier 1 derived for both Hermitian and non-Hermitian H, universally in N.**

Hermitian H (this proof, original Welle): five-step proof.

- Step 1: reduction to ‖M_{+i}‖² = ‖M_{-i}‖² (Π-eigenspace decomposition of M_plus_half / M_minus_half).
- Step 2: bit_b-homogeneous c implies `np.kron(c, c.conj())` is Π²-conj +1 (via F38/F63 Π² eigenvalue formula on Pauli strings).
- Step 3: Π²-conj +1 eigenspace = Π-conj {+1, −1}, hence no +i or −i content.
- Step 4: M_{+i} and M_{-i} come entirely from L_H, with norms 2 · ‖L_{H,±i}‖².
- Step 5: For Hermitian H, L_H is anti-Hermitian as superoperator (Lemma B), and dagger maps Π +i ↔ Π −i bijectively while preserving Frobenius (Lemma A). Combining gives ‖L_{H,+i}‖² = ‖L_{H,-i}‖².

Non-Hermitian H (Welle 11, 2026-05-27): two-lemma structural proof in [PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md). Reduces the equality to F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0, then closes the per-pair identity F(σ_α, σ_β) = 0 structurally via Pauli-basis matrix-support disjointness (Lemmas A and B). Welle 10 numerical anchor at N ≤ 5 (559,912 pairs all bit-exact 0) preserved as historical empirical validation.

The theorem is typed as a single Tier1Derived claim in C# Core (`LindbladBitBPiBalance`) covering both H cases universally in N.

## Significance

If formalized as Tier1Derived, F112 becomes the structural identity behind the polarity_coordinates_from_L diagnostic:

- For ANY standard Lindblad system (Hermitian H + single-Pauli c, the physically meaningful case), the diagnostic asymmetry is 0 by Steps 1-5.
- Asymmetry ≠ 0 is the precise witness for L outside the bit_b-homogeneous-c regime (either non-standard c, or L not in Lindblad form at all).

Connections:
- **F38**: Π² = (-1)^{bit_b} on Pauli strings (foundational input).
- **F63**: [L, Π²] = 0 for Z-deph (foundational input via Π²-eigenvalue commutation).
- **F108 Part 1/2/3**: the bilinear set {XX, YY, YZ, ZY, ZZ} that F108 palindromizes is exactly the bit_b=0 (Π²-Z-even) family. F108's closure mechanism and F112's balance mechanism are both consequences of the bit_b Z₂ grading on the Pauli group.
- **F87 dissipator-resonance law**: orthogonal axis, **derived 2026-06-10** (previously empirical via `_polarity_probe_f87_connection.py`). F87 lives in M's spectrum-palindrome structure; F112 lives in M_anti's Π +i/-i split. The derivation (scope inclusion + mechanism separation + the scoped F113 one-way bridge) is in [The F87 orthogonality, derived](#the-f87-orthogonality-derived-2026-06-10); committed verifier `simulations/f112_f87_orthogonality.py`.
- **`polarity_coordinates_from_L`**: F112 makes the primitive's diagnostic value precise. Asymmetry ≠ 0 detects c with cross-bit_b Pauli support, which is OUTSIDE the F108-closure regime.

## Open

- ~~**Step 5 extension to non-Hermitian H**~~: **CLOSED 2026-05-27 (Welle 11)** in [PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md) via the two-lemma structural proof (Lemma N-A: ‖L_{σ,−i}‖² = 4^N for BitB-odd σ; Lemma N-B: ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0 for σ_α ≠ σ_β both BitB-odd; both via Pauli-basis matrix-support disjointness). F112 non-Hermitian extension is now Tier1Derived universal N.
- ~~**F87 ↔ F112 orthogonality**~~: **CLOSED 2026-06-10** in [The F87 orthogonality, derived](#the-f87-orthogonality-derived-2026-06-10) above: (a) scope inclusion (F87's entire domain satisfies this proof's hypotheses, so the F112 asymmetry is identically zero where F87 lives), (b) mechanism separation (shared dagger involution M_rec† = 𝓕 M_rec 𝓕 on the diagonal cell; F112 reads it at degree 2, F87 hardness at odd degree), (c) the scoped one-way F113 bridge (balance-broken ⟹ F87-hard via the shared moment t₁^(l) = Tr(Z_l H)). Committed verifier `simulations/f112_f87_orthogonality.py`.
- **C# Core typing**: Tier1Derived for both Hermitian and non-Hermitian H (universal N).
- **Connection to F104, F105, F106** (F87 Z₂³-cubed refinements at various N, k): potential bridge to the F112 balance via shared bit_b structure.
