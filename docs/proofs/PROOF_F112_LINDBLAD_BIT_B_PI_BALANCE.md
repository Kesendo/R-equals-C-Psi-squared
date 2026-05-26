# PROOF F112: Lindblad Π-Eigenvalue Balance under bit_b Homogeneity

**Status:** Tier 1 derived for Hermitian H (rigorous proof via dagger + anti-Hermitian L_H). Tier 1 candidate for non-Hermitian H extension (empirical bit-exact across 20+ random configs; structural proof open).
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- F38 (Π² = (-1)^{w_YZ} on Pauli strings; `docs/ANALYTICAL_FORMULAS.md` F38 entry)
- F63 ([L, Π²] = 0 for Z-dephasing; `docs/ANALYTICAL_FORMULAS.md` F63 entry)
- F108 Part 1 ([PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md)) — shared bit_b Z₂ grading
- `polarity_coordinates_from_L` primitive (`simulations/framework/diagnostics/polarity_coordinates.py`, added 2026-05-25)
- F87 dissipator-resonance law (orthogonal axis; empirically established via `simulations/_polarity_probe_f87_connection.py`)

## Theorem

Let L = -i[H, ·] + Σ_k γ_k · `np.kron(c_k, c_k^*)` be a Lindblad-form Liouvillian on N qubits with Hermitian H, arbitrary γ_k ∈ ℂ, and operators c_k. If each c_k has Pauli-string support entirely within a single Π²-Z parity sector (every Pauli string P in c_k's expansion satisfies bit_b(P) = (#Y(P) + #Z(P)) mod 2 = const), then the `polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I satisfies

    ‖M_plus_half‖² = ‖M_minus_half‖²

bit-exactly (machine precision), for any choice of complex coefficients in each c_k's Pauli expansion.

**Empirical extension:** the equality is also observed bit-exact for non-Hermitian H across 20 random configurations at N=2, N=3 (probe 14 `_polarity_step5_stress.py` Tests 2-3). The structural proof for non-Hermitian H is open; the rigorous proof in Step 5 below covers Hermitian H only.

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

## Non-Hermitian H extension (empirical, structural proof open)

Probe 14 Tests 2 and 3 (non-Hermitian Pauli sums + random complex matrix H, 20 configurations at N=2, 3): all 20 also give ‖L_{H,+i}‖² = ‖L_{H,-i}‖² bit-exact.

The structural reason is non-obvious. Writing H = H_re + i H_im (Hermitian decomposition with H_re = (H + H^†)/2, H_im = (H − H^†)/(2i) both Hermitian), the equality reduces to a specific identity:

    Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0

This identity holds bit-exact empirically for arbitrary Hermitian H_re, H_im. Its rigorous proof is open; candidate routes include an explicit Π-eigenbasis construction or a refined dagger-based argument that exploits both H_re's and H_im's anti-Hermitian L_H structure jointly.

For the typed F112 claim, the Hermitian-H scope (covering all standard Lindblad systems) is sufficient. The non-Hermitian extension is documented as empirical with bit-exact anchor.

## Status

**F112 is rigorously proven for the Hermitian-H case** (the physically relevant scope for standard Lindblad systems). All five steps in the Hermitian case are rigorous:

- Step 1: reduction to ‖M_{+i}‖² = ‖M_{-i}‖² (Π-eigenspace decomposition of M_plus_half / M_minus_half).
- Step 2: bit_b-homogeneous c implies `np.kron(c, c.conj())` is Π²-conj +1 (via F38/F63 Π² eigenvalue formula on Pauli strings).
- Step 3: Π²-conj +1 eigenspace = Π-conj {+1, −1}, hence no +i or −i content.
- Step 4: M_{+i} and M_{-i} come entirely from L_H, with norms 2 · ‖L_{H,±i}‖².
- Step 5: For Hermitian H, L_H is anti-Hermitian as superoperator (Lemma B), and dagger maps Π +i ↔ Π −i bijectively while preserving Frobenius (Lemma A). Combining gives ‖L_{H,+i}‖² = ‖L_{H,-i}‖².

**Non-Hermitian H extension is empirical** (bit-exact across 20 random configurations); the structural proof reduces to a specific identity Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im, whose rigorous derivation is open.

The theorem now stands ready for typing as a Tier1Derived claim in C# Core (`StandardLindbladBitBPiBalance` or similar F112) restricted to Hermitian H; the non-Hermitian extension can be typed as a separate Tier1Candidate corollary.

## Significance

If formalized as Tier1Derived, F112 becomes the structural identity behind the polarity_coordinates_from_L diagnostic:

- For ANY standard Lindblad system (Hermitian H + single-Pauli c, the physically meaningful case), the diagnostic asymmetry is 0 by Steps 1-5.
- Asymmetry ≠ 0 is the precise witness for L outside the bit_b-homogeneous-c regime (either non-standard c, or L not in Lindblad form at all).

Connections:
- **F38**: Π² = (-1)^{bit_b} on Pauli strings (foundational input).
- **F63**: [L, Π²] = 0 for Z-deph (foundational input via Π²-eigenvalue commutation).
- **F108 Part 1/2/3**: the bilinear set {XX, YY, YZ, ZY, ZZ} that F108 palindromizes is exactly the bit_b=0 (Π²-Z-even) family. F108's closure mechanism and F112's balance mechanism are both consequences of the bit_b Z₂ grading on the Pauli group.
- **F87 dissipator-resonance law**: empirically established as orthogonal axis via `_polarity_probe_f87_connection.py` (all three F87 trichotomy classes give F112 balance = 0 bit-exact at N=3). F87 lives in M's spectrum-palindrome structure; F112 lives in M_anti's Π +i/-i split. Both projections of the same bit_b Z₂-grading.
- **`polarity_coordinates_from_L`**: F112 makes the primitive's diagnostic value precise. Asymmetry ≠ 0 detects c with cross-bit_b Pauli support, which is OUTSIDE the F108-closure regime.

## Open

- **Step 5 extension to non-Hermitian H**: reduces to proving Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im. Candidate routes: explicit Π-eigenbasis construction, or refined dagger argument using both H_re and H_im anti-Hermitian L_H structure jointly.
- **F87 ↔ F112 orthogonality**: empirically confirmed via the F87-connection probe; structural derivation that both axes are projections of the bit_b Z₂-grading on the Pauli group is deferred.
- **C# Core typing**: Tier1Derived for Hermitian H; Tier1Candidate corollary for non-Hermitian H extension.
- **Connection to F104, F105, F106** (F87 Z₂³-cubed refinements at various N, k): potential bridge to the F112 balance via shared bit_b structure.
