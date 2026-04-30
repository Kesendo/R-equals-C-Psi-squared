# Proof of F81: Π-Conjugation of M Decomposes into the Π²-Odd Hamiltonian Commutator

**Tier:** 1 (closed-form algebraic proof + numerical verification at machine precision).
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_SVD_CLUSTER_STRUCTURE.md](PROOF_SVD_CLUSTER_STRUCTURE.md) (Master Lemma: M is γ-independent for pure Z-dephasing; Π² acts on Pauli string σ_α as (-1)^{bit_b(α)})
- [`framework/symmetry.py`](../../simulations/framework/symmetry.py) (`build_pi_full`, `pi_squared_eigenvalue`)
- [`framework/lindblad.py`](../../simulations/framework/lindblad.py) (`palindrome_residual`, `lindbladian_z_dephasing`)

**Statement (Theorem F81):** For any 2-bilinear Hamiltonian H = H_even + H_odd (Π²-decomposition into bit_b-even and bit_b-odd Pauli bilinears) under uniform Z-dephasing,

    Π · M · Π⁻¹ = M − 2 · L_{H_odd}

where L_{H_odd} = -i[H_odd, ·] is the unitary commutator induced by the Π²-odd part of H. Equivalently the Π-conjugation symmetric/antisymmetric decomposition of M reads

    M_sym  = (M + Π·M·Π⁻¹) / 2 = Π·L·Π⁻¹ + L_diss + L_{H_even} + 2Σγ·I
    M_anti = (M − Π·M·Π⁻¹) / 2 = L_{H_odd}

with M_sym ⊥_F M_anti (Frobenius-orthogonal) and ‖M‖²_F = ‖M_sym‖²_F + ‖M_anti‖²_F.

---

## Numerical verification (N=3, γ_Z=0.1, Σγ=0.3, all residuals at machine precision 1e-16)

| Hamiltonian | trichotomy | H_odd | ‖Π·M·Π⁻¹ − (M − 2·L_{H_odd})‖ | ‖M_anti − L_{H_odd}‖ |
|-------------|------------|-------|--------------------------------|----------------------|
| XX+YY | truly | 0 | 1.5e-16 (M=0 trivially) | 0 |
| YZ+ZY | soft (Π²-even non-truly) | 0 | 1.5e-16 | 0 |
| XY+YX | soft (Π²-odd) | XY+YX | 1.5e-16 | 7.6e-17 |
| XX+XY | hard (mixed) | XY only | 1.5e-16 | 7.6e-17 |
| pure XY | (Π²-odd) | XY | 1.5e-16 | (analogous) |
| pure XZ | (Π²-odd) | XZ | 1.5e-16 | (analogous) |

For pure Π²-odd 2-body chain Hamiltonians at any N and any γ, the Π-decomposition has ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 exactly (50/50 split). Verified numerically at N=3, 4, 5 with γ_Z ∈ {0, 0.05, 0.1, 0.5, 1.0}; analytical reason given in Step 8 below. The hard mixed case (XX+XY) at N=3 also gives 50/50 by a separate numerical coincidence specific to that combination, not a general structural theorem.

---

## Proof

### Step 1: Π² acts on Pauli strings as a Π²-parity sign

For a Pauli string σ_α with bit_b-encoding (a₁b₁)(a₂b₂)...(aNbN), the framework's Π² operator is diagonal in the Pauli-string basis with eigenvalue

    Π² σ_α = (-1)^{bit_b(α)} σ_α       where bit_b(α) = (b₁ + b₂ + ... + bN) mod 2.

In other words, in the Pauli basis Π² is the diagonal operator D = diag(±1) with sign on each basis vector determined by its total bit_b-parity. This is a standard property of the framework's Π construction (see [`pi_squared_eigenvalue`](../../simulations/framework/symmetry.py)).

### Step 2: Π² conjugation acts on the Liouville superoperator by Pauli-product sign

Since Π² is diagonal in the Pauli basis, conjugation of a 4^N × 4^N matrix M by Π² acts entry-wise:

    (Π² · M · Π⁻²)_{γβ} = D_γγ · M_{γβ} · D⁻¹_ββ = (-1)^{bit_b(γ)+bit_b(β)} · M_{γβ}.

The matrix element M_{γβ} is non-zero only when M takes |σ_β⟩ → |σ_γ⟩ under its action.

### Step 3: For L_H = -i[H, ·], the matrix-element sign factor reduces to bit_b(α)

For H = Σ_α h_α σ_α and a Pauli basis vector |σ_β⟩, compute

    L_H |σ_β⟩ = -i Σ_α h_α |[σ_α, σ_β]⟩.

The commutator of two Pauli strings is either zero (if they commute) or proportional to ±2i · σ_{α·β} (if they anticommute), where α·β denotes the Pauli-product index. Using the identity

    bit_b(α · β) = bit_b(α) + bit_b(β) (mod 2)

(Pauli multiplication is bit-wise XOR on (a, b) indices), the sign factor in Step 2 simplifies:

    (-1)^{bit_b(γ)+bit_b(β)} = (-1)^{bit_b(α·β)+bit_b(β)} = (-1)^{bit_b(α)}.

Hence each term -i · h_α · [σ_α, ·] of L_H gets multiplied by (-1)^{bit_b(α)} under Π² conjugation. The Hamiltonian matrix elements with σ_α Π²-even (bit_b(α) = 0) are preserved; those with σ_α Π²-odd (bit_b(α) = 1) flip sign. Splitting H = H_even + H_odd by bit_b-parity:

    Π² · L_H · Π⁻² = L_{H_even} − L_{H_odd}.

### Step 4: Z-dephasing dissipator is Π²-symmetric

For the Z-dephasing dissipator L_diss = Σ_l γ_l (Z_l · ρ · Z_l − ρ), the channel ρ → Z_l ρ Z_l is diagonal in the Pauli basis: (Z_l σ_β Z_l) = (-1)^{bit_a(β at site l)} · σ_β. Since both L_diss and Π² are diagonal in this basis, they commute: [Π², L_diss] = 0. Therefore

    Π² · L_diss · Π⁻² = L_diss.

(For non-Z dissipators, this step changes; F81 in its current form is Z-dephasing-specific.)

### Step 5: Combining to get the F81 identity

Splitting L = L_H + L_diss = L_{H_even} + L_{H_odd} + L_diss, Steps 3 and 4 give

    Π² · L · Π⁻² = L_{H_even} − L_{H_odd} + L_diss = L − 2 · L_{H_odd}.

Now apply Π conjugation to the palindrome equation M = Π·L·Π⁻¹ + L + 2Σγ·I:

    Π · M · Π⁻¹ = Π² · L · Π⁻² + Π · L · Π⁻¹ + 2Σγ·I
                = (L − 2·L_{H_odd}) + Π·L·Π⁻¹ + 2Σγ·I
                = (Π·L·Π⁻¹ + L + 2Σγ·I) − 2·L_{H_odd}
                = M − 2 · L_{H_odd}.    ∎

### Step 6: Π-symmetric/antisymmetric decomposition

Define M_sym = (M + Π·M·Π⁻¹)/2 and M_anti = (M − Π·M·Π⁻¹)/2. From Step 5:

    M_sym  = (M + M − 2·L_{H_odd}) / 2 = M − L_{H_odd}
           = Π·L·Π⁻¹ + L_{H_even} + L_diss + 2Σγ·I,
    M_anti = (M − M + 2·L_{H_odd}) / 2 = L_{H_odd}.

The Frobenius orthogonality M_sym ⊥_F M_anti follows because they live in disjoint Π²-eigenspaces (M_sym is in the +1 eigenspace of the Π-conjugation operation on operators; M_anti is in the −1 eigenspace), and ‖M‖²_F = ‖M_sym‖²_F + ‖M_anti‖²_F by Pythagoras.

### Step 7: Spectral consequence

Spec(Π·M·Π⁻¹) = Spec(M) by unitary invariance of the spectrum (Π is unitary). Combined with Step 5:

    Spec(M) = Spec(M − 2·L_{H_odd}).

Subtracting 2·L_{H_odd} from M yields a similar matrix (related by Π-conjugation). Combined with F80's structural identity Spec(M) = ±2i · Spec_{many-body}(H_non-truly), this constrains how L_{H_odd}'s spectrum interacts with M's: Spec(L_{H_odd}) need not equal Spec(M) in general (and does not, since M has eigenvalues outside L_{H_odd}'s range when H_even or L_diss is nonzero), but the *similarity* M ~ M − 2·L_{H_odd} via Π-conjugation is exact.

### Step 8: 50/50 split when all non-truly bilinears are Π²-odd

The 50/50 split ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 holds for any 2-body chain Hamiltonian whose **non-truly bilinears are all Π²-odd**. This includes both pure Π²-odd H (H = H_odd) and mixed truly + Π²-odd H (H = H_truly + H_odd, e.g. XX+XY). The condition fails when H contains Π²-even non-truly bilinears (YZ-type), where M_anti shrinks relative to M_sym.

Derivation: Frobenius orthogonality from Step 6 gives ‖M‖² = ‖M_sym‖² + ‖M_anti‖² always. M_anti = L_{H_odd} = -i[H_odd, ·] by F81. So 50/50 ⟺ ‖M_anti‖² = ‖M‖²/2 ⟺ ‖M‖² = 2·‖L_{H_odd}‖².

Two ingredients:

(i) **Master Lemma for the residual:** truly bilinears contribute zero to M (PROOF_SVD_CLUSTER_STRUCTURE.md). Therefore ‖M‖² depends only on the non-truly part of H. If all non-truly bilinears are Π²-odd, then "non-truly" coincides with "Π²-odd," and the Frobenius residual scaling (F49 chain version) reads:

    ‖M‖²_F = 4 · ‖H_odd‖²_F · 2^N.

(ii) **Standard commutator-Frobenius identity:** for traceless Hermitian H_odd acting on a 2^N-dimensional Hilbert space,

    ‖L_{H_odd}‖²_F = 2 · 2^N · ‖H_odd‖²_F.

(Standard: ‖[H, ·]‖²_F = 2·d·‖H‖²_F − 2·|tr(H)|², with tr(H_odd) = 0 since H_odd is a sum of non-identity Pauli strings.)

Substituting: 2·‖L_{H_odd}‖² / ‖M‖² = (2 · 2 · 2^N · ‖H_odd‖²) / (4 · ‖H_odd‖² · 2^N) = 1. Hence ‖M_anti‖² = ‖M‖² / 2 exactly, independent of N, γ, or the relative coefficients of truly and Π²-odd bilinears in H.

**Empirically verified**: pure Π²-odd (XY+YX, XY, XZ), mixed truly + Π²-odd (XX+XY, YY+XY), at N=3, 4, 5 with γ_Z ∈ {0, 0.05, 0.1, 0.5, 1.0}, all to machine precision. Cases that violate the condition (Π²-even non-truly content) give different splits:

  - Pure Π²-even non-truly (YZ alone, YZ+ZY): 100/0 (M_anti = 0 trivially since H_odd = 0).
  - Mixed Π²-odd + Π²-even non-truly (XY+YZ, XX+XY+YZ): 5/6 sym + 1/6 anti at N=3 and N=4 (the 1/6 reflects the smaller fraction of Π²-odd content within H_non-truly).

The 5/6:1/6 split for mixed odd+even cases is empirical at N=3,4; an analytical derivation analogous to the 50/50 case would require a generalized Frobenius-residual scaling for mixed H_non-truly. Open.

**Scope of 50/50**: this analytical result covers truly H (trivial: M=0), pure Π²-odd 2-body chain H, and mixed truly + Π²-odd 2-body chain H. It does not cover Hamiltonians with Π²-even non-truly content.

---

## Reading

F81 is the algebraic backbone of the geometric reading in [reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md](../../reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md). It states precisely how M and its Π-conjugate relate as operators when H has Π²-odd content. The "through-line" of the framework is therefore:

- For Π²-even H: M itself (Π·M·Π⁻¹ = M).
- For Π²-odd H: the spectrum Spec(M), with M and Π·M·Π⁻¹ split by the explicit shift −2·L_{H_odd}.

Companion to F80: F80 states what Spec(M) is (= ±2i · Spec(H_non-truly)). F81 states how M and Π·M·Π⁻¹ relate as operators sharing that spectrum. Together they characterize the Π-action on M completely.

**Verified:** N=3 all listed cases at machine precision (1e-16 residuals on Π·M·Π⁻¹ vs M − 2·L_{H_odd} and on M_anti vs L_{H_odd}). Pytest lock pending in next commit.

**Open generalizations:**
- Non-Z dissipators (T1 amplitude damping, X/Y dephasing): Π² no longer commutes with L_diss. F81 should generalize with an additional dissipator-correction term; analytical form not yet worked out.
- Higher-body Hamiltonians (3-body, 4-body): the Π²-parity classification carries over; the proof structure should generalize verbatim.
