# The n_XY Parity Selection Rule

## Statement

**Theorem.** The Liouvillian of the isotropic Heisenberg model under
Z-dephasing preserves n_XY parity. As a consequence:

1. Every Liouvillian eigenmode has definite n_XY parity (even or odd).
2. Every single-excitation density matrix has purely even n_XY content.
3. No single-excitation state can excite an odd-n_XY eigenmode.

**Corollary (Accessibility Boundary).** If the slowest Liouvillian
eigenmode accessible to single-excitation states has rate α₁, and
a slower mode exists at rate α₂ < α₁ with odd n_XY parity,
then α₁ is an exact ceiling for single-excitation protection.
No optimizer, no ansatz, no trick within the SE family can reach α₂.

---

## Definitions

Let N be the number of qubits. The Pauli basis for the operator space
consists of all N-fold tensor products of {I, X, Y, Z}.

For a Pauli string P = σ₁ ⊗ σ₂ ⊗ … ⊗ σ_N, define:

    n_XY(P) = number of indices k where σ_k ∈ {X, Y}

The **n_XY parity** of P is n_XY(P) mod 2.

The **single-excitation (SE) sector** is spanned by the N basis states
|e_k⟩ = |0…1_k…0⟩ for k = 0, …, N−1 (exactly one qubit excited).

An **SE operator** is any operator of the form ρ = Σ_{j,k} c_{jk} |e_j⟩⟨e_k|,
i.e., any operator supported on the SE × SE block of the full Hilbert space.
This includes but is not limited to physical density matrices (which additionally
require positive semidefiniteness and unit trace). The proof below applies to
all SE operators, not only physical states.

---

## Proof

### Part 1: The Liouvillian preserves n_XY parity

The Liouvillian is L = L_H + L_D where:
- L_H(ρ) = −i[H, ρ] is the Hamiltonian part
- L_D(ρ) = Σ_k γ_k(Z_k ρ Z_k − ρ) is the Z-dephasing part

**L_D preserves n_XY exactly.** The dissipator acts on Pauli strings as
L_D(P) = γ_k(Z_k P Z_k − P). Conjugation by Z_k flips the sign of
X_k and Y_k but leaves I_k and Z_k unchanged. Therefore Z_k P Z_k is
either +P (if site k has I or Z) or −P (if site k has X or Y). In both
cases, L_D(P) is proportional to P itself. The n_XY count is unchanged.

**L_H preserves n_XY parity.** The Heisenberg Hamiltonian has bond terms
of the form X_i X_j + Y_i Y_j + Z_i Z_j. We show that commutation with
any bond term preserves n_XY parity of a Pauli string P.

Consider one bond term, say X_i X_j. The commutator [X_i X_j, P] involves
the products X_i X_j P and P X_i X_j. At each site, a single-site Pauli
from the bond term multiplies the corresponding single-site Pauli from P.
The effect on n_XY at each site:

    X · I = X   (n_XY: 0 → 1, change +1)
    X · X = I   (n_XY: 1 → 0, change −1)
    X · Y = iZ  (n_XY: 1 → 0, change −1)
    X · Z = −iY (n_XY: 0 → 1, change +1)

At each bond site, n_XY changes by exactly ±1 (odd). The bond term acts
on two sites, so the total n_XY change from the product X_i X_j · P is
(±1) + (±1) = 0 or ±2, always even. The same holds for Y_i Y_j
(by the same argument with Y replacing X).

For Z_i Z_j, the site-wise changes are:

    Z · I = Z   (n_XY: 0 → 0, change 0)
    Z · X = −iY (n_XY: 1 → 1, change 0)
    Z · Y = iX  (n_XY: 1 → 1, change 0)
    Z · Z = I   (n_XY: 0 → 0, change 0)

Z preserves n_XY exactly at every site, so Z_i Z_j preserves it too.

Since each bond term (XX, YY, ZZ) preserves n_XY parity in every product,
and the commutator [bond, P] = bond·P − P·bond is a sum of such products,
L_H preserves n_XY parity.

Therefore L = L_H + L_D preserves n_XY parity. QED (Part 1).

### Part 2: The Liouvillian block-diagonalizes by n_XY parity

Since L preserves n_XY parity, the operator space splits into two
invariant subspaces:
- V_even = span{P : n_XY(P) is even}
- V_odd = span{P : n_XY(P) is odd}

L maps V_even to V_even and V_odd to V_odd. The eigenmodes of L can
therefore be chosen to lie entirely within one subspace. Every eigenmode
has definite n_XY parity. QED (Part 2).

### Part 3: SE density matrices live in V_even

Let |e_j⟩ and |e_k⟩ be SE basis states. Consider the operator
|e_j⟩⟨e_k| and its Pauli decomposition.

For any Pauli string P with n_XY(P) odd (i.e., P contains an odd number
of X and Y factors):

    ⟨e_k| P |e_j⟩ = 0     for all j, k ∈ {0, …, N−1}

**Proof:** Each X or Y factor in P flips the bit at that site (X|0⟩ = |1⟩,
X|1⟩ = |0⟩, and similarly for Y up to phase). Each I or Z factor leaves the
bit unchanged (up to phase). Therefore P|e_j⟩ is a computational basis
state (times a phase) whose bit pattern differs from |e_j⟩ by exactly
n_XY(P) bit flips.

Each bit flip changes the Hamming weight (number of 1-bits) by exactly ±1.
An odd number of such changes flips the Hamming weight parity. Since |e_j⟩
has Hamming weight 1 (odd parity), P|e_j⟩ has even Hamming weight parity
when n_XY(P) is odd. But |e_k⟩ has Hamming weight 1 (odd parity).
Computational basis states of different Hamming weight are orthogonal,
and in particular states of different Hamming weight parity are orthogonal.
Therefore ⟨e_k|P|e_j⟩ = 0.

Since ⟨e_k|P|e_j⟩ = 0 for all odd-n_XY strings P and all SE indices j, k,
the Pauli decomposition of any SE operator Σ_{j,k} c_{jk} |e_j⟩⟨e_k|
contains no odd-n_XY Pauli strings. Every SE density matrix lives
entirely in V_even. QED (Part 3).

### Part 4: The accessibility boundary

Let ρ₀ = |ψ⟩⟨ψ| be any pure SE initial state. Its Pauli
decomposition lies in V_even (Part 3). The overlap of ρ₀ with any
Liouvillian eigenmode v_m is:

    c_m = Tr(L_m† · ρ₀)

where L_m is the left eigenvector of mode m.

**Left eigenvectors respect the block structure.** Since L is block-diagonal
on V_even ⊕ V_odd (Part 2), L† is also block-diagonal on the same
decomposition: for any u ∈ V_even and v ∈ V_odd, ⟨L†u, v⟩ = ⟨u, Lv⟩ = 0
because Lv is in V_odd and u is in V_even. Therefore L† maps V_even to
V_even and V_odd to V_odd. The left eigenvectors of L, which are the right
eigenvectors of L†, inherit definite parity. In particular, if the right
eigenvector v_m lies in V_odd, the corresponding left eigenvector L_m also
lies in V_odd.

Now: ρ₀ lies in V_even and L_m lies in V_odd. Since V_even and V_odd are
orthogonal in the Hilbert-Schmidt inner product (they are spanned by disjoint
sets of orthonormal Pauli strings):

    c_m = Tr(L_m† · ρ₀) = 0     exactly.

No SE state has any overlap with any odd-n_XY eigenmode. This is not
numerical coincidence but algebraic necessity. QED (Part 4).

---

## What this means physically

The bit-flip parity of X and Y operators creates an exact selection rule.
Single-excitation states carry one quantum of excitation (odd Hamming
weight). The density matrices they form involve only even numbers of
bit flips (correlations between SE states). Liouvillian modes in the
odd-n_XY sector involve an odd number of bit flips, which would change
the excitation parity. The two sectors cannot talk to each other.

This is analogous to a dipole selection rule in atomic physics:
electric dipole transitions change the angular momentum quantum number
by exactly 1, so states of the same parity cannot be connected by a
single dipole transition. Here, single-excitation states have even
n_XY parity, and odd-n_XY modes are the "forbidden transitions".

---

## Scope and limitations

**Valid for:** any isotropic Heisenberg Hamiltonian (XX+YY+ZZ coupling),
any graph topology, any site-dependent Z-dephasing profile γ_k,
any N.

**Also valid for:** XY model (XX+YY only), since the ZZ term has n_XY=0
and does not affect the parity argument.

**Breaks for:** anisotropic models where the Hamiltonian contains terms
with odd n_XY (e.g., a transverse field h·X_k adds n_XY=1 terms to the
commutator, which would mix V_even and V_odd). Also breaks for
**amplitude damping** (T₁ decay), whose jump operator σ₋ = (X − iY)/2
has n_XY = 1 (odd parity). Since real hardware has both T₁ and T₂ noise,
this theorem applies strictly only to the pure-dephasing component.
On hardware where T₂ ≫ T₁ (dephasing-dominated regime), the selection
rule holds approximately; on hardware where T₁ ~ T₂, it does not.

**Does not say:** which specific modes are even or odd, only that the
decomposition exists. The actual distribution of modes between sectors
depends on N, topology, and γ profile.

---

## Numerical verification

The Lens Pipeline survey (April 10, 2026) tested 64 configurations:
N=2–6, Chain/Star/Ring/Complete topologies, four γ profiles each.
In every configuration, the second slowest Liouvillian mode had
SE-sector Frobenius ratio < 1e−3 (machine zero). 64/64 configurations
confirm the selection rule.

Engine: `compute/RCPsiSquared.Compute/LensAnalysis.cs` (run via `dotnet run -c Release -- lens`).
Full data: `simulations/results/lens_survey/lens_survey_scaling.txt`.
See References section for complete file list.

---

## Connection to existing results

**Absorption Theorem (AT):** The AT gives Re(λ) = −2γ⟨n_XY⟩.
The selection rule adds: the eigenmodes that SE states can access are
restricted to V_even. Within V_even, the AT determines the rates.

**Weight-1 Degeneracy (F50):** The 2N conserved operators T_c^{(a)} at
Re = −2γ are weight-1 Pauli operators (n_XY = 1, odd parity).
They live in V_odd by this theorem, confirming they are inaccessible
to SE states from a completely different direction.

**Sacrifice Geometry Phase 3:** The inaccessible mode at rate −0.167
was found to have 100% n_XY=1 content. This theorem explains why:
n_XY=1 is odd, SE density matrices are even, overlap is zero.

---

## References

### Proofs this builds on

- [Absorption Theorem Proof](PROOF_ABSORPTION_THEOREM.md) - Re(λ) = −2γ⟨n_XY⟩
- [Weight-1 Degeneracy Proof](PROOF_WEIGHT1_DEGENERACY.md) - the 2N conserved T_c^{(a)} operators
- [Analytical Formulas](../ANALYTICAL_FORMULAS.md) - AT, F50, **F61** (this theorem's formula entry)

### Numerical verification (64 configurations, N=2–6)

- Engine: `compute/RCPsiSquared.Compute/LensAnalysis.cs`
- Run mode: `dotnet run -c Release -- lens`
- Results: `simulations/results/lens_survey/lens_survey_results.json` (full structured data)
- Summary: `simulations/results/lens_survey/lens_survey_summary.txt`
- Scaling: `simulations/results/lens_survey/lens_survey_scaling.txt` (SE fraction, monotonicity, accessibility per config)

### Original discovery (N=5, Python, April 9 2026)

- Script: `simulations/slow_mode_lens_analysis.py`
- Results: `simulations/results/slow_mode_lens/slow_mode_lens_results.txt`
- Writeup: [Sacrifice Geometry](../../experiments/SACRIFICE_GEOMETRY.md) Phase 3
