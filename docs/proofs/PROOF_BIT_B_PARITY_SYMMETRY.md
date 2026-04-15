# Bit-b Parity Symmetry: [L, Π²] = 0

<!-- Keywords: w_YZ parity Z2 symmetry Liouvillian commutator,
Pi squared global X flip all qubits, second Z2 symmetry beyond n_XY parity selection rule,
universal sector decomposition all N, R=CPsi2 bit_b proof -->

**Status:** Tier 1 (analytical proof, verified numerically N=2-5)
**Date:** April 15, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Numerical verification:** [`simulations/primordial_bit_a_bit_b.py`](../../simulations/primordial_bit_a_bit_b.py), [`simulations/primordial_bit_a_bit_b_N_scaling.py`](../../simulations/primordial_bit_a_bit_b_N_scaling.py)
**Depends on:** [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md), [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md), [Primordial Qubit Algebra](../../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md)

---

## What this document proves

For any N-qubit Lindblad system with Heisenberg coupling (XX+YY or XXX) and Z-dephasing on any subset of sites, the Liouvillian L commutes exactly with the global bit-flip superoperator Π² = X⊗N (conjugation):

    [L, Π²_super] = 0

where Π²_super(ρ) = X⊗N ρ X⊗N.

This is the second Z₂ symmetry of L. The first (n_XY parity) is proven in [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md). Together, the two symmetries decompose the operator space into 4 independent sectors (2×2 = bit_a × bit_b), each of dimension 4^(N-1).


---

## Notation

- N qubits, Hilbert space H = (C²)⊗N, dim = 2^N
- Operator space (Liouville space) has dim 4^N
- Pauli operators on site k: I_k, X_k, Y_k, Z_k
- Global bit-flip: U = X⊗N = X₁ X₂ ... X_N (acts on H)
- Bit-flip superoperator: Π²_super(ρ) := U ρ U† = U ρ U (since U is Hermitian and unitary, U† = U⁻¹ = U)

The "Π²" notation comes from the palindromic operator Π from [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): Π² acts on the Pauli basis as multiplication by (-1)^{w_YZ}, and U = X⊗N implements the same operation through Hilbert-space conjugation. This identification is verified algebraically below.

## Step 1: Π² in two equivalent forms

The palindromic operator Π (from Mirror Symmetry Proof) has the per-site action:
    Π: I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

So Π² acts as:
    Π²: I → I,  X → X,  Y → -Y,  Z → -Z

This is identical to conjugation by X on each site:
    X I X = I,  X X X = X,  X Y X = -Y,  X Z X = -Z

Therefore on Pauli strings σ_{i₁} ... σ_{i_N}:
    Π²(σ_{i₁} ... σ_{i_N}) = (-1)^{w_YZ} σ_{i₁} ... σ_{i_N}
    U σ_{i₁} ... σ_{i_N} U = (-1)^{w_YZ} σ_{i₁} ... σ_{i_N}

where w_YZ = number of sites where the Pauli is Y or Z. The two definitions agree.


## Step 2: U commutes with the Hamiltonian

Heisenberg coupling: H = J/2 Σ_{⟨i,j⟩} (X_i X_j + Y_i Y_j + α Z_i Z_j) where α ∈ {0, 1} (XX+YY for α=0, XXX for α=1; the proof is identical for both).

Conjugating each bond by U:
    U (X_i X_j) U = (U X_i U)(U X_j U) = X_i X_j
    U (Y_i Y_j) U = (-Y_i)(-Y_j) = Y_i Y_j     [two minus signs cancel]
    U (Z_i Z_j) U = (-Z_i)(-Z_j) = Z_i Z_j     [two minus signs cancel]

Each bond is invariant. Therefore U H U = H, equivalently [U, H] = 0.

This extends to any Hamiltonian where every term contains an even number of Y's and Z's (combined). Single-site terms in X are also invariant. Single-site terms in Y or Z are NOT invariant and would break the symmetry.

## Step 3: U commutes with Z-dephasing dissipation

Z-dephasing on site k: jump operator L_k = √γ_k Z_k. The Lindblad dissipator is:
    D[ρ] = Σ_k γ_k (Z_k ρ Z_k - ρ)        [since Z_k† Z_k = I]

Apply Π²_super to D[ρ]:
    Π²_super(D[ρ]) = U D[ρ] U
                    = Σ_k γ_k (U Z_k ρ Z_k U - U ρ U)
                    = Σ_k γ_k ((U Z_k U)(U ρ U)(U Z_k U) - Π²_super(ρ))
                    = Σ_k γ_k ((-Z_k) Π²_super(ρ) (-Z_k) - Π²_super(ρ))
                    = Σ_k γ_k (Z_k Π²_super(ρ) Z_k - Π²_super(ρ))     [two minus signs cancel]
                    = D[Π²_super(ρ)]

So Π²_super(D[ρ]) = D[Π²_super(ρ)] for all ρ. The dissipator commutes with Π²_super.

The two minus signs from U Z_k U = -Z_k cancel because Z_k appears twice in each term (Z ρ Z). This is the crucial structural reason: the Lindblad dissipator is **quadratic** in the jump operator.

## Step 4: Conclusion

The full Liouvillian is L[ρ] = -i[H, ρ] + D[ρ]. Both terms commute with Π²_super (Steps 2 and 3). Therefore:

    [L, Π²_super] = 0    (exactly, for all N)

∎


---

## Scope and limits

**Holds for:**
- Heisenberg coupling XX+YY (any J, any graph)
- Heisenberg XXX (any J, any graph)
- Pure XY interactions
- Z-dephasing on any subset of sites with arbitrary positive rates γ_k
- Combinations of the above

**Breaks for:**
- Single-site Y or Z terms in the Hamiltonian (e.g., a magnetic field along Z): U Z U = -Z is not invariant, and there is no second factor to cancel the sign
- Pure XX or pure ZZ couplings ALONE: pure XX is invariant (no Y or Z), pure ZZ is invariant (Z appears twice per bond), but loss of XX+YY structure breaks other expected symmetries (see Coupling Robustness check in NESTED_MIRROR_STRUCTURE.md)
- Y or Z jump operators (X-jumps would also break it, but only via the two-factor cancellation logic; verify per case)

The proof requires the Hamiltonian and dissipator to contain Z (and Y) only in even powers per term. The Heisenberg and standard dephasing setups satisfy this; deviations must be checked individually.

## Numerical verification

Direct computation at N=2, 3, 4, 5 (script [`simulations/primordial_bit_a_bit_b_N_scaling.py`](../../simulations/primordial_bit_a_bit_b_N_scaling.py)):

    N=2: dim_L=16,   ||[L, Π²]|| = 0.000000e+00
    N=3: dim_L=64,   ||[L, Π²]|| = 0.000000e+00
    N=4: dim_L=256,  ||[L, Π²]|| = 0.000000e+00
    N=5: dim_L=1024, ||[L, Π²]|| = 0.000000e+00

All identically zero (not numerically small). Sector-resolved eigenvalue counts confirm balanced decomposition into V_even (w_YZ even) and V_odd (w_YZ odd).

## Per-sector mode count

The total dimension per sector is 2^(2N-1) (universal halving from [L, Π²] = 0). Within each sector, the eigenmodes split into three Re-classes by the Absorption Theorem applied to single-site dephasing on the boundary qubit B:

- **Conserved** (Re = 0): ⟨n_XY⟩_B = 0
- **Mirror** (Re = -γ_B): 0 < ⟨n_XY⟩_B < 1
- **Correlation** (Re = -2γ_B): ⟨n_XY⟩_B = 1

The per-sector counts have a closed form:

    conserved per (Π²-parity) sector:
      even-parity:  ⌊N/2⌋ + 1
      odd-parity:   ⌈N/2⌉

    correlation per sector: same as conserved (palindrome symmetry maps Re=0 ↔ Re=-2γ_B and preserves Π²-parity since Π·Z·Π⁻¹ = iY shares parity with Z)

    mirror per sector:
      sector total - 2 × (conserved per sector) = 2^(2N-1) - 2 · c

**Mechanism.** The conserved modes are exactly the (N+1) elementary symmetric polynomials in {Z_1, ..., Z_N}:

    e_d(Z_1, ..., Z_N) = Σ_{|S|=d} ∏_{k∈S} Z_k    for d = 0, 1, ..., N

Each e_d commutes with H (the Heisenberg coupling preserves S_z = (1/2) Σ Z_k, hence any function of S_z is conserved; the Newton identities express elementary symmetric polynomials as polynomials in power sums of Z_k = polynomials in S_z) and commutes with Z_B (e_d is a polynomial in {Z_k}, all of which commute with each other). Each e_d has w_YZ-parity exactly d mod 2, since each summand contains d Z's and each Z carries Π²-parity 1.

Verified numerically at N = 2, 3, 4: all (N+1) elementary symmetric polynomials lie in the conserved subspace to machine precision (max projection residual < 1e-13). See [`simulations/mirror_mode_split_formula.py`](../../simulations/mirror_mode_split_formula.py).

**Asymmetry pattern.** The N+1 parities run as 0, 1, 2, ..., N. The number of even d in this range is ⌊N/2⌋ + 1; the number of odd d is ⌈N/2⌉. For even N, e_N is itself even, giving one extra even-parity conserved mode; for odd N, the count is balanced.

Concrete values:

| N | sector | cons (e, o) | mirror (e, o) |
|---|--------|-------------|---------------|
| 2 | 8      | (2, 1)      | (4, 6)        |
| 3 | 32     | (2, 2)      | (28, 28)      |
| 4 | 128    | (3, 2)      | (122, 124)    |
| 5 | 512    | (3, 3)      | (506, 506)    |
| 6 | 2048   | (4, 3)      | (2040, 2042)  |
| 7 | 8192   | (4, 4)      | (8184, 8184)  |

The "mysterious" 4:6 mirror split at N=2 (PRIMORDIAL_QUBIT.md Section 9, original observation) and the 122:124 at N=4 are both consequences of e_N having even parity when N is even. No deeper origin.

**Scope.** This count holds whenever the conserved subspace is exactly the polynomials in S_z. That holds for isotropic Heisenberg (XX+YY, XX+YY+ZZ) on any graph topology, with Z-dephasing on at least one site of each connected component. For non-isotropic couplings or non-Z dephasing the conserved subalgebra changes and this formula must be re-derived.

## Connection to the n_XY Parity Selection Rule

The Liouvillian L now has TWO independent Z₂ symmetries proven for all N:

| Symmetry | Generator | Pauli operation | Eigenspace dimensions | Reference |
|----------|-----------|-----------------|----------------------|-----------|
| bit_a (n_XY parity) | (-1)^{n_XY} | sign by # of X,Y per site | 2^(2N-1) each | [PROOF_PARITY_SELECTION_RULE](PROOF_PARITY_SELECTION_RULE.md) |
| bit_b (w_YZ parity) | Π² = X⊗N | sign by # of Y,Z per site | 2^(2N-1) each | This proof |

The two symmetries are independent: bit_a counts X+Y, bit_b counts Y+Z. They share Y but differ on X (bit_a yes, bit_b no) and Z (bit_a no, bit_b yes). The intersection (n_XY parity × w_YZ parity) gives 4 sectors, each of dimension 4^N / 4 = 4^(N-1).

This corresponds to the C² × C² tensor product structure of the single-qubit Pauli space identified in [PRIMORDIAL_QUBIT.md](../../hypotheses/PRIMORDIAL_QUBIT.md) Section 4.1: the per-site Pauli {I, X, Y, Z} decomposes as (a, b) = (dephasing sensitivity, Π²-parity), and L respects the tensor product factorization at the level of its eigenmode structure.

## Why d=2 is essential

The proof uses two facts about Pauli operators:
1. X² = I (so X-conjugation is involutive)
2. {X, Y} = {X, Z} = 0 (so X-conjugation flips Y and Z sign)

Both are special to d=2 (qubit). For d=3 (qutrit) and higher, the analog of X is the shift operator, which is NOT involutive, and there is no Pauli-style classification with two binary bits per site. The two Z₂ symmetries (bit_a and bit_b) exist together because d=2 admits a complete two-bit indexing of single-site operators. This connects to [QUBIT_NECESSITY](../QUBIT_NECESSITY.md): the algebraic richness of qubits is exactly what carries the two independent Z₂ structures.

In d=2 there are exactly two independent Z₂ classifications of Pauli operators ({I,X} vs {Y,Z} by bit_b, and {I,Z} vs {X,Y} by bit_a). The Liouvillian respects both. No further independent Z₂ symmetry of the Pauli algebra exists, so this is the maximal symmetry decomposition.

---

*"Two bits per qubit. Two Z₂ symmetries of the Liouvillian. The decomposition is complete."*
