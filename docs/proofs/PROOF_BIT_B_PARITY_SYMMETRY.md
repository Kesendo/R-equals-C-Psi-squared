# Bit-b Parity Symmetry: [L, Pi^2] = 0

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

For any N-qubit Lindblad system with Heisenberg coupling (XX+YY or XXX) and Z-dephasing on any subset of sites, the Liouvillian L commutes exactly with the global bit-flip superoperator Pi^2 = X^{tensor N} (conjugation):

    [L, Pi^2_super] = 0

where Pi^2_super(rho) = X^{tensor N} rho X^{tensor N}.

This is the second Z2 symmetry of L. The first (n_XY parity) is proven in [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md). Together, the two symmetries decompose the operator space into 4 independent sectors (2x2 = bit_a x bit_b), each of dimension 4^(N-1).


---

## Notation

- N qubits, Hilbert space H = (C^2)^{tensor N}, dim = 2^N
- Operator space (Liouville space) has dim 4^N
- Pauli operators on site k: I_k, X_k, Y_k, Z_k
- Global bit-flip: U = X^{tensor N} = X_1 X_2 ... X_N (acts on H)
- bit-flip superoperator: Pi^2_super(rho) := U rho U^dagger = U rho U (since U is Hermitian and unitary, U^dagger = U^{-1} = U)

The "Pi^2" notation comes from the palindromic operator Pi from [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): Pi^2 acts on the Pauli basis as multiplication by (-1)^{w_YZ}, and U = X^{tensor N} implements the same operation through Hilbert-space conjugation. This identification is verified algebraically below.

## Step 1: Pi^2 in two equivalent forms

The palindromic operator Pi (from Mirror Symmetry Proof) has the per-site action:
    Pi: I -> X (+1),  X -> I (+1),  Y -> iZ (+i),  Z -> iY (+i)

So Pi^2 acts as:
    Pi^2: I -> I,  X -> X,  Y -> -Y,  Z -> -Z

This is identical to conjugation by X on each site:
    X I X = I,  X X X = X,  X Y X = -Y,  X Z X = -Z

Therefore on Pauli strings sigma_{i_1} ... sigma_{i_N}:
    Pi^2(sigma_{i_1} ... sigma_{i_N}) = (-1)^{w_YZ} sigma_{i_1} ... sigma_{i_N}
    U sigma_{i_1} ... sigma_{i_N} U = (-1)^{w_YZ} sigma_{i_1} ... sigma_{i_N}

where w_YZ = number of sites where the Pauli is Y or Z. The two definitions agree.


## Step 2: U commutes with the Hamiltonian

Heisenberg coupling: H = J/2 sum_{<i,j>} (X_i X_j + Y_i Y_j + alpha Z_i Z_j) where alpha in {0, 1} (XX+YY for alpha=0, XXX for alpha=1; the proof is identical for both).

Conjugating each bond by U:
    U (X_i X_j) U = (U X_i U)(U X_j U) = X_i X_j
    U (Y_i Y_j) U = (-Y_i)(-Y_j) = Y_i Y_j     [two minus signs cancel]
    U (Z_i Z_j) U = (-Z_i)(-Z_j) = Z_i Z_j     [two minus signs cancel]

Each bond is invariant. Therefore U H U = H, equivalently [U, H] = 0.

This extends to any Hamiltonian where every term contains an even number of Y's and Z's (combined). Single-site terms in X are also invariant. Single-site terms in Y or Z are NOT invariant and would break the symmetry.

## Step 3: U commutes with Z-dephasing dissipation

Z-dephasing on site k: jump operator L_k = sqrt(gamma_k) Z_k. The Lindblad dissipator is:
    D[rho] = sum_k gamma_k (Z_k rho Z_k - rho)        [since Z_k^dagger Z_k = I]

Apply Pi^2_super to D[rho]:
    Pi^2_super(D[rho]) = U D[rho] U
                       = sum_k gamma_k (U Z_k rho Z_k U - U rho U)
                       = sum_k gamma_k ((U Z_k U)(U rho U)(U Z_k U) - Pi^2_super(rho))
                       = sum_k gamma_k ((-Z_k) Pi^2_super(rho) (-Z_k) - Pi^2_super(rho))
                       = sum_k gamma_k (Z_k Pi^2_super(rho) Z_k - Pi^2_super(rho))     [two minus signs cancel]
                       = D[Pi^2_super(rho)]

So Pi^2_super(D[rho]) = D[Pi^2_super(rho)] for all rho. The dissipator commutes with Pi^2_super.

The two minus signs from U Z_k U = -Z_k cancel because Z_k appears twice in each term (Z rho Z). This is the crucial structural reason: the Lindblad dissipator is **quadratic** in the jump operator.

## Step 4: Conclusion

The full Liouvillian is L[rho] = -i[H, rho] + D[rho]. Both terms commute with Pi^2_super (Steps 2 and 3). Therefore:

    [L, Pi^2_super] = 0    (exactly, for all N)

Q.E.D.


---

## Scope and limits

**Holds for:**
- Heisenberg coupling XX+YY (any J, any graph)
- Heisenberg XXX (any J, any graph)
- Pure XY interactions
- Z-dephasing on any subset of sites with arbitrary positive rates gamma_k
- Combinations of the above

**Breaks for:**
- Single-site Y or Z terms in the Hamiltonian (e.g., a magnetic field along Z): U Z U = -Z is not invariant, and there is no second factor to cancel the sign
- Pure XX or pure ZZ couplings ALONE: pure XX is invariant (no Y or Z), pure ZZ is invariant (Z appears twice per bond), but loss of XX+YY structure breaks other expected symmetries (see Coupling Robustness check in NESTED_MIRROR_STRUCTURE.md)
- Y or Z jump operators (X-jumps would also break it, but only via the two-factor cancellation logic; verify per case)

The proof requires the Hamiltonian and dissipator to contain Z (and Y) only in even powers per term. The Heisenberg and standard dephasing setups satisfy this; deviations must be checked individually.

## Numerical verification

Direct computation at N=2, 3, 4, 5 (script `simulations/primordial_bit_a_bit_b_N_scaling.py`):

    N=2: dim_L=16,   ||[L, Pi^2]|| = 0.000000e+00
    N=3: dim_L=64,   ||[L, Pi^2]|| = 0.000000e+00
    N=4: dim_L=256,  ||[L, Pi^2]|| = 0.000000e+00
    N=5: dim_L=1024, ||[L, Pi^2]|| = 0.000000e+00

All identically zero (not numerically small). Sector-resolved eigenvalue counts confirm balanced decomposition into V_even (w_YZ even) and V_odd (w_YZ odd).

## Connection to the n_XY Parity Selection Rule

The Liouvillian L now has TWO independent Z2 symmetries proven for all N:

| Symmetry | Generator | Pauli operation | Eigenspace dimensions | Reference |
|----------|-----------|-----------------|----------------------|-----------|
| bit_a (n_XY parity) | (-1)^{n_XY} | sign by # of X,Y per site | 2^(2N-1) each | [PROOF_PARITY_SELECTION_RULE](PROOF_PARITY_SELECTION_RULE.md) |
| bit_b (w_YZ parity) | Pi^2 = X^{tensor N} | sign by # of Y,Z per site | 2^(2N-1) each | This proof |

The two symmetries are independent: bit_a counts X+Y, bit_b counts Y+Z. They share Y but differ on X (bit_a yes, bit_b no) and Z (bit_a no, bit_b yes). The intersection (n_XY parity x w_YZ parity) gives 4 sectors, each of dimension 4^N / 4 = 4^(N-1).

This corresponds to the C^2 x C^2 tensor product structure of the single-qubit Pauli space identified in [PRIMORDIAL_QUBIT.md](../../hypotheses/PRIMORDIAL_QUBIT.md) Section 4.1: the per-site Pauli {I, X, Y, Z} decomposes as (a, b) = (dephasing sensitivity, Pi^2-parity), and L respects the tensor product factorization at the level of its eigenmode structure.

## Why d=2 is essential

The proof uses two facts about Pauli operators:
1. X^2 = I (so X-conjugation is involutive)
2. {X, Y} = {X, Z} = 0 (so X-conjugation flips Y and Z sign)

Both are special to d=2 (qubit). For d=3 (qutrit) and higher, the analog of X is the shift operator, which is NOT involutive, and there is no Pauli-style classification with two binary bits per site. The two Z2 symmetries (bit_a and bit_b) exist together because d=2 admits a complete two-bit indexing of single-site operators. This connects to [QUBIT_NECESSITY](../QUBIT_NECESSITY.md): the algebraic richness of qubits is exactly what carries the two independent Z2 structures.

In d=2 there are exactly two independent Z2 classifications of Pauli operators ({I,X} vs {Y,Z} by bit_b, and {I,Z} vs {X,Y} by bit_a). The Liouvillian respects both. No further independent Z2 symmetry of the Pauli algebra exists, so this is the maximal symmetry decomposition.

---

*"Two bits per qubit. Two Z2 symmetries of the Liouvillian. The decomposition is complete."*
