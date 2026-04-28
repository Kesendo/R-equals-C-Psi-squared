"""Framework primitives for R=CΨ²: ur-Pauli, ur-Heisenberg, ur-eigenvalues.

Built from the framework's foundational structure:
- d²−2d=0 selects d=2 (uniquely solvable for the qubit immune-decaying balance)
- Pauli space at d=2 is C²⊗C², indexed by two Z₂ parities:
    bit_a (n_XY): dephasing axis. 0=immune (I, Z), 1=decaying (X, Y).
    bit_b (n_YZ): palindromic axis. 0=Π²-even (I, X), 1=Π²-odd (Y, Z).
- Pauli operators are labeled by (bit_a, bit_b) ∈ {0,1}² in Weyl form:
    P_{a,b} = i^{a·b} · X^a · Z^b
    (0,0)=I, (1,0)=X, (0,1)=Z, (1,1)=Y

From these primitives:
- The Heisenberg/XXZ form is the unique both-parity-even 2-body bilinear:
    {II, XX, YY, ZZ} = the (a,a)·(a,a)·... diagonal sum on bonds.
- The Pauli identity (σ·σ)² = 3I − 2(σ·σ) follows from σ²_a = I and Clifford
  anti-commutation, giving the eigenvalue equation x² + 2x − 3 = 0
  with roots +1 (triplet, 3-fold) and −3 (singlet, 1-fold).
- The Π conjugation per site flips bit_a, preserves bit_b (with phase i if
  bit_b=1).

Use this module instead of redefining σ_x, σ_y, σ_z for every simulation.
"""
import math
import sys
from itertools import product as iproduct

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# ----------------------------------------------------------------------
# Section 1: ur-Pauli via (bit_a, bit_b) indexing
# ----------------------------------------------------------------------

# The four Pauli operators at d=2, indexed by (bit_a, bit_b)
# (0,0) = I (identity, immune to Z-dephasing, Π²-even)
# (1,0) = X (decaying, Π²-even)
# (0,1) = Z (immune, Π²-odd)
# (1,1) = Y (decaying, Π²-odd, picks up i in Weyl form)

_PAULI_MATRICES = {
    (0, 0): np.eye(2, dtype=complex),
    (1, 0): np.array([[0, 1], [1, 0]], dtype=complex),
    (0, 1): np.array([[1, 0], [0, -1]], dtype=complex),
    (1, 1): np.array([[0, -1j], [1j, 0]], dtype=complex),
}

PAULI_LABELS = {(0, 0): 'I', (1, 0): 'X', (0, 1): 'Z', (1, 1): 'Y'}
LABEL_TO_INDEX = {v: k for k, v in PAULI_LABELS.items()}


def ur_pauli(a, b=None):
    """Return the 2×2 Pauli matrix at index (a, b) or by label string.

    Examples:
        ur_pauli(0, 0) → I
        ur_pauli(1, 0) → X
        ur_pauli('Y') → Y (looked up by label)
    """
    if b is None and isinstance(a, str):
        return _PAULI_MATRICES[LABEL_TO_INDEX[a]]
    return _PAULI_MATRICES[(a, b)]


def bit_a(idx):
    """bit_a (dephasing parity, n_XY) of a Pauli (idx is (a,b) tuple or label)."""
    if isinstance(idx, str):
        idx = LABEL_TO_INDEX[idx]
    return idx[0]


def bit_b(idx):
    """bit_b (Π²-parity, n_YZ) of a Pauli (idx is (a,b) tuple or label)."""
    if isinstance(idx, str):
        idx = LABEL_TO_INDEX[idx]
    return idx[1]


def total_bit_a(indices):
    """XY-weight of a Pauli string (total n_XY across sites)."""
    return sum(bit_a(idx) for idx in indices)


def total_bit_b_parity(indices):
    """Π²-parity of a Pauli string (total n_YZ mod 2)."""
    return sum(bit_b(idx) for idx in indices) % 2


# ----------------------------------------------------------------------
# Section 2: Multi-qubit Pauli strings
# ----------------------------------------------------------------------

def pauli_string(indices_or_labels):
    """Build a 2^N × 2^N Pauli string from a list of N (a,b) tuples or labels.

    Examples:
        pauli_string([(1,0), (1,0)]) → X⊗X
        pauli_string(['X', 'X']) → X⊗X
        pauli_string(['I', 'Z', 'Z']) → I⊗Z⊗Z (a w=0 string at N=3)
    """
    out = ur_pauli(*_resolve(indices_or_labels[0]))
    for idx in indices_or_labels[1:]:
        out = np.kron(out, ur_pauli(*_resolve(idx)))
    return out


def _resolve(idx):
    """Convert label string or (a,b) tuple into (a,b) tuple."""
    if isinstance(idx, str):
        return LABEL_TO_INDEX[idx]
    return idx


def site_op(N, site, idx_or_label):
    """N-qubit operator with single-Pauli idx on the given site, I elsewhere."""
    ops = [ur_pauli('I')] * N
    ops[site] = ur_pauli(*_resolve(idx_or_label))
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


# ----------------------------------------------------------------------
# Section 3: Π conjugation (per site and full)
# ----------------------------------------------------------------------

# Π per site: I↔X (sign 1), Y↔Z (sign i). In (bit_a, bit_b) language:
# bit_a flips, bit_b preserved, phase = i if bit_b=1.

def pi_action(idx):
    """Return ((new_a, new_b), phase) for Π acting on Pauli at idx."""
    a, b = _resolve(idx)
    new_a = 1 - a   # bit_a flips
    new_b = b       # bit_b preserved
    phase = 1j if b == 1 else 1
    return (new_a, new_b), phase


def pi_squared_eigenvalue(indices):
    """Π² eigenvalue on a Pauli string = (-1)^{total_bit_b}."""
    return (-1) ** total_bit_b_parity(indices)


def build_pi_full(N):
    """Π in the 4^N Pauli-string basis: 4^N × 4^N matrix."""
    d2 = 4 ** N
    Pi = np.zeros((d2, d2), dtype=complex)
    for k in range(d2):
        indices = _k_to_indices(k, N)
        new_indices = []
        sign = 1
        for idx in indices:
            (na, nb), phase = pi_action(idx)
            new_indices.append((na, nb))
            sign *= phase
        new_k = _indices_to_k(new_indices)
        Pi[new_k, k] = sign
    return Pi


def _k_to_indices(k, N):
    """Convert flat index k ∈ [0, 4^N) to N-tuple of (a,b) tuples."""
    out = []
    kk = k
    for _ in range(N):
        i = kk % 4
        kk //= 4
        # i ∈ {0,1,2,3}: I, X, Z, Y in some order; use (a, b) encoding
        # We use: 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
        out.append((i & 1, (i >> 1) & 1))
    return tuple(reversed(out))


def _indices_to_k(indices):
    """Inverse of _k_to_indices."""
    k = 0
    for (a, b) in indices:
        i = a + 2 * b
        k = k * 4 + i
    return k


# ----------------------------------------------------------------------
# Section 4: Parity selection: which 2-body bilinears respect framework
# ----------------------------------------------------------------------

def respects_bit_a_parity(idx_pair):
    """Total bit_a even (a₁ + a₂ ≡ 0 mod 2)."""
    return (bit_a(idx_pair[0]) + bit_a(idx_pair[1])) % 2 == 0


def respects_bit_b_parity(idx_pair):
    """Total bit_b even (b₁ + b₂ ≡ 0 mod 2)."""
    return (bit_b(idx_pair[0]) + bit_b(idx_pair[1])) % 2 == 0


def is_both_parity_even(idx_pair):
    """The Heisenberg-form selection: both Z₂ parities even."""
    return respects_bit_a_parity(idx_pair) and respects_bit_b_parity(idx_pair)


def both_parity_even_terms():
    """Return all 2-body bilinears that pass the both-parity-even filter.

    Result: [('I','I'), ('X','X'), ('Y','Y'), ('Z','Z')] = the {II, XX, YY, ZZ}
    selection. This is the Heisenberg/XXZ family (modulo identity).
    """
    out = []
    for a, b in iproduct(['I', 'X', 'Y', 'Z'], repeat=2):
        if is_both_parity_even((a, b)):
            out.append((a, b))
    return out


# ----------------------------------------------------------------------
# Section 5: ur-Heisenberg and ur-XXZ
# ----------------------------------------------------------------------

def ur_heisenberg(N, J=1.0, bonds=None):
    """N-qubit Heisenberg H = J Σ_bond (XX + YY + ZZ).

    By construction: the unique SU(2)-invariant bilinear in
    {II, XX, YY, ZZ} (the both-parity-even set, modulo identity).
    Default bonds: nearest-neighbor open chain.
    """
    if bonds is None:
        bonds = [(i, i + 1) for i in range(N - 1)]
    return _build_bilinear(N, bonds, [('X', 'X', J), ('Y', 'Y', J), ('Z', 'Z', J)])


def ur_xxz(N, J=1.0, delta=1.0, bonds=None):
    """N-qubit XXZ H = J Σ_bond (XX + YY + Δ ZZ).

    Both-parity-even (preserves both Z₂'s). Δ=1 → Heisenberg. Δ=0 → XY-model.
    """
    if bonds is None:
        bonds = [(i, i + 1) for i in range(N - 1)]
    return _build_bilinear(N, bonds, [('X', 'X', J), ('Y', 'Y', J), ('Z', 'Z', J * delta)])


def ur_xyz(N, Jx=1.0, Jy=1.0, Jz=1.0, bonds=None):
    """N-qubit XYZ H = Jx XX + Jy YY + Jz ZZ. Most general both-parity-even."""
    if bonds is None:
        bonds = [(i, i + 1) for i in range(N - 1)]
    return _build_bilinear(N, bonds, [('X', 'X', Jx), ('Y', 'Y', Jy), ('Z', 'Z', Jz)])


def _build_bilinear(N, bonds, terms):
    """Build H = Σ_bond Σ_term coeff · σ_a^l σ_b^{l+1}.

    terms = list of (label1, label2, coeff).
    """
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for (la, lb, coeff) in terms:
            H = H + coeff * site_op(N, i, la) @ site_op(N, j, lb)
    return H


# ----------------------------------------------------------------------
# Section 6: Lindbladian and palindrome residual
# ----------------------------------------------------------------------

def lindbladian_general(H, c_ops):
    """General Lindbladian L(ρ) = -i[H, ρ] + Σ_k (c_k ρ c_k† − ½{c_k†c_k, ρ}).

    c_ops: list of Lindblad jump operators (each a d×d matrix). Returned in the
    operator-vec convention used throughout this repository (flatten('F')):
        L = -i(H⊗I - I⊗Hᵀ) + Σ_k [ c_k ⊗ c_k* − ½(c_k†c_k ⊗ I) − ½(I ⊗ (c_k†c_k)ᵀ) ].
    Compatible with `_vec_to_pauli_basis_transform` and `palindrome_residual`.
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c in c_ops:
        c_dag_c = c.conj().T @ c
        L = L + (np.kron(c, c.conj())
                 - 0.5 * np.kron(c_dag_c, Id)
                 - 0.5 * np.kron(Id, c_dag_c.T))
    return L


def lindbladian_z_dephasing(H, gamma_l):
    """Lindbladian L(ρ) = -i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l - ρ).

    Specific case of `lindbladian_general` with c_l = √γ_l · Z_l. The
    γ·(Z·ρ·Z − ρ) form is the pure-dephasing limit; the framework's
    palindrome (Π·L·Π⁻¹ + L + 2Σγ·I = 0) is derived for this Lindbladian.
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    N = int(round(math.log2(d)))
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        Zl = site_op(N, l, 'Z')
        L = L + gamma * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l):
    """Z-dephasing + T1 amplitude damping.

    L(ρ) = -i[H, ρ]
           + Σ_l γ_l · (Z_l ρ Z_l − ρ)        (pure Z-dephasing)
           + Σ_l γ^{T1}_l · (σ⁻_l ρ σ⁺_l − ½{σ⁺_l σ⁻_l, ρ})  (amplitude damping)

    σ⁻ = (X − iY)/2 is the lowering operator (|1⟩→|0⟩). When γ^{T1}_l = 0
    on every site, this reduces to `lindbladian_z_dephasing`. The framework's
    palindrome holds for the Z-dephasing piece alone; the T1 piece introduces
    breaking that breaks the operator equation. This function lets the caller
    measure how `pi_protected_observables` shifts with increasing T1.
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    N = int(round(math.log2(d)))
    c_ops = []
    # Z-dephasing jump operators
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        c_ops.append(np.sqrt(gamma) * site_op(N, l, 'Z'))
    # T1 amplitude damping (σ⁻ on each site)
    sigma_minus_2 = np.array([[0, 1], [0, 0]], dtype=complex)  # |0⟩⟨1|
    for l, gamma_t1 in enumerate(gamma_t1_l):
        if gamma_t1 == 0:
            continue
        # site_op with σ⁻ — direct construction
        ops = [np.eye(2, dtype=complex)] * N
        ops[l] = sigma_minus_2
        sigma_minus_l = ops[0]
        for op in ops[1:]:
            sigma_minus_l = np.kron(sigma_minus_l, op)
        c_ops.append(np.sqrt(gamma_t1) * sigma_minus_l)
    return lindbladian_general(H, c_ops)


def palindrome_residual(L, Sigma_gamma, N):
    """Compute Π·L·Π⁻¹ + L + 2Σγ·I in Pauli-string basis.

    Returns the residual matrix (4^N × 4^N) for sector-block analysis.
    Uses the orthogonality M†M = 2^N · I of the Pauli basis, so the
    inverse transform is M†/2^N (no solve needed).
    """
    Pi = build_pi_full(N)
    M = _vec_to_pauli_basis_transform(N)
    L_pauli = (M.conj().T @ L @ M) / (2 ** N)
    Pi_inv = np.linalg.inv(Pi)
    return Pi @ L_pauli @ Pi_inv + L_pauli + 2 * Sigma_gamma * np.eye(4 ** N)


def _vec_to_pauli_basis_transform(N):
    """Transformation matrix from Pauli-string basis to column-stack vec."""
    d2 = 4 ** N
    d = 2 ** N
    M = np.zeros((d * d, d2), dtype=complex)
    for k in range(d2):
        indices = _k_to_indices(k, N)
        sigma = pauli_string(list(indices))
        M[:, k] = sigma.flatten('F')
    return M


def palindrome_residual_norm_squared_factor_graph(N, B, D2, hamiltonian_class='main'):
    """Universal N-scaling factor F(N, G) = ‖M(N, G)‖² / c_H for any graph G.

    Generalises the chain formula to arbitrary connectivity:

        main class         ‖M(N, G)‖² = c_H · B(G) · 4^(N − 2)
        single-body class  ‖M(N, G)‖² = c_H · (D2(G) / 2) · 4^(N − 2)

    where:
      B(G)  = number of bonds in the graph
      D2(G) = Σ_i deg_G(i)² = sum of squared vertex degrees

    The chain case recovers (N − 1) and (2N − 3) from B = N − 1 and
    D2 = 4N − 6 (degree sequence 1, 2, 2, ..., 2, 1).

    Verified to machine precision (m/p = 1.000000) on:
      chain (B = N−1)             N = 4, 5
      ring  (B = N)               N = 4, 5
      star  (B = N−1, hub)        N = 4, 5
      K_N   (B = N(N−1)/2)        N = 4, 5

    The c_H anchor is graph-independent: at N=2 chain and any topology
    with one bond it equals ‖M(2, single bond)‖². See
    experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md.

    Args:
        N: number of sites in the graph (>= 2)
        B: bond count of G (must be >= 1)
        D2: Σ_i deg_G(i)² (must equal 2B for bond-count consistency check)
        hamiltonian_class: 'main' or 'single_body'

    Returns:
        F(N, G) such that ‖M(N, G)‖² = c_H · F(N, G), where c_H = ‖M(2)‖²
        for any one-bond Hamiltonian with the same Pauli-pair structure.
    """
    if N < 2:
        raise ValueError(f"N must be >= 2; got {N}")
    if B < 1:
        raise ValueError(f"B must be >= 1 (graph needs at least one bond); got {B}")
    deg_sum = 2 * B  # handshake lemma; cheap consistency check on D2
    if D2 < deg_sum:
        # D2 >= sum(deg) by Cauchy-Schwarz only when each deg >= 1, but
        # we can at least require D2 >= 2B (each end of each bond contributes >=1)
        raise ValueError(
            f"D2 ({D2}) inconsistent with B ({B}); "
            f"D2 must be at least 2B = {deg_sum}"
        )
    if hamiltonian_class == 'main':
        return B * (4 ** (N - 2))
    if hamiltonian_class == 'single_body':
        return (D2 / 2) * (4 ** (N - 2))
    raise ValueError(
        f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}"
    )


def palindrome_residual_norm_squared_factor(N, hamiltonian_class='main'):
    """Universal N-scaling factor F(N) such that ‖M(N)‖_F² = c_H · F(N).

    For chain Pauli-pair Hamiltonians under uniform Z-dephasing, the squared
    Frobenius norm of the palindrome residual M = Π·L·Π⁻¹ + L + 2Σγ·I obeys
    a closed-form *absolute* scaling, not just a ratio law:

        main class         ‖M(N)‖² = c_H · (N − 1) · 4^(N − 2)
        single-body class  ‖M(N)‖² = c_H · (2N − 3) · 4^(N − 2)

    where c_H = ‖M(2)‖² is a per-Hamiltonian anchor constant fixed at the
    minimum chain length (one bond). Verified to machine precision (std
    < 3·10⁻¹⁶) at N ∈ {2, 3, 4, 5} on representative Hamiltonians from each
    class; see experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md.

    The (N − 1) factor counts bonds; the 4^(N − 2) factor is the operator-
    space extension when adding qubits beyond the minimal 2-site bond. The
    single-body class replaces (N − 1) by (2N − 3) = (N − 1) + (N − 2),
    bond count plus interior-site count, reflecting the doubled middle-site
    weight of bilinears of the form (Iσ, σI).

    For non-chain topologies use `palindrome_residual_norm_squared_factor_graph`,
    which takes the bond count and the second moment of the degree sequence
    directly. The chain formulas recover from B = N − 1 and D2 = 4N − 6.

    Args:
        N: chain length (>= 2)
        hamiltonian_class: 'main' or 'single_body'

    Returns:
        F(N) such that ‖M(N)‖_F² = c_H · F(N), where c_H = ‖M(2)‖².
        Caller multiplies by c_H (computed once at N=2) to predict all
        larger N without re-eigenvaluing.

    Example:
        >>> N2_norm_sq = ...  # measure at N=2 once
        >>> for N in [3, 4, 5, 6, 7]:
        ...     predicted = N2_norm_sq * fw.palindrome_residual_norm_squared_factor(N)
    """
    if N < 2:
        raise ValueError(f"N must be >= 2 (chain needs at least one bond); got {N}")
    if hamiltonian_class == 'main':
        return (N - 1) * (4 ** (N - 2))
    if hamiltonian_class == 'single_body':
        return (2 * N - 3) * (4 ** (N - 2))
    raise ValueError(
        f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}"
    )


def palindrome_residual_norm_ratio_squared(N1, N2, hamiltonian_class='main'):
    """Closed-form ‖M(N2)‖_F² / ‖M(N1)‖_F² for chain Pauli-pair Hamiltonians.

    For a chain with uniform-coefficient bilinear Pauli-pair bond Hamiltonians
    under uniform Z-dephasing, the Frobenius norm of the palindrome residual
    M = Π·L·Π⁻¹ + L + 2Σγ·I follows two rational scaling laws between
    adjacent N (N2 = N1 + 1):

        main class         ratio² = 4·k / (k − 1)
        single-body class  ratio² = 4·(2k − 1) / (2k − 3)

    where k = N1. Both → 4 as k → ∞ (the universal d²-extension factor:
    one extra qubit quadruples the operator-space dimension).

    The **main class** consists of Pauli pairs whose bond bilinear has at
    least one genuine 2-body term. For these, ‖M(N)‖² ∝ (N − 1)·4^(N − 2)
    (bond count × per-bond Liouvillian extension), so the ratio is
    4·(N − 1)/(N − 2) evaluated at source N = k.

    The **single-body class** consists of pairs of the form (Iσ, σI) for
    σ ∈ {X, Y, Z}, where the bond bilinear σ·I + I·σ reduces to a sum of
    single-site Pauli terms. For σ ∈ {X, Y} the H-commutator part of M
    cancels (Π σ Π⁻¹ = −σ, so [ΠHΠ⁻¹, ·] + [H, ·] = 0); only the
    dissipator + trace remainder contributes, giving the modified law
    4·(2k − 1)/(2k − 3) (slightly above the main class, converging to 4
    from above).

    The IZ+ZI case lands in **hard** rather than soft (Π Z Π⁻¹ = +Z, so
    H-commutator does not cancel), but its residual still respects the
    same single-body algebraic structure.

    Verified at N=3 → N=4 and N=4 → N=5 across all 120 unordered
    Pauli-pair chain Hamiltonians at γ = 0.1, J = 1.0; see
    experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md.

    Args:
        N1: source chain length (>= 2, i.e. >= 1 bond)
        N2: target chain length (must equal N1 + 1)
        hamiltonian_class: 'main' or 'single_body'

    Returns:
        Predicted ratio of squared Frobenius norms.
    """
    if N1 < 2:
        raise ValueError(f"N1 must be >= 2 (chain needs at least one bond); got {N1}")
    if N2 != N1 + 1:
        raise NotImplementedError(
            f"closed-form scaling verified for adjacent N (N2 = N1 + 1) only; "
            f"got N1={N1}, N2={N2}"
        )
    k = N1
    if hamiltonian_class == 'main':
        return 4.0 * k / (k - 1)
    if hamiltonian_class == 'single_body':
        if 2 * k - 3 <= 0:
            raise ValueError(
                f"single-body formula requires N1 >= 2; got N1={N1} "
                f"(denominator 2N1 − 3 = {2*k - 3})"
            )
        return 4.0 * (2 * k - 1) / (2 * k - 3)
    raise ValueError(
        f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}"
    )


# ----------------------------------------------------------------------
# Section 7: ur-Eigenvalues from Pauli closure
# ----------------------------------------------------------------------

def pauli_dot_pauli_identity():
    """Return the (σ·σ)² identity as a string description.

    The identity (σ_1·σ_2)² = 3I − 2(σ_1·σ_2) follows from:
    1. σ²_a = I for each Pauli (σ²=I is the C=1/2 axiom realized).
    2. {σ_a, σ_b} = 2δ_{ab}I (Clifford anticommutation).
    3. ε_{abc}ε_{abd} = 2δ_{cd} when summed over a, b.

    Combining: (σ·σ)² = 3I + Σ_{a≠b} cross terms = 3I − 2(σ·σ).
    """
    return "(σ_1·σ_2)² = 3I − 2(σ_1·σ_2),  i.e., x² + 2x − 3 = 0,  x ∈ {+1, −3}"


def pauli_dot_pauli_eigenvalues():
    """Return eigenvalues of σ_1·σ_2: +1 (triplet, 3-fold) and −3 (singlet, 1-fold).

    Derived analytically from x² + 2x − 3 = (x − 1)(x + 3) = 0.
    Multiplicities from trace(σ_1·σ_2) = 0: 3·(+1) + 1·(−3) = 0. ✓
    """
    return [(+1.0, 3, 'triplet'), (-3.0, 1, 'singlet')]


def heisenberg_pair_energies(J=1.0):
    """Energy levels of a single Heisenberg pair H = J σ_1·σ_2.

    Triplet: +J (3-fold). Singlet: −3J (1-fold). Splitting 4J.
    """
    return [(J, 3, 'triplet'), (-3 * J, 1, 'singlet')]


def heisenberg_singlet_triplet_gap(J=1.0):
    """The 4J splitting that controls all Level-1 chemistry."""
    return 4 * J


# ----------------------------------------------------------------------
# Section 8: ur-states (singlet, triplet, X-Néel, etc.)
# ----------------------------------------------------------------------

def ur_singlet_2qubit():
    """The 2-qubit singlet state (|01⟩ − |10⟩)/√2."""
    return np.array([0, 1, -1, 0], dtype=complex) / math.sqrt(2)


def ur_triplet_2qubit():
    """The three triplet states: |T_+⟩=|00⟩, |T_0⟩=(|01⟩+|10⟩)/√2, |T_−⟩=|11⟩."""
    T_plus = np.array([1, 0, 0, 0], dtype=complex)
    T_zero = np.array([0, 1, 1, 0], dtype=complex) / math.sqrt(2)
    T_minus = np.array([0, 0, 0, 1], dtype=complex)
    return T_plus, T_zero, T_minus


def ur_xneel(N, sign_pattern=None):
    """X-basis Néel state ⊗_l |s_l⟩ where s_l ∈ {+, −}.

    Default: |+−+−...⟩ (alternating from +). Use sign_pattern to override.
    """
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    if sign_pattern is None:
        sign_pattern = [+1 if l % 2 == 0 else -1 for l in range(N)]
    state = (plus if sign_pattern[0] > 0 else minus).copy()
    for s in sign_pattern[1:]:
        state = np.kron(state, plus if s > 0 else minus)
    return state


# ----------------------------------------------------------------------
# Section 9: V-Effect bridge analysis
# ----------------------------------------------------------------------

def v_effect_emergent_exchange(alpha, J_intra=1.0):
    """Predicted Level-1 effective exchange J_eff from V-Effect bridge α.

    Two singlet pairs each with intra-coupling J_intra, bonded by Heisenberg
    bridge α. Second-order PT gives:
        δE_GS = −(3/8) α² / J_intra
    The 3/8 prefactor: Pauli identity gives ⟨(σ·σ)²⟩_singlet-singlet = 3
    (since ⟨σ·σ⟩ on bridge = 0 between two singlets), divided by gap 8J_intra
    (both pairs flip singlet→triplet costing 4J each).
    """
    return -(3.0 / 8.0) * alpha ** 2 / J_intra


# ----------------------------------------------------------------------
# Section 10: Π-protected observables  (operator → state-level bridge)
# ----------------------------------------------------------------------

def pauli_basis_vector(rho, N):
    """Express a 2^N × 2^N density matrix as a 4^N Pauli-basis coefficient vector.

    Returns vec[k] = (1/2^N) · Tr(σ_k · ρ), so that ρ = (1/2^N) · Σ_k vec[k] · σ_k.
    The Pauli-string indexing matches `_k_to_indices` / `_indices_to_k`.
    """
    d2 = 4 ** N
    vec = np.zeros(d2, dtype=complex)
    for k in range(d2):
        sigma_k = pauli_string(list(_k_to_indices(k, N)))
        vec[k] = np.trace(sigma_k @ rho) / (2 ** N)
    return vec


def _pauli_label(k, N):
    """Convert flat Pauli index k to label string like 'XIZ' (left-most-first)."""
    return ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(k, N))


def pi_protected_observables(H, gamma_l, rho_0, N, threshold=1e-9, cluster_tol=1e-8):
    """Identify which Pauli-string observables σ_α have ⟨σ_α(t)⟩ = 0 for all t,
    given an initial state `rho_0` under Lindblad evolution L = -i[H, ·] +
    Σ γ_l (Z_l ρ Z_l − ρ).

    Algebraic content (no fit, no time-evolution):

        ⟨σ_α(t)⟩ = 2^N · Σ_k V[α, k] · c[k] · exp(λ_k t)
                 = 2^N · Σ_λ S_λ(α) · exp(λ t)

    where (λ_k, V[:,k]) are eigenvalues / right-eigenvectors of L_pauli,
    c = V⁻¹ · ρ_0_pauli are the left-projections of the initial state,
    and S_λ(α) = Σ_{k : λ_k = λ} V[α, k] · c[k] is the total contribution
    at decay rate λ (cluster of degenerate eigenvalues counts once).

    Sums of exponentials with distinct rates are identically zero iff the
    coefficient at each rate vanishes:

        σ_α is **Π-protected** iff S_λ(α) = 0 for every cluster λ.

    This is strictly weaker than "each individual V[α, k]·c[k] vanishes":
    contributions can cancel within a degenerate cluster (Heisenberg
    + Z-dephasing on |+−+⟩: ⟨X₀IZ₂⟩ = 0 for all t, despite individual
    eigenmodes contributing — they cancel pairwise within SU(2)
    multiplets).

    Returns:
        {
          'protected': list of {'k', 'pauli', 'max_cluster_contribution'},
          'active':    list of {'k', 'pauli', 'max_cluster_contribution',
                                'dominant_eigenvalue'},
          'eigenvalues': L_pauli eigenvalues,
          'n_clusters': number of distinct eigenvalue clusters.
        }

    Identity (k=0) is excluded; ⟨I⟩ = 1 trivially.
    """
    L = lindbladian_z_dephasing(H, gamma_l)
    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)

    rho_pauli = pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli

    # Cluster eigenvalues by approximate equality (degeneracies cancel inside)
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cluster = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cluster.append(j)
                used[j] = True
        clusters.append(cluster)

    d2 = 4 ** N
    protected, active = [], []
    for alpha in range(1, d2):  # skip identity
        # Sum V[α, k] · c[k] within each degenerate cluster
        cluster_sums = []
        for cluster in clusters:
            S = sum(V[alpha, k] * c[k] for k in cluster)
            cluster_sums.append((S, cluster[0]))
        max_S = max((abs(S) for S, _ in cluster_sums), default=0.0)
        entry = {
            'k': alpha,
            'pauli': _pauli_label(alpha, N),
            'max_cluster_contribution': float(max_S),
        }
        if max_S < threshold:
            protected.append(entry)
        else:
            # Largest-cluster eigenvalue (the dominant rate of departure from 0)
            dom_S, dom_idx = max(cluster_sums, key=lambda x: abs(x[0]))
            entry['dominant_eigenvalue'] = complex(evals[dom_idx])
            active.append(entry)
    return {'protected': protected, 'active': active,
            'eigenvalues': evals, 'n_clusters': len(clusters)}


# ----------------------------------------------------------------------
# Section 11: Cockpit panel — Lebensader skeleton + cusp-pattern classifier
# ----------------------------------------------------------------------
#
# Adds two predictive structures to the Decoherence Cockpit's three measured
# observables (Purity, Concurrence, L1-coherence). The Cockpit measures what
# the trajectory does; these primitives predict, from the framework's
# operator structure, what the trajectory WILL do before running it.
#
#   pi_lebensader_panel:    skeleton (Π-protected count under pure-Z and +T1)
#                           + trace (θ-trajectory metrics on the full state)
#                           + Lebensader rating combining the two.
#
#   cusp_pattern_classifier: number of CΨ=1/4 crossings, dominant L-eigenmode
#                            at the first crossing, classification
#                            (monotonic / heartbeat / factorising).
#
# Both share the eigendecomposition of L (in Pauli basis) so they're
# combined into one dispatch, `cockpit_panel()`, that returns both views.


def _purity_psi_norm_cpsi(rho):
    """Tom's CΨ glossary: Purity·Ψ-norm. Ψ-norm = L1(ρ)/(d-1)."""
    p = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi_n = l1 / (d - 1)
    return p, psi_n, p * psi_n


def _theta_from_cpsi(c):
    if c <= 0.25:
        return 0.0
    return math.degrees(math.atan(math.sqrt(4 * c - 1)))


def _cluster_eigenvalues(evals, tol=1e-8):
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)
    return clusters


def _trajectory_via_eigendecomp(L_pauli, V, Vinv, rho_0, N, times):
    """ρ(t) at sample times, via the precomputed L eigendecomposition.
    Avoids re-running expm at each time. Returns list of full d×d matrices."""
    evals = np.diag(Vinv @ L_pauli @ V).copy()  # = eigenvalues
    M_basis = _vec_to_pauli_basis_transform(N)
    rho_pauli_0 = pauli_basis_vector(rho_0, N)
    c0 = Vinv @ rho_pauli_0
    out = []
    d = 2 ** N
    for t in times:
        rho_pauli_t = V @ (c0 * np.exp(evals * t))
        rho_vec = M_basis @ rho_pauli_t
        out.append(rho_vec.reshape(d, d).T)
    return out


def _find_cpsi_crossings(times, cpsi_t, threshold=0.25):
    """Linear-interpolate CΨ - threshold zeros. Returns list of (t_cross, sign_change)."""
    out = []
    for i in range(len(times) - 1):
        a, b = cpsi_t[i] - threshold, cpsi_t[i + 1] - threshold
        if a * b < 0:
            frac = a / (a - b)
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            out.append((float(t_cross), -1 if a > 0 else +1, i))
    return out


def cockpit_panel(H, gamma_l, rho_0, N,
                  gamma_t1_l=None,
                  t_max=10.0, dt=0.005,
                  threshold=1e-9, cluster_tol=1e-8):
    """Cockpit predictive panel: Lebensader (skeleton + trace) + cusp pattern.

    Returns a dict with two top-level keys:

      'lebensader':  skeleton (Π-protected counts under pure-Z and +T1),
                     trace (θ-trajectory metrics: max, tail duration, α
                     descent exponent), and a qualitative rating.
      'cusp':        number of CΨ=1/4 crossings, dominant L-eigenvalue
                     at first crossing, classification.

    The two views share the L eigendecomposition (computed once per channel).
    `gamma_t1_l = None` disables T1; both pure-Z and +T1 panels are returned
    when γ_T1 ≠ 0 is supplied.
    """
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N
    times = np.linspace(0.0, t_max, int(t_max / dt) + 1)

    # ---- Skeleton: protected counts under pure-Z and +T1 ----
    pi_pure = pi_protected_observables(
        H, gamma_l, rho_0, N, threshold=threshold, cluster_tol=cluster_tol,
    )
    n_pure = len(pi_pure['protected'])

    if any(g != 0 for g in gamma_t1_l):
        L_t1 = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
        # Replicate pi_protected_observables logic for L_t1 directly
        M_basis = _vec_to_pauli_basis_transform(N)
        L_t1_pauli = (M_basis.conj().T @ L_t1 @ M_basis) / (2 ** N)
        evals_t1, V_t1 = np.linalg.eig(L_t1_pauli)
        Vinv_t1 = np.linalg.inv(V_t1)
        c_t1 = Vinv_t1 @ pauli_basis_vector(rho_0, N)
        clusters_t1 = _cluster_eigenvalues(evals_t1, tol=cluster_tol)
        n_t1 = 0
        for alpha in range(1, 4 ** N):
            max_S = 0.0
            for cl in clusters_t1:
                S = sum(V_t1[alpha, k] * c_t1[k] for k in cl)
                max_S = max(max_S, abs(S))
            if max_S < threshold:
                n_t1 += 1
    else:
        n_t1 = n_pure

    skeleton = {'n_protected_pure': n_pure, 'n_protected_t1': n_t1,
                'drop': n_pure - n_t1}

    # ---- Trace: θ trajectory on full state under L_t1 (if T1) else L_pure ----
    if any(g != 0 for g in gamma_t1_l):
        L_active = L_t1
        M_basis = _vec_to_pauli_basis_transform(N)
        L_active_pauli = (M_basis.conj().T @ L_active @ M_basis) / (2 ** N)
        evals_active, V_active = np.linalg.eig(L_active_pauli)
        Vinv_active = np.linalg.inv(V_active)
    else:
        L_active = lindbladian_z_dephasing(H, gamma_l)
        M_basis = _vec_to_pauli_basis_transform(N)
        L_active_pauli = (M_basis.conj().T @ L_active @ M_basis) / (2 ** N)
        evals_active, V_active = np.linalg.eig(L_active_pauli)
        Vinv_active = np.linalg.inv(V_active)

    traj = _trajectory_via_eigendecomp(L_active_pauli, V_active, Vinv_active,
                                        rho_0, N, times)
    cpsi_t = np.array([_purity_psi_norm_cpsi(r)[2] for r in traj])
    theta_t = np.array([_theta_from_cpsi(c) for c in cpsi_t])

    crossings = _find_cpsi_crossings(times, cpsi_t)
    n_crossings = len(crossings)
    theta_max = float(theta_t.max())

    # Last permanent boundary departure: find last time θ > 0
    above = theta_t > 0
    if above.any():
        last_above_idx = len(times) - 1 - int(np.argmax(above[::-1]))
        t_last_above = float(times[last_above_idx])
    else:
        t_last_above = 0.0

    # Tail duration: time θ spends in (0°, 5°) regime (the "letzte Erinnerung")
    tail_mask = (theta_t > 0) & (theta_t < 5.0)
    if tail_mask.any():
        tail_idx = np.where(tail_mask)[0]
        tail_duration = float(times[tail_idx[-1]] - times[tail_idx[0]])
    else:
        tail_duration = 0.0

    # α descent exponent: fit log θ ~ α log(t* - t) on the last 0.3 unit
    # before the LAST crossing (most informative cusp).
    alpha = None
    if crossings:
        t_cross_last = crossings[-1][0]
        win_mask = ((times >= t_cross_last - 0.3) & (times < t_cross_last)
                    & (theta_t > 0.01) & (theta_t < 10.0))
        if win_mask.sum() >= 4:
            x = np.log(np.maximum(t_cross_last - times[win_mask], 1e-9))
            y = np.log(theta_t[win_mask])
            coef = np.polyfit(x, y, 1)
            alpha = float(coef[0])

    # ---- Cusp-pattern classifier ----
    # Two orthogonal axes: pattern (n_crossings) and dominant-mode TYPE.
    # The combined label is (pattern, mode-type) so users can read both.
    if not crossings:
        pattern = 'never crosses'
        dom_eigval = None
        dom_eigval_nonzero = None
        mode_type = 'n/a'
    else:
        t_first, _, idx_first = crossings[0]
        rho_at_cross = traj[idx_first]
        rho_pauli_at = pauli_basis_vector(rho_at_cross, N)
        c_at = Vinv_active @ rho_pauli_at
        clusters = _cluster_eigenvalues(evals_active, tol=cluster_tol)
        cluster_norms = []
        for cl in clusters:
            norm_cl = float(np.linalg.norm([c_at[k] for k in cl]))
            cluster_norms.append((norm_cl, cl))
        cluster_norms.sort(key=lambda x: -x[0])
        dom_eigval = complex(evals_active[cluster_norms[0][1][0]])
        # Dominant NON-ZERO eigenvalue cluster (excludes steady-state if it's
        # the absolute largest — the actual driving mode of the dynamics).
        nonzero_norms = [(n, cl) for n, cl in cluster_norms
                          if abs(complex(evals_active[cl[0]])) > 1e-6]
        if nonzero_norms:
            dom_eigval_nonzero = complex(evals_active[nonzero_norms[0][1][0]])
        else:
            dom_eigval_nonzero = None

        pattern = 'monotonic' if n_crossings == 1 else 'heartbeat'

        # Classify the dominant-mode TYPE (uses the non-zero mode if steady
        # state is the absolute largest, otherwise the absolute largest).
        ref_eigval = (dom_eigval_nonzero if abs(dom_eigval) < 1e-6
                       else dom_eigval)
        if ref_eigval is None:
            mode_type = 'pure steady state'
        elif abs(dom_eigval) < 1e-6:
            # Steady state is the largest cluster → "settling" regime
            if abs(ref_eigval.imag) > 1e-3:
                mode_type = 'steady-state + oscillatory sub-mode'
            else:
                mode_type = 'steady-state + real-decay sub-mode'
        elif abs(ref_eigval.imag) > 1e-3:
            mode_type = 'oscillatory'
        else:
            mode_type = 'real-decay'

    classification = (f"{pattern} | {mode_type}"
                      if pattern != 'never crosses' else 'never crosses')

    trace = {'theta_max': theta_max,
             'n_crossings': n_crossings,
             't_last_theta_positive': t_last_above,
             'tail_duration_sub5deg': tail_duration,
             'alpha_descent': alpha}

    cusp = {'n_crossings': n_crossings,
            'dominant_eigenvalue_at_first_crossing': dom_eigval,
            'dominant_nonzero_eigenvalue': dom_eigval_nonzero,
            'pattern': pattern,
            'mode_type': mode_type,
            'classification': classification}

    # Lebensader rating: report skeleton and trace as SEPARATE statuses,
    # plus a combined reading. They decouple at N ≥ 4, so the user should
    # read both axes independently.
    skel_status = 'intact' if skeleton['drop'] <= 1 else (
        'partial' if skeleton['drop'] <= 5 else 'broken'
    )
    # Trace threshold: N-scaled. Empirically max tail among skeleton-intact
    # cases is ~0.86 at N=3, 0.010 at N=4, 0.005 at N=5. A threshold of
    # max(0.05, 0.86 / 2^(N-3)) catches the same proportion at each N.
    tail_threshold = max(0.005, 0.86 / (2 ** max(0, N - 3)))
    trace_status = 'long' if tail_duration > tail_threshold else (
        'short' if tail_duration > tail_threshold / 10 else 'absent'
    )
    if skel_status == 'intact' and trace_status == 'long':
        leb_rating = 'intact (skeleton holds, trace lives)'
    elif skel_status == 'intact' and trace_status != 'long':
        leb_rating = 'partial (skeleton holds, trace short)'
    elif skel_status != 'intact' and trace_status == 'long':
        leb_rating = 'partial (skeleton bleeds, trace persists)'
    else:
        leb_rating = 'collapsed (both skeleton and trace gone)'

    # ---- Chiral panel: K_full sublattice symmetry ----
    chiral = chiral_panel(H, rho_0, N)

    # ---- Y-parity panel: bit_a · bit_b grading ----
    y_parity = y_parity_panel(H, gamma_l, rho_0, N, gamma_t1_l=gamma_t1_l)

    return {
        'lebensader': {'skeleton': skeleton, 'trace': trace,
                       'rating': leb_rating,
                       'skeleton_status': skel_status,
                       'trace_status': trace_status,
                       'tail_threshold_for_N': tail_threshold},
        'cusp': cusp,
        'chiral': chiral,
        'y_parity': y_parity,
        '_trajectory_for_inspection': {'times': times, 'cpsi': cpsi_t,
                                        'theta': theta_t},
    }


# ----------------------------------------------------------------------
# Section 12: Chiral panel — K_full = ⊗_{odd i} Z_i sublattice symmetry
# ----------------------------------------------------------------------
#
# K is the chiral / sublattice symmetry of the XY chain (AZ class BDI).
# K_full = ⊗_{odd i} Z_i (Z on odd-indexed sites) anti-commutes with the
# bilinear bond Hamiltonian H_xy = (X_a X_{a+1} + Y_a Y_{a+1})/2:
#     K_full · H_xy · K_full = −H_xy
# Z-dephasing commutes with K_full trivially (Z operators on any site
# commute with each Z on odd sites). So when K_full anti-commutes with
# H, the full Lindbladian L = −i[H, ·] + L_dephasing has K as a
# super-operator symmetry: K · L[ρ] · K = L[K ρ K]. This block-
# diagonalises L into K-symmetric and K-antisymmetric sectors of ρ.
#
# The chiral panel reports:
#   - K_full structure
#   - K-symmetry status of H ([K,H] = 0 or {K,H} = 0 or neither)
#   - Pauli classification: each P_α is K-even (K P K = +P) or K-odd
#   - ρ_0's K-decomposition: ρ_+ + ρ_- with weights w_+, w_- = Tr(ρ_±²)
#   - chiral skeleton: Pauli observables that are K-protected
#     (⟨P⟩(0) = 0 by K-mismatch and stay zero forever if K is a symmetry)


def chiral_K_full(N):
    """K_full = ⊗_{odd i} Z_i, the chain chiral / sublattice operator."""
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    out = np.array([[1.0]], dtype=complex)
    for i in range(N):
        out = np.kron(out, Z if (i % 2 == 1) else I2)
    return out


def k_classify_pauli(N, K_full=None):
    """For each Pauli string α, return its K-eigenvalue (+1 or −1).

    Returns a dict {alpha_int: ±1} for α = 0..4^N-1, plus a list of
    K-even labels and K-odd labels.
    """
    if K_full is None:
        K_full = chiral_K_full(N)
    out = {}
    for alpha in range(4 ** N):
        idx_tuple = _k_to_indices(alpha, N)
        P = np.array([[1.0]], dtype=complex)
        for idx in idx_tuple:
            P = np.kron(P, _PAULI_MATRICES[idx])
        KPK = K_full @ P @ K_full
        if np.allclose(KPK, P, atol=1e-10):
            out[alpha] = +1
        elif np.allclose(KPK, -P, atol=1e-10):
            out[alpha] = -1
        else:
            out[alpha] = 0  # mixed (shouldn't happen for K_full = ⊗Z's)
    return out


def k_classify_hamiltonian(H, N, K_full=None):
    """Test whether H is K-even ([K,H]=0), K-odd ({K,H}=0), or mixed.

    Returns 'K-even', 'K-odd', or 'K-mixed'.
    """
    if K_full is None:
        K_full = chiral_K_full(N)
    KHK = K_full @ H @ K_full
    if np.allclose(KHK, H, atol=1e-10):
        return 'K-even'
    if np.allclose(KHK, -H, atol=1e-10):
        return 'K-odd'
    return 'K-mixed'


def chiral_panel(H, rho_0, N, K_full=None):
    """Chiral structure panel for (H, ρ_0). Reports K-symmetry status of L
    and the K-character of the initial state — a STATIC structural reading,
    not dynamical Pauli-protection.

    K dynamical protection of ⟨P⟩(t) = 0 requires:
      (1) K is a symmetry of L (K-status of H in {'K-even', 'K-odd'})
      (2) ρ_0 is in a K-eigenstate (w_+ = 0 OR w_- = 0)

    Both conditions are reported. The 'chiral_protected_dynamic' set is
    populated only when both hold.

    Separately, 'static_zero_at_t0' reports Pauli observables with
    ⟨P_α⟩(0) = 0 by virtue of the (K-class of P, K-projection of ρ_0)
    contraction. These are NOT guaranteed to stay zero under L unless
    the dynamic conditions also hold.

    Returns dict with:
      'K_status':        'K-even' / 'K-odd' / 'K-mixed'
      'K_symmetric_L':   bool — is K a symmetry of L?
                          True iff K_status ∈ {'K-even', 'K-odd'}
      'rho_decomposition': {'w_plus': Tr(ρ_+²), 'w_minus': Tr(ρ_-²),
                             'rho_plus': ρ_+, 'rho_minus': ρ_-}
      'rho_is_K_eigenstate': bool — w_+ < 1e-12 OR w_- < 1e-12
      'pauli_classification': {'n_K_even': count, 'n_K_odd': count}
      'static_zero_at_t0':  Pauli labels with ⟨P⟩(0) = 0 from K-structure
                             (does NOT imply dynamical protection)
      'chiral_protected_dynamic': Pauli labels GUARANTEED to stay zero
                                    (only nonempty when both conditions hold)
    """
    if K_full is None:
        K_full = chiral_K_full(N)

    K_status = k_classify_hamiltonian(H, N, K_full=K_full)
    K_symmetric_L = (K_status in ('K-even', 'K-odd'))

    rho_KKt = K_full @ rho_0 @ K_full
    rho_plus = (rho_0 + rho_KKt) / 2
    rho_minus = (rho_0 - rho_KKt) / 2
    w_plus = float(np.real(np.trace(rho_plus @ rho_plus)))
    w_minus = float(np.real(np.trace(rho_minus @ rho_minus)))
    rho_is_K_eigenstate = (w_plus < 1e-12) or (w_minus < 1e-12)

    k_class = k_classify_pauli(N, K_full=K_full)
    n_even = sum(1 for v in k_class.values() if v == +1)
    n_odd = sum(1 for v in k_class.values() if v == -1)

    static_zero = []
    for alpha in range(1, 4 ** N):
        # ⟨P_α⟩(0) = Tr(P_α ρ_0). If P is K-even, only ρ_+ contributes;
        # if K-odd, only ρ_-.
        P_alpha_idx = _k_to_indices(alpha, N)
        P_alpha = np.array([[1.0]], dtype=complex)
        for idx in P_alpha_idx:
            P_alpha = np.kron(P_alpha, _PAULI_MATRICES[idx])
        if k_class[alpha] == +1:
            exp_val = float(np.real(np.trace(P_alpha @ rho_plus)))
        elif k_class[alpha] == -1:
            exp_val = float(np.real(np.trace(P_alpha @ rho_minus)))
        else:
            exp_val = 1.0  # K-mixed P shouldn't appear at all
        if abs(exp_val) < 1e-10:
            static_zero.append(alpha)

    # Dynamical protection: only when BOTH conditions hold.
    chiral_protected_dynamic = []
    if K_symmetric_L and rho_is_K_eigenstate:
        for alpha in range(1, 4 ** N):
            if (k_class[alpha] == +1 and w_minus < 1e-12) or \
               (k_class[alpha] == -1 and w_plus < 1e-12):
                chiral_protected_dynamic.append(alpha)

    static_zero_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(a, N))
        for a in static_zero
    ]
    chiral_dynamic_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(a, N))
        for a in chiral_protected_dynamic
    ]

    return {
        'K_status': K_status,
        'K_symmetric_L': K_symmetric_L,
        'rho_decomposition': {
            'w_plus': w_plus, 'w_minus': w_minus,
            'rho_plus': rho_plus, 'rho_minus': rho_minus,
        },
        'rho_is_K_eigenstate': rho_is_K_eigenstate,
        'pauli_classification': {
            'n_K_even': n_even, 'n_K_odd': n_odd,
        },
        'static_zero_at_t0': static_zero_labels,
        'chiral_protected_dynamic': chiral_dynamic_labels,
    }


# ----------------------------------------------------------------------
# Section 13: Y-parity panel — Z₂ grading on Pauli-Y count
# ----------------------------------------------------------------------
#
# Y-parity = parity of the number of Y operators in a Pauli string.
# Equivalently: parity of sum over sites of (bit_a · bit_b) since Y is
# the unique Pauli with both bit_a = 1 (decaying) AND bit_b = 1 (Π²-odd).
#
# Unlike Π and K, Y-parity is NOT a conjugation symmetry of the Pauli
# algebra (no Pauli operator U has U P U = (-1)^{#Y(P)} P). Pauli
# multiplication does not preserve Y-count in general (e.g., X·Z = ±iY
# adds a Y).
#
# Y-parity preservation is therefore a system-specific algebraic
# property: it requires that L (as a super-operator) maps Y-parity-even
# operators to Y-parity-even operators, equivalently that L is
# block-diagonal in the Y-parity decomposition of the Pauli basis.
#
# This panel checks this empirically by projecting L_pauli onto the
# Y-even/Y-odd Pauli subspaces and measuring the off-diagonal weight.
#
# Empirical finding (from cockpit-120 enum + intact_hards_analysis):
# at N=3, |+−+⟩, the 6 Lebensader-intact cases (XY+YX, IY+YI, XY+YZ,
# XY+ZY, YX+YZ, YX+ZY) all have L preserving Y-parity, and ρ_0 = |+−+⟩
# has Y-parity-0 (no Y components). So all 28 Y-odd observables stay
# strictly zero forever for these 6 cases.


def _y_parity_classify_paulis(N):
    """Return Y-parity (0 or 1) for each Pauli index α = 0..4^N-1."""
    out = np.zeros(4 ** N, dtype=int)
    for alpha in range(4 ** N):
        idx_tuple = _k_to_indices(alpha, N)
        n_y = sum(1 for idx in idx_tuple if idx == (1, 1))  # Y = (1, 1)
        out[alpha] = n_y % 2
    return out


def y_parity_panel(H, gamma_l, rho_0, N, gamma_t1_l=None):
    """Y-parity panel: checks whether L preserves Y-parity, and reports
    ρ_0's Y-parity decomposition + the protected Y-odd Pauli set.

    Returns dict with:
      'L_preserves_Y_parity':  bool — is L block-diagonal in Y-parity?
      'L_offdiag_weight':      Frobenius norm of off-diagonal Y-parity blocks
                                (in the Pauli basis). Small = Y-parity preserved.
      'rho_0_Y_decomposition': {'w_even': ..., 'w_odd': ...,
                                 'rho_even_pauli_norm': ...}
                                Pauli-vector norms of Y-even and Y-odd parts
                                of ρ_0 in the Pauli basis.
      'pauli_classification':  {'n_Y_even': count, 'n_Y_odd': count}
      'Y_odd_observables':     List of Pauli labels in the Y-odd sector
                                (28 at N=3).
      'Y_parity_protected':    Y-odd Pauli labels that stay zero forever
                                (requires L_preserves_Y_parity = True AND
                                ρ_0 has no Y-odd content).
    """
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N

    # Build L (with T1 if active)
    if any(g != 0 for g in gamma_t1_l):
        L = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = lindbladian_z_dephasing(H, gamma_l)

    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

    y_parity = _y_parity_classify_paulis(N)
    even_idx = np.where(y_parity == 0)[0]
    odd_idx = np.where(y_parity == 1)[0]
    n_even = len(even_idx)
    n_odd = len(odd_idx)

    # Off-diagonal blocks: L_pauli[odd, even] and L_pauli[even, odd].
    # If L preserves Y-parity, both should be ≈ 0.
    L_oe = L_pauli[np.ix_(odd_idx, even_idx)]
    L_eo = L_pauli[np.ix_(even_idx, odd_idx)]
    offdiag_norm = float(np.linalg.norm(L_oe)) + float(np.linalg.norm(L_eo))
    diag_norm = float(np.linalg.norm(L_pauli[np.ix_(even_idx, even_idx)])) \
                 + float(np.linalg.norm(L_pauli[np.ix_(odd_idx, odd_idx)]))
    relative_offdiag = offdiag_norm / (diag_norm + 1e-15)
    L_preserves_Y_parity = (relative_offdiag < 1e-10)

    # ρ_0 Y-parity decomposition (in Pauli basis)
    rho_pauli_vec = pauli_basis_vector(rho_0, N)
    rho_even_pauli = rho_pauli_vec.copy()
    rho_even_pauli[odd_idx] = 0
    rho_odd_pauli = rho_pauli_vec.copy()
    rho_odd_pauli[even_idx] = 0
    w_even = float(np.linalg.norm(rho_even_pauli))
    w_odd = float(np.linalg.norm(rho_odd_pauli))

    # Y-odd Pauli labels
    y_odd_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
        for a in odd_idx
    ]

    # Y-parity-protected: Y-odd Paulis that stay 0 forever
    Y_parity_protected = []
    if L_preserves_Y_parity and w_odd < 1e-10:
        # ρ_0 is Y-even, L preserves Y-parity → all Y-odd ⟨P⟩ = 0 forever
        Y_parity_protected = list(y_odd_labels)
    elif L_preserves_Y_parity and w_even < 1e-10:
        # ρ_0 is Y-odd, L preserves Y-parity → all Y-even ⟨P⟩ = 0 forever
        # (rare case; would protect 27 Y-even non-identity strings + identity)
        even_labels = [
            ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
            for a in even_idx if a != 0
        ]
        Y_parity_protected = even_labels

    return {
        'L_preserves_Y_parity': L_preserves_Y_parity,
        'L_offdiag_weight': offdiag_norm,
        'L_relative_offdiag': relative_offdiag,
        'rho_0_Y_decomposition': {
            'w_even': w_even, 'w_odd': w_odd,
        },
        'pauli_classification': {
            'n_Y_even': int(n_even), 'n_Y_odd': int(n_odd),
        },
        'Y_odd_observables': y_odd_labels,
        'Y_parity_protected': Y_parity_protected,
    }


# ----------------------------------------------------------------------
# Section 14: bit_a / bit_b parity panels — multiplicative Z₂ gradings
# ----------------------------------------------------------------------
#
# The Pauli algebra has TWO underlying Z₂ gradings (multiplicative under
# Pauli multiplication, since (a₁,b₁) · (a₂,b₂) = (a₁⊕a₂, b₁⊕b₂)):
#
#   bit_a-parity (#X + #Y mod 2):  parity of "decaying" sites in P
#   bit_b-parity (#Y + #Z mod 2):  parity of "Π²-odd" sites in P
#
# Their XOR (bit_a XOR bit_b) gives a third multiplicative grading:
#   #X + #Z mod 2  (parity of "Π²-even non-I" sites, excluding Y).
#
# Y-parity (#Y mod 2) is NOT a multiplicative grading (e.g., X·Z = ±iY
# adds a Y) but IS preserved empirically by certain L's (Section 13).
#
# Z-dephasing preserves ALL three multiplicative parities (σ_z has
# bit_a=0, bit_b=1, fixed values; conjugation only flips Pauli signs).
# T1 amplitude damping σ⁻ has bit_a=1 (pure) but mixed bit_b → T1
# preserves bit_a-parity but NOT bit_b-parity in general.
#
# A multiplicative parity is preserved by L iff [H, ·] preserves it,
# which happens iff H has parity-even content in its Pauli decomposition.


def _bit_a_b_classify_paulis(N):
    """Return (bit_a_parity, bit_b_parity) ∈ {0,1}² for each Pauli α."""
    bit_a_arr = np.zeros(4 ** N, dtype=int)
    bit_b_arr = np.zeros(4 ** N, dtype=int)
    for alpha in range(4 ** N):
        idx_tuple = _k_to_indices(alpha, N)
        a_total = sum(idx[0] for idx in idx_tuple) % 2
        b_total = sum(idx[1] for idx in idx_tuple) % 2
        bit_a_arr[alpha] = a_total
        bit_b_arr[alpha] = b_total
    return bit_a_arr, bit_b_arr


def parity_panel(H, gamma_l, rho_0, N, parity_kind, gamma_t1_l=None):
    """Multiplicative Z₂ parity panel: 'bit_a', 'bit_b', or 'bit_a_xor_b'.

    Reports:
      L_preserves_parity:  is L block-diagonal in this parity?
      L_offdiag_weight, L_relative_offdiag
      rho_0_decomposition: w_even, w_odd norms in Pauli basis
      pauli_classification: n_even, n_odd
      protected:  Pauli labels with ⟨P⟩ guaranteed zero (parity-mismatch
                   on a parity-eigenstate ρ_0)
    """
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N
    if any(g != 0 for g in gamma_t1_l):
        L = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = lindbladian_z_dephasing(H, gamma_l)

    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

    bit_a_arr, bit_b_arr = _bit_a_b_classify_paulis(N)
    if parity_kind == 'bit_a':
        parity = bit_a_arr
    elif parity_kind == 'bit_b':
        parity = bit_b_arr
    elif parity_kind == 'bit_a_xor_b':
        parity = (bit_a_arr ^ bit_b_arr).astype(int)
    else:
        raise ValueError(f"Unknown parity_kind: {parity_kind}")

    even_idx = np.where(parity == 0)[0]
    odd_idx = np.where(parity == 1)[0]

    L_oe = L_pauli[np.ix_(odd_idx, even_idx)]
    L_eo = L_pauli[np.ix_(even_idx, odd_idx)]
    offdiag_norm = float(np.linalg.norm(L_oe)) + float(np.linalg.norm(L_eo))
    diag_norm = float(np.linalg.norm(L_pauli[np.ix_(even_idx, even_idx)])) \
                 + float(np.linalg.norm(L_pauli[np.ix_(odd_idx, odd_idx)]))
    relative_offdiag = offdiag_norm / (diag_norm + 1e-15)
    L_preserves = (relative_offdiag < 1e-10)

    rho_pauli_vec = pauli_basis_vector(rho_0, N)
    rho_even_vec = rho_pauli_vec.copy()
    rho_even_vec[odd_idx] = 0
    rho_odd_vec = rho_pauli_vec.copy()
    rho_odd_vec[even_idx] = 0
    w_even = float(np.linalg.norm(rho_even_vec))
    w_odd = float(np.linalg.norm(rho_odd_vec))

    odd_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
        for a in odd_idx
    ]

    protected = []
    if L_preserves and w_odd < 1e-10:
        protected = list(odd_labels)
    elif L_preserves and w_even < 1e-10:
        even_labels = [
            ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
            for a in even_idx if a != 0
        ]
        protected = even_labels

    return {
        'parity_kind': parity_kind,
        'L_preserves_parity': L_preserves,
        'L_offdiag_weight': offdiag_norm,
        'L_relative_offdiag': relative_offdiag,
        'rho_0_decomposition': {'w_even': w_even, 'w_odd': w_odd},
        'pauli_classification': {
            'n_even': int(len(even_idx)), 'n_odd': int(len(odd_idx)),
        },
        'protected': protected,
    }


# ----------------------------------------------------------------------
# Section 15: recommend_initial_state — best ρ_0 for given H
# ----------------------------------------------------------------------
#
# Given a Hamiltonian H + dissipator (γ_dephasing, γ_T1), test a catalog
# of standard initial states and return the one that activates the most
# Z₂ / U(1) protections.
#
# Protection counts come from pi_protected_observables under L_full:
# they include cluster cancellation from Y-parity, bit_a, K_full, and
# any other symmetries that happen to apply for the given (H, ρ_0).
#
# The recommender doesn't classify which symmetry kicks in (use the
# 4-panel cockpit for that diagnosis); it just returns the maximum-
# protection initial state from a catalog of computational and
# superposition states. All catalog members are operationally trivial
# to prepare on a real quantum device.


def _compute_n_protected(H, gamma_l, gamma_t1_l, rho_0, N,
                          threshold=1e-9, cluster_tol=1e-8):
    """Count pi_protected Pauli observables under L_full = L_dephasing + L_T1."""
    if any(g != 0 for g in gamma_t1_l):
        L = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = lindbladian_z_dephasing(H, gamma_l)
    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)
    rho_pauli = pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli

    n_eig = len(evals)
    used = np.zeros(n_eig, dtype=bool)
    clusters = []
    for i in range(n_eig):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n_eig):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    n_protected = 0
    for alpha in range(1, 4 ** N):
        max_S = max(
            abs(sum(V[alpha, k] * c[k] for k in cl)) for cl in clusters
        )
        if max_S < threshold:
            n_protected += 1
    return n_protected


def standard_initial_state_catalog(N):
    """Catalog of ρ_0 candidates for the recommender.

    Returns a list of (label, ρ_0_matrix). All candidates are pure states
    that can be prepared on a real quantum device with a constant number
    of single-qubit + nearest-neighbour gates.

    Includes:
      - All 2^N computational basis states |b_{N-1}…b_1 b_0⟩
      - |+⟩^N (all-Hadamard)
      - |+−+−…⟩ (alternating Hadamard + Z)
      - |GHZ⟩ = (|0…0⟩ + |1…1⟩)/√2
      - |W⟩ = symmetric single-excitation
      - (|0⟩ + i|1⟩)^N / 2^(N/2) (Y-eigenstate, bit_b-pure)
    """
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    one = np.array([0, 1], dtype=complex)
    plus_y = np.array([1, 1j], dtype=complex) / math.sqrt(2)  # +Y eigenstate

    out = []
    d = 2 ** N

    # Computational basis states
    for b in range(d):
        psi = np.zeros(d, dtype=complex)
        psi[b] = 1.0
        bits = format(b, f'0{N}b')
        out.append((f"|{bits}⟩", np.outer(psi, psi.conj())))

    # |+⟩^N
    psi_plus = plus.copy()
    for _ in range(N - 1):
        psi_plus = np.kron(psi_plus, plus)
    out.append(("|+⟩^N", np.outer(psi_plus, psi_plus.conj())))

    # |+−+−…⟩
    psi_alt = plus.copy()
    for k in range(1, N):
        psi_alt = np.kron(psi_alt, minus if k % 2 == 1 else plus)
    out.append(("|+−+…⟩", np.outer(psi_alt, psi_alt.conj())))

    # |GHZ⟩
    psi_ghz = np.zeros(d, dtype=complex)
    psi_ghz[0] = psi_ghz[d - 1] = 1.0 / math.sqrt(2)
    out.append(("|GHZ⟩", np.outer(psi_ghz, psi_ghz.conj())))

    # |W⟩ = symmetric single-excitation
    psi_w = np.zeros(d, dtype=complex)
    for k in range(N):
        psi_w[1 << k] = 1.0 / math.sqrt(N)
    out.append(("|W⟩", np.outer(psi_w, psi_w.conj())))

    # Y-eigenstate: (|0⟩+i|1⟩)/√2 on each qubit (bit_b-pure)
    psi_yeig = plus_y.copy()
    for _ in range(N - 1):
        psi_yeig = np.kron(psi_yeig, plus_y)
    out.append(("|+y⟩^N", np.outer(psi_yeig, psi_yeig.conj())))

    # |+y -y +y …⟩ alternating Y-eigenstates
    minus_y = np.array([1, -1j], dtype=complex) / math.sqrt(2)
    psi_yalt = plus_y.copy()
    for k in range(1, N):
        psi_yalt = np.kron(psi_yalt, minus_y if k % 2 == 1 else plus_y)
    out.append(("|+y -y +y…⟩", np.outer(psi_yalt, psi_yalt.conj())))

    return out


def _bit_a_even_fraction(rho_0, N):
    """Fraction of ρ_0's Pauli-basis coefficient norm on bit_a-EVEN sector.

    bit_a-even Pauli operators (composed of I and Z only on each site)
    are immune to Z-dephasing (the framework's 'slow / immune' sector,
    F61). bit_a-odd operators (containing X or Y) decay under Z-dephasing.

    README Section 10 Rule 1: GHZ projects 100% onto bit_a-odd / XOR modes
    → bit_a_even_fraction ≈ 0. W distributes across modes → higher.

    Returns ‖c on bit_a-even Paulis‖ / ‖c on non-identity Paulis‖.
    """
    rho_pauli = pauli_basis_vector(rho_0, N)
    bit_a_arr, _ = _bit_a_b_classify_paulis(N)
    # Exclude identity (alpha=0) since it always has weight 1 trivially.
    mask_even = (bit_a_arr == 0)
    mask_even[0] = False
    mask_total = np.ones(4 ** N, dtype=bool)
    mask_total[0] = False
    norm_even = float(np.linalg.norm(rho_pauli[mask_even]))
    norm_total = float(np.linalg.norm(rho_pauli[mask_total]))
    if norm_total < 1e-12:
        return 1.0  # only identity content; trivially "stable"
    return norm_even / norm_total


def theoretical_max_protected_pure(N):
    """Theoretical maximum n_protected for any pure-state ρ_0 at N qubits.

    For pure ρ_0, Tr(ρ²) = 1 ⇒ Σ_α ⟨P_α⟩² = 2^N. Excluding identity
    (⟨I⟩=1), Σ_{α≠0} ⟨P_α⟩² = 2^N − 1. Since |⟨P⟩| ≤ 1, the minimum
    number of non-zero non-identity Pauli expectations is 2^N − 1
    (achieved iff each non-zero ⟨P⟩ = ±1, i.e., ρ_0 is a STABILIZER STATE).

    Therefore max n_protected = (4^N − 1) − (2^N − 1) = 4^N − 2^N.

    At N=3: 56. At N=4: 240. At N=5: 992.

    Saturation conditions: ρ_0 is a stabilizer state AND the Hamiltonian
    leaves the (N stabilizer-eigenvalue) constraints invariant under
    L = -i[H, ·] + dephasing + T1.

    Computational basis states (|0…0⟩, |010⟩, etc.) are Z-stabilizer
    states: their stabilizer is generated by {Z_0, Z_1, …, Z_{N-1}}.
    They saturate for any Hamiltonian conserving total Z (U(1)).

    |GHZ⟩ is a stabilizer state with stabilizer ⟨X⊗X⊗…⊗X, Z_iZ_j⟩.
    |+⟩^N is the all-X stabilizer.
    """
    return (2 ** N) * (2 ** N) - 2 ** N  # = 4^N - 2^N


def recommend_initial_state(H, gamma_l, N, gamma_t1_l=None,
                              candidates=None, top_k=5,
                              threshold=1e-9, cluster_tol=1e-8,
                              score_weight_slow=0.5):
    """For a given H + dissipator, recommend ρ_0 maximising COMBINED score.

    Combines two metrics:
      n_protected:    number of Pauli observables strictly zero forever
                       (Π-protection count)
      slow_fraction:  fraction of ρ_0's L-eigenmode coefficient norm on
                       slow modes (= long signal lifetime for the
                       non-protected observables)

    README Section 10 Rule 1 says 'use W states, not GHZ': GHZ has
    n_protected high but slow_fraction = 0 (all weight on fast XOR modes,
    so the experimental signal decays before measurement). The combined
    score is:
      score = n_protected / max_n_protected
              + score_weight_slow * slow_fraction

    Default score_weight_slow = 0.5 gives equal weight to (normalised)
    protection count and slow-mode fraction. Set score_weight_slow = 0
    to recover pure-protection ranking.

    Returns dict:
      'best':       {rho_0, label, n_protected, slow_fraction, score}
      'top_k':      list sorted by combined score, with warnings
      'pure_protection_best':  for comparison, the n_protected-only winner
      'catalog_size': int
    """
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N
    if candidates is None:
        candidates = standard_initial_state_catalog(N)

    results = []
    for label, rho_0 in candidates:
        n_prot = _compute_n_protected(
            H, gamma_l, gamma_t1_l, rho_0, N,
            threshold=threshold, cluster_tol=cluster_tol,
        )
        slow_frac = _bit_a_even_fraction(rho_0, N)
        # Detect "trivially classical" states (eigenstate of Z⊗N basis):
        # only I and Z components, no X or Y. Pure diagonals have
        # bit_a_even_fraction = 1.0 AND no off-diagonal coherence.
        rho_pauli = pauli_basis_vector(rho_0, N)
        bit_a_arr, _ = _bit_a_b_classify_paulis(N)
        # Coherence content = norm of bit_a-odd part (anything with X or Y)
        bit_a_odd_norm = float(np.linalg.norm(rho_pauli[bit_a_arr == 1]))
        is_classical_diagonal = bit_a_odd_norm < 1e-10

        results.append({
            'label': label,
            'n_protected': n_prot,
            'bit_a_even_fraction': slow_frac,
            'is_classical_diagonal': is_classical_diagonal,
            'rho_0': rho_0,
        })

    if results:
        max_n = max(r['n_protected'] for r in results)
        max_n = max(max_n, 1)
        for r in results:
            # Combined score: protection + slow-mode weight - classical penalty
            classical_penalty = 0.5 if r['is_classical_diagonal'] else 0.0
            r['score'] = (r['n_protected'] / max_n
                           + score_weight_slow * r['bit_a_even_fraction']
                           - classical_penalty)

    # Sort by combined score
    by_score = sorted(results, key=lambda r: -r['score'])
    # Sort by pure n_protected (for comparison)
    by_n = sorted(results, key=lambda r: -r['n_protected'])
    # Filter out classical-diagonal, then by score (for "meaningful demo" rec)
    quantum_only = [r for r in results if not r['is_classical_diagonal']]
    by_quantum_score = sorted(quantum_only, key=lambda r: -r['score'])

    best = by_score[0]
    pure_best = by_n[0]
    quantum_best = by_quantum_score[0] if quantum_only else None

    warnings = []
    if pure_best['is_classical_diagonal']:
        warnings.append(
            f"Pure-protection winner is {pure_best['label']}, a classical "
            f"Z-eigenstate (no X/Y coherence). Protection is mostly trivial "
            f"conservation — no quantum demonstration value. Consider "
            f"{quantum_best['label']} for a meaningful experiment."
        )
    # GHZ caveat (README Section 10 Rule 1) — always warn if GHZ in top 3
    ghz_in_top = any('GHZ' in r['label'] for r in by_score[:3])
    if ghz_in_top:
        ghz_record = next((r for r in results if r['label'] == '|GHZ⟩'), None)
        w_record = next((r for r in results if r['label'] == '|W⟩'), None)
        msg = (
            f"|GHZ⟩ in top recommendations. README Section 10 Rule 1: "
            f"'Use W states, not GHZ. GHZ excites only the fastest-absorbing "
            f"modes (all light, maximum absorption). W distributes across modes.' "
            f"For experimentally-meaningful demos, the non-protected "
            f"observables decay slower with W."
        )
        if ghz_record and w_record:
            msg += (f" |GHZ⟩: n_prot={ghz_record['n_protected']}, "
                    f"|W⟩: n_prot={w_record['n_protected']}.")
        warnings.append(msg)

    # Theoretical maximum bound + saturation
    max_pure = theoretical_max_protected_pure(N)
    if pure_best['n_protected'] >= max_pure:
        saturation_status = 'SATURATED — stabilizer-state regime'
    else:
        saturation_frac = pure_best['n_protected'] / max_pure
        saturation_status = (f"{saturation_frac:.0%} of theoretical max "
                              f"({pure_best['n_protected']}/{max_pure})")

    return {
        'best': best,
        'pure_protection_best': pure_best,
        'quantum_best': quantum_best,
        'top_k': by_score[:top_k],
        'top_k_by_protection': by_n[:top_k],
        'catalog_size': len(candidates),
        'warnings': warnings,
        'theoretical_max_pure': max_pure,
        'saturation_status': saturation_status,
    }


# ----------------------------------------------------------------------
# Section 16: Cavity-mode-exposure (F64) — single-excitation coherence
# ----------------------------------------------------------------------

def single_excitation_h1(N, bonds, J=1.0):
    """Single-excitation block H_1 (N×N) for an XX+YY (hopping) Hamiltonian.

    H_1[i, j] = J_{ij} if (i, j) ∈ bonds, else 0. Encodes the |i⟩→|j⟩ amplitude
    in the basis where |i⟩ has a single excitation at site i.

    For Heisenberg/XXZ Hamiltonians, the XX+YY part contributes the off-diagonal
    hopping; the ZZ part contributes a constant shift on the diagonal that does
    not affect eigenvectors. So this primitive — and the F64 formula derived
    from it — apply identically to XY, Heisenberg, XXZ in the single-excitation
    sector.

    `J` may be a scalar (uniform) or a list of length len(bonds) (per-bond).
    """
    H1 = np.zeros((N, N), dtype=complex)
    if np.isscalar(J):
        J_per_bond = [float(J)] * len(bonds)
    else:
        J_per_bond = list(J)
        if len(J_per_bond) != len(bonds):
            raise ValueError("len(J) must match len(bonds) when J is a list")
    for (i, j), J_ij in zip(bonds, J_per_bond):
        H1[i, j] += J_ij
        H1[j, i] += J_ij
    return H1


def cavity_coh_liouvillian(H1, B, gamma_B):
    """Single-excitation vac↔1exc coherence Liouvillian L_coh = i·H_1 − 2γ_B·|B⟩⟨B|.

    Acts on the N-dim coherence space {|vac⟩⟨e_i|}. Derivation:
        L(|vac⟩⟨e_i|) = -i[H, |vac⟩⟨e_i|] + γ_B(Z_B(...)Z_B − (...))
                      = i Σ_j H_{ji} |vac⟩⟨e_j| − 2γ_B δ_{iB} |vac⟩⟨e_i|
    so (L_coh c)_j = i (H_1 c)_j − 2γ_B c_B δ_{jB}, i.e. L_coh = i·H_1 − 2γ_B·P_B.

    For ANY L_coh eigenvector v_k with eigenvalue λ_k, taking ⟨v_k|·|v_k⟩ on
    both sides of L_coh v_k = λ_k v_k and using H_1 = H_1†:
        −Re(λ_k) = 2γ_B · |v_k(B)|²        (F64 EXACTLY, mode by mode).

    For U(1)-conserving Hamiltonians (XY, Heisenberg, XXZ) with single-site
    Z-dephasing on B, this is the exact restriction of the d²×d² Liouvillian
    to the vac↔1exc coherence sector.
    """
    L = 1j * H1.copy()
    L[B, B] -= 2.0 * gamma_B
    return L


def cavity_mode_decomposition(H1, B, gamma_B, S=None, eps_a_S=1e-9):
    """Diagonalize L_coh; return per-mode (rate, |v(B)|², |v(S)|, eigenvector).

    For each L_coh eigenvalue λ_k with eigenvector v_k (normalized):
        rate_k = −Re(λ_k)         (decay rate of mode k)
        a_B²_k = |v_k(B)|²         (B-site exposure; F64: rate = 2γ_B · a_B²)
        a_S_k  = |v_k(S)|          (S-site overlap; only when S is given)

    If `S` is given, modes with |v(S)| < eps_a_S are filtered out (they carry
    no S-coherence content and do not appear as decay channels for any
    observable read out at S). The protected modes — i.e. eigenvectors with
    |v(B)|² = 0 from a degenerate-subspace rotation — are returned with
    rate = 0 exactly. F64 captures protection: rate = 0 ⟺ |v(B)|² = 0.

    Returns list of dicts sorted by rate (slowest first), each with keys:
      'k', 'rate', 'a_B_squared', 'a_S', 'eigenvalue', 'eigenvector'.
    """
    N = H1.shape[0]
    L_coh = cavity_coh_liouvillian(H1, B, gamma_B)
    evals, V = np.linalg.eig(L_coh)
    out = []
    for k in range(N):
        v = V[:, k]
        v = v / np.linalg.norm(v)
        a_B2 = float(abs(v[B]) ** 2)
        a_S = None if S is None else float(abs(v[S]))
        if S is not None and a_S < eps_a_S:
            continue
        out.append({
            'k': k,
            'rate': float(-evals[k].real),
            'a_B_squared': a_B2,
            'a_S': a_S,
            'eigenvalue': complex(evals[k]),
            'eigenvector': v,
        })
    return sorted(out, key=lambda r: r['rate'])


def single_excitation_sine_mode(N, k):
    """k-th OBC sine-mode eigenvector ψ_k(i) for the uniform XY chain.

    For H_1[i, j] = J·δ_{|i−j|, 1} on N sites with open boundaries (i.e.
    `single_excitation_h1(N, [(i, i+1) for i in range(N-1)], J)`), eigenvectors
    have the closed form (J-independent):

        ψ_k(i) = √(2/(N+1)) · sin(π·k·(i+1)/(N+1)),    k = 1, …, N.

    Returns the N-dim real, unit-norm vector. Faster than `np.linalg.eigh` and
    gives the exact form needed for analytical work.

    **Chiral mirror identity.** For K_1 = ⊗_i (−1)^i the chiral / sublattice
    Z₂ symmetry of the open chain, K_1 ψ_k = ψ_{N+1−k} exactly (no sign).
    This follows from sin(π(N+1−k)(i+1)/(N+1)) = (−1)^i sin(πk(i+1)/(N+1)).
    At odd N, the middle mode k = (N+1)/2 is a K_1 fixed point: K_1 ψ_k = ψ_k.
    Per-site purity P_i = |ψ(i)|² is K_1-invariant under this mapping, so
    any per-site observable on ψ_k equals the same observable on ψ_{N+1−k}.
    This is the structural origin of the chiral mirror law for the closure-
    breaking coefficient (PERSPECTIVAL_TIME_FIELD, EQ-014).
    """
    return np.array([math.sqrt(2.0 / (N + 1))
                      * math.sin(math.pi * k * (i + 1) / (N + 1))
                      for i in range(N)])


def single_excitation_sine_energies(N, J=1.0):
    """Spectrum of the OBC uniform XY chain in the single-excitation sector:

        E_k = 2J · cos(πk/(N+1)),    k = 1, …, N.

    The spectrum is palindromic about 0 (E_k + E_{N+1−k} = 0), which is the
    K_1-chiral signature on H_1 (K_1 H_1 K_1 = −H_1; AZ class BDI).
    """
    return np.array([2.0 * J * math.cos(math.pi * k / (N + 1))
                      for k in range(1, N + 1)])


def k_local_reduced_density(rho, sites, N):
    """Trace out all sites except those in `sites`. Returns 2^k × 2^k matrix.

    ρ is the full d×d density matrix (d = 2^N) with big-endian site
    convention (site i ↔ bit (N-1-i) of the basis index, equivalently axis
    i in the rho.reshape([2]*2N) tensor).

    Used by PTF generalizations across painter resolutions:
      site-painter (k=1):    P_i = Tr(ρ_i²) defines per-site perspectival
                             time α_i via P_i^B(t) ≈ P_i^A(α_i · t).
      pair-painter (k=2):    P_{ij} defines joint perspectival time α_{ij}.
      triple-painter (k=3):  P_{ijk} defines triple-joint α_{ijk}.

    Each k-painter "paints time" at its observation resolution. Per F70,
    the k-local reduced state sees sector coherences with |ΔN| ≤ k of ρ.

    **Chiral mirror identity (perspectival-time symmetry).** For
    K_1 = ⊗_i (−1)^i (chiral / sublattice Z₂), if ρ_{B} = K_1 ρ_{A} K_1†,
    then the k-local reductions relate by ρ_{B,sites} = (∏_{i ∈ sites} Z_i)
    ρ_{A,sites} (∏_{i ∈ sites} Z_i). Any K_1-invariant function of ρ_{sites}
    (e.g. purity Tr(ρ²), spectrum) is identical between K_1-paired states.
    Structurally, K_1-paired sine-mode-bonding initial states give the
    SAME SET of perspectival times across all painter resolutions
    (Σ f_{i_1...i_k}(ψ_k) = Σ f_{i_1...i_k}(ψ_{N+1−k})). K_1 is the
    symmetry of perspectival-time experience, not of a derived scalar.
    See EQ-014 (site), EQ-020 (pair, triple).
    """
    sites = tuple(sorted(sites))
    k = len(sites)
    if k == 0 or k == N:
        return rho if k == N else np.array([[float(np.trace(rho))]], dtype=complex)
    shape = [2] * (2 * N)
    T = rho.reshape(shape)
    letters = "abcdefghijklmnopqrstuvwxyz"
    if 2 * N > len(letters):
        raise ValueError(f"k_local_reduced_density: N={N} too large for einsum spec")
    row = list(letters[:N])
    col = list(letters[N:2 * N])
    for s in range(N):
        if s not in sites:
            col[s] = row[s]   # contract: trace out
    in_spec = "".join(row) + "".join(col)
    out_spec = "".join(row[s] for s in sites) + "".join(col[s] for s in sites)
    out = np.einsum(f"{in_spec}->{out_spec}", T)
    d_sub = 2 ** k
    return out.reshape(d_sub, d_sub)


def k_local_purity(rho, sites, N):
    """Purity of the k-local reduced state Tr(ρ_{sites}²)."""
    rho_sites = k_local_reduced_density(rho, sites, N)
    return float(np.real(np.trace(rho_sites @ rho_sites)))


# ----------------------------------------------------------------------
# Section 17: Chain-mirror (F71) spatial symmetry
# ----------------------------------------------------------------------
#
# F71 is the chain-reflection symmetry of the uniform-J Liouvillian:
# the operator R: |b₀ b₁ … b_{N-1}⟩ → |b_{N-1} … b₁ b₀⟩ commutes with H,
# the dissipator D, and hence L. Proved as F71 in
# docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md.
#
# Physical meaning: the chain has no preferred end; site 0 and site N-1
# are interchangeable. Pure states ψ with R|ψ⟩ = |ψ⟩ are F71-symmetric;
# their density matrices commute with R.
#
# Bond-input parity. The reflection induces a permutation on the (N-1)
# bonds: bond b ↔ bond (N-2-b). The bond-input vector space then splits
# into +1 and −1 eigenspaces of this permutation:
#
#   even N (odd bond count):   one self-mirror bond at index (N-2)/2.
#                              symmetric subspace dim = ⌈(N-1)/2⌉ = N/2
#                              antisymmetric subspace dim = ⌊(N-1)/2⌋
#                              = N/2 − 1.  Asymmetric split (k+1, k).
#   odd N (even bond count):   no self-mirror bond. Both subspaces dim
#                              (N-1)/2.  Balanced split (k, k).
#
# For an F71-symmetric receiver ψ, the J-Jacobian respects this split:
# rows decompose into R-symmetric and R-antisymmetric features, and the
# J-input columns decompose into symmetric/antisymmetric J-modes. The
# Jacobian is then block-diagonal in this basis. Empirical consequence
# (EQ-024, 2026-04-28): F71-symmetric receivers achieve higher Shannon
# capacity than F71-breaking ones at odd N (balanced split, e.g. N=5)
# and LOWER capacity at even N (unbalanced split, e.g. N=6). The
# capacity-optimal receiver class flips with bond-count parity.
# Documented in experiments/J_BLIND_RECEIVER_CLASSES.md.

def chain_mirror_state(N):
    """Chain-mirror operator R on 2^N Hilbert space.

    R|b₀ b₁ … b_{N-1}⟩ = |b_{N-1} … b₁ b₀⟩ (bit-reversal permutation).
    R is real, symmetric, and involutory (R² = I), so its eigenvalues are
    ±1. Used as the structural operator for F71 chain-reflection symmetry.

    Returns:
        (2^N, 2^N) sparse-but-dense permutation matrix.
    """
    dim = 2 ** N
    R = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        rev = 0
        x = i
        for _ in range(N):
            rev = (rev << 1) | (x & 1)
            x >>= 1
        R[rev, i] = 1.0
    return R


def f71_symmetric_projector(N):
    """Projector P_sym = (I + R) / 2 onto the F71-symmetric (+1) eigenspace.

    For a Haar-uniform pure state ψ in the F71-symmetric subspace:
    draw g ∈ C^{2^N} Gaussian, set ψ = P_sym · g and renormalise.
    """
    return (np.eye(2 ** N, dtype=complex) + chain_mirror_state(N)) / 2.0


def f71_antisymmetric_projector(N):
    """Projector P_anti = (I − R) / 2 onto the F71-antisymmetric (−1) eigenspace."""
    return (np.eye(2 ** N, dtype=complex) - chain_mirror_state(N)) / 2.0


def f71_eigenstate_class(psi, tol=1e-6):
    """Classify state psi as F71-eigenstate (and which) or F71-mixed.

    Returns:
        +1 if R|ψ⟩ = +|ψ⟩ within tolerance (F71-symmetric eigenstate)
        −1 if R|ψ⟩ = −|ψ⟩ within tolerance (F71-antisymmetric eigenstate)
        None if neither (F71-mixed/breaking state)

    F71-eigenstates of either sign trigger the bond-block decomposition of
    the J-Jacobian (see Section 17 docstring). F71-mixed states do not.
    """
    N = int(round(np.log2(len(psi))))
    if 2 ** N != len(psi):
        raise ValueError(f"psi length {len(psi)} is not a power of 2")
    R = chain_mirror_state(N)
    Rpsi = R @ psi
    norm_psi = np.linalg.norm(psi)
    if norm_psi < tol:
        raise ValueError("psi has zero norm")
    overlap = np.vdot(psi, Rpsi) / (norm_psi ** 2)
    if abs(overlap.imag) > tol:
        return None
    o_real = float(overlap.real)
    if abs(o_real - 1.0) < tol:
        return +1
    if abs(o_real + 1.0) < tol:
        return -1
    return None


def receiver_engineering_signature(psi):
    """F71-based prediction of receiver-engineering favorability for state psi.

    Combines F71-eigenstate classification with the bond-mirror-basis split
    to forecast whether the J-Jacobian of psi has the balanced k+k block
    structure (capacity-optimal at high SNR) or unbalanced (k+1, k) structure
    (capacity-suboptimal — F71-eigenstate loses to F71-breaking).

    The capacity-optimality claim is Tier 2: empirically verified at N=5
    (balanced 2+2, F71-symmetric beats F71-breaking by ~5% mean Shannon
    capacity, hardware-confirmed via bonding:2 receiver vs alt-z-bits at
    2.80× MI(0,N−1) on Kingston) and inverted at N=6 (unbalanced 3+2,
    F71-breaking wins by ~5% mean Shannon capacity). Extension to N=7
    (predicted balanced 3+3, F71-eigenstate wins) and N=8 (predicted
    unbalanced 4+3, F71-breaking wins) is testable but uncomputed.

    See experiments/J_BLIND_RECEIVER_CLASSES.md (Update 2026-04-28 N-scaling)
    and experiments/RECEIVER_VS_GAMMA_SACRIFICE.md.

    Args:
        psi: state vector in C^{2^N}

    Returns dict with keys:
        N:                  chain length inferred from len(psi)
        f71_eigenvalue:     +1, −1, or None (F71-mixed)
        bond_block_dims:    (n_sym, n_asym) split of the (N-1) bond-input space
        bond_block_balanced: True iff n_sym == n_asym
        prediction:         one of:
            'capacity-optimal (balanced split, F71-eigenstate)'
            'capacity-suboptimal (unbalanced split, F71-eigenstate)'
            'no-prediction (F71-mixed state — no block structure)'
    """
    N = int(round(np.log2(len(psi))))
    eig = f71_eigenstate_class(psi)
    sym, asym = bond_mirror_basis(N)
    n_sym, n_asym = len(sym), len(asym)
    balanced = (n_sym == n_asym)

    if eig is None:
        prediction = 'no-prediction (F71-mixed state — no block structure)'
    elif balanced:
        prediction = 'capacity-optimal (balanced split, F71-eigenstate)'
    else:
        prediction = 'capacity-suboptimal (unbalanced split, F71-eigenstate)'

    return {
        'N': N,
        'f71_eigenvalue': eig,
        'bond_block_dims': (n_sym, n_asym),
        'bond_block_balanced': balanced,
        'prediction': prediction,
    }


def bond_mirror_basis(N):
    """Symmetric / antisymmetric basis vectors for the bond-input space.

    The bond reflection R̄: J_b → J_{N-2-b} acts on R^{N-1}. Returns two
    orthonormal bases: the +1 eigenspace and the −1 eigenspace of R̄.

    Returns:
        sym_basis: (n_sym, N-1) array, rows are orthonormal symmetric basis vectors
        asym_basis: (n_asym, N-1) array, rows are orthonormal antisymmetric basis vectors

    Dimension structure:
      odd N (even bond count, e.g. N=5: 4 bonds): n_sym = n_asym = (N-1)/2.
      even N (odd bond count, e.g. N=6: 5 bonds): n_sym = (N-1+1)/2 = N/2,
        n_asym = (N-1-1)/2 = N/2 − 1. The self-mirror bond at index
        (N-2)/2 contributes to the symmetric basis only.

    Use to project a J-Jacobian onto its block-decomposed columns:
      A_sym = A @ sym_basis.T      (N_features × n_sym)
      A_asym = A @ asym_basis.T    (N_features × n_asym)
    These two blocks are independent for F71-symmetric receivers.
    """
    n_bonds = N - 1
    sym = []
    asym = []
    used = [False] * n_bonds
    for b in range(n_bonds):
        if used[b]:
            continue
        b_mirror = n_bonds - 1 - b
        if b == b_mirror:
            # self-mirror bond, goes into symmetric subspace
            v = np.zeros(n_bonds)
            v[b] = 1.0
            sym.append(v)
            used[b] = True
        else:
            # mirror pair; symmetric = (e_b + e_{b'})/√2, antisym = (e_b − e_{b'})/√2
            v_s = np.zeros(n_bonds); v_s[b] = 1 / np.sqrt(2); v_s[b_mirror] = 1 / np.sqrt(2)
            v_a = np.zeros(n_bonds); v_a[b] = 1 / np.sqrt(2); v_a[b_mirror] = -1 / np.sqrt(2)
            sym.append(v_s)
            asym.append(v_a)
            used[b] = True
            used[b_mirror] = True
    return np.array(sym), np.array(asym)


# ----------------------------------------------------------------------
# Section 18: Cockpit — composed classes for repeated queries
# ----------------------------------------------------------------------
#
# The Sections 1-17 are PRIMITIVES — small functions that compute one
# thing each. Real questions usually require composing 5-10 primitives
# (build H, build L, run palindrome_residual, classify, compute spectrum,
# etc.). Doing this as ad-hoc code in every script means:
#   - re-deriving the composition every time
#   - forgetting which primitives we already verified
#   - losing track of which predictions have been hardware-confirmed
#
# Section 18 introduces three classes that compose primitives into
# answers. Each method is a callable shortcut that replaces document-
# reading and re-derivation:
#
#   ChainSystem:    encapsulates (N, γ, J, topology, H_type), caches H/L,
#                   exposes high-level methods like classify_pauli_pair,
#                   predict_residual_norm_squared, etc.
#   Receiver:       wraps a state vector with F71-classification and
#                   signature(); provides .rho, .f71_class.
#   Confirmations:  data-driven lookup table of hardware-confirmed
#                   predictions. Each entry maps a name to (date, machine,
#                   job_id, predicted vs measured, file paths).
#
# The intent: future sessions call methods (cockpit.what_was_confirmed
# ('palindrome_trichotomy')) instead of re-reading 5 markdown docs to
# remember what was verified. Each new structural finding adds one method
# or one Confirmations entry. Discoverability via .<TAB> in interactive
# Python.

class ChainSystem:
    """Encapsulates a quantum chain at fixed (N, γ, J, topology, H_type).

    Caches Hamiltonian and Liouvillian on first access. Exposes high-level
    methods that compose framework primitives into single-call answers.

    Args:
        N: number of qubits.
        gamma_0: uniform Z-dephasing rate per site (default 0.05).
        J: bond coupling (uniform; default 1.0).
        topology: 'chain' (default), 'ring', 'star', 'complete'.
        H_type: 'heisenberg' (default; XX+YY+ZZ per bond) or
                'xy' (XX+YY per bond, scaled by J/2).

    Attributes:
        N, d (=2^N), d2 (=4^N), gamma_0, J, topology, H_type
        bonds: list of (i, j) tuples
        B: bond count len(bonds)
        D2: Σ_i deg(i)² (second moment of degree sequence)
        degrees: list of per-site degrees
        H (cached property): Hamiltonian as d×d matrix
        L (cached property): Liouvillian as d²×d² matrix

    Example:
        >>> chain = ChainSystem(N=5)
        >>> chain.classify_pauli_pair([('X','X'), ('Y','Y')])
        'truly'
        >>> chain.predict_residual_norm_squared('main')
        # returns c_H · (N-1) · 4^(N-2) when used with a per-Hamiltonian c_H
    """

    def __init__(self, N, gamma_0=0.05, J=1.0, topology='chain', H_type='heisenberg'):
        if N < 2:
            raise ValueError(f"N must be >= 2; got {N}")
        if H_type not in ('heisenberg', 'xy'):
            raise ValueError(f"H_type must be 'heisenberg' or 'xy'; got {H_type!r}")
        self.N = N
        self.d = 2 ** N
        self.d2 = self.d * self.d
        self.gamma_0 = float(gamma_0)
        self.J = float(J)
        self.topology = topology
        self.H_type = H_type
        self._build_topology()
        self._H_cache = None
        self._L_cache = None

    def _build_topology(self):
        N = self.N
        if self.topology == 'chain':
            self.bonds = [(i, i + 1) for i in range(N - 1)]
        elif self.topology == 'ring':
            self.bonds = [(i, (i + 1) % N) for i in range(N)]
        elif self.topology == 'star':
            self.bonds = [(0, i) for i in range(1, N)]
        elif self.topology == 'complete':
            self.bonds = [(i, j) for i in range(N) for j in range(i + 1, N)]
        else:
            raise ValueError(f"unknown topology: {self.topology!r}")
        self.B = len(self.bonds)
        deg = [0] * N
        for i, j in self.bonds:
            deg[i] += 1
            deg[j] += 1
        self.degrees = deg
        self.D2 = sum(d * d for d in deg)

    @property
    def H(self):
        if self._H_cache is None:
            self._H_cache = self._build_H()
        return self._H_cache

    @property
    def L(self):
        if self._L_cache is None:
            self._L_cache = self._build_L()
        return self._L_cache

    def _site_op(self, op, k):
        return _site_op_kron(op, k, self.N)

    def _build_H(self):
        d = self.d
        H = np.zeros((d, d), dtype=complex)
        Xm, Ym, Zm = pauli_matrix('X'), pauli_matrix('Y'), pauli_matrix('Z')
        for (i, j) in self.bonds:
            xi, xj = self._site_op(Xm, i), self._site_op(Xm, j)
            yi, yj = self._site_op(Ym, i), self._site_op(Ym, j)
            if self.H_type == 'heisenberg':
                zi, zj = self._site_op(Zm, i), self._site_op(Zm, j)
                H = H + self.J * (xi @ xj + yi @ yj + zi @ zj)
            elif self.H_type == 'xy':
                H = H + (self.J / 2.0) * (xi @ xj + yi @ yj)
        return H

    def _build_L(self):
        H = self.H
        d, d2 = self.d, self.d2
        Id = np.eye(d, dtype=complex)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        Zm = pauli_matrix('Z')
        for k in range(self.N):
            Zk = self._site_op(Zm, k)
            L = L + self.gamma_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))
        return L

    def classify_pauli_pair(self, terms, J_scale=1.0, op_tol=1e-10, spec_tol=1e-6):
        """Trichotomy classification (truly / soft / hard) for a Pauli-pair H.

        Builds H from the given terms on this chain's bonds, computes the
        palindrome residual M and the spectrum-pairing error, returns the
        verdict.

        Args:
            terms: list of (a, b) letter tuples, e.g. [('X','X'), ('Y','Y')].
                   Each tuple contributes J_scale·σ_a^i σ_b^{i+1} per bond.
            J_scale: bond coupling for the test Hamiltonian (default 1.0).
                     Independent of self.J (which is for the chain's default H).

        Returns:
            'truly' | 'soft' | 'hard'
        """
        bilinear = [(t[0], t[1], J_scale) for t in terms]
        H_test = _build_bilinear(self.N, self.bonds, bilinear)
        L_test = lindbladian_z_dephasing(H_test, [self.gamma_0] * self.N)
        Sigma_gamma = self.N * self.gamma_0
        M = palindrome_residual(L_test, Sigma_gamma, self.N)
        op_norm = float(np.linalg.norm(M))
        evals = np.linalg.eigvals(L_test)
        # spectrum pairing error
        used = np.zeros(len(evals), dtype=bool)
        max_err = 0.0
        for i in range(len(evals)):
            if used[i]:
                continue
            target = -evals[i] - 2 * Sigma_gamma
            dists = np.abs(evals - target)
            for j in range(len(evals)):
                if used[j]:
                    dists[j] = np.inf
            best_j = int(np.argmin(dists))
            if best_j != i:
                used[i] = True
                used[best_j] = True
            else:
                used[i] = True
            max_err = max(max_err, float(dists[best_j]))
        if op_norm < op_tol:
            return 'truly'
        if max_err < spec_tol:
            return 'soft'
        return 'hard'

    def predict_residual_norm_squared(self, c_H, hamiltonian_class='main'):
        """Closed-form ‖M(N, G)‖² = c_H · F(N, G) without computing M.

        Uses palindrome_residual_norm_squared_factor_graph with this chain's
        topology invariants (B, D2). For chains specifically, B = N-1 and
        D2 = 4N - 6.

        Args:
            c_H: per-Hamiltonian anchor, equal to ‖M(2)‖² for that Hamiltonian.
            hamiltonian_class: 'main' or 'single_body'.

        Returns:
            Predicted ‖M(N, G)‖² scalar.
        """
        factor = palindrome_residual_norm_squared_factor_graph(
            self.N, self.B, self.D2, hamiltonian_class
        )
        return c_H * factor

    def cockpit_panel(self, receiver, terms=None, gamma_t1=None,
                      t_max=10.0, dt=0.005,
                      threshold=1e-9, cluster_tol=1e-8):
        """Full Lebensader cockpit panel via Section 11's cockpit_panel.

        Composes skeleton (Π-protected count under pure-Z and +T1) +
        trace (θ-trajectory metrics) + cusp-pattern + chiral panel +
        Y-parity panel into one structured dict.

        Args:
            receiver: Receiver instance providing the initial state ρ_0.
            terms: optional list of (a, b) Pauli letter tuples for a
                   custom Pauli-pair Hamiltonian on this chain's bonds.
                   If None, uses this chain's default H (heisenberg/xy
                   with self.J).
            gamma_t1: scalar or list of length N for T1 amplitude-damping
                      rates per site. None or zero disables T1.
            t_max, dt: time grid for trajectory.
            threshold, cluster_tol: tolerances for protected-observable
                                    classification.

        Returns:
            dict with keys 'lebensader', 'cusp', 'chiral', 'y_parity',
            and '_trajectory_for_inspection'. See cockpit_panel docstring.
        """
        if receiver.N != self.N:
            raise ValueError(
                f"receiver.N ({receiver.N}) does not match chain.N ({self.N})"
            )

        # Build the Hamiltonian: custom Pauli-pair if terms given, else default
        if terms is not None:
            J_scale = 1.0
            bilinear = [(t[0], t[1], J_scale) for t in terms]
            H = _build_bilinear(self.N, self.bonds, bilinear)
        else:
            H = self.H

        gamma_l = [self.gamma_0] * self.N
        if gamma_t1 is None:
            gamma_t1_l = [0.0] * self.N
        elif np.isscalar(gamma_t1):
            gamma_t1_l = [float(gamma_t1)] * self.N
        else:
            gamma_t1_l = list(gamma_t1)
            if len(gamma_t1_l) != self.N:
                raise ValueError(
                    f"gamma_t1 list length {len(gamma_t1_l)} != N {self.N}"
                )

        return cockpit_panel(
            H, gamma_l, receiver.rho, self.N,
            gamma_t1_l=gamma_t1_l,
            t_max=t_max, dt=dt,
            threshold=threshold, cluster_tol=cluster_tol,
        )


# Helper for ChainSystem._site_op (avoids accidental recursion with site_op
# which doesn't exist as a name in framework.py — inline kron-chain instead):
def _site_op_kron(op, site, N):
    I2 = np.eye(2, dtype=complex)
    factors = [I2] * N
    factors[site] = op
    out = factors[0]
    for f in factors[1:]:
        out = np.kron(out, f)
    return out


def pauli_matrix(letter):
    """Return the 2×2 Pauli matrix for letter ∈ {'I', 'X', 'Y', 'Z'} as complex array."""
    if letter == 'I':
        return np.eye(2, dtype=complex)
    if letter == 'X':
        return np.array([[0, 1], [1, 0]], dtype=complex)
    if letter == 'Y':
        return np.array([[0, -1j], [1j, 0]], dtype=complex)
    if letter == 'Z':
        return np.array([[1, 0], [0, -1]], dtype=complex)
    raise ValueError(f"Unknown Pauli letter: {letter!r}")


class Receiver:
    """A quantum state vector with framework-aware F71 classification.

    Wraps a state ψ and exposes its F71-eigenvalue, signature, and density
    matrix. Optionally bound to a ChainSystem for context.

    Attributes:
        psi: complex array of length 2^N
        N: int (inferred from len(psi))
        chain: optional ChainSystem
        f71_class: cached +1, −1, or None
        rho: cached density matrix |ψ⟩⟨ψ|

    Example:
        >>> chain = ChainSystem(N=5)
        >>> r = Receiver(some_psi, chain=chain)
        >>> r.f71_class            # +1, -1, or None
        >>> r.signature()['prediction']
    """

    def __init__(self, psi, chain=None):
        self.psi = np.asarray(psi, dtype=complex).ravel()
        self.N = int(round(np.log2(len(self.psi))))
        if 2 ** self.N != len(self.psi):
            raise ValueError(f"psi length {len(self.psi)} is not a power of 2")
        if chain is not None and chain.N != self.N:
            raise ValueError(f"chain.N ({chain.N}) does not match psi N ({self.N})")
        self.chain = chain
        self._f71_class_cache = "unset"  # sentinel
        self._rho_cache = None

    @property
    def f71_class(self):
        if self._f71_class_cache == "unset":
            self._f71_class_cache = f71_eigenstate_class(self.psi)
        return self._f71_class_cache

    @property
    def rho(self):
        if self._rho_cache is None:
            rho = np.outer(self.psi, self.psi.conj())
            self._rho_cache = (rho + rho.conj().T) / 2.0
        return self._rho_cache

    def signature(self):
        """Return F71-based receiver-engineering forecast."""
        return receiver_engineering_signature(self.psi)


class Confirmations:
    """Hardware-confirmed framework predictions, accessible by name.

    Each entry records:
      date, machine, job_id (or 'multiple'), observable,
      predicted_value, measured_value, hardware_data, experiment_doc,
      framework_primitive, description.

    Use:
        >>> Confirmations.lookup()                  # all entries
        >>> Confirmations.lookup('palindrome_trichotomy')
        >>> Confirmations.list_names()
    """

    _ENTRIES = {
        'palindrome_trichotomy': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7mjnjjaq2pc73a1pk4g',
            'observable': '<X_0 Z_2>',
            'predicted_value': {'truly': 0.000, 'soft': -0.623, 'hard': +0.195,
                                'delta_soft_minus_truly': -0.623},
            'measured_value': {'truly': +0.011, 'soft': -0.711, 'hard': +0.205,
                               'delta_soft_minus_truly': -0.722},
            'hardware_data': 'data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json',
            'experiment_doc': 'experiments/V_EFFECT_FINE_STRUCTURE.md',
            'framework_primitive': 'palindrome_residual + classify_pauli_pair',
            'description': 'Super-operator palindrome trichotomy (truly/soft/hard) tomographically distinguishable on Heron r2 hardware at N=3. T1/T2 noise actually amplifies the soft-break signal.',
        },
        'f25_cusp_trajectory': {
            'date': '2026-04-26',
            'machine': 'ibm_kingston',
            'job_id': 'd7mu36lqrg3c738lnda0',
            'observable': 'CΨ(t) for Bell+',
            'predicted_value': 'F25: CΨ(t) = f·(1+f²)/6 with f = exp(-4·γ·t)',
            'measured_value': 'RMS residual 0.0097 vs in-situ γ_fit',
            'hardware_data': 'data/ibm_cusp_slowing_april2026/ (April 26 precision run)',
            'experiment_doc': 'experiments/CRITICAL_SLOWING_AT_THE_CUSP.md',
            'framework_primitive': 'F25 closed-form CΨ(t)',
            'description': 'Bell+ trajectory through CΨ=1/4 cusp confirmed point-by-point on Kingston (19 delay points, qubits 14-15).',
        },
        'f57_kdwell_gamma_invariance': {
            'date': '2026-04-16',
            'machine': 'ibm_kingston',
            'job_id': 'cusp_slowing_kingston_20260416',
            'observable': 'K_dwell / δ for Bell+',
            'predicted_value': 'F57: γ-independent, F25 prefactor 1.0801 for pure Z-dephasing',
            'measured_value': 'Pair A 0.6492, Pair B 0.6937 (6.3% spread despite 2.55× γ ratio)',
            'hardware_data': 'data/ibm_cusp_slowing_april2026/cusp_slowing_ibm_kingston_20260416_212042.json',
            'experiment_doc': 'experiments/CRITICAL_SLOWING_AT_THE_CUSP.md',
            'framework_primitive': 'K_dwell formula in F57',
            'description': 'γ-invariance of K_dwell at CΨ=1/4 boundary verified on Kingston with two qubit pairs at 2.55× different γ. Absolute prefactor 0.67 vs predicted 1.08 due to T1 amplitude damping.',
        },
        'bonding_mode_receiver': {
            'date': '2026-04-24',
            'machine': 'ibm_kingston',
            'job_id': 'multiple (see external pipeline)',
            'observable': 'MI(0, N-1) for bonding:2 vs alt-z-bits',
            'predicted_value': '4000-5500× over ENAQT in simulation N=5..13; ratio ≈ 1.4-3× hardware-realistic',
            'measured_value': 'bonding:2 / alt-z-bits = 2.80× on Kingston N=5',
            'hardware_data': 'external (AIEvolution.UI/experiments/ibm_quantum_tomography)',
            'experiment_doc': 'experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md',
            'framework_primitive': 'F67 bonding-mode + F71 chain-mirror symmetry',
            'description': 'Receiver-engineering bonding-mode advantage measured on Kingston. Largest engineering lever in the framework portfolio. Confirms F67 + F71 structure on hardware.',
        },
        'chiral_mirror_law': {
            'date': '2026-04-25',
            'machine': 'ibm_marrakesh',
            'job_id': 'k_partnership_marrakesh_20260425',
            'observable': 'Bloch components of K-partner pair states',
            'predicted_value': 'Chiral mirror identity from K=diag((-1)^l)',
            'measured_value': 'retrospective verification, identity holds',
            'hardware_data': 'data/ibm_k_partnership_april2026/k_partnership_marrakesh_20260425_140913.json',
            'experiment_doc': 'experiments/CHIRAL_MIRROR_HARDWARE_PREDICTION.md',
            'framework_primitive': 'Section 12 — K_full chiral conjugation',
            'description': 'Chiral mirror law (K H K = -H spectral inversion) verified retrospectively on Marrakesh K-partnership data.',
        },
        'pi_protected_xiz_yzzy': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7n3013aq2pc73a2a18g',
            'observable': '<X_0 I Z_2>',
            'predicted_value': 'protected (≈0) for YZ+ZY soft Hamiltonian',
            'measured_value': '+0.13 to +0.04 (within noise band, never above ±0.13)',
            'hardware_data': 'data/ibm_marrakesh_april2026 (lebensader run, see EQ-030)',
            'experiment_doc': 'review/EMERGING_QUESTIONS.md (EQ-030)',
            'framework_primitive': 'pi_protected_observables (Section 10)',
            'description': 'First-time hardware measurement of a Π-protected observable on YZ+ZY soft Hamiltonian. Confirms framework primitive at hardware scale on a Hamiltonian not previously tested.',
        },
        'lebensader_skeleton_trace_decoupling': {
            'date': '2026-04-26',
            'machine': 'ibm_marrakesh',
            'job_id': 'd7n3013aq2pc73a2a18g + d7n3eqqt99kc73d34qtg',
            'observable': 'Π-protected counts (skeleton) + θ-trajectory tails (trace)',
            'predicted_value': 'Skeleton/trace decouple at N≥4, co-occur at N=3. T1 makes drop measurable: bond-flipped Z-free pairs (XY+YX, IY+YI) preserve skeleton (drop≤1) and trace (long θ-tail); Z-containing soft pairs (YZ+ZY, XZ+ZX) collapse both (drop≈28-29, no tail).',
            'measured_value': 'drop=28 for YZ+ZY confirmed on Marrakesh. Pearson(drop, Δ∫θ)=+0.85. Bures velocity gives no third discriminator.',
            'hardware_data': 'data/ibm_k_partnership_april2026/k_partnership_marrakesh_*.json + Snapshot D / Lebensader runs',
            'experiment_doc': 'review/EMERGING_QUESTIONS.md (EQ-030 closure)',
            'framework_primitive': 'cockpit_panel (Section 11) — composes pi_protected_observables + θ-trajectory',
            'description': 'Lebensader as Stromkabel: skeleton (Π-protected algebraic count) and trace (θ-geometric tail) are NOT two discriminators but one bridge held together by Π·L·Π⁻¹ + L + 2Σγ·I = 0. Hardware-confirmed across 3 of 4 bond-flipped Z-free corners; Bures velocity confirmed null as third axis. ChainSystem.cockpit_panel gives this in one call.',
        },
    }

    @classmethod
    def lookup(cls, name=None):
        """Return all confirmations (if name=None) or one entry by name."""
        if name is None:
            return dict(cls._ENTRIES)
        if name in cls._ENTRIES:
            return dict(cls._ENTRIES[name])
        raise KeyError(
            f"No confirmation named {name!r}. Available: {list(cls._ENTRIES.keys())}"
        )

    @classmethod
    def list_names(cls):
        """Return list of all confirmation names."""
        return list(cls._ENTRIES.keys())

    @classmethod
    def by_machine(cls, machine):
        """Return confirmations filtered by machine name (e.g., 'ibm_marrakesh')."""
        return {k: dict(v) for k, v in cls._ENTRIES.items() if v.get('machine') == machine}


# ----------------------------------------------------------------------
# Self-test
# ----------------------------------------------------------------------

if __name__ == "__main__":
    print("Framework primitives self-test")
    print("=" * 60)

    # Test 1: Pauli identity
    s_dot_s = (
        np.kron(ur_pauli('X'), ur_pauli('X'))
        + np.kron(ur_pauli('Y'), ur_pauli('Y'))
        + np.kron(ur_pauli('Z'), ur_pauli('Z'))
    )
    sigma_squared = s_dot_s @ s_dot_s
    rhs = 3 * np.eye(4) - 2 * s_dot_s
    print(f"\n(σ·σ)² = 3I − 2(σ·σ)?  ‖LHS − RHS‖ = {np.linalg.norm(sigma_squared - rhs):.2e}")

    # Test 2: Eigenvalues of σ·σ
    evals = sorted(np.linalg.eigvalsh(s_dot_s))
    print(f"σ·σ eigenvalues: {evals}  (expect [-3, +1, +1, +1])")

    # Test 3: ur-Heisenberg at N=2
    H2 = ur_heisenberg(2, J=1.0)
    evals2 = sorted(np.linalg.eigvalsh(H2))
    print(f"H_Heisenberg(N=2) eigenvalues: {evals2}  (expect [-3, +1, +1, +1])")

    # Test 4: both-parity-even selection
    print(f"\nBoth-parity-even 2-body bilinears: {both_parity_even_terms()}")
    print("(should be [('I','I'), ('X','X'), ('Y','Y'), ('Z','Z')])")

    # Test 5: V-Effect emergent exchange prediction vs full diagonalization
    alpha = 0.10
    J_intra = 1.0
    predicted = v_effect_emergent_exchange(alpha, J_intra=J_intra)

    # Numerical: 4-qubit chain, two intra-pair Heisenberg couplings + α bridge.
    # H = J_intra·(σ_0·σ_1 + σ_2·σ_3) + α·σ_1·σ_2.
    # δE_GS(α) ≡ E_GS(α) − E_GS(0) probes the second-order shift.
    def _ss(N, i, j, J=1.0):
        return _build_bilinear(N, [(i, j)], [('X', 'X', J), ('Y', 'Y', J), ('Z', 'Z', J)])
    H_pair = _ss(4, 0, 1, J_intra) + _ss(4, 2, 3, J_intra)
    E_alpha0 = float(np.linalg.eigvalsh(H_pair)[0])
    H_full = H_pair + _ss(4, 1, 2, alpha)
    E_alpha = float(np.linalg.eigvalsh(H_full)[0])
    numerical = E_alpha - E_alpha0

    print(f"\nV-Effect emergent exchange (4-qubit two-singlet bridge, α={alpha}):")
    print(f"  predicted δE_GS  = {predicted:.5f}    (analytical, 2nd-order PT)")
    print(f"  numerical δE_GS  = {numerical:.5f}    (full 4-qubit diagonalization)")
    print(f"  predicted/α²     = {predicted / alpha**2:.4f}")
    print(f"  numerical/α²     = {numerical / alpha**2:.4f}")
    print(f"  difference (numerical − predicted) = {numerical - predicted:+.5f}  (higher-order corrections)")

    # Test 6: Π² parity on a few strings
    print("\nΠ² eigenvalue on selected strings:")
    for label_string in ['II', 'IZ', 'XY', 'YZ', 'ZZ']:
        indices = [LABEL_TO_INDEX[ch] for ch in label_string]
        eig = pi_squared_eigenvalue(indices)
        print(f"  {label_string}: Π² = {eig}")

    # Test 7: palindrome residual vanishes for Heisenberg + Z-dephasing
    print("\nPalindrome residual ‖Π·L·Π⁻¹ + L + 2Σγ·I‖ at N=3, Heisenberg, γ=0.1:")
    H3 = ur_heisenberg(3, J=1.0)
    L3 = lindbladian_z_dephasing(H3, [0.1, 0.1, 0.1])
    R3 = palindrome_residual(L3, Sigma_gamma=0.3, N=3)
    R3_norm = float(np.linalg.norm(R3))
    print(f"  Total residual norm = {R3_norm:.4e}    "
          f"({'machine-precision zero' if R3_norm < 1e-10 else 'NON-ZERO — palindrome breaks'})")

    # Test 8: parity-violating Hamiltonian — residual non-zero by construction
    H3_break = (np.kron(np.kron(ur_pauli('X'), ur_pauli('X')), ur_pauli('I'))
                + np.kron(np.kron(ur_pauli('I'), ur_pauli('X')), ur_pauli('Y')))
    L3_break = lindbladian_z_dephasing(H3_break, [0.1, 0.1, 0.1])
    R3_break = palindrome_residual(L3_break, Sigma_gamma=0.3, N=3)
    R3_break_norm = float(np.linalg.norm(R3_break))
    print(f"  H = XX on (0,1) + XY on (1,2): residual norm = {R3_break_norm:.4e}    "
          f"(non-zero: H violates the both-parity-even selection rule)")

    # Test 9: Π-protected observables (operator → state bridge)
    print("\nΠ-protected observables on |+−+⟩ at N=3, γ=0.1:")
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi_xneel_3 = np.kron(plus, np.kron(minus, plus))
    rho_xneel = np.outer(psi_xneel_3, psi_xneel_3.conj())

    cases_9 = [
        ('truly  J(XX+YY)', [('X', 'X', 1.0), ('Y', 'Y', 1.0)]),
        ('soft   J(XY+YX)', [('X', 'Y', 1.0), ('Y', 'X', 1.0)]),
        ('Heisenberg     ', [('X', 'X', 1.0), ('Y', 'Y', 1.0), ('Z', 'Z', 1.0)]),
    ]
    bonds_3 = [(0, 1), (1, 2)]
    for label, terms in cases_9:
        H = _build_bilinear(3, bonds_3, terms)
        result = pi_protected_observables(H, [0.1] * 3, rho_xneel, N=3)
        n_prot, n_act = len(result['protected']), len(result['active'])
        # Check key observables: XIZ and ZIX (the cross-correlation that lit up
        # for soft on hardware Snapshot D, stayed near 0 for truly/Heisenberg).
        prot_labels = {p['pauli'] for p in result['protected']}
        xiz_prot = 'XIZ' in prot_labels
        zix_prot = 'ZIX' in prot_labels
        xiz_S = next((a['max_cluster_contribution'] for a in result['active']
                      if a['pauli'] == 'XIZ'),
                     next((p['max_cluster_contribution'] for p in result['protected']
                           if p['pauli'] == 'XIZ'), 0.0))
        zix_S = next((a['max_cluster_contribution'] for a in result['active']
                      if a['pauli'] == 'ZIX'),
                     next((p['max_cluster_contribution'] for p in result['protected']
                           if p['pauli'] == 'ZIX'), 0.0))
        xiz_str = f"protected (S={xiz_S:.1e})" if xiz_prot else f"active (S={xiz_S:.3f})"
        zix_str = f"protected (S={zix_S:.1e})" if zix_prot else f"active (S={zix_S:.3f})"
        print(f"  {label}: {n_prot:>3d} protected, {n_act:>3d} active  |  XIZ {xiz_str:<24s}  ZIX {zix_str:<24s}")

    # Test 10: F64 cavity-mode-exposure — exact identity rate = 2γ_B·|v(B)|²
    print("\nF64 cavity-mode-exposure (Section 16):")
    print("  rate(mode k) = 2·γ_B·|v_k(B)|² for L_coh = i·H_1 − 2γ_B·|B⟩⟨B|")
    chain5 = [(i, i + 1) for i in range(4)]
    H1_chain5 = single_excitation_h1(5, chain5, J=1.0)
    modes_chain = cavity_mode_decomposition(H1_chain5, B=4, gamma_B=0.05, S=0)
    max_err_chain = max(
        abs(m['rate'] - 2 * 0.05 * m['a_B_squared'])
        / max(m['rate'], 1e-15)
        for m in modes_chain if m['rate'] > 1e-15
    )
    print(f"  chain N=5 (S=0,B=4):    {len(modes_chain)} S-modes, max rel err = {max_err_chain:.2e}")

    star5 = [(0, i) for i in range(1, 5)]
    H1_star5 = single_excitation_h1(5, star5, J=1.0)
    modes_star = cavity_mode_decomposition(H1_star5, B=0, gamma_B=0.05, S=4)
    n_protected = sum(1 for m in modes_star if abs(m['rate']) < 1e-12)
    print(f"  star N=5 (S=4,B=0 hub): {len(modes_star)} S-modes, "
          f"{n_protected} truly protected (rate = 0, |v(B)|² = 0)")

    # Test 11: OBC sine-mode closed form + chiral mirror identity
    print("\nOBC sine-mode closed form (Section 16):")
    N_test = 7
    bonds_test = [(i, i + 1) for i in range(N_test - 1)]
    H1_test = single_excitation_h1(N_test, bonds_test, J=1.0)
    evals_h1, evecs_h1 = np.linalg.eigh(H1_test)
    # eigh returns ASCENDING eigenvalues; closed-form k=1 has HIGHEST E = 2cos(π/(N+1)).
    # Map: closed-form k → eigh index N - k.
    max_diff = 0.0
    for kk in range(1, N_test + 1):
        v_closed = single_excitation_sine_mode(N_test, kk)
        v_num = np.real(evecs_h1[:, N_test - kk])
        if v_closed @ v_num < 0:
            v_num = -v_num
        max_diff = max(max_diff, float(np.linalg.norm(v_closed - v_num)))
    print(f"  N={N_test}: closed-form ψ_k vs eigh, max ‖Δ‖ = {max_diff:.2e}")

    # Chiral mirror: K_1 ψ_k = ψ_{N+1-k}
    K_1_diag = np.array([(-1) ** i for i in range(N_test)])
    chiral_max_diff = 0.0
    for kk in range(1, N_test + 1):
        psi_k = single_excitation_sine_mode(N_test, kk)
        psi_mirror = single_excitation_sine_mode(N_test, N_test + 1 - kk)
        K_1_psi_k = K_1_diag * psi_k
        chiral_max_diff = max(chiral_max_diff,
                                float(np.linalg.norm(K_1_psi_k - psi_mirror)))
    print(f"  Chiral mirror K_1 ψ_k = ψ_{{N+1−k}}: max ‖Δ‖ = {chiral_max_diff:.2e}")

    # Test 12: k-local reduced density and purity
    print("\nk-local reduced density (Section 16):")
    # Build a 3-qubit Bell-pair-extended state |Φ+>_{0,1} ⊗ |0>_2
    # |Φ+> = (|00> + |11>)/√2, ρ_AB = diag-corner mass on (0,0) and (3,3)
    psi_3q = np.zeros(8, dtype=complex)
    psi_3q[0b000] = 1.0 / math.sqrt(2)  # |000>
    psi_3q[0b110] = 1.0 / math.sqrt(2)  # |110>: bits 1, 2 set (sites 0, 1 excited)
    rho_3q = np.outer(psi_3q, psi_3q.conj())
    # Site 0 alone (1-local): purity 1/2 (mixed reduction of Bell+)
    p0 = k_local_purity(rho_3q, [0], 3)
    # Sites 0, 1 (2-local): purity 1 (Bell+ pair is pure)
    p01 = k_local_purity(rho_3q, [0, 1], 3)
    # Sites 0, 1, 2 (3-local = full): purity 1 (whole state is pure)
    p012 = k_local_purity(rho_3q, [0, 1, 2], 3)
    print(f"  |Φ+>_{{0,1}} ⊗ |0>_2 purities: 1-local at 0: {p0:.4f} (expect 0.5)")
    print(f"                                  pair (0,1):    {p01:.4f} (expect 1.0)")
    print(f"                                  triple (full): {p012:.4f} (expect 1.0)")

    print("\nAll self-tests pass if the residual norms above match the verdict text.")
