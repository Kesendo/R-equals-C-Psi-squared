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
# Section 4: Parity selection — which 2-body bilinears respect framework
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

def lindbladian_z_dephasing(H, gamma_l):
    """Lindbladian L(ρ) = -i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l - ρ).

    Returned as a d²×d² matrix in column-stack vec convention.
    """
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


def palindrome_residual(L, Sigma_gamma, N):
    """Compute ‖Π·L·Π⁻¹ + L + 2Σγ·I‖ in Pauli-string basis.

    Returns the residual matrix (4^N × 4^N) for sector-block analysis.
    """
    Pi = build_pi_full(N)
    M = _vec_to_pauli_basis_transform(N)
    L_pauli = np.linalg.solve(M.conj().T @ M, M.conj().T) @ L @ M
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

    # Test 5: V-Effect emergent exchange prediction
    alpha = 0.10
    print(f"\nV-Effect at α={alpha}, J_intra=1: predicted δE_GS = "
          f"{v_effect_emergent_exchange(alpha):.5f}")
    print(f"  (compare to numerical δE/α² ≈ -0.384, i.e., δE ≈ -0.00384)")

    # Test 6: Π² parity on a few strings
    print("\nΠ² eigenvalue on selected strings:")
    for label_string in ['II', 'IZ', 'XY', 'YZ', 'ZZ']:
        indices = [LABEL_TO_INDEX[ch] for ch in label_string]
        eig = pi_squared_eigenvalue(indices)
        print(f"  {label_string}: Π² = {eig}")

    print("\nAll self-tests pass if numerical results match expectations.")
