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

    # Lebensader rating: simple qualitative reading from skeleton + trace
    if skeleton['drop'] <= 1 and tail_duration > 0.05:
        leb_rating = 'intact (skeleton holds, trace lives)'
    elif skeleton['drop'] <= 1 and tail_duration <= 0.05:
        leb_rating = 'partial (skeleton holds, trace short)'
    elif skeleton['drop'] > 1 and tail_duration > 0.05:
        leb_rating = 'partial (skeleton bleeds, trace persists)'
    else:
        leb_rating = 'collapsed (both skeleton and trace gone)'

    # ---- Chiral panel: K_full sublattice symmetry ----
    chiral = chiral_panel(H, rho_0, N)

    # ---- Y-parity panel: bit_a · bit_b grading ----
    y_parity = y_parity_panel(H, gamma_l, rho_0, N, gamma_t1_l=gamma_t1_l)

    return {
        'lebensader': {'skeleton': skeleton, 'trace': trace,
                       'rating': leb_rating},
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

    print("\nAll self-tests pass if the residual norms above match the verdict text.")
