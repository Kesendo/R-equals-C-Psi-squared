"""
Q3: Entanglement Structure of the Pi Operator

Is the palindromic conjugation operator Pi entangled in operator space?

Pi acts on Liouville space C^{d^2}. For N qubits, Liouville space has
the tensor product structure C^{d_A^2} x C^{d_B^2}. If Pi cannot be
written as P_A x P_B, it is "operator-entangled": the palindromic
symmetry is inherently non-local.

Method: build Pi from eigenvalue pairing, compute operator Schmidt
decomposition (SVD of reshaped operator matrix), report rank and
coefficients.

Cases:
  1. N=2 Heisenberg (XX+YY+ZZ) + Z-deph  (should be product, rank 1)
  2. N=2 XZ+ZY + Z-deph                   (known non-local Pi)
  3. N=2 YZ+ZX + Z-deph                   (known non-local Pi)
  4. N=3 Heisenberg chain + Z-deph         (should be product across A|BC)
  5. N=3 XY + Z-deph chain                 (alternating Pi family)
"""

import numpy as np
from itertools import product as iproduct

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS_MAT = [I2, sx, sy, sz]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']


def kron_list(mats):
    """Kronecker product of a list of matrices."""
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result


def build_hamiltonian(N, terms, J=1.0, topology='chain'):
    """Build N-qubit Hamiltonian from Pauli coupling terms.

    terms: list of (pauli_A_matrix, pauli_B_matrix) for each bond
    topology: 'chain' (bonds 0-1, 1-2, ...) or 'single' (bond 0-1 only)
    """
    d = 2**N
    H = np.zeros((d, d), dtype=complex)

    if topology == 'single' or N == 2:
        bonds = [(0, 1)]
    else:
        bonds = [(i, i + 1) for i in range(N - 1)]

    for site_a, site_b in bonds:
        for pA, pB in terms:
            ops = [I2] * N
            ops[site_a] = pA
            ops[site_b] = pB
            H += J * kron_list(ops)
    return H


def build_z_dephasing(N, gamma=0.05):
    """Z-dephasing Lindblad operators on each qubit."""
    c_ops = []
    for k in range(N):
        ops = [I2] * N
        ops[k] = np.sqrt(gamma) * sz
        c_ops.append(kron_list(ops))
    return c_ops, N * gamma


def build_liouvillian(H, c_ops):
    """Build Liouvillian superoperator L. vec(L(rho)) = L @ vec(rho)."""
    d = H.shape[0]
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for c in c_ops:
        cd = c.conj().T
        cdc = cd @ c
        L += (np.kron(c, c.conj())
              - 0.5 * np.kron(cdc, eye)
              - 0.5 * np.kron(eye, cdc.T))
    return L


def find_pi_operator(L, sum_gamma, tol=1e-6):
    """Construct Pi by pairing eigenvalues lambda <-> -lambda - 2*sum_gamma."""
    d2 = L.shape[0]
    evals, R = np.linalg.eig(L)

    # Pair eigenvalues
    target = -evals - 2 * sum_gamma
    pairs = []
    used = set()

    for i in range(d2):
        if i in used:
            continue
        dists = np.abs(evals - target[i])
        # Find closest unpaired eigenvalue to target[i]
        best_j, best_d = None, np.inf
        for j in range(d2):
            if j not in used and dists[j] < best_d:
                best_j = j
                best_d = dists[j]
        if best_d < tol:
            pairs.append((i, best_j))
            used.add(i)
            used.add(best_j)

    if len(used) != d2:
        print(f"  WARNING: {d2 - len(used)} eigenvalues unpaired!")
        return None, pairs, evals

    # Construct Pi from eigenvector pairing
    L_left = np.linalg.inv(R)
    Pi = np.zeros((d2, d2), dtype=complex)
    for i, j in pairs:
        Pi += np.outer(R[:, j], L_left[i, :])
        if i != j:
            Pi += np.outer(R[:, i], L_left[j, :])

    return Pi, pairs, evals


def verify_pi(Pi, L, sum_gamma):
    """Check Pi L Pi^{-1} = -L - 2*sum_gamma*I. Return relative error."""
    d2 = L.shape[0]
    Pi_inv = np.linalg.inv(Pi)
    LHS = Pi @ L @ Pi_inv
    RHS = -L - 2 * sum_gamma * np.eye(d2, dtype=complex)
    err = np.linalg.norm(LHS - RHS) / np.linalg.norm(RHS)
    return err


def build_bipartition_perm_2q():
    """Permutation from standard vec ordering to A|B grouped ordering for N=2.

    Standard: k = 4*i + j, i = 2*i_A + i_B, j = 2*j_A + j_B
    Grouped:  m = 4*alpha + beta, alpha = 2*i_A + j_A, beta = 2*i_B + j_B
    """
    perm = np.zeros(16, dtype=int)
    for i_A in range(2):
        for j_A in range(2):
            for i_B in range(2):
                for j_B in range(2):
                    k = 8 * i_A + 4 * i_B + 2 * j_A + j_B
                    m = 8 * i_A + 4 * j_A + 2 * i_B + j_B
                    perm[m] = k
    return perm


def build_bipartition_perm_3q_A_BC():
    """Permutation from standard vec to A|BC grouped ordering for N=3.

    d_A=2, d_BC=4, d=8, d^2=64.
    Standard: k = 8*i + j, i = 4*i_A + i_BC, j = 4*j_A + j_BC
    Grouped:  m = 16*alpha + beta, alpha = 2*i_A + j_A, beta = 4*i_BC + j_BC
    """
    d_A, d_BC = 2, 4
    d = d_A * d_BC  # 8
    d2 = d * d  # 64
    perm = np.zeros(d2, dtype=int)
    for i_A in range(d_A):
        for j_A in range(d_A):
            for i_BC in range(d_BC):
                for j_BC in range(d_BC):
                    i = d_BC * i_A + i_BC
                    j = d_BC * j_A + j_BC
                    k = d * i + j
                    alpha = d_A * i_A + j_A
                    beta = d_BC * i_BC + j_BC
                    m = d_BC**2 * alpha + beta
                    perm[m] = k
    return perm


def operator_schmidt(Pi, perm, d_A2, d_B2):
    """Operator Schmidt decomposition of Pi across bipartition.

    perm: permutation from standard to grouped ordering
    d_A2: dimension of A Liouville space (d_A^2)
    d_B2: dimension of B Liouville space (d_B^2)

    Returns sorted Schmidt coefficients (singular values of reshaped operator).
    """
    Pi_r = Pi[np.ix_(perm, perm)]

    # Reshape to 4-index tensor: T[alpha, beta, alpha', beta']
    T = Pi_r.reshape(d_A2, d_B2, d_A2, d_B2)

    # Regroup: M[(alpha, alpha'), (beta, beta')]
    M = T.transpose(0, 2, 1, 3).reshape(d_A2**2, d_B2**2)

    # SVD
    U, S, Vh = np.linalg.svd(M, full_matrices=False)

    return S


def analyze_case(name, N, terms, topology='chain', J=1.0, gamma=0.05):
    """Full analysis pipeline for one case."""
    print(f"\n{'=' * 60}")
    print(f"Case: {name}")
    print(f"N={N}, J={J}, gamma={gamma}, topology={topology}")
    print('=' * 60)

    # Build system
    H = build_hamiltonian(N, terms, J=J, topology=topology)
    c_ops, sum_gamma = build_z_dephasing(N, gamma)
    L = build_liouvillian(H, c_ops)

    # Find Pi
    Pi, pairs, evals = find_pi_operator(L, sum_gamma)
    if Pi is None:
        print("  FAILED: Could not construct Pi (incomplete pairing)")
        return None

    # Verify
    err = verify_pi(Pi, L, sum_gamma)
    print(f"  Conjugation equation error: {err:.2e}")
    if err > 1e-6:
        print("  WARNING: Pi does not satisfy conjugation equation!")
        return None

    print(f"  Eigenvalue pairs: {len(pairs)}")
    print(f"  Pi shape: {Pi.shape}")
    print(f"  Pi is {'real' if np.allclose(Pi.imag, 0) else 'complex'}")

    # Schmidt decomposition
    if N == 2:
        perm = build_bipartition_perm_2q()
        S = operator_schmidt(Pi, perm, d_A2=4, d_B2=4)
        cut_label = "A|B"
    elif N == 3:
        perm = build_bipartition_perm_3q_A_BC()
        S = operator_schmidt(Pi, perm, d_A2=4, d_B2=16)
        cut_label = "A|BC"
    else:
        print("  Schmidt analysis not implemented for N > 3")
        return None

    # Report
    S_normalized = S / S[0] if S[0] > 0 else S
    rank = np.sum(S > 1e-10 * S[0])

    print(f"\n  Operator Schmidt decomposition ({cut_label} cut):")
    print(f"  Schmidt rank: {rank}")
    print(f"  Schmidt coefficients (normalized to largest):")
    for i, (s, sn) in enumerate(zip(S[:min(rank + 2, len(S))], S_normalized)):
        marker = " *" if s > 1e-10 * S[0] else ""
        print(f"    sigma_{i} = {s:.6f}  (relative: {sn:.6f}){marker}")

    if rank == 1:
        print(f"\n  >>> Pi is a PRODUCT operator (rank 1)")
        print(f"  >>> The palindromic symmetry is LOCAL: Pi = Pi_A x Pi_B")
    else:
        print(f"\n  >>> Pi is ENTANGLED in operator space (rank {rank})")
        print(f"  >>> The palindromic symmetry is NON-LOCAL")
        # Entanglement entropy
        S_prob = (S[:rank]**2) / np.sum(S[:rank]**2)
        ent_entropy = -np.sum(S_prob * np.log2(S_prob + 1e-30))
        print(f"  >>> Operator entanglement entropy: {ent_entropy:.4f} bits")
        print(f"  >>> Max possible: {np.log2(min(4**2, (2**(N-1))**4)):.1f} bits")

    return {'rank': rank, 'S': S, 'Pi': Pi}


def build_single_site_map(pi_perm, pi_sign):
    """Build the 4x4 vec-basis superoperator for one per-site Pauli map."""
    PAULI_MATS = [I2, sx, sy, sz]
    M_vec = np.zeros((4, 4), dtype=complex)
    for col_idx in range(4):
        i, j = col_idx // 2, col_idx % 2
        result = np.zeros((2, 2), dtype=complex)
        for a in range(4):
            coeff = PAULI_MATS[a][j, i]
            result += pi_sign[a] * coeff * 0.5 * PAULI_MATS[pi_perm[a]]
        for row_idx in range(4):
            k, l = row_idx // 2, row_idx % 2
            M_vec[row_idx, col_idx] = result[k, l]
    return M_vec


def grouped_to_standard_perm(N):
    """Build permutation matrix P: grouped vec -> standard vec for N qubits."""
    d = 2**N
    d2 = d * d
    perm_to_std = np.zeros(d2, dtype=int)
    for m in range(d2):
        alphas = []
        temp = m
        for _ in range(N):
            alphas.append(temp % 4)
            temp //= 4
        alphas.reverse()
        i_bits = [a // 2 for a in alphas]
        j_bits = [a % 2 for a in alphas]
        i_val = sum(b * (2**(N - 1 - k)) for k, b in enumerate(i_bits))
        j_val = sum(b * (2**(N - 1 - k)) for k, b in enumerate(j_bits))
        perm_to_std[m] = d * i_val + j_val
    P = np.zeros((d2, d2), dtype=complex)
    for m_idx in range(d2):
        P[perm_to_std[m_idx], m_idx] = 1.0
    return P


def build_analytical_pi_multi(N, site_maps):
    """Build Pi from per-site maps (may differ per site) in standard vec basis.

    site_maps: list of N (perm, sign) tuples, one per site.
    """
    M_vecs = [build_single_site_map(p, s) for p, s in site_maps]
    Pi_grouped = M_vecs[0].copy()
    for m in M_vecs[1:]:
        Pi_grouped = np.kron(Pi_grouped, m)
    P = grouped_to_standard_perm(N)
    return P @ Pi_grouped @ P.conj().T


def commutant_dimension(L, tol=1e-8):
    """Dimension of {S : [S, L] = 0} = kernel of the adjoint map ad_L.

    ad_L(S) = LS - SL. Vectorize: ad_L = L (x) I - I (x) L^T.
    dim(ker(ad_L)) = number of singular values below tol.
    """
    d2 = L.shape[0]
    eye = np.eye(d2, dtype=complex)
    ad_L = np.kron(L, eye) - np.kron(eye, L.T)
    sv = np.linalg.svd(ad_L, compute_uv=False)
    dim = np.sum(sv < tol * sv[0]) if sv[0] > 0 else d2 * d2
    return dim, sv


def build_analytical_pi_product(N, pi_perm, pi_sign):
    """Build analytical product Pi = M_1 x M_2 x ... x M_N in standard vec basis.

    pi_perm: dict mapping Pauli index a -> perm(a) (per site)
    pi_sign: dict mapping Pauli index a -> phase factor (per site)

    Works by computing the per-site map in the 4x4 vec basis,
    then tensoring and permuting to standard vec ordering.
    """
    # Per-site map in vec basis: P1 maps |i><j| to P1(|i><j|)
    # Using Pauli decomposition: |i><j| = (1/2) sum_a (sigma_a)_{ji} sigma_a
    # P1(|i><j|) = (1/2) sum_a sign(a) (sigma_a)_{ji} sigma_{perm(a)}

    PAULI_MATS = [I2, sx, sy, sz]
    d_single = 2

    # Build per-site superoperator M_vec (4x4)
    M_vec = np.zeros((4, 4), dtype=complex)
    for col_idx in range(4):  # input: |i><j| with col_idx = 2*i + j
        i, j = col_idx // 2, col_idx % 2
        # Compute P1(|i><j|) using Pauli decomposition
        result = np.zeros((2, 2), dtype=complex)
        for a in range(4):
            coeff = PAULI_MATS[a][j, i]  # (sigma_a)_{ji}
            result += pi_sign[a] * coeff * 0.5 * PAULI_MATS[pi_perm[a]]
        # Convert result to vec
        for row_idx in range(4):
            k, l = row_idx // 2, row_idx % 2
            M_vec[row_idx, col_idx] = result[k, l]

    # For N qubits: Pi_grouped = M_vec^{otimes N} in A1 x A2 x ... grouped basis
    Pi_grouped = M_vec.copy()
    for _ in range(N - 1):
        Pi_grouped = np.kron(Pi_grouped, M_vec)

    # Permute to standard vec basis
    d = 2**N
    d2 = d * d

    # Build permutation: grouped -> standard
    # grouped index: (alpha_1, alpha_2, ..., alpha_N) where alpha_k = 2*i_k + j_k
    # standard index: k = d * sum(i_k * 2^(N-1-k)) + sum(j_k * 2^(N-1-k))
    perm_to_std = np.zeros(d2, dtype=int)
    for m in range(d2):
        # Decode grouped index into (alpha_1, ..., alpha_N)
        alphas = []
        temp = m
        for _ in range(N):
            alphas.append(temp % 4)
            temp //= 4
        alphas.reverse()

        # Decode alphas into (i_k, j_k) pairs
        i_bits = [a // 2 for a in alphas]
        j_bits = [a % 2 for a in alphas]

        # Compute standard indices
        i_val = sum(b * (2**(N - 1 - k)) for k, b in enumerate(i_bits))
        j_val = sum(b * (2**(N - 1 - k)) for k, b in enumerate(j_bits))
        k_std = d * i_val + j_val
        perm_to_std[m] = k_std

    # Build permutation matrix P: P[k_std, m_grouped] = 1
    P = np.zeros((d2, d2), dtype=complex)
    for m_idx in range(d2):
        P[perm_to_std[m_idx], m_idx] = 1.0

    # Pi in standard vec basis
    Pi_std = P @ Pi_grouped @ P.conj().T

    return Pi_std, M_vec


def analyze_with_analytical(name, N, terms, pi_perm, pi_sign,
                            topology='chain', J=1.0, gamma=0.05):
    """Analysis with both analytical (product) and numerical Pi."""
    print(f"\n{'=' * 60}")
    print(f"Case: {name}")
    print(f"N={N}, J={J}, gamma={gamma}")
    print('=' * 60)

    # Build system
    H = build_hamiltonian(N, terms, J=J, topology=topology)
    c_ops, sum_gamma = build_z_dephasing(N, gamma)
    L = build_liouvillian(H, c_ops)

    # --- Analytical product Pi ---
    Pi_anal, M_vec = build_analytical_pi_product(N, pi_perm, pi_sign)
    err_anal = verify_pi(Pi_anal, L, sum_gamma)
    print(f"\n  ANALYTICAL Pi (product M^{{xN}}):")
    print(f"    Conjugation equation error: {err_anal:.2e}")

    if err_anal < 1e-6:
        print(f"    VERIFIED: product Pi satisfies the conjugation equation")
        if N == 2:
            perm = build_bipartition_perm_2q()
            S_anal = operator_schmidt(Pi_anal, perm, 4, 4)
        elif N == 3:
            perm = build_bipartition_perm_3q_A_BC()
            S_anal = operator_schmidt(Pi_anal, perm, 4, 16)
        rank_anal = np.sum(S_anal > 1e-10 * S_anal[0])
        print(f"    Schmidt rank: {rank_anal}")
        if rank_anal == 1:
            print(f"    >>> Pi is a PRODUCT operator (rank 1)")
            print(f"    >>> The palindromic symmetry is LOCAL")
        else:
            print(f"    >>> Pi has rank {rank_anal} (unexpected for product!)")
            for i in range(min(rank_anal + 1, len(S_anal))):
                print(f"      sigma_{i} = {S_anal[i]:.6f}")
    else:
        print(f"    FAILED: product Pi does NOT work for this Hamiltonian")
        print(f"    The palindromic symmetry REQUIRES non-local Pi")

        # Fall back to numerical Pi
        Pi_num, pairs, evals = find_pi_operator(L, sum_gamma)
        if Pi_num is not None:
            err_num = verify_pi(Pi_num, L, sum_gamma)
            print(f"\n  NUMERICAL Pi (from eigenvalue pairing):")
            print(f"    Conjugation equation error: {err_num:.2e}")
            if N == 2:
                perm = build_bipartition_perm_2q()
                S_num = operator_schmidt(Pi_num, perm, 4, 4)
            elif N == 3:
                perm = build_bipartition_perm_3q_A_BC()
                S_num = operator_schmidt(Pi_num, perm, 4, 16)
            rank_num = np.sum(S_num > 1e-10 * S_num[0])
            print(f"    Schmidt rank: {rank_num} (upper bound on minimum)")
            print(f"    Schmidt coefficients:")
            for i in range(min(rank_num, len(S_num))):
                sn = S_num[i] / S_num[0]
                print(f"      sigma_{i} = {S_num[i]:.6f}  (relative: {sn:.6f})")

            S_prob = (S_num[:rank_num]**2) / np.sum(S_num[:rank_num]**2)
            ent = -np.sum(S_prob * np.log2(S_prob + 1e-30))
            print(f"    Operator entanglement entropy: {ent:.4f} bits")
            print(f"    >>> Pi is ENTANGLED in operator space")
            print(f"    >>> The palindromic symmetry is genuinely NON-LOCAL")


# === Main ===
if __name__ == '__main__':
    print("Q3: Entanglement Structure of the Pi Operator")
    print("Is the palindromic mirror entangled in operator space?")
    print()
    print("Method: build ANALYTICAL product Pi (per-site map tensored)")
    print("If it satisfies Pi*L*Pi^{-1} = -L - 2Sg*I, the symmetry is LOCAL.")
    print("If not, no per-site Pi exists and the symmetry is NON-LOCAL.")

    # Known per-site maps (from pi_time_reversal_verify.py)
    P1_PERM = {0: 1, 1: 0, 2: 3, 3: 2}   # I<->X, Y<->Z
    P1_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

    P4_PERM = {0: 2, 1: 3, 2: 0, 3: 1}   # I<->Y, X<->Z
    P4_SIGN = {0: 1, 1: 1j, 2: 1, 3: 1j}

    # --- Case 1: Heisenberg N=2 (P1 family, should be product) ---
    analyze_with_analytical(
        "N=2 Heisenberg (XX+YY+ZZ) + Z-deph [P1 family]",
        N=2, terms=[(sx, sx), (sy, sy), (sz, sz)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 2: N=2 Ising ZZ (P1 family) ---
    analyze_with_analytical(
        "N=2 Ising (ZZ) + Z-deph [P1 family]",
        N=2, terms=[(sz, sz)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 3: N=2 XX (P1 family) ---
    analyze_with_analytical(
        "N=2 XX + Z-deph [P1 family]",
        N=2, terms=[(sx, sx)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 4: N=2 XZ+ZY (known non-local, test P1 -> should FAIL) ---
    analyze_with_analytical(
        "N=2 XZ+ZY + Z-deph [test P1 - should fail]",
        N=2, terms=[(sx, sz), (sz, sy)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 5: N=2 YZ+ZX (known non-local, test P1 -> should FAIL) ---
    analyze_with_analytical(
        "N=2 YZ+ZX + Z-deph [test P1 - should fail]",
        N=2, terms=[(sy, sz), (sz, sx)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 6: N=2 XY (alternating Pi, test P1 -> might fail) ---
    analyze_with_analytical(
        "N=2 XY + Z-deph [test P1]",
        N=2, terms=[(sx, sy)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
    )

    # --- Case 7: N=3 Heisenberg chain (P1 family) ---
    analyze_with_analytical(
        "N=3 Heisenberg chain (XX+YY+ZZ) + Z-deph [P1 family]",
        N=3, terms=[(sx, sx), (sy, sy), (sz, sz)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
        topology='chain',
    )

    # --- Case 8: N=3 Ising chain (P1 family) ---
    analyze_with_analytical(
        "N=3 Ising chain (ZZ) + Z-deph [P1 family]",
        N=3, terms=[(sz, sz)],
        pi_perm=P1_PERM, pi_sign=P1_SIGN,
        topology='chain',
    )

    # ==================================================================
    # Q3.3: Alternating Pi for XY coupling
    # ==================================================================
    print("\n" + "=" * 60)
    print("Q3.3: ALTERNATING Pi FOR XY COUPLING")
    print("=" * 60)

    # M2 from algebraic_pi_search: I->+1Y, X->-iZ, Y->+1I, Z->-iX
    M2_PERM = {0: 2, 1: 3, 2: 0, 3: 1}   # I->Y, X->Z, Y->I, Z->X
    M2_SIGN = {0: 1, 1: -1j, 2: 1, 3: -1j}

    # N=2 XY: site 0 = P1, site 1 = M2
    print("\n  N=2 XY + Z-deph [alternating P1 x M2]:")
    H = build_hamiltonian(2, [(sx, sy)], J=1.0)
    c_ops, sum_gamma = build_z_dephasing(2, 0.05)
    L = build_liouvillian(H, c_ops)
    Pi_alt = build_analytical_pi_multi(2, [
        (P1_PERM, P1_SIGN),
        (M2_PERM, M2_SIGN),
    ])
    err = verify_pi(Pi_alt, L, sum_gamma)
    print(f"    Conjugation equation error: {err:.2e}")
    if err < 1e-6:
        perm = build_bipartition_perm_2q()
        S = operator_schmidt(Pi_alt, perm, 4, 4)
        rank = np.sum(S > 1e-10 * S[0])
        print(f"    Schmidt rank: {rank}")
        print(f"    >>> VERIFIED: alternating Pi is {'PRODUCT (LOCAL)' if rank == 1 else f'rank {rank}'}")
    else:
        print(f"    FAILED (error {err:.2e})")

    # Also test M2 x P1 (reversed)
    print("\n  N=2 XY + Z-deph [alternating M2 x P1]:")
    Pi_alt2 = build_analytical_pi_multi(2, [
        (M2_PERM, M2_SIGN),
        (P1_PERM, P1_SIGN),
    ])
    err2 = verify_pi(Pi_alt2, L, sum_gamma)
    print(f"    Conjugation equation error: {err2:.2e}")
    if err2 < 1e-6:
        S2 = operator_schmidt(Pi_alt2, perm, 4, 4)
        rank2 = np.sum(S2 > 1e-10 * S2[0])
        print(f"    Schmidt rank: {rank2}")
        print(f"    >>> VERIFIED: reversed alternating Pi is {'PRODUCT (LOCAL)' if rank2 == 1 else f'rank {rank2}'}")
    else:
        print(f"    FAILED (error {err2:.2e})")

    # Also test P4 for XZ+ZY (should also fail)
    print("\n  N=2 XZ+ZY + Z-deph [test P4 - should also fail]:")
    H_xz_zy = build_hamiltonian(2, [(sx, sz), (sz, sy)], J=1.0)
    L_xz_zy = build_liouvillian(H_xz_zy, c_ops)
    Pi_p4 = build_analytical_pi_multi(2, [
        (P4_PERM, P4_SIGN),
        (P4_PERM, P4_SIGN),
    ])
    err_p4 = verify_pi(Pi_p4, L_xz_zy, sum_gamma)
    print(f"    P4 x P4 error: {err_p4:.2e} {'FAIL' if err_p4 > 1e-6 else 'OK'}")
    Pi_p1p4 = build_analytical_pi_multi(2, [(P1_PERM, P1_SIGN), (P4_PERM, P4_SIGN)])
    err_p1p4 = verify_pi(Pi_p1p4, L_xz_zy, sum_gamma)
    print(f"    P1 x P4 error: {err_p1p4:.2e} {'FAIL' if err_p1p4 > 1e-6 else 'OK'}")
    Pi_p4p1 = build_analytical_pi_multi(2, [(P4_PERM, P4_SIGN), (P1_PERM, P1_SIGN)])
    err_p4p1 = verify_pi(Pi_p4p1, L_xz_zy, sum_gamma)
    print(f"    P4 x P1 error: {err_p4p1:.2e} {'FAIL' if err_p4p1 > 1e-6 else 'OK'}")
    Pi_m2m2 = build_analytical_pi_multi(2, [(M2_PERM, M2_SIGN), (M2_PERM, M2_SIGN)])
    err_m2m2 = verify_pi(Pi_m2m2, L_xz_zy, sum_gamma)
    print(f"    M2 x M2 error: {err_m2m2:.2e} {'FAIL' if err_m2m2 > 1e-6 else 'OK'}")
    if err_p4p1 < 1e-6:
        S_p4p1 = operator_schmidt(Pi_p4p1, perm, 4, 4)
        rank_p4p1 = np.sum(S_p4p1 > 1e-10 * S_p4p1[0])
        print(f"    !!! SURPRISE: P4 x P1 WORKS for XZ+ZY!")
        print(f"    !!! Schmidt rank: {rank_p4p1}")
        print(f"    !!! XZ+ZY is LOCAL with non-uniform per-site maps!")
    else:
        print(f"    >>> ALL per-site products fail for XZ+ZY: genuinely non-local")

    # Test YZ+ZX with all per-site combinations
    print(f"\n  N=2 YZ+ZX + Z-deph [testing all per-site combos]:")
    H_yz_zx = build_hamiltonian(2, [(sy, sz), (sz, sx)], J=1.0)
    L_yz_zx = build_liouvillian(H_yz_zx, c_ops)
    combos = [
        ("P1 x P1", P1_PERM, P1_SIGN, P1_PERM, P1_SIGN),
        ("P4 x P4", P4_PERM, P4_SIGN, P4_PERM, P4_SIGN),
        ("P1 x P4", P1_PERM, P1_SIGN, P4_PERM, P4_SIGN),
        ("P4 x P1", P4_PERM, P4_SIGN, P1_PERM, P1_SIGN),
        ("P1 x M2", P1_PERM, P1_SIGN, M2_PERM, M2_SIGN),
        ("M2 x P1", M2_PERM, M2_SIGN, P1_PERM, P1_SIGN),
        ("P4 x M2", P4_PERM, P4_SIGN, M2_PERM, M2_SIGN),
        ("M2 x P4", M2_PERM, M2_SIGN, P4_PERM, P4_SIGN),
    ]
    found_product = False
    for label, p1, s1, p2, s2 in combos:
        Pi_test = build_analytical_pi_multi(2, [(p1, s1), (p2, s2)])
        err_test = verify_pi(Pi_test, L_yz_zx, sum_gamma)
        status = "OK" if err_test < 1e-6 else "FAIL"
        print(f"    {label:10s} error: {err_test:.2e}  {status}")
        if err_test < 1e-6:
            S_test = operator_schmidt(Pi_test, perm, 4, 4)
            rank_test = np.sum(S_test > 1e-10 * S_test[0])
            print(f"               Schmidt rank: {rank_test} -> {'LOCAL' if rank_test == 1 else 'ENTANGLED'}")
            found_product = True
    if found_product:
        print(f"    >>> SURPRISE: YZ+ZX also has a product Pi!")
    else:
        print(f"    >>> YZ+ZX: no product found among tested combinations")

    # ==================================================================
    # CRITICAL CHECK: Test the ACTUAL cases from NON_HEISENBERG_PALINDROME
    # XZ+YZ = sigma_X x sigma_Z + sigma_Y x sigma_Z  (X and Y on SAME site)
    # ZX+ZY = sigma_Z x sigma_X + sigma_Z x sigma_Y  (X and Y on SAME site)
    # These are DIFFERENT from XZ+ZY and YZ+ZX tested above!
    # ==================================================================
    print(f"\n  --- CRITICAL: Testing actual NON_HEISENBERG cases ---")

    for h_name, h_terms in [
        ("XZ+YZ", [(sx, sz), (sy, sz)]),
        ("ZX+ZY", [(sz, sx), (sz, sy)]),
    ]:
        print(f"\n  N=2 {h_name} + Z-deph [exhaustive per-site test]:")
        H_test = build_hamiltonian(2, h_terms, J=1.0)
        L_test = build_liouvillian(H_test, c_ops)
        found = False
        for label, p1, s1, p2, s2 in combos:
            Pi_t = build_analytical_pi_multi(2, [(p1, s1), (p2, s2)])
            err_t = verify_pi(Pi_t, L_test, sum_gamma)
            status = "OK" if err_t < 1e-6 else ""
            if err_t < 1e-6:
                S_t = operator_schmidt(Pi_t, perm, 4, 4)
                r_t = np.sum(S_t > 1e-10 * S_t[0])
                print(f"    {label:10s} error: {err_t:.2e}  OK  rank={r_t}")
                found = True
            # Only print failures in compact form
        if not found:
            print(f"    All 8 per-site combinations FAIL.")
            print(f"    >>> {h_name} is genuinely NON-LOCAL (confirmed)")
        else:
            print(f"    >>> {h_name} has a product Pi (LOCAL)")

    # ==================================================================
    # Q3.4: Commutant dimension (space of valid Pi operators)
    # ==================================================================
    print("\n" + "=" * 60)
    print("Q3.4: DIMENSION OF VALID-Pi SPACE")
    print("=" * 60)
    print("\n  Space of valid Pi = Pi_0 * {S : [S,L]=0}")
    print("  dim(valid Pi space) = dim(commutant of L)")

    test_cases = [
        ("N=2 Heisenberg", [(sx, sx), (sy, sy), (sz, sz)], 2),
        ("N=2 Ising ZZ", [(sz, sz)], 2),
        ("N=2 XX", [(sx, sx)], 2),
        ("N=2 XY", [(sx, sy)], 2),
        ("N=2 XZ+ZY", [(sx, sz), (sz, sy)], 2),
        ("N=2 YZ+ZX", [(sy, sz), (sz, sx)], 2),
    ]

    for name, terms, N in test_cases:
        H = build_hamiltonian(N, terms, J=1.0)
        c_ops, sg = build_z_dephasing(N, 0.05)
        L = build_liouvillian(H, c_ops)
        dim, sv = commutant_dimension(L)
        # Count distinct eigenvalues to understand degeneracy
        evals = np.linalg.eigvals(L)
        n_distinct = len(set(np.round(evals, 8)))
        print(f"  {name:25s}  commutant dim = {dim:3d}  "
              f"(L has {n_distinct} distinct eigenvalues of {len(evals)})")

    # --- Summary ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("""
The question is not "is the numerical Pi entangled?" (it always can be,
because the space of valid Pi operators includes entangled ones when
eigenspaces are degenerate). The question is: "does a PRODUCT Pi exist?"

If a product Pi exists (passes conjugation equation with per-site M^{xN}):
  The palindromic symmetry is LOCAL. Each qubit has its own
  independent palindromic structure.

If NO product Pi exists (all per-site M fail):
  The palindromic symmetry is genuinely NON-LOCAL. It cannot be
  understood as independent per-site symmetries. Implications:
  - Error correction using the palindrome requires entangled
    measurements across qubits
  - The palindrome has a genuinely multi-body character
  - A topological invariant (Q2) becomes plausible
""")
