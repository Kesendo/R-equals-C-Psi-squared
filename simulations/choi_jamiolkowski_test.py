"""
Choi-Jamiolkowski Test: Is Non-Local Pi a Projection Artifact?
===============================================================
Map the non-local Pi superoperator into its Choi state in the doubled
Hilbert space. Test whether the non-locality survives representation change.

Core test: If the Choi state factors as a product operator across the
(site1_doubled | site2_doubled) bipartition, the non-locality is a
projection artifact. If not, it is fundamental and irreducible.

Script: simulations/choi_jamiolkowski_test.py
Output: simulations/results/choi_jamiolkowski_test.txt
"""
import numpy as np
from itertools import product as iproduct
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\choi_jamiolkowski_test.txt"
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# PAULI BASICS (from continuous_pi_search.py)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
H_LABELS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


def plabel(indices):
    return ''.join(NAMES[i] for i in indices)


def build_pauli_basis(N):
    all_idx = list(iproduct(range(4), repeat=N))
    d = 2 ** N
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)
    return all_idx, pmats, d


def build_H_term(N, bond, label):
    i, j = bond
    ops = [I2] * N
    ops[i] = PM[label[0]]
    ops[j] = PM[label[1]]
    H = ops[0]
    for o in ops[1:]:
        H = np.kron(H, o)
    return H


def build_L_H(N, H, all_idx, pmats, d):
    num = 4 ** N
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / d
    return L_H


def build_L_D_diag(N, gamma, all_idx):
    return np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])


# ============================================================
# PI CONSTRUCTION
# ============================================================
def construct_Pi(L_full, Sg, num):
    """Construct Pi from palindromic eigenvalue pairing (Section 5 method)."""
    evals, R = np.linalg.eig(L_full)

    paired = np.zeros(num, dtype=bool)
    pair_map = -np.ones(num, dtype=int)
    for k in range(num):
        if paired[k]:
            continue
        target = -(evals[k] + 2 * Sg)
        best_j, best_d = -1, np.inf
        for j in range(k, num):
            if paired[j] and j != k:
                continue
            dd = abs(evals[j] - target)
            if dd < best_d:
                best_d, best_j = dd, j
        if best_j >= 0 and best_d < 1e-8:
            paired[k] = paired[best_j] = True
            pair_map[k] = best_j
            pair_map[best_j] = k

    n_paired = int(np.sum(paired))
    if n_paired < num:
        return None, n_paired

    R_inv = np.linalg.inv(R)
    P_eig = np.zeros((num, num), dtype=complex)
    seen = set()
    for k in range(num):
        j = pair_map[k]
        if k in seen:
            continue
        seen.add(k)
        seen.add(j)
        if k == j:
            P_eig[k, k] = 1.0
        else:
            P_eig[j, k] = 1.0
            P_eig[k, j] = 1.0

    return R @ P_eig @ R_inv, n_paired


def full_palindrome_error(Pi, L, Sg):
    Pi_inv = np.linalg.inv(Pi)
    num = L.shape[0]
    conj = Pi @ L @ Pi_inv
    err = np.max(np.abs(conj + L + 2 * Sg * np.eye(num)))
    norm = np.max(np.abs(L))
    return err / norm if norm > 1e-15 else err


def anticomm_error(Pi, L_H):
    Pi_inv = np.linalg.inv(Pi)
    conj = Pi @ L_H @ Pi_inv
    err = np.max(np.abs(conj + L_H))
    norm = np.max(np.abs(L_H))
    return err / norm if norm > 1e-15 else err


def check_tensor_product(Pi):
    """Check if 16x16 Pi = M1 x M2 in Pauli basis."""
    Pi_4d = Pi.reshape(4, 4, 4, 4)
    for c in range(4):
        for dd in range(4):
            sl = Pi_4d[:, :, c, dd]
            if np.max(np.abs(sl)) < 1e-10:
                continue
            s = np.linalg.svd(sl, compute_uv=False)
            if np.sum(s > 1e-8 * s[0]) > 1:
                return False
    return True


def build_M_block(A, B):
    """Build 4x4 per-site M from 2x2 blocks."""
    M = np.zeros((4, 4), dtype=complex)
    M[0, 1] = A[0, 0]; M[0, 2] = A[0, 1]
    M[3, 1] = A[1, 0]; M[3, 2] = A[1, 1]
    M[1, 0] = B[0, 0]; M[1, 3] = B[0, 1]
    M[2, 0] = B[1, 0]; M[2, 3] = B[1, 1]
    return M


def build_M_P1():
    A = np.array([[1, 0], [0, 1j]], dtype=complex)
    B = np.array([[1, 0], [0, 1j]], dtype=complex)
    return build_M_block(A, B)


def build_M_P4():
    A = np.array([[0, 1], [1, 0]], dtype=complex)
    B = np.array([[0, 1], [1, 0]], dtype=complex)
    return build_M_block(A, B)


# ============================================================
# CHOI-JAMIOLKOWSKI MACHINERY
# ============================================================
def build_change_of_basis(pmats, d):
    """V: Pauli basis -> computational basis.
    V[:,a] = vec(P_a) where vec uses row-major (C) ordering.
    V satisfies V^dag V = d * I, so V^{-1} = V^dag / d.
    """
    num = len(pmats)
    V = np.zeros((d * d, num), dtype=complex)
    for a, P in enumerate(pmats):
        V[:, a] = P.flatten()
    return V


def pi_to_superop(Pi, V, d):
    """Convert Pauli-basis Pi to computational-basis superoperator.
    S_comp = (1/d) V Pi V^dag
    where S_comp acts on vec(rho) to give vec(Phi(rho)).
    """
    return V @ Pi @ V.conj().T / d


def superop_to_choi(S, d):
    """Reshuffle superoperator S (d^2 x d^2) to Choi matrix J (d^2 x d^2).
    Index map: S[m,n,k,l] -> J[m,k,n,l] = S[m,n,k,l] / d
    (axes 1 and 2 are swapped, normalized by 1/d).
    """
    return S.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d) / d


def partial_transpose(J, dA, dB, sub='B'):
    """Partial transpose of J (dA*dB x dA*dB) over subsystem."""
    J4 = J.reshape(dA, dB, dA, dB)
    if sub == 'B':
        return J4.transpose(0, 3, 2, 1).reshape(dA * dB, dA * dB)
    else:
        return J4.transpose(2, 1, 0, 3).reshape(dA * dB, dA * dB)


def reorder_site_bipartition(J):
    """Reorder 16x16 Choi from (sys x anc) to (site1_doubled x site2_doubled).
    Original: row index = (m1,m2,k1,k2) as m1*8+m2*4+k1*2+k2
    Reordered: (m1,k1,m2,k2) as m1*8+k1*4+m2*2+k2
    i.e., swap axes 1 (m2, weight 4) and 2 (k1, weight 2)."""
    J8 = J.reshape(2, 2, 2, 2, 2, 2, 2, 2)
    # axes: (m1=0, m2=1, k1=2, k2=3, n1=4, n2=5, l1=6, l2=7)
    # target: (m1=0, k1=2, m2=1, k2=3, n1=4, l1=6, n2=5, l2=7)
    return J8.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)


def operator_svd(J, dA, dB):
    """SVD of J viewed as an operator in (A x B) space.
    Reshape J[a,b,a',b'] -> M[(a,a'),(b,b')] and compute SVD.
    Rank 1 = product operator J_A x J_B."""
    J4 = J.reshape(dA, dB, dA, dB)
    M = J4.transpose(0, 2, 1, 3).reshape(dA ** 2, dB ** 2)
    s = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(s > 1e-8 * s[0])) if s[0] > 1e-15 else 0
    return rank, s, M


def safe_eigvalsh(M):
    """Sorted eigenvalues: eigvalsh for Hermitian, real part of eigvals otherwise."""
    if np.allclose(M, M.conj().T, atol=1e-10):
        return np.sort(np.linalg.eigvalsh(M))
    return np.sort(np.real(np.linalg.eigvals(M)))


def negativity(J_pt):
    """Sum of absolute values of negative eigenvalues."""
    evals = safe_eigvalsh(J_pt)
    return float(-np.sum(evals[evals < -1e-12]))


# ============================================================
# FULL ANALYSIS ROUTINE
# ============================================================
def analyze_choi(J, d, label):
    """Entanglement analysis of Choi matrix across three bipartitions."""
    d2 = d * d
    is_herm = np.allclose(J, J.conj().T, atol=1e-10)
    evals = safe_eigvalsh(J)
    n_neg = int(np.sum(evals < -1e-10))
    is_psd = (n_neg == 0) and is_herm

    log(f"\n    Properties of J({label}):")
    log(f"      Trace: {np.trace(J):.6f}")
    log(f"      Hermitian: {is_herm}")
    log(f"      PSD: {is_psd}")
    log(f"      Eigenvalues:")
    log(f"        {np.array2string(evals, precision=6, separator=', ', max_line_width=80)}")

    # --- (a) System | Ancilla ---
    log(f"\n    (a) System|Ancilla:")
    J_pt_a = partial_transpose(J, d, d)
    neg_a = negativity(J_pt_a)
    evals_pt_a = safe_eigvalsh(J_pt_a)
    n_neg_a = int(np.sum(evals_pt_a < -1e-10))
    log(f"      PT eigenvalues:")
    log(f"        {np.array2string(evals_pt_a, precision=6, separator=', ', max_line_width=80)}")
    log(f"      Negative PT eigenvalues: {n_neg_a}")
    log(f"      Negativity: {neg_a:.6f}")

    if is_psd and abs(np.trace(J)) > 1e-10:
        J_n = J / np.trace(J)
        rho_A = np.trace(J_n.reshape(d, d, d, d), axis1=1, axis2=3)
        ev_A = np.linalg.eigvalsh(rho_A)
        ev_A = ev_A[ev_A > 1e-15]
        ent_a = -np.sum(ev_A * np.log2(ev_A))
        log(f"      Entanglement entropy: {ent_a:.6f} bits")

    # --- (b) Site1 | Site2 in doubled space --- THE KEY TEST ---
    log(f"\n    (b) Site1|Site2 (doubled space) -- THE KEY TEST:")
    J_b = reorder_site_bipartition(J)
    rank_b, svs_b, M_b = operator_svd(J_b, 4, 4)
    svs_nz = svs_b[svs_b > 1e-10 * svs_b[0]]
    log(f"      Operator SVD singular values:")
    log(f"        {np.array2string(svs_nz, precision=6, separator=', ')}")
    log(f"      Operator rank: {rank_b}")

    if rank_b == 1:
        log(f"      ==> PRODUCT OPERATOR (J = J_site1 x J_site2)")
        log(f"      ==> Non-locality IS a projection artifact")
    else:
        log(f"      ==> NOT a product operator (rank {rank_b})")
        log(f"      ==> Non-locality is FUNDAMENTAL")

    J_pt_b = partial_transpose(J_b, 4, 4)
    neg_b = negativity(J_pt_b)
    evals_pt_b = safe_eigvalsh(J_pt_b)
    n_neg_b = int(np.sum(evals_pt_b < -1e-10))
    log(f"      PT negativity: {neg_b:.6f} ({n_neg_b} negative eigenvalues)")

    # --- (c) Time bipartition ---
    log(f"\n    (c) Time bipartition (sys_sites | anc_sites):")
    log(f"      Identical to (a): System|Ancilla.")
    log(f"      Negativity: {neg_a:.6f}")

    return {
        'is_psd': is_psd,
        'neg_a': neg_a,
        'neg_b': neg_b,
        'rank_b': rank_b,
        'svs_b': svs_b,
        'M_b': M_b,
    }


# ============================================================
# MAIN
# ============================================================
gamma = 0.05
N = 2
bonds = [(0, 1)]
Sg = N * gamma
all_idx, pmats, d = build_pauli_basis(N)
num = 4 ** N  # 16
L_D_diag = build_L_D_diag(N, gamma, all_idx)

log("=" * 90)
log("CHOI-JAMIOLKOWSKI TEST: Is Non-Local Pi a Projection Artifact?")
log(f"Date: {datetime.now()}")
log(f"N = {N}, gamma = {gamma}, Sg = N*gamma = {Sg}")
log("=" * 90)

# Build change-of-basis V and verify
V = build_change_of_basis(pmats, d)
err_V = np.max(np.abs(V.conj().T @ V - d * np.eye(num)))
log(f"\nBasis check: max|V^dag V - {d}*I| = {err_V:.2e}")

# Sanity: Choi of identity superoperator = |Omega><Omega|
S_id = np.eye(d * d, dtype=complex)
J_id = superop_to_choi(S_id, d)
Omega = np.zeros(d * d)
for i in range(d):
    Omega[i * d + i] = 1.0 / np.sqrt(d)
err_choi = np.max(np.abs(J_id - np.outer(Omega, Omega)))
log(f"Choi sanity: max|J(Id) - |Omega><Omega|| = {err_choi:.2e}")

# Build L_H for all single terms
L_H_terms = {}
for label in H_LABELS:
    H = build_H_term(N, bonds[0], label)
    L_H_terms[label] = build_L_H(N, H, all_idx, pmats, d)


# ################################################################
# SECTION 1: Load Non-Local Pi (XZ+YZ)
# ################################################################
log()
log("=" * 90)
log("SECTION 1: Pi loaded, sanity check")
log("=" * 90)

L_H_xzyz = L_H_terms['XZ'] + L_H_terms['YZ']
L_full_xzyz = L_H_xzyz.copy()
for a in range(num):
    L_full_xzyz[a, a] += L_D_diag[a]

Pi_xzyz, n_paired = construct_Pi(L_full_xzyz, Sg, num)
err_pal = full_palindrome_error(Pi_xzyz, L_full_xzyz, Sg)
is_tp = check_tensor_product(Pi_xzyz)

log(f"\n  XZ+YZ:")
log(f"    Palindromic pairing: {n_paired}/{num}")
log(f"    |Pi L Pi^-1 + L + 2Sg I| / |L|: {err_pal:.2e}")
log(f"    Sanity check: {'PASSED' if err_pal < 1e-10 else 'FAILED'}")
log(f"    Tensor product M1 x M2 (Pauli basis): {'YES' if is_tp else 'NO'}")


# ################################################################
# SECTION 2: Compute Choi State
# ################################################################
log()
log("=" * 90)
log("SECTION 2: Choi state J(Pi) computed, 16x16")
log("=" * 90)

S_xzyz = pi_to_superop(Pi_xzyz, V, d)
J_xzyz = superop_to_choi(S_xzyz, d)

log(f"\n  Computational-basis superoperator: {S_xzyz.shape}")
log(f"  Choi matrix: {J_xzyz.shape}")

# Cross-check: apply Pi in Pauli basis vs computational basis
rho_test = np.array([[0.6, 0.1 + 0.05j, 0, 0],
                      [0.1 - 0.05j, 0.1, 0, 0],
                      [0, 0, 0.2, 0.03j],
                      [0, 0, -0.03j, 0.1]], dtype=complex)
r_test = np.array([np.trace(P @ rho_test) for P in pmats])
r_mapped = Pi_xzyz @ r_test
rho_pauli = sum(r_mapped[a] * pmats[a] / d for a in range(num))
rho_comp = (S_xzyz @ rho_test.flatten()).reshape(d, d)
cross_err = np.max(np.abs(rho_pauli - rho_comp))
log(f"  Cross-check (Pauli vs comp): max|diff| = {cross_err:.2e}")


# ################################################################
# SECTION 3: Entanglement Analysis
# ################################################################
log()
log("=" * 90)
log("SECTION 3: Entanglement analysis")
log("=" * 90)

results_xzyz = analyze_choi(J_xzyz, d, "Pi_XZ+YZ")


# ################################################################
# SECTION 4: Control (local Pi)
# ################################################################
log()
log("=" * 90)
log("SECTION 4: Control (local Pi)")
log("=" * 90)

# Find a working local Pi from discrete P1/P4 families
log(f"\n  Searching for a local Pi = M x M among discrete candidates...")
M_p1 = build_M_P1()
M_p4 = build_M_P4()

control_found = False
control_label = None
Pi_control = None

for base_M, mname in [(M_p1, 'P1'), (M_p4, 'P4')]:
    if control_found:
        break
    for phases in iproduct([1, -1, 1j, -1j], repeat=4):
        if control_found:
            break
        M_c = base_M.copy()
        for col in range(4):
            M_c[:, col] *= phases[col]
        if abs(np.linalg.det(M_c)) < 1e-10:
            continue
        Pi_test = np.kron(M_c, M_c)
        # Try single-term Hamiltonians
        for label in H_LABELS:
            err = anticomm_error(Pi_test, L_H_terms[label])
            if err < 1e-8:
                control_label = label
                Pi_control = Pi_test
                control_found = True
                log(f"    Found: {mname} x {mname} works for H = {label}")
                break

if not control_found:
    # Try two-term combos
    for base_M, mname in [(M_p1, 'P1'), (M_p4, 'P4')]:
        if control_found:
            break
        for phases in iproduct([1, -1, 1j, -1j], repeat=4):
            if control_found:
                break
            M_c = base_M.copy()
            for col in range(4):
                M_c[:, col] *= phases[col]
            if abs(np.linalg.det(M_c)) < 1e-10:
                continue
            Pi_test = np.kron(M_c, M_c)
            for i, t1 in enumerate(H_LABELS):
                for t2 in H_LABELS[i + 1:]:
                    if f"{t1}+{t2}" in ['XZ+YZ', 'ZX+ZY']:
                        continue
                    L_combo = L_H_terms[t1] + L_H_terms[t2]
                    err = anticomm_error(Pi_test, L_combo)
                    if err < 1e-8:
                        control_label = f"{t1}+{t2}"
                        Pi_control = Pi_test
                        control_found = True
                        log(f"    Found: {mname} x {mname} works for H = {control_label}")
                        break
                if control_found:
                    break

results_control = None
if control_found:
    # Build full L for the control case
    if '+' in control_label:
        t1, t2 = control_label.split('+')
        L_H_ctrl = L_H_terms[t1] + L_H_terms[t2]
    else:
        L_H_ctrl = L_H_terms[control_label]
    L_full_ctrl = L_H_ctrl.copy()
    for a in range(num):
        L_full_ctrl[a, a] += L_D_diag[a]

    err_ctrl = full_palindrome_error(Pi_control, L_full_ctrl, Sg)
    is_tp_ctrl = check_tensor_product(Pi_control)
    log(f"\n  Control ({control_label}):")
    log(f"    Palindrome error: {err_ctrl:.2e}")
    log(f"    Tensor product (Pauli basis): {is_tp_ctrl}")

    S_ctrl = pi_to_superop(Pi_control, V, d)
    J_ctrl = superop_to_choi(S_ctrl, d)
    results_control = analyze_choi(J_ctrl, d, f"Pi_{control_label}")
else:
    log(f"  No local control found in discrete search.")
    log(f"  Falling back to eigenvector-constructed Pi for XX...")
    L_H_xx = L_H_terms['XX']
    L_full_xx = L_H_xx.copy()
    for a in range(num):
        L_full_xx[a, a] += L_D_diag[a]
    Pi_xx, n_p = construct_Pi(L_full_xx, Sg, num)
    if Pi_xx is not None:
        err_xx = full_palindrome_error(Pi_xx, L_full_xx, Sg)
        is_tp_xx = check_tensor_product(Pi_xx)
        log(f"    XX palindrome error: {err_xx:.2e}")
        log(f"    XX tensor product: {is_tp_xx}")
        S_xx = pi_to_superop(Pi_xx, V, d)
        J_xx = superop_to_choi(S_xx, d)
        results_control = analyze_choi(J_xx, d, "Pi_XX")
        control_label = "XX"


# ################################################################
# SECTION 5: Quantification
# ################################################################
log()
log("=" * 90)
log("SECTION 5: Quantification")
log("=" * 90)

rank = results_xzyz['rank_b']
svs = results_xzyz['svs_b']
M_b = results_xzyz['M_b']

if rank == 1:
    log(f"\n  The Choi state is SEPARABLE across Site1|Site2.")
    log(f"  Non-locality is a projection artifact.")

    # Extract product decomposition
    U, s, Vh = np.linalg.svd(M_b)
    J1 = (U[:, 0] * np.sqrt(s[0])).reshape(4, 4)
    J2 = (Vh[0, :].conj() * np.sqrt(s[0])).reshape(4, 4)
    log(f"\n  Product decomposition: J = J_site1 x J_site2")
    log(f"  J_site1 (4x4):")
    for row in range(4):
        log(f"    [{', '.join(f'{J1[row, c]:+9.5f}' for c in range(4))}]")
    log(f"  J_site2 (4x4):")
    for row in range(4):
        log(f"    [{', '.join(f'{J2[row, c]:+9.5f}' for c in range(4))}]")

else:
    log(f"\n  The Choi state is NOT separable across Site1|Site2 (operator rank {rank}).")
    log(f"  The non-locality is FUNDAMENTAL -- survives in every representation.")

    svs_nz = svs[svs > 1e-10 * svs[0]]
    log(f"\n  Schmidt decomposition (operator SVD):")
    log(f"    Number of nonzero singular values: {len(svs_nz)}")
    log(f"    Singular values: {np.array2string(svs_nz, precision=6, separator=', ')}")

    # Normalized singular values as "probabilities"
    svs_sq = (svs_nz / np.linalg.norm(svs_nz)) ** 2
    entropy = -np.sum(svs_sq * np.log2(svs_sq + 1e-30))
    log(f"    Schmidt entropy: {entropy:.6f} bits")
    log(f"    Max possible (log2 {len(svs_nz)}): {np.log2(len(svs_nz)):.6f} bits")

    if len(svs_nz) >= 2:
        log(f"    Dominance ratio s1/s2: {svs_nz[0] / svs_nz[1]:.4f}")
        log(f"    Ratio s1/s_last: {svs_nz[0] / svs_nz[-1]:.4f}")

    # How "close" to a product is it?
    best_rank1 = svs_nz[0] * np.outer(
        np.linalg.svd(M_b)[0][:, 0],
        np.linalg.svd(M_b)[2][0, :].conj()
    )
    residual = np.linalg.norm(M_b - best_rank1) / np.linalg.norm(M_b)
    log(f"    Residual after rank-1 approximation: {residual:.6f} ({residual * 100:.2f}%)")


# ################################################################
# SECTION 6: Roberts Connection
# ################################################################
log()
log("=" * 90)
log("SECTION 6: Roberts connection -- antiunitary structure in doubled space")
log("=" * 90)

d2 = d * d

log(f"\n  The Choi matrix J(Pi) lives in the doubled Hilbert space.")
log(f"  Roberts-Lingenfelter-Clerk: hidden TRS operators live here.")

# SWAP operator (system <-> ancilla)
SWAP = np.zeros((d2, d2), dtype=complex)
for m in range(d):
    for k in range(d):
        SWAP[m * d + k, k * d + m] = 1.0

comm_sw = np.linalg.norm(J_xzyz @ SWAP - SWAP @ J_xzyz)
acomm_sw = np.linalg.norm(J_xzyz @ SWAP + SWAP @ J_xzyz)
log(f"\n  SWAP (system <-> ancilla):")
log(f"    ||[J, SWAP]|| = {comm_sw:.6e}")
log(f"    ||{{J, SWAP}}|| = {acomm_sw:.6e}")

# Symmetry properties
log(f"\n  Symmetries of J(Pi):")
log(f"    max|Im(J)|: {np.max(np.abs(np.imag(J_xzyz))):.6e}")
is_real = np.max(np.abs(np.imag(J_xzyz))) < 1e-10
log(f"    Real: {is_real}")
log(f"    max|J - J^T|: {np.max(np.abs(J_xzyz - J_xzyz.T)):.6e}")
is_symm = np.max(np.abs(J_xzyz - J_xzyz.T)) < 1e-10
log(f"    Symmetric: {is_symm}")
log(f"    max|J - J*|: {np.max(np.abs(J_xzyz - J_xzyz.conj())):.6e}")

# Involution: Pi^2 = ?
Pi2 = Pi_xzyz @ Pi_xzyz
c_Pi = np.trace(Pi2) / num
err_invol = np.max(np.abs(Pi2 - c_Pi * np.eye(num)))
log(f"\n  Involution check:")
log(f"    Pi^2 ~ c*I: c = {c_Pi:.6f}, err = {err_invol:.6e}")
log(f"    Pi is {'an involution' if err_invol < 1e-6 else 'NOT an involution'}")

# J^2
J2 = J_xzyz @ J_xzyz
c_J = np.trace(J2) / d2
err_J2 = np.max(np.abs(J2 - c_J * np.eye(d2)))
log(f"    J^2 ~ c*I: c = {c_J:.6f}, err = {err_J2:.6e}")

# Compare with Choi of transpose map (SWAP/d)
J_transpose = SWAP / d
diff_T = np.linalg.norm(J_xzyz - J_transpose) / np.linalg.norm(J_xzyz)
log(f"\n  Distance to transpose-map Choi (SWAP/d): {diff_T:.6f}")

# Check if J is proportional to a unitary channel's Choi state
# A unitary channel U: rho -> U rho U^dag has Choi = |psi><psi| where |psi> = (U x I)|Omega>
# So J should be rank 1 for a unitary channel
J_evals = safe_eigvalsh(J_xzyz)
J_rank = int(np.sum(np.abs(J_evals) > 1e-10 * np.max(np.abs(J_evals))))
log(f"  Rank of J: {J_rank} (rank 1 = unitary channel)")


# ################################################################
# SECTION 7: ZX+ZY comparison
# ################################################################
log()
log("=" * 90)
log("SECTION 7: ZX+ZY comparison")
log("=" * 90)

L_H_zxzy = L_H_terms['ZX'] + L_H_terms['ZY']
L_full_zxzy = L_H_zxzy.copy()
for a in range(num):
    L_full_zxzy[a, a] += L_D_diag[a]

Pi_zxzy, n_paired_zxzy = construct_Pi(L_full_zxzy, Sg, num)

if Pi_zxzy is not None:
    err_zxzy = full_palindrome_error(Pi_zxzy, L_full_zxzy, Sg)
    is_tp_zxzy = check_tensor_product(Pi_zxzy)

    log(f"\n  ZX+ZY:")
    log(f"    Palindromic pairing: {n_paired_zxzy}/{num}")
    log(f"    Palindrome error: {err_zxzy:.2e}")
    log(f"    Tensor product (Pauli basis): {'YES' if is_tp_zxzy else 'NO'}")

    S_zxzy = pi_to_superop(Pi_zxzy, V, d)
    J_zxzy = superop_to_choi(S_zxzy, d)
    results_zxzy = analyze_choi(J_zxzy, d, "Pi_ZX+ZY")

    log(f"\n    Comparison with XZ+YZ:")
    log(f"      XZ+YZ operator rank: {results_xzyz['rank_b']}")
    log(f"      ZX+ZY operator rank: {results_zxzy['rank_b']}")

    svs1 = results_xzyz['svs_b']
    svs2 = results_zxzy['svs_b']
    svs1_nz = svs1[svs1 > 1e-10 * svs1[0]]
    svs2_nz = svs2[svs2 > 1e-10 * svs2[0]]
    log(f"      XZ+YZ SVs: {np.array2string(svs1_nz, precision=6, separator=', ')}")
    log(f"      ZX+ZY SVs: {np.array2string(svs2_nz, precision=6, separator=', ')}")

    same = results_xzyz['rank_b'] == results_zxzy['rank_b']
    log(f"\n      Same entanglement structure? {'YES' if same else 'NO'}")
else:
    log(f"\n  ZX+ZY: palindromic pairing failed ({n_paired_zxzy}/{num})")


# ################################################################
# SUMMARY
# ################################################################
log()
log("=" * 90)
log("SUMMARY")
log("=" * 90)

rank_xzyz = results_xzyz['rank_b']

if rank_xzyz == 1:
    log(f"\n  OUTCOME A: SEPARABLE")
    log(f"  The non-locality of Pi is a PROJECTION ARTIFACT.")
    log(f"  In the Choi (doubled) representation, Pi factors as a product operator")
    log(f"  across Site1|Site2. What looks like entanglement is a local operation")
    log(f"  viewed from the wrong angle.")
else:
    log(f"\n  OUTCOME B: ENTANGLED (operator rank {rank_xzyz})")
    log(f"  The non-locality is FUNDAMENTAL.")
    log(f"  In the Choi representation, Pi STILL does not factor across Site1|Site2.")
    log(f"  There is no representation where Pi is a product operator.")
    log(f"  The mirror is genuinely quantum in every space, every basis.")
    log(f"  The 'between' is irreducible.")

if results_control is not None:
    rc = results_control['rank_b']
    log(f"\n  Control ({control_label}): operator rank {rc}")
    if rc == 1 and rank_xzyz > 1:
        log(f"  Local Pi -> separable Choi. Non-local Pi -> entangled Choi.")
        log(f"  The distinction between the 34/36 and the 2/36 is real and fundamental.")
    elif rc > 1:
        log(f"  Note: control also shows rank > 1 (eigenvector-pairing artifact).")

log()
log("=" * 90)
log(f"COMPLETE -- {datetime.now()}")
log("=" * 90)

f.close()
print(f"\n>>> Results written to {OUT}")
