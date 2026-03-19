"""
Continuous Pi Search -- Find Pi for the Last 2/36 Cases (XZ+YZ, ZX+ZY)
=======================================================================
The discrete Pauli permutation search found 34/36 two-term Hamiltonians
palindromic under Z-dephasing. XZ+YZ and ZX+ZY are palindromic numerically
but have no discrete Pi. This script searches for CONTINUOUS per-site
operators (rotations, general unitaries, general invertible maps).

Script: simulations/continuous_pi_search.py
Output: simulations/results/continuous_pi_search.txt
"""
import numpy as np
from itertools import product as iproduct
from datetime import datetime
from scipy.optimize import minimize
from scipy.linalg import expm

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\continuous_pi_search.txt"
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# PAULI BASICS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}

# H term labels: all 9 Pauli-pair terms
H_LABELS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


def plabel(indices):
    return ''.join(NAMES[i] for i in indices)


# ============================================================
# BUILD SYSTEM IN PAULI BASIS
# ============================================================
def build_pauli_basis(N):
    """Build all Pauli string matrices and index tuples for N sites."""
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
    """Build Hamiltonian for a single Pauli-pair term on a bond."""
    d = 2 ** N
    i, j = bond
    P1 = PM[label[0]]
    P2 = PM[label[1]]
    ops = [I2] * N
    ops[i] = P1
    ops[j] = P2
    H = ops[0]
    for o in ops[1:]:
        H = np.kron(H, o)
    return H


def build_L_H(N, H, all_idx, pmats, d):
    """Build Hamiltonian superoperator in Pauli basis."""
    num = 4 ** N
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / d
    return L_H


def build_L_D_diag(N, gamma, all_idx):
    """Build Z-dephasing diagonal in Pauli basis."""
    return np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])


def build_Pi_from_M(M, N):
    """Build N-site Pi = M (x) M (x) ... (x) M using Kronecker products."""
    Pi = M.copy()
    for _ in range(N - 1):
        Pi = np.kron(Pi, M)
    return Pi


def build_Pi_nonuniform(Ms, N):
    """Build N-site Pi = M1 (x) M2 (x) ... (x) MN."""
    Pi = Ms[0].copy()
    for i in range(1, N):
        Pi = np.kron(Pi, Ms[i])
    return Pi


def anticomm_error(Pi, L_H):
    """Compute ||Pi L_H Pi^{-1} + L_H|| / ||L_H||."""
    try:
        Pi_inv = np.linalg.inv(Pi)
    except np.linalg.LinAlgError:
        return 1e10
    conj = Pi @ L_H @ Pi_inv
    err = np.max(np.abs(conj + L_H))
    norm = np.max(np.abs(L_H))
    return err / norm if norm > 1e-15 else err


def full_palindrome_error(Pi, L, Sg):
    """Compute ||Pi L Pi^{-1} + L + 2Sg I|| / ||L||."""
    try:
        Pi_inv = np.linalg.inv(Pi)
    except np.linalg.LinAlgError:
        return 1e10
    num = L.shape[0]
    conj = Pi @ L @ Pi_inv
    err = np.max(np.abs(conj + L + 2 * Sg * np.eye(num)))
    norm = np.max(np.abs(L))
    return err / norm if norm > 1e-15 else err


# ============================================================
# PER-SITE M BUILDERS (block off-diagonal structure)
# ============================================================
def build_M_block(A, B):
    """Build 4x4 per-site M from 2x2 blocks A ({X,Y}->{I,Z}) and B ({I,Z}->{X,Y}).

    In standard basis (I=0, X=1, Y=2, Z=3):
    M = | 0     A[0,0] A[0,1] 0      |
        | B[0,0] 0      0      B[0,1] |
        | B[1,0] 0      0      B[1,1] |
        | 0     A[1,0] A[1,1] 0      |
    """
    M = np.zeros((4, 4), dtype=complex)
    # A maps {X,Y} -> {I,Z}: rows 0,3 cols 1,2
    M[0, 1] = A[0, 0]
    M[0, 2] = A[0, 1]
    M[3, 1] = A[1, 0]
    M[3, 2] = A[1, 1]
    # B maps {I,Z} -> {X,Y}: rows 1,2 cols 0,3
    M[1, 0] = B[0, 0]
    M[1, 3] = B[0, 1]
    M[2, 0] = B[1, 0]
    M[2, 3] = B[1, 1]
    return M


def M_P1():
    """P1 family: I->X, X->I, Y->iZ, Z->iY."""
    A = np.array([[1, 0], [0, 1j]], dtype=complex)
    B = np.array([[1, 0], [0, 1j]], dtype=complex)
    return build_M_block(A, B)


def M_P4():
    """P4 family: I->Y, X->Z, Y->I, Z->X."""
    A = np.array([[0, 1], [1, 0]], dtype=complex)
    B = np.array([[0, 1], [1, 0]], dtype=complex)
    return build_M_block(A, B)


def M_rotation(theta):
    """Single-parameter interpolation between P1 and P4.
    theta=0 -> P1, theta=pi/2 -> P4."""
    A_P1 = np.array([[1, 0], [0, 1j]], dtype=complex)
    A_P4 = np.array([[0, 1], [1, 0]], dtype=complex)
    A = np.cos(theta) * A_P1 + np.sin(theta) * A_P4
    B = A.copy()  # symmetric: A = B
    return build_M_block(A, B)


def M_rotation_AB(theta_A, theta_B):
    """Two-parameter: separate rotations for A and B blocks."""
    A_P1 = np.array([[1, 0], [0, 1j]], dtype=complex)
    A_P4 = np.array([[0, 1], [1, 0]], dtype=complex)
    A = np.cos(theta_A) * A_P1 + np.sin(theta_A) * A_P4
    B = np.cos(theta_B) * A_P1 + np.sin(theta_B) * A_P4
    return build_M_block(A, B)


def M_from_params(params):
    """Build M from 16 real parameters (8 complex entries of A and B)."""
    # params[0:4] = Re,Im of A[0,0], A[0,1]
    # params[4:8] = Re,Im of A[1,0], A[1,1]
    # params[8:12] = Re,Im of B[0,0], B[0,1]
    # params[12:16] = Re,Im of B[1,0], B[1,1]
    A = np.array([
        [params[0] + 1j * params[1], params[2] + 1j * params[3]],
        [params[4] + 1j * params[5], params[6] + 1j * params[7]],
    ])
    B = np.array([
        [params[8] + 1j * params[9], params[10] + 1j * params[11]],
        [params[12] + 1j * params[13], params[14] + 1j * params[15]],
    ])
    return build_M_block(A, B)


# ============================================================
# MAIN
# ============================================================
gamma = 0.05

log("=" * 90)
log("CONTINUOUS PI SEARCH: Find Pi for XZ+YZ and ZX+ZY")
log(f"Date: {datetime.now()}")
log(f"gamma = {gamma}")
log("=" * 90)


# ############################################################
# SECTION 1: Verify the 2/36 problem
# ############################################################
log()
log("=" * 90)
log("SECTION 1: Verify palindromic status and discrete Pi availability")
log("=" * 90)

N = 2
bonds = [(0, 1)]
Sg = N * gamma
all_idx, pmats, d = build_pauli_basis(N)
num = 4 ** N
L_D_diag = build_L_D_diag(N, gamma, all_idx)

# Build L_H for each of the 9 terms
L_H_terms = {}
for label in H_LABELS:
    H = build_H_term(N, bonds[0], label)
    L_H_terms[label] = build_L_H(N, H, all_idx, pmats, d)

# Build discrete Pi candidates: P1 and P4 with phases
discrete_Ms = []
for base_M, name in [(M_P1(), 'P1'), (M_P4(), 'P4')]:
    for phases in iproduct([1, -1, 1j, -1j], repeat=4):
        M = base_M.copy()
        for col in range(4):
            M[:, col] *= phases[col]
        discrete_Ms.append((M, name))

log(f"\n  {len(discrete_Ms)} discrete per-site maps (P1 + P4 with all phase combos)")

# Check each term and pair
target_combos = ['XZ+YZ', 'ZX+ZY']
all_combos = []
for i, t1 in enumerate(H_LABELS):
    for t2 in H_LABELS[i + 1:]:
        all_combos.append(f"{t1}+{t2}")

log(f"\n  {'Combo':>10}  {'Palindromic':>12}  {'#Discrete Pi':>13}  {'Status':>10}")
log(f"  {'-' * 50}")

for combo in target_combos:
    t1, t2 = combo.split('+')
    L_H_combo = L_H_terms[t1] + L_H_terms[t2]
    L_full = L_H_combo.copy()
    for a in range(num):
        L_full[a, a] += L_D_diag[a]

    # Check palindrome
    evals = np.linalg.eigvals(L_full)
    paired = 0
    for k in range(num):
        target = -(evals[k] + 2 * Sg)
        if np.min(np.abs(evals - target)) < 1e-8:
            paired += 1
    is_palindromic = paired == num

    # Count discrete Pi that work
    n_valid = 0
    for M_cand, _ in discrete_Ms:
        Pi_cand = build_Pi_from_M(M_cand, N)
        err = anticomm_error(Pi_cand, L_H_combo)
        if err < 1e-8:
            n_valid += 1

    status = "PROBLEM" if is_palindromic and n_valid == 0 else "OK"
    log(f"  {combo:>10}  {'YES' if is_palindromic else 'NO':>12}  {n_valid:>13}  {status:>10}")


# ############################################################
# SECTION 2: Single-parameter rotation sweep
# ############################################################
log()
log("=" * 90)
log("SECTION 2: Single-parameter rotation sweep")
log("  M(theta) interpolates between P1 (theta=0) and P4 (theta=pi/2)")
log("  A(theta) = B(theta) = cos(theta)*A_P1 + sin(theta)*A_P4")
log("=" * 90)

for combo in target_combos:
    t1, t2 = combo.split('+')
    L_H_combo = L_H_terms[t1] + L_H_terms[t2]

    log(f"\n  {combo}:")
    log(f"  {'theta':>10}  {'theta/pi':>10}  {'error':>12}  {'det(M)':>12}")
    log(f"  {'-' * 48}")

    best_theta, best_err = 0, np.inf
    thetas = np.linspace(0, np.pi, 181)  # 1-degree resolution
    for theta in thetas:
        M = M_rotation(theta)
        if abs(np.linalg.det(M)) < 1e-10:
            continue
        Pi = build_Pi_from_M(M, N)
        err = anticomm_error(Pi, L_H_combo)
        det_M = np.linalg.det(M)
        if theta in [0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2,
                     5 * np.pi / 8, 3 * np.pi / 4, 7 * np.pi / 8, np.pi]:
            log(f"  {theta:>10.4f}  {theta / np.pi:>10.4f}  "
                f"{err:>12.4e}  {abs(det_M):>12.4e}")
        if err < best_err:
            best_err, best_theta = err, theta

    log(f"\n  Best theta: {best_theta:.6f} ({best_theta / np.pi:.6f}*pi)")
    log(f"  Best error: {best_err:.6e}")
    log(f"  {'FOUND' if best_err < 1e-8 else 'NOT FOUND'}")


# ############################################################
# SECTION 3: Two-parameter rotation (separate A, B)
# ############################################################
log()
log("=" * 90)
log("SECTION 3: Two-parameter rotation (theta_A, theta_B)")
log("=" * 90)

for combo in target_combos:
    t1, t2 = combo.split('+')
    L_H_combo = L_H_terms[t1] + L_H_terms[t2]

    best_err = np.inf
    best_params = (0, 0)
    thetas = np.linspace(0, np.pi, 91)
    for tA in thetas:
        for tB in thetas:
            M = M_rotation_AB(tA, tB)
            if abs(np.linalg.det(M)) < 1e-10:
                continue
            Pi = build_Pi_from_M(M, N)
            err = anticomm_error(Pi, L_H_combo)
            if err < best_err:
                best_err = err
                best_params = (tA, tB)

    log(f"\n  {combo}:")
    log(f"    Best (theta_A, theta_B) = ({best_params[0]:.4f}, {best_params[1]:.4f})")
    log(f"    = ({best_params[0] / np.pi:.4f}*pi, {best_params[1] / np.pi:.4f}*pi)")
    log(f"    Best error: {best_err:.6e}")
    log(f"    {'FOUND' if best_err < 1e-8 else 'NOT FOUND'}")


# ############################################################
# SECTION 4: Full 16-parameter optimization
# ############################################################
log()
log("=" * 90)
log("SECTION 4: Full block optimization (16 real parameters)")
log("  A and B are arbitrary 2x2 complex matrices (not necessarily unitary)")
log("=" * 90)

for combo in target_combos:
    t1, t2 = combo.split('+')
    L_H_combo = L_H_terms[t1] + L_H_terms[t2]

    def objective(params):
        M = M_from_params(params)
        if abs(np.linalg.det(M)) < 1e-12:
            return 1e10
        Pi = build_Pi_from_M(M, N)
        Pi_inv = np.linalg.inv(Pi)
        conj = Pi @ L_H_combo @ Pi_inv
        return np.sum(np.abs(conj + L_H_combo) ** 2)

    # Try many random starting points
    best_result = None
    best_fun = np.inf
    np.random.seed(42)

    # Also try starting from known discrete solutions
    starts = [np.random.randn(16) * 0.5 for _ in range(200)]
    # Add P1 and P4 as starting points
    for M_start in [M_P1(), M_P4()]:
        A_s = np.array([[M_start[0, 1], M_start[0, 2]],
                        [M_start[3, 1], M_start[3, 2]]])
        B_s = np.array([[M_start[1, 0], M_start[1, 3]],
                        [M_start[2, 0], M_start[2, 3]]])
        p = np.zeros(16)
        p[0], p[1] = A_s[0, 0].real, A_s[0, 0].imag
        p[2], p[3] = A_s[0, 1].real, A_s[0, 1].imag
        p[4], p[5] = A_s[1, 0].real, A_s[1, 0].imag
        p[6], p[7] = A_s[1, 1].real, A_s[1, 1].imag
        p[8], p[9] = B_s[0, 0].real, B_s[0, 0].imag
        p[10], p[11] = B_s[0, 1].real, B_s[0, 1].imag
        p[12], p[13] = B_s[1, 0].real, B_s[1, 0].imag
        p[14], p[15] = B_s[1, 1].real, B_s[1, 1].imag
        starts.append(p)
        # Also add perturbations
        for _ in range(20):
            starts.append(p + np.random.randn(16) * 0.3)

    for s_idx, x0 in enumerate(starts):
        res = minimize(objective, x0, method='Nelder-Mead',
                       options={'maxiter': 5000, 'xatol': 1e-12, 'fatol': 1e-14})
        if res.fun < best_fun:
            best_fun = res.fun
            best_result = res

    log(f"\n  {combo}:")
    log(f"    Optimization over {len(starts)} starting points")
    log(f"    Best objective: {best_fun:.6e}")
    M_opt = M_from_params(best_result.x)
    Pi_opt = build_Pi_from_M(M_opt, N)
    err_H = anticomm_error(Pi_opt, L_H_combo)
    log(f"    Anti-commutation error: {err_H:.6e}")

    if err_H < 1e-6:
        log(f"    FOUND per-site continuous Pi!")
        # Show the M matrix
        log(f"    M =")
        for row in range(4):
            entries = []
            for col in range(4):
                v = M_opt[row, col]
                if abs(v) < 1e-10:
                    entries.append("   0   ")
                elif abs(v.imag) < 1e-10:
                    entries.append(f"{v.real:+7.4f}")
                elif abs(v.real) < 1e-10:
                    entries.append(f"{v.imag:+6.4f}i")
                else:
                    entries.append(f"{v:.3f}")
                    entries[-1] = entries[-1][:7]
            log(f"      [{', '.join(entries)}]")

        # Check full palindrome
        L_full = L_H_combo.copy()
        for a in range(num):
            L_full[a, a] += L_D_diag[a]
        err_full = full_palindrome_error(Pi_opt, L_full, Sg)
        log(f"    Full palindrome error: {err_full:.6e}")
    else:
        log(f"    NOT FOUND with per-site M (x) M")
        log(f"    Minimum achievable error: {err_H:.6e}")


# ############################################################
# SECTION 5: Direct eigenvector construction (N=2)
# ############################################################
log()
log("=" * 90)
log("SECTION 5: Direct eigenvector construction")
log("  Build Pi from palindromic eigenvalue pairing (no per-site assumption)")
log("=" * 90)

for combo in target_combos:
    t1, t2 = combo.split('+')
    L_H_combo = L_H_terms[t1] + L_H_terms[t2]
    L_full = L_H_combo.copy()
    for a in range(num):
        L_full[a, a] += L_D_diag[a]

    evals, R = np.linalg.eig(L_full)

    # Pair eigenvalues palindromically
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
    log(f"\n  {combo}:")
    log(f"    Palindromic pairing: {n_paired}/{num}")

    if n_paired < num:
        log(f"    Incomplete pairing. Cannot construct Pi.")
        continue

    # Build Pi from eigenvector pairing: Pi |r_k> = alpha_k |r_{k'}>
    # Pi = R D R^{-1} where D maps column k to column pair_map[k]
    R_inv = np.linalg.inv(R)

    # Build the permutation+scaling in eigenbasis
    # For each k, Pi maps r_k -> r_{k'}. We need Pi = R @ P @ R_inv
    # where P is the pairing matrix in eigenbasis
    P_eig = np.zeros((num, num), dtype=complex)
    seen = set()
    for k in range(num):
        j = pair_map[k]
        if k in seen:
            continue
        seen.add(k)
        seen.add(j)
        if k == j:
            # Self-paired: Pi maps this to itself (up to phase)
            P_eig[k, k] = 1.0
        else:
            P_eig[j, k] = 1.0
            P_eig[k, j] = 1.0

    Pi_constructed = R @ P_eig @ R_inv

    # Verify it works
    err_full = full_palindrome_error(Pi_constructed, L_full, Sg)
    log(f"    Full palindrome error of constructed Pi: {err_full:.6e}")

    if err_full > 1e-4:
        log(f"    Construction failed (error too large)")
        # Try with scaling factors
        log(f"    Trying with optimal scaling per pair...")

        def opt_scaling(params):
            n_pairs = len(params)
            P_e = np.zeros((num, num), dtype=complex)
            pi_idx = 0
            seen2 = set()
            for kk in range(num):
                jj = pair_map[kk]
                if kk in seen2:
                    continue
                seen2.add(kk)
                seen2.add(jj)
                if kk == jj:
                    P_e[kk, kk] = 1.0
                else:
                    alpha = params[pi_idx] if pi_idx < len(params) else 1.0
                    P_e[jj, kk] = alpha
                    P_e[kk, jj] = 1.0 / alpha if abs(alpha) > 1e-15 else 1.0
                    pi_idx += 1
            Pi_c = R @ P_e @ R_inv
            try:
                Pi_c_inv = np.linalg.inv(Pi_c)
                return np.sum(np.abs(Pi_c @ L_full @ Pi_c_inv + L_full + 2 * Sg * np.eye(num)) ** 2)
            except Exception:
                return 1e10

        n_free_pairs = sum(1 for k in range(num) if pair_map[k] > k)
        x0 = np.ones(n_free_pairs)
        res = minimize(opt_scaling, x0, method='Nelder-Mead',
                       options={'maxiter': 10000, 'xatol': 1e-14, 'fatol': 1e-16})
        log(f"    Optimized scaling error: {res.fun:.6e}")
    else:
        log(f"    Pi constructed successfully!")

    # Check if Pi is a tensor product M1 (x) M2
    log(f"\n    Checking tensor product structure:")
    # For N=2: Pi is 16x16. If Pi = M1 (x) M2 where M1, M2 are 4x4,
    # then Pi[4*a+b, 4*c+d] = M1[a,c] * M2[b,d]
    # Extract M1 and M2 by looking at specific entries
    # M1[a,c] = Pi[4*a+0, 4*c+0] / M2[0,0] (if M2[0,0] != 0)
    # This is fragile. Try a different approach: SVD-based factorization.

    # Reshape Pi into (4,4,4,4) and check if it factors
    Pi_4d = Pi_constructed.reshape(4, 4, 4, 4)

    # If Pi = M1 (x) M2, then Pi_4d[a,b,c,d] = M1[a,c] * M2[b,d]
    # This means Pi_4d[:,:,c,d] = M2[b,d] * M1[:,c] is a rank-1 matrix
    # Check rank of various "slices"
    is_product = True
    for c in range(4):
        for dd in range(4):
            slice_mat = Pi_4d[:, :, c, dd]
            if np.max(np.abs(slice_mat)) < 1e-10:
                continue
            s = np.linalg.svd(slice_mat, compute_uv=False)
            rank = np.sum(s > 1e-8 * s[0])
            if rank > 1:
                is_product = False
                break
        if not is_product:
            break

    log(f"    Is tensor product M1 (x) M2? {'YES' if is_product else 'NO'}")

    if is_product:
        # Extract M1 and M2
        # Find a nonzero slice to use as reference
        for c0 in range(4):
            for d0 in range(4):
                ref = Pi_4d[:, :, c0, d0]
                if np.max(np.abs(ref)) > 1e-8:
                    U, s, Vh = np.linalg.svd(ref)
                    # ref = s[0] * U[:,0] @ Vh[0,:]
                    # This means M1[:,c0] ~ U[:,0] and M2[:,d0] ~ Vh[0,:]
                    break
            else:
                continue
            break
        log(f"    (Extraction of M1, M2 from tensor product)")
        # M1[a,c] * M2[b,d] = Pi_4d[a,b,c,d]
        # From the reference slice (c0,d0): M1[a,c0] * M2[b,d0] = Pi_4d[a,b,c0,d0]
        # So M1[:,c0] (x) M2[:,d0] = vec(Pi_4d[:,:,c0,d0])
        # Use SVD to extract
        M1_cols = {}
        M2_cols = {}
        for c in range(4):
            for dd in range(4):
                sl = Pi_4d[:, :, c, dd]
                if np.max(np.abs(sl)) < 1e-10:
                    continue
                U, s, Vh = np.linalg.svd(sl)
                m1_col = U[:, 0] * np.sqrt(s[0])
                m2_col = Vh[0, :].conj() * np.sqrt(s[0])
                # Normalize relative to reference
                if c not in M1_cols:
                    M1_cols[c] = m1_col
                    M2_cols[dd] = m2_col
        if len(M1_cols) == 4 and len(M2_cols) == 4:
            M1_est = np.column_stack([M1_cols[c] for c in range(4)])
            M2_est = np.column_stack([M2_cols[d] for d in range(4)])
            log(f"    M1 =")
            for row in range(4):
                log(f"      [{', '.join(f'{M1_est[row,c]:+8.4f}' for c in range(4))}]")
            log(f"    M2 =")
            for row in range(4):
                log(f"      [{', '.join(f'{M2_est[row,c]:+8.4f}' for c in range(4))}]")
    else:
        log(f"    Pi is GENUINELY NON-LOCAL (not a tensor product of per-site maps)")
        log(f"    This is a major structural finding: XZ+YZ requires entangled Pi")

        # Show the 16x16 Pi structure
        log(f"\n    Pi matrix (nonzero entries):")
        for a in range(num):
            for b in range(num):
                if abs(Pi_constructed[a, b]) > 1e-6:
                    log(f"      Pi[{plabel(all_idx[a])},{plabel(all_idx[b])}] = "
                        f"{Pi_constructed[a, b]:.4f}")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("SEARCH COMPLETE")
log(f"Date: {datetime.now()}")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
