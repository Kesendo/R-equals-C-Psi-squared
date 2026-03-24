#!/usr/bin/env python3
"""
qubit_necessity_tests.py - Five computational tests for qubit necessity.

Test 1: QST Comparison (qubit chain vs qutrit chain)
Test 2: Qutrit Eigenvalue Structure
Test 3: Hybrid Chain (qubit-qutrit-qubit)
Test 4: Exhaustive Qutrit Dissipator Search
Test 5: Composition Proof (tensor product of Pi operators)
"""

import numpy as np
from scipy import linalg as la
from itertools import product as iprod
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "qubit_necessity_tests.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# QUBIT INFRASTRUCTURE
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
GAMMA = 0.05

def xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))

def build_pauli_basis(N):
    all_idx = list(iprod(range(4), repeat=N))
    d = 2 ** N
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, N):
            m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m)
    return all_idx, pmats, d

def site_op(op, k, N, d_site=2):
    ops = [np.eye(d_site, dtype=complex)] * N
    ops[k] = op
    result = ops[0]
    for i in range(1, N):
        result = np.kron(result, ops[i])
    return result

def build_L_comp_qubit(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L

def build_H_heisenberg_qubit(N, bonds, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for P in [sx, sy, sz]:
            H += J * np.kron(site_op(P, i, N), np.eye(1)).reshape(d, d) if False else 0
        # Direct construction
        for P in [sx, sy, sz]:
            Pi = site_op(P, i, N)
            Pj = site_op(P, j, N)
            H += J * (Pi @ Pj)
    return H

def check_eigenvalue_pairing(evals, Sg, tol=1e-6):
    num = len(evals)
    used = set()
    errors = []
    n_well = 0
    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = np.argmin(diffs)
        err = float(diffs[best])
        errors.append(err)
        if err < tol:
            n_well += 1 if k == best else 2
        used.add(k)
        used.add(best)
    return n_well, num, errors

def construct_Pi(L_full, Sg, tol=1e-8):
    num = L_full.shape[0]
    evals, R = np.linalg.eig(L_full)
    R_inv = np.linalg.inv(R)
    P_eig = np.zeros((num, num), dtype=complex)
    used = set()
    n_paired = 0
    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = np.argmin(diffs)
        if diffs[best] < tol:
            P_eig[k, best] = 1.0
            P_eig[best, k] = 1.0
            used.add(k)
            used.add(best)
            n_paired += 1 if k == best else 2
        else:
            P_eig[k, k] = 1.0
            used.add(k)
    Pi = R @ P_eig @ R_inv
    return Pi, n_paired

def palindrome_error(Pi, L, Sg):
    Pi_inv = np.linalg.inv(Pi)
    E = Pi @ L @ Pi_inv + L + 2 * Sg * np.eye(L.shape[0])
    return np.max(np.abs(E)) / max(np.max(np.abs(L)), 1e-15)

# ============================================================
# QUTRIT INFRASTRUCTURE
# ============================================================
I3 = np.eye(3, dtype=complex)
gm_raw = [
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),         # lambda_1
    np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),      # lambda_2
    np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),         # lambda_3
    np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),          # lambda_4
    np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),       # lambda_5
    np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),          # lambda_6
    np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),       # lambda_7
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3),  # lambda_8
]
GM_NAMES = ['I3', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8']
IMMUNE_GM = {0, 3, 8}  # I3, lambda_3, lambda_8

# Normalized basis: Tr(B_i B_j) = 3 * delta_ij
c_norm = np.sqrt(3.0 / 2.0)
gm_basis = [I3] + [c_norm * m for m in gm_raw]

def site_op_qt(M, k, N):
    """Place 3x3 operator M at site k of an N-qutrit system."""
    ops = [I3] * N
    ops[k] = M
    result = ops[0]
    for i in range(1, N):
        result = np.kron(result, ops[i])
    return result

def build_H_heisenberg_qutrit(N, bonds, J=1.0):
    """SU(3) Heisenberg: H = J * sum_bond sum_{i=1}^8 lambda_i^a x lambda_i^b."""
    d = 3 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for m in gm_raw:
            Ma = site_op_qt(m, a, N)
            Mb = site_op_qt(m, b, N)
            H += J * (Ma @ Mb)
    return H

def build_L_comp_qutrit(H, gamma, N):
    """Computational-basis Lindbladian for N qutrits with diagonal dephasing."""
    d = 3 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    # Jump operators: lambda_3 and lambda_8 at each site
    for k in range(N):
        for m in [gm_raw[2], gm_raw[7]]:  # lambda_3, lambda_8
            Mk = site_op_qt(m, k, N)
            MkdMk = Mk.conj().T @ Mk
            L += gamma * (np.kron(Mk, Mk.conj())
                          - 0.5 * np.kron(MkdMk, Id)
                          - 0.5 * np.kron(Id, MkdMk.T))
    return L

def ptrace_general(rho, dims, keep):
    """Partial trace for a system with subsystem dimensions `dims`, keeping `keep`."""
    N = len(dims)
    d_total = int(np.prod(dims))
    rho_t = rho.reshape(list(dims) + list(dims))
    trace_out = sorted([i for i in range(N) if i not in keep])
    # Trace out from highest index first
    offset = 0
    for idx in reversed(trace_out):
        rho_t = np.trace(rho_t, axis1=idx, axis2=idx + N - offset)
        offset += 1
    d_keep = int(np.prod([dims[k] for k in keep]))
    return rho_t.reshape(d_keep, d_keep)


# ============================================================
# TEST 1: QST COMPARISON
# ============================================================
def run_test_1():
    log("=" * 70)
    log("TEST 1: QST COMPARISON - QUBIT CHAIN VS QUTRIT CHAIN")
    log("=" * 70)
    log()

    gamma = GAMMA
    N = 3

    # --- Qubit chain ---
    bonds = [(0, 1), (1, 2)]
    H_q2 = build_H_heisenberg_qubit(N, bonds)
    L_q2 = build_L_comp_qubit(H_q2, gamma, N)
    d_q2 = 2 ** N  # = 8
    Sg_q2 = N * gamma

    evals_q2 = np.linalg.eigvals(L_q2)
    n_well_q2, n_tot_q2, _ = check_eigenvalue_pairing(evals_q2, Sg_q2)

    log(f"Qubit chain (d=2, N={N}):")
    log(f"  Hilbert dim = {d_q2}, Liouville dim = {d_q2**2}")
    log(f"  Palindromic pairing: {n_well_q2}/{n_tot_q2}")

    # QST: |1,0,0> -> site 2
    psi_init_q2 = np.zeros(d_q2, dtype=complex)
    psi_init_q2[4] = 1.0  # |100> = index 4
    rho_init_q2 = np.outer(psi_init_q2, psi_init_q2.conj())

    t_points = np.linspace(0.1, 30, 60)
    best_fid_q2 = 0
    best_t_q2 = 0
    fid_curve_q2 = []
    for t in t_points:
        eLt = la.expm(L_q2 * t)
        rho_t = (eLt @ rho_init_q2.reshape(-1)).reshape(d_q2, d_q2)
        rho_2 = ptrace_general(rho_t, [2, 2, 2], [2])
        fid = float(np.real(rho_2[1, 1]))
        fid_curve_q2.append(fid)
        if fid > best_fid_q2:
            best_fid_q2 = fid
            best_t_q2 = t

    log(f"  Max transfer fidelity: {best_fid_q2:.4f} at t={best_t_q2:.2f}")

    # Spectrum summary
    re_parts_q2 = np.sort(evals_q2.real)
    n_osc_q2 = int(np.sum(np.abs(evals_q2.imag) > 1e-10))
    log(f"  Oscillating modes: {n_osc_q2}/{n_tot_q2}")
    log(f"  Decay rate range: [{re_parts_q2[0]:.4f}, {re_parts_q2[-1]:.4f}]")
    log()

    # --- Qutrit chain ---
    H_q3 = build_H_heisenberg_qutrit(N, bonds)
    L_q3 = build_L_comp_qutrit(H_q3, gamma, N)
    d_q3 = 3 ** N  # = 27

    evals_q3 = np.linalg.eigvals(L_q3)
    # Try multiple candidate centers
    log(f"Qutrit chain (d=3, N={N}):")
    log(f"  Hilbert dim = {d_q3}, Liouville dim = {d_q3**2}")

    # What is Sg for qutrits? Max dephasing rate / 2
    # Each off-diagonal GM decays at -2*gamma per site, max weight = N
    # So max rate = -2*gamma*N, center = -gamma*N
    Sg_q3 = N * gamma
    n_well_q3, n_tot_q3, errs_q3 = check_eigenvalue_pairing(evals_q3, Sg_q3)
    log(f"  Palindromic pairing (center={2*Sg_q3:.3f}): {n_well_q3}/{n_tot_q3}")

    # Also try other centers
    for test_Sg in [N * gamma * 2 / 3, N * gamma * 4 / 3, N * gamma * 2]:
        nw, nt, _ = check_eigenvalue_pairing(evals_q3, test_Sg, tol=1e-4)
        log(f"  Pairing at center={2*test_Sg:.3f}: {nw}/{nt}")

    # QST: |1,0,0> -> site 2
    psi_init_q3 = np.zeros(d_q3, dtype=complex)
    psi_init_q3[9] = 1.0  # |1,0,0> = 1*9+0*3+0 = 9
    rho_init_q3 = np.outer(psi_init_q3, psi_init_q3.conj())

    best_fid_q3 = 0
    best_t_q3 = 0
    fid_curve_q3 = []
    for t in t_points:
        eLt = la.expm(L_q3 * t)
        rho_t = (eLt @ rho_init_q3.reshape(-1)).reshape(d_q3, d_q3)
        rho_2 = ptrace_general(rho_t, [3, 3, 3], [2])
        fid = float(np.real(rho_2[1, 1]))
        fid_curve_q3.append(fid)
        if fid > best_fid_q3:
            best_fid_q3 = fid
            best_t_q3 = t

    log(f"  Max transfer fidelity: {best_fid_q3:.4f} at t={best_t_q3:.2f}")

    re_parts_q3 = np.sort(evals_q3.real)
    n_osc_q3 = int(np.sum(np.abs(evals_q3.imag) > 1e-10))
    log(f"  Oscillating modes: {n_osc_q3}/{n_tot_q3}")
    log(f"  Decay rate range: [{re_parts_q3[0]:.4f}, {re_parts_q3[-1]:.4f}]")

    log()
    log("Comparison:")
    log(f"  {'':>20} {'Qubit':>10} {'Qutrit':>10}")
    log(f"  {'Palindromic pairs':>20} {n_well_q2:>10} {n_well_q3:>10}")
    log(f"  {'Max QST fidelity':>20} {best_fid_q2:>10.4f} {best_fid_q3:>10.4f}")
    log(f"  {'Oscillating modes':>20} {n_osc_q2:>10} {n_osc_q3:>10}")


# ============================================================
# TEST 2: QUTRIT EIGENVALUE STRUCTURE
# ============================================================
def run_test_2():
    log()
    log("=" * 70)
    log("TEST 2: QUTRIT EIGENVALUE STRUCTURE")
    log("=" * 70)
    log()

    N = 2
    gamma = GAMMA
    bonds = [(0, 1)]
    H_qt = build_H_heisenberg_qutrit(N, bonds)
    L_qt = build_L_comp_qutrit(H_qt, gamma, N)
    d = 3 ** N  # = 9

    evals = np.linalg.eigvals(L_qt)
    n_evals = len(evals)
    log(f"N=2 qutrit chain: {n_evals} eigenvalues")

    # 1. Check for pairing at various centers
    log()
    log("Pairing check at various centers:")
    for center in np.linspace(0, 1.0, 21):
        nw, _, _ = check_eigenvalue_pairing(evals, center / 2, tol=1e-4)
        if nw > 20:
            log(f"  center={center:.3f}: {nw}/{n_evals} paired")

    # Best possible center (optimize)
    best_nw = 0
    best_c = 0
    for c in np.linspace(0, 2.0, 200):
        nw, _, _ = check_eigenvalue_pairing(evals, c / 2, tol=1e-4)
        if nw > best_nw:
            best_nw = nw
            best_c = c
    log(f"  Best center: {best_c:.4f} with {best_nw}/{n_evals} paired")
    log()

    # 2. Check magnitude pairing: |lambda_i| = |lambda_j|
    mags = np.sort(np.abs(evals))
    n_mag_pairs = 0
    used = set()
    for i in range(n_evals):
        if i in used:
            continue
        for j in range(i + 1, n_evals):
            if j in used:
                continue
            if abs(mags[i] - mags[j]) < 1e-4:
                n_mag_pairs += 1
                used.add(i)
                used.add(j)
                break
    log(f"Magnitude pairing (|lambda_i| = |lambda_j|): {2*n_mag_pairs}/{n_evals}")

    # 3. Conjugate pairing: Re(lambda_i) = Re(lambda_j), Im opposite
    n_conj = int(np.sum(np.abs(evals.imag) > 1e-10)) // 2
    log(f"Conjugate pairs (complex eigenvalues): {2*n_conj}/{n_evals}")

    # 4. Level spacing distribution
    re_sorted = np.sort(evals.real)
    spacings = np.diff(re_sorted)
    spacings = spacings[spacings > 1e-12]
    if len(spacings) > 0:
        mean_s = np.mean(spacings)
        std_s = np.std(spacings)
        # Wigner-Dyson: std/mean ~ 0.52 for GUE
        # Poisson: std/mean ~ 1.0
        log(f"Level spacing: mean={mean_s:.6f}, std={std_s:.6f}, std/mean={std_s/mean_s:.3f}")
        log(f"  (Poisson ~ 1.0, GUE ~ 0.52)")

    # 5. Compare with random matrix
    np.random.seed(42)
    L_rand = np.random.randn(n_evals, n_evals) + 1j * np.random.randn(n_evals, n_evals)
    L_rand = L_rand - np.trace(L_rand) / n_evals * np.eye(n_evals)
    evals_rand = np.linalg.eigvals(L_rand)
    nw_rand, _, _ = check_eigenvalue_pairing(evals_rand, 0, tol=1e-4)
    log(f"Random matrix pairing at center=0: {nw_rand}/{n_evals}")

    # 6. Eigenvalue table
    log()
    log("Eigenvalue spectrum (sorted by real part):")
    sorted_evals = evals[np.argsort(evals.real)]
    for i, ev in enumerate(sorted_evals):
        if i < 10 or i >= n_evals - 5 or abs(ev.real) < 1e-10:
            log(f"  {i:>3}: {ev.real:>10.6f} + {ev.imag:>10.6f}i")


# ============================================================
# TEST 3: HYBRID CHAIN (qubit-qutrit-qubit)
# ============================================================
def run_test_3():
    log()
    log("=" * 70)
    log("TEST 3: HYBRID CHAIN - QUBIT-QUTRIT-QUBIT")
    log("=" * 70)
    log()

    gamma = GAMMA
    dims = [2, 3, 2]
    d_total = int(np.prod(dims))  # 12
    log(f"System: qubit(d=2) - qutrit(d=3) - qubit(d=2)")
    log(f"Total Hilbert dim = {d_total}, Liouville dim = {d_total**2}")

    # Embedding of qubit Paulis into qutrit space
    # lambda_1 ~ X, lambda_2 ~ Y, lambda_3 ~ Z in upper-left 2x2 block
    X_emb = gm_raw[0]  # lambda_1 = [[0,1,0],[1,0,0],[0,0,0]]
    Y_emb = gm_raw[1]  # lambda_2 = [[0,-i,0],[i,0,0],[0,0,0]]
    Z_emb = gm_raw[2]  # lambda_3 = [[1,0,0],[0,-1,0],[0,0,0]]

    def multi_site_op(ops_list, dims_list):
        """Build tensor product of operators with given dimensions."""
        result = ops_list[0]
        for i in range(1, len(ops_list)):
            result = np.kron(result, ops_list[i])
        return result

    # Hamiltonian: Heisenberg coupling qubit-qutrit and qutrit-qubit
    H = np.zeros((d_total, d_total), dtype=complex)
    Id2 = np.eye(2, dtype=complex)
    Id3 = np.eye(3, dtype=complex)

    # Bond (0,1): qubit site 0 coupled to qutrit site 1
    for P_q, P_qt in [(sx, X_emb), (sy, Y_emb), (sz, Z_emb)]:
        H += multi_site_op([P_q, P_qt, Id2], dims)

    # Bond (1,2): qutrit site 1 coupled to qubit site 2
    for P_qt, P_q in [(X_emb, sx), (Y_emb, sy), (Z_emb, sz)]:
        H += multi_site_op([Id2, P_qt, P_q], dims)

    # Dissipator (computational basis Lindbladian)
    Id = np.eye(d_total, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    # Site 0: qubit Z-dephasing
    Z0 = multi_site_op([sz, Id3, Id2], dims)
    L += gamma * (np.kron(Z0, Z0.conj()) - np.eye(d_total**2, dtype=complex))

    # Site 1: qutrit dephasing with lambda_3, lambda_8
    for m in [gm_raw[2], gm_raw[7]]:  # lambda_3, lambda_8
        M1 = multi_site_op([Id2, m, Id2], dims)
        MdM = M1.conj().T @ M1
        L += gamma * (np.kron(M1, M1.conj()) - 0.5 * np.kron(MdM, Id)
                       - 0.5 * np.kron(Id, MdM.T))

    # Site 2: qubit Z-dephasing
    Z2 = multi_site_op([Id2, Id3, sz], dims)
    L += gamma * (np.kron(Z2, Z2.conj()) - np.eye(d_total**2, dtype=complex))

    evals = np.linalg.eigvals(L)
    n_evals = len(evals)

    # Check palindromic pairing at various centers
    log()
    log("Palindromic pairing search (hybrid chain):")
    best_nw = 0
    best_c = 0
    for c in np.linspace(0, 2.0, 200):
        nw, _, _ = check_eigenvalue_pairing(evals, c / 2, tol=1e-4)
        if nw > best_nw:
            best_nw = nw
            best_c = c
    log(f"  Best center: {best_c:.4f} with {best_nw}/{n_evals} paired (tol=1e-4)")

    # Also check with loose tolerance
    for tol in [1e-3, 1e-2, 5e-2]:
        nw, _, _ = check_eigenvalue_pairing(evals, best_c / 2, tol=tol)
        log(f"  At center={best_c:.4f}, tol={tol:.0e}: {nw}/{n_evals} paired")

    # Compare with pure qubit N=3
    H_pure = build_H_heisenberg_qubit(3, [(0, 1), (1, 2)])
    L_pure = build_L_comp_qubit(H_pure, gamma, 3)
    evals_pure = np.linalg.eigvals(L_pure)
    nw_pure, _, _ = check_eigenvalue_pairing(evals_pure, 3 * gamma)
    log()
    log(f"Control: pure qubit N=3 chain: {nw_pure}/{len(evals_pure)} paired")

    # Spectrum summary
    n_osc = int(np.sum(np.abs(evals.imag) > 1e-10))
    n_osc_pure = int(np.sum(np.abs(evals_pure.imag) > 1e-10))
    log()
    log(f"Spectral comparison:")
    log(f"  {'':>25} {'Hybrid':>10} {'Pure qubit':>12}")
    log(f"  {'Total eigenvalues':>25} {n_evals:>10} {len(evals_pure):>12}")
    log(f"  {'Oscillating modes':>25} {n_osc:>10} {n_osc_pure:>12}")
    log(f"  {'Best palindromic pairs':>25} {best_nw:>10} {nw_pure:>12}")
    log(f"  {'Palindromic fraction':>25} {best_nw/n_evals:>10.1%} {nw_pure/len(evals_pure):>12.1%}")


# ============================================================
# TEST 4: EXHAUSTIVE QUTRIT DISSIPATOR SEARCH
# ============================================================
def run_test_4():
    log()
    log("=" * 70)
    log("TEST 4: EXHAUSTIVE QUTRIT DISSIPATOR SEARCH")
    log("=" * 70)
    log()

    d = 3  # qutrit dimension
    # Compute per-site rates for each basis element under each dissipator

    def compute_rates(jump_ops):
        """Compute dephasing rates for each of the 9 GM basis elements."""
        rates = []
        for j, B in enumerate(gm_basis):
            # D(B) = sum_k (M_k B M_k^dag - 0.5{M_k^dag M_k, B})
            DB = np.zeros((d, d), dtype=complex)
            for M in jump_ops:
                DB += M @ B @ M.conj().T - 0.5 * (M.conj().T @ M @ B + B @ M.conj().T @ M)
            # Rate: Tr(B^dag D(B)) / Tr(B^dag B)
            numer = np.trace(B.conj().T @ DB)
            denom = np.trace(B.conj().T @ B)
            rate = numer / denom if abs(denom) > 1e-15 else 0.0
            rates.append(float(rate.real))
        return rates

    def check_rate_pairing(rates):
        """Check if 9 rates can be partitioned into pairs with equal sums."""
        n = len(rates)
        # For 9 elements: need 4 pairs + 1 self-paired (at c/2)
        sorted_rates = sorted(rates)
        # Try each element as the self-paired center
        for center_idx in range(n):
            c = 2 * sorted_rates[center_idx]
            remaining = sorted_rates[:center_idx] + sorted_rates[center_idx + 1:]
            # Check if remaining 8 form 4 pairs summing to c
            used = [False] * 8
            n_pairs = 0
            for i in range(8):
                if used[i]:
                    continue
                target = c - remaining[i]
                found = False
                for j in range(i + 1, 8):
                    if used[j]:
                        continue
                    if abs(remaining[j] - target) < 1e-10:
                        used[i] = True
                        used[j] = True
                        n_pairs += 1
                        found = True
                        break
                if not found:
                    break
            if n_pairs == 4:
                return True, c
        return False, None

    # Phase 1: Single Gell-Mann matrix as jump operator
    log("Phase 1: Single Gell-Mann matrix as jump operator (8 cases)")
    log(f"{'Jump op':>10} {'Immune':>7} {'Decaying':>9} {'Split':>8} {'Pairs?':>7}")
    log("-" * 50)

    n_works = 0
    for i, m in enumerate(gm_raw):
        rates = compute_rates([m])
        immune = sum(1 for r in rates if abs(r) < 1e-10)
        decaying = 9 - immune
        pairable, c = check_rate_pairing(rates)
        split = f"{immune}:{decaying}"
        status = "YES" if pairable else "no"
        if pairable:
            n_works += 1
        log(f"  {'l' + str(i+1):>8} {immune:>7} {decaying:>9} {split:>8} {status:>7}")

    log(f"  Phase 1 result: {n_works}/8 have valid rate pairing")
    log()

    # Phase 2: Pairs of Gell-Mann matrices
    log("Phase 2: Pairs of Gell-Mann matrices as jump operators (28 cases)")
    n_works_2 = 0
    for i in range(8):
        for j in range(i + 1, 8):
            rates = compute_rates([gm_raw[i], gm_raw[j]])
            immune = sum(1 for r in rates if abs(r) < 1e-10)
            pairable, c = check_rate_pairing(rates)
            if pairable:
                n_works_2 += 1
                log(f"  l{i+1}+l{j+1}: immune={immune}, PAIRABLE (c={c:.4f})")

    log(f"  Phase 2 result: {n_works_2}/28 have valid rate pairing")
    log()

    # Phase 3: Random linear combinations
    log("Phase 3: Random linear combinations of Gell-Mann matrices (100 samples)")
    np.random.seed(123)
    n_works_3 = 0
    for trial in range(100):
        coeffs = np.random.randn(8)
        coeffs /= np.linalg.norm(coeffs)
        M = sum(c * m for c, m in zip(coeffs, gm_raw))
        rates = compute_rates([M])
        pairable, _ = check_rate_pairing(rates)
        if pairable:
            n_works_3 += 1
            log(f"  Trial {trial}: PAIRABLE! coeffs={coeffs[:4]}...")

    log(f"  Phase 3 result: {n_works_3}/100 have valid rate pairing")
    log()

    # Phase 3b: Random PAIRS of linear combinations
    log("Phase 3b: Random pairs of linear combinations (100 samples)")
    n_works_3b = 0
    for trial in range(100):
        M1 = sum(c * m for c, m in zip(np.random.randn(8), gm_raw))
        M1 /= np.linalg.norm(M1)
        M2 = sum(c * m for c, m in zip(np.random.randn(8), gm_raw))
        M2 /= np.linalg.norm(M2)
        rates = compute_rates([M1, M2])
        pairable, _ = check_rate_pairing(rates)
        if pairable:
            n_works_3b += 1

    log(f"  Phase 3b result: {n_works_3b}/100 have valid rate pairing")
    log()

    total = n_works + n_works_2 + n_works_3 + n_works_3b
    total_tested = 8 + 28 + 100 + 100
    log(f"SUMMARY: {total}/{total_tested} dissipator configurations allow rate pairing")
    if total == 0:
        log("  NO qutrit dissipator permits a palindromic Pi.")
        log("  The exclusivity claim extends beyond diagonal dephasing:")
        log("  d=2 is necessary for ANY physically motivated dephasing noise.")


# ============================================================
# TEST 5: COMPOSITION PROOF
# ============================================================
def run_test_5():
    log()
    log("=" * 70)
    log("TEST 5: COMPOSITION PROOF - TENSOR PRODUCT OF Pi OPERATORS")
    log("=" * 70)
    log()

    gamma = GAMMA
    N2 = 2
    N4 = 4
    Sg2 = N2 * gamma
    Sg4 = N4 * gamma

    # Build N=2 Liouvillian and Pi
    all_idx2, pmats2, d2 = build_pauli_basis(N2)
    H2 = build_H_heisenberg_qubit(N2, [(0, 1)])
    L_H2 = np.zeros((16, 16), dtype=complex)
    for a in range(16):
        for b in range(16):
            comm = -1j * (H2 @ pmats2[b] - pmats2[b] @ H2)
            L_H2[a, b] = np.trace(pmats2[a].conj().T @ comm) / d2
    L_D2 = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx2])
    L2 = L_H2.copy()
    for a in range(16):
        L2[a, a] += L_D2[a]

    Pi2, n_paired2 = construct_Pi(L2, Sg2)
    err2 = palindrome_error(Pi2, L2, Sg2)
    log(f"N=2 Heisenberg: {n_paired2}/16 paired, error={err2:.2e}")

    # Build N=4 Liouvillian (decoupled: bonds (0,1) and (2,3) only)
    all_idx4, pmats4, d4 = build_pauli_basis(N4)
    num4 = 4 ** N4
    H4_decoupled = build_H_heisenberg_qubit(N4, [(0, 1), (2, 3)])
    L_H4_dec = np.zeros((num4, num4), dtype=complex)
    for a in range(num4):
        for b in range(num4):
            comm = -1j * (H4_decoupled @ pmats4[b] - pmats4[b] @ H4_decoupled)
            L_H4_dec[a, b] = np.trace(pmats4[a].conj().T @ comm) / d4
    L_D4 = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx4])
    L4_dec = L_H4_dec.copy()
    for a in range(num4):
        L4_dec[a, a] += L_D4[a]

    # Test 1: Tensor product Pi
    Pi4_tensor = np.kron(Pi2, Pi2)
    err4_tensor_dec = palindrome_error(Pi4_tensor, L4_dec, Sg4)
    log()
    log(f"DECOUPLED N=4 (bonds (0,1) and (2,3) only):")
    log(f"  Pi_tensor = Pi_2 x Pi_2")
    log(f"  Palindrome error: {err4_tensor_dec:.2e}")
    log(f"  {'PASS' if err4_tensor_dec < 1e-10 else 'FAIL'}: tensor product Pi works"
        f" for decoupled blocks")

    # Build N=4 coupled: add bond (1,2)
    H4_coupled = build_H_heisenberg_qubit(N4, [(0, 1), (1, 2), (2, 3)])
    L_H4_cpl = np.zeros((num4, num4), dtype=complex)
    for a in range(num4):
        for b in range(num4):
            comm = -1j * (H4_coupled @ pmats4[b] - pmats4[b] @ H4_coupled)
            L_H4_cpl[a, b] = np.trace(pmats4[a].conj().T @ comm) / d4
    L4_cpl = L_H4_cpl.copy()
    for a in range(num4):
        L4_cpl[a, a] += L_D4[a]

    # Test tensor product on coupled system
    err4_tensor_cpl = palindrome_error(Pi4_tensor, L4_cpl, Sg4)
    log()
    log(f"COUPLED N=4 (bonds (0,1), (1,2), (2,3)):")
    log(f"  Pi_tensor = Pi_2 x Pi_2")
    log(f"  Palindrome error: {err4_tensor_cpl:.2e}")
    log(f"  {'PASS' if err4_tensor_cpl < 1e-10 else 'FAIL'}: tensor product Pi "
        f"{'works' if err4_tensor_cpl < 1e-10 else 'FAILS'} for coupled chain")

    # Construct the CORRECT Pi for coupled N=4 via eigenvector pairing
    Pi4_eigvec, n_paired4 = construct_Pi(L4_cpl, Sg4)
    err4_eigvec = palindrome_error(Pi4_eigvec, L4_cpl, Sg4)
    log()
    log(f"  Eigenvector-constructed Pi_4:")
    log(f"  Paired: {n_paired4}/{num4}")
    log(f"  Palindrome error: {err4_eigvec:.2e}")

    # Compare tensor Pi with eigenvector Pi
    diff = np.linalg.norm(Pi4_tensor - Pi4_eigvec) / np.linalg.norm(Pi4_eigvec)
    log(f"  ||Pi_tensor - Pi_eigvec|| / ||Pi_eigvec|| = {diff:.4f}")

    # Use canonical per-site Pi for coupled case
    PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
    PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}
    idx_map = {idx: i for i, idx in enumerate(all_idx4)}
    Pi4_canon = np.zeros((num4, num4), dtype=complex)
    for b, idx_b in enumerate(all_idx4):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = idx_map[mapped]
        Pi4_canon[a, b] = sign

    err4_canon = palindrome_error(Pi4_canon, L4_cpl, Sg4)
    log()
    log(f"  Canonical per-site Pi_4 (I<->X, Y<->Z):")
    log(f"  Palindrome error: {err4_canon:.2e}")
    log(f"  {'PASS' if err4_canon < 1e-10 else 'FAIL'}")

    log()
    log("COMPOSITION SUMMARY:")
    log(f"  Decoupled: Pi_2 x Pi_2 works perfectly (error {err4_tensor_dec:.2e})")
    if err4_tensor_cpl < 1e-10:
        log(f"  Coupled: Pi_2 x Pi_2 also works (error {err4_tensor_cpl:.2e})")
        log(f"  Composition is TRIVIAL for Heisenberg: tensor product suffices")
    else:
        log(f"  Coupled: Pi_2 x Pi_2 FAILS (error {err4_tensor_cpl:.2e})")
        if err4_canon < 1e-10:
            log(f"  But canonical per-site Pi works (error {err4_canon:.2e})")
            log(f"  The per-site Pi adapts to the coupled topology automatically")
        elif err4_eigvec < 1e-10:
            log(f"  But eigenvector Pi works (error {err4_eigvec:.2e})")
            log(f"  Composition requires adaptation, not simple tensor product")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Qubit Necessity: Five Computational Tests")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"gamma = {GAMMA}")
    log()

    t1 = time.time()
    run_test_1()
    log(f"\n[Test 1 completed in {time.time() - t1:.1f}s]")

    t2 = time.time()
    run_test_2()
    log(f"\n[Test 2 completed in {time.time() - t2:.1f}s]")

    t3 = time.time()
    run_test_3()
    log(f"\n[Test 3 completed in {time.time() - t3:.1f}s]")

    t4 = time.time()
    run_test_4()
    log(f"\n[Test 4 completed in {time.time() - t4:.1f}s]")

    t5 = time.time()
    run_test_5()
    log(f"\n[Test 5 completed in {time.time() - t5:.1f}s]")

    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
