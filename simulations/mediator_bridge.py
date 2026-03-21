#!/usr/bin/env python3
"""
The Mediator Bridge: S Was Always the Answer
==============================================
Test 1: Palindrome survives (theorem says yes, verify)
Test 2: Cross-pair information flow (the central question)
Test 3: Three conditions scale up
Test 4: Standing wave across the bridge
Test 5: Decoder across the bridge
Test 6: Robustness against direct leakage

Topology:  1 -- 2 -- M -- 4 -- 5
          [Pair A]   ^   [Pair B]
                   bridge

Script: simulations/mediator_bridge.py
Output: simulations/results/mediator_bridge.txt
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "mediator_bridge.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# INFRASTRUCTURE (N=5)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
PNAMES = ['I', 'X', 'Y', 'Z']
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 5
d = 2 ** N   # 32
d2 = d * d   # 1024


def site_op(op, k, nq=N):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H_heisenberg(bonds, Js, nq=N):
    """Heisenberg Hamiltonian. bonds: list of (i,j), Js: list of coupling strengths."""
    dim = 2 ** nq
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j), J in zip(bonds, Js):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, nq) @ site_op(P, j, nq)
    return H


def build_L_H(H, nq=N):
    dim = 2 ** nq
    Id = np.eye(dim)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_L_D_zdeph(gammas, nq=N):
    dim = 2 ** nq
    dim2 = dim * dim
    L_D = np.zeros((dim2, dim2), dtype=complex)
    for k in range(nq):
        Zk = site_op(sz, k, nq)
        L_D += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(dim2))
    return L_D


def add_lindblad_jump(L, F, kappa, nq=N):
    dim = 2 ** nq
    Id = np.eye(dim)
    FdF = F.conj().T @ F
    L += kappa * (np.kron(F, F.conj())
                  - 0.5 * np.kron(FdF, Id)
                  - 0.5 * np.kron(Id, FdF.T))
    return L


def palindrome_err(evals, Sg, tol=1e-6):
    n = len(evals); max_err = 0; n_paired = 0
    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        best = np.min(np.abs(evals - target))
        if best < tol: n_paired += 1
        if best > max_err: max_err = best
    return max_err, n_paired, n


def best_palindrome(evals, Sg_base, tol=1e-10):
    """Two-stage Sg search."""
    best_np = 0; best_Sg = Sg_base; best_err = 1e30
    candidates = np.concatenate([[Sg_base], np.linspace(0, Sg_base * 5, 60)])
    for Sg_try in candidates:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=1e-4)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try
    fine = np.linspace(best_Sg - 0.05, best_Sg + 0.05, 200)
    for Sg_try in fine:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=tol)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try
    return best_err, best_np, best_Sg


# ============================================================
# PARTIAL TRACE & QUANTUM METRICS
# ============================================================
def ptrace_keep(rho_full, nq, keep):
    """Partial trace keeping qubits in 'keep' list. General N-qubit."""
    dim = 2 ** nq
    nk = len(keep)
    dk = 2 ** nk
    rho_r = np.zeros((dk, dk), dtype=complex)
    traced = [k for k in range(nq) if k not in keep]
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (nq - 1 - k)) & 1 for k in range(nq)]
            bj = [(j >> (nq - 1 - k)) & 1 for k in range(nq)]
            if all(bi[k] == bj[k] for k in traced):
                ki = sum(bi[keep[m]] << (nk - 1 - m) for m in range(nk))
                kj = sum(bj[keep[m]] << (nk - 1 - m) for m in range(nk))
                rho_r[ki, kj] += rho_full[i, j]
    return rho_r


def purity(rho):
    return np.real(np.trace(rho @ rho))


def l1_coherence(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def psi_norm(rho):
    dd = rho.shape[0]
    if dd <= 1: return 0.0
    return l1_coherence(rho) / (dd - 1)


def compute_CPsi(rho):
    """CΨ = purity * normalized L1 coherence."""
    C = purity(rho)
    Psi = psi_norm(rho)
    return C * Psi, C, Psi


def compute_theta(CPsi):
    if CPsi > 0.25:
        return np.degrees(np.arctan(np.sqrt(4 * CPsi - 1)))
    return 0.0


def concurrence_2q(rho):
    """Wootters concurrence for 2-qubit density matrix."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho.conj() @ sy2
    R = rho @ rho_tilde
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])


def von_neumann(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


def mutual_info(rho_full, nq, keep_A, keep_B):
    """I(A:B) = S(A) + S(B) - S(AB)."""
    keep_AB = sorted(set(keep_A) | set(keep_B))
    rho_A = ptrace_keep(rho_full, nq, keep_A)
    rho_B = ptrace_keep(rho_full, nq, keep_B)
    rho_AB = ptrace_keep(rho_full, nq, keep_AB)
    rho_A = (rho_A + rho_A.conj().T) / 2
    rho_B = (rho_B + rho_B.conj().T) / 2
    rho_AB = (rho_AB + rho_AB.conj().T) / 2
    return von_neumann(rho_A) + von_neumann(rho_B) - von_neumann(rho_AB)


# ============================================================
# SYSTEM CONSTRUCTION
# ============================================================
def build_mediator_system(J_A=1.0, J_AM=1.0, J_MB=1.0, J_B=1.0,
                          gammas=None):
    """Build 5-qubit mediator bridge: 0-1-2-3-4.
    Pair A: 0,1. Mediator: 2. Pair B: 3,4.
    """
    if gammas is None:
        gammas = [0.05] * N

    bonds = [(0, 1), (1, 2), (2, 3), (3, 4)]
    Js = [J_A, J_AM, J_MB, J_B]

    H = build_H_heisenberg(bonds, Js)
    L = build_L_H(H) + build_L_D_zdeph(gammas)
    Sg = sum(gammas)
    return L, Sg, H


def make_initial_state(name='bell_A_0M_pp_B'):
    """Build initial state for 5-qubit system."""
    if name == 'bell_A_0M_pp_B':
        # Bell(0,1) x |0>_2 x |+>_3 x |+>_4
        psi_A = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
        psi = np.kron(np.kron(psi_A, up), np.kron(plus, plus))
    elif name == 'bell_A_0M_bell_B':
        psi_A = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
        psi_B = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
        psi = np.kron(np.kron(psi_A, up), psi_B)
    elif name == '00000':
        psi = np.zeros(d, dtype=complex); psi[0] = 1
    elif name == '01010':
        idx = 0b01010
        psi = np.zeros(d, dtype=complex); psi[idx] = 1
    elif name == 'W5':
        psi = np.zeros(d, dtype=complex)
        for k in range(N):
            idx = 1 << (N - 1 - k)
            psi[idx] = 1 / np.sqrt(N)
    elif name == 'GHZ5':
        psi = np.zeros(d, dtype=complex)
        psi[0] = 1 / np.sqrt(2)
        psi[d - 1] = 1 / np.sqrt(2)
    else:
        raise ValueError(f"Unknown state: {name}")
    return np.outer(psi, psi.conj())


def evolve_rho(L, rho0, t):
    """Evolve density matrix to time t using matrix exponential."""
    rho0_vec = rho0.flatten()
    rho_t_vec = expm(L * t) @ rho0_vec
    rho_t = rho_t_vec.reshape(d, d)
    return (rho_t + rho_t.conj().T) / 2


# ============================================================
# PAULI BASIS (for Pi operator)
# ============================================================
all_idx = list(iprod(range(4), repeat=N))
num = 4 ** N  # 1024

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def build_Pi_pauli(idx_list, nq=N):
    n = len(idx_list)
    Pi = np.zeros((n, n), dtype=complex)
    idx_map = {idx: i for i, idx in enumerate(idx_list)}
    for b, idx_b in enumerate(idx_list):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = idx_map[mapped]
        Pi[a, b] = sign
    return Pi


def xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))


# ============================================================
# TEST 1: PALINDROME SURVIVES
# ============================================================
def run_test_1(J_A=1.0, J_AM=1.0, J_MB=1.0, J_B=1.0, gamma=0.05):
    log("=" * 70)
    log("TEST 1: DOES THE PALINDROME SURVIVE?")
    log("=" * 70)
    log()

    gammas = [gamma] * N
    L, Sg, H = build_mediator_system(J_A, J_AM, J_MB, J_B, gammas)

    log(f"N={N}, d={d}, d^2={d2}")
    log(f"J_A={J_A}, J_AM={J_AM}, J_MB={J_MB}, J_B={J_B}")
    log(f"gamma={gamma} (uniform), Sg={Sg:.4f}")
    log(f"Topology: 0 -- 1 -- M(2) -- 3 -- 4")
    log()

    t_eig = time.time()
    evals = np.linalg.eigvals(L)
    log(f"Eigendecomposition: {time.time() - t_eig:.1f}s")

    # Palindrome check
    b_err, b_np, b_Sg = best_palindrome(evals, Sg)
    log()
    log(f"Palindromic pairs: {b_np}/{d2}")
    log(f"Pairing error: {b_err:.2e}")
    log(f"Best Sg: {b_Sg:.6f} (expected {Sg:.6f})")
    log()

    if b_np == d2:
        log("PALINDROME CONFIRMED. All 1024 eigenvalues paired.")
    else:
        log(f"WARNING: Only {b_np}/{d2} pairs found!")
    log()

    # X^N parity
    XN = sx.copy()
    for k in range(1, N): XN = np.kron(XN, sx)
    XN_super = np.kron(XN, XN.conj())
    comm_xn = np.max(np.abs(XN_super @ L - L @ XN_super))
    log(f"[X^{N}, L] = {comm_xn:.2e}")
    log(f"X^{N} parity: {'CONSERVED' if comm_xn < 1e-10 else 'BROKEN'}")
    log()

    # Pi operator check
    log("--- Pi operator conjugation ---")
    Pi_pauli = build_Pi_pauli(all_idx)

    # Build Pauli basis matrices
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, N):
            m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m)

    V = np.zeros((d2, num), dtype=complex)
    for a in range(num):
        V[:, a] = pmats[a].flatten()
    V_inv = np.linalg.inv(V)

    Pi_comp = V @ Pi_pauli @ V_inv
    Pi_inv = np.linalg.inv(Pi_comp)

    Sg_candidates = np.concatenate([[Sg], np.linspace(Sg - 0.02, Sg + 0.02, 100)])
    L_norm = max(np.max(np.abs(L)), 1)
    best_pi_err = 1e30; best_pi_Sg = Sg
    for Sg_try in Sg_candidates:
        E = Pi_comp @ L @ Pi_inv + L + 2 * Sg_try * np.eye(d2)
        err = np.max(np.abs(E)) / L_norm
        if err < best_pi_err:
            best_pi_err = err; best_pi_Sg = Sg_try

    log(f"Pi conjugation error: {best_pi_err:.2e} (Sg={best_pi_Sg:.6f})")
    if best_pi_err < 1e-6:
        log("Pi WORKS for the mediator bridge.")
    else:
        log("Pi does NOT work.")
    log()

    # Spectral statistics
    rates = -evals.real
    n_ss = int(np.sum(np.abs(rates) < 1e-8))
    unique_rates = len(np.unique(np.round(rates, 6)))
    log(f"Steady states: {n_ss}")
    log(f"Distinct rates: {unique_rates}")
    log(f"Rate range: [{np.min(rates):.6f}, {np.max(rates):.6f}]")
    log()

    # Compare with N=4 chain and isolated pairs
    log("--- Comparison systems ---")
    log()

    # N=4 chain (no mediator)
    bonds4 = [(0, 1), (1, 2), (2, 3)]
    H4 = build_H_heisenberg(bonds4, [J_A, J_AM, J_MB], nq=4)
    d4 = 16; d4_2 = 256
    L4 = build_L_H(H4, 4) + build_L_D_zdeph([gamma] * 4, 4)
    evals4 = np.linalg.eigvals(L4)
    b4_err, b4_np, b4_Sg = best_palindrome(evals4, 4 * gamma)
    log(f"N=4 chain (0-1-2-3): {b4_np}/{d4_2} palindromic, err={b4_err:.2e}")

    # Isolated pairs
    bonds_iso = [(0, 1), (2, 3)]
    H_iso = build_H_heisenberg(bonds_iso, [J_A, J_B], nq=4)
    L_iso = build_L_H(H_iso, 4) + build_L_D_zdeph([gamma] * 4, 4)
    evals_iso = np.linalg.eigvals(L_iso)
    bi_err, bi_np, bi_Sg = best_palindrome(evals_iso, 4 * gamma)
    log(f"Isolated pairs (0-1, 2-3): {bi_np}/{d4_2} palindromic, err={bi_err:.2e}")
    log()

    return L, Sg, H, evals, Pi_comp, V, V_inv


# ============================================================
# TEST 2: CROSS-PAIR INFORMATION FLOW
# ============================================================
def run_test_2(L, Sg, H):
    log("=" * 70)
    log("TEST 2: CROSS-PAIR INFORMATION FLOW")
    log("=" * 70)
    log()

    t_points = np.linspace(0, 30, 60)

    initial_states = {
        'Bell_A|0>M|++>B': 'bell_A_0M_pp_B',
        'Bell_A|0>M|Bell>B': 'bell_A_0M_bell_B',
        '|00000>': '00000',
        '|01010>': '01010',
        'W5': 'W5',
        'GHZ5': 'GHZ5',
    }

    # Pair A = qubits 0,1. Pair B = qubits 3,4. M = qubit 2.
    log("Metrics: CPsi_15 (edge-edge), CPsi_24 (boundary-boundary),")
    log("         Conc_15, Conc_24, MI(A:B)")
    log()

    for sname, skey in initial_states.items():
        log(f"--- Initial state: {sname} ---")
        log()

        rho0 = make_initial_state(skey)

        # Time series
        cpsi_15 = []; cpsi_24 = []; conc_15 = []; conc_24 = []
        mi_ab = []; theta_15 = []; theta_24 = []

        for t in t_points:
            if t == 0:
                rho_t = rho0.copy()
            else:
                rho_t = evolve_rho(L, rho0, t)

            # Reduced states
            rho_15 = ptrace_keep(rho_t, N, [0, 4])   # edge pair
            rho_24 = ptrace_keep(rho_t, N, [1, 3])   # boundary pair
            rho_A = ptrace_keep(rho_t, N, [0, 1])    # pair A
            rho_B = ptrace_keep(rho_t, N, [3, 4])    # pair B

            cp15, _, _ = compute_CPsi(rho_15)
            cp24, _, _ = compute_CPsi(rho_24)
            c15 = concurrence_2q(rho_15)
            c24 = concurrence_2q(rho_24)
            mi = mutual_info(rho_t, N, [0, 1], [3, 4])

            cpsi_15.append(cp15)
            cpsi_24.append(cp24)
            conc_15.append(c15)
            conc_24.append(c24)
            mi_ab.append(mi)
            theta_15.append(compute_theta(cp15))
            theta_24.append(compute_theta(cp24))

        cpsi_15 = np.array(cpsi_15); cpsi_24 = np.array(cpsi_24)
        conc_15 = np.array(conc_15); conc_24 = np.array(conc_24)
        mi_ab = np.array(mi_ab)
        theta_15 = np.array(theta_15); theta_24 = np.array(theta_24)

        # Summary
        max_cp15 = np.max(cpsi_15)
        max_cp24 = np.max(cpsi_24)
        max_c15 = np.max(conc_15)
        max_c24 = np.max(conc_24)
        max_mi = np.max(mi_ab)

        cross_15 = np.any(cpsi_15 > 0.25)
        cross_24 = np.any(cpsi_24 > 0.25)

        log(f"  CPsi_15 max: {max_cp15:.6f} {'(> 1/4!)' if cross_15 else '(< 1/4)'}")
        log(f"  CPsi_24 max: {max_cp24:.6f} {'(> 1/4!)' if cross_24 else '(< 1/4)'}")
        log(f"  Conc_15 max: {max_c15:.6f}")
        log(f"  Conc_24 max: {max_c24:.6f}")
        log(f"  MI(A:B) max: {max_mi:.6f}")

        if cross_15 or cross_24:
            log(f"  *** CROSS-PAIR CPsi EXCEEDS 1/4 ***")
            if cross_15:
                t_cross = t_points[np.argmax(cpsi_15 > 0.25)]
                log(f"  First crossing (edge-edge) at t = {t_cross:.2f}")
            if cross_24:
                t_cross = t_points[np.argmax(cpsi_24 > 0.25)]
                log(f"  First crossing (boundary) at t = {t_cross:.2f}")

        # Time series at key points
        log()
        log(f"  {'t':>6}  {'CPsi15':>8}  {'CPsi24':>8}  {'C15':>8}  "
            f"{'C24':>8}  {'MI':>8}  {'th15':>6}  {'th24':>6}")
        log(f"  {'-'*65}")
        for idx in range(0, len(t_points), 6):
            t = t_points[idx]
            log(f"  {t:>6.1f}  {cpsi_15[idx]:>8.4f}  {cpsi_24[idx]:>8.4f}  "
                f"{conc_15[idx]:>8.4f}  {conc_24[idx]:>8.4f}  "
                f"{mi_ab[idx]:>8.4f}  {theta_15[idx]:>6.1f}  {theta_24[idx]:>6.1f}")
        log()

    # QST Fidelity
    log("--- QST Fidelity: qubit 0 -> qubit 4 ---")
    log()

    mub_states = [up, dn, plus, (up - dn) / np.sqrt(2),
                  (up + 1j * dn) / np.sqrt(2), (up - 1j * dn) / np.sqrt(2)]
    mub_names = ['|0>', '|1>', '|+>', '|->', '|+i>', '|-i>']

    log(f"  {'Input':>6}  {'t_max_F':>8}  {'max_F':>8}  {'F(t=1)':>8}  "
        f"{'F(t=5)':>8}")
    log(f"  {'-'*42}")

    best_avg_F = 0; best_t = 0
    fid_vs_t = np.zeros(len(t_points))

    for mi, psi_in in enumerate(mub_states):
        # Prepare: psi_in on qubit 0, rest |0>
        psi_full = psi_in.copy()
        for k in range(1, N):
            psi_full = np.kron(psi_full, up)
        rho0 = np.outer(psi_full, psi_full.conj())

        fids = []
        for t in t_points:
            if t == 0:
                rho_t = rho0.copy()
            else:
                rho_t = evolve_rho(L, rho0, t)
            rho_5 = ptrace_keep(rho_t, N, [4])  # qubit 4 (far end)
            F = np.real(psi_in.conj() @ rho_5 @ psi_in)
            fids.append(F)

        fids = np.array(fids)
        fid_vs_t += fids
        t_max = t_points[np.argmax(fids)]
        F_max = np.max(fids)
        F_1 = fids[np.argmin(np.abs(t_points - 1))]
        F_5 = fids[np.argmin(np.abs(t_points - 5))]

        log(f"  {mub_names[mi]:>6}  {t_max:>8.2f}  {F_max:>8.4f}  "
            f"{F_1:>8.4f}  {F_5:>8.4f}")

    fid_vs_t /= len(mub_states)
    t_best = t_points[np.argmax(fid_vs_t)]
    F_best = np.max(fid_vs_t)
    log()
    log(f"  Average QST fidelity: {F_best:.4f} at t = {t_best:.2f}")
    log(f"  (Random baseline: 0.5)")
    log()


# ============================================================
# TEST 3: THREE CONDITIONS SCALE UP
# ============================================================
def run_test_3(gamma=0.05):
    log("=" * 70)
    log("TEST 3: THREE CONDITIONS SCALE UP")
    log("=" * 70)
    log()

    rho0 = make_initial_state('bell_A_0M_pp_B')
    t_measure = 5.0

    # (a) Coupling asymmetry: sweep J_AM and J_MB
    log("--- (a) Coupling asymmetry: J_AM vs J_MB ---")
    log()
    log(f"  {'J_AM':>6}  {'J_MB':>6}  {'CPsi15':>8}  {'CPsi24':>8}  "
        f"{'Conc15':>8}  {'MI':>8}")
    log(f"  {'-'*50}")

    J_range = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    for J_AM in J_range:
        for J_MB in [1.0]:  # fix J_MB, sweep J_AM
            L, Sg, H = build_mediator_system(1.0, J_AM, J_MB, 1.0,
                                              [gamma] * N)
            rho_t = evolve_rho(L, rho0, t_measure)
            rho_15 = ptrace_keep(rho_t, N, [0, 4])
            rho_24 = ptrace_keep(rho_t, N, [1, 3])
            cp15, _, _ = compute_CPsi(rho_15)
            cp24, _, _ = compute_CPsi(rho_24)
            c15 = concurrence_2q(rho_15)
            mi = mutual_info(rho_t, N, [0, 1], [3, 4])
            log(f"  {J_AM:>6.2f}  {J_MB:>6.2f}  {cp15:>8.4f}  {cp24:>8.4f}  "
                f"{c15:>8.4f}  {mi:>8.4f}")
    log()

    # Same but sweep J_MB
    log(f"  {'J_AM':>6}  {'J_MB':>6}  {'CPsi15':>8}  {'CPsi24':>8}  "
        f"{'Conc15':>8}  {'MI':>8}")
    log(f"  {'-'*50}")
    for J_MB in J_range:
        L, Sg, H = build_mediator_system(1.0, 1.0, J_MB, 1.0,
                                          [gamma] * N)
        rho_t = evolve_rho(L, rho0, t_measure)
        rho_15 = ptrace_keep(rho_t, N, [0, 4])
        rho_24 = ptrace_keep(rho_t, N, [1, 3])
        cp15, _, _ = compute_CPsi(rho_15)
        cp24, _, _ = compute_CPsi(rho_24)
        c15 = concurrence_2q(rho_15)
        mi = mutual_info(rho_t, N, [0, 1], [3, 4])
        log(f"  {1.0:>6.2f}  {J_MB:>6.2f}  {cp15:>8.4f}  {cp24:>8.4f}  "
            f"{c15:>8.4f}  {mi:>8.4f}")
    log()

    # (b) Mediator noise sweep
    log("--- (b) Mediator noise: gamma_M sweep ---")
    log()
    log(f"  {'gM':>6}  {'CPsi15':>8}  {'CPsi24':>8}  {'Conc15':>8}  "
        f"{'MI':>8}  {'Paired':>8}")
    log(f"  {'-'*52}")

    for gM in np.linspace(0, 0.5, 15):
        gammas = [gamma, gamma, gM, gamma, gamma]
        L, Sg, H = build_mediator_system(gammas=gammas)
        evals = np.linalg.eigvals(L)
        b_err, b_np, b_Sg = best_palindrome(evals, Sg)

        rho_t = evolve_rho(L, rho0, t_measure)
        rho_15 = ptrace_keep(rho_t, N, [0, 4])
        rho_24 = ptrace_keep(rho_t, N, [1, 3])
        cp15, _, _ = compute_CPsi(rho_15)
        cp24, _, _ = compute_CPsi(rho_24)
        c15 = concurrence_2q(rho_15)
        mi = mutual_info(rho_t, N, [0, 1], [3, 4])
        log(f"  {gM:>6.3f}  {cp15:>8.4f}  {cp24:>8.4f}  {c15:>8.4f}  "
            f"{mi:>8.4f}  {b_np:>4}/{d2}")
    log()


# ============================================================
# TEST 4: STANDING WAVE
# ============================================================
def run_test_4(L, H, gamma=0.05):
    log("=" * 70)
    log("TEST 4: STANDING WAVE ACROSS THE BRIDGE")
    log("=" * 70)
    log()

    rho0 = make_initial_state('bell_A_0M_pp_B')
    t_points = np.linspace(0, 20, 40)

    # Two-point correlators
    pairs = [(0, 4), (1, 3), (0, 2), (2, 4), (0, 1), (3, 4)]
    pair_names = ['15(edge)', '24(bnd)', '1M', 'M4', '12(A)', '45(B)']

    hdr = f"  {'t':>6}"
    for pn in pair_names:
        hdr += f"  {'XX_'+pn:>10}  {'ZZ_'+pn:>10}"
    log(hdr)
    log(f"  {'-'*130}")

    for t in t_points[::4]:
        rho_t = evolve_rho(L, rho0, t) if t > 0 else rho0.copy()
        row = f"  {t:>6.1f}"
        for i, j in pairs:
            XX = site_op(sx, i) @ site_op(sx, j)
            xx_val = np.real(np.trace(XX @ rho_t))
            ZZ = site_op(sz, i) @ site_op(sz, j)
            zz_val = np.real(np.trace(ZZ @ rho_t))
            row += f"  {xx_val:>10.4f}  {zz_val:>10.4f}"
        log(row)
    log()


# ============================================================
# TEST 5: DECODER
# ============================================================
def run_test_5(L, Sg, gamma=0.05):
    log("=" * 70)
    log("TEST 5: DECODER ACROSS THE BRIDGE")
    log("=" * 70)
    log()

    # Fisher information for J_AM sensitivity
    log("--- Fisher information for J_AM at t=5 ---")
    log()

    dJ = 1e-5
    L_ref = L
    L_pert, _, _ = build_mediator_system(J_AM=1.0 + dJ,
                                          gammas=[gamma] * N)
    t_fish = 5.0

    test_states = {
        '|00000>': '00000',
        '|01010>': '01010',
        'Bell_A|0>M|++>B': 'bell_A_0M_pp_B',
        'W5': 'W5',
        'GHZ5': 'GHZ5',
    }

    log(f"  {'State':>20}  {'F(J_AM)':>12}  {'Sensitivity':>12}")
    log(f"  {'-'*48}")

    for sname, skey in test_states.items():
        rho0 = make_initial_state(skey)
        rho0_vec = rho0.flatten()

        rho_ref = (expm(L_ref * t_fish) @ rho0_vec).reshape(d, d)
        rho_pert = (expm(L_pert * t_fish) @ rho0_vec).reshape(d, d)
        rho_ref = (rho_ref + rho_ref.conj().T) / 2
        rho_pert = (rho_pert + rho_pert.conj().T) / 2

        p0 = np.real(np.diag(rho_ref))
        pp = np.real(np.diag(rho_pert))
        dp = (pp - p0) / dJ

        fisher = 0
        for k in range(d):
            if p0[k] > 1e-12:
                fisher += dp[k] ** 2 / p0[k]

        sensitivity = np.linalg.norm(rho_pert - rho_ref) / dJ

        log(f"  {sname:>20}  {fisher:>12.4f}  {sensitivity:>12.4f}")

    log()

    # Response matrix for per-site gamma
    log("--- Response matrix: sensitivity to per-site gamma ---")
    log()

    evals_ref, R_ref = np.linalg.eig(L)
    gammas_ref = [gamma] * N

    rho0 = make_initial_state('bell_A_0M_pp_B')
    rho0_vec = rho0.flatten()
    R_inv_ref = np.linalg.inv(R_ref)
    coeffs_ref = R_inv_ref @ rho0_vec

    # Build palindromic pairs
    b_err, b_np, b_Sg = best_palindrome(evals_ref, Sg, tol=1e-6)
    pairs = []
    used = set()
    for k in range(d2):
        if k in used: continue
        target = -(evals_ref[k] + 2 * b_Sg)
        dists = np.abs(evals_ref - target)
        for u in used: dists[u] = 1e30
        j = np.argmin(dists)
        if dists[j] < 1e-4:
            pairs.append((k, j))
            used.add(k); used.add(j)

    n_pairs = len(pairs)
    log(f"Palindromic pairs for response matrix: {n_pairs}")

    if n_pairs > 0:
        dg = 1e-4
        response = np.zeros((n_pairs, N))

        for site in range(N):
            gammas_p = list(gammas_ref)
            gammas_p[site] += dg
            L_ps, _, _ = build_mediator_system(gammas=gammas_p)
            evals_ps, R_ps = np.linalg.eig(L_ps)
            R_inv_ps = np.linalg.inv(R_ps)
            coeffs_ps = R_inv_ps @ rho0_vec

            for pi, (k, j) in enumerate(pairs):
                amp_ref = abs(coeffs_ref[k]) + (abs(coeffs_ref[j]) if j != k else 0)
                dk = np.argmin(np.abs(evals_ps - evals_ref[k]))
                amp_pert = abs(coeffs_ps[dk])
                if j != k:
                    dj = np.argmin(np.abs(evals_ps - evals_ref[j]))
                    amp_pert += abs(coeffs_ps[dj])
                response[pi, site] = (amp_pert - amp_ref) / dg

        U, sv, Vt = np.linalg.svd(response, full_matrices=False)
        n_indep = int(np.sum(sv > 1e-6 * sv[0])) if len(sv) > 0 else 0

        log(f"Response matrix: {n_pairs} pairs x {N} sites")
        log(f"SVD singular values: {', '.join(f'{s:.4f}' for s in sv[:min(N+2, len(sv))])}")
        log(f"Independent noise directions: {n_indep}/{N}")

        if n_indep == N:
            log("FULL RANK: All per-site gammas readable through the bridge.")
        else:
            log(f"RANK DEFICIENT: Only {n_indep}/{N} gamma parameters visible.")
    log()


# ============================================================
# TEST 6: ROBUSTNESS AGAINST DIRECT LEAKAGE
# ============================================================
def run_test_6(gamma=0.05):
    log("=" * 70)
    log("TEST 6: ROBUSTNESS AGAINST DIRECT LEAKAGE")
    log("=" * 70)
    log()

    log("Adding XZ cross-dissipation between qubits 1 and 3 (boundary)")
    log("on top of the mediator bridge.")
    log()

    log(f"  {'epsilon':>10}  {'Paired':>8}  {'Err':>10}  {'MI(t=5)':>10}  "
        f"{'CPsi24':>10}")
    log(f"  {'-'*55}")

    rho0 = make_initial_state('bell_A_0M_pp_B')

    epsilons = np.concatenate([[0], np.logspace(-5, -1, 20)])

    for eps in epsilons:
        gammas = [gamma] * N
        L, Sg, H = build_mediator_system(gammas=gammas)

        if eps > 0:
            # Add cross-boundary dissipation between qubits 1 and 3
            F1 = site_op(sx, 1) @ site_op(sz, 3)
            F2 = site_op(sy, 1) @ site_op(sz, 3)
            F3 = site_op(sz, 1) @ site_op(sx, 3)
            F4 = site_op(sz, 1) @ site_op(sy, 3)
            L = add_lindblad_jump(L, F1, eps)
            L = add_lindblad_jump(L, F2, eps)
            L = add_lindblad_jump(L, F3, eps)
            L = add_lindblad_jump(L, F4, eps)

        evals = np.linalg.eigvals(L)
        b_err, b_np, b_Sg = best_palindrome(evals, Sg)

        rho_t = evolve_rho(L, rho0, 5.0)
        mi = mutual_info(rho_t, N, [0, 1], [3, 4])
        rho_24 = ptrace_keep(rho_t, N, [1, 3])
        cp24, _, _ = compute_CPsi(rho_24)

        log(f"  {eps:>10.5f}  {b_np:>4}/{d2}  {b_err:>10.2e}  "
            f"{mi:>10.4f}  {cp24:>10.4f}")

    log()


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("The Mediator Bridge: S Was Always the Answer")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, d={d}, d^2={d2}")
    log()

    # Test 1
    t1 = time.time()
    L, Sg, H, evals, Pi_comp, V, V_inv = run_test_1()
    log(f"[Test 1 completed in {time.time() - t1:.1f}s]")
    log()

    # Test 2
    t2 = time.time()
    run_test_2(L, Sg, H)
    log(f"[Test 2 completed in {time.time() - t2:.1f}s]")
    log()

    # Test 3
    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time() - t3:.1f}s]")
    log()

    # Test 4
    t4 = time.time()
    run_test_4(L, H)
    log(f"[Test 4 completed in {time.time() - t4:.1f}s]")
    log()

    # Test 5
    t5 = time.time()
    run_test_5(L, Sg)
    log(f"[Test 5 completed in {time.time() - t5:.1f}s]")
    log()

    # Test 6
    t6 = time.time()
    run_test_6()
    log(f"[Test 6 completed in {time.time() - t6:.1f}s]")
    log()

    # Summary
    log("=" * 70)
    log("OVERALL SUMMARY")
    log("=" * 70)
    log()
    log("Test 1: Palindrome verification (theorem check)")
    log("Test 2: Cross-pair information flow (CPsi, concurrence, MI, QST)")
    log("Test 3: Three conditions (coupling asymmetry, mediator noise)")
    log("Test 4: Standing wave (correlators across bridge)")
    log("Test 5: Decoder (Fisher info, response matrix)")
    log("Test 6: Robustness (direct leakage tolerance)")
    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
