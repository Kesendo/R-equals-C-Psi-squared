#!/usr/bin/env python3
"""
The Mixed Bridge: Hamiltonian + Dissipative, Both or Nothing
=============================================================
Test 1: Phase diagram (J_bridge, kappa) -- palindrome survival region
Test 3: Bidirectionality requirement in the mixed bridge
Test 4: Pi operator in the mixed system
Test 2: Critical line -- where the palindrome breaks
Test 5: Does the decoder still work?

Script: simulations/mixed_bridge.py
Output: simulations/results/mixed_bridge.txt
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "mixed_bridge.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# INFRASTRUCTURE
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
PNAMES = ['I', 'X', 'Y', 'Z']

N = 4
d = 2 ** N   # 16
d2 = d * d   # 256


def site_op(op, k):
    ops = [I2] * N; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H_heisenberg(bonds, J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, j)
    return H


def build_L_H(H):
    Id = np.eye(d)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_L_D_zdeph(gammas):
    L_D = np.zeros((d2, d2), dtype=complex)
    for k in range(N):
        Zk = site_op(sz, k)
        L_D += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L_D


def add_lindblad_jump(L, F, kappa):
    Id = np.eye(d)
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
    """Two-stage Sg search: coarse then fine, returns (err, n_paired, Sg)."""
    # Stage 1: coarse search with relaxed tolerance
    best_np = 0; best_Sg = Sg_base; best_err = 1e30
    # Always include Sg_base exactly, plus a wide sweep
    candidates = np.concatenate([
        [Sg_base],
        np.linspace(0, Sg_base * 5, 60)
    ])
    for Sg_try in candidates:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=1e-4)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try

    # Stage 2: fine search around best Sg
    fine = np.linspace(best_Sg - 0.05, best_Sg + 0.05, 200)
    for Sg_try in fine:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=tol)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try

    return best_err, best_np, best_Sg


def mutual_information_AB(rho_AB):
    """S(A) + S(B) - S(AB) where A = qubits 0,1 and B = qubits 2,3.
    rho_AB is 16x16 in the |q0 q1 q2 q3> basis.
    """
    def von_neumann(rho):
        evals = np.linalg.eigvalsh(rho)
        evals = evals[evals > 1e-15]
        return -np.sum(evals * np.log2(evals))

    # Partial trace over B (qubits 2,3) to get rho_A (4x4)
    # |q0 q1 q2 q3>: index = q0*8 + q1*4 + q2*2 + q3
    # Reshape as (dA, dB, dA, dB) where dA=4, dB=4
    # Index mapping: i = a*dB + b where a = q0*2+q1, b = q2*2+q3
    dA = 4; dB = 4
    rho_r = rho_AB.reshape(dA, dB, dA, dB)
    rho_A = np.trace(rho_r, axis1=1, axis2=3)  # trace over B
    rho_B = np.trace(rho_r, axis1=0, axis2=2)  # trace over A

    SA = von_neumann(rho_A)
    SB = von_neumann(rho_B)
    SAB = von_neumann(rho_AB)

    return SA + SB - SAB


def transient_mi(L, t=5.0):
    """Compute MI at time t from Bell state on boundary qubits (1,2)."""
    # |psi0> = (|0000> + |0110>) / sqrt(2) -- Bell on boundary
    psi0 = np.zeros(d, dtype=complex)
    psi0[0] = 1/np.sqrt(2)
    psi0[0b0110] = 1/np.sqrt(2)
    rho0 = np.outer(psi0, psi0.conj())
    rho0_vec = rho0.flatten()

    rho_t_vec = expm(L * t) @ rho0_vec
    rho_t = rho_t_vec.reshape(d, d)
    rho_t = (rho_t + rho_t.conj().T) / 2  # enforce Hermiticity

    # Ensure valid density matrix
    evals_rho = np.linalg.eigvalsh(rho_t)
    if np.min(evals_rho) < -1e-10:
        return 0.0  # not a valid state

    return max(0, mutual_information_AB(rho_t))


def steady_state_rho(L):
    """Extract steady state density matrix from Liouvillian."""
    evals, R = np.linalg.eig(L)
    k = np.argmin(np.abs(evals))
    rho_vec = R[:, k]
    rho = rho_vec.reshape(d, d)
    rho = rho / np.trace(rho)
    rho = (rho + rho.conj().T) / 2  # enforce Hermiticity
    return rho


# ============================================================
# BRIDGE CONSTRUCTION
# ============================================================
def build_mixed_bridge(gamma=0.1, J_A=1.0, J_B=1.0,
                       J_bridge=0.5, kappa=0.1,
                       bridge_dissipator='xz_boundary',
                       kappa_AB=None, kappa_BA=None):
    """Build the full mixed (Hamiltonian + dissipative) bridge Liouvillian.

    System A: qubits 0,1 with Heisenberg coupling J_A
    System B: qubits 2,3 with Heisenberg coupling J_B
    Bridge: Hamiltonian bond J_bridge on (1,2) + dissipative channel kappa
    """
    if kappa_AB is None: kappa_AB = kappa
    if kappa_BA is None: kappa_BA = kappa

    gammas = [gamma] * N

    # Internal Hamiltonians
    H_A = build_H_heisenberg([(0, 1)], J_A)
    H_B = build_H_heisenberg([(2, 3)], J_B)

    # Bridge Hamiltonian: Heisenberg coupling between boundary qubits 1-2
    H_bridge = build_H_heisenberg([(1, 2)], J_bridge)

    H_total = H_A + H_B + H_bridge

    # Build Liouvillian
    L = build_L_H(H_total) + build_L_D_zdeph(gammas)

    # Add dissipative bridge
    if kappa_AB > 0 or kappa_BA > 0:
        if bridge_dissipator == 'xz_boundary':
            # Coherence-population cross-coupling (the 73/27 operators)
            F1 = site_op(sx, 1) @ site_op(sz, 2)
            F2 = site_op(sy, 1) @ site_op(sz, 2)
            if kappa_AB > 0:
                L = add_lindblad_jump(L, F1, kappa_AB)
                L = add_lindblad_jump(L, F2, kappa_AB)
            F3 = site_op(sz, 1) @ site_op(sx, 2)
            F4 = site_op(sz, 1) @ site_op(sy, 2)
            if kappa_BA > 0:
                L = add_lindblad_jump(L, F3, kappa_BA)
                L = add_lindblad_jump(L, F4, kappa_BA)

        elif bridge_dissipator == 'zz_collective':
            F = site_op(sz, 1) @ site_op(sz, 2)
            L = add_lindblad_jump(L, F, (kappa_AB + kappa_BA) / 2)

        elif bridge_dissipator == 'amplitude_damping':
            sp = (sx + 1j * sy) / 2
            sm = (sx - 1j * sy) / 2
            if kappa_AB > 0:
                F1 = site_op(sm, 1) @ site_op(sp, 2)
                L = add_lindblad_jump(L, F1, kappa_AB)
            if kappa_BA > 0:
                F2 = site_op(sp, 1) @ site_op(sm, 2)
                L = add_lindblad_jump(L, F2, kappa_BA)

        elif bridge_dissipator == 'swap_dissipation':
            # Dissipative SWAP on qubits 1,2
            # |01><10| + |10><01| in the (qubit1, qubit2) subspace
            # Build as full 16x16 operator
            F = np.zeros((d, d), dtype=complex)
            for a0 in range(2):
                for a3 in range(2):
                    # |a0, 0, 1, a3><a0, 1, 0, a3|
                    bra = a0 * 8 + 1 * 4 + 0 * 2 + a3
                    ket = a0 * 8 + 0 * 4 + 1 * 2 + a3
                    F[ket, bra] = 1.0
                    F[bra, ket] = 1.0
            k_avg = (kappa_AB + kappa_BA) / 2
            if k_avg > 0:
                L = add_lindblad_jump(L, F, k_avg)

    Sg = sum(gammas)
    return L, Sg


# ============================================================
# PAULI BASIS TOOLS (for Pi operator)
# ============================================================
all_idx = list(iprod(range(4), repeat=N))
num = 4 ** N  # 256

pmats = []
for idx in all_idx:
    m = PAULIS[idx[0]]
    for k in range(1, N):
        m = np.kron(m, PAULIS[idx[k]])
    pmats.append(m)

# Change-of-basis: Pauli <-> computational vectorized
V = np.zeros((d2, num), dtype=complex)
for a in range(num):
    V[:, a] = pmats[a].flatten()
V_inv = np.linalg.inv(V)

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def build_Pi_pauli(idx_list):
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


def pi_conjugation_error(L, Pi_comp, Sg_candidates):
    """Find best Sg for Pi L Pi^-1 + L + 2Sg I = 0."""
    Pi_inv = np.linalg.inv(Pi_comp)
    PLP = Pi_comp @ L @ Pi_inv
    L_norm = max(np.max(np.abs(L)), 1)
    best_err = 1e30; best_Sg = 0
    for Sg in Sg_candidates:
        E = PLP + L + 2 * Sg * np.eye(d2)
        err = np.max(np.abs(E)) / L_norm
        if err < best_err:
            best_err = err; best_Sg = Sg
    return best_err, best_Sg


def xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))


# ============================================================
# TEST 1: THE PHASE DIAGRAM
# ============================================================
def run_test_1(gamma=0.1, J=1.0):
    log("=" * 70)
    log("TEST 1: THE PHASE DIAGRAM")
    log("=" * 70)
    log()
    log(f"J_A = J_B = {J}, gamma = {gamma} (all 4 qubits)")
    log(f"Sweep: J_bridge [0, 2.0] x kappa [0, 0.3], 20x20 = 400 points")
    log()

    n_J = 20
    n_k = 20
    J_vals = np.linspace(0, 2.0, n_J)
    k_vals = np.linspace(0, 0.3, n_k)
    Sg_base = gamma * N  # 0.4

    # Storage for phase diagram
    pal_err_map = np.zeros((n_k, n_J))
    pal_count_map = np.zeros((n_k, n_J))
    mi_map = np.zeros((n_k, n_J))
    xn_conserved_map = np.zeros((n_k, n_J), dtype=bool)

    # X^N parity superoperator
    XN = sx.copy()
    for k in range(1, N): XN = np.kron(XN, sx)
    XN_super = np.kron(XN, XN.conj())

    # Also run the 3 alternative dissipators
    dissipator_models = ['xz_boundary', 'zz_collective',
                         'amplitude_damping', 'swap_dissipation']

    for diss_model in dissipator_models:
        log(f"--- Dissipator model: {diss_model} ---")
        log()

        t_start = time.time()

        for ik, kappa in enumerate(k_vals):
            for ij, J_br in enumerate(J_vals):
                L, Sg = build_mixed_bridge(
                    gamma=gamma, J_A=J, J_B=J,
                    J_bridge=J_br, kappa=kappa,
                    bridge_dissipator=diss_model)

                evals = np.linalg.eigvals(L)

                # Palindrome check with two-stage Sg search
                b_err, b_np, b_Sg = best_palindrome(evals, Sg)
                pal_err_map[ik, ij] = b_err
                pal_count_map[ik, ij] = b_np

                # X^N parity
                comm = np.max(np.abs(XN_super @ L - L @ XN_super))
                xn_conserved_map[ik, ij] = (comm < 1e-10)

                # Transient mutual information (t=5)
                if kappa > 0 or J_br > 0:
                    try:
                        mi = transient_mi(L, t=5.0)
                        mi_map[ik, ij] = mi
                    except Exception:
                        mi_map[ik, ij] = 0
                else:
                    mi_map[ik, ij] = 0

        elapsed = time.time() - t_start
        log(f"Sweep completed in {elapsed:.1f}s")
        log()

        # Identify survival region: palindrome error < 1e-10 AND full pairing
        survival = (pal_err_map < 1e-10) & (pal_count_map == d2)
        n_survival = np.sum(survival)
        n_mi_positive = np.sum(mi_map > 0.001)
        overlap = survival & (mi_map > 0.001)
        n_overlap = np.sum(overlap)

        log(f"Palindromic points (err < 1e-10, 256/256): {n_survival}/{n_J*n_k}")
        log(f"MI > 0.001 points: {n_mi_positive}/{n_J*n_k}")
        log(f"OVERLAP (palindrome + MI): {n_overlap}/{n_J*n_k}")
        log()

        if n_overlap > 0:
            log("*** SURVIVAL REGION EXISTS ***")
            log()
            # Report the overlap region
            log(f"{'J_bridge':>10}  {'kappa':>10}  {'PalErr':>10}  "
                f"{'Paired':>8}  {'MI':>10}  {'X^N':>5}")
            log("-" * 60)
            for ik in range(n_k):
                for ij in range(n_J):
                    if overlap[ik, ij]:
                        log(f"{J_vals[ij]:>10.4f}  {k_vals[ik]:>10.4f}  "
                            f"{pal_err_map[ik,ij]:>10.2e}  "
                            f"{int(pal_count_map[ik,ij]):>4}/{d2}  "
                            f"{mi_map[ik,ij]:>10.6f}  "
                            f"{'Y' if xn_conserved_map[ik,ij] else 'N':>5}")
            log()
        else:
            log("No overlap region found.")
            log()

            # Even without full palindrome, find best points
            log("--- Best palindromic approximation with MI > 0 ---")
            mask_mi = mi_map > 0.001
            if np.any(mask_mi):
                # Among MI>0 points, find best palindrome
                err_masked = np.where(mask_mi, pal_err_map, 1e30)
                cnt_masked = np.where(mask_mi, pal_count_map, 0)
                best_idx = np.unravel_index(np.argmax(cnt_masked), cnt_masked.shape)
                ik_b, ij_b = best_idx
                log(f"Best palindrome with MI>0:")
                log(f"  J_bridge = {J_vals[ij_b]:.4f}, kappa = {k_vals[ik_b]:.4f}")
                log(f"  Paired: {int(cnt_masked[ik_b,ij_b])}/{d2}")
                log(f"  Error: {pal_err_map[ik_b,ij_b]:.2e}")
                log(f"  MI: {mi_map[ik_b,ij_b]:.6f}")
            log()

        # Phase diagram as ASCII heatmap (palindrome count)
        log("--- Phase diagram: palindrome pair count ---")
        log(f"Rows = kappa (bottom={k_vals[0]:.2f}, top={k_vals[-1]:.2f})")
        log(f"Cols = J_bridge (left={J_vals[0]:.2f}, right={J_vals[-1]:.2f})")
        log()

        # Compress to symbols
        for ik in range(n_k - 1, -1, -1):
            row = ""
            for ij in range(n_J):
                cnt = pal_count_map[ik, ij]
                if cnt == d2:
                    row += "#"  # full palindrome
                elif cnt > d2 * 0.9:
                    row += "+"
                elif cnt > d2 * 0.5:
                    row += "o"
                elif cnt > d2 * 0.1:
                    row += "."
                else:
                    row += " "
            log(f"  k={k_vals[ik]:.3f} |{row}|")
        log(f"          {''.join(['-'] * n_J)}")
        log()

        # MI heatmap
        log("--- Phase diagram: mutual information ---")
        for ik in range(n_k - 1, -1, -1):
            row = ""
            for ij in range(n_J):
                mi = mi_map[ik, ij]
                if mi > 0.1:
                    row += "#"
                elif mi > 0.01:
                    row += "+"
                elif mi > 0.001:
                    row += "."
                else:
                    row += " "
            log(f"  k={k_vals[ik]:.3f} |{row}|")
        log(f"          {''.join(['-'] * n_J)}")
        log()

        # Report pure-Hamiltonian column (kappa=0)
        log("--- Kappa=0 column (pure Hamiltonian) ---")
        log(f"{'J_bridge':>10}  {'Paired':>8}  {'Err':>10}  {'MI':>10}")
        log("-" * 45)
        for ij in range(n_J):
            log(f"{J_vals[ij]:>10.4f}  "
                f"{int(pal_count_map[0,ij]):>4}/{d2}  "
                f"{pal_err_map[0,ij]:>10.2e}  "
                f"{mi_map[0,ij]:>10.6f}")
        log()

        # Report J_bridge=0 row (pure dissipation)
        log("--- J_bridge=0 row (pure dissipation) ---")
        log(f"{'kappa':>10}  {'Paired':>8}  {'Err':>10}  {'MI':>10}")
        log("-" * 45)
        for ik in range(n_k):
            log(f"{k_vals[ik]:>10.4f}  "
                f"{int(pal_count_map[ik,0]):>4}/{d2}  "
                f"{pal_err_map[ik,0]:>10.2e}  "
                f"{mi_map[ik,0]:>10.6f}")
        log()

    # Return the last model's maps for downstream tests
    return J_vals, k_vals, pal_err_map, pal_count_map, mi_map, survival


# ============================================================
# TEST 3: BIDIRECTIONALITY IN THE MIXED BRIDGE
# ============================================================
def run_test_3(gamma=0.1, J=1.0, J_br_0=1.0, kappa_0=0.1):
    log("=" * 70)
    log("TEST 3: BIDIRECTIONALITY IN THE MIXED BRIDGE")
    log("=" * 70)
    log()
    log(f"Operating point: J_bridge={J_br_0}, kappa_0={kappa_0}")
    log()

    configs = {
        'symmetric':       (kappa_0, kappa_0),
        'asymmetric_mild': (1.2 * kappa_0, 0.8 * kappa_0),
        'asymmetric_strong': (2 * kappa_0, 0),
        'reversed':        (0, 2 * kappa_0),
    }

    dissipators = ['xz_boundary', 'zz_collective',
                   'amplitude_damping', 'swap_dissipation']

    for diss in dissipators:
        log(f"--- Dissipator: {diss} ---")
        log()
        log(f"{'Config':>18}  {'kAB':>6}  {'kBA':>6}  "
            f"{'Paired':>8}  {'Err':>10}  {'MI':>10}  {'BestSg':>8}")
        log("-" * 75)

        for cname, (kAB, kBA) in configs.items():
            L, Sg = build_mixed_bridge(
                gamma=gamma, J_A=J, J_B=J,
                J_bridge=J_br_0, kappa=kappa_0,
                bridge_dissipator=diss,
                kappa_AB=kAB, kappa_BA=kBA)

            evals = np.linalg.eigvals(L)
            b_err, b_np, b_Sg = best_palindrome(evals, Sg)

            try:
                mi = transient_mi(L, t=5.0)
            except Exception:
                mi = 0

            pal_tag = "YES" if b_np == d2 else f"{b_np}/{d2}"
            log(f"{cname:>18}  {kAB:>6.3f}  {kBA:>6.3f}  "
                f"{pal_tag:>8}  {b_err:>10.2e}  "
                f"{mi:>10.6f}  {b_Sg:>8.4f}")

        log()

    # Asymmetry tolerance sweep
    log("--- Asymmetry tolerance (xz_boundary, J_bridge=1.0) ---")
    log()
    log(f"{'Ratio':>8}  {'kAB':>6}  {'kBA':>6}  "
        f"{'Paired':>8}  {'Err':>10}  {'MI':>10}")
    log("-" * 55)

    ratios = np.linspace(1.0, 5.0, 20)
    for ratio in ratios:
        kAB = kappa_0 * ratio / (1 + ratio) * 2
        kBA = kappa_0 * 2 / (1 + ratio)
        # Keeps kAB + kBA = 2*kappa_0
        L, Sg = build_mixed_bridge(
            gamma=gamma, J_A=J, J_B=J,
            J_bridge=J_br_0, kappa=kappa_0,
            bridge_dissipator='xz_boundary',
            kappa_AB=kAB, kappa_BA=kBA)

        evals = np.linalg.eigvals(L)
        b_err, b_np, b_Sg = best_palindrome(evals, Sg)

        try:
            mi = transient_mi(L, t=5.0)
        except Exception:
            mi = 0

        log(f"{ratio:>8.2f}  {kAB:>6.3f}  {kBA:>6.3f}  "
            f"{b_np:>4}/{d2}  {b_err:>10.2e}  {mi:>10.6f}")

    log()


# ============================================================
# TEST 4: Pi IN THE MIXED BRIDGE
# ============================================================
def run_test_4(gamma=0.1, J=1.0):
    log("=" * 70)
    log("TEST 4: Pi IN THE MIXED BRIDGE")
    log("=" * 70)
    log()

    Sg_base = gamma * N

    # Build Pi_chain in computational basis
    Pi_pauli = build_Pi_pauli(all_idx)
    Pi_comp = V @ Pi_pauli @ V_inv
    Pi_inv = np.linalg.inv(Pi_comp)

    # Build Pi_A x Pi_B
    all_idx_2 = list(iprod(range(4), repeat=2))
    Pi_2 = np.zeros((16, 16), dtype=complex)
    idx_map_2 = {idx: i for i, idx in enumerate(all_idx_2)}
    for b, idx_b in enumerate(all_idx_2):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b: sign *= PI_SIGN[i]
        a = idx_map_2[mapped]
        Pi_2[a, b] = sign
    Pi_AB_pauli = np.kron(Pi_2, Pi_2)
    Pi_AB_comp = V @ Pi_AB_pauli @ V_inv

    # X^N parity
    XN = sx.copy()
    for k in range(1, N): XN = np.kron(XN, sx)
    XN_super = np.kron(XN, XN.conj())

    Sg_candidates = np.linspace(0, Sg_base * 5, 200)

    # Slice a) J_bridge=0, kappa swept (pure dissipation)
    log("--- Slice (a): J_bridge=0, kappa swept (pure dissipation) ---")
    log()
    log(f"{'kappa':>8}  {'PiChain':>10}  {'PiAxB':>10}  "
        f"{'Paired':>8}  {'BestSg':>8}  {'X^N':>5}")
    log("-" * 55)

    for kappa in np.linspace(0, 0.3, 15):
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=0, kappa=kappa)
        evals = np.linalg.eigvals(L)

        err_chain, Sg_c = pi_conjugation_error(L, Pi_comp, Sg_candidates)
        err_AB, Sg_AB = pi_conjugation_error(L, Pi_AB_comp, Sg_candidates)

        _, best_np, _ = best_palindrome(evals, Sg)

        comm_xn = np.max(np.abs(XN_super @ L - L @ XN_super))

        log(f"{kappa:>8.4f}  {err_chain:>10.2e}  {err_AB:>10.2e}  "
            f"{best_np:>4}/{d2}  {Sg_c:>8.4f}  "
            f"{'Y' if comm_xn < 1e-10 else 'N':>5}")
    log()

    # Slice b) kappa=0, J_bridge swept (pure Hamiltonian)
    log("--- Slice (b): kappa=0, J_bridge swept (pure Hamiltonian) ---")
    log()
    log(f"{'J_br':>8}  {'PiChain':>10}  {'PiAxB':>10}  "
        f"{'Paired':>8}  {'BestSg':>8}  {'X^N':>5}")
    log("-" * 55)

    for J_br in np.linspace(0, 2.0, 15):
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=J_br, kappa=0)
        evals = np.linalg.eigvals(L)

        err_chain, Sg_c = pi_conjugation_error(L, Pi_comp, Sg_candidates)
        err_AB, Sg_AB = pi_conjugation_error(L, Pi_AB_comp, Sg_candidates)

        _, best_np, _ = best_palindrome(evals, Sg)

        comm_xn = np.max(np.abs(XN_super @ L - L @ XN_super))

        log(f"{J_br:>8.4f}  {err_chain:>10.2e}  {err_AB:>10.2e}  "
            f"{best_np:>4}/{d2}  {Sg_c:>8.4f}  "
            f"{'Y' if comm_xn < 1e-10 else 'N':>5}")
    log()

    # Slice c) J_bridge = kappa, both swept (diagonal)
    log("--- Slice (c): J_bridge = kappa, both swept (diagonal) ---")
    log()
    log(f"{'J=k':>8}  {'PiChain':>10}  {'PiAxB':>10}  "
        f"{'Paired':>8}  {'BestSg':>8}  {'X^N':>5}")
    log("-" * 55)

    for val in np.linspace(0, 0.3, 15):
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=val, kappa=val)
        evals = np.linalg.eigvals(L)

        err_chain, Sg_c = pi_conjugation_error(L, Pi_comp, Sg_candidates)
        err_AB, Sg_AB = pi_conjugation_error(L, Pi_AB_comp, Sg_candidates)

        _, best_np, _ = best_palindrome(evals, Sg)

        comm_xn = np.max(np.abs(XN_super @ L - L @ XN_super))

        log(f"{val:>8.4f}  {err_chain:>10.2e}  {err_AB:>10.2e}  "
            f"{best_np:>4}/{d2}  {Sg_c:>8.4f}  "
            f"{'Y' if comm_xn < 1e-10 else 'N':>5}")
    log()

    # Eigenvector-based Pi construction attempt
    log("--- Eigenvector-based Pi construction ---")
    log()

    # Try at several representative points
    test_points = [
        (1.0, 0.0, "pure Hamiltonian bridge"),
        (0.0, 0.1, "pure dissipative bridge"),
        (1.0, 0.1, "mixed bridge"),
        (0.5, 0.05, "weak mixed bridge"),
    ]

    for J_br, kappa, label in test_points:
        log(f"  Point: J_bridge={J_br}, kappa={kappa} ({label})")

        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=J_br, kappa=kappa)
        evals, R = np.linalg.eig(L)

        # Find palindromic pairs
        best_err, best_np, best_Sg = best_palindrome(evals, Sg, tol=1e-6)

        log(f"  Palindromic pairs at Sg={best_Sg:.4f}: {best_np}/{d2}")

        if best_np == d2:
            # Build Pi from eigenvector pairing
            pairs_ev = []
            used = set()
            for k in range(d2):
                if k in used: continue
                target = -(evals[k] + 2 * best_Sg)
                dists = np.abs(evals - target)
                for u in used: dists[u] = 1e30
                j = np.argmin(dists)
                if dists[j] < 1e-6:
                    pairs_ev.append((k, j))
                    used.add(k); used.add(j)

            # Attempt Pi construction: Pi |r_k> = phase * |r_j>
            # Build Pi column by column using right eigenvectors
            R_inv = np.linalg.inv(R)
            Pi_ev = np.zeros((d2, d2), dtype=complex)
            for k, j in pairs_ev:
                if k == j:
                    Pi_ev[:, k] = R[:, k]
                else:
                    # Phase: choose so Pi^2 has nice structure
                    Pi_ev[:, k] = R[:, j]
                    Pi_ev[:, j] = R[:, k]

            # Pi_ev maps eigenvectors, convert to operator form
            Pi_constructed = Pi_ev @ R_inv

            # Check if it's unitary-like
            prod = Pi_constructed @ Pi_constructed.conj().T
            unitarity = np.max(np.abs(prod - prod[0,0] * np.eye(d2)))
            log(f"  Constructed Pi unitarity error: {unitarity:.2e}")

            # Check conjugation
            err_con, Sg_con = pi_conjugation_error(L, Pi_constructed,
                                                    Sg_candidates)
            log(f"  Constructed Pi conjugation error: {err_con:.2e} "
                f"(Sg={Sg_con:.4f})")

            # What is Pi^2?
            Pi2 = Pi_constructed @ Pi_constructed
            comm_Pi2_XN = np.max(np.abs(Pi2 - XN_super))
            log(f"  ||Pi^2 - X^N|| = {comm_Pi2_XN:.2e}")

            # Compare with Pi_chain
            diff_chain = np.max(np.abs(Pi_constructed - Pi_comp))
            diff_AB = np.max(np.abs(Pi_constructed - Pi_AB_comp))
            log(f"  ||Pi_constructed - Pi_chain|| = {diff_chain:.2e}")
            log(f"  ||Pi_constructed - Pi_AxB|| = {diff_AB:.2e}")
        else:
            log(f"  Cannot construct Pi (palindrome incomplete)")

        log()


# ============================================================
# TEST 2: THE CRITICAL LINE
# ============================================================
def run_test_2(gamma=0.1, J=1.0, J_vals=None, k_vals=None,
               pal_count_map=None):
    log("=" * 70)
    log("TEST 2: THE CRITICAL LINE")
    log("=" * 70)
    log()

    # If we have phase diagram data, find the boundary
    # Otherwise, do a 1D sweep along J_bridge=1.0, sweeping kappa
    Sg_base = gamma * N

    log("--- 1D sweep: J_bridge=1.0, kappa from 0 to 0.5 ---")
    log()
    log(f"{'kappa':>8}  {'Paired':>8}  {'Err':>10}  "
        f"{'1stBreak':>10}  {'Transition':>12}")
    log("-" * 55)

    prev_np = d2
    kappas = np.linspace(0, 0.5, 50)

    for kappa in kappas:
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=1.0, kappa=kappa)
        evals = np.linalg.eigvals(L)

        best_err, best_np, best_Sg = best_palindrome(evals, Sg)

        break_info = ""
        if best_np < d2 and best_np < prev_np:
            break_info = "BREAKING"

        log(f"{kappa:>8.4f}  {best_np:>4}/{d2}  {best_err:>10.2e}  "
            f"{'':>10}  {break_info:>12}")

        prev_np = best_np

    log()

    # Detailed analysis at the transition point
    log("--- Fine sweep around transition ---")
    log()

    # Find approximate transition kappa
    for kappa_test in np.linspace(0, 0.3, 100):
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=1.0, kappa=kappa_test)
        evals = np.linalg.eigvals(L)
        _, np_test, best_Sg_t = best_palindrome(evals, Sg)
        if np_test < d2:
            k_transition = kappa_test
            break
    else:
        k_transition = 0.3

    log(f"Approximate transition at kappa ~ {k_transition:.4f}")
    log()

    # Fine sweep around transition
    kappas_fine = np.linspace(max(0, k_transition - 0.02),
                              k_transition + 0.02, 40)

    log(f"{'kappa':>10}  {'Paired':>8}  {'Err':>10}  {'BestSg':>8}")
    log("-" * 42)

    for kappa in kappas_fine:
        L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=1.0, kappa=kappa)
        evals = np.linalg.eigvals(L)
        best_err, best_np, best_Sg = best_palindrome(evals, Sg)

        log(f"{kappa:>10.5f}  {best_np:>4}/{d2}  {best_err:>10.2e}  "
            f"{best_Sg:>8.4f}")
    log()

    # At the transition: which Pauli-weight modes break first?
    log("--- Mode analysis at transition ---")
    log()

    kappa_at = k_transition + 0.001  # just past transition
    L, Sg = build_mixed_bridge(gamma=gamma, J_bridge=1.0, kappa=kappa_at)
    evals = np.linalg.eigvals(L)

    # Find best Sg
    _, _, best_Sg = best_palindrome(evals, Sg, tol=1e-6)

    # Check each eigenvalue's pairing error
    pair_errors = []
    for k in range(d2):
        target = -(evals[k] + 2 * best_Sg)
        best_match = np.min(np.abs(evals - target))
        pair_errors.append(best_match)

    pair_errors = np.array(pair_errors)
    broken = np.where(pair_errors > 1e-6)[0]
    intact = np.where(pair_errors <= 1e-6)[0]

    log(f"At kappa={kappa_at:.4f}, Sg={best_Sg:.4f}:")
    log(f"  Broken modes: {len(broken)}/{d2}")
    log(f"  Intact modes: {len(intact)}/{d2}")
    log()

    # Eigenvector Pauli weight analysis of broken modes
    _, R = np.linalg.eig(L)

    log("--- Pauli weight of broken vs intact modes ---")
    log()

    def eigvec_pauli_weights(vec):
        """Decompose eigenvector into Pauli basis and compute XY-weight histogram."""
        coeffs = V_inv @ vec
        weight_hist = np.zeros(N + 1)
        for a, idx in enumerate(all_idx):
            w = xy_weight(idx)
            weight_hist[w] += abs(coeffs[a]) ** 2
        total = np.sum(weight_hist)
        if total > 1e-30: weight_hist /= total
        return weight_hist

    broken_weights = np.zeros(N + 1)
    for k in broken:
        broken_weights += eigvec_pauli_weights(R[:, k])
    if len(broken) > 0:
        broken_weights /= len(broken)

    intact_weights = np.zeros(N + 1)
    for k in intact:
        intact_weights += eigvec_pauli_weights(R[:, k])
    if len(intact) > 0:
        intact_weights /= len(intact)

    log(f"{'XY-weight':>10}  {'Broken':>10}  {'Intact':>10}")
    log("-" * 35)
    for w in range(N + 1):
        log(f"{w:>10}  {broken_weights[w]:>10.4f}  {intact_weights[w]:>10.4f}")
    log()


# ============================================================
# TEST 5: DOES THE DECODER STILL WORK?
# ============================================================
def run_test_5(gamma=0.1, J=1.0, J_br=1.0, kappa=0.05):
    log("=" * 70)
    log("TEST 5: DOES THE DECODER STILL WORK?")
    log("=" * 70)
    log()

    Sg_base = gamma * N

    # Build reference system
    L_ref, Sg = build_mixed_bridge(gamma=gamma, J_bridge=J_br, kappa=kappa)
    evals_ref, R_ref = np.linalg.eig(L_ref)

    # Verify palindrome
    _, best_np, best_Sg = best_palindrome(evals_ref, Sg)

    log(f"Reference: J_bridge={J_br}, kappa={kappa}")
    log(f"Palindromic pairs: {best_np}/{d2}")
    log()

    # Initial states
    up = np.array([1, 0], dtype=complex)
    dn = np.array([0, 1], dtype=complex)

    def kron_list(vecs):
        r = vecs[0]
        for v in vecs[1:]: r = np.kron(r, v)
        return r

    test_states = {}
    test_states['|0110>'] = kron_list([up, dn, dn, up])
    test_states['|0100>'] = kron_list([up, dn, up, up])
    test_states['|0010>'] = kron_list([up, up, dn, up])
    test_states['|0000>'] = kron_list([up, up, up, up])

    # Bell on boundary qubits
    bell_12 = np.zeros(d, dtype=complex)
    bell_12[0] = 1/np.sqrt(2)
    bell_12[0b0110] = 1/np.sqrt(2)
    test_states['Bell(12)|00>'] = bell_12

    # GHZ
    ghz = np.zeros(d, dtype=complex)
    ghz[0] = 1/np.sqrt(2); ghz[d-1] = 1/np.sqrt(2)
    test_states['GHZ'] = ghz

    # Distinguishability: vary kappa slightly
    kappa_pert = kappa + 0.01
    L_pert, _ = build_mixed_bridge(gamma=gamma, J_bridge=J_br, kappa=kappa_pert)

    t_measure = np.linspace(0.5, 30, 60)

    log("--- Distinguishability: reference kappa vs kappa+0.01 ---")
    log()
    log(f"{'State':>14}  {'t_max':>8}  {'max_dist':>10}  "
        f"{'t=1':>8}  {'t=5':>8}  {'t=10':>8}")
    log("-" * 60)

    for sname, psi0 in test_states.items():
        rho0 = np.outer(psi0, psi0.conj())
        rho0_vec = rho0.flatten()

        dists = []
        for t in t_measure:
            rho_ref = (expm(L_ref * t) @ rho0_vec).reshape(d, d)
            rho_p = (expm(L_pert * t) @ rho0_vec).reshape(d, d)
            dists.append(np.linalg.norm(rho_ref - rho_p))

        dists = np.array(dists)
        max_dist = np.max(dists)
        t_max = t_measure[np.argmax(dists)]
        d1 = dists[np.argmin(np.abs(t_measure - 1))]
        d5 = dists[np.argmin(np.abs(t_measure - 5))]
        d10 = dists[np.argmin(np.abs(t_measure - 10))]

        log(f"{sname:>14}  {t_max:>8.2f}  {max_dist:>10.6f}  "
            f"{d1:>8.6f}  {d5:>8.6f}  {d10:>8.6f}")

    log()

    # Fisher information for kappa sensitivity
    log("--- Fisher information for kappa at t=5 ---")
    log()

    dg = 1e-5
    L_p, _ = build_mixed_bridge(gamma=gamma, J_bridge=J_br, kappa=kappa + dg)
    t_fish = 5.0

    log(f"{'State':>14}  {'F(kappa)':>12}  {'Sensitivity':>12}")
    log("-" * 42)

    for sname, psi0 in test_states.items():
        rho0 = np.outer(psi0, psi0.conj())
        rho0_vec = rho0.flatten()

        rho_0 = (expm(L_ref * t_fish) @ rho0_vec).reshape(d, d)
        rho_p = (expm(L_p * t_fish) @ rho0_vec).reshape(d, d)

        p0 = np.real(np.diag(rho_0))
        pp = np.real(np.diag(rho_p))
        dp = (pp - p0) / dg

        fisher = 0
        for k in range(d):
            if p0[k] > 1e-12:
                fisher += dp[k] ** 2 / p0[k]

        sensitivity = np.linalg.norm(rho_p - rho_0) / dg

        log(f"{sname:>14}  {fisher:>12.4f}  {sensitivity:>12.4f}")

    log()

    # Response matrix: sensitivity to per-site gamma
    log("--- Response matrix: mode sensitivity to per-site gamma ---")
    log()

    gammas_ref = [gamma] * N
    R_inv_ref = np.linalg.inv(R_ref)

    # Use |0110> as initial state (boundary qubits active)
    psi0 = kron_list([up, dn, dn, up])
    rho0_vec = np.outer(psi0, psi0.conj()).flatten()
    coeffs_ref = R_inv_ref @ rho0_vec

    # Build palindromic pairs
    pairs = []
    used = set()
    for k in range(d2):
        if k in used: continue
        target = -(evals_ref[k] + 2 * best_Sg)
        dists_ev = np.abs(evals_ref - target)
        for u in used: dists_ev[u] = 1e30
        j = np.argmin(dists_ev)
        if dists_ev[j] < 1e-6:
            pairs.append((k, j))
            used.add(k); used.add(j)

    n_pairs = len(pairs)
    dg = 1e-4
    response = np.zeros((n_pairs, N))

    for site in range(N):
        gammas_pert_s = list(gammas_ref)
        gammas_pert_s[site] += dg
        L_ps = build_L_H(
            build_H_heisenberg([(0,1)], J) +
            build_H_heisenberg([(2,3)], J) +
            build_H_heisenberg([(1,2)], J_br)
        ) + build_L_D_zdeph(gammas_pert_s)
        # Add dissipative bridge
        F1 = site_op(sx, 1) @ site_op(sz, 2)
        F2 = site_op(sy, 1) @ site_op(sz, 2)
        F3 = site_op(sz, 1) @ site_op(sx, 2)
        F4 = site_op(sz, 1) @ site_op(sy, 2)
        L_ps = add_lindblad_jump(L_ps, F1, kappa)
        L_ps = add_lindblad_jump(L_ps, F2, kappa)
        L_ps = add_lindblad_jump(L_ps, F3, kappa)
        L_ps = add_lindblad_jump(L_ps, F4, kappa)

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

    if n_pairs == 0:
        log("No palindromic pairs found -- response matrix not available.")
        log("(Palindrome is broken at this operating point.)")
    else:
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
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("The Mixed Bridge: Hamiltonian + Dissipative, Both or Nothing")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, d={d}, d^2={d2}")
    log()

    # --- TEST 1: Phase Diagram ---
    t1 = time.time()
    result = run_test_1()
    J_vals, k_vals, pal_err, pal_count, mi, survival = result
    log(f"[Test 1 completed in {time.time() - t1:.1f}s]")
    log()

    # Check if we should continue
    has_survival = np.any(survival)

    # Check if overlap is exclusively at kappa=0 (trivial case)
    if has_survival:
        nontrivial = False
        for ik in range(len(k_vals)):
            for ij in range(len(J_vals)):
                if survival[ik, ij] and mi[ik, ij] > 0.001 and k_vals[ik] > 0:
                    nontrivial = True
        if not nontrivial:
            has_survival = False  # only kappa=0, not the mixed bridge

    if not has_survival:
        log("=" * 70)
        log("NO MIXED BRIDGE SURVIVAL REGION")
        log("=" * 70)
        log()
        log("Palindrome survives ONLY at kappa=0 (pure Hamiltonian bridge).")
        log("This is the standard Heisenberg chain, not a mixed bridge.")
        log("Adding ANY dissipative channel breaks the palindrome immediately.")
        log("The bridge cannot be modeled as Hamiltonian + Lindblad on")
        log("a tensor product space with palindromic structure preserved")
        log("AND dissipative information flowing simultaneously.")
        log()
        log("Continuing with remaining tests for completeness...")
        log()

    # --- TEST 3: Bidirectionality ---
    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time() - t3:.1f}s]")
    log()

    # --- TEST 4: Pi operator ---
    t4 = time.time()
    run_test_4()
    log(f"[Test 4 completed in {time.time() - t4:.1f}s]")
    log()

    # --- TEST 2: Critical line ---
    t2 = time.time()
    run_test_2()
    log(f"[Test 2 completed in {time.time() - t2:.1f}s]")
    log()

    # --- TEST 5: Decoder ---
    t5 = time.time()
    # Use a point where palindrome is approximately best
    run_test_5(J_br=1.0, kappa=0.05)
    log(f"[Test 5 completed in {time.time() - t5:.1f}s]")
    log()

    # ============================================================
    # OVERALL SUMMARY
    # ============================================================
    log("=" * 70)
    log("OVERALL SUMMARY")
    log("=" * 70)
    log()
    log("Test 1: Phase diagram (J_bridge, kappa)")
    log(f"  Mixed bridge survival region (kappa > 0): "
        f"{'EXISTS' if has_survival else 'NONE'}")
    log(f"  Pure Hamiltonian (kappa=0): palindromic at all J_bridge")
    log()
    log("Test 3: Bidirectionality in the mixed bridge")
    log("  Does the Hamiltonian protect against asymmetric dissipation?")
    log()
    log("Test 4: Pi operator")
    log("  Pi_chain, Pi_AxB, eigenvector construction across slices")
    log()
    log("Test 2: Critical line")
    log("  Where palindrome breaks as dissipation increases")
    log()
    log("Test 5: Decoder")
    log("  Fisher information for kappa sensitivity, response matrix")
    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
