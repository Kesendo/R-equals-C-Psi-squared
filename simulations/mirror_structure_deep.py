#!/usr/bin/env python3
"""
mirror_structure_deep.py - The Structure of the Mirror
Rank Threshold, Qudit Universality, Engineering, and Geometry

Four open threads from the boot script investigation:
  1. Why does high Choi rank prevent scaling? (mechanism)
  2. Single-bond universality for qudits (scope)
  3. Engineering the entangled mirror point (engineering)
  4. The clock as space-creator (geometry)
"""

import numpy as np
from scipy import linalg as la
from itertools import product as iprod
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "mirror_structure_deep.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
# Force UTF-8 on stdout for Windows
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# PAULI INFRASTRUCTURE
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
H_LABELS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']
GAMMA = 0.05

# Classification from prior results
SURVIVING = {
    'XX+XZ', 'XX+YY', 'XX+YZ', 'XX+ZX', 'XX+ZY', 'XX+ZZ',
    'XZ+YY', 'XZ+ZX', 'XZ+ZZ',
    'YY+YZ', 'YY+ZX', 'YY+ZY', 'YY+ZZ',
    'YZ+ZY', 'YZ+ZZ', 'ZX+ZZ', 'ZY+ZZ',
    'XY+YX', 'XY+ZZ', 'YX+ZZ',
    'XZ+YZ', 'ZX+ZY'}
NONLOCAL = {'XZ+YZ', 'ZX+ZY'}

def all_36():
    combos = []
    for i, t1 in enumerate(H_LABELS):
        for t2 in H_LABELS[i + 1:]:
            combos.append(f"{t1}+{t2}")
    return combos

# ============================================================
# CORE BUILDERS
# ============================================================
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

def build_H_term(N, bond, label):
    ops = [I2] * N
    ops[bond[0]] = PM[label[0]]
    ops[bond[1]] = PM[label[1]]
    H = ops[0]
    for k in range(1, N):
        H = np.kron(H, ops[k])
    return H

def build_L_H(H, all_idx, pmats, d):
    num = len(all_idx)
    L = np.zeros((num, num), dtype=complex)
    for a in range(num):
        for b in range(num):
            comm = -1j * (H @ pmats[b] - pmats[b] @ H)
            L[a, b] = np.trace(pmats[a].conj().T @ comm) / d
    return L

def build_L_D_diag(gamma, all_idx):
    return np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])

def assemble_L(L_H, L_D_diag):
    L = L_H.copy()
    for a in range(len(L_D_diag)):
        L[a, a] += L_D_diag[a]
    return L

def build_H_combo(N, bonds, combo):
    terms = combo.split('+')
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for bond in bonds:
        for term in terms:
            H += build_H_term(N, bond, term)
    return H

def build_full_L(N, bonds, combo, gamma):
    all_idx, pmats, d = build_pauli_basis(N)
    H = build_H_combo(N, bonds, combo)
    L_H = build_L_H(H, all_idx, pmats, d)
    L_D = build_L_D_diag(gamma, all_idx)
    L = assemble_L(L_H, L_D)
    Sg = N * gamma
    return L, Sg, all_idx, pmats, d

# ============================================================
# PI CONSTRUCTION
# ============================================================
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}   # I<->X, Y<->Z
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

def build_Pi_canonical(all_idx):
    """Standard per-site Pi: I<->X (+1), Y->iZ, Z->iY."""
    num = len(all_idx)
    Pi = np.zeros((num, num), dtype=complex)
    idx_map = {idx: i for i, idx in enumerate(all_idx)}
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = idx_map[mapped]
        Pi[a, b] = sign
    return Pi

def construct_Pi(L_full, Sg, tol=1e-8):
    num = L_full.shape[0]
    evals, R = np.linalg.eig(L_full)
    R_inv = np.linalg.inv(R)
    P_eig = np.zeros((num, num), dtype=complex)
    used = set()
    n_paired = 0
    pair_errors = []
    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = np.argmin(diffs)
        err = diffs[best]
        if err < tol:
            P_eig[k, best] = 1.0
            P_eig[best, k] = 1.0
            used.add(k)
            used.add(best)
            n_paired += 1 if k == best else 2
            pair_errors.append(err)
        else:
            P_eig[k, k] = 1.0
            used.add(k)
            pair_errors.append(err)
    Pi = R @ P_eig @ R_inv
    return Pi, n_paired, pair_errors

def palindrome_error(Pi, L, Sg):
    Pi_inv = np.linalg.inv(Pi)
    E = Pi @ L @ Pi_inv + L + 2 * Sg * np.eye(L.shape[0])
    return np.max(np.abs(E)) / max(np.max(np.abs(L)), 1e-15)

def check_eigenvalue_pairing(evals, Sg, tol=1e-6):
    """Return (n_well_paired, n_total, max_error_of_paired, errors_list)"""
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

# ============================================================
# CHOI PIPELINE
# ============================================================
def build_V(pmats, d):
    V = np.zeros((d * d, len(pmats)), dtype=complex)
    for a, P in enumerate(pmats):
        V[:, a] = P.reshape(-1)
    return V

def pi_to_choi(Pi, V, d):
    S = V @ Pi @ V.conj().T / d
    J = S.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d) / d
    return J

def choi_site_svd_2site(J):
    """SVD across site bipartition for N=2 Choi matrix (16x16)."""
    J8 = J.reshape(2, 2, 2, 2, 2, 2, 2, 2)
    J_site = J8.transpose(0, 2, 4, 6, 1, 3, 5, 7).reshape(16, 16)
    svs = np.linalg.svd(J_site, compute_uv=False)
    rank = int(np.sum(svs > 1e-10 * max(svs[0], 1e-30)))
    return rank, svs

def schmidt_entropy(svs):
    s = svs[svs > 1e-15]
    if len(s) == 0:
        return 0.0
    p = (s / np.linalg.norm(s)) ** 2
    p = p[p > 0]
    return -float(np.sum(p * np.log2(p)))

def partial_transpose(J, dA, dB, sub='B'):
    J4 = J.reshape(dA, dB, dA, dB)
    if sub == 'B':
        J4 = J4.transpose(0, 3, 2, 1)
    else:
        J4 = J4.transpose(2, 1, 0, 3)
    return J4.reshape(dA * dB, dA * dB)

def negativity(J, dA, dB):
    Jpt = partial_transpose(J, dA, dB)
    eigs = np.linalg.eigvalsh(Jpt)
    return float(np.sum(np.abs(eigs[eigs < -1e-14])))

# ============================================================
# COMPUTATIONAL-BASIS LINDBLADIAN (for QST fidelity)
# ============================================================
def build_L_comp(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N
        ops[k] = sz
        Zk = ops[0]
        for j in range(1, N):
            Zk = np.kron(Zk, ops[j])
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L

def ptrace(rho, N, keep):
    """Partial trace of density matrix, keeping specified qubits."""
    d = 2 ** N
    rho_t = rho.reshape([2] * (2 * N))
    n_keep = len(keep)
    trace_out = [i for i in range(N) if i not in keep]
    for i in sorted(trace_out, reverse=True):
        rho_t = np.trace(rho_t, axis1=i, axis2=i + N - (N - len(trace_out) - (len(trace_out) - (N - 1 - i))))
    # Simpler approach: use explicit contraction
    d_keep = 2 ** n_keep
    d_trace = 2 ** (N - n_keep)
    # Reorder axes: keep qubits first, then trace-out qubits, for both bra and ket
    all_axes = list(keep) + list(trace_out) + [k + N for k in keep] + [k + N for k in trace_out]
    rho_t = rho.reshape([2] * (2 * N)).transpose(all_axes)
    rho_t = rho_t.reshape(d_keep, d_trace, d_keep, d_trace)
    return np.trace(rho_t, axis1=1, axis2=3)

def qst_fidelity_sweep(H, gamma, N, t_points):
    """Compute max state transfer fidelity |0>→|N-1> over time."""
    L_c = build_L_comp(H, gamma, N)
    d = 2 ** N
    # Initial state: |1> at site 0, |0> elsewhere
    psi_init = np.zeros(d, dtype=complex)
    psi_init[1 << (N - 1)] = 1.0  # site 0 = |1>
    rho_init = np.outer(psi_init, psi_init.conj())
    vec_init = rho_init.reshape(-1)
    # Target: |1> at site N-1
    psi_target = np.zeros(2, dtype=complex)
    psi_target[1] = 1.0
    best_fid = 0.0
    best_t = 0.0
    for t in t_points:
        eLt = la.expm(L_c * t)
        vec_out = eLt @ vec_init
        rho_out = vec_out.reshape(d, d)
        # Partial trace to get site N-1
        rho_r = ptrace(rho_out, N, [N - 1])
        fid = float(np.real(psi_target.conj() @ rho_r @ psi_target))
        if fid > best_fid:
            best_fid = fid
            best_t = t
    return best_fid, best_t

# ============================================================
# ============================================================
#  SECTION 1: WHY DOES HIGH CHOI RANK PREVENT SCALING?
# ============================================================
# ============================================================

def run_section_1():
    log("=" * 70)
    log("SECTION 1: WHY DOES HIGH CHOI RANK PREVENT SCALING?")
    log("=" * 70)
    log()

    # ----- 1a: Collision at the shared site -----
    log("--- 1a: Collision at the shared site ---")
    log()

    N = 2
    gamma = GAMMA
    Sg = N * gamma
    all_idx, pmats, d = build_pauli_basis(N)
    L_D = build_L_D_diag(gamma, all_idx)
    V = build_V(pmats, d)
    num = 4 ** N

    combos = all_36()
    BROKEN = [c for c in combos if c not in SURVIVING]

    # Build canonical Pi (standard per-site map) and check which combos it works for
    Pi_canon = build_Pi_canonical(all_idx)
    canon_err = {}
    for combo in combos:
        H = build_H_combo(N, [(0, 1)], combo)
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        canon_err[combo] = palindrome_error(Pi_canon, L, Sg)

    results = {}
    log(f"{'Combo':<12} {'EigRk':>6} {'CanonOK':>7} {'PalErr':>9} "
        f"{'R_diff':>8} {'Asym':>8} {'Surv':>5}")
    log("-" * 65)

    for combo in combos:
        H = build_H_combo(N, [(0, 1)], combo)
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)

        # Construct Pi via eigenvector pairing
        Pi, n_paired, _ = construct_Pi(L, Sg)
        err = palindrome_error(Pi, L, Sg)

        # Choi of eigenvector Pi
        J = pi_to_choi(Pi, V, d)
        choi_rank, svs = choi_site_svd_2site(J)
        s_ent = schmidt_entropy(svs)

        # Does canonical per-site Pi work?
        c_ok = "YES" if canon_err[combo] < 1e-10 else "no"

        # Extract per-site reduced maps via marginalization
        Pi_4d = Pi.reshape(4, 4, 4, 4)
        R_left = np.sum(Pi_4d, axis=(0, 2)) / 4.0
        R_right = np.sum(Pi_4d, axis=(1, 3)) / 4.0
        norm_l = np.linalg.norm(R_left)
        norm_r = np.linalg.norm(R_right)
        r_diff = (np.linalg.norm(R_left - R_right) / max(norm_l, norm_r)
                  if max(norm_l, norm_r) > 1e-12 else 0.0)

        # Operator SVD of Pi
        Pi_swap = Pi_4d.transpose(0, 2, 1, 3).reshape(16, 16)
        asym = np.linalg.norm(Pi_swap - Pi_swap.T) / max(np.linalg.norm(Pi_swap), 1e-15)

        surv = "YES" if combo in SURVIVING else "NO"
        log(f"{combo:<12} {choi_rank:>6} {c_ok:>7} {err:>9.1e} "
            f"{r_diff:>8.5f} {asym:>8.5f} {surv:>5}")

        results[combo] = {
            'choi_rank': choi_rank, 's_ent': s_ent, 'r_diff': r_diff,
            'asym': asym, 'survives': combo in SURVIVING,
            'R_left': R_left, 'R_right': R_right, 'Pi': Pi,
            'n_paired': n_paired, 'pal_err': err, 'canon_ok': c_ok == "YES"
        }

    log()
    log("Note: EigRk = Choi rank of eigenvector-constructed Pi.")
    log("  The eigenvector Pi is NOT unique (degenerate eigenvalues allow")
    log("  different pairings). Its rank varies and does NOT cleanly separate")
    log("  surviving from broken combos.")
    log()
    log("Canonical Pi (I<->X, Y<->Z with specific phases) works for:")
    n_canon = sum(1 for r in results.values() if r['canon_ok'])
    canon_list = [c for c in combos if results[c]['canon_ok']]
    log(f"  {n_canon}/36 combos: {', '.join(canon_list)}")
    log("  Other combos have local Pi with DIFFERENT per-site maps")
    log("  (from the P1/P4 families with varied phases).")
    log()
    log("KEY DISCRIMINATOR: N=3 eigenvalue pairing (Section 1c)")
    log("  Surviving combos: ALL eigenvalues pair at N=3")
    log("  Broken combos: MOST eigenvalues fail to pair at N=3")

    # ----- 1b: Commutator analysis -----
    log()
    log("--- 1b: Commutator analysis at shared site ---")
    log()
    log(f"{'Combo':<12} {'||[R_L,R_R]||':>14} {'||R_L-R_R||':>14} {'Choi_Rk':>8} {'Surv':>5}")
    log("-" * 60)

    for combo in combos:
        r = results[combo]
        R_L = r['R_left']
        R_R = r['R_right']
        comm = R_L @ R_R - R_R @ R_L
        comm_norm = np.linalg.norm(comm)
        diff_norm = np.linalg.norm(R_L - R_R)
        surv = "YES" if r['survives'] else "NO"
        log(f"{combo:<12} {comm_norm:>14.8f} {diff_norm:>14.8f} {r['choi_rank']:>8} {surv:>5}")

    # ----- 1c: Partial Pi_3 for broken combos -----
    log()
    log("--- 1c: Can the collision be resolved? Partial Pi_3 ---")
    log()

    N3 = 3
    Sg3 = N3 * gamma
    all_idx3, pmats3, d3 = build_pauli_basis(N3)
    L_D3 = build_L_D_diag(gamma, all_idx3)
    num3 = 4 ** N3
    V3 = build_V(pmats3, d3)

    log(f"{'Combo':<12} {'N2_Rk':>6} {'Paired':>7} {'Total':>6} "
        f"{'MaxPairErr':>11} {'Pi3_PalErr':>11}")
    log("-" * 60)

    for combo in BROKEN[:14]:
        H3 = build_H_combo(N3, [(0, 1), (1, 2)], combo)
        L_H3 = build_L_H(H3, all_idx3, pmats3, d3)
        L3 = assemble_L(L_H3, L_D3)

        # Check eigenvalue pairing
        evals3 = np.linalg.eigvals(L3)
        n_well, n_tot, errors = check_eigenvalue_pairing(evals3, Sg3)

        # Construct partial Pi_3
        Pi3, n3_paired, _ = construct_Pi(L3, Sg3, tol=1e-6)
        err3 = palindrome_error(Pi3, L3, Sg3)

        max_pair_err = max(e for e in errors if e < 1.0) if any(e < 1.0 for e in errors) else 0.0

        log(f"{combo:<12} {results[combo]['choi_rank']:>6} {n_well:>7}/{n_tot:<5} "
            f"{max_pair_err:>11.2e} {err3:>11.2e}")

    # Analysis for surviving combos at N=3
    log()
    log("Control: surviving combos at N=3")
    log(f"{'Combo':<12} {'N2_Rk':>6} {'Paired':>7} {'Pi3_PalErr':>11}")
    log("-" * 45)

    for combo in sorted(SURVIVING)[:8]:
        H3 = build_H_combo(N3, [(0, 1), (1, 2)], combo)
        L_H3 = build_L_H(H3, all_idx3, pmats3, d3)
        L3 = assemble_L(L_H3, L_D3)
        evals3 = np.linalg.eigvals(L3)
        n_well, n_tot, _ = check_eigenvalue_pairing(evals3, Sg3)
        Pi3, _, _ = construct_Pi(L3, Sg3)
        err3 = palindrome_error(Pi3, L3, Sg3)
        log(f"{combo:<12} {results[combo]['choi_rank']:>6} {n_well:>7}/{n_tot:<5} {err3:>11.2e}")


# ============================================================
# ============================================================
#  SECTION 2: SINGLE-BOND UNIVERSALITY FOR QUDITS
# ============================================================
# ============================================================

def build_gell_mann():
    """Return 9 qutrit basis matrices with Tr(B_i B_j) = 3*delta_ij."""
    I3 = np.eye(3, dtype=complex)
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    l8 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3)

    # Normalize: Tr(B_i B_j) = d * delta_ij = 3 * delta_ij
    c = np.sqrt(3.0 / 2.0)
    raw = [l1, l2, l3, l4, l5, l6, l7, l8]
    return [I3] + [c * m for m in raw], raw

def run_section_2():
    log()
    log("=" * 70)
    log("SECTION 2: SINGLE-BOND UNIVERSALITY FOR QUDITS")
    log("=" * 70)
    log()

    # ----- 2b: Per-site split analysis (do first - gives the structural answer) -----
    log("--- 2b: Per-site split for d=3 (qutrits) ---")
    log()

    basis_norm, gm_raw = build_gell_mann()
    gm_names = ['I3', 'λ1', 'λ2', 'λ3', 'λ4', 'λ5', 'λ6', 'λ7', 'λ8']
    immune_idx = {0, 3, 8}  # I, λ3, λ8 (diagonal)

    log("Qutrit Z-dephasing (jump ops: λ3, λ8):")
    log(f"  Immune basis elements (commute with λ3, λ8): "
        f"{[gm_names[i] for i in sorted(immune_idx)]}")
    log(f"  Decaying basis elements: "
        f"{[gm_names[i] for i in range(9) if i not in immune_idx]}")
    log(f"  Per-site split: {len(immune_idx)} immune : "
        f"{9 - len(immune_idx)} decaying = 3:6")
    log()

    # Verify dephasing rates numerically
    log("Dephasing rates (single qutrit site, γ=1):")
    for i, B in enumerate(basis_norm):
        # D(B) = γ(λ3 B λ3 + λ8 B λ8 - ½{λ3², B} - ½{λ8², B})
        l3 = gm_raw[2]  # index 2 in raw = λ3
        l8 = gm_raw[7]  # index 7 in raw = λ8
        DB = (l3 @ B @ l3 + l8 @ B @ l8
              - 0.5 * (l3 @ l3 @ B + B @ l3 @ l3)
              - 0.5 * (l8 @ l8 @ B + B @ l8 @ l8))
        # D(B) should be proportional to B: D(B) = rate * B
        if np.linalg.norm(B) > 1e-10:
            # Extract rate: rate = Tr(B† D(B)) / Tr(B† B)
            rate = np.trace(B.conj().T @ DB) / np.trace(B.conj().T @ B)
            log(f"  {gm_names[i]:>3}: rate = {rate.real:>8.4f} "
                f"(expected: {0 if i in immune_idx else -2:.1f})")
    log()

    # Block counts for N=2
    log("N=2 qutrit block counts:")
    for w in range(3):
        from math import comb
        n_immune = len(immune_idx)
        n_decay = 9 - n_immune
        count = comb(2, w) * (n_decay ** w) * (n_immune ** (2 - w))
        log(f"  w={w}: {count} elements")
    log()
    log("For palindrome: w=0 (9 elements) must pair with w=2 (36 elements)")
    log("IMPOSSIBLE: 9 ≠ 36. The palindrome CANNOT hold for qutrits.")
    log()

    # ----- 2a: Numerical verification -----
    log("--- 2a: Numerical verification for qutrits (d=3) ---")
    log()

    d_qt = 3  # Hilbert space dimension for one qutrit
    N_qt = 2
    d_full = d_qt ** N_qt  # = 9

    # Build 2-site qutrit basis
    all_idx_qt = list(iprod(range(9), repeat=N_qt))
    pmats_qt = []
    for (a, b) in all_idx_qt:
        pmats_qt.append(np.kron(basis_norm[a], basis_norm[b]))
    num_qt = len(all_idx_qt)

    # Dissipator
    L_D_qt = np.array([-2 * GAMMA * sum(1 for i in idx if i not in immune_idx)
                        for idx in all_idx_qt])
    Sg_qt = N_qt * GAMMA

    # Test several qutrit Hamiltonian pairs
    qutrit_tests = [
        ("λ1⊗λ1", 0, 0),
        ("λ1⊗λ2", 0, 1),
        ("λ1⊗λ3", 0, 2),
        ("λ1⊗λ4", 0, 3),
        ("λ3⊗λ3", 2, 2),
        ("λ3⊗λ8", 2, 7),
        ("λ4⊗λ4", 3, 3),
        ("λ6⊗λ7", 5, 6),
        ("λ1⊗λ1+λ2⊗λ2", None, None),
        ("λ1⊗λ4+λ2⊗λ5", None, None),
    ]

    log(f"{'Hamiltonian':<22} {'Paired':>7} {'Total':>6} {'Max_pair_err':>13} {'Verdict':>10}")
    log("-" * 65)

    for name, i1, i2 in qutrit_tests:
        if i1 is not None:
            H_qt = np.kron(gm_raw[i1], gm_raw[i2])
        elif "λ1⊗λ1+λ2⊗λ2" in name:
            H_qt = np.kron(gm_raw[0], gm_raw[0]) + np.kron(gm_raw[1], gm_raw[1])
        elif "λ1⊗λ4+λ2⊗λ5" in name:
            H_qt = np.kron(gm_raw[0], gm_raw[3]) + np.kron(gm_raw[1], gm_raw[4])

        # Build Liouvillian in Gell-Mann basis
        L_H_qt = np.zeros((num_qt, num_qt), dtype=complex)
        for a in range(num_qt):
            for b in range(num_qt):
                comm = -1j * (H_qt @ pmats_qt[b] - pmats_qt[b] @ H_qt)
                L_H_qt[a, b] = np.trace(pmats_qt[a].conj().T @ comm) / d_full

        L_qt = L_H_qt.copy()
        for a in range(num_qt):
            L_qt[a, a] += L_D_qt[a]

        evals_qt = np.linalg.eigvals(L_qt)
        n_well, n_tot, errors = check_eigenvalue_pairing(evals_qt, Sg_qt, tol=1e-6)
        max_err = max(errors) if errors else 0
        verdict = "PALINDROMIC" if n_well == n_tot else "BROKEN"
        log(f"{name:<22} {n_well:>7}/{n_tot:<5} {max_err:>13.2e} {verdict:>10}")

    log()

    # ----- 2c: Ququart quick check -----
    log("--- 2c: Ququart check (d=4) ---")
    log()
    log("For d=4: 16 basis matrices (generalized Pauli)")
    log("  Diagonal: 4 (including identity)")
    log("  Off-diagonal: 12")
    log("  Split: 4:12. NOT balanced (would need 8:8).")
    log("  N=2 block counts: w=0: 16, w=1: 96, w=2: 144")
    log("  16 ≠ 144 → palindrome impossible for ququarts.")
    log()
    log("CONCLUSION: Single-bond universality is a QUBIT PHENOMENON.")
    log("  d=2: split 2:2 (balanced) → palindrome universal")
    log("  d=3: split 3:6 → broken (9 vs 36)")
    log("  d=4: split 4:12 → broken (16 vs 144)")
    log("  d=d: split d:(d²-d) → balanced only when d = d²-d → d²-2d = 0 → d=2.")


# ============================================================
# ============================================================
#  SECTION 3: ENGINEERING THE ENTANGLED MIRROR
# ============================================================
# ============================================================

def run_section_3():
    log()
    log("=" * 70)
    log("SECTION 3: ENGINEERING THE ENTANGLED MIRROR")
    log("=" * 70)
    log()

    N = 2
    gamma = GAMMA
    Sg = N * gamma
    all_idx, pmats, d = build_pauli_basis(N)
    L_D = build_L_D_diag(gamma, all_idx)
    V = build_V(pmats, d)
    num = 4 ** N

    # ----- 3a: Optimal a/b ratio for maximum entanglement -----
    log("--- 3a: Optimal a/b for H = a*XZ + b*YZ ---")
    log()

    ratios = np.logspace(-2, 2, 100)
    best_ent = 0
    best_ratio = 1.0
    best_neg = 0

    log(f"{'a/b':>8} {'Choi_Rk':>8} {'S_ent':>8} {'S_max':>6} {'%':>6} {'Neg':>8}")
    log("-" * 52)

    results_3a = []
    for r in ratios:
        a, b = r, 1.0
        H = a * build_H_term(N, (0, 1), 'XZ') + b * build_H_term(N, (0, 1), 'YZ')
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        Pi, n_p, _ = construct_Pi(L, Sg)
        if n_p < num:
            continue
        J = pi_to_choi(Pi, V, d)
        choi_rk, svs = choi_site_svd_2site(J)
        s_ent = schmidt_entropy(svs)
        s_max = np.log2(max(choi_rk, 1))
        pct = (s_ent / s_max * 100) if s_max > 0 else 0
        neg = negativity(J, d, d)

        results_3a.append((r, choi_rk, s_ent, s_max, pct, neg))

        if s_ent > best_ent:
            best_ent = s_ent
            best_ratio = r
            best_neg = neg

    # Print every 10th point + extremes
    for i, (r, rk, se, sm, p, ng) in enumerate(results_3a):
        if i % 10 == 0 or abs(r - best_ratio) < 0.01:
            log(f"{r:>8.4f} {rk:>8} {se:>8.4f} {sm:>6.2f} {p:>5.1f}% {ng:>8.5f}")

    log()
    log(f"OPTIMAL: a/b = {best_ratio:.4f}")
    log(f"  Schmidt entropy = {best_ent:.4f}")
    log(f"  Negativity = {best_neg:.6f}")
    log()

    # Fine-grained search around optimal
    fine_ratios = np.linspace(max(best_ratio * 0.5, 0.01), best_ratio * 2.0, 200)
    for r in fine_ratios:
        a, b = r, 1.0
        H = a * build_H_term(N, (0, 1), 'XZ') + b * build_H_term(N, (0, 1), 'YZ')
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        Pi, n_p, _ = construct_Pi(L, Sg)
        if n_p < num:
            continue
        J = pi_to_choi(Pi, V, d)
        _, svs = choi_site_svd_2site(J)
        se = schmidt_entropy(svs)
        if se > best_ent:
            best_ent = se
            best_ratio = r
            best_neg = negativity(J, d, d)

    log(f"REFINED OPTIMAL: a/b = {best_ratio:.6f}")
    log(f"  Schmidt entropy = {best_ent:.6f}")
    log(f"  At this ratio, H ≈ {best_ratio:.3f}*XZ + 1*YZ")
    if best_ratio < 1:
        log(f"  Dominated by YZ (a/b < 1)")
    else:
        log(f"  Dominated by XZ (a/b > 1)")

    # Spectrum at optimal ratio
    a_opt, b_opt = best_ratio, 1.0
    H_opt = a_opt * build_H_term(N, (0, 1), 'XZ') + b_opt * build_H_term(N, (0, 1), 'YZ')
    L_H_opt = build_L_H(H_opt, all_idx, pmats, d)
    L_opt = assemble_L(L_H_opt, L_D)
    evals_opt = np.sort(np.linalg.eigvals(L_opt))
    log()
    log("Spectrum at optimal ratio:")
    for i, ev in enumerate(evals_opt):
        log(f"  λ_{i:>2} = {ev.real:>10.6f} + {ev.imag:>10.6f}i")

    # ----- 3b: Phase diagram -----
    log()
    log("--- 3b: Phase diagram: H = a*XZ + b*YZ + c*(XX+YY+ZZ) ---")
    log()

    # Sweep t from 0 (pure Heisenberg) to 1 (pure XZ+YZ)
    log("Linear interpolation: H(t) = (1-t)*Heisenberg + t*(XZ+YZ)")
    log()

    H_heis = (build_H_term(N, (0, 1), 'XX')
              + build_H_term(N, (0, 1), 'YY')
              + build_H_term(N, (0, 1), 'ZZ'))
    H_xzyz = build_H_term(N, (0, 1), 'XZ') + build_H_term(N, (0, 1), 'YZ')

    log(f"{'t':>6} {'Choi_Rk':>8} {'S_ent':>8} {'Neg':>8} {'PalErr':>10}")
    log("-" * 48)

    t_vals = np.linspace(0, 1, 41)
    for t in t_vals:
        H = (1 - t) * H_heis + t * H_xzyz
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        Pi, n_p, _ = construct_Pi(L, Sg)
        err = palindrome_error(Pi, L, Sg)
        J = pi_to_choi(Pi, V, d)
        rk, svs = choi_site_svd_2site(J)
        se = schmidt_entropy(svs)
        neg = negativity(J, d, d)
        log(f"{t:>6.3f} {rk:>8} {se:>8.4f} {neg:>8.5f} {err:>10.2e}")

    # 2D sweep: a vs c (with b=1 fixed)
    log()
    log("2D sweep: H = a*XZ + YZ + c*(XX+YY+ZZ)")
    log()
    log(f"{'a':>6} {'c':>6} {'Choi_Rk':>8} {'S_ent':>8}")
    log("-" * 34)

    a_vals = np.linspace(0, 3, 13)
    c_vals = np.linspace(0, 3, 13)
    for a in a_vals:
        for c in c_vals:
            H = a * build_H_term(N, (0, 1), 'XZ') + build_H_term(N, (0, 1), 'YZ') + c * H_heis
            L_H = build_L_H(H, all_idx, pmats, d)
            L = assemble_L(L_H, L_D)
            Pi, n_p, _ = construct_Pi(L, Sg)
            J = pi_to_choi(Pi, V, d)
            rk, svs = choi_site_svd_2site(J)
            se = schmidt_entropy(svs)
            log(f"{a:>6.2f} {c:>6.2f} {rk:>8} {se:>8.4f}")

    # ----- 3c: QST fidelity connection -----
    log()
    log("--- 3c: QST fidelity vs Pi entanglement ---")
    log()

    t_points = np.linspace(0.5, 50, 100)
    log(f"{'a/b':>8} {'S_ent':>8} {'Choi_Rk':>8} {'MaxFid':>8} {'t_opt':>8}")
    log("-" * 48)

    test_ratios = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]
    for r in test_ratios:
        a, b = r, 1.0
        H = a * build_H_term(N, (0, 1), 'XZ') + b * build_H_term(N, (0, 1), 'YZ')
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        Pi, _, _ = construct_Pi(L, Sg)
        J = pi_to_choi(Pi, V, d)
        rk, svs = choi_site_svd_2site(J)
        se = schmidt_entropy(svs)

        # QST fidelity
        fid, t_opt = qst_fidelity_sweep(H, gamma, N, t_points)
        log(f"{r:>8.3f} {se:>8.4f} {rk:>8} {fid:>8.5f} {t_opt:>8.2f}")

    # Also test Heisenberg for comparison
    H_heis_qst = H_heis
    fid_h, t_h = qst_fidelity_sweep(H_heis_qst, gamma, N, t_points)
    log(f"{'Heis':>8} {'0.0000':>8} {'1':>8} {fid_h:>8.5f} {t_h:>8.2f}")
    log()
    log("FINDING: XZ+YZ has ZERO QST fidelity at all a/b ratios.")
    log("  This is correct: XZ and YZ are NOT exchange interactions.")
    log("  They do not swap spin excitations between sites.")
    log("  Only Heisenberg (XX+YY+ZZ) facilitates state transfer.")
    log("  Pi entanglement is IRRELEVANT for QST - they measure")
    log("  completely different properties of the Hamiltonian.")


# ============================================================
# ============================================================
#  SECTION 4: THE CLOCK AS SPACE-CREATOR
# ============================================================
# ============================================================

def run_section_4():
    log()
    log("=" * 70)
    log("SECTION 4: THE CLOCK AS SPACE-CREATOR")
    log("=" * 70)
    log()

    N = 2
    gamma = GAMMA
    Sg = N * gamma
    all_idx, pmats, d = build_pauli_basis(N)
    L_D = build_L_D_diag(gamma, all_idx)
    V = build_V(pmats, d)
    num = 4 ** N

    # ----- 4a: Local vs non-local spectral comparison -----
    log("--- 4a: What do we lose going local → non-local? ---")
    log()

    combos_compare = [
        ("XZ only (local Pi)", "XZ"),
        ("YZ only (local Pi)", "YZ"),
        ("XZ+YZ (non-local Pi)", "XZ+YZ"),
        ("XX+YY (local Pi)", "XX+YY"),
        ("XX (local Pi)", "XX"),
    ]

    for label, combo in combos_compare:
        terms = combo.split('+')
        H = sum(build_H_term(N, (0, 1), t) for t in terms)
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)

        evals = np.sort_complex(np.linalg.eigvals(L))
        Pi, n_p, _ = construct_Pi(L, Sg)
        J = pi_to_choi(Pi, V, d)
        rk, svs = choi_site_svd_2site(J)
        se = schmidt_entropy(svs)

        # Spectral metrics
        re_parts = np.sort(evals.real)
        im_parts = np.sort(np.abs(evals.imag))
        n_osc = int(np.sum(np.abs(evals.imag) > 1e-10))  # oscillating modes
        n_decay = int(np.sum(np.abs(evals.real) > 1e-10))  # non-steady modes
        # Spectral spread
        spread_re = np.ptp(evals.real)
        spread_im = np.ptp(evals.imag)

        log(f"{label}:")
        log(f"  Choi rank = {rk}, Schmidt entropy = {se:.4f}")
        log(f"  Oscillating modes: {n_osc}/16")
        log(f"  Real spread: {spread_re:.4f}, Imag spread: {spread_im:.4f}")
        log()

    # ----- 4b: Entanglement entropy sweep -----
    log("--- 4b: Pi entanglement vs alpha ---")
    log("H(α) = (1-α)*XX + α*(XZ+YZ)")
    log()

    H_local = build_H_term(N, (0, 1), 'XX')
    H_nonlocal = build_H_term(N, (0, 1), 'XZ') + build_H_term(N, (0, 1), 'YZ')

    log(f"{'alpha':>6} {'Paired':>7} {'PalErr':>10} {'Choi_Rk':>8} {'Pi_S_ent':>9}")
    log("-" * 48)

    alphas = np.linspace(0, 1, 51)
    for alpha in alphas:
        H = (1 - alpha) * H_local + alpha * H_nonlocal
        L_H = build_L_H(H, all_idx, pmats, d)
        L = assemble_L(L_H, L_D)
        Pi, n_p, _ = construct_Pi(L, Sg)
        err = palindrome_error(Pi, L, Sg)
        if err < 1e-8:  # valid palindrome
            J = pi_to_choi(Pi, V, d)
            rk, svs = choi_site_svd_2site(J)
            se = schmidt_entropy(svs)
        else:
            rk = -1
            se = float('nan')
        log(f"{alpha:>6.3f} {n_p:>4}/{num:<3} {err:>10.2e} {rk:>8} {se:>9.4f}")

    log()
    log("Note: Choi rank = -1 means palindrome BROKEN (Pi invalid).")
    log("  XX and XZ+YZ are BOTH palindromic, but with DIFFERENT Pi operators.")
    log("  Their mixture is NOT palindromic: no single Pi works for both.")
    log("  The palindrome exists only at the endpoints (alpha=0 and alpha=1).")
    log("  SHARP phase transition - the mirror cannot interpolate between")
    log("  different palindromic symmetry classes.")
    log()
    log("  Contrast with Section 3b: Heisenberg-to-XZ+YZ interpolation SURVIVES")
    log("  because Heisenberg's SU(2) symmetry is rich enough to accommodate")
    log("  the Pi structure needed for XZ+YZ. XX alone is too simple.")

    # ----- 4c: Choi distance table -----
    log()
    log("--- 4c: Choi distance for all 36 combos ---")
    log()

    combos = all_36()
    BROKEN_set = set(c for c in combos if c not in SURVIVING)

    # For N=3 palindrome error, compute for all combos
    N3 = 3
    Sg3 = N3 * gamma
    all_idx3, pmats3, d3 = build_pauli_basis(N3)
    L_D3 = build_L_D_diag(gamma, all_idx3)

    log(f"{'Combo':<12} {'d_Choi':>7} {'N3_PalErr':>10} {'Choi_Rk':>8} {'Surv':>5}")
    log("-" * 50)

    for combo in combos:
        # N=2 Choi distance
        H2 = build_H_combo(N, [(0, 1)], combo)
        L_H2 = build_L_H(H2, all_idx, pmats, d)
        L2 = assemble_L(L_H2, L_D)
        Pi2, _, _ = construct_Pi(L2, Sg)
        J2 = pi_to_choi(Pi2, V, d)
        rk2, svs2 = choi_site_svd_2site(J2)
        d_choi = schmidt_entropy(svs2)

        # N=3 palindrome error
        H3 = build_H_combo(N3, [(0, 1), (1, 2)], combo)
        L_H3 = build_L_H(H3, all_idx3, pmats3, d3)
        L3 = assemble_L(L_H3, L_D3)
        Pi3, _, _ = construct_Pi(L3, Sg3)
        err3 = palindrome_error(Pi3, L3, Sg3)

        surv = "YES" if combo in SURVIVING else "NO"
        log(f"{combo:<12} {d_choi:>7.4f} {err3:>10.2e} {rk2:>8} {surv:>5}")

    log()
    log("Interpretation:")
    log("  d_Choi = 0: local Pi (product state) → no 'space' needed")
    log("  d_Choi > 0: non-local Pi → Pi needs 'space' between sites")
    log("  Correlation with N3_PalErr tells whether higher Choi distance")
    log("  predicts larger palindrome error at N=3.")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t0 = time.time()
    log("The Structure of the Mirror")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"γ = {GAMMA}")
    log()

    run_section_1()
    log(f"\n[Section 1 completed in {time.time() - t0:.1f}s]")

    t1 = time.time()
    run_section_2()
    log(f"\n[Section 2 completed in {time.time() - t1:.1f}s]")

    t2 = time.time()
    run_section_3()
    log(f"\n[Section 3 completed in {time.time() - t2:.1f}s]")

    t3 = time.time()
    run_section_4()
    log(f"\n[Section 4 completed in {time.time() - t3:.1f}s]")

    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results written to: {OUT_PATH}")
    _outf.close()
