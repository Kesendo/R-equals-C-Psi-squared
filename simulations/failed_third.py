#!/usr/bin/env python3
"""
The Failed Third: Origin of Noise from Qubit Decay
=====================================================
Test 1: Does Q3 decay create effective dephasing on Q1-Q2?
Test 2: Is the effective noise Z-like?
Test 3: Does a palindrome form in Q1-Q2?
Test 4: Signature of the failed third (Markovian vs non-Markovian)
Test 5: The 2+1 split

Script: simulations/failed_third.py
Output: simulations/results/failed_third.txt
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "failed_third.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# INFRASTRUCTURE (N=3)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sp = (sx + 1j * sy) / 2  # sigma_plus = |1><0|
sm = (sx - 1j * sy) / 2  # sigma_minus = |0><1|
PAULIS = [I2, sx, sy, sz]
PNAMES = ['I', 'X', 'Y', 'Z']

N = 3
d = 8
d2 = 64


def site_op(op, k, nq=N):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H(J=1.0):
    """Heisenberg coupling between all neighbors: 0-1, 1-2, 0-2."""
    H = np.zeros((d, d), dtype=complex)
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, j)
    return H


def build_L_H(H):
    Id = np.eye(d)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def add_lindblad_jump(L, F, kappa):
    Id = np.eye(d)
    FdF = F.conj().T @ F
    L += kappa * (np.kron(F, F.conj())
                  - 0.5 * np.kron(FdF, Id)
                  - 0.5 * np.kron(Id, FdF.T))
    return L


def build_L_D_zdeph(gammas):
    L_D = np.zeros((d2, d2), dtype=complex)
    for k in range(N):
        Zk = site_op(sz, k)
        L_D += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L_D


def ptrace_keep(rho, nq, keep):
    dim = 2 ** nq; nk = len(keep); dk = 2 ** nk
    traced = [k for k in range(nq) if k not in keep]
    rho_r = np.zeros((dk, dk), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (nq - 1 - k)) & 1 for k in range(nq)]
            bj = [(j >> (nq - 1 - k)) & 1 for k in range(nq)]
            if all(bi[k] == bj[k] for k in traced):
                ki = sum(bi[keep[m]] << (nk - 1 - m) for m in range(nk))
                kj = sum(bj[keep[m]] << (nk - 1 - m) for m in range(nk))
                rho_r[ki, kj] += rho[i, j]
    return rho_r


def palindrome_err(evals, Sg, tol=1e-6):
    n = len(evals); max_err = 0; n_paired = 0
    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        best = np.min(np.abs(evals - target))
        if best < tol: n_paired += 1
        if best > max_err: max_err = best
    return max_err, n_paired, n


def best_palindrome(evals, Sg_base):
    best_np = 0; best_Sg = Sg_base; best_err = 1e30
    candidates = np.concatenate([[Sg_base], np.linspace(0, max(Sg_base * 5, 1.0), 100)])
    for Sg_try in candidates:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=1e-4)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try
    fine = np.linspace(best_Sg - 0.1, best_Sg + 0.1, 300)
    for Sg_try in fine:
        err_try, np_try, _ = palindrome_err(evals, Sg_try, tol=1e-10)
        if np_try > best_np or (np_try == best_np and err_try < best_err):
            best_err = err_try; best_np = np_try; best_Sg = Sg_try
    return best_err, best_np, best_Sg


def evolve(L, rho0, t):
    rho_vec = expm(L * t) @ rho0.flatten()
    rho = rho_vec.reshape(d, d)
    return (rho + rho.conj().T) / 2


# ============================================================
# BUILD ORIGIN MODELS (Q3 unstable, Q1-Q2 no external noise)
# ============================================================
def build_origin_model(option, J=1.0, kappa=0.1, h=10.0):
    """Build Liouvillian for the origin model.
    Q1-Q2 have NO dephasing. Q3 has an instability."""

    H = build_H(J)

    if option == 'A':
        # Amplitude damping on Q3
        H_total = H
        L = build_L_H(H_total)
        F = site_op(sm, 2)  # sigma_minus on Q3
        L = add_lindblad_jump(L, F, kappa)

    elif option == 'B':
        # Strong detuning on Q3
        H_total = H + h * site_op(sz, 2)
        L = build_L_H(H_total)

    elif option == 'C':
        # Thermal bath on Q3
        H_total = H
        L = build_L_H(H_total)
        F_down = site_op(sm, 2)
        F_up = site_op(sp, 2)
        L = add_lindblad_jump(L, F_down, kappa)
        L = add_lindblad_jump(L, F_up, kappa * 0.5)  # detailed balance

    elif option == 'D':
        # Q3 with strong X+Y dephasing (destroys Z-coherence)
        H_total = H
        L = build_L_H(H_total)
        Fx = site_op(sx, 2)
        Fy = site_op(sy, 2)
        L = add_lindblad_jump(L, Fx, kappa)
        L = add_lindblad_jump(L, Fy, kappa)

    return L


# ============================================================
# TEST 1: Effective dephasing on Q1-Q2
# ============================================================
def run_test_1():
    log("=" * 70)
    log("TEST 1: DOES Q3 DECAY CREATE EFFECTIVE DEPHASING ON Q1-Q2?")
    log("=" * 70)
    log()

    # Initial state: Bell(Q1,Q2) x |0>_Q3
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)       # |000>
    psi[0b110] = 1 / np.sqrt(2)   # |110>
    rho0 = np.outer(psi, psi.conj())

    t_points = np.linspace(0, 20, 100)
    options = ['A', 'B', 'C', 'D']
    option_names = {
        'A': 'Amplitude damping on Q3',
        'B': 'Strong detuning on Q3 (h=10)',
        'C': 'Thermal bath on Q3',
        'D': 'X+Y dephasing on Q3',
    }

    results = {}

    for opt in options:
        log(f"--- Option {opt}: {option_names[opt]} ---")
        log()

        L = build_origin_model(opt)

        purities = []
        offdiag_01_10 = []

        log(f"  {'t':>6}  {'Purity12':>10}  {'|rho01_10|':>12}  {'|rho00_11|':>12}")
        log(f"  {'-'*45}")

        for t in t_points:
            rho_t = evolve(L, rho0, t)
            rho_12 = ptrace_keep(rho_t, N, [0, 1])

            pur = np.real(np.trace(rho_12 @ rho_12))
            purities.append(pur)
            offdiag_01_10.append(abs(rho_12[1, 2]))  # |01><10|

            if t < 0.01 or t % 2 < 0.02:
                log(f"  {t:>6.1f}  {pur:>10.6f}  {abs(rho_12[1,2]):>12.6f}  "
                    f"{abs(rho_12[0,3]):>12.6f}")

        purities = np.array(purities)
        offdiag = np.array(offdiag_01_10)

        # Fit exponential decay: |rho_01_10|(t) ~ A * exp(-gamma_eff * t)
        mask = offdiag > 1e-6
        if np.sum(mask) > 5 and offdiag[0] > 1e-6:
            log_offdiag = np.log(offdiag[mask] / offdiag[0])
            t_masked = t_points[mask]
            # Linear fit: log(offdiag) = -gamma_eff * t
            coeffs = np.polyfit(t_masked, log_offdiag, 1)
            gamma_eff = -coeffs[0]
        else:
            gamma_eff = 0

        purity_drop = purities[0] - purities[-1]
        results[opt] = gamma_eff

        log()
        log(f"  Purity drop: {purity_drop:.6f}")
        log(f"  Effective dephasing rate gamma_eff: {gamma_eff:.6f}")
        if gamma_eff > 0.001:
            log(f"  -> Q3 decay CREATES decoherence on Q1-Q2")
        else:
            log(f"  -> No significant decoherence")
        log()

    log("--- Summary ---")
    log()
    for opt in options:
        log(f"  Option {opt}: gamma_eff = {results[opt]:.6f}")
    log()

    return results


# ============================================================
# TEST 2: Noise type characterization
# ============================================================
def run_test_2():
    log("=" * 70)
    log("TEST 2: IS THE EFFECTIVE NOISE Z-LIKE?")
    log("=" * 70)
    log()

    # Build 2-qubit Pauli basis
    pauli_2q = []
    pauli_names_2q = []
    for i in range(4):
        for j in range(4):
            pauli_2q.append(np.kron(PAULIS[i], PAULIS[j]))
            pauli_names_2q.append(PNAMES[i] + PNAMES[j])

    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[0b110] = 1 / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    t_measure = 10.0
    options = ['A', 'B', 'C', 'D']

    for opt in options:
        L = build_origin_model(opt)
        rho_t = evolve(L, rho0, t_measure)
        rho_12 = ptrace_keep(rho_t, N, [0, 1])
        rho_12_0 = ptrace_keep(rho0, N, [0, 1])

        log(f"--- Option {opt} at t={t_measure} ---")
        log()
        log(f"  {'Pauli':>6}  {'t=0':>10}  {'t={}'.format(int(t_measure)):>10}  {'Ratio':>10}  {'Status':>10}")
        log(f"  {'-'*50}")

        survived = []
        decayed = []

        for idx, (P, name) in enumerate(zip(pauli_2q, pauli_names_2q)):
            val_0 = np.real(np.trace(rho_12_0 @ P)) / 4
            val_t = np.real(np.trace(rho_12 @ P)) / 4
            ratio = abs(val_t / val_0) if abs(val_0) > 1e-8 else 0

            if abs(val_0) > 1e-4:
                status = "SURVIVES" if ratio > 0.5 else "DECAYS"
                if ratio > 0.5:
                    survived.append(name)
                else:
                    decayed.append(name)
                log(f"  {name:>6}  {val_0:>10.4f}  {val_t:>10.4f}  {ratio:>10.4f}  {status:>10}")

        log()
        # Determine noise type
        surv_set = set(survived)
        if all(n[0] in ('I', 'Z') and n[1] in ('I', 'Z') for n in survived if len(survived) > 0):
            noise_type = "Z-DEPHASING (palindrome possible)"
        elif all(n[0] in ('I', 'X') and n[1] in ('I', 'X') for n in survived if len(survived) > 0):
            noise_type = "X-DEPHASING"
        else:
            noise_type = "MIXED/OTHER"
        log(f"  Survived: {', '.join(survived)}")
        log(f"  Decayed: {', '.join(decayed)}")
        log(f"  Noise type: {noise_type}")
        log()


# ============================================================
# TEST 3: Palindrome in Q1-Q2 effective dynamics
# ============================================================
def run_test_3():
    log("=" * 70)
    log("TEST 3: DOES A PALINDROME FORM IN Q1-Q2?")
    log("=" * 70)
    log()

    # Process tomography: prepare Q1-Q2 in each Pauli basis state,
    # propagate full system, trace Q3, measure output
    pauli_2q = []
    for i in range(4):
        for j in range(4):
            pauli_2q.append(np.kron(PAULIS[i], PAULIS[j]))

    options = ['A', 'B', 'C', 'D']

    for opt in options:
        log(f"--- Option {opt}: Process tomography ---")
        log()

        L = build_origin_model(opt)
        t_eff = 5.0  # measure at this time

        # Build effective transfer matrix in Pauli basis
        # T[a, b] = Tr(P_a @ Lambda(P_b)) / 4
        # where Lambda is the effective Q1-Q2 channel
        T = np.zeros((16, 16), dtype=complex)

        for b in range(16):
            # Input: P_b / 4 (normalized Pauli state)
            rho_in_12 = pauli_2q[b] / 4.0
            # Embed in 3-qubit space: rho_12 ⊗ |0><0|_Q3
            q3_state = np.array([[1, 0], [0, 0]], dtype=complex)
            rho_in = np.kron(rho_in_12, q3_state)
            # Propagate
            rho_out = evolve(L, rho_in, t_eff)
            # Trace Q3
            rho_out_12 = ptrace_keep(rho_out, N, [0, 1])
            # Measure in Pauli basis
            for a in range(16):
                T[a, b] = np.trace(pauli_2q[a] @ rho_out_12) / 4.0

        # T is the effective transfer matrix = exp(L_eff * t)
        # L_eff = log(T) / t (if T is invertible)
        T_real = np.real(T)  # should be real for physical channel

        try:
            from scipy.linalg import logm
            L_eff = logm(T_real) / t_eff
            evals_eff = np.linalg.eigvals(L_eff)

            # Check palindrome
            b_err, b_np, b_Sg = best_palindrome(evals_eff, 0.1)

            log(f"  L_eff eigenvalues (16):")
            rates = -np.real(evals_eff)
            for i, ev in enumerate(sorted(evals_eff, key=lambda x: x.real)):
                log(f"    {i:>3}  Re={ev.real:>10.6f}  Im={ev.imag:>10.6f}")

            log()
            log(f"  Palindromic pairs: {b_np}/16")
            log(f"  Pairing error: {b_err:.2e}")
            log(f"  Best Sg: {b_Sg:.6f}")

            if b_np == 16:
                log(f"  *** PALINDROME FORMS ***")
            elif b_np > 10:
                log(f"  Partial palindrome ({b_np}/16)")
            else:
                log(f"  No palindrome")
        except Exception as e:
            log(f"  Process tomography failed: {e}")

        log()

    # Also check: standard Q1-Q2 with hand-tuned Z-dephasing for comparison
    log("--- Reference: standard Z-dephasing on Q1-Q2 (no Q3) ---")
    log()
    nq2 = 2; d2q = 4; d2q2 = 16
    H2 = np.zeros((d2q, d2q), dtype=complex)
    for P in [sx, sy, sz]:
        H2 += np.kron(P, P)
    Id2 = np.eye(d2q)
    L2 = -1j * (np.kron(H2, Id2) - np.kron(Id2, H2.T))
    gamma_ref = 0.1
    for k in range(2):
        ops = [I2, I2]; ops[k] = sz
        Zk = ops[0]
        for o in ops[1:]: Zk = np.kron(Zk, o)
        L2 += gamma_ref * (np.kron(Zk, Zk.conj()) - np.eye(d2q2))

    evals_ref = np.linalg.eigvals(L2)
    b_err, b_np, b_Sg = best_palindrome(evals_ref, 2 * gamma_ref)
    log(f"  Standard N=2, gamma={gamma_ref}: {b_np}/16 palindromic, "
        f"err={b_err:.2e}, Sg={b_Sg:.6f}")
    log()


# ============================================================
# TEST 4: Markovian vs non-Markovian signature
# ============================================================
def run_test_4():
    log("=" * 70)
    log("TEST 4: MARKOVIAN VS NON-MARKOVIAN SIGNATURE")
    log("=" * 70)
    log()

    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[0b110] = 1 / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    t_points = np.linspace(0, 20, 200)

    for opt in ['A', 'B', 'C', 'D']:
        L = build_origin_model(opt)

        # Track trace distance between rho_12(t) and Markovian prediction
        offdiag = []
        for t in t_points:
            rho_t = evolve(L, rho0, t)
            rho_12 = ptrace_keep(rho_t, N, [0, 1])
            offdiag.append(abs(rho_12[1, 2]))

        offdiag = np.array(offdiag)

        # Fit exponential to first 20% of data
        n_fit = len(t_points) // 5
        if offdiag[0] > 1e-6 and np.all(offdiag[:n_fit] > 1e-10):
            log_od = np.log(offdiag[:n_fit] / offdiag[0])
            coeffs = np.polyfit(t_points[:n_fit], log_od, 1)
            gamma_eff = -coeffs[0]
            markov_pred = offdiag[0] * np.exp(-gamma_eff * t_points)
        else:
            gamma_eff = 0
            markov_pred = np.full_like(t_points, offdiag[0])

        # BLP non-Markovianity: count increases in trace distance
        increases = 0
        for i in range(1, len(offdiag)):
            if offdiag[i] > offdiag[i-1] + 1e-8:
                increases += 1

        # Deviation from Markovian
        if gamma_eff > 0:
            deviation = np.max(np.abs(offdiag - markov_pred))
        else:
            deviation = 0

        log(f"  Option {opt}: gamma_eff={gamma_eff:.6f}, "
            f"non-Markov increases={increases}/{len(offdiag)-1}, "
            f"max deviation={deviation:.6f}")

    log()


# ============================================================
# TEST 5: The 2+1 split
# ============================================================
def run_test_5():
    log("=" * 70)
    log("TEST 5: THE 2+1 SPLIT")
    log("=" * 70)
    log()

    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[0b110] = 1 / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    t_measure = 10.0

    # (a) 1 fails (Q3), 2 survive (Q1-Q2)
    log("--- (a) 1 fails, 2 survive ---")
    L_a = build_origin_model('A')
    rho_a = evolve(L_a, rho0, t_measure)
    rho_12_a = ptrace_keep(rho_a, N, [0, 1])
    pur_a = np.real(np.trace(rho_12_a @ rho_12_a))
    log(f"  Q1-Q2 purity at t={t_measure}: {pur_a:.6f}")
    log(f"  Q1-Q2 off-diag |01><10|: {abs(rho_12_a[1,2]):.6f}")
    log()

    # (b) 2 fail (Q2, Q3), 1 survives (Q1)
    log("--- (b) 2 fail, 1 survives ---")
    H = build_H()
    L_b = build_L_H(H)
    L_b = add_lindblad_jump(L_b, site_op(sm, 1), 0.1)  # Q2 decays
    L_b = add_lindblad_jump(L_b, site_op(sm, 2), 0.1)  # Q3 decays
    rho_b = evolve(L_b, rho0, t_measure)
    rho_1 = ptrace_keep(rho_b, N, [0])  # single qubit
    pur_1 = np.real(np.trace(rho_1 @ rho_1))
    log(f"  Q1 state at t={t_measure}:")
    log(f"    rho = [[{rho_1[0,0].real:.4f}, {rho_1[0,1]:.4f}],")
    log(f"           [{rho_1[1,0]:.4f}, {rho_1[1,1].real:.4f}]]")
    log(f"  Purity: {pur_1:.6f}")
    log(f"  (N=1 has no palindrome — trivial)")
    log()

    # (c) 0 fail, 3 survive (pure Hamiltonian, no noise)
    log("--- (c) 0 fail, 3 survive (no noise) ---")
    L_c = build_L_H(build_H())
    rho_c = evolve(L_c, rho0, t_measure)
    rho_12_c = ptrace_keep(rho_c, N, [0, 1])
    pur_c = np.real(np.trace(rho_12_c @ rho_12_c))
    log(f"  Q1-Q2 purity at t={t_measure}: {pur_c:.6f}")
    log(f"  (No dephasing → no palindrome, just oscillation)")
    log()

    log("--- Summary: which split creates structure? ---")
    log(f"  0 fail, 3 survive: purity {pur_c:.4f} (no dephasing, no palindrome)")
    log(f"  1 fail, 2 survive: purity {pur_a:.4f} (effective dephasing, palindrome?)")
    log(f"  2 fail, 1 survive: purity {pur_1:.4f} (single qubit, trivial)")
    log()


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("The Failed Third: Origin of Noise from Qubit Decay")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, d={d}, d^2={d2}")
    log()

    t1 = time.time()
    results_1 = run_test_1()
    log(f"[Test 1 completed in {time.time() - t1:.1f}s]")
    log()

    t2 = time.time()
    run_test_2()
    log(f"[Test 2 completed in {time.time() - t2:.1f}s]")
    log()

    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time() - t3:.1f}s]")
    log()

    t4 = time.time()
    run_test_4()
    log(f"[Test 4 completed in {time.time() - t4:.1f}s]")
    log()

    t5 = time.time()
    run_test_5()
    log(f"[Test 5 completed in {time.time() - t5:.1f}s]")
    log()

    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
