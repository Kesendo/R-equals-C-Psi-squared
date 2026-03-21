#!/usr/bin/env python3
"""
Close the Proof Roadmap: Remaining Computational Gaps
======================================================
Test 1: Generalized Pauli channels (X, Y, Z, depolarizing, asymmetric)
Test 2: Direct amplitude damping CΨ trajectory
Test 3: Non-Markovian revival test
Test 4: Analytical δ(N) for GHZ

Script: simulations/proof_roadmap_close.py
Output: simulations/results/proof_roadmap_close.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "proof_roadmap_close.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sm = np.array([[0, 0], [1, 0]], dtype=complex)  # |0><1| = sigma_minus

N = 2
d = 4
d2 = 16


def site_op(op, k, nq=N):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_L_H(H):
    Id = np.eye(d)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def add_jump(L, F, kappa):
    Id = np.eye(d)
    FdF = F.conj().T @ F
    L += kappa * (np.kron(F, F.conj())
                  - 0.5 * np.kron(FdF, Id)
                  - 0.5 * np.kron(Id, FdF.T))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def purity(rho):
    return np.real(np.trace(rho @ rho))


def l1_coh(rho):
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def cpsi(rho):
    C = purity(rho)
    Psi = l1_coh(rho) / (d - 1)
    return C * Psi


def build_H():
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += site_op(P, 0) @ site_op(P, 1)
    return H


def bell_plus():
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1/np.sqrt(2)
    psi[3] = 1/np.sqrt(2)
    return np.outer(psi, psi.conj())


# ============================================================
# TEST 1: Generalized Pauli Channels
# ============================================================
def run_test_1():
    log("=" * 70)
    log("TEST 1: GENERALIZED PAULI CHANNELS")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_plus()
    t_points = np.linspace(0, 5, 500)

    channels = {
        'Z-dephasing': [(sz, 0.05, 0), (sz, 0.05, 1)],
        'X-noise':     [(sx, 0.05, 0), (sx, 0.05, 1)],
        'Y-noise':     [(sy, 0.05, 0), (sy, 0.05, 1)],
        'Depolarizing': [(sx, 0.05/3, 0), (sy, 0.05/3, 0), (sz, 0.05/3, 0),
                         (sx, 0.05/3, 1), (sy, 0.05/3, 1), (sz, 0.05/3, 1)],
        'Asymmetric':  [(sx, 0.01, 0), (sy, 0.02, 0), (sz, 0.03, 0),
                        (sx, 0.01, 1), (sy, 0.02, 1), (sz, 0.03, 1)],
    }

    log(f"  {'Channel':>14}  {'CPsi(0)':>8}  {'t_cross':>8}  {'CPsi_cross':>10}  {'CPsi(5)':>8}  {'Below?':>7}")
    log(f"  {'-'*60}")

    for name, jumps in channels.items():
        L = build_L_H(H)
        for op, kappa, qubit in jumps:
            F = site_op(op, qubit)
            L = add_jump(L, F, kappa)

        cpsi_traj = []
        for t in t_points:
            rho = evolve(L, rho0, t) if t > 0 else rho0.copy()
            cpsi_traj.append(cpsi(rho))

        cpsi_traj = np.array(cpsi_traj)
        cpsi0 = cpsi_traj[0]
        cpsi_end = cpsi_traj[-1]

        # Find crossing
        t_cross = None
        cpsi_at_cross = None
        for i in range(1, len(cpsi_traj)):
            if cpsi_traj[i-1] > 0.25 and cpsi_traj[i] <= 0.25:
                frac = (cpsi_traj[i-1] - 0.25) / (cpsi_traj[i-1] - cpsi_traj[i])
                t_cross = t_points[i-1] + frac * (t_points[i] - t_points[i-1])
                cpsi_at_cross = 0.25
                break

        # Check if stays below after crossing
        stays_below = True
        if t_cross is not None:
            idx_cross = np.argmin(np.abs(t_points - t_cross))
            if np.any(cpsi_traj[idx_cross:] > 0.251):
                stays_below = False

        tc_str = f"{t_cross:.4f}" if t_cross else "NEVER"
        cc_str = f"{cpsi_at_cross:.4f}" if cpsi_at_cross else "N/A"
        below = "YES" if stays_below and t_cross else ("NO CROSS" if not t_cross else "NO")

        log(f"  {name:>14}  {cpsi0:>8.4f}  {tc_str:>8}  {cc_str:>10}  {cpsi_end:>8.4f}  {below:>7}")

    log()


# ============================================================
# TEST 2: Direct Amplitude Damping
# ============================================================
def run_test_2():
    log("=" * 70)
    log("TEST 2: DIRECT AMPLITUDE DAMPING")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_plus()
    t_points = np.linspace(0, 10, 500)

    log(f"  {'gamma':>8}  {'CPsi(0)':>8}  {'t_cross':>8}  {'CPsi(10)':>8}  {'Z-deph t_cross':>15}")
    log(f"  {'-'*55}")

    for gamma in [0.01, 0.05, 0.10]:
        # Amplitude damping
        L_ad = build_L_H(H)
        for q in range(N):
            F = site_op(sm, q)
            L_ad = add_jump(L_ad, F, gamma)

        # Z-dephasing reference
        L_zd = build_L_H(H)
        for q in range(N):
            F = site_op(sz, q)
            L_zd = add_jump(L_zd, F, gamma)

        t_cross_ad = None
        t_cross_zd = None

        for i, t in enumerate(t_points):
            if t == 0: continue
            rho_ad = evolve(L_ad, rho0, t)
            rho_zd = evolve(L_zd, rho0, t)
            cp_ad = cpsi(rho_ad)
            cp_zd = cpsi(rho_zd)

            if t_cross_ad is None and i > 0:
                rho_prev = evolve(L_ad, rho0, t_points[i-1])
                cp_prev = cpsi(rho_prev)
                if cp_prev > 0.25 and cp_ad <= 0.25:
                    frac = (cp_prev - 0.25) / (cp_prev - cp_ad)
                    t_cross_ad = t_points[i-1] + frac * (t - t_points[i-1])

            if t_cross_zd is None and i > 0:
                rho_prev = evolve(L_zd, rho0, t_points[i-1])
                cp_prev = cpsi(rho_prev)
                if cp_prev > 0.25 and cp_zd <= 0.25:
                    frac = (cp_prev - 0.25) / (cp_prev - cp_zd)
                    t_cross_zd = t_points[i-1] + frac * (t - t_points[i-1])

        rho_end = evolve(L_ad, rho0, 10.0)
        cp_end = cpsi(rho_end)

        tc_ad = f"{t_cross_ad:.4f}" if t_cross_ad else "NEVER"
        tc_zd = f"{t_cross_zd:.4f}" if t_cross_zd else "NEVER"

        log(f"  {gamma:>8.2f}  {cpsi(rho0):>8.4f}  {tc_ad:>8}  {cp_end:>8.6f}  {tc_zd:>15}")

    log()


# ============================================================
# TEST 3: Non-Markovian Revival
# ============================================================
def run_test_3():
    log("=" * 70)
    log("TEST 3: NON-MARKOVIAN REVIVAL")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_plus()
    gamma = 0.05

    # Phase 1: dephase until well below 1/4
    L_deph = build_L_H(H)
    for q in range(N):
        L_deph = add_jump(L_deph, site_op(sz, q), gamma)

    # Find crossing time
    t_fine = np.linspace(0, 2, 200)
    t_cross = None
    for i, t in enumerate(t_fine[1:], 1):
        rho = evolve(L_deph, rho0, t)
        if cpsi(rho) <= 0.25:
            t_cross = t
            break

    if t_cross is None:
        log("  CΨ never crossed 1/4. Cannot test revival.")
        return

    t1 = 2 * t_cross  # well below 1/4
    rho_t1 = evolve(L_deph, rho0, t1)
    cp_before = cpsi(rho_t1)
    log(f"  Phase 1: dephased to t={t1:.3f}, CΨ = {cp_before:.6f} (below 1/4)")
    log()

    # Phase 2: unitary revival pulses
    log(f"  {'theta':>8}  {'CPsi_after':>12}  {'Above 1/4?':>12}  {'Stays above?':>13}")
    log(f"  {'-'*50}")

    XX_YY = site_op(sx, 0) @ site_op(sx, 1) + site_op(sy, 0) @ site_op(sy, 1)

    for theta in np.linspace(0, np.pi, 20):
        U = expm(-1j * theta * XX_YY)
        rho_revived = U @ rho_t1 @ U.conj().T
        rho_revived = (rho_revived + rho_revived.conj().T) / 2
        cp_revived = cpsi(rho_revived)
        above = cp_revived > 0.25

        # Phase 3: resume dephasing, check if stays above
        stays = False
        if above:
            stays = True
            for t3 in np.linspace(0.01, 2.0, 100):
                rho_p3 = evolve(L_deph, rho_revived, t3)
                if cpsi(rho_p3) <= 0.25:
                    stays = False
                    break

        log(f"  {theta:>8.4f}  {cp_revived:>12.6f}  {'YES' if above else 'no':>12}  "
            f"{'YES' if stays else 'no':>13}")

    log()


# ============================================================
# TEST 4: Analytical δ(N) for GHZ
# ============================================================
def run_test_4():
    log("=" * 70)
    log("TEST 4: ANALYTICAL delta(N) FOR GHZ")
    log("=" * 70)
    log()

    gamma = 0.1
    t_eval = 1.0

    log(f"  {'N':>3}  {'CPsi_analytic':>14}  {'CPsi_Lindblad':>14}  {'delta':>10}  {'t_cross':>10}  {'Match?':>8}")
    log(f"  {'-'*65}")

    for nq in range(2, 7):  # N=7 Liouvillian is 8GB, N=8 is 73GB
        dq = 2 ** nq

        # Analytical: GHZ under Z-dephasing
        # rho has diagonal 1/2^N each, off-diag rho[0, 2^N-1] = (1/2)*exp(-N*gamma*t)
        # Purity = 1/2^N + (1/2)*exp(-2*N*gamma*t)
        # l1 = 2 * (1/2) * exp(-N*gamma*t) = exp(-N*gamma*t)
        # Psi = l1 / (d^2 - 1) = exp(-N*gamma*t) / (4^N - 1)
        # CPsi = Purity * Psi

        C_anal = 1.0 / dq + 0.5 * np.exp(-2 * nq * gamma * t_eval)
        l1_anal = np.exp(-nq * gamma * t_eval)
        Psi_anal = l1_anal / (dq**2 - 1)
        CPsi_anal = C_anal * Psi_anal

        # Lindblad numerical
        dq2 = dq * dq
        H_nq = np.zeros((dq, dq), dtype=complex)
        for i in range(nq):
            j = (i + 1) % nq
            for P in [sx, sy, sz]:
                opi = [I2] * nq; opi[i] = P
                opj = [I2] * nq; opj[j] = P
                Pi = opi[0]
                for o in opi[1:]: Pi = np.kron(Pi, o)
                Pj = opj[0]
                for o in opj[1:]: Pj = np.kron(Pj, o)
                H_nq += Pi @ Pj

        Id_nq = np.eye(dq)
        L_nq = -1j * (np.kron(H_nq, Id_nq) - np.kron(Id_nq, H_nq.T))
        for k in range(nq):
            ops = [I2] * nq; ops[k] = sz
            Zk = ops[0]
            for o in ops[1:]: Zk = np.kron(Zk, o)
            L_nq += gamma * (np.kron(Zk, Zk.conj()) - np.eye(dq2))

        psi_ghz = np.zeros(dq, dtype=complex)
        psi_ghz[0] = 1/np.sqrt(2)
        psi_ghz[dq-1] = 1/np.sqrt(2)
        rho0_ghz = np.outer(psi_ghz, psi_ghz.conj())

        rho_t = (expm(L_nq * t_eval) @ rho0_ghz.flatten()).reshape(dq, dq)
        rho_t = (rho_t + rho_t.conj().T) / 2

        C_num = np.real(np.trace(rho_t @ rho_t))
        l1_num = np.sum(np.abs(rho_t)) - np.sum(np.abs(np.diag(rho_t)))
        Psi_num = l1_num / (dq**2 - 1)
        CPsi_num = C_num * Psi_num

        delta = CPsi_num - CPsi_anal
        match = abs(delta) < 1e-6

        # Crossing time (analytical)
        # CPsi(t) = [1/2^N + (1/2)*exp(-2Ngt)] * exp(-Ngt) / (4^N - 1) = 1/4
        # Solve numerically
        from scipy.optimize import brentq
        def cpsi_of_t(t):
            C = 1.0/dq + 0.5*np.exp(-2*nq*gamma*t)
            Psi = np.exp(-nq*gamma*t) / (dq**2 - 1)
            return C * Psi - 0.25

        try:
            if cpsi_of_t(0) > 0:
                tc = brentq(cpsi_of_t, 0, 50)
            else:
                tc = 0
        except ValueError:
            tc = float('inf')

        log(f"  {nq:>3}  {CPsi_anal:>14.8f}  {CPsi_num:>14.8f}  {delta:>10.2e}  {tc:>10.4f}  "
            f"{'YES' if match else 'NO':>8}")

    log()
    log("  Analytical formula:")
    log("    C(t) = 1/2^N + (1/2)*exp(-2*N*gamma*t)")
    log("    Psi(t) = exp(-N*gamma*t) / (4^N - 1)")
    log("    CPsi(t) = C(t) * Psi(t)")
    log()


# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Close the Proof Roadmap: Remaining Computational Gaps")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    t1 = time.time()
    run_test_1()
    log(f"[Test 1 completed in {time.time()-t1:.1f}s]")
    log()

    t2 = time.time()
    run_test_2()
    log(f"[Test 2 completed in {time.time()-t2:.1f}s]")
    log()

    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time()-t3:.1f}s]")
    log()

    t4 = time.time()
    run_test_4()
    log(f"[Test 4 completed in {time.time()-t4:.1f}s]")
    log()

    log(f"Total runtime: {time.time()-t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
