"""
K-Dosimetry: the exposure number K = gamma*t of a state the Hamiltonian cannot touch
===================================================================================
|+>^N under the isotropic Heisenberg chain with equal Z-dephasing on every site is
Hamiltonian-dead: rho(t) = ((I + g X)/2)^{(x)N} with g = exp(-2 gamma t), so the purity
is ((1 + g^2)/2)^N and the time to reach a purity target P is K/gamma with

    K = -(1/4) * ln(2 * P^(1/N) - 1),

at every gamma and every J. This script checks the closed form against the full
Lindblad propagation (dense expm, root-found crossing), reads the Bell+ and |01>
floors, and runs the one profile that is not Hamiltonian-dead (the N=4 edge
concentrator) with the same exact propagator.

Output: simulations/results/k_dosimetry.txt
"""

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

RESULTS_DIR = Path(__file__).parent / "results"
J = 1.0

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_liouvillian(N, gammas, coupling=J):
    """Row-major vec: L = -i(H(x)I - I(x)H^T) + sum_k D[sqrt(gamma_k) Z_k]."""
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in (Xm, Ym, Zm):
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += coupling * kron_chain(ops)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return L


def plus_state(N):
    p = np.ones(2, dtype=complex) / np.sqrt(2)
    v = p
    for _ in range(N - 1):
        v = np.kron(v, p)
    return np.outer(v, v.conj())


def purity_at(L, rho0, t):
    d = rho0.shape[0]
    rho = (expm(L * t) @ rho0.ravel()).reshape(d, d)
    return float(np.real(np.trace(rho @ rho)))


def crossing_time(L, rho0, target, t_lo, t_hi, n_scan=400):
    """First time the purity falls below target: scan, then brentq on the bracket."""
    ts = np.linspace(t_lo, t_hi, n_scan)
    prev = purity_at(L, rho0, ts[0])
    for a, b in zip(ts[:-1], ts[1:]):
        cur = purity_at(L, rho0, b)
        if prev >= target > cur:
            return brentq(lambda t: purity_at(L, rho0, t) - target, a, b, xtol=1e-14, rtol=1e-14)
        prev = cur
    return None


def k_closed(N, target):
    return -0.25 * np.log(2.0 * target ** (1.0 / N) - 1.0)


def main():
    out = []

    def log(msg=""):
        print(msg)
        out.append(msg)

    log("=" * 75)
    log("K-DOSIMETRY: K = gamma*t FOR A HAMILTONIAN-DEAD STATE")
    log("=" * 75)
    log("closed form (|+>^N, uniform gamma, any J): K = -(1/4) ln(2 P^(1/N) - 1)")
    log()

    log("STEP 1: N=2, |++>, target purity 0.26, gamma sweep (exact propagation)")
    target = 0.26
    K0 = k_closed(2, target)
    log(f"closed form K = {K0:.12f}")
    log(f"{'gamma':>8s} {'t_cross':>14s} {'K = gamma*t':>16s} {'rel. dev.':>12s}")
    rho0 = plus_state(2)
    worst = 0.0
    for gamma in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        L = build_liouvillian(2, [gamma] * 2)
        t = crossing_time(L, rho0, target, 0.0, 3.0 / gamma)
        K = gamma * t
        dev = abs(K - K0) / K0
        worst = max(worst, dev)
        log(f"{gamma:8.3f} {t:14.6f} {K:16.12f} {dev:12.2e}")
    log(f"largest relative deviation from the closed form: {worst:.2e}")
    log()

    log("STEP 2: N=2, gamma=0.05, Bell+ and |01>: purity floor (populations survive)")
    L = build_liouvillian(2, [0.05] * 2)
    for name, v in [("Bell+", np.array([1, 0, 0, 1]) / np.sqrt(2)), ("|01>", np.array([0, 1, 0, 0]))]:
        r = np.outer(v, v.conj()).astype(complex)
        ps = [purity_at(L, r, t) for t in (50.0, 200.0, 800.0)]
        log(f"{name:>6s}: purity at t = 50, 200, 800: " + ", ".join(f"{p:.10f}" for p in ps))
    log()

    log("STEP 3: N=2, |++>, gamma=0.05, purity-target sweep")
    log(f"{'target':>8s} {'K closed':>12s} {'t_cross':>12s} {'K numeric':>12s}")
    L = build_liouvillian(2, [0.05] * 2)
    for tg in [0.90, 0.75, 0.50, 0.40, 0.30, 0.27, 0.26, 0.255, 0.251]:
        t = crossing_time(L, rho0, tg, 0.0, 40.0)
        log(f"{tg:8.3f} {k_closed(2, tg):12.6f} {t:12.6f} {0.05 * t:12.6f}")
    log()

    log("STEP 4: N=2..5, |+>^N, gamma=0.05, target 1/2^N + 0.01")
    log(f"{'N':>3s} {'target':>9s} {'K closed':>12s} {'K numeric':>12s} {'K_system':>12s}")
    t_uniform4 = None
    for N in range(2, 6):
        tg = 1.0 / 2**N + 0.01
        L = build_liouvillian(N, [0.05] * N)
        t = crossing_time(L, plus_state(N), tg, 0.0, 40.0, n_scan=160)
        if N == 4:
            t_uniform4 = t
        log(f"{N:3d} {tg:9.5f} {k_closed(N, tg):12.6f} {0.05 * t:12.6f} {N * 0.05 * t:12.6f}")
    log()

    log("STEP 5: N=4 edge-concentrator profile, same total sum(gamma) = 0.2, target 1/16 + 0.01")
    eps = 0.001
    gs = [4 * 0.05 - 3 * eps] + [eps] * 3
    tg = 1.0 / 16 + 0.01
    L = build_liouvillian(4, gs)
    t_s = crossing_time(L, plus_state(4), tg, 0.0, 40.0, n_scan=160)
    t_free = brentq(lambda t: (1 + np.exp(-4 * gs[0] * t)) / 2 * ((1 + np.exp(-4 * eps * t)) / 2) ** 3 - tg, 1.0, 5000.0)
    log(f"profile {gs}")
    log(f"uniform:   t = {t_uniform4:.6f}, K_total = {0.2 * t_uniform4:.6f}")
    log(f"edge conc.: t = {t_s:.6f}, K_total = {0.2 * t_s:.6f}")
    log(f"K_total ratio (edge concentrator/uniform) = {t_s / t_uniform4:.6f}")
    log(f"same profile without coupling (J = 0, closed form): t = {t_free:.4f}")
    log()

    RESULTS_DIR.mkdir(exist_ok=True)
    (RESULTS_DIR / "k_dosimetry.txt").write_text("\n".join(out) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8")
    main()
