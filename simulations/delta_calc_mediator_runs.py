"""
The mediator page's tool numbers, in the retired tool's reading and exactly
==========================================================================
hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md reports runs of the retired delta_calc MCP tool (its
source is kept outside the repo and is transcribed here, not imported): the C_int / C_ext
comparison of §2.3 and Appendix A.4, the local-noise GHZ_3 mediator of Appendix A.1, and the
N-scaling of §7.1 and Appendix A.5. The Claude Desktop log of the chat's MCP holds one of the
comparison calls (2026-02-11 17:17 UTC, gamma 0.1); the others are regenerated from the page's
settings and the tool's defaults (J = 1, h = 0).

compute_delta_cint evolves Bell+ under H = J (XX + YY + ZZ) with sigma_z on site A at gamma and on
site B at gamma (C_int) or 0 (C_ext), by an Euler step of dt 0.01 without clipping, to t = 1, and
reports delta = purity - (diag + offdiag * exp(-2 * total_gamma * t)), total_gamma = 2 gamma for
C_int and gamma for C_ext. Bell+ and Bell- are eigenstates of the bond with the same energy and
dephasing keeps the state in their span, so H does nothing: the coherence |00><11| decays at
2 gamma per dephased site, the purity is 1/2 + 1/2 exp(-8 gamma t) (C_int) or 1/2 + 1/2 exp(-4 gamma t)
(C_ext), and the subtracted curve is the purity that dephasing at half that rate would give. So
exactly

    delta_int = (exp(-8 gamma t) - exp(-4 gamma t)) / 2,   delta_ext = (exp(-4 gamma t) - exp(-2 gamma t)) / 2,

and both vanish against a prediction at each run's own rate; the rates of the two sites add.

Appendix A.1 is simulate_dynamic_lindblad with noise_type local (the scalar law, rate gamma * C with
the mutual-purity bridge C, sigma_z on every site, Euler dt 0.01 with the clipping), GHZ_3 on the
Heisenberg ring, gamma 0.05, t_max 10. Appendix A.5 is GHZ_N on the ring under operator feedback,
which is blind on GHZ_N (A.2), so delta is the GHZ coherence's dephasing against the tool's
half-rate curve, (exp(-4 N gamma t) - exp(-2 N gamma t)) / 2 at t = 5.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

from functools import reduce
import sys

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
EPS = np.finfo(float).eps

# The logged C_int / C_ext call (2026-02-11 17:17 UTC, Bell+, Heisenberg, t 1, gamma 0.1).
LOG_CINT = {"c_int purity": 0.724304, "c_int delta": -0.110856, "c_ext purity": 0.835026,
            "c_ext delta": -0.07434, "difference": -0.036516}

CLAIMS, CHECKS = [], []


def site(P, k, n):
    return reduce(np.kron, [P if i == k else I2 for i in range(n)])


def heisenberg(n, ring):
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i in range(n if ring else n - 1):
        j = (i + 1) % n
        for P in (X, Y, Z):
            H += site(P, i, n) @ site(P, j, n)
    return H


def pure(v):
    v = np.asarray(v, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def bell():
    return pure([1, 0, 0, 1])


def ghz(n):
    v = np.zeros(2**n)
    v[0] = v[-1] = 1
    return pure(v)


def purity(r):
    return float(np.real(np.trace(r @ r)))


def one_site(rho, n, k):
    t = rho.reshape([2] * (2 * n))
    keep = [k, n + k]
    rest = [a for a in range(2 * n) if a not in keep]
    m = 2 ** (n - 1)
    return np.einsum("abii->ab", np.transpose(t, keep + rest).reshape(2, 2, m, m))


def mutual_purity(rho, n):
    return float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]) ** (1.0 / n))


def psi(r):
    return float((np.abs(r).sum() - np.abs(np.diag(r)).sum()) / (r.shape[0] - 1))


def prediction(rho0, t, total_gamma):
    diag = float(np.sum(np.abs(np.diag(rho0)) ** 2))
    return diag + (purity(rho0) - diag) * np.exp(-2 * total_gamma * t)


def tool_cint(gamma, gammas, t=1.0, dt=0.01):
    """compute_delta_cint, transcribed: plain Euler, no clipping."""
    rho0, H = bell(), heisenberg(2, False)
    L = [(g, site(Z, k, 2)) for k, g in enumerate(gammas) if g > 0]
    rho = rho0.copy()
    for _ in range(int(t / dt)):
        d = -1j * (H @ rho - rho @ H)
        for g, Lk in L:
            d += g * (Lk @ rho @ Lk - rho)
        rho = rho + dt * d
    rho = 0.5 * (rho + rho.conj().T)
    p = purity(rho)
    return p, p - prediction(rho0, t, sum(gammas))


def exact_purity(rho0, H, gammas, t):
    n = int(np.log2(rho0.shape[0]))
    d = rho0.shape[0]
    Lv = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d)))
    for k, g in enumerate(gammas):
        Lk = site(Z, k, n)
        Lv += g * (np.kron(Lk.conj(), Lk) - np.kron(np.eye(d), np.eye(d)))
    return purity((expm(Lv * t) @ rho0.reshape(-1, order="F")).reshape(d, d, order="F"))


def tool_local(rho0, n, H, gamma, t_max, dt=0.01):
    """simulate_dynamic_lindblad, noise_type local (the scalar law, mutual-purity bridge), transcribed."""
    L = [site(Z, k, n) for k in range(n)]
    rho, rows = rho0.copy(), []
    for step in range(int(t_max / dt) + 1):
        if step % 10 == 0:
            rows.append((round(purity(rho), 6), round(mutual_purity(rho, n), 6), round(psi(rho), 6)))
        if step == int(t_max / dt):
            break
        g = gamma * mutual_purity(rho, n)
        d = -1j * (H @ rho - rho @ H)
        for Lk in L:
            d += g * (Lk @ rho @ Lk - rho)
        rho = rho + dt * d
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
    return rows


def claim(label, ok, detail=""):
    CLAIMS.append((label, ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    CHECKS.append((label, ok))
    print(f"  {label:<62s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")


def main():
    print("Sec. 2.3 / A.4, the C_int / C_ext comparison (Bell+, Heisenberg, t = 1):")
    pi, di = tool_cint(0.1, [0.1, 0.1])
    pe, de = tool_cint(0.1, [0.1, 0.0])
    claim("the logged gamma 0.1 call regenerates (purities, deltas, difference, six digits; four shown)",
          [round(pi, 6), round(di, 6), round(pe, 6), round(de, 6), round(di - de, 6)]
          == [LOG_CINT[k] for k in ("c_int purity", "c_int delta", "c_ext purity", "c_ext delta", "difference")],
          f"({pi:.6f}, {di:.6f}, {pe:.6f}, {de:.6f})")
    pi, di = tool_cint(0.05, [0.05, 0.05])
    pe, de = tool_cint(0.05, [0.05, 0.0])
    check("gamma 0.05, C_int purity (the tool)", 0.835, pi)
    check("gamma 0.05, C_ext purity (the tool)", 0.909, pe)
    rows = {0.05: (-0.074, -0.043, -0.031), 0.10: (-0.111, -0.074, -0.037), 0.20: (-0.124, -0.111, -0.014)}
    for g, (pdi, pde, pdd) in rows.items():
        _, di = tool_cint(g, [g, g])
        _, de = tool_cint(g, [g, 0.0])
        check(f"gamma {g:.2f}: delta_int (the tool)", pdi, di)
        check(f"gamma {g:.2f}: delta_ext (the tool)", pde, de)
        check(f"gamma {g:.2f}: delta_int - delta_ext (the tool)", pdd, di - de)
    print("  exactly (expm), the closed forms and the matched prediction:")
    for g in (0.05, 0.10, 0.20):
        p_int = exact_purity(bell(), heisenberg(2, False), [g, g], 1.0)
        p_ext = exact_purity(bell(), heisenberg(2, False), [g, 0.0], 1.0)
        ci = (np.exp(-8 * g) - np.exp(-4 * g)) / 2
        ce = (np.exp(-4 * g) - np.exp(-2 * g)) / 2
        dev = max(abs(p_int - prediction(bell(), 1.0, 2 * g) - ci), abs(p_ext - prediction(bell(), 1.0, g) - ce))
        claim(f"gamma {g:.2f}: delta_int = {ci:.6f}, delta_ext = {ce:.6f} in closed form (<= 100 eps)",
              dev <= 100 * EPS, f"(largest deviation {dev:.1e})")
        check(f"gamma {g:.2f}: delta_int - delta_ext exactly", {0.05: -0.031, 0.10: -0.036, 0.20: -0.013}[g], ci - ce)
        matched = max(abs(p_int - prediction(bell(), 1.0, 4 * g)), abs(p_ext - prediction(bell(), 1.0, 2 * g)))
        claim(f"gamma {g:.2f}: against a prediction at each run's own rate both deltas vanish (<= 100 eps)",
              matched <= 100 * EPS, f"(largest {matched:.1e})")

    print("  the sign change of delta_int - delta_ext at gamma t = ln(phi)/2:")
    phi = (1 + 5 ** 0.5) / 2
    grid = np.linspace(0.01, 3.0, 300)
    us = np.exp(-2 * grid)
    dd = (np.exp(-8 * grid) - np.exp(-4 * grid)) / 2 - (np.exp(-4 * grid) - np.exp(-2 * grid)) / 2
    fact = us * (us - 1) * (us * us + us - 1) / 2
    claim("  Delta delta = u(u-1)(u^2+u-1)/2 with u = exp(-2 gamma t), on a grid of gamma t in (0, 3] (<= 100 eps)",
          np.abs(dd - fact).max() <= 100 * EPS, f"(largest {np.abs(dd - fact).max():.1e})")
    zero = np.log(phi) / 2
    claim(f"  its only zero in (0, 3] is at gamma t = ln(phi)/2 = {zero:.4f}: |delta_int| > |delta_ext| below, < above",
          all((g < zero) == (abs((np.exp(-8 * g) - np.exp(-4 * g)) / 2) > abs((np.exp(-4 * g) - np.exp(-2 * g)) / 2))
              for g in grid if abs(g - zero) > 1e-3))
    print("  Bell+ and Bell- are degenerate eigenstates of the bond, and dephasing keeps their span (integers, exact):")
    Hb = np.rint(np.real(heisenberg(2, False))).astype(int)
    P = np.zeros((4, 4), dtype=int)
    P[0, 0] = P[3, 3] = 1
    P[0, 3] = P[3, 0] = 1  # 2 |Bell+><Bell+|
    M = np.zeros((4, 4), dtype=int)
    M[0, 0] = M[3, 3] = 1
    M[0, 3] = M[3, 0] = -1  # 2 |Bell-><Bell-|
    claim("  H (2 Bell+) = 2 Bell+ and H (2 Bell-) = 2 Bell-: the bond acts as 1 on both", (Hb @ P == P).all()
          and (Hb @ M == M).all() and not (Hb @ np.diag([0, 1, 0, 0]) == np.diag([0, 1, 0, 0])).all())
    ZA = np.rint(np.real(site(Z, 0, 2))).astype(int)
    claim("  and dephasing maps Bell+ to Bell- (Z_A (2 Bell+) Z_A == 2 Bell-), so the span is kept",
          (ZA @ P @ ZA == M).all())

    print("\nA.1, GHZ_3 on the ring, local noise (the scalar law), gamma 0.05, t = 10:")
    rows = tool_local(ghz(3), 3, heisenberg(3, True), 0.05, 10.0)
    check("purity at t = 10", 0.525, rows[-1][0])
    claim("the bridge C stays 1/2 at every recorded time (six digits)", all(c == 0.5 for _, c, _ in rows))
    check("Psi at t = 0", 0.143, rows[0][2])
    check("Psi at t = 10", 0.032, rows[-1][2])
    cpsi = [c * p for _, c, p in rows]
    check("C Psi at t = 0", 0.071, cpsi[0])
    check("C Psi at t = 10", 0.016, cpsi[-1])
    claim("C Psi stays below 1/4 at every recorded time", max(cpsi) < 0.25, f"(largest {max(cpsi):.4f})")
    for a, b in ((0, 1), (0, 2), (1, 2)):
        r = ghz(3)
        keep = [a, b]
        t = r.reshape([2] * 6)
        other = [k for k in range(3) if k not in keep][0]
        pair = np.trace(np.moveaxis(t, [other, 3 + other], [0, 1]), axis1=0, axis2=1).reshape(4, 4)
        claim(f"  the GHZ_3 pair ({a},{b}) carries no coherence (l1 == 0: |000><111| differs on the traced site)",
              psi(pair) == 0.0)
    print("\nThe Layer-3 GHZ table of the quarter-boundary roadmap (gamma 0.1 on every site, t = 1):")
    for n in range(2, 7):
        Hn = np.rint(np.real(heisenberg(n, True))).astype(int)
        e0 = np.zeros(2**n, dtype=int)
        e1 = np.zeros(2**n, dtype=int)
        e0[0] = 1
        e1[-1] = 1
        lam = int((Hn @ e0)[0])
        claim(f"  N = {n}: |0..0> and |1..1> are degenerate eigenstates of the ring (H e = {lam} e, integers)",
              (Hn @ e0 == lam * e0).all() and (Hn @ e1 == lam * e1).all())
    for n, page in ((3, -0.1244), (4, -0.1244)):
        rho0, Hn = ghz(n), heisenberg(n, True)
        rho = rho0.copy()
        for _ in range(100):
            d = -1j * (Hn @ rho - rho @ Hn)
            for k in range(n):
                Zk = site(Z, k, n)
                d += 0.1 * (Zk @ rho @ Zk - rho)
            rho = rho + 0.01 * d
        rho = 0.5 * (rho + rho.conj().T)
        check(f"  N = {n}: the tool's delta (Euler dt 0.01), the table's four-digit tie", page,
              purity(rho) - prediction(rho0, 1.0, 0.1 * n), digits=4)
    print("  delta(N) = (e^(-2a) - e^(-a))/2, a = 2 N gamma t,")
    aa = np.linspace(0.01, 5, 49901)
    amin = aa[np.argmin((np.exp(-2 * aa) - np.exp(-aa)) / 2)]
    check("  is most negative at a = ln 2 (grid argmin)", np.log(2), amin, digits=3)
    check("  so at N, gamma t = 0.1", 3.5, amin / 0.2, digits=1)
    check("  N = 3 exactly", -0.1238, (np.exp(-1.2) - np.exp(-0.6)) / 2, digits=4)
    check("  N = 4 exactly", -0.1237, (np.exp(-1.6) - np.exp(-0.8)) / 2, digits=4)

    print("\nSec. 7.1 / A.5, GHZ_N on the ring, gamma 0.05, t = 5 (closed forms); the feedback is blind:")
    for n in (3, 4, 5):
        O = site(X, 0, n) @ site(X, 1, n)
        deph = 0.5 * (np.diag(np.eye(2**n)[0]) + np.diag(np.eye(2**n)[-1]))
        claim(f"  N = {n}: <X0 X1> == 0 on GHZ_N and on its dephased image", np.trace(ghz(n) @ O) == 0
              and np.trace(deph @ O) == 0)
    for n, pd, pp in ((3, -0.087, 0.525), (4, -0.059, 0.509), (5, -0.038, 0.503)):
        check(f"N = {n}: delta_final", pd, (np.exp(-4 * n * 0.05 * 5) - np.exp(-2 * n * 0.05 * 5)) / 2)
        check(f"N = {n}: purity_final", pp, 0.5 + 0.5 * np.exp(-4 * n * 0.05 * 5))

    held = sum(ok for _, ok in CLAIMS)
    matched = sum(ok for _, ok in CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; "
          f"{held} of {len(CLAIMS)} regenerations, identities and claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
