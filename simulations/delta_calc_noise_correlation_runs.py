"""
The correlation bridge of Noise Robustness §3, in the retired tool's step and exactly
====================================================================================
experiments/NOISE_ROBUSTNESS.md §3 reports the correlation bridge C of Bell+ on the two-site
Heisenberg chain (J = 1, h = 0) under local sigma_z dephasing at gamma 0.05, the tool's scalar
law: the rate gamma * C(rho) in front of sigma_z on both sites, with
C = max(0, min(1, 2 (P_AB - P_A P_B))). The retired delta_calc tool's source is kept outside the
repo and is transcribed here, not imported: an Euler step of dt 0.01, then the Hermitian part,
the negative eigenvalues set to zero and the trace renormalized, recorded every tenth step and
rounded to six digits. The page's table holds C = 1 until t = 1.7 and 0.987 and 0.950 at t = 1.8
and 2.0; this script gives the tool's step and exact propagation of the same law (DOP853) at those
times. On this trajectory rho = (|00><00| + |11><11| + f (|00><11| + |11><00|))/2, so
C = min(1, 1/2 + f^2) and the step reduces to f <- f (1 - 4 gamma C dt).

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import sys

import numpy as np
from scipy.integrate import solve_ivp

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
GAMMA, DT = 0.05, 0.01
H = sum(np.kron(P, P) for P in (X, Y, Z))
LS = [np.kron(Z, I2), np.kron(I2, Z)]
CLAIMS, CHECKS = [], []


def purity(r):
    return float(np.real(np.trace(r @ r)))


def bridge(rho):
    t = rho.reshape(2, 2, 2, 2)
    pa = purity(np.einsum("ijkj->ik", t))
    pb = purity(np.einsum("ijil->jl", t))
    return max(0.0, min(1.0, 2 * (purity(rho) - pa * pb)))


def rhs(rho):
    d = -1j * (H @ rho - rho @ H)
    g = GAMMA * bridge(rho)
    for L in LS:
        d += g * (L @ rho @ L - rho)
    return d


def bell():
    v = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    return np.outer(v, v.conj())


def tool_records(t_max=2.0):
    rho, recs = bell(), {}
    for step in range(int(round(t_max / DT)) + 1):
        if step % 10 == 0:
            recs[round(step * DT, 4)] = round(bridge(rho), 6)
        rho = rho + DT * rhs(rho)
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
    return recs


def claim(label, ok, detail=""):
    CLAIMS.append((label, ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    CHECKS.append((label, ok))
    print(f"  {label:<52s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")


def main():
    recs = tool_records()
    print("The tool's step (dt 0.01, records every 0.1):")
    claim("C = 1 at every record up to t = 1.7", all(recs[t] == 1.0 for t in recs if t <= 1.7))
    check("C at t = 1.8", 0.986, recs[1.8])
    check("C at t = 2.0", 0.950, recs[2.0])
    f, scalar = 1.0, {}
    for step in range(201):
        if step % 10 == 0:
            scalar[round(step * DT, 4)] = round(min(1.0, 0.5 + f * f), 6)
        f *= 1 - 4 * GAMMA * min(1.0, 0.5 + f * f) * DT
    claim("the scalar recursion f <- f (1 - 4 gamma C dt) gives every record", scalar == recs)

    print("Exact propagation of the same law (DOP853):")
    sol = solve_ivp(lambda _, v: rhs(v.reshape(4, 4)).reshape(-1), (0, 2.0), bell().reshape(-1),
                    method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    c = lambda t: bridge(sol.sol(t).reshape(4, 4))
    claim("C = 1 up to t = 1.7", all(c(t) == 1.0 for t in np.arange(0, 1.7001, 0.1)))
    check("C at t = 1.8", 0.987, c(1.8))
    check("C at t = 2.0", 0.951, c(2.0))

    held = sum(ok for _, ok in CLAIMS)
    matched = sum(ok for _, ok in CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; "
          f"{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
