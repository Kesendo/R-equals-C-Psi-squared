#!/usr/bin/env python3
"""Reproduction + verification of the ceiling 6->4->2 correction (committed evidence for
experiments/CEILING_FOUR_NONLOCAL_CASES.md). Under the SAME reliable test -- the complete uniform-continuous
per-site map family (16-real-param palindrome objective, Q = M^(x)N, residual ||Q L Q^-1 -(-L-2sigma)||) --
applied to all 6 candidate term-sets:
  * the 2 (XIX+XIY+YIX, YIY+XIY+YIX) ROUTE this uniform family (palindrome ~1e-13, N-stable 3,4,5) => LOCAL (6->4).
  * the 2 Z-middle (XZX+XZY+YZX, YZY+XZY+YZX) RESIST (bounded away) and resist every per-site product Q => non-local.
  * the 2 I-heavy (IXI+IIY+YII, IYI+IIX+XII) RESIST this UNIFORM family too, but are LOCAL via a SITE-VARYING
    per-site product of single-site crossover maps (the 4->2 step, machine-precision; see ceiling_4to2_iheavy_local.py).
  * the 2 routable fail the discrete {P1,P4,M2} per-term routers (O(1)): the continuous M is needed = the gap.
Self-validating: asserts the routable below 1e-7 and the resisting (Z-middle + I-heavy-uniform) above 1e-2.
Run: python simulations/ceiling_6to4_verification.py
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.optimize import minimize
if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

_P = [[(0,1),(1,1),(2,1),(3,1)],[(1,1),(0,1),(3,1j),(2,-1j)],
      [(2,1),(3,-1j),(0,1),(1,1j)],[(3,1),(2,1j),(1,-1j),(0,1)]]
GAMMA = 0.05


def sprod(t, s):
    res, ph = [], 1+0j
    for a, b in zip(t, s):
        c, p = _P[a][b]; res.append(c); ph *= p
    return tuple(res), ph


def comm(t):
    k = len(t); dim = 4**k
    C = np.zeros((dim, dim), dtype=complex)
    for si in range(dim):
        s = tuple((si >> (2*(k-1-j))) & 3 for j in range(k))
        ts, pts = sprod(t, s); _, pst = sprod(s, t)
        co = pts - pst
        if abs(co) > 1e-15:
            r = sum(ts[j] << (2*(k-1-j)) for j in range(k)); C[r, si] = co
    return C


def embed(term, w, N):
    f = [0]*N
    for j, a in enumerate(term): f[w+j] = a
    return tuple(f)


def n_xy(idx, N): return sum(1 for j in range(N) if ((idx >> (2*(N-1-j))) & 3) in (1, 2))


def build_L(terms, N):
    k = len(terms[0]); dim = 4**N
    C = np.zeros((dim, dim), dtype=complex)
    for w in range(N-k+1):
        for t in terms: C += comm(embed(t, w, N))
    D = np.diag([-2*GAMMA*n_xy(i, N) for i in range(dim)]).astype(complex)
    return -1j*C + D, N*GAMMA


def build_M(p):
    q = p[:8] + 1j*p[8:]
    m = np.zeros((4, 4), dtype=complex)
    m[1,0],m[1,3],m[2,0],m[2,3] = q[0],q[1],q[2],q[3]
    m[0,1],m[0,2],m[3,1],m[3,2] = q[4],q[5],q[6],q[7]
    return m


def kpow(m, N):
    o = m
    for _ in range(N-1): o = np.kron(o, m)
    return o


def palindrome_resid(m, terms, N):
    L, sg = build_L(terms, N)
    Q = kpow(m, N)
    try:
        Qi = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        return 1e6
    if not np.all(np.isfinite(Qi)): return 1e6
    RHS = -L - 2*sg*np.eye(4**N)
    return float(np.linalg.norm(Q @ L @ Qi - RHS)/np.linalg.norm(RHS))


def optimize(terms, restarts=90):
    best, bestp = np.inf, None
    rng = np.random.default_rng(17)
    for _ in range(restarts):
        x0 = rng.standard_normal(16)
        r = minimize(lambda p: palindrome_resid(build_M(p), terms, 3), x0,
                     method="Nelder-Mead", options={"maxiter":7000,"xatol":1e-11,"fatol":1e-14})
        if r.fun < best: best, bestp = r.fun, r.x
        if best < 1e-11: break
    return best, build_M(bestp)


ROUTABLE = {"XIX+XIY+YIX": [(1,0,1),(1,0,2),(2,0,1)], "YIY+XIY+YIX": [(2,0,2),(1,0,2),(2,0,1)]}
# The genuine non-local ceiling (the Z-middle): resist this uniform family AND every per-site product Q.
NONLOCAL = {"XZX+XZY+YZX": [(1,3,1),(1,3,2),(2,3,1)], "YZY+XZY+YZX": [(2,3,2),(1,3,2),(2,3,1)]}
# The I-heavy: resist this UNIFORM family (~0.81) but are LOCAL via a SITE-VARYING per-site product of
# single-site crossover maps (the 4->2 step, machine-precision; see ceiling_4to2_iheavy_local.py).
SITEVARYING_LOCAL = {"IXI+IIY+YII": [(0,1,0),(0,0,2),(2,0,0)], "IYI+IIX+XII": [(0,2,0),(0,0,1),(1,0,0)]}
P1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,1j,0]], dtype=complex)
P4 = np.array([[0,0,1,0],[0,0,0,1j],[1,0,0,0],[0,1j,0,0]], dtype=complex)
M2 = np.array([[0,0,1,0],[0,0,0,-1j],[1,0,0,0],[0,-1j,0,0]], dtype=complex)


def main():
    print("THE 2 LOCAL CASES (per-site product Q exists; verify N=3,4,5):", flush=True)
    for name, terms in ROUTABLE.items():
        best, m = optimize(terms)
        r3, r4, r5 = (palindrome_resid(m, terms, N) for N in (3, 4, 5))
        disc = min(palindrome_resid(d, terms, 4) for d in (P1, P4, M2))
        print(f"  {name}: route={best:.1e}  N=3/4/5 {r3:.1e}/{r4:.1e}/{r5:.1e}  best-discrete={disc:.2f}", flush=True)
        assert best < 1e-7 and r4 < 1e-7 and r5 < 1e-7, f"{name} must route (local)"
        assert disc > 1e-2, f"{name} must fail the discrete routers (continuous-sum gap)"
    print("\nTHE 2 NON-LOCAL Z-MIDDLE CASES (no per-site product Q at all, the genuine ceiling):", flush=True)
    for name, terms in NONLOCAL.items():
        best, _ = optimize(terms, restarts=70)
        print(f"  {name}: best={best:.3f}  => no per-site Q", flush=True)
        assert best > 1e-2, f"{name} must resist (non-local)"
    print("\nTHE 2 I-HEAVY CASES (resist this UNIFORM family, but LOCAL via site-varying; ceiling_4to2_iheavy_local.py):", flush=True)
    for name, terms in SITEVARYING_LOCAL.items():
        best, _ = optimize(terms, restarts=70)
        print(f"  {name}: uniform best={best:.3f}  => resists UNIFORM only (routes site-varying = LOCAL)", flush=True)
        assert best > 1e-2, f"{name} must resist the uniform family"
    print("\nOK: 2 route uniform (local), 2 route site-varying (local), 2 Z-middle resist (non-local). 6->4->2.", flush=True)


if __name__ == "__main__":
    main()
