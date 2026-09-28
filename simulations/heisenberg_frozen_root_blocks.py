"""Locate the Heisenberg profiles whose frozen root carries the largest Jordan block.

Producer for docs/proofs/PROOF_R90_FROZEN_DIVISOR.md Section 9.3. On the R90
locus at gbar != 0 the corner block Mtilde = L_(1,1) + 4 gbar has the root 0 with
multiplicity floor(N/2) at generic coupling; k extra eigenvalues arrive together
where the cofactor q(0) = c_m (m = floor(N/2)) and the next k - 1 coefficients of
det(eps - Mtilde) vanish at one coupling. With d1 = 1 (Mtilde is homogeneous in
(J, gamma)) the unknowns are w = J^2, G = gbar^2 and the floor(N/2) - 1 remaining
offset ratios: floor(N/2) + 1 of them. Setting c_m .. c_(2m) to zero is therefore a
square system, whose isolated solutions carry floor(N/2) + 1 arrivals, a block of
size floor(N/2) + 2 where the geometric count stays floor(N/2).

What this script does:
  1. the coefficients c_m .. c_(2m) symbolically over ZZ[z, g, d2, ...], z = -iJ
     (sympy DomainMatrix; N = 6 takes minutes, N = 7 about eleven hours);
  2. a Newton search from random real starts on the square system, in floats;
  3. discards the degenerate components (w ~ 0, gbar ~ 0, and the lines where two
     offsets coincide in size or one is 0 or 1, on which q(0) loses its constant
     term), refines what is left at 60 digits (deduplicated, failures counted), and prints each real carrier with
     its Jacobian determinant, the next coefficient c_(2m+1) there (nonzero means
     the block stops at floor(N/2) + 2), and its margin to the box of non-negative
     rates (gbar^2 - max offset^2).

It locates carriers; it does not read block sizes. The gate
simulations/r90_frozen_divisor_gate.py (G19) reads the chain at the carriers this
found, by a precision law, and certifies the smaller cases exactly.

Usage:
  python simulations/heisenberg_frozen_root_blocks.py N [starts] [seed]
      N = 4..7; starts per run (default 4000); the coefficient tables are cached as
      simulations/results/heisenberg_frozen_root_blocks/coeffs_N{N}.json (N = 6, 7
      committed, so the symbolic step is not repeated)
  python simulations/heisenberg_frozen_root_blocks.py N refine w G s2 [s3 ...]
      refines one carrier from a start (the 16 digits printed by a search) to 170
      digits and prints it to 140: this is how the gate's stored points were made
"""
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results",
                   "heisenberg_frozen_root_blocks")


def se_h(N):
    """Single-excitation Heisenberg h, Pauli convention: hop 2, ZZ diagonal."""
    h = [[0] * N for _ in range(N)]
    for a in range(N - 1):
        h[a][a + 1] = h[a + 1][a] = 2
    for a in range(N):
        h[a][a] = sum(-1 if (a == b or a == b + 1) else 1 for b in range(N - 1))
    return h


def coefficients(N):
    """c_m .. c_(2m+1) of det(eps - Mtilde) over ZZ[z, g, d2, ...], d1 = 1."""
    m = N // 2
    z, g = sp.symbols('z g')
    ds = [sp.Integer(1)] + list(sp.symbols(f'd2:{m + 1}'))
    gl = [g + d for d in ds] + ([g] if N % 2 else []) + [g - d for d in reversed(ds)]
    h, n = se_h(N), N * N
    rows = [[sp.Integer(0)] * n for _ in range(n)]
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                if h[a][c]:
                    rows[r][c * N + b] += z * h[a][c]
                if h[c][b]:
                    rows[r][a * N + c] -= z * h[c][b]
            if a != b:
                rows[r][r] += -2 * (gl[a] + gl[b])
            rows[r][r] += 4 * g
    gens = [z, g] + ds[1:]
    dom = sp.ZZ[tuple(gens)]
    cp = DomainMatrix.from_list_sympy(n, n, rows).convert_to(dom).charpoly()[::-1]
    assert all(cp[i] == dom.zero for i in range(m)), "root multiplicity below floor(N/2)"
    out = []
    for k in range(m, 2 * m + 2):
        P = sp.Poly(dom.to_sympy(cp[k]), *gens)
        mz = min(e[0] for e in P.monoms())
        mg = min(e[1] for e in P.monoms())
        terms = {}
        for e, c in P.terms():   # strip the monomial; z^2 = -w, g^2 = G
            a_, b_ = e[0] - mz, e[1] - mg
            assert a_ % 2 == 0 and b_ % 2 == 0
            key = (a_ // 2, b_ // 2) + tuple(e[2:])
            terms[key] = terms.get(key, 0) + int(c) * (-1) ** (a_ // 2)
        out.append(terms)
    return out


class Poly:
    """One coefficient as a dense term table, evaluated in numpy and in mpmath."""

    def __init__(self, terms):
        self.E = np.array(list(terms.keys()), dtype=np.int64)
        self.C = [int(v) for v in terms.values()]
        self.scale = max(abs(c) for c in self.C)
        self.Cf = np.array([c / self.scale for c in self.C])

    def num(self, T):
        out = np.ones((T[0].shape[0], len(self.Cf)))
        for j, Tj in enumerate(T):
            out *= Tj[:, self.E[:, j]]
        return out @ self.Cf

    def grad(self, T):
        cols = []
        for j in range(len(T)):
            E2 = self.E.copy()
            c = self.Cf * E2[:, j]
            E2[:, j] = np.maximum(E2[:, j] - 1, 0)
            out = np.ones((T[0].shape[0], len(c)))
            for k, Tk in enumerate(T):
                out *= Tk[:, E2[:, k]]
            cols.append(out @ c)
        return np.array(cols).T

    def mpv(self, x):
        pw = [[x[j] ** k for k in range(int(self.E[:, j].max()) + 1)] for j in range(len(x))]
        return mp.fsum(mp.mpf(c) * mp.fprod(pw[j][e[j]] for j in range(len(x)))
                       for c, e in zip(self.C, self.E)) / self.scale

    def mpg(self, x):
        pw = [[x[j] ** k for k in range(int(self.E[:, j].max()) + 1)] for j in range(len(x))]
        r = []
        for j in range(len(x)):
            acc = [mp.mpf(c) * int(e[j]) * mp.fprod(pw[k][e[k] - (1 if k == j else 0)]
                                                    for k in range(len(x)))
                   for c, e in zip(self.C, self.E) if e[j] > 0]
            r.append(mp.fsum(acc) / self.scale)
        return r


def tables(N):
    """The coefficient tables c_m .. c_(2m+1), from the cache or computed once."""
    os.makedirs(OUT, exist_ok=True)
    cache = os.path.join(OUT, f"coeffs_N{N}.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return [{tuple(k): v for k, v in t} for t in json.load(f)]
    t = time.time()
    terms = coefficients(N)
    with open(cache, "w") as f:
        json.dump([[[list(k), v] for k, v in t.items()] for t in terms], f)
    print(f"N={N}: coefficients in {time.time() - t:.0f} s", flush=True)
    return terms


def degenerate(x):
    """The components where q(0) loses its constant term or the profile collapses, and
    the solutions with no real coupling: w = J^2 below 10^-3 (J ~ 0, or J imaginary),
    gbar ~ 0, an offset ratio ~ 0 or of size ~ 1, two ratios of equal size."""
    offs = [1.0] + [abs(float(v)) for v in x[2:]]
    pairs = [abs(a - b) for i, a in enumerate(offs) for b in offs[i + 1:]]
    return float(x[0]) < 1e-3 or abs(float(x[1])) < 1e-5 or min(offs[1:] + pairs, default=1) < 1e-3


def refine(system, start, dps):
    """Newton at dps digits from start; returns the point or None."""
    mp.mp.dps = dps + 20
    try:
        r = mp.findroot(lambda *v: [p.mpv(v) for p in system], [mp.mpf(v) for v in start],
                        J=lambda *v: [p.mpg(v) for p in system],
                        tol=mp.mpf(10) ** -(2 * dps), maxsteps=200)
    except (ZeroDivisionError, ValueError):
        return None
    return [r[i] for i in range(len(start))]


def report(system, nxt, v, digits):
    res = max(abs(p.mpv(v)) for p in system)
    dj = mp.det(mp.matrix([p.mpg(v) for p in system]))
    bound = max([mp.mpf(1)] + [s ** 2 for s in v[2:]])
    print(f"  w={mp.nstr(v[0], digits)} G={mp.nstr(v[1], digits)} "
          f"s={[mp.nstr(s, digits) for s in v[2:]]} residual={mp.nstr(res, 3)} "
          f"detJac={mp.nstr(dj, 3)} next={mp.nstr(nxt.mpv(v), 4)} "
          f"margin={mp.nstr(v[1] - bound, 4)}", flush=True)


def search(N, starts, seed):
    polys = [Poly(t) for t in tables(N)]
    system, nxt = polys[:-1], polys[-1]
    nv = N // 2 + 1                                   # w, G, s2 .. s_m
    deg = max(int(p.E.max()) for p in polys) + 1
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.exp(rng.uniform(-3, 4.5, starts)), np.exp(rng.uniform(-3, 5, starts))]
                        + [rng.uniform(-10, 10, starts) for _ in range(nv - 2)])
    for _ in range(80):
        T = [X[:, j:j + 1] ** np.arange(deg)[None, :] for j in range(nv)]
        Fv = np.array([p.num(T) for p in system]).T
        Jv = np.array([p.grad(T) for p in system]).transpose(1, 0, 2)
        try:
            dx = np.linalg.solve(Jv, Fv[..., None])[..., 0]
        except np.linalg.LinAlgError:
            dx = np.array([np.linalg.lstsq(a, b, rcond=None)[0] for a, b in zip(Jv, Fv)])
        X = X - np.minimum(1, 5 / np.maximum(1e-300, np.abs(dx).max(1)))[:, None] * dx
        X[~np.isfinite(X)] = 0
    T = [X[:, j:j + 1] ** np.arange(deg)[None, :] for j in range(nv)]
    good = X[np.abs(np.array([p.num(T) for p in system]).T).max(1) < 1e-10]
    cands = [x for x in good if not degenerate(x)]
    print(f"N={N} seed {seed}: {len(good)}/{starts} converged in floats, "
          f"{len(cands)} off the degenerate components", flush=True)
    found, failed = [], 0
    for x in cands:
        v = refine(system, x, 60)
        if v is None or max(abs(p.mpv(v)) for p in system) > mp.mpf(10) ** -50:
            failed += 1
            continue
        if degenerate(v) or any(max(abs(a - b) for a, b in zip(v, u)) < mp.mpf(10) ** -30 for u in found):
            continue
        found.append(v)
        report(system, nxt, v, 16)
    print(f"N={N} seed {seed}: {len(found)} distinct carriers at 60 digits, {failed} starts "
          "that did not refine", flush=True)


def refine_point(N, start):
    """Refine one carrier to 170 digits and print it to 140: the stored points of the gate."""
    polys = [Poly(t) for t in tables(N)]
    v = refine(polys[:-1], start, 150)
    if v is None:
        print("did not refine")
        return
    report(polys[:-1], polys[-1], v, 140)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "refine":
        refine_point(int(sys.argv[1]), [float(x) for x in sys.argv[3:]])
    else:
        search(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 4000,
               int(sys.argv[3]) if len(sys.argv) > 3 else 1)
