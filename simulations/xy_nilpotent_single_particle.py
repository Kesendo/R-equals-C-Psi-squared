"""Find the real profiles that make the XY chain's single-particle matrix m nilpotent.

Producer for docs/proofs/PROOF_R90_FROZEN_DIVISOR.md Section 9.2. On the R90 locus
m = -iJh - 2 Delta (h the XY hopping, 2 on every bond; Delta = diag(gamma_l - gbar)).
The similarity S m S^-1, S = diag(i^a), and the scale -1/(2J) turn m into the real tridiagonal matrix
with diagonal d = delta / J, superdiagonal +1 and subdiagonal -1, whose characteristic
polynomial is the continuant p_k = (x - d_k) p_(k-1) + p_(k-2). On the locus d is odd
under the chain reflection, so the odd coefficients vanish and m is nilpotent iff the
floor(N/2) even coefficients do: a square system F(u) = 0 in the floor(N/2) free
entries u of d. Where it holds, the corner block's frozen root carries a Jordan block
of size 2N - 3 (the Clebsch-Gordan image of one N x N block).

What this script does:
  1. a Newton search from random real starts (floats), closed under u -> -u (F is
     even) and deduplicated;
  2. each root refined to 60 digits and certified by the Krawczyk test on a real box
     of radius 10^-30 (K(X) inside the interior of X: exactly one root in X, hence
     real), and the certified boxes checked pairwise disjoint;
  3. prints per N the number found and certified, and the root of smallest max |d|
     (the smallest gbar / J that keeps every rate non-negative).

A random-start search gives a LOWER bound on the number of real roots; from N = 10 on
the default run and a run with seed 2 find different sets, so neither is complete. The gate
simulations/r90_frozen_divisor_gate.py (G20) certifies one root per N = 3..12 from
the seeds this prints (to 17 significant digits), and reads the N = 6 integer root exactly.

Usage:
  python simulations/xy_nilpotent_single_particle.py [Nmax] [starts] [seed]
      N = 3..Nmax (default 12), starts per N (default 300), RNG seed (default 1)
"""
import sys

import mpmath as mp
import numpy as np
from numpy.polynomial import polynomial as P


def charpoly(d):
    p0, p1 = np.array([1.0]), np.array([-d[0], 1.0])
    for k in range(1, len(d)):
        p0, p1 = p1, P.polyadd(P.polymul(np.array([-d[k], 1.0]), p1), p0)
    return p1


def full(u, N):
    h = N // 2
    d = np.zeros(N)
    d[:h] = u
    d[N - h:] = -u[::-1]
    return d


def F_float(u, N):
    c = charpoly(full(u, N))
    return np.array([c[N - 2 * k] for k in range(1, N // 2 + 1)])


def newton(u, N):
    for _ in range(200):
        f = F_float(u, N)
        if not np.all(np.isfinite(f)) or np.max(abs(u)) > 1e3:
            return None
        if np.max(abs(f)) < 1e-12:
            return u
        e = 1e-7
        Jm = np.array([(F_float(u + e * np.eye(len(u))[j], N) - f) / e for j in range(len(u))]).T
        try:
            step = np.linalg.solve(Jm, -f)
        except np.linalg.LinAlgError:
            return None
        lam = 1.0
        while lam > 1e-4 and np.max(abs(F_float(u + lam * step, N))) >= np.max(abs(f)):
            lam /= 2
        u = u + lam * step
    return None


def dual_F(u, N, one, zero):
    """Even continuant coefficients and their gradients (forward mode); the same
    routine as the gate's nm_dual_F."""
    h = N // 2
    d = []
    for a in range(N):
        if a < h:
            d.append((u[a], [one if j == a else zero for j in range(h)]))
        elif a >= N - h:
            b = N - 1 - a
            d.append((-u[b], [-one if j == b else zero for j in range(h)]))
        else:
            d.append((zero, [zero] * h))

    def times_x_minus(p, dk):
        out = [(zero, [zero] * h) for _ in range(len(p) + 1)]
        dv, dg = dk
        for i, (v, g) in enumerate(p):
            ov, og = out[i + 1]
            out[i + 1] = (ov + v, [a + b for a, b in zip(og, g)])
            ov, og = out[i]
            out[i] = (ov - dv * v, [a - (dv * b + v * c) for a, b, c in zip(og, g, dg)])
        return out

    p0 = [(one, [zero] * h)]
    p1 = times_x_minus(p0, d[0])
    for k in range(1, N):
        p2 = times_x_minus(p1, d[k])
        for i, (v, g) in enumerate(p0):
            ov, og = p2[i]
            p2[i] = (ov + v, [a + b for a, b in zip(og, g)])
        p0, p1 = p1, p2
    rows = [p1[N - 2 * k] for k in range(1, h + 1)]
    return [r[0] for r in rows], [r[1] for r in rows]


def refine(u, N):
    x = mp.matrix([mp.mpf(t) for t in u])
    for _ in range(80):
        F, Jm = dual_F(list(x), N, mp.mpf(1), mp.mpf(0))
        dx = mp.lu_solve(mp.matrix(Jm), -mp.matrix(F))
        x = x + dx
        if max(abs(t) for t in dx) < mp.mpf(10) ** -(mp.mp.dps - 5):
            break
    return [x[i] for i in range(N // 2)]


def krawczyk(y, N, r):
    h, iv = N // 2, mp.iv
    iv.dps = mp.mp.dps
    _, Jy = dual_F(y, N, mp.mpf(1), mp.mpf(0))
    Y = mp.matrix(Jy) ** -1
    X = [iv.mpf([t - r, t + r]) for t in y]
    yi = [iv.mpf(t) for t in y]
    Fy, _ = dual_F(yi, N, iv.mpf(1), iv.mpf(0))
    _, JX = dual_F(X, N, iv.mpf(1), iv.mpf(0))
    for i in range(h):
        s = yi[i] - sum(iv.mpf(Y[i, j]) * Fy[j] for j in range(h))
        for k in range(h):
            Mik = iv.mpf(1 if i == k else 0) - sum(iv.mpf(Y[i, j]) * JX[j][k] for j in range(h))
            s += Mik * (X[k] - yi[k])
        if not (s.a > X[i].a and s.b < X[i].b):
            return False
    return True


def main():
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    starts = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    mp.mp.dps = 60
    r = mp.mpf(10) ** -30
    for N in range(3, Nmax + 1):
        found = []
        for _ in range(starts):
            u = newton(rng.normal(scale=rng.uniform(0.5, 3), size=N // 2), N)
            # F is even (F(-u) = F(u), the mirrored profile) and u = 0 is never a root,
            # so the roots come in pairs: keep both members of every pair found
            for v in ((u, -u) if u is not None else ()):
                if all(np.max(abs(v - k)) > 1e-5 for k in found):
                    found.append(v)
        cert = [y for y in (refine(u, N) for u in found) if krawczyk(y, N, r)]
        disjoint = all(max(abs(a - b) for a, b in zip(p, q)) > 2 * r
                       for i, p in enumerate(cert) for q in cert[i + 1:])
        print(f"N={N}: found {len(found)}, certified {len(cert)}, boxes disjoint {disjoint}")
        if cert:
            best = min(cert, key=lambda y: max(abs(t) for t in y))
            print(f"  smallest max|d| = {mp.nstr(max(abs(t) for t in best), 8)} at u = "
                  f"[{', '.join(mp.nstr(t, 17) for t in best)}]", flush=True)


if __name__ == "__main__":
    main()
