"""Holonomies of the complement connection for Heisenberg bonds under dephasing axes off the letters.

The second route of docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md, last section. Every site is turned
so its axis n_l is z; a basis state is then a product of spinors along +-n_l, and Theorem 2 of
PROOF_PALINDROME_COMPLEMENT_CONNECTION reads the far end through the ratios H_xy / H_{xbar,ybar}. Four loops:

1. the triangle of pair flips: ((d + i det)/(d - i det))^2, d = 1 + a.b + b.c + c.a, (a, b, c) = (s0 n0, -s1 n1,
   -s2 n2), s_l = (-1)^{x_l}; exact (sympy) at four triples and all eight vertices;
2. the chain's two pair flips: e^{-+4i phi}, phi the angle at n1 between the planes (n0, n1) and (n1, n2), read as
   the area of the lune n0, n1, n2, -n1; floating point, random triples;
3. the chain's two single flips, at site 2 and then site 1: e^{2i arg(u_perp.v_perp + i a.(u x v))}, a = alpha_1,
   u = J01 alpha_0 - J12 alpha_2, v = J01 alpha_0 + J12 alpha_2, alpha_l = s_l n_l; floating point, random
   triples and signed weights (loops 2 and 3 at random weights too);
4. two fields at the sites 0 and 2, all axes along z, at weights 0.5 and -2, the loop 000 -> 100 -> 010 -> 001 -> 000 (the field at 0,
   two hops, the field at 2): e^{2i(theta_2 - theta_0)}; floating point.
The triangle loop is x -> x^e01 -> x^e02 -> x; the reverse orientation gives the inverse.

Run: python simulations/heisenberg_colour_holonomies.py (about a minute; the exact part dominates).
"""
import itertools
import sys

import numpy as np
import sympy as sp

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)
I2 = np.eye(2)
PAULI = [X, Y, Z]


def site(n, l, m):
    r = np.array([[1]], complex)
    for k in range(n):
        r = np.kron(r, m if k == l else I2)
    return r


def rot_to_z(v, phi):
    """A proper rotation carrying the unit direction of v to z, followed by a turn phi about z (a gauge)."""
    v = np.array(v, float)
    v /= np.linalg.norm(v)
    z = np.array([0, 0, 1.0])
    c = v @ z
    ax = np.cross(v, z)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        r = np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
    else:
        k = ax / s
        kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        r = np.eye(3) + s * kx + (1 - c) * kx @ kx
    rz = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    return rz @ r


def turned_h(axes, edges, fields, phis):
    """H = sum J_ij sigma_i.sigma_j + sum h_l.sigma_l, edges (i, j, J_ij), in the frame where every axis is z."""
    n = len(axes)
    rs = [rot_to_z(a, p) for a, p in zip(axes, phis)]
    h = np.zeros((2 ** n, 2 ** n), complex)
    for i, j, w in edges:
        o = w * rs[i] @ rs[j].T
        for a in range(3):
            for b in range(3):
                h += o[a, b] * site(n, i, PAULI[a]) @ site(n, j, PAULI[b])
    for l, f in enumerate(fields):
        fp = rs[l] @ np.array(f, float)
        for a in range(3):
            h += fp[a] * site(n, l, PAULI[a])
    return h


def loop(h, n, path):
    m = 2 ** n - 1
    r = 1
    for u, v in zip(path, path[1:]):
        r *= h[u, v] / h[m ^ u, m ^ v]
    return r


def solid(a, b, c):
    return 2 * np.arctan2(np.linalg.det(np.array([a, b, c])), 1 + a @ b + b @ c + c @ a)


def bit(n, l):
    return 1 << (n - 1 - l)


def triangle_exact():
    """Loop 1 exactly: spinor bases from sympy eigenvectors, the holonomy against the closed form."""
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])
    kron = sp.kronecker_product

    def spinor(v, s):
        u = sp.Matrix(v) / sp.sqrt(sum(c ** 2 for c in v))
        m = u[0] * sx + u[1] * sy + u[2] * sz
        for val, _, vecs in m.eigenvects():
            if sp.simplify(val - s) == 0:
                w = vecs[0]
                return w / sp.sqrt(sp.simplify((w.H * w)[0]))
        raise ValueError(v)

    def heis(i, j):
        def st(l, m):
            r = sp.Matrix([[1]])
            for k in range(3):
                r = kron(r, m if k == l else sp.eye(2))
            return r
        return sum((st(i, p) * st(j, p) for p in (sx, sy, sz)), sp.zeros(8, 8))

    h = heis(0, 1) + heis(1, 2) + heis(0, 2)
    ok = True
    for axes in [[(1, 2, 0), (0, 1, 2), (2, 0, 1)], [(1, 1, 0), (0, 1, 1), (1, 0, 1)],
                 [(2, 2, 1), (1, 2, -1), (1, 0, 2)], [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]:
        def basis(x):
            v = sp.Matrix([[1]])
            for l in range(3):
                v = kron(v, spinor(axes[l], (-1) ** ((x >> (2 - l)) & 1)))
            return v

        def entry(x, y):
            return sp.simplify((basis(x).H * h * basis(y))[0])

        units = [sp.Matrix(a) / sp.sqrt(sum(c ** 2 for c in a)) for a in axes]
        for x in range(8):
            path = [x, x ^ bit(3, 0) ^ bit(3, 1), x ^ bit(3, 0) ^ bit(3, 2), x]
            num = sp.prod([entry(u, v) for u, v in zip(path, path[1:])])
            den = sp.prod([entry(7 ^ u, 7 ^ v) for u, v in zip(path, path[1:])])
            hol = sp.simplify(num / den)
            s = [(-1) ** ((x >> (2 - l)) & 1) for l in range(3)]
            a, b, c = s[0] * units[0], -s[1] * units[1], -s[2] * units[2]
            d = 1 + a.dot(b) + b.dot(c) + c.dot(a)
            det = a.dot(b.cross(c))
            pred = ((d + sp.I * det) / (d - sp.I * det)) ** 2
            same = sp.simplify(sp.expand_complex(hol - pred)) == 0
            ok &= same
        print(f"  triangle, axes {axes}: closed form met at all eight vertices so far: {ok}")
    return ok


def chain_lune(rng):
    worst = 0.0
    for _ in range(40):
        a = rng.normal(size=(3, 3))
        a /= np.linalg.norm(a, axis=1)[:, None]
        j01, j12 = rng.choice([-1, 1], 2) * rng.uniform(0.3, 3, 2)
        h = turned_h(list(a), [(0, 1, j01), (1, 2, j12)], [(0, 0, 0)] * 3, tuple(rng.uniform(0, 6, 3)))
        p, q = bit(3, 0) ^ bit(3, 1), bit(3, 1) ^ bit(3, 2)
        for x in range(8):
            hol = loop(h, 3, [x, x ^ p, x ^ p ^ q, x ^ q, x])
            lune = solid(a[0], a[1], a[2]) + solid(a[0], a[2], -a[1])
            k = -1 if (x >> 1) & 1 == 0 else 1
            worst = max(worst, abs(hol - np.exp(k * 2j * lune)))
            # the lune's area is twice the angle at n1 between the planes (n0, n1) and (n1, n2)
            m1, m2 = np.cross(a[1], a[0]), np.cross(a[1], a[2])
            phi = np.arctan2(np.cross(m1, m2) @ a[1], m1 @ m2)
            worst = max(worst, abs(np.exp(2j * lune) - np.exp(-4j * phi)))
    return worst


def chain_square(rng):
    worst = 0.0
    for _ in range(200):
        a = rng.normal(size=(3, 3))
        a /= np.linalg.norm(a, axis=1)[:, None]
        j01, j12 = rng.choice([-1, 1], 2) * rng.uniform(0.3, 3, 2)
        h = turned_h(list(a), [(0, 1, j01), (1, 2, j12)], [(0, 0, 0)] * 3, tuple(rng.uniform(0, 6, 3)))
        for x in range(8):
            s = np.array([(-1) ** ((x >> (2 - l)) & 1) for l in range(3)])
            al = a * s[:, None]
            av, u, v = al[1], j01 * al[0] - j12 * al[2], j01 * al[0] + j12 * al[2]
            up, vp = u - (u @ av) * av, v - (v @ av) * av
            zz = up @ vp + 1j * (av @ np.cross(u, v))
            hol = loop(h, 3, [x, x ^ bit(3, 2), x ^ bit(3, 2) ^ bit(3, 1), x ^ bit(3, 1), x])
            worst = max(worst, abs(hol - np.exp(2j * np.angle(zz))))
    return worst


def fields_angle():
    worst = 0.0
    for t0 in [0.0, 0.4, 1.1, 2.5]:
        for t in [0.0, 0.3, np.pi / 4, np.pi / 2, 2.0, np.pi]:
            h = turned_h([(0, 0, 1)] * 3, [(0, 1, 0.5), (1, 2, -2.0)],
                         [(0.8 * np.cos(t0), 0.8 * np.sin(t0), 0), (0, 0, 0), (1.7 * np.cos(t), 1.7 * np.sin(t), 0)],
                         (0, 0, 0))
            hol = loop(h, 3, [0b000, 0b100, 0b010, 0b001, 0b000])
            worst = max(worst, abs(hol - np.exp(2j * (t - t0))))
    return worst


def triangle_weights(rng):
    """Loop 1 in floating point at random signed weights: the weights cancel in the ratio."""
    worst = 0.0
    for _ in range(40):
        a = rng.normal(size=(3, 3))
        a /= np.linalg.norm(a, axis=1)[:, None]
        w = rng.choice([-1, 1], 3) * rng.uniform(0.3, 3, 3)
        h = turned_h(list(a), [(0, 1, w[0]), (1, 2, w[1]), (0, 2, w[2])], [(0, 0, 0)] * 3, tuple(rng.uniform(0, 6, 3)))
        for x in range(8):
            s = np.array([(-1) ** ((x >> (2 - l)) & 1) for l in range(3)])
            aa, bb, cc = s[0] * a[0], -s[1] * a[1], -s[2] * a[2]
            d = 1 + aa @ bb + bb @ cc + cc @ aa
            det = np.linalg.det(np.array([aa, bb, cc]))
            pred = ((d + 1j * det) / (d - 1j * det)) ** 2
            hol = loop(h, 3, [x, x ^ bit(3, 0) ^ bit(3, 1), x ^ bit(3, 0) ^ bit(3, 2), x])
            worst = max(worst, abs(hol - pred))
    return worst


def main():
    rng = np.random.default_rng(20261027)
    print("1. triangle of pair flips, exact:")
    exact = triangle_exact()
    tri_w = triangle_weights(rng)
    lune = chain_lune(rng)
    square = chain_square(rng)
    field = fields_angle()
    print(f"1'. triangle at random signed weights, 40 random triples: worst deviation {tri_w:.2e}")
    print(f"2. chain, two pair flips, 40 random triples: worst deviation from the lune {lune:.2e}")
    print(f"3. chain, two single flips, 200 random triples: worst deviation from the triple-product form {square:.2e}")
    print(f"4. two fields, all axes z, weights 0.5 and -2, 24 angle pairs: worst deviation from e^(2i dtheta) {field:.2e}")
    ok = exact and tri_w < 1e-9 and lune < 1e-9 and square < 1e-9 and field < 1e-9
    print("ALL HOLD" if ok else "A LOOP DEVIATES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
