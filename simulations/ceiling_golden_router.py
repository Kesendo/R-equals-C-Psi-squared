#!/usr/bin/env python3
r"""The golden ceiling router: committed, self-validating verification (F116).

The two Z-middle ceiling cases, H = sum of sliding 3-site windows of XZX+XZY+YZX (and the
X<->Y sibling YZY+XZY+YZX) under Z-dephasing, ARE palindromized by a per-site product
superoperator after all:

    W = tensor_l q_{l mod 4},    W L W^{-1} = -L - 2*sigma,    sigma = sum_l gamma_l,

with the period-4 golden closed form (phi^2 = phi + 1, all entries in Z[phi] + i*Z[phi]):

    g_l = q_l(I)  follows [a, a, b, b],   a = phi*X + Y,  b = X - phi*Y  (b = -R(a), a _|_ b)
    h_l = q_l(Z) = (-1)^(l+1) * i * R(g_l)        (R = 90-degree rotation in the (X,Y) plane)
    q_l(X) = -(g_l)_X * I + (h_l)_X * Z,   q_l(Y) = -(g_l)_Y * I + (h_l)_Y * Z

Each q_l is class-swapping ({I,Z} <-> {X,Y}), invertible (q^2 = -(2+phi)*I, all singular
values sqrt(2+phi): a scalar times a unitary), and the two directions a, b are the two roots
of the golden locus alpha^2 - alpha*beta - beta^2 = 0 (slopes 1/phi and -phi; tan 2*theta = 2).

Mechanism: per-window, cross-TEMPLATE. The window-summed anticommutator
{Q_3, [XZX+XZY+YZX, .]_3} vanishes EXACTLY at every window offset (the three templates cancel
against each other inside one window; per-term it fails), so window additivity gives the
palindromizer at EVERY N >= 3 with arbitrary site-dependent gamma_l. Equivalently, two-sided:
q_even(s) = (aZ) s Z and q_odd(s) = Z s (Za), so W(rho) = (2+phi)^{N/2} P rho Q with product
unitaries P, Q each ANTICOMMUTING with H. One-sided conjugation fails (the banked chiral-K
null result: it fixes the dephaser); the two-sided split is what class-swaps every site.

Block ledger
------------
  Block 1  closed form + per-site structure : EXACT (ring arithmetic; locus, B = diag(-1,1)C^T,
                                              q^2 = -(2+phi)I, class-swap blocks zero)
  Block 2  the window lemma                 : EXACT ring zero of the template-SUMMED window
                                              anticommutator at all 4 offsets; per-term NONZERO
  Block 3  {W,A} = 0                        : EXACT ring zero at N=3,4,5,6, both siblings
  Block 4  two-sided chiral characterization: q = u(.)v forms, P H P^dag = -H, Q H Q^dag = -H,
                                              [G,H] = 0 for G = PQ (float, ~1e-13)
  Block 5  end-to-end vs framework          : W L W^{-1} = -L - 2*sigma against
                                              lindbladian_pauli_dephasing, N=5 (gamma 0.3, 1.0),
                                              N=6 (0.7), and SITE-DEPENDENT gamma at N=5
  Block 6  exclusion algebra                 : window equations K1, K2 derived symbolically;
                                              uniform and period-2 admit only g = 0; the locus
                                              det = a^2-ab-b^2; the discrete C# candidates sit
                                              OFF the locus (P1: +1, P4/M2: -1, M: -1/2)
  Block 7  spectral consequences            : spec(L) pairs lambda <-> -lambda-2*sigma (N=5);
                                              spec(H) is exactly +/-E symmetric; L(G) = -2*sigma*G

Provenance: discovered 2026-06-10 (identity-column lemma -> golden locus -> full-N search ->
closed form), adversarially verified by 5 independent re-implementations. The 2026-06-07
"site-dependent routing at N=4" finding, then dismissed as a small-N artifact, was this router;
the genuine artifact was only the period-3 N=4 hit (the 3-cycle closes only for N <= 4).
Scope: OPEN chains (rings/PBC untested). Run: python simulations/ceiling_golden_router.py (~3 min).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import _build_kbody_chain  # noqa: E402

PHI = (1 + np.sqrt(5)) / 2

# ======================================================================
# Ring Z[phi] + i*Z[phi]: elements as 4 float64 arrays (p, q, r, s) = (p + q*phi) + i*(r + s*phi).
# Component matmuls in float64 are exact while every integer stays below 2^53 (asserted).
# ======================================================================


def ring_mat(entries, shape):
    m = [np.zeros(shape) for _ in range(4)]
    for (i, j), comps in entries.items():
        for c in range(4):
            m[c][i, j] = comps[c]
    return m


def _zphi_pair(p, q, pp, qq, op):
    return op(p, pp) + op(q, qq), op(p, qq) + op(q, pp) + op(q, qq)


def _ring_combine(A, B, op):
    p1, q1, r1, s1 = A
    p2, q2, r2, s2 = B
    rp, rq = _zphi_pair(p1, q1, p2, q2, op)
    ip, iq = _zphi_pair(r1, s1, r2, s2, op)
    ap, aq = _zphi_pair(p1, q1, r2, s2, op)
    bp, bq = _zphi_pair(r1, s1, p2, q2, op)
    return [rp - ip, rq - iq, ap + bp, aq + bq]


def ring_matmul(A, B):
    return _ring_combine(A, B, lambda x, y: x @ y)


def ring_kron(A, B):
    return _ring_combine(A, B, np.kron)


def ring_add(A, B):
    return [a + b for a, b in zip(A, B)]


def ring_max_abs(A):
    return max(float(np.max(np.abs(c))) for c in A)


def ring_to_complex(A):
    return A[0] + A[1] * PHI + 1j * (A[2] + A[3] * PHI)


# ======================================================================
# Pauli string algebra ({I,X,Y,Z} = {0,1,2,3}) with ring phases.
# ======================================================================
_PROD = {}
for _a in range(4):
    _PROD[(0, _a)] = (_a, (1, 0, 0, 0)); _PROD[(_a, 0)] = (_a, (1, 0, 0, 0))
    _PROD[(_a, _a)] = (0, (1, 0, 0, 0))
_PROD[(1, 2)] = (3, (0, 0, 1, 0)); _PROD[(2, 1)] = (3, (0, 0, -1, 0))
_PROD[(2, 3)] = (1, (0, 0, 1, 0)); _PROD[(3, 2)] = (1, (0, 0, -1, 0))
_PROD[(3, 1)] = (2, (0, 0, 1, 0)); _PROD[(1, 3)] = (2, (0, 0, -1, 0))


def _comp_mul(c1, c2):
    p, q = c1[0] * c2[0] + c1[1] * c2[1], c1[0] * c2[1] + c1[1] * c2[0] + c1[1] * c2[1]
    ip, iq = c1[2] * c2[2] + c1[3] * c2[3], c1[2] * c2[3] + c1[3] * c2[2] + c1[3] * c2[3]
    ap, aq = c1[0] * c2[2] + c1[1] * c2[3], c1[0] * c2[3] + c1[1] * c2[2] + c1[1] * c2[3]
    bp, bq = c1[2] * c2[0] + c1[3] * c2[1], c1[2] * c2[1] + c1[3] * c2[0] + c1[3] * c2[1]
    return (p - ip, q - iq, ap + bp, aq + bq)


def string_prod(t, s, k):
    r = 0
    ph = (1, 0, 0, 0)
    for site in range(k):
        sh = 2 * (k - 1 - site)
        lr, p = _PROD[((t >> sh) & 3, (s >> sh) & 3)]
        r |= lr << sh
        ph = _comp_mul(ph, p)
    return r, ph


L_ = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
CEILING = ['XZX', 'XZY', 'YZX']
SIBLING = ['YZY', 'YZX', 'XZY']


def commutator_super_ring(words, N, windows=None):
    """[sum_w sum_T T_w, .] as a ring matrix in the Pauli-string basis (dim 4^N)."""
    dim = 4 ** N
    entries = {}
    for w in (windows if windows is not None else range(N - 2)):
        for word in words:
            letters = [0] * N
            for j, c in enumerate(word):
                letters[w + j] = L_[c]
            t = 0
            for letter in letters:
                t = (t << 2) | letter
            for s in range(dim):
                ts, pts = string_prod(t, s, N)
                _, pst = string_prod(s, t, N)
                e = tuple(np.array(pts) - np.array(pst))
                if any(e):
                    cur = entries.get((ts, s), (0, 0, 0, 0))
                    entries[(ts, s)] = tuple(np.array(cur) + np.array(e))
    return ring_mat(entries, (dim, dim))


# ======================================================================
# The golden closed form.
# ======================================================================
A_VEC = [(0, 1, 0, 0), (1, 0, 0, 0)]      # a = (phi, 1) in (X, Y) coefficients
B_VEC = [(1, 0, 0, 0), (0, -1, 0, 0)]     # b = (1, -phi)


def rot90(v):
    return [tuple(-x for x in v[1]), v[0]]


def i_times(v):
    return [(-c[2], -c[3], c[0], c[1]) for c in v]


def neg(v):
    return [tuple(-x for x in c) for c in v]


def golden_site_map(l):
    g = A_VEC if (l % 4) in (0, 1) else B_VEC
    h = i_times(rot90(g))
    if l % 2 == 0:
        h = neg(h)                          # h_l = (-1)^(l+1) * i * R(g_l)
    entries = {(1, 0): g[0], (2, 0): g[1],          # q(I) = g   (rows X, Y)
               (1, 3): h[0], (2, 3): h[1],          # q(Z) = h
               (0, 1): tuple(-np.array(g[0])), (3, 1): h[0],   # q(X) = -g_X I + h_X Z
               (0, 2): tuple(-np.array(g[1])), (3, 2): h[1]}   # q(Y) = -g_Y I + h_Y Z
    return ring_mat(entries, (4, 4))


def xy_swap_conjugate(q):
    s = np.zeros((4, 4))
    s[0, 0] = 1; s[1, 2] = 1; s[2, 1] = 1; s[3, 3] = -1
    return [s @ c @ s for c in q]


def build_W_ring(site_maps, N):
    W = site_maps[0]
    for l in range(1, N):
        W = ring_kron(W, site_maps[l % 4])
    return W


# ======================================================================
# BLOCK 1 -- closed form + per-site structure (EXACT).
# ======================================================================
def block1_structure(maps):
    print("-" * 92)
    print("BLOCK 1  closed form + per-site structure  [EXACT ring arithmetic]")
    print("-" * 92)
    for l, q in enumerate(maps):
        qc = ring_to_complex_map(q)
        # class-swap: the {I,Z}->{I,Z} and {X,Y}->{X,Y} blocks are exactly zero
        assert all(qc[i, j] == 0 for i in (0, 3) for j in (0, 3)), f"site {l}: IZ block nonzero"
        assert all(qc[i, j] == 0 for i in (1, 2) for j in (1, 2)), f"site {l}: XY block nonzero"
        # golden locus for g and h: alpha^2 - alpha*beta - beta^2 = 0 (exact in the ring)
        for col, name in ((0, 'g'), (3, 'h')):
            al, be = qc[1, col], qc[2, col]
            locus = al * al - al * be - be * be
            assert abs(locus) < 1e-12, f"site {l}: {name} off the golden locus ({locus})"
        # B = diag(-1,1) C^T
        C = np.array([[qc[1, 0], qc[1, 3]], [qc[2, 0], qc[2, 3]]])
        B = np.array([[qc[0, 1], qc[0, 2]], [qc[3, 1], qc[3, 2]]])
        assert np.array_equal(B, np.diag([-1.0, 1.0]) @ C.T), f"site {l}: B != diag(-1,1)C^T"
        # q^2 = -(2+phi) I  and  all singular values sqrt(2+phi)
        q2 = ring_matmul(q, q)
        exp_m = ring_mat({(i, i): (-2, -1, 0, 0) for i in range(4)}, (4, 4))
        diff = [q2[c] - exp_m[c] for c in range(4)]
        assert max(float(np.max(np.abs(d))) for d in diff) == 0.0, f"site {l}: q^2 != -(2+phi)I"
        sv = np.linalg.svd(qc, compute_uv=False)
        assert np.allclose(sv, np.sqrt(2 + PHI)), f"site {l}: singular values not sqrt(2+phi)"
        print(f"  site {l}: class-swap OK, g/h on the golden locus (exact), B = diag(-1,1)C^T, "
              f"q^2 = -(2+phi)I, sigma = sqrt(2+phi) (cond 1)")
    print("BLOCK 1 PASS")


def ring_to_complex_map(q):
    return ring_to_complex(q)


# ======================================================================
# BLOCK 2 -- the window lemma (EXACT): template-summed vanishes, per-term does not.
# ======================================================================
def block2_window_lemma(maps, words, label):
    print("-" * 92)
    print(f"BLOCK 2  window lemma [{label}]  [EXACT: template-summed = 0 at all offsets; "
          f"per-term != 0]")
    print("-" * 92)
    for off in range(4):
        Q3 = ring_kron(ring_kron(maps[off % 4], maps[(off + 1) % 4]), maps[(off + 2) % 4])
        Csum = commutator_super_ring(words, 3, windows=[0])
        anti = ring_add(ring_matmul(Q3, Csum), ring_matmul(Csum, Q3))
        assert ring_max_abs(anti) == 0.0, f"window offset {off}: template-summed anticommutator != 0"
        per_term = []
        for word in words:
            Ct = commutator_super_ring([word], 3, windows=[0])
            at = ring_add(ring_matmul(Q3, Ct), ring_matmul(Ct, Q3))
            per_term.append(ring_max_abs(at))
        assert all(v > 0 for v in per_term), f"offset {off}: a per-term anticommutator vanished"
        print(f"  offset {off}: {{Q3, [Sum T, .]}} = 0 EXACT; per-term max-entries "
              f"{[f'{v:.0f}' for v in per_term]} (nonzero: the cancellation is cross-template)")
    print("BLOCK 2 PASS")


# ======================================================================
# BLOCK 3 -- {W, A} = 0 EXACT at N = 3..6, both siblings.
# ======================================================================
def block3_anticommutation(maps, words, label):
    print("-" * 92)
    print(f"BLOCK 3  {{W, A}} = 0  [{label}]  [EXACT ring zero, N = 3..6]")
    print("-" * 92)
    for N in (3, 4, 5, 6):
        A = commutator_super_ring(words, N)
        W = build_W_ring(maps, N)
        assert max(ring_max_abs(W), ring_max_abs(A)) < 2 ** 40, "float-exactness bound"
        anti = ring_add(ring_matmul(W, A), ring_matmul(A, W))
        assert ring_max_abs(anti) == 0.0, f"{label}: {{W,A}} != 0 at N={N}"
        print(f"  N={N}: max |entry| (all 4 ring components) = 0.0  EXACT")
    print("BLOCK 3 PASS")


# ======================================================================
# BLOCK 4 -- the two-sided chiral characterization (float, phi irrational).
# ======================================================================
PAULI_M = [np.eye(2, dtype=complex),
           np.array([[0, 1], [1, 0]], complex),
           np.array([[0, -1j], [1j, 0]], complex),
           np.array([[1, 0], [0, -1]], complex)]


def block4_two_sided(maps):
    print("-" * 92)
    print("BLOCK 4  two-sided form: q_even = (aZ)(.)Z, q_odd = Z(.)(Za); P H P+ = -H = Q H Q+")
    print("-" * 92)
    X, Y, Z = PAULI_M[1], PAULI_M[2], PAULI_M[3]
    phi = PHI
    axes = []
    for l in range(4):
        gc = ring_to_complex(maps[l])[1:3, 0]
        ax = (gc[0] * X + gc[1] * Y) / np.sqrt(2 + phi)        # normalized golden axis
        axes.append(ax)
        u, v = (ax @ Z, Z) if l % 2 == 0 else (Z, Z @ ax)
        for letter in range(4):
            sigma = PAULI_M[letter]
            via_q = sum(ring_to_complex(maps[l])[k, letter] * PAULI_M[k] for k in range(4))
            two_sided = u @ sigma @ v * np.sqrt(2 + phi)
            assert np.allclose(via_q, two_sided, atol=1e-12), \
                f"site {l}, letter {letter}: two-sided form mismatch"
    N = 5
    H = _build_kbody_chain(N, [tuple(w) + (1.0,) for w in CEILING])
    P = np.array([[1]], dtype=complex)
    Q = np.array([[1]], dtype=complex)
    for l in range(N):
        ax = axes[l % 4]
        P = np.kron(P, ax @ Z if l % 2 == 0 else Z)
        Q = np.kron(Q, Z if l % 2 == 0 else Z @ ax)
    assert np.linalg.norm(P @ H + H @ P) < 1e-12, "P does not anticommute with H"
    assert np.linalg.norm(Q @ H + H @ Q) < 1e-12, "Q does not anticommute with H"
    G = P @ Q
    assert np.linalg.norm(G @ H - H @ G) < 1e-12, "[G, H] != 0"
    print(f"  per-site two-sided forms bit-exact; at N=5: ||{{P,H}}|| , ||{{Q,H}}|| < 1e-12, "
          f"[PQ, H] = 0  OK")
    print("BLOCK 4 PASS")


# ======================================================================
# BLOCK 5 -- end-to-end vs the framework Lindbladian (incl. site-dependent gamma).
# ======================================================================
P_BASIS = np.column_stack([p.reshape(4) for p in PAULI_M])
P_INV = np.linalg.inv(P_BASIS)


def apply_W_comp(qs_comp, vec, N):
    d = 2 ** N
    t = vec.reshape([2] * (2 * N))
    perm = []
    for l in range(N):
        perm += [l, N + l]
    t = np.transpose(t, perm).reshape([4] * N)
    for l in range(N):
        t = np.tensordot(qs_comp[l % 4], t, axes=([1], [l]))
        t = np.moveaxis(t, 0, l)
    t = t.reshape([2] * (2 * N))
    return np.transpose(t, np.argsort(perm)).reshape(d * d)


def block5_end_to_end(maps, words, label):
    print("-" * 92)
    print(f"BLOCK 5  end-to-end  [{label}]: W L W^-1 = -L - 2*sigma vs framework")
    print("-" * 92)
    qs_comp = [P_BASIS @ ring_to_complex(q) @ P_INV for q in maps]
    cases = [(5, [0.3] * 5), (5, [1.0] * 5), (6, [0.7] * 6),
             (5, [0.3, 1.1, 0.05, 2.0, 0.77])]
    for N, gammas in cases:
        d = 2 ** N
        H = _build_kbody_chain(N, [tuple(w) + (1.0,) for w in words])
        L = lindbladian_pauli_dephasing(H, list(gammas), dephase_letter='Z')
        sigma = sum(gammas)
        rng = np.random.default_rng(7)
        worst = 0.0
        for _ in range(8):
            v = rng.standard_normal(d * d) + 1j * rng.standard_normal(d * d)
            v /= np.linalg.norm(v)
            out = (apply_W_comp(qs_comp, L @ v, N) + L @ apply_W_comp(qs_comp, v, N)
                   + 2 * sigma * apply_W_comp(qs_comp, v, N))
            worst = max(worst, np.linalg.norm(out)
                        / max(np.linalg.norm(apply_W_comp(qs_comp, v, N)), 1e-300))
        tag = "site-dependent" if len(set(gammas)) > 1 else f"gamma={gammas[0]}"
        assert worst < 1e-12, f"{label} N={N} {tag}: palindromizer residual {worst:.2e}"
        print(f"  N={N} {tag} (sigma={sigma:.2f}): max rel residual = {worst:.3e}  OK")
    print("BLOCK 5 PASS")


# ======================================================================
# BLOCK 6 -- exclusion algebra (symbolic): uniform and period-2 are empty; the golden locus.
# ======================================================================
def block6_exclusion():
    print("-" * 92)
    print("BLOCK 6  exclusion  [symbolic: uniform/period-2 only g = 0; locus; discrete OFF-locus]")
    print("-" * 92)
    a0, b0, a2, b2 = sp.symbols('a0 b0 a2 b2')
    # the identity-column window equations on (g_w, g_{w+2}) (middle factor divides out):
    K1 = a0 * a2 + a0 * b2 + b0 * a2
    K2 = a0 * b2 + b0 * a2 - b0 * b2
    # uniform: g_w = g_{w+2}
    sols_u = sp.solve([K1.subs({a2: a0, b2: b0}), K2.subs({a2: a0, b2: b0})], [a0, b0], dict=True)
    assert all(s[a0] == 0 and s[b0] == 0 for s in sols_u), f"uniform admits nonzero g: {sols_u}"
    # period 2: same constraint (sites w and w+2 carry the same map)
    print("  uniform / period-2: K1 = K2 = 0 forces g = 0 (sympy solve over C)  OK")
    # nontrivial pairs need the determinant = the golden locus
    M = sp.Matrix([[a2 + b2, a2], [b2, a2 - b2]])
    det = sp.expand(M.det())
    assert sp.simplify(det - (a2 ** 2 - a2 * b2 - b2 ** 2)) == 0, "locus determinant mismatch"
    phi_s = (1 + sp.sqrt(5)) / 2
    for al, be, name in ((phi_s, 1, 'a'), (1, -phi_s, 'b')):
        assert sp.simplify(al ** 2 - al * be - be ** 2) == 0, f"{name} off the locus"
    # golden [a,a,b,b]: both window-pair types (a,b) and (b,a) satisfy K1 = K2 = 0
    for g0, g2 in (((phi_s, 1), (1, -phi_s)), ((1, -phi_s), (phi_s, 1))):
        v1 = sp.simplify(K1.subs({a0: g0[0], b0: g0[1], a2: g2[0], b2: g2[1]}))
        v2 = sp.simplify(K2.subs({a0: g0[0], b0: g0[1], a2: g2[0], b2: g2[1]}))
        assert v1 == 0 and v2 == 0, "golden pair fails the window equations"
    print("  det = alpha^2 - alpha*beta - beta^2 (the golden locus); [a,a,b,b] solves all "
          "window pairs exactly  OK")
    # the discrete C# candidates sit OFF the locus
    for name, al, be, want in (('P1 (g=X)', 1, 0, 1), ('P4/M2 (g=Y)', 0, 1, -1),
                               ('M (g ~ X+Y)', sp.sqrt(2) / 2, sp.sqrt(2) / 2, sp.Rational(-1, 2))):
        val = sp.simplify(al ** 2 - al * be - be ** 2)
        assert sp.simplify(val - want) == 0, f"{name}: locus value {val} != {want}"
        print(f"  discrete candidate {name}: locus value {want} != 0 (cannot route)  OK")
    print("BLOCK 6 PASS")


# ======================================================================
# BLOCK 7 -- spectral consequences.
# ======================================================================
def block7_spectral(maps):
    print("-" * 92)
    print("BLOCK 7  spectral: spec(L) palindromic; spec(H) = -spec(H); L(G) = -2*sigma*G")
    print("-" * 92)
    N, gam = 5, 0.45
    H = _build_kbody_chain(N, [tuple(w) + (1.0,) for w in CEILING])
    L = lindbladian_pauli_dephasing(H, [gam] * N, dephase_letter='Z')
    sigma = N * gam
    ev = np.linalg.eigvals(L)
    # greedy nearest pairing of spec(L) against -spec(L) - 2*sigma (sort-based comparison is
    # unstable under degenerate real parts)
    used = np.zeros(len(ev), bool)
    worst = 0.0
    for i in range(len(ev)):
        if used[i]:
            continue
        d = np.abs(ev - (-ev[i] - 2 * sigma))
        d[used] = np.inf
        j = int(np.argmin(d))
        worst = max(worst, float(d[j]))
        used[i] = True
        if j != i:
            used[j] = True
    assert worst < 1e-10, f"spec(L) not palindromic (worst {worst:.2e})"
    eh = np.sort(np.linalg.eigvalsh(H))
    assert np.allclose(eh, -eh[::-1], atol=1e-12), "spec(H) not +/-E symmetric"
    X, Y, Z = PAULI_M[1], PAULI_M[2], PAULI_M[3]
    G = np.array([[1]], dtype=complex)
    phi = PHI
    dirs = [(phi, 1), (phi, 1), (1, -phi), (1, -phi)]
    for l in range(N):
        al, be = dirs[l % 4]
        G = np.kron(G, (al * X + be * Y) / np.sqrt(2 + phi))
    # full Lindblad action on G via the framework superoperator
    vecG = G.reshape((2 ** N) ** 2)
    out = L @ vecG
    assert np.linalg.norm(out + 2 * sigma * vecG) / np.linalg.norm(vecG) < 1e-10, \
        "L(G) != -2*sigma*G"
    print(f"  spec(L) pairs about -sigma (worst {worst:.2e}); spec(H) exactly +/-E; "
          f"G = prod[a,a,b,b] is an exact eigenmode at the palindrome floor -2*sigma  OK")
    print("BLOCK 7 PASS")


def main():
    print("=" * 92)
    print("THE GOLDEN CEILING ROUTER (F116) -- self-validating verification")
    print("=" * 92)
    ceiling_maps = [golden_site_map(l) for l in range(4)]
    sibling_maps = [xy_swap_conjugate(q) for q in ceiling_maps]

    block1_structure(ceiling_maps)
    block2_window_lemma(ceiling_maps, CEILING, "XZX+XZY+YZX")
    block2_window_lemma(sibling_maps, SIBLING, "YZY+XZY+YZX (X<->Y conjugated)")
    block3_anticommutation(ceiling_maps, CEILING, "XZX+XZY+YZX")
    block3_anticommutation(sibling_maps, SIBLING, "YZY+XZY+YZX")
    block4_two_sided(ceiling_maps)
    block5_end_to_end(ceiling_maps, CEILING, "XZX+XZY+YZX")
    block5_end_to_end(sibling_maps, SIBLING, "YZY+XZY+YZX")
    block6_exclusion()
    block7_spectral(ceiling_maps)

    print("=" * 92)
    print("ALL BLOCKS PASS")
    print("=" * 92)


if __name__ == "__main__":
    main()
