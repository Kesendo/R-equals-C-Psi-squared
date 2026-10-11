"""
The shared field's odd moments: where the third moment is blind, the fifth sees the break
=========================================================================================
experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md, the global break. Z dephasing with one covariance kappa common to
every pair of sites (a_kk = gamma, a_kl = kappa) on the open XY chain H = J sum (XX + YY) and on the open
Heisenberg chain H = J sum (XX + YY + ZZ). A palindrome about the only candidate centre -sigma (sigma = N gamma, the
mean eigenvalue) makes every odd power sum of (lambda + sigma) vanish. The script computes the third and fifth
moments M_n = Tr((L + sigma)^n) as polynomials in kappa and J at gamma = 1, every coefficient an integer, from the
cell sums alone (no eigenvalue): with B_c = sigma - r_c on the cells c = (i, j) (r_c = 1/2 u^T a u) and
S = ad_H = H x 1 - 1 x H^T (real symmetric; L + sigma = -i S + B),

    M3 = Tr B^3 - 3 Tr(S^2 B),
    M5 = Tr B^5 - 5 Tr(S^2 B^3) - 5 sum_{c,c'} S_cc'^2 B_c' B_c^2 + 5 Tr(S^4 B)

(the terms odd in S vanish because H is real, which the direct-trace check below confirms). It checks both against
the direct trace at a rational point, forms the resultant in kappa of M3/kappa and M5, and for every positive real
root J0^2 of the resultant locates the common kappa exactly: the first subresultant gives it as a rational function
of J^2, and its sign against the anti-correlated end -1/(N - 1) is decided on an isolating interval of J0^2 by an
exact root count (no float decides anything; a vanishing resultant, an identically zero sign polynomial, a root on a
rational point, or a common factor of degree two that no interval isolates raises instead of passing). The cell sums are
checked against direct traces: exact sympy matrices at N = 3, sparse integer matrix powers at N = 2 and 4, both with a
Hamiltonian built from Pauli products rather than from the cell sums' bit flips; and the third moment against the page's
closed forms on both chains at every N. At N = 3 on the XY chain the two moments are read at the
anti-correlated end with J = gamma, where the third vanishes and the fifth is -11520 gamma^5, the number the page's
verifier meets by matrix powers. The claim it checks: for N = 2..7 on both chains, M3 and M5 have no
common zero with -1/(N - 1) <= kappa < 0 and J > 0, so at every coupling and every kappa in that range one of them
is nonzero and the spectrum does not pair about -sigma (for kappa > 0 the third moment alone is nonzero).

Run: python simulations/shared_field_odd_moments.py   (seconds)
Prints only; exits with status 1 if any check fails.
"""

import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import sympy as sym

kap, J, x = sym.symbols("kappa J x")
CLAIMS = []


def claim(label, ok):
    CLAIMS.append(bool(ok))
    print(f"  {label}  {'holds' if ok else 'FAILS'}")
    sys.stdout.flush()


def zvals(N):
    return np.array([[1 - 2 * ((i >> (N - 1 - k)) & 1) for k in range(N)] for i in range(2 ** N)], dtype=np.int64)


def chain_h1(N, jz):
    """the open chain at J = 1: XX + YY moves an antiparallel pair with amplitude 2, jz ZZ on the diagonal"""
    d = 2 ** N
    zv = zvals(N)
    rows, cols, vals = [], [], []
    diag = np.zeros(d, dtype=np.int64)
    for k in range(N - 1):
        for i in range(d):
            if zv[i, k] != zv[i, k + 1]:
                rows.append(i ^ (1 << (N - 1 - k)) ^ (1 << (N - 2 - k)))
                cols.append(i)
                vals.append(2)
            diag[i] += jz * zv[i, k] * zv[i, k + 1]
    return sp.coo_matrix((vals, (rows, cols)), shape=(d, d), dtype=np.int64).tocsr() + sp.diags(diag)


def moment_polys(N, jz):
    d = 2 ** N
    zv = zvals(N)
    U = zv[:, None, :] - zv[None, :, :]
    k_dis = (U != 0).sum(axis=2).reshape(-1)
    s = (U.sum(axis=2) // 2).reshape(-1)
    # r = 2 gamma k + 2 kappa (s^2 - k) at gamma = 1, so B = sigma - r = (N - 2k) - 2 kappa (s^2 - k)
    beta0, beta1 = N - 2 * k_dis, -2 * (s ** 2 - k_dis)
    H1 = chain_h1(N, jz)
    Id = sp.identity(d, dtype=np.int64, format="csr")
    S1 = (sp.kron(H1, Id) - sp.kron(Id, H1.T)).tocsr()
    S1sq = (S1 @ S1).tocsr()
    diag2 = S1sq.diagonal()
    diag4 = np.asarray(S1sq.multiply(S1sq).sum(axis=1)).reshape(-1)
    count, t2, t4 = defaultdict(int), defaultdict(int), defaultdict(int)
    for c in range(d * d):
        key = (int(beta0[c]), int(beta1[c]))
        count[key] += 1
        t2[key] += int(diag2[c])
        t4[key] += int(diag4[c])
    pair = defaultdict(int)
    S1c = S1.tocoo()
    for c, cp, v in zip(S1c.row, S1c.col, S1c.data):
        pair[((int(beta0[c]), int(beta1[c])), (int(beta0[cp]), int(beta1[cp])))] += int(v) ** 2
    Bk = lambda key: key[0] + kap * key[1]
    TrB3 = sum(n * Bk(k) ** 3 for k, n in count.items())
    TrB5 = sum(n * Bk(k) ** 5 for k, n in count.items())
    M3 = sym.expand(TrB3 - 3 * J ** 2 * sum(v * Bk(k) for k, v in t2.items()))
    cross = sum(n * Bk(kc) ** 2 * Bk(kcp) for (kc, kcp), n in pair.items())
    M5 = sym.expand(TrB5 - 5 * J ** 2 * (sum(v * Bk(k) ** 3 for k, v in t2.items()) + cross)
                    + 5 * J ** 4 * sum(v * Bk(k) for k, v in t4.items()))
    return M3, M5


def pauli_chain(N, jz, Jv):
    """J sum (XX + YY + jz ZZ) over the open chain from the Pauli matrices themselves, a route independent of
    chain_h1's bit flips (so the direct check below also checks the Hamiltonian)"""
    X = sym.Matrix([[0, 1], [1, 0]])
    Y = sym.Matrix([[0, -sym.I], [sym.I, 0]])
    Zm = sym.Matrix([[1, 0], [0, -1]])

    def at(ops):
        out = sym.Matrix([[1]])
        for k in range(N):
            out = sym.kronecker_product(out, ops.get(k, sym.eye(2)))
        return out

    d = 2 ** N
    H = sym.zeros(d, d)
    for k in range(N - 1):
        H += at({k: X, k + 1: X}) + at({k: Y, k + 1: Y}) + jz * at({k: Zm, k + 1: Zm})
    return H * Jv


def direct_moments(N, jz, kv, Jv):
    """Tr((L + sigma)^3) and Tr((L + sigma)^5) at rational kappa and J, exact over Q(i) via sympy matrices (N <= 3)"""
    d = 2 ** N
    zv = zvals(N)
    H = pauli_chain(N, jz, Jv)
    Id = sym.eye(d)
    S = sym.Matrix(sym.kronecker_product(H, Id) - sym.kronecker_product(Id, H.T))
    r = []
    for i in range(d):
        for j in range(d):
            u = [int(zv[i, k] - zv[j, k]) for k in range(N)]
            r.append(sym.Rational(1, 2) * sum(u[k] * u[l] * (1 if k == l else kv) for k in range(N) for l in range(N)))
    M = -sym.I * S + sym.diag(*[N - rr for rr in r])
    M2 = M * M
    return sym.expand((M2 * M).trace()), sym.expand((M2 * M2 * M).trace())


def integer_chain(N, jz):
    """the open chain at J = 1 as an integer matrix built from Pauli products (Y x Y = -(YR x YR), Y = i YR),
    independent of chain_h1's bit flips"""
    X = np.array([[0, 1], [1, 0]], dtype=np.int64)
    YR = np.array([[0, -1], [1, 0]], dtype=np.int64)
    Zm = np.array([[1, 0], [0, -1]], dtype=np.int64)

    def at(ops):
        out = np.ones((1, 1), dtype=np.int64)
        for k in range(N):
            out = np.kron(out, ops.get(k, np.eye(2, dtype=np.int64)))
        return out

    d = 2 ** N
    H = np.zeros((d, d), dtype=np.int64)
    for k in range(N - 1):
        H += at({k: X, k + 1: X}) - at({k: YR, k + 1: YR}) + jz * at({k: Zm, k + 1: Zm})
    return H


def direct_integer(N, jz):
    """Tr((L + sigma)^3) and Tr((L + sigma)^5) at kappa = -2/15, J = 7/10 (gamma = 1) by sparse integer matrix powers
    of M = 30 (L + sigma): the real part is the cell diagonal 30 (sigma - r), the imaginary part -21 ad_H1; returned
    divided by 30^3 and 30^5 as Fractions"""
    from fractions import Fraction
    d = 2 ** N
    zv = zvals(N)
    U = zv[:, None, :] - zv[None, :, :]
    k_dis = (U != 0).sum(axis=2).reshape(-1)
    s = (U.sum(axis=2) // 2).reshape(-1)
    # 30 (sigma - r) with r = 2k + 2 kappa (s^2 - k), kappa = -2/15: 30 N - 60 k + 8 (s^2 - k)
    R = sp.diags((30 * N - 60 * k_dis + 8 * (s ** 2 - k_dis)).astype(np.int64)).tocsr()
    H1 = sp.csr_matrix(integer_chain(N, jz))
    Id = sp.identity(d, dtype=np.int64, format="csr")
    Im = (-21 * (sp.kron(H1, Id) - sp.kron(Id, H1.T))).tocsr()
    rho = int(max(np.abs(R).sum(axis=1).max(), 1)) + int(np.abs(Im).sum(axis=1).max())
    assert d * d * rho ** 5 < 2 ** 62, "int64 headroom"
    R2, I2 = R @ R - Im @ Im, R @ Im + Im @ R
    R4, I4 = R2 @ R2 - I2 @ I2, R2 @ I2 + I2 @ R2
    tr3 = int((R2 @ R - I2 @ Im).diagonal().sum())
    tr5 = int((R4 @ R - I4 @ Im).diagonal().sum())
    return Fraction(tr3, 30 ** 3), Fraction(tr5, 30 ** 5)


def sign_at_root(poly, root_interval, f):
    """the sign of the polynomial poly(x) at the root of f in the open rational interval (a, b): exact, by refining the
    interval until poly has no root in it, then reading poly at a rational point inside"""
    a, b = root_interval
    pf = sym.Poly(f, x)
    pp = sym.Poly(poly, x)
    for _ in range(200):
        if pp.is_zero:
            raise ValueError("the polynomial whose sign is asked for is identically zero")
        if pp.count_roots(a, b) == 0:
            mid = (a + b) / 2
            v = pp.eval(mid)
            if v == 0:
                raise ValueError("a zero inside an interval free of roots")
            return 1 if v > 0 else -1
        mid = (a + b) / 2
        fa, fm = pf.eval(a), pf.eval(mid)
        if fm == 0:
            raise ValueError("the root sits on a rational point")
        a, b = (a, mid) if fa * fm < 0 else (mid, b)
    raise ValueError("no isolation reached: the polynomial shares the root (a common factor of degree two or more, "
                     "or a degenerate subresultant)")


def main():
    print("the shared field's third and fifth moments on two open chains, gamma = 1")
    for jz in (0, 1):
        M3, M5 = moment_polys(3, jz)
        kv, Jv = sym.Rational(-2, 15), sym.Rational(7, 10)
        d3, d5 = direct_moments(3, jz, kv, Jv)
        claim(f"N=3 {'XY' if jz == 0 else 'Heisenberg'}: the cell-sum polynomials equal the direct traces at kappa = -2/15, "
              f"J = 7/10 (third {M3.subs({kap: kv, J: Jv})}, fifth {M5.subs({kap: kv, J: Jv})})",
              M3.subs({kap: kv, J: Jv}) == d3 and M5.subs({kap: kv, J: Jv}) == d5)
    for jz in (0, 1):
        for N in (2, 4):
            M3, M5 = moment_polys(N, jz)
            kv, Jv = sym.Rational(-2, 15), sym.Rational(7, 10)
            d3, d5 = direct_integer(N, jz)
            claim(f"N={N} {'XY' if jz == 0 else 'Heisenberg'}: the cell-sum polynomials equal sparse integer matrix powers "
                  f"of a Pauli-built generator at kappa = -2/15, J = 7/10",
                  M3.subs({kap: kv, J: Jv}) == sym.Rational(d3.numerator, d3.denominator)
                  and M5.subs({kap: kv, J: Jv}) == sym.Rational(d5.numerator, d5.denominator))
    M3, M5 = moment_polys(3, 0)
    at_end = {kap: sym.Rational(-1, 2), J: 1}
    claim(f"N=3 XY at the anti-correlated end kappa = -1/2 with J = 1: the third moment is {M3.subs(at_end)}, the fifth "
          f"{M5.subs(at_end)} (the page's verifier meets the same number by matrix powers)",
          M3.subs(at_end) == 0 and M5.subs(at_end) == -11520)
    for jz, name in ((0, "XY"), (1, "Heisenberg")):
        for N in range(2, 8):
            M3, M5 = moment_polys(N, jz)
            if jz == 0:
                closed = -4 ** N * 4 * (N - 1) * kap * (2 * N * (N - 2) * kap ** 2 + 3 * N * kap + 3 * J ** 2)
            else:
                closed = -4 ** N * 4 * kap * (2 * N * (N - 1) * (N - 2) * kap ** 2 + 3 * N * (N - 1) * kap + 3 * J ** 2)
            claim(f"{name} N={N}: the third moment is the page's closed form", sym.expand(M3 - closed) == 0)
            P3 = sym.cancel(M3 / kap)
            sub = sym.subresultants(sym.Poly(P3, kap), sym.Poly(M5, kap))
            lin = [p for p in sub if p.degree(kap) == 1]
            res = sym.expand(sym.resultant(P3, M5, kap))
            if res == 0:
                raise ValueError(f"{name} N={N}: the resultant vanishes identically, so P3 and M5 share a factor")
            res_x = sym.factor_list(res.subs(J, sym.sqrt(x)))
            roots = []
            for fac, _ in res_x[1]:
                if sym.degree(fac, x) < 1:
                    continue
                for (a_, b_), _m in sym.Poly(fac, x).intervals():
                    if a_ < 0 < b_:
                        raise ValueError("an isolating interval straddles zero")
                    if b_ > 0 and a_ >= 0 and not (a_ == b_ == 0):
                        roots.append((fac, (sym.Rational(a_), sym.Rational(b_))))
            outside_all = True
            details = []
            for fac, rt in roots:
                # the common kappa: the linear subresultant s1(J) kappa + s0(J) = 0 at J^2 = rt
                L1 = sym.Poly(lin[-1].as_expr(), kap)
                s1, s0 = L1.coeff_monomial(kap), L1.coeff_monomial(1)
                s1x, s0x = sym.expand(s1.subs(J, sym.sqrt(x))), sym.expand(s0.subs(J, sym.sqrt(x)))
                a_, b_ = rt
                if a_ == b_:
                    raise ValueError("a rational root of the resultant: decide it by hand")
                # kappa_c + 1/(N - 1) = (s1 - (N - 1) s0) / ((N - 1) s1): the common root lies below the end when < 0
                num = sym.expand(s1x - (N - 1) * s0x)
                sn = sign_at_root(num, (a_, b_), fac)
                sd = sign_at_root(s1x, (a_, b_), fac)
                below = sn * sd < 0
                outside_all &= below
                details.append(f"J0^2 in ({float(a_):.6f}, {float(b_):.6f}): the common kappa lies "
                               f"{'below' if below else 'NOT below'} -1/(N-1)")
            tag = "; ".join(details) if details else "no common zero with kappa != 0 at any J > 0"
            claim(f"{name} N={N}: M3 = {sym.factor(M3)}; {tag}", outside_all)
            claim(f"{name} N={N}: -M3/kappa has only positive coefficients in kappa and J, so M3 < 0 for every kappa > 0",
                  all(c > 0 for c in sym.Poly(sym.expand(-M3 / kap), kap, J).coeffs()))
            if N > 2:
                q0 = sym.Poly(sym.cancel(-M3.subs(J, 0) / kap ** 2), kap)
                claim(f"{name} N={N}: at J = 0 the third moment is -kappa^2 times a polynomial positive on the PSD range",
                      q0.count_roots(sym.Rational(-1, N - 1), 1) == 0 and q0.eval(0) > 0)
    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
