"""The fifth centred trace under loss and gain: flux times the area of the rate triangle (F137).

For H = sum over bonds of J(XX + YY) + D_z(X_a Y_b - Y_a X_b) + Delta Z_a Z_b on any graph, loss gamma_down and gain
gamma_up per site, no fields, no dephasing, c = sum(gamma_down + gamma_up):

    Tr((L + c/2)^5) = 60 * 4^N * sum over triangles Im(t_ab t_bc t_ca) * det[(1,1,1), gamma_down, gamma_up],

t_ab = J + i D_z the amplitude of the hop a -> b. The gate works with M = 2L + c and s = gamma_down + gamma_up,
a = gamma_down - gamma_up, where the same statement reads

    Tr(M^5) = -960 * 4^N * sum over ordered pairs l != k of s_l a_k * sum_j Im(t_lk t_kj t_jl).

numpy alone, integer inputs, every comparison exact (all values are Gaussian integers far below 2^53, see (E);
the settings with rational rates are run at 120 or 72 times the rates, the fifth trace being of degree two in
them). A bond may carry ZZ without hopping, and the random rates take either sign, the identities being polynomial.

  (M) the three maps as identities of generators, with Hc the complex conjugate of H:
      rho -> Z^N rho, then the transpose of the generator:  M(H; s, a) -> M(Hc; -s, a)
      the spin flip rho -> X^N rho X^N:                      M(H; s, a) -> M(Hc; s, -a)
      rho -> rho^T:                                          M(H; s, a) -> M(-Hc; s, a)
  (S) Tr M^m(H; -s, -a) = Tr M^m(H; s, a): the degree in the rates is even
  (T) Tr M^m(H; s, a) = Tr M^m(-H; s, -a): the degree in H plus the degree in a is even
  (A) the two second-order traces of the derivation, with k = 2^(N-2), t_xy the amplitude of the hop x -> y,
      P = 4 sum_j t_lj t_jk (the detour l -> j -> k) and tau = 2 t_lk (the direct hop):
      Tr(H^2 s+_l s-_k) = k (P - 2 Delta_lk tau),  Tr(s+_l H s-_k H) = k (P + 2 Delta_lk tau)
  (F) Tr M^5 against the closed form, on random graphs at N = 3, 4, 5
  (R) the values F137 records: -656 sin(Phi) at Phi = pi/2, gate M11's control -3936 (run at Delta = 7, the
      trace not depending on Delta); a fifth trace of zero at collinear rate points with three ratios, on M11's
      triangle with a pendant site at another ratio, and on a ring of four with flux pi/2 and four ratios (not
      M11's own ring, whose rates would leave the exact range of floats; no ring of four has a triangle);
      and a temperature gradient at one coupling, gamma_down = g (nbar + 1), gamma_up = g nbar, a line of rate
      points that misses the origin: third and fifth trace zero, Tr((L + c/2)^7) = 430080 on the unit triangle
      at flux pi/2 with g = 2, nbar = (1, 2, 3)
  (E) every entry of every power is an integer, and the sum of |P_ij| |M_jk| over j, which bounds every partial
      sum of the products, stays below 2^50
  (C) controls that the closed form can fail: a longitudinal field, and a bond that does not conserve the
      excitation number, each move Tr M^5 off it
"""
import itertools
import random
import sys
from fractions import Fraction

import numpy as np

SM = np.array([[0, 1], [0, 0]], dtype=complex)      # sigma^- = |0><1|, bit 1 excited
SP = SM.T.copy()
Z = np.diag([1, -1]).astype(complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)


def site(op, l, n):
    m = np.eye(1, dtype=complex)
    for q in range(n):
        m = np.kron(m, op if q == l else np.eye(2))
    return m


def sup(A, B):
    """rho -> A rho B on the row-major vectorisation."""
    return np.kron(A, B.T)


def t_of(hop, x, y):
    """Amplitude of the hop x -> y."""
    if (x, y) in hop:
        return hop[(x, y)]
    if (y, x) in hop:
        return hop[(y, x)].conjugate()
    return 0


def hamiltonian(n, hop, zz):
    H = np.zeros((1 << n, 1 << n), dtype=complex)
    for (a, b), t in hop.items():
        h = 2 * t * site(SP, b, n) @ site(SM, a, n)
        H += h + h.conj().T
    for (a, b), d in zz.items():
        H += d * site(Z, a, n) @ site(Z, b, n)
    return H


def generator(n, H, s, a):
    """M = 2L + c with loss (s + a)/2 and gain (s - a)/2 per site."""
    I = np.eye(1 << n)
    M = -2j * (sup(H, I) - sup(I, H))
    for l in range(n):
        p, m = site(SP, l, n), site(SM, l, n)
        d, u = (s[l] + a[l]) / 2, (s[l] - a[l]) / 2
        nn = p @ m
        M += 2 * d * (sup(m, p) - (sup(nn, I) + sup(I, nn)) / 2)
        M += 2 * u * (sup(p, m) - (sup(I - nn, I) + sup(I, I - nn)) / 2)
        M += (d + u) * np.eye(M.shape[0])
    return M


def closed_form(n, hop, s, a):
    total = 0
    for l, k in itertools.permutations(range(n), 2):
        loop = sum((t_of(hop, l, k) * t_of(hop, k, j) * t_of(hop, j, l)).imag for j in range(n) if j not in (l, k))
        total += -960 * 4 ** n * s[l] * a[k] * loop
    return total


BIGGEST = [0.0]
EXACT = [True]


def trace_powers(M, top):
    out, P = {}, np.eye(M.shape[0])
    for m in range(1, top + 1):
        if m > 1:
            BIGGEST[0] = max(BIGGEST[0], float((np.abs(P) @ np.abs(M)).max()))
        P = P @ M
        out[m] = np.trace(P)
        BIGGEST[0] = max(BIGGEST[0], float(np.abs(P).sum(axis=1).max()))
        EXACT[0] = EXACT[0] and bool(np.all(P.real == np.round(P.real))) and bool(np.all(P.imag == np.round(P.imag)))
    return out


def maps_hold(n, H, s, a):
    dim = 1 << n
    zn, xn = np.eye(1, dtype=complex), np.eye(1, dtype=complex)
    for _ in range(n):
        zn, xn = np.kron(zn, Z), np.kron(xn, X)
    I = np.eye(dim)
    M = generator(n, H, s, a)
    neg_s, neg_a = [-x for x in s], [-x for x in a]
    K = sup(zn, I)
    first = np.array_equal((K @ M @ K).T, generator(n, H.conj(), neg_s, a))
    F = sup(xn, xn)
    second = np.array_equal(F @ M @ F, generator(n, H.conj(), s, neg_a))
    T = np.zeros((dim * dim, dim * dim))
    for i in range(dim):
        for j in range(dim):
            T[j * dim + i, i * dim + j] = 1
    third = np.array_equal(T @ M @ T, generator(n, -H.conj(), s, a))
    return first, second, third


FAILED = []


def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILED.append(label)


def random_setting(rnd, n):
    hop, zz = {}, {}
    for a, b in itertools.combinations(range(n), 2):
        if rnd.random() < 0.8:
            hop[(a, b)] = complex(rnd.randint(-3, 3), rnd.randint(-3, 3))
        if rnd.random() < 0.8:
            zz[(a, b)] = rnd.randint(-2, 2)
    s = [2 * rnd.randint(1, 4) for _ in range(n)]
    a = [2 * rnd.randint(-3, 3) for _ in range(n)]
    return hop, zz, s, a


def main():
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    nonzero = 0
    for n, trials, top in ((3, 6, 7), (4, 6, 5), (5, 2, 5)):
        for _ in range(trials):
            hop, zz, s, a = random_setting(rnd, n)
            H = hamiltonian(n, hop, zz)
            tr = trace_powers(generator(n, H, s, a), top)
            neg = [-x for x in s], [-x for x in a]
            tr_s = trace_powers(generator(n, H, neg[0], neg[1]), top)
            tr_t = trace_powers(generator(n, -H, s, neg[1]), top)
            sym_s = all(tr[m] == tr_s[m] for m in tr)
            sym_t = all(tr[m] == tr_t[m] for m in tr)
            k2 = 1 << (n - 2)
            lem = True
            for l, k in itertools.permutations(range(n), 2):
                tau = 2 * t_of(hop, l, k)
                P = sum(4 * t_of(hop, l, j) * t_of(hop, j, k) for j in range(n) if j not in (l, k))
                D = zz.get((min(l, k), max(l, k)), 0)
                a20 = np.trace(H @ H @ site(SP, l, n) @ site(SM, k, n))
                a11 = np.trace(site(SP, l, n) @ H @ site(SM, k, n) @ H)
                lem = lem and a20 == k2 * (P - 2 * D * tau) and a11 == k2 * (P + 2 * D * tau)
            want = closed_form(n, hop, s, a)
            nonzero += want != 0
            parts = dict(zip(("M1", "M2", "M3"), maps_hold(n, H, s, a)))
            parts.update(S=sym_s, T=sym_t, A=bool(lem), F=bool(tr[5] == want))
            bad = [name for name, ok in parts.items() if not ok]
            check(f"N={n}, {len(hop)} bonds: (M), (S), (T), (A), (F) Tr M^5 = {tr[5].real:.0f}",
                  not bad, f"closed form {want:.0f}" + (f"; failed: {bad}" if bad else ""))
    check("the random settings are not all silent (at least six with a nonzero fifth trace)", nonzero >= 6, f"{nonzero}")

    # (R) F137's recorded values, rates times 120: gamma_down = (1/2, 1, 3/2), gamma_up = (1/8, 1/3, 1/5)
    d, u = [60, 120, 180], [15, 40, 24]
    s, a = [x + y for x, y in zip(d, u)], [x - y for x, y in zip(d, u)]
    scale = 32 * 120 ** 2       # (2L + c)^5 = 32 (L + c/2)^5, degree two in the rates
    hop = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j}
    tr5 = trace_powers(generator(3, hamiltonian(3, hop, {}), s, a), 5)[5]
    check("F137's triangle at flux pi/2: Tr((L + c/2)^5) = -656", tr5 == -656 * scale and closed_form(3, hop, s, a) == -656 * scale)
    tr5x = trace_powers(generator(3, hamiltonian(3, hop, {}), s, [-x for x in a]), 5)[5]
    check("the same with loss and gain exchanged: +656", tr5x == 656 * scale)
    hop = {(0, 1): 1 + 0j, (1, 2): 2 + 0j, (2, 0): 3j}
    zz = {(0, 1): 7, (1, 2): 7, (0, 2): 7}     # Delta = 7/10 in gate M11; any Delta gives the same trace
    tr5 = trace_powers(generator(3, hamiltonian(3, hop, zz), s, a), 5)[5]
    check("gate M11's control (hoppings 1, 2, 3i, ZZ on every bond): Tr((L + c/2)^5) = -3936",
          tr5 == -3936 * scale and closed_form(3, hop, s, a) == -3936 * scale)
    line = [x // 3 + 20 for x in d]            # gamma_up = gamma_down/3 + 1/6: collinear, three different ratios
    sl, al = [x + y for x, y in zip(d, line)], [x - y for x, y in zip(d, line)]
    check("collinear rate points at three different ratios: the fifth trace is zero",
          trace_powers(generator(3, hamiltonian(3, hop, zz), sl, al), 5)[5] == 0 and len({Fraction(x, y) for x, y in zip(line, d)}) == 3)

    dg, ug = [4, 6, 8], [2, 4, 6]                # g = 2, nbar = (1, 2, 3)
    trg = trace_powers(generator(3, hamiltonian(3, {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j}, {}),
                                 [x + y for x, y in zip(dg, ug)], [x - y for x, y in zip(dg, ug)]), 7)
    check("a temperature gradient at one coupling on the unit triangle at flux pi/2: third and fifth trace zero, "
          "Tr((L + c/2)^7) = 430080", trg[3] == 0 and trg[5] == 0 and trg[7] == 430080 * 128, f"{trg[7].real:.0f}")

    # two broken settings at N = 4, rates times 72: the breaks sit above the fifth trace
    d4, u4 = [36, 72, 108, 72], [9, 18, 27, 8]          # ratio 1/4 on the triangle, 1/9 on the fourth site
    s4, a4 = [x + y for x, y in zip(d4, u4)], [x - y for x, y in zip(d4, u4)]
    hop4 = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j, (3, 0): 1 + 0j}
    check("M11's triangle at one ratio with a fourth site at another, joined by a hopping bond: the fifth trace is zero",
          trace_powers(generator(4, hamiltonian(4, hop4, {(0, 3): 1}), s4, a4), 5)[5] == 0)
    ring = {(0, 1): 1 + 0j, (1, 2): 2 + 0j, (2, 3): 3 + 0j, (3, 0): 17j}
    d4, u4 = [36, 72, 108, 72], [9, 24, 20, 8]          # four ratios
    s4, a4 = [x + y for x, y in zip(d4, u4)], [x - y for x, y in zip(d4, u4)]
    check("a ring of four with flux pi/2, ZZ on every bond and four ratios: the fifth trace is zero",
          trace_powers(generator(4, hamiltonian(4, ring, {(0, 1): 7, (1, 2): 7, (2, 3): 7, (0, 3): 7}), s4, a4), 5)[5] == 0)

    # (C) controls
    Hf = hamiltonian(3, hop, zz) + 3 * site(Z, 0, 3)
    check("control, a longitudinal field on one site: the closed form fails",
          trace_powers(generator(3, Hf, s, a), 5)[5] != closed_form(3, hop, s, a))
    Hx = hamiltonian(3, hop, zz) + 2 * site(X, 0, 3) @ site(X, 1, 3)
    check("control, an XX term on one bond (the excitation number not conserved): the closed form fails",
          trace_powers(generator(3, Hx, s, a), 5)[5] != closed_form(3, hop, s, a))

    check("every entry of the powers is an integer and every partial sum of the products stays below 2^50",
          EXACT[0] and BIGGEST[0] < 2 ** 50, f"largest {BIGGEST[0]:.3g}")

    print("ALL PASS" if not FAILED else f"FAILED: {FAILED}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
