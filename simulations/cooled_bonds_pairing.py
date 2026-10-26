"""Which two-site bonds pair under cooling alone: exact checks with controls.

experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md. Every site carries cooling sigma^- = |0><1| at its own rate and nothing
else; H is a chain of one bond repeated, J_ab·sigma_a(x)sigma_b over a, b in {X, Y} plus Δ·ZZ, or a bond with a mixed
letter (XZ, ZX). Pairing means the spectrum of L is symmetric about −c/2, c = Σ rates. Rows print PASS or FAIL; the
script exits nonzero on any FAIL. With --map it first prints the float survey the document quotes (a reading, not a
gate): random integer couplings over every subset of the nine two-letter bonds, N = 2, 3, 4.

  A  The turn lemma: for every real 2x2 J there is a rotation R with R·Jᵀ·Rᵀ = F·J·F, F = diag(1, −1) (Jᵀ and FJF share
     trace, antisymmetric part and the length of the traceless symmetric part). Checked on six J.
  B  The carrier at rates the chain's reversal P keeps: E = (Ad X^N ∘ Ad(U·P))⁻¹ ∘ S_IZ, S_IZ: rho -> rho·Z^N, U the
     uniform turn about Z realizing R, gives E·L·E⁻¹ = −L† − c exactly (F137's Ad_V∘S_XY with V = (U·P)†); controls:
     rates P does not keep.
  C  The Theodolite's elements (exact, read off the letters and verified): X/Y letters alone at any rates (a staggered
     A·rhoᵀ·B); with ZZ, a real or an imaginary X/Y part (A·rho·B), J symmetric (a = 0) or traceless (s = 0) at
     unequal rates, w = 0 (XXZ with an axial DM term), and any J on one bond
     (turned); silent where s, a and w are all nonzero at rates the reversal
     does not keep. J = s·1 + a·ε + T, ε = [[0, 1], [−1, 0]], T traceless symmetric, w = T₁₁ + i·T₁₂.
  D  Breaks, exact: an odd centred power sum Tr((L − tr L/dim)^k) that is not zero rules out every element; the first
     such k is asserted: 7 where s, a, w are all nonzero (XX + 2XY + ZZ at (2, 1, 1), (1, 2, 3) and, N = 4,
     (1, 2, 3, 4); 2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3)), 5 for a mixed letter (XZ; ZX + XX). Control: zero up to
     k = 7 at rates the reversal keeps.
  E  The seventh sum in closed form: Tr((L + c/2)^7) = −7·2¹⁴·4^(N−3)·Δ·s·a·|w|²·(γ₁γ₂ − γ_{N−1}γ_N), exact as a
     polynomial in the rates at N = 2 (zero, one J), 3 and 4 for three J (interpolated modulo primes, three surplus
     points asserted to fit, the CRT value unchanged by a fifth prime), at random rates at N = 5 and 6 (modulo a
     prime, every entry an integer); controls a = 0 and Δ = 0 give zero; where the seventh is silent off the mirror,
     rates (1, 6, 2, 3) and (1, 0, 2), the ninth is not.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding='utf-8')

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import cube_letters as cl  # noqa: E402

X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Z = sp.diag(1, -1); I2 = sp.eye(2)
LOW = ct.LOWER
LET = {'X': X, 'Y': Y, 'Z': Z}
F = sp.diag(1, -1)
FAILURES = []


def check(name, ok, detail=''):
    print(f"  {'PASS' if ok else 'FAIL'}: {name}" + (f"  [{detail}]" if detail else ''))
    if not ok:
        FAILURES.append(name)


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.kronecker_product(out, m)
    return out


def emb(ops, n):
    return kron(*[ops.get(l, I2) for l in range(n)])


def chain_terms(n, bond):
    """bond: {two-letter string: coefficient}, repeated on every bond of the open chain."""
    return [(s, (b, b + 1), c) for b in range(n - 1) for s, c in bond.items()]


def hamiltonian(n, bond):
    H = sp.zeros(2 ** n, 2 ** n)
    for s, (a, b), c in chain_terms(n, bond):
        H += c * emb({a: LET[s[0]], b: LET[s[1]]}, n)
    return H


def generator(H, rates):
    n = len(rates); d = 2 ** n; Id = sp.eye(d)
    L = -sp.I * (kron(H, Id) - kron(Id, H.T))
    for l, r in enumerate(rates):
        A = emb({l: LOW}, n); AdA = A.H * A
        L += r * (kron(A, A.conjugate()) - sp.Rational(1, 2) * (kron(AdA, Id) + kron(Id, AdA.T)))
    return L


def bond_of(J, zz):
    out = {a + b: J[i, j] for i, a in enumerate('XY') for j, b in enumerate('XY') if J[i, j] != 0}
    if zz:
        out['ZZ'] = zz
    return out


def turn_for(J):
    c, s = sp.symbols('c s', real=True)
    R = sp.Matrix([[c, -s], [s, c]])
    sols = sp.solve(list(R * J.T * R.T - F * J * F) + [c ** 2 + s ** 2 - 1], [c, s], dict=True)
    if not sols:
        return None
    sol = sols[0]
    if s not in sol:                       # a free angle (J commutes with every rotation's action): take s = 0
        sol = {c: sp.sympify(sol.get(c, 1)).subs(s, 0), s: 0}
        if sol[c] ** 2 != 1:
            sol[c] = 1
    return sol[c], sol[s]


def reversal(n):
    P = sp.zeros(2 ** n, 2 ** n)
    for k in range(2 ** n):
        bits = [(k >> (n - 1 - l)) & 1 for l in range(n)]
        P[sum(b << (n - 1 - l) for l, b in enumerate(reversed(bits))), k] = 1
    return P


def carrier_reflects(n, J, zz, rates):
    cth, sth = turn_for(J)
    W = emb({l: sp.diag(1, cth + sp.I * sth) for l in range(n)}, n) * reversal(n)
    Xn, Zn = emb({l: X for l in range(n)}, n), emb({l: Z for l in range(n)}, n)
    Ad = lambda A: kron(A, A.conjugate())                       # rho -> A·rho·A†, row-major
    E = (Ad(Xn) * Ad(W)).inv() * kron(sp.eye(2 ** n), Zn.T)    # S_IZ: rho -> rho·Z^N
    L = generator(hamiltonian(n, bond_of(J, zz)), rates)
    res = E * L + (L.H + sum(rates) * sp.eye(4 ** n)) * E
    return res.applyfunc(lambda e: sp.simplify(sp.expand(e))).is_zero_matrix


def setting(n, bond, rates):
    return ct.Setting(n, chain_terms(n, bond), [((l,), LOW, r) for l, r in enumerate(rates)])


def letters(s, rel):
    rd = cl.read(s, rel)
    if rd is not None and cl.verify(s, rel, rd.maps):
        return 'transposed' if 'ᵀ' in cl.describe(rd)[0] else 'plain'
    return None


def site_turns(s, rel):
    """The Theodolite's left turn x_l per site (site factor [[0, x], [1, 0]]), exact, or None."""
    rd = cl.read(s, rel)
    if rd is None or not cl.verify(s, rel, rd.maps):
        return None
    out = []
    for line in cl.describe(rd)[1:]:
        if line.startswith('site ') and 'A = [[0, x], [1, 0]] with x = ' in line:
            out.append(sp.sympify(line.split('with x = ')[1].split(';')[0], locals={'I': sp.I}))
    return out


def first_odd_trace(n, bond, rates, kmax=7):
    L = generator(hamiltonian(n, bond), rates)
    dim = L.shape[0]
    M = DomainMatrix.from_Matrix((dim * L - L.trace() * sp.eye(dim)).applyfunc(sp.expand)).convert_to(sp.QQ_I)
    P = M
    for k in range(1, kmax + 1):
        if k > 1:
            P = P * M
        if k % 2:
            tr = sum(P.to_Matrix()[i, i] for i in range(dim))
            if tr != 0:
                return k
    return None


def survey():
    """The float reading the document quotes: per (N, rates), how many random bonds pair, by letter class."""
    P = {k: np.array(v, dtype=complex) for k, v in ct.PAULI.items()}
    low = np.array(LOW, dtype=complex)
    letters9 = [a + b for a in 'XYZ' for b in 'XYZ']
    rng = np.random.default_rng(1)
    for n, rates in [(2, [1, 2]), (3, [1, 1, 1]), (3, [1, 2, 1]), (3, [1, 2, 3]), (4, [1, 2, 2, 1]), (4, [1, 2, 3, 5])]:
        d = 2 ** n; Id = np.eye(d); tally = {}; margin = {True: 0.0, False: np.inf}
        for k in range(1, 10):
            for sub in itertools.combinations(letters9, k):
                H = np.zeros((d, d), complex)
                coef = {s: rng.integers(1, 4) * rng.choice([-1, 1]) for s in sub}   # one bond, repeated
                for b in range(n - 1):
                    for s in sub:
                        ops = [np.eye(2)] * n
                        ops[b], ops[b + 1] = P[s[0]], P[s[1]]
                        term = ops[0]
                        for o in ops[1:]:
                            term = np.kron(term, o)
                        H += coef[s] * term
                L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
                for l, r in enumerate(rates):
                    ops = [np.eye(2)] * n; ops[l] = low
                    A = ops[0]
                    for o in ops[1:]:
                        A = np.kron(A, o)
                    AdA = A.conj().T @ A
                    L += r * (np.kron(A, A.conj()) - 0.5 * np.kron(AdA, Id) - 0.5 * np.kron(Id, AdA.T))
                ev = np.linalg.eigvals(L)
                dist = max(np.min(np.abs(ev + sum(rates) + e)) for e in ev)
                mixed = any(set(s) & {'Z'} and s != 'ZZ' for s in sub)
                Jc = {k: coef.get(k, 0) for k in ('XX', 'XY', 'YX', 'YY')}
                s_, a_ = (Jc['XX'] + Jc['YY']) / 2, (Jc['XY'] - Jc['YX']) / 2
                w_ = complex((Jc['XX'] - Jc['YY']) / 2, (Jc['XY'] + Jc['YX']) / 2)
                cls = ('mixed letter' if mixed else 'Δ = 0' if 'ZZ' not in sub
                       else 'Δ ≠ 0, s·a·w = 0' if s_ * a_ * abs(w_) == 0 else 'Δ ≠ 0, s, a, w ≠ 0')
                key = (cls, dist < 1e-6)
                tally[key] = tally.get(key, 0) + 1
                margin[key[1]] = (max if key[1] else min)(margin[key[1]], dist)
        print(f'  N = {n}, rates {rates}:', {f'{c} {"pairs" if p else "breaks"}': v for (c, p), v in sorted(tally.items())},
              f'(largest distance among the pairing {margin[True]:.1e}, smallest among the breaking {margin[False]:.2f})')


if '--map' in sys.argv:
    print("map (float reading, one random bond in ±{1, 2, 3} repeated along the chain, every subset of the nine letters)")
    survey()

print("A  the turn lemma")
for J in [sp.Matrix([[1, 2], [0, 0]]), sp.Matrix([[3, -1], [2, 5]]), sp.Matrix([[1, 0], [0, -4]]),
          sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[2, 3], [-3, 2]]), sp.Matrix([[0, 2], [3, 0]])]:
    t = turn_for(J)
    R = sp.Matrix([[t[0], -t[1]], [t[1], t[0]]]) if t else None
    check(f"J = {J.tolist()}: a rotation R with R·Jᵀ·Rᵀ = F·J·F exists",
          t is not None and sp.simplify(t[0] ** 2 + t[1] ** 2) == 1
          and (R * J.T * R.T - F * J * F).applyfunc(sp.simplify).is_zero_matrix, f'cos, sin = {t}')

print("B  the carrier: chain reversal, uniform turn, S_IZ")
J1, J2 = sp.Matrix([[1, 2], [0, 0]]), sp.Matrix([[3, -1], [2, 5]])
check("XX + 2XY + ZZ, N = 3, rates (1, 2, 1): E·L·E⁻¹ = −L† − c", carrier_reflects(3, J1, 1, [1, 2, 1]))
check("XX + 2XY + ZZ, N = 3, rates (1, 1, 1)", carrier_reflects(3, J1, 1, [1, 1, 1]))
check("J = [[3, −1], [2, 5]] − 2ZZ, N = 3, rates (2, 7, 2)", carrier_reflects(3, J2, -2, [2, 7, 2]))
check("XX + 2XY + ZZ, N = 4, rates (1, 2, 2, 1)", carrier_reflects(4, J1, 1, [1, 2, 2, 1]))
check("controls: rates the reversal does not keep, (2, 1, 1) and (2, 7, 3): it is not",
      not carrier_reflects(3, J1, 1, [2, 1, 1]) and not carrier_reflects(3, J2, -2, [2, 7, 3]))

print("C  the Theodolite's elements")
xy_only = {'XX': 1, 'XY': 2, 'YX': 3, 'YY': 5}
check("X/Y letters alone, XX + 2XY + 3YX + 5YY, rates (1, 2, 3): a transposed sandwich reflects ('refl')",
      letters(setting(3, xy_only, [1, 2, 3]), 'refl') == 'transposed')
check("real part with ZZ, XX + 3YY + ZZ, rates (1, 2, 3): a sandwich reflects",
      letters(setting(3, {'XX': 1, 'YY': 3, 'ZZ': 1}, [1, 2, 3]), 'refl') == 'plain')
check("imaginary part with ZZ, XY + 3YX + ZZ, rates (1, 2, 3): a sandwich reflects",
      letters(setting(3, {'XY': 1, 'YX': 3, 'ZZ': 1}, [1, 2, 3]), 'refl') == 'plain')
xs = site_turns(setting(3, {'XX': 1, 'XY': 2, 'YX': 2, 'ZZ': 1}, [1, 2, 3]), 'refl')
check("J symmetric (a = 0), XX + 2XY + 2YX + ZZ, rates (1, 2, 3): one turn on every site, x = ±(−1 + 4i)/√17",
      xs is not None and len(xs) == 3 and all(sp.simplify(x - xs[0]) == 0 for x in xs)
      and sp.simplify(xs[0] ** 2 - (-1 + 4 * sp.I) ** 2 / 17) == 0, f'x = {xs}')
xs = site_turns(setting(3, {'XY': 1, 'YX': 3, 'ZZ': 1}, [1, 2, 3]), 'refl')
check("imaginary part with ZZ, XY + 3YX + ZZ: the turn alternates by a half turn, x = (−1, 1, −1) up to sign",
      xs is not None and [sp.simplify(x / xs[0]) for x in xs] == [1, -1, 1], f'x = {xs}')
check("J traceless (s = 0), 2XX + 3XY − YX − 2YY + ZZ, rates (1, 2, 3): a turned sandwich reflects",
      letters(setting(3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': -2, 'ZZ': 1}, [1, 2, 3]), 'refl') is not None)
check("w = 0 (XXZ with an axial DM term), XX + YY + XY − YX + ZZ, rates (1, 2, 3): a turned sandwich reflects",
      letters(setting(3, {'XX': 1, 'YY': 1, 'XY': 1, 'YX': -1, 'ZZ': 1}, [1, 2, 3]), 'refl') is not None)
check("s, a, w all nonzero with ZZ on one bond, XX + 2XY + ZZ, N = 2, rates (1, 2): a turned sandwich reflects",
      letters(setting(2, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2]), 'refl') is not None)
check("s, a, w all nonzero with ZZ, N = 3, rates the reversal keeps (1, 2, 1): the letters are silent (no site map "
      "carries the reversal)", letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [1, 2, 1]), 'refl') is None)
check("s, a, w all nonzero with ZZ, N = 3: the letters are silent (XX + 2XY + ZZ at (2, 1, 1), "
      "2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3))",
      letters(setting(3, {'XX': 1, 'XY': 2, 'ZZ': 1}, [2, 1, 1]), 'refl') is None
      and letters(setting(3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3]), 'refl') is None)

print("D  breaks, exact")
mixed = {'XX': 1, 'XY': 2, 'ZZ': 1}
for n, bond, rates, name, first in [
        (3, mixed, [2, 1, 1], 'XX + 2XY + ZZ, rates (2, 1, 1)', 7),
        (3, mixed, [1, 2, 3], 'XX + 2XY + ZZ, rates (1, 2, 3)', 7),
        (4, mixed, [1, 2, 3, 4], 'XX + 2XY + ZZ, N = 4, rates (1, 2, 3, 4)', 7),
        (3, {'XX': 2, 'XY': 3, 'YX': -1, 'YY': 2, 'ZZ': 1}, [1, 2, 3], '2XX + 3XY − YX + 2YY + ZZ, rates (1, 2, 3)', 7),
        (3, {'XZ': 1}, [1, 1, 1], 'XZ, rates (1, 1, 1)', 5),
        (3, {'ZX': 1, 'XX': 1}, [1, 2, 3], 'ZX + XX, rates (1, 2, 3)', 5)]:
    k = first_odd_trace(n, bond, rates)
    check(f"{name}: the first odd centred power sum that is not zero is the {first}th", k == first, f'k = {k}')
check("control: XX + 2XY + ZZ at rates (1, 2, 1) has every odd power sum zero up to k = 7",
      first_odd_trace(3, mixed, [1, 2, 1]) is None)

print("E  the seventh sum in closed form")
PRIMES = [p_ for p_ in sp.primerange(900000, 1048576) if p_ % 4 == 1][-5:]   # float BLAS stays exact below 2^20


def modp_generator(n, bond, rates, p):
    """2L + c over GF(p), i the square root of −1 (p = 1 mod 4); row-major vec, rho -> A rho B is A ⊗ Bᵀ."""
    i = next(r for r in (pow(g, (p - 1) // 4, p) for g in range(2, p)) if r * r % p == p - 1)
    let = {'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, p - i], [i, 0]]), 'Z': np.array([[1, 0], [0, p - 1]])}
    low = np.array([[0, 1], [0, 0]])

    def emb_(ops):
        out = np.array([[1]], dtype=np.int64)
        for l in range(n):
            out = np.kron(out, ops.get(l, np.eye(2, dtype=np.int64))) % p
        return out
    d = 2 ** n
    H = np.zeros((d, d), dtype=np.int64)
    for b in range(n - 1):
        for lab, v in bond.items():
            H = (H + v * emb_({b: let[lab[0]], b + 1: let[lab[1]]})) % p
    Id = np.eye(d, dtype=np.int64)
    M = (-2 * i) * ((np.kron(H, Id) - np.kron(Id, H.T)) % p) % p
    for l, g in enumerate(rates):
        A = emb_({l: low}); AdA = A.T @ A
        M = (M + 2 * g * np.kron(A, A) - g * (np.kron(AdA, Id) + np.kron(Id, AdA.T)) + g * np.eye(d * d, dtype=np.int64)) % p
    return M


def modp_trace_power(M, k, p):
    R = M
    for _ in range(k - 1):
        R = (R.astype(np.float64) @ M.astype(np.float64) % p).astype(np.int64)
    return int(np.trace(R) % p)


def modp_solve(A, b, p):
    A = np.concatenate([np.array(A, dtype=np.int64) % p, np.array(b, dtype=np.int64)[:, None] % p], axis=1)
    n = A.shape[1] - 1
    for c in range(n):
        k = c + int(np.nonzero(A[c:, c])[0][0])
        A[[c, k]] = A[[k, c]]
        A[c] = A[c] * pow(int(A[c, c]), p - 2, p) % p
        f = A[:, c].copy(); f[c] = 0
        A = (A - (f[:, None] * A[c][None, :]) % p) % p
    assert not A[n:].any(), 'the surplus points do not fit the polynomial'
    return [int(x) for x in A[:n, n]]


def seventh_sum_polynomial(n, bond, k=7):
    """Tr((2L + c)^k) as a polynomial in the rates (degree at most k), interpolated modulo each prime; returns the
    CRT polynomial from four primes and whether a fifth prime changes it."""
    gs = sp.symbols(f'g1:{n + 1}')
    monos = sorted(e for e in itertools.product(range(k + 1), repeat=n) if sum(e) <= k)
    rng = np.random.default_rng(1)
    pts = [tuple(int(x) for x in rng.integers(1, 10 ** 5, n)) for _ in range(len(monos) + 3)]
    per_p = []
    for p in PRIMES:
        A = [[int(np.prod([pow(x, e_, p) for x, e_ in zip(pt, e)], dtype=object) % p) for e in monos] for pt in pts]
        per_p.append(modp_solve(A, [modp_trace_power(modp_generator(n, bond, pt, p), k, p) for pt in pts], p))

    def crt(m):
        out = 0
        for j, e in enumerate(monos):
            r, mod = sp.ntheory.modular.crt(PRIMES[:m], [c[j] for c in per_p[:m]])
            r = int(r) - (int(mod) if int(r) > int(mod) // 2 else 0)
            out += r * sp.Mul(*[g ** e_ for g, e_ in zip(gs, e)])
        return sp.expand(out)
    four = crt(4)
    return four, gs, sp.expand(four - crt(5)) == 0


def closed_form(n, bond, gs):
    J = {lab: sp.Integer(bond.get(lab, 0)) for lab in ('XX', 'XY', 'YX', 'YY', 'ZZ')}
    s_, a_ = (J['XX'] + J['YY']) / 2, (J['XY'] - J['YX']) / 2
    w2 = ((J['XX'] - J['YY']) / 2) ** 2 + ((J['XY'] + J['YX']) / 2) ** 2
    return -7 * 2 ** 14 * 4 ** (n - 3) * J['ZZ'] * s_ * a_ * w2 * (gs[0] * gs[1] - gs[-2] * gs[-1]) * 2 ** 7


bonds_E = [({'XX': 1, 'XY': 2, 'ZZ': 1}, 'XX + 2XY + ZZ'), ({'XX': 3, 'XY': -1, 'YX': 2, 'YY': 5, 'ZZ': -2},
           '3XX − XY + 2YX + 5YY − 2ZZ'), ({'XX': -2, 'XY': 1, 'YX': 3, 'YY': 1, 'ZZ': 3}, '−2XX + XY + 3YX + YY + 3ZZ')]
for n in (2, 3, 4):
    for bond, name in (bonds_E if n > 2 else bonds_E[:1]):
        P, gs, stable = seventh_sum_polynomial(n, bond)
        target = 0 if n == 2 else closed_form(n, bond, gs)
        check(f"N = {n}, {name}: the seventh sum is the closed form as a polynomial in the rates"
              + (" (zero)" if n == 2 else ""), stable and sp.expand(P - target) == 0, f'{sp.factor(P)}')
for bond, name in [({'XX': 1, 'XY': 2, 'YX': 2, 'ZZ': 1}, 'a = 0, XX + 2XY + 2YX + ZZ'),
                   ({'XX': 1, 'XY': 2}, 'Δ = 0, XX + 2XY')]:
    P, gs, stable = seventh_sum_polynomial(3, bond)
    check(f"control, N = 3, {name}: the seventh sum is zero at every rate", stable and P == 0)
rng = np.random.default_rng(3)
for n in (5, 6):
    for bond, name in bonds_E[:2]:
        g = [int(x) for x in rng.integers(1, 1000, n)]
        p = PRIMES[0]
        v = closed_form(n, bond, sp.symbols(f'g1:{n + 1}')).subs(dict(zip(sp.symbols(f'g1:{n + 1}'), g)))
        check(f"N = {n}, {name}, rates {tuple(g)}: the seventh sum equals the closed form modulo {p}",
              modp_trace_power(modp_generator(n, bond, g, p), 7, p) == int(v) % p)
M = modp_generator(4, bonds_E[0][0], [1, 6, 2, 3], PRIMES[0])
check("N = 4, XX + 2XY + ZZ, rates (1, 6, 2, 3) (γ₁γ₂ = γ₃γ₄, not mirror-symmetric): the seventh sum is zero and "
      "the ninth is not", modp_trace_power(M, 7, PRIMES[0]) == 0 and modp_trace_power(M, 9, PRIMES[0]) != 0)
M = modp_generator(3, bonds_E[0][0], [1, 0, 2], PRIMES[0])
check("N = 3, XX + 2XY + ZZ, rates (1, 0, 2) (no cooling on the middle site): the seventh sum is zero and the ninth "
      "is not", modp_trace_power(M, 7, PRIMES[0]) == 0 and modp_trace_power(M, 9, PRIMES[0]) != 0)

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILED")
    sys.exit(1)
print("ALL PASS")
